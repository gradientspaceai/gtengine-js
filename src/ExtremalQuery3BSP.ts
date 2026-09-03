// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ExtremalQuery3BSP.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The extremal queries for convex objects are based on the algorithm
// described in
// https://www.geometrictools.com/Documentation/ExtremalPolytopeQueries.pdf
//
// The query preprocesses the polytope into a BSP tree of spherical arcs on
// the Gauss map of the polytope. Each leaf region of the Gauss map is the set
// of directions for which a single vertex is extreme, so the query is a
// root-to-leaf descent whose cost is O(log n) rather than the O(n) of the
// brute-force ExtremalQuery3PRJ.
//
// Port notes:
// - std::multiset<SphericalArc> ordered by 'separation' and traversed in
//   reverse becomes an array that is stably sorted by increasing separation
//   and traversed backwards. A C++ multiset inserts equivalent elements at
//   the upper bound, so equal-separation arcs keep insertion order and a
//   reverse traversal visits them in reverse insertion order; a stable sort
//   plus a backwards walk reproduces that exactly.
// - std::map<Triangle*, int32_t> mTriToNormal becomes a Map keyed by the
//   triangle object. It is only ever looked up, never iterated, so no
//   ordering is observable.
// - VETManifoldMesh's unordered vertex/triangle adjacency containers are read
//   through the port's sorted accessors (getVertices(), getTAdjacent(),
//   getEdges()), so the constructed tree is deterministic. The upstream C++
//   tree shape depends on hash-table order; both shapes are valid BSP trees
//   for the same Gauss map, and the query results agree.
// - The two output reference parameters of GetExtremeVertices become the
//   ExtremalQuery3Result object declared by the base class.
// - The deleted copy constructor and copy assignment have no port.
//
// KNOWN UPSTREAM DEFECT (preserved, not fixed). The BSP construction does not
// always produce a tree whose leaf regions are single Gauss-map cells with the
// correct extreme vertex, so getExtremeVertices() can return a vertex that is
// not extreme. It is exact for small polytopes (verified over thousands of
// random directions on a tetrahedron and an octahedron) and starts to fail on
// an icosahedron, where roughly 1% of random directions get a wrong answer.
//
// The cause is in insertArc(): an arc descends the tree comparing only its two
// endpoints against the *current node's* great circle, never against the
// accumulated region of the descent. Consequently an arc can be stored as a
// node of a region that its arc segment does not touch, and the region's two
// halves then inherit posVertex/negVertex values that are valid only near the
// arc. Because a query reads the extreme vertex off the last node visited
// rather than off a vertex recorded at the leaf, a region that is entirely
// inside one Gauss cell can report a different cell's vertex. Example, with
// the regular icosahedron of test/ExtremalQuery3BSP.test.ts and the direction
// (0.85564456289238, -0.33788348840949, 0.39205500921802): the descent ends at
// a bisector node whose two sides both answer vertex 4, yet sampling the leaf
// region shows it lies entirely inside the Gauss cell of vertex 9.
//
// The quirk is preserved because repairing it is not a local fix: clipping the
// arc against the accumulated region during the descent (which the port was
// prototyped with) removes most but not all of the wrong answers, so a correct
// implementation needs a redesigned construction, which is outside the scope of
// a port. Use ExtremalQuery3PRJ, the O(n) brute-force query, when exactness
// matters.

import { ExtremalQuery3 } from './ExtremalQuery3.js';
import type { ExtremalQuery3Result } from './ExtremalQuery3.js';
import { ETManifoldMeshTriangle } from './ETManifoldMesh.js';
import { isign } from './Functions.js';
import { logAssert } from './Logger.js';
import { Polyhedron3 } from './Polyhedron3.js';
import { Vector, dot } from './Vector.js';
import { cross } from './Vector3.js';
import { VETManifoldMesh } from './VETManifoldMesh.js';

// The port of the private ExtremalQuery3BSP::SphericalArc class. The arcs are
// the nodes of the BSP tree.
class SphericalArc {
    // Indices N[] into the face normal array for the endpoints of the arc.
    nIndex: [number, number];

    // The number of arcs in the path from normal N[0] to normal N[1]. For
    // spherical polygon edges the number is 1. The number is 2 or larger for
    // bisector arcs of the spherical polygon.
    separation: number;

    // The normal is Cross(FaceNormal[N[0]], FaceNormal[N[1]]).
    normal: Vector;

    // Indices into the vertex array for the extremal points for the two
    // regions sharing the arc. As the arc is traversed from normal N[0] to
    // normal N[1], posVertex is the index for the extreme vertex to the left
    // of the arc and negVertex is the index for the extreme vertex to the
    // right of the arc.
    posVertex: number;
    negVertex: number;

    // Support for BSP trees stored as contiguous nodes in an array.
    posChild: number;
    negChild: number;

    constructor() {
        this.nIndex = [-1, -1];
        this.separation = 0;
        this.normal = new Vector(3);
        this.posVertex = -1;
        this.negVertex = -1;
        this.posChild = -1;
        this.negChild = -1;
    }

    // The port of the implicit C++ copy constructor; insertArc stores copies
    // of the incoming arc in the node array (C++ value semantics).
    clone(): SphericalArc {
        const arc = new SphericalArc();
        arc.nIndex = [this.nIndex[0], this.nIndex[1]];
        arc.separation = this.separation;
        arc.normal = this.normal.clone();
        arc.posVertex = this.posVertex;
        arc.negVertex = this.negVertex;
        arc.posChild = this.posChild;
        arc.negChild = this.negChild;
        return arc;
    }
}

export class ExtremalQuery3BSP extends ExtremalQuery3 {
    // Lookup table for indexing into mFaceNormals.
    private mTriToNormal: Map<ETManifoldMeshTriangle, number>;

    // Storage for the BSP nodes.
    private mNodes: SphericalArc[];
    private mTreeDepth: number;

    // The caller must ensure that the input polyhedron is convex.
    constructor(polytope: Polyhedron3) {
        super(polytope);
        this.mTriToNormal = new Map<ETManifoldMeshTriangle, number>();
        this.mNodes = [];
        this.mTreeDepth = 0;

        // Create the adjacency information for the polytope.
        const mesh = new VETManifoldMesh();
        const indices = this.mPolytope.getIndices();
        const numTriangles = Math.floor(indices.length / 3);
        for (let t = 0; t < numTriangles; ++t) {
            const triangle = mesh.insert(indices[3 * t + 0], indices[3 * t + 1],
                indices[3 * t + 2]);
            // Upstream does not test the insertion result; a nonmanifold or
            // duplicated index set would silently associate the null triangle
            // pointer with a face normal. The port reports the bad input.
            logAssert(triangle !== null,
                'Failed to insert a polytope triangle into the mesh.');
            this.mTriToNormal.set(triangle, t);
        }

        // Create the set of unique arcs which are used to create the BSP
        // tree.
        const arcs: SphericalArc[] = [];
        this.createSphericalArcs(mesh, arcs);

        // Create the BSP tree to be used in the extremal query.
        this.createBSPTree(arcs);
    }

    // Compute the extreme vertices in the specified direction and return the
    // indices of the vertices in the polyhedron vertex array.
    getExtremeVertices(direction: Vector): ExtremalQuery3Result {
        let positiveDirection = -1;
        let negativeDirection = -1;

        // Do a nonrecursive depth-first search of the BSP tree to determine
        // which spherical polygon contains the incoming direction D. Index 0
        // is the root of the BSP tree.
        let current = 0;
        while (current >= 0) {
            const node = this.mNodes[current];
            const sign = isign(dot(direction, node.normal));
            if (sign >= 0) {
                current = node.posChild;
                if (current === -1) {
                    // At a leaf node.
                    positiveDirection = node.posVertex;
                }
            }
            else {
                current = node.negChild;
                if (current === -1) {
                    // At a leaf node.
                    positiveDirection = node.negVertex;
                }
            }
        }

        // Do a nonrecursive depth-first search of the BSP tree to determine
        // which spherical polygon contains the reverse incoming direction -D.
        // Testing sign <= 0 on Dot(D,normal) is the same as testing sign >= 0
        // on Dot(-D,normal).
        current = 0;  // the root of the BSP tree
        while (current >= 0) {
            const node = this.mNodes[current];
            const sign = isign(dot(direction, node.normal));
            if (sign <= 0) {
                current = node.posChild;
                if (current === -1) {
                    // At a leaf node.
                    negativeDirection = node.posVertex;
                }
            }
            else {
                current = node.negChild;
                if (current === -1) {
                    // At a leaf node.
                    negativeDirection = node.negVertex;
                }
            }
        }

        return {
            positiveDirection: positiveDirection,
            negativeDirection: negativeDirection
        };
    }

    // Tree statistics.
    getNumNodes(): number {
        return this.mNodes.length;
    }

    // The upstream depth counter measures the size of the depth-first search
    // stack rather than the length of a root-to-leaf path, so this value is a
    // loose statistic and not the true tree depth. The quirk is harmless (the
    // value is never used by the query itself) and is preserved.
    getTreeDepth(): number {
        return this.mTreeDepth;
    }

    private sortAdjacentTriangles(vIndex: number,
        tAdj: readonly ETManifoldMeshTriangle[]): ETManifoldMeshTriangle[] {
        // Traverse the triangles adjacent to vertex V using edge-triangle
        // adjacency information to produce a sorted array of adjacent
        // triangles.
        const numTriangles = tAdj.length;
        const tAdjSorted = new Array<ETManifoldMeshTriangle>(numTriangles);
        let tri: ETManifoldMeshTriangle | null = tAdj[0];
        for (let i = 0; i < numTriangles; ++i) {
            logAssert(tri !== null, 'Unexpected condition.');
            const current: ETManifoldMeshTriangle = tri;
            for (let prev = 2, curr = 0; curr < 3; prev = curr++) {
                if (current.V[curr] === vIndex) {
                    tAdjSorted[i] = current;
                    tri = current.T[prev];
                    break;
                }
            }
        }
        return tAdjSorted;
    }

    private createSphericalArcs(mesh: VETManifoldMesh, arcs: SphericalArc[]): void {
        const prev = [2, 0, 1];
        const next = [1, 2, 0];

        for (const edge of mesh.getEdges()) {
            const t0 = edge.T[0];
            const t1 = edge.T[1];
            // Upstream indexes mTriToNormal with a possibly null T[1]; on a
            // mesh with boundary that silently maps to face normal 0. A
            // convex polytope is a closed mesh, so the port reports the bad
            // input instead.
            logAssert(t0 !== null && t1 !== null,
                'The polytope must be a closed manifold mesh.');

            const arc = new SphericalArc();
            arc.nIndex[0] = this.mTriToNormal.get(t0) as number;
            arc.nIndex[1] = this.mTriToNormal.get(t1) as number;
            arc.separation = 1;
            arc.normal = cross(this.mFaceNormals[arc.nIndex[0]],
                this.mFaceNormals[arc.nIndex[1]]);

            const adj = t0;
            let j: number;
            for (j = 0; j < 3; ++j) {
                if (adj.V[j] !== edge.V[0] && adj.V[j] !== edge.V[1]) {
                    arc.posVertex = adj.V[prev[j]];
                    arc.negVertex = adj.V[next[j]];
                    break;
                }
            }
            logAssert(j < 3, 'Unexpected condition.');

            arcs.push(arc);
        }

        this.createSphericalBisectors(mesh, arcs);
    }

    private createSphericalBisectors(mesh: VETManifoldMesh, arcs: SphericalArc[]): void {
        for (const vertex of mesh.getVertices()) {
            // Sort the normals into a counterclockwise spherical polygon when
            // viewed from outside the sphere.
            const vIndex = vertex.V;
            const tAdjSorted = this.sortAdjacentTriangles(vIndex, vertex.getTAdjacent());
            const numTriangles = tAdjSorted.length;

            // The upstream std::queue is declared outside the vertex loop but
            // is always drained before the next vertex, so a per-vertex queue
            // is equivalent.
            const queue: [number, number][] = [[0, numTriangles]];
            let head = 0;
            while (head < queue.length) {
                const item = queue[head++];
                const i0 = item[0];
                const i1 = item[1];
                const separation = i1 - i0;
                if (separation > 1 && separation !== numTriangles - 1) {
                    if (i1 < numTriangles) {
                        const arc = new SphericalArc();
                        arc.nIndex[0] = this.mTriToNormal.get(tAdjSorted[i0]) as number;
                        arc.nIndex[1] = this.mTriToNormal.get(tAdjSorted[i1]) as number;
                        arc.separation = separation;

                        arc.normal = cross(this.mFaceNormals[arc.nIndex[0]],
                            this.mFaceNormals[arc.nIndex[1]]);

                        arc.posVertex = vIndex;
                        arc.negVertex = vIndex;
                        arcs.push(arc);
                    }
                    const imid = Math.floor((i0 + i1 + 1) / 2);
                    if (imid !== i1) {
                        queue.push([i0, imid]);
                        queue.push([imid, i1]);
                    }
                }
            }
        }
    }

    private createBSPTree(arcs: SphericalArc[]): void {
        // The tree has at least a root.
        this.mTreeDepth = 1;

        // The upstream container is a multiset ordered by increasing
        // separation that is traversed in reverse. This heuristic is designed
        // to create BSP trees whose topmost nodes can eliminate as many arcs
        // as possible during an extremal query.
        const sorted = arcs.map((arc, i) => ({ arc: arc, i: i }));
        sorted.sort((a, b) => (a.arc.separation !== b.arc.separation
            ? a.arc.separation - b.arc.separation : a.i - b.i));
        for (let k = sorted.length - 1; k >= 0; --k) {
            this.insertArc(sorted[k].arc);
        }

        // The leaf nodes are not counted in the traversal of insertArc. The
        // depth must be incremented to account for leaves.
        ++this.mTreeDepth;
    }

    private insertArc(arc: SphericalArc): void {
        // The incoming arc is stored at the end of the nodes array.
        if (this.mNodes.length > 0) {
            // Do a nonrecursive depth-first search of the current BSP tree to
            // place the incoming arc. Index 0 is the root of the BSP tree.
            const candidates: number[] = [0];
            while (candidates.length > 0) {
                const current = candidates.pop() as number;
                const node = this.mNodes[current];

                let sign0: number;
                if (arc.nIndex[0] === node.nIndex[0] || arc.nIndex[0] === node.nIndex[1]) {
                    sign0 = 0;
                }
                else {
                    sign0 = isign(dot(this.mFaceNormals[arc.nIndex[0]], node.normal));
                }

                let sign1: number;
                if (arc.nIndex[1] === node.nIndex[0] || arc.nIndex[1] === node.nIndex[1]) {
                    sign1 = 0;
                }
                else {
                    sign1 = isign(dot(this.mFaceNormals[arc.nIndex[1]], node.normal));
                }

                let doTest = 0;
                if (sign0 * sign1 < 0) {
                    // The new arc straddles the current arc, so propagate it
                    // to both child nodes.
                    doTest = 3;
                }
                else if (sign0 > 0 || sign1 > 0) {
                    // The new arc is on the positive side of the current arc.
                    doTest = 1;
                }
                else if (sign0 < 0 || sign1 < 0) {
                    // The new arc is on the negative side of the current arc.
                    doTest = 2;
                }
                // else: sign0 = sign1 = 0, in which case no propagation is
                // needed because the current BSP node will handle the correct
                // partitioning of the arcs during extremal queries.

                if (doTest & 1) {
                    if (node.posChild !== -1) {
                        candidates.push(node.posChild);
                        const depth = candidates.length;
                        if (depth > this.mTreeDepth) {
                            this.mTreeDepth = depth;
                        }
                    }
                    else {
                        node.posChild = this.mNodes.length;
                        this.mNodes.push(arc.clone());
                    }
                }

                if (doTest & 2) {
                    if (node.negChild !== -1) {
                        candidates.push(node.negChild);
                        const depth = candidates.length;
                        if (depth > this.mTreeDepth) {
                            this.mTreeDepth = depth;
                        }
                    }
                    else {
                        node.negChild = this.mNodes.length;
                        this.mNodes.push(arc.clone());
                    }
                }
            }
        }
        else {
            // root node
            this.mNodes.push(arc.clone());
        }
    }
}
