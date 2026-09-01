// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ETManifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The ETManifoldMesh class represents an edge-triangle manifold mesh. The 'E'
// stands for edge and the 'T' stands for triangle. It is general purpose,
// allowing insertion and removal of triangles at any time. However, the
// performance is limited because of the use of the container classes. If your
// application requires an edge-triangle manifold mesh for which no triangles
// will be removed, a better choice is StaticVETManifoldMesh.
//
// Port notes:
//   - Upstream stores the edges in an std::unordered_map keyed by
//     EdgeKey<false> and the triangles in an std::unordered_map keyed by
//     TriangleKey<true>. A JavaScript Map compares object keys by identity,
//     so the port keys its Maps by FeatureKey.mapKey() (the ETNonmanifoldMesh
//     and BoxManager precedent). Because the upstream containers are
//     *unordered*, the iteration order of the mesh contents is
//     implementation-defined in C++; the port instead iterates in increasing
//     feature-key order (FeatureKey.compare) everywhere the order is
//     observable, which is deterministic and is what the ETNonmanifoldMesh
//     port does for the ordered std::map containers. The one genuinely
//     ordered upstream container, the std::map of boundary edges in
//     GetBoundaryPolygons, is likewise iterated in increasing key order.
//   - std::unique_ptr<Edge>/std::unique_ptr<Triangle> map values and the raw
//     Edge*/Triangle* back references become plain object references and
//     null. Every pointer that upstream nulls out is nulled here too, so the
//     observable behavior matches.
//   - The ECreator/TCreator function pointers, which let derived meshes
//     attach their own edge and triangle data, become optional constructor
//     callbacks (the VEManifoldMesh/ETNonmanifoldMesh precedent). The copy
//     constructor and operator=, which deep copy by reinserting the
//     triangles, become clone() and assign().
//   - std::numeric_limits<size_t>::max(), the "no adjacent triangle" sentinel
//     of CreateCompactGraph, becomes -1. The upstream comments describe the
//     sentinel as -1 already; JavaScript numbers are signed, so -1 is the
//     natural representation.
//   - The GetOppositeVertexOfEdge output reference parameter becomes a
//     returned { found, uOpposite } object.
//
// Upstream quirk, preserved: when insert() rejects a triangle because the
// mesh would become nonmanifold and throwOnNonmanifoldInsertion(false) has
// been set, upstream returns nullptr but leaves behind the state built by the
// loop iterations that already ran: edges created by this call stay in the
// edge map with T[0] referencing a triangle that was never added to the
// triangle map, and an already-present adjacent triangle may have had its
// T[j] pointed at that same phantom triangle. This is the ETManifoldMesh form
// of the VEManifoldMesh::Insert leak (upstream issue #73), and the port
// reproduces it rather than rolling the insertion back, because it matches
// the sibling port and callers of the graceful path may rely on the upstream
// state. Prefer the default throwing behavior.

import { logAssert, logError } from './Logger';
import { EdgeKey } from './EdgeKey';
import { TriangleKey } from './TriangleKey';
import { FeatureKey } from './FeatureKey';

// The port of ETManifoldMesh::Edge.
export class ETManifoldMeshEdge {
    // Vertices of the edge. These are stored in the order in which the edge
    // was first encountered by insert(), so V[0] < V[1] is not guaranteed;
    // the unordered EdgeKey of (V[0],V[1]) is the map key.
    V: [number, number];

    // Triangles sharing the edge.
    T: [ETManifoldMeshTriangle | null, ETManifoldMeshTriangle | null];

    constructor(v0: number, v1: number) {
        this.V = [v0, v1];
        this.T = [null, null];
    }
}

// The port of ETManifoldMesh::Triangle.
export class ETManifoldMeshTriangle {
    // Vertices, listed in counterclockwise order (V[0],V[1],V[2]).
    V: [number, number, number];

    // Adjacent edges. E[i] points to edge (V[i],V[(i+1)%3]).
    E: [ETManifoldMeshEdge | null, ETManifoldMeshEdge | null, ETManifoldMeshEdge | null];

    // Adjacent triangles. T[i] points to the adjacent triangle sharing edge
    // E[i].
    T: [ETManifoldMeshTriangle | null, ETManifoldMeshTriangle | null,
        ETManifoldMeshTriangle | null];

    constructor(v0: number, v1: number, v2: number) {
        this.V = [v0, v1, v2];
        this.E = [null, null, null];
        this.T = [null, null, null];
    }

    // The edge <u0,u1> is directed. Determine whether the triangle has an
    // edge <V[i],V[(i+1)%3]> = <u0,u1> (return +1) or an edge
    // <V[i],V[(i+1)%3]> = <u1,u0> (return -1) or does not have an edge
    // meeting either condition (return 0).
    whichSideOfEdge(u0: number, u1: number): number {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            if (this.V[i0] === u0 && this.V[i1] === u1) {
                return +1;
            }
            if (this.V[i0] === u1 && this.V[i1] === u0) {
                return -1;
            }
        }
        return 0;
    }

    getAdjacentOfEdge(u0: number, u1: number): ETManifoldMeshTriangle | null {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            if ((this.V[i0] === u0 && this.V[i1] === u1) ||
                (this.V[i0] === u1 && this.V[i1] === u0)) {
                return this.T[i0];
            }
        }
        return null;
    }

    // The port of GetOppositeVertexOfEdge; the uOpposite output reference
    // parameter is returned in the result object. When found is false,
    // uOpposite is the value the upstream caller's variable would keep,
    // which is unspecified; the port reports -1.
    getOppositeVertexOfEdge(u0: number, u1: number): { found: boolean; uOpposite: number } {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            if ((this.V[i0] === u0 && this.V[i1] === u1) ||
                (this.V[i0] === u1 && this.V[i1] === u0)) {
                return { found: true, uOpposite: this.V[(i1 + 1) % 3] };
            }
        }
        return { found: false, uOpposite: -1 };
    }
}

export type ETManifoldMeshECreator = (v0: number, v1: number) => ETManifoldMeshEdge;
export type ETManifoldMeshTCreator =
    (v0: number, v1: number, v2: number) => ETManifoldMeshTriangle;

// The port of the ETManifoldMesh::BoundaryEdge helper struct. The 'directed'
// member is a port addition; it stores the (v0,v1) key of the std::map entry
// so that the map can be iterated in increasing key order.
export class ETManifoldMeshBoundaryEdge {
    triangle: ETManifoldMeshTriangle | null;
    index: number;
    visited: boolean;
    directed: [number, number];

    constructor(directed: [number, number], triangle: ETManifoldMeshTriangle | null = null,
        index = -1, visited = false) {
        this.triangle = triangle;
        this.index = index;
        this.visited = visited;
        this.directed = directed;
    }
}

// The unordered key of an edge, EdgeKey<false>(V[0],V[1]).
function edgeKeyOf(edge: ETManifoldMeshEdge): EdgeKey {
    return new EdgeKey(false, edge.V[0], edge.V[1]);
}

// The ordered key of a triangle, TriangleKey<true>(V[0],V[1],V[2]).
function triangleKeyOf(tri: ETManifoldMeshTriangle): TriangleKey {
    return new TriangleKey(true, tri.V[0], tri.V[1], tri.V[2]);
}

// The key of the boundary-edge map, the directed pair <v0,v1>.
function boundaryKey(v0: number, v1: number): string {
    return `${v0},${v1}`;
}

// The port of std::map::operator[] on the boundary-edge map: a missing key is
// inserted with a default-constructed BoundaryEdge.
function boundaryEdgeAt(boundaryEdges: Map<string, ETManifoldMeshBoundaryEdge>,
    v0: number, v1: number): ETManifoldMeshBoundaryEdge {
    const key = boundaryKey(v0, v1);
    let edge = boundaryEdges.get(key);
    if (edge === undefined) {
        edge = new ETManifoldMeshBoundaryEdge([v0, v1]);
        boundaryEdges.set(key, edge);
    }
    return edge;
}

export class ETManifoldMesh {
    protected mECreator: ETManifoldMeshECreator;
    protected mEMap: Map<string, ETManifoldMeshEdge>;
    protected mTCreator: ETManifoldMeshTCreator;
    protected mTMap: Map<string, ETManifoldMeshTriangle>;
    protected mThrowOnNonmanifoldInsertion: boolean;  // default: true

    constructor(eCreator?: ETManifoldMeshECreator, tCreator?: ETManifoldMeshTCreator) {
        this.mECreator = eCreator ?? ETManifoldMesh.createEdge;
        this.mEMap = new Map<string, ETManifoldMeshEdge>();
        this.mTCreator = tCreator ?? ETManifoldMesh.createTriangle;
        this.mTMap = new Map<string, ETManifoldMeshTriangle>();
        this.mThrowOnNonmanifoldInsertion = true;
    }

    // Support for a deep copy of the mesh; the port of operator=. The edge
    // and triangle objects are not shared between the meshes. Note that the
    // triangles are reinserted using the vertices of the triangle keys rather
    // than the vertices of the triangle objects, so a copied triangle may
    // have its vertices cyclically rotated relative to the original (the
    // winding order is preserved). This is the upstream behavior.
    assign(mesh: ETManifoldMesh): this {
        this.clear();

        this.mECreator = mesh.mECreator;
        this.mTCreator = mesh.mTCreator;
        this.mThrowOnNonmanifoldInsertion = mesh.mThrowOnNonmanifoldInsertion;
        for (const tkey of mesh.getTriangleKeys()) {
            this.insert(tkey.V[0], tkey.V[1], tkey.V[2]);
        }

        return this;
    }

    // The port of the upstream copy constructor.
    clone(): ETManifoldMesh {
        return new ETManifoldMesh().assign(this);
    }

    // Member access. Upstream returns the containers themselves; the port
    // returns arrays of the values in increasing feature-key order.
    getEdges(): ETManifoldMeshEdge[] {
        const edges = Array.from(this.mEMap.values());
        edges.sort((e0, e1) => FeatureKey.compare(edgeKeyOf(e0), edgeKeyOf(e1)));
        return edges;
    }

    getTriangles(): ETManifoldMeshTriangle[] {
        const triangles = Array.from(this.mTMap.values());
        triangles.sort((t0, t1) => FeatureKey.compare(triangleKeyOf(t0), triangleKeyOf(t1)));
        return triangles;
    }

    getEdgeKeys(): EdgeKey[] {
        return this.getEdges().map(edge => edgeKeyOf(edge));
    }

    getTriangleKeys(): TriangleKey[] {
        return this.getTriangles().map(tri => triangleKeyOf(tri));
    }

    // The port of the upstream map find operations.
    getEdge(v0: number, v1: number): ETManifoldMeshEdge | null {
        return this.mEMap.get(new EdgeKey(false, v0, v1).mapKey()) ?? null;
    }

    getTriangle(v0: number, v1: number, v2: number): ETManifoldMeshTriangle | null {
        return this.mTMap.get(new TriangleKey(true, v0, v1, v2).mapKey()) ?? null;
    }

    getNumEdges(): number {
        return this.mEMap.size;
    }

    getNumTriangles(): number {
        return this.mTMap.size;
    }

    // If the insertion of a triangle fails because the mesh would become
    // nonmanifold, the default behavior is to throw an exception. You can
    // disable this behavior and continue gracefully without an exception.
    // The return value is the previous value of the internal state.
    throwOnNonmanifoldInsertion(doException: boolean): boolean {
        const previous = this.mThrowOnNonmanifoldInsertion;
        this.mThrowOnNonmanifoldInsertion = doException;
        return previous;
    }

    // If <v0,v1,v2> is not in the mesh, a Triangle object is created and
    // returned; otherwise, <v0,v1,v2> is in the mesh and null is returned. If
    // the insertion leads to a nonmanifold mesh, the call fails with null
    // returned (or an exception, see throwOnNonmanifoldInsertion).
    insert(v0: number, v1: number, v2: number): ETManifoldMeshTriangle | null {
        const tkey = new TriangleKey(true, v0, v1, v2).mapKey();
        if (this.mTMap.has(tkey)) {
            // The triangle already exists. Return null as a signal to the
            // caller that the insertion failed.
            return null;
        }

        // Create the new triangle. It will be added to mTMap at the end of
        // the function so that if an assertion is triggered and the function
        // returns early, the (bad) triangle will not be part of the mesh.
        const tri = this.mTCreator(v0, v1, v2);

        // Add the edges to the mesh if they do not already exist.
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            const ekey = new EdgeKey(false, tri.V[i0], tri.V[i1]).mapKey();
            let edge = this.mEMap.get(ekey);
            if (edge === undefined) {
                // This is the first time the edge is encountered.
                edge = this.mECreator(tri.V[i0], tri.V[i1]);
                this.mEMap.set(ekey, edge);

                // Update the edge and triangle.
                edge.T[0] = tri;
                tri.E[i0] = edge;
            } else {
                // This is the second time the edge is encountered.
                logAssert(edge !== null, 'Unexpected condition.');

                if (this.mThrowOnNonmanifoldInsertion) {
                    // tri and edge.T[0] must have a shared edge
                    // (tri.V[i0],tri.V[i1]). For tri, the directed edge is
                    // <tri.V[i0],tri.V[i1]>. For edge.T[0], the directed edge
                    // must be <tri.V[i1],tri.V[i0]>.
                    const first = edge.T[0];
                    logAssert(first !== null, 'Unexpected condition.');
                    for (let j = 0; j < 3; ++j) {
                        if (first.V[j] === tri.V[i0]) {
                            logAssert(first.V[(j + 2) % 3] === tri.V[i1],
                                'Attempt to create nonmanifold mesh.');
                        }
                    }
                }

                // Update the edge.
                if (edge.T[1]) {
                    if (this.mThrowOnNonmanifoldInsertion) {
                        logError('Attempt to create nonmanifold mesh.');
                    } else {
                        return null;
                    }
                }
                edge.T[1] = tri;

                // Update the adjacent triangles.
                const adjacent = edge.T[0];
                logAssert(adjacent !== null, 'Unexpected condition.');
                for (let j = 0; j < 3; ++j) {
                    if (adjacent.E[j] === edge) {
                        adjacent.T[j] = tri;
                        break;
                    }
                }

                // Update the triangle.
                tri.E[i0] = edge;
                tri.T[i0] = adjacent;
            }
        }

        this.mTMap.set(tkey, tri);
        return tri;
    }

    // If <v0,v1,v2> is in the mesh, it is removed and 'true' is returned;
    // otherwise, <v0,v1,v2> is not in the mesh and 'false' is returned.
    remove(v0: number, v1: number, v2: number): boolean {
        const tkey = new TriangleKey(true, v0, v1, v2).mapKey();
        const tri = this.mTMap.get(tkey);
        if (tri === undefined) {
            // The triangle does not exist.
            return false;
        }

        // Remove the edges and update adjacent triangles if necessary.
        for (let i = 0; i < 3; ++i) {
            // Inform the edges the triangle is being deleted.
            const edge: ETManifoldMeshEdge | null = tri.E[i];
            logAssert(edge !== null, 'Unexpected condition.');

            if (edge.T[0] === tri) {
                // One-triangle edges always have the reference at index zero.
                edge.T[0] = edge.T[1];
                edge.T[1] = null;
            } else if (edge.T[1] === tri) {
                edge.T[1] = null;
            } else {
                logError('Unexpected condition.');
            }

            // Remove the edge if you have the last reference to it.
            if (!edge.T[0] && !edge.T[1]) {
                this.mEMap.delete(edgeKeyOf(edge).mapKey());
            }

            // Inform adjacent triangles the triangle is being deleted.
            const adjacent: ETManifoldMeshTriangle | null = tri.T[i];
            if (adjacent) {
                for (let j = 0; j < 3; ++j) {
                    if (adjacent.T[j] === tri) {
                        adjacent.T[j] = null;
                        break;
                    }
                }
            }
        }

        this.mTMap.delete(tkey);
        return true;
    }

    // Destroy the edges and triangles to obtain an empty mesh.
    clear(): void {
        this.mEMap.clear();
        this.mTMap.clear();
    }

    // A manifold mesh is closed if each edge is shared twice. A closed mesh
    // is not necessarily oriented. For example, you could have a mesh with
    // spherical topology. The upper hemisphere has outer-facing normals and
    // the lower hemisphere has inner-facing normals. The discontinuity in
    // orientation occurs on the circle shared by the hemispheres.
    isClosed(): boolean {
        for (const edge of this.mEMap.values()) {
            if (!edge.T[0] || !edge.T[1]) {
                return false;
            }
        }
        return true;
    }

    // Test whether all triangles in the mesh are oriented consistently and
    // that no two triangles are coincident. The latter means that you cannot
    // have both triangles <v0,v1,v2> and <v0,v2,v1> in the mesh to be
    // considered oriented.
    isOriented(): boolean {
        for (const edge of this.mEMap.values()) {
            const ekey = edgeKeyOf(edge);
            const t0 = edge.T[0];
            const t1 = edge.T[1];
            if (t0 && t1) {
                // In each triangle, find the ordered edge that corresponds to
                // the unordered edge key. Also find the vertex opposite that
                // edge.
                const edgePositive: [boolean, boolean] = [false, false];
                const vOpposite: [number, number] = [-1, -1];
                const adjacent: [ETManifoldMeshTriangle, ETManifoldMeshTriangle] = [t0, t1];
                for (let j = 0; j < 2; ++j) {
                    const tri = adjacent[j];
                    for (let i = 0; i < 3; ++i) {
                        if (tri.V[i] === ekey.V[0]) {
                            const vNext = tri.V[(i + 1) % 3];
                            if (vNext === ekey.V[1]) {
                                edgePositive[j] = true;
                                vOpposite[j] = tri.V[(i + 2) % 3];
                            } else {
                                edgePositive[j] = false;
                                vOpposite[j] = vNext;
                            }
                            break;
                        }
                    }
                }

                // To be oriented consistently, the edges must have reversed
                // ordering and the opposite vertices cannot match.
                if (edgePositive[0] === edgePositive[1] || vOpposite[0] === vOpposite[1]) {
                    return false;
                }
            }
        }
        return true;
    }

    // Compute the connected components of the edge-triangle graph that the
    // mesh represents. The first function returns references to the triangle
    // objects of 'this' mesh, so you must consume the components before
    // clearing 'this'. The second function returns triangle keys, which
    // allows you to clear 'this' before consuming the components.
    getComponents(): ETManifoldMeshTriangle[][] {
        const components: ETManifoldMeshTriangle[][] = [];

        // visited: 0 (unvisited), 1 (discovered), 2 (finished)
        const visited = new Map<ETManifoldMeshTriangle, number>();
        for (const tri of this.mTMap.values()) {
            visited.set(tri, 0);
        }

        for (const tri of this.getTriangles()) {
            if (visited.get(tri) === 0) {
                const component: ETManifoldMeshTriangle[] = [];
                this.depthFirstSearch(tri, visited, component);
                components.push(component);
            }
        }

        return components;
    }

    // The port of the GetComponents overload that returns triangle keys.
    getComponentKeys(): TriangleKey[][] {
        const components: TriangleKey[][] = [];

        // visited: 0 (unvisited), 1 (discovered), 2 (finished)
        const visited = new Map<ETManifoldMeshTriangle, number>();
        for (const tri of this.mTMap.values()) {
            visited.set(tri, 0);
        }

        for (const tri of this.getTriangles()) {
            if (visited.get(tri) === 0) {
                const component: ETManifoldMeshTriangle[] = [];
                this.depthFirstSearch(tri, visited, component);
                components.push(component.map(t => triangleKeyOf(t)));
            }
        }

        return components;
    }

    // Create a compact edge-triangle graph. The vertex indices are those
    // integers passed to an insert(...) call. These have no meaning to the
    // semantics of maintaining an edge-triangle manifold mesh, so this class
    // makes no assumption about them. The vertex indices do not necessarily
    // start at 0 and they are not necessarily contiguous numbers. The
    // triangles are represented by triples of vertex indices. The compact
    // graph stores these in an array of N triples, say,
    //   T[0] = (v0,v1,v2), T[1] = (v3,v4,v5), ...
    // Each triangle has up to 3 adjacent triangles. The compact graph stores
    // the adjacency information in an array of N triples, say,
    //   A[0] = (t0,t1,t2), A[1] = (t3,t4,t5), ...
    // where the ti are indices into the array of triangles. For example, the
    // triangle T[0] has edges E[0] = (v0,v1), E[1] = (v1,v2) and
    // E[2] = (v2,v0). The edge E[0] has adjacent triangle T[0]. If E[0] has
    // another adjacent triangle, it is T[A[0][0]]. If it does not have another
    // adjacent triangle, then A[0][0] = -1. Similar assignments are made for
    // the other two edges, which produces A[0][1] for E[1] and A[0][2] for
    // E[2].
    createCompactGraph(): { triangles: number[][]; adjacents: number[][] } {
        const numTriangles = this.mTMap.size;
        logAssert(numTriangles > 0, 'Invalid input.');

        const triangles: number[][] = new Array<number[]>(numTriangles);
        const adjacents: number[][] = new Array<number[]>(numTriangles);

        // The triangles are enumerated in increasing triangle-key order; see
        // the port notes about the unordered upstream container.
        const ordered = this.getTriangles();
        const triIndexMap = new Map<ETManifoldMeshTriangle | null, number>();
        triIndexMap.set(null, -1);
        let index = 0;
        for (const tri of ordered) {
            triIndexMap.set(tri, index++);
        }

        index = 0;
        for (const triPtr of ordered) {
            const tri: number[] = [0, 0, 0];
            const adj: number[] = [0, 0, 0];
            for (let j = 0; j < 3; ++j) {
                tri[j] = triPtr.V[j];
                adj[j] = triIndexMap.get(triPtr.T[j]) as number;
            }
            triangles[index] = tri;
            adjacents[index] = adj;
            ++index;
        }

        return { triangles, adjacents };
    }

    // The output of createCompactGraph can be used to compute the connected
    // components of the graph, each component having triangles with the same
    // chirality (winding order). Using only the mesh topology, it is not
    // possible to ensure that the chirality is the same for all the
    // components. Additional application-specific geometric information is
    // required.
    //
    // The returned 'components' contains indices into the 'triangles' array
    // and is partitioned into C subarrays, each representing a connected
    // component. The lengths of the subarrays are stored in
    // 'numComponentTriangles'. The number of elements of that array is C. It
    // is the case that the number of triangles in the mesh is
    // sum_{i=0}^{C-1} numComponentTriangles[i].
    //
    // On return, 'triangles' and 'adjacents' have been modified in place and
    // have the correct chirality.
    static getComponentsConsistentChirality(triangles: number[][], adjacents: number[][]):
        { components: number[]; numComponentTriangles: number[] } {
        logAssert(triangles.length > 0 && triangles.length === adjacents.length,
            'Invalid inputs.');

        // Use a breadth-first search to process the chirality of the
        // triangles. Keep track of the connected components.
        const numTriangles = triangles.length;
        const visited = new Array<boolean>(numTriangles).fill(false);
        const components: number[] = [];
        const numComponentTriangles: number[] = [];

        // The 'firstUnvisited' index is that of the first triangle to process
        // in a connected component of the mesh.
        for (; ;) {
            // Let n[i] be the number of elements of the i-th connected
            // component. Let C be the number of components. During the
            // execution of this loop, the array numComponentTriangles stores
            //   {0, n[0], n[0]+n[1], ..., n[0]+...+n[C-1]=numTriangles}
            // At the end of this function, the array is modified to
            //   {n[0], n[1], ..., n[C-1]}
            numComponentTriangles.push(components.length);

            // Find the starting index of a connected component.
            let firstUnvisited = numTriangles;
            for (let i = 0; i < numTriangles; ++i) {
                if (!visited[i]) {
                    firstUnvisited = i;
                    break;
                }
            }
            if (firstUnvisited === numTriangles) {
                // All connected components have been found.
                break;
            }

            // Initialize the queue to start at the first unvisited triangle
            // of a connected component.
            const triQueue: number[] = [firstUnvisited];
            let queueFront = 0;
            visited[firstUnvisited] = true;
            components.push(firstUnvisited);

            // Perform the breadth-first search.
            while (queueFront < triQueue.length) {
                const curIndex = triQueue[queueFront++];

                const curTriangle = triangles[curIndex];
                for (let i0 = 0; i0 < 3; ++i0) {
                    const adjIndex = adjacents[curIndex][i0];
                    if (adjIndex !== -1 && !visited[adjIndex]) {
                        // The current triangle has a directed edge
                        // <curTriangle[i0],curTriangle[i1]> and there is a
                        // triangle adjacent to it.
                        const i1 = (i0 + 1) % 3;
                        const adjTriangle = triangles[adjIndex];
                        const tv0 = curTriangle[i0];
                        const tv1 = curTriangle[i1];

                        // To have the same chirality, it is required that the
                        // adjacent triangle have the directed edge
                        // <curTriangle[i1],curTriangle[i0]>.
                        let sameChirality = true;
                        let j0 = 0;
                        let j1 = 0;
                        for (j0 = 0; j0 < 3; ++j0) {
                            j1 = (j0 + 1) % 3;
                            if (adjTriangle[j0] === tv0) {
                                if (adjTriangle[j1] === tv1) {
                                    // The adjacent triangle has the same
                                    // directed edge as the current triangle,
                                    // so the chiralities do not match.
                                    sameChirality = false;
                                }
                                break;
                            }
                        }
                        logAssert(j0 < 3, 'Unexpected condition.');

                        if (!sameChirality) {
                            // Swap the vertices of the adjacent triangle that
                            // form the shared directed edge of the current
                            // triangle. This requires that the adjacency
                            // information for the other two edges of the
                            // adjacent triangle be swapped.
                            const adjAdjacent = adjacents[adjIndex];
                            const j2 = (j1 + 1) % 3;
                            const tmpV = adjTriangle[j0];
                            adjTriangle[j0] = adjTriangle[j1];
                            adjTriangle[j1] = tmpV;
                            const tmpA = adjAdjacent[j1];
                            adjAdjacent[j1] = adjAdjacent[j2];
                            adjAdjacent[j2] = tmpA;
                        }

                        // The adjacent triangle has been processed, but it
                        // might have neighbors that need to be processed. Push
                        // the adjacent triangle into the queue to ensure this
                        // happens. Insert the adjacent triangle into the
                        // active connected component.
                        triQueue.push(adjIndex);
                        visited[adjIndex] = true;
                        components.push(adjIndex);
                    }
                }
            }
        }

        // Read the comments at the beginning of this function.
        const numSizes = numComponentTriangles.length;
        logAssert(numSizes > 1, 'Unexpected condition.');
        for (let i0 = 0, i1 = 1; i1 < numComponentTriangles.length; i0 = i1++) {
            numComponentTriangles[i0] = numComponentTriangles[i1] - numComponentTriangles[i0];
        }
        numComponentTriangles.length = numSizes - 1;

        return { components, numComponentTriangles };
    }

    // This is a simple wrapper around createCompactGraph(...) and
    // getComponentsConsistentChirality(...), in particular when you do not
    // need to work directly with the connected components. The mesh is
    // reconstructed, because the bookkeeping details of trying to modify the
    // mesh in-place are horrendous. NOTE: If your mesh has more than 1
    // connected component, you should read the comments for
    // getComponentsConsistentChirality(...) about the potential for different
    // chiralities between components.
    makeConsistentChirality(): void {
        const { triangles, adjacents } = this.createCompactGraph();

        ETManifoldMesh.getComponentsConsistentChirality(triangles, adjacents);

        // Only the 'triangles' array is needed to reconstruct the mesh. The
        // other arrays are discarded.
        this.clear();

        for (const triangle of triangles) {
            this.insert(triangle[0], triangle[1], triangle[2]);
        }
    }

    // Compute the boundary-edge components of the mesh. These are polygons
    // that are simple for the strict definition of manifold mesh that
    // disallows bow-tie configurations. The GTE mesh implementations do allow
    // bow-tie configurations, in which case some polygons might not be simple.
    // If you select duplicateEndpoints to be false, a component has
    // consecutive vertices (v[0], v[1], ..., v[n-1]) and the polygon has edges
    // (v[0],v[1]), (v[1],v[2]), ..., (v[n-2],v[n-1]), (v[n-1],v[0]). If
    // duplicateEndpoints is set to true, a component has consecutive vertices
    // (v[0], v[1], ..., v[n-1], v[0]), emphasizing that the component is
    // closed.
    getBoundaryPolygons(duplicateEndpoints: boolean): number[][] {
        const polygons: number[][] = [];

        // Get the boundary edges. The index into the Triangle.T[] adjacency
        // array is also stored to help with the traversal of polygons.
        const boundaryEdges = new Map<string, ETManifoldMeshBoundaryEdge>();
        for (const tri of this.getTriangles()) {
            for (let i = 0; i < 3; ++i) {
                if (tri.T[i] === null) {
                    const directed: [number, number] = [tri.V[i], tri.V[(i + 1) % 3]];
                    const key = boundaryKey(directed[0], directed[1]);
                    if (!boundaryEdges.has(key)) {
                        // std::map::insert does not overwrite an existing
                        // entry.
                        boundaryEdges.set(key,
                            new ETManifoldMeshBoundaryEdge(directed, tri, i, false));
                    }
                }
            }
        }

        // The upstream container is an std::map, so it is visited in
        // increasing order of the directed pair.
        const initialEdges = Array.from(boundaryEdges.values());
        initialEdges.sort((e0, e1) =>
            (e0.directed[0] - e1.directed[0]) || (e0.directed[1] - e1.directed[1]));

        for (const initialEdge of initialEdges) {
            if (!initialEdge.visited) {
                logAssert(initialEdge.triangle !== null, 'Unexpected condition.');
                const polygon: number[] = [];
                this.getBoundaryPolygon(initialEdge.triangle, initialEdge.index,
                    boundaryEdges, polygon);
                polygons.push(polygon);
            }
        }

        if (!duplicateEndpoints) {
            for (const polygon of polygons) {
                polygon.length = polygon.length - 1;
            }
        }

        return polygons;
    }

    // The edge data and default edge creation.
    protected static createEdge(v0: number, v1: number): ETManifoldMeshEdge {
        return new ETManifoldMeshEdge(v0, v1);
    }

    // The triangle data and default triangle creation.
    protected static createTriangle(v0: number, v1: number, v2: number): ETManifoldMeshTriangle {
        return new ETManifoldMeshTriangle(v0, v1, v2);
    }

    // Support for computing connected components. This is a straightforward
    // depth-first search of the graph but uses a preallocated stack rather
    // than a recursive function that could possibly overflow the call stack.
    protected depthFirstSearch(tInitial: ETManifoldMeshTriangle,
        visited: Map<ETManifoldMeshTriangle, number>,
        component: ETManifoldMeshTriangle[]): void {
        // Allocate the maximum-size stack that can occur in the depth-first
        // search. The stack is empty when the index top is -1.
        const tStack = new Array<ETManifoldMeshTriangle | null>(this.mTMap.size).fill(null);
        let top = -1;
        tStack[++top] = tInitial;
        while (top >= 0) {
            const tri = tStack[top];
            logAssert(tri !== null, 'Unexpected condition.');
            visited.set(tri, 1);
            let i: number;
            for (i = 0; i < 3; ++i) {
                const adj = tri.T[i];
                if (adj && visited.get(adj) === 0) {
                    tStack[++top] = adj;
                    break;
                }
            }
            if (i === 3) {
                visited.set(tri, 2);
                component.push(tri);
                --top;
            }
        }
    }

    protected getBoundaryPolygon(initialTriangle: ETManifoldMeshTriangle,
        initialIndex: number, boundaryEdges: Map<string, ETManifoldMeshBoundaryEdge>,
        polygon: number[]): void {
        let tri: ETManifoldMeshTriangle = initialTriangle;
        let i0 = initialIndex;
        let i1 = (i0 + 1) % 3;
        const vEdge: [number, number] = [tri.V[i0], tri.V[i1]];
        polygon.push(vEdge[0]);
        while (!boundaryEdgeAt(boundaryEdges, vEdge[0], vEdge[1]).visited) {
            polygon.push(vEdge[1]);
            boundaryEdgeAt(boundaryEdges, vEdge[0], vEdge[1]).visited = true;

            // Traverse the triangle strip with vertex at vEdge[1] until the
            // last triangle is encountered. The final edge of the last
            // triangle is the next boundary edge and starts at vEdge[1].
            const visited = new Set<ETManifoldMeshTriangle>();
            visited.add(tri);
            while (tri.T[i1] !== null) {
                tri = tri.T[i1] as ETManifoldMeshTriangle;

                // If this assertion is triggered, try calling isOriented()
                // before calling getBoundaryPolygons(). If isOriented()
                // returns false, the call to getBoundaryPolygons() will fail.
                logAssert(!visited.has(tri),
                    'Triangle already visited. Is the mesh orientable?');
                visited.add(tri);

                let j: number;
                for (j = 0; j < 3; ++j) {
                    if (vEdge[1] === tri.V[j]) {
                        i1 = j;
                        break;
                    }
                }
                logAssert(j < 3, 'Unexpected condition.');
            }

            const i2 = (i1 + 1) % 3;
            vEdge[0] = vEdge[1];
            vEdge[1] = tri.V[i2];
            i0 = i1;
            i1 = i2;
        }
    }
}
