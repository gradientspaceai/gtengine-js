// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ETNonmanifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// ETNonmanifoldMesh is a dynamic edge-triangle mesh that allows an edge to be
// shared by any number of triangles. The mesh stores the edges in a map keyed
// by the unordered edge key <min(v0,v1),max(v0,v1)> and the triangles in a
// map keyed by the ordered triangle key of <v0,v1,v2>.
//
// Port notes: upstream uses std::shared_ptr<Edge>/std::shared_ptr<Triangle>
// for the map values and std::weak_ptr with the WeakPtrLT comparator of
// WeakPtrCompare.h (a C++-runtime file that is intentionally not ported) for
// the back references. The port uses plain object references and null; every
// weak reference that upstream would let expire is explicitly dropped by
// remove() and clear(), so the two behave the same.
//
// The upstream std::map<EdgeKey<false>, ...> and std::map<TriangleKey<true>,
// ...> become JavaScript Maps keyed by FeatureKey.mapKey(). Where upstream
// relies on the ordered iteration of those containers, the port sorts
// explicitly with FeatureKey.compare:
//   - getEdges() and getTriangles() return arrays in std::map order.
//   - assign() copies the triangles in triangle-key order, as upstream
//     operator= does.
//   - getComponents()/getComponentKeys() visit the triangles in triangle-key
//     order and, within depthFirstSearch, visit the triangles adjacent to an
//     edge in triangle-key order; the latter is the order in which the
//     upstream std::set<std::weak_ptr<Triangle>, WeakPtrLT<Triangle>>
//     iterates, because WeakPtrLT forwards to Triangle::operator<, which
//     compares TriangleKey<true> values. Both orderings determine the order
//     of the components and of the triangles within a component.
// The std::map<std::shared_ptr<Triangle>, int32_t> that records the visited
// state is keyed by pointer and is never iterated, so it becomes a Map keyed
// by the Triangle object.
//
// The upstream ECreator/TCreator function pointers, which let derived meshes
// attach their own edge and triangle data, become the optional constructor
// callbacks of the same purpose. The upstream copy constructor and
// operator=, which perform a deep copy by reinserting the triangles, become
// clone() and assign().

import { logAssert } from './Logger.js';
import { EdgeKey } from './EdgeKey.js';
import { TriangleKey } from './TriangleKey.js';
import { FeatureKey } from './FeatureKey.js';

// The port of ETNonmanifoldMesh::Edge.
export class ETNonmanifoldMeshEdge {
    // Vertices of the edge. Note that these are stored in the order in which
    // the edge was first encountered by insert(), so V[0] < V[1] is not
    // guaranteed; the unordered EdgeKey of (V[0],V[1]) is the map key.
    V: [number, number];

    // Triangles sharing the edge, keyed by the ordered triangle key. The
    // port of the upstream std::set of weak pointers.
    T: Map<string, ETNonmanifoldMeshTriangle>;

    constructor(v0: number, v1: number) {
        this.V = [v0, v1];
        this.T = new Map<string, ETNonmanifoldMeshTriangle>();
    }

    // The triangles sharing the edge, in the order in which the upstream
    // std::set of weak pointers iterates (increasing triangle key).
    getTriangles(): ETNonmanifoldMeshTriangle[] {
        return sortTriangles(Array.from(this.T.values()));
    }

    // The port of Edge::operator<, the comparison of the unordered edge keys.
    lessThan(other: ETNonmanifoldMeshEdge): boolean {
        return edgeKeyOf(this).lessThan(edgeKeyOf(other));
    }
}

// The port of ETNonmanifoldMesh::Triangle.
export class ETNonmanifoldMeshTriangle {
    // Vertices listed in counterclockwise order (V[0],V[1],V[2]).
    V: [number, number, number];

    // Adjacent edges. E[i] points to edge (V[i],V[(i+1)%3]).
    E: [ETNonmanifoldMeshEdge | null, ETNonmanifoldMeshEdge | null, ETNonmanifoldMeshEdge | null];

    constructor(v0: number, v1: number, v2: number) {
        this.V = [v0, v1, v2];
        this.E = [null, null, null];
    }

    // The port of Triangle::operator<, the comparison of the ordered
    // triangle keys.
    lessThan(other: ETNonmanifoldMeshTriangle): boolean {
        return triangleKeyOf(this).lessThan(triangleKeyOf(other));
    }
}

export type ETNonmanifoldMeshECreator = (v0: number, v1: number) => ETNonmanifoldMeshEdge;
export type ETNonmanifoldMeshTCreator =
    (v0: number, v1: number, v2: number) => ETNonmanifoldMeshTriangle;

// The unordered key of an edge, EdgeKey<false>(V[0],V[1]).
function edgeKeyOf(edge: ETNonmanifoldMeshEdge): EdgeKey {
    return new EdgeKey(false, edge.V[0], edge.V[1]);
}

// The ordered key of a triangle, TriangleKey<true>(V[0],V[1],V[2]).
function triangleKeyOf(tri: ETNonmanifoldMeshTriangle): TriangleKey {
    return new TriangleKey(true, tri.V[0], tri.V[1], tri.V[2]);
}

// Sort triangles the way the upstream std::set<std::weak_ptr<Triangle>,
// WeakPtrLT<Triangle>> and std::map<TriangleKey<true>, ...> iterate.
function sortTriangles(triangles: ETNonmanifoldMeshTriangle[]): ETNonmanifoldMeshTriangle[] {
    triangles.sort((t0, t1) => FeatureKey.compare(triangleKeyOf(t0), triangleKeyOf(t1)));
    return triangles;
}

export class ETNonmanifoldMesh {
    protected mECreator: ETNonmanifoldMeshECreator;
    protected mEMap: Map<string, ETNonmanifoldMeshEdge>;
    protected mTCreator: ETNonmanifoldMeshTCreator;
    protected mTMap: Map<string, ETNonmanifoldMeshTriangle>;

    constructor(eCreator?: ETNonmanifoldMeshECreator, tCreator?: ETNonmanifoldMeshTCreator) {
        this.mECreator = eCreator ?? ETNonmanifoldMesh.createEdge;
        this.mEMap = new Map<string, ETNonmanifoldMeshEdge>();
        this.mTCreator = tCreator ?? ETNonmanifoldMesh.createTriangle;
        this.mTMap = new Map<string, ETNonmanifoldMeshTriangle>();
    }

    // Support for a deep copy of the mesh. The edge and triangle objects are
    // not shared between the meshes; the copy reinserts the triangles of the
    // source mesh, which is the upstream operator=. Note that the triangles
    // are reinserted using the vertices of the triangle keys rather than the
    // vertices of the triangle objects, so a copied triangle may have its
    // vertices cyclically rotated relative to the original (the winding
    // order is preserved). This is the upstream behavior.
    assign(mesh: ETNonmanifoldMesh): this {
        this.clear();

        this.mECreator = mesh.mECreator;
        this.mTCreator = mesh.mTCreator;
        for (const tkey of mesh.getTriangleKeys()) {
            this.insert(tkey.V[0], tkey.V[1], tkey.V[2]);
        }

        return this;
    }

    // The port of the upstream copy constructor.
    clone(): ETNonmanifoldMesh {
        return new ETNonmanifoldMesh().assign(this);
    }

    // The edges in increasing order of the unordered edge key, the order in
    // which the upstream std::map iterates.
    getEdges(): ETNonmanifoldMeshEdge[] {
        const edges = Array.from(this.mEMap.values());
        edges.sort((e0, e1) => FeatureKey.compare(edgeKeyOf(e0), edgeKeyOf(e1)));
        return edges;
    }

    // The triangles in increasing order of the ordered triangle key, the
    // order in which the upstream std::map iterates.
    getTriangles(): ETNonmanifoldMeshTriangle[] {
        return sortTriangles(Array.from(this.mTMap.values()));
    }

    // The keys of the edge map, in the order the upstream std::map iterates.
    getEdgeKeys(): EdgeKey[] {
        return this.getEdges().map(edge => edgeKeyOf(edge));
    }

    // The keys of the triangle map, in the order the upstream std::map
    // iterates.
    getTriangleKeys(): TriangleKey[] {
        return this.getTriangles().map(tri => triangleKeyOf(tri));
    }

    // The port of the upstream std::map find operations.
    getEdge(v0: number, v1: number): ETNonmanifoldMeshEdge | null {
        return this.mEMap.get(new EdgeKey(false, v0, v1).mapKey()) ?? null;
    }

    getTriangle(v0: number, v1: number, v2: number): ETNonmanifoldMeshTriangle | null {
        return this.mTMap.get(new TriangleKey(true, v0, v1, v2).mapKey()) ?? null;
    }

    getNumEdges(): number {
        return this.mEMap.size;
    }

    getNumTriangles(): number {
        return this.mTMap.size;
    }

    // If <v0,v1,v2> is not in the mesh, a Triangle object is created and
    // returned; otherwise, <v0,v1,v2> is in the mesh and null is returned.
    insert(v0: number, v1: number, v2: number): ETNonmanifoldMeshTriangle | null {
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
            } else {
                // The edge was previously encountered and created.
                logAssert(edge !== null, 'Unexpected condition.');
            }

            // Associate the edge with the triangle.
            tri.E[i0] = edge;

            // Update the adjacent set of triangles for the edge.
            edge.T.set(tkey, tri);
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
            const edge = tri.E[i];
            logAssert(edge !== null, 'Unexpected condition.');

            // Remove the triangle from the edge's set of adjacent triangles.
            const numRemoved = edge.T.delete(tkey) ? 1 : 0;
            logAssert(numRemoved > 0, 'Unexpected condition.');

            // Remove the edge if you have the last reference to it.
            if (edge.T.size === 0) {
                this.mEMap.delete(edgeKeyOf(edge).mapKey());
            }
        }

        // Remove the triangle from the graph.
        this.mTMap.delete(tkey);
        return true;
    }

    // Destroy the edges and triangles to obtain an empty mesh.
    clear(): void {
        this.mEMap.clear();
        this.mTMap.clear();
    }

    // A manifold mesh has the property that an edge is shared by at most two
    // triangles.
    isManifold(): boolean {
        for (const edge of this.mEMap.values()) {
            if (edge.T.size > 2) {
                return false;
            }
        }
        return true;
    }

    // A manifold mesh is closed if each edge is shared twice. A closed mesh
    // is not necessarily oriented. For example, you could have a mesh with
    // spherical topology. The upper hemisphere has outer-facing normals and
    // the lower hemisphere has inner-facing normals. The discontinuity in
    // orientation occurs on the circle shared by the hemispheres.
    isClosed(): boolean {
        for (const edge of this.mEMap.values()) {
            if (edge.T.size !== 2) {
                return false;
            }
        }
        return true;
    }

    // Compute the connected components of the edge-triangle graph that the
    // mesh represents. The first function returns references to the triangle
    // objects of 'this' mesh, so you must consume the components before
    // clearing 'this'. The second function returns triangle keys, which
    // allows you to clear 'this' before consuming the components.
    getComponents(): ETNonmanifoldMeshTriangle[][] {
        const components: ETNonmanifoldMeshTriangle[][] = [];

        // visited: 0 (unvisited), 1 (discovered), 2 (finished)
        const visited = new Map<ETNonmanifoldMeshTriangle, number>();
        for (const tri of this.mTMap.values()) {
            visited.set(tri, 0);
        }

        for (const tri of this.getTriangles()) {
            if (visited.get(tri) === 0) {
                const component: ETNonmanifoldMeshTriangle[] = [];
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
        const visited = new Map<ETNonmanifoldMeshTriangle, number>();
        for (const tri of this.mTMap.values()) {
            visited.set(tri, 0);
        }

        for (const tri of this.getTriangles()) {
            if (visited.get(tri) === 0) {
                const component: ETNonmanifoldMeshTriangle[] = [];
                this.depthFirstSearch(tri, visited, component);
                components.push(component.map(t => triangleKeyOf(t)));
            }
        }

        return components;
    }

    // The edge data and default edge creation.
    protected static createEdge(v0: number, v1: number): ETNonmanifoldMeshEdge {
        return new ETNonmanifoldMeshEdge(v0, v1);
    }

    // The triangle data and default triangle creation.
    protected static createTriangle(v0: number, v1: number, v2: number): ETNonmanifoldMeshTriangle {
        return new ETNonmanifoldMeshTriangle(v0, v1, v2);
    }

    // Support for computing connected components. This is a straightforward
    // depth-first search of the graph but uses a preallocated stack rather
    // than a recursive function that could possibly overflow the call stack.
    protected depthFirstSearch(tInitial: ETNonmanifoldMeshTriangle,
        visited: Map<ETNonmanifoldMeshTriangle, number>,
        component: ETNonmanifoldMeshTriangle[]): void {
        // Allocate the maximum-size stack that can occur in the depth-first
        // search. The stack is empty when the index top is -1.
        const tStack = new Array<ETNonmanifoldMeshTriangle | null>(this.mTMap.size).fill(null);
        let top = -1;
        tStack[++top] = tInitial;
        while (top >= 0) {
            const tri = tStack[top];
            logAssert(tri !== null, 'Unexpected condition.');
            visited.set(tri, 1);
            let i: number;
            for (i = 0; i < 3; ++i) {
                const edge = tri.E[i];
                logAssert(edge !== null, 'Unexpected condition.');

                let foundUnvisited = false;
                for (const adj of edge.getTriangles()) {
                    if (visited.get(adj) === 0) {
                        tStack[++top] = adj;
                        foundUnvisited = true;
                        break;
                    }
                }

                if (foundUnvisited) {
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
}
