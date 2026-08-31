// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VEManifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// VEManifoldMesh is a dynamic vertex-edge manifold mesh: a graph of directed
// edges in which each vertex is shared by at most two edges. The mesh is a
// union of open and closed polylines.
//
// Port notes: upstream stores the vertices in a std::map keyed by the vertex
// index and the edges in a std::map keyed by the directed pair (v0,v1) of the
// edge; note that the key is not sorted, so <v0,v1> and <v1,v0> are distinct
// edges. The port uses JavaScript Maps with the pair key encoded as a string
// and exposes the contents through getVertices() and getEdges(), which return
// arrays in the same sorted order that the upstream std::map iteration
// produces (vertices by index, edges lexicographically by (V[0],V[1])). The
// upstream std::shared_ptr/std::weak_ptr become plain references and null:
// every weak reference that upstream would let expire is explicitly cleared
// by remove(), so the two behave the same. The upstream VCreator/ECreator
// function pointers, which let derived meshes attach their own vertex and
// edge data, become the optional constructor callbacks of the same purpose.
//
// Upstream quirk, preserved: insert() adds the new edge to the edge map
// before it checks the vertices for a nonmanifold configuration. When
// throwOnNonmanifoldInsertion(false) has been set and the insertion fails,
// upstream returns nullptr but leaves the rejected edge in the edge map (and
// leaves any vertex already updated by the partial loop pointing at it). The
// port reproduces that behavior rather than rolling the insertion back,
// because callers of the graceful path may rely on the upstream state.

import { logAssert, logError } from './Logger';

// The port of VEManifoldMesh::Vertex.
export class VEManifoldMeshVertex {
    // The unique vertex index.
    V: number;

    // The edges (if any) sharing the vertex.
    E: [VEManifoldMeshEdge | null, VEManifoldMeshEdge | null];

    constructor(v: number) {
        this.V = v;
        this.E = [null, null];
    }
}

// The port of VEManifoldMesh::Edge.
export class VEManifoldMeshEdge {
    // Vertices, listed as a directed edge <V[0],V[1]>.
    V: [number, number];

    // Adjacent edges. E[i] points to the edge sharing V[i].
    E: [VEManifoldMeshEdge | null, VEManifoldMeshEdge | null];

    constructor(v0: number, v1: number) {
        this.V = [v0, v1];
        this.E = [null, null];
    }
}

export type VEManifoldMeshVCreator = (v: number) => VEManifoldMeshVertex;
export type VEManifoldMeshECreator = (v0: number, v1: number) => VEManifoldMeshEdge;

export class VEManifoldMesh {
    protected mVCreator: VEManifoldMeshVCreator;
    protected mVMap: Map<number, VEManifoldMeshVertex>;
    protected mECreator: VEManifoldMeshECreator;
    protected mEMap: Map<string, VEManifoldMeshEdge>;
    protected mThrowOnNonmanifoldInsertion: boolean;  // default: true

    constructor(vCreator?: VEManifoldMeshVCreator, eCreator?: VEManifoldMeshECreator) {
        this.mVCreator = vCreator ?? VEManifoldMesh.createVertex;
        this.mVMap = new Map<number, VEManifoldMeshVertex>();
        this.mECreator = eCreator ?? VEManifoldMesh.createEdge;
        this.mEMap = new Map<string, VEManifoldMeshEdge>();
        this.mThrowOnNonmanifoldInsertion = true;
    }

    // The vertices in increasing order of the vertex index, the order in
    // which the upstream std::map iterates.
    getVertices(): VEManifoldMeshVertex[] {
        const vertices = Array.from(this.mVMap.values());
        vertices.sort((v0, v1) => v0.V - v1.V);
        return vertices;
    }

    // The edges in lexicographically increasing order of the directed pair
    // (V[0],V[1]), the order in which the upstream std::map iterates.
    getEdges(): VEManifoldMeshEdge[] {
        const edges = Array.from(this.mEMap.values());
        edges.sort((e0, e1) =>
            (e0.V[0] !== e1.V[0] ? e0.V[0] - e1.V[0] : e0.V[1] - e1.V[1]));
        return edges;
    }

    // The port of the upstream std::map find operations.
    getVertex(v: number): VEManifoldMeshVertex | null {
        return this.mVMap.get(v) ?? null;
    }

    getEdge(v0: number, v1: number): VEManifoldMeshEdge | null {
        return this.mEMap.get(VEManifoldMesh.edgeKey(v0, v1)) ?? null;
    }

    getNumVertices(): number {
        return this.mVMap.size;
    }

    getNumEdges(): number {
        return this.mEMap.size;
    }

    // If the insertion of an edge fails because the mesh would become
    // nonmanifold, the default behavior is to throw an exception. You can
    // disable this behavior and continue gracefully without an exception.
    throwOnNonmanifoldInsertion(doException: boolean): void {
        this.mThrowOnNonmanifoldInsertion = doException;
    }

    // If <v0,v1> is not in the mesh, an Edge object is created and returned;
    // otherwise, <v0,v1> is in the mesh and null is returned. If the
    // insertion leads to a nonmanifold mesh, the call fails with null
    // returned.
    insert(v0: number, v1: number): VEManifoldMeshEdge | null {
        const ekey = VEManifoldMesh.edgeKey(v0, v1);
        if (this.mEMap.has(ekey)) {
            // The edge already exists. Return null as a signal to the caller
            // that the insertion failed.
            return null;
        }

        // Add the new edge.
        const edge = this.mECreator(v0, v1);
        this.mEMap.set(ekey, edge);

        // Add the vertices if they do not already exist.
        for (let i = 0; i < 2; ++i) {
            const v = edge.V[i];
            let vertex = this.mVMap.get(v);
            if (vertex === undefined) {
                // This is the first time the vertex is encountered.
                vertex = this.mVCreator(v);
                this.mVMap.set(v, vertex);

                // Update the vertex.
                vertex.E[0] = edge;
            } else {
                // This is the second time the vertex is encountered.

                // Update the vertex.
                if (vertex.E[1]) {
                    if (this.mThrowOnNonmanifoldInsertion) {
                        logError('The mesh must be manifold.');
                    } else {
                        return null;
                    }
                }
                vertex.E[1] = edge;

                // Update the adjacent edge.
                const adjacent = vertex.E[0];
                logAssert(adjacent !== null, 'Unexpected condition.');
                for (let j = 0; j < 2; ++j) {
                    if (adjacent.V[j] === v) {
                        adjacent.E[j] = edge;
                        break;
                    }
                }

                // Update the edge.
                edge.E[i] = adjacent;
            }
        }

        return edge;
    }

    // If <v0,v1> is in the mesh, it is removed and 'true' is returned;
    // otherwise, <v0,v1> is not in the mesh and 'false' is returned.
    remove(v0: number, v1: number): boolean {
        const ekey = VEManifoldMesh.edgeKey(v0, v1);
        const edge = this.mEMap.get(ekey);
        if (edge === undefined) {
            // The edge does not exist.
            return false;
        }

        // Remove the vertices if necessary (when they are not shared).
        for (let i = 0; i < 2; ++i) {
            // Inform the vertices the edge is being deleted.
            const vertex = this.mVMap.get(edge.V[i]);
            logAssert(vertex !== undefined, 'Unexpected condition.');

            if (vertex.E[0] === edge) {
                // One-edge vertices always have the edge at index zero.
                vertex.E[0] = vertex.E[1];
                vertex.E[1] = null;
            } else if (vertex.E[1] === edge) {
                vertex.E[1] = null;
            } else {
                logError('Unexpected condition.');
            }

            // Remove the vertex if you have the last reference to it.
            if (!vertex.E[0] && !vertex.E[1]) {
                this.mVMap.delete(vertex.V);
            }

            // Inform adjacent edges the edge is being deleted.
            const adjacent: VEManifoldMeshEdge | null = edge.E[i];
            if (adjacent) {
                for (let j = 0; j < 2; ++j) {
                    if (adjacent.E[j] === edge) {
                        adjacent.E[j] = null;
                        break;
                    }
                }
            }
        }

        this.mEMap.delete(ekey);
        return true;
    }

    // A manifold mesh is closed if each vertex is shared twice.
    isClosed(): boolean {
        for (const vertex of this.mVMap.values()) {
            if (!vertex.E[0] || !vertex.E[1]) {
                return false;
            }
        }
        return true;
    }

    // The key of the upstream std::pair<int32_t,int32_t> edge map key. The
    // pair is the directed edge, so the vertex indices are not sorted.
    protected static edgeKey(v0: number, v1: number): string {
        return v0 + ',' + v1;
    }

    // The vertex data and default vertex creation.
    protected static createVertex(v0: number): VEManifoldMeshVertex {
        return new VEManifoldMeshVertex(v0);
    }

    // The edge data and default edge creation.
    protected static createEdge(v0: number, v1: number): VEManifoldMeshEdge {
        return new VEManifoldMeshEdge(v0, v1);
    }
}
