// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VETManifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The VETManifoldMesh class represents an edge-triangle manifold mesh but
// additionally stores vertex adjacency information. It is general purpose,
// allowing insertion and removal of triangles at any time. However, the
// performance is limited because of the use of the container classes. If your
// application requires a vertex-edge-triangle manifold mesh for which no
// triangles will be removed, a better choice is StaticVETManifoldMesh.
//
// Port notes (following the B73 ETManifoldMesh and B84 VETNonmanifoldMesh
// precedent):
//   - The upstream std::unordered_map<int32_t, std::unique_ptr<Vertex>>
//     becomes a Map keyed by the vertex index; getVertices() sorts by index
//     so that the observable iteration order is deterministic (the upstream
//     container is unordered, so its order is implementation-defined).
//   - The upstream std::unordered_set members VAdjacent, EAdjacent and
//     TAdjacent become Sets, with sorted accessors getVAdjacent(),
//     getEAdjacent() and getTAdjacent() for deterministic iteration.
//   - The VCreator function pointer becomes an optional constructor callback,
//     and the copy constructor and operator= become clone() and assign().

import { EdgeKey } from './EdgeKey.js';
import { FeatureKey } from './FeatureKey.js';
import {
    ETManifoldMesh, ETManifoldMeshEdge, ETManifoldMeshTriangle
} from './ETManifoldMesh.js';
import type { ETManifoldMeshECreator, ETManifoldMeshTCreator } from './ETManifoldMesh.js';
import { logAssert } from './Logger.js';
import { TriangleKey } from './TriangleKey.js';

// The port of VETManifoldMesh::Vertex.
export class VETManifoldMeshVertex {
    // The index into the vertex pool of the mesh.
    V: number;

    // Adjacent objects.
    VAdjacent: Set<number>;
    EAdjacent: Set<ETManifoldMeshEdge>;
    TAdjacent: Set<ETManifoldMeshTriangle>;

    constructor(vIndex: number) {
        this.V = vIndex;
        this.VAdjacent = new Set<number>();
        this.EAdjacent = new Set<ETManifoldMeshEdge>();
        this.TAdjacent = new Set<ETManifoldMeshTriangle>();
    }

    // The adjacent vertex indices in increasing order.
    getVAdjacent(): number[] {
        return Array.from(this.VAdjacent).sort((a, b) => a - b);
    }

    // The adjacent edges in increasing order of the unordered edge key.
    getEAdjacent(): ETManifoldMeshEdge[] {
        const edges = Array.from(this.EAdjacent);
        edges.sort((e0, e1) => FeatureKey.compare(
            new EdgeKey(false, e0.V[0], e0.V[1]),
            new EdgeKey(false, e1.V[0], e1.V[1])));
        return edges;
    }

    // The adjacent triangles in increasing order of the ordered triangle key.
    getTAdjacent(): ETManifoldMeshTriangle[] {
        const triangles = Array.from(this.TAdjacent);
        triangles.sort((t0, t1) => FeatureKey.compare(
            new TriangleKey(true, t0.V[0], t0.V[1], t0.V[2]),
            new TriangleKey(true, t1.V[0], t1.V[1], t1.V[2])));
        return triangles;
    }
}

export type VETManifoldMeshVCreator = (vIndex: number) => VETManifoldMeshVertex;

export class VETManifoldMesh extends ETManifoldMesh {
    protected mVCreator: VETManifoldMeshVCreator;
    protected mVMap: Map<number, VETManifoldMeshVertex>;

    constructor(vCreator?: VETManifoldMeshVCreator, eCreator?: ETManifoldMeshECreator,
        tCreator?: ETManifoldMeshTCreator) {
        super(eCreator, tCreator);
        this.mVCreator = vCreator ?? VETManifoldMesh.createVertex;
        this.mVMap = new Map<number, VETManifoldMeshVertex>();
    }

    // Support for a deep copy of the mesh; the port of operator=. The vertex,
    // edge and triangle objects are not shared between the meshes; the copy
    // reinserts the triangles of the source mesh.
    override assign(mesh: ETManifoldMesh): this {
        this.clear();
        if (mesh instanceof VETManifoldMesh) {
            this.mVCreator = mesh.mVCreator;
        }
        super.assign(mesh);
        return this;
    }

    // The port of the upstream copy constructor.
    override clone(): VETManifoldMesh {
        return new VETManifoldMesh().assign(this);
    }

    // Member access. The vertices are in increasing order of the vertex
    // index.
    getVertices(): VETManifoldMeshVertex[] {
        const vertices = Array.from(this.mVMap.values());
        vertices.sort((v0, v1) => v0.V - v1.V);
        return vertices;
    }

    getVertex(vIndex: number): VETManifoldMeshVertex | null {
        return this.mVMap.get(vIndex) ?? null;
    }

    getNumVertices(): number {
        return this.mVMap.size;
    }

    // If the triangle is not in the mesh, a Triangle object is created and
    // returned; otherwise, the triangle is in the mesh and null is returned.
    // If the insertion leads to a nonmanifold mesh, the call fails with null
    // returned (or an exception, see throwOnNonmanifoldInsertion).
    override insert(v0: number, v1: number, v2: number): ETManifoldMeshTriangle | null {
        const tri = super.insert(v0, v1, v2);
        if (!tri) {
            return null;
        }

        for (let i = 0; i < 3; ++i) {
            const vIndex = tri.V[i];
            let vertex = this.mVMap.get(vIndex);
            if (vertex === undefined) {
                vertex = this.mVCreator(vIndex);
                this.mVMap.set(vIndex, vertex);
            }

            vertex.TAdjacent.add(tri);

            for (let j = 0; j < 3; ++j) {
                const edge = tri.E[j];
                logAssert(edge !== null, 'Unexpected condition.');
                if (edge.V[0] === vIndex) {
                    vertex.VAdjacent.add(edge.V[1]);
                    vertex.EAdjacent.add(edge);
                }
                else if (edge.V[1] === vIndex) {
                    vertex.VAdjacent.add(edge.V[0]);
                    vertex.EAdjacent.add(edge);
                }
            }
        }

        return tri;
    }

    // If the triangle is in the mesh, it is removed and 'true' is returned;
    // otherwise, the triangle is not in the mesh and 'false' is returned.
    override remove(v0: number, v1: number, v2: number): boolean {
        const tri = this.getTriangle(v0, v1, v2);
        if (tri === null) {
            return false;
        }

        for (let i = 0; i < 3; ++i) {
            const vIndex = tri.V[i];
            const vertex = this.mVMap.get(vIndex);
            logAssert(vertex !== undefined, 'Unexpected condition.');
            for (let j = 0; j < 3; ++j) {
                const edge = tri.E[j];
                logAssert(edge !== null, 'Unexpected condition.');

                // When the edge is shared only by the triangle being removed,
                // ETManifoldMesh.remove will destroy it, so remove the vertex
                // references to it. In a manifold mesh, an unordered vertex
                // pair has at most one edge, so VAdjacent and EAdjacent stay
                // in agreement with the surviving edges.
                if (edge.T[0] && !edge.T[1]) {
                    if (edge.V[0] === vIndex) {
                        vertex.VAdjacent.delete(edge.V[1]);
                        vertex.EAdjacent.delete(edge);
                    }
                    else if (edge.V[1] === vIndex) {
                        vertex.VAdjacent.delete(edge.V[0]);
                        vertex.EAdjacent.delete(edge);
                    }
                }
            }

            vertex.TAdjacent.delete(tri);

            // If the vertex is no longer shared by any triangle, remove it.
            if (vertex.TAdjacent.size === 0) {
                logAssert(vertex.VAdjacent.size === 0 && vertex.EAdjacent.size === 0,
                    'Unexpected condition.');

                this.mVMap.delete(vIndex);
            }
        }

        return super.remove(v0, v1, v2);
    }

    // Destroy the vertices, edges and triangles to obtain an empty mesh.
    override clear(): void {
        this.mVMap.clear();
        super.clear();
    }

    // The vertex data and default vertex creation.
    protected static createVertex(vIndex: number): VETManifoldMeshVertex {
        return new VETManifoldMeshVertex(vIndex);
    }
}
