// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VETNonmanifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The VETNonmanifoldMesh class represents an edge-triangle nonmanifold mesh
// but additionally stores vertex adjacency information.
//
// Port notes (following the B63/B73 mesh precedent):
//   - The upstream std::map<int32_t, std::shared_ptr<Vertex>> becomes a Map
//     keyed by the vertex index; getVertices() sorts by index to reproduce
//     the ordered iteration of std::map.
//   - The upstream std::set<int32_t> VAdjacent and the std::set of shared
//     pointers EAdjacent/TAdjacent (ordered by SharedPtrLT, which forwards to
//     Edge::operator< and Triangle::operator<, i.e. to the feature keys)
//     become Sets with sorted accessors getVAdjacent()/getEAdjacent()/
//     getTAdjacent().
//   - The upstream VCreator function pointer becomes an optional constructor
//     callback, and the copy constructor and operator= become clone() and
//     assign().

import { EdgeKey } from './EdgeKey';
import { FeatureKey } from './FeatureKey';
import {
    ETNonmanifoldMesh, ETNonmanifoldMeshEdge, ETNonmanifoldMeshTriangle
} from './ETNonmanifoldMesh';
import type {
    ETNonmanifoldMeshECreator, ETNonmanifoldMeshTCreator
} from './ETNonmanifoldMesh';
import { logAssert } from './Logger';
import { TriangleKey } from './TriangleKey';

// The port of VETNonmanifoldMesh::Vertex.
export class VETNonmanifoldMeshVertex {
    // The index into the vertex pool of the mesh.
    V: number;

    // Adjacent objects.
    VAdjacent: Set<number>;
    EAdjacent: Set<ETNonmanifoldMeshEdge>;
    TAdjacent: Set<ETNonmanifoldMeshTriangle>;

    constructor(vIndex: number) {
        this.V = vIndex;
        this.VAdjacent = new Set<number>();
        this.EAdjacent = new Set<ETNonmanifoldMeshEdge>();
        this.TAdjacent = new Set<ETNonmanifoldMeshTriangle>();
    }

    // The port of Vertex::operator<.
    lessThan(other: VETNonmanifoldMeshVertex): boolean {
        return this.V < other.V;
    }

    // The adjacent vertex indices in increasing order, the order in which the
    // upstream std::set<int32_t> iterates.
    getVAdjacent(): number[] {
        return Array.from(this.VAdjacent).sort((a, b) => a - b);
    }

    // The adjacent edges in increasing order of the unordered edge key, the
    // order in which the upstream std::set with SharedPtrLT iterates.
    getEAdjacent(): ETNonmanifoldMeshEdge[] {
        const edges = Array.from(this.EAdjacent);
        edges.sort((e0, e1) => FeatureKey.compare(
            new EdgeKey(false, e0.V[0], e0.V[1]),
            new EdgeKey(false, e1.V[0], e1.V[1])));
        return edges;
    }

    // The adjacent triangles in increasing order of the ordered triangle key,
    // the order in which the upstream std::set with SharedPtrLT iterates.
    getTAdjacent(): ETNonmanifoldMeshTriangle[] {
        const triangles = Array.from(this.TAdjacent);
        triangles.sort((t0, t1) => FeatureKey.compare(
            new TriangleKey(true, t0.V[0], t0.V[1], t0.V[2]),
            new TriangleKey(true, t1.V[0], t1.V[1], t1.V[2])));
        return triangles;
    }
}

export type VETNonmanifoldMeshVCreator = (vIndex: number) => VETNonmanifoldMeshVertex;

export class VETNonmanifoldMesh extends ETNonmanifoldMesh {
    protected mVCreator: VETNonmanifoldMeshVCreator;
    protected mVMap: Map<number, VETNonmanifoldMeshVertex>;

    constructor(vCreator?: VETNonmanifoldMeshVCreator,
        eCreator?: ETNonmanifoldMeshECreator, tCreator?: ETNonmanifoldMeshTCreator) {
        super(eCreator, tCreator);
        this.mVCreator = vCreator ?? VETNonmanifoldMesh.createVertex;
        this.mVMap = new Map<number, VETNonmanifoldMeshVertex>();
    }

    // Support for a deep copy of the mesh. The vertex, edge, and triangle
    // objects are not shared between the meshes; the copy reinserts the
    // triangles of the source mesh, which is the upstream operator=.
    override assign(mesh: ETNonmanifoldMesh): this {
        this.clear();
        if (mesh instanceof VETNonmanifoldMesh) {
            this.mVCreator = mesh.mVCreator;
        }
        super.assign(mesh);
        return this;
    }

    // The port of the upstream copy constructor.
    override clone(): VETNonmanifoldMesh {
        return new VETNonmanifoldMesh().assign(this);
    }

    // Member access. The vertices are in increasing order of the vertex
    // index, the order in which the upstream std::map iterates.
    getVertices(): VETNonmanifoldMeshVertex[] {
        const vertices = Array.from(this.mVMap.values());
        vertices.sort((v0, v1) => v0.V - v1.V);
        return vertices;
    }

    getVertex(vIndex: number): VETNonmanifoldMeshVertex | null {
        return this.mVMap.get(vIndex) ?? null;
    }

    getNumVertices(): number {
        return this.mVMap.size;
    }

    // If <v0,v1,v2> is not in the mesh, a Triangle object is created and
    // returned; otherwise, <v0,v1,v2> is in the mesh and null is returned.
    override insert(v0: number, v1: number, v2: number): ETNonmanifoldMeshTriangle | null {
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

    // If <v0,v1,v2> is in the mesh, it is removed and 'true' is returned;
    // otherwise, <v0,v1,v2> is not in the mesh and 'false' is returned.
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

                // If the edge will be removed by ETNonmanifoldMesh::remove,
                // remove the vertex references to it. (Upstream wraps the
                // body in a loop over edge->T; that loop always executes
                // exactly once here because edge->T.size() == 1, so the port
                // drops it.)
                if (edge.T.size === 1) {
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
                // Upstream asserts
                //   VAdjacent.size() != 0 || EAdjacent.size() != 0
                // which is inverted: when the last triangle sharing the
                // vertex is removed, every edge of that triangle containing
                // the vertex has been dropped from VAdjacent and EAdjacent
                // above, so both sets are empty and the upstream assertion
                // always fires. The port asserts that both sets are empty,
                // which is the condition the message ("Malformed mesh.")
                // describes.
                logAssert(vertex.VAdjacent.size === 0 && vertex.EAdjacent.size === 0,
                    'Malformed mesh.');

                this.mVMap.delete(vIndex);
            }
        }

        return super.remove(v0, v1, v2);
    }

    // Destroy the vertices, edges, and triangles to obtain an empty mesh.
    override clear(): void {
        this.mVMap.clear();
        super.clear();
    }

    // The vertex data and default vertex creation.
    protected static createVertex(vIndex: number): VETNonmanifoldMeshVertex {
        return new VETNonmanifoldMeshVertex(vIndex);
    }
}
