// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VTSManifoldMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The VTSManifoldMesh class represents a triangle-tetrahedron manifold mesh
// but additionally stores vertex adjacency information. The 'V' stands for
// vertex, the 'T' stands for triangle (face) and the 'S' stands for simplex
// (tetrahedron). It is general purpose, allowing insertion and removal of
// tetrahedra at any time. However, the performance is limited because of the
// use of the container classes. If your application requires a static
// vertex-triangle-simplex manifold mesh for which no modifications will
// occur, a better choice is StaticVTSManifoldMesh.
//
// Port notes (following the B73 TSManifoldMesh and B84 VETNonmanifoldMesh
// precedent):
//   - The upstream std::unordered_map<int32_t, std::unique_ptr<Vertex>>
//     becomes a Map keyed by the vertex index; getVertices() sorts by index
//     so that the observable iteration order is deterministic (the upstream
//     container is unordered, so its order is implementation-defined).
//   - The upstream std::unordered_set members VAdjacent, TAdjacent and
//     SAdjacent become Sets, with sorted accessors getVAdjacent(),
//     getTAdjacent() and getSAdjacent() for deterministic iteration.
//   - The VCreator function pointer becomes an optional constructor callback,
//     and the copy constructor and operator= become clone() and assign().
//   - Upstream bug, FIXED here: see the comment in remove().

import { FeatureKey } from './FeatureKey';
import { logAssert } from './Logger';
import { TetrahedronKey } from './TetrahedronKey';
import { TriangleKey } from './TriangleKey';
import {
    TSManifoldMesh, TSManifoldMeshTetrahedron, TSManifoldMeshTriangle
} from './TSManifoldMesh';
import type { TSManifoldMeshSCreator, TSManifoldMeshTCreator } from './TSManifoldMesh';

// The port of VTSManifoldMesh::Vertex.
export class VTSManifoldMeshVertex {
    // The index into the vertex pool of the mesh.
    V: number;

    // Adjacent objects.
    VAdjacent: Set<number>;
    TAdjacent: Set<TSManifoldMeshTriangle>;
    SAdjacent: Set<TSManifoldMeshTetrahedron>;

    constructor(vIndex: number) {
        this.V = vIndex;
        this.VAdjacent = new Set<number>();
        this.TAdjacent = new Set<TSManifoldMeshTriangle>();
        this.SAdjacent = new Set<TSManifoldMeshTetrahedron>();
    }

    // The adjacent vertex indices in increasing order.
    getVAdjacent(): number[] {
        return Array.from(this.VAdjacent).sort((a, b) => a - b);
    }

    // The adjacent faces in increasing order of the unordered triangle key.
    getTAdjacent(): TSManifoldMeshTriangle[] {
        const faces = Array.from(this.TAdjacent);
        faces.sort((t0, t1) => FeatureKey.compare(
            new TriangleKey(false, t0.V[0], t0.V[1], t0.V[2]),
            new TriangleKey(false, t1.V[0], t1.V[1], t1.V[2])));
        return faces;
    }

    // The adjacent tetrahedra in increasing order of the ordered tetrahedron
    // key.
    getSAdjacent(): TSManifoldMeshTetrahedron[] {
        const tetras = Array.from(this.SAdjacent);
        tetras.sort((s0, s1) => FeatureKey.compare(
            new TetrahedronKey(true, s0.V[0], s0.V[1], s0.V[2], s0.V[3]),
            new TetrahedronKey(true, s1.V[0], s1.V[1], s1.V[2], s1.V[3])));
        return tetras;
    }
}

export type VTSManifoldMeshVCreator = (vIndex: number) => VTSManifoldMeshVertex;

export class VTSManifoldMesh extends TSManifoldMesh {
    protected mVCreator: VTSManifoldMeshVCreator;
    protected mVMap: Map<number, VTSManifoldMeshVertex>;

    constructor(vCreator?: VTSManifoldMeshVCreator, tCreator?: TSManifoldMeshTCreator,
        sCreator?: TSManifoldMeshSCreator) {
        super(tCreator, sCreator);
        this.mVCreator = vCreator ?? VTSManifoldMesh.createVertex;
        this.mVMap = new Map<number, VTSManifoldMeshVertex>();
    }

    // Support for a deep copy of the mesh; the port of operator=. The vertex,
    // face and tetrahedron objects are not shared between the meshes; the
    // copy reinserts the tetrahedra of the source mesh.
    override assign(mesh: TSManifoldMesh): this {
        this.clear();
        if (mesh instanceof VTSManifoldMesh) {
            this.mVCreator = mesh.mVCreator;
        }
        super.assign(mesh);
        return this;
    }

    // The port of the upstream copy constructor.
    override clone(): VTSManifoldMesh {
        return new VTSManifoldMesh().assign(this);
    }

    // Member access. The vertices are in increasing order of the vertex
    // index.
    getVertices(): VTSManifoldMeshVertex[] {
        const vertices = Array.from(this.mVMap.values());
        vertices.sort((v0, v1) => v0.V - v1.V);
        return vertices;
    }

    getVertex(vIndex: number): VTSManifoldMeshVertex | null {
        return this.mVMap.get(vIndex) ?? null;
    }

    getNumVertices(): number {
        return this.mVMap.size;
    }

    // If the tetrahedron is not in the mesh, a Tetrahedron object is created
    // and returned; otherwise, the tetrahedron is in the mesh and null is
    // returned. If the insertion leads to a nonmanifold mesh, the call fails
    // with null returned (or an exception, see throwOnNonmanifoldInsertion).
    override insert(v0: number, v1: number, v2: number, v3: number):
        TSManifoldMeshTetrahedron | null {
        const tetra = super.insert(v0, v1, v2, v3);
        if (!tetra) {
            return null;
        }

        for (let i = 0; i < 4; ++i) {
            const vIndex = tetra.V[i];
            let vertex = this.mVMap.get(vIndex);
            if (vertex === undefined) {
                vertex = this.mVCreator(vIndex);
                this.mVMap.set(vIndex, vertex);
            }

            vertex.SAdjacent.add(tetra);

            for (let j = 0; j < 4; ++j) {
                const tri = tetra.T[j];
                logAssert(tri !== null, 'Unexpected condition.');
                if (tri.V[0] === vIndex) {
                    vertex.VAdjacent.add(tri.V[1]);
                    vertex.VAdjacent.add(tri.V[2]);
                    vertex.TAdjacent.add(tri);
                }
                else if (tri.V[1] === vIndex) {
                    vertex.VAdjacent.add(tri.V[0]);
                    vertex.VAdjacent.add(tri.V[2]);
                    vertex.TAdjacent.add(tri);
                }
                else if (tri.V[2] === vIndex) {
                    vertex.VAdjacent.add(tri.V[0]);
                    vertex.VAdjacent.add(tri.V[1]);
                    vertex.TAdjacent.add(tri);
                }
            }
        }

        return tetra;
    }

    // If the tetrahedron is in the mesh, it is removed and 'true' is
    // returned; otherwise, the tetrahedron is not in the mesh and 'false' is
    // returned.
    override remove(v0: number, v1: number, v2: number, v3: number): boolean {
        const tetra = this.getTetrahedron(v0, v1, v2, v3);
        if (tetra === null) {
            return false;
        }

        for (let i = 0; i < 4; ++i) {
            const vIndex = tetra.V[i];
            const vertex = this.mVMap.get(vIndex);
            logAssert(vertex !== undefined, 'Unexpected condition.');
            for (let j = 0; j < 4; ++j) {
                const tri = tetra.T[j];
                logAssert(tri !== null, 'Unexpected condition.');

                // When the face is shared only by the tetrahedron being
                // removed, TSManifoldMesh.remove will destroy it, so remove
                // the vertex reference to it.
                if (tri.S[0] && !tri.S[1]) {
                    if (tri.V[0] === vIndex || tri.V[1] === vIndex ||
                        tri.V[2] === vIndex) {
                        vertex.TAdjacent.delete(tri);
                    }
                }
            }

            vertex.SAdjacent.delete(tetra);

            // Upstream bug, fixed here. Upstream erases the two other
            // vertices of each destroyed face from VAdjacent. Unlike the
            // edge-based VETManifoldMesh, where an unordered vertex pair
            // belongs to exactly one edge, a vertex pair can belong to
            // several faces, so erasing on behalf of one destroyed face drops
            // adjacencies that surviving faces still support. For example,
            // for tetrahedra <0,1,2,3> and <0,1,2,4> sharing the face
            // <0,1,2>, removing <0,1,2,3> erases 1 and 2 from the VAdjacent
            // set of vertex 0 even though the faces of <0,1,2,4> still make 0
            // adjacent to both. The invariant maintained by insert() is that
            // VAdjacent is the union, over the faces in TAdjacent, of the
            // other two vertices of each face; the port restores that
            // invariant after the destroyed faces have been removed from
            // TAdjacent, which agrees with upstream whenever upstream is
            // correct.
            vertex.VAdjacent.clear();
            for (const tri of vertex.TAdjacent) {
                for (let k = 0; k < 3; ++k) {
                    if (tri.V[k] === vIndex) {
                        vertex.VAdjacent.add(tri.V[(k + 1) % 3]);
                        vertex.VAdjacent.add(tri.V[(k + 2) % 3]);
                        break;
                    }
                }
            }

            // If the vertex is no longer shared by any tetrahedron, remove it.
            if (vertex.SAdjacent.size === 0) {
                logAssert(vertex.VAdjacent.size === 0 && vertex.TAdjacent.size === 0,
                    'Unexpected condition.');

                this.mVMap.delete(vIndex);
            }
        }

        return super.remove(v0, v1, v2, v3);
    }

    // Destroy the vertices, faces and tetrahedra to obtain an empty mesh.
    override clear(): void {
        this.mVMap.clear();
        super.clear();
    }

    // The vertex data and default vertex creation.
    protected static createVertex(vIndex: number): VTSManifoldMeshVertex {
        return new VTSManifoldMeshVertex(vIndex);
    }
}
