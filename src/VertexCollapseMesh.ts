// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VertexCollapseMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Decimate a triangle mesh by repeatedly removing interior vertices. The
// triangles that share the vertex are removed and the polygon boundary of
// that triangle fan (the "link" of the vertex) is retriangulated. The
// vertices are visited in order of a weight that measures the local area and
// the deviation of the neighbors from the tangent plane at the vertex, so
// nearly planar neighborhoods are collapsed first.
//
// Port notes:
// - Upstream 'template <typename Real>' becomes 'number'. The positions are
//   3D 'Vector' objects (the array is referenced, not copied, as upstream).
// - The nested struct 'Record' and the 'bool DoCollapse(Record&)' output
//   parameter become the returned 'VertexCollapseMeshResult' object with the
//   boolean in its 'collapsed' field.
// - The private nested 'VCVertex' derived from VETManifoldMesh::Vertex
//   becomes the exported class VertexCollapseMeshVertex, supplied to the mesh
//   through the VETManifoldMesh vertex-creator callback.
// - The upstream private enumerates VCM_* become the module-private
//   'CollapseStatus'. VCM_NO_MORE_ALLOWED is never returned by upstream; it
//   is retained for documentation only.
// - std::map<int32_t,int32_t> edgeMap becomes a Map with an explicit sorted
//   start key, so the traversal matches the upstream std::map ordering.

import type { ETManifoldMeshTriangle } from './ETManifoldMesh';
import { logAssert } from './Logger';
import { MinHeap } from './MinHeap';
import type { MinHeapRecord } from './MinHeap';
import { Polygon2 } from './Polygon2';
import { TriangleKey } from './TriangleKey';
import { TriangulateEC } from './TriangulateEC';
import { Vector, dot, length as vectorLength, normalize, sub } from './Vector';
import { computeOrthogonalComplement3, cross } from './Vector3';
import { VETManifoldMesh, VETManifoldMeshVertex } from './VETManifoldMesh';

// The upstream invalid-vertex marker 0x80000000, which is INT32_MIN when
// stored in the int32_t 'record.vertex'.
export const VERTEX_COLLAPSE_MESH_INVALID_VERTEX = -0x80000000;

// The port of VertexCollapseMesh::VCVertex.
export class VertexCollapseMeshVertex extends VETManifoldMeshVertex {
    // The area-weighted average of the normals of the triangles sharing the
    // vertex, normalized. It is computed as a side effect of computeWeight.
    normal: Vector;
    isBoundary: boolean;

    constructor(vIndex: number) {
        super(vIndex);
        this.normal = new Vector(3);
        this.isBoundary = false;
    }

    // The weight depends on the area of the triangles sharing the vertex and
    // the lengths of the projections of the adjacent vertices onto the vertex
    // normal line. A side effect of the call is that the vertex normal is
    // computed and stored.
    computeWeight(positions: readonly Vector[]): number {
        let weight = 0;

        const normal = new Vector(3);
        for (const tri of this.getTAdjacent()) {
            const E0 = sub(positions[tri.V[1]], positions[tri.V[0]]);
            const E1 = sub(positions[tri.V[2]], positions[tri.V[0]]);
            const N = cross(E0, E1);
            normal.values[0] += N.values[0];
            normal.values[1] += N.values[1];
            normal.values[2] += N.values[2];
            weight += vectorLength(N);
        }
        normalize(normal);
        this.normal = normal;

        for (const index of this.getVAdjacent()) {
            const diff = sub(positions[index], positions[this.V]);
            weight += Math.abs(dot(normal, diff));
        }

        return weight;
    }
}

// The port of VertexCollapseMesh::Record, extended with the boolean that
// upstream returns from DoCollapse.
export interface VertexCollapseMeshResult {
    // True when a vertex collapse occurred. Once it is false, no more vertex
    // collapses are allowed, so you may then stop calling doCollapse.
    collapsed: boolean;

    // The index of the interior vertex that is removed from the mesh. The
    // triangles adjacent to the vertex are 'removed' from the mesh. The
    // polygon boundary of the adjacent triangles is triangulated and the new
    // triangles are 'inserted' into the mesh. When 'collapsed' is false, the
    // vertex is VERTEX_COLLAPSE_MESH_INVALID_VERTEX.
    vertex: number;
    removed: TriangleKey[];
    inserted: TriangleKey[];
}

// The functions triangulateLink and collapsed return one of these:
//
// NO_MORE_ALLOWED:
//     Either the mesh has no more interior vertices or a collapse will lead
//     to a mesh fold-over or to a nonmanifold mesh. (Upstream declares this
//     enumerate but never returns it.)
//
// ALLOWED:
//     An interior vertex v may be removed. The vertex normal is the weighted
//     average of non-unit-length normals of triangles sharing v; the weights
//     are the triangle areas. The adjacent vertices are projected onto a
//     plane containing v and having normal equal to the vertex normal. If the
//     projection is a simple polygon in the plane, the collapse is allowed.
//
// DEFERRED:
//     The projection polygon is not simple (at least one pair of edges
//     overlaps at some edge-interior point), so the collapse would produce a
//     fold-over in the mesh, or the retriangulation would create a
//     nonmanifold edge. It is possible that such a vertex may be collapsed
//     later as its neighbors are adjusted by other collapses.
//
// UNEXPECTED_ERROR:
//     A condition occurred that is not expected for a theoretically correct
//     implementation.
enum CollapseStatus {
    NO_MORE_ALLOWED,
    ALLOWED,
    DEFERRED,
    UNEXPECTED_ERROR
}

// The port of the intermediate state shared by triangulateLink and
// collapsed (upstream passes them as reference parameters).
interface LinkInfo {
    status: CollapseStatus;
    removed: TriangleKey[];
    inserted: TriangleKey[];
    linkVertices: number[];
}

export class VertexCollapseMesh {
    private mNumPositions: number;
    private mPositions: readonly Vector[];
    private mMesh: VETManifoldMesh;
    private mMinHeap: MinHeap<number, number>;
    private mHeapRecords: Map<number, MinHeapRecord<number, number>>;

    // Construction. The 'positions' array is referenced, not copied (as
    // upstream, which stores a raw pointer). The 'indices' array is a list of
    // triples of indices into 'positions'. If the input is invalid, the
    // object is constructed in a state for which doCollapse always reports
    // 'collapsed: false'.
    constructor(positions: readonly Vector[] | null,
        indices: readonly number[] | null) {
        this.mNumPositions = 0;
        this.mPositions = [];
        this.mMesh = new VETManifoldMesh(
            (v: number) => new VertexCollapseMeshVertex(v));
        this.mMinHeap = new MinHeap<number, number>(0);
        this.mHeapRecords = new Map<number, MinHeapRecord<number, number>>();

        if (positions === null || positions === undefined ||
            positions.length <= 0 || indices === null ||
            indices === undefined || indices.length < 3) {
            return;
        }

        this.mNumPositions = positions.length;
        this.mPositions = positions;

        // Build the manifold mesh from the inputs.
        const numTriangles = Math.floor(indices.length / 3);
        for (let t = 0, j = 0; t < numTriangles; ++t) {
            const v0 = indices[j++];
            const v1 = indices[j++];
            const v2 = indices[j++];
            this.mMesh.insert(v0, v1, v2);
        }

        // Locate the vertices (if any) on the mesh boundary.
        for (const edge of this.mMesh.getEdges()) {
            if (!edge.T[1]) {
                for (let i = 0; i < 2; ++i) {
                    const vertex = this.getVCVertex(edge.V[i]);
                    logAssert(vertex !== null, 'Unexpected condition.');
                    vertex.isBoundary = true;
                }
            }
        }

        // Build the priority queue of weights for the interior vertices.
        const vertices = this.mMesh.getVertices();
        this.mMinHeap.reset(vertices.length);
        for (const velement of vertices) {
            const vertex = velement as VertexCollapseMeshVertex;
            const weight = (vertex.isBoundary ? Number.MAX_VALUE :
                vertex.computeWeight(this.mPositions));
            const record = this.mMinHeap.insert(vertex.V, weight);
            if (record !== null) {
                this.mHeapRecords.set(vertex.V, record);
            }
        }
    }

    // Access the current state of the mesh, whether the original built in the
    // constructor or a decimated mesh during doCollapse calls.
    getMesh(): VETManifoldMesh {
        return this.mMesh;
    }

    // Decimate the mesh using vertex collapses. The 'collapsed' field of the
    // result is true when a vertex collapse occurs. Once it is false, no more
    // vertex collapses are allowed so you may then stop calling the function.
    // The implementation has several consistency tests that should not fail
    // with a theoretically correct implementation. If a test fails, the
    // function reports 'collapsed: false' and the vertex is set to the
    // invalid integer VERTEX_COLLAPSE_MESH_INVALID_VERTEX.
    doCollapse(): VertexCollapseMeshResult {
        const failure: VertexCollapseMeshResult = {
            collapsed: false,
            vertex: VERTEX_COLLAPSE_MESH_INVALID_VERTEX,
            removed: [],
            inserted: []
        };

        if (this.mNumPositions === 0) {
            // The constructor failed, so there is nothing to collapse.
            return failure;
        }

        while (this.mMinHeap.getNumElements() > 0) {
            const minimum = this.mMinHeap.getMinimum();
            if (minimum === null) {
                // Unexpected condition.
                return failure;
            }

            const v = minimum.key;
            if (minimum.value === Number.MAX_VALUE) {
                // There are no more interior vertices to collapse.
                return failure;
            }

            const vertex = this.getVCVertex(v);
            if (vertex === null) {
                // Unexpected condition.
                return failure;
            }

            const info = this.triangulateLink(vertex);
            if (info.status === CollapseStatus.UNEXPECTED_ERROR) {
                return failure;
            }

            if (info.status === CollapseStatus.ALLOWED) {
                const status = this.collapse(info.removed, info.inserted,
                    info.linkVertices);
                if (status === CollapseStatus.UNEXPECTED_ERROR) {
                    return failure;
                }

                if (status === CollapseStatus.ALLOWED) {
                    // Remove the vertex and associated weight.
                    this.mMinHeap.remove();
                    this.mHeapRecords.delete(v);

                    // Update the weights of the link vertices.
                    for (const vlink of info.linkVertices) {
                        const linkVertex = this.getVCVertex(vlink);
                        if (linkVertex === null) {
                            // Unexpected condition.
                            return failure;
                        }

                        if (!linkVertex.isBoundary) {
                            const record = this.mHeapRecords.get(vlink);
                            if (record === undefined) {
                                // Unexpected condition.
                                return failure;
                            }

                            this.mMinHeap.update(record,
                                linkVertex.computeWeight(this.mPositions));
                        }
                    }

                    return {
                        collapsed: true,
                        vertex: v,
                        removed: info.removed,
                        inserted: info.inserted
                    };
                }
                // else: status === CollapseStatus.DEFERRED
            }

            // To get here, the status must be DEFERRED. The vertex collapse
            // would cause mesh fold-over. Temporarily set the vertex weight
            // to infinity. After removal of other triangles, the vertex
            // weight will be updated to a finite value and the vertex
            // possibly can be removed at that time.
            const record = this.mHeapRecords.get(v);
            if (record === undefined) {
                // Unexpected condition.
                return failure;
            }
            this.mMinHeap.update(record, Number.MAX_VALUE);
        }

        // We do not expect to reach this line of code, even for a closed
        // mesh.
        return failure;
    }

    private getVCVertex(vIndex: number): VertexCollapseMeshVertex | null {
        const vertex = this.mMesh.getVertex(vIndex);
        return (vertex !== null ? vertex as VertexCollapseMeshVertex : null);
    }

    private triangulateLink(vertex: VertexCollapseMeshVertex): LinkInfo {
        const result: LinkInfo = {
            status: CollapseStatus.UNEXPECTED_ERROR,
            removed: [],
            inserted: [],
            linkVertices: []
        };

        // Create the (CCW) polygon boundary of the link of the vertex. The
        // incoming vertex is interior, so the number of triangles sharing the
        // vertex is equal to the number of vertices of the polygon. A
        // precondition of the function call is that the vertex normal has
        // already been computed.

        // Get the edges of the link that are opposite the incoming vertex.
        const adjacent: ETManifoldMeshTriangle[] = vertex.getTAdjacent();
        const numVertices = adjacent.length;
        const removed = new Array<TriangleKey>(numVertices);
        const edgeMap = new Map<number, number>();
        let j = 0;
        for (const tri of adjacent) {
            for (let i = 0; i < 3; ++i) {
                if (tri.V[i] === vertex.V) {
                    const key = tri.V[(i + 1) % 3];
                    // std::map::insert does not overwrite an existing key.
                    if (!edgeMap.has(key)) {
                        edgeMap.set(key, tri.V[(i + 2) % 3]);
                    }
                    break;
                }
            }
            removed[j++] = new TriangleKey(true, tri.V[0], tri.V[1], tri.V[2]);
        }
        if (edgeMap.size !== numVertices) {
            return result;
        }

        // Connect the edges into a polygon. The upstream traversal starts at
        // the smallest key of the std::map.
        const linkVertices = new Array<number>(numVertices);
        let current = Math.min(...edgeMap.keys());
        for (let i = 0; i < numVertices; ++i) {
            linkVertices[i] = current;
            const next = edgeMap.get(current);
            if (next === undefined || !edgeMap.has(next)) {
                return result;
            }
            current = next;
        }
        if (current !== linkVertices[0]) {
            return result;
        }

        // Project the polygon onto the plane containing the incoming vertex
        // and having the vertex normal. The projected polygon is computed so
        // that the incoming vertex is projected to (0,0).
        const center = this.mPositions[vertex.V];
        const basis: Vector[] = [vertex.normal.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const projected = new Array<Vector>(numVertices);
        const indices = new Array<number>(numVertices);
        for (let i = 0; i < numVertices; ++i) {
            const diff = sub(this.mPositions[linkVertices[i]], center);
            projected[i] = Vector.fromArray(
                [dot(basis[1], diff), dot(basis[2], diff)]);
            indices[i] = i;
        }

        // The polygon must be simple in order to triangulate it.
        const polygon = new Polygon2(projected, indices, true);
        if (polygon.isValid() && polygon.isSimple()) {
            const triangulator = new TriangulateEC(projected, numVertices);
            triangulator.triangulate();
            const triangles = triangulator.getTriangles();
            if (triangles.length === 0) {
                return result;
            }

            const inserted = new Array<TriangleKey>(triangles.length);
            for (let t = 0; t < triangles.length; ++t) {
                inserted[t] = new TriangleKey(true,
                    linkVertices[triangles[t][0]],
                    linkVertices[triangles[t][1]],
                    linkVertices[triangles[t][2]]);
            }

            result.status = CollapseStatus.ALLOWED;
            result.removed = removed;
            result.inserted = inserted;
            result.linkVertices = linkVertices;
            return result;
        }

        result.status = CollapseStatus.DEFERRED;
        result.linkVertices = linkVertices;
        return result;
    }

    private collapse(removed: readonly TriangleKey[],
        inserted: readonly TriangleKey[],
        linkVertices: readonly number[]): CollapseStatus {
        // The triangles that were disconnected from the link edges are
        // guaranteed to allow manifold reconnection to the 'inserted'
        // triangles. On the insertion, each diagonal of the link becomes a
        // mesh edge and shares two (link) triangles. It is possible that the
        // mesh already contains the (diagonal) edge, which will lead to a
        // nonmanifold connection, which we cannot allow. The following code
        // traps this condition before any triangle is removed.
        const edges = new Set<string>();
        for (const tri of inserted) {
            for (let k0 = 2, k1 = 0; k1 < 3; k0 = k1++) {
                const v0 = Math.min(tri.V[k0], tri.V[k1]);
                const v1 = Math.max(tri.V[k0], tri.V[k1]);
                const key = `${v0},${v1}`;
                if (!edges.has(key)) {
                    edges.add(key);
                }
                else {
                    // The edge has been visited twice, so it is a diagonal of
                    // the link.
                    const edge = this.mMesh.getEdge(v0, v1);
                    if (edge !== null && edge.T[1]) {
                        // The edge will not allow a manifold connection.
                        return CollapseStatus.DEFERRED;
                    }

                    edges.delete(key);
                }
            }
        }

        // Remove the old triangle neighborhood, which will lead to the vertex
        // itself being removed from the mesh.
        for (const tri of removed) {
            this.mMesh.remove(tri.V[0], tri.V[1], tri.V[2]);
        }

        // Insert the new triangulation.
        for (const tri of inserted) {
            this.mMesh.insert(tri.V[0], tri.V[1], tri.V[2]);
        }

        // If the remove(...) calls remove a boundary vertex that is in the
        // link vertices, the insert(...) calls will insert the boundary
        // vertex again. We must re-tag those boundary vertices.
        const numVertices = linkVertices.length;
        for (let i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++) {
            const edge = this.mMesh.getEdge(linkVertices[i0],
                linkVertices[i1]);
            if (edge === null) {
                return CollapseStatus.UNEXPECTED_ERROR;
            }

            if (edge.T[0] && !edge.T[1]) {
                for (let k = 0; k < 2; ++k) {
                    const vertex = this.getVCVertex(edge.V[k]);
                    if (vertex === null) {
                        return CollapseStatus.UNEXPECTED_ERROR;
                    }

                    vertex.isBoundary = true;
                }
            }
        }

        return CollapseStatus.ALLOWED;
    }
}
