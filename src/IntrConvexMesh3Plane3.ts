// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrConvexMesh3Plane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection of a convex polyhedron (a ConvexMesh3) and a
// plane. The query reports how the two objects are configured relative to
// each other, and on request it returns the polygon (and a triangulation of
// it) in which they intersect and the two convex polyhedra into which the
// plane splits the input polyhedron.
//
// Port notes:
// - Upstream static_asserts that Real is an arbitrary-precision type with a
//   division operator, because the sign tests Dot(N,X) - c are exact only
//   then. The port instantiates the query for 'number' like every other
//   query in the library; with floating-point input the sign classification
//   is subject to rounding, exactly as it would be for a C++ instantiation
//   with double if the static_assert were removed.
// - The CFG_* and REQ_* constants are static members of the query class, as
//   upstream, so that the flat module index does not need 20 more global
//   names.
// - The upstream 'requested' argument of operator() is the third argument of
//   find(); it defaults to REQ_ALL so that the class satisfies the two-
//   argument FIQuery interface.
// - std::map/std::set of EdgeKey<false> becomes a Map keyed by a string built
//   from the ordered index pair, with explicit sorting wherever upstream
//   relies on the std::map iteration order.
// - The two upstream overloads named GetIntersectionPolygon become the
//   module-private functions getIntersectionPolygonTangential (the coplanar
//   face case) and getIntersectionPolygonSplit (the transverse case).

import { ConvexMesh3 } from './ConvexMesh3';
import type { ConvexMesh3Triangle, ConvexMesh3Vertex } from './ConvexMesh3';
import type { FIQuery } from './FIQuery';
import type { Hyperplane } from './Hyperplane';
import { logAssert, logError } from './Logger';
import { UniqueVerticesSimplices } from './UniqueVerticesSimplices';
import { Vector, add, dot, mul } from './Vector';

// The result of the find-intersection query.
export interface IntrConvexMesh3Plane3FIResult {
    // The configuration describes geometrically how the input convex
    // polyhedron and the plane intersect. It is one of the CFG_* constants.
    configuration: number;

    // The information that was requested from the query, one of the REQ_*
    // constants (or a bitwise combination of them).
    requested: number;

    // The intersection of the convex polyhedron and the plane is either
    // empty, a single vertex, a single edge or a convex polygon. The
    // intersection members have the properties:
    //   empty:    vertices has 0 elements, the mesh is empty
    //   vertex:   vertices has 1 element, the mesh is empty
    //   edge:     vertices has 2 elements, the mesh is empty
    //   polygon:  vertices has 3 or more elements, the mesh is a
    //             triangulation of the convex polygon
    // The convex polygon vertices are listed in the order consistent with
    // that of the positive polyhedron triangles.
    intersectionMesh: ConvexMesh3;
    intersectionPolygon: ConvexMesh3Vertex[];

    // If the configuration is CFG_POS_SIDE* or CFG_SPLIT, this convex
    // polyhedron is the portion of the input convex polyhedron on the
    // positive side of the plane with possibly a vertex or edge on the plane.
    positivePolyhedron: ConvexMesh3;

    // If the configuration is CFG_NEG_SIDE* or CFG_SPLIT, this convex
    // polyhedron is the portion of the input convex polyhedron on the
    // negative side of the plane with possibly a vertex or edge on the plane.
    negativePolyhedron: ConvexMesh3;
}

// The port of the upstream Result default constructor.
export function defaultIntrConvexMesh3Plane3FIResult():
    IntrConvexMesh3Plane3FIResult {
    return {
        configuration: IntrConvexMesh3Plane3FI.CFG_EMPTY,
        requested: IntrConvexMesh3Plane3FI.REQ_CONFIGURATION_ONLY,
        intersectionMesh: new ConvexMesh3(),
        intersectionPolygon: [],
        positivePolyhedron: new ConvexMesh3(),
        negativePolyhedron: new ConvexMesh3()
    };
}

// The unordered-edge key of the port's std::set<EdgeKey<false>> replacement.
function edgeKey(v0: number, v1: number): string {
    return v0 < v1 ? v0 + '_' + v1 : v1 + '_' + v0;
}

// Sort unordered edges the way std::set<EdgeKey<false>> orders them, by
// (V[0], V[1]) with V[0] = min and V[1] = max.
function sortedEdges(edges: Iterable<[number, number]>): [number, number][] {
    const sorted = Array.from(edges);
    sorted.sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]));
    return sorted;
}

export class IntrConvexMesh3Plane3FI implements
    FIQuery<ConvexMesh3, Hyperplane, IntrConvexMesh3Plane3FIResult> {

    // The configuration describes geometrically how the input convex
    // polyhedron and the plane intersect.
    static readonly CFG_EMPTY = 0x00000000;
    static readonly CFG_POS_SIDE = 0x00000010;
    static readonly CFG_NEG_SIDE = 0x00000020;

    // The plane intersects the convex polyhedron transversely. The set of
    // intersection is a convex polygon. The convex polyhedron is split into
    // two convex polyhedra, one on the positive side of the plane and one on
    // the negative side of the plane, both polyhedra sharing the convex
    // polygon of intersection.
    static readonly CFG_SPLIT =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE |
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE;  // 48

    // The convex polyhedron is strictly on the positive side of the plane.
    static readonly CFG_POS_SIDE_STRICT =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE;  // 16

    // The convex polyhedron is on the positive side of the plane with one
    // vertex in the plane.
    static readonly CFG_POS_SIDE_VERTEX =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE | 1;  // 17

    // The convex polyhedron is on the positive side of the plane with one
    // edge in the plane.
    static readonly CFG_POS_SIDE_EDGE =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE | 2;  // 18

    // The convex polyhedron is on the positive side of the plane with a
    // polygonal face in the plane. The face can consist of multiple
    // triangles.
    static readonly CFG_POS_SIDE_POLYGON =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE | 4;  // 20

    // Flags for any of the tangential cases (vertex touching, edge touching,
    // face touching).
    static readonly CFG_POS_SIDE_TANGENT =
        IntrConvexMesh3Plane3FI.CFG_POS_SIDE | 7;  // 23

    // The convex polyhedron is strictly on the negative side of the plane.
    static readonly CFG_NEG_SIDE_STRICT =
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE;  // 32

    // The convex polyhedron is on the negative side of the plane with one
    // vertex in the plane.
    static readonly CFG_NEG_SIDE_VERTEX =
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE | 1;  // 33

    // The convex polyhedron is on the negative side of the plane with one
    // edge in the plane.
    static readonly CFG_NEG_SIDE_EDGE =
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE | 2;  // 34

    // The convex polyhedron is on the negative side of the plane with a
    // polygonal face in the plane. The face can consist of multiple
    // triangles.
    static readonly CFG_NEG_SIDE_POLYGON =
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE | 4;  // 36

    // Flags for any of the tangential cases (vertex touching, edge touching,
    // face touching).
    static readonly CFG_NEG_SIDE_TANGENT =
        IntrConvexMesh3Plane3FI.CFG_NEG_SIDE | 7;  // 39

    // Requested information for the query to compute.
    static readonly REQ_CONFIGURATION_ONLY = 0x00000000;
    static readonly REQ_INTR_MESH = 0x00000001;
    static readonly REQ_INTR_POLYGON = 0x00000002;
    static readonly REQ_INTR_BOTH =
        IntrConvexMesh3Plane3FI.REQ_INTR_MESH |
        IntrConvexMesh3Plane3FI.REQ_INTR_POLYGON;
    static readonly REQ_POLYHEDRON_POS = 0x00000004;
    static readonly REQ_POLYHEDRON_NEG = 0x00000008;
    static readonly REQ_POLYHEDRON_BOTH =
        IntrConvexMesh3Plane3FI.REQ_POLYHEDRON_POS |
        IntrConvexMesh3Plane3FI.REQ_POLYHEDRON_NEG;
    static readonly REQ_ALL = 0x0000000F;

    find(polyhedron: ConvexMesh3, plane: Hyperplane,
        requested: number = IntrConvexMesh3Plane3FI.REQ_ALL):
        IntrConvexMesh3Plane3FIResult {
        const Q = IntrConvexMesh3Plane3FI;
        const result = defaultIntrConvexMesh3Plane3FIResult();
        result.requested = requested;

        // Storage for (Dot(N,X) - c) for each vertex X, where N is a plane
        // normal (not necessarily unit length) and c is the corresponding
        // plane constant.
        let numPositive = 0;
        let numNegative = 0;
        let numZero = 0;
        const numVertices = polyhedron.vertices.length;
        const dotValues: number[] = new Array<number>(numVertices);
        const sign: number[] = new Array<number>(numVertices);
        for (let i = 0; i < numVertices; ++i) {
            dotValues[i] = dot(plane.normal, polyhedron.vertices[i]) -
                plane.constant;
            if (dotValues[i] > 0) {
                sign[i] = +1;
                ++numPositive;
            }
            else if (dotValues[i] < 0) {
                sign[i] = -1;
                ++numNegative;
            }
            else {
                sign[i] = 0;
                ++numZero;
            }
        }

        if (numPositive === 0) {
            result.configuration =
                Q.CFG_NEG_SIDE | (numZero < 3 ? numZero : 4);

            if ((requested & Q.REQ_POLYHEDRON_NEG) !== 0) {
                result.negativePolyhedron = cloneMesh(polyhedron);
            }

            if (numZero > 0 && ((requested & Q.REQ_INTR_BOTH) !== 0)) {
                getIntersection(polyhedron, numZero, sign, result);
            }
        }
        else if (numNegative === 0) {
            result.configuration =
                Q.CFG_POS_SIDE | (numZero < 3 ? numZero : 4);

            if ((requested & Q.REQ_POLYHEDRON_POS) !== 0) {
                result.positivePolyhedron = cloneMesh(polyhedron);
            }

            if (numZero > 0 && ((requested & Q.REQ_INTR_BOTH) !== 0)) {
                getIntersection(polyhedron, numZero, sign, result);
            }
        }
        else {
            result.configuration = Q.CFG_SPLIT;

            if (requested !== Q.REQ_CONFIGURATION_ONLY) {
                splitPolyhedron(polyhedron, dotValues, sign, result);
            }
        }
        return result;
    }
}

// The port of C++ copy assignment of a ConvexMesh3 (value semantics).
function cloneMesh(mesh: ConvexMesh3): ConvexMesh3 {
    const copy = new ConvexMesh3();
    copy.configuration = mesh.configuration;
    copy.vertices = mesh.vertices.map(vertex => vertex.clone());
    copy.triangles = mesh.triangles.map(triangle =>
        [triangle[0], triangle[1], triangle[2]] as ConvexMesh3Triangle);
    return copy;
}

function getIntersection(polyhedron: ConvexMesh3, numZero: number,
    sign: readonly number[], result: IntrConvexMesh3Plane3FIResult): void {
    const Q = IntrConvexMesh3Plane3FI;
    const wantIntrMesh = (result.requested & Q.REQ_INTR_MESH) !== 0;
    const wantIntrPolygon = (result.requested & Q.REQ_INTR_POLYGON) !== 0;

    if (numZero === 1) {
        getIntersectionVertex(polyhedron, sign, wantIntrMesh, wantIntrPolygon,
            result);
    }
    else if (numZero === 2) {
        getIntersectionEdge(polyhedron, sign, wantIntrMesh, wantIntrPolygon,
            result);
    }
    else {  // numZero >= 3
        getIntersectionPolygonTangential(polyhedron, sign, wantIntrMesh,
            wantIntrPolygon, result);
    }
}

function getIntersectionVertex(polyhedron: ConvexMesh3,
    sign: readonly number[], wantIntrMesh: boolean, wantIntrPolygon: boolean,
    result: IntrConvexMesh3Plane3FIResult): void {
    result.intersectionMesh.configuration = ConvexMesh3.CFG_POINT;

    const numVertices = polyhedron.vertices.length;
    for (let i = 0; i < numVertices; ++i) {
        if (sign[i] === 0) {
            if (wantIntrMesh) {
                result.intersectionMesh.vertices = [
                    polyhedron.vertices[i].clone()
                ];
            }
            if (wantIntrPolygon) {
                result.intersectionPolygon = [polyhedron.vertices[i].clone()];
            }
            return;
        }
    }
}

function getIntersectionEdge(polyhedron: ConvexMesh3,
    sign: readonly number[], wantIntrMesh: boolean, wantIntrPolygon: boolean,
    result: IntrConvexMesh3Plane3FIResult): void {
    result.intersectionMesh.configuration = ConvexMesh3.CFG_SEGMENT;

    const found: ConvexMesh3Vertex[] = [];
    const numVertices = polyhedron.vertices.length;
    for (let i = 0; i < numVertices; ++i) {
        if (sign[i] === 0) {
            found.push(polyhedron.vertices[i].clone());
            if (found.length === 2) {
                break;
            }
        }
    }

    if (wantIntrMesh) {
        result.intersectionMesh.vertices = found.map(v => v.clone());
    }
    if (wantIntrPolygon) {
        result.intersectionPolygon = found.map(v => v.clone());
    }
}

function getIntersectionPolygonTangential(polyhedron: ConvexMesh3,
    sign: readonly number[], wantIntrMesh: boolean, wantIntrPolygon: boolean,
    result: IntrConvexMesh3Plane3FIResult): void {
    result.intersectionMesh.configuration = ConvexMesh3.CFG_POLYGON;

    const intersectionMeshTriangles: ConvexMesh3Triangle[] = [];
    for (const triangle of polyhedron.triangles) {
        if (sign[triangle[0]] === 0 && sign[triangle[1]] === 0 &&
            sign[triangle[2]] === 0) {
            intersectionMeshTriangles.push(
                [triangle[0], triangle[1], triangle[2]]);
        }
    }

    if (intersectionMeshTriangles.length === 0) {
        // The polyhedron has three or more vertices in the plane but no
        // triangle face is contained in the plane; this can happen only for a
        // degenerate polyhedron (for example, three collinear vertices in the
        // plane). There is no polygon of intersection to report. Upstream
        // would pass an empty triangle list to UniqueVerticesSimplices.
        return;
    }

    const uvt = new UniqueVerticesSimplices<ConvexMesh3Vertex>(3);
    const { vertices: outVertices, simplices: outTriangles } =
        uvt.removeDuplicateAndUnusedVerticesGrouped(polyhedron.vertices,
            intersectionMeshTriangles);

    if (wantIntrPolygon) {
        // Get the boundary edges with ordering consistent with the triangle
        // face chirality. An edge visited twice is interior, so it is removed.
        const edgeMap = new Map<string, [number, number]>();
        for (const triangle of outTriangles) {
            for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
                const key = edgeKey(triangle[j0], triangle[j1]);
                if (edgeMap.has(key)) {
                    // The edge is now visited twice, so it cannot be a
                    // boundary edge.
                    edgeMap.delete(key);
                }
                else {
                    // The edge is visited the first time, so it might be a
                    // boundary edge. The value is the directed edge in the
                    // chirality of the triangle that produced it.
                    edgeMap.set(key, [triangle[j0], triangle[j1]]);
                }
            }
        }

        // Construct the boundary polygon by walking the directed boundary
        // edges.
        //
        // Port fix for an upstream bug. Upstream builds
        //   std::vector<int32_t> polygonIndices(edgeMap.size(), -1);
        //   for (element : edgeMap)
        //       polygonIndices[element.second[0]] = element.second[1];
        // and then reads polygonIndices[0], polygonIndices[1], ... in index
        // order. That array is a predecessor map indexed by vertex index; it
        // is never traversed as a cycle, so the reported polygon vertices are
        // an arbitrary permutation of the boundary rather than an ordered
        // polygon. Moreover the array is sized by the number of boundary
        // edges while it is indexed by vertex indices, which writes out of
        // bounds whenever the coplanar face has a vertex that is not on its
        // boundary. The port walks the cycle instead, in the chirality of the
        // triangles, which is what the Result documentation promises.
        const nextOf = new Map<number, number>();
        for (const [from, to] of sortedEdges(edgeMap.values())) {
            nextOf.set(from, to);
        }

        const numBoundary = nextOf.size;
        const polygonIndices: number[] = [];
        if (numBoundary > 0) {
            // Start at the smallest vertex index so the traversal is
            // deterministic.
            let current = Math.min(...nextOf.keys());
            for (let i = 0; i < numBoundary; ++i) {
                polygonIndices.push(current);
                const next = nextOf.get(current);
                logAssert(next !== undefined,
                    'The coplanar face boundary is not a closed polygon.');
                current = next as number;
            }
        }

        result.intersectionPolygon = polygonIndices.map(
            index => outVertices[index].clone());
    }

    if (wantIntrMesh) {
        result.intersectionMesh.vertices = outVertices;
        result.intersectionMesh.triangles = outTriangles.map(triangle =>
            [triangle[0], triangle[1], triangle[2]] as ConvexMesh3Triangle);
    }
}

function splitPolyhedron(polyhedron: ConvexMesh3,
    dotValues: readonly number[], sign: readonly number[],
    result: IntrConvexMesh3Plane3FIResult): void {
    const Q = IntrConvexMesh3Plane3FI;
    const wantPosMesh = (result.requested & Q.REQ_POLYHEDRON_POS) !== 0;
    const wantNegMesh = (result.requested & Q.REQ_POLYHEDRON_NEG) !== 0;
    const wantIntrMesh = (result.requested & Q.REQ_INTR_MESH) !== 0;
    const wantIntrPolygon = (result.requested & Q.REQ_INTR_POLYGON) !== 0;

    // The split polyhedra use the input polyhedron's vertices and any
    // edge-interior intersections between the plane and the mesh edges. The
    // center point of the polygon of intersection (if any) is also used as a
    // vertex.
    const { splitVertices, eiVMap } =
        getVertexCandidates(polyhedron, dotValues, sign);

    // Split each triangle face of the polyhedron by the plane.
    const posMesh: ConvexMesh3Triangle[] = [];
    const negMesh: ConvexMesh3Triangle[] = [];
    const posIntersection = new Map<number, number>();
    doSplit(polyhedron, sign, eiVMap, wantPosMesh, posMesh, wantNegMesh,
        negMesh, posIntersection);

    // Get the polygon of intersection. This is used by all of the requested
    // features.
    const polygon = getIntersectionPolygonSplit(posIntersection, splitVertices,
        wantIntrPolygon, result);

    if (wantPosMesh || wantNegMesh || wantIntrMesh) {
        // Get the polyhedra split by the plane. The polygon of intersection is
        // also computed and used to close the polyhedra.
        getSplitPolyhedra(splitVertices, polygon, wantIntrMesh, wantPosMesh,
            posMesh, wantNegMesh, negMesh, result);
    }
}

function getVertexCandidates(polyhedron: ConvexMesh3,
    dotValues: readonly number[], sign: readonly number[]):
    { splitVertices: ConvexMesh3Vertex[], eiVMap: Map<string, number> } {
    // Get the edges of the polyhedron.
    const edgeSet = new Map<string, [number, number]>();
    for (const triangle of polyhedron.triangles) {
        for (const [a, b] of [
            [triangle[0], triangle[1]],
            [triangle[1], triangle[2]],
            [triangle[2], triangle[0]]
        ]) {
            const key = edgeKey(a, b);
            if (!edgeSet.has(key)) {
                edgeSet.set(key, a < b ? [a, b] : [b, a]);
            }
        }
    }

    // The vertex candidates include the original vertices, any edge-interior
    // intersections between the plane and polyhedron, and the average of the
    // convex-polygon intersection (if there is such an intersection).
    const splitVertices: ConvexMesh3Vertex[] =
        polyhedron.vertices.map(vertex => vertex.clone());

    // Compute edge-interior points of intersection between the plane and the
    // mesh edges. The eiVMap container allows accessing the edge-interior
    // vertices when each triangle face of the polyhedron is processed for
    // intersection with the plane.
    const eiVMap = new Map<string, number>();
    for (const [v0, v1] of sortedEdges(edgeSet.values())) {
        if (sign[v0] * sign[v1] < 0) {
            const denom = dotValues[v1] - dotValues[v0];
            const w0 = dotValues[v1] / denom;
            const w1 = -dotValues[v0] / denom;
            const eiVertex = add(mul(w0, polyhedron.vertices[v0]),
                mul(w1, polyhedron.vertices[v1]));
            const eiIndex = splitVertices.length;
            eiVMap.set(edgeKey(v0, v1), eiIndex);
            splitVertices.push(eiVertex);
        }
    }

    // The average point will be appended to splitVertices later when
    // necessary.
    return { splitVertices, eiVMap };
}

// Look up the edge-interior intersection vertex for the mesh edge (v0,v1).
function eiVertexIndex(eiVMap: Map<string, number>, v0: number,
    v1: number): number {
    const index = eiVMap.get(edgeKey(v0, v1));
    logAssert(index !== undefined,
        'Missing edge-interior intersection vertex.');
    return index as number;
}

function doSplit(polyhedron: ConvexMesh3, sign: readonly number[],
    eiVMap: Map<string, number>,
    wantPosMesh: boolean, posMesh: ConvexMesh3Triangle[],
    wantNegMesh: boolean, negMesh: ConvexMesh3Triangle[],
    posIntersection: Map<number, number>): void {
    // The port of std::map::insert, which does not overwrite an existing key.
    const insertPos = (key: number, value: number): void => {
        if (!posIntersection.has(key)) {
            posIntersection.set(key, value);
        }
    };

    for (const triangle of polyhedron.triangles) {
        const v0 = triangle[0];
        const v1 = triangle[1];
        const v2 = triangle[2];
        let v01 = -1;
        let v12 = -1;
        let v20 = -1;

        if (sign[v0] > 0) {
            if (sign[v1] > 0) {
                if (sign[v2] > 0) {
                    // +++
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                }
                else if (sign[v2] < 0) {
                    // ++-
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v0, v12, v20]);
                        posMesh.push([v0, v1, v12]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v2, v20, v12]);
                    }
                    insertPos(v20, v12);
                }
                else {
                    // ++0
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                }
            }
            else if (sign[v1] < 0) {
                if (sign[v2] > 0) {
                    // +-+
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    if (wantPosMesh) {
                        posMesh.push([v0, v01, v12]);
                        posMesh.push([v0, v12, v2]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v1, v12, v01]);
                    }
                    insertPos(v12, v01);
                }
                else if (sign[v2] < 0) {
                    // +--
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v0, v01, v20]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v1, v20, v01]);
                        negMesh.push([v1, v2, v20]);
                    }
                    insertPos(v20, v01);
                }
                else {
                    // +-0
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    if (wantPosMesh) {
                        posMesh.push([v2, v0, v01]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v2, v01, v1]);
                    }
                    insertPos(v2, v01);
                }
            }
            else {
                if (sign[v2] > 0) {
                    // +0+
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                }
                else if (sign[v2] < 0) {
                    // +0-
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v1, v20, v0]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v1, v2, v20]);
                    }
                    insertPos(v20, v1);
                }
                else {
                    // +00
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                    insertPos(v2, v1);
                }
            }
        }
        else if (sign[v0] < 0) {
            if (sign[v1] > 0) {
                if (sign[v2] > 0) {
                    // -++
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v1, v20, v01]);
                        posMesh.push([v1, v2, v20]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v0, v01, v20]);
                    }
                    insertPos(v01, v20);
                }
                else if (sign[v2] < 0) {
                    // -+-
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    if (wantPosMesh) {
                        posMesh.push([v1, v12, v01]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v0, v01, v12]);
                        negMesh.push([v0, v12, v2]);
                    }
                    insertPos(v01, v12);
                }
                else {
                    // -+0
                    v01 = eiVertexIndex(eiVMap, v0, v1);
                    if (wantPosMesh) {
                        posMesh.push([v1, v2, v01]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v2, v0, v01]);
                    }
                    insertPos(v01, v2);
                }
            }
            else if (sign[v1] < 0) {
                if (sign[v2] > 0) {
                    // --+
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v2, v20, v12]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v12]);
                        negMesh.push([v0, v12, v20]);
                    }
                    insertPos(v12, v20);
                }
                else if (sign[v2] < 0) {
                    // ---
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
                else {
                    // --0
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
            }
            else {
                if (sign[v2] > 0) {
                    // -0+
                    v20 = eiVertexIndex(eiVMap, v2, v0);
                    if (wantPosMesh) {
                        posMesh.push([v2, v20, v1]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v20]);
                    }
                    insertPos(v1, v20);
                }
                else if (sign[v2] < 0) {
                    // -0-
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
                else {
                    // -00
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
            }
        }
        else {
            if (sign[v1] > 0) {
                if (sign[v2] > 0) {
                    // 0++
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                }
                else if (sign[v2] < 0) {
                    // 0+-
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    if (wantPosMesh) {
                        posMesh.push([v1, v12, v0]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v2, v0, v12]);
                    }
                    insertPos(v0, v12);
                }
                else {
                    // 0+0
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                    insertPos(v0, v2);
                }
            }
            else if (sign[v1] < 0) {
                if (sign[v2] > 0) {
                    // 0-+
                    v12 = eiVertexIndex(eiVMap, v1, v2);
                    if (wantPosMesh) {
                        posMesh.push([v2, v0, v12]);
                    }
                    if (wantNegMesh) {
                        negMesh.push([v1, v12, v0]);
                    }
                    insertPos(v12, v0);
                }
                else if (sign[v2] < 0) {
                    // 0--
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
                else {
                    // 0-0
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
            }
            else {
                if (sign[v2] > 0) {
                    // 00+
                    if (wantPosMesh) {
                        posMesh.push([v0, v1, v2]);
                    }
                    insertPos(v1, v0);
                }
                else if (sign[v2] < 0) {
                    // 00-
                    if (wantNegMesh) {
                        negMesh.push([v0, v1, v2]);
                    }
                }
                else {
                    // 000
                    // This case cannot occur with exact arithmetic, because it
                    // would have been trapped previously by the tests
                    // numPositive == 0 or numNegative == 0.
                    logError('This case cannot occur with exact arithmetic.');
                }
            }
        }
    }
}

function getIntersectionPolygonSplit(
    posIntersection: ReadonlyMap<number, number>,
    splitVertices: readonly ConvexMesh3Vertex[], wantIntrPolygon: boolean,
    result: IntrConvexMesh3Plane3FIResult): number[] {
    const numVertices = posIntersection.size;
    const polygon: number[] = new Array<number>(numVertices);

    // The upstream traversal starts at posIntersection.begin(), the smallest
    // key of the std::map, and follows the successor chain.
    if (numVertices > 0) {
        let current = Math.min(...posIntersection.keys());
        for (let i = 0; i < numVertices; ++i) {
            polygon[i] = current;
            const next = posIntersection.get(current);
            logAssert(next !== undefined,
                'The polygon of intersection is not a closed polygon.');
            current = next as number;
        }
    }

    if (wantIntrPolygon) {
        result.intersectionPolygon = polygon.map(
            index => splitVertices[index].clone());
    }

    return polygon;
}

function getSplitPolyhedra(splitVertices: ConvexMesh3Vertex[],
    polygon: readonly number[], wantIntrMesh: boolean, wantPosMesh: boolean,
    posMesh: ConvexMesh3Triangle[], wantNegMesh: boolean,
    negMesh: ConvexMesh3Triangle[],
    result: IntrConvexMesh3Plane3FIResult): void {
    // Triangulate the polygon for use by the positive polyhedron. A triangle
    // fan will not always work when the polygon has collinear vertices. The
    // average of the polygon vertices is inserted as an extra vertex. The
    // triangulation includes each triangle that is formed by the average
    // point and an edge of the polygon. The negative polyhedron uses the same
    // triangulation but with opposite chirality. Upstream NOTE: to avoid
    // biases in the average due to vertex distribution, use the center of
    // mass of the polygon instead.
    const numVertices = polygon.length;
    logAssert(numVertices > 0,
        'The transverse split has an empty polygon of intersection.');
    let average = Vector.zero(3);
    for (const i of polygon) {
        average = add(average, splitVertices[i]);
    }
    for (let d = 0; d < 3; ++d) {
        average.values[d] /= numVertices;
    }
    const iAvrIndex = splitVertices.length;
    splitVertices.push(average);

    const intrMesh: ConvexMesh3Triangle[] = [];
    for (let i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++) {
        if (wantPosMesh) {
            posMesh.push([iAvrIndex, polygon[i0], polygon[i1]]);
        }

        if (wantNegMesh) {
            negMesh.push([iAvrIndex, polygon[i1], polygon[i0]]);
        }

        if (wantIntrMesh) {
            intrMesh.push([iAvrIndex, polygon[i0], polygon[i1]]);
        }
    }

    const uvt = new UniqueVerticesSimplices<ConvexMesh3Vertex>(3);

    if (wantPosMesh) {
        result.positivePolyhedron.configuration = ConvexMesh3.CFG_POLYHEDRON;
        const out = uvt.removeDuplicateAndUnusedVerticesGrouped(splitVertices,
            posMesh);
        result.positivePolyhedron.vertices = out.vertices;
        result.positivePolyhedron.triangles = out.simplices.map(triangle =>
            [triangle[0], triangle[1], triangle[2]] as ConvexMesh3Triangle);
    }

    if (wantNegMesh) {
        result.negativePolyhedron.configuration = ConvexMesh3.CFG_POLYHEDRON;
        const out = uvt.removeDuplicateAndUnusedVerticesGrouped(splitVertices,
            negMesh);
        result.negativePolyhedron.vertices = out.vertices;
        result.negativePolyhedron.triangles = out.simplices.map(triangle =>
            [triangle[0], triangle[1], triangle[2]] as ConvexMesh3Triangle);
    }

    if (wantIntrMesh) {
        result.intersectionMesh.configuration = ConvexMesh3.CFG_POLYGON;
        const out = uvt.removeDuplicateAndUnusedVerticesGrouped(splitVertices,
            intrMesh);
        result.intersectionMesh.vertices = out.vertices;
        result.intersectionMesh.triangles = out.simplices.map(triangle =>
            [triangle[0], triangle[1], triangle[2]] as ConvexMesh3Triangle);
    }
}
