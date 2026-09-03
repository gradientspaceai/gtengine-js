// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IncrementalDelaunay2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Incremental insertion and removal of vertices in a Delaunay triangulation.
// The triangles are counterclockwise ordered. The insertion code is that of
// Delaunay2 (src/Delaunay2.ts); upstream copy-pasted it into this header and
// asks that any change be made in both places, so the port keeps the two
// implementations textually parallel as well.
//
// The removal code is an implementation of the algorithm in
//     Olivier Devillers,
//     "On Deletion in Delaunay Triangulations",
//     International Journal of Computational Geometry and Applications,
//     World Scientific Publishing, 2002, 12, pp. 193-205.
//     https://hal.inria.fr/inria-00167201/document
// The weight function for the priority queue, implemented as a min-heap, is
// the negative of the power(p, circle(q0,q1,q2)) function described in the
// paper.
//
// The paper appears to assume that the removal point is an interior point of
// the triangulation. Just as the insertion algorithms are different for
// interior points and for boundary points, the removal algorithms are
// different for interior points and for boundary points.
//
// The paper mentions that degeneracies (collinear points, cocircular points)
// are handled by jittering. Although one hopes that jittering prevents
// degeneracies -- and perhaps probabilistically this is acceptable, the only
// guarantee for a correct result is to use exact arithmetic on the input
// points. The implementation here uses a blend of interval and rational
// arithmetic for exactness; the input points are not jittered.
//
// The details of the algorithms and implementation are provided in
// https://www.geometrictools.com/Documentation/IncrementalDelaunayTriangulation.pdf
//
// Port notes:
// * Following the Delaunay2 precedent, upstream's InputRational =
//   BSNumber<UIntegerFP32<2 or 4>>, ComputeRational =
//   BSNumber<UIntegerFP32<36 or 264>>, the preallocated mCRPool of scratch
//   compute-rationals and the Copy() widening helper are all dropped. The
//   port's BSNumber is bigint-backed and grows as needed, so a single
//   rational type suffices and the exact-arithmetic results are identical.
// * std::set<EdgeKey<true>> becomes DirectedEdgeKeySet, a Map iterated in
//   sorted key order. std::map<int32_t,int32_t> and std::map<size_t,size_t>
//   likewise become Maps whose observable iteration is in sorted key order.
//   std::unordered_set<Triangle*> containers are read through the port's
//   sorted VETManifoldMesh accessors.
// * std::map<Vector2<T>, size_t> mVertexIndexMap is only ever searched and
//   erased, never iterated, so it becomes a Map keyed by the string form of
//   the coordinate pair. C++ compares -0.0 and 0.0 as equivalent keys and so
//   does the string form, since String(-0) is "0".
// * size_t 'invalid' (std::numeric_limits<size_t>::max()) becomes -1.
// * Output-reference parameters become returned objects: GetTriangulation()
//   returns { vertices, triangles }, GetTriangle(t)/GetAdjacent(t) return the
//   tuple or null, and GetHull() returns the index array.
// * The unused private ToTriangle() helper (dead code upstream) is not
//   ported.

import { BSNumber } from './BSNumber.js';
import { ETManifoldMesh, ETManifoldMeshTriangle } from './ETManifoldMesh.js';
import { logAssert, logError } from './Logger.js';
import { MinHeap } from './MinHeap.js';
import type { MinHeapRecord } from './MinHeap.js';
import { SWInterval } from './SWInterval.js';
import { Vector } from './Vector.js';
import { VETManifoldMesh } from './VETManifoldMesh.js';

// A rational 2D point, the port of Vector2<InputRational>.
type RationalPoint2 = [BSNumber, BSNumber];

// The sentinel for "no such index"; the port of the upstream size_t constant
// 'invalid'. It also selects the stored query point in toLine().
const invalid = -1;

// Indexing for the vertices of the triangle adjacent to a vertex. The edge
// adjacent to vertex j is <mIndex[j][0], mIndex[j][1]> and is listed so that
// the triangle interior is to your left as you walk around the edges.
const mIndex: readonly (readonly [number, number])[] =
    [[0, 1], [1, 2], [2, 0]];

// The port of std::set<EdgeKey<true>>; iterated in increasing key order.
class DirectedEdgeKeySet {
    private mMap: Map<string, [number, number]>;

    constructor() {
        this.mMap = new Map<string, [number, number]>();
    }

    insert(v0: number, v1: number): void {
        this.mMap.set(`${v0},${v1}`, [v0, v1]);
    }

    get size(): number {
        return this.mMap.size;
    }

    keys(): [number, number][] {
        const keys = Array.from(this.mMap.values());
        keys.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]));
        return keys;
    }
}

// The port of the SearchInfo struct used by getContainingTriangle(). The
// first triangle searched is 'initialTriangle'. On return, 'path' stores the
// (ordered) triangle indices visited during the search, with 'numPath' valid
// entries. The last visited triangle has index 'finalTriangle' and vertex
// indices finalV[0,1,2] stored in counterclockwise order. The last edge of
// the search is <finalV[0],finalV[1]>. For spatially coherent query points in
// numerous calls, specify 'finalTriangle' of the previous call as
// 'initialTriangle' of the next call to reduce search times.
export class IncrementalDelaunay2SearchInfo {
    initialTriangle: number;
    finalTriangle: number;
    finalV: [number, number, number];
    numPath: number;
    path: number[];

    constructor() {
        this.initialTriangle = invalid;
        this.finalTriangle = invalid;
        this.finalV = [invalid, invalid, invalid];
        this.numPath = 0;
        this.path = [];
    }
}

// The weight of a removal-polygon vertex. See the class comments below.
enum RPWeightType {
    finite,
    infinite,
    unmodifiable
}

// Let Vc be a vertex in the removal polygon. If Vc is not an ear, its weight
// is +infinity. If Vc is an ear, let Vp be its predecessor and let Vn be its
// successor when traversing counterclockwise. Let P be the removal point. The
// weight is
//   W = -H(Vp, Vc, Vn, P) / D(Vp, Vc, Vn) = WNumer / WDenom
// where
//           +-                -+
//           | Vp.x  Vc.x  Vn.x |       +-                    -+
//   D = det | Vp.y  Vc.y  Vn.y | = det | Vc.x-Vp.x  Vn.x-Vp.x |
//           |   1     1     1  |       | Vc.y-Vp.y  Vn.y-Vp.y |
//           +-                -+       +-                    -+
// and
//           +-                             -+
//           | Vp.x    Vc.x    Vn.x    P.x   |
//   H = det | Vp.y    Vc.y    Vn.y    P.y   |
//           | |Vp|^2  |Vc|^2  |Vn|^2  |P|^2 |
//           |   1       1       1       1   |
//           +-                             -+
//            +-                                          -+
//            | Vc.x-Vp.x      Vn.x-Vp.x      P.x-Vp.x     |
//     = -det | Vc.y-Vp.y      Vn.y-Vp.y      P.y-Vp.y     |
//            | |Vc|^2-|Vp|^2  |Vn|^2-|Vp|^2  |P|^2-|Vp|^2 |
//            +-                                          -+
//
// To use BSNumber-based rationals, the weight is a ratio stored as a pair
// (numerator, denominator) with denominator > 0. The comparison of weights is
// WN0/WD0 < WN1/WD1, implemented as WN0*WD1 < WN1*WD0. Additionally, a type
// flag distinguishes a finite ratio (for convex vertices) from a weight that
// is infinite (for reflex vertices) and from a vertex that must not be
// clipped at all.
class RPWeight {
    // The finite weight is numerator/denominator with a nonnegative numerator
    // and a positive denominator. If the weight is not finite, the numerator
    // and denominator values are invalid.
    numerator: BSNumber;
    denominator: BSNumber;
    type: RPWeightType;

    constructor(inType: RPWeightType = RPWeightType.unmodifiable) {
        this.numerator = BSNumber.fromNumber(0);
        this.denominator = BSNumber.fromNumber(
            inType === RPWeightType.finite ? 1 : 0);
        this.type = inType;
    }

    // The port of operator<: finite < infinite < unmodifiable.
    static lessThan(lhs: RPWeight, rhs: RPWeight): boolean {
        if (lhs.type === RPWeightType.finite) {
            if (rhs.type === RPWeightType.finite) {
                return lhs.numerator.mul(rhs.denominator).lessThan(
                    rhs.numerator.mul(lhs.denominator));
            }
            else {
                // rhs.type is infinite or unmodifiable
                return true;
            }
        }
        else if (lhs.type === RPWeightType.infinite) {
            if (rhs.type === RPWeightType.finite) {
                return false;
            }
            else {
                // rhs.type is infinite or unmodifiable
                return rhs.type === RPWeightType.unmodifiable;
            }
        }
        else {
            // lhs.type is unmodifiable
            return false;
        }
    }
}

// The port of the private RPVertex class.
class RPVertex {
    // The index relative to IncrementalDelaunay2 mVertices[].
    vIndex: number;

    // A vertex is either convex or reflex.
    isConvex: boolean;

    // Vertex indices for the polygon. These are indices relative to
    // RPPolygon mVertices[].
    iPrev: number;
    iNext: number;

    // Support for the priority queue of ears.
    record: MinHeapRecord<number, RPWeight> | null;

    constructor() {
        this.vIndex = invalid;
        this.isConvex = false;
        this.iPrev = invalid;
        this.iNext = invalid;
        this.record = null;
    }
}

type ToLineFunction = (vPrev: number, vCurr: number, vNext: number) => number;

// The port of the private RPPolygon class, a circular list of the removal
// polygon vertices that supports dynamic removal.
class RPPolygon {
    private mNumActive: number;
    private mVertices: RPVertex[];

    constructor(polygon: readonly number[], toLine: ToLineFunction) {
        this.mNumActive = polygon.length;
        this.mVertices = new Array<RPVertex>(polygon.length);

        // Create a circular list of the polygon vertices for dynamic removal
        // of vertices.
        const numVertices = polygon.length;
        for (let i = 0; i < numVertices; ++i) {
            const vertex = new RPVertex();
            vertex.vIndex = polygon[i];
            vertex.iPrev = (i > 0 ? i - 1 : numVertices - 1);
            vertex.iNext = (i < numVertices - 1 ? i + 1 : 0);
            this.mVertices[i] = vertex;
        }

        // Determine whether each vertex is convex or reflex.
        for (let i = 0; i < numVertices; ++i) {
            const tri = this.getTriangle(i);
            this.mVertices[i].isConvex = (toLine(tri.vPrev, tri.vCurr, tri.vNext) < 0);
        }
    }

    vertex(i: number): RPVertex {
        return this.mVertices[i];
    }

    getTriangle(i: number): { vPrev: number; vCurr: number; vNext: number } {
        const vertex = this.mVertices[i];
        return {
            vPrev: this.mVertices[vertex.iPrev].vIndex,
            vCurr: vertex.vIndex,
            vNext: this.mVertices[vertex.iNext].vIndex
        };
    }

    classify(i: number, toLine: ToLineFunction): void {
        const tri = this.getTriangle(i);
        this.mVertices[i].isConvex = (toLine(tri.vPrev, tri.vCurr, tri.vNext) < 0);
    }

    getNumActive(): number {
        return this.mNumActive;
    }

    getActive(): number {
        for (let i = 0; i < this.mVertices.length; ++i) {
            if (this.mVertices[i].iPrev !== invalid) {
                return i;
            }
        }

        logError('Expecting to find an active vertex.');
        return invalid;
    }

    remove(i: number): void {
        // Remove the vertex from the polygon.
        const vertex = this.mVertices[i];
        const iPrev = vertex.iPrev;
        const iNext = vertex.iNext;
        this.mVertices[iPrev].iNext = iNext;
        this.mVertices[iNext].iPrev = iPrev;

        vertex.vIndex = invalid;
        vertex.isConvex = false;
        vertex.iPrev = invalid;
        vertex.iNext = invalid;
        vertex.record = null;

        --this.mNumActive;
    }
}

export class IncrementalDelaunay2 {
    // The rectangular domain in which all input points live.
    private mXMin: number;
    private mYMin: number;
    private mXMax: number;
    private mYMax: number;

    // The rectangular domain is always inserted into the triangulation first.
    // After all your insert and remove calls, if you remove the rectangle via
    // finalizeTriangulation(), you can no longer insert or remove points. The
    // values of this member are
    //   0: rectangle has not been removed
    //   1: rectangle is in the process of being removed
    //   2: rectangle is removed
    // The 3-valued member allows an exception to be thrown by remove() when
    // the state is 2, while finalization can use remove() without exceptions
    // when the state is 1.
    private mRectangleRemoved: number;

    // The current vertices.
    private mVertexIndexMap: Map<string, number>;
    private mVertices: Vector[];
    private mIRVertices: RationalPoint2[];

    // The graph is the current triangulation.
    private mGraph: VETManifoldMesh;

    // Support for queries associated with the mesh of Delaunay triangles.
    private mTriangles: [number, number, number][];
    private mAdjacencies: [number, number, number][];
    private mTrianglesAndAdjacenciesNeedUpdate: boolean;
    private mQueryPoint: Vector;
    private mIRQueryPoint: RationalPoint2;

    // Wrap the toLine function for use in retriangulating the removal
    // polygon.
    private mToLineWrapper: ToLineFunction;

    // A bounding rectangle for the input points must be specified. NOTE: The
    // bounding rectangle is inserted automatically into the triangulation as
    // two triangles. Once you are finished inserting and removing points,
    // call finalizeTriangulation(). After that call, you cannot insert or
    // remove points.
    constructor(xMin: number, yMin: number, xMax: number, yMax: number) {
        this.mXMin = xMin;
        this.mYMin = yMin;
        this.mXMax = xMax;
        this.mYMax = yMax;
        this.mRectangleRemoved = 0;
        this.mVertexIndexMap = new Map<string, number>();
        this.mVertices = [];
        this.mIRVertices = [];
        this.mGraph = new VETManifoldMesh();
        this.mTriangles = [];
        this.mAdjacencies = [];
        this.mTrianglesAndAdjacenciesNeedUpdate = true;
        this.mQueryPoint = new Vector(2);
        this.mIRQueryPoint = [BSNumber.fromNumber(0), BSNumber.fromNumber(0)];
        this.mToLineWrapper = (vPrev: number, vCurr: number, vNext: number) =>
            this.toLine(vPrev, vCurr, vNext);

        logAssert(this.mXMin < this.mXMax && this.mYMin < this.mYMax,
            'Invalid bounding rectangle.');

        // Create the vertices for a supertriangle that contains the input
        // rectangle
        //   V[0] = (x0,y0) = (xmin - dx, ymin - dy)
        //   V[1] = (x1,y1) = (xmin + 5 * dx, ymin - dy)
        //   V[2] = (x2,y2) = (xmin - dx, ymin + 5 * dy)
        // Create the vertices for the input rectangle
        //   V[3] = (x3,y3) = (xmin, ymin)
        //   V[4] = (x4,y4) = (xmax, ymin)
        //   V[5] = (x5,y5) = (xmin, ymax)
        //   V[6] = (x6,y6) = (xmax, ymax)
        const xDelta = this.mXMax - this.mXMin;
        const yDelta = this.mYMax - this.mYMin;
        const x0 = this.mXMin - xDelta;
        const y0 = this.mYMin - yDelta;
        const x1 = this.mXMin + 5 * xDelta;
        const y1 = y0;
        const x2 = x0;
        const y2 = this.mYMin + 5 * yDelta;
        const vertices: readonly (readonly [number, number])[] = [
            [x0, y0],
            [x1, y1],
            [x2, y2],
            [this.mXMin, this.mYMin],
            [this.mXMax, this.mYMin],
            [this.mXMin, this.mYMax],
            [this.mXMax, this.mYMax]
        ];

        // Insert the vertices into the vertex storage.
        for (let i = 0; i < vertices.length; ++i) {
            this.addVertex(vertices[i][0], vertices[i][1], i);
        }

        // Create the triangles formed by the supervertices and the input
        // rectangle vertices.
        const triangles: readonly (readonly [number, number, number])[] = [
            [0, 5, 2], [0, 3, 5], [0, 4, 3], [0, 1, 4], [1, 6, 4],
            [1, 2, 6], [2, 5, 6], [3, 4, 6], [3, 6, 5]
        ];

        // Insert the triangles into the triangulation.
        for (const tri of triangles) {
            const inserted = this.mGraph.insert(tri[0], tri[1], tri[2]);
            logAssert(inserted !== null, 'Failed to insert initial triangle.');
        }
    }

    // Insert a point into the triangulation. It is required that the point be
    // strictly inside the input rectangle; if it is not, an exception is
    // thrown. If the input point already exists, its vertex map index is
    // returned; otherwise, the point is inserted into the vertex map and the
    // index associated with the insertion is returned. When the input
    // rectangle has been removed, no insertion occurs and 'invalid' (-1) is
    // returned.
    insert(position: Vector): number {
        logAssert(position.size === 2, 'IncrementalDelaunay2 requires 2D points.');
        logAssert(
            this.mXMin < position.get(0) && position.get(0) < this.mXMax &&
            this.mYMin < position.get(1) && position.get(1) < this.mYMax,
            'The position must be strictly inside the domain specified in the constructor.');

        if (this.mRectangleRemoved === 2) {
            // You cannot insert points after the input rectangle is removed.
            return invalid;
        }

        this.mTrianglesAndAdjacenciesNeedUpdate = true;

        const key = IncrementalDelaunay2.vertexKey(position.get(0), position.get(1));
        const existing = this.mVertexIndexMap.get(key);
        if (existing !== undefined) {
            // The vertex already exists.
            return existing;
        }

        // Store the position in the various pools.
        const posIndex = this.mVertices.length;
        this.addVertex(position.get(0), position.get(1), posIndex);

        this.update(posIndex);
        return posIndex;
    }

    // Remove a point from the triangulation. The return value is the index
    // associated with the vertex in the vertex map when that vertex exists.
    // If the vertex does not exist, the return value is 'invalid' (-1).
    remove(position: Vector): number {
        logAssert(position.size === 2, 'IncrementalDelaunay2 requires 2D points.');
        if (this.mRectangleRemoved === 0) {
            logAssert(
                this.mXMin < position.get(0) && position.get(0) < this.mXMax &&
                this.mYMin < position.get(1) && position.get(1) < this.mYMax,
                'The position must be strictly inside the domain specified in the constructor.');
        }

        if (this.mRectangleRemoved === 2) {
            // You cannot remove points after the input rectangle is removed.
            return invalid;
        }

        this.mTrianglesAndAdjacenciesNeedUpdate = true;

        const key = IncrementalDelaunay2.vertexKey(position.get(0), position.get(1));
        const found = this.mVertexIndexMap.get(key);
        if (found === undefined) {
            // The position is not a vertex of the triangulation.
            return invalid;
        }
        const vRemovalIndex = found;

        if (this.mVertexIndexMap.size === 4) {
            // The last vertex of the input rectangle is to be removed.
            for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
                const removed = this.mGraph.remove(vRemovalIndex, i0, i1);
                logAssert(removed, 'Unexpected removal failure.');
            }

            const inserted = this.mGraph.insert(0, 1, 2);
            logAssert(inserted !== null, 'Failed to insert supertriangle.');

            this.mVertexIndexMap.delete(key);
            return vRemovalIndex;
        }

        // Locate the position in the vertices of the graph.
        const vertex = this.mGraph.getVertex(vRemovalIndex);
        logAssert(vertex !== null,
            'Expecting to find the to-be-removed vertex in the triangulation.');

        let removalPointOnBoundary = false;
        for (const vIndex of vertex.getVAdjacent()) {
            if (IncrementalDelaunay2.isSupervertex(vIndex)) {
                // The triangle has a supervertex, so the removal point is on
                // the boundary of the Delaunay triangulation.
                removalPointOnBoundary = true;
                break;
            }
        }

        const adjacents = vertex.getTAdjacent();
        const polygon = this.deleteRemovalPolygon(vRemovalIndex, adjacents);

        if (removalPointOnBoundary) {
            this.retriangulateBoundaryRemovalPolygon(vRemovalIndex, polygon);
        }
        else {
            this.retriangulateInteriorRemovalPolygon(vRemovalIndex, polygon);
        }

        this.mVertexIndexMap.delete(key);
        return vRemovalIndex;
    }

    // Call this only after you are finished inserting points into or removing
    // points from the triangulation.
    finalizeTriangulation(): boolean {
        if (this.mRectangleRemoved === 2) {
            // You cannot remove the input rectangle more than once.
            return false;
        }

        // Remove the input rectangle vertices. The triangles strictly
        // interior to the input rectangle form the Delaunay triangulation.
        // However, the triangles sharing a supervertex still exist in the
        // graph.
        const vertex: readonly (readonly [number, number])[] = [
            [this.mXMin, this.mYMin],
            [this.mXMax, this.mYMin],
            [this.mXMin, this.mYMax],
            [this.mXMax, this.mYMax]
        ];

        this.mRectangleRemoved = 1;
        for (let i = 0; i < vertex.length; ++i) {
            const index = this.remove(Vector.fromArray([vertex[i][0], vertex[i][1]]));
            logAssert(index === i + 3, 'Incorrect index for vertex.');
        }
        this.mRectangleRemoved = 2;
        return true;
    }

    // Get the current triangulation including the supervertices and the
    // triangles containing supervertices.
    getTriangulation(): { vertices: Vector[]; triangles: [number, number, number][] } {
        const vertices = this.mVertices.map(v => v.clone());
        const triangles: [number, number, number][] = [];
        for (const tKey of this.mGraph.getTriangleKeys()) {
            triangles.push([tKey.V[0], tKey.V[1], tKey.V[2]]);
        }
        return { vertices: vertices, triangles: triangles };
    }

    // Get the current graph, which includes all triangles whether Delaunay or
    // those containing a supervertex.
    getGraph(): VETManifoldMesh {
        return this.mGraph;
    }

    // Queries associated with the mesh of Delaunay triangles. The triangles
    // containing a supervertex are not included in the triangle queries.
    //
    // Note that getNumVertices()/getVertices() report the full vertex pool,
    // which always begins with the 3 supervertices followed by the 4 input
    // rectangle vertices; the triangle indices are relative to that pool.
    // Upstream documents these accessors as excluding the supervertices,
    // which is inaccurate, but the indices stored in getTriangles() are pool
    // indices, so the behavior (not the comment) is the correct one.
    getNumVertices(): number {
        return this.mVertices.length;
    }

    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getNumTriangles(): number {
        this.updateTrianglesAndAdjacenciesIfNeeded();
        return this.mTriangles.length;
    }

    getTriangles(): readonly [number, number, number][] {
        this.updateTrianglesAndAdjacenciesIfNeeded();
        return this.mTriangles;
    }

    getAdjacencies(): readonly [number, number, number][] {
        this.updateTrianglesAndAdjacenciesIfNeeded();
        return this.mAdjacencies;
    }

    // Get the vertex indices for triangle t. The return value is null when t
    // is not a valid triangle index.
    getTriangle(t: number): [number, number, number] | null {
        this.updateTrianglesAndAdjacenciesIfNeeded();
        if (0 <= t && t < this.mTriangles.length) {
            const tri = this.mTriangles[t];
            return [tri[0], tri[1], tri[2]];
        }
        return null;
    }

    // Get the indices for triangles adjacent to triangle t. The return value
    // is null when t is not a valid triangle index. When valid, triangle t
    // has ordered vertices <V[0], V[1], V[2]>. The value adjacent[0] is the
    // index for the triangle adjacent to edge <V[0], V[1]>, adjacent[1] is
    // the index for the triangle adjacent to edge <V[1], V[2]>, and
    // adjacent[2] is the index for the triangle adjacent to edge <V[2], V[0]>.
    getAdjacent(t: number): [number, number, number] | null {
        this.updateTrianglesAndAdjacenciesIfNeeded();
        if (0 <= t && t < this.mAdjacencies.length) {
            const adj = this.mAdjacencies[t];
            return [adj[0], adj[1], adj[2]];
        }
        return null;
    }

    // Get the convex polygon that is the hull of the Delaunay triangles. The
    // polygon is counterclockwise ordered with vertices V[hull[0]],
    // V[hull[1]], ..., V[hull[hull.length-1]].
    getHull(): number[] {
        this.updateTrianglesAndAdjacenciesIfNeeded();

        // The hull edges are shared by the triangles with exactly one
        // supervertex.
        const edges = new Map<number, number>();
        for (let v = 0; v < 3; ++v) {
            const vertex = this.mGraph.getVertex(v);
            logAssert(vertex !== null,
                'Expecting the supervertices to exist in the graph.');

            for (const adj of vertex.getTAdjacent()) {
                for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2, ++i2) {
                    if (adj.V[i0] === v) {
                        if (IncrementalDelaunay2.isDelaunayVertex(adj.V[i1]) &&
                            IncrementalDelaunay2.isDelaunayVertex(adj.V[i2])) {
                            // std::map::insert does not overwrite an existing
                            // key, so neither does this.
                            if (!edges.has(adj.V[i2])) {
                                edges.set(adj.V[i2], adj.V[i1]);
                            }
                            break;
                        }
                    }
                }
            }
        }

        // Upstream dereferences edges.begin() unconditionally, which is
        // undefined behavior when the Delaunay triangulation is empty (for
        // example, before any point is inserted). The port returns an empty
        // hull instead.
        if (edges.size === 0) {
            return [];
        }

        // Repackage the edges into a convex polygon with vertices ordered
        // counterclockwise. The upstream std::map is traversed from its
        // smallest key, so the port starts at the smallest key as well.
        const numEdges = edges.size;
        const hull = new Array<number>(numEdges);
        const vStart = Math.min(...Array.from(edges.keys()));
        let vNext = edges.get(vStart) as number;
        let i = 0;
        hull[0] = vStart;
        while (vNext !== vStart) {
            // Upstream bug (fixed in the port): the traversal is unbounded. If
            // the edge map is not a single cycle through the smallest key, the
            // upstream loop runs forever and writes past the end of 'hull'.
            // That happens for a degenerate triangulation: for collinear input
            // points there are no Delaunay triangles at all, yet the triangles
            // sharing a supervertex still contribute edges, and those edges
            // form a path with a 2-cycle at its far end rather than a closed
            // polygon. The port detects the condition and reports it.
            logAssert(i + 1 < numEdges,
                'The Delaunay triangulation is degenerate (the input points are '
                + 'collinear), so it has no convex hull polygon.');
            hull[++i] = vNext;
            const next = edges.get(vNext);
            logAssert(next !== undefined, 'Expecting to find a hull edge.');
            vNext = next;
        }
        return hull;
    }

    // Support for searching the Delaunay triangles that contain the point p.
    // If there is a containing triangle, the returned value is a triangle
    // index i with 0 <= i < getNumTriangles(). If there is not a containing
    // triangle, 'invalid' (-1) is returned. The computations are performed
    // using exact rational arithmetic.
    getContainingTriangle(p: Vector, info: IncrementalDelaunay2SearchInfo): number {
        logAssert(p.size === 2, 'IncrementalDelaunay2 requires 2D points.');
        this.updateTrianglesAndAdjacenciesIfNeeded();

        this.mQueryPoint = p.clone();
        this.mIRQueryPoint = [BSNumber.fromNumber(p.get(0)), BSNumber.fromNumber(p.get(1))];

        const numTriangles = this.mTriangles.length;
        info.path = new Array<number>(numTriangles).fill(invalid);
        info.numPath = 0;
        let tIndex: number;
        if (0 <= info.initialTriangle && info.initialTriangle < numTriangles) {
            tIndex = info.initialTriangle;
        }
        else {
            info.initialTriangle = 0;
            tIndex = 0;
        }

        for (let t = 0; t < numTriangles; ++t) {
            const v = this.mTriangles[tIndex];
            const adj = this.mAdjacencies[tIndex];

            info.finalTriangle = tIndex;
            info.finalV = [v[0], v[1], v[2]];
            info.path[info.numPath++] = tIndex;

            let i0 = 1, i1 = 2, i2 = 0;
            for (; i2 < 3; i0 = i1, i1 = i2++) {
                // toLine(pIndex, v0Index, v1Index) uses mQueryPoint when
                // pIndex is 'invalid'.
                if (this.toLine(invalid, v[i0], v[i1]) > 0) {
                    tIndex = adj[i0];
                    if (tIndex === invalid) {
                        info.finalV[0] = v[i0];
                        info.finalV[1] = v[i1];
                        info.finalV[2] = v[i2];
                        return invalid;
                    }
                    break;
                }
            }
            if (i2 === 3) {
                return tIndex;
            }
        }
        return invalid;
    }

    // The string key for the map of vertex positions to indices.
    private static vertexKey(x: number, y: number): string {
        return `${x},${y}`;
    }

    private addVertex(x: number, y: number, index: number): void {
        this.mVertexIndexMap.set(IncrementalDelaunay2.vertexKey(x, y), index);
        this.mVertices.push(Vector.fromArray([x, y]));
        this.mIRVertices.push([BSNumber.fromNumber(x), BSNumber.fromNumber(y)]);
    }

    // Vertices 0, 1 and 2 are the supervertices of the enclosing
    // supertriangle; every other vertex belongs to the Delaunay
    // triangulation.
    private static isDelaunayVertex(vIndex: number): boolean {
        return vIndex >= 3;
    }

    private static isDelaunayTriangle(v0: number, v1: number, v2: number): boolean {
        return v0 >= 3 && v1 >= 3 && v2 >= 3;
    }

    private static isSupervertex(vIndex: number): boolean {
        return vIndex < 3;
    }

    // The floating-point query point for the given index, where 'invalid'
    // selects the stored query point.
    private inPoint(index: number): readonly number[] {
        return (index !== invalid ? this.mVertices[index].values : this.mQueryPoint.values);
    }

    // The rational query point for the given index, where 'invalid' selects
    // the stored query point.
    private irPoint(index: number): RationalPoint2 {
        return (index !== invalid ? this.mIRVertices[index] : this.mIRQueryPoint);
    }

    // Given a line with origin V0 and direction <V0,V1> and a query point P,
    // toLine returns
    //   +1, P on right of line
    //   -1, P on left of line
    //    0, P on the line
    private toLine(pIndex: number, v0Index: number, v1Index: number): number {
        // The expression tree has 13 nodes consisting of 6 input leaves and
        // 7 compute nodes.

        // Use interval arithmetic to determine the sign if possible.
        const inP = this.inPoint(pIndex);
        const inV0 = this.mVertices[v0Index].values;
        const inV1 = this.mVertices[v1Index].values;

        const x0 = SWInterval.sub(inP[0], inV0[0]);
        const y0 = SWInterval.sub(inP[1], inV0[1]);
        const x1 = SWInterval.sub(inV1[0], inV0[0]);
        const y1 = SWInterval.sub(inV1[1], inV0[1]);
        const x0y1 = x0.mul(y1);
        const x1y0 = x1.mul(y0);
        const det = x0y1.sub(x1y0);

        const zero = 0;
        if (det.get(0) > zero) {
            return +1;
        }
        else if (det.get(1) < zero) {
            return -1;
        }

        // The exact sign of the determinant is not known, so compute the
        // determinant using rational arithmetic.
        const irP = this.irPoint(pIndex);
        const irV0 = this.mIRVertices[v0Index];
        const irV1 = this.mIRVertices[v1Index];

        const crX0 = irP[0].sub(irV0[0]);
        const crY0 = irP[1].sub(irV0[1]);
        const crX1 = irV1[0].sub(irV0[0]);
        const crY1 = irV1[1].sub(irV0[1]);
        const crX0Y1 = crX0.mul(crY1);
        const crX1Y0 = crX1.mul(crY0);
        const crDet = crX0Y1.sub(crX1Y0);
        return crDet.getSign();
    }

    // For a triangle with counterclockwise vertices V0, V1 and V2 and a query
    // point P, toCircumcircle returns
    //   +1, P outside circumcircle of triangle
    //   -1, P inside circumcircle of triangle
    //    0, P on circumcircle of triangle
    private toCircumcircle(pIndex: number, v0Index: number, v1Index: number,
        v2Index: number): number {
        // The expression tree has 43 nodes consisting of 8 input leaves and
        // 35 compute nodes.

        // Use interval arithmetic to determine the sign if possible.
        const inP = this.inPoint(pIndex);
        const inV0 = this.mVertices[v0Index].values;
        const inV1 = this.mVertices[v1Index].values;
        const inV2 = this.mVertices[v2Index].values;

        const x0 = SWInterval.sub(inV0[0], inP[0]);
        const y0 = SWInterval.sub(inV0[1], inP[1]);
        const s00 = SWInterval.add(inV0[0], inP[0]);
        const s01 = SWInterval.add(inV0[1], inP[1]);
        const x1 = SWInterval.sub(inV1[0], inP[0]);
        const y1 = SWInterval.sub(inV1[1], inP[1]);
        const s10 = SWInterval.add(inV1[0], inP[0]);
        const s11 = SWInterval.add(inV1[1], inP[1]);
        const x2 = SWInterval.sub(inV2[0], inP[0]);
        const y2 = SWInterval.sub(inV2[1], inP[1]);
        const s20 = SWInterval.add(inV2[0], inP[0]);
        const s21 = SWInterval.add(inV2[1], inP[1]);
        const t00 = s00.mul(x0);
        const t01 = s01.mul(y0);
        const t10 = s10.mul(x1);
        const t11 = s11.mul(y1);
        const t20 = s20.mul(x2);
        const t21 = s21.mul(y2);
        const z0 = t00.add(t01);
        const z1 = t10.add(t11);
        const z2 = t20.add(t21);
        const y0z1 = y0.mul(z1);
        const y0z2 = y0.mul(z2);
        const y1z0 = y1.mul(z0);
        const y1z2 = y1.mul(z2);
        const y2z0 = y2.mul(z0);
        const y2z1 = y2.mul(z1);
        const c0 = y1z2.sub(y2z1);
        const c1 = y2z0.sub(y0z2);
        const c2 = y0z1.sub(y1z0);
        const x0c0 = x0.mul(c0);
        const x1c1 = x1.mul(c1);
        const x2c2 = x2.mul(c2);
        const det = x0c0.add(x1c1).add(x2c2);

        const zero = 0;
        if (det.get(0) > zero) {
            return -1;
        }
        else if (det.get(1) < zero) {
            return +1;
        }

        // The exact sign of the determinant is not known, so compute the
        // determinant using rational arithmetic.
        const irP = this.irPoint(pIndex);
        const irV0 = this.mIRVertices[v0Index];
        const irV1 = this.mIRVertices[v1Index];
        const irV2 = this.mIRVertices[v2Index];

        const crP0 = irP[0], crP1 = irP[1];
        const crV00 = irV0[0], crV01 = irV0[1];
        const crV10 = irV1[0], crV11 = irV1[1];
        const crV20 = irV2[0], crV21 = irV2[1];

        // Evaluate the expression tree.
        const crX0 = crV00.sub(crP0);
        const crY0 = crV01.sub(crP1);
        const crS00 = crV00.add(crP0);
        const crS01 = crV01.add(crP1);
        const crT00 = crS00.mul(crX0);
        const crT01 = crS01.mul(crY0);
        const crZ0 = crT00.add(crT01);

        const crX1 = crV10.sub(crP0);
        const crY1 = crV11.sub(crP1);
        const crS10 = crV10.add(crP0);
        const crS11 = crV11.add(crP1);
        const crT10 = crS10.mul(crX1);
        const crT11 = crS11.mul(crY1);
        const crZ1 = crT10.add(crT11);

        const crX2 = crV20.sub(crP0);
        const crY2 = crV21.sub(crP1);
        const crS20 = crV20.add(crP0);
        const crS21 = crV21.add(crP1);
        const crT20 = crS20.mul(crX2);
        const crT21 = crS21.mul(crY2);
        const crZ2 = crT20.add(crT21);

        const crY0Z1 = crY0.mul(crZ1);
        const crY0Z2 = crY0.mul(crZ2);
        const crY1Z0 = crY1.mul(crZ0);
        const crY1Z2 = crY1.mul(crZ2);
        const crY2Z0 = crY2.mul(crZ0);
        const crY2Z1 = crY2.mul(crZ1);

        const crC0 = crY1Z2.sub(crY2Z1);
        const crC1 = crY2Z0.sub(crY0Z2);
        const crC2 = crY0Z1.sub(crY1Z0);
        const crX0C0 = crX0.mul(crC0);
        const crX1C1 = crX1.mul(crC1);
        const crX2C2 = crX2.mul(crC2);
        const crTerm = crX0C0.add(crX1C1);
        const crDet = crTerm.add(crX2C2);
        return -crDet.getSign();
    }

    // The port of the upstream 'bool GetContainingTriangle(pIndex, tri)' with
    // its Triangle*& in/out parameter; the port returns both.
    private getContainingTriangleOfVertex(pIndex: number, tri: ETManifoldMeshTriangle):
        { found: boolean; tri: ETManifoldMeshTriangle } {
        const numTriangles = this.mGraph.getNumTriangles();
        for (let t = 0; t < numTriangles; ++t) {
            let j: number;
            for (j = 0; j < 3; ++j) {
                const v0Index = tri.V[mIndex[j][0]];
                const v1Index = tri.V[mIndex[j][1]];
                if (this.toLine(pIndex, v0Index, v1Index) > 0) {
                    // The point sees edge <v0,v1> from outside the triangle.
                    const adjTri = tri.T[j];
                    if (adjTri) {
                        // Traverse to the triangle sharing the edge.
                        tri = adjTri;
                        break;
                    }
                    else {
                        // We reached a hull edge, so the point is outside the
                        // hull.
                        return { found: false, tri };
                    }
                }
            }

            if (j === 3) {
                // The point is inside all three edges, so the point is inside
                // a triangle.
                return { found: true, tri };
            }
        }

        logError('Unexpected termination of loop while searching for a triangle.');
        return { found: false, tri };
    }

    private getAndRemoveInsertionPolygon(pIndex: number,
        candidates: Set<ETManifoldMeshTriangle>, boundary: DirectedEdgeKeySet): void {
        // Locate the triangles that make up the insertion polygon.
        const polygon = new ETManifoldMesh();
        while (candidates.size > 0) {
            const tri = candidates.values().next().value as ETManifoldMeshTriangle;
            candidates.delete(tri);

            for (let j = 0; j < 3; ++j) {
                const adj = tri.T[j];
                if (adj && !candidates.has(adj)) {
                    const v0Index = adj.V[0];
                    const v1Index = adj.V[1];
                    const v2Index = adj.V[2];
                    if (IncrementalDelaunay2.isDelaunayTriangle(v0Index, v1Index, v2Index) &&
                        this.toCircumcircle(pIndex, v0Index, v1Index, v2Index) <= 0) {
                        // Point P is in the circumcircle.
                        candidates.add(adj);
                    }
                }
            }

            const inserted = polygon.insert(tri.V[0], tri.V[1], tri.V[2]);
            logAssert(inserted !== null, 'Unexpected insertion failure.');
            const removed = this.mGraph.remove(tri.V[0], tri.V[1], tri.V[2]);
            logAssert(removed, 'Unexpected removal failure.');
        }

        // Get the boundary edges of the insertion polygon.
        for (const tri of polygon.getTriangles()) {
            for (let j = 0; j < 3; ++j) {
                if (!tri.T[j]) {
                    boundary.insert(tri.V[mIndex[j][0]], tri.V[mIndex[j][1]]);
                }
            }
        }
    }

    private update(pIndex: number): void {
        const tmap = this.mGraph.getTriangles();
        const containing = this.getContainingTriangleOfVertex(pIndex, tmap[0]);
        if (containing.found) {
            // The point is inside the convex hull. The insertion polygon
            // contains only triangles in the current triangulation; the hull
            // does not change.

            // Use a depth-first search for those triangles whose
            // circumcircles contain point P.
            const candidates = new Set<ETManifoldMeshTriangle>();
            candidates.add(containing.tri);

            // Get the boundary of the insertion polygon C that contains the
            // triangles whose circumcircles contain point P. Polygon C
            // contains this point.
            const boundary = new DirectedEdgeKeySet();
            this.getAndRemoveInsertionPolygon(pIndex, candidates, boundary);

            // The insertion polygon consists of the triangles formed by point
            // P and the faces of C.
            for (const key of boundary.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                if (this.toLine(pIndex, v0Index, v1Index) < 0) {
                    const inserted = this.mGraph.insert(pIndex, v0Index, v1Index);
                    logAssert(inserted !== null, 'Unexpected insertion failure.');
                }
            }
        }
        else {
            // The point is outside the convex hull. The insertion polygon is
            // formed by point P and any triangles in the current
            // triangulation whose circumcircles contain point P.

            // Locate the convex hull of the triangles.
            const hull = new DirectedEdgeKeySet();
            for (const t of tmap) {
                for (let j = 0; j < 3; ++j) {
                    if (!t.T[j]) {
                        hull.insert(t.V[mIndex[j][0]], t.V[mIndex[j][1]]);
                    }
                }
            }

            // Iterate over all the hull edges and use the ones visible to
            // point P to locate the insertion polygon.
            const candidates = new Set<ETManifoldMeshTriangle>();
            const visible = new DirectedEdgeKeySet();
            for (const key of hull.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                if (this.toLine(pIndex, v0Index, v1Index) > 0) {
                    const edge = this.mGraph.getEdge(v0Index, v1Index);
                    if (edge !== null && edge.T[1] === null) {
                        const adj = edge.T[0];
                        if (adj && !candidates.has(adj)) {
                            const a0Index = adj.V[0];
                            const a1Index = adj.V[1];
                            const a2Index = adj.V[2];
                            if (this.toCircumcircle(pIndex, a0Index, a1Index, a2Index) <= 0) {
                                // Point P is in the circumcircle.
                                candidates.add(adj);
                            }
                            else {
                                // Point P is not in the circumcircle but the
                                // hull edge is visible.
                                visible.insert(key[0], key[1]);
                            }
                        }
                    }
                    else {
                        logError('This condition should not occur for rational arithmetic.');
                    }
                }
            }

            // Get the boundary of the insertion subpolygon C that contains
            // the triangles whose circumcircles contain point P.
            const boundary = new DirectedEdgeKeySet();
            this.getAndRemoveInsertionPolygon(pIndex, candidates, boundary);

            // The insertion polygon P consists of the triangles formed by
            // point P and the back edges of C and by the visible edges of
            // mGraph-C.
            for (const key of boundary.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                if (this.toLine(pIndex, v0Index, v1Index) < 0) {
                    // This is a back edge of the boundary.
                    const inserted = this.mGraph.insert(pIndex, v0Index, v1Index);
                    logAssert(inserted !== null, 'Unexpected insertion failure.');
                }
            }
            for (const key of visible.keys()) {
                const inserted = this.mGraph.insert(pIndex, key[1], key[0]);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }
        }
    }

    // Support for triangulating the removal polygon.
    private computeWeight(iConvexIndex: number, vRemovalIndex: number,
        rpPolygon: RPPolygon): RPWeight {
        // Get the triangle <VP,VC,VN> with convex vertex VC.
        const tri = rpPolygon.getTriangle(iConvexIndex);

        const VP = this.mIRVertices[tri.vPrev];
        const VC = this.mIRVertices[tri.vCurr];
        const VN = this.mIRVertices[tri.vNext];
        const PR = this.mIRVertices[vRemovalIndex];

        const subVCVP: RationalPoint2 = [VC[0].sub(VP[0]), VC[1].sub(VP[1])];
        const subVNVP: RationalPoint2 = [VN[0].sub(VP[0]), VN[1].sub(VP[1])];
        const subPRVP: RationalPoint2 = [PR[0].sub(VP[0]), PR[1].sub(VP[1])];
        const addVCVP: RationalPoint2 = [VC[0].add(VP[0]), VC[1].add(VP[1])];
        const addVNVP: RationalPoint2 = [VN[0].add(VP[0]), VN[1].add(VP[1])];
        const addPRVP: RationalPoint2 = [PR[0].add(VP[0]), PR[1].add(VP[1])];
        const c20 = dotPerpR(subVNVP, subPRVP);
        const c21 = dotPerpR(subPRVP, subVCVP);
        const c22 = dotPerpR(subVCVP, subVNVP);
        const a20 = dotR(subVCVP, addVCVP);
        const a21 = dotR(subVNVP, addVNVP);
        const a22 = dotR(subPRVP, addPRVP);

        const weight = new RPWeight(RPWeightType.finite);
        weight.numerator = a20.mul(c20).add(a21.mul(c21)).add(a22.mul(c22)).negated();
        weight.denominator = c22;
        if (weight.denominator.getSign() < 0) {
            weight.numerator = weight.numerator.negated();
            weight.denominator = weight.denominator.negated();
        }
        return weight;
    }

    private doEarClipping(earHeap: MinHeap<number, RPWeight>,
        weightFunction: (i: number) => RPWeight, rpPolygon: RPPolygon): void {
        // Remove the finite-weight vertices from the priority queue, one at a
        // time.
        while (earHeap.getNumElements() >= 3) {
            // Get the ear of minimum weight. The vertex at index i must be
            // convex.
            const minimum = earHeap.getMinimum();
            logAssert(minimum !== null, 'Unexpected condition.');
            if (minimum.value.type !== RPWeightType.finite) {
                break;
            }
            const removedMin = earHeap.remove();
            logAssert(removedMin !== null, 'Unexpected condition.');
            const i = removedMin.key;

            // Get the triangle associated with the ear.
            const tri = rpPolygon.getTriangle(i);

            // Insert the triangle into the graph.
            const inserted = this.mGraph.insert(tri.vPrev, tri.vCurr, tri.vNext);
            logAssert(inserted !== null, 'Unexpected insertion failure.');
            if (earHeap.getNumElements() < 3) {
                earHeap.reset(0);
                break;
            }

            // Remove the vertex from the polygon. The previous and next
            // neighbor indices are required to update the adjacent vertices
            // after the removal.
            const vertex = rpPolygon.vertex(i);
            const iPrev = vertex.iPrev;
            const iNext = vertex.iNext;
            rpPolygon.remove(i);

            // Removal of the ear can cause an adjacent vertex to become an
            // ear or to stop being an ear.
            const vertexP = rpPolygon.vertex(iPrev);
            const wasConvexP = vertexP.isConvex;
            rpPolygon.classify(iPrev, this.mToLineWrapper);
            const nowConvexP = vertexP.isConvex;
            if (wasConvexP) {
                // The clipped vertex is convex. If 'vertexP' was convex, it
                // cannot become reflex after the ear is clipped.
                logAssert(nowConvexP, 'Unexpected condition.');

                if (vertexP.record !== null &&
                    vertexP.record.value.type !== RPWeightType.unmodifiable) {
                    earHeap.update(vertexP.record, weightFunction(iPrev));
                }
            }
            else if (nowConvexP) {
                // 'vertexP' was reflex and is now convex.
                if (vertexP.record !== null &&
                    vertexP.record.value.type !== RPWeightType.unmodifiable) {
                    earHeap.update(vertexP.record, weightFunction(iPrev));
                }
            }

            const vertexN = rpPolygon.vertex(iNext);
            const wasConvexN = vertexN.isConvex;
            rpPolygon.classify(iNext, this.mToLineWrapper);
            const nowConvexN = vertexN.isConvex;
            if (wasConvexN) {
                logAssert(nowConvexN, 'Unexpected condition.');

                if (vertexN.record !== null &&
                    vertexN.record.value.type !== RPWeightType.unmodifiable) {
                    earHeap.update(vertexN.record, weightFunction(iNext));
                }
            }
            else if (nowConvexN) {
                if (vertexN.record !== null &&
                    vertexN.record.value.type !== RPWeightType.unmodifiable) {
                    earHeap.update(vertexN.record, weightFunction(iNext));
                }
            }
        }
    }

    private deleteRemovalPolygon(vRemovalIndex: number,
        adjacents: readonly ETManifoldMeshTriangle[]): number[] {
        // Get the edges of the removal polygon. The polygon is star shaped
        // relative to the removal position.
        const edges = new Map<number, number>();
        for (const adj of adjacents) {
            let i: number;
            for (i = 0; i < 3; ++i) {
                if (vRemovalIndex === adj.V[i]) {
                    break;
                }
            }
            logAssert(i < 3, 'Unexpected condition.');

            const opposite1 = adj.V[(i + 1) % 3];
            const opposite2 = adj.V[(i + 2) % 3];
            if (!edges.has(opposite1)) {
                edges.set(opposite1, opposite2);
            }
        }

        // Remove the triangles. The upstream std::map is traversed in
        // increasing key order.
        const sortedKeys = Array.from(edges.keys()).sort((a, b) => a - b);
        for (const first of sortedKeys) {
            const second = edges.get(first) as number;
            const removed = this.mGraph.remove(vRemovalIndex, first, second);
            logAssert(removed, 'Unexpected removal failure.');
        }

        // Create the removal polygon; its vertices are counterclockwise
        // ordered.
        const polygon: number[] = [];
        const vStart = sortedKeys[0];
        let vCurr = edges.get(vStart) as number;
        polygon.push(vStart);
        while (vCurr !== vStart) {
            polygon.push(vCurr);
            const next = edges.get(vCurr);
            logAssert(next !== undefined, 'Unexpected condition.');
            vCurr = next;
        }
        return polygon;
    }

    private retriangulateInteriorRemovalPolygon(vRemovalIndex: number,
        polygon: readonly number[]): void {
        // Create a representation of 'polygon' that can be processed using a
        // priority queue.
        const rpPolygon = new RPPolygon(polygon, this.mToLineWrapper);

        const weightFunction = (i: number) =>
            this.computeWeight(i, vRemovalIndex, rpPolygon);

        // Create a priority queue of vertices. Convex vertices have a finite
        // and positive weight. Reflex vertices have a weight of +infinity.
        const earHeap = new MinHeap<number, RPWeight>(polygon.length, RPWeight.lessThan);
        for (let i = 0; i < polygon.length; ++i) {
            const vertex = rpPolygon.vertex(i);
            const weight = (vertex.isConvex ? weightFunction(i)
                : new RPWeight(RPWeightType.infinite));
            vertex.record = earHeap.insert(i, weight);
        }

        // Remove the finite-weight vertices from the priority queue, one at a
        // time.
        this.doEarClipping(earHeap, weightFunction, rpPolygon);
        logAssert(earHeap.getNumElements() === 0,
            'Expecting the hole to be completely filled.');
    }

    private retriangulateBoundaryRemovalPolygon(vRemovalIndex: number,
        polygon: readonly number[]): void {
        const numPolygon = polygon.length;
        if (numPolygon >= 3) {
            // Create a representation of 'polygon' that can be processed
            // using a priority queue.
            const rpPolygon = new RPPolygon(polygon, this.mToLineWrapper);

            const weightFunction = (i: number) =>
                this.computeWeight(i, vRemovalIndex, rpPolygon);

            const zeroWeightFunction = (_i: number) => new RPWeight(RPWeightType.finite);

            // Create a priority queue of vertices. The vertices adjacent to a
            // supervertex are unmodifiable. Of the other vertices, convex
            // vertices have a finite and positive weight and reflex vertices
            // have a weight of +infinity.
            const earHeap = new MinHeap<number, RPWeight>(numPolygon, RPWeight.lessThan);

            let iPrev = numPolygon - 2;
            let iCurr = iPrev + 1;
            let iNext = 0;
            for (; iNext < numPolygon; iPrev = iCurr, iCurr = iNext, ++iNext) {
                const vertexPrev = rpPolygon.vertex(iPrev);
                const vertexCurr = rpPolygon.vertex(iCurr);
                const vertexNext = rpPolygon.vertex(iNext);
                if (IncrementalDelaunay2.isSupervertex(vertexPrev.vIndex) ||
                    IncrementalDelaunay2.isSupervertex(vertexCurr.vIndex) ||
                    IncrementalDelaunay2.isSupervertex(vertexNext.vIndex)) {
                    vertexCurr.record = earHeap.insert(iCurr,
                        new RPWeight(RPWeightType.unmodifiable));
                }
                else if (vertexCurr.isConvex) {
                    vertexCurr.record = earHeap.insert(iCurr, weightFunction(iCurr));
                }
                else {
                    vertexCurr.record = earHeap.insert(iCurr,
                        new RPWeight(RPWeightType.infinite));
                }
            }

            // Remove the finite-weight vertices from the priority queue, one
            // at a time. This process fills in the subpolygon of the removal
            // polygon that is contained by the Delaunay triangulation.
            this.doEarClipping(earHeap, weightFunction, rpPolygon);

            // Get the subpolygon of the removal polygon that is external to
            // the Delaunay triangulation.
            let numExternal = rpPolygon.getNumActive();
            let external = new Array<number>(numExternal);
            iCurr = rpPolygon.getActive();
            for (let i = 0; i < numExternal; ++i) {
                external[i] = iCurr;
                rpPolygon.classify(iCurr, this.mToLineWrapper);
                iCurr = rpPolygon.vertex(iCurr).iNext;
            }

            earHeap.reset(numExternal);
            for (let i = 0; i < numExternal; ++i) {
                const index = external[i];
                const vertex = rpPolygon.vertex(index);
                if (IncrementalDelaunay2.isSupervertex(vertex.vIndex)) {
                    vertex.record = earHeap.insert(index,
                        new RPWeight(RPWeightType.unmodifiable));
                }
                else if (vertex.isConvex) {
                    vertex.record = earHeap.insert(index, zeroWeightFunction(index));
                }
                else {
                    vertex.record = earHeap.insert(index,
                        new RPWeight(RPWeightType.infinite));
                }
            }

            // Remove the finite-weight vertices from the priority queue, one
            // at a time. This process fills in a portion or all of the
            // subpolygon of the removal polygon that is external to the
            // Delaunay triangulation.
            this.doEarClipping(earHeap, zeroWeightFunction, rpPolygon);
            if (earHeap.getNumElements() === 0) {
                // The external polygon contained only 1 supervertex.
                return;
            }

            // The remaining external polygon is a triangle fan with 2 or 3
            // supervertices.
            numExternal = rpPolygon.getNumActive();
            external = new Array<number>(numExternal);
            iCurr = rpPolygon.getActive();
            for (let i = 0; i < numExternal; ++i) {
                external[i] = iCurr;
                rpPolygon.classify(iCurr, this.mToLineWrapper);
                iCurr = rpPolygon.vertex(iCurr).iNext;
            }

            earHeap.reset(numExternal);
            iPrev = numExternal - 2;
            iCurr = iPrev + 1;
            iNext = 0;
            for (; iNext < numExternal; iPrev = iCurr, iCurr = iNext, ++iNext) {
                const index = external[iCurr];
                const vertexPrev = rpPolygon.vertex(external[iPrev]);
                const vertexCurr = rpPolygon.vertex(index);
                const vertexNext = rpPolygon.vertex(external[iNext]);
                if (IncrementalDelaunay2.isSupervertex(vertexCurr.vIndex)) {
                    if (IncrementalDelaunay2.isDelaunayVertex(vertexPrev.vIndex) ||
                        IncrementalDelaunay2.isDelaunayVertex(vertexNext.vIndex)) {
                        logAssert(vertexCurr.isConvex, 'Unexpected condition.');

                        vertexCurr.record = earHeap.insert(index,
                            zeroWeightFunction(index));
                    }
                    else {
                        vertexCurr.record = earHeap.insert(index,
                            new RPWeight(RPWeightType.unmodifiable));
                    }
                }
                else {
                    vertexCurr.record = earHeap.insert(index,
                        new RPWeight(RPWeightType.infinite));
                }
            }

            // Remove the finite-weight vertices from the priority queue, one
            // at a time. This process fills in the triangle fan of the
            // subpolygon of the removal polygon that is external to the
            // Delaunay triangulation.
            this.doEarClipping(earHeap, zeroWeightFunction, rpPolygon);
            logAssert(earHeap.getNumElements() === 0,
                'Expecting the hole to be completely filled.');
        }
        else {
            // numPolygon == 2. The removal polygon degenerates to a 2-gon,
            // which happens only when a single Delaunay vertex remains after
            // the removal. The graph is rebuilt as the fan of that vertex
            // with the supertriangle.
            //
            // Upstream bug (fixed in the port): upstream selects the surviving
            // vertex with
            //   vOtherIndex = (polygon[0] == vRemovalIndex ? polygon[1] : polygon[0]);
            // but 'polygon' is built from the vertices *opposite* the removal
            // vertex in the removal triangles, so it never contains
            // vRemovalIndex; the test is always false and polygon[0] is
            // always chosen. When polygon[0] is a supervertex, the rebuilt fan
            // is centered on a supervertex and the triangulation is corrupted.
            // The port selects the Delaunay vertex of the 2-gon instead, which
            // is the evident intent of the name 'vOtherIndex'. (The branch
            // appears to be unreachable in practice: the single-remaining-
            // vertex case is handled earlier by the mVertexIndexMap.size == 4
            // path in remove(). The fix is defensive.)
            const vOtherIndex = (IncrementalDelaunay2.isDelaunayVertex(polygon[0])
                ? polygon[0] : polygon[1]);

            this.mGraph.clear();
            for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
                const inserted = this.mGraph.insert(vOtherIndex, i0, i1);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }
        }
    }

    private updateTrianglesAndAdjacenciesIfNeeded(): void {
        if (this.mTrianglesAndAdjacenciesNeedUpdate) {
            this.updateTrianglesAndAdjacencies();
            this.mTrianglesAndAdjacenciesNeedUpdate = false;
        }
    }

    private updateTrianglesAndAdjacencies(): void {
        // Assign integer values to the triangles.
        const tmap = this.mGraph.getTriangles();
        if (tmap.length === 0) {
            this.mTriangles = [];
            this.mAdjacencies = [];
            return;
        }

        const permute = new Map<ETManifoldMeshTriangle | null, number>();
        permute.set(null, invalid);
        let numTriangles = 0;
        for (const tri of tmap) {
            if (IncrementalDelaunay2.isDelaunayTriangle(tri.V[0], tri.V[1], tri.V[2])) {
                permute.set(tri, numTriangles++);
            }
            else {
                permute.set(tri, invalid);
            }
        }

        this.mTriangles = new Array<[number, number, number]>(numTriangles);
        this.mAdjacencies = new Array<[number, number, number]>(numTriangles);
        let t = 0;
        for (const tri of tmap) {
            if (permute.get(tri) !== invalid) {
                this.mTriangles[t] = [tri.V[0], tri.V[1], tri.V[2]];
                this.mAdjacencies[t] = [
                    permute.get(tri.T[0]) as number,
                    permute.get(tri.T[1]) as number,
                    permute.get(tri.T[2]) as number
                ];
                ++t;
            }
        }
    }
}

// Dot(u, v) for rational 2D points.
function dotR(u: RationalPoint2, v: RationalPoint2): BSNumber {
    return u[0].mul(v[0]).add(u[1].mul(v[1]));
}

// DotPerp(u, v) = Dot(Perp(u), v) = u.x * v.y - u.y * v.x for rational 2D
// points.
function dotPerpR(u: RationalPoint2, v: RationalPoint2): BSNumber {
    return u[0].mul(v[1]).sub(u[1].mul(v[0]));
}
