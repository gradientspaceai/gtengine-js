// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Delaunay2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Delaunay triangulation of points (intrinsic dimensionality 2).
//   VQ = number of vertices
//   V  = array of vertices
//   TQ = number of triangles
//   I  = Array of 3-tuples of indices into V that represent the triangles
//        (3*TQ total elements). Access via getIndices().
//   A  = Array of 3-tuples of indices into I that represent the adjacent
//        triangles (3*TQ total elements). Access via getAdjacencies().
// The i-th triangle has vertices
//   vertex[0] = V[I[3*i+0]]
//   vertex[1] = V[I[3*i+1]]
//   vertex[2] = V[I[3*i+2]]
// and edge index pairs
//   edge[0] = <I[3*i+0],I[3*i+1]>
//   edge[1] = <I[3*i+1],I[3*i+2]>
//   edge[2] = <I[3*i+2],I[3*i+0]>
// The triangles adjacent to these edges have indices
//   adjacent[0] = A[3*i+0] is the triangle sharing edge[0]
//   adjacent[1] = A[3*i+1] is the triangle sharing edge[1]
//   adjacent[2] = A[3*i+2] is the triangle sharing edge[2]
// If there is no adjacent triangle, the A[*] value is set to -1.
//
// The class uses a blend of interval arithmetic (the fast path) and exact
// rational arithmetic (the fallback when the interval straddles zero) for the
// ToLine and ToCircumcircle sign classifications, which is what makes the
// triangulation robust.
//
// Port notes:
// * Upstream Delaunay2.h declares a variadic class template with two
//   specializations: the deprecated Delaunay2<InputType, ComputeType>, which
//   upstream states "will be removed in a future release", and the
//   replacement Delaunay2<T>. Only the replacement is ported. The deprecated
//   class routes its predicates through PrimalQuery2<ComputeType>; the
//   replacement embeds its own interval/rational predicates, which is the
//   behavior ported here (src/PrimalQuery2.ts is number-only and is not used
//   by this file).
// * Upstream selects InputRational = BSNumber<UIntegerFP32<2 or 4>> and
//   ComputeRational = BSNumber<UIntegerFP32<36 or 264>> with a preallocated
//   pool of compute-rational scratch values (mCRPool) and a Copy() helper
//   that widens an input rational into a compute rational. The port's
//   BSNumber is bigint-backed and grows as needed, so the fixed word counts,
//   the two distinct rational types, the pool and Copy() are all unnecessary
//   and are dropped. The exact-arithmetic results are identical.
// * std::unordered_set of directed EdgeKey<true> becomes DirectedEdgeKeySet,
//   a Map keyed by the directed pair and iterated in sorted key order so the
//   triangulation is deterministic.
// * size_t 'negOne' (std::numeric_limits<size_t>::max()) becomes -1.
// * The 'bool GetHull(std::vector<size_t>&)' output-parameter form becomes
//   'getHull(): number[]'; the bool was 'true' on every path that returns.
// * The upstream operator() becomes compute(vertices).

import { logAssert, logError } from './Logger.js';
import { Line } from './Line.js';
import type { Line2 } from './Line.js';
import { Vector } from './Vector.js';
import { IntrinsicsVector2 } from './Vector2.js';
import { BSNumber } from './BSNumber.js';
import { SWInterval } from './SWInterval.js';
import { VETManifoldMesh } from './VETManifoldMesh.js';
import { ETManifoldMesh, ETManifoldMeshTriangle } from './ETManifoldMesh.js';

// A rational 2D point, the port of Vector2<InputRational>.
type RationalPoint2 = [BSNumber, BSNumber];

// The sentinel for "the query point is not one of the input vertices"; the
// port of the upstream size_t constant negOne.
const negOne = -1;

// Indexing for the vertices of the triangle adjacent to a vertex. The edge
// adjacent to vertex j is <mIndex[j][0], mIndex[j][1]> and is listed so that
// the triangle interior is to your left as you walk around the edges.
const mIndex: readonly (readonly [number, number])[] =
    [[0, 1], [1, 2], [2, 0]];

// The port of std::unordered_set<EdgeKey<true>, ...>. The port iterates in
// increasing key order so that the constructed triangulation does not depend
// on hash-table ordering.
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
export class Delaunay2SearchInfo {
    initialTriangle: number;
    numPath: number;
    finalTriangle: number;
    finalV: [number, number, number];
    path: number[];

    constructor() {
        this.initialTriangle = negOne;
        this.numPath = 0;
        this.finalTriangle = 0;
        this.finalV = [0, 0, 0];
        this.path = [];
    }
}

export class Delaunay2 {
    // The vertices used for geometric queries. The input vertices are
    // read-only, so they can be represented exactly by rational numbers.
    protected mNumVertices: number;
    protected mVertices: readonly Vector[];
    protected mIRVertices: RationalPoint2[];

    protected mGraph: VETManifoldMesh;

    // If a vertex occurs multiple times in the 'vertices' input, the first
    // processed occurrence of that vertex has an index stored in this array.
    // If there are no duplicates, then mDuplicates[i] = i for all i.
    protected mDuplicates: number[];
    protected mNumUniqueVertices: number;

    // If the intrinsic dimension of the input vertices is 0 or 1, compute()
    // returns early. The caller is responsible for retrieving the dimension
    // and taking an alternate path should the dimension be smaller than 2. If
    // the dimension is 0, all vertices are the same. If the dimension is 1,
    // the vertices lie on a line, in which case the caller can project
    // vertices[] onto the line for further processing.
    protected mDimension: number;
    protected mLine: Line2;

    // These are computed by updateIndicesAdjacencies(). They are used for
    // point-containment queries in the triangle mesh.
    protected mNumTriangles: number;
    protected mIndices: number[];
    protected mAdjacencies: number[];

    // The query point for update(), getContainingTriangle() and
    // getAndRemoveInsertionPolygon() when the point is not an input vertex.
    // toLine() and toCircumcircle() are passed indices into the vertex array.
    // When the index is valid, mVertices[] and mIRVertices[] are used for
    // lookups. When the index is negOne, the query point is used.
    private mQueryPoint: Vector;
    private mIRQueryPoint: RationalPoint2;

    constructor() {
        this.mNumVertices = 0;
        this.mVertices = [];
        this.mIRVertices = [];
        this.mGraph = new VETManifoldMesh();
        this.mDuplicates = [];
        this.mNumUniqueVertices = 0;
        this.mDimension = 0;
        this.mLine = new Line(2);
        this.mLine.direction.makeZero();
        this.mNumTriangles = 0;
        this.mIndices = [];
        this.mAdjacencies = [];
        this.mQueryPoint = new Vector(2);
        this.mIRQueryPoint = [BSNumber.fromNumber(0), BSNumber.fromNumber(0)];
    }

    // The input is the array of vertices whose Delaunay triangulation is
    // required. The return value is 'true' if and only if the intrinsic
    // dimension of the points is 2. If the intrinsic dimension is 1, the
    // points lie exactly on a line which is then accessible via getLine().
    // If the intrinsic dimension is 0, the points are all the same point.
    compute(vertices: readonly Vector[]): boolean {
        // Initialize values in case they were set by a previous call to
        // compute().
        logAssert(vertices.length > 0, 'Invalid argument.');
        for (const vertex of vertices) {
            logAssert(vertex.size === 2, 'Delaunay2 requires 2D vertices.');
        }

        this.mNumVertices = vertices.length;
        this.mVertices = vertices;
        this.mIRVertices = [];
        this.mDuplicates = [];
        this.mLine.origin.makeZero();
        this.mLine.direction.makeZero();
        this.mNumUniqueVertices = 0;
        this.mNumTriangles = 0;
        this.mGraph = new VETManifoldMesh();
        this.mIndices = [];
        this.mAdjacencies = [];
        this.mQueryPoint = new Vector(2);
        this.mIRQueryPoint = [BSNumber.fromNumber(0), BSNumber.fromNumber(0)];

        // Compute the intrinsic dimension and return early if that dimension
        // is 0 or 1.
        const info = new IntrinsicsVector2(vertices, 0);
        if (info.dimension === 0) {
            // The vertices are the same point.
            this.mDimension = 0;
            this.mLine.origin = info.origin.clone();
            return false;
        }

        // Convert the floating-point inputs to rational type. Upstream does
        // this only after it has decided the intrinsic dimension is 2; the
        // port needs the rational vertices earlier, for the exact
        // classification below.
        this.mIRVertices = new Array<RationalPoint2>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            const v = vertices[i].values;
            this.mIRVertices[i] = [BSNumber.fromNumber(v[0]), BSNumber.fromNumber(v[1])];
        }

        // Classify the intrinsic dimension exactly and select the seed
        // triangle from the extreme vertices found by the intrinsics
        // computation.
        //
        // Upstream bug (fixed in the port): IntrinsicsVector2 determines the
        // intrinsic dimension in floating-point arithmetic and Delaunay2<T>
        // hardcodes its epsilon to 0. The direction of the line through the
        // first two extremes is normalized, so it carries roundoff, and the
        // perpendicular distances of an exactly collinear input set are then
        // nonzero at the 1-ulp level. Such a set is reported as dimension 2
        // with an exactly degenerate extreme triangle. For example, the two
        // vertices (-3,-9) and (0,0) yield info.dimension = 2 with
        // info.extreme = [0,1,1]: upstream inserts the degenerate triangle
        // <0,1,1> (its LogAssert only tests that the mesh insertion returned
        // a triangle, not that the triangle has area) and returns true with
        // getIndices() = [0,1,1]. Four collinear vertices such as (0,0),
        // (1,3), (2,6), (3,9) instead throw from deep inside the incremental
        // update. The class's own exact ToLine predicate settles the
        // question: keep upstream's extreme[2] when the seed triangle really
        // is nondegenerate; otherwise look for a vertex that is exactly off
        // the line through the first two extremes; and if there is no such
        // vertex the input is exactly collinear, which is dimension 1 with
        // the line upstream's dimension-1 branch would have reported. The
        // exact sign also replaces info.extremeCCW, which is the same
        // roundoff-prone quantity, so the seed triangle is counterclockwise
        // as the circumcircle-visibility algorithm requires.
        //
        // The dimension-0 case above needs no such check: with epsilon = 0
        // its test is 'maxRange == 0', which holds exactly when the bounding
        // box of the input is a single point.
        const e0 = info.extreme[0], e1 = info.extreme[1];
        let e2 = info.extreme[2];
        let toLineSign = (e2 !== e0 && e2 !== e1 ? this.toLine(e2, e0, e1) : 0);
        if (toLineSign === 0) {
            for (let i = 0; i < this.mNumVertices; ++i) {
                if (i === e0 || i === e1) {
                    continue;
                }
                const sign = this.toLine(i, e0, e1);
                if (sign !== 0) {
                    e2 = i;
                    toLineSign = sign;
                    break;
                }
            }
        }
        if (toLineSign === 0) {
            // The vertices are exactly collinear.
            this.mDimension = 1;
            this.mIRVertices = [];
            this.mLine.origin = info.origin.clone();
            this.mLine.direction = info.direction[0].clone();
            return false;
        }

        // The vertices necessarily will have a triangulation.
        this.mDimension = 2;

        // Assume initially the vertices are unique. If duplicates are found
        // during the Delaunay update, mDuplicates[] will be modified
        // accordingly.
        this.mDuplicates = new Array<number>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mDuplicates[i] = i;
        }

        // Insert the nondegenerate triangle constructed by the intrinsics
        // computation. This is necessary for the circumcircle-visibility
        // algorithm to work correctly. toLine returns -1 when the point is to
        // the left of the directed line, which is the counterclockwise
        // orientation; upstream swaps extreme[1] and extreme[2] in the
        // clockwise case.
        const extreme: [number, number, number] =
            (toLineSign < 0 ? [e0, e1, e2] : [e0, e2, e1]);

        const inserted = this.mGraph.insert(extreme[0], extreme[1], extreme[2]);
        logAssert(inserted !== null, 'The triangle should not be degenerate.');

        // Incrementally update the triangulation. The set of processed points
        // is maintained to eliminate duplicates. The upstream
        // std::unordered_set<ProcessedVertex> is a Map keyed by the exact
        // floating-point coordinate pair; the mapped value is the 'location'
        // member.
        const processed = new Map<string, number>();
        const keyOf = (i: number): string => {
            const v = vertices[i].values;
            return `${v[0]},${v[1]}`;
        };
        for (let i = 0; i < 3; ++i) {
            const j = extreme[i];
            processed.set(keyOf(j), j);
            this.mDuplicates[j] = j;
        }
        for (let i = 0; i < this.mNumVertices; ++i) {
            const key = keyOf(i);
            const location = processed.get(key);
            if (location === undefined) {
                this.update(i);
                processed.set(key, i);
                this.mDuplicates[i] = i;
            }
            else {
                this.mDuplicates[i] = location;
            }
        }
        this.mNumUniqueVertices = processed.size;

        // Assign integer values to the triangles for use by the caller and
        // copy the triangle information to compact arrays mIndices and
        // mAdjacencies.
        this.updateIndicesAdjacencies();

        return true;
    }

    // Dimensional information. If getDimension() returns 1, the points lie on
    // a line P+t*D. You can sort these if you need a polyline output by
    // projecting onto the line each vertex X = P+t*D, where t = Dot(D,X-P).
    getDimension(): number {
        return this.mDimension;
    }

    getLine(): Line2 {
        return this.mLine;
    }

    // Member access.
    //
    // Upstream bug (fixed in the port): GetNumVertices() returns
    // mIRVertices.size(), which is 0 whenever the intrinsic dimension is 0 or
    // 1 because operator() clears mIRVertices and returns before resizing it.
    // The port returns mNumVertices, which is the number of input vertices on
    // every path.
    getNumVertices(): number {
        return this.mNumVertices;
    }

    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getNumUniqueVertices(): number {
        return this.mNumUniqueVertices;
    }

    // If 'vertices' has no duplicates, getDuplicates()[i] = i for all i. If
    // vertices[i] is the first occurrence of a vertex and if vertices[j] is
    // found later, then getDuplicates()[j] = i.
    getDuplicates(): readonly number[] {
        return this.mDuplicates;
    }

    getNumTriangles(): number {
        return this.mNumTriangles;
    }

    getGraph(): ETManifoldMesh {
        return this.mGraph;
    }

    getIndices(): readonly number[] {
        return this.mIndices;
    }

    getAdjacencies(): readonly number[] {
        return this.mAdjacencies;
    }

    // Locate those triangle edges that do not share other triangles. The
    // returned array has hull.length = 2*numEdges, each pair representing an
    // edge. The edges are not ordered, but the pair of vertices for an edge
    // is ordered so that they conform to a counterclockwise traversal of the
    // hull. The dimension must be 2.
    getHull(): number[] {
        if (this.mDimension === 2) {
            // Count the number of edges that are not shared by two triangles.
            let numEdges = 0;
            for (const adj of this.mAdjacencies) {
                if (adj === -1) {
                    ++numEdges;
                }
            }

            if (numEdges > 0) {
                // Enumerate the edges.
                const hull = new Array<number>(2 * numEdges);
                let current = 0, i = 0;
                for (const adj of this.mAdjacencies) {
                    if (adj === -1) {
                        const tri = Math.floor(i / 3), j = i % 3;
                        hull[current++] = this.mIndices[3 * tri + j];
                        hull[current++] = this.mIndices[3 * tri + ((j + 1) % 3)];
                    }
                    ++i;
                }
                return hull;
            }
            else {
                logError('Unexpected condition. There must be at least one triangle.');
            }
        }
        else {
            logError('The dimension must be 2.');
        }
        // Unreachable; logError throws.
        return [];
    }

    // Copy Delaunay triangles to compact arrays mIndices and mAdjacencies.
    updateIndicesAdjacencies(): void {
        // Assign integer values to the triangles.
        const tmap = this.mGraph.getTriangles();
        const permute = new Map<ETManifoldMeshTriangle, number>();
        for (let i = 0; i < tmap.length; ++i) {
            permute.set(tmap[i], i);
        }

        this.mNumTriangles = tmap.length;
        const numindices = 3 * this.mNumTriangles;
        if (numindices > 0) {
            this.mIndices = new Array<number>(numindices);
            this.mAdjacencies = new Array<number>(numindices);
            let i = 0;
            for (const tri of tmap) {
                for (let j = 0; j < 3; ++j, ++i) {
                    this.mIndices[i] = tri.V[j];
                    const adj = tri.T[j];
                    this.mAdjacencies[i] = (adj !== null ? permute.get(adj) as number : -1);
                }
            }
        }
    }

    // Get the vertex indices for triangle t. The function returns the
    // indices when the dimension is 2 and t is a valid triangle index;
    // otherwise, it returns null.
    getTriangleIndices(t: number): [number, number, number] | null {
        if (this.mDimension === 2) {
            const numTriangles = this.mIndices.length / 3;
            if (0 <= t && t < numTriangles) {
                return [this.mIndices[3 * t], this.mIndices[3 * t + 1],
                    this.mIndices[3 * t + 2]];
            }
        }
        return null;
    }

    // Get the indices for triangles adjacent to triangle t. The function
    // returns the adjacencies when the dimension is 2 and t is a valid
    // triangle index; otherwise, it returns null.
    getTriangleAdjacencies(t: number): [number, number, number] | null {
        if (this.mDimension === 2) {
            const numTriangles = this.mIndices.length / 3;
            if (0 <= t && t < numTriangles) {
                return [this.mAdjacencies[3 * t], this.mAdjacencies[3 * t + 1],
                    this.mAdjacencies[3 * t + 2]];
            }
        }
        return null;
    }

    // Support for searching the triangulation for a triangle that contains a
    // point. If there is a containing triangle, the returned value is a
    // triangle index t with 0 <= t < getNumTriangles(). If there is not a
    // containing triangle, -1 is returned. The computations are performed
    // using exact rational arithmetic.
    getContainingTriangle(inP: Vector, info: Delaunay2SearchInfo): number {
        logAssert(this.mDimension === 2, 'Invalid dimension for triangle search.');
        logAssert(inP.size === 2, 'Delaunay2 requires 2D vertices.');

        this.mQueryPoint = inP.clone();
        this.mIRQueryPoint = [BSNumber.fromNumber(inP.values[0]),
            BSNumber.fromNumber(inP.values[1])];

        const numTriangles = this.mIndices.length / 3;
        info.path = new Array<number>(numTriangles).fill(0);
        info.numPath = 0;
        let triangle: number;
        if (0 <= info.initialTriangle && info.initialTriangle < numTriangles) {
            triangle = info.initialTriangle;
        }
        else {
            info.initialTriangle = 0;
            triangle = 0;
        }

        // Use triangle edges as binary separating lines.
        for (let i = 0; i < numTriangles; ++i) {
            const ibase = 3 * triangle;
            const v0 = this.mIndices[ibase];
            const v1 = this.mIndices[ibase + 1];
            const v2 = this.mIndices[ibase + 2];

            info.path[info.numPath++] = triangle;
            info.finalTriangle = triangle;
            info.finalV[0] = v0;
            info.finalV[1] = v1;
            info.finalV[2] = v2;

            if (this.toLine(negOne, v0, v1) > 0) {
                const adjacent = this.mAdjacencies[ibase];
                if (adjacent === -1) {
                    info.finalV[0] = v0;
                    info.finalV[1] = v1;
                    info.finalV[2] = v2;
                    return negOne;
                }
                triangle = adjacent;
                continue;
            }

            if (this.toLine(negOne, v1, v2) > 0) {
                const adjacent = this.mAdjacencies[ibase + 1];
                if (adjacent === -1) {
                    info.finalV[0] = v1;
                    info.finalV[1] = v2;
                    info.finalV[2] = v0;
                    return negOne;
                }
                triangle = adjacent;
                continue;
            }

            if (this.toLine(negOne, v2, v0) > 0) {
                const adjacent = this.mAdjacencies[ibase + 2];
                if (adjacent === -1) {
                    info.finalV[0] = v2;
                    info.finalV[1] = v0;
                    info.finalV[2] = v1;
                    return negOne;
                }
                triangle = adjacent;
                continue;
            }

            return triangle;
        }

        logError('Unexpected termination of loop while searching for a triangle.');
        return negOne;
    }

    // The floating-point query point for the given index, where negOne
    // selects the stored query point.
    private inPoint(index: number): readonly number[] {
        return (index !== negOne ? this.mVertices[index].values : this.mQueryPoint.values);
    }

    // The rational query point for the given index, where negOne selects the
    // stored query point.
    private irPoint(index: number): RationalPoint2 {
        return (index !== negOne ? this.mIRVertices[index] : this.mIRQueryPoint);
    }

    // Given a line with origin V0 and direction <V0,V1> and a query point P,
    // toLine returns
    //   +1, P on right of line
    //   -1, P on left of line
    //    0, P on the line
    protected toLine(pIndex: number, v0Index: number, v1Index: number): number {
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
    protected toCircumcircle(pIndex: number, v0Index: number, v1Index: number,
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
        // The '| 0' reproduces the C++ negation of an int32_t: JavaScript's
        // unary minus turns the sign 0 into -0, which is a distinct value
        // under Object.is and would leak out of a predicate documented to
        // return -1, 0 or +1.
        return -crDet.getSign() | 0;
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
                // a triangle. (Upstream's comment says "four edges", a stale
                // copy from the 3D code.)
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
                    if (this.toCircumcircle(pIndex, v0Index, v1Index, v2Index) <= 0) {
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
        const first = tmap[0];
        const containing = this.getContainingTriangleOfVertex(pIndex, first);
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
                // else: Point P is on an edge of the boundary, so the
                // subdivision would have degenerate triangles. Ignore these.
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
}
