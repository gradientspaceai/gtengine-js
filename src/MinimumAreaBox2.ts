// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumAreaBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute a minimum-area oriented box containing the specified points. The
// algorithm uses the rotating calipers method, but with a dual pair of
// calipers. For details, see
// http://www-cgrl.cs.mcgill.ca/~godfried/research/calipers.html
// https://web.archive.org/web/20150330010154/http://cgm.cs.mcgill.ca/~orm/rotcal.html
// The box is supported by the convex hull of the points, so the algorithm is
// really about computing the minimum-area box containing a convex polygon.
// The rotating calipers approach is O(n) in time for n polygon edges.
//
// A detailed description of the algorithm and implementation is found in
// https://www.geometrictools.com/Documentation/MinimumAreaRectangle.pdf
//
// Port notes:
// * Upstream templates the class on <InputType, ComputeType> and warns that
//   a correct output is guaranteed only when ComputeType is an exact
//   arithmetic type that supports division (the recommended choice is
//   BSRational<UIntegerAP32>). The port instantiates only that exact path,
//   using the bigint-backed BSRational, so the results are always the
//   guaranteed-correct ones. This follows the ConvexHull2 (B100) and
//   RotatingCalipers (B103) precedent of doing the predicate arithmetic
//   exactly. Consequently 'useRotatingCalipers' defaults to true, which is
//   upstream's default when ComputeType is not a floating-point type.
// * The four operator() overloads become two methods:
//     compute(points, useRotatingCalipers)
//       the port of overloads 1 and 2, which compute the convex hull of
//       arbitrary points;
//     computeConvexPolygon(points, indices, useRotatingCalipers)
//       the port of overload 3, where the caller guarantees the input is a
//       counterclockwise nondegenerate convex polygon; an empty 'indices'
//       means the points themselves are the polygon.
//   Upstream overload 4 (the std::vector form of overload 3) is not ported
//   separately because it is inconsistent with overload 3: for an empty
//   indices vector it forwards to the convex-hull path rather than treating
//   the points as the polygon. Callers wanting the hull path call compute().
// * The internal Box struct is a module-private interface; C++ value-copy
//   assignment of it becomes an explicit cloneBox().
// * Output reference parameters become returned values; GetHull(),
//   GetSupportIndices() and GetArea() become getHull(), getSupportIndices()
//   and getArea().

import { ConvexHull2 } from './ConvexHull2.js';
import { logAssert } from './Logger.js';
import { BSRational } from './BSRational.js';
import { OrientedBox, type OrientedBox2 } from './OrientedBox.js';
import { Vector, dot, normalize, sub } from './Vector.js';

// The port of Vector2<ComputeType>.
type RationalPoint2 = [BSRational, BSRational];

const R_ZERO = BSRational.fromNumber(0);
const R_ONE = BSRational.fromNumber(1);
const R_HALF = BSRational.fromNumber(0.5);

function rAdd(v0: RationalPoint2, v1: RationalPoint2): RationalPoint2 {
    return [v0[0].add(v1[0]), v0[1].add(v1[1])];
}

function rSub(v0: RationalPoint2, v1: RationalPoint2): RationalPoint2 {
    return [v0[0].sub(v1[0]), v0[1].sub(v1[1])];
}

function rDot(v0: RationalPoint2, v1: RationalPoint2): BSRational {
    return v0[0].mul(v1[0]).add(v0[1].mul(v1[1]));
}

function rDotPerp(v0: RationalPoint2, v1: RationalPoint2): BSRational {
    return v0[0].mul(v1[1]).sub(v0[1].mul(v1[0]));
}

function rNegate(v: RationalPoint2): RationalPoint2 {
    return [v[0].negated(), v[1].negated()];
}

// The port of -Perp(v) = -(v1,-v0) = (-v1,v0).
function rNegPerp(v: RationalPoint2): RationalPoint2 {
    return [v[1].negated(), v[0]];
}

// The box axes are U[i] and are usually not unit-length in order to allow
// exact arithmetic. The box is supported by the polygon vertices
// vertices[index[i]]. The box axes are not necessarily unit length, but they
// have the same length. They need to be normalized for conversion back to
// floating point.
interface Box {
    U: [RationalPoint2, RationalPoint2];
    // The support order is bottom, right, top, left.
    index: [number, number, number, number];
    sqrLenU0: BSRational;
    area: BSRational;
}

// The port of C++ value-copy assignment of a Box. The RationalPoint2 values
// are never mutated in place, so sharing them is safe.
function cloneBox(box: Box): Box {
    return {
        U: [box.U[0], box.U[1]],
        index: [box.index[0], box.index[1], box.index[2], box.index[3]],
        sqrLenU0: box.sqrLenU0,
        area: box.area
    };
}

// An (angle, box-edge) pair of the A[] array of ComputeAngles.
interface AnglePair {
    sinThetaSqr: BSRational;
    edge: number;
}

// The rotating calipers algorithm has a loop invariant that requires the
// convex polygon not to have collinear points. Any such points must be
// removed first. The code is also executed for the O(n^2) algorithm to
// reduce the number of processed edges.
//
// Upstream bug (MinimumAreaBox2.h, RemoveCollinearPoints): the "arriving"
// edge used in the DotPerp test is the immediately preceding edge of the
// input array, which may be the zero-length edge produced by a duplicate
// point. DotPerp is then zero and a genuine polygon corner is discarded
// along with the duplicate. This is the same defect that was found in
// RotatingCalipers::CreatePolygon (see gtengine-js issue #286). It is
// unreachable through compute(), because ConvexHull2 removes duplicates and
// guarantees no three consecutive collinear hull points, but it is reachable
// through computeConvexPolygon(). The port compares against the most recent
// *nonzero* edge, which is identical to upstream whenever the input has no
// duplicate points.
function removeCollinearPoints(vertices: readonly RationalPoint2[]):
    RationalPoint2[] {
    const numVertices = vertices.length;

    const edges: RationalPoint2[] = [];
    const isZeroEdge: boolean[] = [];
    for (let i = 0; i < numVertices; ++i) {
        const edge = rSub(vertices[(i + 1) % numVertices], vertices[i]);
        edges.push(edge);
        isZeroEdge.push(edge[0].equals(R_ZERO) && edge[1].equals(R_ZERO));
    }

    // The most recent nonzero edge before index 0, cyclically.
    let prevNonzero = -1;
    for (let i = numVertices - 1; i >= 0; --i) {
        if (!isZeroEdge[i]) {
            prevNonzero = i;
            break;
        }
    }

    const result: RationalPoint2[] = [];
    if (prevNonzero < 0) {
        // All the input points are identical; there are no edges at all.
        return result;
    }

    for (let i0 = 0; i0 < numVertices; ++i0) {
        if (isZeroEdge[i0]) {
            continue;
        }
        if (rDotPerp(edges[prevNonzero], edges[i0]).notEquals(R_ZERO)) {
            result.push(vertices[i0]);
        }
        prevNonzero = i0;
    }

    return result;
}

// Compute the smallest box for the polygon edge <V[i0],V[i1]>.
function smallestBox(i0: number, i1: number,
    vertices: readonly RationalPoint2[]): Box {
    const u0 = rSub(vertices[i1], vertices[i0]);
    const box: Box = {
        U: [u0, rNegPerp(u0)],
        index: [i1, i1, i1, i1],
        sqrLenU0: rDot(u0, u0),
        area: R_ZERO
    };

    const origin = vertices[i1];
    const support: RationalPoint2[] = [
        [R_ZERO, R_ZERO], [R_ZERO, R_ZERO], [R_ZERO, R_ZERO], [R_ZERO, R_ZERO]
    ];

    for (let i = 0; i < vertices.length; ++i) {
        const diff = rSub(vertices[i], origin);
        const v: RationalPoint2 = [rDot(box.U[0], diff), rDot(box.U[1], diff)];

        // The right-most vertex of the bottom edge is vertices[i1]. The
        // assumption of no triple of collinear vertices guarantees that
        // box.index[0] is i1, which is the initial value assigned at the
        // beginning of this function. Therefore, there is no need to test
        // for other vertices farther to the right than vertices[i1].

        if (v[0].greaterThan(support[1][0]) ||
            (v[0].equals(support[1][0]) && v[1].greaterThan(support[1][1]))) {
            // New right maximum OR same right maximum but closer to top.
            box.index[1] = i;
            support[1] = v;
        }

        if (v[1].greaterThan(support[2][1]) ||
            (v[1].equals(support[2][1]) && v[0].lessThan(support[2][0]))) {
            // New top maximum OR same top maximum but closer to left.
            box.index[2] = i;
            support[2] = v;
        }

        if (v[0].lessThan(support[3][0]) ||
            (v[0].equals(support[3][0]) && v[1].lessThan(support[3][1]))) {
            // New left minimum OR same left minimum but closer to bottom.
            box.index[3] = i;
            support[3] = v;
        }
    }

    // The comment in the loop has the implication that support[0] = {0,0},
    // so the scaled height (support[2][1] - support[0][1]) is simply
    // support[2][1].
    const scaledWidth = support[1][0].sub(support[3][0]);
    const scaledHeight = support[2][1];
    box.area = scaledWidth.mul(scaledHeight).div(box.sqrLenU0);
    return box;
}

// Compute (sin(angle))^2 for the polygon edges emanating from the support
// vertices of the box. The returned array is empty when no angle is in
// [0,pi/2), in which case the original polygon must be a rectangle.
function computeAngles(vertices: readonly RationalPoint2[],
    box: Box): AnglePair[] {
    const numVertices = vertices.length;
    const A: AnglePair[] = [];
    for (let k0 = 3, k1 = 0; k1 < 4; k0 = k1++) {
        if (box.index[k0] !== box.index[k1]) {
            // The box edges are ordered in k0 as U[0], U[1], -U[0], -U[1].
            const D = ((k0 & 2) !== 0 ? rNegate(box.U[k0 & 1]) : box.U[k0 & 1]);
            const j0 = box.index[k0];
            let j1 = j0 + 1;
            if (j1 === numVertices) {
                j1 = 0;
            }
            const E = rSub(vertices[j1], vertices[j0]);
            const dp = rDotPerp(D, E);
            const esqrlen = rDot(E, E);
            A.push({ sinThetaSqr: dp.mul(dp).div(esqrlen), edge: k0 });
        }
    }
    return A;
}

// Sort the angles indirectly. The sorted indices are returned. This avoids
// swapping elements of A[], which is expensive for an exact rational type.
// The upstream sorting networks are replicated exactly, because the choice
// among tied minima determines which box edge becomes the new bottom edge.
function sortAngles(A: readonly AnglePair[]): [number, number, number, number] {
    const sort: [number, number, number, number] = [0, 1, 2, 3];
    const swap = (a: number, b: number): void => {
        if (A[sort[a]].sinThetaSqr.greaterThan(A[sort[b]].sinThetaSqr)) {
            const temp = sort[a];
            sort[a] = sort[b];
            sort[b] = temp;
        }
    };

    const numA = A.length;
    if (numA === 2) {
        swap(0, 1);
    }
    else if (numA === 3) {
        swap(0, 1);
        swap(0, 2);
        swap(1, 2);
    }
    else if (numA === 4) {
        swap(0, 1);
        swap(2, 3);
        swap(0, 2);
        swap(1, 3);
        swap(1, 2);
    }
    return sort;
}

// Update the supporting indices (box.index[]) and the box axis directions
// (box.U[]). The return value is false when the new bottom polygon edge has
// already been processed, at which point the search is over.
function updateSupport(A: readonly AnglePair[],
    sort: readonly [number, number, number, number],
    vertices: readonly RationalPoint2[], visited: boolean[], box: Box): boolean {
    // Replace the support vertices of those edges attaining minimum angle
    // with the other endpoints of the edges.
    const numVertices = vertices.length;
    const amin = A[sort[0]];
    for (let k = 0; k < A.length; ++k) {
        const a = A[sort[k]];
        if (a.sinThetaSqr.equals(amin.sinThetaSqr)) {
            if (++box.index[a.edge] === numVertices) {
                box.index[a.edge] = 0;
            }
        }
    }

    const bottom = box.index[amin.edge];
    if (visited[bottom]) {
        // This polygon edge has already been processed.
        return false;
    }
    visited[bottom] = true;

    // Cycle the vertices so that the bottom support occurs first.
    const nextIndex: [number, number, number, number] = [0, 0, 0, 0];
    for (let k = 0; k < 4; ++k) {
        nextIndex[k] = box.index[(amin.edge + k) % 4];
    }
    box.index = nextIndex;

    // Compute the box axis directions.
    const j1 = box.index[0];
    let j0 = j1 - 1;
    if (j0 < 0) {
        j0 = numVertices - 1;
    }
    box.U[0] = rSub(vertices[j1], vertices[j0]);
    box.U[1] = rNegPerp(box.U[0]);
    box.sqrLenU0 = rDot(box.U[0], box.U[0]);

    // Compute the box area.
    const diff0 = rSub(vertices[box.index[1]], vertices[box.index[3]]);
    const diff1 = rSub(vertices[box.index[2]], vertices[box.index[0]]);
    box.area = rDot(box.U[0], diff0).mul(rDot(box.U[1], diff1))
        .div(box.sqrLenU0);
    return true;
}

// The slow O(n^2) search.
function computeBoxForEdgeOrderNSqr(vertices: readonly RationalPoint2[]): Box {
    let minBox: Box | null = null;
    const numIndices = vertices.length;
    for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
        const box = smallestBox(i0, i1, vertices);
        if (minBox === null || box.area.lessThan(minBox.area)) {
            minBox = box;
        }
    }
    logAssert(minBox !== null, 'MinimumAreaBox2: the polygon has no edges.');
    return minBox;
}

// The fast O(n) search.
function computeBoxForEdgeOrderN(vertices: readonly RationalPoint2[]): Box {
    // The inputs are assumed to be the vertices of a convex polygon that is
    // counterclockwise ordered. The input points must not contain three
    // consecutive collinear points.

    // When the bounding box corresponding to a polygon edge is computed, the
    // edge is marked as visited. If the edge is encountered later, the
    // algorithm terminates.
    const visited: boolean[] = new Array<boolean>(vertices.length).fill(false);

    // Start the minimum-area rectangle search with the edge from the last
    // polygon vertex to the first. When updating the extremes, we want the
    // bottom-most point on the left edge, the top-most point on the right
    // edge, the left-most point on the top edge, and the right-most point on
    // the bottom edge. The polygon edges starting at these points are then
    // guaranteed not to coincide with a box edge except when an extreme
    // point is shared by two box edges (at a corner).
    let minBox = smallestBox(vertices.length - 1, 0, vertices);
    visited[minBox.index[0]] = true;

    // Execute the rotating calipers algorithm.
    const box = cloneBox(minBox);
    for (let i = 0; i < vertices.length; ++i) {
        const A = computeAngles(vertices, box);
        if (A.length === 0) {
            // The polygon is a rectangle, so the search is over.
            break;
        }

        // Indirectly sort the A-array.
        const sort = sortAngles(A);

        if (!updateSupport(A, sort, vertices, visited, box)) {
            // This box polygon edge has already been processed, so the
            // search is over.
            break;
        }

        if (box.area.lessThan(minBox.area)) {
            minBox = cloneBox(box);
        }
    }

    return minBox;
}

export class MinimumAreaBox2 {
    // The input points to be bound.
    private mNumPoints: number;
    private mPoints: readonly Vector[];

    // The indices into mPoints for the convex hull vertices.
    private mHull: number[];

    // The support indices for the minimum-area box. NOTE: these are lookups
    // into the polygon that the box search processed, which is the
    // collinear-point-free polygon derived from the hull (compute) or from
    // the caller's polygon (computeConvexPolygon). Because ConvexHull2
    // guarantees no three consecutive collinear hull points, for compute()
    // they are lookups into getHull(), so the supporting input points are
    // getPoints()[getHull()[getSupportIndices()[i]]].
    private mSupportIndices: [number, number, number, number];

    // The area of the minimum-area box. The rational value is exact, so the
    // only rounding error occurs in the conversion to floating point
    // (round-to-nearest-ties-to-even).
    private mArea: number;

    // The class is a functor to support computing the minimum-area box of
    // multiple data sets using the same class object.
    constructor() {
        this.mNumPoints = 0;
        this.mPoints = [];
        this.mHull = [];
        this.mSupportIndices = [0, 0, 0, 0];
        this.mArea = 0;
    }

    // The points are arbitrary, so the convex hull must be computed from
    // them in order to compute the minimum-area box. NOTE: ConvexHull2
    // guarantees that the hull does not have three consecutive collinear
    // points.
    compute(points: readonly Vector[],
        useRotatingCalipers: boolean = true): OrientedBox2 {
        // Upstream does not reset the support indices and the area on the
        // degenerate (dimension 0 or 1) paths, so getSupportIndices() and
        // getArea() then report stale values from an earlier call on the
        // same object. The port resets them for every query.
        this.mNumPoints = points.length;
        this.mPoints = points;
        this.mHull = [];
        this.mSupportIndices = [0, 0, 0, 0];
        this.mArea = 0;

        // Get the convex hull of the points.
        const hullQuery = new ConvexHull2();
        hullQuery.compute(points);
        const dimension = hullQuery.getDimension();

        const minBox = new OrientedBox(2);

        if (dimension === 0) {
            // The points are all the same.
            minBox.center = points[0].clone();
            minBox.axis[0] = Vector.unit(2, 0);
            minBox.axis[1] = Vector.unit(2, 1);
            minBox.extent = Vector.zero(2);
            this.mHull = [0];
            return minBox;
        }

        if (dimension === 1) {
            // The points lie on a line. Determine the extreme t-values for
            // the points represented as P = origin + t*direction. The
            // 'origin' is an input vertex, so t = 0 is attained.
            //
            // Upstream bug (MinimumAreaBox2.h, operator(), dimension == 1):
            // tmin and tmax start at zero with imin = imax = 0, on the
            // grounds that the line origin is an input vertex whose t-value
            // is zero. But the line origin is points[hull[0]], not
            // necessarily points[0], so the reported extreme indices are
            // wrong (and mHull is then wrong) unless hull[0] == 0. The port
            // seeds both extremes from point 0 instead, which yields the
            // same tmin and tmax and the correct indices.
            const line = hullQuery.getLine();
            let tmin = dot(sub(points[0], line.origin), line.direction);
            let tmax = tmin;
            let imin = 0;
            let imax = 0;
            for (let i = 1; i < points.length; ++i) {
                const t = dot(sub(points[i], line.origin), line.direction);
                if (t > tmax) {
                    tmax = t;
                    imax = i;
                }
                else if (t < tmin) {
                    tmin = t;
                    imin = i;
                }
            }

            const half = 0.5;
            minBox.center = Vector.fromArray([
                line.origin.get(0) + half * (tmin + tmax) * line.direction.get(0),
                line.origin.get(1) + half * (tmin + tmax) * line.direction.get(1)
            ]);
            minBox.extent = Vector.fromArray([half * (tmax - tmin), 0]);
            minBox.axis[0] = line.direction.clone();
            minBox.axis[1] = Vector.fromArray([
                -line.direction.get(1), line.direction.get(0)
            ]);
            this.mHull = [imin, imax];
            return minBox;
        }

        this.mHull = [...hullQuery.getHull()];
        const vertices = hullQuery.getPoints();
        const computePoints: RationalPoint2[] = this.mHull.map(h => [
            BSRational.fromNumber(vertices[h].get(0)),
            BSRational.fromNumber(vertices[h].get(1))
        ] as RationalPoint2);

        return this.computeFromPolygon(computePoints, useRotatingCalipers,
            minBox);
    }

    // The points already form a counterclockwise, nondegenerate convex
    // polygon. If the points directly are the convex polygon, pass an empty
    // 'indices' array. If the polygon vertices are a subset of the incoming
    // points, that subset is identified by an 'indices' array of at least 3
    // elements.
    computeConvexPolygon(points: readonly Vector[],
        indices: readonly number[],
        useRotatingCalipers: boolean = true): OrientedBox2 {
        this.mNumPoints = points.length;
        this.mPoints = points;
        this.mHull = [];
        this.mSupportIndices = [0, 0, 0, 0];
        this.mArea = 0;

        const minBox = new OrientedBox(2);

        if (points.length < 3 ||
            (indices.length > 0 && indices.length < 3)) {
            minBox.center = Vector.zero(2);
            minBox.axis[0] = Vector.unit(2, 0);
            minBox.axis[1] = Vector.unit(2, 1);
            minBox.extent = Vector.zero(2);
            return minBox;
        }

        if (indices.length > 0) {
            this.mHull = [...indices];
        }
        else {
            this.mHull = [];
            for (let i = 0; i < points.length; ++i) {
                this.mHull.push(i);
            }
        }

        const computePoints: RationalPoint2[] = this.mHull.map(h => [
            BSRational.fromNumber(points[h].get(0)),
            BSRational.fromNumber(points[h].get(1))
        ] as RationalPoint2);

        return this.computeFromPolygon(computePoints, useRotatingCalipers,
            minBox);
    }

    // Member access.
    getNumPoints(): number {
        return this.mNumPoints;
    }

    getPoints(): readonly Vector[] {
        return this.mPoints;
    }

    getHull(): readonly number[] {
        return this.mHull;
    }

    getSupportIndices(): readonly [number, number, number, number] {
        return this.mSupportIndices;
    }

    getArea(): number {
        return this.mArea;
    }

    private computeFromPolygon(computePoints: RationalPoint2[],
        useRotatingCalipers: boolean, minBox: OrientedBox2): OrientedBox2 {
        const vertices = removeCollinearPoints(computePoints);
        logAssert(vertices.length >= 3,
            'MinimumAreaBox2: the polygon must have at least 3 noncollinear vertices.');

        const box = useRotatingCalipers
            ? computeBoxForEdgeOrderN(vertices)
            : computeBoxForEdgeOrderNSqr(vertices);

        this.convertTo(box, vertices, minBox);
        return minBox;
    }

    // Convert the rational box to the floating-point box. The conversions
    // are deferred until the last step to avoid precision loss.
    private convertTo(box: Box, computePoints: readonly RationalPoint2[],
        itMinBox: OrientedBox2): void {
        // The sum, difference and center are all computed exactly.
        const sum: [RationalPoint2, RationalPoint2] = [
            rAdd(computePoints[box.index[1]], computePoints[box.index[3]]),
            rAdd(computePoints[box.index[2]], computePoints[box.index[0]])
        ];

        const difference: [RationalPoint2, RationalPoint2] = [
            rSub(computePoints[box.index[1]], computePoints[box.index[3]]),
            rSub(computePoints[box.index[2]], computePoints[box.index[0]])
        ];

        const dot0 = rDot(box.U[0], sum[0]);
        const dot1 = rDot(box.U[1], sum[1]);
        const center: RationalPoint2 = [
            R_HALF.mul(dot0.mul(box.U[0][0]).add(dot1.mul(box.U[1][0])))
                .div(box.sqrLenU0),
            R_HALF.mul(dot0.mul(box.U[0][1]).add(dot1.mul(box.U[1][1])))
                .div(box.sqrLenU0)
        ];

        // Calculate the squared extent using rational arithmetic to avoid
        // loss of precision before computing a square root.
        const sqrExtent: BSRational[] = [];
        for (let i = 0; i < 2; ++i) {
            const e = R_HALF.mul(rDot(box.U[i], difference[i]));
            sqrExtent.push(e.mul(e).div(box.sqrLenU0));
        }

        for (let i = 0; i < 2; ++i) {
            itMinBox.center.set(i, center[i].toNumber());
            itMinBox.extent.set(i, Math.sqrt(sqrExtent[i].toNumber()));

            // Before converting to floating point, factor out the maximum
            // component using rational arithmetic to generate rational
            // numbers in a range that avoids loss of precision during the
            // conversion and normalization.
            const axis = box.U[i];
            const abs0 = BSRational.fabs(axis[0]);
            const abs1 = BSRational.fabs(axis[1]);
            const cmax = abs0.lessThan(abs1) ? abs1 : abs0;
            const invCMax = R_ONE.div(cmax);
            itMinBox.axis[i] = Vector.fromArray([
                axis[0].mul(invCMax).toNumber(),
                axis[1].mul(invCMax).toNumber()
            ]);
            normalize(itMinBox.axis[i]);
        }

        this.mSupportIndices = [
            box.index[0], box.index[1], box.index[2], box.index[3]
        ];
        this.mArea = box.area.toNumber();
    }
}
