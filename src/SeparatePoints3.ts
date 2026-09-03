// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SeparatePoints3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Separate two point sets, if possible, by computing a plane for which the
// point sets lie on opposite sides. The algorithm computes the convex hull
// of the point sets, then uses the method of separating axes to determine
// whether the two convex polyhedra are disjoint.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
//
// Port notes:
// * Upstream templates the class on <Real, ComputeType>, but ComputeType is
//   never used: the current ConvexHull3 is templated only on the
//   floating-point type and performs its own exact predicates internally.
//   The port drops the vestigial parameter, matching the SeparatePoints2
//   port.
// * operator() returns a bool and writes the separating plane through an
//   output reference. The port's compute(...) returns
//   { separated, separatingPlane }. Upstream leaves 'separatingPlane'
//   holding the last candidate plane tested when the return value is false;
//   the port returns a plane in the same state, but the field is meaningful
//   only when 'separated' is true.
// * The code assumes each point set has at least 4 noncoplanar points. When
//   either hull has dimension other than 3, the query reports "not
//   separated", as upstream does.
// * Upstream builds std::set<std::pair<size_t,size_t>> edge sets whose
//   iteration order is lexicographic on the ordered pair (the reversed pair
//   is a distinct element). The port replicates that with sorted arrays of
//   ordered pairs so the reported separating plane does not depend on
//   container ordering.
//
// Upstream bug (SeparatePoints3.h, OnSameSide and WhichSide): this is the
// 3D instance of the SeparatePoints2 defect (see that file and its issue).
// The side classification compares Dot(plane.normal, P) to plane.constant
// with exact floating-point comparisons, where the plane is built by
// Plane3<Real>({P0,P1,P2}), whose normal is UnitCross(P1-P0, P2-P0) and
// whose constant is Dot(normal, P0); both are rounded. The vertices of the
// hull that owns the candidate face should classify as "on the plane" or
// "negative side" (the normal of a counterclockwise hull face points
// outward), but round-off makes vertices of that very face classify as
// strictly positive. WhichSide then returns +1 instead of -1 for that hull,
// so 'side0 * side1 <= 0' succeeds and the query reports a separation for
// point sets that plainly overlap. The same round-off affects the
// cross-product plane loop, where the plane normal UnitCross(diff0, diff1)
// and the constant Dot(normal, P) are independently rounded so that P
// itself need not satisfy Dot(normal, P) == constant.
//
// The port fixes this by evaluating the side of a point with an exact
// orientation predicate: for the plane through P0 with (unnormalized)
// normal N, the side of Q is the sign of Dot(N, Q - P0), which has the same
// sign as the upstream expression in exact arithmetic. This follows the
// ConvexHull3 (B118) precedent of computing the geometric predicates
// exactly. The separating plane that is returned is still the normalized
// floating-point plane of upstream.
//
// Upstream quirk (fixed in the port): in the cross-product loop upstream
// assigns separatingPlane.normal and separatingPlane.constant directly,
// which leaves separatingPlane.origin holding the origin of whichever face
// plane was constructed last. Hyperplane documents that the members cannot
// be set independently. The port builds the plane with
// Hyperplane.fromNormalConstant so that origin = constant * normal is
// consistent with the other two members.

import { BSNumber } from './BSNumber.js';
import { ConvexHull3 } from './ConvexHull3.js';
import { Hyperplane, type Plane3 } from './Hyperplane.js';
import { Vector, dot, sub } from './Vector.js';
import { unitCross } from './Vector3.js';

// The result of the separation query. 'separated' is upstream's boolean
// return value: true if and only if 'separatingPlane' separates the two
// point sets.
export interface SeparatePoints3Result {
    separated: boolean;
    separatingPlane: Plane3;
}

// The exact port of Vector3<Real> coordinates for the orientation predicate.
type ExactPoint3 = [BSNumber, BSNumber, BSNumber];

// A lazily filled cache of the exact coordinates of an input point set. Only
// the hull vertices are ever needed, so the conversions are memoized in the
// style of ConvexHull3::mRPoints.
class ExactPoints {
    private mPoints: readonly Vector[];
    private mExact: (ExactPoint3 | null)[];

    constructor(points: readonly Vector[]) {
        this.mPoints = points;
        this.mExact = new Array<ExactPoint3 | null>(points.length).fill(null);
    }

    get(index: number): ExactPoint3 {
        let exact = this.mExact[index];
        if (exact === null) {
            const v = this.mPoints[index].values;
            exact = [
                BSNumber.fromNumber(v[0]),
                BSNumber.fromNumber(v[1]),
                BSNumber.fromNumber(v[2])
            ];
            this.mExact[index] = exact;
        }
        return exact;
    }
}

function exactSub(p1: ExactPoint3, p0: ExactPoint3): ExactPoint3 {
    return [p1[0].sub(p0[0]), p1[1].sub(p0[1]), p1[2].sub(p0[2])];
}

function exactCross(u: ExactPoint3, v: ExactPoint3): ExactPoint3 {
    return [
        u[1].mul(v[2]).sub(u[2].mul(v[1])),
        u[2].mul(v[0]).sub(u[0].mul(v[2])),
        u[0].mul(v[1]).sub(u[1].mul(v[0]))
    ];
}

// A candidate plane in exact arithmetic: the plane through 'origin' with the
// (not necessarily unit-length) normal 'normal'.
interface ExactPlane {
    normal: ExactPoint3;
    origin: ExactPoint3;
}

// The sign of Dot(normal, q - origin), which is negative when q is on the
// inner side of an outward-pointing hull face normal, positive when q is on
// the outer side, and zero when q is on the plane.
function sideSign(plane: ExactPlane, q: ExactPoint3): number {
    const d = exactSub(q, plane.origin);
    const n = plane.normal;
    return n[0].mul(d[0]).add(n[1].mul(d[1])).add(n[2].mul(d[2])).getSign();
}

// Test whether all hull vertices are on the same side of the plane. The
// return value is 0 when the plane splits the point set, +1 when the points
// are on the positive side and -1 otherwise. A hull whose vertices all lie
// exactly on the plane reports -1, matching upstream.
function onSameSide(plane: ExactPlane, hull: readonly number[],
    points: ExactPoints): number {
    let posSide = 0;
    let negSide = 0;

    for (let t = 0; t < hull.length; ++t) {
        const s = sideSign(plane, points.get(hull[t]));
        if (s > 0) {
            ++posSide;
        }
        else if (s < 0) {
            ++negSide;
        }

        if (posSide !== 0 && negSide !== 0) {
            // The plane splits the point set.
            return 0;
        }
    }

    return posSide !== 0 ? +1 : -1;
}

// Establish which side of the plane the hull is on. The return value is 0
// when the hull is effectively coplanar with the plane.
function whichSide(plane: ExactPlane, hull: readonly number[],
    points: ExactPoints): number {
    for (let t = 0; t < hull.length; ++t) {
        const s = sideSign(plane, points.get(hull[t]));
        if (s > 0) {
            // The hull is on the positive side.
            return +1;
        }
        if (s < 0) {
            // The hull is on the negative side.
            return -1;
        }
    }

    // The hull is effectively coplanar.
    return 0;
}

// The port of std::set<std::pair<size_t,size_t>>: the distinct ordered index
// pairs of the hull triangles, in lexicographic order.
function buildEdgeSet(hull: readonly number[]): [number, number][] {
    const map = new Map<string, [number, number]>();
    const numTriangles = Math.floor(hull.length / 3);
    for (let i = 0; i < numTriangles; ++i) {
        const i0 = hull[3 * i];
        const i1 = hull[3 * i + 1];
        const i2 = hull[3 * i + 2];
        for (const pair of [[i0, i1], [i0, i2], [i1, i2]] as [number, number][]) {
            map.set(`${pair[0]},${pair[1]}`, pair);
        }
    }

    const edges = Array.from(map.values());
    edges.sort((e0, e1) => (e0[0] !== e1[0] ? e0[0] - e1[0] : e0[1] - e1[1]));
    return edges;
}

export class SeparatePoints3 {
    // The return value 'separated' is true if and only if there is a
    // separation. If true, the returned plane is a separating plane. The
    // code assumes that each point set has at least 4 noncoplanar points.
    compute(points0: readonly Vector[], points1: readonly Vector[]):
        SeparatePoints3Result {
        let separatingPlane = new Hyperplane(3);

        // Construct the convex hull of point set 0.
        const ch0 = new ConvexHull3();
        ch0.compute(points0);
        if (ch0.getDimension() !== 3) {
            return { separated: false, separatingPlane };
        }

        // Construct the convex hull of point set 1.
        const ch1 = new ConvexHull3();
        ch1.compute(points1);
        if (ch1.getDimension() !== 3) {
            return { separated: false, separatingPlane };
        }

        const hull0 = ch0.getHull();
        const hull1 = ch1.getHull();
        const numTriangles0 = Math.floor(hull0.length / 3);
        const numTriangles1 = Math.floor(hull1.length / 3);
        const exact0 = new ExactPoints(points0);
        const exact1 = new ExactPoints(points1);

        // Test the faces of hull 0 for possible separation of the points.
        for (let i = 0; i < numTriangles0; ++i) {
            // Look up the face (assert: i0 != i1 && i0 != i2 && i1 != i2).
            const i0 = hull0[3 * i];
            const i1 = hull0[3 * i + 1];
            const i2 = hull0[3 * i + 2];

            // Compute the potential separating plane
            // (assert: normal != (0,0,0)).
            const plane = facePlane(exact0, i0, i1, i2);

            // Determine whether hull 1 is on the same side of the plane.
            const side1 = onSameSide(plane, hull1, exact1);

            if (side1 !== 0) {
                // Determine on which side of the plane hull 0 lies.
                const side0 = whichSide(plane, hull0, exact0);
                if (side0 * side1 <= 0) {
                    // The plane separates the hulls.
                    separatingPlane = Hyperplane.fromPoints(
                        [points0[i0], points0[i1], points0[i2]]);
                    return { separated: true, separatingPlane };
                }
            }
        }

        // Test the faces of hull 1 for possible separation of the points.
        for (let i = 0; i < numTriangles1; ++i) {
            // Look up the face (assert: i0 != i1 && i0 != i2 && i1 != i2).
            const i0 = hull1[3 * i];
            const i1 = hull1[3 * i + 1];
            const i2 = hull1[3 * i + 2];

            // Compute the perpendicular to the face
            // (assert: normal != (0,0,0)).
            const plane = facePlane(exact1, i0, i1, i2);

            // Determine whether hull 0 is on the same side of the plane.
            const side0 = onSameSide(plane, hull0, exact0);
            if (side0 !== 0) {
                // Determine on which side of the plane hull 1 lies.
                const side1 = whichSide(plane, hull1, exact1);
                if (side0 * side1 <= 0) {
                    // The plane separates the hulls.
                    separatingPlane = Hyperplane.fromPoints(
                        [points1[i0], points1[i1], points1[i2]]);
                    return { separated: true, separatingPlane };
                }
            }
        }

        // Build the edge sets for the hulls.
        const edgeSet0 = buildEdgeSet(hull0);
        const edgeSet1 = buildEdgeSet(hull1);

        // Test the planes whose normals are cross products of two edges, one
        // from each hull.
        for (const e0 of edgeSet0) {
            // Get the edge.
            const origin = exact0.get(e0[0]);
            const diff0 = exactSub(exact0.get(e0[1]), origin);

            for (const e1 of edgeSet1) {
                const diff1 = exactSub(exact1.get(e1[1]), exact1.get(e1[0]));

                // Compute the potential separating plane. The exact normal is
                // Cross(diff0, diff1), a positive multiple of the upstream
                // UnitCross(diff0, diff1) whenever the cross product is not
                // the zero vector. When it is the zero vector, every point
                // classifies as "on the plane", which matches the upstream
                // behavior of a zero normal and a zero constant.
                const plane: ExactPlane = {
                    normal: exactCross(diff0, diff1),
                    origin
                };

                // Determine whether the hulls are on the same side of the
                // plane.
                const side0 = onSameSide(plane, hull0, exact0);
                const side1 = onSameSide(plane, hull1, exact1);
                if (side0 * side1 < 0) {
                    // The plane separates the hulls.
                    const normal = unitCross(
                        sub(points0[e0[1]], points0[e0[0]]),
                        sub(points1[e1[1]], points1[e1[0]]));
                    separatingPlane = Hyperplane.fromNormalConstant(
                        normal, dot(normal, points0[e0[0]]));
                    return { separated: true, separatingPlane };
                }
            }
        }

        return { separated: false, separatingPlane };
    }
}

// The exact plane of the hull face <i0,i1,i2>, whose normal is the positive
// multiple Cross(P1 - P0, P2 - P0) of the upstream unit-length normal.
function facePlane(points: ExactPoints, i0: number, i1: number,
    i2: number): ExactPlane {
    const p0 = points.get(i0);
    const edge0 = exactSub(points.get(i1), p0);
    const edge1 = exactSub(points.get(i2), p0);
    return { normal: exactCross(edge0, edge1), origin: p0 };
}
