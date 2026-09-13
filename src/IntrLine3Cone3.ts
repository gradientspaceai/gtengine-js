// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the cone to be single sided and solid. The cone height
// range is [hmin,hmax]. The cone can be infinite where hmin = 0 and
// hmax = +infinity, infinite truncated where hmin > 0 and hmax = +infinity,
// finite where hmin = 0 and hmax < +infinity, or a cone frustum where
// hmin > 0 and hmax < +infinity. The algorithm details are found in
// https://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, which becomes IntrLine3Cone3FI. The static
// 'int32_t constexpr' result kinds (isEmpty, isPoint, ...) become the
// exported enum IntrLine3Cone3FIResultType, and the protected DoQuery,
// DoQuerySpecial, the Case* handlers and the Set* helpers become
// module-private functions that write into the result object.
//
// The query computes with QFNumber quadratic-field values x + y*sqrt(d) so
// that it is error-free for exact arithmetic types. A 3D point with
// quadratic-field components is the tuple type IntrLine3Cone3QFPoint rather
// than a Vector (the port's Vector stores numbers).
//
// The port fixes upstream's handling of the degenerate configuration where
// the line contains the cone vertex, which is result corrupting (gtengine-js
// issue #465 item 1 and its V35 addition). The quadratic then has a double
// root, so the discriminant is zero in exact arithmetic but a cancelling
// difference in floating point, and upstream's case analysis for the rounded
// branches assumes a nondegenerate configuration: it reports a point, the
// empty set, or a segment reaching the cone's maximum height where the answer
// is a ray or the vertex alone. The fix is in four places, each explained at
// its site: the vertex test hoisted to the top of doQuerySpecial, and the
// c2 > 0 cases of caseC2NotZeroDiscrNeg, caseC2NotZeroDiscrPos and
// caseC2NotZeroDiscrZero. A line that passes very close to, but not exactly
// through, the vertex remains ill conditioned: the discriminant has no
// significant digits there, and the port keeps upstream's answer.
//
// Upstream uses FIIntervalInterval<QFN1>, that is, the interval-interval
// find-intersection query instantiated for quadratic-field numbers. The port
// of IntrIntervals.ts is specialized to 'number', so the two interval
// overloads that this query needs are reimplemented here (module-private
// qfFindIntervals and qfFindFiniteSemiInfinite) as literal transcriptions of
// the IntrIntervals.ts logic with QFNumber comparisons.

import type { Cone3 } from './Cone.js';
import type { FIQuery } from './FIQuery.js';
import type { Line3 } from './Line.js';
import { logAssert, logError } from './Logger.js';
import { QFNumber } from './QFNumber.js';
import { Vector, dot, negate, sub } from './Vector.js';

// A 3D point whose components are quadratic-field numbers (the port of
// Vector3<QFNumber<Real,1>>).
export type IntrLine3Cone3QFPoint = [QFNumber, QFNumber, QFNumber];

// Because the intersection of a line and a cone with infinite height can be a
// ray or a line, the result has a 'type' value that says how to interpret the
// t[] and P[] values.
export enum IntrLine3Cone3FIResultType {
    // No intersection.
    isEmpty = 0,

    // t[0] is finite, t[1] is set to t[0], P[0] is the point of intersection,
    // P[1] is set to P[0].
    isPoint = 1,

    // t[0] and t[1] are finite with t[0] < t[1], P[0] and P[1] are the
    // endpoints of the segment of intersection.
    isSegment = 2,

    // Dot(line.direction, cone.ray.direction) > 0: t[0] is finite, t[1] is
    // +infinity (set to +1), P[0] is the ray origin, P[1] is the ray
    // direction (set to line.direction). The ray starts at P[0] and you walk
    // away from it in the line direction.
    isRayPositive = 3,

    // Dot(line.direction, cone.ray.direction) < 0: t[0] is -infinity (set to
    // -1), t[1] is finite, P[0] is the ray endpoint, P[1] is the ray
    // direction (set to line.direction). The ray ends at P[0] and you walk
    // towards it in the line direction.
    isRayNegative = 4
}

// The result of IntrLine3Cone3FI.find.
export interface IntrLine3Cone3FIResult {
    intersect: boolean;
    type: IntrLine3Cone3FIResultType;
    t: [QFNumber, QFNumber];
    P: [IntrLine3Cone3QFPoint, IntrLine3Cone3QFPoint];
}

function zeroQFPoint(): IntrLine3Cone3QFPoint {
    return [new QFNumber(), new QFNumber(), new QFNumber()];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Cone3FIResult(): IntrLine3Cone3FIResult {
    return {
        intersect: false,
        type: IntrLine3Cone3FIResultType.isEmpty,
        t: [new QFNumber(), new QFNumber()],
        P: [zeroQFPoint(), zeroQFPoint()]
    };
}

// The port of the templated 'Result::Convert' helpers, which upstream never
// instantiates (QFNumber has no conversion operator to Real, so the templates
// would not compile if they were used). The obvious intent is the numeric
// value x[0] + x[1] * sqrt(d).
export function intrLine3Cone3Convert(input: QFNumber): number {
    return (input.x[0] as number) + (input.x[1] as number) * Math.sqrt(input.d);
}

export function intrLine3Cone3ConvertPoint(input: IntrLine3Cone3QFPoint):
    Vector {
    return Vector.fromArray([
        intrLine3Cone3Convert(input[0]),
        intrLine3Cone3Convert(input[1]),
        intrLine3Cone3Convert(input[2])
    ]);
}

// The port of 'Result::ComputePoints'.
export function intrLine3Cone3ComputePoints(origin: Vector, direction: Vector,
    result: IntrLine3Cone3FIResult): void {
    switch (result.type) {
        case IntrLine3Cone3FIResultType.isEmpty:
            for (let i = 0; i < 3; ++i) {
                result.P[0][i] = new QFNumber();
                result.P[1][i] = result.P[0][i].clone();
            }
            break;
        case IntrLine3Cone3FIResultType.isPoint:
            for (let i = 0; i < 3; ++i) {
                result.P[0][i] = result.t[0].mul(direction.values[i])
                    .add(origin.values[i]);
                result.P[1][i] = result.P[0][i].clone();
            }
            break;
        case IntrLine3Cone3FIResultType.isSegment:
            for (let i = 0; i < 3; ++i) {
                result.P[0][i] = result.t[0].mul(direction.values[i])
                    .add(origin.values[i]);
                result.P[1][i] = result.t[1].mul(direction.values[i])
                    .add(origin.values[i]);
            }
            break;
        case IntrLine3Cone3FIResultType.isRayPositive:
            for (let i = 0; i < 3; ++i) {
                result.P[0][i] = result.t[0].mul(direction.values[i])
                    .add(origin.values[i]);
                result.P[1][i] = new QFNumber(direction.values[i], 0,
                    result.t[0].d);
            }
            break;
        case IntrLine3Cone3FIResultType.isRayNegative:
            for (let i = 0; i < 3; ++i) {
                result.P[0][i] = result.t[1].mul(direction.values[i])
                    .add(origin.values[i]);
                result.P[1][i] = new QFNumber(direction.values[i], 0,
                    result.t[1].d);
            }
            break;
        default:
            logError('Invalid case.');
            break;
    }
}

// Find-intersection query for a line and a solid cone in 3D.
export class IntrLine3Cone3FI implements
    FIQuery<Line3, Cone3, IntrLine3Cone3FIResult> {

    find(line: Line3, cone: Cone3): IntrLine3Cone3FIResult {
        logAssert(cone.dimension === 3 && line.origin.size === 3,
            'IntrLine3Cone3FI: mismatched sizes.');
        const result = defaultIntrLine3Cone3FIResult();
        intrLine3Cone3FIDoQuery(line.origin, line.direction, cone, result);
        intrLine3Cone3ComputePoints(line.origin, line.direction, result);
        result.intersect = (result.type !== IntrLine3Cone3FIResultType.isEmpty);
        return result;
    }
}

// The port of the protected 'FIQuery::DoQuery'. The result.type and result.t
// values are computed here; the result.P and result.intersect values are
// computed from them by 'find'.
export function intrLine3Cone3FIDoQuery(lineOrigin: Vector,
    lineDirection: Vector, cone: Cone3,
    result: IntrLine3Cone3FIResult): void {
    // The algorithm implemented in doQuerySpecial avoids extra branches if we
    // choose a line whose direction forms an acute angle with the cone
    // direction.
    if (dot(lineDirection, cone.ray.direction) >= 0) {
        doQuerySpecial(lineOrigin, lineDirection, cone, result);
    }
    else {
        doQuerySpecial(lineOrigin, negate(lineDirection), cone, result);
        const t0 = result.t[0].negate();
        const t1 = result.t[1].negate();
        result.t[0] = t1;
        result.t[1] = t0;
        if (result.type === IntrLine3Cone3FIResultType.isRayPositive) {
            result.type = IntrLine3Cone3FIResultType.isRayNegative;
        }
    }
}

function doQuerySpecial(lineOrigin: Vector, lineDirection: Vector,
    cone: Cone3, result: IntrLine3Cone3FIResult): void {
    // Compute the number of real-valued roots and represent them using
    // rational quadratic field elements to support exact rational arithmetic.
    const PmV = sub(lineOrigin, cone.ray.origin);
    const UdU = dot(lineDirection, lineDirection);
    const DdU = dot(cone.ray.direction, lineDirection);  // >= 0
    const DdPmV = dot(cone.ray.direction, PmV);
    const UdPmV = dot(lineDirection, PmV);
    const PmVdPmV = dot(PmV, PmV);
    const c2 = DdU * DdU - cone.cosAngleSqr * UdU;
    const c1 = DdU * DdPmV - cone.cosAngleSqr * UdPmV;
    const c0 = DdPmV * DdPmV - cone.cosAngleSqr * PmVdPmV;

    // Port fix for an upstream defect (gtengine-js issue #465 item 1 and its
    // V35 addition), result corrupting and therefore fixed rather than
    // preserved. When the line contains the cone vertex V, the quadratic
    // Q(t) = c2*t^2 + 2*c1*t + c0 is c2*(t-tv)^2, where tv is the parameter
    // of V: it has a double root, so discr = c1*c1 - c0*c2 is zero. Upstream
    // reaches the two blocks that handle this only when discr evaluates to
    // exactly zero, but discr is a difference of two nearly equal products,
    // so rounding routes the query into the discr > 0 or discr < 0 branches
    // instead and it reports a point, the empty set, or a segment reaching
    // the cone's maximum height in place of the correct set. The test below
    // is on the input itself: tv is the parameter of the point of the line
    // closest to V (no cancellation), and V is on the line exactly when
    // (P - V) + tv * U is the zero vector. The answer then follows from the
    // sign of c2 alone, with no discriminant: c2 > 0 means U is interior to
    // the cone's angular region, so Q >= 0 (the double-sided cone) holds for
    // every t and the positive cone is the ray of nonnegative heights from V;
    // c2 < 0 means U is exterior, so Q >= 0 holds only at V; c2 == 0 means
    // the line is a cone ruling through V, again a ray.
    const tv = -UdPmV / UdU;
    let onVertex = true;
    for (let i = 0; i < 3; ++i) {
        if (PmV.values[i] + tv * lineDirection.values[i] !== 0) {
            onVertex = false;
            break;
        }
    }

    if (onVertex) {
        if (c2 < 0) {
            // Block 4. The line is outside the double-sided cone and
            // intersects it only at V.
            setPointClamp(new QFNumber(tv, 0, 0), new QFNumber(0, 0, 0), cone,
                result);
        }
        else {
            // Block 5. The line is inside the double-sided cone (or on its
            // boundary when c2 == 0), so the intersection is a ray with
            // origin V.
            setRayClamp(new QFNumber(0, 0, 0), DdU, DdPmV, cone, result);
        }
        return;
    }

    if (c2 !== 0) {
        const discr = c1 * c1 - c0 * c2;
        if (discr < 0) {
            caseC2NotZeroDiscrNeg(c1, c2, DdU, DdPmV, cone, result);
        }
        else if (discr > 0) {
            caseC2NotZeroDiscrPos(c1, c2, discr, DdU, DdPmV, cone, result);
        }
        else {  // discr == 0
            caseC2NotZeroDiscrZero(c1, c2, PmV, lineDirection, DdU, DdPmV,
                cone, result);
        }
    }
    else if (c1 !== 0) {
        caseC2ZeroC1NotZero(c0, c1, DdU, DdPmV, cone, result);
    }
    else {
        caseC2ZeroC1Zero(c0, UdU, UdPmV, DdU, DdPmV, cone, result);
    }
}

// The intersection of the line with the positive cone when the quadratic has
// a double root at 't' -- or when its discriminant is the rounding of a
// double root -- and c2 > 0. The whole line is then inside the double-sided
// cone, so the positive cone contributes the ray of heights at or above the
// height of the double root. That height is zero (the vertex) in exact
// arithmetic; 'h' is its rounded value, clamped to be nonnegative.
function setRayClampAtDoubleRoot(t: number, DdU: number, DdPmV: number,
    cone: Cone3, result: IntrLine3Cone3FIResult): void {
    const h = t * DdU + DdPmV;
    setRayClamp(new QFNumber(h > 0 ? h : 0, 0, 0), DdU, DdPmV, cone, result);
}

function caseC2NotZeroDiscrNeg(c1: number, c2: number, DdU: number,
    DdPmV: number, cone: Cone3, result: IntrLine3Cone3FIResult): void {
    if (c2 < 0) {
        // Block 0. The quadratic has no real-valued roots. The line does not
        // intersect the double-sided cone.
        setEmpty(result);
    }
    else {
        // Port fix (see doQuerySpecial). With c2 > 0 the line direction is
        // interior to the cone's angular region, so the line runs into both
        // nappes and the quadratic has real roots: discr >= 0 in exact
        // arithmetic, with equality exactly when the line contains the
        // vertex. A negative discr is therefore the rounding of a double
        // root, and upstream's empty result discards a whole ray -- the cone
        // axis itself is reported as missing the cone.
        setRayClampAtDoubleRoot(-c1 / c2, DdU, DdPmV, cone, result);
    }
}

function caseC2NotZeroDiscrPos(c1: number, c2: number, discr: number,
    DdU: number, DdPmV: number, cone: Cone3,
    result: IntrLine3Cone3FIResult): void {
    // The quadratic has two distinct real-valued roots, t[0] and t[1] with
    // t[0] < t[1].
    const x = -c1 / c2;
    const y = (c2 > 0 ? 1 / c2 : -1 / c2);
    const t: [QFNumber, QFNumber] = [
        new QFNumber(x, -y, discr),
        new QFNumber(x, y, discr)
    ];

    // Compute the signed heights at the intersection points, h[0] and h[1]
    // with h[0] <= h[1]. The ordering is guaranteed because the input line
    // was arranged to satisfy Dot(D,U) >= 0.
    const h: [QFNumber, QFNumber] = [
        t[0].mul(DdU).add(DdPmV),
        t[1].mul(DdU).add(DdPmV)
    ];

    const zero = new QFNumber(0, 0, discr);
    if (h[0].greaterThanEqual(zero)) {
        // Block 1. The line intersects the positive cone in two points.
        setSegmentClamp(t, h, DdU, DdPmV, cone, result);
    }
    else if (h[1].lessThanEqual(zero)) {
        // Block 2. The line intersects the negative cone in two points.
        setEmpty(result);
    }
    else if (c2 > 0) {
        // Block 3. The line direction is interior to the cone's angular
        // region, so Q >= 0 (the double-sided cone) holds on the two rays
        // t <= t[0] and t >= t[1]; the line intersects the positive cone in
        // the ray of heights at or above h[1] and the negative cone in the
        // ray of heights at or below h[0].
        setRayClamp(h[1], DdU, DdPmV, cone, result);
    }
    else {
        // Port fix for an upstream defect (gtengine-js issue #465 item 1),
        // result corrupting and therefore fixed rather than preserved.
        // Upstream reaches SetRayClamp here for either sign of c2, but that
        // analysis holds only for c2 > 0. With c2 < 0 the set Q >= 0 is the
        // segment [t[0],t[1]], not a pair of rays, and the positive cone
        // takes the part of it with nonnegative height: the heights run from
        // zero up to h[1], not from h[1] to the cone's maximum height.
        //
        // Reaching this block needs the height to vanish strictly inside the
        // segment, where Q > 0; but where the height is zero Q is
        // -cosAngleSqr*|P-V|^2 <= 0, so it happens only for a line that
        // passes through (numerically, near) the cone vertex, whose double
        // root has been split by rounding of the discriminant. Upstream then
        // reported a segment reaching the cone's maximum height where the
        // true intersection is the vertex alone.
        setSegmentClamp(t, [zero, h[1]], DdU, DdPmV, cone, result);
    }
}

function caseC2NotZeroDiscrZero(c1: number, c2: number, PmV: Vector,
    lineDirection: Vector, DdU: number, DdPmV: number, cone: Cone3,
    result: IntrLine3Cone3FIResult): void {
    const t = -c1 / c2;

    // The line touches the double-sided cone at the single point
    // X = P + t * U. Determine whether that point is the cone vertex V, that
    // is, whether (P - V) + t * U is the zero vector.
    //
    // Port fix for an upstream bug. Upstream tests 't * UdU + UdPmV == 0',
    // which is only the U-component of that vector equation: it says that X
    // is the point of the line closest to V, not that X is V. When the line
    // is tangent to the cone at a point X != V whose parameter also minimizes
    // the distance to V (for example, the 45-degree cone about +z touched by
    // the line x = 1, z = 1 with direction (0,1,0)), upstream takes the
    // vertex branch and hands SetPointClamp the height 0 instead of the
    // height of X. The height range test then rejects a tangent point that is
    // inside a truncated cone's height range, or accepts one outside it. The
    // port tests the full vector equation, so such a configuration falls
    // through to the tangency branch below, which computes the correct
    // height.
    let onVertex = true;
    for (let i = 0; i < 3; ++i) {
        if (PmV.values[i] + t * lineDirection.values[i] !== 0) {
            onVertex = false;
            break;
        }
    }

    if (onVertex) {
        // To get here, it must be that V = P + (-c1/c2) * U, where U is not
        // necessarily a unit-length vector. The line intersects the cone
        // vertex.
        if (c2 < 0) {
            // Block 4. The line is outside the double-sided cone and
            // intersects it only at V.
            setPointClamp(new QFNumber(t, 0, 0), new QFNumber(0, 0, 0), cone,
                result);
        }
        else {
            // Block 5. The line is inside the double-sided cone, so the
            // intersection is a ray with origin V.
            setRayClamp(new QFNumber(0, 0, 0), DdU, DdPmV, cone, result);
        }
    }
    else if (c2 > 0) {
        // Port fix (see doQuerySpecial). A line with c2 > 0 that meets the
        // double-sided cone at a single point meets it at the vertex: the
        // whole line satisfies Q >= 0, and a tangency at a point X != V would
        // put the line in the tangent plane at X, hence along a ruling, which
        // gives c2 == 0. So this is the vertex configuration with the double
        // root rounded away from the vertex, and the intersection is a ray,
        // not the tangency point of Block 6.
        setRayClampAtDoubleRoot(t, DdU, DdPmV, cone, result);
    }
    else {
        // The line is tangent to the cone at a point different from the
        // vertex.
        const h = t * DdU + DdPmV;
        if (h >= 0) {
            // Block 6. The line is tangent to the positive cone.
            setPointClamp(new QFNumber(t, 0, 0), new QFNumber(h, 0, 0), cone,
                result);
        }
        else {
            // Block 7. The line is tangent to the negative cone.
            setEmpty(result);
        }
    }
}

function caseC2ZeroC1NotZero(c0: number, c1: number, DdU: number,
    DdPmV: number, cone: Cone3, result: IntrLine3Cone3FIResult): void {
    // U is a direction vector on the cone boundary. Compute the t-value for
    // the intersection point and compute the corresponding height h to
    // determine whether that point is on the positive cone or negative cone.
    const t = -0.5 * c0 / c1;
    const h = t * DdU + DdPmV;
    if (h > 0) {
        // Block 8. The line intersects the positive cone and the ray of
        // intersection is interior to the positive cone. The intersection is
        // a ray or segment.
        setRayClamp(new QFNumber(h, 0, 0), DdU, DdPmV, cone, result);
    }
    else {
        // Block 9. The line intersects the negative cone and the ray of
        // intersection is interior to the negative cone.
        setEmpty(result);
    }
}

function caseC2ZeroC1Zero(c0: number, UdU: number, UdPmV: number, DdU: number,
    DdPmV: number, cone: Cone3, result: IntrLine3Cone3FIResult): void {
    if (c0 !== 0) {
        // Block 10. The line does not intersect the double-sided cone.
        setEmpty(result);
    }
    else {
        // Block 11. The line is on the cone boundary. The intersection with
        // the positive cone is a ray that contains the cone vertex. The
        // intersection is either a ray or segment.
        const t = -UdPmV / UdU;
        const h = t * DdU + DdPmV;
        setRayClamp(new QFNumber(h, 0, 0), DdU, DdPmV, cone, result);
    }
}

function setEmpty(result: IntrLine3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isEmpty;
    result.t[0] = new QFNumber();
    result.t[1] = new QFNumber();
}

function setPoint(t: QFNumber, result: IntrLine3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isPoint;
    result.t[0] = t;
    result.t[1] = result.t[0].clone();
}

function setSegment(t0: QFNumber, t1: QFNumber,
    result: IntrLine3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isSegment;
    result.t[0] = t0;
    result.t[1] = t1;
}

function setRayPositive(t: QFNumber, result: IntrLine3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isRayPositive;
    result.t[0] = t;
    result.t[1] = new QFNumber(+1, 0, t.d);  // +infinity
}

function setPointClamp(t: QFNumber, h: QFNumber, cone: Cone3,
    result: IntrLine3Cone3FIResult): void {
    if (cone.heightInRange(h.x[0] as number)) {
        // P0.
        setPoint(t, result);
    }
    else {
        // P1.
        setEmpty(result);
    }
}

function setSegmentClamp(t: readonly [QFNumber, QFNumber],
    h: readonly [QFNumber, QFNumber], DdU: number, DdPmV: number,
    cone: Cone3, result: IntrLine3Cone3FIResult): void {
    const hrange: [QFNumber, QFNumber] = [
        new QFNumber(cone.getMinHeight(), 0, h[0].d),
        new QFNumber(cone.getMaxHeight(), 0, h[0].d)
    ];

    if (h[1].greaterThan(h[0])) {
        const iir = (cone.isFinite()
            ? qfFindIntervals(h, hrange)
            : qfFindFiniteSemiInfinite(h, hrange[0], true));
        if (iir.numIntersections === 2) {
            // S0.
            setSegment(iir.overlap[0].sub(DdPmV).div(DdU),
                iir.overlap[1].sub(DdPmV).div(DdU), result);
        }
        else if (iir.numIntersections === 1) {
            // S1.
            setPoint(iir.overlap[0].sub(DdPmV).div(DdU), result);
        }
        else {  // iir.numIntersections == 0
            // S2.
            setEmpty(result);
        }
    }
    else {  // h[1] == h[0]
        if (hrange[0].lessThanEqual(h[0])
            && (cone.isFinite() ? h[0].lessThanEqual(hrange[1]) : true)) {
            // S3. DdU > 0 and the line is not perpendicular to the cone axis.
            setSegment(t[0], t[1], result);
        }
        else {
            // S4. DdU == 0 and the line is perpendicular to the cone axis.
            setEmpty(result);
        }
    }
}

function setRayClamp(h: QFNumber, DdU: number, DdPmV: number, cone: Cone3,
    result: IntrLine3Cone3FIResult): void {
    const hrange: [QFNumber, QFNumber] = [
        new QFNumber(cone.getMinHeight(), 0, h.d),
        new QFNumber(cone.getMaxHeight(), 0, h.d)
    ];

    if (cone.isFinite()) {
        const iir = qfFindFiniteSemiInfinite(hrange, h, true);
        if (iir.numIntersections === 2) {
            // R0.
            setSegment(iir.overlap[0].sub(DdPmV).div(DdU),
                iir.overlap[1].sub(DdPmV).div(DdU), result);
        }
        else if (iir.numIntersections === 1) {
            // R1.
            setPoint(iir.overlap[0].sub(DdPmV).div(DdU), result);
        }
        else {  // iir.numIntersections == 0
            // R2.
            setEmpty(result);
        }
    }
    else {
        // R3.
        const maxValue = (hrange[0].lessThan(h) ? h : hrange[0]);
        setRayPositive(maxValue.sub(DdPmV).div(DdU), result);
    }
}

// The QFNumber specializations of the two IntrIntervals.ts find-intersection
// overloads used above. Only 'numIntersections' and 'overlap' are consumed by
// this file, so only those are returned.
interface QFIntervalResult {
    numIntersections: number;
    overlap: [QFNumber, QFNumber];
}

function qfFindIntervals(interval0: readonly [QFNumber, QFNumber],
    interval1: readonly [QFNumber, QFNumber]): QFIntervalResult {
    const result: QFIntervalResult = {
        numIntersections: 0,
        overlap: [new QFNumber(), new QFNumber()]
    };

    if (interval0[1].lessThan(interval1[0])
        || interval0[0].greaterThan(interval1[1])) {
        result.numIntersections = 0;
    }
    else if (interval0[1].greaterThan(interval1[0])) {
        if (interval0[0].lessThan(interval1[1])) {
            result.overlap[0] = (interval0[0].lessThan(interval1[0])
                ? interval1[0] : interval0[0]);
            result.overlap[1] = (interval0[1].greaterThan(interval1[1])
                ? interval1[1] : interval0[1]);
            result.numIntersections =
                (result.overlap[0].lessThan(result.overlap[1]) ? 2 : 1);
        }
        else {  // interval0[0] == interval1[1]
            result.numIntersections = 1;
            result.overlap[0] = interval0[0];
            result.overlap[1] = result.overlap[0].clone();
        }
    }
    else {  // interval0[1] == interval1[0]
        result.numIntersections = 1;
        result.overlap[0] = interval0[1];
        result.overlap[1] = result.overlap[0].clone();
    }

    return result;
}

function qfFindFiniteSemiInfinite(finite: readonly [QFNumber, QFNumber],
    a: QFNumber, isPositiveInfinite: boolean): QFIntervalResult {
    const result: QFIntervalResult = {
        numIntersections: 0,
        overlap: [new QFNumber(), new QFNumber()]
    };

    if (isPositiveInfinite) {
        if (finite[1].greaterThan(a)) {
            result.overlap[0] = (finite[0].lessThan(a) ? a : finite[0]);
            result.overlap[1] = finite[1];
            result.numIntersections =
                (result.overlap[0].lessThan(result.overlap[1]) ? 2 : 1);
        }
        else if (finite[1].equals(a)) {
            result.numIntersections = 1;
            result.overlap[0] = a;
            result.overlap[1] = result.overlap[0].clone();
        }
        else {
            result.numIntersections = 0;
        }
    }
    else {  // is negative-infinite
        if (finite[0].lessThan(a)) {
            result.overlap[0] = finite[0];
            result.overlap[1] = (finite[1].greaterThan(a) ? a : finite[1]);
            result.numIntersections =
                (result.overlap[0].lessThan(result.overlap[1]) ? 2 : 1);
        }
        else if (finite[0].equals(a)) {
            result.numIntersections = 1;
            result.overlap[0] = a;
            result.overlap[1] = result.overlap[0].clone();
        }
        else {
            result.numIntersections = 0;
        }
    }

    return result;
}
