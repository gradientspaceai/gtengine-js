// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistCircle3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The 3D circle-circle distance algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
// The notation used in the code matches that of the document.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Circle3<T>, Circle3<T>>' becomes the
// class DistCircle3Circle3 with the result type DistCircle3Circle3Result.
//
// - Upstream builds the polynomials with 'Polynomial1<Rational>', where
//   Rational is BSRational<UIntegerAP32>. The port's Polynomial1 is
//   number-valued only, so the small set of rational polynomial operations
//   needed here (add, subtract, multiply, scalar multiply, evaluate) is
//   implemented as module-private functions on BSRational[] coefficient
//   arrays, listed in increasing order of power. The behavior matches
//   Polynomial1: addition and subtraction eliminate high-order zero
//   coefficients, multiplication does not.
// - The private class SCPolynomial in the upstream header is dead code (no
//   member of the class references it), so it is not ported.
// - The private helpers PrepareCircles and DoQueryParallelPlanes become
//   module-private functions.

import { AxisAngle } from './AxisAngle';
import { BSRational } from './BSRational';
import { Circle3 } from './Circle3';
import type { DCPQuery } from './DCPQuery';
import { logAssert } from './Logger';
import { Matrix, mulMatrix } from './Matrix';
import { RootsGeneralPolynomial } from './RootsGeneralPolynomial';
import { Rotation } from './Rotation';
import {
    Vector, add, div, dot, getOrthogonal, length, mul, negate, normalize, sub
} from './Vector';
import { cross, unitCross } from './Vector3';

export interface DistCircle3Circle3Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, either 1 or 2.
    numClosestPairs: number;

    // circle0Closest[i] and circle1Closest[i] are the i-th pair of closest
    // points, for 0 <= i < numClosestPairs.
    circle0Closest: [Vector, Vector];
    circle1Closest: [Vector, Vector];

    // 'true' when the reported closest pair is only one representative of
    // infinitely many equidistant pairs (concentric circles, for example).
    equidistant: boolean;
}

// ---------------------------------------------------------------------------
// Rational polynomial arithmetic (see the port notes). A polynomial is stored
// as an array of coefficients in increasing order of power; the array is
// never empty.
// ---------------------------------------------------------------------------

type RPoly = BSRational[];

function r(x: number): BSRational {
    return BSRational.fromNumber(x);
}

// The port of Polynomial1::EliminateLeadingZeros.
function rpTrim(p: RPoly): RPoly {
    let leading = p.length - 1;
    while (leading > 0 && p[leading].getSign() === 0) {
        --leading;
    }
    return p.slice(0, leading + 1);
}

function rpAdd(p0: RPoly, p1: RPoly): RPoly {
    const n = Math.max(p0.length, p1.length);
    const result: RPoly = new Array<BSRational>(n);
    for (let i = 0; i < n; ++i) {
        if (i < p0.length && i < p1.length) {
            result[i] = p0[i].add(p1[i]);
        } else {
            result[i] = (i < p0.length ? p0[i] : p1[i]).clone();
        }
    }
    return rpTrim(result);
}

function rpSub(p0: RPoly, p1: RPoly): RPoly {
    const n = Math.max(p0.length, p1.length);
    const result: RPoly = new Array<BSRational>(n);
    for (let i = 0; i < n; ++i) {
        if (i < p0.length && i < p1.length) {
            result[i] = p0[i].sub(p1[i]);
        } else if (i < p0.length) {
            result[i] = p0[i].clone();
        } else {
            result[i] = p1[i].negated();
        }
    }
    return rpTrim(result);
}

// As upstream's Polynomial1::operator*, no high-order zeros are eliminated.
function rpMul(p0: RPoly, p1: RPoly): RPoly {
    const n = p0.length + p1.length - 1;
    const result: RPoly = new Array<BSRational>(n);
    for (let i = 0; i < n; ++i) {
        result[i] = r(0);
    }
    for (let i0 = 0; i0 < p0.length; ++i0) {
        for (let i1 = 0; i1 < p1.length; ++i1) {
            result[i0 + i1] = result[i0 + i1].add(p0[i0].mul(p1[i1]));
        }
    }
    return result;
}

// As upstream's Polynomial1::operator*(scalar), no high-order zeros are
// eliminated.
function rpScale(scalar: BSRational, p: RPoly): RPoly {
    return p.map(c => scalar.mul(c));
}

// Horner's method (the port of Polynomial1::operator()).
function rpEval(p: RPoly, x: BSRational): BSRational {
    let result = p[p.length - 1];
    for (let i = p.length - 2; i >= 0; --i) {
        result = result.mul(x).add(p[i]);
    }
    return result;
}

// The polynomial is identically zero. This is the port of the upstream test
// '!(poly.GetDegree() > 0 || poly[0].GetSign() != 0)'.
function rpIsZero(p: RPoly): boolean {
    return p.length === 1 && p[0].getSign() === 0;
}

// The port of the std::set<Rational> uniqueness filtering of the roots. The
// solver returns the roots in increasing order, so a linear pass suffices and
// the increasing order of std::set is preserved.
function solveUnique(p: RPoly): BSRational[] {
    const roots = RootsGeneralPolynomial.solveRational(p, true);
    const unique: BSRational[] = [];
    for (const root of roots) {
        if (unique.length === 0 ||
            unique[unique.length - 1].notEquals(root)) {
            unique.push(root);
        }
    }
    return unique;
}

// The rational cross product of two 3-tuples.
function rCross(u: readonly BSRational[], v: readonly BSRational[]):
    BSRational[] {
    return [
        u[1].mul(v[2]).sub(u[2].mul(v[1])),
        u[2].mul(v[0]).sub(u[0].mul(v[2])),
        u[0].mul(v[1]).sub(u[1].mul(v[0]))
    ];
}

function rDot(u: readonly BSRational[], v: readonly BSRational[]): BSRational {
    return u[0].mul(v[0]).add(u[1].mul(v[1])).add(u[2].mul(v[2]));
}

// ---------------------------------------------------------------------------

// The port of the private struct ClosestInfo.
interface ClosestInfo {
    sqrDistance: number;
    circle0Closest: Vector;
    circle1Closest: Vector;
    equidistant: boolean;
}

// The port of PrepareCircles. Upstream writes the transformed circles and the
// transformation through reference parameters and returns the swap flag; the
// port returns all of them.
function prepareCircles(inCircle0: Circle3, inCircle1: Circle3): {
    circle0: Circle3, circle1: Circle3, rotate: Matrix, translate: Vector,
    scale: number, swapped: boolean
} {
    // Order the circles so that circle1.radius has the larger radius of the
    // two circles.
    const swapped = (inCircle0.radius > inCircle1.radius);
    const circle0 = (swapped ? inCircle1 : inCircle0).clone();
    const circle1 = (swapped ? inCircle0 : inCircle1).clone();

    // Ensure both circles have normals with z-value in [0,1].
    if (circle0.normal.values[2] < 0) {
        circle0.normal = negate(circle0.normal);
    }
    if (circle1.normal.values[2] < 0) {
        circle1.normal = negate(circle1.normal);
    }

    // Apply a translation, rotation, and uniform scaling so that
    // circle1.center = (0,0,0), circle1.normal = (0,0,1), and
    // circle1.radius = 1. A consequence is that circle0.radius <= 1.
    const aa = new AxisAngle(
        unitCross(circle1.normal, Vector.unit(3, 2)),
        Math.acos(circle1.normal.values[2]));
    let rotate = Rotation.fromAxisAngle(aa).toMatrix();
    const translate = negate(circle1.center);
    const scale = 1 / circle1.radius;

    circle0.center = add(circle0.center, translate);
    circle0.center = mul(scale, mulMatrix(rotate, circle0.center) as Vector);
    circle0.normal = mulMatrix(rotate, circle0.normal) as Vector;
    circle0.radius *= scale;
    circle1.center = new Vector(3);
    circle1.normal = Vector.unit(3, 2);
    circle1.radius = 1;

    // Rotate about circle1.normal to transform circle0.center to (0,k1,k2);
    // that is, the x-component is 0.
    if (circle0.center.values[0] !== 0) {
        const c0 = circle0.center.values[0];
        const c1 = circle0.center.values[1];
        const len = Math.sqrt(c0 * c0 + c1 * c1);
        const sn = c0 / len;
        const cs = c1 / len;
        const rot1 = new Matrix(3, 3);
        rot1.setCol(0, Vector.fromArray([cs, sn, 0]));
        rot1.setCol(1, Vector.fromArray([-sn, cs, 0]));
        rot1.setCol(2, Vector.fromArray([0, 0, 1]));
        circle0.center = mulMatrix(rot1, circle0.center) as Vector;
        circle0.center.values[0] = 0;
        circle0.normal = mulMatrix(rot1, circle0.normal) as Vector;
        rotate = mulMatrix(rot1, rotate) as Matrix;
    }

    return { circle0, circle1, rotate, translate, scale, swapped };
}

// The port of DoQueryParallelPlanes. The two circles are in parallel planes
// where D = C1 - C0, the difference of circle centers.
function doQueryParallelPlanes(circle0: Circle3, circle1: Circle3, D: Vector,
    result: DistCircle3Circle3Result): void {
    const n0dD = dot(circle0.normal, D);
    const normProj = mul(n0dD, circle0.normal);
    const compProj = sub(D, normProj);
    let U = compProj.clone();
    const d = normalize(U);

    // The configuration is determined by the relative location of the
    // intervals of projection of the circles on to the D-line. Circle0
    // projects to [-r0,r0] and circle1 projects to [d-r1,d+r1].
    const r0 = circle0.radius, r1 = circle1.radius;
    const dmr1 = d - r1;
    let distance: number;
    if (dmr1 >= r0) {  // d >= r0 + r1
        // The circles are separated (d > r0 + r1) or tangent with one outside
        // the other (d = r0 + r1).
        distance = dmr1 - r0;
        result.numClosestPairs = 1;
        result.circle0Closest[0] = add(circle0.center, mul(r0, U));
        result.circle1Closest[0] = sub(circle1.center, mul(r1, U));
        result.equidistant = false;
    } else {  // d < r0 + r1
        // The cases implicitly use the knowledge that d >= 0.
        const dpr1 = d + r1;
        if (dpr1 <= r0) {
            // Circle1 is inside circle0.
            distance = r0 - dpr1;
            result.numClosestPairs = 1;
            if (d > 0) {
                result.circle0Closest[0] = add(circle0.center, mul(r0, U));
                result.circle1Closest[0] = add(circle1.center, mul(r1, U));
                result.equidistant = false;
            } else {
                // The circles are concentric, so U = (0,0,0). Construct a
                // vector perpendicular to N0 to use for closest points.
                U = getOrthogonal(circle0.normal, true);
                result.circle0Closest[0] = add(circle0.center, mul(r0, U));
                result.circle1Closest[0] = add(circle1.center, mul(r1, U));
                result.equidistant = true;
            }
        } else if (dmr1 <= -r0) {
            // Circle0 is inside circle1.
            distance = -r0 - dmr1;
            result.numClosestPairs = 1;
            if (d > 0) {
                result.circle0Closest[0] = sub(circle0.center, mul(r0, U));
                result.circle1Closest[0] = sub(circle1.center, mul(r1, U));
                result.equidistant = false;
            } else {
                // The circles are concentric, so U = (0,0,0). Construct a
                // vector perpendicular to N0 to use for closest points.
                U = getOrthogonal(circle0.normal, true);
                result.circle0Closest[0] = add(circle0.center, mul(r0, U));
                result.circle1Closest[0] = add(circle1.center, mul(r1, U));
                result.equidistant = true;
            }
        } else {
            // The circles are overlapping. The two points of intersection are
            // C0 + s*(C1-C0) +/- h*Cross(N,U), where
            // s = (1 + (r0^2 - r1^2)/d^2)/2 and h = sqrt(r0^2 - s^2 * d^2).
            const r0sqr = r0 * r0, r1sqr = r1 * r1, dsqr = d * d;
            const s = (1 + (r0sqr - r1sqr) / dsqr) / 2;
            const arg = Math.max(r0sqr - dsqr * s * s, 0);
            const h = Math.sqrt(arg);
            const midpoint = add(circle0.center, mul(s, compProj));
            const hNxU = mul(h, cross(circle0.normal, U));
            distance = 0;
            result.numClosestPairs = 2;
            result.circle0Closest[0] = add(midpoint, hNxU);
            result.circle0Closest[1] = sub(midpoint, hNxU);
            result.circle1Closest[0] = add(result.circle0Closest[0], normProj);
            result.circle1Closest[1] = add(result.circle0Closest[1], normProj);
            result.equidistant = false;
        }
    }

    result.sqrDistance = distance * distance + n0dD * n0dD;
    result.distance = Math.sqrt(result.sqrDistance);
}

export class DistCircle3Circle3
    implements DCPQuery<Circle3, Circle3, DistCircle3Circle3Result> {
    compute(inCircle0: Circle3, inCircle1: Circle3): DistCircle3Circle3Result {
        const result: DistCircle3Circle3Result = {
            distance: 0,
            sqrDistance: 0,
            numClosestPairs: 0,
            circle0Closest: [new Vector(3), new Vector(3)],
            circle1Closest: [new Vector(3), new Vector(3)],
            equidistant: false
        };

        // Transform the circles by a translation, rotation, and uniform
        // scaling so that circle1.center = (0,0,0), circle1.normal = (0,0,1),
        // circle1.radius = 1, circle0.center[0] = 0, and circle0.radius =
        // min(inCircle0.radius, inCircle1.radius) /
        // max(inCircle0.radius, inCircle1.radius). The transformation is
        // Q = scale*rotate*(P+translate). The inverse transformation is
        // P = (1/scale)*Transpose(rotate)*Q - translate.
        const { circle0, circle1, rotate, translate, scale, swapped } =
            prepareCircles(inCircle0, inCircle1);

        if (circle0.normal.values[2] < 1) {
            // Convert the circle members to rationals. This is required to
            // avoid significant floating-point rounding errors when creating
            // the polynomials. If this is not done, the polynomial root
            // finder produces inaccurate results.
            const rZero = r(0), rOne = r(1), rTwo = r(2);
            const rC0 = [rZero, r(circle0.center.values[1]),
                r(circle0.center.values[2])];
            const rN0 = [r(circle0.normal.values[0]),
                r(circle0.normal.values[1]), r(circle0.normal.values[2])];
            const rR0 = r(circle0.radius);
            // D = -C0, U1 = (1,0,0), V1 = (0,1,0)

            // Construct the polynomial phi(cos(theta)).
            const rR0sqr = rR0.mul(rR0);
            const rN0xD = rCross(rN0, rC0.map(c => c.negated()));
            const rA0 = rC0[0].negated();  // r1 * Dot(D,U1)
            const rA1 = rC0[1].negated();  // r1 * Dot(D,V1)
            const rA2 = rDot(rN0xD, rN0xD);
            // r1 * Dot(N0xD, N0xU1), where N0xU1 = (0, N0[2], -N0[1])
            const rA3 = rN0xD[1].mul(rN0[2]).sub(rN0xD[2].mul(rN0[1]));
            // r1 * Dot(N0xD, N0xV1), where N0xV1 = (-N0[2], 0, N0[0])
            const rA4 = rN0xD[2].mul(rN0[0]).sub(rN0xD[0].mul(rN0[2]));
            // r1^2 * Dot(N0xU1, N0xU1)
            const rA5 = rN0[1].mul(rN0[1]).add(rN0[2].mul(rN0[2]));
            // r1^2 * Dot(N0xU1, N0xV1)
            const rA6 = rN0[0].mul(rN0[1]).negated();
            // r1^2 * Dot(N0xV1, N0xV1)
            const rA7 = rN0[0].mul(rN0[0]).add(rN0[2].mul(rN0[2]));

            const rP0: RPoly = [rA2.add(rA7), rTwo.mul(rA3), rA5.sub(rA7)];
            const rP1: RPoly = [rTwo.mul(rA4), rTwo.mul(rA6)];
            const rP2: RPoly = [rZero, rA1];
            const rP3: RPoly = [rA0.negated()];
            const rP4: RPoly = [rA6.negated(), rA4, rTwo.mul(rA6)];
            const rP5: RPoly = [rA3.negated(), rA7.sub(rA5)];
            const rTmp0: RPoly = [rOne, rZero, rOne.negated()];  // 1 - c^2
            const rTmp1 = rpAdd(rpMul(rP2, rP2),
                rpMul(rpMul(rTmp0, rP3), rP3));
            const rTmp2 = rpMul(rpScale(rTwo, rP2), rP3);
            const rTmp3 = rpAdd(rpMul(rP4, rP4),
                rpMul(rpMul(rTmp0, rP5), rP5));
            const rTmp4 = rpMul(rpScale(rTwo, rP4), rP5);
            const rP6 = rpSub(rpAdd(rpMul(rP0, rTmp1),
                rpMul(rpMul(rTmp0, rP1), rTmp2)), rpScale(rR0sqr, rTmp3));
            const rP7 = rpSub(rpAdd(rpMul(rP0, rTmp2), rpMul(rP1, rTmp1)),
                rpScale(rR0sqr, rTmp4));

            // Parameters for polynomial root finding. Only the unique roots
            // are needed. The pairs[] array stores the (cosine,sine)
            // information mentioned in the PDF.
            const pairs: Array<[number, number]> = [];

            if (!rpIsZero(rP7)) {
                // H(cs,sn) = p6(cs) + sn * p7(cs)
                const rPhi = rpSub(rpMul(rP6, rP6),
                    rpMul(rpMul(rTmp0, rP7), rP7));
                logAssert(rPhi.length - 1 > 0, 'Unexpected degree for phi.');

                for (const rCos of solveUnique(rPhi)) {
                    if (!BSRational.fabs(rCos).lessThanOrEqual(rOne)) {
                        continue;
                    }
                    const rValue = rpEval(rP7, rCos);
                    if (rValue.getSign() !== 0) {
                        // Because phi(cs) = 0, the sine below satisfies
                        // sn^2 = 1 - cs^2, so (cs,sn) is a point of the unit
                        // circle as required.
                        const rSin = rpEval(rP6, rCos).negated().div(rValue);
                        pairs.push([rCos.toNumber(), rSin.toNumber()]);
                    } else {
                        const rSin = BSRational.sqrt(rOne.sub(rCos.mul(rCos)));
                        pairs.push([rCos.toNumber(), rSin.toNumber()]);
                        if (rSin.getSign() !== 0) {
                            pairs.push([rCos.toNumber(), -rSin.toNumber()]);
                        }
                    }
                }

                // Upstream uses only the roots of phi. The roots of p6 are
                // added here for robustness. Squaring H to eliminate the sine
                // gives phi = p6^2 - (1-cs^2)*p7^2, which has a double root
                // wherever the configuration is mirror symmetric. Rounding
                // errors in the transformed circles split such a double root
                // into two simple roots that can be closer together than the
                // resolution of the floating-point bisection, so the root
                // finder reports neither and upstream then reports a critical
                // point that is not the closest pair. In that configuration p7
                // is numerically negligible and H reduces to p6, whose roots
                // pair with sn = +/-sqrt(1-cs^2). Adding those candidates is
                // safe: each is a point of circle1 whose distance to circle0
                // is computed exactly below, so an extra candidate can never
                // produce a distance smaller than the true minimum. See the
                // PR body.
                for (const rCos of solveUnique(rP6)) {
                    if (!BSRational.fabs(rCos).lessThanOrEqual(rOne)) {
                        continue;
                    }
                    const rSin = BSRational.sqrt(rOne.sub(rCos.mul(rCos)));
                    pairs.push([rCos.toNumber(), rSin.toNumber()]);
                    if (rSin.getSign() !== 0) {
                        pairs.push([rCos.toNumber(), -rSin.toNumber()]);
                    }
                }
            } else {
                // H(cs,sn) = p6(cs)
                logAssert(rP6.length - 1 > 0, 'Unexpected degree for p6.');

                for (const rCos of solveUnique(rP6)) {
                    if (!BSRational.fabs(rCos).lessThanOrEqual(rOne)) {
                        continue;
                    }
                    const rSin = BSRational.sqrt(rOne.sub(rCos.mul(rCos)));
                    pairs.push([rCos.toNumber(), rSin.toNumber()]);
                    if (rSin.getSign() !== 0) {
                        pairs.push([rCos.toNumber(), -rSin.toNumber()]);
                    }
                }
            }

            // Upstream indexes a default-constructed std::array of 16
            // ClosestInfo objects. When no (cosine,sine) pair survives the
            // filtering, candidates[0] is that default object and the query
            // silently reports distance 0 at the origin. The port traps the
            // configuration instead; see the PR body.
            logAssert(pairs.length > 0,
                'DistCircle3Circle3: no closest-point candidates.');

            // Convert the rational values to floating-point values for fast
            // computation of the closest-point candidates.
            const candidates: ClosestInfo[] = [];
            for (let i = 0; i < pairs.length; ++i) {
                let delta = add(sub(circle1.center, circle0.center),
                    mul(circle1.radius,
                        Vector.fromArray([pairs[i][0], pairs[i][1], 0])));
                const circle1Closest = add(circle0.center, delta);

                const n0xDelta = cross(circle0.normal, delta);
                const lenN0xDelta = length(n0xDelta);
                let sqrDistance: number;
                let circle0Closest: Vector;
                let equidistant: boolean;
                if (lenN0xDelta > 0) {
                    const n0dDelta = dot(circle0.normal, delta);
                    const diff = lenN0xDelta - circle0.radius;
                    sqrDistance = n0dDelta * n0dDelta + diff * diff;
                    delta = sub(delta, mul(n0dDelta, circle0.normal));
                    normalize(delta);
                    circle0Closest = add(circle0.center,
                        mul(circle0.radius, delta));
                    equidistant = false;
                } else {
                    // Delta is parallel to the normal of circle0, so every
                    // point of circle0 is equidistant from circle1Closest.
                    let U0: Vector;
                    if (Math.abs(circle0.normal.values[0]) >
                        Math.abs(circle0.normal.values[1])) {
                        U0 = Vector.fromArray([-circle0.normal.values[2], 0,
                            circle0.normal.values[0]]);
                    } else {
                        U0 = Vector.fromArray([0, circle0.normal.values[2],
                            -circle0.normal.values[1]]);
                    }
                    normalize(U0);

                    const r0U0 = mul(circle0.radius, U0);
                    const diff = sub(delta, r0U0);
                    sqrDistance = dot(diff, diff);
                    circle0Closest = add(circle0.center, r0U0);
                    equidistant = true;
                }

                candidates.push({
                    sqrDistance, circle0Closest, circle1Closest, equidistant
                });
            }

            candidates.sort((a, b) => a.sqrDistance - b.sqrDistance);

            result.numClosestPairs = 1;
            result.sqrDistance = candidates[0].sqrDistance;
            result.distance = Math.sqrt(result.sqrDistance);
            result.circle0Closest[0] = candidates[0].circle0Closest;
            result.circle1Closest[0] = candidates[0].circle1Closest;
            result.equidistant = candidates[0].equidistant;
            // Upstream tests rRoots.size() > 1, the number of roots before
            // the uniqueness and |cos| <= 1 filtering, rather than the number
            // of candidates. When a single candidate survives, candidates[1]
            // is a default-constructed object with sqrDistance 0 and closest
            // points at the origin, so a zero-distance query reports a bogus
            // second pair. The port tests the candidate count; see the PR
            // body.
            if (candidates.length > 1 &&
                candidates[1].sqrDistance === candidates[0].sqrDistance) {
                result.numClosestPairs = 2;
                result.circle0Closest[1] = candidates[1].circle0Closest;
                result.circle1Closest[1] = candidates[1].circle1Closest;
            }
        } else {
            // The planes of the circles are parallel. Whether the planes are
            // the same or different, the problem reduces to determining how
            // two circles in the same plane are separated, tangent with one
            // circle outside the other, overlapping, or one circle contained
            // inside the other circle.
            const D = sub(circle1.center, circle0.center);
            doQueryParallelPlanes(circle0, circle1, D, result);
        }

        result.distance /= scale;
        result.sqrDistance = result.distance * result.distance;
        for (let i = 0; i < result.numClosestPairs; ++i) {
            result.circle0Closest[i] = sub(
                div(mulMatrix(result.circle0Closest[i], rotate) as Vector,
                    scale), translate);
            result.circle1Closest[i] = sub(
                div(mulMatrix(result.circle1Closest[i], rotate) as Vector,
                    scale), translate);
        }

        if (swapped) {
            const temp = result.circle0Closest;
            result.circle0Closest = result.circle1Closest;
            result.circle1Closest = temp;
        }
        return result;
    }
}
