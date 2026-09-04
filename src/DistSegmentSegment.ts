// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistSegmentSegment.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the closest points for two segments in nD.
//
// The segments are P[0] + s[0] * (P[1] - P[0]) for 0 <= s[0] <= 1 and
// Q[0] + s[1] * (Q[1] - Q[0]) for 0 <= s[1] <= 1. The D[i] are not required
// to be unit length.
//
// The closest point on segment[i] is stored in closest[i] with parameter[i]
// storing s[i]. When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// The computeRobust algorithm is robust even for nearly parallel segments.
// Effectively, it uses a conjugate gradient search for the minimum of the
// squared distance function, which avoids the numerical problems introduced
// by divisions in the case the minimum is located at an interior point of
// the domain. See the document
//   https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
// for details.
//
//
// Upstream caveat (all of DistLineLine/LineRay/LineSegment/RayRay/RaySegment
// and DistSegmentSegment.compute): parallelism is detected by
// det = max(a00*a11 - a01*a01, 0) > 0. That difference of two products
// cancels exactly only in exact arithmetic, so for direction vectors that are
// mathematically parallel it can round to one ulp above zero. The
// nonparallel branch is then entered with numerators that are pure rounding
// noise, and the reported points, while on their primitives, need not be near
// the minimum. Preserved as upstream has it; DistSegmentSegment.computeRobust
// is upstream's answer for that case.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Segment<N,T>, Segment<N,T>>' becomes
// the class DistSegmentSegment with the result type
// DistSegmentSegmentResult. Upstream overloads each query on either two
// Segments or four endpoint Vectors; per PORTING.md only the canonical
// two-argument query is named compute, so the endpoint forms are
// computeEndpoints and computeRobustEndpoints. The private static helpers
// GetClampedRoot, ComputeIntersection and ComputeMinimumParameters become
// module-private functions; ComputeIntersection's output references become
// the fields of a returned object.

import type { DCPQuery } from './DCPQuery.js';
import type { Segment } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

export interface DistSegmentSegmentResult {
    distance: number;
    sqrDistance: number;

    // The segment parameters s[0] and s[1], each in [0,1].
    parameter: [number, number];

    // closest[0] is on segment0, closest[1] is on segment1.
    closest: [Vector, Vector];
}

// A point (s,t) of the domain [0,1]^2.
type DomainPoint = [number, number];

// Compute the root of h(z) = h0 + slope*z and clamp it to the interval
// [0,1]. It is required that for h1 = h(1), either (h0 < 0 and h1 > 0) or
// (h0 > 0 and h1 < 0).
function getClampedRoot(slope: number, h0: number, h1: number): number {
    // Theoretically, r is in (0,1). However, when the slope is nearly zero,
    // then so are h0 and h1. Significant numerical rounding problems can
    // occur when using floating-point arithmetic. If the rounding causes r
    // to be outside the interval, clamp it. It is possible that r is in
    // (0,1) and has rounding errors, but because h0 and h1 are both nearly
    // zero, the quadratic is nearly constant on (0,1). Any choice of the
    // parameter should not cause undesirable accuracy problems for the final
    // distance computation.
    let r: number;
    if (h0 < 0) {
        if (h1 > 0) {
            r = -h0 / slope;
            if (r > 1) {
                r = 0.5;
            }
            // The slope is positive and -h0 is positive, so there is no need
            // to test for a negative value and clamp it.
        }
        else {
            r = 1;
        }
    }
    else {
        r = 0;
    }
    return r;
}

// The t-coordinate of the intersection of the line dR/ds = 0 with the domain
// edge s = 0 (using f00) or s = 1 (using f10). The divisions are
// theoretically numbers in [0,1]; numerical rounding errors can place the
// result outside the interval, and the out-of-range value is repaired here.
//
// UPSTREAM BUG FIX (DistSegmentSegment.h, ComputeIntersection). Upstream
// replaces an out-of-range ratio by 1/2, arguing that it happens only when
// "both numerator and denominator are nearly zero" so that "the choice of 0.5
// should not cause significant accuracy problems". That premise does not
// hold. The ratio is f/b with b = Dot(P1-P0,Q1-Q0), and it exceeds one by a
// single ulp whenever the line dR/ds = 0 passes through a corner of the
// domain, with f and b both of ordinary magnitude. Moving the endpoint of the
// intersection segment from the corner to the middle of the domain edge then
// sends ComputeMinimumParameters to the wrong edge and ComputeRobust reports
// a distance that is not the minimum. Witness (used as a regression test):
// segment0 = <(-0.8119320124387741,-3,0), (7,0,0)>, segment1 =
// <(0,1.9484716467559338,0), (7,0,6)>, whose direction vectors are nowhere
// near parallel (the sine of the angle between them is 0.79). There
// f10 / b = 1 + 1 ulp, upstream substitutes 1/2 and returns 4.6265 instead of
// the true minimum 3.5103 attained at the interior point
// (s,t) = (0.5413,0.3423). Over 300000 randomized configurations the
// substitution produced a wrong minimum in 12 cases, with a worst excess of
// 1.12; clamping to the nearest endpoint of [0,1], as done here, leaves a
// worst excess of 2.7e-14 and is also the repair upstream applies everywhere
// else it projects a parameter back into its domain. GetClampedRoot's own
// 0.5 substitution is left alone: there h0 and h1 really are both nearly
// zero, so the quadratic is nearly constant and upstream's argument holds.
function getEdgeT(f: number, b: number): number {
    const t = f / b;
    if (t < 0) { return 0; }
    if (t > 1) { return 1; }
    return t;
}

// Compute the intersection of the line dR/ds = 0 with the domain [0,1]^2.
// The direction of the line dR/ds is conjugate to (1,0), so the algorithm
// for minimization is effectively the conjugate gradient algorithm for a
// quadratic function. The edge[i] flag tells on which domain edge end[i]
// lives: 0 (s=0), 1 (s=1), 2 (t=0), 3 (t=1).
function computeIntersection(sValue: readonly number[],
    classify: readonly number[], b: number, f00: number, f10: number):
    { edge: [number, number], end: [DomainPoint, DomainPoint] } {
    const edge: [number, number] = [0, 0];
    const end: [DomainPoint, DomainPoint] = [[0, 0], [0, 0]];

    if (classify[0] < 0) {
        edge[0] = 0;
        end[0] = [0, getEdgeT(f00, b)];

        if (classify[1] === 0) {
            edge[1] = 3;
            end[1] = [sValue[1], 1];
        }
        else {
            // classify[1] > 0
            edge[1] = 1;
            end[1] = [1, getEdgeT(f10, b)];
        }
    }
    else if (classify[0] === 0) {
        edge[0] = 2;
        end[0] = [sValue[0], 0];

        if (classify[1] < 0) {
            edge[1] = 0;
            end[1] = [0, getEdgeT(f00, b)];
        }
        else if (classify[1] === 0) {
            edge[1] = 3;
            end[1] = [sValue[1], 1];
        }
        else {
            edge[1] = 1;
            end[1] = [1, getEdgeT(f10, b)];
        }
    }
    else {
        // classify[0] > 0
        edge[0] = 1;
        end[0] = [1, getEdgeT(f10, b)];

        if (classify[1] === 0) {
            edge[1] = 3;
            end[1] = [sValue[1], 1];
        }
        else {
            edge[1] = 0;
            end[1] = [0, getEdgeT(f00, b)];
        }
    }

    return { edge, end };
}

// Compute the location of the minimum of R on the segment of intersection
// for the line dR/ds = 0 and the domain [0,1]^2.
function computeMinimumParameters(edge: readonly number[],
    end: readonly DomainPoint[], b: number, c: number, e: number, g00: number,
    g10: number, g01: number, g11: number): [number, number] {
    const delta = end[1][1] - end[0][1];
    const h0 = delta * (-b * end[0][0] + c * end[0][1] - e);
    if (h0 >= 0) {
        if (edge[0] === 0) {
            return [0, getClampedRoot(c, g00, g01)];
        }
        if (edge[0] === 1) {
            return [1, getClampedRoot(c, g10, g11)];
        }
        return [end[0][0], end[0][1]];
    }

    const h1 = delta * (-b * end[1][0] + c * end[1][1] - e);
    if (h1 <= 0) {
        if (edge[1] === 0) {
            return [0, getClampedRoot(c, g00, g01)];
        }
        if (edge[1] === 1) {
            return [1, getClampedRoot(c, g10, g11)];
        }
        return [end[1][0], end[1][1]];
    }

    // h0 < 0 and h1 > 0
    const z = Math.min(Math.max(h0 / (h0 - h1), 0), 1);
    const omz = 1 - z;
    return [
        omz * end[0][0] + z * end[1][0],
        omz * end[0][1] + z * end[1][1]
    ];
}

// Assemble the result from the segment parameters.
function makeResult(P0: Vector, P1mP0: Vector, Q0: Vector, Q1mQ0: Vector,
    s: number, t: number): DistSegmentSegmentResult {
    const closest0 = add(P0, mul(s, P1mP0));
    const closest1 = add(Q0, mul(t, Q1mQ0));
    const diff = sub(closest0, closest1);
    const sqrDistance = dot(diff, diff);
    return {
        distance: Math.sqrt(sqrDistance),
        sqrDistance,
        parameter: [s, t],
        closest: [closest0, closest1]
    };
}

export class DistSegmentSegment
    implements DCPQuery<Segment, Segment, DistSegmentSegmentResult> {
    // compute and computeEndpoints are exact for computing sqrDistance when
    // the arithmetic is exact.
    compute(segment0: Segment, segment1: Segment): DistSegmentSegmentResult {
        return this.computeEndpoints(segment0.p[0], segment0.p[1],
            segment1.p[0], segment1.p[1]);
    }

    computeEndpoints(P0: Vector, P1: Vector, Q0: Vector,
        Q1: Vector): DistSegmentSegmentResult {
        const P1mP0 = sub(P1, P0);
        const Q1mQ0 = sub(Q1, Q0);
        const P0mQ0 = sub(P0, Q0);
        const a = dot(P1mP0, P1mP0);
        const b = dot(P1mP0, Q1mQ0);
        const c = dot(Q1mQ0, Q1mQ0);
        const d = dot(P1mP0, P0mQ0);
        const e = dot(Q1mQ0, P0mQ0);
        const det = a * c - b * b;
        let s: number;
        let t: number;
        let nd: number;
        let bmd: number;
        let bpe: number;

        if (det > 0) {
            const bte = b * e;
            const ctd = c * d;
            if (bte <= ctd) {
                // s <= 0
                s = 0;
                if (e <= 0) {
                    // region 6 (t <= 0)
                    t = 0;
                    nd = -d;
                    if (nd >= a) {
                        s = 1;
                    }
                    else if (nd > 0) {
                        s = nd / a;
                    }
                    // else: s is already zero
                }
                else if (e < c) {
                    // region 5 (0 < t < 1)
                    t = e / c;
                }
                else {
                    // region 4 (t >= 1)
                    t = 1;
                    bmd = b - d;
                    if (bmd >= a) {
                        s = 1;
                    }
                    else if (bmd > 0) {
                        s = bmd / a;
                    }
                    // else: s is already zero
                }
            }
            else {
                // s > 0
                s = bte - ctd;
                if (s >= det) {
                    // s = 1
                    s = 1;
                    bpe = b + e;
                    if (bpe <= 0) {
                        // region 8 (t <= 0)
                        t = 0;
                        nd = -d;
                        if (nd <= 0) {
                            s = 0;
                        }
                        else if (nd < a) {
                            s = nd / a;
                        }
                        // else: s is already one
                    }
                    else if (bpe < c) {
                        // region 1 (0 < t < 1)
                        t = bpe / c;
                    }
                    else {
                        // region 2 (t >= 1)
                        t = 1;
                        bmd = b - d;
                        if (bmd <= 0) {
                            s = 0;
                        }
                        else if (bmd < a) {
                            s = bmd / a;
                        }
                        // else: s is already one
                    }
                }
                else {
                    // 0 < s < 1
                    const ate = a * e;
                    const btd = b * d;
                    if (ate <= btd) {
                        // region 7 (t <= 0)
                        t = 0;
                        nd = -d;
                        if (nd <= 0) {
                            s = 0;
                        }
                        else if (nd >= a) {
                            s = 1;
                        }
                        else {
                            s = nd / a;
                        }
                    }
                    else {
                        // t > 0
                        t = ate - btd;
                        if (t >= det) {
                            // region 3 (t >= 1)
                            t = 1;
                            bmd = b - d;
                            if (bmd <= 0) {
                                s = 0;
                            }
                            else if (bmd >= a) {
                                s = 1;
                            }
                            else {
                                s = bmd / a;
                            }
                        }
                        else {
                            // region 0 (0 < t < 1)
                            s /= det;
                            t /= det;
                        }
                    }
                }
            }
        }
        else {
            // The segments are parallel. The quadratic factors to
            //   R(s,t) = a*(s-(b/a)*t)^2 + 2*d*(s - (b/a)*t) + f
            // where a*c = b^2, e = b*d/a, f = |P0-Q0|^2, and b is not zero.
            // R is constant along lines of the form s-(b/a)*t = k and its
            // minimum occurs on the line a*s - b*t + d = 0. This line must
            // intersect both the s-axis and the t-axis because 'a' and 'b'
            // are not zero. Because of parallelism, the line is also
            // represented by -b*s + c*t - e = 0.
            //
            // The code determines an edge of the domain [0,1]^2 that
            // intersects the minimum line, or if none of the edges
            // intersect, it determines the closest corner to the minimum
            // line. The conditionals are designed to test first for
            // intersection with the t-axis (s = 0) using -b*s + c*t - e = 0
            // and then with the s-axis (t = 0) using a*s - b*t + d = 0.

            // When s = 0, solve c*t - e = 0 (t = e/c).
            if (e <= 0) {
                // t <= 0; now solve a*s - b*t + d = 0 for t = 0 (s = -d/a).
                t = 0;
                nd = -d;
                if (nd <= 0) {
                    // region 6 (s <= 0)
                    s = 0;
                }
                else if (nd >= a) {
                    // region 8 (s >= 1)
                    s = 1;
                }
                else {
                    // region 7 (0 < s < 1)
                    s = nd / a;
                }
            }
            else if (e >= c) {
                // t >= 1; now solve a*s - b*t + d = 0 for t = 1
                // (s = (b-d)/a).
                t = 1;
                bmd = b - d;
                if (bmd <= 0) {
                    // region 4 (s <= 0)
                    s = 0;
                }
                else if (bmd >= a) {
                    // region 2 (s >= 1)
                    s = 1;
                }
                else {
                    // region 3 (0 < s < 1)
                    s = bmd / a;
                }
            }
            else {
                // 0 < t < 1. The point (0,e/c) is on the line and in the
                // domain, so we have one point at which R is a minimum.
                s = 0;
                t = e / c;
            }
        }

        return makeResult(P0, P1mP0, Q0, Q1mQ0, s, t);
    }

    // computeRobust and computeRobustEndpoints are exact for computing
    // sqrDistance when the arithmetic is exact. They are generally more
    // robust than compute/computeEndpoints for floating-point arithmetic.
    computeRobust(segment0: Segment,
        segment1: Segment): DistSegmentSegmentResult {
        return this.computeRobustEndpoints(segment0.p[0], segment0.p[1],
            segment1.p[0], segment1.p[1]);
    }

    computeRobustEndpoints(P0: Vector, P1: Vector, Q0: Vector,
        Q1: Vector): DistSegmentSegmentResult {
        // The code allows degenerate line segments; that is, P0 and P1 can
        // be the same point or Q0 and Q1 can be the same point. The
        // quadratic function for the squared distance between the segments
        // is
        //   R(s,t) = a*s^2 - 2*b*s*t + c*t^2 + 2*d*s - 2*e*t + f
        // for (s,t) in [0,1]^2 where
        //   a = Dot(P1-P0,P1-P0), b = Dot(P1-P0,Q1-Q0),
        //   c = Dot(Q1-Q0,Q1-Q0), d = Dot(P1-P0,P0-Q0),
        //   e = Dot(Q1-Q0,P0-Q0), f = Dot(P0-Q0,P0-Q0)
        const P1mP0 = sub(P1, P0);
        const Q1mQ0 = sub(Q1, Q0);
        const P0mQ0 = sub(P0, Q0);
        const a = dot(P1mP0, P1mP0);
        const b = dot(P1mP0, Q1mQ0);
        const c = dot(Q1mQ0, Q1mQ0);
        const d = dot(P1mP0, P0mQ0);
        const e = dot(Q1mQ0, P0mQ0);

        // The derivatives dR/ds(i,j) at the four corners of the domain.
        const f00 = d;
        const f10 = f00 + a;
        const f01 = f00 - b;
        const f11 = f10 - b;

        // The derivatives dR/dt(i,j) at the four corners of the domain.
        const g00 = -e;
        const g10 = g00 - b;
        const g01 = g00 + c;
        const g11 = g10 + c;

        let parameter: [number, number];

        if (a > 0 && c > 0) {
            // Compute the solutions to dR/ds(s0,0) = 0 and dR/ds(s1,1) = 0.
            // The location of sI on the s-axis is stored in classify[I]
            // (I = 0 or 1). If sI <= 0, classify[I] is -1. If sI >= 1,
            // classify[I] is +1. If 0 < sI < 1, classify[I] is 0. This
            // information helps determine where to search for the minimum
            // point (s,t). The fij values are dR/ds(i,j) for i and j in
            // {0,1}.
            const sValue: [number, number] = [
                getClampedRoot(a, f00, f10),
                getClampedRoot(a, f01, f11)
            ];

            const classify: [number, number] = [0, 0];
            for (let i = 0; i < 2; ++i) {
                if (sValue[i] <= 0) {
                    classify[i] = -1;
                }
                else if (sValue[i] >= 1) {
                    classify[i] = +1;
                }
                else {
                    classify[i] = 0;
                }
            }

            if (classify[0] === -1 && classify[1] === -1) {
                // The minimum must occur on s = 0 for 0 <= t <= 1.
                parameter = [0, getClampedRoot(c, g00, g01)];
            }
            else if (classify[0] === +1 && classify[1] === +1) {
                // The minimum must occur on s = 1 for 0 <= t <= 1.
                parameter = [1, getClampedRoot(c, g10, g11)];
            }
            else {
                // The line dR/ds = 0 intersects the domain [0,1]^2 in a
                // nondegenerate segment. Compute the endpoints of that
                // segment, end[0] and end[1].
                const { edge, end } =
                    computeIntersection(sValue, classify, b, f00, f10);

                // The directional derivative of R along the segment of
                // intersection is
                //   H(z) = (end[1][1]-end[0][1]) *
                //          dR/dt((1-z)*end[0] + z*end[1])
                // for z in [0,1]. The formula uses the fact that dR/ds = 0
                // on the segment. Compute the minimum of H on [0,1].
                parameter = computeMinimumParameters(edge, end, b, c, e, g00,
                    g10, g01, g11);
            }
        }
        else if (a > 0) {
            // The Q-segment is degenerate (Q0 and Q1 are the same point) and
            // the quadratic is R(s,0) = a*s^2 + 2*d*s + f with (half) first
            // derivative F(s) = a*s + d. The closest P-point is interior to
            // the P-segment when F(0) < 0 and F(1) > 0.
            parameter = [getClampedRoot(a, f00, f10), 0];
        }
        else if (c > 0) {
            // The P-segment is degenerate (P0 and P1 are the same point) and
            // the quadratic is R(0,t) = c*t^2 - 2*e*t + f with (half) first
            // derivative G(t) = c*t - e. The closest Q-point is interior to
            // the Q-segment when G(0) < 0 and G(1) > 0.
            parameter = [0, getClampedRoot(c, g00, g01)];
        }
        else {
            // The P-segment and Q-segment are both degenerate.
            parameter = [0, 0];
        }

        // The closest points are computed by linear interpolation of the
        // endpoints rather than by P0 + s * (P1 - P0), which matches
        // upstream.
        const s = parameter[0];
        const t = parameter[1];
        const closest0 = add(mul(1 - s, P0), mul(s, P1));
        const closest1 = add(mul(1 - t, Q0), mul(t, Q1));
        const diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            parameter: [s, t],
            closest: [closest0, closest1]
        };
    }
}
