// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointTriangle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a solid triangle in nD.
//
// The triangle has vertices <V[0],V[1],V[2]>. A triangle point is
// X = sum_{i=0}^2 b[i] * V[i], where 0 <= b[i] <= 1 for all i and
// sum_{i=0}^2 b[i] = 1.
//
// The input point is stored in closest[0]. The closest point on the triangle
// is stored in closest[1] with barycentric coordinates (b[0],b[1],b[2]).
//
// For a description of the algebraic details of the quadratic minimization
// approach used by compute(), see
//   https://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
// Although the document describes the 3D case, the construction applies in
// general dimensions N. The useConjugateGradient function uses conjugate
// gradient minimization.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Triangle<N,T>>' becomes
// the class DistPointTriangle with the result type DistPointTriangleResult.
// The private helpers GetMinEdge02, GetMinEdge12 and GetMinInterior become
// module-private functions. The second query 'UseConjugateGradient' keeps a
// descriptive name, useConjugateGradient, because only the canonical query
// may be named compute.

import type { DCPQuery } from './DCPQuery';
import type { Triangle } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistPointTriangleResult {
    distance: number;
    sqrDistance: number;

    // The barycentric coordinates (b[0],b[1],b[2]) of closest[1] relative to
    // the triangle vertices; they are nonnegative and sum to 1.
    barycentric: [number, number, number];

    // closest[0] is the input point, closest[1] is the closest triangle
    // point.
    closest: [Vector, Vector];
}

// A point in the (s,t) parameter domain of the triangle; the triangle point
// is V[0] + s * (V[1] - V[0]) + t * (V[2] - V[0]).
type Param = { s: number, t: number };

// The port of GetMinEdge02: minimize on the edge <V[0],V[2]>, that is, on
// s = 0.
function getMinEdge02(a11: number, b1: number): Param {
    if (b1 >= 0) {
        return { s: 0, t: 0 };
    }
    if (a11 + b1 <= 0) {
        return { s: 0, t: 1 };
    }
    return { s: 0, t: -b1 / a11 };
}

// The port of GetMinEdge12: minimize on the edge <V[1],V[2]>, that is, on
// s + t = 1.
function getMinEdge12(a01: number, a11: number, b1: number, f10: number,
    f01: number): Param {
    let t: number;
    const h0 = a01 + b1 - f10;
    if (h0 >= 0) {
        t = 0;
    }
    else {
        const h1 = a11 + b1 - f01;
        if (h1 <= 0) {
            t = 1;
        }
        else {
            t = h0 / (h0 - h1);
        }
    }
    return { s: 1 - t, t };
}

// The port of GetMinInterior: the minimum is at an interior point of the
// segment from p0 to p1, located by the linear interpolation of the
// directional derivatives h0 and h1.
function getMinInterior(p0: Param, h0: number, p1: Param,
    h1: number): Param {
    const z = h0 / (h0 - h1);
    const omz = 1 - z;
    return {
        s: omz * p0.s + z * p1.s,
        t: omz * p0.t + z * p1.t
    };
}

// Assemble the result common to both queries.
function makeResult(point: Vector, triangle: Triangle, edge0: Vector,
    edge1: Vector, s: number, t: number): DistPointTriangleResult {
    const closest0 = point.clone();
    const closest1 = add(triangle.v[0], add(mul(s, edge0), mul(t, edge1)));
    const diff = sub(closest0, closest1);
    const sqrDistance = dot(diff, diff);
    return {
        distance: Math.sqrt(sqrDistance),
        sqrDistance,
        barycentric: [1 - s - t, s, t],
        closest: [closest0, closest1]
    };
}

export class DistPointTriangle
    implements DCPQuery<Vector, Triangle, DistPointTriangleResult> {
    // This query is exact when using arbitrary-precision arithmetic. It can
    // be used also for floating-point arithmetic, but rounding errors can
    // sometimes lead to an inaccurate result. For floating-point, consider
    // useConjugateGradient(...) which is more robust.
    compute(point: Vector, triangle: Triangle): DistPointTriangleResult {
        const diff = sub(triangle.v[0], point);
        const edge0 = sub(triangle.v[1], triangle.v[0]);
        const edge1 = sub(triangle.v[2], triangle.v[0]);
        const a00 = dot(edge0, edge0);
        const a01 = dot(edge0, edge1);
        const a11 = dot(edge1, edge1);
        const b0 = dot(diff, edge0);
        const b1 = dot(diff, edge1);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s = a01 * b1 - a11 * b0;
        let t = a01 * b0 - a00 * b1;

        if (s + t <= det) {
            if (s < 0) {
                if (t < 0) {
                    // region 4
                    if (b0 < 0) {
                        t = 0;
                        if (-b0 >= a00) {
                            s = 1;
                        }
                        else {
                            s = -b0 / a00;
                        }
                    }
                    else {
                        s = 0;
                        if (b1 >= 0) {
                            t = 0;
                        }
                        else if (-b1 >= a11) {
                            t = 1;
                        }
                        else {
                            t = -b1 / a11;
                        }
                    }
                }
                else {
                    // region 3
                    s = 0;
                    if (b1 >= 0) {
                        t = 0;
                    }
                    else if (-b1 >= a11) {
                        t = 1;
                    }
                    else {
                        t = -b1 / a11;
                    }
                }
            }
            else if (t < 0) {
                // region 5
                t = 0;
                if (b0 >= 0) {
                    s = 0;
                }
                else if (-b0 >= a00) {
                    s = 1;
                }
                else {
                    s = -b0 / a00;
                }
            }
            else {
                // region 0: the minimum is at an interior point.
                s /= det;
                t /= det;
            }
        }
        else {
            let tmp0: number;
            let tmp1: number;
            let numer: number;
            let denom: number;

            if (s < 0) {
                // region 2
                tmp0 = a01 + b0;
                tmp1 = a11 + b1;
                if (tmp1 > tmp0) {
                    numer = tmp1 - tmp0;
                    denom = a00 - 2 * a01 + a11;
                    if (numer >= denom) {
                        s = 1;
                        t = 0;
                    }
                    else {
                        s = numer / denom;
                        t = 1 - s;
                    }
                }
                else {
                    s = 0;
                    if (tmp1 <= 0) {
                        t = 1;
                    }
                    else if (b1 >= 0) {
                        t = 0;
                    }
                    else {
                        t = -b1 / a11;
                    }
                }
            }
            else if (t < 0) {
                // region 6
                tmp0 = a01 + b1;
                tmp1 = a00 + b0;
                if (tmp1 > tmp0) {
                    numer = tmp1 - tmp0;
                    denom = a00 - 2 * a01 + a11;
                    if (numer >= denom) {
                        t = 1;
                        s = 0;
                    }
                    else {
                        t = numer / denom;
                        s = 1 - t;
                    }
                }
                else {
                    t = 0;
                    if (tmp1 <= 0) {
                        s = 1;
                    }
                    else if (b0 >= 0) {
                        s = 0;
                    }
                    else {
                        s = -b0 / a00;
                    }
                }
            }
            else {
                // region 1
                numer = a11 + b1 - a01 - b0;
                if (numer <= 0) {
                    s = 0;
                    t = 1;
                }
                else {
                    denom = a00 - 2 * a01 + a11;
                    if (numer >= denom) {
                        s = 1;
                        t = 0;
                    }
                    else {
                        s = numer / denom;
                        t = 1 - s;
                    }
                }
            }
        }

        return makeResult(point, triangle, edge0, edge1, s, t);
    }

    // The query is designed to be robust when using floating-point
    // arithmetic. For arbitrary-precision arithmetic, use compute(...).
    useConjugateGradient(point: Vector,
        triangle: Triangle): DistPointTriangleResult {
        const diff = sub(point, triangle.v[0]);
        const edge0 = sub(triangle.v[1], triangle.v[0]);
        const edge1 = sub(triangle.v[2], triangle.v[0]);
        const a00 = dot(edge0, edge0);
        const a01 = dot(edge0, edge1);
        const a11 = dot(edge1, edge1);
        const b0 = -dot(diff, edge0);
        const b1 = -dot(diff, edge1);

        const f00 = b0;
        const f10 = b0 + a00;
        const f01 = b0 + a01;

        let p0: Param;
        let p1: Param;
        let p: Param;
        let dt1: number;
        let h0: number;
        let h1: number;

        // Compute the endpoints p0 and p1 of the segment. The segment is
        // parameterized by L(z) = (1-z)*p0 + z*p1 for z in [0,1] and the
        // directional derivative of half the quadratic on the segment is
        // H(z) = Dot(p1-p0,gradient[Q](L(z))/2), where gradient[Q]/2 =
        // (F,G). By design, F(L(z)) = 0 for cases (2), (4), (5) and (6).
        // Cases (1) and (3) can correspond to no-intersection or
        // intersection of F = 0 with the triangle.
        if (f00 >= 0) {
            if (f01 >= 0) {
                // (1) p0 = (0,0), p1 = (0,1), H(z) = G(L(z))
                p = getMinEdge02(a11, b1);
            }
            else {
                // (2) p0 = (0,t10), p1 = (t01,1-t01),
                // H(z) = (t11 - t10)*G(L(z))
                p0 = { s: 0, t: f00 / (f00 - f01) };
                p1 = { s: f01 / (f01 - f10), t: 0 };
                p1.t = 1 - p1.s;
                dt1 = p1.t - p0.t;
                h0 = dt1 * (a11 * p0.t + b1);
                if (h0 >= 0) {
                    p = getMinEdge02(a11, b1);
                }
                else {
                    h1 = dt1 * (a01 * p1.s + a11 * p1.t + b1);
                    if (h1 <= 0) {
                        p = getMinEdge12(a01, a11, b1, f10, f01);
                    }
                    else {
                        p = getMinInterior(p0, h0, p1, h1);
                    }
                }
            }
        }
        else if (f01 <= 0) {
            if (f10 <= 0) {
                // (3) p0 = (1,0), p1 = (0,1), H(z) = G(L(z)) - F(L(z))
                p = getMinEdge12(a01, a11, b1, f10, f01);
            }
            else {
                // (4) p0 = (t00,0), p1 = (t01,1-t01), H(z) = t11*G(L(z))
                p0 = { s: f00 / (f00 - f10), t: 0 };
                p1 = { s: f01 / (f01 - f10), t: 0 };
                p1.t = 1 - p1.s;
                h0 = p1.t * (a01 * p0.s + b1);
                if (h0 >= 0) {
                    p = p0;  // getMinEdge01
                }
                else {
                    h1 = p1.t * (a01 * p1.s + a11 * p1.t + b1);
                    if (h1 <= 0) {
                        p = getMinEdge12(a01, a11, b1, f10, f01);
                    }
                    else {
                        p = getMinInterior(p0, h0, p1, h1);
                    }
                }
            }
        }
        else if (f10 <= 0) {
            // (5) p0 = (0,t10), p1 = (t01,1-t01),
            // H(z) = (t11 - t10)*G(L(z))
            p0 = { s: 0, t: f00 / (f00 - f01) };
            p1 = { s: f01 / (f01 - f10), t: 0 };
            p1.t = 1 - p1.s;
            dt1 = p1.t - p0.t;
            h0 = dt1 * (a11 * p0.t + b1);
            if (h0 >= 0) {
                p = getMinEdge02(a11, b1);
            }
            else {
                h1 = dt1 * (a01 * p1.s + a11 * p1.t + b1);
                if (h1 <= 0) {
                    p = getMinEdge12(a01, a11, b1, f10, f01);
                }
                else {
                    p = getMinInterior(p0, h0, p1, h1);
                }
            }
        }
        else {
            // (6) p0 = (t00,0), p1 = (0,t11), H(z) = t11*G(L(z))
            p0 = { s: f00 / (f00 - f10), t: 0 };
            p1 = { s: 0, t: f00 / (f00 - f01) };
            h0 = p1.t * (a01 * p0.s + b1);
            if (h0 >= 0) {
                p = p0;  // getMinEdge01
            }
            else {
                h1 = p1.t * (a11 * p1.t + b1);
                if (h1 <= 0) {
                    p = getMinEdge02(a11, b1);
                }
                else {
                    p = getMinInterior(p0, h0, p1, h1);
                }
            }
        }

        return makeResult(point, triangle, edge0, edge1, p.s, p.t);
    }
}
