import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrRay3Ellipsoid3TI,
    IntrRay3Ellipsoid3FI,
    defaultIntrRay3Ellipsoid3FIResult,
    intrRay3Ellipsoid3FIDoQuery
} from '../src/IntrRay3Ellipsoid3.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { mulMatrix } from '../src/Matrix.js';
import { intrLine3Ellipsoid3FIDoQuery } from '../src/IntrLine3Ellipsoid3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, ray as arbRay,
    rotationFrame, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function ellipsoid(center: number[], extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)],
        Vector.fromArray(extent));
}

// The value of (X-C)^T M (X-C) - 1: zero on the ellipsoid, negative inside.
function level(E: Hyperellipsoid, X: Vector): number {
    const d = sub(X, E.center);
    let sum = -1;
    for (let i = 0; i < 3; ++i) {
        const t = d.values[i] / E.extent.values[i];
        sum += t * t;
    }
    return sum;
}

const ti = new IntrRay3Ellipsoid3TI();
const fi = new IntrRay3Ellipsoid3FI();

describe('IntrRay3Ellipsoid3', () => {
    it('has an empty default result', () => {
        const result = defaultIntrRay3Ellipsoid3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds both crossings when the ray starts outside', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const result = fi.find(ray([-10, 0, 0], [1, 0, 0]), E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(7, 12);
        expect(result.parameter[1]).toBeCloseTo(13, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-3, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(ray([-10, 0, 0], [1, 0, 0]), E).intersect).toBe(true);
    });

    it('clips the near crossing when the ray starts inside', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([0, 0, 0], [1, 0, 0]);
        const result = fi.find(R, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(ti.test(R, E).intersect).toBe(true);
    });

    it('reports a grazing ray at the tangent point', () => {
        // The ray y = 2 grazes the ellipsoid at (0,2,0). Rounding decides
        // whether the discriminant is zero or a tiny positive number, so the
        // query reports one or two nearly identical crossings.
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([-10, 2, 0], [1, 0, 0]);
        const result = fi.find(R, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBeGreaterThanOrEqual(1);
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.point[i].values[0]).toBeCloseTo(0, 6);
            expect(result.point[i].values[1]).toBeCloseTo(2, 12);
        }
        expect(ti.test(R, E).intersect).toBe(true);
    });

    it('reports no intersection for a ray pointing away', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const R = ray([-10, 0, 0], [-1, 0, 0]);
        expect(fi.find(R, E).intersect).toBe(false);
        expect(ti.test(R, E).intersect).toBe(false);
        // A ray that misses the ellipsoid entirely.
        const M = ray([-10, 5, 0], [1, 0, 0]);
        expect(fi.find(M, E).intersect).toBe(false);
        expect(ti.test(M, E).intersect).toBe(false);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrRay3Ellipsoid3FIResult();
        intrRay3Ellipsoid3FIDoQuery(vec(-10, 0, 0), vec(1, 0, 0),
            ellipsoid([0, 0, 0], [3, 2, 1]), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(7, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 112233;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        let numInside = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const E = ellipsoid([rand(), rand(), rand()],
                [0.5 + Math.abs(rand()) * 2, 0.5 + Math.abs(rand()) * 2,
                    0.5 + Math.abs(rand()) * 2]);
            const R = ray([rand() * 4, rand() * 4, rand() * 4],
                [rand(), rand(), rand() + 0.001]);
            const result = fi.find(R, E);
            expect(ti.test(R, E).intersect).toBe(result.intersect);

            if (result.intersect) {
                ++numHits;
                const originInside = level(E, R.origin) <= 0;
                if (originInside) {
                    ++numInside;
                    expect(result.parameter[0]).toBe(0);
                }
                for (let i = 0; i < result.numIntersections; ++i) {
                    expect(result.parameter[i]).toBeGreaterThanOrEqual(0);
                    const onRay = add(R.origin,
                        mul(result.parameter[i], R.direction));
                    expect(length(sub(result.point[i], onRay)))
                        .toBeCloseTo(0, 10);
                    // The point is on the ellipsoid unless it is the clipped
                    // ray origin.
                    if (!(i === 0 && result.parameter[0] === 0
                        && originInside)) {
                        expect(level(E, result.point[i])).toBeCloseTo(0, 8);
                    }
                }
            }
        }
        expect(numHits).toBeGreaterThan(15);
        expect(numInside).toBeGreaterThan(0);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrRay3Ellipsoid3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrRay3Ellipsoid3 verification', () => {
    const tiQ = new IntrRay3Ellipsoid3TI();
    const fiQ = new IntrRay3Ellipsoid3FI();

    const arbEllipsoid = fc.tuple(wellScaledVector(3), rotationFrame(3),
        fc.tuple(positive(4, 0.05), positive(4, 0.05), positive(4, 0.05)))
        .map(([c, axis, e]) => Hyperellipsoid.fromCenterAxisExtent(c, axis,
            Vector.fromArray([e[0], e[1], e[2]])));

    // Q(X) = (X-C)^T M (X-C) - 1; negative inside, zero on the surface.
    function quadratic(E: Hyperellipsoid, X: Vector): number {
        const d = sub(X, E.center);
        return dot(d, mulMatrix(E.getM(), d) as Vector) - 1;
    }

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            expect(tiQ.test(R, E).intersect).toBe(fiQ.find(R, E).intersect);
        });
    });

    it('the FI parameters are nonnegative and their points are on the ray', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            const r = fiQ.find(R, E);
            if (!r.intersect) {
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.numIntersections).toBeGreaterThan(0);
            for (let i = 0; i < r.numIntersections; ++i) {
                expect(r.parameter[i]).toBeGreaterThanOrEqual(0);
                expectVectorClose(r.point[i],
                    add(R.origin, mul(r.parameter[i], R.direction)),
                    1e-12, 1e-12);
            }
            if (r.numIntersections === 2) {
                expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            }
        });
    });

    it('each reported point is on the ellipsoid or is the ray origin', () => {
        // The FI query clips the line t-interval to [0, +infinity), so the
        // first endpoint is either a surface point or the ray origin (when
        // the origin is inside the solid ellipsoid).
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            const r = fiQ.find(R, E);
            if (!r.intersect) { return; }
            const scale = 1 + length(E.extent) + length(sub(R.origin,
                E.center));
            for (let i = 0; i < r.numIntersections; ++i) {
                const q = quadratic(E, r.point[i]);
                const atOrigin = r.parameter[i] === 0;
                if (!atOrigin) {
                    expectClose(q, 0, 1e-7 * scale, 1e-7);
                } else {
                    // The clipped endpoint is inside or on the surface.
                    expect(q).toBeLessThanOrEqual(1e-9 * scale);
                }
            }
        });
    });

    it('agrees with a dense sampling of the ray against the solid', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            const maxExtent = Math.max(E.extent.values[0], E.extent.values[1],
                E.extent.values[2]);
            const reach = length(sub(E.center, R.origin)) + maxExtent + 1;
            let inside = false;
            for (let i = 0; i <= 4000; ++i) {
                const t = (reach * i) / 4000;
                if (quadratic(E, add(R.origin, mul(t, R.direction))) < 0) {
                    inside = true;
                    break;
                }
            }
            if (inside) { expect(tiQ.test(R, E).intersect).toBe(true); }
        }, 60);
    });

    it('is the line query clipped to t >= 0', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            const line = defaultIntrRay3Ellipsoid3FIResult();
            intrLine3Ellipsoid3FIDoQuery(R.origin, R.direction, E, line);
            const r = fiQ.find(R, E);
            if (!line.intersect) {
                expect(r.intersect).toBe(false);
                return;
            }
            const t0 = line.parameter[0];
            const t1 = line.parameter[1];
            if (t1 < 0) {
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.intersect).toBe(true);
            expect(r.parameter[0]).toBe(Math.max(t0, 0));
            expect(r.parameter[1]).toBe(t1);
            expect(r.numIntersections).toBe(
                Math.max(t0, 0) < t1 ? 2 : 1);
        });
    });

    it('reports a ray starting inside the ellipsoid', () => {
        check(fc.tuple(arbEllipsoid, unitVector(3),
            fc.double({ min: 0, max: 0.8, noNaN: true }), unitVector(3)),
            ([E, u, frac, d]) => {
                // A point strictly inside: C + frac * (extent-scaled u).
                const inside = add(E.center, add(
                    mul(frac * E.extent.values[0] * u.values[0], E.axis[0]),
                    add(mul(frac * E.extent.values[1] * u.values[1], E.axis[1]),
                        mul(frac * E.extent.values[2] * u.values[2],
                            E.axis[2]))));
                expect(quadratic(E, inside)).toBeLessThan(0);
                const R = Ray.fromOriginDirection(inside, d);
                expect(tiQ.test(R, E).intersect).toBe(true);
                const r = fiQ.find(R, E);
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(2);
                expect(r.parameter[0]).toBe(0);
                expectVectorClose(r.point[0], inside, 1e-12, 1e-12);
                const scale = 1 + length(E.extent);
                expectClose(quadratic(E, r.point[1]), 0, 1e-7 * scale, 1e-7);
            });
    });

    it('reports a ray pointing away from the ellipsoid as empty', () => {
        check(fc.tuple(arbEllipsoid, unitVector(3), positive(5, 1)),
            ([E, d, extra]) => {
                const maxExtent = Math.max(E.extent.values[0],
                    E.extent.values[1], E.extent.values[2]);
                const origin = add(E.center, mul(maxExtent + extra + 1, d));
                const R = Ray.fromOriginDirection(origin, d);
                expect(tiQ.test(R, E).intersect).toBe(false);
                const r = fiQ.find(R, E);
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
            });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid, rotationFrame(3),
            wellScaledVector(3)),
            ([R, E, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v), dot(Rot[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const R2 = Ray.fromOriginDirection(map(R.origin),
                    rot(R.direction));
                const E2 = Hyperellipsoid.fromCenterAxisExtent(map(E.center),
                    E.axis.map(rot), E.extent);
                const r0 = fiQ.find(R, E);
                const r1 = fiQ.find(R2, E2);
                if (r0.numIntersections === 1
                    || r1.numIntersections === 1) {
                    return;  // tangency or clipped endpoint
                }
                expect(r1.intersect).toBe(r0.intersect);
                expect(r1.numIntersections).toBe(r0.numIntersections);
                for (let i = 0; i < r0.numIntersections; ++i) {
                    expectClose(r1.parameter[i], r0.parameter[i], 1e-6, 1e-6);
                }
            });
    });

    it('the exported DoQuery agrees with find on the parameters', () => {
        check(fc.tuple(arbRay(3), arbEllipsoid), ([R, E]) => {
            const d = defaultIntrRay3Ellipsoid3FIResult();
            intrRay3Ellipsoid3FIDoQuery(R.origin, R.direction, E, d);
            const f = fiQ.find(R, E);
            expect(d.intersect).toBe(f.intersect);
            expect(d.numIntersections).toBe(f.numIntersections);
            expect(d.parameter[0]).toBe(f.parameter[0]);
            expect(d.parameter[1]).toBe(f.parameter[1]);
        });
    });
});
