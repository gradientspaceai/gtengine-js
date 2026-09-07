import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrRay3Plane3TI,
    IntrRay3Plane3FI,
    defaultIntrRay3Plane3FIResult,
    defaultIntrRay3Plane3TIResult,
    intrRay3Plane3FIDoQuery
} from '../src/IntrRay3Plane3.js';
import { IntrLine3Plane3FI } from '../src/IntrLine3Plane3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

const ti = new IntrRay3Plane3TI();
const fi = new IntrRay3Plane3FI();

describe('IntrRay3Plane3', () => {
    it('default-constructs results as no intersection', () => {
        expect(defaultIntrRay3Plane3TIResult()).toEqual({ intersect: false });
        const r = defaultIntrRay3Plane3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
        expect(r.parameter).toBe(0);
        expect(r.point.values).toEqual([0, 0, 0]);
    });

    it('finds a transverse intersection at a known point', () => {
        const R = ray([0, 0, 5], [0, 0, -1]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(R, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.point.values).toEqual([0, 0, 2]);
        expect(ti.test(R, P).intersect).toBe(true);
    });

    it('rejects the intersection behind the ray origin', () => {
        // The line hits the plane at t = -3, which is not on the ray.
        const R = ray([0, 0, 5], [0, 0, 1]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        const result = fi.find(R, P);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(R, P).intersect).toBe(false);
    });

    it('accepts an origin exactly on the plane with t = 0', () => {
        const P = plane([0, 0, 1], [0, 0, 2]);
        for (const dir of [[0, 0, 1], [0, 0, -1], [1, 0, 0]]) {
            const R = ray([1, -1, 2], dir);
            const result = fi.find(R, P);
            expect(result.intersect).toBe(true);
            expect(ti.test(R, P).intersect).toBe(true);
            if (dir[2] !== 0) {
                expect(result.numIntersections).toBe(1);
                expect(result.parameter).toBeCloseTo(0, 12);
                expect(result.point.values[0]).toBeCloseTo(1, 12);
                expect(result.point.values[1]).toBeCloseTo(-1, 12);
                expect(result.point.values[2]).toBeCloseTo(2, 12);
            }
            else {
                // The ray is on the plane.
                expect(result.numIntersections).toBe(2147483647);
            }
        }
    });

    it('reports a parallel disjoint ray as no intersection', () => {
        const R = ray([0, 0, 5], [1, 0, 0]);
        const P = plane([0, 0, 1], [0, 0, 2]);
        expect(fi.find(R, P).intersect).toBe(false);
        expect(ti.test(R, P).intersect).toBe(false);
    });

    it('matches upstream TI logic on the sign of the ray direction', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // Origin above the plane (signedDistance > 0).
        expect(ti.test(ray([0, 0, 3], [0, 0, -1]), P).intersect).toBe(true);
        expect(ti.test(ray([0, 0, 3], [0, 0, 1]), P).intersect).toBe(false);
        // Origin below the plane (signedDistance < 0).
        expect(ti.test(ray([0, 0, -3], [0, 0, 1]), P).intersect).toBe(true);
        expect(ti.test(ray([0, 0, -3], [0, 0, -1]), P).intersect).toBe(false);
    });

    it('exposes the DoQuery helper used by derived queries', () => {
        const P = plane([1, 1, 1], [0, 0, 0]);
        const result = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(vec(1, 1, 1), vec(-1, -1, -1), P, result);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(1, 12);

        const behind = defaultIntrRay3Plane3FIResult();
        intrRay3Plane3FIDoQuery(vec(1, 1, 1), vec(1, 1, 1), P, behind);
        expect(behind.intersect).toBe(false);
        expect(behind.numIntersections).toBe(0);
    });

    it('agrees with TI and puts the point on both primitives (randomized)', () => {
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 6 * rand() - 3;

        let hits = 0;
        for (let trial = 0; trial < 4000; ++trial) {
            const R = ray([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const P = plane([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const fiResult = fi.find(R, P);
            const tiResult = ti.test(R, P);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            if (fiResult.intersect && fiResult.numIntersections === 1) {
                ++hits;
                expect(fiResult.parameter).toBeGreaterThanOrEqual(0);
                // The point is on the ray.
                const onRay = add(R.origin, mul(fiResult.parameter, R.direction));
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.point.values[i]).toBeCloseTo(
                        onRay.values[i], 10);
                }
                // The point is on the plane.
                expect(dot(P.normal, fiResult.point) - P.constant)
                    .toBeCloseTo(0, 8);
            }
        }
        expect(hits).toBeGreaterThan(1000);
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): property-based cross-checks against the line query and
// against the geometry of the reported point.
// ---------------------------------------------------------------------------

describe('IntrRay3Plane3 verification', () => {
    const ti = new IntrRay3Plane3TI();
    const fiv = new IntrRay3Plane3FI();
    const lineFI = new IntrLine3Plane3FI();
    const INT32_MAX = 2147483647;

    // A ray and a plane that are both moderately scaled, so that a
    // catastrophic cancellation in Dot(N, O) - c cannot make the sign tests
    // meaningless.
    const rayPlane = fc.tuple(unitVector(3), wellScaledVector(3, -5, 5),
        unitVector(3), wellScaledVector(3, -5, 5))
        .map(([n, po, d, ro]) => ({
            plane: Hyperplane.fromNormalOrigin(n, po),
            ray: Ray.fromOriginDirection(ro, d)
        }));

    it('TI and FI agree on intersect', () => {
        check(rayPlane, ({ ray: R, plane: P }) => {
            // The two queries are algebraically the same predicate: the FI
            // parameter is -signedDistance/DdN, which is nonnegative exactly
            // when the TI sign test succeeds. No tolerance is needed.
            const a = ti.test(R, P).intersect;
            const b = fiv.find(R, P).intersect;
            expect(b).toBe(a);
        });
    });

    it('equals the line query restricted to t >= 0', () => {
        check(rayPlane, ({ ray: R, plane: P }) => {
            const lr = lineFI.find(
                Line.fromOriginDirection(R.origin, R.direction), P);
            const rr = fiv.find(R, P);
            const expected = lr.intersect && !(lr.parameter < 0);
            expect(rr.intersect).toBe(expected);
            if (expected) {
                expect(rr.numIntersections).toBe(lr.numIntersections);
                expect(rr.parameter).toBe(lr.parameter);
            }
            else {
                expect(rr.numIntersections).toBe(0);
            }
        });
    });

    it('reports a point on the plane and on the ray', () => {
        check(rayPlane, ({ ray: R, plane: P }) => {
            const r = fiv.find(R, P);
            if (!r.intersect) {
                return;
            }
            for (let i = 0; i < 3; ++i) {
                expect(Number.isNaN(r.point.get(i))).toBe(false);
            }
            expect(r.parameter >= 0).toBe(true);
            if (r.numIntersections === 1) {
                // |N| = 1, so Dot(N,X) - c is the signed distance. The point
                // is O + t*D with |t| <= |signedDistance(O)|/|DdN|; the
                // residual is bounded by the rounding of that product.
                const scale = 1 + Math.abs(r.parameter);
                expectClose(dot(P.normal, r.point) - P.constant, 0,
                    1e-12 * scale, 1e-12);
                expectVectorClose(r.point,
                    add(R.origin, mul(r.parameter, R.direction)), 1e-12, 1e-12);
            }
            else {
                // The ray lies on the plane.
                expect(r.numIntersections).toBe(INT32_MAX);
                expect(r.parameter).toBe(0);
                expectVectorClose(r.point, R.origin);
            }
        });
    });

    it('reports the whole ray when the ray lies exactly in the plane', () => {
        // The ray must lie in the plane exactly, not merely to within
        // rounding: for a nearly coplanar ray the FI query divides a tiny
        // signed distance by a tiny Dot(D,N) and the resulting parameter is
        // arbitrary, so the classification is a knife edge. An axis-aligned
        // plane and points sharing that coordinate give exact zeros.
        check(fc.tuple(fc.integer({ min: 0, max: 2 }),
            wellScaledVector(3, -5, 5), wellScaledVector(3, -5, 5),
            wellScaledVector(3, -5, 5)),
            ([k, po, a, b]) => {
                const o = a.clone();
                const q = b.clone();
                o.values[k] = po.values[k];
                q.values[k] = po.values[k];
                const d = sub(q, o);
                if (length(d) < 1e-2) {
                    return;
                }
                normalize(d);
                // The k-th component of d is exactly zero, so Dot(D,N) = 0.
                const P = Hyperplane.fromNormalOrigin(Vector.unit(3, k), po);
                const R = Ray.fromOriginDirection(o, d);
                const r = fiv.find(R, P);
                expect(ti.test(R, P).intersect).toBe(true);
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(INT32_MAX);
                expect(r.parameter).toBe(0);
                expectVectorClose(r.point, R.origin);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(rayPlane, rotationFrame(3), wellScaledVector(3, -5, 5)),
            ([{ ray: R, plane: P }, frame, t]) => {
                const xf = (v: Vector): Vector => {
                    const w = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        w.values[i] = frame[0].values[i] * v.values[0]
                            + frame[1].values[i] * v.values[1]
                            + frame[2].values[i] * v.values[2]
                            + t.values[i];
                    }
                    return w;
                };
                const xfDir = (v: Vector): Vector => {
                    const w = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        w.values[i] = frame[0].values[i] * v.values[0]
                            + frame[1].values[i] * v.values[1]
                            + frame[2].values[i] * v.values[2];
                    }
                    return w;
                };
                const r0 = fiv.find(R, P);
                const r1 = fiv.find(
                    Ray.fromOriginDirection(xf(R.origin), xfDir(R.direction)),
                    Hyperplane.fromNormalOrigin(xfDir(P.normal),
                        xf(add(mul(P.constant, P.normal), new Vector(3)))));
                // Grazing configurations (the transformed Dot(D,N) landing on
                // the other side of zero) are knife edges; skip them.
                if (Math.abs(dot(R.direction, P.normal)) < 1e-6) {
                    return;
                }
                expect(r1.intersect).toBe(r0.intersect);
                if (r0.intersect && r0.numIntersections === 1) {
                    expectClose(r1.parameter, r0.parameter, 1e-8, 1e-8);
                    expectVectorClose(r1.point, xf(r0.point), 1e-8, 1e-8);
                }
            });
    });
});
