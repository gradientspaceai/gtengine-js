import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrPlane3Sphere3TI,
    IntrPlane3Sphere3FI,
    defaultIntrPlane3Sphere3FIResult
} from '../src/IntrPlane3Sphere3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, fc, plane as arbPlane, rotationFrame,
    sphere as arbSphere, wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function sphere(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

const ti = new IntrPlane3Sphere3TI();
const fi = new IntrPlane3Sphere3FI();

describe('IntrPlane3Sphere3', () => {
    it('has an all-zero invalid circle in the default result', () => {
        const result = defaultIntrPlane3Sphere3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.isCircle).toBe(false);
        expect(result.circle.center.values).toEqual([0, 0, 0]);
        expect(result.circle.normal.values).toEqual([0, 0, 0]);
        expect(result.circle.radius).toBe(0);
    });

    it('cuts a sphere in a circle of the expected radius', () => {
        // The plane z = 3 cuts the sphere of radius 5 centered at the origin
        // in a circle of radius 4.
        const P = plane([0, 0, 1], [0, 0, 3]);
        const S = sphere([0, 0, 0], 5);
        const result = fi.find(P, S);
        expect(result.intersect).toBe(true);
        expect(result.isCircle).toBe(true);
        expect(result.circle.radius).toBeCloseTo(4, 12);
        expect(result.circle.center.values[2]).toBeCloseTo(3, 12);
        expect(result.circle.normal.values).toEqual([0, 0, 1]);
        expect(ti.test(P, S).intersect).toBe(true);
    });

    it('reports a tangent plane as a single point', () => {
        const P = plane([0, 0, 1], [0, 0, 5]);
        const S = sphere([0, 0, 0], 5);
        const result = fi.find(P, S);
        expect(result.intersect).toBe(true);
        expect(result.isCircle).toBe(false);
        expect(result.point.values[0]).toBeCloseTo(0, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(5, 12);
        expect(ti.test(P, S).intersect).toBe(true);
    });

    it('reports separation when the plane misses the sphere', () => {
        const P = plane([0, 0, 1], [0, 0, 6]);
        const S = sphere([0, 0, 0], 5);
        expect(fi.find(P, S).intersect).toBe(false);
        expect(ti.test(P, S).intersect).toBe(false);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 987654321;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numCircles = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const S = sphere([rand() * 3, rand() * 3, rand() * 3],
                0.5 + 2 * Math.abs(rand()));
            const result = fi.find(P, S);
            expect(ti.test(P, S).intersect).toBe(result.intersect);

            const signedDistance = dot(P.normal, S.center) - P.constant;
            expect(result.intersect).toBe(Math.abs(signedDistance) <= S.radius);

            if (result.isCircle) {
                ++numCircles;
                // The circle center is on the plane and its radius satisfies
                // r^2 + d^2 = R^2.
                expect(dot(P.normal, result.circle.center) - P.constant)
                    .toBeCloseTo(0, 10);
                expect(result.circle.radius * result.circle.radius
                    + signedDistance * signedDistance)
                    .toBeCloseTo(S.radius * S.radius, 8);
                // A point on the circle is on the sphere.
                const toCenter = sub(result.circle.center, S.center);
                expect(dot(toCenter, toCenter)
                    + result.circle.radius * result.circle.radius)
                    .toBeCloseTo(S.radius * S.radius, 8);
            }
        }
        expect(numCircles).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3Sphere3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3Sphere3 verification', () => {
    const tiQ = new IntrPlane3Sphere3TI();
    const fiQ = new IntrPlane3Sphere3FI();

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbPlane(3), arbSphere(3)), ([P, S]) => {
            expect(fiQ.find(P, S).intersect).toBe(tiQ.test(P, S).intersect);
        });
    });

    it('the FI circle lies on the plane and on the sphere', () => {
        check(fc.tuple(arbPlane(3), arbSphere(3)), ([P, S]) => {
            const r = fiQ.find(P, S);
            if (!r.intersect || !r.isCircle) { return; }
            const scale = 1 + S.radius + length(S.center);
            // The circle normal is the plane normal and the circle center is
            // on the plane.
            expect(r.circle.normal.values).toEqual(P.normal.values);
            expectClose(dot(P.normal, r.circle.center) - P.constant, 0,
                1e-13 * scale, 0);
            // Every circle point is on the sphere: sample a few using an
            // orthonormal basis of the plane.
            const helper = Math.abs(P.normal.values[0]) < 0.9
                ? Vector.fromArray([1, 0, 0]) : Vector.fromArray([0, 1, 0]);
            const u = cross(P.normal, helper);
            normalize(u);
            const v = cross(P.normal, u);
            for (let i = 0; i < 8; ++i) {
                const a = (2 * Math.PI * i) / 8;
                const X = add(r.circle.center,
                    add(mul(r.circle.radius * Math.cos(a), u),
                        mul(r.circle.radius * Math.sin(a), v)));
                expectClose(length(sub(X, S.center)), S.radius,
                    1e-11 * scale, 1e-11);
                expectClose(dot(P.normal, X) - P.constant, 0,
                    1e-11 * scale, 0);
            }
        });
    });

    it('the FI tangent point is on the plane and on the sphere', () => {
        check(fc.tuple(arbPlane(3), wellScaledVector(3)), ([P, c]) => {
            // Build a sphere whose radius is exactly the point-plane
            // distance, so the '==' branch is taken.
            const radius = Math.abs(dot(P.normal, c) - P.constant);
            if (radius < 1e-6) { return; }
            const S = Hypersphere.fromCenterRadius(c, radius);
            const r = fiQ.find(P, S);
            expect(r.intersect).toBe(true);
            expect(r.isCircle).toBe(false);
            const scale = 1 + radius + length(c);
            expectClose(dot(P.normal, r.point) - P.constant, 0,
                1e-12 * scale, 0);
            expectClose(length(sub(r.point, c)), radius, 1e-12 * scale, 1e-12);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbSphere(3), rotationFrame(3),
            wellScaledVector(3)),
            ([P, S, R, t]) => {
                const rot = (x: Vector) => Vector.fromArray([
                    dot(R[0], x), dot(R[1], x), dot(R[2], x)]);
                const map = (x: Vector) => add(rot(x), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const S2 = Hypersphere.fromCenterRadius(map(S.center),
                    S.radius);
                const r0 = fiQ.find(P, S);
                const r1 = fiQ.find(P2, S2);
                const d = Math.abs(dot(P.normal, S.center) - P.constant);
                if (Math.abs(d - S.radius) < 1e-9 * (1 + S.radius)) { return; }
                expect(r1.intersect).toBe(r0.intersect);
                expect(r1.isCircle).toBe(r0.isCircle);
                if (r0.isCircle) {
                    expectClose(r1.circle.radius, r0.circle.radius, 1e-9,
                        1e-9);
                }
            });
    });

    it('reports a sphere entirely on one side as no intersection', () => {
        check(fc.tuple(arbPlane(3), arbSphere(3)), ([P, S]) => {
            const d = Math.abs(dot(P.normal, S.center) - P.constant);
            if (d <= S.radius) { return; }
            const r = fiQ.find(P, S);
            expect(r.intersect).toBe(false);
            expect(r.isCircle).toBe(false);
            expect(r.circle.radius).toBe(0);
            expect(r.point.values).toEqual([0, 0, 0]);
        });
    });

    it('reports a plane through the sphere center as a great circle', () => {
        check(fc.tuple(arbSphere(3), wellScaledVector(3)), ([S, n]) => {
            if (dot(n, n) < 1e-6) { return; }
            const N = n.clone();
            normalize(N);
            const P = Hyperplane.fromNormalOrigin(N, S.center);
            const r = fiQ.find(P, S);
            expect(r.intersect).toBe(true);
            expect(r.isCircle).toBe(true);
            expectClose(r.circle.radius, S.radius, 1e-12, 1e-12);
        });
    });
});
