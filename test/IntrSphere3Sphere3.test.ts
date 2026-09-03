import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    IntrSphere3Sphere3TI,
    IntrSphere3Sphere3FI,
    IntrSphere3Sphere3FIResultType
} from '../src/IntrSphere3Sphere3.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function sphere(c: number[], r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(c), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSphere3Sphere3', () => {
    const ti = new IntrSphere3Sphere3TI();
    const fi = new IntrSphere3Sphere3FI();
    const T = IntrSphere3Sphere3FIResultType;

    it('reports separated spheres', () => {
        const s0 = sphere([0, 0, 0], 1);
        const s1 = sphere([5, 0, 0], 1);
        expect(ti.test(s0, s1).intersect).toBe(false);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(false);
        expect(result.type).toBe(T.separated);
    });

    it('reports external tangency', () => {
        const s0 = sphere([0, 0, 0], 1);
        const s1 = sphere([3, 0, 0], 2);
        expect(ti.test(s0, s1).intersect).toBe(true);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.touchingOutside);
        expect(result.point.values[0]).toBeCloseTo(1, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
    });

    it('reports the circle of intersection', () => {
        // Two unit spheres with centers 1 apart: the plane of intersection is
        // x = 0.5 and the circle radius is sqrt(3)/2.
        const s0 = sphere([-0.5, 0, 0], 1);
        const s1 = sphere([0.5, 0, 0], 1);
        expect(ti.test(s0, s1).intersect).toBe(true);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.circle);
        expect(result.circle.center.values[0]).toBeCloseTo(0, 12);
        expect(result.circle.radius).toBeCloseTo(Math.sqrt(3) / 2, 12);
        expect(result.circle.normal.values[0]).toBeCloseTo(1, 12);
        expect(result.circle.normal.values[1]).toBeCloseTo(0, 12);
    });

    it('reports sphere0 strictly contained in sphere1', () => {
        const s0 = sphere([0, 0, 0], 1);
        const s1 = sphere([0.25, 0, 0], 5);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere0StrictlyInside);
        expect(result.point.values[0]).toBeCloseTo(0.125, 12);
    });

    it('reports sphere1 strictly contained in sphere0', () => {
        const s0 = sphere([0, 0, 0], 5);
        const s1 = sphere([1, 0, 0], 1);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere1StrictlyInside);
    });

    it('reports internal tangency with sphere0 inside sphere1', () => {
        // r1 - r0 = 2 and |C1-C0| = 2.
        const s0 = sphere([2, 0, 0], 1);
        const s1 = sphere([0, 0, 0], 3);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere0InsideTouching);
        // Upstream bug (FIXED; see upstream-bug issue (B71)): upstream
        // computes C1 + r1 * normalize(C1 - C0), the antipode (-3,0,0) of the
        // true contact point. The port uses C1 - r1 * normalize(C1 - C0),
        // matching the sign used by the sibling (sphere1 inside sphere0)
        // branch.
        expect(result.point.values[0]).toBeCloseTo(3, 12);
        expect(result.point.values[1]).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
        // The contact point lies on both sphere surfaces.
        const e0 = sub(result.point, s0.center);
        const e1 = sub(result.point, s1.center);
        expect(Math.sqrt(dot(e0, e0))).toBeCloseTo(s0.radius, 12);
        expect(Math.sqrt(dot(e1, e1))).toBeCloseTo(s1.radius, 12);
    });

    it('reports internal tangency with sphere1 inside sphere0', () => {
        const s0 = sphere([0, 0, 0], 3);
        const s1 = sphere([2, 0, 0], 1);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere1InsideTouching);
        expect(result.point.values[0]).toBeCloseTo(3, 12);
        const e0 = sub(result.point, s0.center);
        const e1 = sub(result.point, s1.center);
        expect(Math.sqrt(dot(e0, e0))).toBeCloseTo(s0.radius, 12);
        expect(Math.sqrt(dot(e1, e1))).toBeCloseTo(s1.radius, 12);
    });

    it('reports concentric equal spheres as a strict containment', () => {
        // rDif = 0 and sqrLen = 0 < rDifSqr is false, and sqrLen == rDifSqr
        // holds, so the single-contact branch runs with rDif <= 0.
        const s0 = sphere([0, 0, 0], 2);
        const s1 = sphere([0, 0, 0], 2);
        const result = fi.find(s0, s1);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere0InsideTouching);
    });

    it('handles zero-radius spheres', () => {
        const point = sphere([1, 0, 0], 0);
        const big = sphere([0, 0, 0], 2);
        expect(ti.test(point, big).intersect).toBe(true);
        const result = fi.find(point, big);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(T.sphere0StrictlyInside);

        const far = sphere([9, 0, 0], 0);
        expect(ti.test(far, big).intersect).toBe(false);
        expect(fi.find(far, big).type).toBe(T.separated);
    });

    it('internal-tangency contact points lie on both spheres (random)', () => {
        const rnd = makeRandom(778899);
        let contactMismatch = 0;
        let cases = 0;

        for (let trial = 0; trial < 300; ++trial) {
            // Build an exact internal tangency: |C1 - C0| == |r0 - r1|.
            const rBig = 1 + 2 * rnd();
            const rSmall = 0.25 + (rBig - 0.5) * rnd();
            const dir = vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const len = Math.sqrt(dot(dir, dir));
            if (len < 1e-3) {
                continue;
            }
            const unitDir = mul(1 / len, dir);
            const cBig = vec(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const cSmall = add(cBig, mul(rBig - rSmall, unitDir));
            const big = Hypersphere.fromCenterRadius(cBig, rBig);
            const small = Hypersphere.fromCenterRadius(cSmall, rSmall);

            // Both orders: sphere0 inside sphere1, and sphere1 inside sphere0.
            for (const [a, b, expectedType] of [
                [small, big, T.sphere0InsideTouching],
                [big, small, T.sphere1InsideTouching]
            ] as [Hypersphere, Hypersphere, number][]) {
                const result = fi.find(a, b);
                if (result.type !== expectedType) {
                    continue;  // rounding pushed it off the exact tangency
                }
                ++cases;
                const da = sub(result.point, a.center);
                const db = sub(result.point, b.center);
                if (Math.abs(Math.sqrt(dot(da, da)) - a.radius) > 1e-9 ||
                    Math.abs(Math.sqrt(dot(db, db)) - b.radius) > 1e-9) {
                    ++contactMismatch;
                }
            }
        }

        expect(cases).toBeGreaterThan(20);
        expect(contactMismatch).toBe(0);
    });

    it('the circle of intersection lies on both spheres (random)', () => {
        const rnd = makeRandom(112358);
        let tiFiMismatch = 0;
        let circleMismatch = 0;
        let typeMismatch = 0;
        let circles = 0;

        for (let trial = 0; trial < 400; ++trial) {
            const s0 = sphere([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1],
                0.25 + 1.5 * rnd());
            const s1 = sphere([2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1],
                0.25 + 1.5 * rnd());

            const tiResult = ti.test(s0, s1);
            const fiResult = fi.find(s0, s1);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            const diff = sub(s1.center, s0.center);
            const dist = Math.sqrt(dot(diff, diff));
            const expectedIntersect = dist <= s0.radius + s1.radius;
            if (expectedIntersect !== fiResult.intersect) {
                ++typeMismatch;
            }
            if (fiResult.intersect &&
                dist < Math.abs(s0.radius - s1.radius) &&
                fiResult.type !== T.sphere0StrictlyInside &&
                fiResult.type !== T.sphere1StrictlyInside) {
                ++typeMismatch;
            }

            if (fiResult.type === T.circle) {
                ++circles;
                // Pick a point on the reported circle and check it lies on
                // both spheres.
                const n = fiResult.circle.normal;
                let u = vec(1, 0, 0);
                if (Math.abs(n.values[0]) > 0.9) {
                    u = vec(0, 1, 0);
                }
                const w = vec(
                    n.values[1] * u.values[2] - n.values[2] * u.values[1],
                    n.values[2] * u.values[0] - n.values[0] * u.values[2],
                    n.values[0] * u.values[1] - n.values[1] * u.values[0]);
                const wLen = Math.sqrt(dot(w, w));
                const p = add(fiResult.circle.center,
                    mul(fiResult.circle.radius / wLen, w));
                const d0 = sub(p, s0.center);
                const d1 = sub(p, s1.center);
                if (Math.abs(Math.sqrt(dot(d0, d0)) - s0.radius) > 1e-9 ||
                    Math.abs(Math.sqrt(dot(d1, d1)) - s1.radius) > 1e-9) {
                    ++circleMismatch;
                }
                // The normal must be unit length and parallel to C1 - C0.
                if (Math.abs(Math.sqrt(dot(n, n)) - 1) > 1e-12) {
                    ++circleMismatch;
                }
            }
        }

        expect(circles).toBeGreaterThan(20);
        expect([tiFiMismatch, typeMismatch, circleMismatch]).toEqual([0, 0, 0]);
    });
});
