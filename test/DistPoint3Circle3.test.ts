import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistPoint3Circle3 } from '../src/DistPoint3Circle3.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPoint3Circle3', () => {
    const query = new DistPoint3Circle3();

    it('measures a point in the plane of the circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 2);
        const result = query.compute(v(5, 0, 0), c);
        expect(result.equidistant).toBe(false);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a point inside the circle in its plane', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 5);
        const result = query.compute(v(1, 0, 0), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(5, 12);
    });

    it('measures a point off the plane (3-4-5 configuration)', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 1);
        // Radial distance 4 from the circle, height 3.
        const result = query.compute(v(5, 0, 3), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports zero distance for a point on the circle', () => {
        const c = circle([1, 2, 3], [0, 1, 0], 4);
        const result = query.compute(v(5, 2, 3), c);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('reports equidistance for a point on the normal line', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 3);
        const result = query.compute(v(0, 0, 4), c);
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBeCloseTo(5, 12);
        // The closest point is on the circle.
        const d = sub(result.closest[1], c.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(3, 10);
        expect(Math.abs(dot(d, c.normal))).toBeLessThan(1e-12);
    });

    it('reports equidistance at the circle center', () => {
        const c = circle([2, -1, 4], [1, 1, 1], 2);
        const result = query.compute(v(2, -1, 4), c);
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBeCloseTo(2, 12);
    });

    it('handles a zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const result = query.compute(v(3, 4, 12), c);
        expect(result.distance).toBeCloseTo(13, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 10);
    });

    it('is accurate for a point nearly on the normal line', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 1);
        const result = query.compute(v(1e-14, 0, 1), c);
        expect(result.equidistant).toBe(false);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 10);
        const d = sub(result.closest[1], c.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(1, 10);
    });

    it('agrees with a dense sampling of the circle', () => {
        const rnd = makeRandom(31415);
        const center = v(0.5, -1, 2);
        const normal = v(1, 2, -0.5);
        normalize(normal);
        const radius = 1.5;
        const c = Circle3.fromCenterNormalRadius(center, normal, radius);

        // An orthonormal basis of the plane of the circle.
        const U = v(-normal.values[1], normal.values[0], 0);
        normalize(U);
        const W = cross(normal, U);

        for (let trial = 0; trial < 50; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, c);

            // The reported closest point is on the circle.
            const d = sub(result.closest[1], center);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 8);
            expect(Math.abs(dot(d, normal))).toBeLessThan(1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 8);

            // No sampled circle point is closer.
            const n = 5000;
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, add(mul(radius * Math.cos(t), U),
                    mul(radius * Math.sin(t), W)));
                const f = sub(p, q);
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
