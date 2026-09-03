import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2.js';
import { DistRay2Circle2 } from '../src/DistRay2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistRay2Circle2', () => {
    const query = new DistRay2Circle2();
    const unit = circle([0, 0], 1);

    it('reports two closest pairs when the ray crosses the circle', () => {
        const result = query.compute(ray([-5, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(-1, 10);
        expect(result.closest[1][1].values[0]).toBeCloseTo(1, 10);
    });

    it('drops the t0-point when the ray origin is inside the circle', () => {
        const result = query.compute(ray([0, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 10);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('uses the ray origin when the ray points away from the circle', () => {
        // The ray origin at (5,0) pointing to +x; the whole line crosses the
        // circle only at negative parameters.
        const result = query.compute(ray([5, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0][0].values).toEqual([5, 0]);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('matches the point-circle query when the ray points away', () => {
        const pcQuery = new DistPoint2Circle2();
        const origin = v(3, 4);
        const result = query.compute(
            Ray.fromOriginDirection(origin, v(1, 1)), unit);
        const pcResult = pcQuery.compute(origin, unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(pcResult.distance, 12);
        for (let i = 0; i < 2; ++i) {
            expect(result.closest[0][1].values[i]).toBeCloseTo(
                pcResult.closest[1].values[i], 12);
        }
    });

    it('handles a line that misses the circle with a positive parameter',
        () => {
            const result = query.compute(ray([-5, 3], [1, 0]), unit);
            expect(result.numClosestPairs).toBe(1);
            expect(result.parameter[0]).toBeGreaterThan(0);
            expect(result.distance).toBeCloseTo(2, 10);
        });

    it('handles a tangent ray', () => {
        const result = query.compute(ray([-5, 1], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('agrees with a dense sampling of the ray and circle', () => {
        const rnd = makeRandom(70001);
        const center = v(0.5, -0.75);
        const radius = 1.25;
        const c = Hypersphere.fromCenterRadius(center, radius);

        for (let trial = 0; trial < 60; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r = Ray.fromOriginDirection(origin, dir);
            const result = query.compute(r, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                expect(result.parameter[j]).toBeGreaterThanOrEqual(-1e-12);
                const p = add(r.origin, mul(result.parameter[j], r.direction));
                for (let i = 0; i < 2; ++i) {
                    expect(p.values[i]).toBeCloseTo(
                        result.closest[j][0].values[i], 7);
                }
                const d = sub(result.closest[j][1], center);
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 7);
                const e = sub(result.closest[j][0], result.closest[j][1]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 6);
            }

            // No sampled (ray point, circle point) pair is closer.
            const n = 3000;
            const dd = dot(r.direction, r.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, v(radius * Math.cos(t),
                    radius * Math.sin(t)));
                const w = sub(q, r.origin);
                let s = dot(w, r.direction) / dd;
                if (s < 0) {
                    s = 0;
                }
                const f = sub(w, mul(s, r.direction));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});
