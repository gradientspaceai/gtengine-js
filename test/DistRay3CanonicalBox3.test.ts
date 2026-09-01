import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox';
import { DistRay3CanonicalBox3 } from '../src/DistRay3CanonicalBox3';
import { Line } from '../src/Line';
import { Ray } from '../src/Ray';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistRay3CanonicalBox3', () => {
    const query = new DistRay3CanonicalBox3();
    const unitBox = CanonicalBox.fromExtent(v(1, 1, 1));

    it('reports zero distance for a ray that enters the box', () => {
        const result = query.compute(ray([-5, 0, 0], [1, 0, 0]), unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeGreaterThanOrEqual(0);
    });

    it('reports zero distance for a ray originating inside the box', () => {
        const result = query.compute(ray([0, 0, 0], [1, 2, 3]), unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('uses the ray origin when the ray points away from the box', () => {
        // The origin is at (5,0,0) and the ray points to +x.
        const result = query.compute(ray([5, 0, 0], [1, 0, 0]), unitBox);
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0].values).toEqual([5, 0, 0]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('matches the point-box query when the ray points away', () => {
        const pbQuery = new DistPointCanonicalBox();
        const origin = v(3, -4, 5);
        const result = query.compute(
            Ray.fromOriginDirection(origin, v(1, -1, 1)), unitBox);
        const pbResult = pbQuery.compute(origin, unitBox);
        expect(result.distance).toBeCloseTo(pbResult.distance, 12);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(
                pbResult.closest[1].values[i], 12);
        }
    });

    it('matches the line query when the closest line point is on the ray',
        () => {
            const lbQuery = new DistLine3CanonicalBox3();
            const r = ray([-5, 3, 0], [1, 0, 0]);
            const result = query.compute(r, unitBox);
            const lbResult = lbQuery.compute(
                Line.fromOriginDirection(r.origin, r.direction), unitBox);
            expect(lbResult.parameter).toBeGreaterThanOrEqual(0);
            expect(result.distance).toBeCloseTo(lbResult.distance, 12);
        });

    it('handles a degenerate box with zero extents', () => {
        const box = CanonicalBox.fromExtent(v(0, 0, 0));
        const result = query.compute(ray([0, 3, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the ray and box', () => {
        const rnd = makeRandom(51515);
        const extent = v(1, 0.75, 1.5);
        const box = CanonicalBox.fromExtent(extent);

        for (let trial = 0; trial < 25; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r = Ray.fromOriginDirection(origin, dir);
            const result = query.compute(r, box);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onRay.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
                expect(Math.abs(result.closest[1].values[i]))
                    .toBeLessThanOrEqual(extent.values[i] + 1e-8);
            }
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled (ray point, box point) pair is closer.
            const nb = 10;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= nb; ++i) {
                for (let j = 0; j <= nb; ++j) {
                    for (let k = 0; k <= nb; ++k) {
                        const q = v(
                            (2 * i / nb - 1) * extent.values[0],
                            (2 * j / nb - 1) * extent.values[1],
                            (2 * k / nb - 1) * extent.values[2]);
                        // The closest ray point to q, clamped to t >= 0.
                        const w = sub(q, r.origin);
                        let t = dot(w, r.direction)
                            / dot(r.direction, r.direction);
                        if (t < 0) {
                            t = 0;
                        }
                        const f = sub(w, mul(t, r.direction));
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-8);
        }
    });
});
