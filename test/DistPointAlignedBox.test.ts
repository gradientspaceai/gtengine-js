import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { Vector, sub, dot } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// A simple deterministic PRNG so the randomized checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPointAlignedBox', () => {
    const query = new DistPointAlignedBox();

    it('reports zero distance for a point inside the box', () => {
        const result = query.compute(v(0.25, -0.5, 0.75),
            box([-1, -1, -1], [1, 1, 1]));
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(result.closest[1].values).toEqual([0.25, -0.5, 0.75]);
    });

    it('reports zero distance for a point on the box boundary', () => {
        const result = query.compute(v(1, 0, 0), box([-1, -1, -1], [1, 1, 1]));
        expect(result.distance).toBe(0);
    });

    it('measures a face-region point', () => {
        // The box [0,2]^3 and a point above the +z face.
        const result = query.compute(v(1, 1, 5), box([0, 0, 0], [2, 2, 2]));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([1, 1, 2]);
    });

    it('measures an edge-region point', () => {
        const result = query.compute(v(5, 1, 6), box([0, 0, 0], [2, 2, 2]));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([2, 1, 2]);
    });

    it('measures a corner-region point', () => {
        const result = query.compute(v(-3, -4, -12), box([0, 0, 0], [2, 2, 2]));
        expect(result.distance).toBeCloseTo(13, 12);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('handles a degenerate box with zero extents (a single point)', () => {
        const result = query.compute(v(3, 4, 0), box([1, 1, 1], [1, 1, 1]));
        expect(result.closest[1].values).toEqual([1, 1, 1]);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9 + 1), 12);
    });

    it('works in 2D', () => {
        const result = query.compute(v(4, 0), box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([1, 0]);
    });

    it('is translation invariant', () => {
        const p = v(3, -2, 7);
        const shift = v(10, -5, 2);
        const r0 = query.compute(p, box([-1, -2, -3], [4, 5, 6]));
        const r1 = query.compute(Vector.fromArray([13, -7, 9]),
            box([9, -7, -1], [14, 0, 8]));
        expect(r1.distance).toBeCloseTo(r0.distance, 12);
        for (let i = 0; i < 3; ++i) {
            expect(r1.closest[1].values[i]).toBeCloseTo(
                r0.closest[1].values[i] + shift.values[i], 12);
        }
    });

    it('agrees with a dense brute-force sampling of the box', () => {
        const rnd = makeRandom(12345);
        const b = box([-1, -0.5, -2], [2, 1.5, 0.5]);
        for (let trial = 0; trial < 40; ++trial) {
            const p = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const result = query.compute(p, b);

            // The reported closest point is on the box.
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeGreaterThanOrEqual(
                    b.min.values[i] - 1e-12);
                expect(result.closest[1].values[i]).toBeLessThanOrEqual(
                    b.max.values[i] + 1e-12);
            }

            // The reported closest point realizes the reported distance.
            const d = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(result.distance, 10);

            // No sampled box point is closer.
            const n = 12;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    for (let k = 0; k <= n; ++k) {
                        const q = v(
                            b.min.values[0] + (i / n)
                                * (b.max.values[0] - b.min.values[0]),
                            b.min.values[1] + (j / n)
                                * (b.max.values[1] - b.min.values[1]),
                            b.min.values[2] + (k / n)
                                * (b.max.values[2] - b.min.values[2]));
                        const e = sub(p, q);
                        best = Math.min(best, dot(e, e));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
