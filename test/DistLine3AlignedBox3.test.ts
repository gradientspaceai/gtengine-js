import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistLine3AlignedBox3 } from '../src/DistLine3AlignedBox3.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3AlignedBox3', () => {
    const query = new DistLine3AlignedBox3();

    it('reports zero distance for a line through the box', () => {
        const box = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
        const result = query.compute(line([-5, 0, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to a box face', () => {
        const box = AlignedBox.fromMinMax(v(0, 0, 0), v(2, 2, 2));
        const result = query.compute(line([1, 1, 7], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(2, 12);
    });

    it('measures a line offset from a box edge', () => {
        const box = AlignedBox.fromMinMax(v(0, 0, 0), v(1, 1, 1));
        // A line parallel to the x-axis passing through (0, 4, 4) offsets.
        const result = query.compute(line([5, 4, 5], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('matches the canonical-box query for a centered box', () => {
        const rnd = makeRandom(97531);
        const cbQuery = new DistLine3CanonicalBox3();
        const extent = v(1, 2, 0.5);
        const cbox = CanonicalBox.fromExtent(extent);
        const abox = AlignedBox.fromMinMax(v(-1, -2, -0.5), v(1, 2, 0.5));
        for (let trial = 0; trial < 60; ++trial) {
            const ln = line([10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5],
                [2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1]);
            if (dot(ln.direction, ln.direction) < 1e-4) {
                continue;
            }
            const r0 = cbQuery.compute(ln, cbox);
            const r1 = query.compute(ln, abox);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });

    it('is translation invariant', () => {
        const rnd = makeRandom(1357);
        const shift = v(3, -4, 5);
        const box0 = AlignedBox.fromMinMax(v(-1, -2, -3), v(2, 1, 0));
        const box1 = AlignedBox.fromMinMax(add(v(-1, -2, -3), shift),
            add(v(2, 1, 0), shift));
        for (let trial = 0; trial < 40; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r0 = query.compute(Line.fromOriginDirection(origin, dir),
                box0);
            const r1 = query.compute(
                Line.fromOriginDirection(add(origin, shift), dir), box1);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });

    it('agrees with a dense sampling of the box', () => {
        const rnd = makeRandom(80808);
        const box = AlignedBox.fromMinMax(v(-1, -0.5, 0), v(2, 1.5, 1));

        for (let trial = 0; trial < 25; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, box);

            // The reported line point matches the reported parameter.
            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onLine.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
                // The reported box point is in the box.
                expect(result.closest[1].values[i]).toBeGreaterThanOrEqual(
                    box.min.values[i] - 1e-8);
                expect(result.closest[1].values[i]).toBeLessThanOrEqual(
                    box.max.values[i] + 1e-8);
            }
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled box point is closer to the line.
            const n = 16;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    for (let k = 0; k <= n; ++k) {
                        const q = v(
                            box.min.values[0] + (i / n)
                                * (box.max.values[0] - box.min.values[0]),
                            box.min.values[1] + (j / n)
                                * (box.max.values[1] - box.min.values[1]),
                            box.min.values[2] + (k / n)
                                * (box.max.values[2] - box.min.values[2]));
                        const w = sub(q, ln.origin);
                        const t = dot(w, ln.direction) / dd;
                        const f = sub(w, mul(t, ln.direction));
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-8);
        }
    });
});
