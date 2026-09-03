import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine3AlignedBox3 } from '../src/DistLine3AlignedBox3.js';
import { DistLine3OrientedBox3 } from '../src/DistLine3OrientedBox3.js';
import { Line } from '../src/Line.js';
import { OrientedBox } from '../src/OrientedBox.js';
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

describe('DistLine3OrientedBox3', () => {
    const query = new DistLine3OrientedBox3();
    const axisAligned = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];

    it('reports zero distance for a line through the box', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0), axisAligned,
            v(1, 1, 1));
        const result = query.compute(line([-5, 0, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to a box face', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0), axisAligned,
            v(1, 1, 1));
        const result = query.compute(line([0, 0, 6], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('measures a line against a box rotated 45 degrees about z', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(1, 1, 1));
        // A line parallel to the z-axis at (3,0).
        const result = query.compute(line([3, 0, 0], [0, 0, 1]), box);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 9);
    });

    it('handles a degenerate box with zero extents', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 5),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(0, 0, 0));
        const result = query.compute(line([0, 0, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 9);
    });

    it('matches the aligned-box query for an axis-aligned oriented box',
        () => {
            const rnd = makeRandom(60606);
            const abQuery = new DistLine3AlignedBox3();
            const abox = AlignedBox.fromMinMax(v(-1, -2, -3), v(2, 1, 0));
            const { center, extent } = abox.getCenteredForm();
            const obox = OrientedBox.fromCenterAxisExtent(center, axisAligned,
                extent);
            for (let trial = 0; trial < 60; ++trial) {
                const ln = line(
                    [10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5],
                    [2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1]);
                if (dot(ln.direction, ln.direction) < 1e-4) {
                    continue;
                }
                const r0 = abQuery.compute(ln, abox);
                const r1 = query.compute(ln, obox);
                expect(r1.distance).toBeCloseTo(r0.distance, 9);
            }
        });

    it('agrees with a dense sampling of the box', () => {
        const rnd = makeRandom(70707);
        const c = Math.cos(0.5), s = Math.sin(0.5);
        const axis = [v(c, s, 0), v(-s, c, 0), v(0, 0, 1)];
        const extent = v(1.25, 0.75, 2);
        const box = OrientedBox.fromCenterAxisExtent(v(0.5, -1, 0.25), axis,
            extent);

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
            }

            // The reported box point is in the box.
            const delta = sub(result.closest[1], box.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(delta, axis[i]))).toBeLessThanOrEqual(
                    extent.values[i] + 1e-8);
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
                        const q = add(box.center, add(
                            mul((2 * i / n - 1) * extent.values[0], axis[0]),
                            add(
                                mul((2 * j / n - 1) * extent.values[1],
                                    axis[1]),
                                mul((2 * k / n - 1) * extent.values[2],
                                    axis[2]))));
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
