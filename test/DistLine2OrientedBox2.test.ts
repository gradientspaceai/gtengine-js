import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2.js';
import { DistLine2OrientedBox2 } from '../src/DistLine2OrientedBox2.js';
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

describe('DistLine2OrientedBox2', () => {
    const query = new DistLine2OrientedBox2();
    const axisAligned = [v(1, 0), v(0, 1)];

    it('reports zero distance for a line through the box', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0), axisAligned,
            v(1, 1));
        const result = query.compute(line([-5, 0], [1, 0]), box);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.sqrDistance).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to a box face', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0), axisAligned,
            v(1, 1));
        const result = query.compute(line([0, 4], [1, 0]), box);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
    });

    it('measures a line against a 45-degree rotated square', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0),
            [v(c, c), v(-c, c)], v(1, 1));
        // The horizontal line y = 3; the rotated square reaches y = sqrt(2).
        const result = query.compute(line([0, 3], [1, 0]), box);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
    });

    it('handles a degenerate box with zero extents', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 5),
            [v(c, c), v(-c, c)], v(0, 0));
        const result = query.compute(line([0, 0], [1, 0]), box);
        expect(result.distance).toBeCloseTo(5, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(5, 10);
    });

    it('matches the aligned-box query for an axis-aligned oriented box',
        () => {
            const rnd = makeRandom(24680);
            const abQuery = new DistLine2AlignedBox2();
            const abox = AlignedBox.fromMinMax(v(-1, -2), v(2, 1));
            const { center, extent } = abox.getCenteredForm();
            const obox = OrientedBox.fromCenterAxisExtent(center, axisAligned,
                extent);
            for (let trial = 0; trial < 60; ++trial) {
                const ln = line([10 * rnd() - 5, 10 * rnd() - 5],
                    [2 * rnd() - 1, 2 * rnd() - 1]);
                if (dot(ln.direction, ln.direction) < 1e-6) {
                    continue;
                }
                const r0 = abQuery.compute(ln, abox);
                const r1 = query.compute(ln, obox);
                expect(r1.distance).toBeCloseTo(r0.distance, 9);
            }
        });

    it('agrees with a dense sampling of line and box points', () => {
        const rnd = makeRandom(1122);
        const c = Math.cos(0.6), s = Math.sin(0.6);
        const axis = [v(c, s), v(-s, c)];
        const extent = v(1.5, 0.75);
        const box = OrientedBox.fromCenterAxisExtent(v(0.25, -0.5), axis,
            extent);

        for (let trial = 0; trial < 40; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, box);

            // The reported closest points realize the reported distance and
            // lie on their primitives.
            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            expect(onLine.values[0]).toBeCloseTo(result.closest[0].values[0],
                8);
            expect(onLine.values[1]).toBeCloseTo(result.closest[0].values[1],
                8);
            const delta = sub(result.closest[1], box.center);
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(dot(delta, axis[i]))).toBeLessThanOrEqual(
                    extent.values[i] + 1e-8);
            }
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 8);

            // No sampled box point is farther than the reported distance from
            // the line.
            const n = 60;
            let best = Number.MAX_VALUE;
            const dd = dot(ln.direction, ln.direction);
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    const s0 = (2 * i / n - 1) * extent.values[0];
                    const s1 = (2 * j / n - 1) * extent.values[1];
                    const q = add(box.center,
                        add(mul(s0, axis[0]), mul(s1, axis[1])));
                    // The squared distance from q to the line.
                    const w = sub(q, ln.origin);
                    const t = dot(w, ln.direction) / dd;
                    const f = sub(w, mul(t, ln.direction));
                    best = Math.min(best, dot(f, f));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-8);
        }
    });
});
