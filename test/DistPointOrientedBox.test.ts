import { describe, expect, it } from 'vitest';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox';
import { AlignedBox } from '../src/AlignedBox';
import { OrientedBox } from '../src/OrientedBox';
import { Vector, add, dot, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPointOrientedBox', () => {
    const query = new DistPointOrientedBox();

    it('reports zero distance for a point inside the box', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 2, 3));
        const result = query.compute(v(0.5, -1, 2), box);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values).toEqual([0.5, -1, 2]);
    });

    it('measures a point outside an axis-aligned oriented box', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = query.compute(v(4, 0, 0), box);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('measures a point against a 45-degree rotated 2D box', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0),
            [v(c, c), v(-c, c)], v(1, 1));
        // The rotated square has vertices at distance sqrt(2) along the axes.
        const result = query.compute(v(3, 0), box);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
    });

    it('handles a degenerate box with zero extents', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(1, 2),
            [v(c, c), v(-c, c)], v(0, 0));
        const result = query.compute(v(4, 6), box);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(2, 12);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('matches the aligned-box query for an axis-aligned oriented box', () => {
        const rnd = makeRandom(7331);
        const abQuery = new DistPointAlignedBox();
        const abox = AlignedBox.fromMinMax(v(-1, -2, -3), v(2, 1, 4));
        const { center, extent } = abox.getCenteredForm();
        const obox = OrientedBox.fromCenterAxisExtent(center,
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], extent);
        for (let trial = 0; trial < 50; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = abQuery.compute(p, abox);
            const r1 = query.compute(p, obox);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i]).toBeCloseTo(
                    r0.closest[1].values[i], 10);
            }
        }
    });

    it('is invariant to rotating the point and the box together', () => {
        const rnd = makeRandom(4242);
        // A rotation by angle a about the z-axis.
        const a = 0.7;
        const ca = Math.cos(a), sa = Math.sin(a);
        const rot = (p: Vector): Vector => v(
            ca * p.values[0] - sa * p.values[1],
            sa * p.values[0] + ca * p.values[1],
            p.values[2]);

        const axis = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        const box0 = OrientedBox.fromCenterAxisExtent(v(1, 2, 3), axis,
            v(1, 2, 0.5));
        const box1 = OrientedBox.fromCenterAxisExtent(rot(v(1, 2, 3)),
            [rot(axis[0]), rot(axis[1]), rot(axis[2])], v(1, 2, 0.5));
        for (let trial = 0; trial < 30; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = query.compute(p, box0);
            const r1 = query.compute(rot(p), box1);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            const expected = rot(r0.closest[1]);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i]).toBeCloseTo(
                    expected.values[i], 10);
            }
        }
    });

    it('agrees with a dense brute-force sampling of the box', () => {
        const rnd = makeRandom(999);
        const c = Math.SQRT1_2;
        const axis = [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)];
        const extent = v(1.5, 0.5, 2);
        const box = OrientedBox.fromCenterAxisExtent(v(0.5, -1, 0.25), axis,
            extent);

        for (let trial = 0; trial < 30; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, box);

            // The closest point is on the box.
            const delta = sub(result.closest[1], box.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(delta, axis[i]))).toBeLessThanOrEqual(
                    extent.values[i] + 1e-9);
            }

            // The closest point realizes the reported distance.
            const d = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(result.distance, 10);

            // No sampled box point is closer.
            const n = 10;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    for (let k = 0; k <= n; ++k) {
                        const s0 = (2 * i / n - 1) * extent.values[0];
                        const s1 = (2 * j / n - 1) * extent.values[1];
                        const s2 = (2 * k / n - 1) * extent.values[2];
                        const q = add(box.center, add(mul(s0, axis[0]),
                            add(mul(s1, axis[1]), mul(s2, axis[2]))));
                        const e = sub(p, q);
                        best = Math.min(best, dot(e, e));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});
