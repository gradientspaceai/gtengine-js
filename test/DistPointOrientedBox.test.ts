import { describe, expect, it } from 'vitest';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistPointOrientedBox.h.
// ---------------------------------------------------------------------------

describe('DistPointOrientedBox verification', () => {
    const query = new DistPointOrientedBox();
    const alignedQuery = new DistPointAlignedBox();

    const boxArb = (n: number): fc.Arbitrary<OrientedBox> =>
        fc.tuple(wellScaledVector(n, -5, 5), rotationFrame(n),
            fc.array(finite(0, 4), { minLength: n, maxLength: n }))
            .map(([c, axis, e]) => OrientedBox.fromCenterAxisExtent(c, axis,
                Vector.fromArray(e)));

    // Independent closed form: clamp the box-frame coordinates.
    function closestInBox(p: Vector, b: OrientedBox): Vector {
        const n = b.extent.size;
        const delta = sub(p, b.center);
        let q = b.center.clone();
        for (let i = 0; i < n; ++i) {
            const e = b.extent.get(i);
            const y = Math.min(Math.max(dot(b.axis[i], delta), -e), e);
            q = add(q, mul(y, b.axis[i]));
        }
        return q;
    }

    for (const n of [2, 3]) {
        it(`matches the box-frame clamp in ${n}D`, () => {
            check(fc.tuple(wellScaledVector(n, -8, 8), boxArb(n)),
                ([p, b]) => {
                    const r = query.compute(p, b);
                    const q = closestInBox(p, b);
                    expectVectorClose(r.closest[1], q, 1e-9, 1e-9);
                    expectClose(r.distance, length(sub(p, q)), 1e-9, 1e-9);
                    expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                        1e-12);
                    expectVectorClose(r.closest[0], p, 0, 0);
                    expect(r.closest[0]).not.toBe(p);
                });
        });
    }

    it('agrees with the aligned-box query for an identity frame', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        check(fc.tuple(wellScaledVector(3, -8, 8), wellScaledVector(3, -5, 5),
            fc.array(finite(0, 4), { minLength: 3, maxLength: 3 })),
            ([p, c, e]) => {
                const ext = Vector.fromArray(e);
                const ob = OrientedBox.fromCenterAxisExtent(c, axes, ext);
                const ab = AlignedBox.fromMinMax(sub(c, ext), add(c, ext));
                const r0 = query.compute(p, ob);
                const r1 = alignedQuery.compute(p, ab);
                expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-12, 1e-12);
            });
    });

    it('reports zero distance for points inside the box', () => {
        check(fc.tuple(boxArb(3), fc.array(finite(-1, 1),
            { minLength: 3, maxLength: 3 })), ([b, u]) => {
            let p = b.center.clone();
            for (let i = 0; i < 3; ++i) {
                p = add(p, mul(u[i] * b.extent.get(i), b.axis[i]));
            }
            const r = query.compute(p, b);
            expectClose(r.distance, 0, 1e-9, 1e-9);
            expectVectorClose(r.closest[1], p, 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), boxArb(3),
            rotationFrame(3), wellScaledVector(3, -6, 6)),
            ([p, b, frame, shift]) => {
                const rot = (q: Vector): Vector =>
                    add(add(mul(q.values[0], frame[0]),
                        mul(q.values[1], frame[1])),
                        mul(q.values[2], frame[2]));
                const moved = OrientedBox.fromCenterAxisExtent(
                    add(shift, rot(b.center)),
                    [rot(b.axis[0]), rot(b.axis[1]), rot(b.axis[2])],
                    b.extent);
                const r0 = query.compute(p, b);
                const r1 = query.compute(add(shift, rot(p)), moved);
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectVectorClose(add(shift, rot(r0.closest[1])),
                    r1.closest[1], 1e-7, 1e-7);
            });
    });

    it('is invariant to negating a box axis', () => {
        // Negating an axis describes the same solid box.
        check(fc.tuple(wellScaledVector(3, -8, 8), boxArb(3), fc.nat(2)),
            ([p, b, k]) => {
                const axes = [b.axis[0].clone(), b.axis[1].clone(),
                    b.axis[2].clone()];
                axes[k] = mul(-1, axes[k]);
                const flipped = OrientedBox.fromCenterAxisExtent(b.center,
                    axes, b.extent);
                const r0 = query.compute(p, b);
                const r1 = query.compute(p, flipped);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-9, 1e-9);
            });
    });
});
