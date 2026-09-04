import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2.js';
import { DistLine2OrientedBox2 } from '../src/DistLine2OrientedBox2.js';
import { Line } from '../src/Line.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine2OrientedBox2.h.
// ---------------------------------------------------------------------------

describe('DistLine2OrientedBox2 verification', () => {
    const query = new DistLine2OrientedBox2();
    const alignedQuery = new DistLine2AlignedBox2();

    // A well-scaled oriented box: rotation-frame axes and extents bounded
    // away from zero keep the box-frame coordinates non-subnormal.
    const boxArb = fc.tuple(wellScaledVector(2, -5, 5), rotationFrame(2),
        fc.tuple(finite(0, 4), finite(0, 4)))
        .map(([c, axis, e]) => OrientedBox.fromCenterAxisExtent(c, axis,
            v(e[0], e[1])));

    const lineArb = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Independent closed-form distance from a point to a solid oriented box.
    function pointBoxDistance(p: Vector, b: OrientedBox): number {
        const delta = sub(p, b.center);
        let sum = 0;
        for (let i = 0; i < 2; ++i) {
            const y = dot(b.axis[i], delta);
            const over = Math.abs(y) - b.extent.values[i];
            if (over > 0) { sum += over * over; }
        }
        return Math.sqrt(sum);
    }

    // Minimize the convex function f over [lo,hi] by ternary search. The
    // distance from a point moving along a line to a convex set is convex in
    // the line parameter, so this converges to the true minimum.
    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const r = query.compute(ln, b);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-9, 1e-9);
            expectVectorClose(r.closest[0],
                add(ln.origin, mul(r.parameter, ln.direction)), 1e-9, 1e-9);
            // closest[1] lies in the box.
            const delta = sub(r.closest[1], b.center);
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(dot(b.axis[i], delta)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
            }
        });
    });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const r = query.compute(ln, b);
            const best = ternaryMin(
                t => pointBoxDistance(add(ln.origin, mul(t, ln.direction)), b),
                -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });

    it('agrees with the aligned-box query for an identity frame', () => {
        const axes = [v(1, 0), v(0, 1)];
        check(fc.tuple(lineArb, wellScaledVector(2, -5, 5),
            fc.tuple(finite(0, 4), finite(0, 4))), ([ln, c, e]) => {
            const ext = v(e[0], e[1]);
            const ob = OrientedBox.fromCenterAxisExtent(c, axes, ext);
            const ab = AlignedBox.fromMinMax(sub(c, ext), add(c, ext));
            const r0 = query.compute(ln, ob);
            const r1 = alignedQuery.compute(ln, ab);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(r0.closest[0], r1.closest[0], 1e-9, 1e-9);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-9, 1e-9);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, boxArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([ln, b, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                v(frame[0].values[0] * p.values[0]
                    + frame[1].values[0] * p.values[1],
                    frame[0].values[1] * p.values[0]
                    + frame[1].values[1] * p.values[1]);
            const movedLine = Line.fromOriginDirection(
                add(shift, rot(ln.origin)), rot(ln.direction));
            const movedBox = OrientedBox.fromCenterAxisExtent(
                add(shift, rot(b.center)),
                [rot(b.axis[0]), rot(b.axis[1])], b.extent);
            const r0 = query.compute(ln, b);
            const r1 = query.compute(movedLine, movedBox);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(add(shift, rot(r0.closest[1])), r1.closest[1],
                1e-7, 1e-7);
        });
    });

    it('handles a degenerate box with zero extents as a point query', () => {
        check(fc.tuple(lineArb, wellScaledVector(2, -5, 5), rotationFrame(2)),
            ([ln, c, frame]) => {
                const b = OrientedBox.fromCenterAxisExtent(c, frame, v(0, 0));
                const r = query.compute(ln, b);
                const diff = sub(c, ln.origin);
                const t = dot(diff, ln.direction)
                    / dot(ln.direction, ln.direction);
                expectClose(r.distance,
                    length(sub(diff, mul(t, ln.direction))), 1e-9, 1e-9);
                expectVectorClose(r.closest[1], c, 1e-9, 1e-9);
            });
    });
});
