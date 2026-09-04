import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistLine3AlignedBox3 } from '../src/DistLine3AlignedBox3.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, unitVector,
    wellScaledVector
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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine3AlignedBox3.h.
// ---------------------------------------------------------------------------

describe('DistLine3AlignedBox3 verification', () => {
    const query = new DistLine3AlignedBox3();
    const canonicalQuery = new DistLine3CanonicalBox3();

    const boxArb = fc.tuple(wellScaledVector(3, -5, 5),
        fc.array(finite(0, 4), { minLength: 3, maxLength: 3 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            v(c.values[0] - e[0], c.values[1] - e[1], c.values[2] - e[2]),
            v(c.values[0] + e[0], c.values[1] + e[1], c.values[2] + e[2])));

    const lineArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Independent closed-form distance from a point to a solid aligned box.
    function pointBoxDistance(p: Vector, b: AlignedBox): number {
        let sum = 0;
        for (let i = 0; i < 3; ++i) {
            const under = b.min.values[i] - p.values[i];
            const over = p.values[i] - b.max.values[i];
            const d = Math.max(under, over, 0);
            sum += d * d;
        }
        return Math.sqrt(sum);
    }

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
            for (let i = 0; i < 3; ++i) {
                expect(r.closest[1].values[i])
                    .toBeGreaterThanOrEqual(b.min.values[i] - 1e-9);
                expect(r.closest[1].values[i])
                    .toBeLessThanOrEqual(b.max.values[i] + 1e-9);
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

    it('equals the canonical-box query on the translated line', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const { center, extent } = b.getCenteredForm();
            const cbox = CanonicalBox.fromExtent(extent);
            const shifted = Line.fromOriginDirection(sub(ln.origin, center),
                ln.direction);
            const r0 = query.compute(ln, b);
            const r1 = canonicalQuery.compute(shifted, cbox);
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectClose(r0.parameter, r1.parameter, 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], add(r1.closest[1], center),
                1e-9, 1e-9);
        });
    });

    it('is equivariant under translation', () => {
        check(fc.tuple(lineArb, boxArb, wellScaledVector(3, -6, 6)),
            ([ln, b, shift]) => {
                const movedLine = Line.fromOriginDirection(
                    add(ln.origin, shift), ln.direction);
                const movedBox = AlignedBox.fromMinMax(add(b.min, shift),
                    add(b.max, shift));
                const r0 = query.compute(ln, b);
                const r1 = query.compute(movedLine, movedBox);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(add(r0.closest[1], shift), r1.closest[1],
                    1e-7, 1e-7);
            });
    });

    it('reduces to the point-line distance for a degenerate box', () => {
        check(fc.tuple(lineArb, wellScaledVector(3, -5, 5)), ([ln, c]) => {
            const b = AlignedBox.fromMinMax(c, c);
            const r = query.compute(ln, b);
            const diff = sub(c, ln.origin);
            const t = dot(diff, ln.direction) / dot(ln.direction, ln.direction);
            expectClose(r.distance, length(sub(diff, mul(t, ln.direction))),
                1e-9, 1e-9);
            expectVectorClose(r.closest[1], c, 1e-12, 1e-12);
        });
    });
});
