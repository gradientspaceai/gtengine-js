import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { Vector, add, dot, length, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistPointAlignedBox.h.
// ---------------------------------------------------------------------------

describe('DistPointAlignedBox verification', () => {
    const query = new DistPointAlignedBox();

    const boxArb = (n: number): fc.Arbitrary<AlignedBox> =>
        fc.tuple(wellScaledVector(n, -5, 5),
            fc.array(finite(0, 4), { minLength: n, maxLength: n }))
            .map(([c, e]) => {
                const lo = new Vector(n);
                const hi = new Vector(n);
                for (let i = 0; i < n; ++i) {
                    lo.set(i, c.get(i) - e[i]);
                    hi.set(i, c.get(i) + e[i]);
                }
                return AlignedBox.fromMinMax(lo, hi);
            });

    // The closest point of an aligned box is the componentwise clamp, which
    // is an independent closed form for the whole query.
    function clampToBox(p: Vector, b: AlignedBox): Vector {
        const q = new Vector(p.size);
        for (let i = 0; i < p.size; ++i) {
            q.set(i, Math.min(Math.max(p.get(i), b.min.get(i)),
                b.max.get(i)));
        }
        return q;
    }

    for (const n of [2, 3, 4]) {
        it(`matches the componentwise clamp in ${n}D`, () => {
            check(fc.tuple(wellScaledVector(n, -8, 8), boxArb(n)),
                ([p, b]) => {
                    const r = query.compute(p, b);
                    const q = clampToBox(p, b);
                    expectVectorClose(r.closest[1], q, 1e-12, 1e-12);
                    let sqr = 0;
                    for (let i = 0; i < n; ++i) {
                        sqr += (p.get(i) - q.get(i)) ** 2;
                    }
                    expectClose(r.sqrDistance, sqr, 1e-12, 1e-12);
                    expectClose(r.distance, Math.sqrt(sqr), 1e-12, 1e-12);
                    // closest[0] is a copy of the input point.
                    expectVectorClose(r.closest[0], p, 0, 0);
                    expect(r.closest[0]).not.toBe(p);
                });
        });
    }

    it('reports zero distance for points in the box', () => {
        check(fc.tuple(boxArb(3), fc.array(finite(0, 1),
            { minLength: 3, maxLength: 3 })), ([b, u]) => {
            const p = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                p.set(i, b.min.get(i) + u[i] * (b.max.get(i) - b.min.get(i)));
            }
            const r = query.compute(p, b);
            // Not exactly zero: the query works in the centered form, and
            // min - 0.5*(min+max) differs from -0.5*(max-min) by a rounding
            // error, so a point on the boundary clamps a half-ulp outside.
            expectClose(r.distance, 0, 1e-12, 1e-12);
            expectVectorClose(r.closest[1], p, 1e-12, 1e-12);
        });
    });

    it('is equivariant under translation', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), boxArb(3),
            wellScaledVector(3, -6, 6)), ([p, b, shift]) => {
            const moved = AlignedBox.fromMinMax(add(b.min, shift),
                add(b.max, shift));
            const r0 = query.compute(p, b);
            const r1 = query.compute(add(p, shift), moved);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(add(r0.closest[1], shift), r1.closest[1],
                1e-9, 1e-9);
        });
    });

    it('agrees with the canonical-box query on the centered form', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), boxArb(3)), ([p, b]) => {
            const { center, extent } = b.getCenteredForm();
            const r0 = query.compute(p, b);
            const r1 = new DistPointCanonicalBox().compute(sub(p, center),
                CanonicalBox.fromExtent(extent));
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], add(r1.closest[1], center),
                1e-12, 1e-12);
        });
    });

    it('handles a degenerate box as a point query', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), wellScaledVector(3, -5, 5)),
            ([p, c]) => {
                const r = query.compute(p, AlignedBox.fromMinMax(c, c));
                expectVectorClose(r.closest[1], c, 1e-12, 1e-12);
                expectClose(r.distance, length(sub(p, c)), 1e-12, 1e-12);
            });
    });
});
