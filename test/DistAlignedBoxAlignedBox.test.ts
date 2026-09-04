import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistAlignedBoxAlignedBox } from '../src/DistAlignedBoxAlignedBox.js';
import { Vector, add, dot, sub } from '../src/Vector.js';
import { check, expectClose, fc, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The reference distance: for aligned boxes the squared distance is the sum
// of the per-axis interval gaps.
function referenceSqrDistance(b0: AlignedBox, b1: AlignedBox): number {
    let sum = 0;
    for (let i = 0; i < b0.min.size; ++i) {
        const gap = Math.max(0,
            Math.max(b0.min.values[i] - b1.max.values[i],
                b1.min.values[i] - b0.max.values[i]));
        sum += gap * gap;
    }
    return sum;
}

describe('DistAlignedBoxAlignedBox', () => {
    const query = new DistAlignedBoxAlignedBox();

    it('reports zero distance for overlapping boxes', () => {
        const b0 = box([0, 0], [2, 2]);
        const b1 = box([1, 1], [3, 3]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        // The closest sets are the overlap box for both inputs.
        expect(result.closest[0].min.values).toEqual([1, 1]);
        expect(result.closest[0].max.values).toEqual([2, 2]);
        expect(result.closest[1].min.values).toEqual([1, 1]);
        expect(result.closest[1].max.values).toEqual([2, 2]);
    });

    it('reports zero distance for a contained box', () => {
        const b0 = box([0, 0, 0], [10, 10, 10]);
        const b1 = box([2, 3, 4], [5, 6, 7]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        expect(result.closest[0].min.values).toEqual([2, 3, 4]);
        expect(result.closest[0].max.values).toEqual([5, 6, 7]);
    });

    it('measures a one-axis gap', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([4, 0], [5, 1]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBeCloseTo(3, 12);
        // The closest sets degenerate in x and span the overlap in y.
        expect(result.closest[0].min.values).toEqual([1, 0]);
        expect(result.closest[0].max.values).toEqual([1, 1]);
        expect(result.closest[1].min.values).toEqual([4, 0]);
        expect(result.closest[1].max.values).toEqual([4, 1]);
    });

    it('measures a corner-to-corner gap', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([4, 5], [6, 7]);
        const result = query.compute(b0, b1);
        expect(result.sqrDistance).toBeCloseTo(9 + 16, 12);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[0].min.values).toEqual([1, 1]);
        expect(result.closest[0].max.values).toEqual([1, 1]);
        expect(result.closest[1].min.values).toEqual([4, 5]);
        expect(result.closest[1].max.values).toEqual([4, 5]);
    });

    it('reports zero distance for boxes that just touch', () => {
        const b0 = box([0, 0], [1, 1]);
        const b1 = box([1, 0], [2, 1]);
        const result = query.compute(b0, b1);
        expect(result.distance).toBe(0);
        expect(result.closest[0].min.values[0]).toBe(1);
        expect(result.closest[0].max.values[0]).toBe(1);
    });

    it('is symmetric', () => {
        const b0 = box([-1, 2, 0], [3, 4, 1]);
        const b1 = box([7, -5, 0.5], [9, -1, 2]);
        const a = query.compute(b0, b1);
        const b = query.compute(b1, b0);
        expect(b.distance).toBeCloseTo(a.distance, 12);
        expect(b.closest[1].min.values).toEqual(a.closest[0].min.values);
        expect(b.closest[0].min.values).toEqual(a.closest[1].min.values);
    });

    it('agrees with the per-axis gap reference and the closest sets realize '
        + 'the distance', () => {
        let seed = 2468;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 10 - 5;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const lo0 = [rand(), rand(), rand()];
            const lo1 = [rand(), rand(), rand()];
            const b0 = box(lo0, lo0.map(x => x + Math.abs(rand()) + 0.1));
            const b1 = box(lo1, lo1.map(x => x + Math.abs(rand()) + 0.1));
            const result = query.compute(b0, b1);
            expect(result.sqrDistance).toBeCloseTo(
                referenceSqrDistance(b0, b1), 9);

            // The centers of the closest sets realize the distance and lie
            // inside their respective input boxes.
            let sum = 0;
            for (let i = 0; i < 3; ++i) {
                const p0 = 0.5 * (result.closest[0].min.values[i]
                    + result.closest[0].max.values[i]);
                const p1 = 0.5 * (result.closest[1].min.values[i]
                    + result.closest[1].max.values[i]);
                expect(p0).toBeGreaterThanOrEqual(b0.min.values[i] - 1e-12);
                expect(p0).toBeLessThanOrEqual(b0.max.values[i] + 1e-12);
                expect(p1).toBeGreaterThanOrEqual(b1.min.values[i] - 1e-12);
                expect(p1).toBeLessThanOrEqual(b1.max.values[i] + 1e-12);
                sum += (p0 - p1) * (p0 - p1);
            }
            expect(sum).toBeCloseTo(result.sqrDistance, 9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistAlignedBoxAlignedBox.ts against the upstream header
// DistAlignedBoxAlignedBox.h.
// ---------------------------------------------------------------------------

const boxOfDim = (n: number): fc.Arbitrary<AlignedBox> =>
    fc.tuple(wellScaledVector(n, -8, 8), wellScaledVector(n, -8, 8))
        .map(([a, b]) => {
            const lo = new Vector(n);
            const hi = new Vector(n);
            for (let i = 0; i < n; ++i) {
                lo.values[i] = Math.min(a.values[i], b.values[i]);
                hi.values[i] = Math.max(a.values[i], b.values[i]);
            }
            return AlignedBox.fromMinMax(lo, hi);
        });

const boxPair = (n: number): fc.Arbitrary<[AlignedBox, AlignedBox]> =>
    fc.tuple(boxOfDim(n), boxOfDim(n));

describe('DistAlignedBoxAlignedBox verification', () => {
    const query = new DistAlignedBoxAlignedBox();

    it('matches the per-axis gap formula', () => {
        check(fc.integer({ min: 1, max: 5 }).chain(n => boxPair(n)),
            ([b0, b1]) => {
                const r = query.compute(b0, b1);
                let expected = 0;
                for (let i = 0; i < b0.min.size; ++i) {
                    const gap = Math.max(0,
                        b0.min.values[i] - b1.max.values[i],
                        b1.min.values[i] - b0.max.values[i]);
                    expected += gap * gap;
                }
                expectClose(r.sqrDistance, expected, 1e-12, 1e-12);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
            });
    });

    it('reports closest boxes contained in their inputs', () => {
        check(fc.integer({ min: 1, max: 4 }).chain(n => boxPair(n)),
            ([b0, b1]) => {
                const r = query.compute(b0, b1);
                const inputs = [b0, b1];
                for (let j = 0; j < 2; ++j) {
                    for (let i = 0; i < b0.min.size; ++i) {
                        expect(r.closest[j].min.values[i])
                            .toBeLessThanOrEqual(r.closest[j].max.values[i]);
                        expect(r.closest[j].min.values[i])
                            .toBeGreaterThanOrEqual(inputs[j].min.values[i]);
                        expect(r.closest[j].max.values[i])
                            .toBeLessThanOrEqual(inputs[j].max.values[i]);
                    }
                }
            });
    });

    // The upstream file comment claims that "any choice of P0 in closest[0]
    // and any choice of P1 in closest[1]" is a closest pair. That is an
    // overstatement: on an axis where the input boxes overlap, both closest
    // boxes carry the same nondegenerate interval, and only pairs that agree
    // on that coordinate realize the distance (picking opposite ends of the
    // overlap gives a strictly larger separation). The correct statement,
    // tested here, uses a common interpolation parameter per axis; upstream's
    // own suggestion of the two box centers is the special case t = 1/2.
    it('realizes the distance for matched points of the closest boxes', () => {
        const rand = seededRandom(0x51db);
        check(boxPair(3), ([b0, b1]) => {
            const r = query.compute(b0, b1);
            for (let k = 0; k < 12; ++k) {
                const p = new Vector(3);
                const q = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    const t = rand();
                    p.values[i] = (1 - t) * r.closest[0].min.values[i]
                        + t * r.closest[0].max.values[i];
                    q.values[i] = (1 - t) * r.closest[1].min.values[i]
                        + t * r.closest[1].max.values[i];
                }
                const diff = sub(p, q);
                expectClose(dot(diff, diff), r.sqrDistance, 1e-9, 1e-9);
            }
        }, 60);
    });

    it('is minimal over sampled box points', () => {
        const rand = seededRandom(0x51dc);
        check(boxPair(3), ([b0, b1]) => {
            const r = query.compute(b0, b1);
            const p = new Vector(3);
            const q = new Vector(3);
            for (let k = 0; k < 24; ++k) {
                for (let i = 0; i < 3; ++i) {
                    const t = rand();
                    const u = rand();
                    p.values[i] = (1 - t) * b0.min.values[i]
                        + t * b0.max.values[i];
                    q.values[i] = (1 - u) * b1.min.values[i]
                        + u * b1.max.values[i];
                }
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance)
                    .toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is symmetric under argument swap', () => {
        check(fc.integer({ min: 1, max: 4 }).chain(n => boxPair(n)),
            ([b0, b1]) => {
                const r0 = query.compute(b0, b1);
                const r1 = query.compute(b1, b0);
                expectClose(r0.sqrDistance, r1.sqrDistance, 0, 0);
                for (let i = 0; i < b0.min.size; ++i) {
                    expectClose(r0.closest[0].min.values[i] + 0,
                        r1.closest[1].min.values[i] + 0, 0, 0);
                    expectClose(r0.closest[0].max.values[i] + 0,
                        r1.closest[1].max.values[i] + 0, 0, 0);
                    expectClose(r0.closest[1].min.values[i] + 0,
                        r1.closest[0].min.values[i] + 0, 0, 0);
                    expectClose(r0.closest[1].max.values[i] + 0,
                        r1.closest[0].max.values[i] + 0, 0, 0);
                }
            });
    });

    it('is invariant under a common translation', () => {
        check(fc.tuple(boxPair(3), wellScaledVector(3, -5, 5)),
            ([pair, tr]) => {
                const move = (b: AlignedBox): AlignedBox =>
                    AlignedBox.fromMinMax(add(b.min, tr), add(b.max, tr));
                expectClose(query.compute(pair[0], pair[1]).sqrDistance,
                    query.compute(move(pair[0]), move(pair[1])).sqrDistance,
                    1e-9, 1e-9);
            });
    });

    it('reports the intersection box for overlapping boxes', () => {
        check(fc.integer({ min: 1, max: 4 }).chain(n =>
            fc.tuple(boxOfDim(n), wellScaledVector(n, -1, 1)))
            .map(([b, shift]) => {
                // Shift by less than the extent on each axis, so the boxes
                // always overlap on every axis.
                const n = b.min.size;
                const lo = new Vector(n);
                const hi = new Vector(n);
                for (let i = 0; i < n; ++i) {
                    const half = 0.4 * (b.max.values[i] - b.min.values[i]);
                    lo.values[i] = b.min.values[i] + half * shift.values[i];
                    hi.values[i] = b.max.values[i] + half * shift.values[i];
                }
                return [b, AlignedBox.fromMinMax(lo, hi)] as
                    [AlignedBox, AlignedBox];
            }),
        ([b0, b1]) => {
            const r = query.compute(b0, b1);
            expect(r.sqrDistance).toBe(0);
            expect(r.distance).toBe(0);
            for (let j = 0; j < 2; ++j) {
                for (let i = 0; i < b0.min.size; ++i) {
                    expectClose(r.closest[j].min.values[i],
                        Math.max(b0.min.values[i], b1.min.values[i]), 0, 0);
                    expectClose(r.closest[j].max.values[i],
                        Math.min(b0.max.values[i], b1.max.values[i]), 0, 0);
                }
            }
        });
    });
});
