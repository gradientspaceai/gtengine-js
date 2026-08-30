import { describe, it, expect } from 'vitest';
import {
    IntrIntervalsTI,
    IntrIntervalsFI,
    IntrIntervalsFIResultType
} from '../src/IntrIntervals';

const ti = new IntrIntervalsTI();
const fi = new IntrIntervalsFI();

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Independent (brute force) static overlap of two finite intervals: the
// intersection of [a0,b0] and [a1,b1] is [max(a0,a1), min(b0,b1)] when that
// is nonempty.
function bruteForceOverlap(interval0: number[], interval1: number[]): number[] | null {
    const lo = Math.max(interval0[0], interval1[0]);
    const hi = Math.min(interval0[1], interval1[1]);
    return lo <= hi ? [lo, hi] : null;
}

// The interval translated by 'speed * t'.
function moved(interval: number[], speed: number, t: number): number[] {
    return [interval[0] + speed * t, interval[1] + speed * t];
}

describe('IntrIntervalsTI static query of two finite intervals', () => {
    it('reports overlapping intervals', () => {
        expect(ti.test([0, 2], [1, 3]).intersect).toBe(true);
    });

    it('reports intervals touching at a single endpoint', () => {
        expect(ti.test([0, 2], [2, 5]).intersect).toBe(true);
        expect(ti.test([2, 5], [0, 2]).intersect).toBe(true);
    });

    it('reports disjoint intervals', () => {
        expect(ti.test([0, 2], [3, 5]).intersect).toBe(false);
        expect(ti.test([3, 5], [0, 2]).intersect).toBe(false);
    });

    it('reports nested intervals', () => {
        expect(ti.test([0, 10], [3, 4]).intersect).toBe(true);
        expect(ti.test([3, 4], [0, 10]).intersect).toBe(true);
    });

    it('handles degenerate (point) intervals', () => {
        // The point is interior to the other interval.
        expect(ti.test([1, 1], [0, 2]).intersect).toBe(true);
        // The point is an endpoint of the other interval.
        expect(ti.test([2, 2], [0, 2]).intersect).toBe(true);
        expect(ti.test([0, 0], [0, 2]).intersect).toBe(true);
        // The point is outside the other interval.
        expect(ti.test([3, 3], [0, 2]).intersect).toBe(false);
        // Two identical points and two distinct points.
        expect(ti.test([7, 7], [7, 7]).intersect).toBe(true);
        expect(ti.test([7, 7], [8, 8]).intersect).toBe(false);
    });

    it('leaves the contact times at zero for the static query', () => {
        const result = ti.test([0, 2], [1, 3]);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBe(0);
    });

    it('is symmetric in its arguments', () => {
        const random = makeRandom(12345);
        for (let trial = 0; trial < 200; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const forward = ti.test([a0, b0], [a1, b1]).intersect;
            const reverse = ti.test([a1, b1], [a0, b0]).intersect;
            expect(forward).toBe(reverse);
            // Cross-check against the brute-force overlap.
            expect(forward).toBe(bruteForceOverlap([a0, b0], [a1, b1]) !== null);
        }
    });
});

describe('IntrIntervalsTI static query with one semiinfinite interval', () => {
    it('handles [a,+infinity)', () => {
        // finite = [0,2] versus [a,+infinity).
        expect(ti.testFiniteSemiInfinite([0, 2], 1, true).intersect).toBe(true);
        expect(ti.testFiniteSemiInfinite([0, 2], 2, true).intersect).toBe(true);
        expect(ti.testFiniteSemiInfinite([0, 2], 3, true).intersect).toBe(false);
        expect(ti.testFiniteSemiInfinite([0, 2], -5, true).intersect).toBe(true);
    });

    it('handles (-infinity,a]', () => {
        expect(ti.testFiniteSemiInfinite([0, 2], 1, false).intersect).toBe(true);
        expect(ti.testFiniteSemiInfinite([0, 2], 0, false).intersect).toBe(true);
        expect(ti.testFiniteSemiInfinite([0, 2], -1, false).intersect).toBe(false);
        expect(ti.testFiniteSemiInfinite([0, 2], 5, false).intersect).toBe(true);
    });

    it('agrees with a large-but-finite stand-in for the infinite endpoint', () => {
        const random = makeRandom(999);
        const big = 1000;
        for (let trial = 0; trial < 100; ++trial) {
            const f0 = Math.round(20 * random() - 10);
            const f1 = f0 + Math.round(5 * random());
            const a = Math.round(20 * random() - 10);
            expect(ti.testFiniteSemiInfinite([f0, f1], a, true).intersect)
                .toBe(ti.test([f0, f1], [a, big]).intersect);
            expect(ti.testFiniteSemiInfinite([f0, f1], a, false).intersect)
                .toBe(ti.test([f0, f1], [-big, a]).intersect);
        }
    });
});

describe('IntrIntervalsTI static query with two semiinfinite intervals', () => {
    it('two positive-infinite intervals always intersect', () => {
        expect(ti.testSemiInfiniteSemiInfinite(0, true, 100, true).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(100, true, 0, true).intersect).toBe(true);
    });

    it('two negative-infinite intervals always intersect', () => {
        expect(ti.testSemiInfiniteSemiInfinite(0, false, 100, false).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(100, false, 0, false).intersect).toBe(true);
    });

    it('positive-infinite versus negative-infinite intersects when a0 <= a1', () => {
        expect(ti.testSemiInfiniteSemiInfinite(1, true, 3, false).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(3, true, 3, false).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(4, true, 3, false).intersect).toBe(false);
    });

    it('negative-infinite versus positive-infinite intersects when a0 >= a1', () => {
        expect(ti.testSemiInfiniteSemiInfinite(3, false, 1, true).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(3, false, 3, true).intersect).toBe(true);
        expect(ti.testSemiInfiniteSemiInfinite(3, false, 4, true).intersect).toBe(false);
    });

    it('is symmetric under swapping the two semiinfinite intervals', () => {
        const values = [-2, 0, 1, 3];
        for (const a0 of values) {
            for (const a1 of values) {
                for (const p0 of [true, false]) {
                    for (const p1 of [true, false]) {
                        expect(ti.testSemiInfiniteSemiInfinite(a0, p0, a1, p1).intersect)
                            .toBe(ti.testSemiInfiniteSemiInfinite(a1, p1, a0, p0).intersect);
                    }
                }
            }
        }
    });
});

describe('IntrIntervalsTI dynamic query', () => {
    it('computes first and last contact when interval0 is on the left', () => {
        // interval0 = [0,1] moving right at speed 2, interval1 = [5,7] at
        // rest. The right endpoint of interval0 reaches 5 at t = 4/2 = 2.
        // The left endpoint of interval0 passes 7 at t = 7/2 = 3.5.
        const result = ti.testDynamic(10, [0, 1], 2, [5, 7], 0);
        expect(result.intersect).toBe(true);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
    });

    it('computes first and last contact when interval0 is on the right', () => {
        // Mirror image of the previous case.
        const result = ti.testDynamic(10, [5, 7], -2, [0, 1], 0);
        expect(result.intersect).toBe(true);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
    });

    it('reports contact times beyond maxTime with intersect false', () => {
        const result = ti.testDynamic(1, [0, 1], 2, [5, 7], 0);
        expect(result.intersect).toBe(false);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
    });

    it('reports intersect exactly at maxTime', () => {
        const result = ti.testDynamic(2, [0, 1], 2, [5, 7], 0);
        expect(result.intersect).toBe(true);
        expect(result.firstTime).toBeCloseTo(2, 12);
    });

    it('reports no intersection when the intervals separate', () => {
        const result = ti.testDynamic(10, [0, 1], -1, [5, 7], 1);
        expect(result.intersect).toBe(false);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBe(0);
    });

    it('reports no intersection when disjoint intervals have equal speeds', () => {
        const result = ti.testDynamic(10, [0, 1], 3, [5, 7], 3);
        expect(result.intersect).toBe(false);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBe(0);
    });

    it('handles intervals that initially intersect and later separate', () => {
        // interval0 = [0,4] at rest, interval1 = [2,6] moving right at 1.
        // The left endpoint of interval1 passes 4 at t = (4-2)/1 = 2.
        const faster1 = ti.testDynamic(10, [0, 4], 0, [2, 6], 1);
        expect(faster1.intersect).toBe(true);
        expect(faster1.firstTime).toBe(0);
        expect(faster1.lastTime).toBeCloseTo(2, 12);

        // interval0 = [0,4] moving right at 1, interval1 = [2,6] at rest.
        // The left endpoint of interval0 passes 6 at t = (6-0)/1 = 6.
        const faster0 = ti.testDynamic(10, [0, 4], 1, [2, 6], 0);
        expect(faster0.intersect).toBe(true);
        expect(faster0.firstTime).toBe(0);
        expect(faster0.lastTime).toBeCloseTo(6, 12);
    });

    it('reports an unbounded last contact time for equal speeds while overlapping', () => {
        const result = ti.testDynamic(10, [0, 4], 5, [2, 6], 5);
        expect(result.intersect).toBe(true);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBe(Number.MAX_VALUE);
    });

    it('is symmetric under swapping the two moving intervals', () => {
        const random = makeRandom(4242);
        for (let trial = 0; trial < 200; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const s0 = Math.round(10 * random() - 5);
            const s1 = Math.round(10 * random() - 5);
            const forward = ti.testDynamic(3, [a0, b0], s0, [a1, b1], s1);
            const reverse = ti.testDynamic(3, [a1, b1], s1, [a0, b0], s0);
            expect(forward.intersect).toBe(reverse.intersect);
            expect(forward.firstTime).toBeCloseTo(reverse.firstTime, 12);
            expect(forward.lastTime).toBeCloseTo(reverse.lastTime, 12);
        }
    });

    it('agrees with an independent simulation of the motion', () => {
        const random = makeRandom(777);
        const maxTime = 100;
        const tolerance = 1e-9;
        for (let trial = 0; trial < 200; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + 1 + Math.round(4 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + 1 + Math.round(4 * random());
            const s0 = Math.round(10 * random() - 5);
            const s1 = Math.round(10 * random() - 5);
            const interval0 = [a0, b0];
            const interval1 = [a1, b1];
            const result = ti.testDynamic(maxTime, interval0, s0, interval1, s1);

            // Separation of the moved intervals: negative when they overlap,
            // zero when they touch and positive when they are disjoint.
            const separation = (t: number): number => {
                const p0 = moved(interval0, s0, t);
                const p1 = moved(interval1, s1, t);
                return Math.max(p0[0] - p1[1], p1[0] - p0[1]);
            };

            if (!result.intersect && result.firstTime === 0) {
                // The intervals never meet, so sample the motion and check
                // that they stay disjoint.
                for (let i = 0; i <= 20; ++i) {
                    expect(separation((maxTime * i) / 20)).toBeGreaterThan(-tolerance);
                }
                continue;
            }

            const first = result.firstTime;
            if (first > 0) {
                // The intervals were initially disjoint, so they exactly
                // touch at the first contact time and are disjoint just
                // before it.
                expect(Math.abs(separation(first))).toBeLessThan(tolerance);
                expect(separation(first * (1 - 1e-3))).toBeGreaterThan(tolerance);
                // Halfway between first and last contact they overlap.
                const last = Math.min(result.lastTime, first + 1e6);
                expect(separation(0.5 * (first + last))).toBeLessThan(tolerance);
            }
            else {
                // The intervals intersect at time zero.
                expect(separation(0)).toBeLessThan(tolerance);
            }

            // The intervals touch at the last contact time and are disjoint
            // just after it (the unbounded case is skipped).
            if (result.lastTime < Number.MAX_VALUE) {
                const last = result.lastTime;
                expect(Math.abs(separation(last))).toBeLessThan(tolerance);
                expect(separation(last + Math.max(1e-3, last * 1e-3))).toBeGreaterThan(tolerance);
            }
        }
    });
});

describe('IntrIntervalsFI static query of two finite intervals', () => {
    it('finds the overlap of properly overlapping intervals', () => {
        const result = fi.find([0, 2], [1, 3]);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.overlap).toEqual([1, 2]);
        expect(result.type).toBe(IntrIntervalsFIResultType.isFinite);
    });

    it('finds the overlap of nested intervals', () => {
        const result = fi.find([0, 10], [3, 4]);
        expect(result.numIntersections).toBe(2);
        expect(result.overlap).toEqual([3, 4]);
        expect(result.type).toBe(IntrIntervalsFIResultType.isFinite);

        const swapped = fi.find([3, 4], [0, 10]);
        expect(swapped.overlap).toEqual([3, 4]);
        expect(swapped.type).toBe(IntrIntervalsFIResultType.isFinite);
    });

    it('finds a single point when the intervals touch', () => {
        const right = fi.find([0, 2], [2, 5]);
        expect(right.intersect).toBe(true);
        expect(right.numIntersections).toBe(1);
        expect(right.overlap).toEqual([2, 2]);
        expect(right.type).toBe(IntrIntervalsFIResultType.isPoint);

        const left = fi.find([2, 5], [0, 2]);
        expect(left.numIntersections).toBe(1);
        expect(left.overlap).toEqual([2, 2]);
        expect(left.type).toBe(IntrIntervalsFIResultType.isPoint);
    });

    it('reports an empty intersection for disjoint intervals', () => {
        const result = fi.find([0, 2], [3, 5]);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.overlap).toEqual([0, 0]);
        expect(result.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('handles degenerate (point) intervals', () => {
        const inside = fi.find([1, 1], [0, 2]);
        expect(inside.numIntersections).toBe(1);
        expect(inside.overlap).toEqual([1, 1]);
        expect(inside.type).toBe(IntrIntervalsFIResultType.isPoint);

        const atEnd = fi.find([2, 2], [0, 2]);
        expect(atEnd.numIntersections).toBe(1);
        expect(atEnd.overlap).toEqual([2, 2]);

        const same = fi.find([7, 7], [7, 7]);
        expect(same.numIntersections).toBe(1);
        expect(same.overlap).toEqual([7, 7]);

        const apart = fi.find([7, 7], [8, 8]);
        expect(apart.intersect).toBe(false);
        expect(apart.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('agrees with the brute-force overlap and is symmetric', () => {
        const random = makeRandom(2024);
        for (let trial = 0; trial < 300; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const result = fi.find([a0, b0], [a1, b1]);
            const swapped = fi.find([a1, b1], [a0, b0]);
            const expected = bruteForceOverlap([a0, b0], [a1, b1]);

            expect(result.intersect).toBe(expected !== null);
            expect(swapped.intersect).toBe(result.intersect);
            expect(swapped.numIntersections).toBe(result.numIntersections);
            expect(swapped.overlap).toEqual(result.overlap);
            expect(swapped.type).toBe(result.type);

            if (expected === null) {
                expect(result.numIntersections).toBe(0);
                expect(result.overlap).toEqual([0, 0]);
                expect(result.type).toBe(IntrIntervalsFIResultType.isEmpty);
            }
            else {
                expect(result.overlap).toEqual(expected);
                if (expected[0] === expected[1]) {
                    expect(result.numIntersections).toBe(1);
                    expect(result.type).toBe(IntrIntervalsFIResultType.isPoint);
                }
                else {
                    expect(result.numIntersections).toBe(2);
                    expect(result.type).toBe(IntrIntervalsFIResultType.isFinite);
                }
            }
            // The find and test queries agree on the Boolean result.
            expect(ti.test([a0, b0], [a1, b1]).intersect).toBe(result.intersect);
        }
    });
});

describe('IntrIntervalsFI static query with one semiinfinite interval', () => {
    it('intersects a finite interval with [a,+infinity)', () => {
        const overlap = fi.findFiniteSemiInfinite([0, 2], 1, true);
        expect(overlap.numIntersections).toBe(2);
        expect(overlap.overlap).toEqual([1, 2]);
        expect(overlap.type).toBe(IntrIntervalsFIResultType.isFinite);

        const contains = fi.findFiniteSemiInfinite([0, 2], -5, true);
        expect(contains.overlap).toEqual([0, 2]);
        expect(contains.type).toBe(IntrIntervalsFIResultType.isFinite);

        const touch = fi.findFiniteSemiInfinite([0, 2], 2, true);
        expect(touch.numIntersections).toBe(1);
        expect(touch.overlap).toEqual([2, 2]);
        expect(touch.type).toBe(IntrIntervalsFIResultType.isPoint);

        const empty = fi.findFiniteSemiInfinite([0, 2], 3, true);
        expect(empty.intersect).toBe(false);
        expect(empty.overlap).toEqual([0, 0]);
        expect(empty.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('intersects a finite interval with (-infinity,a]', () => {
        const overlap = fi.findFiniteSemiInfinite([0, 2], 1, false);
        expect(overlap.numIntersections).toBe(2);
        expect(overlap.overlap).toEqual([0, 1]);
        expect(overlap.type).toBe(IntrIntervalsFIResultType.isFinite);

        const contains = fi.findFiniteSemiInfinite([0, 2], 5, false);
        expect(contains.overlap).toEqual([0, 2]);

        const touch = fi.findFiniteSemiInfinite([0, 2], 0, false);
        expect(touch.numIntersections).toBe(1);
        expect(touch.overlap).toEqual([0, 0]);
        expect(touch.type).toBe(IntrIntervalsFIResultType.isPoint);

        const empty = fi.findFiniteSemiInfinite([0, 2], -1, false);
        expect(empty.intersect).toBe(false);
        expect(empty.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('handles a degenerate finite interval', () => {
        // The point 1 is in [1,+infinity) and in (-infinity,1], and the
        // degenerate overlap is reported as a point.
        const onBoundary = fi.findFiniteSemiInfinite([1, 1], 1, true);
        expect(onBoundary.numIntersections).toBe(1);
        expect(onBoundary.overlap).toEqual([1, 1]);
        expect(onBoundary.type).toBe(IntrIntervalsFIResultType.isPoint);

        const inside = fi.findFiniteSemiInfinite([1, 1], 0, true);
        expect(inside.numIntersections).toBe(1);
        expect(inside.overlap).toEqual([1, 1]);
        expect(inside.type).toBe(IntrIntervalsFIResultType.isPoint);

        const outside = fi.findFiniteSemiInfinite([1, 1], 2, true);
        expect(outside.intersect).toBe(false);
    });

    it('agrees with a large-but-finite stand-in for the infinite endpoint', () => {
        const random = makeRandom(31337);
        const big = 1000;
        for (let trial = 0; trial < 200; ++trial) {
            const f0 = Math.round(20 * random() - 10);
            const f1 = f0 + Math.round(5 * random());
            const a = Math.round(20 * random() - 10);

            const positive = fi.findFiniteSemiInfinite([f0, f1], a, true);
            const positiveRef = fi.find([f0, f1], [a, big]);
            expect(positive.intersect).toBe(positiveRef.intersect);
            expect(positive.numIntersections).toBe(positiveRef.numIntersections);
            expect(positive.overlap).toEqual(positiveRef.overlap);
            expect(positive.type).toBe(positiveRef.type);

            const negative = fi.findFiniteSemiInfinite([f0, f1], a, false);
            const negativeRef = fi.find([f0, f1], [-big, a]);
            expect(negative.intersect).toBe(negativeRef.intersect);
            expect(negative.numIntersections).toBe(negativeRef.numIntersections);
            expect(negative.overlap).toEqual(negativeRef.overlap);
            expect(negative.type).toBe(negativeRef.type);

            // The TI query agrees with the FI query.
            expect(ti.testFiniteSemiInfinite([f0, f1], a, true).intersect).toBe(positive.intersect);
            expect(ti.testFiniteSemiInfinite([f0, f1], a, false).intersect).toBe(negative.intersect);
        }
    });
});

describe('IntrIntervalsFI static query with two semiinfinite intervals', () => {
    it('intersects two positive-infinite intervals', () => {
        const result = fi.findSemiInfiniteSemiInfinite(1, true, 3, true);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        // overlap[1] is the +1 marker for the +infinity endpoint.
        expect(result.overlap).toEqual([3, 1]);
        expect(result.type).toBe(IntrIntervalsFIResultType.isPositiveInfinite);
    });

    it('intersects two negative-infinite intervals', () => {
        const result = fi.findSemiInfiniteSemiInfinite(1, false, 3, false);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        // overlap[0] is the -1 marker for the -infinity endpoint.
        expect(result.overlap).toEqual([-1, 1]);
        expect(result.type).toBe(IntrIntervalsFIResultType.isNegativeInfinite);
    });

    it('intersects [a0,+infinity) with (-infinity,a1]', () => {
        const finite = fi.findSemiInfiniteSemiInfinite(1, true, 3, false);
        expect(finite.numIntersections).toBe(2);
        expect(finite.overlap).toEqual([1, 3]);
        expect(finite.type).toBe(IntrIntervalsFIResultType.isFinite);

        const point = fi.findSemiInfiniteSemiInfinite(3, true, 3, false);
        expect(point.numIntersections).toBe(1);
        expect(point.overlap).toEqual([3, 3]);
        expect(point.type).toBe(IntrIntervalsFIResultType.isPoint);

        const empty = fi.findSemiInfiniteSemiInfinite(4, true, 3, false);
        expect(empty.intersect).toBe(false);
        expect(empty.overlap).toEqual([0, 0]);
        expect(empty.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('intersects (-infinity,a0] with [a1,+infinity)', () => {
        const finite = fi.findSemiInfiniteSemiInfinite(3, false, 1, true);
        expect(finite.numIntersections).toBe(2);
        expect(finite.overlap).toEqual([1, 3]);
        expect(finite.type).toBe(IntrIntervalsFIResultType.isFinite);

        const point = fi.findSemiInfiniteSemiInfinite(3, false, 3, true);
        expect(point.numIntersections).toBe(1);
        expect(point.overlap).toEqual([3, 3]);
        expect(point.type).toBe(IntrIntervalsFIResultType.isPoint);

        const empty = fi.findSemiInfiniteSemiInfinite(1, false, 3, true);
        expect(empty.intersect).toBe(false);
        expect(empty.overlap).toEqual([0, 0]);
        expect(empty.type).toBe(IntrIntervalsFIResultType.isEmpty);
    });

    it('is symmetric and agrees with the TI query', () => {
        const values = [-2, 0, 1, 3];
        for (const a0 of values) {
            for (const a1 of values) {
                for (const p0 of [true, false]) {
                    for (const p1 of [true, false]) {
                        const forward = fi.findSemiInfiniteSemiInfinite(a0, p0, a1, p1);
                        const reverse = fi.findSemiInfiniteSemiInfinite(a1, p1, a0, p0);
                        expect(forward.intersect).toBe(reverse.intersect);
                        expect(forward.numIntersections).toBe(reverse.numIntersections);
                        expect(forward.overlap).toEqual(reverse.overlap);
                        expect(forward.type).toBe(reverse.type);
                        expect(ti.testSemiInfiniteSemiInfinite(a0, p0, a1, p1).intersect)
                            .toBe(forward.intersect);
                    }
                }
            }
        }
    });
});

describe('IntrIntervalsFI dynamic query', () => {
    it('computes first and last contact when interval0 is on the left', () => {
        // interval0 = [0,1] moving right at speed 2, interval1 = [5,7] at
        // rest. First contact at t = (5-1)/2 = 2, where the right endpoint
        // of interval0 is at 1 + 2*2 = 5. Last contact at t = (7-0)/2 = 3.5.
        const result = fi.findDynamic(10, [0, 1], 2, [5, 7], 0);
        expect(result.intersect).toBe(true);
        expect(result.type).toBe(IntrIntervalsFIResultType.isDynamicQuery);
        expect(result.numIntersections).toBe(1);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
        expect(result.overlap[0]).toBeCloseTo(5, 12);
        expect(result.overlap[1]).toBeCloseTo(5, 12);
    });

    it('computes first and last contact when interval0 is on the right', () => {
        // Mirror image: first contact at t = (5-1)/2 = 2, where the left
        // endpoint of interval0 is at 5 - 2*2 = 1.
        const result = fi.findDynamic(10, [5, 7], -2, [0, 1], 0);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
        expect(result.overlap[0]).toBeCloseTo(1, 12);
        expect(result.overlap[1]).toBeCloseTo(1, 12);
    });

    it('reports contact times beyond maxTime with intersect false', () => {
        const result = fi.findDynamic(1, [0, 1], 2, [5, 7], 0);
        expect(result.intersect).toBe(false);
        expect(result.firstTime).toBeCloseTo(2, 12);
        expect(result.lastTime).toBeCloseTo(3.5, 12);
        expect(result.type).toBe(IntrIntervalsFIResultType.isDynamicQuery);
    });

    it('reports no intersection when the intervals separate', () => {
        const result = fi.findDynamic(10, [0, 1], -1, [5, 7], 1);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.overlap).toEqual([0, 0]);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBe(0);
        expect(result.type).toBe(IntrIntervalsFIResultType.isDynamicQuery);
    });

    it('reports the current overlap for initially intersecting intervals', () => {
        const result = fi.findDynamic(10, [0, 4], 0, [2, 6], 1);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.overlap).toEqual([2, 4]);
        expect(result.firstTime).toBe(0);
        expect(result.lastTime).toBeCloseTo(2, 12);

        const other = fi.findDynamic(10, [0, 4], 1, [2, 6], 0);
        expect(other.overlap).toEqual([2, 4]);
        expect(other.lastTime).toBeCloseTo(6, 12);
    });

    it('reports a degenerate overlap for initially touching intervals', () => {
        const right = fi.findDynamic(10, [0, 2], 0, [2, 5], 1);
        expect(right.numIntersections).toBe(1);
        expect(right.overlap).toEqual([2, 2]);
        expect(right.firstTime).toBe(0);
        expect(right.lastTime).toBe(0);

        const left = fi.findDynamic(10, [2, 5], 0, [0, 2], -1);
        expect(left.numIntersections).toBe(1);
        expect(left.overlap).toEqual([2, 2]);
        expect(left.firstTime).toBe(0);
    });

    it('reports an unbounded last contact time for equal speeds while overlapping', () => {
        const result = fi.findDynamic(10, [0, 4], 5, [2, 6], 5);
        expect(result.intersect).toBe(true);
        expect(result.lastTime).toBe(Number.MAX_VALUE);
        expect(result.overlap).toEqual([2, 4]);
    });

    it('reports the same intersection state as the TI dynamic query', () => {
        const random = makeRandom(8675309);
        for (let trial = 0; trial < 300; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const s0 = Math.round(10 * random() - 5);
            const s1 = Math.round(10 * random() - 5);
            const testResult = ti.testDynamic(4, [a0, b0], s0, [a1, b1], s1);
            const findResult = fi.findDynamic(4, [a0, b0], s0, [a1, b1], s1);
            expect(findResult.intersect).toBe(testResult.intersect);
            expect(findResult.firstTime).toBe(testResult.firstTime);
            expect(findResult.lastTime).toBe(testResult.lastTime);
            expect(findResult.type).toBe(IntrIntervalsFIResultType.isDynamicQuery);
        }
    });

    it('is symmetric under swapping the two moving intervals', () => {
        const random = makeRandom(555);
        for (let trial = 0; trial < 300; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const s0 = Math.round(10 * random() - 5);
            const s1 = Math.round(10 * random() - 5);
            const forward = fi.findDynamic(4, [a0, b0], s0, [a1, b1], s1);
            const reverse = fi.findDynamic(4, [a1, b1], s1, [a0, b0], s0);
            expect(forward.intersect).toBe(reverse.intersect);
            expect(forward.numIntersections).toBe(reverse.numIntersections);
            expect(forward.firstTime).toBeCloseTo(reverse.firstTime, 12);
            expect(forward.lastTime).toBeCloseTo(reverse.lastTime, 12);
            // The reported overlap is the same set, so the port fixes the
            // upstream asymmetry in the contact point of the separated case.
            expect(forward.overlap[0]).toBeCloseTo(reverse.overlap[0], 12);
            expect(forward.overlap[1]).toBeCloseTo(reverse.overlap[1], 12);
        }
    });

    it('reports a contact point that lies in both intervals at the first contact time', () => {
        const random = makeRandom(20260828);
        let separatedContacts = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const a0 = Math.round(20 * random() - 10);
            const b0 = a0 + Math.round(5 * random());
            const a1 = Math.round(20 * random() - 10);
            const b1 = a1 + Math.round(5 * random());
            const s0 = Math.round(10 * random() - 5);
            const s1 = Math.round(10 * random() - 5);
            const interval0 = [a0, b0];
            const interval1 = [a1, b1];
            const result = fi.findDynamic(100, interval0, s0, interval1, s1);
            if (result.numIntersections === 0) {
                continue;
            }

            // The reported overlap must be inside both intervals at the
            // first contact time. This is the independent cross-check that
            // detects the upstream contact-point error.
            const t = result.firstTime;
            const p0 = moved(interval0, s0, t);
            const p1 = moved(interval1, s1, t);
            const eps = 1e-9;
            for (const x of result.overlap) {
                expect(x).toBeGreaterThanOrEqual(p0[0] - eps);
                expect(x).toBeLessThanOrEqual(p0[1] + eps);
                expect(x).toBeGreaterThanOrEqual(p1[0] - eps);
                expect(x).toBeLessThanOrEqual(p1[1] + eps);
            }

            if (t > 0) {
                ++separatedContacts;
                // Separated intervals meet at a single point.
                expect(result.numIntersections).toBe(1);
                expect(result.overlap[0]).toBe(result.overlap[1]);
            }
        }
        // The random sample really exercised the separated-contact branch.
        expect(separatedContacts).toBeGreaterThan(20);
    });
});
