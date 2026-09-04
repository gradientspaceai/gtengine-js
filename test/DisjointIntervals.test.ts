import { describe, it, expect } from 'vitest';
import { DisjointIntervals } from '../src/DisjointIntervals.js';
import { check, fc, finite, scaled } from './helpers/arbitraries.js';

// Collect the intervals of a set as [xmin, xmax] pairs for comparison.
function intervalsOf(set: DisjointIntervals): [number, number][] {
    const result: [number, number][] = [];
    for (let i = 0; i < set.getNumIntervals(); ++i) {
        const interval = set.getInterval(i)!;
        result.push([interval.xmin, interval.xmax]);
    }
    return result;
}

// Build a set from a list of already-disjoint, sorted [xmin, xmax) pairs.
function makeSet(intervals: [number, number][]): DisjointIntervals {
    const set = new DisjointIntervals();
    for (const [xmin, xmax] of intervals) {
        set.insert(xmin, xmax);
    }
    return set;
}

// Independent membership predicate: x is in the set when some interval
// [xmin, xmax) contains it.
function contains(set: DisjointIntervals, x: number): boolean {
    for (const [xmin, xmax] of intervalsOf(set)) {
        if (xmin <= x && x < xmax) {
            return true;
        }
    }
    return false;
}

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// A random disjoint-interval set with integer endpoints in [0, range).
function randomSet(rand: () => number, range: number): DisjointIntervals {
    const endpoints = new Set<number>();
    const count = 2 * (1 + Math.floor(rand() * 4));
    while (endpoints.size < count) {
        endpoints.add(Math.floor(rand() * range));
    }
    const sorted = [...endpoints].sort((a, b) => a - b);
    const set = new DisjointIntervals();
    for (let i = 0; i < sorted.length; i += 2) {
        set.insert(sorted[i], sorted[i + 1]);
    }
    return set;
}

describe('DisjointIntervals', () => {
    it('default-constructs empty', () => {
        const set = new DisjointIntervals();
        expect(set.getNumIntervals()).toBe(0);
        expect(set.getInterval(0)).toBeNull();
    });

    it('constructs a single interval when xmin < xmax', () => {
        const set = new DisjointIntervals(1, 3);
        expect(set.getNumIntervals()).toBe(1);
        expect(set.getInterval(0)).toEqual({ xmin: 1, xmax: 3 });
    });

    it('constructs empty when xmin >= xmax', () => {
        expect(new DisjointIntervals(3, 1).getNumIntervals()).toBe(0);
        expect(new DisjointIntervals(2, 2).getNumIntervals()).toBe(0);
    });

    it('getInterval returns null out of range', () => {
        const set = new DisjointIntervals(0, 1);
        expect(set.getInterval(-1)).toBeNull();
        expect(set.getInterval(1)).toBeNull();
        expect(set.getInterval(0)).toEqual({ xmin: 0, xmax: 1 });
    });

    it('clear empties the set', () => {
        const set = new DisjointIntervals(0, 1);
        set.clear();
        expect(set.getNumIntervals()).toBe(0);
    });

    it('clone makes an independent copy', () => {
        const set = new DisjointIntervals(0, 1);
        const copy = set.clone();
        set.insert(5, 6);
        expect(intervalsOf(copy)).toEqual([[0, 1]]);
        expect(intervalsOf(set)).toEqual([[0, 1], [5, 6]]);
    });

    describe('insert', () => {
        it('rejects xmin >= xmax', () => {
            const set = new DisjointIntervals();
            expect(set.insert(2, 2)).toBe(false);
            expect(set.insert(3, 1)).toBe(false);
            expect(set.getNumIntervals()).toBe(0);
        });

        it('keeps disjoint intervals separate and sorted', () => {
            const set = makeSet([[5, 6], [1, 2]]);
            expect(intervalsOf(set)).toEqual([[1, 2], [5, 6]]);
        });

        it('merges overlapping intervals', () => {
            const set = makeSet([[1, 3], [2, 5]]);
            expect(intervalsOf(set)).toEqual([[1, 5]]);
        });

        it('merges abutting half-open intervals: [1,3) union [3,5) = [1,5)', () => {
            const set = makeSet([[1, 3], [3, 5]]);
            expect(intervalsOf(set)).toEqual([[1, 5]]);
        });

        it('merges an interval spanning several existing intervals', () => {
            const set = makeSet([[0, 1], [2, 3], [4, 5], [8, 9]]);
            set.insert(0.5, 4.5);
            expect(intervalsOf(set)).toEqual([[0, 5], [8, 9]]);
        });
    });

    describe('remove', () => {
        it('rejects xmin >= xmax', () => {
            const set = new DisjointIntervals(0, 10);
            expect(set.remove(4, 4)).toBe(false);
            expect(intervalsOf(set)).toEqual([[0, 10]]);
        });

        it('splits an interval', () => {
            const set = new DisjointIntervals(0, 10);
            expect(set.remove(3, 5)).toBe(true);
            expect(intervalsOf(set)).toEqual([[0, 3], [5, 10]]);
        });

        it('trims interval ends and removes whole intervals', () => {
            const set = makeSet([[0, 2], [3, 5], [6, 8]]);
            set.remove(1, 7);
            expect(intervalsOf(set)).toEqual([[0, 1], [7, 8]]);
        });

        it('removing a superset empties the set', () => {
            const set = makeSet([[1, 2], [3, 4]]);
            set.remove(0, 10);
            expect(set.getNumIntervals()).toBe(0);
        });
    });

    describe('static set operations (hand-computed)', () => {
        const A = makeSet([[0, 4], [6, 9]]);
        const B = makeSet([[2, 7], [8, 11]]);

        it('union', () => {
            expect(intervalsOf(DisjointIntervals.union(A, B)))
                .toEqual([[0, 11]]);
        });

        it('intersection', () => {
            expect(intervalsOf(DisjointIntervals.intersection(A, B)))
                .toEqual([[2, 4], [6, 7], [8, 9]]);
        });

        it('difference', () => {
            expect(intervalsOf(DisjointIntervals.difference(A, B)))
                .toEqual([[0, 2], [7, 8]]);
            expect(intervalsOf(DisjointIntervals.difference(B, A)))
                .toEqual([[4, 6], [9, 11]]);
        });

        it('exclusiveOr', () => {
            expect(intervalsOf(DisjointIntervals.exclusiveOr(A, B)))
                .toEqual([[0, 2], [4, 6], [7, 8], [9, 11]]);
        });

        it('operations with an empty set', () => {
            const empty = new DisjointIntervals();
            expect(intervalsOf(DisjointIntervals.union(A, empty))).toEqual(intervalsOf(A));
            expect(intervalsOf(DisjointIntervals.intersection(A, empty))).toEqual([]);
            expect(intervalsOf(DisjointIntervals.difference(A, empty))).toEqual(intervalsOf(A));
            expect(intervalsOf(DisjointIntervals.difference(empty, A))).toEqual([]);
            expect(intervalsOf(DisjointIntervals.exclusiveOr(A, empty))).toEqual(intervalsOf(A));
        });

        it('shared endpoints collapse: [0,2) xor [2,4) has endpoint 2 removed', () => {
            const left = new DisjointIntervals(0, 2);
            const right = new DisjointIntervals(2, 4);
            expect(intervalsOf(DisjointIntervals.exclusiveOr(left, right)))
                .toEqual([[0, 4]]);
            expect(intervalsOf(DisjointIntervals.intersection(left, right)))
                .toEqual([]);
        });
    });

    describe('set-algebra identities (randomized cross-check)', () => {
        it('operations agree with pointwise membership over random sets', () => {
            const rand = makeRandom(0x2b04);
            const range = 24;
            for (let trial = 0; trial < 50; ++trial) {
                const A = randomSet(rand, range);
                const B = randomSet(rand, range);
                const union = DisjointIntervals.union(A, B);
                const inter = DisjointIntervals.intersection(A, B);
                const diff = DisjointIntervals.difference(A, B);
                const xor = DisjointIntervals.exclusiveOr(A, B);
                // Sample at half-integers to probe every unit cell strictly
                // inside [0, range), avoiding endpoint ambiguity.
                for (let k = 0; k < range; ++k) {
                    const x = k + 0.5;
                    const inA = contains(A, x);
                    const inB = contains(B, x);
                    expect(contains(union, x)).toBe(inA || inB);
                    expect(contains(inter, x)).toBe(inA && inB);
                    expect(contains(diff, x)).toBe(inA && !inB);
                    expect(contains(xor, x)).toBe(inA !== inB);
                }
            }
        });

        it('structural identities hold on random sets', () => {
            const rand = makeRandom(0xb04b04);
            for (let trial = 0; trial < 50; ++trial) {
                const A = randomSet(rand, 24);
                const B = randomSet(rand, 24);
                // Commutativity.
                expect(intervalsOf(DisjointIntervals.union(A, B)))
                    .toEqual(intervalsOf(DisjointIntervals.union(B, A)));
                expect(intervalsOf(DisjointIntervals.intersection(A, B)))
                    .toEqual(intervalsOf(DisjointIntervals.intersection(B, A)));
                // A xor B = (A - B) union (B - A).
                expect(intervalsOf(DisjointIntervals.exclusiveOr(A, B)))
                    .toEqual(intervalsOf(DisjointIntervals.union(
                        DisjointIntervals.difference(A, B),
                        DisjointIntervals.difference(B, A))));
                // (A - B) union (A intersect B) = A.
                expect(intervalsOf(DisjointIntervals.union(
                    DisjointIntervals.difference(A, B),
                    DisjointIntervals.intersection(A, B))))
                    .toEqual(intervalsOf(A));
                // (A - B) intersect B is empty.
                expect(DisjointIntervals.intersection(
                    DisjointIntervals.difference(A, B), B).getNumIntervals())
                    .toBe(0);
                // Idempotence: A union A = A, A intersect A = A, A - A = {}.
                expect(intervalsOf(DisjointIntervals.union(A, A))).toEqual(intervalsOf(A));
                expect(intervalsOf(DisjointIntervals.intersection(A, A))).toEqual(intervalsOf(A));
                expect(DisjointIntervals.difference(A, A).getNumIntervals()).toBe(0);
                expect(DisjointIntervals.exclusiveOr(A, A).getNumIntervals()).toBe(0);
            }
        });
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). The Boolean operators are pure
// comparison-and-copy algorithms on the endpoint arrays: no arithmetic is
// performed on the scalars, so every property below holds *exactly* (no
// tolerances), including at the interval endpoints themselves.
// ---------------------------------------------------------------------------

// A set of half-open intervals built by inserting random [x, x+len) pairs. The
// endpoints are drawn from a coarse grid so that overlaps, abutments and
// shared endpoints occur often; a second generator uses raw doubles.
const gridSet = (n: number, cells: number): fc.Arbitrary<DisjointIntervals> =>
    fc.array(fc.tuple(fc.integer({ min: 0, max: cells - 1 }),
        fc.integer({ min: 1, max: 4 })), { minLength: 0, maxLength: n })
        .map(pairs => {
            const set = new DisjointIntervals();
            for (const [x, len] of pairs) { set.insert(x, x + len); }
            return set;
        });

const doubleSet = (n: number): fc.Arbitrary<DisjointIntervals> =>
    fc.array(fc.tuple(scaled(-8, 8, 32), scaled(0.5, 4, 8)),
        { minLength: 0, maxLength: n })
        .map(pairs => {
            const set = new DisjointIntervals();
            for (const [x, len] of pairs) { set.insert(x, x + len); }
            return set;
        });

// Every endpoint of the inputs plus the midpoint of each consecutive pair and
// one probe outside the hull: a finite sample set that distinguishes any two
// different unions of half-open intervals with these endpoints.
function probePoints(...sets: DisjointIntervals[]): number[] {
    const ends = new Set<number>();
    for (const set of sets) {
        for (const [xmin, xmax] of intervalsOf(set)) { ends.add(xmin); ends.add(xmax); }
    }
    const sorted = [...ends].sort((a, b) => a - b);
    const probes: number[] = [];
    for (let i = 0; i < sorted.length; ++i) {
        probes.push(sorted[i]);
        if (i + 1 < sorted.length) { probes.push(0.5 * (sorted[i] + sorted[i + 1])); }
    }
    probes.push((sorted[0] ?? 0) - 1, (sorted[sorted.length - 1] ?? 0) + 1);
    return probes;
}

// A well-formed set has strictly increasing endpoints (sorted, disjoint and
// non-abutting), i.e. the representation is canonical.
function expectWellFormed(set: DisjointIntervals): void {
    const intervals = intervalsOf(set);
    for (let i = 0; i < intervals.length; ++i) {
        expect(intervals[i][0]).toBeLessThan(intervals[i][1]);
        if (i > 0) { expect(intervals[i - 1][1]).toBeLessThan(intervals[i][0]); }
    }
}

describe('DisjointIntervals verification', () => {
    it('the four operators agree with pointwise membership (grid endpoints)', () => {
        check(fc.tuple(gridSet(5, 12), gridSet(5, 12)), ([A, B]) => {
            const union = DisjointIntervals.union(A, B);
            const inter = DisjointIntervals.intersection(A, B);
            const diff = DisjointIntervals.difference(A, B);
            const xor = DisjointIntervals.exclusiveOr(A, B);
            for (const x of probePoints(A, B)) {
                const a = contains(A, x), b = contains(B, x);
                expect(contains(union, x)).toBe(a || b);
                expect(contains(inter, x)).toBe(a && b);
                expect(contains(diff, x)).toBe(a && !b);
                expect(contains(xor, x)).toBe(a !== b);
            }
        });
    });

    it('the four operators agree with pointwise membership (double endpoints)', () => {
        check(fc.tuple(doubleSet(4), doubleSet(4)), ([A, B]) => {
            for (const x of probePoints(A, B)) {
                const a = contains(A, x), b = contains(B, x);
                expect(contains(DisjointIntervals.union(A, B), x)).toBe(a || b);
                expect(contains(DisjointIntervals.intersection(A, B), x)).toBe(a && b);
                expect(contains(DisjointIntervals.difference(A, B), x)).toBe(a && !b);
                expect(contains(DisjointIntervals.exclusiveOr(A, B), x)).toBe(a !== b);
            }
        });
    });

    it('every operator returns a canonical (strictly increasing) endpoint list', () => {
        check(fc.tuple(gridSet(5, 12), gridSet(5, 12)), ([A, B]) => {
            expectWellFormed(A);
            expectWellFormed(B);
            expectWellFormed(DisjointIntervals.union(A, B));
            expectWellFormed(DisjointIntervals.intersection(A, B));
            expectWellFormed(DisjointIntervals.difference(A, B));
            expectWellFormed(DisjointIntervals.exclusiveOr(A, B));
        });
    });

    it('insert/remove equal union/difference with a single interval', () => {
        check(fc.tuple(gridSet(4, 12), fc.integer({ min: 0, max: 11 }),
            fc.integer({ min: 1, max: 5 })), ([A, x, len]) => {
            const single = new DisjointIntervals(x, x + len);
            const inserted = A.clone();
            expect(inserted.insert(x, x + len)).toBe(true);
            expect(intervalsOf(inserted))
                .toEqual(intervalsOf(DisjointIntervals.union(A, single)));
            const removed = A.clone();
            expect(removed.remove(x, x + len)).toBe(true);
            expect(intervalsOf(removed))
                .toEqual(intervalsOf(DisjointIntervals.difference(A, single)));
        });
    });

    it('set-algebra laws hold (associativity, distributivity, De Morgan)', () => {
        check(fc.tuple(gridSet(4, 10), gridSet(4, 10), gridSet(4, 10)),
            ([A, B, C]) => {
                const u = DisjointIntervals.union;
                const n = DisjointIntervals.intersection;
                const d = DisjointIntervals.difference;
                expect(intervalsOf(u(u(A, B), C))).toEqual(intervalsOf(u(A, u(B, C))));
                expect(intervalsOf(n(n(A, B), C))).toEqual(intervalsOf(n(A, n(B, C))));
                expect(intervalsOf(n(A, u(B, C))))
                    .toEqual(intervalsOf(u(n(A, B), n(A, C))));
                // De Morgan relative to a universe U that contains every set.
                const U = new DisjointIntervals(-4, 20);
                expect(intervalsOf(d(U, u(A, B))))
                    .toEqual(intervalsOf(n(d(U, A), d(U, B))));
                expect(intervalsOf(d(U, n(A, B))))
                    .toEqual(intervalsOf(u(d(U, A), d(U, B))));
            });
    });

    it('operators only copy endpoints of their inputs (no arithmetic)', () => {
        // Every endpoint of a result must be an endpoint of one of the inputs;
        // a translation error that recomputed values would break this.
        check(fc.tuple(doubleSet(4), doubleSet(4)), ([A, B]) => {
            const inputEnds = new Set<number>();
            for (const set of [A, B]) {
                for (const [lo, hi] of intervalsOf(set)) {
                    inputEnds.add(lo);
                    inputEnds.add(hi);
                }
            }
            const ops = [DisjointIntervals.union, DisjointIntervals.intersection,
                DisjointIntervals.difference, DisjointIntervals.exclusiveOr];
            for (const op of ops) {
                for (const [lo, hi] of intervalsOf(op(A, B))) {
                    expect(inputEnds.has(lo)).toBe(true);
                    expect(inputEnds.has(hi)).toBe(true);
                }
            }
        });
    });

    it('degenerate and out-of-contract inputs behave as documented', () => {
        check(fc.tuple(finite(-5, 5), finite(-5, 5)), ([a, b]) => {
            const set = new DisjointIntervals(a, b);
            expect(set.getNumIntervals()).toBe(a < b ? 1 : 0);
            const other = new DisjointIntervals();
            expect(other.insert(a, b)).toBe(a < b);
            expect(other.remove(a, b)).toBe(a < b);
        });
        // NaN endpoints never satisfy xmin < xmax, so the set stays empty.
        expect(new DisjointIntervals(NaN, 1).getNumIntervals()).toBe(0);
        expect(new DisjointIntervals(0, NaN).getNumIntervals()).toBe(0);
        expect(new DisjointIntervals(0, 1).insert(NaN, NaN)).toBe(false);
        // Negative and too-large indices return null (upstream returns false).
        const single = new DisjointIntervals(0, 1);
        expect(single.getInterval(-1)).toBeNull();
        expect(single.getInterval(1)).toBeNull();
        // Empty operands.
        const empty = new DisjointIntervals();
        expect(intervalsOf(DisjointIntervals.union(empty, single))).toEqual([[0, 1]]);
        expect(intervalsOf(DisjointIntervals.union(single, empty))).toEqual([[0, 1]]);
        expect(DisjointIntervals.intersection(single, empty).getNumIntervals()).toBe(0);
        expect(intervalsOf(DisjointIntervals.difference(single, empty))).toEqual([[0, 1]]);
        expect(DisjointIntervals.difference(empty, single).getNumIntervals()).toBe(0);
        expect(intervalsOf(DisjointIntervals.exclusiveOr(empty, single))).toEqual([[0, 1]]);
    });
});
