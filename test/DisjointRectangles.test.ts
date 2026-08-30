import { describe, it, expect } from 'vitest';
import { DisjointRectangles } from '../src/DisjointRectangles';
import { DisjointIntervals } from '../src/DisjointIntervals';

type Rect = [number, number, number, number];  // [xmin, xmax, ymin, ymax]

// Collect the rectangles of a set as [xmin, xmax, ymin, ymax] tuples.
function rectanglesOf(set: DisjointRectangles): Rect[] {
    const result: Rect[] = [];
    for (let i = 0; i < set.getNumRectangles(); ++i) {
        const r = set.getRectangle(i)!;
        result.push([r.xmin, r.xmax, r.ymin, r.ymax]);
    }
    return result;
}

// Build a set from a list of rectangles by Boolean-union insertion.
function makeSet(rects: Rect[]): DisjointRectangles {
    const set = new DisjointRectangles();
    for (const [xmin, xmax, ymin, ymax] of rects) {
        set.insert(xmin, xmax, ymin, ymax);
    }
    return set;
}

// Independent membership predicate: (x,y) is in the set when some rectangle
// [xmin,xmax)x[ymin,ymax) contains it.
function contains(set: DisjointRectangles, x: number, y: number): boolean {
    for (const [xmin, xmax, ymin, ymax] of rectanglesOf(set)) {
        if (xmin <= x && x < xmax && ymin <= y && y < ymax) {
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

// A random rectangle set built from a few random integer rectangles in
// [0, range) x [0, range).
function randomSet(rand: () => number, range: number): DisjointRectangles {
    const set = new DisjointRectangles();
    const count = 1 + Math.floor(rand() * 4);
    for (let k = 0; k < count; ++k) {
        const x0 = Math.floor(rand() * (range - 1));
        const x1 = x0 + 1 + Math.floor(rand() * (range - x0 - 1));
        const y0 = Math.floor(rand() * (range - 1));
        const y1 = y0 + 1 + Math.floor(rand() * (range - y0 - 1));
        set.insert(x0, x1, y0, y1);
    }
    return set;
}

// Verify the structural invariants of a disjoint rectangle set: strips are
// sorted in y and non-overlapping, and rectangles within each strip come
// from a disjoint interval set.
function verifyDisjoint(set: DisjointRectangles): void {
    let previousYMax = -Infinity;
    let rectangleCount = 0;
    for (let s = 0; s < set.getNumStrips(); ++s) {
        const strip = set.getStrip(s)!;
        expect(strip.ymin).toBeLessThan(strip.ymax);
        expect(strip.ymin).toBeGreaterThanOrEqual(previousYMax);
        previousYMax = strip.ymax;
        let previousXMax = -Infinity;
        for (let i = 0; i < strip.intervalSet.getNumIntervals(); ++i) {
            const interval = strip.intervalSet.getInterval(i)!;
            expect(interval.xmin).toBeLessThan(interval.xmax);
            expect(interval.xmin).toBeGreaterThanOrEqual(previousXMax);
            previousXMax = interval.xmax;
            ++rectangleCount;
        }
    }
    expect(set.getNumRectangles()).toBe(rectangleCount);
}

describe('DisjointRectangles', () => {
    it('default-constructs empty', () => {
        const set = new DisjointRectangles();
        expect(set.getNumRectangles()).toBe(0);
        expect(set.getNumStrips()).toBe(0);
        expect(set.getRectangle(0)).toBeNull();
        expect(set.getStrip(0)).toBeNull();
    });

    it('constructs a single rectangle when xmin < xmax and ymin < ymax', () => {
        const set = new DisjointRectangles(1, 3, 2, 5);
        expect(set.getNumRectangles()).toBe(1);
        expect(set.getNumStrips()).toBe(1);
        expect(set.getRectangle(0)).toEqual({ xmin: 1, xmax: 3, ymin: 2, ymax: 5 });
    });

    it('constructs empty for degenerate inputs', () => {
        expect(new DisjointRectangles(3, 1, 0, 1).getNumRectangles()).toBe(0);
        expect(new DisjointRectangles(0, 1, 5, 5).getNumRectangles()).toBe(0);
        expect(new DisjointRectangles(2, 2, 2, 2).getNumRectangles()).toBe(0);
    });

    it('getRectangle enumerates across strips and returns null out of range', () => {
        // Two strips: y in [0,1) has x-intervals [0,1) and [2,3); y in [5,6)
        // has x-interval [0,4).
        const set = makeSet([[0, 1, 0, 1], [2, 3, 0, 1], [0, 4, 5, 6]]);
        expect(set.getNumRectangles()).toBe(3);
        expect(rectanglesOf(set)).toEqual([
            [0, 1, 0, 1],
            [2, 3, 0, 1],
            [0, 4, 5, 6]
        ]);
        expect(set.getRectangle(3)).toBeNull();
        expect(set.getRectangle(-1)).toBeNull();
    });

    it('getStrip returns strip data and null out of range', () => {
        const set = makeSet([[0, 1, 0, 1], [2, 3, 0, 1], [0, 4, 5, 6]]);
        expect(set.getNumStrips()).toBe(2);
        const strip0 = set.getStrip(0)!;
        expect(strip0.ymin).toBe(0);
        expect(strip0.ymax).toBe(1);
        expect(strip0.intervalSet.getNumIntervals()).toBe(2);
        expect(strip0.intervalSet.getInterval(0)).toEqual({ xmin: 0, xmax: 1 });
        expect(strip0.intervalSet.getInterval(1)).toEqual({ xmin: 2, xmax: 3 });
        expect(set.getStrip(-1)).toBeNull();
        expect(set.getStrip(2)).toBeNull();
    });

    it('getStrip returns a copy of the interval set (value semantics)', () => {
        const set = new DisjointRectangles(0, 4, 0, 1);
        const strip = set.getStrip(0)!;
        strip.intervalSet.insert(10, 20);
        expect(set.getRectangle(0)).toEqual({ xmin: 0, xmax: 4, ymin: 0, ymax: 1 });
        expect(set.getNumRectangles()).toBe(1);
    });

    it('clear empties the set', () => {
        const set = new DisjointRectangles(0, 1, 0, 1);
        set.clear();
        expect(set.getNumRectangles()).toBe(0);
        expect(set.getNumStrips()).toBe(0);
    });

    it('clone makes an independent deep copy', () => {
        const set = new DisjointRectangles(0, 2, 0, 2);
        const copy = set.clone();
        set.insert(5, 6, 5, 6);
        expect(rectanglesOf(copy)).toEqual([[0, 2, 0, 2]]);
        expect(rectanglesOf(set)).toEqual([[0, 2, 0, 2], [5, 6, 5, 6]]);
    });

    describe('insert', () => {
        it('rejects degenerate rectangles', () => {
            const set = new DisjointRectangles();
            expect(set.insert(2, 2, 0, 1)).toBe(false);
            expect(set.insert(0, 1, 3, 1)).toBe(false);
            expect(set.getNumRectangles()).toBe(0);
        });

        it('merges horizontally abutting rectangles in the same strip', () => {
            const set = makeSet([[0, 1, 0, 1], [1, 2, 0, 1]]);
            expect(rectanglesOf(set)).toEqual([[0, 2, 0, 1]]);
        });

        it('keeps vertically abutting rectangles as separate strips', () => {
            const set = makeSet([[0, 1, 0, 1], [0, 1, 1, 2]]);
            expect(set.getNumStrips()).toBe(2);
            expect(rectanglesOf(set)).toEqual([[0, 1, 0, 1], [0, 1, 1, 2]]);
        });

        it('splits strips on partial y-overlap (cross shape)', () => {
            // Vertical bar [2,4)x[0,6) union horizontal bar [0,6)x[2,4).
            const set = makeSet([[2, 4, 0, 6], [0, 6, 2, 4]]);
            expect(rectanglesOf(set)).toEqual([
                [2, 4, 0, 2],
                [0, 6, 2, 4],
                [2, 4, 4, 6]
            ]);
            verifyDisjoint(set);
        });
    });

    describe('remove', () => {
        it('rejects degenerate rectangles', () => {
            const set = new DisjointRectangles(0, 10, 0, 10);
            expect(set.remove(4, 4, 0, 1)).toBe(false);
            expect(rectanglesOf(set)).toEqual([[0, 10, 0, 10]]);
        });

        it('punches a hole, leaving a frame of four rectangles', () => {
            const set = new DisjointRectangles(0, 10, 0, 10);
            expect(set.remove(3, 7, 4, 6)).toBe(true);
            expect(rectanglesOf(set)).toEqual([
                [0, 10, 0, 4],
                [0, 3, 4, 6],
                [7, 10, 4, 6],
                [0, 10, 6, 10]
            ]);
            verifyDisjoint(set);
        });

        it('removing a superset empties the set', () => {
            const set = makeSet([[1, 2, 1, 2], [3, 4, 3, 4]]);
            set.remove(0, 10, 0, 10);
            expect(set.getNumRectangles()).toBe(0);
        });
    });

    describe('static set operations (hand-computed)', () => {
        // Two overlapping unit-aligned squares: A = [0,4)x[0,4),
        // B = [2,6)x[2,6).
        const A = new DisjointRectangles(0, 4, 0, 4);
        const B = new DisjointRectangles(2, 6, 2, 6);

        it('union', () => {
            const result = DisjointRectangles.union(A, B);
            expect(rectanglesOf(result)).toEqual([
                [0, 4, 0, 2],
                [0, 6, 2, 4],
                [2, 6, 4, 6]
            ]);
            verifyDisjoint(result);
        });

        it('intersection', () => {
            const result = DisjointRectangles.intersection(A, B);
            expect(rectanglesOf(result)).toEqual([[2, 4, 2, 4]]);
        });

        it('difference', () => {
            expect(rectanglesOf(DisjointRectangles.difference(A, B))).toEqual([
                [0, 4, 0, 2],
                [0, 2, 2, 4]
            ]);
            expect(rectanglesOf(DisjointRectangles.difference(B, A))).toEqual([
                [4, 6, 2, 4],
                [2, 6, 4, 6]
            ]);
        });

        it('exclusiveOr', () => {
            const result = DisjointRectangles.exclusiveOr(A, B);
            expect(rectanglesOf(result)).toEqual([
                [0, 4, 0, 2],
                [0, 2, 2, 4],
                [4, 6, 2, 4],
                [2, 6, 4, 6]
            ]);
            verifyDisjoint(result);
        });

        it('operations with an empty set', () => {
            const empty = new DisjointRectangles();
            expect(rectanglesOf(DisjointRectangles.union(A, empty))).toEqual(rectanglesOf(A));
            expect(rectanglesOf(DisjointRectangles.union(empty, B))).toEqual(rectanglesOf(B));
            expect(DisjointRectangles.intersection(A, empty).getNumRectangles()).toBe(0);
            expect(rectanglesOf(DisjointRectangles.difference(A, empty))).toEqual(rectanglesOf(A));
            expect(DisjointRectangles.difference(empty, A).getNumRectangles()).toBe(0);
            expect(rectanglesOf(DisjointRectangles.exclusiveOr(A, empty))).toEqual(rectanglesOf(A));
        });

        it('xor of identical sets is empty', () => {
            const C = makeSet([[0, 3, 0, 3], [5, 7, 1, 2]]);
            expect(DisjointRectangles.exclusiveOr(C, C).getNumRectangles()).toBe(0);
            expect(DisjointRectangles.difference(C, C).getNumRectangles()).toBe(0);
        });

        it('disjoint sets: union keeps both, intersection empty', () => {
            const L = new DisjointRectangles(0, 1, 0, 1);
            const R = new DisjointRectangles(5, 6, 5, 6);
            expect(rectanglesOf(DisjointRectangles.union(L, R)))
                .toEqual([[0, 1, 0, 1], [5, 6, 5, 6]]);
            expect(DisjointRectangles.intersection(L, R).getNumRectangles()).toBe(0);
        });
    });

    describe('randomized cross-checks', () => {
        it('operations agree with pointwise membership over random sets', () => {
            const rand = makeRandom(0xb14);
            const range = 12;
            for (let trial = 0; trial < 40; ++trial) {
                const A = randomSet(rand, range);
                const B = randomSet(rand, range);
                const union = DisjointRectangles.union(A, B);
                const inter = DisjointRectangles.intersection(A, B);
                const diff = DisjointRectangles.difference(A, B);
                const xor = DisjointRectangles.exclusiveOr(A, B);
                verifyDisjoint(union);
                verifyDisjoint(inter);
                verifyDisjoint(diff);
                verifyDisjoint(xor);
                // Sample at half-integers to probe every unit cell strictly
                // inside [0, range)^2, avoiding endpoint ambiguity.
                for (let ky = 0; ky < range; ++ky) {
                    for (let kx = 0; kx < range; ++kx) {
                        const x = kx + 0.5;
                        const y = ky + 0.5;
                        const inA = contains(A, x, y);
                        const inB = contains(B, x, y);
                        expect(contains(union, x, y)).toBe(inA || inB);
                        expect(contains(inter, x, y)).toBe(inA && inB);
                        expect(contains(diff, x, y)).toBe(inA && !inB);
                        expect(contains(xor, x, y)).toBe(inA !== inB);
                    }
                }
            }
        });

        it('structural identities hold on random sets', () => {
            const rand = makeRandom(0xd15c0);
            for (let trial = 0; trial < 40; ++trial) {
                const A = randomSet(rand, 12);
                const B = randomSet(rand, 12);
                // Commutativity.
                expect(rectanglesOf(DisjointRectangles.union(A, B)))
                    .toEqual(rectanglesOf(DisjointRectangles.union(B, A)));
                expect(rectanglesOf(DisjointRectangles.intersection(A, B)))
                    .toEqual(rectanglesOf(DisjointRectangles.intersection(B, A)));
                // Idempotence.
                expect(rectanglesOf(DisjointRectangles.union(A, A))).toEqual(rectanglesOf(A));
                expect(rectanglesOf(DisjointRectangles.intersection(A, A))).toEqual(rectanglesOf(A));
                expect(DisjointRectangles.difference(A, A).getNumRectangles()).toBe(0);
                // (A - B) intersect B is empty.
                expect(DisjointRectangles.intersection(
                    DisjointRectangles.difference(A, B), B).getNumRectangles()).toBe(0);
            }
        });
    });

    it('interoperates with DisjointIntervals strip contents', () => {
        const set = makeSet([[0, 2, 0, 1], [4, 6, 0, 1]]);
        const strip = set.getStrip(0)!;
        const other = new DisjointIntervals(1, 5);
        const merged = DisjointIntervals.union(strip.intervalSet, other);
        expect(merged.getNumIntervals()).toBe(1);
        expect(merged.getInterval(0)).toEqual({ xmin: 0, xmax: 6 });
    });
});
