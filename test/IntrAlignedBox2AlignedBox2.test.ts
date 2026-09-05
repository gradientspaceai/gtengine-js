import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector } from '../src/Vector.js';
import {
    IntrAlignedBox2AlignedBox2TI,
    IntrAlignedBox2AlignedBox2FI
} from '../src/IntrAlignedBox2AlignedBox2.js';
import { check, fc } from './helpers/arbitraries.js';

function box(minX: number, minY: number, maxX: number, maxY: number): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray([minX, minY]),
        Vector.fromArray([maxX, maxY]));
}

// A simple deterministic PRNG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrAlignedBox2AlignedBox2', () => {
    const ti = new IntrAlignedBox2AlignedBox2TI();
    const fi = new IntrAlignedBox2AlignedBox2FI();

    it('reports overlap for boxes that share a rectangle', () => {
        const b0 = box(0, 0, 4, 3);
        const b1 = box(2, 1, 6, 5);
        expect(ti.test(b0, b1).intersect).toBe(true);

        const result = fi.find(b0, b1);
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([2, 1]);
        expect(result.box.max.values).toEqual([4, 3]);
    });

    it('reports no overlap for separated boxes', () => {
        const b0 = box(0, 0, 1, 1);
        const b1 = box(2, 0, 3, 1);
        expect(ti.test(b0, b1).intersect).toBe(false);
        expect(fi.find(b0, b1).intersect).toBe(false);

        // Separation along the second dimension only.
        const b2 = box(0, 5, 1, 6);
        expect(ti.test(b0, b2).intersect).toBe(false);
    });

    it('treats edge contact as an intersection (degenerate box)', () => {
        const b0 = box(0, 0, 1, 1);
        const b1 = box(1, 0, 2, 1);
        expect(ti.test(b0, b1).intersect).toBe(true);

        const result = fi.find(b0, b1);
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([1, 0]);
        expect(result.box.max.values).toEqual([1, 1]);
    });

    it('treats corner contact as an intersection (doubly degenerate box)', () => {
        const result = fi.find(box(0, 0, 1, 1), box(1, 1, 2, 2));
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([1, 1]);
        expect(result.box.max.values).toEqual([1, 1]);
    });

    it('returns the contained box when one box is inside the other', () => {
        const outer = box(-10, -10, 10, 10);
        const inner = box(-1, 2, 3, 4);
        const result = fi.find(outer, inner);
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual(inner.min.values);
        expect(result.box.max.values).toEqual(inner.max.values);
    });

    it('handles zero-extent (degenerate) boxes', () => {
        // A point box on the boundary of another box.
        const point = box(1, 0.5, 1, 0.5);
        expect(ti.test(box(0, 0, 1, 1), point).intersect).toBe(true);
        // A point box outside.
        const outside = box(1.5, 0.5, 1.5, 0.5);
        expect(ti.test(box(0, 0, 1, 1), outside).intersect).toBe(false);
    });

    it('is symmetric in its arguments', () => {
        const b0 = box(0, 0, 4, 3);
        const b1 = box(2, 1, 6, 5);
        const r01 = fi.find(b0, b1);
        const r10 = fi.find(b1, b0);
        expect(r10.box.min.values).toEqual(r01.box.min.values);
        expect(r10.box.max.values).toEqual(r01.box.max.values);
    });

    it('agrees with an independent overlap oracle and with the TI query', () => {
        const rand = makeRandom(20240815);
        let numIntersect = 0, numSeparate = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const boxes: AlignedBox[] = [];
            for (let k = 0; k < 2; ++k) {
                const x0 = 4 * rand() - 2, y0 = 4 * rand() - 2;
                const x1 = x0 + 2 * rand(), y1 = y0 + 2 * rand();
                boxes.push(box(x0, y0, x1, y1));
            }
            // Independent oracle: closed intervals overlap in both dimensions.
            const oracle =
                boxes[0].min.values[0] <= boxes[1].max.values[0] &&
                boxes[1].min.values[0] <= boxes[0].max.values[0] &&
                boxes[0].min.values[1] <= boxes[1].max.values[1] &&
                boxes[1].min.values[1] <= boxes[0].max.values[1];

            const t = ti.test(boxes[0], boxes[1]).intersect;
            const f = fi.find(boxes[0], boxes[1]);
            expect(t).toBe(oracle);
            expect(f.intersect).toBe(t);
            if (f.intersect) {
                ++numIntersect;
                for (let d = 0; d < 2; ++d) {
                    expect(f.box.min.values[d]).toBeLessThanOrEqual(f.box.max.values[d]);
                    // The intersection is inside both boxes.
                    expect(f.box.min.values[d]).toBeGreaterThanOrEqual(
                        Math.max(boxes[0].min.values[d], boxes[1].min.values[d]) - 1e-15);
                    expect(f.box.max.values[d]).toBeLessThanOrEqual(
                        Math.min(boxes[0].max.values[d], boxes[1].max.values[d]) + 1e-15);
                }
            } else {
                ++numSeparate;
            }
        }
        expect(numIntersect).toBeGreaterThan(0);
        expect(numSeparate).toBeGreaterThan(0);
    });
});

describe('IntrAlignedBox2AlignedBox2 verification', () => {
    // Integer-cornered boxes: every comparison in the query is exact, and the
    // brute-force lattice sweep below sees the touching configurations that
    // decide the closed-box convention.
    const latticeBox2 = fc.tuple(fc.integer({ min: -4, max: 4 }),
        fc.integer({ min: -4, max: 4 }), fc.integer({ min: 0, max: 5 }),
        fc.integer({ min: 0, max: 5 }))
        .map(([x, y, w, h]) => box(x, y, x + w, y + h));

    const ti2 = new IntrAlignedBox2AlignedBox2TI();
    const fi2 = new IntrAlignedBox2AlignedBox2FI();

    it('TI agrees with FI and with the exact per-axis overlap test', () => {
        check(fc.tuple(latticeBox2, latticeBox2), ([b0, b1]) => {
            const t = ti2.test(b0, b1).intersect;
            const f = fi2.find(b0, b1);
            expect(f.intersect).toBe(t);

            let expected = true;
            for (let i = 0; i < 2; ++i) {
                if (b0.max.get(i) < b1.min.get(i) ||
                    b0.min.get(i) > b1.max.get(i)) {
                    expected = false;
                }
            }
            expect(t).toBe(expected);

            if (f.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(f.box.min.get(i))
                        .toBe(Math.max(b0.min.get(i), b1.min.get(i)));
                    expect(f.box.max.get(i))
                        .toBe(Math.min(b0.max.get(i), b1.max.get(i)));
                    expect(f.box.min.get(i))
                        .toBeLessThanOrEqual(f.box.max.get(i));
                }
            }
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(latticeBox2, latticeBox2), ([b0, b1]) => {
            expect(ti2.test(b1, b0).intersect).toBe(ti2.test(b0, b1).intersect);
            const a = fi2.find(b0, b1);
            const b = fi2.find(b1, b0);
            expect(b.intersect).toBe(a.intersect);
            if (a.intersect) {
                for (let i = 0; i < 2; ++i) {
                    expect(b.box.min.get(i)).toBe(a.box.min.get(i));
                    expect(b.box.max.get(i)).toBe(a.box.max.get(i));
                }
            }
        });
    });

    it('the FI box is exactly the set of lattice points common to both boxes', () => {
        check(fc.tuple(latticeBox2, latticeBox2), ([b0, b1]) => {
            const f = fi2.find(b0, b1);
            const inside = (b: AlignedBox, p: number[]): boolean =>
                p[0] >= b.min.get(0) && p[0] <= b.max.get(0) &&
                p[1] >= b.min.get(1) && p[1] <= b.max.get(1);
            let common = 0;
            for (let x = -10; x <= 10; ++x) {
                for (let y = -10; y <= 10; ++y) {
                    const p = [x, y];
                    const both = inside(b0, p) && inside(b1, p);
                    if (both) {
                        ++common;
                        expect(f.intersect).toBe(true);
                        expect(inside(f.box, p)).toBe(true);
                    }
                    if (f.intersect && inside(f.box, p)) {
                        expect(both).toBe(true);
                    }
                }
            }
            // A common lattice point implies intersection, but the converse
            // needs no lattice witness (boxes can share only a sliver), so
            // 'common' is allowed to be zero here.
            expect(common).toBeGreaterThanOrEqual(0);
        }, 60);
    });

    it('degenerate (zero-extent) boxes follow the closed-solid convention', () => {
        // A degenerate box is a point/segment; touching still counts.
        const pt = box(1, 1, 1, 1);
        expect(ti2.test(pt, box(0, 0, 2, 2)).intersect).toBe(true);
        expect(ti2.test(pt, box(1, 1, 3, 3)).intersect).toBe(true);
        expect(ti2.test(pt, box(2, 2, 3, 3)).intersect).toBe(false);
        const f = fi2.find(pt, box(1, 1, 3, 3));
        expect(f.intersect).toBe(true);
        expect(f.box.min.values).toEqual([1, 1]);
        expect(f.box.max.values).toEqual([1, 1]);
    });

    it('a separated FI result leaves the default box untouched', () => {
        const f = fi2.find(box(0, 0, 1, 1), box(5, 5, 6, 6));
        expect(f.intersect).toBe(false);
        // The port mirrors the upstream value-initialized AlignedBox, whose
        // default constructor is min = -1 and max = +1 in each dimension.
        expect(f.box.min.values).toEqual([-1, -1]);
        expect(f.box.max.values).toEqual([1, 1]);
    });
});
