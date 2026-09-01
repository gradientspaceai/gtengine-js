import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Vector } from '../src/Vector';
import {
    IntrAlignedBox2AlignedBox2TI,
    IntrAlignedBox2AlignedBox2FI
} from '../src/IntrAlignedBox2AlignedBox2';

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
