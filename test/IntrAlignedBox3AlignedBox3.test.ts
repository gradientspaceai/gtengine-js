import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector } from '../src/Vector.js';
import {
    IntrAlignedBox3AlignedBox3TI,
    IntrAlignedBox3AlignedBox3FI
} from '../src/IntrAlignedBox3AlignedBox3.js';

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(Vector.fromArray(min), Vector.fromArray(max));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrAlignedBox3AlignedBox3', () => {
    const ti = new IntrAlignedBox3AlignedBox3TI();
    const fi = new IntrAlignedBox3AlignedBox3FI();

    it('reports overlap and the intersection box', () => {
        const b0 = box([0, 0, 0], [4, 3, 2]);
        const b1 = box([2, 1, -1], [6, 5, 1]);
        expect(ti.test(b0, b1).intersect).toBe(true);

        const result = fi.find(b0, b1);
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([2, 1, 0]);
        expect(result.box.max.values).toEqual([4, 3, 1]);
    });

    it('reports no overlap when separated in any single dimension', () => {
        const b0 = box([0, 0, 0], [1, 1, 1]);
        expect(ti.test(b0, box([2, 0, 0], [3, 1, 1])).intersect).toBe(false);
        expect(ti.test(b0, box([0, -3, 0], [1, -2, 1])).intersect).toBe(false);
        expect(ti.test(b0, box([0, 0, 5], [1, 1, 6])).intersect).toBe(false);
        expect(fi.find(b0, box([0, 0, 5], [1, 1, 6])).intersect).toBe(false);
    });

    it('treats face contact as an intersection with a flat box', () => {
        const result = fi.find(box([0, 0, 0], [1, 1, 1]), box([1, 0, 0], [2, 1, 1]));
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([1, 0, 0]);
        expect(result.box.max.values).toEqual([1, 1, 1]);
    });

    it('treats corner contact as an intersection with a point box', () => {
        const result = fi.find(box([0, 0, 0], [1, 1, 1]), box([1, 1, 1], [2, 2, 2]));
        expect(result.intersect).toBe(true);
        expect(result.box.min.values).toEqual([1, 1, 1]);
        expect(result.box.max.values).toEqual([1, 1, 1]);
    });

    it('returns the contained box for nested boxes', () => {
        const inner = box([-1, 0, 0.5], [1, 2, 0.75]);
        const result = fi.find(box([-5, -5, -5], [5, 5, 5]), inner);
        expect(result.box.min.values).toEqual(inner.min.values);
        expect(result.box.max.values).toEqual(inner.max.values);
    });

    it('handles zero-extent boxes', () => {
        const point = box([0.5, 0.5, 1], [0.5, 0.5, 1]);
        expect(ti.test(box([0, 0, 0], [1, 1, 1]), point).intersect).toBe(true);
        const away = box([0.5, 0.5, 1.25], [0.5, 0.5, 1.25]);
        expect(ti.test(box([0, 0, 0], [1, 1, 1]), away).intersect).toBe(false);
    });

    it('agrees with an independent overlap oracle and with the TI query', () => {
        const rand = makeRandom(775533);
        let numIntersect = 0, numSeparate = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const boxes: AlignedBox[] = [];
            for (let k = 0; k < 2; ++k) {
                const lo: number[] = [], hi: number[] = [];
                for (let d = 0; d < 3; ++d) {
                    const a = 4 * rand() - 2;
                    lo.push(a);
                    hi.push(a + 2 * rand());
                }
                boxes.push(box(lo, hi));
            }
            let oracle = true;
            for (let d = 0; d < 3; ++d) {
                if (boxes[0].min.values[d] > boxes[1].max.values[d] ||
                    boxes[1].min.values[d] > boxes[0].max.values[d]) {
                    oracle = false;
                }
            }

            const t = ti.test(boxes[0], boxes[1]).intersect;
            const f = fi.find(boxes[0], boxes[1]);
            expect(t).toBe(oracle);
            expect(f.intersect).toBe(t);
            if (f.intersect) {
                ++numIntersect;
                // The center of the reported box is inside both boxes.
                for (let d = 0; d < 3; ++d) {
                    const c = 0.5 * (f.box.min.values[d] + f.box.max.values[d]);
                    for (let k = 0; k < 2; ++k) {
                        expect(c).toBeGreaterThanOrEqual(boxes[k].min.values[d] - 1e-15);
                        expect(c).toBeLessThanOrEqual(boxes[k].max.values[d] + 1e-15);
                    }
                }
            } else {
                ++numSeparate;
            }
        }
        expect(numIntersect).toBeGreaterThan(0);
        expect(numSeparate).toBeGreaterThan(0);
    });
});
