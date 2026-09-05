import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Vector } from '../src/Vector.js';
import {
    IntrAlignedBox3AlignedBox3TI,
    IntrAlignedBox3AlignedBox3FI
} from '../src/IntrAlignedBox3AlignedBox3.js';
import { check, fc } from './helpers/arbitraries.js';

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

describe('IntrAlignedBox3AlignedBox3 verification', () => {
    // Integer-cornered boxes: every comparison in the query is exact, and the
    // touching configurations that decide the closed-box convention are hit
    // often.
    const latticeBox3 = fc.tuple(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        fc.array(fc.integer({ min: 0, max: 5 }), { minLength: 3, maxLength: 3 }))
        .map(([lo, ext]) => box(lo, [lo[0] + ext[0], lo[1] + ext[1], lo[2] + ext[2]]));

    const ti3 = new IntrAlignedBox3AlignedBox3TI();
    const fi3 = new IntrAlignedBox3AlignedBox3FI();

    it('TI agrees with FI and with the exact per-axis overlap test', () => {
        check(fc.tuple(latticeBox3, latticeBox3), ([b0, b1]) => {
            const t = ti3.test(b0, b1).intersect;
            const f = fi3.find(b0, b1);
            expect(f.intersect).toBe(t);

            let expected = true;
            for (let i = 0; i < 3; ++i) {
                if (b0.max.get(i) < b1.min.get(i) ||
                    b0.min.get(i) > b1.max.get(i)) {
                    expected = false;
                }
            }
            expect(t).toBe(expected);

            if (f.intersect) {
                for (let i = 0; i < 3; ++i) {
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
        check(fc.tuple(latticeBox3, latticeBox3), ([b0, b1]) => {
            expect(ti3.test(b1, b0).intersect).toBe(ti3.test(b0, b1).intersect);
            const a = fi3.find(b0, b1);
            const b = fi3.find(b1, b0);
            expect(b.intersect).toBe(a.intersect);
            if (a.intersect) {
                for (let i = 0; i < 3; ++i) {
                    expect(b.box.min.get(i)).toBe(a.box.min.get(i));
                    expect(b.box.max.get(i)).toBe(a.box.max.get(i));
                }
            }
        });
    });

    it('the FI box is exactly the set of lattice points common to both boxes', () => {
        check(fc.tuple(latticeBox3, latticeBox3), ([b0, b1]) => {
            const f = fi3.find(b0, b1);
            const inside = (b: AlignedBox, p: number[]): boolean => {
                for (let i = 0; i < 3; ++i) {
                    if (p[i] < b.min.get(i) || p[i] > b.max.get(i)) {
                        return false;
                    }
                }
                return true;
            };
            for (let x = -6; x <= 10; ++x) {
                for (let y = -6; y <= 10; ++y) {
                    for (let z = -6; z <= 10; ++z) {
                        const p = [x, y, z];
                        const both = inside(b0, p) && inside(b1, p);
                        if (both) {
                            expect(f.intersect).toBe(true);
                            expect(inside(f.box, p)).toBe(true);
                        }
                        if (f.intersect && inside(f.box, p)) {
                            expect(both).toBe(true);
                        }
                    }
                }
            }
        }, 40);
    }, 30000);

    it('degenerate (zero-extent) boxes follow the closed-solid convention', () => {
        const pt = box([1, 1, 1], [1, 1, 1]);
        expect(ti3.test(pt, box([0, 0, 0], [2, 2, 2])).intersect).toBe(true);
        expect(ti3.test(pt, box([1, 1, 1], [3, 3, 3])).intersect).toBe(true);
        expect(ti3.test(pt, box([2, 2, 2], [3, 3, 3])).intersect).toBe(false);
        const f = fi3.find(pt, box([1, 1, 1], [3, 3, 3]));
        expect(f.intersect).toBe(true);
        expect(f.box.min.values).toEqual([1, 1, 1]);
        expect(f.box.max.values).toEqual([1, 1, 1]);
    });

    it('a separated FI result leaves the default box untouched', () => {
        const f = fi3.find(box([0, 0, 0], [1, 1, 1]), box([5, 5, 5], [6, 6, 6]));
        expect(f.intersect).toBe(false);
        // The port mirrors the upstream value-initialized AlignedBox, whose
        // default constructor is min = -1 and max = +1 in each dimension.
        expect(f.box.min.values).toEqual([-1, -1, -1]);
        expect(f.box.max.values).toEqual([1, 1, 1]);
    });
});
