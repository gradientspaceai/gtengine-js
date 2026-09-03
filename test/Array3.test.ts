import { describe, it, expect } from 'vitest';
import { Array3 } from '../src/Array3.js';

describe('Array3', () => {
    it('reports its bounds', () => {
        const a = new Array3<number>(4, 3, 2);
        expect(a.getBound0()).toBe(4);
        expect(a.getBound1()).toBe(3);
        expect(a.getBound2()).toBe(2);
        expect(a.data().length).toBe(24);
    });

    it('stores elements in lexicographical order (i0 fastest, i2 slowest)', () => {
        const a = new Array3<number>(2, 2, 2);
        let flat = 0;
        for (let i2 = 0; i2 < 2; ++i2) {
            for (let i1 = 0; i1 < 2; ++i1) {
                for (let i0 = 0; i0 < 2; ++i0) {
                    a.set(i0, i1, i2, flat++);
                }
            }
        }
        // Writing in lexicographical loop order must produce 0..7 in the
        // flat storage: flat index is i0 + bound0 * (i1 + bound1 * i2).
        expect(a.data()).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);
        expect(a.get(1, 1, 1)).toBe(7);
        expect(a.get(0, 1, 0)).toBe(2);
        expect(a.get(0, 0, 1)).toBe(4);
    });

    it('round trips get/set over the full index space', () => {
        const a = new Array3<string>(3, 2, 4);
        for (let i2 = 0; i2 < 4; ++i2) {
            for (let i1 = 0; i1 < 2; ++i1) {
                for (let i0 = 0; i0 < 3; ++i0) {
                    a.set(i0, i1, i2, `${i0},${i1},${i2}`);
                }
            }
        }
        for (let i2 = 0; i2 < 4; ++i2) {
            for (let i1 = 0; i1 < 2; ++i1) {
                for (let i0 = 0; i0 < 3; ++i0) {
                    expect(a.get(i0, i1, i2)).toBe(`${i0},${i1},${i2}`);
                }
            }
        }
    });

    it('aliases caller-owned storage', () => {
        const objects = Array.from({ length: 8 }, (_, i) => i);
        const a = new Array3<number>(2, 2, 2, objects);
        expect(a.get(1, 0, 1)).toBe(5);
        a.set(1, 0, 1, -5);
        expect(objects[5]).toBe(-5);
    });

    it('rejects caller-owned storage of the wrong length', () => {
        expect(() => new Array3<number>(2, 2, 2, [1, 2, 3])).toThrow();
    });

    it('fill sets every element', () => {
        const a = new Array3<number>(2, 1, 2);
        a.fill(3);
        expect(a.data()).toEqual([3, 3, 3, 3]);
    });
});
