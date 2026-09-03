import { describe, it, expect } from 'vitest';
import { Array2 } from '../src/Array2.js';

describe('Array2', () => {
    it('reports its bounds', () => {
        const a = new Array2<number>(3, 2);
        expect(a.getBound0()).toBe(3);
        expect(a.getBound1()).toBe(2);
        expect(a.data().length).toBe(6);
    });

    it('defaults to empty bounds', () => {
        const a = new Array2<number>();
        expect(a.getBound0()).toBe(0);
        expect(a.getBound1()).toBe(0);
        expect(a.data().length).toBe(0);
    });

    it('stores elements in lexicographical order (i0 fastest)', () => {
        const a = new Array2<number>(3, 2);
        for (let i1 = 0; i1 < 2; ++i1) {
            for (let i0 = 0; i0 < 3; ++i0) {
                a.set(i0, i1, 10 * i1 + i0);
            }
        }
        // Flat index is i0 + bound0 * i1.
        expect(a.data()).toEqual([0, 1, 2, 10, 11, 12]);
        expect(a.get(2, 1)).toBe(12);
        expect(a.get(0, 1)).toBe(10);
    });

    it('round trips get/set over the full index space', () => {
        const a = new Array2<string>(4, 3);
        for (let i1 = 0; i1 < 3; ++i1) {
            for (let i0 = 0; i0 < 4; ++i0) {
                a.set(i0, i1, `${i0},${i1}`);
            }
        }
        for (let i1 = 0; i1 < 3; ++i1) {
            for (let i0 = 0; i0 < 4; ++i0) {
                expect(a.get(i0, i1)).toBe(`${i0},${i1}`);
            }
        }
    });

    it('aliases caller-owned storage (the second upstream constructor)', () => {
        const objects = [1, 2, 3, 4, 5, 6];
        const a = new Array2<number>(3, 2, objects);
        expect(a.get(0, 0)).toBe(1);
        expect(a.get(1, 1)).toBe(5);
        a.set(1, 1, -5);
        expect(objects[1 + 3 * 1]).toBe(-5);
        objects[0] = 99;
        expect(a.get(0, 0)).toBe(99);
    });

    it('rejects caller-owned storage of the wrong length', () => {
        expect(() => new Array2<number>(3, 2, [1, 2, 3])).toThrow();
    });

    it('fill sets every element', () => {
        const a = new Array2<number>(2, 2);
        a.fill(7);
        expect(a.data()).toEqual([7, 7, 7, 7]);
    });
});
