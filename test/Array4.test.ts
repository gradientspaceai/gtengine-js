import { describe, it, expect } from 'vitest';
import { Array4 } from '../src/Array4';

describe('Array4', () => {
    it('reports its bounds', () => {
        const a = new Array4<number>(5, 4, 3, 2);
        expect(a.getBound0()).toBe(5);
        expect(a.getBound1()).toBe(4);
        expect(a.getBound2()).toBe(3);
        expect(a.getBound3()).toBe(2);
        expect(a.data().length).toBe(120);
    });

    it('stores elements in lexicographical order (i0 fastest, i3 slowest)', () => {
        const a = new Array4<number>(2, 2, 2, 2);
        let flat = 0;
        for (let i3 = 0; i3 < 2; ++i3) {
            for (let i2 = 0; i2 < 2; ++i2) {
                for (let i1 = 0; i1 < 2; ++i1) {
                    for (let i0 = 0; i0 < 2; ++i0) {
                        a.set(i0, i1, i2, i3, flat++);
                    }
                }
            }
        }
        expect(a.data()).toEqual([...Array(16).keys()]);
        // Flat index is i0 + bound0 * (i1 + bound1 * (i2 + bound2 * i3)).
        expect(a.get(1, 1, 1, 1)).toBe(15);
        expect(a.get(0, 0, 0, 1)).toBe(8);
        expect(a.get(1, 0, 1, 0)).toBe(5);
    });

    it('round trips get/set over the full index space', () => {
        const a = new Array4<string>(2, 3, 2, 2);
        for (let i3 = 0; i3 < 2; ++i3) {
            for (let i2 = 0; i2 < 2; ++i2) {
                for (let i1 = 0; i1 < 3; ++i1) {
                    for (let i0 = 0; i0 < 2; ++i0) {
                        a.set(i0, i1, i2, i3, `${i0},${i1},${i2},${i3}`);
                    }
                }
            }
        }
        for (let i3 = 0; i3 < 2; ++i3) {
            for (let i2 = 0; i2 < 2; ++i2) {
                for (let i1 = 0; i1 < 3; ++i1) {
                    for (let i0 = 0; i0 < 2; ++i0) {
                        expect(a.get(i0, i1, i2, i3)).toBe(`${i0},${i1},${i2},${i3}`);
                    }
                }
            }
        }
    });

    it('aliases caller-owned storage', () => {
        const objects = Array.from({ length: 16 }, (_, i) => i);
        const a = new Array4<number>(2, 2, 2, 2, objects);
        expect(a.get(0, 1, 0, 1)).toBe(10);
        a.set(0, 1, 0, 1, -10);
        expect(objects[10]).toBe(-10);
    });

    it('rejects caller-owned storage of the wrong length', () => {
        expect(() => new Array4<number>(2, 2, 2, 2, [1, 2, 3])).toThrow();
    });

    it('fill sets every element', () => {
        const a = new Array4<number>(1, 2, 1, 2);
        a.fill(-1);
        expect(a.data()).toEqual([-1, -1, -1, -1]);
    });
});
