import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
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

describe('Array2 verification', () => {
    it('flat index equals i0 + bound0 * i1 for random bounds', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 8 }), fc.integer({ min: 1, max: 8 })),
            ([b0, b1]) => {
                const a = new Array2<number>(b0, b1);
                a.fill(0);
                for (let i1 = 0; i1 < b1; ++i1) {
                    for (let i0 = 0; i0 < b0; ++i0) {
                        a.set(i0, i1, i0 + b0 * i1);
                    }
                }
                const flat = a.data();
                for (let k = 0; k < b0 * b1; ++k) {
                    if (flat[k] !== k) { return false; }
                }
                return true;
            });
    });

    it('get/set round trips; every cell keeps the value written to it', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 6 }), fc.integer({ min: 1, max: 6 }),
            fc.array(finite(), { minLength: 36, maxLength: 36 })),
            ([b0, b1, vals]) => {
                const a = new Array2<number>(b0, b1);
                a.fill(Number.NaN);
                const ref = new Map<string, number>();
                for (let i1 = 0; i1 < b1; ++i1) {
                    for (let i0 = 0; i0 < b0; ++i0) {
                        const v = vals[i0 + b0 * i1]!;
                        a.set(i0, i1, v);
                        ref.set(`${i0},${i1}`, v);
                    }
                }
                for (const [key, v] of ref) {
                    const parts = key.split(',').map(Number);
                    if (!Object.is(a.get(parts[0]!, parts[1]!), v)) { return false; }
                }
                return true;
            });
    });

    it('caller-owned storage stays aliased in both directions', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 5 }), fc.integer({ min: 1, max: 5 }), finite()),
            ([b0, b1, v]) => {
                const objects = new Array<number>(b0 * b1).fill(0);
                const a = new Array2<number>(b0, b1, objects);
                const i0 = b0 - 1, i1 = b1 - 1;
                a.set(i0, i1, v);
                if (!Object.is(objects[i0 + b0 * i1], v)) { return false; }
                objects[0] = v + 1;
                return Object.is(a.get(0, 0), v + 1);
            });
    });

    it('fill covers exactly bound0 * bound1 entries', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 7 }), fc.integer({ min: 0, max: 7 })),
            ([b0, b1]) => {
                const a = new Array2<number>(b0, b1);
                a.fill(5);
                const flat = a.data();
                if (flat.length !== b0 * b1) { return false; }
                return flat.every(x => x === 5);
            });
    });
});
