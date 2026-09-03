import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { Array2 } from '../src/Array2.js';
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

describe('Array3 verification', () => {
    it('flat index equals i0 + b0 * (i1 + b1 * i2) for random bounds', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 5 }), fc.integer({ min: 1, max: 5 }),
            fc.integer({ min: 1, max: 5 })),
            ([b0, b1, b2]) => {
                const a = new Array3<number>(b0, b1, b2);
                a.fill(0);
                for (let i2 = 0; i2 < b2; ++i2) {
                    for (let i1 = 0; i1 < b1; ++i1) {
                        for (let i0 = 0; i0 < b0; ++i0) {
                            a.set(i0, i1, i2, i0 + b0 * (i1 + b1 * i2));
                        }
                    }
                }
                const flat = a.data();
                for (let k = 0; k < b0 * b1 * b2; ++k) {
                    if (flat[k] !== k) { return false; }
                }
                return true;
            });
    });

    it('get/set round trips over the whole index space', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }), fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 }), fc.array(finite(), { minLength: 64, maxLength: 64 })),
            ([b0, b1, b2, vals]) => {
                const a = new Array3<number>(b0, b1, b2);
                a.fill(Number.NaN);
                for (let i2 = 0; i2 < b2; ++i2) {
                    for (let i1 = 0; i1 < b1; ++i1) {
                        for (let i0 = 0; i0 < b0; ++i0) {
                            a.set(i0, i1, i2, vals[i0 + b0 * (i1 + b1 * i2)]!);
                        }
                    }
                }
                for (let i2 = 0; i2 < b2; ++i2) {
                    for (let i1 = 0; i1 < b1; ++i1) {
                        for (let i0 = 0; i0 < b0; ++i0) {
                            if (!Object.is(a.get(i0, i1, i2), vals[i0 + b0 * (i1 + b1 * i2)])) {
                                return false;
                            }
                        }
                    }
                }
                return true;
            });
    });

    it('a slice of an Array3 matches the corresponding Array2', () => {
        // Fixing i2 gives a contiguous bound0 x bound1 block, which upstream's
        // pointer indirection also exposes as a plain 2D array.
        check(fc.tuple(fc.integer({ min: 1, max: 5 }), fc.integer({ min: 1, max: 5 }),
            fc.integer({ min: 1, max: 5 })),
            ([b0, b1, b2]) => {
                const a = new Array3<number>(b0, b1, b2);
                a.fill(0);
                for (let i2 = 0; i2 < b2; ++i2) {
                    for (let i1 = 0; i1 < b1; ++i1) {
                        for (let i0 = 0; i0 < b0; ++i0) {
                            a.set(i0, i1, i2, 100 * i2 + 10 * i1 + i0);
                        }
                    }
                }
                const flat = a.data();
                for (let i2 = 0; i2 < b2; ++i2) {
                    const slice = flat.slice(b0 * b1 * i2, b0 * b1 * (i2 + 1));
                    const two = new Array2<number>(b0, b1, slice);
                    for (let i1 = 0; i1 < b1; ++i1) {
                        for (let i0 = 0; i0 < b0; ++i0) {
                            if (two.get(i0, i1) !== a.get(i0, i1, i2)) { return false; }
                        }
                    }
                }
                return true;
            });
    });

    it('caller-owned storage stays aliased in both directions', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }), fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 }), finite()),
            ([b0, b1, b2, v]) => {
                const objects = new Array<number>(b0 * b1 * b2).fill(0);
                const a = new Array3<number>(b0, b1, b2, objects);
                const i0 = b0 - 1, i1 = b1 - 1, i2 = b2 - 1;
                a.set(i0, i1, i2, v);
                if (!Object.is(objects[i0 + b0 * (i1 + b1 * i2)], v)) { return false; }
                objects[0] = v + 1;
                return Object.is(a.get(0, 0, 0), v + 1);
            });
    });

    it('fill covers exactly bound0 * bound1 * bound2 entries', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }), fc.integer({ min: 0, max: 5 }),
            fc.integer({ min: 0, max: 5 })),
            ([b0, b1, b2]) => {
                const a = new Array3<number>(b0, b1, b2);
                a.fill(-2);
                const flat = a.data();
                if (flat.length !== b0 * b1 * b2) { return false; }
                return flat.every(x => x === -2);
            });
    });
});
