import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { Array3 } from '../src/Array3.js';
import { Array4 } from '../src/Array4.js';

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

describe('Array4 verification', () => {
    it('flat index equals i0 + b0 * (i1 + b1 * (i2 + b2 * i3))', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 4 }), fc.integer({ min: 1, max: 4 }),
            fc.integer({ min: 1, max: 4 }), fc.integer({ min: 1, max: 4 })),
            ([b0, b1, b2, b3]) => {
                const a = new Array4<number>(b0, b1, b2, b3);
                a.fill(0);
                for (let i3 = 0; i3 < b3; ++i3) {
                    for (let i2 = 0; i2 < b2; ++i2) {
                        for (let i1 = 0; i1 < b1; ++i1) {
                            for (let i0 = 0; i0 < b0; ++i0) {
                                a.set(i0, i1, i2, i3, i0 + b0 * (i1 + b1 * (i2 + b2 * i3)));
                            }
                        }
                    }
                }
                const flat = a.data();
                for (let k = 0; k < b0 * b1 * b2 * b3; ++k) {
                    if (flat[k] !== k) { return false; }
                }
                return true;
            });
    });

    it('get/set round trips over the whole index space', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }),
            fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }),
            fc.array(finite(), { minLength: 81, maxLength: 81 })),
            ([b0, b1, b2, b3, vals]) => {
                const a = new Array4<number>(b0, b1, b2, b3);
                a.fill(Number.NaN);
                const idx = (i0: number, i1: number, i2: number, i3: number) =>
                    i0 + b0 * (i1 + b1 * (i2 + b2 * i3));
                for (let i3 = 0; i3 < b3; ++i3) {
                    for (let i2 = 0; i2 < b2; ++i2) {
                        for (let i1 = 0; i1 < b1; ++i1) {
                            for (let i0 = 0; i0 < b0; ++i0) {
                                a.set(i0, i1, i2, i3, vals[idx(i0, i1, i2, i3)]!);
                            }
                        }
                    }
                }
                for (let i3 = 0; i3 < b3; ++i3) {
                    for (let i2 = 0; i2 < b2; ++i2) {
                        for (let i1 = 0; i1 < b1; ++i1) {
                            for (let i0 = 0; i0 < b0; ++i0) {
                                if (!Object.is(a.get(i0, i1, i2, i3), vals[idx(i0, i1, i2, i3)])) {
                                    return false;
                                }
                            }
                        }
                    }
                }
                return true;
            });
    });

    it('a cuboid of an Array4 matches the corresponding Array3', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }),
            fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 })),
            ([b0, b1, b2, b3]) => {
                const a = new Array4<number>(b0, b1, b2, b3);
                a.fill(0);
                for (let i3 = 0; i3 < b3; ++i3) {
                    for (let i2 = 0; i2 < b2; ++i2) {
                        for (let i1 = 0; i1 < b1; ++i1) {
                            for (let i0 = 0; i0 < b0; ++i0) {
                                a.set(i0, i1, i2, i3, 1000 * i3 + 100 * i2 + 10 * i1 + i0);
                            }
                        }
                    }
                }
                const flat = a.data();
                const stride = b0 * b1 * b2;
                for (let i3 = 0; i3 < b3; ++i3) {
                    const cuboid = new Array3<number>(b0, b1, b2,
                        flat.slice(stride * i3, stride * (i3 + 1)));
                    for (let i2 = 0; i2 < b2; ++i2) {
                        for (let i1 = 0; i1 < b1; ++i1) {
                            for (let i0 = 0; i0 < b0; ++i0) {
                                if (cuboid.get(i0, i1, i2) !== a.get(i0, i1, i2, i3)) {
                                    return false;
                                }
                            }
                        }
                    }
                }
                return true;
            });
    });

    it('caller-owned storage stays aliased in both directions', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }),
            fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }), finite()),
            ([b0, b1, b2, b3, v]) => {
                const objects = new Array<number>(b0 * b1 * b2 * b3).fill(0);
                const a = new Array4<number>(b0, b1, b2, b3, objects);
                const i0 = b0 - 1, i1 = b1 - 1, i2 = b2 - 1, i3 = b3 - 1;
                a.set(i0, i1, i2, i3, v);
                if (!Object.is(objects[i0 + b0 * (i1 + b1 * (i2 + b2 * i3))], v)) { return false; }
                objects[0] = v + 1;
                return Object.is(a.get(0, 0, 0, 0), v + 1);
            });
    });

    it('fill covers exactly the product of the four bounds', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 4 }), fc.integer({ min: 0, max: 4 }),
            fc.integer({ min: 0, max: 4 }), fc.integer({ min: 0, max: 4 })),
            ([b0, b1, b2, b3]) => {
                const a = new Array4<number>(b0, b1, b2, b3);
                a.fill(9);
                const flat = a.data();
                if (flat.length !== b0 * b1 * b2 * b3) { return false; }
                return flat.every(x => x === 9);
            });
    });
});
