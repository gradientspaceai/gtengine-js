import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { LexicoArray2 } from '../src/LexicoArray2.js';

describe('LexicoArray2', () => {
    // A 2x3 matrix
    //   [ 1 2 3 ]
    //   [ 4 5 6 ]
    // in row-major flat storage.
    const rowMajorData = [1, 2, 3, 4, 5, 6];
    // The same matrix in column-major flat storage.
    const colMajorData = [1, 4, 2, 5, 3, 6];

    it('reports its shape and ordering', () => {
        const a = new LexicoArray2(true, 2, 3, [...rowMajorData]);
        expect(a.isRowMajor()).toBe(true);
        expect(a.getNumRows()).toBe(2);
        expect(a.getNumCols()).toBe(3);
    });

    it('row-major get(r, c) reads element c + numCols * r', () => {
        const a = new LexicoArray2(true, 2, 3, [...rowMajorData]);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(a.get(r, c)).toBe(1 + c + 3 * r);
            }
        }
    });

    it('column-major get(r, c) reads element r + numRows * c', () => {
        const a = new LexicoArray2(false, 2, 3, [...colMajorData]);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(a.get(r, c)).toBe(1 + c + 3 * r);
            }
        }
    });

    it('both orderings agree on the logical matrix', () => {
        const rm = new LexicoArray2(true, 2, 3, [...rowMajorData]);
        const cm = new LexicoArray2(false, 2, 3, [...colMajorData]);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(rm.get(r, c)).toBe(cm.get(r, c));
            }
        }
    });

    it('set writes through to the caller-owned flat array (aliasing)', () => {
        const backing = [0, 0, 0, 0, 0, 0];
        const a = new LexicoArray2(true, 2, 3, backing);
        a.set(1, 2, 42);
        expect(backing[2 + 3 * 1]).toBe(42);
        expect(a.get(1, 2)).toBe(42);

        const b = new LexicoArray2(false, 2, 3, backing);
        b.set(0, 1, 7);
        expect(backing[0 + 2 * 1]).toBe(7);
    });
});

describe('LexicoArray2 verification', () => {
    const shape = fc.tuple(fc.integer({ min: 1, max: 7 }), fc.integer({ min: 1, max: 7 }));

    it('row-major index is c + numCols * r for random shapes', () => {
        check(shape, ([rows, cols]) => {
            const data = new Array<number>(rows * cols).fill(0);
            const a = new LexicoArray2(true, rows, cols, data);
            for (let r = 0; r < rows; ++r) {
                for (let c = 0; c < cols; ++c) { a.set(r, c, c + cols * r); }
            }
            for (let k = 0; k < rows * cols; ++k) {
                if (data[k] !== k) { return false; }
            }
            return true;
        });
    });

    it('column-major index is r + numRows * c for random shapes', () => {
        check(shape, ([rows, cols]) => {
            const data = new Array<number>(rows * cols).fill(0);
            const a = new LexicoArray2(false, rows, cols, data);
            for (let r = 0; r < rows; ++r) {
                for (let c = 0; c < cols; ++c) { a.set(r, c, r + rows * c); }
            }
            for (let k = 0; k < rows * cols; ++k) {
                if (data[k] !== k) { return false; }
            }
            return true;
        });
    });

    it('the two orderings describe the same logical matrix', () => {
        check(fc.tuple(shape, fc.array(finite(), { minLength: 49, maxLength: 49 })),
            ([[rows, cols], vals]) => {
                const rmData = new Array<number>(rows * cols).fill(0);
                const cmData = new Array<number>(rows * cols).fill(0);
                const rm = new LexicoArray2(true, rows, cols, rmData);
                const cm = new LexicoArray2(false, rows, cols, cmData);
                for (let r = 0; r < rows; ++r) {
                    for (let c = 0; c < cols; ++c) {
                        const v = vals[c + cols * r]!;
                        rm.set(r, c, v);
                        cm.set(r, c, v);
                    }
                }
                // Same logical entries, different flat layouts.
                for (let r = 0; r < rows; ++r) {
                    for (let c = 0; c < cols; ++c) {
                        if (!Object.is(rm.get(r, c), cm.get(r, c))) { return false; }
                    }
                }
                // Column-major storage of A equals row-major storage of A^T.
                for (let r = 0; r < rows; ++r) {
                    for (let c = 0; c < cols; ++c) {
                        if (!Object.is(cmData[r + rows * c], rmData[c + cols * r])) {
                            return false;
                        }
                    }
                }
                return true;
            });
    });

    it('get/set round trip and the backing array is aliased, not copied', () => {
        check(fc.tuple(shape, fc.boolean(), finite()), ([[rows, cols], rowMajor, v]) => {
            const data = new Array<number>(rows * cols).fill(0);
            const a = new LexicoArray2(rowMajor, rows, cols, data);
            const r = rows - 1, c = cols - 1;
            a.set(r, c, v);
            const flat = rowMajor ? c + cols * r : r + rows * c;
            if (!Object.is(data[flat], v)) { return false; }
            data[0] = v + 1;
            return Object.is(a.get(0, 0), v + 1);
        });
    });

    it('shape accessors report the constructor arguments unchanged', () => {
        check(fc.tuple(shape, fc.boolean()), ([[rows, cols], rowMajor]) => {
            const a = new LexicoArray2(rowMajor, rows, cols,
                new Array<number>(rows * cols).fill(0));
            return a.getNumRows() === rows && a.getNumCols() === cols
                && a.isRowMajor() === rowMajor;
        });
    });

    it('every (r, c) maps to a distinct flat index (no aliasing of cells)', () => {
        check(fc.tuple(shape, fc.boolean()), ([[rows, cols], rowMajor]) => {
            const data = new Array<number>(rows * cols).fill(-1);
            const a = new LexicoArray2(rowMajor, rows, cols, data);
            let k = 0;
            for (let r = 0; r < rows; ++r) {
                for (let c = 0; c < cols; ++c) { a.set(r, c, k++); }
            }
            return new Set(data).size === rows * cols;
        });
    });
});
