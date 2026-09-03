import { describe, it, expect } from 'vitest';
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
