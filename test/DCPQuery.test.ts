import { describe, it, expect } from 'vitest';
import { type DCPQuery } from '../src/DCPQuery.js';

// A miniature distance query in the style of the Dist* files: point-to-point
// distance in 1D, with a Result shaped like upstream's nested structs.
interface Result {
    distance: number;
    sqrDistance: number;
    closest: [number, number];
}

class DistPointPoint1 implements DCPQuery<number, number, Result> {
    compute(primitive0: number, primitive1: number): Result {
        const diff = primitive0 - primitive1;
        return {
            distance: Math.abs(diff),
            sqrDistance: diff * diff,
            closest: [primitive0, primitive1]
        };
    }
}

describe('DCPQuery', () => {
    it('is implementable by concrete distance queries with compute()', () => {
        const query = new DistPointPoint1();
        const result = query.compute(1, 4);
        expect(result.distance).toBe(3);
        expect(result.sqrDistance).toBe(9);
        expect(result.closest).toEqual([1, 4]);
    });

    it('supports use through the interface type', () => {
        const query: DCPQuery<number, number, Result> = new DistPointPoint1();
        expect(query.compute(-2, 2).distance).toBe(4);
        expect(query.compute(5, 5).distance).toBe(0);
    });
});
