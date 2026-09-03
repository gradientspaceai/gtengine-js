import { describe, it, expect } from 'vitest';
import { check, finite, expectClose, fc } from './helpers/arbitraries.js';
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

describe('DCPQuery verification', () => {
    // The interface carries no behaviour of its own (upstream's primary
    // template is empty), so the properties exercise the contract every Dist*
    // implementation must satisfy through it: compute() returns a result with
    // the distance invariants and never mutates its inputs.
    const query: DCPQuery<number, number, Result> = new DistPointPoint1();

    it('distance equals |closest[0] - closest[1]| and sqrDistance its square', () => {
        check(fc.tuple(finite(-1e3, 1e3), finite(-1e3, 1e3)), ([a, b]) => {
            const r = query.compute(a, b);
            expectClose(r.distance, Math.abs(r.closest[0] - r.closest[1]), 1e-12, 1e-12);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-9, 1e-12);
        });
    });

    it('the query is symmetric under argument swap', () => {
        check(fc.tuple(finite(-1e3, 1e3), finite(-1e3, 1e3)), ([a, b]) => {
            const r0 = query.compute(a, b);
            const r1 = query.compute(b, a);
            expect(r1.distance).toBe(r0.distance);
            expect(r1.closest).toEqual([r0.closest[1], r0.closest[0]]);
        });
    });

    it('successive results do not alias each other', () => {
        check(fc.tuple(finite(), finite(), finite(), finite()), ([a, b, c, d]) => {
            const r0 = query.compute(a, b);
            const snapshot = [...r0.closest];
            const r1 = query.compute(c, d);
            expect(r0.closest).toEqual(snapshot);
            expect(r1).not.toBe(r0);
        });
    });

    it('distance is zero exactly for coincident inputs', () => {
        check(finite(), a => {
            expect(query.compute(a, a).distance).toBe(0);
            expect(query.compute(a, a).sqrDistance).toBe(0);
        });
    });
});
