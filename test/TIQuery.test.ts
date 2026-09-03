import { describe, it, expect } from 'vitest';
import { type TIQuery } from '../src/TIQuery.js';

// A miniature test-intersection query in the style of the Intr* files:
// interval-interval overlap in 1D.
interface Interval {
    min: number;
    max: number;
}

interface Result {
    intersect: boolean;
}

class IntrIntervalIntervalTI implements TIQuery<Interval, Interval, Result> {
    test(primitive0: Interval, primitive1: Interval): Result {
        return {
            intersect: primitive0.min <= primitive1.max && primitive1.min <= primitive0.max
        };
    }
}

describe('TIQuery', () => {
    it('is implementable by concrete intersection queries with test()', () => {
        const query = new IntrIntervalIntervalTI();
        expect(query.test({ min: 0, max: 2 }, { min: 1, max: 3 }).intersect).toBe(true);
        expect(query.test({ min: 0, max: 2 }, { min: 3, max: 4 }).intersect).toBe(false);
        // Touching endpoints intersect.
        expect(query.test({ min: 0, max: 2 }, { min: 2, max: 4 }).intersect).toBe(true);
    });

    it('supports use through the interface type', () => {
        const query: TIQuery<Interval, Interval, Result> = new IntrIntervalIntervalTI();
        expect(query.test({ min: -1, max: 1 }, { min: 0, max: 0 }).intersect).toBe(true);
    });
});
