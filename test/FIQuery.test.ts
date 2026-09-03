import { describe, it, expect } from 'vitest';
import { type FIQuery } from '../src/FIQuery.js';

// A miniature find-intersection query in the style of the Intr* files:
// interval-interval intersection in 1D.
interface Interval {
    min: number;
    max: number;
}

interface Result {
    intersect: boolean;
    overlap: Interval | null;
}

class IntrIntervalIntervalFI implements FIQuery<Interval, Interval, Result> {
    find(primitive0: Interval, primitive1: Interval): Result {
        const min = Math.max(primitive0.min, primitive1.min);
        const max = Math.min(primitive0.max, primitive1.max);
        if (min <= max) {
            return { intersect: true, overlap: { min, max } };
        }
        return { intersect: false, overlap: null };
    }
}

describe('FIQuery', () => {
    it('is implementable by concrete intersection queries with find()', () => {
        const query = new IntrIntervalIntervalFI();
        const result = query.find({ min: 0, max: 2 }, { min: 1, max: 3 });
        expect(result.intersect).toBe(true);
        expect(result.overlap).toEqual({ min: 1, max: 2 });
    });

    it('reports no intersection for disjoint inputs', () => {
        const query: FIQuery<Interval, Interval, Result> = new IntrIntervalIntervalFI();
        const result = query.find({ min: 0, max: 1 }, { min: 2, max: 3 });
        expect(result.intersect).toBe(false);
        expect(result.overlap).toBeNull();
    });
});
