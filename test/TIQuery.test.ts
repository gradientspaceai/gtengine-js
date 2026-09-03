import { describe, it, expect } from 'vitest';
import { check, finite, positive, seededRandom, fc } from './helpers/arbitraries.js';
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

describe('TIQuery verification', () => {
    const query: TIQuery<Interval, Interval, Result> = new IntrIntervalIntervalTI();
    const interval = fc.tuple(finite(), finite())
        .map(([a, b]): Interval => ({ min: Math.min(a, b), max: Math.max(a, b) }));

    it('test() is symmetric under argument swap', () => {
        check(fc.tuple(interval, interval), ([i0, i1]) =>
            query.test(i0, i1).intersect === query.test(i1, i0).intersect);
    });

    it('test() agrees with brute-force sampling of the two intervals', () => {
        const rand = seededRandom(0x71D3);
        check(fc.tuple(interval, interval), ([i0, i1]) => {
            const intersect = query.test(i0, i1).intersect;
            const lo = Math.min(i0.min, i1.min), hi = Math.max(i0.max, i1.max);
            for (let k = 0; k < 64; ++k) {
                const t = lo + (hi - lo) * rand();
                if (t >= i0.min && t <= i0.max && t >= i1.min && t <= i1.max && !intersect) {
                    return false;
                }
            }
            return true;
        });
    });

    it('an interval always intersects itself', () => {
        check(interval, i => query.test(i, i).intersect);
    });

    it('disjoint intervals never intersect', () => {
        check(fc.tuple(interval, positive()), ([i, gap]) => {
            const shifted: Interval = { min: i.max + gap, max: i.max + gap + 1 };
            return !query.test(i, shifted).intersect;
        });
    });
});
