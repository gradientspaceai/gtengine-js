import { describe, it, expect } from 'vitest';
import { check, finite, seededRandom, fc } from './helpers/arbitraries.js';
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

describe('FIQuery verification', () => {
    const query: FIQuery<Interval, Interval, Result> = new IntrIntervalIntervalFI();
    const interval = fc.tuple(finite(), finite())
        .map(([a, b]): Interval => ({ min: Math.min(a, b), max: Math.max(a, b) }));

    it('find() is symmetric under argument swap', () => {
        check(fc.tuple(interval, interval), ([i0, i1]) => {
            const r0 = query.find(i0, i1);
            const r1 = query.find(i1, i0);
            expect(r1.intersect).toBe(r0.intersect);
            expect(r1.overlap).toEqual(r0.overlap);
        });
    });

    it('the overlap is contained in both inputs and contains their common points', () => {
        check(fc.tuple(interval, interval), ([i0, i1]) => {
            const r = query.find(i0, i1);
            if (!r.intersect) {
                // No point can lie in both intervals.
                expect(r.overlap).toBeNull();
                return Math.max(i0.min, i1.min) > Math.min(i0.max, i1.max);
            }
            const o = r.overlap!;
            return o.min >= i0.min && o.min >= i1.min
                && o.max <= i0.max && o.max <= i1.max
                && o.min <= o.max;
        });
    });

    it('brute-force sampling agrees with the reported overlap', () => {
        const rand = seededRandom(0x51D3);
        check(fc.tuple(interval, interval), ([i0, i1]) => {
            const r = query.find(i0, i1);
            const lo = Math.min(i0.min, i1.min), hi = Math.max(i0.max, i1.max);
            for (let k = 0; k < 64; ++k) {
                const t = lo + (hi - lo) * rand();
                const inBoth = t >= i0.min && t <= i0.max && t >= i1.min && t <= i1.max;
                if (inBoth && !r.intersect) { return false; }
                if (inBoth && r.overlap !== null
                    && (t < r.overlap.min || t > r.overlap.max)) { return false; }
            }
            return true;
        });
    });

    it('the unset output is null, never a stale object from a previous call', () => {
        const r0 = query.find({ min: 0, max: 1 }, { min: 0.5, max: 2 });
        expect(r0.overlap).toEqual({ min: 0.5, max: 1 });
        const r1 = query.find({ min: 0, max: 1 }, { min: 5, max: 6 });
        expect(r1.overlap).toBeNull();
        expect(r0.overlap).toEqual({ min: 0.5, max: 1 });
    });
});
