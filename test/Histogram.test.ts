import { describe, it, expect } from 'vitest';
import { Histogram } from '../src/Histogram';

describe('Histogram', () => {
    it('constructs empty for incremental updates', () => {
        const histogram = new Histogram(4);
        expect(histogram.getBuckets()).toEqual([0, 0, 0, 0]);
        expect(histogram.getExcessLess()).toBe(0);
        expect(histogram.getExcessGreater()).toBe(0);
    });

    it('throws on nonpositive bucket counts and empty samples', () => {
        expect(() => new Histogram(0)).toThrow('Invalid input.');
        expect(() => Histogram.fromIntegerSamples(4, [], true)).toThrow('Invalid input.');
        expect(() => Histogram.fromRealSamples(0, [1])).toThrow('Invalid input.');
    });

    it('insert increments buckets without bounds checking', () => {
        const histogram = new Histogram(3);
        histogram.insert(0);
        histogram.insert(2);
        histogram.insert(2);
        expect(histogram.getBuckets()).toEqual([1, 0, 2]);
    });

    it('insertCheck counts out-of-range values as excess', () => {
        const histogram = new Histogram(3);
        histogram.insertCheck(-1);
        histogram.insertCheck(0);
        histogram.insertCheck(2);
        histogram.insertCheck(3);
        histogram.insertCheck(100);
        expect(histogram.getBuckets()).toEqual([1, 0, 1]);
        expect(histogram.getExcessLess()).toBe(1);
        expect(histogram.getExcessGreater()).toBe(2);
    });

    it('maps integer samples directly with noRescaling', () => {
        const samples = [0, 1, 1, 3, 3, 3, -2, 4, 7];
        const histogram = Histogram.fromIntegerSamples(4, samples, true);
        expect(histogram.getBuckets()).toEqual([1, 2, 0, 3]);
        expect(histogram.getExcessLess()).toBe(1);      // -2
        expect(histogram.getExcessGreater()).toBe(2);   // 4, 7
    });

    it('rescales integer samples so the extremes map to the end buckets', () => {
        // min = 10 maps to bucket 0, max = 50 maps to bucket 4 (numBuckets-1).
        const samples = [10, 50, 30, 20, 40];
        const histogram = Histogram.fromIntegerSamples(5, samples, false);
        // index = trunc((5-1)/(50-10) * (v-10)) = trunc((v-10)/10)
        expect(histogram.getBuckets()).toEqual([1, 1, 1, 1, 1]);
        expect(histogram.getExcessLess()).toBe(0);
        expect(histogram.getExcessGreater()).toBe(0);
    });

    it('places constant integer samples in bucket 0', () => {
        const histogram = Histogram.fromIntegerSamples(8, [5, 5, 5], false);
        expect(histogram.getBuckets()).toEqual([3, 0, 0, 0, 0, 0, 0, 0]);
    });

    it('rescales real samples with correct bin placement and boundaries', () => {
        // min = 0 -> bucket 0, max = 1 -> bucket 3; mult = 3.
        const samples = [0, 0.1, 0.34, 0.5, 0.67, 0.99, 1];
        const histogram = Histogram.fromRealSamples(4, samples);
        // indices: trunc(3*v) = 0, 0, 1, 1, 2, 2, 3
        expect(histogram.getBuckets()).toEqual([2, 2, 2, 1]);
    });

    it('handles boundary values that land exactly on bin edges', () => {
        // mult = (4-1)/(4-0) = 0.75; values 0,1,2,3,4 -> trunc(0.75*v) =
        // 0, 0, 1, 2, 3.
        const histogram = Histogram.fromRealSamples(4, [0, 1, 2, 3, 4]);
        expect(histogram.getBuckets()).toEqual([2, 1, 1, 1]);
    });

    it('places constant real samples in bucket 0', () => {
        const histogram = Histogram.fromRealSamples(4, [2.5, 2.5]);
        expect(histogram.getBuckets()).toEqual([2, 0, 0, 0]);
    });

    it('computes the lower tail index', () => {
        const histogram = new Histogram(5);
        // buckets = [10, 20, 40, 20, 10], N = 100
        const counts = [10, 20, 40, 20, 10];
        for (let b = 0; b < 5; ++b) {
            for (let i = 0; i < counts[b]; ++i) {
                histogram.insert(b);
            }
        }

        // tailAmount 0.05 -> hTailSum = 5: cdf(0) = 10 >= 5 -> L = 0.
        expect(histogram.getLowerTail(0.05)).toBe(0);
        // tailAmount 0.25 -> hTailSum = 25: cdf(0) = 10 < 25,
        // cdf(1) = 30 >= 25 -> L = 1.
        expect(histogram.getLowerTail(0.25)).toBe(1);
        // Verify the documented cdf properties for L = getLowerTail(t).
        const cdf = [10, 30, 70, 90, 100];
        for (const t of [0.05, 0.1, 0.25, 0.5, 0.9]) {
            const lower = histogram.getLowerTail(t);
            const hTailSum = Math.trunc(t * 100);
            expect(cdf[lower]).toBeGreaterThanOrEqual(hTailSum);
            if (lower > 0) {
                expect(cdf[lower - 1]).toBeLessThan(hTailSum);
            }
        }
    });

    it('computes the upper tail index', () => {
        const histogram = new Histogram(5);
        const counts = [10, 20, 40, 20, 10];
        for (let b = 0; b < 5; ++b) {
            for (let i = 0; i < counts[b]; ++i) {
                histogram.insert(b);
            }
        }

        // tailAmount 0.05 -> hTailSum = 5: bucket 4 = 10 >= 5 -> U = 4.
        expect(histogram.getUpperTail(0.05)).toBe(4);
        // tailAmount 0.25 -> hTailSum = 25: bucket 4 = 10 < 25,
        // buckets 3+4 = 30 >= 25 -> U = 3.
        expect(histogram.getUpperTail(0.25)).toBe(3);
    });

    it('computes both tails with half the tail amount on each side', () => {
        const histogram = new Histogram(5);
        const counts = [10, 20, 40, 20, 10];
        for (let b = 0; b < 5; ++b) {
            for (let i = 0; i < counts[b]; ++i) {
                histogram.insert(b);
            }
        }

        // tailAmount 0.5 -> hTailSum = 25 on each side: L = 1, U = 3.
        const { lower, upper } = histogram.getTails(0.5);
        expect(lower).toBe(1);
        expect(upper).toBe(3);
        expect(lower).toBe(histogram.getLowerTail(0.25));
        expect(upper).toBe(histogram.getUpperTail(0.25));
    });

    it('cumulative tail queries are consistent on a skewed histogram', () => {
        const histogram = Histogram.fromIntegerSamples(
            4, [0, 0, 0, 0, 0, 0, 0, 0, 1, 3], true);
        // buckets = [8, 1, 0, 1], N = 10.
        expect(histogram.getBuckets()).toEqual([8, 1, 0, 1]);
        // tailAmount 0.5 -> hTailSum = 5: cdf(0) = 8 >= 5 -> L = 0; from
        // above, 1 + 0 + 1 = 2 < 5, adding bucket 0 gives 10 >= 5 -> U = 0.
        expect(histogram.getLowerTail(0.5)).toBe(0);
        expect(histogram.getUpperTail(0.5)).toBe(0);
        // A smaller tail amount stays in the true upper tail:
        // hTailSum = 1 -> bucket 3 alone suffices -> U = 3.
        expect(histogram.getUpperTail(0.1)).toBe(3);
    });
});
