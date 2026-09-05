import { describe, it, expect } from 'vitest';
import { Histogram } from '../src/Histogram.js';
import { check, fc, finite, scaled } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream Histogram.h semantics.
// ---------------------------------------------------------------------------

describe('Histogram verification', () => {
    const numBuckets = fc.integer({ min: 1, max: 24 });
    const intSamples = fc.array(fc.integer({ min: -30, max: 60 }),
        { minLength: 1, maxLength: 60 });
    const realSamples = fc.array(finite(-50, 50), { minLength: 1, maxLength: 60 });

    const sum = (a: readonly number[]) => a.reduce((p, q) => p + q, 0);

    it('noRescaling: every sample is counted exactly once', () => {
        check(fc.tuple(numBuckets, intSamples), ([b, samples]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, true);
            const buckets = histogram.getBuckets();
            expect(buckets.length).toBe(b);
            expect(sum(buckets) + histogram.getExcessLess() + histogram.getExcessGreater())
                .toBe(samples.length);

            // Brute-force bucket assignment.
            const brute = new Array<number>(b).fill(0);
            let less = 0, greater = 0;
            for (const s of samples) {
                if (s < 0) { ++less; }
                else if (s >= b) { ++greater; }
                else { ++brute[s]; }
            }
            expect([...buckets]).toEqual(brute);
            expect(histogram.getExcessLess()).toBe(less);
            expect(histogram.getExcessGreater()).toBe(greater);
        });
    });

    it('rescaling integer samples: counts sum to the sample count and match the formula', () => {
        check(fc.tuple(numBuckets, intSamples), ([b, samples]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, false);
            const buckets = histogram.getBuckets();
            expect(sum(buckets)).toBe(samples.length);
            // Rescaling never produces excess counts.
            expect(histogram.getExcessLess()).toBe(0);
            expect(histogram.getExcessGreater()).toBe(0);

            const minValue = Math.min(...samples);
            const maxValue = Math.max(...samples);
            const brute = new Array<number>(b).fill(0);
            if (minValue < maxValue) {
                const mult = (b - 1) / (maxValue - minValue);
                for (const s of samples) {
                    ++brute[Math.trunc(mult * (s - minValue))];
                }
            } else {
                brute[0] = samples.length;
            }
            expect([...buckets]).toEqual(brute);
        });
    });

    it('rescaling real samples: bin indices stay in range and match the formula', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 24 }), realSamples), ([b, samples]) => {
            const histogram = Histogram.fromRealSamples(b, samples);
            const buckets = histogram.getBuckets();
            expect(sum(buckets)).toBe(samples.length);
            expect(buckets.length).toBe(b);

            const minValue = Math.min(...samples);
            const maxValue = Math.max(...samples);
            const brute = new Array<number>(b).fill(0);
            if (minValue < maxValue) {
                const mult = (b - 1) / (maxValue - minValue);
                for (const s of samples) {
                    const index = Math.trunc(mult * (s - minValue));
                    // No sample may fall outside the bucket array; upstream
                    // does not bounds check the rescaled index.
                    expect(index).toBeGreaterThanOrEqual(0);
                    expect(index).toBeLessThan(b);
                    ++brute[index];
                }
                // The minimum always maps to bucket 0. The maximum maps to
                // bucket b - 1 only when ((b-1)/d)*d rounds up to b-1; see
                // the upstream rounding quirk noted in the PR.
                expect(Math.trunc(mult * (minValue - minValue))).toBe(0);
            } else {
                brute[0] = samples.length;
            }
            expect([...buckets]).toEqual(brute);
        });
    });

    it('rescaling assigns nondecreasing bin indices to sorted samples', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 24 }),
            fc.array(scaled(-50, 50), { minLength: 2, maxLength: 40 })),
            ([b, samples]) => {
                const sorted = [...samples].sort((p, q) => p - q);
                const minValue = sorted[0], maxValue = sorted[sorted.length - 1];
                if (minValue >= maxValue) { return; }
                const mult = (b - 1) / (maxValue - minValue);
                let previous = -1;
                for (const s of sorted) {
                    const index = Math.trunc(mult * (s - minValue));
                    expect(index).toBeGreaterThanOrEqual(previous);
                    previous = index;
                }
                // The histogram is the tally of exactly those indices.
                const histogram = Histogram.fromRealSamples(b, samples);
                expect(sum(histogram.getBuckets())).toBe(samples.length);
            });
    });

    it('insert and insertCheck agree for in-range values', () => {
        check(fc.tuple(numBuckets, fc.array(fc.integer({ min: -5, max: 30 }),
            { minLength: 1, maxLength: 40 })), ([b, values]) => {
                const checked = new Histogram(b);
                for (const v of values) {
                    checked.insertCheck(v);
                }
                const unchecked = new Histogram(b);
                for (const v of values) {
                    if (0 <= v && v < b) {
                        unchecked.insert(v);
                    }
                }
                expect([...checked.getBuckets()]).toEqual([...unchecked.getBuckets()]);
                expect(sum(checked.getBuckets()) + checked.getExcessLess()
                    + checked.getExcessGreater()).toBe(values.length);
            });
    });

    it('getLowerTail returns the first index whose cdf reaches the tail count', () => {
        check(fc.tuple(numBuckets, intSamples, finite(0, 1)), ([b, samples, tail]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, true);
            const buckets = histogram.getBuckets();
            const hSum = sum(buckets);
            const hTailSum = Math.trunc(tail * hSum);
            const lower = histogram.getLowerTail(tail);

            // Independent scan of the cumulative distribution.
            let cdf = 0, expected = b;
            for (let i = 0; i < b; ++i) {
                cdf += buckets[i];
                if (cdf >= hTailSum) { expected = i; break; }
            }
            expect(lower).toBe(expected);
            // Upstream returns numBuckets when the loop never breaks, which
            // happens only if every partial sum stays below the tail count.
            expect(lower <= b).toBe(true);
        });
    });

    it('getUpperTail scans from the last bucket downward', () => {
        check(fc.tuple(numBuckets, intSamples, finite(0, 1)), ([b, samples, tail]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, true);
            const buckets = histogram.getBuckets();
            const hTailSum = Math.trunc(tail * sum(buckets));
            const upper = histogram.getUpperTail(tail);

            let cdf = 0, expected = -1;
            for (let i = b - 1; i >= 0; --i) {
                cdf += buckets[i];
                if (cdf >= hTailSum) { expected = i; break; }
            }
            expect(upper).toBe(expected);
        });
    });

    it('getTails equals the single-sided queries at half the tail amount', () => {
        check(fc.tuple(numBuckets, intSamples, finite(0, 1)), ([b, samples, tail]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, true);
            const { lower, upper } = histogram.getTails(tail);
            // getTails uses trunc(0.5 * tail * hSum) for both scans; the
            // single-sided queries use trunc(tail * hSum), so halving the
            // input reproduces them exactly only when the products agree.
            const hSum = sum(histogram.getBuckets());
            if (Math.trunc(0.5 * tail * hSum) === Math.trunc((0.5 * tail) * hSum)) {
                expect(lower).toBe(histogram.getLowerTail(0.5 * tail));
                expect(upper).toBe(histogram.getUpperTail(0.5 * tail));
            }
        });
    });

    it('a zero tail amount selects bucket 0 and the last bucket', () => {
        check(fc.tuple(numBuckets, intSamples), ([b, samples]) => {
            const histogram = Histogram.fromIntegerSamples(b, samples, true);
            // hTailSum is 0, so both loops break on their first iteration.
            expect(histogram.getLowerTail(0)).toBe(0);
            expect(histogram.getUpperTail(0)).toBe(b - 1);
            const { lower, upper } = histogram.getTails(0);
            expect(lower).toBe(0);
            expect(upper).toBe(b - 1);
        });
    });
});
