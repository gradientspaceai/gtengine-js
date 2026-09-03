// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Histogram.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: the four upstream constructors become the constructor for the
// incremental-update mode plus two static factories, following the ported
// convention that ambiguous C++ constructor overloads become factories:
//   Histogram(numBuckets, numSamples, int32_t const*, noRescaling)
//     -> Histogram.fromIntegerSamples(numBuckets, samples, noRescaling)
//   Histogram(numBuckets, numSamples, float const*) and
//   Histogram(numBuckets, numSamples, double const*)
//     -> Histogram.fromRealSamples(numBuckets, samples)
//   Histogram(numBuckets) -> new Histogram(numBuckets)
// The float and double sample constructors have identical logic once Real is
// number, so they are one factory. The (numSamples, pointer) pairs become
// arrays whose length is the sample count. The out-parameters of GetTails
// become the return value { lower, upper }.

import { logAssert } from './Logger.js';

export class Histogram {
    private mBuckets: number[];
    private mExcessLess: number;
    private mExcessGreater: number;

    // Construction for the case when you plan on updating the histogram
    // incrementally. The incremental update is implemented only for integer
    // samples and no rescaling; see insert() and insertCheck().
    constructor(numBuckets: number) {
        logAssert(numBuckets > 0, 'Invalid input.');
        this.mBuckets = new Array<number>(numBuckets).fill(0);
        this.mExcessLess = 0;
        this.mExcessGreater = 0;
    }

    // Construction from integer-valued samples. Set noRescaling to 'true'
    // when you want the sample values mapped directly to the buckets.
    // Typically, you know that the sample values are in the set of numbers
    // {0,1,...,numBuckets-1}, but in the event of out-of-range values, the
    // histogram stores a count for those numbers smaller than 0 and those
    // numbers larger or equal to numBuckets.
    static fromIntegerSamples(numBuckets: number, samples: readonly number[],
        noRescaling: boolean): Histogram {
        logAssert(numBuckets > 0 && samples.length > 0, 'Invalid input.');
        const histogram = new Histogram(numBuckets);
        const numSamples = samples.length;

        if (noRescaling) {
            // Map to the buckets, also counting out-of-range pixels.
            for (let i = 0; i < numSamples; ++i) {
                const value = samples[i];
                if (0 <= value) {
                    if (value < numBuckets) {
                        ++histogram.mBuckets[value];
                    } else {
                        ++histogram.mExcessGreater;
                    }
                } else {
                    ++histogram.mExcessLess;
                }
            }
        } else {
            // Compute the extremes.
            let minValue = samples[0], maxValue = minValue;
            for (let i = 1; i < numSamples; ++i) {
                const value = samples[i];
                if (value < minValue) {
                    minValue = value;
                } else if (value > maxValue) {
                    maxValue = value;
                }
            }

            // Map to the buckets.
            if (minValue < maxValue) {
                // The image is not constant.
                const numer = numBuckets - 1;
                const denom = maxValue - minValue;
                const mult = numer / denom;
                for (let i = 0; i < numSamples; ++i) {
                    const index = Math.trunc(mult * (samples[i] - minValue));
                    ++histogram.mBuckets[index];
                }
            } else {
                // The image is constant.
                histogram.mBuckets[0] = numSamples;
            }
        }

        return histogram;
    }

    // Construction from real-valued samples. The samples are rescaled so
    // the extreme values map to the first and last buckets.
    static fromRealSamples(numBuckets: number, samples: readonly number[]): Histogram {
        logAssert(numBuckets > 0 && samples.length > 0, 'Invalid input.');
        const histogram = new Histogram(numBuckets);
        const numSamples = samples.length;

        // Compute the extremes.
        let minValue = samples[0], maxValue = minValue;
        for (let i = 1; i < numSamples; ++i) {
            const value = samples[i];
            if (value < minValue) {
                minValue = value;
            } else if (value > maxValue) {
                maxValue = value;
            }
        }

        // Map to the buckets.
        if (minValue < maxValue) {
            // The image is not constant.
            const numer = numBuckets - 1;
            const denom = maxValue - minValue;
            const mult = numer / denom;
            for (let i = 0; i < numSamples; ++i) {
                const index = Math.trunc(mult * (samples[i] - minValue));
                ++histogram.mBuckets[index];
            }
        } else {
            // The image is constant.
            histogram.mBuckets[0] = numSamples;
        }

        return histogram;
    }

    // This function is called when you have used the Histogram(numBuckets)
    // constructor. No bounds checking is used; you must ensure that the
    // input value is in {0,...,numBuckets-1}.
    insert(value: number): void {
        ++this.mBuckets[value];
    }

    // This function is called when you have used the Histogram(numBuckets)
    // constructor. Bounds checking is used.
    insertCheck(value: number): void {
        if (0 <= value) {
            if (value < this.mBuckets.length) {
                ++this.mBuckets[value];
            } else {
                ++this.mExcessGreater;
            }
        } else {
            ++this.mExcessLess;
        }
    }

    // Member access.
    getBuckets(): readonly number[] {
        return this.mBuckets;
    }

    getExcessLess(): number {
        return this.mExcessLess;
    }

    getExcessGreater(): number {
        return this.mExcessGreater;
    }

    // In the following, define cdf(V) = sum_{i=0}^{V} bucket[i], where
    // 0 <= V < B and B is the number of buckets. Define N = cdf(B-1),
    // which must be the number of pixels in the image.

    // Get the lower tail of the histogram. The returned index L has the
    // properties: cdf(L-1)/N < tailAmount and cdf(L)/N >= tailAmount.
    getLowerTail(tailAmount: number): number {
        const numBuckets = this.mBuckets.length;
        let hSum = 0;
        for (let i = 0; i < numBuckets; ++i) {
            hSum += this.mBuckets[i];
        }

        const hTailSum = Math.trunc(tailAmount * hSum);
        let hLowerSum = 0;
        let lower: number;
        for (lower = 0; lower < numBuckets; ++lower) {
            hLowerSum += this.mBuckets[lower];
            if (hLowerSum >= hTailSum) {
                break;
            }
        }
        return lower;
    }

    // Get the upper tail of the histogram. The returned index U has the
    // properties: cdf(U)/N >= 1-tailAmount and cdf(U+1) < 1-tailAmount.
    getUpperTail(tailAmount: number): number {
        const numBuckets = this.mBuckets.length;
        let hSum = 0;
        for (let i = 0; i < numBuckets; ++i) {
            hSum += this.mBuckets[i];
        }

        const hTailSum = Math.trunc(tailAmount * hSum);
        let hUpperSum = 0;
        let upper: number;
        for (upper = numBuckets - 1; upper >= 0; --upper) {
            hUpperSum += this.mBuckets[upper];
            if (hUpperSum >= hTailSum) {
                break;
            }
        }
        return upper;
    }

    // Get the lower and upper tails of the histogram. The returned indices
    // are L and U and have the properties:
    //   cdf(L-1)/N < tailAmount/2, cdf(L)/N >= tailAmount/2,
    //   cdf(U)/N >= 1-tailAmount/2, and cdf(U+1) < 1-tailAmount/2.
    getTails(tailAmount: number): { lower: number, upper: number } {
        const numBuckets = this.mBuckets.length;
        let hSum = 0;
        for (let i = 0; i < numBuckets; ++i) {
            hSum += this.mBuckets[i];
        }

        const hTailSum = Math.trunc(0.5 * tailAmount * hSum);
        let hLowerSum = 0;
        let lower: number;
        for (lower = 0; lower < numBuckets; ++lower) {
            hLowerSum += this.mBuckets[lower];
            if (hLowerSum >= hTailSum) {
                break;
            }
        }

        let hUpperSum = 0;
        let upper: number;
        for (upper = numBuckets - 1; upper >= 0; --upper) {
            hUpperSum += this.mBuckets[upper];
            if (hUpperSum >= hTailSum) {
                break;
            }
        }

        return { lower, upper };
    }
}
