// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FastMarch.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The topic of fast marching methods are discussed in the book
//   Level Set Methods and Fast Marching Methods:
//     Evolving Interfaces in Computational Geometry, Fluid Mechanics,
//     Computer Vision, and Materials Science
//   J.A. Sethian,
//   Cambridge University Press, 1999
//
// Port notes: the two upstream protected constructors (per-pixel speeds
// versus a single constant speed) are merged into one whose 'speeds'
// parameter is either an array or a single number. The upstream
// std::numeric_limits<Real>::max() sentinel becomes Number.MAX_VALUE. The
// out-parameters of GetTimeExtremes and GetInterior become return values
// ({ minValue, maxValue } and number[]). The mTrials records are
// MinHeapRecord<number, number> | null, matching the upstream
// MinHeap<size_t, Real>::Record* array whose null entries mark non-trial
// pixels.

import { MinHeap } from './MinHeap';
import type { MinHeapRecord } from './MinHeap';

export abstract class FastMarch {
    // The seed points have a crossing time of 0. As the iterations occur,
    // some of the non-seed points are visited by the moving front. Define
    // maxReal to be Number.MAX_VALUE. The valid crossing times are
    // 0 <= t < maxReal. A value of maxReal indicates the pixel has not yet
    // been reached by the moving front. If the speed value at a pixel is 0,
    // the pixel is marked with a time of -maxReal. Such pixels can never be
    // visited; the minus sign distinguishes these from pixels not yet
    // reached during iteration.
    //
    // Trial pixels are identified by having min-heap records associated
    // with them. Known or far pixels have no associated record.
    //
    // The speeds must be nonnegative and are inverted because the
    // reciprocals are all that are needed in the numerical method.

    protected mQuantity: number;
    protected mTimes: number[];
    protected mInvSpeeds: number[];
    protected mHeap: MinHeap<number, number>;
    protected mTrials: (MinHeapRecord<number, number> | null)[];

    // The 'speeds' input is either a per-pixel array of length 'quantity'
    // or a single constant speed for all pixels.
    protected constructor(quantity: number, seeds: readonly number[],
        speeds: readonly number[] | number) {
        this.mQuantity = quantity;
        this.mTimes = new Array<number>(quantity).fill(Number.MAX_VALUE);
        this.mHeap = new MinHeap<number, number>(quantity);
        this.mTrials = new Array<MinHeapRecord<number, number> | null>(quantity).fill(null);

        if (typeof speeds === 'number') {
            this.mInvSpeeds = new Array<number>(quantity).fill(1 / speeds);
            for (const seed of seeds) {
                this.mTimes[seed] = 0;
            }
        } else {
            this.mInvSpeeds = new Array<number>(quantity);
            for (const seed of seeds) {
                this.mTimes[seed] = 0;
            }

            for (let i = 0; i < this.mQuantity; ++i) {
                if (speeds[i] > 0) {
                    this.mInvSpeeds[i] = 1 / speeds[i];
                } else {
                    this.mInvSpeeds[i] = Number.MAX_VALUE;
                    this.mTimes[i] = -Number.MAX_VALUE;
                }
            }
        }
    }

    // Member access.
    getQuantity(): number {
        return this.mQuantity;
    }

    setTime(i: number, time: number): void {
        this.mTimes[i] = time;
    }

    getTime(i: number): number {
        return this.mTimes[i];
    }

    getTimeExtremes(): { minValue: number, maxValue: number } {
        let minValue = Number.MAX_VALUE;
        let maxValue = -Number.MAX_VALUE;
        let i: number;
        for (i = 0; i < this.mQuantity; ++i) {
            if (this.isValid(i)) {
                minValue = this.mTimes[i];
                maxValue = minValue;
                break;
            }
        }

        // Assert: At least one time must be valid, in which case
        // i < mQuantity at this point. If all times are invalid,
        // minValue = +maxReal and maxValue = -maxReal on exit.

        for (/**/; i < this.mQuantity; ++i) {
            if (this.isValid(i)) {
                if (this.mTimes[i] < minValue) {
                    minValue = this.mTimes[i];
                } else if (this.mTimes[i] > maxValue) {
                    maxValue = this.mTimes[i];
                }
            }
        }

        return { minValue, maxValue };
    }

    // Image element classification.
    isValid(i: number): boolean {
        return 0 <= this.mTimes[i] && this.mTimes[i] < Number.MAX_VALUE;
    }

    isTrial(i: number): boolean {
        return this.mTrials[i] !== null;
    }

    isFar(i: number): boolean {
        return this.mTimes[i] === Number.MAX_VALUE;
    }

    isZeroSpeed(i: number): boolean {
        return this.mTimes[i] === -Number.MAX_VALUE;
    }

    isInterior(i: number): boolean {
        return this.isValid(i) && !this.isTrial(i);
    }

    getInterior(): number[] {
        const interior: number[] = [];
        for (let i = 0; i < this.mQuantity; ++i) {
            if (this.isValid(i) && !this.isTrial(i)) {
                interior.push(i);
            }
        }
        return interior;
    }

    abstract getBoundary(): number[];
    abstract isBoundary(i: number): boolean;

    // Run one step of the fast marching algorithm.
    abstract iterate(): void;
}
