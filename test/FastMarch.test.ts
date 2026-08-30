import { describe, it, expect } from 'vitest';
import { FastMarch } from '../src/FastMarch';

// FastMarch is an abstract base class; exercise it through a minimal
// concrete 1-dimensional subclass that exposes the protected trial/heap
// machinery the way FastMarch2/FastMarch3 use it.
class TestFastMarch extends FastMarch {
    constructor(quantity: number, seeds: readonly number[],
        speeds: readonly number[] | number) {
        super(quantity, seeds, speeds);
    }

    getInvSpeed(i: number): number {
        return this.mInvSpeeds[i];
    }

    addTrial(i: number, time: number): void {
        this.mTimes[i] = time;
        this.mTrials[i] = this.mHeap.insert(i, time);
    }

    updateTrial(i: number, time: number): void {
        this.mTimes[i] = time;
        this.mHeap.update(this.mTrials[i], time);
    }

    override getBoundary(): number[] {
        const boundary: number[] = [];
        for (let i = 0; i < this.mQuantity; ++i) {
            if (this.isBoundary(i)) {
                boundary.push(i);
            }
        }
        return boundary;
    }

    override isBoundary(i: number): boolean {
        // A valid non-trial pixel with a 1-dimensional neighbor that is not
        // valid.
        if (!this.isValid(i) || this.isTrial(i)) {
            return false;
        }
        const left = i > 0 ? this.isValid(i - 1) : true;
        const right = i < this.mQuantity - 1 ? this.isValid(i + 1) : true;
        return !left || !right;
    }

    override iterate(): void {
        // Remove the minimum-time trial pixel and finalize it.
        const minimum = this.mHeap.remove();
        if (minimum !== null) {
            this.mTrials[minimum.key] = null;
        }
    }
}

describe('FastMarch', () => {
    it('constructs with a constant speed and marks seeds with time zero', () => {
        const march = new TestFastMarch(8, [2, 5], 2);
        expect(march.getQuantity()).toBe(8);
        for (let i = 0; i < 8; ++i) {
            expect(march.getInvSpeed(i)).toBe(0.5);
            if (i === 2 || i === 5) {
                expect(march.getTime(i)).toBe(0);
                expect(march.isValid(i)).toBe(true);
                expect(march.isFar(i)).toBe(false);
                expect(march.isInterior(i)).toBe(true);
            } else {
                expect(march.getTime(i)).toBe(Number.MAX_VALUE);
                expect(march.isValid(i)).toBe(false);
                expect(march.isFar(i)).toBe(true);
                expect(march.isInterior(i)).toBe(false);
            }
            expect(march.isTrial(i)).toBe(false);
            expect(march.isZeroSpeed(i)).toBe(false);
        }
    });

    it('constructs with per-pixel speeds and marks zero-speed pixels', () => {
        const speeds = [1, 2, 0, 4, 0.5];
        const march = new TestFastMarch(5, [0], speeds);
        expect(march.getInvSpeed(0)).toBe(1);
        expect(march.getInvSpeed(1)).toBe(0.5);
        expect(march.getInvSpeed(2)).toBe(Number.MAX_VALUE);
        expect(march.getInvSpeed(3)).toBe(0.25);
        expect(march.getInvSpeed(4)).toBe(2);

        // The zero-speed pixel has time -maxReal and is neither valid nor
        // far.
        expect(march.getTime(2)).toBe(-Number.MAX_VALUE);
        expect(march.isZeroSpeed(2)).toBe(true);
        expect(march.isValid(2)).toBe(false);
        expect(march.isFar(2)).toBe(false);

        // A zero-speed seed would have been overwritten to -maxReal only if
        // it had zero speed; seed 0 has speed 1 and remains time 0.
        expect(march.getTime(0)).toBe(0);
    });

    it('sets and gets times', () => {
        const march = new TestFastMarch(4, [0], 1);
        march.setTime(3, 1.5);
        expect(march.getTime(3)).toBe(1.5);
        expect(march.isValid(3)).toBe(true);
        expect(march.isFar(3)).toBe(false);
    });

    it('computes time extremes over valid pixels only', () => {
        const march = new TestFastMarch(6, [1], [1, 1, 0, 1, 1, 1]);
        march.setTime(3, 4);
        march.setTime(4, 2.5);
        // Pixel 0 and 5 are far, pixel 2 is zero speed; valid times are
        // {0, 4, 2.5}.
        const { minValue, maxValue } = march.getTimeExtremes();
        expect(minValue).toBe(0);
        expect(maxValue).toBe(4);
    });

    it('returns [+maxReal, -maxReal] extremes when no time is valid', () => {
        const march = new TestFastMarch(3, [], 1);
        const { minValue, maxValue } = march.getTimeExtremes();
        expect(minValue).toBe(Number.MAX_VALUE);
        expect(maxValue).toBe(-Number.MAX_VALUE);
    });

    it('classifies trial pixels through the min-heap records', () => {
        const march = new TestFastMarch(6, [0], 1);
        march.addTrial(1, 3);
        march.addTrial(2, 1);

        expect(march.isTrial(1)).toBe(true);
        expect(march.isTrial(2)).toBe(true);
        expect(march.isTrial(3)).toBe(false);

        // Trial pixels are valid but not interior.
        expect(march.isValid(1)).toBe(true);
        expect(march.isInterior(1)).toBe(false);

        // getInterior reports valid non-trial pixels.
        expect(march.getInterior()).toEqual([0]);

        // Iterating finalizes the minimum-time trial pixel first.
        march.iterate();
        expect(march.isTrial(2)).toBe(false);
        expect(march.isTrial(1)).toBe(true);
        expect(march.getInterior()).toEqual([0, 2]);
    });

    it('supports updating trial times through the heap records', () => {
        const march = new TestFastMarch(6, [0], 1);
        march.addTrial(1, 3);
        march.addTrial(4, 5);
        march.updateTrial(4, 0.5);

        // After the update, pixel 4 has the minimum time.
        march.iterate();
        expect(march.isTrial(4)).toBe(false);
        expect(march.isTrial(1)).toBe(true);
    });

    it('exposes the boundary through the concrete subclass', () => {
        const march = new TestFastMarch(6, [2, 3], 1);
        // Pixels 2 and 3 are interior/valid; both neighbor invalid pixels.
        expect(march.getBoundary()).toEqual([2, 3]);
        expect(march.isBoundary(2)).toBe(true);
        expect(march.isBoundary(0)).toBe(false);
    });
});
