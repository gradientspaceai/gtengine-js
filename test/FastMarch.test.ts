import { describe, it, expect } from 'vitest';
import { FastMarch } from '../src/FastMarch.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream FastMarch.h base class.
// ---------------------------------------------------------------------------

// A complete 1-dimensional fast-marching solver built on the base class, so
// the trial/known/far bookkeeping, the min-heap records and the inverted
// speeds are all exercised. The upwind update in one dimension is
// T(i) = min(T(i-1), T(i+1)) + invSpeed(i), which makes the exact solution a
// shortest-path problem that the test cross-checks with Dijkstra's algorithm.
class Eikonal1D extends FastMarch {
    constructor(quantity: number, seeds: readonly number[],
        speeds: readonly number[] | number) {
        super(quantity, seeds, speeds);
        for (const seed of seeds) {
            if (!this.isZeroSpeed(seed)) {
                this.relaxNeighbors(seed);
            }
        }
    }

    private relaxNeighbors(i: number): void {
        for (const j of [i - 1, i + 1]) {
            if (j < 0 || j >= this.mQuantity || this.isZeroSpeed(j)) {
                continue;
            }
            if (this.isValid(j) && !this.isTrial(j)) {
                continue;   // already known
            }
            const time = this.mTimes[i] + this.mInvSpeeds[j];
            if (this.isFar(j)) {
                this.mTimes[j] = time;
                this.mTrials[j] = this.mHeap.insert(j, time);
            } else if (this.isTrial(j) && time < this.mTimes[j]) {
                this.mTimes[j] = time;
                this.mHeap.update(this.mTrials[j], time);
            }
        }
    }

    override iterate(): void {
        const minimum = this.mHeap.remove();
        if (minimum === null) {
            return;
        }
        const i = minimum.key;
        this.mTrials[i] = null;
        this.mTimes[i] = minimum.value;
        this.relaxNeighbors(i);
    }

    getInvSpeedValue(i: number): number {
        return this.mInvSpeeds[i];
    }

    numTrials(): number {
        return this.mHeap.getNumElements();
    }

    minimumTrialTime(): number | null {
        const minimum = this.mHeap.getMinimum();
        return minimum === null ? null : minimum.value;
    }

    override isBoundary(i: number): boolean {
        if (!this.isValid(i) || this.isTrial(i)) {
            return false;
        }
        const left = i > 0 ? this.isValid(i - 1) : true;
        const right = i < this.mQuantity - 1 ? this.isValid(i + 1) : true;
        return !left || !right;
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
}

// Independent shortest-path solution: the cost of entering pixel j is
// invSpeeds[j]; seeds start at time 0; zero-speed pixels are impassable.
function dijkstra(quantity: number, seeds: readonly number[],
    invSpeeds: readonly number[], blocked: readonly boolean[]): number[] {
    const times = new Array<number>(quantity).fill(Number.POSITIVE_INFINITY);
    const done = new Array<boolean>(quantity).fill(false);
    for (const seed of seeds) {
        if (!blocked[seed]) {
            times[seed] = 0;
        }
    }
    for (; ;) {
        let best = -1;
        for (let i = 0; i < quantity; ++i) {
            if (!done[i] && times[i] < Number.POSITIVE_INFINITY
                && (best < 0 || times[i] < times[best])) {
                best = i;
            }
        }
        if (best < 0) {
            break;
        }
        done[best] = true;
        for (const j of [best - 1, best + 1]) {
            if (j < 0 || j >= quantity || blocked[j]) {
                continue;
            }
            const time = times[best] + invSpeeds[j];
            if (time < times[j]) {
                times[j] = time;
            }
        }
    }
    return times;
}

describe('FastMarch verification', () => {
    const quantity = fc.integer({ min: 2, max: 24 });
    const configuration = quantity.chain(n => fc.tuple(
        fc.constant(n),
        fc.uniqueArray(fc.integer({ min: 0, max: n - 1 }),
            { minLength: 1, maxLength: Math.max(1, Math.min(n, 4)) }),
        // Speeds on a coarse positive grid, with zeros allowed so the
        // zero-speed (impassable) marking is exercised.
        fc.array(fc.constantFrom(0, 0.5, 1, 2, 4), { minLength: n, maxLength: n })));

    it('classifies every pixel as exactly one of valid, far or zero speed', () => {
        check(configuration, ([n, seeds, speeds]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            for (let i = 0; i < n; ++i) {
                const flags = [march.isValid(i), march.isFar(i), march.isZeroSpeed(i)];
                expect(flags.filter(f => f).length).toBe(1);
                // isInterior is isValid && !isTrial by definition.
                expect(march.isInterior(i)).toBe(march.isValid(i) && !march.isTrial(i));
            }
        });
    });

    it('seeds start at time zero unless their speed is zero', () => {
        check(configuration, ([n, seeds, speeds]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            for (const seed of seeds) {
                if (speeds[seed] > 0) {
                    expect(march.getTime(seed)).toBe(0);
                } else {
                    // The speed loop runs after the seed loop upstream, so a
                    // zero-speed seed is overwritten with -maxReal.
                    expect(march.getTime(seed)).toBe(-Number.MAX_VALUE);
                    expect(march.isZeroSpeed(seed)).toBe(true);
                }
            }
            for (let i = 0; i < n; ++i) {
                if (speeds[i] > 0) {
                    expect(march.getInvSpeedValue(i)).toBe(1 / speeds[i]);
                } else {
                    expect(march.getInvSpeedValue(i)).toBe(Number.MAX_VALUE);
                }
            }
        });
    });

    it('marching to completion reproduces the Dijkstra arrival times', () => {
        check(configuration, ([n, seeds, speeds]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            let guard = 4 * n + 4;
            while (march.numTrials() > 0 && guard-- > 0) {
                march.iterate();
            }
            expect(march.numTrials()).toBe(0);

            const invSpeeds = speeds.map(s => (s > 0 ? 1 / s : Number.MAX_VALUE));
            const blocked = speeds.map(s => !(s > 0));
            const reference = dijkstra(n, seeds, invSpeeds, blocked);
            for (let i = 0; i < n; ++i) {
                if (blocked[i]) {
                    expect(march.getTime(i)).toBe(-Number.MAX_VALUE);
                } else if (reference[i] === Number.POSITIVE_INFINITY) {
                    // Unreachable behind a zero-speed pixel: still far.
                    expect(march.isFar(i)).toBe(true);
                } else {
                    expectClose(march.getTime(i), reference[i], 1e-12, 1e-12);
                    expect(march.isValid(i)).toBe(true);
                    expect(march.isTrial(i)).toBe(false);
                }
            }
        });
    });

    it('the times removed from the front are nondecreasing', () => {
        check(configuration, ([n, seeds, speeds]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            let previous = 0;
            let guard = 4 * n + 4;
            while (march.numTrials() > 0 && guard-- > 0) {
                const next = march.minimumTrialTime();
                expect(next).not.toBeNull();
                expect(next as number).toBeGreaterThanOrEqual(previous - 1e-12);
                previous = next as number;
                march.iterate();
            }
        });
    });

    it('getInterior lists exactly the valid non-trial pixels at any stage', () => {
        check(fc.tuple(configuration, fc.nat({ max: 40 })), ([[n, seeds, speeds], steps]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            for (let s = 0; s < steps && march.numTrials() > 0; ++s) {
                march.iterate();
            }
            const brute: number[] = [];
            for (let i = 0; i < n; ++i) {
                if (march.isValid(i) && !march.isTrial(i)) {
                    brute.push(i);
                }
            }
            expect(march.getInterior()).toEqual(brute);
        });
    });

    it('getTimeExtremes equals the min and max over the valid times', () => {
        check(fc.tuple(configuration, fc.nat({ max: 40 })), ([[n, seeds, speeds], steps]) => {
            const march = new Eikonal1D(n, seeds, speeds);
            for (let s = 0; s < steps && march.numTrials() > 0; ++s) {
                march.iterate();
            }
            const valid: number[] = [];
            for (let i = 0; i < n; ++i) {
                if (march.isValid(i)) {
                    valid.push(march.getTime(i));
                }
            }
            const { minValue, maxValue } = march.getTimeExtremes();
            if (valid.length === 0) {
                expect(minValue).toBe(Number.MAX_VALUE);
                expect(maxValue).toBe(-Number.MAX_VALUE);
            } else {
                expect(minValue).toBe(Math.min(...valid));
                expect(maxValue).toBe(Math.max(...valid));
            }
        });
    });

    it('setTime drives the classification predicates', () => {
        const sentinels = fc.constantFrom(0, 0.5, 1, 17, Number.MAX_VALUE,
            -Number.MAX_VALUE, -1);
        check(fc.tuple(fc.integer({ min: 1, max: 12 }),
            fc.array(sentinels, { minLength: 1, maxLength: 12 })),
            ([n, times]) => {
                const march = new Eikonal1D(n, [0], 1);
                for (let i = 0; i < n; ++i) {
                    const time = times[i % times.length];
                    march.setTime(i, time);
                    expect(march.getTime(i)).toBe(time);
                    expect(march.isValid(i)).toBe(0 <= time && time < Number.MAX_VALUE);
                    expect(march.isFar(i)).toBe(time === Number.MAX_VALUE);
                    expect(march.isZeroSpeed(i)).toBe(time === -Number.MAX_VALUE);
                }
            });
    });

    it('a constant speed gives arrival times proportional to the seed distance', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 24 }),
            fc.constantFrom(0.5, 1, 2, 4)).chain(([n, speed]) =>
                fc.tuple(fc.constant(n), fc.constant(speed),
                    fc.uniqueArray(fc.integer({ min: 0, max: n - 1 }),
                        { minLength: 1, maxLength: Math.min(n, 3) }))),
            ([n, speed, seeds]) => {
                const march = new Eikonal1D(n, seeds, speed);
                let guard = 4 * n + 4;
                while (march.numTrials() > 0 && guard-- > 0) {
                    march.iterate();
                }
                for (let i = 0; i < n; ++i) {
                    const distance = Math.min(...seeds.map(s => Math.abs(i - s)));
                    expectClose(march.getTime(i), distance / speed, 1e-12, 1e-12);
                }
            });
    });
});
