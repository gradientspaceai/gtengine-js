import { describe, it, expect } from 'vitest';
import { GaussianBlur2 } from '../src/GaussianBlur2.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

// A spike image: all zeros except a single 1 at (cx,cy).
function spike(xBound: number, yBound: number, cx: number, cy: number): number[] {
    const data = new Array<number>(xBound * yBound).fill(0);
    data[cx + xBound * cy] = 1;
    return data;
}

function image(filter: GaussianBlur2, xBound: number, yBound: number): number[] {
    const u: number[] = [];
    for (let y = 0; y < yBound; ++y) {
        for (let x = 0; x < xBound; ++x) {
            u.push(filter.getU(x, y));
        }
    }
    return u;
}

function total(values: readonly number[]): number {
    return values.reduce((a, b) => a + b, 0);
}

describe('GaussianBlur2', () => {
    it('computes the stable maximum time step from the spacings', () => {
        const data = spike(5, 5, 2, 2);
        const unit = new GaussianBlur2(5, 5, 1, 1, data, null, 0, PdeFilterScaleType.NONE);
        // 0.5 / (1/dx^2 + 1/dy^2) = 0.5/2 = 0.25.
        expect(unit.getMaximumTimeStep()).toBeCloseTo(0.25, 15);

        const scaled = new GaussianBlur2(5, 5, 2, 0.5, data, null, 0, PdeFilterScaleType.NONE);
        // 0.5 / (1/4 + 4) = 0.5/4.25.
        expect(scaled.getMaximumTimeStep()).toBeCloseTo(0.5 / 4.25, 15);
    });

    it('takes the known explicit-Euler step on a spike at the maximum time step', () => {
        const filter = new GaussianBlur2(5, 5, 1, 1, spike(5, 5, 2, 2), null, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(filter.getMaximumTimeStep());
        filter.update();

        // u <- u + dt*(uxx + uyy) with dt = 1/4. The center loses all of its
        // mass, each 4-neighbor gains 1/4 of it, and the diagonals are still
        // untouched after one step.
        expect(filter.getU(2, 2)).toBeCloseTo(0, 15);
        expect(filter.getU(1, 2)).toBeCloseTo(0.25, 15);
        expect(filter.getU(3, 2)).toBeCloseTo(0.25, 15);
        expect(filter.getU(2, 1)).toBeCloseTo(0.25, 15);
        expect(filter.getU(2, 3)).toBeCloseTo(0.25, 15);
        expect(filter.getU(1, 1)).toBeCloseTo(0, 15);
        expect(filter.getU(3, 3)).toBeCloseTo(0, 15);
    });

    it('conserves total mass while the front has not reached the border', () => {
        const filter = new GaussianBlur2(7, 7, 1, 1, spike(7, 7, 3, 3), null, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.2);
        for (let i = 0; i < 3; ++i) {
            filter.update();
            expect(total(image(filter, 7, 7))).toBeCloseTo(1, 12);
        }
    });

    it('smooths a spike: the peak decreases and the neighborhood spreads', () => {
        const filter = new GaussianBlur2(9, 9, 1, 1, spike(9, 9, 4, 4), null, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.1);

        let peak = filter.getU(4, 4);
        for (let i = 0; i < 8; ++i) {
            filter.update();
            const newPeak = filter.getU(4, 4);
            expect(newPeak).toBeLessThan(peak);
            peak = newPeak;
        }

        // The result is still nonnegative, still peaked at the center, and
        // symmetric under the reflections that fix the spike.
        const u = image(filter, 9, 9);
        for (const value of u) {
            expect(value).toBeGreaterThanOrEqual(0);
            expect(Number.isFinite(value)).toBe(true);
        }
        for (let y = 0; y < 9; ++y) {
            for (let x = 0; x < 9; ++x) {
                expect(filter.getU(x, y)).toBeCloseTo(filter.getU(8 - x, y), 12);
                expect(filter.getU(x, y)).toBeCloseTo(filter.getU(y, x), 12);
                if (x !== 4 || y !== 4) {
                    expect(filter.getU(x, y)).toBeLessThan(peak + 1e-15);
                }
            }
        }
    });

    it('does not diffuse an image that is already constant', () => {
        const data = new Array<number>(25).fill(3);
        const filter = new GaussianBlur2(5, 5, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.25);
        filter.update();
        // The base PdeFilter stores mOffset + (d - mMin)*mScale even for
        // ScaleType.NONE, so a constant image is shifted to zero (upstream
        // issue #60). Whatever the offset, the Laplacian vanishes and the
        // image is unchanged by the update.
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                expect(filter.getU(x, y)).toBeCloseTo(0, 12);
            }
        }
    });

    it('keeps the Neumann image border at its constructor-time duplicates', () => {
        // Upstream recomputes only the mask border in OnPreUpdate, so the
        // 1-pixel image border stays at the values duplicated in the
        // constructor (see upstream issue #60, Neumann border staleness).
        // Pin that: the border-adjacent pixel evolves as if its outside
        // neighbor were frozen at the original value.
        const data = [
            0, 0, 0,
            0, 1, 0,
            0, 0, 0
        ];
        const filter = new GaussianBlur2(3, 3, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.25);
        filter.update();

        // Pixel (1,1) is the spike: uxx = uyy = -2, so u -> 1 - 1 = 0.
        expect(filter.getU(1, 1)).toBeCloseTo(0, 15);
        // Pixel (0,1): the padded neighbor at x = -1 duplicates u(0,1) = 0,
        // so uxx = 1 - 0 + 0 = 1 and uyy = 0, giving 0.25.
        expect(filter.getU(0, 1)).toBeCloseTo(0.25, 15);
        // Total mass grows because the frozen border injects flux.
        expect(total(image(filter, 3, 3))).toBeCloseTo(1, 12);
    });

    it('applies the UNIT scale type to the stored image', () => {
        const data = [
            2, 2, 2,
            2, 6, 2,
            2, 2, 2
        ];
        const filter = new GaussianBlur2(3, 3, 1, 1, data, null, 0,
            PdeFilterScaleType.UNIT);
        // d' = (d - min)/(max - min) in [0,1].
        expect(filter.getU(1, 1)).toBeCloseTo(1, 15);
        expect(filter.getU(0, 0)).toBeCloseTo(0, 15);
    });

    it('respects a mask by leaving masked-out pixels unchanged', () => {
        const data = spike(5, 5, 2, 2);
        const mask = new Array<number>(25).fill(1);
        mask[1 + 5 * 2] = 0;    // mask out the left neighbor of the spike
        const filter = new GaussianBlur2(5, 5, 1, 1, data, mask, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.25);
        filter.update();
        expect(filter.getMask(1, 2)).toBe(0);
        // The masked pixel is never visited by onUpdate, so it keeps the
        // Dirichlet mask-border value assigned in the constructor.
        expect(filter.getU(1, 2)).toBeCloseTo(0, 15);
        expect(filter.getU(3, 2)).toBeCloseTo(0.25, 15);
    });
});
