import { describe, it, expect } from 'vitest';
import { GaussianBlur3 } from '../src/GaussianBlur3.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

function index(x: number, y: number, z: number, xBound: number, yBound: number): number {
    return x + xBound * (y + yBound * z);
}

// A spike image: all zeros except a single 1 at (cx,cy,cz).
function spike(bound: number, cx: number, cy: number, cz: number): number[] {
    const data = new Array<number>(bound * bound * bound).fill(0);
    data[index(cx, cy, cz, bound, bound)] = 1;
    return data;
}

function total(filter: GaussianBlur3, bound: number): number {
    let sum = 0;
    for (let z = 0; z < bound; ++z) {
        for (let y = 0; y < bound; ++y) {
            for (let x = 0; x < bound; ++x) {
                sum += filter.getU(x, y, z);
            }
        }
    }
    return sum;
}

describe('GaussianBlur3', () => {
    it('computes the stable maximum time step from the spacings', () => {
        const data = spike(5, 2, 2, 2);
        const unit = new GaussianBlur3(5, 5, 5, 1, 1, 1, data, null, 0,
            PdeFilterScaleType.NONE);
        expect(unit.getMaximumTimeStep()).toBeCloseTo(0.5 / 3, 15);

        const scaled = new GaussianBlur3(5, 5, 5, 2, 0.5, 1, data, null, 0,
            PdeFilterScaleType.NONE);
        // 0.5 / (1/4 + 4 + 1).
        expect(scaled.getMaximumTimeStep()).toBeCloseTo(0.5 / 5.25, 15);
    });

    it('takes the known explicit-Euler step on a spike at the maximum time step', () => {
        const filter = new GaussianBlur3(5, 5, 5, 1, 1, 1, spike(5, 2, 2, 2), null, 0,
            PdeFilterScaleType.NONE);
        const dt = filter.getMaximumTimeStep();   // 1/6
        filter.setTimeStep(dt);
        filter.update();

        // The center loses all of its mass to the six face neighbors.
        expect(filter.getU(2, 2, 2)).toBeCloseTo(0, 15);
        expect(filter.getU(1, 2, 2)).toBeCloseTo(1 / 6, 15);
        expect(filter.getU(3, 2, 2)).toBeCloseTo(1 / 6, 15);
        expect(filter.getU(2, 1, 2)).toBeCloseTo(1 / 6, 15);
        expect(filter.getU(2, 3, 2)).toBeCloseTo(1 / 6, 15);
        expect(filter.getU(2, 2, 1)).toBeCloseTo(1 / 6, 15);
        expect(filter.getU(2, 2, 3)).toBeCloseTo(1 / 6, 15);
        // Edge and corner neighbors are untouched after one step.
        expect(filter.getU(1, 1, 2)).toBeCloseTo(0, 15);
        expect(filter.getU(1, 1, 1)).toBeCloseTo(0, 15);
    });

    it('conserves total mass while the front has not reached the border', () => {
        const filter = new GaussianBlur3(7, 7, 7, 1, 1, 1, spike(7, 3, 3, 3), null, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(0.1);
        for (let i = 0; i < 3; ++i) {
            filter.update();
            expect(total(filter, 7)).toBeCloseTo(1, 12);
        }
    });

    it('smooths a spike and preserves the octahedral symmetry', () => {
        const bound = 7;
        const filter = new GaussianBlur3(bound, bound, bound, 1, 1, 1,
            spike(bound, 3, 3, 3), null, 0, PdeFilterScaleType.NONE);
        filter.setTimeStep(0.1);

        let peak = filter.getU(3, 3, 3);
        for (let i = 0; i < 4; ++i) {
            filter.update();
            const newPeak = filter.getU(3, 3, 3);
            expect(newPeak).toBeLessThan(peak);
            peak = newPeak;
        }

        for (let z = 0; z < bound; ++z) {
            for (let y = 0; y < bound; ++y) {
                for (let x = 0; x < bound; ++x) {
                    const value = filter.getU(x, y, z);
                    expect(Number.isFinite(value)).toBe(true);
                    expect(value).toBeGreaterThanOrEqual(0);
                    expect(value).toBeLessThanOrEqual(peak + 1e-15);
                    // Reflection and coordinate-permutation symmetry.
                    expect(value).toBeCloseTo(filter.getU(6 - x, y, z), 12);
                    expect(value).toBeCloseTo(filter.getU(z, y, x), 12);
                    expect(value).toBeCloseTo(filter.getU(x, z, y), 12);
                }
            }
        }
    });

    it('anisotropic spacings weight the axes differently', () => {
        // With dz much smaller than dx and dy, diffusion along z dominates.
        const filter = new GaussianBlur3(5, 5, 5, 1, 1, 0.5, spike(5, 2, 2, 2), null, 0,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(filter.getMaximumTimeStep());
        filter.update();
        expect(filter.getU(2, 2, 1)).toBeGreaterThan(filter.getU(1, 2, 2));
        expect(filter.getU(2, 2, 1)).toBeCloseTo(filter.getU(2, 2, 3), 15);
        expect(filter.getU(1, 2, 2)).toBeCloseTo(filter.getU(2, 1, 2), 15);
    });

    it('does not diffuse an image that is already constant', () => {
        const data = new Array<number>(125).fill(-2);
        const filter = new GaussianBlur3(5, 5, 5, 1, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE);
        filter.setTimeStep(filter.getMaximumTimeStep());
        filter.update();
        // A constant image is shifted to zero by the base-class scaling
        // (upstream issue #60) and the Laplacian vanishes.
        for (let z = 0; z < 5; ++z) {
            for (let y = 0; y < 5; ++y) {
                for (let x = 0; x < 5; ++x) {
                    expect(filter.getU(x, y, z)).toBeCloseTo(0, 12);
                }
            }
        }
    });

    it('applies the SYMMETRIC scale type to the stored image', () => {
        const data = new Array<number>(27).fill(4);
        data[index(1, 1, 1, 3, 3)] = 10;
        const filter = new GaussianBlur3(3, 3, 3, 1, 1, 1, data, null, 0,
            PdeFilterScaleType.SYMMETRIC);
        // d' = -1 + 2*(d - min)/(max - min) in [-1,1].
        expect(filter.getU(1, 1, 1)).toBeCloseTo(1, 15);
        expect(filter.getU(0, 0, 0)).toBeCloseTo(-1, 15);
    });
});
