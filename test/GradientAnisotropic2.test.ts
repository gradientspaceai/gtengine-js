import { describe, it, expect } from 'vitest';
import { GradientAnisotropic2 } from '../src/GradientAnisotropic2.js';
import { GaussianBlur2 } from '../src/GaussianBlur2.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

// Exposes the protected conductance parameters for verification.
class TestGradientAnisotropic2 extends GradientAnisotropic2 {
    parameter(): number {
        return this.mParameter;
    }

    mHalf(): number {
        return this.mMHalfParameter;
    }
}

function build(xBound: number, yBound: number,
    f: (x: number, y: number) => number): number[] {
    const data: number[] = [];
    for (let y = 0; y < yBound; ++y) {
        for (let x = 0; x < xBound; ++x) {
            data.push(f(x, y));
        }
    }
    return data;
}

describe('GradientAnisotropic2', () => {
    it('computes the conductance parameter from the average squared gradient', () => {
        // u = x on a 4x4 image with Neumann borders. The centered x-derivative
        // is 1 at the interior columns and 0.5 at the two border columns
        // (where the padded neighbor duplicates the column itself); uy = 0
        // everywhere. The average of |grad u|^2 is
        // (0.25 + 1 + 1 + 0.25)*4/16 = 0.625.
        const filter = new TestGradientAnisotropic2(4, 4, 1, 1,
            build(4, 4, (x) => x), null, NEUMANN, PdeFilterScaleType.NONE, 1);
        expect(filter.parameter()).toBeCloseTo(1 / 0.625, 12);
        expect(filter.mHalf()).toBeCloseTo(-0.5 / 0.625, 12);
    });

    it('scales the parameter by 1/K^2', () => {
        const data = build(4, 4, (x) => x);
        const one = new TestGradientAnisotropic2(4, 4, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE, 1);
        const two = new TestGradientAnisotropic2(4, 4, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE, 2);
        expect(two.parameter()).toBeCloseTo(one.parameter() / 4, 12);
    });

    it('produces finite values on every pixel, including the last one', () => {
        // Regression for the upstream ComputeParameter off-by-one: iterating
        // over padded coordinates reads one element past the end of the padded
        // buffer, which is undefined -> NaN in TypeScript and would poison the
        // whole image.
        const filter = new TestGradientAnisotropic2(6, 5, 1, 1,
            build(6, 5, (x, y) => x * x + 2 * y), null, NEUMANN,
            PdeFilterScaleType.NONE, 1.5);
        expect(Number.isFinite(filter.parameter())).toBe(true);
        filter.setTimeStep(0.05);
        for (let i = 0; i < 5; ++i) {
            filter.update();
        }
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 6; ++x) {
                expect(Number.isFinite(filter.getU(x, y))).toBe(true);
            }
        }
    });

    it('reproduces the upstream update formula on a hand-computed pixel', () => {
        const data = [
            0, 1, 4,
            2, 3, 5,
            6, 7, 8
        ];
        const K = 1.25;
        const filter = new TestGradientAnisotropic2(3, 3, 1, 1, data, null, 0,
            PdeFilterScaleType.NONE, K);
        const dt = 0.05;
        filter.setTimeStep(dt);

        const mHalf = filter.mHalf();
        const u = (x: number, y: number) => data[x + 3 * y];
        const uxFwd = u(2, 1) - u(1, 1);
        const uxBwd = u(1, 1) - u(0, 1);
        const uyFwd = u(1, 2) - u(1, 1);
        const uyBwd = u(1, 1) - u(1, 0);
        const uxCenM = 0.5 * (u(2, 0) - u(0, 0));
        const uxCenZ = 0.5 * (u(2, 1) - u(0, 1));
        const uxCenP = 0.5 * (u(2, 2) - u(0, 2));
        const uyCenM = 0.5 * (u(0, 2) - u(0, 0));
        const uyCenZ = 0.5 * (u(1, 2) - u(1, 0));
        const uyCenP = 0.5 * (u(2, 2) - u(2, 0));
        const conductance = (s: number) => Math.exp(mHalf * s);
        const uyEstP = 0.5 * (uyCenZ + uyCenP);
        const cxp = conductance(uxCenZ * uxCenZ + uyEstP * uyEstP);
        const uyEstM = 0.5 * (uyCenZ + uyCenM);
        const cxm = conductance(uxCenZ * uxCenZ + uyEstM * uyEstM);
        const uxEstP = 0.5 * (uxCenZ + uxCenP);
        const cyp = conductance(uyCenZ * uyCenZ + uxEstP * uxEstP);
        const uxEstM = 0.5 * (uxCenZ + uxCenM);
        const cym = conductance(uyCenZ * uyCenZ + uxEstM * uxEstM);
        const expected = u(1, 1) + dt * (cxp * uxFwd - cxm * uxBwd
            + cyp * uyFwd - cym * uyBwd);

        filter.update();
        expect(filter.getU(1, 1)).toBeCloseTo(expected, 12);
    });

    it('recomputes the parameter on every iteration', () => {
        const filter = new TestGradientAnisotropic2(7, 7, 1, 1,
            build(7, 7, (x, y) => (x < 3 ? 0 : 1) + 0.1 * y), null, NEUMANN,
            PdeFilterScaleType.NONE, 1);
        filter.setTimeStep(0.05);
        const before = filter.parameter();
        for (let i = 0; i < 10; ++i) {
            filter.update();
        }
        const after = filter.parameter();
        // Diffusion lowers the average squared gradient, so the parameter,
        // which is its reciprocal, grows.
        expect(after).toBeGreaterThan(before);
    });

    it('preserves an edge better than an isotropic Gaussian blur', () => {
        const bound = 11;
        const step = (x: number) => (x < 5 ? 0 : 1);
        const data = build(bound, bound, (x) => step(x));
        const dt = 0.05;
        const iterations = 30;
        const K = 0.05;   // small K keeps the conductance tiny at the edge

        const aniso = new GradientAnisotropic2(bound, bound, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE, K);
        aniso.setTimeStep(dt);
        const blur = new GaussianBlur2(bound, bound, 1, 1, data, null, NEUMANN,
            PdeFilterScaleType.NONE);
        blur.setTimeStep(dt);
        for (let i = 0; i < iterations; ++i) {
            aniso.update();
            blur.update();
        }

        const y = 5;
        const anisoJump = aniso.getU(5, y) - aniso.getU(4, y);
        const blurJump = blur.getU(5, y) - blur.getU(4, y);
        expect(anisoJump).toBeGreaterThan(blurJump);
        // The anisotropic filter keeps nearly the whole unit step.
        expect(anisoJump).toBeGreaterThan(0.9);
        expect(blurJump).toBeLessThan(0.6);

        // Far from the edge, both filters leave the flat regions flat.
        expect(aniso.getU(1, y)).toBeCloseTo(0, 6);
        expect(aniso.getU(9, y)).toBeCloseTo(1, 6);
    });

    it('smooths within a region while the edge is preserved', () => {
        const bound = 11;
        const data = build(bound, bound, (x, y) =>
            (x < 5 ? 0 : 1) + (x === 2 && y === 5 ? 0.2 : 0));
        const filter = new GradientAnisotropic2(bound, bound, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE, 0.5);
        filter.setTimeStep(0.05);
        const noise = filter.getU(2, 5) - filter.getU(2, 4);
        for (let i = 0; i < 20; ++i) {
            filter.update();
        }
        const residual = filter.getU(2, 5) - filter.getU(2, 4);
        expect(Math.abs(residual)).toBeLessThan(Math.abs(noise));
    });
});
