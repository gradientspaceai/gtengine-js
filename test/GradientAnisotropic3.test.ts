import { describe, it, expect } from 'vitest';
import { GradientAnisotropic3 } from '../src/GradientAnisotropic3.js';
import { GaussianBlur3 } from '../src/GaussianBlur3.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';

const NEUMANN = Number.MAX_VALUE;

// Exposes the protected conductance parameters for verification.
class TestGradientAnisotropic3 extends GradientAnisotropic3 {
    parameter(): number {
        return this.mParameter;
    }

    mHalf(): number {
        return this.mMHalfParameter;
    }
}

function buildBounds(xBound: number, yBound: number, zBound: number,
    f: (x: number, y: number, z: number) => number): number[] {
    const data: number[] = [];
    for (let z = 0; z < zBound; ++z) {
        for (let y = 0; y < yBound; ++y) {
            for (let x = 0; x < xBound; ++x) {
                data.push(f(x, y, z));
            }
        }
    }
    return data;
}

function build(bound: number, f: (x: number, y: number, z: number) => number): number[] {
    return buildBounds(bound, bound, bound, f);
}

describe('GradientAnisotropic3', () => {
    it('computes the conductance parameter from the average squared gradient', () => {
        // u = x on a 4x4x4 image with Neumann borders. Along each row the
        // centered x-derivatives are 0.5, 1, 1, 0.5 and uy = uz = 0, so the
        // average of |grad u|^2 is 2.5/4 = 0.625.
        const filter = new TestGradientAnisotropic3(4, 4, 4, 1, 1, 1,
            build(4, (x) => x), null, NEUMANN, PdeFilterScaleType.NONE, 1);
        expect(filter.parameter()).toBeCloseTo(1 / 0.625, 12);
        expect(filter.mHalf()).toBeCloseTo(-0.5 / 0.625, 12);
    });

    it('scales the parameter by 1/K^2', () => {
        const data = build(4, (x, y) => x + y);
        const one = new TestGradientAnisotropic3(4, 4, 4, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE, 1);
        const half = new TestGradientAnisotropic3(4, 4, 4, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE, 0.5);
        expect(half.parameter()).toBeCloseTo(4 * one.parameter(), 10);
    });

    it('produces finite values on every voxel, including the last one', () => {
        // Regression for the upstream ComputeParameter off-by-one, which reads
        // past the end of the padded buffer on the final voxel.
        const filter = new TestGradientAnisotropic3(5, 4, 6, 1, 1, 1,
            buildBounds(5, 4, 6, (x, y, z) => x * x + 2 * y - z),
            null, NEUMANN, PdeFilterScaleType.NONE, 1.5);
        expect(Number.isFinite(filter.parameter())).toBe(true);
        filter.setTimeStep(0.05);
        for (let i = 0; i < 4; ++i) {
            filter.update();
        }
        for (let z = 0; z < 6; ++z) {
            for (let y = 0; y < 4; ++y) {
                for (let x = 0; x < 5; ++x) {
                    expect(Number.isFinite(filter.getU(x, y, z))).toBe(true);
                }
            }
        }
    });

    it('reproduces the upstream update formula on a hand-computed voxel', () => {
        const bound = 3;
        const value = (x: number, y: number, z: number) =>
            1 + x + 2 * y + 3 * z + 0.5 * x * y - 0.25 * y * z + 0.75 * x * z;
        const data = build(bound, value);
        const K = 0.8;
        const filter = new TestGradientAnisotropic3(bound, bound, bound, 1, 1, 1,
            data, null, 0, PdeFilterScaleType.NONE, K);
        const dt = 0.02;
        filter.setTimeStep(dt);

        const mHalf = filter.mHalf();
        const shift = -Math.min(...data);
        const u = (x: number, y: number, z: number) =>
            data[x + bound * (y + bound * z)] + shift;
        const c = (s: number) => Math.exp(mHalf * s);

        const uxFwd = u(2, 1, 1) - u(1, 1, 1);
        const uxBwd = u(1, 1, 1) - u(0, 1, 1);
        const uyFwd = u(1, 2, 1) - u(1, 1, 1);
        const uyBwd = u(1, 1, 1) - u(1, 0, 1);
        const uzFwd = u(1, 1, 2) - u(1, 1, 1);
        const uzBwd = u(1, 1, 1) - u(1, 1, 0);

        const duvzz = 0.5 * (u(2, 1, 1) - u(0, 1, 1));
        const duvpz = 0.5 * (u(2, 2, 1) - u(0, 2, 1));
        const duvmz = 0.5 * (u(2, 0, 1) - u(0, 0, 1));
        const duvzp = 0.5 * (u(2, 1, 2) - u(0, 1, 2));
        const duvzm = 0.5 * (u(2, 1, 0) - u(0, 1, 0));

        const duzvz = 0.5 * (u(1, 2, 1) - u(1, 0, 1));
        const dupvz = 0.5 * (u(2, 2, 1) - u(2, 0, 1));
        const dumvz = 0.5 * (u(0, 2, 1) - u(0, 0, 1));
        const duzvp = 0.5 * (u(1, 2, 2) - u(1, 0, 2));
        const duzvm = 0.5 * (u(1, 2, 0) - u(1, 0, 0));

        const duzzv = 0.5 * (u(1, 1, 2) - u(1, 1, 0));
        const dupzv = 0.5 * (u(2, 1, 2) - u(2, 1, 0));
        const dumzv = 0.5 * (u(0, 1, 2) - u(0, 1, 0));
        const duzpv = 0.5 * (u(1, 2, 2) - u(1, 2, 0));
        const duzmv = 0.5 * (u(1, 0, 2) - u(1, 0, 0));

        const uxCenSqr = duvzz * duvzz;
        const uyCenSqr = duzvz * duzvz;
        const uzCenSqr = duzzv * duzzv;

        let uxEst: number, uyEst: number, uzEst: number;
        uyEst = 0.5 * (duzvz + dupvz);
        uzEst = 0.5 * (duzzv + dupzv);
        const cxp = c(uxCenSqr + uyEst * uyEst + uzEst * uzEst);
        uyEst = 0.5 * (duzvz + dumvz);
        uzEst = 0.5 * (duzzv + dumzv);
        const cxm = c(uxCenSqr + uyEst * uyEst + uzEst * uzEst);
        uxEst = 0.5 * (duvzz + duvpz);
        uzEst = 0.5 * (duzzv + duzpv);
        const cyp = c(uxEst * uxEst + uyCenSqr + uzEst * uzEst);
        uxEst = 0.5 * (duvzz + duvmz);
        uzEst = 0.5 * (duzzv + duzmv);
        const cym = c(uxEst * uxEst + uyCenSqr + uzEst * uzEst);
        uxEst = 0.5 * (duvzz + duvzp);
        uyEst = 0.5 * (duzvz + duzvp);
        const czp = c(uxEst * uxEst + uyEst * uyEst + uzCenSqr);
        uxEst = 0.5 * (duvzz + duvzm);
        uyEst = 0.5 * (duzvz + duzvm);
        const czm = c(uxEst * uxEst + uyEst * uyEst + uzCenSqr);

        const expected = u(1, 1, 1) + dt * (cxp * uxFwd - cxm * uxBwd
            + cyp * uyFwd - cym * uyBwd + czp * uzFwd - czm * uzBwd);

        filter.update();
        expect(filter.getU(1, 1, 1)).toBeCloseTo(expected, 12);
    });

    it('recomputes the parameter on every iteration', () => {
        const filter = new TestGradientAnisotropic3(7, 7, 7, 1, 1, 1,
            build(7, (x, y) => (x < 3 ? 0 : 1) + 0.1 * y), null, NEUMANN,
            PdeFilterScaleType.NONE, 1);
        filter.setTimeStep(0.05);
        const before = filter.parameter();
        for (let i = 0; i < 8; ++i) {
            filter.update();
        }
        expect(filter.parameter()).toBeGreaterThan(before);
    });

    it('preserves an edge better than an isotropic Gaussian blur', () => {
        const bound = 9;
        const data = build(bound, (x) => (x < 4 ? 0 : 1));
        const dt = 0.05;
        const iterations = 25;

        const aniso = new GradientAnisotropic3(bound, bound, bound, 1, 1, 1, data,
            null, NEUMANN, PdeFilterScaleType.NONE, 0.05);
        aniso.setTimeStep(dt);
        const blur = new GaussianBlur3(bound, bound, bound, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE);
        blur.setTimeStep(dt);
        for (let i = 0; i < iterations; ++i) {
            aniso.update();
            blur.update();
        }

        const anisoJump = aniso.getU(4, 4, 4) - aniso.getU(3, 4, 4);
        const blurJump = blur.getU(4, 4, 4) - blur.getU(3, 4, 4);
        expect(anisoJump).toBeGreaterThan(blurJump);
        expect(anisoJump).toBeGreaterThan(0.9);
        expect(blurJump).toBeLessThan(0.6);
    });

    it('is invariant under permuting the coordinate axes of a symmetric image', () => {
        const bound = 7, c = 3;
        const radial = (x: number, y: number, z: number) =>
            Math.sqrt((x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2);
        const filter = new GradientAnisotropic3(bound, bound, bound, 1, 1, 1,
            build(bound, radial), null, NEUMANN, PdeFilterScaleType.NONE, 1);
        filter.setTimeStep(0.05);
        filter.update();
        for (let z = 0; z < bound; ++z) {
            for (let y = 0; y < bound; ++y) {
                for (let x = 0; x < bound; ++x) {
                    expect(filter.getU(x, y, z)).toBeCloseTo(filter.getU(y, x, z), 12);
                    expect(filter.getU(x, y, z)).toBeCloseTo(filter.getU(x, z, y), 12);
                    expect(filter.getU(x, y, z)).toBeCloseTo(
                        filter.getU(bound - 1 - x, y, z), 12);
                }
            }
        }
    });
});
