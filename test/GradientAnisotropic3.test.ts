import { describe, it, expect } from 'vitest';
import { GradientAnisotropic3 } from '../src/GradientAnisotropic3.js';
import { GaussianBlur3 } from '../src/GaussianBlur3.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

describe('GradientAnisotropic3 verification', () => {
    // The padded image the constructor holds under Neumann conditions and the
    // NONE scale type (which still subtracts the data minimum).
    function neumannPadded(xB: number, yB: number, zB: number,
        data: readonly number[]): number[] {
        const min = Math.min(...data);
        const max = Math.max(...data);
        const shift = (min === max ? () => 0 : (d: number) => d - min);
        const w = xB + 2, h = yB + 2;
        const p = new Array<number>(w * h * (zB + 2)).fill(0);
        const at = (px: number, py: number, pz: number) => px + w * (py + h * pz);
        for (let z = 0; z < zB; ++z) {
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    p[at(x + 1, y + 1, z + 1)] = shift(data[x + xB * (y + yB * z)]);
                }
            }
        }
        for (let pz = 0; pz < zB + 2; ++pz) {
            for (let py = 0; py < h; ++py) {
                for (let px = 0; px < w; ++px) {
                    if (px === 0 || px === xB + 1 || py === 0 || py === yB + 1
                        || pz === 0 || pz === zB + 1) {
                        p[at(px, py, pz)] = p[at(
                            Math.min(Math.max(px, 1), xB),
                            Math.min(Math.max(py, 1), yB),
                            Math.min(Math.max(pz, 1), zB))];
                    }
                }
            }
        }
        return p;
    }

    function values(filter: GradientAnisotropic3, xB: number, yB: number,
        zB: number): number[] {
        const u: number[] = [];
        for (let z = 0; z < zB; ++z) {
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    u.push(filter.getU(x, y, z));
                }
            }
        }
        return u;
    }

    // Nonconstant images only: a constant image has zero average squared
    // gradient, so mParameter is infinite and every conductance is NaN (an
    // upstream quirk with no guard, preserved by the port).
    const config = fc.tuple(
        fc.integer({ min: 3, max: 5 }),
        fc.integer({ min: 3, max: 5 }),
        fc.integer({ min: 3, max: 5 }),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2, 5),
        fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 125, maxLength: 125 }))
        .map(([xB, yB, zB, hx, hy, hz, K, pool]) => ({
            xB, yB, zB, hx, hy, hz, K, data: pool.slice(0, xB * yB * zB)
        }))
        .filter(({ data }) => new Set(data).size > 1);

    it('the conductance parameter averages over the unpadded coordinates', () => {
        // The average must be taken over the unpadded domain the accessors
        // document and divided by mQuantity; upstream feeds padded
        // coordinates and reads past the padded buffer (upstream issue #122).
        check(config, ({ xB, yB, zB, hx, hy, hz, K, data }) => {
            const filter = new TestGradientAnisotropic3(xB, yB, zB, hx, hy, hz,
                data, null, NEUMANN, PdeFilterScaleType.NONE, K);
            let sum = 0;
            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        const ux = filter.getUx(x, y, z);
                        const uy = filter.getUy(x, y, z);
                        const uz = filter.getUz(x, y, z);
                        sum += ux * ux + uy * uy + uz * uz;
                    }
                }
            }
            const avg = sum / (xB * yB * zB);
            expect(Number.isFinite(filter.parameter())).toBe(true);
            expectClose(filter.parameter(), 1 / (K * K * avg), 1e-9, 1e-9);
            expectClose(filter.mHalf(), -0.5 * filter.parameter(), 1e-12, 1e-12);
        });
    });

    it('one step reproduces the upstream update from the padded image', () => {
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xB, yB, zB, hx, hy, hz, K, data }, dt]) => {
                const filter = new TestGradientAnisotropic3(xB, yB, zB, hx, hy,
                    hz, data, null, NEUMANN, PdeFilterScaleType.NONE, K);
                filter.setTimeStep(dt);
                const mHalf = filter.mHalf();

                const p = neumannPadded(xB, yB, zB, data);
                const w = xB + 2, h = yB + 2;
                const at = (px: number, py: number, pz: number) =>
                    px + w * (py + h * pz);
                const expected: number[] = [];
                for (let z = 1; z <= zB; ++z) {
                    for (let y = 1; y <= yB; ++y) {
                        for (let x = 1; x <= xB; ++x) {
                            const c = p[at(x, y, z)];
                            const uxFwd = (p[at(x + 1, y, z)] - c) / hx;
                            const uxBwd = (c - p[at(x - 1, y, z)]) / hx;
                            const uyFwd = (p[at(x, y + 1, z)] - c) / hy;
                            const uyBwd = (c - p[at(x, y - 1, z)]) / hy;
                            const uzFwd = (p[at(x, y, z + 1)] - c) / hz;
                            const uzBwd = (c - p[at(x, y, z - 1)]) / hz;

                            const dx = (a: number, b: number) =>
                                (p[a] - p[b]) / (2 * hx);
                            const dy = (a: number, b: number) =>
                                (p[a] - p[b]) / (2 * hy);
                            const dz = (a: number, b: number) =>
                                (p[a] - p[b]) / (2 * hz);

                            const duvzz = dx(at(x + 1, y, z), at(x - 1, y, z));
                            const duvpz = dx(at(x + 1, y + 1, z), at(x - 1, y + 1, z));
                            const duvmz = dx(at(x + 1, y - 1, z), at(x - 1, y - 1, z));
                            const duvzp = dx(at(x + 1, y, z + 1), at(x - 1, y, z + 1));
                            const duvzm = dx(at(x + 1, y, z - 1), at(x - 1, y, z - 1));

                            const duzvz = dy(at(x, y + 1, z), at(x, y - 1, z));
                            const dupvz = dy(at(x + 1, y + 1, z), at(x + 1, y - 1, z));
                            const dumvz = dy(at(x - 1, y + 1, z), at(x - 1, y - 1, z));
                            const duzvp = dy(at(x, y + 1, z + 1), at(x, y - 1, z + 1));
                            const duzvm = dy(at(x, y + 1, z - 1), at(x, y - 1, z - 1));

                            const duzzv = dz(at(x, y, z + 1), at(x, y, z - 1));
                            const dupzv = dz(at(x + 1, y, z + 1), at(x + 1, y, z - 1));
                            const dumzv = dz(at(x - 1, y, z + 1), at(x - 1, y, z - 1));
                            const duzpv = dz(at(x, y + 1, z + 1), at(x, y + 1, z - 1));
                            const duzmv = dz(at(x, y - 1, z + 1), at(x, y - 1, z - 1));

                            const conduct = (g: number) => Math.exp(mHalf * g);
                            const sq = (v: number) => v * v;
                            const uxCenSqr = sq(duvzz);
                            const uyCenSqr = sq(duzvz);
                            const uzCenSqr = sq(duzzv);

                            const cxp = conduct(uxCenSqr
                                + sq(0.5 * (duzvz + dupvz)) + sq(0.5 * (duzzv + dupzv)));
                            const cxm = conduct(uxCenSqr
                                + sq(0.5 * (duzvz + dumvz)) + sq(0.5 * (duzzv + dumzv)));
                            const cyp = conduct(sq(0.5 * (duvzz + duvpz))
                                + uyCenSqr + sq(0.5 * (duzzv + duzpv)));
                            const cym = conduct(sq(0.5 * (duvzz + duvmz))
                                + uyCenSqr + sq(0.5 * (duzzv + duzmv)));
                            const czp = conduct(sq(0.5 * (duvzz + duvzp))
                                + sq(0.5 * (duzvz + duzvp)) + uzCenSqr);
                            const czm = conduct(sq(0.5 * (duvzz + duvzm))
                                + sq(0.5 * (duzvz + duzvm)) + uzCenSqr);

                            expected.push(c + dt * (cxp * uxFwd - cxm * uxBwd
                                + cyp * uyFwd - cym * uyBwd
                                + czp * uzFwd - czm * uzBwd));
                        }
                    }
                }

                filter.update();
                const actual = values(filter, xB, yB, zB);
                for (let i = 0; i < expected.length; ++i) {
                    expect(Number.isFinite(actual[i])).toBe(true);
                    expectClose(actual[i], expected[i], 1e-9, 1e-9);
                }
            });
    });

    it('reduces to the isotropic Gaussian blur as K grows without bound', () => {
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xB, yB, zB, data }, dt]) => {
                const anis = new GradientAnisotropic3(xB, yB, zB, 1, 1, 1, data,
                    null, NEUMANN, PdeFilterScaleType.NONE, 1e6);
                anis.setTimeStep(dt);
                anis.update();

                const iso = new GaussianBlur3(xB, yB, zB, 1, 1, 1, data, null,
                    NEUMANN, PdeFilterScaleType.NONE);
                iso.setTimeStep(dt);
                iso.update();

                for (let z = 0; z < zB; ++z) {
                    for (let y = 0; y < yB; ++y) {
                        for (let x = 0; x < xB; ++x) {
                            expectClose(anis.getU(x, y, z), iso.getU(x, y, z),
                                1e-6, 1e-6);
                        }
                    }
                }
            });
    });

    it('is equivariant under permuting the coordinate axes', () => {
        check(config, ({ xB, yB, zB, hx, hy, hz, K, data }) => {
            const filter = new GradientAnisotropic3(xB, yB, zB, hx, hy, hz, data,
                null, NEUMANN, PdeFilterScaleType.NONE, K);
            filter.setTimeStep(0.1);
            filter.update();

            const permuted = new Array<number>(xB * yB * zB).fill(0);
            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        permuted[y + yB * (z + zB * x)] = data[x + xB * (y + yB * z)];
                    }
                }
            }
            const other = new GradientAnisotropic3(yB, zB, xB, hy, hz, hx,
                permuted, null, NEUMANN, PdeFilterScaleType.NONE, K);
            other.setTimeStep(0.1);
            other.update();

            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        expectClose(other.getU(y, z, x), filter.getU(x, y, z),
                            1e-9, 1e-9);
                    }
                }
            }
        });
    });
});
