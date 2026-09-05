import { describe, it, expect } from 'vitest';
import { GradientAnisotropic2 } from '../src/GradientAnisotropic2.js';
import { GaussianBlur2 } from '../src/GaussianBlur2.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

describe('GradientAnisotropic2 verification', () => {
    // The padded image the constructor holds under Neumann conditions and the
    // NONE scale type (which still subtracts the data minimum).
    function neumannPadded(xB: number, yB: number, data: readonly number[]): number[] {
        const min = Math.min(...data);
        const max = Math.max(...data);
        const shift = (min === max ? () => 0 : (d: number) => d - min);
        const w = xB + 2;
        const p = new Array<number>(w * (yB + 2)).fill(0);
        const at = (px: number, py: number) => px + w * py;
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                p[at(x + 1, y + 1)] = shift(data[x + xB * y]);
            }
        }
        for (let py = 0; py < yB + 2; ++py) {
            for (let px = 0; px < w; ++px) {
                if (px === 0 || px === xB + 1 || py === 0 || py === yB + 1) {
                    p[at(px, py)] = p[at(
                        Math.min(Math.max(px, 1), xB),
                        Math.min(Math.max(py, 1), yB))];
                }
            }
        }
        return p;
    }

    function values(filter: GradientAnisotropic2, xB: number, yB: number): number[] {
        const u: number[] = [];
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                u.push(filter.getU(x, y));
            }
        }
        return u;
    }

    // Nonconstant images only: a constant image has zero average squared
    // gradient, so mParameter is infinite and every conductance is NaN (an
    // upstream quirk with no guard, preserved by the port).
    const config = fc.tuple(
        fc.integer({ min: 3, max: 6 }),
        fc.integer({ min: 3, max: 6 }),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2, 5),
        fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 36, maxLength: 36 }))
        .map(([xB, yB, hx, hy, K, pool]) => ({
            xB, yB, hx, hy, K, data: pool.slice(0, xB * yB)
        }))
        .filter(({ data }) => new Set(data).size > 1);

    it('the conductance parameter averages over the unpadded coordinates', () => {
        // The average must be taken over 0 <= x < xBound and 0 <= y < yBound,
        // the domain GetUx/GetUy document, and divided by mQuantity. Upstream
        // feeds padded coordinates to the unpadded accessors, which shifts
        // the window by one pixel and reads past the padded buffer; the port
        // fixes that (upstream issue #122).
        check(config, ({ xB, yB, hx, hy, K, data }) => {
            const filter = new TestGradientAnisotropic2(xB, yB, hx, hy, data,
                null, NEUMANN, PdeFilterScaleType.NONE, K);
            let sum = 0;
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    const ux = filter.getUx(x, y);
                    const uy = filter.getUy(x, y);
                    sum += ux * ux + uy * uy;
                }
            }
            const avg = sum / (xB * yB);
            expect(Number.isFinite(filter.parameter())).toBe(true);
            expectClose(filter.parameter(), 1 / (K * K * avg), 1e-9, 1e-9);
            expectClose(filter.mHalf(), -0.5 * filter.parameter(), 1e-12, 1e-12);
        });
    });

    it('one step reproduces the upstream update from the padded image', () => {
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xB, yB, hx, hy, K, data }, dt]) => {
                const filter = new TestGradientAnisotropic2(xB, yB, hx, hy, data,
                    null, NEUMANN, PdeFilterScaleType.NONE, K);
                filter.setTimeStep(dt);
                const mHalf = filter.mHalf();

                const p = neumannPadded(xB, yB, data);
                const w = xB + 2;
                const at = (px: number, py: number) => px + w * py;
                const expected: number[] = [];
                for (let y = 1; y <= yB; ++y) {
                    for (let x = 1; x <= xB; ++x) {
                        const uzz = p[at(x, y)];
                        const umz = p[at(x - 1, y)], upz = p[at(x + 1, y)];
                        const uzm = p[at(x, y - 1)], uzp = p[at(x, y + 1)];
                        const umm = p[at(x - 1, y - 1)], upm = p[at(x + 1, y - 1)];
                        const ump = p[at(x - 1, y + 1)], upp = p[at(x + 1, y + 1)];

                        const uxFwd = (upz - uzz) / hx;
                        const uxBwd = (uzz - umz) / hx;
                        const uyFwd = (uzp - uzz) / hy;
                        const uyBwd = (uzz - uzm) / hy;

                        const uxCenM = (upm - umm) / (2 * hx);
                        const uxCenZ = (upz - umz) / (2 * hx);
                        const uxCenP = (upp - ump) / (2 * hx);
                        const uyCenM = (ump - umm) / (2 * hy);
                        const uyCenZ = (uzp - uzm) / (2 * hy);
                        const uyCenP = (upp - upm) / (2 * hy);

                        const conduct = (g: number) => Math.exp(mHalf * g);
                        const uyEstP = 0.5 * (uyCenZ + uyCenP);
                        const cxp = conduct(uxCenZ * uxCenZ + uyEstP * uyEstP);
                        const uyEstM = 0.5 * (uyCenZ + uyCenM);
                        const cxm = conduct(uxCenZ * uxCenZ + uyEstM * uyEstM);
                        const uxEstP = 0.5 * (uxCenZ + uxCenP);
                        const cyp = conduct(uyCenZ * uyCenZ + uxEstP * uxEstP);
                        const uxEstM = 0.5 * (uxCenZ + uxCenM);
                        const cym = conduct(uyCenZ * uyCenZ + uxEstM * uxEstM);

                        expected.push(uzz + dt * (cxp * uxFwd - cxm * uxBwd
                            + cyp * uyFwd - cym * uyBwd));
                    }
                }

                filter.update();
                const actual = values(filter, xB, yB);
                for (let i = 0; i < expected.length; ++i) {
                    expect(Number.isFinite(actual[i])).toBe(true);
                    expectClose(actual[i], expected[i], 1e-9, 1e-9);
                }
            });
    });

    it('reduces to the isotropic Gaussian blur as K grows without bound', () => {
        // A huge K makes every conductance 1, and with unit spacing the
        // one-sided terms telescope into the ordinary discrete Laplacian, so
        // the step becomes the GaussianBlur2 step.
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xB, yB, data }, dt]) => {
                const anis = new GradientAnisotropic2(xB, yB, 1, 1, data, null,
                    NEUMANN, PdeFilterScaleType.NONE, 1e6);
                anis.setTimeStep(dt);
                anis.update();

                const iso = new GaussianBlur2(xB, yB, 1, 1, data, null, NEUMANN,
                    PdeFilterScaleType.NONE);
                iso.setTimeStep(dt);
                iso.update();

                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        expectClose(anis.getU(x, y), iso.getU(x, y), 1e-6, 1e-6);
                    }
                }
            });
    });

    it('is equivariant under reflecting the x axis', () => {
        check(config, ({ xB, yB, hx, hy, K, data }) => {
            const filter = new GradientAnisotropic2(xB, yB, hx, hy, data, null,
                NEUMANN, PdeFilterScaleType.NONE, K);
            filter.setTimeStep(0.1);
            filter.update();

            const mirrored = new Array<number>(xB * yB).fill(0);
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    mirrored[(xB - 1 - x) + xB * y] = data[x + xB * y];
                }
            }
            const other = new GradientAnisotropic2(xB, yB, hx, hy, mirrored, null,
                NEUMANN, PdeFilterScaleType.NONE, K);
            other.setTimeStep(0.1);
            other.update();

            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    expectClose(other.getU(xB - 1 - x, y), filter.getU(x, y),
                        1e-9, 1e-9);
                }
            }
        });
    });
});
