import { describe, it, expect } from 'vitest';
import { GaussianBlur3 } from '../src/GaussianBlur3.js';
import { PdeFilterScaleType } from '../src/PdeFilter.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

describe('GaussianBlur3 verification', () => {
    // The padded image the constructor holds: interior samples shifted by the
    // data minimum (the NONE scale type still subtracts it) surrounded by the
    // Neumann ghost shell, which duplicates the nearest interior sample.
    function neumannPadded(xB: number, yB: number, zB: number,
        data: readonly number[]): number[] {
        const min = Math.min(...data);
        const max = Math.max(...data);
        const shift = (min === max ? (d: number) => 0 * d : (d: number) => d - min);
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

    const config = fc.tuple(
        fc.integer({ min: 3, max: 5 }),
        fc.integer({ min: 3, max: 5 }),
        fc.integer({ min: 3, max: 5 }),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 125, maxLength: 125 }))
        .map(([xB, yB, zB, hx, hy, hz, pool]) => ({
            xB, yB, zB, hx, hy, hz, data: pool.slice(0, xB * yB * zB)
        }));

    function values(filter: GaussianBlur3, xB: number, yB: number, zB: number): number[] {
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

    it('one step is the explicit heat step with the Neumann ghost shell', () => {
        check(fc.tuple(config, fc.constantFrom(0, 0.25, 0.5, 1)),
            ([{ xB, yB, zB, hx, hy, hz, data }, dtFraction]) => {
                const filter = new GaussianBlur3(xB, yB, zB, hx, hy, hz, data,
                    null, NEUMANN, PdeFilterScaleType.NONE);
                expectClose(filter.getMaximumTimeStep(),
                    0.5 / (1 / (hx * hx) + 1 / (hy * hy) + 1 / (hz * hz)),
                    1e-15, 1e-15);
                const dt = dtFraction * filter.getMaximumTimeStep();
                filter.setTimeStep(dt);

                const p = neumannPadded(xB, yB, zB, data);
                const w = xB + 2, h = yB + 2;
                const at = (px: number, py: number, pz: number) =>
                    px + w * (py + h * pz);
                const expected: number[] = [];
                for (let z = 1; z <= zB; ++z) {
                    for (let y = 1; y <= yB; ++y) {
                        for (let x = 1; x <= xB; ++x) {
                            const uxx = (p[at(x + 1, y, z)] - 2 * p[at(x, y, z)]
                                + p[at(x - 1, y, z)]) / (hx * hx);
                            const uyy = (p[at(x, y + 1, z)] - 2 * p[at(x, y, z)]
                                + p[at(x, y - 1, z)]) / (hy * hy);
                            const uzz = (p[at(x, y, z + 1)] - 2 * p[at(x, y, z)]
                                + p[at(x, y, z - 1)]) / (hz * hz);
                            expected.push(p[at(x, y, z)] + dt * (uxx + uyy + uzz));
                        }
                    }
                }

                filter.update();
                const actual = values(filter, xB, yB, zB);
                for (let i = 0; i < expected.length; ++i) {
                    expectClose(actual[i], expected[i], 1e-10, 1e-10);
                }
            });
    });

    it('preserves the total mass under Neumann conditions', () => {
        check(config, ({ xB, yB, zB, hx, hy, hz, data }) => {
            const filter = new GaussianBlur3(xB, yB, zB, hx, hy, hz, data, null,
                NEUMANN, PdeFilterScaleType.NONE);
            filter.setTimeStep(filter.getMaximumTimeStep());
            const before = values(filter, xB, yB, zB).reduce((a, b) => a + b, 0);
            filter.update();
            const after = values(filter, xB, yB, zB).reduce((a, b) => a + b, 0);
            expectClose(after, before, 1e-9, 1e-9);
        });
    });

    it('leaves an affine image unchanged away from the border', () => {
        const affine = fc.tuple(
            fc.integer({ min: 4, max: 6 }),
            fc.constantFrom(0.5, 1, 2),
            fc.array(fc.integer({ min: -5, max: 5 }), { minLength: 4, maxLength: 4 }));
        check(affine, ([bound, h, k]) => {
            const [a, b, c, d] = k;
            const data: number[] = [];
            for (let z = 0; z < bound; ++z) {
                for (let y = 0; y < bound; ++y) {
                    for (let x = 0; x < bound; ++x) {
                        data.push(a + b * x + c * y + d * z);
                    }
                }
            }
            const filter = new GaussianBlur3(bound, bound, bound, h, h, h, data,
                null, NEUMANN, PdeFilterScaleType.NONE);
            const before = values(filter, bound, bound, bound);
            filter.setTimeStep(filter.getMaximumTimeStep());
            filter.update();
            const after = values(filter, bound, bound, bound);
            for (let z = 1; z + 1 < bound; ++z) {
                for (let y = 1; y + 1 < bound; ++y) {
                    for (let x = 1; x + 1 < bound; ++x) {
                        const i = x + bound * (y + bound * z);
                        expectClose(after[i], before[i], 1e-9, 1e-9);
                    }
                }
            }
        });
    });

    it('is equivariant under permuting the coordinate axes', () => {
        // Relabel (x,y,z) as (y,z,x) together with the spacings; the result
        // must be relabelled the same way. A swapped axis in lookUp7 or in
        // the stencil breaks this.
        check(config, ({ xB, yB, zB, hx, hy, hz, data }) => {
            const filter = new GaussianBlur3(xB, yB, zB, hx, hy, hz, data, null,
                NEUMANN, PdeFilterScaleType.NONE);
            filter.setTimeStep(0.5 * filter.getMaximumTimeStep());
            filter.update();

            // The permuted image has bounds (yB, zB, xB) and holds
            // permuted[y + yB*(z + zB*x)] = data[x + xB*(y + yB*z)].
            const permuted = new Array<number>(xB * yB * zB).fill(0);
            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        permuted[y + yB * (z + zB * x)] = data[x + xB * (y + yB * z)];
                    }
                }
            }
            const other = new GaussianBlur3(yB, zB, xB, hy, hz, hx, permuted,
                null, NEUMANN, PdeFilterScaleType.NONE);
            other.setTimeStep(0.5 * other.getMaximumTimeStep());
            other.update();

            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        expectClose(other.getU(y, z, x), filter.getU(x, y, z),
                            1e-10, 1e-10);
                    }
                }
            }
        });
    });
});
