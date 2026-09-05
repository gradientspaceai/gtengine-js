import { describe, it, expect } from 'vitest';
import { PdeFilterScaleType } from '../src/PdeFilter.js';
import { PdeFilter1 } from '../src/PdeFilter1.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

// A concrete 1D filter solving the heat equation u_t = u_xx with an explicit
// Euler step, which is the 1D analogue of the upstream GaussianBlur2 filter.
class HeatFilter1 extends PdeFilter1 {
    constructor(xBound: number, xSpacing: number, data: ArrayLike<number>,
        mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType = PdeFilterScaleType.NONE) {
        super(xBound, xSpacing, data, mask, borderValue, scaleType);
    }

    protected override onUpdateSingle(x: number): void {
        this.lookUp3(x);
        const uxx = this.mInvDxDx * (this.mUp - 2 * this.mUz + this.mUm);
        this.mBuffer[this.mDst][x] = this.mUz + this.mTimeStep * uxx;
    }

    // Test access to the padded source buffer and the neighborhood cache.
    padded(): number[] {
        return Array.from(this.mBuffer[this.mSrc]);
    }

    neighborhood(): [number, number, number] {
        return [this.mUm, this.mUz, this.mUp];
    }

    lookUp(x: number): void {
        this.lookUp3(x);
    }

    values(): number[] {
        const u: number[] = [];
        for (let x = 0; x < this.mXBound; ++x) {
            u.push(this.getU(x));
        }
        return u;
    }
}

const NEUMANN = Number.MAX_VALUE;

function closeTo(actual: number[], expected: number[]): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(actual[i]).toBeCloseTo(expected[i], 12);
    }
}

describe('PdeFilter1', () => {
    it('exposes the image parameters', () => {
        const filter = new HeatFilter1(4, 0.5, [0, 1, 2, 3], null, NEUMANN);
        expect(filter.getXBound()).toBe(4);
        expect(filter.getXSpacing()).toBe(0.5);
        expect(filter.getQuantity()).toBe(4);
        expect(filter.getBorderValue()).toBe(NEUMANN);
    });

    it('duplicates the end samples for a Neumann image border', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 3, 7], null, NEUMANN);
        // The padded buffer is [u(0), u(0), u(1), u(2), u(3), u(3)].
        expect(filter.padded()).toEqual([0, 0, 1, 3, 7, 7]);
        expect(filter.values()).toEqual([0, 1, 3, 7]);
    });

    it('assigns the border value for a Dirichlet image border', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 2, 3], null, -5);
        expect(filter.padded()).toEqual([-5, 0, 1, 2, 3, -5]);
    });

    it('estimates first and second derivatives with a Neumann border', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 3, 7], null, NEUMANN);
        // F = [0,0,1,3,7,7]; Ux(x) = (F[x+2]-F[x])/2, Uxx(x) = F[x+2]-2F[x+1]+F[x].
        closeTo([0, 1, 2, 3].map(x => filter.getUx(x)), [0.5, 1.5, 3, 2]);
        closeTo([0, 1, 2, 3].map(x => filter.getUxx(x)), [1, 1, 2, -4]);
    });

    it('scales the derivatives by the sample spacing', () => {
        const filter = new HeatFilter1(4, 2, [0, 1, 3, 7], null, NEUMANN);
        // dx = 2 gives 1/(2*dx) = 0.25 and 1/(dx*dx) = 0.25.
        closeTo([0, 1, 2, 3].map(x => filter.getUx(x)), [0.25, 0.75, 1.5, 1]);
        closeTo([0, 1, 2, 3].map(x => filter.getUxx(x)), [0.25, 0.25, 0.5, -1]);
    });

    it('caches the 3-tuple neighborhood in lookUp3', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 3, 7], null, NEUMANN);
        // Padded coordinate 2 is image coordinate 1.
        filter.lookUp(2);
        expect(filter.neighborhood()).toEqual([0, 1, 3]);
        filter.lookUp(1);
        expect(filter.neighborhood()).toEqual([0, 0, 1]);
    });

    it('takes Dirichlet heat-equation steps matching hand computation', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 2, 3], null, 0);
        filter.setTimeStep(0.25);
        // F = [0, 0,1,2,3, 0]; Uxx = [1, 0, 0, -4].
        filter.update();
        expect(filter.padded()).toEqual([0, 0.25, 1, 2, 2, 0]);
        // F = [0, 0.25,1,2,2, 0]; Uxx = [0.5, 0.25, -1, -2].
        filter.update();
        expect(filter.padded()).toEqual([0, 0.375, 1.0625, 1.75, 1.5, 0]);
    });

    it('takes Neumann heat-equation steps matching hand computation', () => {
        const filter = new HeatFilter1(4, 1, [0, 1, 3, 7], null, NEUMANN);
        filter.setTimeStep(0.1);
        filter.update();
        // Uxx = [1, 1, 2, -4] on F = [0,0,1,3,7,7].
        closeTo(filter.values(), [0.1, 1.1, 3.2, 6.6]);
        // Upstream never re-runs AssignNeumannImageBorder after construction,
        // so the padded border keeps its original values.
        closeTo(filter.padded(), [0, 0.1, 1.1, 3.2, 6.6, 7]);
        filter.update();
        // Uxx = [0.9, 1.1, 1.3, -3] on F = [0, 0.1,1.1,3.2,6.6, 7].
        closeTo(filter.values(), [0.19, 1.21, 3.33, 6.3]);
    });

    it('applies the requested scale type to the stored samples', () => {
        const unit = new HeatFilter1(3, 1, [2, 6, 4], null, 0, PdeFilterScaleType.UNIT);
        closeTo(unit.values(), [0, 1, 0.5]);

        const symmetric = new HeatFilter1(3, 1, [2, 6, 4], null, 0,
            PdeFilterScaleType.SYMMETRIC);
        closeTo(symmetric.values(), [-1, 1, 0]);

        // ScaleType.NONE has offset 0 and scale 1 but leaves mMin at the data
        // minimum, so the samples are shifted by -min. Upstream behavior; see
        // the PR notes on PdeFilter.
        const none = new HeatFilter1(3, 1, [2, 6, 4], null, 0, PdeFilterScaleType.NONE);
        closeTo(none.values(), [0, 4, 2]);
    });

    it('averages the unmasked neighbors for a Neumann mask border', () => {
        const filter = new HeatFilter1(5, 1, [0, 10, 20, 30, 40], [0, 1, 1, 1, 0],
            NEUMANN);
        // The masked-out end samples take the average of their unmasked
        // neighbors: sample 0 becomes 10 and sample 4 becomes 30.
        expect(filter.padded()).toEqual([0, 10, 10, 20, 30, 30, 40]);
        expect([0, 1, 2, 3, 4].map(x => filter.getMask(x))).toEqual([0, 1, 1, 1, 0]);
    });

    it('updates only unmasked samples and refreshes the Neumann mask border', () => {
        const filter = new HeatFilter1(5, 1, [0, 10, 20, 30, 40], [0, 1, 1, 1, 0],
            NEUMANN);
        filter.setTimeStep(0.5);
        filter.update();
        // Uxx at padded 2,3,4 is 10, 0, -10 on F = [0,10,10,20,30,30,40].
        closeTo(filter.values(), [10, 15, 20, 25, 30]);
        filter.update();
        // onPreUpdate first pulls the masked ends to 15 and 25 on
        // F = [0,15,15,20,25,25,40], where Uxx at padded 2,3,4 is 5, 0, -5.
        closeTo(filter.values(), [15, 17.5, 20, 22.5, 25]);
    });

    it('assigns the border value on a Dirichlet mask border', () => {
        const filter = new HeatFilter1(5, 1, [0, 10, 20, 30, 40], [0, 1, 1, 1, 0], 0);
        // The masked-out samples adjacent to the mask take the border value,
        // replacing 0 and 40.
        expect(filter.padded()).toEqual([0, 0, 10, 20, 30, 0, 0]);
        filter.setTimeStep(0.5);
        filter.update();
        // Uxx at padded 2,3,4 is 0, 0, -40; the Dirichlet border is not
        // recomputed by onPreUpdate.
        closeTo(filter.values(), [0, 10, 20, 10, 0]);
    });

    it('leaves a constant image unchanged (Neumann steady state)', () => {
        const filter = new HeatFilter1(6, 1, [3, 3, 3, 3, 3, 3], null, NEUMANN);
        filter.setTimeStep(0.4);
        for (let i = 0; i < 5; ++i) {
            filter.update();
        }
        // For constant data the base class takes the min == max branch, which
        // sets offset 0 and scale 1 but leaves mMin at the data value, so the
        // stored samples are the data shifted by -min, i.e. all zero (see the
        // PR notes on PdeFilter). A constant image is a steady state of the
        // heat equation with Neumann conditions, so the samples do not change.
        closeTo(filter.values(), [0, 0, 0, 0, 0, 0]);
    });

    it('agrees with an independent explicit heat-equation solver (Dirichlet)', () => {
        const data = [0, 2, 5, 9, 4, 1, 7, 3];
        const dx = 0.5;
        const dt = 0.05;
        const border = 0;
        const filter = new HeatFilter1(data.length, dx, data, null, border);
        filter.setTimeStep(dt);

        // Reference: an explicit-Euler heat solver on a padded array whose
        // border is held at the Dirichlet value.
        let u = [border, ...data, border];
        for (let step = 0; step < 6; ++step) {
            const next = u.slice();
            for (let i = 1; i <= data.length; ++i) {
                next[i] = u[i] + dt * (u[i + 1] - 2 * u[i] + u[i - 1]) / (dx * dx);
            }
            u = next;
            filter.update();
        }
        closeTo(filter.values(), u.slice(1, data.length + 1));
    });
});

describe('PdeFilter1 verification', () => {
    const scaleTypes = [
        PdeFilterScaleType.NONE,
        PdeFilterScaleType.UNIT,
        PdeFilterScaleType.SYMMETRIC,
        PdeFilterScaleType.PRESERVE_ZERO
    ];

    // An independent implementation of the documented PdeFilter scaling
    // (upstream PdeFilter.h comments), written from the definitions rather
    // than from the ported code. Note the NONE branch: upstream still
    // subtracts the minimum, so "as is" actually shifts the data.
    function referenceScale(data: readonly number[],
        scaleType: PdeFilterScaleType): number[] {
        const min = Math.min(...data);
        const max = Math.max(...data);
        if (min === max) {
            return data.map(() => 0);
        }
        switch (scaleType) {
            case PdeFilterScaleType.UNIT:
                return data.map(d => (d - min) / (max - min));
            case PdeFilterScaleType.SYMMETRIC:
                return data.map(d => -1 + (2 / (max - min)) * (d - min));
            case PdeFilterScaleType.PRESERVE_ZERO: {
                const scale = (max >= -min ? 1 / max : -1 / min);
                return data.map(d => d * scale);
            }
            default:
                return data.map(d => d - min);
        }
    }

    // The padded (xBound+2)-element buffer that the constructor builds.
    function referencePadded(data: readonly number[], borderValue: number,
        scaleType: PdeFilterScaleType): number[] {
        const scaled = referenceScale(data, scaleType);
        const n = data.length;
        const p = new Array<number>(n + 2).fill(0);
        for (let x = 0; x < n; ++x) {
            p[x + 1] = scaled[x];
        }
        if (borderValue !== NEUMANN) {
            p[0] = borderValue;
            p[n + 1] = borderValue;
        } else {
            p[0] = p[1];
            p[n + 1] = p[n];
        }
        return p;
    }

    const config = fc.tuple(
        fc.integer({ min: 3, max: 8 }),
        fc.constantFrom(0.25, 0.5, 1, 2),
        fc.constantFrom(...scaleTypes),
        fc.constantFrom(NEUMANN, 0, -3, 7),
        fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 8, maxLength: 8 }))
        .map(([xBound, h, scaleType, borderValue, pool]) => ({
            xBound, h, scaleType, borderValue,
            data: pool.slice(0, xBound)
        }));

    it('builds the padded buffer the border conditions prescribe', () => {
        check(config, ({ xBound, h, scaleType, borderValue, data }) => {
            const filter = new HeatFilter1(xBound, h, data, null, borderValue,
                scaleType);
            const expected = referencePadded(data, borderValue, scaleType);
            const actual = filter.padded();
            expect(actual.length).toBe(expected.length);
            for (let i = 0; i < expected.length; ++i) {
                expectClose(actual[i], expected[i], 1e-12, 1e-12);
            }
        });
    });

    it('getUx and getUxx are exact for a quadratic sample function', () => {
        // Central differences reproduce the first and second derivatives of a
        // quadratic exactly, so any wrong stencil coefficient or padding
        // offset shows up immediately. Only interior samples are used; the
        // border samples see the padded ghost values instead of the
        // quadratic.
        const quadratic = fc.tuple(
            fc.integer({ min: 4, max: 9 }),
            fc.constantFrom(0.25, 0.5, 1, 2),
            fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 }));
        check(quadratic, ([xBound, h, a, b, c]) => {
            const f = (x: number) => a + b * (x * h) + c * (x * h) * (x * h);
            const data: number[] = [];
            for (let x = 0; x < xBound; ++x) {
                data.push(f(x));
            }
            // PRESERVE_ZERO keeps derivatives up to the known scale factor;
            // NONE only shifts by a constant, so derivatives are unchanged.
            const filter = new HeatFilter1(xBound, h, data, null, NEUMANN,
                PdeFilterScaleType.NONE);
            for (let x = 1; x + 1 < xBound; ++x) {
                const scale = Math.abs(a) + Math.abs(b) + Math.abs(c) + 1;
                expectClose(filter.getUx(x), b + 2 * c * (x * h),
                    1e-9 * scale / h, 1e-9);
                expectClose(filter.getUxx(x), 2 * c, 1e-9 * scale / (h * h), 1e-9);
            }
        });
    });

    it('one update step agrees with an independent explicit solver', () => {
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xBound, h, scaleType, borderValue, data }, dtFactor]) => {
                const filter = new HeatFilter1(xBound, h, data, null, borderValue,
                    scaleType);
                const dt = dtFactor * h * h;
                filter.setTimeStep(dt);

                const p = referencePadded(data, borderValue, scaleType);
                const expected: number[] = [];
                for (let x = 1; x <= xBound; ++x) {
                    expected.push(p[x] + (dt / (h * h))
                        * (p[x + 1] - 2 * p[x] + p[x - 1]));
                }

                filter.update();
                const actual = filter.values();
                for (let x = 0; x < xBound; ++x) {
                    expectClose(actual[x], expected[x], 1e-10, 1e-10);
                }
            });
    });

    it('a mask leaves masked samples out of the update', () => {
        const masked = fc.tuple(
            fc.integer({ min: 4, max: 8 }),
            fc.array(fc.integer({ min: -5, max: 5 }), { minLength: 8, maxLength: 8 }),
            fc.array(fc.integer({ min: 0, max: 1 }), { minLength: 8, maxLength: 8 }));
        check(masked, ([xBound, pool, maskPool]) => {
            const data = pool.slice(0, xBound);
            const mask = maskPool.slice(0, xBound);
            const filter = new HeatFilter1(xBound, 1, data, mask, NEUMANN,
                PdeFilterScaleType.NONE);
            // getMask must report exactly the caller's mask on the unpadded
            // samples; the padded border is zeroed but is not visible here.
            for (let x = 0; x < xBound; ++x) {
                expect(filter.getMask(x)).toBe(mask[x]);
            }

            // A masked-out sample is never handed to onUpdateSingle, so the
            // destination buffer keeps whatever the constructor left there.
            const before = filter.values();
            filter.setTimeStep(0);
            filter.update();
            const after = filter.values();
            for (let x = 0; x < xBound; ++x) {
                if (mask[x] !== 0) {
                    // A zero time step makes the update the identity.
                    expectClose(after[x], before[x], 1e-12, 1e-12);
                }
            }
        });
    });
});
