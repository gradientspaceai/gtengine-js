import { describe, it, expect } from 'vitest';
import { PdeFilterScaleType } from '../src/PdeFilter.js';
import { PdeFilter2 } from '../src/PdeFilter2.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

// Concrete subclass solving the linear heat equation u_t = u_xx + u_yy with
// an explicit Euler step, the canonical use of the PdeFilter2 plumbing.
class HeatFilter2 extends PdeFilter2 {
    constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType = PdeFilterScaleType.NONE) {
        super(xBound, yBound, xSpacing, ySpacing, data, mask, borderValue, scaleType);
    }

    protected override onUpdateSingle(x: number, y: number): void {
        this.lookUp5(x, y);
        const uxx = this.mInvDxDx * (this.mUpz - 2 * this.mUzz + this.mUmz);
        const uyy = this.mInvDyDy * (this.mUzp - 2 * this.mUzz + this.mUzm);
        this.mBuffer[this.mDst].set(x, y, this.mUzz + this.mTimeStep * (uxx + uyy));
    }

    // The padded source buffer in row-major order, (xBound+2)-by-(yBound+2).
    padded(): number[] {
        const rows: number[] = [];
        for (let y = 0; y < this.mYBound + 2; ++y) {
            for (let x = 0; x < this.mXBound + 2; ++x) {
                rows.push(this.mBuffer[this.mSrc].get(x, y));
            }
        }
        return rows;
    }

    // The unpadded image in row-major order.
    values(): number[] {
        const u: number[] = [];
        for (let y = 0; y < this.mYBound; ++y) {
            for (let x = 0; x < this.mXBound; ++x) {
                u.push(this.getU(x, y));
            }
        }
        return u;
    }

    maskValues(): number[] {
        const m: number[] = [];
        for (let y = 0; y < this.mYBound; ++y) {
            for (let x = 0; x < this.mXBound; ++x) {
                m.push(this.getMask(x, y));
            }
        }
        return m;
    }

    // Expose the neighborhood caches for the padded pixel (x,y).
    lookUp5At(x: number, y: number): number[] {
        this.lookUp5(x, y);
        return [this.mUzm, this.mUmz, this.mUzz, this.mUpz, this.mUzp];
    }

    lookUp9At(x: number, y: number): number[] {
        this.lookUp9(x, y);
        return [
            this.mUmm, this.mUzm, this.mUpm,
            this.mUmz, this.mUzz, this.mUpz,
            this.mUmp, this.mUzp, this.mUpp
        ];
    }
}

// Subclass that records the pixels visited by onUpdate, in order, and copies
// the source through to the destination unchanged.
class RecordingFilter2 extends PdeFilter2 {
    visited: [number, number][] = [];

    constructor(xBound: number, yBound: number, data: ArrayLike<number>,
        mask: ArrayLike<number> | null, borderValue: number) {
        super(xBound, yBound, 1, 1, data, mask, borderValue, PdeFilterScaleType.NONE);
    }

    protected override onUpdateSingle(x: number, y: number): void {
        this.visited.push([x, y]);
        this.lookUp5(x, y);
        this.mBuffer[this.mDst].set(x, y, this.mUzz);
    }

    src(): number {
        return this.mSrc;
    }

    dst(): number {
        return this.mDst;
    }
}

const NEUMANN = Number.MAX_VALUE;

// Build data[i] = f(x,y) for i = x + xBound * y.
function makeImage(xB: number, yB: number, f: (x: number, y: number) => number): number[] {
    const data = new Array<number>(xB * yB);
    let i = 0;
    for (let y = 0; y < yB; ++y) {
        for (let x = 0; x < xB; ++x, ++i) {
            data[i] = f(x, y);
        }
    }
    return data;
}

function closeTo(actual: number[], expected: number[], digits = 12): void {
    expect(actual.length).toBe(expected.length);
    for (let i = 0; i < expected.length; ++i) {
        expect(actual[i]).toBeCloseTo(expected[i], digits);
    }
}

describe('PdeFilter2', () => {
    it('provides member access', () => {
        const data = makeImage(4, 3, (x, y) => x + y);
        const filter = new HeatFilter2(4, 3, 0.5, 2, data, null, NEUMANN);
        expect(filter.getXBound()).toBe(4);
        expect(filter.getYBound()).toBe(3);
        expect(filter.getXSpacing()).toBe(0.5);
        expect(filter.getYSpacing()).toBe(2);
        expect(filter.getQuantity()).toBe(12);
        expect(filter.getBorderValue()).toBe(NEUMANN);
        expect(filter.getScaleType()).toBe(PdeFilterScaleType.NONE);
    });

    it('stores the data shifted by the minimum with NONE scaling', () => {
        // The base class stores offset + (d - min) * scale; for NONE the
        // offset is 0 and the scale is 1, so the stored image is d - min
        // (upstream semantics).
        const f = (x: number, y: number) => 3 + x + 10 * y;
        const data = makeImage(4, 3, f);
        const filter = new HeatFilter2(4, 3, 1, 1, data, null, NEUMANN);
        for (let y = 0; y < 3; ++y) {
            for (let x = 0; x < 4; ++x) {
                expect(filter.getU(x, y)).toBe(f(x, y) - 3);
            }
        }
    });

    it('applies UNIT, SYMMETRIC and PRESERVE_ZERO scaling', () => {
        const data = [2, 6, 4, 10];  // min 2, max 10

        const unit = new HeatFilter2(2, 2, 1, 1, data, null, 0, PdeFilterScaleType.UNIT);
        closeTo(unit.values(), [0, 0.5, 0.25, 1]);

        const symmetric = new HeatFilter2(2, 2, 1, 1, data, null, 0,
            PdeFilterScaleType.SYMMETRIC);
        closeTo(symmetric.values(), [-1, 0, -0.5, 1]);

        // max = 10 > -min = -2, so d' = d / max and the zero level is kept.
        const preserve = new HeatFilter2(2, 2, 1, 1, data, null, 0,
            PdeFilterScaleType.PRESERVE_ZERO);
        closeTo(preserve.values(), [0.2, 0.6, 0.4, 1]);
    });

    it('duplicates the boundary samples for a Neumann image border', () => {
        const data = [1, 2, 3, 4, 5, 6];  // 3-by-2, min 1
        const filter = new HeatFilter2(3, 2, 1, 1, data, null, NEUMANN);
        // Stored image (shifted by the minimum):
        //     0 1 2
        //     3 4 5
        // The padded image duplicates the nearest interior sample on the
        // edges and the nearest corner sample at the corners.
        expect(filter.padded()).toEqual([
            0, 0, 1, 2, 2,
            0, 0, 1, 2, 2,
            3, 3, 4, 5, 5,
            3, 3, 4, 5, 5
        ]);
    });

    it('assigns the border value for a Dirichlet image border', () => {
        const data = [1, 2, 3, 4, 5, 6];
        const filter = new HeatFilter2(3, 2, 1, 1, data, null, -9);
        expect(filter.padded()).toEqual([
            -9, -9, -9, -9, -9,
            -9, 0, 1, 2, -9,
            -9, 3, 4, 5, -9,
            -9, -9, -9, -9, -9
        ]);
    });

    it('estimates the derivatives of a quadratic exactly', () => {
        // u(i,j) = F(i*dx, j*dy) with
        //   F(X,Y) = 1 + 2X + 3Y + 4X^2 + 5XY + 6Y^2,
        // whose central differences are exact:
        //   Fx = 2 + 8X + 5Y, Fy = 3 + 5X + 12Y, Fxx = 8, Fxy = 5, Fyy = 12.
        const dx = 0.5, dy = 2;
        const xB = 5, yB = 4;
        const F = (X: number, Y: number) =>
            1 + 2 * X + 3 * Y + 4 * X * X + 5 * X * Y + 6 * Y * Y;
        const data = makeImage(xB, yB, (x, y) => F(x * dx, y * dy));
        const filter = new HeatFilter2(xB, yB, dx, dy, data, null, NEUMANN);

        // Only interior pixels avoid the padded border.
        for (let y = 1; y < yB - 1; ++y) {
            for (let x = 1; x < xB - 1; ++x) {
                const X = x * dx, Y = y * dy;
                expect(filter.getUx(x, y)).toBeCloseTo(2 + 8 * X + 5 * Y, 10);
                expect(filter.getUy(x, y)).toBeCloseTo(3 + 5 * X + 12 * Y, 10);
                expect(filter.getUxx(x, y)).toBeCloseTo(8, 10);
                expect(filter.getUxy(x, y)).toBeCloseTo(5, 10);
                expect(filter.getUyy(x, y)).toBeCloseTo(12, 10);
            }
        }
    });

    it('estimates derivatives from the padded border on the image boundary', () => {
        // With a Neumann border the duplicated samples give a zero-valued
        // one-sided estimate at the boundary of a linear image.
        const data = makeImage(3, 3, (x, y) => x + 2 * y);
        const filter = new HeatFilter2(3, 3, 1, 1, data, null, NEUMANN);
        // At x = 0: (u(1,y) - u(0,y)) / 2 = 0.5, not 1.
        expect(filter.getUx(0, 1)).toBeCloseTo(0.5, 12);
        expect(filter.getUx(1, 1)).toBeCloseTo(1, 12);
        expect(filter.getUy(1, 0)).toBeCloseTo(1, 12);
        expect(filter.getUy(1, 1)).toBeCloseTo(2, 12);
        // The duplicated border makes the second derivative nonzero there.
        expect(filter.getUxx(0, 1)).toBeCloseTo(1, 12);
        expect(filter.getUxx(1, 1)).toBeCloseTo(0, 12);
    });

    it('caches the 5- and 9-pixel neighborhoods', () => {
        // Stored image (min 0):
        //     0 1 2
        //     3 4 5
        //     6 7 8
        const data = makeImage(3, 3, (x, y) => x + 3 * y);
        const filter = new HeatFilter2(3, 3, 1, 1, data, null, NEUMANN);
        // Padded (2,2) is image (1,1), the center sample.
        expect(filter.lookUp5At(2, 2)).toEqual([1, 3, 4, 5, 7]);
        expect(filter.lookUp9At(2, 2)).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8]);
        // Padded (1,1) is image (0,0); its neighbors come from the border.
        expect(filter.lookUp9At(1, 1)).toEqual([0, 0, 1, 0, 0, 1, 3, 3, 4]);
    });

    it('takes a heat-equation step matching hand computation', () => {
        // Stored image (min 0):
        //     0 1 2
        //     3 4 5
        //     6 7 8
        // Dirichlet border 0, dx = dy = 1, dt = 0.1.
        const data = makeImage(3, 3, (x, y) => x + 3 * y);
        const filter = new HeatFilter2(3, 3, 1, 1, data, null, 0);
        filter.setTimeStep(0.1);
        filter.update();
        // Laplacian of the padded image at each pixel:
        //   (0,0): 1 + 0 - 4*0 + 3 + 0 = 4      (0,1): 4 + 3 - ...
        const lap = (x: number, y: number): number => {
            const u = (i: number, j: number): number =>
                (0 <= i && i < 3 && 0 <= j && j < 3 ? i + 3 * j : 0);
            return u(x + 1, y) + u(x - 1, y) + u(x, y + 1) + u(x, y - 1) - 4 * u(x, y);
        };
        const expected: number[] = [];
        for (let y = 0; y < 3; ++y) {
            for (let x = 0; x < 3; ++x) {
                expected.push(x + 3 * y + 0.1 * lap(x, y));
            }
        }
        closeTo(filter.values(), expected);
    });

    it('agrees with an independent explicit heat solver (Dirichlet)', () => {
        const xB = 6, yB = 5;
        const dx = 0.5, dy = 0.75, dt = 0.01, border = 0;
        // A deterministic pseudorandom image.
        let seed = 12345;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const data = makeImage(xB, yB, () => 10 * rand());
        const min = Math.min(...data);

        const filter = new HeatFilter2(xB, yB, dx, dy, data, null, border);
        filter.setTimeStep(dt);

        // Reference solver on a padded grid whose border is held at the
        // Dirichlet value. The stored image is shifted by the minimum.
        const w = xB + 2;
        let u = new Array<number>(w * (yB + 2)).fill(border);
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                u[(x + 1) + w * (y + 1)] = data[x + xB * y] - min;
            }
        }

        for (let step = 0; step < 8; ++step) {
            const next = u.slice();
            for (let y = 1; y <= yB; ++y) {
                for (let x = 1; x <= xB; ++x) {
                    const i = x + w * y;
                    next[i] = u[i] + dt * (
                        (u[i + 1] - 2 * u[i] + u[i - 1]) / (dx * dx) +
                        (u[i + w] - 2 * u[i] + u[i - w]) / (dy * dy));
                }
            }
            u = next;
            filter.update();
        }

        const expected: number[] = [];
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                expected.push(u[(x + 1) + w * (y + 1)]);
            }
        }
        closeTo(filter.values(), expected);
    });

    it('keeps a Neumann image border stale across iterations', () => {
        // Upstream never re-runs AssignNeumannImageBorder after construction
        // (onPreUpdate only refreshes a mask border), so the padded border
        // keeps the values duplicated from the original image.
        const data = [1, 2, 3, 4];
        const filter = new HeatFilter2(2, 2, 1, 1, data, null, NEUMANN);
        filter.setTimeStep(0.1);
        const before = filter.padded();
        filter.update();
        const after = filter.padded();
        // The four corners of the padded image are untouched by onUpdate.
        expect(after[0]).toBe(before[0]);
        expect(after[3]).toBe(before[3]);
        expect(after[12]).toBe(before[12]);
        expect(after[15]).toBe(before[15]);
    });

    it('averages the unmasked neighbors for a Neumann mask border', () => {
        // The center pixel is masked out and has eight unmasked neighbors,
        // whose average replaces its value.
        const data = [0, 1, 2, 3, 100, 5, 6, 7, 8];
        const mask = [1, 1, 1, 1, 0, 1, 1, 1, 1];
        const filter = new HeatFilter2(3, 3, 1, 1, data, mask, NEUMANN);
        expect(filter.maskValues()).toEqual(mask);
        // (0+1+2+3+5+6+7+8)/8 = 4
        expect(filter.getU(1, 1)).toBeCloseTo(4, 12);
        // The unmasked pixels are untouched.
        closeTo([filter.getU(0, 0), filter.getU(2, 2)], [0, 8]);
    });

    it('assigns the border value on a Dirichlet mask border', () => {
        const data = [0, 1, 2, 3, 100, 5, 6, 7, 8];
        const mask = [1, 1, 1, 1, 0, 1, 1, 1, 1];
        const filter = new HeatFilter2(3, 3, 1, 1, data, mask, -7);
        expect(filter.getU(1, 1)).toBe(-7);
        closeTo([filter.getU(0, 0), filter.getU(2, 2)], [0, 8]);
    });

    it('leaves a masked pixel with no unmasked neighbors alone', () => {
        // A 5-by-1 image where pixel 0 is masked out and is not an
        // 8-neighbor of any unmasked pixel.
        const data = [100, 200, 0, 10, 20];
        const mask = [0, 0, 1, 1, 1];
        const filter = new HeatFilter2(5, 1, 1, 1, data, mask, NEUMANN);
        // Pixel 0 has no unmasked neighbor, so it keeps its stored value.
        expect(filter.getU(0, 0)).toBe(100);
        // Pixel 1 neighbors the unmasked pixel 2, so it takes its value.
        expect(filter.getU(1, 0)).toBeCloseTo(0, 12);
    });

    it('zeroes the mask on the image border when a mask is present', () => {
        const data = [1, 2, 3, 4];
        const mask = [1, 1, 1, 1];
        const filter = new HeatFilter2(2, 2, 1, 1, data, mask, NEUMANN);
        // getMask is only defined on the unpadded image, which is all ones;
        // the padded border is zeroed so it is never updated.
        expect(filter.maskValues()).toEqual([1, 1, 1, 1]);
    });

    it('visits the unmasked pixels in row-major padded order', () => {
        const data = makeImage(3, 2, (x, y) => x + 3 * y);
        const noMask = new RecordingFilter2(3, 2, data, null, 0);
        noMask.update();
        expect(noMask.visited).toEqual([[1, 1], [2, 1], [3, 1], [1, 2], [2, 2], [3, 2]]);

        const mask = [1, 0, 1, 0, 1, 1];
        const masked = new RecordingFilter2(3, 2, data, mask, 0);
        masked.update();
        expect(masked.visited).toEqual([[1, 1], [3, 1], [2, 2], [3, 2]]);
    });

    it('swaps the ping-pong buffers on every update', () => {
        const data = makeImage(2, 2, (x, y) => x + 2 * y);
        const filter = new RecordingFilter2(2, 2, data, null, 0);
        expect([filter.src(), filter.dst()]).toEqual([0, 1]);
        filter.update();
        expect([filter.src(), filter.dst()]).toEqual([1, 0]);
        filter.update();
        expect([filter.src(), filter.dst()]).toEqual([0, 1]);
    });
});

describe('PdeFilter2 verification', () => {
    const scaleTypes = [
        PdeFilterScaleType.NONE,
        PdeFilterScaleType.UNIT,
        PdeFilterScaleType.SYMMETRIC,
        PdeFilterScaleType.PRESERVE_ZERO
    ];

    // An independent implementation of the documented PdeFilter scaling; the
    // NONE branch still subtracts the minimum, as upstream does.
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

    // The padded (xB+2)-by-(yB+2) buffer the constructor builds, indexed
    // p[px + (xB + 2) * py]. Dirichlet fills the ring with the border value;
    // Neumann duplicates the nearest interior sample (the diagonal neighbor
    // at each corner).
    function referencePadded(xB: number, yB: number, data: readonly number[],
        borderValue: number, scaleType: PdeFilterScaleType): number[] {
        const scaled = referenceScale(data, scaleType);
        const w = xB + 2;
        const p = new Array<number>(w * (yB + 2)).fill(0);
        const at = (px: number, py: number) => px + w * py;
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                p[at(x + 1, y + 1)] = scaled[x + xB * y];
            }
        }
        const dirichlet = (borderValue !== NEUMANN);
        const ghost = (px: number, py: number) => {
            if (dirichlet) {
                return borderValue;
            }
            const sx = Math.min(Math.max(px, 1), xB);
            const sy = Math.min(Math.max(py, 1), yB);
            return p[at(sx, sy)];
        };
        for (let px = 0; px < w; ++px) {
            p[at(px, 0)] = ghost(px, 0);
            p[at(px, yB + 1)] = ghost(px, yB + 1);
        }
        for (let py = 0; py < yB + 2; ++py) {
            p[at(0, py)] = ghost(0, py);
            p[at(xB + 1, py)] = ghost(xB + 1, py);
        }
        return p;
    }

    const config = fc.tuple(
        fc.integer({ min: 3, max: 6 }),
        fc.integer({ min: 3, max: 6 }),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(0.5, 1, 2),
        fc.constantFrom(...scaleTypes),
        fc.constantFrom(NEUMANN, 0, -3, 5),
        fc.array(fc.integer({ min: -6, max: 6 }), { minLength: 36, maxLength: 36 }))
        .map(([xB, yB, hx, hy, scaleType, borderValue, pool]) => ({
            xB, yB, hx, hy, scaleType, borderValue,
            data: pool.slice(0, xB * yB)
        }));

    it('builds the padded buffer the border conditions prescribe', () => {
        check(config, ({ xB, yB, hx, hy, scaleType, borderValue, data }) => {
            const filter = new HeatFilter2(xB, yB, hx, hy, data, null,
                borderValue, scaleType);
            const expected = referencePadded(xB, yB, data, borderValue, scaleType);
            const actual = filter.padded();
            expect(actual.length).toBe(expected.length);
            for (let i = 0; i < expected.length; ++i) {
                expectClose(actual[i], expected[i], 1e-12, 1e-12);
            }
        });
    });

    it('every derivative accessor is exact for a bivariate quadratic', () => {
        // The 3x3 central-difference stencils reproduce all first and second
        // derivatives of a quadratic exactly, so a wrong coefficient, a wrong
        // padding offset or a swapped x/y index is caught here. The mixed
        // derivative uxy is the interesting one: it must be
        // (f(+1,+1) - f(-1,+1) - f(+1,-1) + f(-1,-1)) / (4 dx dy).
        const quadratic = fc.tuple(
            fc.integer({ min: 4, max: 7 }),
            fc.integer({ min: 4, max: 7 }),
            fc.constantFrom(0.5, 1, 2),
            fc.constantFrom(0.5, 1, 2),
            fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 6, maxLength: 6 }));
        check(quadratic, ([xB, yB, hx, hy, k]) => {
            const [a, b, c, d, e, f] = k;
            const F = (X: number, Y: number) =>
                a + b * X + c * Y + d * X * X + e * X * Y + f * Y * Y;
            const data = makeImage(xB, yB, (x, y) => F(x * hx, y * hy));
            const filter = new HeatFilter2(xB, yB, hx, hy, data, null, NEUMANN,
                PdeFilterScaleType.NONE);
            const mag = k.reduce((s, v) => s + Math.abs(v), 1)
                * (1 + xB * hx) * (1 + yB * hy);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    const X = x * hx, Y = y * hy;
                    expectClose(filter.getU(x, y), F(X, Y) - Math.min(...data),
                        1e-9 * mag, 1e-9);
                    expectClose(filter.getUx(x, y), b + 2 * d * X + e * Y,
                        1e-9 * mag / hx, 1e-9);
                    expectClose(filter.getUy(x, y), c + e * X + 2 * f * Y,
                        1e-9 * mag / hy, 1e-9);
                    expectClose(filter.getUxx(x, y), 2 * d,
                        1e-9 * mag / (hx * hx), 1e-9);
                    expectClose(filter.getUyy(x, y), 2 * f,
                        1e-9 * mag / (hy * hy), 1e-9);
                    expectClose(filter.getUxy(x, y), e,
                        1e-9 * mag / (hx * hy), 1e-9);
                }
            }
        });
    });

    it('one update step agrees with an independent explicit heat solver', () => {
        check(fc.tuple(config, fc.constantFrom(0.05, 0.1, 0.2)),
            ([{ xB, yB, hx, hy, scaleType, borderValue, data }, dtFactor]) => {
                const filter = new HeatFilter2(xB, yB, hx, hy, data, null,
                    borderValue, scaleType);
                const dt = dtFactor / (1 / (hx * hx) + 1 / (hy * hy));
                filter.setTimeStep(dt);

                const p = referencePadded(xB, yB, data, borderValue, scaleType);
                const w = xB + 2;
                const at = (px: number, py: number) => px + w * py;
                const expected: number[] = [];
                for (let y = 1; y <= yB; ++y) {
                    for (let x = 1; x <= xB; ++x) {
                        const uxx = (p[at(x + 1, y)] - 2 * p[at(x, y)]
                            + p[at(x - 1, y)]) / (hx * hx);
                        const uyy = (p[at(x, y + 1)] - 2 * p[at(x, y)]
                            + p[at(x, y - 1)]) / (hy * hy);
                        expected.push(p[at(x, y)] + dt * (uxx + uyy));
                    }
                }

                filter.update();
                const actual = filter.values();
                for (let i = 0; i < expected.length; ++i) {
                    expectClose(actual[i], expected[i], 1e-10, 1e-10);
                }
            });
    });

    it('getMask reports the caller mask and onUpdate skips masked pixels', () => {
        const masked = fc.tuple(
            fc.integer({ min: 3, max: 6 }),
            fc.integer({ min: 3, max: 6 }),
            fc.array(fc.integer({ min: -5, max: 5 }), { minLength: 36, maxLength: 36 }),
            fc.array(fc.integer({ min: 0, max: 1 }), { minLength: 36, maxLength: 36 }));
        check(masked, ([xB, yB, pool, maskPool]) => {
            const data = pool.slice(0, xB * yB);
            const mask = maskPool.slice(0, xB * yB);
            const filter = new RecordingFilter2(xB, yB, data, mask, NEUMANN);
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    expect(filter.getMask(x, y)).toBe(mask[x + xB * y]);
                }
            }

            filter.update();
            // Visited pixels are exactly the unmasked ones, in row-major
            // padded order.
            const expected: [number, number][] = [];
            for (let y = 1; y <= yB; ++y) {
                for (let x = 1; x <= xB; ++x) {
                    if (mask[(x - 1) + xB * (y - 1)] !== 0) {
                        expected.push([x, y]);
                    }
                }
            }
            expect(filter.visited).toEqual(expected);
        });
    });
});
