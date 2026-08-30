import { describe, it, expect } from 'vitest';
import { PdeFilter3 } from '../src/PdeFilter3';
import { PdeFilterScaleType } from '../src/PdeFilter';

// Concrete subclass solving the linear heat equation
//   u_t = u_xx + u_yy + u_zz
// with an explicit Euler step, the canonical use of the PdeFilter3 plumbing.
class HeatFilter3 extends PdeFilter3 {
    constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound, yBound, zBound, xSpacing, ySpacing, zSpacing, data,
            mask, borderValue, scaleType);
    }

    protected onUpdateSingle(x: number, y: number, z: number): void {
        this.lookUp7(x, y, z);
        const uxx = this.mInvDxDx * (this.mUpzz - 2 * this.mUzzz + this.mUmzz);
        const uyy = this.mInvDyDy * (this.mUzpz - 2 * this.mUzzz + this.mUzmz);
        const uzz = this.mInvDzDz * (this.mUzzp - 2 * this.mUzzz + this.mUzzm);
        this.mBuffer[this.mDst].set(x, y, z,
            this.mUzzz + this.mTimeStep * (uxx + uyy + uzz));
    }
}

// Subclass that records which voxels were visited by onUpdate and exposes
// the 27-neighborhood lookup for testing.
class RecordingFilter3 extends PdeFilter3 {
    visited: [number, number, number][] = [];

    constructor(xBound: number, yBound: number, zBound: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number) {
        super(xBound, yBound, zBound, 1, 1, 1, data, mask, borderValue,
            PdeFilterScaleType.NONE);
    }

    protected onUpdateSingle(x: number, y: number, z: number): void {
        // Padded coordinates; record and copy through unchanged.
        this.visited.push([x, y, z]);
        this.lookUp7(x, y, z);
        this.mBuffer[this.mDst].set(x, y, z, this.mUzzz);
    }

    // Expose lookUp27 results for the padded voxel (x+1,y+1,z+1).
    neighborhood27(x: number, y: number, z: number): number[] {
        this.lookUp27(x + 1, y + 1, z + 1);
        return [
            this.mUmmm, this.mUzmm, this.mUpmm,
            this.mUmzm, this.mUzzm, this.mUpzm,
            this.mUmpm, this.mUzpm, this.mUppm,
            this.mUmmz, this.mUzmz, this.mUpmz,
            this.mUmzz, this.mUzzz, this.mUpzz,
            this.mUmpz, this.mUzpz, this.mUppz,
            this.mUmmp, this.mUzmp, this.mUpmp,
            this.mUmzp, this.mUzzp, this.mUpzp,
            this.mUmpp, this.mUzpp, this.mUppp
        ];
    }
}

const NEUMANN = Number.MAX_VALUE;

// Build data[i] = f(x,y,z) for i = x + xBound*(y + yBound*z).
function makeImage(xB: number, yB: number, zB: number,
    f: (x: number, y: number, z: number) => number): number[] {
    const data = new Array<number>(xB * yB * zB);
    let i = 0;
    for (let z = 0; z < zB; ++z) {
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x, ++i) {
                data[i] = f(x, y, z);
            }
        }
    }
    return data;
}

describe('PdeFilter3', () => {
    it('provides member access', () => {
        const data = makeImage(3, 4, 5, () => 1);
        // Constant image: min === max path in the base class.
        const filter = new HeatFilter3(3, 4, 5, 0.5, 0.25, 2, data, null,
            NEUMANN, PdeFilterScaleType.NONE);
        expect(filter.getXBound()).toBe(3);
        expect(filter.getYBound()).toBe(4);
        expect(filter.getZBound()).toBe(5);
        expect(filter.getXSpacing()).toBe(0.5);
        expect(filter.getYSpacing()).toBe(0.25);
        expect(filter.getZSpacing()).toBe(2);
        expect(filter.getQuantity()).toBe(60);
        expect(filter.getBorderValue()).toBe(NEUMANN);
    });

    it('getU returns data shifted by the minimum with NONE scaling', () => {
        // The base class stores offset + (d - min) * scale; for NONE the
        // offset is 0 and the scale is 1, so the stored image is d - min
        // (upstream semantics).
        const f = (x: number, y: number, z: number) => 1 + x + 10 * y + 100 * z;
        const data = makeImage(3, 3, 3, f);
        const filter = new HeatFilter3(3, 3, 3, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE);
        const min = 1;  // f(0,0,0)
        for (let z = 0; z < 3; ++z) {
            for (let y = 0; y < 3; ++y) {
                for (let x = 0; x < 3; ++x) {
                    expect(filter.getU(x, y, z)).toBe(f(x, y, z) - min);
                }
            }
        }
    });

    it('getU applies UNIT scaling', () => {
        const data = [2, 6, 4, 10, 2, 6, 4, 10];  // min 2, max 10
        const filter = new HeatFilter3(2, 2, 2, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.UNIT);
        // d' = (d - min) / (max - min)
        expect(filter.getU(0, 0, 0)).toBeCloseTo(0, 14);
        expect(filter.getU(1, 0, 0)).toBeCloseTo(0.5, 14);
        expect(filter.getU(0, 1, 0)).toBeCloseTo(0.25, 14);
        expect(filter.getU(1, 1, 0)).toBeCloseTo(1, 14);
    });

    it('estimates first and second derivatives of a quadratic field', () => {
        // u(x,y,z) = x^2 + 2y^2 + 3z^2 + xy + 2xz + 4yz + x + y + z + 7 in
        // index coordinates. Central differences are exact for quadratics.
        const u = (x: number, y: number, z: number) =>
            x * x + 2 * y * y + 3 * z * z + x * y + 2 * x * z + 4 * y * z
            + x + y + z + 7;
        const dx = 0.5, dy = 1, dz = 2;
        const data = makeImage(5, 5, 5, u);
        const filter = new HeatFilter3(5, 5, 5, dx, dy, dz, data, null,
            NEUMANN, PdeFilterScaleType.NONE);

        for (const [x, y, z] of [[1, 1, 1], [2, 3, 1], [3, 2, 2]] as const) {
            // First derivatives: (index-space central difference) / spacing.
            expect(filter.getUx(x, y, z)).toBeCloseTo((2 * x + y + 2 * z + 1) / dx, 12);
            expect(filter.getUy(x, y, z)).toBeCloseTo((4 * y + x + 4 * z + 1) / dy, 12);
            expect(filter.getUz(x, y, z)).toBeCloseTo((6 * z + 2 * x + 4 * y + 1) / dz, 12);
            // Second derivatives.
            expect(filter.getUxx(x, y, z)).toBeCloseTo(2 / (dx * dx), 12);
            expect(filter.getUyy(x, y, z)).toBeCloseTo(4 / (dy * dy), 12);
            expect(filter.getUzz(x, y, z)).toBeCloseTo(6 / (dz * dz), 12);
            // Mixed derivatives.
            expect(filter.getUxy(x, y, z)).toBeCloseTo(1 / (dx * dy), 12);
            expect(filter.getUxz(x, y, z)).toBeCloseTo(2 / (dx * dz), 12);
            expect(filter.getUyz(x, y, z)).toBeCloseTo(4 / (dy * dz), 12);
        }
    });

    it('Neumann boundary duplicates the nearest interior value', () => {
        // u = 3x + 5y + 7z. With Neumann conditions the padded border
        // duplicates the closest image sample, so the one-sided part of the
        // central difference collapses.
        const u = (x: number, y: number, z: number) => 3 * x + 5 * y + 7 * z;
        const data = makeImage(4, 4, 4, u);
        const filter = new HeatFilter3(4, 4, 4, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE);

        // At x = 0: (u(1) - u(0)) / 2 = 3/2. At x = 3: same by symmetry.
        expect(filter.getUx(0, 1, 1)).toBeCloseTo(1.5, 14);
        expect(filter.getUx(3, 1, 1)).toBeCloseTo(1.5, 14);
        expect(filter.getUy(1, 0, 1)).toBeCloseTo(2.5, 14);
        expect(filter.getUy(1, 3, 1)).toBeCloseTo(2.5, 14);
        expect(filter.getUz(1, 1, 0)).toBeCloseTo(3.5, 14);
        expect(filter.getUz(1, 1, 3)).toBeCloseTo(3.5, 14);

        // Second difference at the border: u(1) - 2u(0) + duplicate(u(0)).
        expect(filter.getUxx(0, 1, 1)).toBeCloseTo(u(1, 1, 1) - u(0, 1, 1), 14);

        // Corner duplicates: getUxy at a corner voxel uses the border
        // corners; on a linear field the interior estimate is 0 and the
        // Neumann padding keeps it finite (no MAX_VALUE contamination).
        expect(Number.isFinite(filter.getUxy(0, 0, 0))).toBe(true);
        expect(Number.isFinite(filter.getUyz(3, 3, 3))).toBe(true);
    });

    it('Dirichlet boundary uses the constant border value', () => {
        // Note: a constant image stores as d - min = 0 everywhere, so the
        // interior values are 0 and the border keeps the raw border value
        // (the border assignment is not rescaled, matching upstream).
        const border = 10;
        const data = makeImage(2, 2, 2, () => 7);
        const filter = new HeatFilter3(2, 2, 2, 1, 1, 1, data, null,
            border, PdeFilterScaleType.NONE);

        expect(filter.getU(0, 0, 0)).toBe(0);
        // getUx(0,0,0) = (u(1,0,0) - border) / 2 = (0 - 10)/2.
        expect(filter.getUx(0, 0, 0)).toBeCloseTo(-5, 14);
        // getUx(1,0,0) = (border - u(0,0,0)) / 2 = (10 - 0)/2.
        expect(filter.getUx(1, 0, 0)).toBeCloseTo(5, 14);
        // getUzz(0,0,0) = u(0,0,1) - 2 u(0,0,0) + border.
        expect(filter.getUzz(0, 0, 0)).toBeCloseTo(10, 14);
    });

    it('performs a hand-computable heat-equation update step', () => {
        // 3x3x3 image, delta impulse at the center, unit spacing, Neumann.
        const data = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1 ? 1 : 0));
        const filter = new HeatFilter3(3, 3, 3, 1, 1, 1, data, null,
            NEUMANN, PdeFilterScaleType.NONE);
        const dt = 0.05;
        filter.setTimeStep(dt);
        expect(filter.getTimeStep()).toBe(dt);
        filter.update();

        // Explicit Euler: u += dt * laplacian(u).
        // Center: 1 + dt * (6*0 - 6*1) = 1 - 6 dt.
        expect(filter.getU(1, 1, 1)).toBeCloseTo(1 - 6 * dt, 14);
        // Face neighbors: 0 + dt * 1 (one neighbor is the center).
        expect(filter.getU(0, 1, 1)).toBeCloseTo(dt, 14);
        expect(filter.getU(2, 1, 1)).toBeCloseTo(dt, 14);
        expect(filter.getU(1, 0, 1)).toBeCloseTo(dt, 14);
        expect(filter.getU(1, 2, 1)).toBeCloseTo(dt, 14);
        expect(filter.getU(1, 1, 0)).toBeCloseTo(dt, 14);
        expect(filter.getU(1, 1, 2)).toBeCloseTo(dt, 14);
        // Edge and corner neighbors are untouched by the 7-point stencil.
        expect(filter.getU(0, 0, 1)).toBe(0);
        expect(filter.getU(0, 0, 0)).toBe(0);
        expect(filter.getU(2, 2, 2)).toBe(0);

        // A second update diffuses further. For the face voxel (0,1,1) the
        // stencil is: border duplicate 0 (assigned at construction and not
        // recomputed per iteration, matching upstream), center 1-6dt, four
        // zero neighbors, self dt:
        //   dt + dt * ((1 - 6dt) - 6dt) = 2dt - 12dt^2.
        // For the center: 6 face neighbors valued dt:
        //   (1 - 6dt) + dt * (6dt - 6(1 - 6dt)) = 1 - 12dt + 42dt^2.
        // The corner voxel's face neighbors were all 0 after step one, so it
        // remains 0.
        filter.update();
        expect(filter.getU(0, 1, 1)).toBeCloseTo(2 * dt - 12 * dt * dt, 14);
        expect(filter.getU(1, 1, 1)).toBeCloseTo(1 - 12 * dt + 42 * dt * dt, 14);
        expect(filter.getU(0, 0, 0)).toBe(0);
    });

    it('onUpdate visits only unmasked voxels and buffers swap correctly', () => {
        const data = makeImage(3, 3, 3, (x, y, z) => x + y + z);
        const mask = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1) || (x === 0 && y === 1 && z === 2) ? 1 : 0);
        const filter = new RecordingFilter3(3, 3, 3, data, mask, NEUMANN);
        filter.update();
        // Padded coordinates of the two masked-in voxels.
        expect(filter.visited).toEqual([[2, 2, 2], [1, 2, 3]]);
        expect(filter.getMask(1, 1, 1)).toBe(1);
        expect(filter.getMask(0, 1, 2)).toBe(1);
        expect(filter.getMask(0, 0, 0)).toBe(0);

        // Without a mask every voxel is visited.
        const nomask = new RecordingFilter3(3, 3, 3, data, null, NEUMANN);
        nomask.update();
        expect(nomask.visited.length).toBe(27);
    });

    it('Neumann mask border averages the masked neighbors', () => {
        // Only the center voxel is masked-in, with value 5. Every other
        // voxel of the 3x3x3 image is a 26-neighbor of the center, so it is
        // assigned the average of its masked neighbors, which is 5.
        const data = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1 ? 5 : 0));
        const mask = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1 ? 1 : 0));
        const filter = new HeatFilter3(3, 3, 3, 1, 1, 1, data, mask,
            NEUMANN, PdeFilterScaleType.NONE);

        for (let z = 0; z < 3; ++z) {
            for (let y = 0; y < 3; ++y) {
                for (let x = 0; x < 3; ++x) {
                    expect(filter.getU(x, y, z)).toBe(5);
                }
            }
        }

        // The center's 7-point neighborhood is constant, so an update leaves
        // it unchanged (and unmasked voxels are never processed).
        filter.setTimeStep(0.05);
        filter.update();
        expect(filter.getU(1, 1, 1)).toBe(5);
        expect(filter.getU(0, 0, 0)).toBe(5);
    });

    it('Dirichlet mask border assigns the border value to mask neighbors', () => {
        const border = 2;
        // The minimum is 0, so the stored image equals the data.
        const data = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1 ? 5 : 0));
        const mask = makeImage(3, 3, 3,
            (x, y, z) => (x === 1 && y === 1 && z === 1 ? 1 : 0));
        const filter = new HeatFilter3(3, 3, 3, 1, 1, 1, data, mask,
            border, PdeFilterScaleType.NONE);

        // All 26 neighbors of the center are masked-out with a masked
        // neighbor, so both buffers hold the border value there.
        expect(filter.getU(0, 0, 0)).toBe(border);
        expect(filter.getU(2, 1, 1)).toBe(border);
        expect(filter.getU(1, 1, 1)).toBe(5);

        // Update: center gets 5 + dt * (6*2 - 6*5) = 5 - 18 dt; the
        // masked-out voxels keep the border value in the swapped buffer.
        const dt = 0.05;
        filter.setTimeStep(dt);
        filter.update();
        expect(filter.getU(1, 1, 1)).toBeCloseTo(5 - 18 * dt, 14);
        expect(filter.getU(0, 0, 0)).toBe(border);
        expect(filter.getU(2, 2, 1)).toBe(border);
    });

    it('lookUp27 loads the full 3x3x3 neighborhood in xyz order', () => {
        const f = (x: number, y: number, z: number) => x + 10 * y + 100 * z;
        const data = makeImage(3, 3, 3, f);
        const filter = new RecordingFilter3(3, 3, 3, data, null, NEUMANN);
        const hood = filter.neighborhood27(1, 1, 1);
        const expected: number[] = [];
        for (let dz = -1; dz <= 1; ++dz) {
            for (let dy = -1; dy <= 1; ++dy) {
                for (let dx = -1; dx <= 1; ++dx) {
                    expected.push(f(1 + dx, 1 + dy, 1 + dz));
                }
            }
        }
        expect(hood).toEqual(expected);
    });
});
