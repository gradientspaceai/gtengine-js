// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PdeFilter1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The abstract base class for PDE-based filtering of 1-dimensional images.
// The internal images for "data" and "mask" are copies of the constructor
// inputs, padded with a 1-pixel border to support filtering on the image
// boundary; they have xBound+2 elements. The lookups into the padded arrays
// are handled internally, so the public accessors take 0 <= x < xBound.
//
// Port notes:
//   - The template parameter Real is number (IEEE double). The ping-pong
//     buffers are Float64Array and the mask is an Int32Array, matching the
//     value-initialization (zero-fill) of the upstream std::vector members.
//   - 'Real const* data' becomes ArrayLike<number> and the nullable
//     'int32_t const* mask' becomes ArrayLike<number> | null.
//   - Upstream declares the per-pixel hook as an overload 'OnUpdate(int32_t)'
//     of the no-argument 'OnUpdate()'. TypeScript cannot give two overloads
//     distinct bodies (and an abstract overload is not expressible), so the
//     per-pixel hook is named onUpdateSingle(x), matching what upstream itself
//     calls the analogous hook in PdeFilter2/PdeFilter3.
//   - Number.MAX_VALUE is the Neumann marker for the border value, matching
//     the upstream std::numeric_limits<Real>::max() test in PdeFilter.

import { PdeFilter, PdeFilterScaleType } from './PdeFilter.js';

export abstract class PdeFilter1 extends PdeFilter {
    // Image parameters.
    protected mXBound: number;
    protected mXSpacing: number;    // dx
    protected mInvDx: number;       // 1/dx
    protected mHalfInvDx: number;   // 1/(2*dx)
    protected mInvDxDx: number;     // 1/(dx*dx)

    // Temporary storage for the 3-tuple neighborhood. In the notation mUx,
    // the x index is in {m,z,p}, referring to subtract 1 (m), no change (z),
    // or add 1 (p) to the index.
    protected mUm: number;
    protected mUz: number;
    protected mUp: number;

    // Successive iterations toggle between two buffers.
    protected mBuffer: [Float64Array, Float64Array];
    protected mSrc: number;
    protected mDst: number;
    protected mMask: Int32Array;
    protected mHasMask: boolean;

    protected constructor(xBound: number, xSpacing: number, data: ArrayLike<number>,
        mask: ArrayLike<number> | null, borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound, data, borderValue, scaleType);

        this.mXBound = xBound;
        this.mXSpacing = xSpacing;
        this.mInvDx = 1 / xSpacing;
        this.mHalfInvDx = 0.5 * this.mInvDx;
        this.mInvDxDx = this.mInvDx * this.mInvDx;
        this.mUm = 0;
        this.mUz = 0;
        this.mUp = 0;
        this.mSrc = 0;
        this.mDst = 1;
        this.mMask = new Int32Array(xBound + 2);
        this.mHasMask = (mask !== null);

        // The mBuffer[] are ping-pong buffers for filtering.
        this.mBuffer = [new Float64Array(xBound + 2), new Float64Array(xBound + 2)];

        for (let x = 0, xp = 1, i = 0; x < this.mXBound; ++x, ++xp, ++i) {
            this.mBuffer[this.mSrc][xp] = this.mOffset + (data[i] - this.mMin) * this.mScale;
            this.mBuffer[this.mDst][xp] = 0;
            this.mMask[xp] = (mask !== null ? mask[i] : 1);
        }

        // Assign values to the 1-pixel image border.
        if (this.mBorderValue !== Number.MAX_VALUE) {
            this.assignDirichletImageBorder();
        } else {
            this.assignNeumannImageBorder();
        }

        // To handle masks that do not cover the entire image, assign values to
        // those pixels that are neighbors of the mask pixels.
        if (this.mHasMask) {
            if (this.mBorderValue !== Number.MAX_VALUE) {
                this.assignDirichletMaskBorder();
            } else {
                this.assignNeumannMaskBorder();
            }
        }
    }

    // Member access.
    getXBound(): number {
        return this.mXBound;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    // Pixel access and derivative estimation. The lookups into the padded data
    // are handled correctly. The estimation involves only the 3-tuple
    // neighborhood of (x), where 0 <= x < xBound. If larger neighborhoods are
    // desired at a later date, the padding and associated code must be
    // adjusted accordingly.
    getU(x: number): number {
        const F = this.mBuffer[this.mSrc];
        return F[x + 1];
    }

    getUx(x: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mHalfInvDx * (F[x + 2] - F[x]);
    }

    getUxx(x: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mInvDxDx * (F[x + 2] - 2 * F[x + 1] + F[x]);
    }

    getMask(x: number): number {
        return this.mMask[x + 1];
    }

    // Assign values to the 1-pixel image border.
    protected assignDirichletImageBorder(): void {
        const xBp1 = this.mXBound + 1;

        // vertex (0)
        this.mBuffer[this.mSrc][0] = this.mBorderValue;
        this.mBuffer[this.mDst][0] = this.mBorderValue;
        if (this.mHasMask) {
            this.mMask[0] = 0;
        }

        // vertex (xmax)
        this.mBuffer[this.mSrc][xBp1] = this.mBorderValue;
        this.mBuffer[this.mDst][xBp1] = this.mBorderValue;
        if (this.mHasMask) {
            this.mMask[xBp1] = 0;
        }
    }

    protected assignNeumannImageBorder(): void {
        const xBp1 = this.mXBound + 1;
        let duplicate: number;

        // vertex (0)
        duplicate = this.mBuffer[this.mSrc][1];
        this.mBuffer[this.mSrc][0] = duplicate;
        this.mBuffer[this.mDst][0] = duplicate;
        if (this.mHasMask) {
            this.mMask[0] = 0;
        }

        // vertex (xmax)
        duplicate = this.mBuffer[this.mSrc][this.mXBound];
        this.mBuffer[this.mSrc][xBp1] = duplicate;
        this.mBuffer[this.mDst][xBp1] = duplicate;
        if (this.mHasMask) {
            this.mMask[xBp1] = 0;
        }
    }

    // Assign values to the 1-pixel mask border.
    protected assignDirichletMaskBorder(): void {
        for (let x = 1; x <= this.mXBound; ++x) {
            if (this.mMask[x] !== 0) {
                continue;
            }

            for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                if (this.mMask[j0] !== 0) {
                    this.mBuffer[this.mSrc][x] = this.mBorderValue;
                    this.mBuffer[this.mDst][x] = this.mBorderValue;
                    break;
                }
            }
        }
    }

    protected assignNeumannMaskBorder(): void {
        // Recompute the values just outside the masked region. This guarantees
        // that derivative estimations use the current values around the
        // boundary.
        for (let x = 1; x <= this.mXBound; ++x) {
            if (this.mMask[x] !== 0) {
                continue;
            }

            let count = 0;
            let average = 0;
            for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                if (this.mMask[j0] !== 0) {
                    average += this.mBuffer[this.mSrc][j0];
                    ++count;
                }
            }

            if (count > 0) {
                average /= count;
                this.mBuffer[this.mSrc][x] = average;
                this.mBuffer[this.mDst][x] = average;
            }
        }
    }

    // This function recomputes the boundary values when Neumann conditions are
    // used. If a derived class overrides this, it must call the base-class
    // onPreUpdate first.
    protected override onPreUpdate(): void {
        if (this.mHasMask && this.mBorderValue === Number.MAX_VALUE) {
            // Neumann boundary conditions are in use, so recompute the mask
            // border.
            this.assignNeumannMaskBorder();
        }
        // else: No mask has been specified or Dirichlet boundary conditions
        // are in use. Nothing to do.
    }

    // Iterate over all the pixels and call onUpdateSingle(x) for each pixel
    // that is not masked out.
    protected override onUpdate(): void {
        for (let x = 1; x <= this.mXBound; ++x) {
            if (!this.mHasMask || this.mMask[x] !== 0) {
                this.onUpdateSingle(x);
            }
        }
    }

    // If a derived class overrides this, it must call the base-class
    // onPostUpdate last. The base-class function swaps the buffers for the
    // next pass.
    protected override onPostUpdate(): void {
        const save = this.mSrc;
        this.mSrc = this.mDst;
        this.mDst = save;
    }

    // The per-pixel processing depends on the PDE algorithm. The x must be in
    // padded coordinates: 1 <= x <= xBound.
    protected abstract onUpdateSingle(x: number): void;

    // Copy source data to temporary storage.
    protected lookUp3(x: number): void {
        const F = this.mBuffer[this.mSrc];
        this.mUm = F[x - 1];
        this.mUz = F[x];
        this.mUp = F[x + 1];
    }
}
