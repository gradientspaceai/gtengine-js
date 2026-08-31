// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PdeFilter2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The abstract base class for PDE-based filtering of 2-dimensional images.
// The internal images for "data" and "mask" are copies of the constructor
// inputs, padded with a 1-pixel thick border to support filtering on the image
// boundary; they are (xBound+2)-by-(yBound+2). The lookups into the padded
// arrays are handled internally, so the public accessors take 0 <= x < xBound
// and 0 <= y < yBound.
//
// Port notes:
//   - The template parameter Real is number (IEEE double). The ping-pong
//     buffers and the mask are Array2 instances, as upstream; following the
//     Array2 port, the C++ 'a[y][x]' becomes 'a.get(x, y)'. Array2 leaves its
//     storage uninitialized, so the buffers and mask are explicitly zero
//     filled to match the value-initialization of the upstream std::vector
//     storage.
//   - 'Real const* data' becomes ArrayLike<number> and the nullable
//     'int32_t const* mask' becomes ArrayLike<number> | null.
//   - Number.MAX_VALUE is the Neumann marker for the border value, matching
//     the upstream std::numeric_limits<Real>::max() test in PdeFilter.

import { Array2 } from './Array2';
import { PdeFilter, PdeFilterScaleType } from './PdeFilter';

export abstract class PdeFilter2 extends PdeFilter {
    // Image parameters.
    protected mXBound: number;
    protected mYBound: number;
    protected mXSpacing: number;        // dx
    protected mYSpacing: number;        // dy
    protected mInvDx: number;           // 1/dx
    protected mInvDy: number;           // 1/dy
    protected mHalfInvDx: number;       // 1/(2*dx)
    protected mHalfInvDy: number;       // 1/(2*dy)
    protected mInvDxDx: number;         // 1/(dx*dx)
    protected mFourthInvDxDy: number;   // 1/(4*dx*dy)
    protected mInvDyDy: number;         // 1/(dy*dy)

    // Temporary storage for the 3x3 neighborhood. In the notation mUxy, the x
    // and y indices are in {m,z,p}, referring to subtract 1 (m), no change
    // (z), or add 1 (p) to the appropriate index.
    protected mUmm: number; protected mUzm: number; protected mUpm: number;
    protected mUmz: number; protected mUzz: number; protected mUpz: number;
    protected mUmp: number; protected mUzp: number; protected mUpp: number;

    // Successive iterations toggle between two buffers.
    protected mBuffer: [Array2<number>, Array2<number>];
    protected mSrc: number;
    protected mDst: number;
    protected mMask: Array2<number>;
    protected mHasMask: boolean;

    protected constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType) {
        super(xBound * yBound, data, borderValue, scaleType);

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mXSpacing = xSpacing;
        this.mYSpacing = ySpacing;
        this.mInvDx = 1 / xSpacing;
        this.mInvDy = 1 / ySpacing;
        this.mHalfInvDx = 0.5 * this.mInvDx;
        this.mHalfInvDy = 0.5 * this.mInvDy;
        this.mInvDxDx = this.mInvDx * this.mInvDx;
        this.mFourthInvDxDy = this.mHalfInvDx * this.mHalfInvDy;
        this.mInvDyDy = this.mInvDy * this.mInvDy;
        this.mUmm = 0; this.mUzm = 0; this.mUpm = 0;
        this.mUmz = 0; this.mUzz = 0; this.mUpz = 0;
        this.mUmp = 0; this.mUzp = 0; this.mUpp = 0;
        this.mSrc = 0;
        this.mDst = 1;
        this.mMask = new Array2<number>(xBound + 2, yBound + 2);
        this.mMask.fill(0);
        this.mHasMask = (mask !== null);

        // The mBuffer[] are ping-pong buffers for filtering.
        this.mBuffer = [
            new Array2<number>(xBound + 2, yBound + 2),
            new Array2<number>(xBound + 2, yBound + 2)
        ];
        this.mBuffer[0].fill(0);
        this.mBuffer[1].fill(0);

        for (let y = 0, yp = 1, i = 0; y < this.mYBound; ++y, ++yp) {
            for (let x = 0, xp = 1; x < this.mXBound; ++x, ++xp, ++i) {
                this.mBuffer[this.mSrc].set(xp, yp, this.mOffset + (data[i] - this.mMin) * this.mScale);
                this.mBuffer[this.mDst].set(xp, yp, 0);
                this.mMask.set(xp, yp, mask !== null ? mask[i] : 1);
            }
        }

        // Assign values to the 1-pixel image border.
        if (this.mBorderValue !== Number.MAX_VALUE) {
            this.assignDirichletImageBorder();
        } else {
            this.assignNeumannImageBorder();
        }

        // To handle masks that do not cover the entire image, assign values to
        // those pixels that are 8-neighbors of the mask pixels.
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

    getYBound(): number {
        return this.mYBound;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    getYSpacing(): number {
        return this.mYSpacing;
    }

    // Pixel access and derivative estimation. The lookups into the padded data
    // are handled correctly. The estimation involves only the 3-by-3
    // neighborhood of (x,y), where 0 <= x < xBound and 0 <= y < yBound. If
    // larger neighborhoods are desired at a later date, the padding and
    // associated code must be adjusted accordingly.
    getU(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return F.get(x + 1, y + 1);
    }

    getUx(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mHalfInvDx * (F.get(x + 2, y + 1) - F.get(x, y + 1));
    }

    getUy(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mHalfInvDy * (F.get(x + 1, y + 2) - F.get(x + 1, y));
    }

    getUxx(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mInvDxDx * (F.get(x + 2, y + 1) - 2 * F.get(x + 1, y + 1) + F.get(x, y + 1));
    }

    getUxy(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mFourthInvDxDy * (F.get(x, y) - F.get(x + 2, y)
            + F.get(x + 2, y + 2) - F.get(x, y + 2));
    }

    getUyy(x: number, y: number): number {
        const F = this.mBuffer[this.mSrc];
        return this.mInvDyDy * (F.get(x + 1, y + 2) - 2 * F.get(x + 1, y + 1) + F.get(x + 1, y));
    }

    getMask(x: number, y: number): number {
        return this.mMask.get(x + 1, y + 1);
    }

    // Assign values to the 1-pixel image border.
    protected assignDirichletImageBorder(): void {
        const xBp1 = this.mXBound + 1;
        const yBp1 = this.mYBound + 1;
        const src = this.mBuffer[this.mSrc];
        const dst = this.mBuffer[this.mDst];
        const value = this.mBorderValue;

        // vertex (0,0)
        src.set(0, 0, value);
        dst.set(0, 0, value);
        if (this.mHasMask) {
            this.mMask.set(0, 0, 0);
        }

        // vertex (xmax,0)
        src.set(xBp1, 0, value);
        dst.set(xBp1, 0, value);
        if (this.mHasMask) {
            this.mMask.set(xBp1, 0, 0);
        }

        // vertex (0,ymax)
        src.set(0, yBp1, value);
        dst.set(0, yBp1, value);
        if (this.mHasMask) {
            this.mMask.set(0, yBp1, 0);
        }

        // vertex (xmax,ymax)
        src.set(xBp1, yBp1, value);
        dst.set(xBp1, yBp1, value);
        if (this.mHasMask) {
            this.mMask.set(xBp1, yBp1, 0);
        }

        // edges (x,0) and (x,ymax)
        for (let x = 1; x <= this.mXBound; ++x) {
            src.set(x, 0, value);
            dst.set(x, 0, value);
            if (this.mHasMask) {
                this.mMask.set(x, 0, 0);
            }

            src.set(x, yBp1, value);
            dst.set(x, yBp1, value);
            if (this.mHasMask) {
                this.mMask.set(x, yBp1, 0);
            }
        }

        // edges (0,y) and (xmax,y)
        for (let y = 1; y <= this.mYBound; ++y) {
            src.set(0, y, value);
            dst.set(0, y, value);
            if (this.mHasMask) {
                this.mMask.set(0, y, 0);
            }

            src.set(xBp1, y, value);
            dst.set(xBp1, y, value);
            if (this.mHasMask) {
                this.mMask.set(xBp1, y, 0);
            }
        }
    }

    protected assignNeumannImageBorder(): void {
        const xBp1 = this.mXBound + 1;
        const yBp1 = this.mYBound + 1;
        const src = this.mBuffer[this.mSrc];
        const dst = this.mBuffer[this.mDst];
        let duplicate: number;

        // vertex (0,0)
        duplicate = src.get(1, 1);
        src.set(0, 0, duplicate);
        dst.set(0, 0, duplicate);
        if (this.mHasMask) {
            this.mMask.set(0, 0, 0);
        }

        // vertex (xmax,0)
        duplicate = src.get(this.mXBound, 1);
        src.set(xBp1, 0, duplicate);
        dst.set(xBp1, 0, duplicate);
        if (this.mHasMask) {
            this.mMask.set(xBp1, 0, 0);
        }

        // vertex (0,ymax)
        duplicate = src.get(1, this.mYBound);
        src.set(0, yBp1, duplicate);
        dst.set(0, yBp1, duplicate);
        if (this.mHasMask) {
            this.mMask.set(0, yBp1, 0);
        }

        // vertex (xmax,ymax)
        duplicate = src.get(this.mXBound, this.mYBound);
        src.set(xBp1, yBp1, duplicate);
        dst.set(xBp1, yBp1, duplicate);
        if (this.mHasMask) {
            this.mMask.set(xBp1, yBp1, 0);
        }

        // edges (x,0) and (x,ymax)
        for (let x = 1; x <= this.mXBound; ++x) {
            duplicate = src.get(x, 1);
            src.set(x, 0, duplicate);
            dst.set(x, 0, duplicate);
            if (this.mHasMask) {
                this.mMask.set(x, 0, 0);
            }

            duplicate = src.get(x, this.mYBound);
            src.set(x, yBp1, duplicate);
            dst.set(x, yBp1, duplicate);
            if (this.mHasMask) {
                this.mMask.set(x, yBp1, 0);
            }
        }

        // edges (0,y) and (xmax,y)
        for (let y = 1; y <= this.mYBound; ++y) {
            duplicate = src.get(1, y);
            src.set(0, y, duplicate);
            dst.set(0, y, duplicate);
            if (this.mHasMask) {
                this.mMask.set(0, y, 0);
            }

            duplicate = src.get(this.mXBound, y);
            src.set(xBp1, y, duplicate);
            dst.set(xBp1, y, duplicate);
            if (this.mHasMask) {
                this.mMask.set(xBp1, y, 0);
            }
        }
    }

    // Assign values to the 1-pixel mask border.
    protected assignDirichletMaskBorder(): void {
        for (let y = 1; y <= this.mYBound; ++y) {
            for (let x = 1; x <= this.mXBound; ++x) {
                if (this.mMask.get(x, y) !== 0) {
                    continue;
                }

                let found = false;
                for (let i1 = 0, j1 = y - 1; i1 < 3 && !found; ++i1, ++j1) {
                    for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                        if (this.mMask.get(j0, j1) !== 0) {
                            this.mBuffer[this.mSrc].set(x, y, this.mBorderValue);
                            this.mBuffer[this.mDst].set(x, y, this.mBorderValue);
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    protected assignNeumannMaskBorder(): void {
        // Recompute the values just outside the masked region. This guarantees
        // that derivative estimations use the current values around the
        // boundary.
        for (let y = 1; y <= this.mYBound; ++y) {
            for (let x = 1; x <= this.mXBound; ++x) {
                if (this.mMask.get(x, y) !== 0) {
                    continue;
                }

                let count = 0;
                let average = 0;
                for (let i1 = 0, j1 = y - 1; i1 < 3; ++i1, ++j1) {
                    for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                        if (this.mMask.get(j0, j1) !== 0) {
                            average += this.mBuffer[this.mSrc].get(j0, j1);
                            ++count;
                        }
                    }
                }

                if (count > 0) {
                    average /= count;
                    this.mBuffer[this.mSrc].set(x, y, average);
                    this.mBuffer[this.mDst].set(x, y, average);
                }
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

    // Iterate over all the pixels and call onUpdateSingle(x,y) for each pixel
    // that is not masked out.
    protected override onUpdate(): void {
        for (let y = 1; y <= this.mYBound; ++y) {
            for (let x = 1; x <= this.mXBound; ++x) {
                if (!this.mHasMask || this.mMask.get(x, y) !== 0) {
                    this.onUpdateSingle(x, y);
                }
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

    // The per-pixel processing depends on the PDE algorithm. The (x,y) must be
    // in padded coordinates: 1 <= x <= xBound and 1 <= y <= yBound.
    protected abstract onUpdateSingle(x: number, y: number): void;

    // Copy source data to temporary storage.
    protected lookUp5(x: number, y: number): void {
        const F = this.mBuffer[this.mSrc];
        const xm = x - 1, xp = x + 1;
        const ym = y - 1, yp = y + 1;
        this.mUzm = F.get(x, ym);
        this.mUmz = F.get(xm, y);
        this.mUzz = F.get(x, y);
        this.mUpz = F.get(xp, y);
        this.mUzp = F.get(x, yp);
    }

    protected lookUp9(x: number, y: number): void {
        const F = this.mBuffer[this.mSrc];
        const xm = x - 1, xp = x + 1;
        const ym = y - 1, yp = y + 1;
        this.mUmm = F.get(xm, ym);
        this.mUzm = F.get(x, ym);
        this.mUpm = F.get(xp, ym);
        this.mUmz = F.get(xm, y);
        this.mUzz = F.get(x, y);
        this.mUpz = F.get(xp, y);
        this.mUmp = F.get(xm, yp);
        this.mUzp = F.get(x, yp);
        this.mUpp = F.get(xp, yp);
    }
}
