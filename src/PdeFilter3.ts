// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PdeFilter3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Abstract base class for 3D PDE-based image filters. The image is stored in
// a linear array with lexicographic order, voxel (x,y,z) at index
// i = x + xBound * (y + yBound * z). Internally the image is copied into
// ping-pong buffers padded with a 1-voxel thick border to support filtering
// on the image boundary.
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// per-voxel border assignments in AssignDirichletImageBorder and
// AssignNeumannImageBorder are factored through the private helper
// assignBorderValue (a behavior-preserving refactor: each site sets the
// source buffer, the destination buffer, and — when a mask is present — the
// mask, exactly as upstream does inline). Upstream Array3 access F[z][y][x]
// maps to F.get(x, y, z) per the Array3 port. Neumann boundary conditions
// are selected by borderValue === Number.MAX_VALUE (upstream
// std::numeric_limits<Real>::max()).

import { Array3 } from './Array3.js';
import { PdeFilter, PdeFilterScaleType } from './PdeFilter.js';

export abstract class PdeFilter3 extends PdeFilter {
    // Image parameters.
    protected mXBound: number;
    protected mYBound: number;
    protected mZBound: number;
    protected mXSpacing: number;       // dx
    protected mYSpacing: number;       // dy
    protected mZSpacing: number;       // dz
    protected mInvDx: number;          // 1/dx
    protected mInvDy: number;          // 1/dy
    protected mInvDz: number;          // 1/dz
    protected mHalfInvDx: number;      // 1/(2*dx)
    protected mHalfInvDy: number;      // 1/(2*dy)
    protected mHalfInvDz: number;      // 1/(2*dz)
    protected mInvDxDx: number;        // 1/(dx*dx)
    protected mFourthInvDxDy: number;  // 1/(4*dx*dy)
    protected mFourthInvDxDz: number;  // 1/(4*dx*dz)
    protected mInvDyDy: number;        // 1/(dy*dy)
    protected mFourthInvDyDz: number;  // 1/(4*dy*dz)
    protected mInvDzDz: number;        // 1/(dz*dz)

    // Temporary storage for the 3x3x3 neighborhood. In the notation mUxyz,
    // the x, y and z indices are in {m,z,p}, referring to subtract 1 (m),
    // no change (z), or add 1 (p) to the appropriate index.
    protected mUmmm: number; protected mUzmm: number; protected mUpmm: number;
    protected mUmzm: number; protected mUzzm: number; protected mUpzm: number;
    protected mUmpm: number; protected mUzpm: number; protected mUppm: number;
    protected mUmmz: number; protected mUzmz: number; protected mUpmz: number;
    protected mUmzz: number; protected mUzzz: number; protected mUpzz: number;
    protected mUmpz: number; protected mUzpz: number; protected mUppz: number;
    protected mUmmp: number; protected mUzmp: number; protected mUpmp: number;
    protected mUmzp: number; protected mUzzp: number; protected mUpzp: number;
    protected mUmpp: number; protected mUzpp: number; protected mUppp: number;

    // Successive iterations toggle between two buffers.
    protected mBuffer: [Array3<number>, Array3<number>];
    protected mSrc: number;
    protected mDst: number;
    protected mMask: Array3<number>;
    protected mHasMask: boolean;

    protected constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound * yBound * zBound, data, borderValue, scaleType);

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mZBound = zBound;
        this.mXSpacing = xSpacing;
        this.mYSpacing = ySpacing;
        this.mZSpacing = zSpacing;
        this.mInvDx = 1 / xSpacing;
        this.mInvDy = 1 / ySpacing;
        this.mInvDz = 1 / zSpacing;
        this.mHalfInvDx = 0.5 * this.mInvDx;
        this.mHalfInvDy = 0.5 * this.mInvDy;
        this.mHalfInvDz = 0.5 * this.mInvDz;
        this.mInvDxDx = this.mInvDx * this.mInvDx;
        this.mFourthInvDxDy = this.mHalfInvDx * this.mHalfInvDy;
        this.mFourthInvDxDz = this.mHalfInvDx * this.mHalfInvDz;
        this.mInvDyDy = this.mInvDy * this.mInvDy;
        this.mFourthInvDyDz = this.mHalfInvDy * this.mHalfInvDz;
        this.mInvDzDz = this.mInvDz * this.mInvDz;

        this.mUmmm = 0; this.mUzmm = 0; this.mUpmm = 0;
        this.mUmzm = 0; this.mUzzm = 0; this.mUpzm = 0;
        this.mUmpm = 0; this.mUzpm = 0; this.mUppm = 0;
        this.mUmmz = 0; this.mUzmz = 0; this.mUpmz = 0;
        this.mUmzz = 0; this.mUzzz = 0; this.mUpzz = 0;
        this.mUmpz = 0; this.mUzpz = 0; this.mUppz = 0;
        this.mUmmp = 0; this.mUzmp = 0; this.mUpmp = 0;
        this.mUmzp = 0; this.mUzzp = 0; this.mUpzp = 0;
        this.mUmpp = 0; this.mUzpp = 0; this.mUppp = 0;

        this.mSrc = 0;
        this.mDst = 1;
        this.mMask = new Array3<number>(xBound + 2, yBound + 2, zBound + 2);
        this.mMask.fill(0);
        this.mHasMask = (mask !== null);

        // The mBuffer[] are ping-pong buffers for filtering.
        this.mBuffer = [
            new Array3<number>(xBound + 2, yBound + 2, zBound + 2),
            new Array3<number>(xBound + 2, yBound + 2, zBound + 2)
        ];
        this.mBuffer[0].fill(0);
        this.mBuffer[1].fill(0);

        for (let z = 0, zp = 1, i = 0; z < this.mZBound; ++z, ++zp) {
            for (let y = 0, yp = 1; y < this.mYBound; ++y, ++yp) {
                for (let x = 0, xp = 1; x < this.mXBound; ++x, ++xp, ++i) {
                    this.mBuffer[this.mSrc].set(xp, yp, zp,
                        this.mOffset + (data[i] - this.mMin) * this.mScale);
                    this.mBuffer[this.mDst].set(xp, yp, zp, 0);
                    this.mMask.set(xp, yp, zp, this.mHasMask ? mask![i] : 1);
                }
            }
        }

        // Assign values to the 1-voxel thick border.
        if (this.mBorderValue !== Number.MAX_VALUE) {
            this.assignDirichletImageBorder();
        } else {
            this.assignNeumannImageBorder();
        }

        // To handle masks that do not cover the entire image, assign values
        // to those voxels that are 26-neighbors of the mask voxels.
        if (this.mHasMask) {
            if (this.mBorderValue !== Number.MAX_VALUE) {
                this.assignDirichletMaskBorder();
            } else {
                this.assignNeumannMaskBorder();
            }
        }
    }

    // Member access. The internal 3D images for "data" and "mask" are copies
    // of the inputs to the constructor but padded with a 1-voxel thick border
    // to support filtering on the image boundary. These images are of size
    // (xbound+2)-by-(ybound+2)-by-(zbound+2). The correct lookups into the
    // padded arrays are handled internally.
    getXBound(): number {
        return this.mXBound;
    }

    getYBound(): number {
        return this.mYBound;
    }

    getZBound(): number {
        return this.mZBound;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    getYSpacing(): number {
        return this.mYSpacing;
    }

    getZSpacing(): number {
        return this.mZSpacing;
    }

    // Voxel access and derivative estimation. The lookups into the padded
    // data are handled correctly. The estimation involves only the
    // 3-by-3-by-3 neighborhood of (x,y,z), where 0 <= x < xbound,
    // 0 <= y < ybound and 0 <= z < zbound.
    getU(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        return F.get(x + 1, y + 1, z + 1);
    }

    getUx(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const yp1 = y + 1, zp1 = z + 1;
        return this.mHalfInvDx * (F.get(x + 2, yp1, zp1) - F.get(x, yp1, zp1));
    }

    getUy(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp1 = x + 1, zp1 = z + 1;
        return this.mHalfInvDy * (F.get(xp1, y + 2, zp1) - F.get(xp1, y, zp1));
    }

    getUz(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp1 = x + 1, yp1 = y + 1;
        return this.mHalfInvDz * (F.get(xp1, yp1, z + 2) - F.get(xp1, yp1, z));
    }

    getUxx(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const yp1 = y + 1, zp1 = z + 1;
        return this.mInvDxDx * (F.get(x + 2, yp1, zp1)
            - 2 * F.get(x + 1, yp1, zp1) + F.get(x, yp1, zp1));
    }

    getUxy(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp2 = x + 2, yp2 = y + 2, zp1 = z + 1;
        return this.mFourthInvDxDy * (F.get(x, y, zp1) - F.get(xp2, y, zp1)
            + F.get(xp2, yp2, zp1) - F.get(x, yp2, zp1));
    }

    getUxz(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp2 = x + 2, yp1 = y + 1, zp2 = z + 2;
        return this.mFourthInvDxDz * (F.get(x, yp1, z) - F.get(xp2, yp1, z)
            + F.get(xp2, yp1, zp2) - F.get(x, yp1, zp2));
    }

    getUyy(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp1 = x + 1, zp1 = z + 1;
        return this.mInvDyDy * (F.get(xp1, y + 2, zp1)
            - 2 * F.get(xp1, y + 1, zp1) + F.get(xp1, y, zp1));
    }

    getUyz(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp1 = x + 1, yp2 = y + 2, zp2 = z + 2;
        return this.mFourthInvDyDz * (F.get(xp1, y, z) - F.get(xp1, yp2, z)
            + F.get(xp1, yp2, zp2) - F.get(xp1, y, zp2));
    }

    getUzz(x: number, y: number, z: number): number {
        const F = this.mBuffer[this.mSrc];
        const xp1 = x + 1, yp1 = y + 1;
        return this.mInvDzDz * (F.get(xp1, yp1, z + 2)
            - 2 * F.get(xp1, yp1, z + 1) + F.get(xp1, yp1, z));
    }

    getMask(x: number, y: number, z: number): number {
        return this.mMask.get(x + 1, y + 1, z + 1);
    }

    // Assign the value to (x,y,z) of both buffers, and zero out the mask
    // there when a mask exists. This is the common per-site operation of the
    // image-border assignment functions.
    private assignBorderValue(x: number, y: number, z: number, value: number): void {
        this.mBuffer[this.mSrc].set(x, y, z, value);
        this.mBuffer[this.mDst].set(x, y, z, value);
        if (this.mHasMask) {
            this.mMask.set(x, y, z, 0);
        }
    }

    // Assign values to the 1-voxel image border.
    protected assignDirichletImageBorder(): void {
        const xBp1 = this.mXBound + 1, yBp1 = this.mYBound + 1, zBp1 = this.mZBound + 1;
        const value = this.mBorderValue;
        let x: number, y: number, z: number;

        // vertices (0,0,0), (xmax,0,0), (0,ymax,0), (xmax,ymax,0),
        // (0,0,zmax), (xmax,0,zmax), (0,ymax,zmax), (xmax,ymax,zmax)
        this.assignBorderValue(0, 0, 0, value);
        this.assignBorderValue(xBp1, 0, 0, value);
        this.assignBorderValue(0, yBp1, 0, value);
        this.assignBorderValue(xBp1, yBp1, 0, value);
        this.assignBorderValue(0, 0, zBp1, value);
        this.assignBorderValue(xBp1, 0, zBp1, value);
        this.assignBorderValue(0, yBp1, zBp1, value);
        this.assignBorderValue(xBp1, yBp1, zBp1, value);

        // edges (x,0,0) and (x,ymax,0)
        for (x = 1; x <= this.mXBound; ++x) {
            this.assignBorderValue(x, 0, 0, value);
            this.assignBorderValue(x, yBp1, 0, value);
        }

        // edges (0,y,0) and (xmax,y,0)
        for (y = 1; y <= this.mYBound; ++y) {
            this.assignBorderValue(0, y, 0, value);
            this.assignBorderValue(xBp1, y, 0, value);
        }

        // edges (x,0,zmax) and (x,ymax,zmax)
        for (x = 1; x <= this.mXBound; ++x) {
            this.assignBorderValue(x, 0, zBp1, value);
            this.assignBorderValue(x, yBp1, zBp1, value);
        }

        // edges (0,y,zmax) and (xmax,y,zmax)
        for (y = 1; y <= this.mYBound; ++y) {
            this.assignBorderValue(0, y, zBp1, value);
            this.assignBorderValue(xBp1, y, zBp1, value);
        }

        // edges (0,0,z) and (xmax,0,z)
        for (z = 1; z <= this.mZBound; ++z) {
            this.assignBorderValue(0, 0, z, value);
            this.assignBorderValue(xBp1, 0, z, value);
        }

        // edges (0,ymax,z) and (xmax,ymax,z)
        for (z = 1; z <= this.mZBound; ++z) {
            this.assignBorderValue(0, yBp1, z, value);
            this.assignBorderValue(xBp1, yBp1, z, value);
        }

        // faces (x,y,0) and (x,y,zmax)
        for (y = 1; y <= this.mYBound; ++y) {
            for (x = 1; x <= this.mXBound; ++x) {
                this.assignBorderValue(x, y, 0, value);
                this.assignBorderValue(x, y, zBp1, value);
            }
        }

        // faces (x,0,z) and (x,ymax,z)
        for (z = 1; z <= this.mZBound; ++z) {
            for (x = 1; x <= this.mXBound; ++x) {
                this.assignBorderValue(x, 0, z, value);
                this.assignBorderValue(x, yBp1, z, value);
            }
        }

        // faces (0,y,z) and (xmax,y,z)
        for (z = 1; z <= this.mZBound; ++z) {
            for (y = 1; y <= this.mYBound; ++y) {
                this.assignBorderValue(0, y, z, value);
                this.assignBorderValue(xBp1, y, z, value);
            }
        }
    }

    protected assignNeumannImageBorder(): void {
        const xB = this.mXBound, yB = this.mYBound, zB = this.mZBound;
        const xBp1 = xB + 1, yBp1 = yB + 1, zBp1 = zB + 1;
        const src = this.mBuffer[this.mSrc];
        let x: number, y: number, z: number;
        let duplicate: number;

        // Duplicate the nearest interior value at each border voxel. The
        // reads are all interior sites and the writes are all border sites,
        // so the order of assignment does not affect the duplicates.

        // vertex (0,0,0)
        duplicate = src.get(1, 1, 1);
        this.assignBorderValue(0, 0, 0, duplicate);

        // vertex (xmax,0,0)
        duplicate = src.get(xB, 1, 1);
        this.assignBorderValue(xBp1, 0, 0, duplicate);

        // vertex (0,ymax,0)
        duplicate = src.get(1, yB, 1);
        this.assignBorderValue(0, yBp1, 0, duplicate);

        // vertex (xmax,ymax,0)
        duplicate = src.get(xB, yB, 1);
        this.assignBorderValue(xBp1, yBp1, 0, duplicate);

        // vertex (0,0,zmax)
        duplicate = src.get(1, 1, zB);
        this.assignBorderValue(0, 0, zBp1, duplicate);

        // vertex (xmax,0,zmax)
        duplicate = src.get(xB, 1, zB);
        this.assignBorderValue(xBp1, 0, zBp1, duplicate);

        // vertex (0,ymax,zmax)
        duplicate = src.get(1, yB, zB);
        this.assignBorderValue(0, yBp1, zBp1, duplicate);

        // vertex (xmax,ymax,zmax)
        duplicate = src.get(xB, yB, zB);
        this.assignBorderValue(xBp1, yBp1, zBp1, duplicate);

        // edges (x,0,0) and (x,ymax,0)
        for (x = 1; x <= xB; ++x) {
            this.assignBorderValue(x, 0, 0, src.get(x, 1, 1));
            this.assignBorderValue(x, yBp1, 0, src.get(x, yB, 1));
        }

        // edges (0,y,0) and (xmax,y,0)
        for (y = 1; y <= yB; ++y) {
            this.assignBorderValue(0, y, 0, src.get(1, y, 1));
            this.assignBorderValue(xBp1, y, 0, src.get(xB, y, 1));
        }

        // edges (x,0,zmax) and (x,ymax,zmax)
        for (x = 1; x <= xB; ++x) {
            this.assignBorderValue(x, 0, zBp1, src.get(x, 1, zB));
            this.assignBorderValue(x, yBp1, zBp1, src.get(x, yB, zB));
        }

        // edges (0,y,zmax) and (xmax,y,zmax)
        for (y = 1; y <= yB; ++y) {
            this.assignBorderValue(0, y, zBp1, src.get(1, y, zB));
            this.assignBorderValue(xBp1, y, zBp1, src.get(xB, y, zB));
        }

        // edges (0,0,z) and (xmax,0,z)
        for (z = 1; z <= zB; ++z) {
            this.assignBorderValue(0, 0, z, src.get(1, 1, z));
            this.assignBorderValue(xBp1, 0, z, src.get(xB, 1, z));
        }

        // edges (0,ymax,z) and (xmax,ymax,z)
        for (z = 1; z <= zB; ++z) {
            this.assignBorderValue(0, yBp1, z, src.get(1, yB, z));
            this.assignBorderValue(xBp1, yBp1, z, src.get(xB, yB, z));
        }

        // faces (x,y,0) and (x,y,zmax)
        for (y = 1; y <= yB; ++y) {
            for (x = 1; x <= xB; ++x) {
                this.assignBorderValue(x, y, 0, src.get(x, y, 1));
                this.assignBorderValue(x, y, zBp1, src.get(x, y, zB));
            }
        }

        // faces (x,0,z) and (x,ymax,z)
        for (z = 1; z <= zB; ++z) {
            for (x = 1; x <= xB; ++x) {
                this.assignBorderValue(x, 0, z, src.get(x, 1, z));
                this.assignBorderValue(x, yBp1, z, src.get(x, yB, z));
            }
        }

        // faces (0,y,z) and (xmax,y,z)
        for (z = 1; z <= zB; ++z) {
            for (y = 1; y <= yB; ++y) {
                this.assignBorderValue(0, y, z, src.get(1, y, z));
                this.assignBorderValue(xBp1, y, z, src.get(xB, y, z));
            }
        }
    }

    // Assign values to the 1-voxel mask border.
    protected assignDirichletMaskBorder(): void {
        for (let z = 1; z <= this.mZBound; ++z) {
            for (let y = 1; y <= this.mYBound; ++y) {
                for (let x = 1; x <= this.mXBound; ++x) {
                    if (this.mMask.get(x, y, z) !== 0) {
                        continue;
                    }

                    let found = false;
                    for (let i2 = 0, j2 = z - 1; i2 < 3 && !found; ++i2, ++j2) {
                        for (let i1 = 0, j1 = y - 1; i1 < 3 && !found; ++i1, ++j1) {
                            for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                                if (this.mMask.get(j0, j1, j2) !== 0) {
                                    this.mBuffer[this.mSrc].set(x, y, z, this.mBorderValue);
                                    this.mBuffer[this.mDst].set(x, y, z, this.mBorderValue);
                                    found = true;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    protected assignNeumannMaskBorder(): void {
        // Recompute the values just outside the masked region. This
        // guarantees that derivative estimations use the current values
        // around the boundary.
        for (let z = 1; z <= this.mZBound; ++z) {
            for (let y = 1; y <= this.mYBound; ++y) {
                for (let x = 1; x <= this.mXBound; ++x) {
                    if (this.mMask.get(x, y, z) !== 0) {
                        continue;
                    }

                    let count = 0;
                    let average = 0;
                    for (let i2 = 0, j2 = z - 1; i2 < 3; ++i2, ++j2) {
                        for (let i1 = 0, j1 = y - 1; i1 < 3; ++i1, ++j1) {
                            for (let i0 = 0, j0 = x - 1; i0 < 3; ++i0, ++j0) {
                                if (this.mMask.get(j0, j1, j2) !== 0) {
                                    average += this.mBuffer[this.mSrc].get(j0, j1, j2);
                                    count++;
                                }
                            }
                        }
                    }

                    if (count > 0) {
                        average /= count;
                        this.mBuffer[this.mSrc].set(x, y, z, average);
                        this.mBuffer[this.mDst].set(x, y, z, average);
                    }
                }
            }
        }
    }

    // This function recomputes the boundary values when Neumann conditions
    // are used. If a derived class overrides this, it must call the
    // base-class onPreUpdate first.
    protected override onPreUpdate(): void {
        if (this.mHasMask && this.mBorderValue === Number.MAX_VALUE) {
            // Neumann boundary conditions are in use, so recompute the mask
            // border.
            this.assignNeumannMaskBorder();
        }
        // else: No mask has been specified or Dirichlet boundary conditions
        // are in use. Nothing to do.
    }

    // Iterate over all the voxels and call onUpdateSingle(x,y,z) for each
    // voxel that is not masked out.
    protected override onUpdate(): void {
        for (let z = 1; z <= this.mZBound; ++z) {
            for (let y = 1; y <= this.mYBound; ++y) {
                for (let x = 1; x <= this.mXBound; ++x) {
                    if (!this.mHasMask || this.mMask.get(x, y, z) !== 0) {
                        this.onUpdateSingle(x, y, z);
                    }
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

    // The per-voxel processing depends on the PDE algorithm. The (x,y,z)
    // must be in padded coordinates: 1 <= x <= xbound, 1 <= y <= ybound and
    // 1 <= z <= zbound.
    protected abstract onUpdateSingle(x: number, y: number, z: number): void;

    // Copy source data to temporary storage.
    protected lookUp7(x: number, y: number, z: number): void {
        const F = this.mBuffer[this.mSrc];
        const xm = x - 1, xp = x + 1;
        const ym = y - 1, yp = y + 1;
        const zm = z - 1, zp = z + 1;
        this.mUzzm = F.get(x, y, zm);
        this.mUzmz = F.get(x, ym, z);
        this.mUmzz = F.get(xm, y, z);
        this.mUzzz = F.get(x, y, z);
        this.mUpzz = F.get(xp, y, z);
        this.mUzpz = F.get(x, yp, z);
        this.mUzzp = F.get(x, y, zp);
    }

    protected lookUp27(x: number, y: number, z: number): void {
        const F = this.mBuffer[this.mSrc];
        const xm = x - 1, xp = x + 1;
        const ym = y - 1, yp = y + 1;
        const zm = z - 1, zp = z + 1;
        this.mUmmm = F.get(xm, ym, zm);
        this.mUzmm = F.get(x, ym, zm);
        this.mUpmm = F.get(xp, ym, zm);
        this.mUmzm = F.get(xm, y, zm);
        this.mUzzm = F.get(x, y, zm);
        this.mUpzm = F.get(xp, y, zm);
        this.mUmpm = F.get(xm, yp, zm);
        this.mUzpm = F.get(x, yp, zm);
        this.mUppm = F.get(xp, yp, zm);
        this.mUmmz = F.get(xm, ym, z);
        this.mUzmz = F.get(x, ym, z);
        this.mUpmz = F.get(xp, ym, z);
        this.mUmzz = F.get(xm, y, z);
        this.mUzzz = F.get(x, y, z);
        this.mUpzz = F.get(xp, y, z);
        this.mUmpz = F.get(xm, yp, z);
        this.mUzpz = F.get(x, yp, z);
        this.mUppz = F.get(xp, yp, z);
        this.mUmmp = F.get(xm, ym, zp);
        this.mUzmp = F.get(x, ym, zp);
        this.mUpmp = F.get(xp, ym, zp);
        this.mUmzp = F.get(xm, y, zp);
        this.mUzzp = F.get(x, y, zp);
        this.mUpzp = F.get(xp, y, zp);
        this.mUmpp = F.get(xm, yp, zp);
        this.mUzpp = F.get(x, yp, zp);
        this.mUppp = F.get(xp, yp, zp);
    }
}
