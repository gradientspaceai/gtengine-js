// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CurvatureFlow3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Curvature flow for a 3D image. The level surfaces of the image evolve in
// their normal direction with speed proportional to the mean curvature. The
// numerator is the sum of the three 2D curvature-flow numerators, one for
// each of the coordinate planes.
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// mBuffer[mDst][z][y][x] becomes mBuffer[mDst].set(x, y, z, ...) per the
// Array3 port.

import { PdeFilter3 } from './PdeFilter3';
import { PdeFilterScaleType } from './PdeFilter';

export class CurvatureFlow3 extends PdeFilter3 {
    constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound, yBound, zBound, xSpacing, ySpacing, zSpacing,
            data, mask, borderValue, scaleType);
    }

    protected override onUpdateSingle(x: number, y: number, z: number): void {
        this.lookUp27(x, y, z);

        const ux = this.mHalfInvDx * (this.mUpzz - this.mUmzz);
        const uy = this.mHalfInvDy * (this.mUzpz - this.mUzmz);
        const uz = this.mHalfInvDz * (this.mUzzp - this.mUzzm);
        const uxx = this.mInvDxDx * (this.mUpzz - 2 * this.mUzzz + this.mUmzz);
        const uxy = this.mFourthInvDxDy * (this.mUmmz + this.mUppz - this.mUpmz - this.mUmpz);
        const uxz = this.mFourthInvDxDz * (this.mUmzm + this.mUpzp - this.mUpzm - this.mUmzp);
        const uyy = this.mInvDyDy * (this.mUzpz - 2 * this.mUzzz + this.mUzmz);
        const uyz = this.mFourthInvDyDz * (this.mUzmm + this.mUzpp - this.mUzpm - this.mUzmp);
        const uzz = this.mInvDzDz * (this.mUzzp - 2 * this.mUzzz + this.mUzzm);

        const denom = ux * ux + uy * uy + uz * uz;
        if (denom > 0) {
            const numer0 = uy * (uxx * uy - uxy * ux) + ux * (uyy * ux - uxy * uy);
            const numer1 = uz * (uxx * uz - uxz * ux) + ux * (uzz * ux - uxz * uz);
            const numer2 = uz * (uyy * uz - uyz * uy) + uy * (uzz * uy - uyz * uz);
            const numer = numer0 + numer1 + numer2;
            this.mBuffer[this.mDst].set(x, y, z, this.mUzzz + this.mTimeStep * numer / denom);
        } else {
            this.mBuffer[this.mDst].set(x, y, z, this.mUzzz);
        }
    }
}
