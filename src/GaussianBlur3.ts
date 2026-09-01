// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GaussianBlur3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Gaussian blurring of a 3D image, implemented as an explicit Euler solver
// for the linear heat equation u_t = u_xx + u_yy + u_zz.
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// mBuffer[mDst][z][y][x] becomes mBuffer[mDst].set(x, y, z, ...) per the
// Array3 port.

import { PdeFilter3 } from './PdeFilter3';
import { PdeFilterScaleType } from './PdeFilter';

export class GaussianBlur3 extends PdeFilter3 {
    // The largest time step for which the explicit Euler update is stable.
    protected mMaximumTimeStep: number;

    constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound, yBound, zBound, xSpacing, ySpacing, zSpacing,
            data, mask, borderValue, scaleType);
        this.mMaximumTimeStep = 0.5 / (this.mInvDxDx + this.mInvDyDy + this.mInvDzDz);
    }

    getMaximumTimeStep(): number {
        return this.mMaximumTimeStep;
    }

    protected override onUpdateSingle(x: number, y: number, z: number): void {
        this.lookUp7(x, y, z);

        const uxx = this.mInvDxDx * (this.mUpzz - 2 * this.mUzzz + this.mUmzz);
        const uyy = this.mInvDyDy * (this.mUzpz - 2 * this.mUzzz + this.mUzmz);
        const uzz = this.mInvDzDz * (this.mUzzp - 2 * this.mUzzz + this.mUzzm);

        this.mBuffer[this.mDst].set(x, y, z, this.mUzzz + this.mTimeStep * (uxx + uyy + uzz));
    }
}
