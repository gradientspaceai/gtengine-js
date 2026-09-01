// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GaussianBlur2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Gaussian blurring of a 2D image, implemented as an explicit Euler solver
// for the linear heat equation u_t = u_xx + u_yy. The solution of the heat
// equation at time t is the convolution of the initial image with a Gaussian
// kernel of variance 2*t, whence the name.
//
// Port notes: the template parameter Real is number (IEEE double). The
// upstream protected member mMaximumTimeStep is exposed by
// getMaximumTimeStep, as upstream. The C++ mBuffer[mDst][y][x] becomes
// mBuffer[mDst].set(x, y, ...) per the Array2 port.

import { PdeFilter2 } from './PdeFilter2';
import { PdeFilterScaleType } from './PdeFilter';

export class GaussianBlur2 extends PdeFilter2 {
    // The largest time step for which the explicit Euler update is stable.
    protected mMaximumTimeStep: number;

    constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType) {
        super(xBound, yBound, xSpacing, ySpacing, data, mask, borderValue, scaleType);
        this.mMaximumTimeStep = 0.5 / (this.mInvDxDx + this.mInvDyDy);
    }

    getMaximumTimeStep(): number {
        return this.mMaximumTimeStep;
    }

    protected override onUpdateSingle(x: number, y: number): void {
        this.lookUp5(x, y);

        const uxx = this.mInvDxDx * (this.mUpz - 2 * this.mUzz + this.mUmz);
        const uyy = this.mInvDyDy * (this.mUzp - 2 * this.mUzz + this.mUzm);

        this.mBuffer[this.mDst].set(x, y, this.mUzz + this.mTimeStep * (uxx + uyy));
    }
}
