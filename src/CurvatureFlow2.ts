// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CurvatureFlow2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Curvature flow for a 2D image. The level curves of the image evolve
// according to u_t = |grad(u)| * div(grad(u)/|grad(u)|), which moves each
// level curve in its normal direction with speed equal to its curvature.
// The right-hand side expands to
//   (u_xx * u_y^2 + u_yy * u_x^2 - 2 * u_xy * u_x * u_y) / (u_x^2 + u_y^2).
// Upstream uses a coefficient of -1/2 rather than -2 on the mixed term; the
// quirk is preserved (see the port notes in the PR).
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// mBuffer[mDst][y][x] becomes mBuffer[mDst].set(x, y, ...) per the Array2
// port.

import { PdeFilter2 } from './PdeFilter2';
import { PdeFilterScaleType } from './PdeFilter';

export class CurvatureFlow2 extends PdeFilter2 {
    constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType) {
        super(xBound, yBound, xSpacing, ySpacing, data, mask, borderValue, scaleType);
    }

    protected override onUpdateSingle(x: number, y: number): void {
        this.lookUp9(x, y);

        const ux = this.mHalfInvDx * (this.mUpz - this.mUmz);
        const uy = this.mHalfInvDy * (this.mUzp - this.mUzm);
        const uxx = this.mInvDxDx * (this.mUpz - 2 * this.mUzz + this.mUmz);
        const uxy = this.mFourthInvDxDy * (this.mUmm + this.mUpp - this.mUmp - this.mUpm);
        const uyy = this.mInvDyDy * (this.mUzp - 2 * this.mUzz + this.mUzm);

        const sqrUx = ux * ux;
        const sqrUy = uy * uy;
        const denom = sqrUx + sqrUy;
        if (denom > 0) {
            const numer = uxx * sqrUy + uyy * sqrUx - 0.5 * uxy * ux * uy;
            this.mBuffer[this.mDst].set(x, y, this.mUzz + this.mTimeStep * numer / denom);
        } else {
            this.mBuffer[this.mDst].set(x, y, this.mUzz);
        }
    }
}
