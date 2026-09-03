// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GradientAnisotropic2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Perona-Malik gradient-magnitude anisotropic diffusion for a 2D image,
//   u_t = div(C(|grad(u)|) * grad(u)),
// where the conductance is C(s) = exp(-s^2/(2*k^2*avg)) and avg is the
// average of the squared gradient magnitudes over the image. The conductance
// is small where the gradient magnitude is large, so edges diffuse much more
// slowly than the interiors of regions.
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// mBuffer[mDst][y][x] becomes mBuffer[mDst].set(x, y, ...) per the Array2
// port.
//
// Upstream bug (fixed here): ComputeParameter iterates
// 1 <= x <= mXBound and 1 <= y <= mYBound and passes those *padded*
// coordinates to GetUx/GetUy, which themselves add the 1-pixel padding
// offset. The average is therefore taken over a window shifted by one pixel,
// and the last sample reads one element past the end of the padded buffer
// (undefined behavior in C++; undefined -> NaN in TypeScript, which would
// poison the whole filter). The loop here uses the unpadded coordinates
// 0 <= x < mXBound and 0 <= y < mYBound that GetUx/GetUy document, which also
// matches the division by mQuantity = xBound*yBound.

import { PdeFilter2 } from './PdeFilter2.js';
import { PdeFilterScaleType } from './PdeFilter.js';

export class GradientAnisotropic2 extends PdeFilter2 {
    // These are updated on each iteration, since they depend on the current
    // average of the squared length of the gradients at the pixels.
    protected mK: number;               // k
    protected mParameter: number;       // 1/(k^2*average(gradMagSqr))
    protected mMHalfParameter: number;  // -0.5*mParameter

    constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null, borderValue: number,
        scaleType: PdeFilterScaleType, K: number) {
        super(xBound, yBound, xSpacing, ySpacing, data, mask, borderValue, scaleType);
        this.mK = K;
        this.mParameter = 0;
        this.mMHalfParameter = 0;
        this.computeParameter();
    }

    protected computeParameter(): void {
        let gradMagSqr = 0;
        for (let y = 0; y < this.mYBound; ++y) {
            for (let x = 0; x < this.mXBound; ++x) {
                const ux = this.getUx(x, y);
                const uy = this.getUy(x, y);
                gradMagSqr += ux * ux + uy * uy;
            }
        }
        gradMagSqr /= this.mQuantity;

        this.mParameter = 1 / (this.mK * this.mK * gradMagSqr);
        this.mMHalfParameter = -0.5 * this.mParameter;
    }

    protected override onPreUpdate(): void {
        this.computeParameter();
    }

    protected override onUpdateSingle(x: number, y: number): void {
        this.lookUp9(x, y);

        // one-sided U-derivative estimates
        const uxFwd = this.mInvDx * (this.mUpz - this.mUzz);
        const uxBwd = this.mInvDx * (this.mUzz - this.mUmz);
        const uyFwd = this.mInvDy * (this.mUzp - this.mUzz);
        const uyBwd = this.mInvDy * (this.mUzz - this.mUzm);

        // centered U-derivative estimates
        const uxCenM = this.mHalfInvDx * (this.mUpm - this.mUmm);
        const uxCenZ = this.mHalfInvDx * (this.mUpz - this.mUmz);
        const uxCenP = this.mHalfInvDx * (this.mUpp - this.mUmp);
        const uyCenM = this.mHalfInvDy * (this.mUmp - this.mUmm);
        const uyCenZ = this.mHalfInvDy * (this.mUzp - this.mUzm);
        const uyCenP = this.mHalfInvDy * (this.mUpp - this.mUpm);

        const uxCenZSqr = uxCenZ * uxCenZ;
        const uyCenZSqr = uyCenZ * uyCenZ;
        let gradMagSqr: number;

        // estimate for C(x+1,y)
        const uyEstP = 0.5 * (uyCenZ + uyCenP);
        gradMagSqr = uxCenZSqr + uyEstP * uyEstP;
        const cxp = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x-1,y)
        const uyEstM = 0.5 * (uyCenZ + uyCenM);
        gradMagSqr = uxCenZSqr + uyEstM * uyEstM;
        const cxm = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y+1)
        const uxEstP = 0.5 * (uxCenZ + uxCenP);
        gradMagSqr = uyCenZSqr + uxEstP * uxEstP;
        const cyp = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y-1)
        const uxEstM = 0.5 * (uxCenZ + uxCenM);
        gradMagSqr = uyCenZSqr + uxEstM * uxEstM;
        const cym = Math.exp(this.mMHalfParameter * gradMagSqr);

        this.mBuffer[this.mDst].set(x, y, this.mUzz + this.mTimeStep * (
            cxp * uxFwd - cxm * uxBwd +
            cyp * uyFwd - cym * uyBwd));
    }
}
