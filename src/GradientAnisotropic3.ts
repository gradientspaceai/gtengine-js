// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GradientAnisotropic3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Perona-Malik gradient-magnitude anisotropic diffusion for a 3D image,
//   u_t = div(C(|grad(u)|) * grad(u)),
// where the conductance is C(s) = exp(-s^2/(2*k^2*avg)) and avg is the
// average of the squared gradient magnitudes over the image.
//
// Port notes: the template parameter Real is number (IEEE double). The C++
// mBuffer[mDst][z][y][x] becomes mBuffer[mDst].set(x, y, z, ...) per the
// Array3 port.
//
// Upstream bug (fixed here, same as GradientAnisotropic2): ComputeParameter
// iterates over the padded coordinates 1 <= x <= mXBound (and likewise for y
// and z) but passes them to GetUx/GetUy/GetUz, which add the padding offset
// themselves. The average is taken over a window shifted by one voxel and the
// final sample reads past the end of the padded buffer (undefined behavior in
// C++; undefined -> NaN in TypeScript). The loop here uses the unpadded
// coordinates the accessors document.

import { PdeFilter3 } from './PdeFilter3.js';
import { PdeFilterScaleType } from './PdeFilter.js';

export class GradientAnisotropic3 extends PdeFilter3 {
    // These are updated on each iteration, since they depend on the current
    // average of the squared length of the gradients at the voxels.
    protected mK: number;               // k
    protected mParameter: number;       // 1/(k^2*average(gradMagSqr))
    protected mMHalfParameter: number;  // -0.5*mParameter

    constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        data: ArrayLike<number>, mask: ArrayLike<number> | null,
        borderValue: number, scaleType: PdeFilterScaleType, K: number) {
        super(xBound, yBound, zBound, xSpacing, ySpacing, zSpacing,
            data, mask, borderValue, scaleType);
        this.mK = K;
        this.mParameter = 0;
        this.mMHalfParameter = 0;
        this.computeParameter();
    }

    protected computeParameter(): void {
        let gradMagSqr = 0;
        for (let z = 0; z < this.mZBound; ++z) {
            for (let y = 0; y < this.mYBound; ++y) {
                for (let x = 0; x < this.mXBound; ++x) {
                    const ux = this.getUx(x, y, z);
                    const uy = this.getUy(x, y, z);
                    const uz = this.getUz(x, y, z);
                    gradMagSqr += ux * ux + uy * uy + uz * uz;
                }
            }
        }
        gradMagSqr /= this.mQuantity;

        this.mParameter = 1 / (this.mK * this.mK * gradMagSqr);
        this.mMHalfParameter = -0.5 * this.mParameter;
    }

    protected override onPreUpdate(): void {
        this.computeParameter();
    }

    protected override onUpdateSingle(x: number, y: number, z: number): void {
        this.lookUp27(x, y, z);

        // one-sided U-derivative estimates
        const uxFwd = this.mInvDx * (this.mUpzz - this.mUzzz);
        const uxBwd = this.mInvDx * (this.mUzzz - this.mUmzz);
        const uyFwd = this.mInvDy * (this.mUzpz - this.mUzzz);
        const uyBwd = this.mInvDy * (this.mUzzz - this.mUzmz);
        const uzFwd = this.mInvDz * (this.mUzzp - this.mUzzz);
        const uzBwd = this.mInvDz * (this.mUzzz - this.mUzzm);

        // centered U-derivative estimates
        const duvzz = this.mHalfInvDx * (this.mUpzz - this.mUmzz);
        const duvpz = this.mHalfInvDx * (this.mUppz - this.mUmpz);
        const duvmz = this.mHalfInvDx * (this.mUpmz - this.mUmmz);
        const duvzp = this.mHalfInvDx * (this.mUpzp - this.mUmzp);
        const duvzm = this.mHalfInvDx * (this.mUpzm - this.mUmzm);

        const duzvz = this.mHalfInvDy * (this.mUzpz - this.mUzmz);
        const dupvz = this.mHalfInvDy * (this.mUppz - this.mUpmz);
        const dumvz = this.mHalfInvDy * (this.mUmpz - this.mUmmz);
        const duzvp = this.mHalfInvDy * (this.mUzpp - this.mUzmp);
        const duzvm = this.mHalfInvDy * (this.mUzpm - this.mUzmm);

        const duzzv = this.mHalfInvDz * (this.mUzzp - this.mUzzm);
        const dupzv = this.mHalfInvDz * (this.mUpzp - this.mUpzm);
        const dumzv = this.mHalfInvDz * (this.mUmzp - this.mUmzm);
        const duzpv = this.mHalfInvDz * (this.mUzpp - this.mUzpm);
        const duzmv = this.mHalfInvDz * (this.mUzmp - this.mUzmm);

        const uxCenSqr = duvzz * duvzz;
        const uyCenSqr = duzvz * duzvz;
        const uzCenSqr = duzzv * duzzv;

        let uxEst: number, uyEst: number, uzEst: number, gradMagSqr: number;

        // estimate for C(x+1,y,z)
        uyEst = 0.5 * (duzvz + dupvz);
        uzEst = 0.5 * (duzzv + dupzv);
        gradMagSqr = uxCenSqr + uyEst * uyEst + uzEst * uzEst;
        const cxp = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x-1,y,z)
        uyEst = 0.5 * (duzvz + dumvz);
        uzEst = 0.5 * (duzzv + dumzv);
        gradMagSqr = uxCenSqr + uyEst * uyEst + uzEst * uzEst;
        const cxm = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y+1,z)
        uxEst = 0.5 * (duvzz + duvpz);
        uzEst = 0.5 * (duzzv + duzpv);
        gradMagSqr = uxEst * uxEst + uyCenSqr + uzEst * uzEst;
        const cyp = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y-1,z)
        uxEst = 0.5 * (duvzz + duvmz);
        uzEst = 0.5 * (duzzv + duzmv);
        gradMagSqr = uxEst * uxEst + uyCenSqr + uzEst * uzEst;
        const cym = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y,z+1)
        uxEst = 0.5 * (duvzz + duvzp);
        uyEst = 0.5 * (duzvz + duzvp);
        gradMagSqr = uxEst * uxEst + uyEst * uyEst + uzCenSqr;
        const czp = Math.exp(this.mMHalfParameter * gradMagSqr);

        // estimate for C(x,y,z-1)
        uxEst = 0.5 * (duvzz + duvzm);
        uyEst = 0.5 * (duzvz + duzvm);
        gradMagSqr = uxEst * uxEst + uyEst * uyEst + uzCenSqr;
        const czm = Math.exp(this.mMHalfParameter * gradMagSqr);

        this.mBuffer[this.mDst].set(x, y, z, this.mUzzz + this.mTimeStep * (
            cxp * uxFwd - cxm * uxBwd +
            cyp * uyFwd - cym * uyBwd +
            czp * uzFwd - czm * uzBwd));
    }
}
