// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpAkimaUniform1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Akima interpolation is described in
// https://en.wikipedia.org/wiki/Akima_spline
//
// Port notes: the interpolator derives from IntpAkima1, so the overloaded
// evaluate(x) and evaluate(order, x) are inherited. The virtual
// 'Lookup(x, index&, dx&)' is ported as 'lookup(x)' returning
// { index, dx }, per the IntpAkima1 base.

import { IntpAkima1 } from './IntpAkima1';
import { logAssert } from './Logger';

export class IntpAkimaUniform1 extends IntpAkima1 {
    protected mXMin: number;
    protected mXMax: number;
    protected mXSpacing: number;

    // Construction. The interpolator is for uniformly spaced x-values. The
    // samples F are aliased, not copied.
    constructor(quantity: number, xMin: number, xSpacing: number, F: readonly number[]) {
        super(quantity, F);

        logAssert(xSpacing > 0, 'Spacing must be positive.');
        this.mXMin = xMin;
        this.mXSpacing = xSpacing;
        this.mXMax = xMin + xSpacing * (quantity - 1);

        // Compute slopes. The slope array has quantity+3 elements. The
        // slopes of the quantity-1 sample intervals are stored at indices 2
        // through quantity. Indices 0, 1, quantity+1 and quantity+2 store
        // the extrapolated slopes that are needed to estimate the
        // derivatives at the boundary samples.
        const invDX = 1 / xSpacing;
        const slope = new Array<number>(quantity + 3).fill(0);
        let i: number, ip1: number, ip2: number;
        for (i = 0, ip1 = 1, ip2 = 2; i < quantity - 1; ++i, ++ip1, ++ip2) {
            slope[ip2] = (this.mF[ip1] - this.mF[i]) * invDX;
        }

        slope[1] = 2 * slope[2] - slope[3];
        slope[0] = 2 * slope[1] - slope[2];
        slope[quantity + 1] = 2 * slope[quantity] - slope[quantity - 1];
        slope[quantity + 2] = 2 * slope[quantity + 1] - slope[quantity];

        // Construct derivatives.
        const FDer = new Array<number>(quantity).fill(0);
        for (i = 0; i < quantity; ++i) {
            FDer[i] = this.computeDerivative(slope, i);
        }

        // Construct polynomials.
        const invDX2 = 1 / (xSpacing * xSpacing);
        const invDX3 = invDX2 / xSpacing;
        for (i = 0, ip1 = 1; i < quantity - 1; ++i, ++ip1) {
            const poly = this.mPoly[i];

            const F0 = F[i];
            const F1 = F[ip1];
            const df = F1 - F0;
            const FDer0 = FDer[i];
            const FDer1 = FDer[ip1];

            poly.setCoeff(0, F0);
            poly.setCoeff(1, FDer0);
            poly.setCoeff(2, (3 * df - xSpacing * (FDer1 + 2 * FDer0)) * invDX2);
            poly.setCoeff(3, (xSpacing * (FDer0 + FDer1) - 2 * df) * invDX3);
        }
    }

    // Member access.
    override getXMin(): number {
        return this.mXMin;
    }

    override getXMax(): number {
        return this.mXMax;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    protected override lookup(x: number): { index: number, dx: number } {
        // The caller has ensured that mXMin <= x <= mXMax.
        let index: number, indexP1: number;
        for (index = 0, indexP1 = 1; indexP1 < this.mQuantity; ++index, ++indexP1) {
            if (x < this.mXMin + this.mXSpacing * indexP1) {
                return { index, dx: x - (this.mXMin + this.mXSpacing * index) };
            }
        }

        --index;
        return { index, dx: x - (this.mXMin + this.mXSpacing * index) };
    }
}
