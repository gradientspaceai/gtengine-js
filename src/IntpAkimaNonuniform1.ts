// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpAkimaNonuniform1.h
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
// { index, dx }, per the IntpAkima1 base. The upstream null-pointer test on
// X becomes a length test.

import { IntpAkima1 } from './IntpAkima1.js';
import { logAssert } from './Logger.js';

export class IntpAkimaNonuniform1 extends IntpAkima1 {
    protected mX: readonly number[];

    // Construction. The interpolator is for arbitrarily spaced x-values.
    // The input arrays must have 'quantity' elements and the X[] array must
    // store increasing values: X[i + 1] > X[i] for all i. The arrays X and
    // F are aliased, not copied.
    constructor(quantity: number, X: readonly number[], F: readonly number[]) {
        super(quantity, F);

        logAssert(X.length >= quantity, 'Invalid input.');
        for (let j0 = 0, j1 = 1; j1 < quantity; ++j0, ++j1) {
            logAssert(X[j1] > X[j0], 'Invalid input.');
        }
        this.mX = X;

        // Compute slopes. The slope array has quantity+3 elements. The
        // slopes of the quantity-1 sample intervals are stored at indices 2
        // through quantity. Indices 0, 1, quantity+1 and quantity+2 store
        // the extrapolated slopes that are needed to estimate the
        // derivatives at the boundary samples.
        const slope = new Array<number>(quantity + 3).fill(0);
        let i: number, ip1: number, ip2: number;
        for (i = 0, ip1 = 1, ip2 = 2; i < quantity - 1; ++i, ++ip1, ++ip2) {
            const dx = X[ip1] - X[i];
            const df = F[ip1] - F[i];
            slope[ip2] = df / dx;
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
        for (i = 0, ip1 = 1; i < quantity - 1; ++i, ++ip1) {
            const poly = this.mPoly[i];

            const F0 = F[i];
            const F1 = F[ip1];
            const FDer0 = FDer[i];
            const FDer1 = FDer[ip1];
            const df = F1 - F0;
            const dx = X[ip1] - X[i];
            const dx2 = dx * dx;
            const dx3 = dx2 * dx;

            poly.setCoeff(0, F0);
            poly.setCoeff(1, FDer0);
            poly.setCoeff(2, (3 * df - dx * (FDer1 + 2 * FDer0)) / dx2);
            poly.setCoeff(3, (dx * (FDer0 + FDer1) - 2 * df) / dx3);
        }
    }

    // Member access.
    getX(): readonly number[] {
        return this.mX;
    }

    override getXMin(): number {
        return this.mX[0];
    }

    override getXMax(): number {
        return this.mX[this.mQuantity - 1];
    }

    protected override lookup(x: number): { index: number, dx: number } {
        // The caller has ensured that mXMin <= x <= mXMax.
        let index: number;
        for (index = 0; index + 1 < this.mQuantity; ++index) {
            if (x < this.mX[index + 1]) {
                return { index, dx: x - this.mX[index] };
            }
        }

        --index;
        return { index, dx: x - this.mX[index] };
    }
}
