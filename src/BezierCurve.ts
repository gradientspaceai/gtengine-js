// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BezierCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A Bezier curve X(t) = sum_{i=0}^{d} Choose(d,i) * (1-t)^{d-i} * t^i * C[i],
// where d is the degree and the C[i] are the d+1 control points. The domain
// is t in [0,1].
//
// Port notes (following BSplineCurve): upstream
// 'template <int32_t N, typename Real>' becomes a runtime dimension passed as
// the first constructor argument and carried by the control points (Vector
// objects). The control points are copied (C++ value semantics). The
// combinatorial table Array2<Real> mChoose stores C++ 'mChoose[r][c]' as
// 'mChoose.get(c, r)', following the Array2 index-order convention.
//
// As upstream, the constructor requires degree >= 2; the first- through
// third-order difference tables that drive the derivative evaluation are not
// built for degree 1.

import { Array2 } from './Array2.js';
import { logAssert } from './Logger.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector } from './Vector.js';

export class BezierCurve extends ParametricCurve {
    protected mDegree: number;
    protected mNumControls: number;

    // mControls[order] holds the control points (order 0) and the
    // first-, second- and third-order differences of them.
    protected mControls: Vector[][];
    protected mChoose: Array2<number>;

    // Construction. The number of control points must be degree + 1. This
    // object copies the input array. The domain is t in [0,1].
    constructor(dimension: number, degree: number,
        controls: readonly Vector[]) {
        super(dimension, 0, 1);

        logAssert(degree >= 2 && controls !== undefined && controls !== null,
            'Invalid input.');

        this.mDegree = degree;
        this.mNumControls = degree + 1;
        this.mChoose = new Array2<number>(this.mNumControls, this.mNumControls);
        this.mChoose.fill(0);
        this.mControls = [[], [], [], []];

        // Copy the controls.
        for (let i = 0; i < this.mNumControls; ++i) {
            this.mControls[0].push(controls[i].clone());
        }

        // Compute first-order differences.
        for (let i = 0, ip1 = 1; ip1 < this.mNumControls; ++i, ++ip1) {
            this.mControls[1].push(
                subtract(this.mControls[0][ip1], this.mControls[0][i]));
        }

        // Compute second-order differences.
        for (let i = 0, ip2 = 2; ip2 < this.mNumControls; ++i, ++ip2) {
            this.mControls[2].push(
                subtract(this.mControls[1][i + 1], this.mControls[1][i]));
        }

        // Compute third-order differences.
        if (degree >= 3) {
            for (let i = 0, ip3 = 3; ip3 < this.mNumControls; ++i, ++ip3) {
                this.mControls[3].push(
                    subtract(this.mControls[2][i + 1], this.mControls[2][i]));
            }
        }

        // Compute combinatorial values Choose(n,k) and store in
        // mChoose[n][k]. The values mChoose[r][c] are invalid for r < c; that
        // is, we use only the entries for r >= c.
        this.mChoose.set(0, 0, 1);
        this.mChoose.set(0, 1, 1);
        this.mChoose.set(1, 1, 1);
        for (let i = 2; i <= this.mDegree; ++i) {
            this.mChoose.set(0, i, 1);
            this.mChoose.set(i, i, 1);
            for (let j = 1; j < i; ++j) {
                this.mChoose.set(j, i,
                    this.mChoose.get(j - 1, i - 1) + this.mChoose.get(j, i - 1));
            }
        }

        this.mConstructed = true;
    }

    // Member access.
    getDegree(): number {
        return this.mDegree;
    }

    getNumControls(): number {
        return this.mNumControls;
    }

    // The returned array aliases internal storage (upstream returns a
    // pointer), so writes through it modify the curve.
    getControls(): Vector[] {
        return this.mControls[0];
    }

    // Evaluation of the curve. The function supports derivative calculation
    // through order 3; that is, order <= 3 is required. If you want only the
    // position, pass in order of 0. If you want the position and first
    // derivative, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order (use createJet()).
    // The values are ordered as: position, first derivative, second
    // derivative, third derivative.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        const supOrder = ParametricCurve.SUP_ORDER;
        if (!this.mConstructed || order >= supOrder) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i < supOrder; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        // Compute position.
        const omt = 1 - t;
        jet[0] = this.compute(t, omt, 0);
        if (order >= 1) {
            // Compute first derivative.
            jet[1] = this.compute(t, omt, 1);
            if (order >= 2) {
                // Compute second derivative.
                jet[2] = this.compute(t, omt, 2);
                if (order >= 3) {
                    // Compute third derivative.
                    if (this.mDegree >= 3) {
                        jet[3] = this.compute(t, omt, 3);
                    }
                    else {
                        jet[3].makeZero();
                    }
                }
            }
        }
    }

    // Support for evaluate(...).
    protected compute(t: number, omt: number, order: number): Vector {
        const controls = this.mControls[order];
        const n = this.mDimension;
        const result = new Vector(n);
        for (let k = 0; k < n; ++k) {
            result.values[k] = omt * controls[0].values[k];
        }

        let tpow = t;
        const isup = this.mDegree - order;
        for (let i = 1; i < isup; ++i) {
            const c = this.mChoose.get(i, isup) * tpow;
            const control = controls[i].values;
            for (let k = 0; k < n; ++k) {
                result.values[k] = (result.values[k] + c * control[k]) * omt;
            }
            tpow *= t;
        }
        const last = controls[isup].values;
        for (let k = 0; k < n; ++k) {
            result.values[k] += tpow * last[k];
        }

        let multiplier = 1;
        for (let i = 0; i < order; ++i) {
            multiplier *= this.mDegree - i;
        }
        for (let k = 0; k < n; ++k) {
            result.values[k] *= multiplier;
        }

        return result;
    }
}

function subtract(v0: Vector, v1: Vector): Vector {
    const result = new Vector(v0.size);
    for (let k = 0; k < v0.size; ++k) {
        result.values[k] = v0.values[k] - v1.values[k];
    }
    return result;
}
