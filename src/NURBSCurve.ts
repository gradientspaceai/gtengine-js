// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NURBSCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A NURBS curve X(t) = (sum_i w[i] * N_i(t) * C[i]) / (sum_i w[i] * N_i(t)),
// where the N_i are the B-spline basis functions of the BasisFunction object,
// the C[i] are the control points and the w[i] are the weights. The domain is
// t in [t[d], t[n]], where t[d] and t[n] are knots with d the degree and n the
// number of control points.
//
// Port notes (following BSplineCurve and NURBSSurface): upstream
// 'template <int32_t N, typename Real>' becomes a runtime dimension passed as
// the first constructor argument and carried by the control points (Vector
// objects). The C++ null pointers used to defer setting the control points and
// weights become optional arguments; when omitted, the control points are zero
// vectors and the weights are zero (as upstream).

import { BasisFunction, BasisFunctionInput } from './BasisFunction';
import { ParametricCurve } from './ParametricCurve';
import { Vector } from './Vector';

export class NURBSCurve extends ParametricCurve {
    protected mBasisFunction: BasisFunction;
    protected mControls: Vector[];
    protected mWeights: number[];

    // Construction. If 'controls' is provided, a copy is made of the control
    // points; similarly for 'weights'. To defer setting them, omit the
    // arguments and later use getControls()/getWeights()/setControl()/
    // setWeight().
    constructor(dimension: number, input: BasisFunctionInput,
        controls?: readonly Vector[], weights?: readonly number[]) {
        super(dimension, 0, 1);

        this.mBasisFunction = new BasisFunction();
        this.mBasisFunction.create(input);

        // The mBasisFunction stores the domain but so does ParametricCurve.
        this.mTime[0] = this.mBasisFunction.getMinDomain();
        this.mTime[this.mTime.length - 1] = this.mBasisFunction.getMaxDomain();

        // The replication of control points for periodic splines is avoided
        // by wrapping the i-loop index in evaluate.
        this.mControls = new Array<Vector>(input.numControls);
        this.mWeights = new Array<number>(input.numControls);
        for (let i = 0; i < input.numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
            this.mWeights[i] = (weights !== undefined ? weights[i] : 0);
        }

        this.mConstructed = true;
    }

    // Member access.
    getBasisFunction(): BasisFunction {
        return this.mBasisFunction;
    }

    getNumControls(): number {
        return this.mControls.length;
    }

    // The returned arrays alias internal storage (upstream returns pointers),
    // so writes through them modify the curve.
    getControls(): Vector[] {
        return this.mControls;
    }

    getWeights(): number[] {
        return this.mWeights;
    }

    setControl(i: number, control: Vector): void {
        if (0 <= i && i < this.getNumControls()) {
            this.mControls[i] = control.clone();
        }
    }

    getControl(i: number): Vector {
        if (0 <= i && i < this.getNumControls()) {
            return this.mControls[i];
        }
        // Invalid index, return something.
        return this.mControls[0];
    }

    setWeight(i: number, weight: number): void {
        if (0 <= i && i < this.getNumControls()) {
            this.mWeights[i] = weight;
        }
    }

    getWeight(i: number): number {
        if (0 <= i && i < this.getNumControls()) {
            return this.mWeights[i];
        }
        // Invalid index, return something.
        return this.mWeights[0];
    }

    // Evaluation of the curve. The function supports derivative calculation
    // through order 3; that is, order <= 3 is required. If you want only the
    // position, pass in order of 0. If you want the position and first
    // derivative, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order (use createJet()). The
    // values are ordered as: position, first derivative, second derivative,
    // third derivative.
    override evaluate(t: number, order: number, jet: Vector[]): void {
        const supOrder = ParametricCurve.SUP_ORDER;
        if (!this.mConstructed || order >= supOrder) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i < supOrder; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        const indices = this.mBasisFunction.evaluate(t, order);
        const imin = indices.minIndex;
        const imax = indices.maxIndex;
        const n = this.mDimension;

        // Compute position.
        const c0 = this.compute(0, imin, imax);
        const invW = 1 / c0.w;
        const jet0 = new Vector(n);
        for (let k = 0; k < n; ++k) {
            jet0.values[k] = invW * c0.x.values[k];
        }
        jet[0] = jet0;

        if (order >= 1) {
            // Compute first derivative.
            const c1 = this.compute(1, imin, imax);
            const jet1 = new Vector(n);
            for (let k = 0; k < n; ++k) {
                jet1.values[k] = invW
                    * (c1.x.values[k] - c1.w * jet0.values[k]);
            }
            jet[1] = jet1;

            if (order >= 2) {
                // Compute second derivative.
                const c2 = this.compute(2, imin, imax);
                const jet2 = new Vector(n);
                for (let k = 0; k < n; ++k) {
                    jet2.values[k] = invW * (c2.x.values[k]
                        - 2 * c1.w * jet1.values[k] - c2.w * jet0.values[k]);
                }
                jet[2] = jet2;

                if (order === 3) {
                    // Compute third derivative.
                    const c3 = this.compute(3, imin, imax);
                    const jet3 = new Vector(n);
                    for (let k = 0; k < n; ++k) {
                        jet3.values[k] = invW * (c3.x.values[k]
                            - 3 * c1.w * jet2.values[k]
                            - 3 * c2.w * jet1.values[k]
                            - c3.w * jet0.values[k]);
                    }
                    jet[3] = jet3;
                }
            }
        }
    }

    // Support for evaluate(...). Upstream writes the weighted sum X and the
    // weight sum w through reference parameters; the port returns them.
    protected compute(order: number, imin: number, imax: number):
        { x: Vector, w: number } {
        // The j-index introduces a tiny amount of overhead in order to handle
        // both aperiodic and periodic splines. For aperiodic splines, j = i
        // always.
        const numControls = this.getNumControls();
        const n = this.mDimension;
        const x = new Vector(n);
        let w = 0;
        for (let i = imin; i <= imax; ++i) {
            const j = (i >= numControls ? i - numControls : i);
            const tmp = this.mBasisFunction.getValue(order, i) * this.mWeights[j];
            const control = this.mControls[j].values;
            for (let k = 0; k < n; ++k) {
                x.values[k] += tmp * control[k];
            }
            w += tmp;
        }
        return { x, w };
    }
}
