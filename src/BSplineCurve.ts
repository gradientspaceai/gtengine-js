// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineCurve.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A B-spline curve X(t) = sum_i N_i(t) * C[i], where the N_i are the B-spline
// basis functions of the BasisFunction object and the C[i] are the control
// points. The domain is t in [t[d], t[n]], where t[d] and t[n] are knots with
// d the degree and n the number of control points.
//
// Port notes (following the BSplineSurface port): upstream
// 'template <int32_t N, typename Real>' becomes a runtime dimension passed to
// the constructor and carried by the control points (Vector objects). The
// C++ 'Vector<N, Real> const* controls' pointer (null to defer setting the
// control points) becomes an optional array; when it is omitted the control
// points are zero vectors of the given dimension. The control points are
// copied (C++ value semantics), and getControl returns a reference into the
// internal storage as upstream does.

import { BasisFunction, BasisFunctionInput } from './BasisFunction.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector } from './Vector.js';

export class BSplineCurve extends ParametricCurve {
    private mBasisFunction: BasisFunction;
    private mControls: Vector[];

    // Construction. If 'controls' is provided, a copy is made of the control
    // points. To defer setting the control points, omit the argument and
    // later access them via getControls() or setControl().
    constructor(dimension: number, input: BasisFunctionInput,
        controls?: readonly Vector[]) {
        super(dimension, 0, 1);

        this.mBasisFunction = new BasisFunction();
        this.mBasisFunction.create(input);

        // The mBasisFunction stores the domain but so does ParametricCurve.
        this.mTime[0] = this.mBasisFunction.getMinDomain();
        this.mTime[this.mTime.length - 1] = this.mBasisFunction.getMaxDomain();

        // The replication of control points for periodic splines is avoided
        // by wrapping the i-loop index in evaluate.
        this.mControls = new Array<Vector>(input.numControls);
        for (let i = 0; i < input.numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
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

    // The returned array aliases internal storage (upstream returns a
    // pointer), so writes through it modify the curve.
    getControls(): Vector[] {
        return this.mControls;
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
        return this.mControls[0];
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

        // Compute position.
        jet[0] = this.compute(0, indices.minIndex, indices.maxIndex);
        if (order >= 1) {
            // Compute first derivative.
            jet[1] = this.compute(1, indices.minIndex, indices.maxIndex);
            if (order >= 2) {
                // Compute second derivative.
                jet[2] = this.compute(2, indices.minIndex, indices.maxIndex);
                if (order === 3) {
                    jet[3] = this.compute(3, indices.minIndex, indices.maxIndex);
                }
            }
        }
    }

    // Support for evaluate(...).
    private compute(order: number, imin: number, imax: number): Vector {
        // The j-index introduces a tiny amount of overhead in order to handle
        // both aperiodic and periodic splines. For aperiodic splines, j = i
        // always.
        const numControls = this.getNumControls();
        const result = new Vector(this.mDimension);
        for (let i = imin; i <= imax; ++i) {
            const tmp = this.mBasisFunction.getValue(order, i);
            const j = (i >= numControls ? i - numControls : i);
            const control = this.mControls[j];
            for (let k = 0; k < this.mDimension; ++k) {
                result.values[k] += tmp * control.values[k];
            }
        }
        return result;
    }
}
