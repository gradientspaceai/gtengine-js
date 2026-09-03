// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineSurface.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A B-spline surface X(u,v) = sum_{i0,i1} N_{i0}(u) * N_{i1}(v) * C[i0][i1],
// where the N are the B-spline basis functions of the two BasisFunction
// objects and C are the control points.
//
// Port notes: upstream 'template <int32_t N, typename Real>' becomes a
// runtime dimension carried by the control points (Vector objects), passed to
// the constructor. The C++ 'Vector<N, Real> const* controls' pointer (null to
// defer setting the control points) becomes an optional array; when it is
// omitted the control points are zero vectors of the given dimension. The
// control points are copied (C++ value semantics), and getControl returns a
// reference into the internal storage as upstream does.

import { BasisFunction, BasisFunctionInput } from './BasisFunction.js';
import { ParametricSurface } from './ParametricSurface.js';
import { Vector } from './Vector.js';

export class BSplineSurface extends ParametricSurface {
    private mBasisFunction: BasisFunction[];
    private mNumControls: number[];
    private mControls: Vector[];

    // Construction. If 'controls' is provided, a copy is made of the control
    // points. To defer setting the control points, omit the argument and
    // later access them via getControls() or setControl(). The input
    // 'controls' must be stored in row-major order,
    // control[i0 + numControls0*i1]. As a 2D array, this corresponds to
    // control2D[i1][i0].
    constructor(dimension: number, input: BasisFunctionInput[],
        controls?: readonly Vector[]) {
        super(dimension, 0, 1, 0, 1, true);

        this.mBasisFunction = [new BasisFunction(), new BasisFunction()];
        this.mNumControls = [0, 0];
        for (let i = 0; i < 2; ++i) {
            this.mNumControls[i] = input[i].numControls;
            this.mBasisFunction[i].create(input[i]);
        }

        // The mBasisFunction stores the domain but so does ParametricSurface.
        this.mUMin = this.mBasisFunction[0].getMinDomain();
        this.mUMax = this.mBasisFunction[0].getMaxDomain();
        this.mVMin = this.mBasisFunction[1].getMinDomain();
        this.mVMax = this.mBasisFunction[1].getMaxDomain();

        // The replication of control points for periodic splines is avoided
        // by wrapping the i-loop index in evaluate.
        const numControls = this.mNumControls[0] * this.mNumControls[1];
        this.mControls = new Array<Vector>(numControls);
        for (let i = 0; i < numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
        }

        this.mConstructed = true;
    }

    // Member access. The index 'dim' must be in {0,1}.
    getBasisFunction(dim: number): BasisFunction {
        return this.mBasisFunction[dim];
    }

    getNumControls(dim: number): number {
        return this.mNumControls[dim];
    }

    // The returned array aliases internal storage (upstream returns a
    // pointer), so writes through it modify the surface.
    getControls(): Vector[] {
        return this.mControls;
    }

    setControl(i0: number, i1: number, control: Vector): void {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)) {
            this.mControls[i0 + this.mNumControls[0] * i1] = control.clone();
        }
    }

    getControl(i0: number, i1: number): Vector {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)) {
            return this.mControls[i0 + this.mNumControls[0] * i1];
        }
        return this.mControls[0];
    }

    // Evaluation of the surface. The function supports derivative calculation
    // through order 2; that is, order <= 2 is required. If you want only the
    // position, pass in order of 0. If you want the position and first-order
    // derivatives, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order (use createJet()). The
    // values are ordered as: position X; first-order derivatives dX/du,
    // dX/dv; second-order derivatives d2X/du2, d2X/dudv, d2X/dv2.
    evaluate(u: number, v: number, order: number, jet: Vector[]): void {
        const supOrder = ParametricSurface.SUP_ORDER;
        if (!this.mConstructed || order >= supOrder) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i < supOrder; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        const u0 = this.mBasisFunction[0].evaluate(u, order);
        const u1 = this.mBasisFunction[1].evaluate(v, order);

        // Compute position.
        jet[0] = this.compute(0, 0, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
        if (order >= 1) {
            // Compute first-order derivatives.
            jet[1] = this.compute(1, 0, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
            jet[2] = this.compute(0, 1, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
            if (order >= 2) {
                // Compute second-order derivatives.
                jet[3] = this.compute(2, 0, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
                jet[4] = this.compute(1, 1, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
                jet[5] = this.compute(0, 2, u0.minIndex, u0.maxIndex, u1.minIndex, u1.maxIndex);
            }
        }
    }

    // Support for evaluate(...).
    private compute(uOrder: number, vOrder: number, iumin: number, iumax: number,
        ivmin: number, ivmax: number): Vector {
        // The j*-indices introduce a tiny amount of overhead in order to
        // handle both aperiodic and periodic splines. For aperiodic splines,
        // j* = i* always.
        const numControls0 = this.mNumControls[0];
        const numControls1 = this.mNumControls[1];
        const result = new Vector(this.mDimension);
        for (let iv = ivmin; iv <= ivmax; ++iv) {
            const tmpv = this.mBasisFunction[1].getValue(vOrder, iv);
            const jv = (iv >= numControls1 ? iv - numControls1 : iv);
            for (let iu = iumin; iu <= iumax; ++iu) {
                const tmpu = this.mBasisFunction[0].getValue(uOrder, iu);
                const ju = (iu >= numControls0 ? iu - numControls0 : iu);
                const control = this.mControls[ju + numControls0 * jv];
                const scalar = tmpu * tmpv;
                for (let d = 0; d < this.mDimension; ++d) {
                    result.values[d] += scalar * control.values[d];
                }
            }
        }
        return result;
    }
}
