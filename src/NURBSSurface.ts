// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NURBSSurface.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A NURBS (nonuniform rational B-spline) surface
//   X(u,v) = [sum_{i0,i1} N_{i0}(u)*N_{i1}(v)*w[i0][i1]*C[i0][i1]] / w(u,v)
//   w(u,v) = sum_{i0,i1} N_{i0}(u)*N_{i1}(v)*w[i0][i1]
// where the N are the B-spline basis functions of the two BasisFunction
// objects, C are the control points and w are the (positive) weights. When
// all the weights are equal, the surface reduces to the B-spline surface with
// the same control net.
//
// Port notes (matching the BSplineSurface precedent): upstream
// 'template <int32_t N, typename Real>' becomes a runtime dimension carried
// by the control points (Vector objects), passed to the constructor. The two
// BasisFunctionInput arguments become an array of two inputs. The C++
// 'Vector<N, Real> const* controls' and 'Real const* weights' pointers (null
// to defer setting them) become optional arrays; when omitted the control
// points are zero vectors and the weights are zero, as upstream. The control
// points are copied (C++ value semantics), and getControl returns a reference
// into the internal storage as upstream does.

import { BasisFunction, BasisFunctionInput } from './BasisFunction.js';
import { ParametricSurface } from './ParametricSurface.js';
import { Vector } from './Vector.js';

export class NURBSSurface extends ParametricSurface {
    protected mBasisFunction: BasisFunction[];
    protected mNumControls: number[];
    protected mControls: Vector[];
    protected mWeights: number[];

    // Construction. If 'controls' is provided, a copy is made of the control
    // points; similarly for 'weights'. To defer setting the control points or
    // weights, omit the arguments and later access them via getControls(),
    // getWeights(), setControl() or setWeight(). The 'controls' and 'weights'
    // must be stored in row-major order, attribute[i0 + numControls0*i1]. As a
    // 2D array, this corresponds to attribute2D[i1][i0].
    constructor(dimension: number, input: BasisFunctionInput[],
        controls?: readonly Vector[], weights?: readonly number[]) {
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
        this.mWeights = new Array<number>(numControls);
        for (let i = 0; i < numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
            this.mWeights[i] = (weights !== undefined ? weights[i] : 0);
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

    // The returned arrays alias internal storage (upstream returns pointers),
    // so writes through them modify the surface.
    getControls(): Vector[] {
        return this.mControls;
    }

    getWeights(): number[] {
        return this.mWeights;
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

    setWeight(i0: number, i1: number, weight: number): void {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)) {
            this.mWeights[i0 + this.mNumControls[0] * i1] = weight;
        }
    }

    getWeight(i0: number, i1: number): number {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)) {
            return this.mWeights[i0 + this.mNumControls[0] * i1];
        }
        return this.mWeights[0];
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

        const ru = this.mBasisFunction[0].evaluate(u, order);
        const rv = this.mBasisFunction[1].evaluate(v, order);
        const iumin = ru.minIndex, iumax = ru.maxIndex;
        const ivmin = rv.minIndex, ivmax = rv.maxIndex;

        // Compute position.
        const p = this.compute(0, 0, iumin, iumax, ivmin, ivmax);
        const invW = 1 / p.w;
        jet[0] = scale(invW, p.x);

        if (order >= 1) {
            // Compute first-order derivatives.
            const du = this.compute(1, 0, iumin, iumax, ivmin, ivmax);
            jet[1] = scale(invW, combine(du.x, -du.w, jet[0]));

            const dv = this.compute(0, 1, iumin, iumax, ivmin, ivmax);
            jet[2] = scale(invW, combine(dv.x, -dv.w, jet[0]));

            if (order >= 2) {
                // Compute second-order derivatives.
                const duu = this.compute(2, 0, iumin, iumax, ivmin, ivmax);
                jet[3] = scale(invW, combine(combine(duu.x, -2 * du.w, jet[1]),
                    -duu.w, jet[0]));

                const duv = this.compute(1, 1, iumin, iumax, ivmin, ivmax);
                jet[4] = scale(invW, combine(combine(combine(duv.x, -du.w, jet[2]),
                    -dv.w, jet[1]), -duv.w, jet[0]));

                const dvv = this.compute(0, 2, iumin, iumax, ivmin, ivmax);
                jet[5] = scale(invW, combine(combine(dvv.x, -2 * dv.w, jet[2]),
                    -dvv.w, jet[0]));
            }
        }
    }

    // Support for evaluate(...). Upstream writes the numerator vector X and
    // the weight w through reference parameters; the port returns them.
    protected compute(uOrder: number, vOrder: number, iumin: number,
        iumax: number, ivmin: number, ivmax: number): { x: Vector, w: number } {
        // The j*-indices introduce a tiny amount of overhead in order to
        // handle both aperiodic and periodic splines. For aperiodic splines,
        // j* = i* always.
        const numControls0 = this.mNumControls[0];
        const numControls1 = this.mNumControls[1];
        const x = new Vector(this.mDimension);
        let w = 0;
        for (let iv = ivmin; iv <= ivmax; ++iv) {
            const tmpv = this.mBasisFunction[1].getValue(vOrder, iv);
            const jv = (iv >= numControls1 ? iv - numControls1 : iv);
            for (let iu = iumin; iu <= iumax; ++iu) {
                const tmpu = this.mBasisFunction[0].getValue(uOrder, iu);
                const ju = (iu >= numControls0 ? iu - numControls0 : iu);
                const index = ju + numControls0 * jv;
                const tmp = tmpu * tmpv * this.mWeights[index];
                const control = this.mControls[index];
                for (let d = 0; d < this.mDimension; ++d) {
                    x.values[d] += tmp * control.values[d];
                }
                w += tmp;
            }
        }
        return { x, w };
    }
}

// s * v, allocating a new vector.
function scale(s: number, v: Vector): Vector {
    const result = new Vector(v.size);
    for (let d = 0; d < result.values.length; ++d) {
        result.values[d] = s * v.values[d];
    }
    return result;
}

// a + s * b, accumulated in place into 'a' (which is always a temporary
// created by compute(...) or by a previous combine(...)).
function combine(a: Vector, s: number, b: Vector): Vector {
    for (let d = 0; d < a.values.length; ++d) {
        a.values[d] += s * b.values[d];
    }
    return a;
}
