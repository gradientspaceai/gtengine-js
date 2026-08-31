// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineVolume.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A B-spline volume
//   X(u,v,w) = sum_{i0,i1,i2} N_{i0}(u)*N_{i1}(v)*N_{i2}(w)*C[i0][i1][i2],
// where the N are the B-spline basis functions of the three BasisFunction
// objects and C are the control points.
//
// Port notes: upstream 'template <int32_t N, typename Real>' becomes a
// runtime dimension carried by the control points (Vector objects), passed to
// the constructor. The C++ 'Vector<N, Real> const* controls' pointer (null to
// defer setting the control points) becomes an optional array. C++
// 'operator bool' has no TS equivalent; it is ported as isConstructed(),
// matching the ParametricSurface precedent.

import { BasisFunction, BasisFunctionInput } from './BasisFunction';
import { Vector } from './Vector';

export class BSplineVolume {
    // The number of entries a 'jet' array must have: position X; first-order
    // derivatives dX/du, dX/dv, dX/dw; second-order derivatives d2X/du2,
    // d2X/dv2, d2X/dw2, d2X/dudv, d2X/dudw, d2X/dvdw.
    static readonly SUP_ORDER = 10;

    private mDimension: number;
    private mBasisFunction: BasisFunction[];
    private mNumControls: number[];
    private mControls: Vector[];
    private mConstructed: boolean;

    // Construction. If 'controls' is provided, a copy is made of the control
    // points. To defer setting the control points, omit the argument and
    // later access them via getControls() or setControl(). The input
    // 'controls' must be stored in lexicographical order,
    // control[i0+numControls0*(i1+numControls1*i2)]. As a 3D array, this
    // corresponds to control3D[i2][i1][i0].
    constructor(dimension: number, input: BasisFunctionInput[],
        controls?: readonly Vector[]) {
        this.mDimension = dimension;
        this.mConstructed = false;

        this.mBasisFunction = [new BasisFunction(), new BasisFunction(), new BasisFunction()];
        this.mNumControls = [0, 0, 0];
        for (let i = 0; i < 3; ++i) {
            this.mNumControls[i] = input[i].numControls;
            this.mBasisFunction[i].create(input[i]);
        }

        // The replication of control points for periodic splines is avoided
        // by wrapping the i-loop index in evaluate.
        const numControls = this.mNumControls[0] * this.mNumControls[1] * this.mNumControls[2];
        this.mControls = new Array<Vector>(numControls);
        for (let i = 0; i < numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
        }

        this.mConstructed = true;
    }

    // To validate construction, create an object as shown:
    //     const volume = new BSplineVolume(parameters);
    //     if (!volume.isConstructed()) { /* constructor failed */ }
    isConstructed(): boolean {
        return this.mConstructed;
    }

    getDimension(): number {
        return this.mDimension;
    }

    // Member access. The index 'dim' must be in {0,1,2}.
    getBasisFunction(dim: number): BasisFunction {
        return this.mBasisFunction[dim];
    }

    getMinDomain(dim: number): number {
        return this.mBasisFunction[dim].getMinDomain();
    }

    getMaxDomain(dim: number): number {
        return this.mBasisFunction[dim].getMaxDomain();
    }

    getNumControls(dim: number): number {
        return this.mNumControls[dim];
    }

    // The returned array aliases internal storage (upstream returns a
    // pointer), so writes through it modify the volume.
    getControls(): Vector[] {
        return this.mControls;
    }

    setControl(i0: number, i1: number, i2: number, control: Vector): void {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            this.mControls[i0 + this.mNumControls[0] * (i1 + this.mNumControls[1] * i2)]
                = control.clone();
        }
    }

    getControl(i0: number, i1: number, i2: number): Vector {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            return this.mControls[i0 + this.mNumControls[0] * (i1 + this.mNumControls[1] * i2)];
        }
        return this.mControls[0];
    }

    // Allocate storage for a jet of this volume: SUP_ORDER zero vectors of
    // the volume dimension. This is the port of the caller-side
    // 'std::array<Vector<N, Real>, SUP_ORDER> jet{}'.
    createJet(): Vector[] {
        const jet = new Array<Vector>(BSplineVolume.SUP_ORDER);
        for (let i = 0; i < BSplineVolume.SUP_ORDER; ++i) {
            jet[i] = new Vector(this.mDimension);
        }
        return jet;
    }

    // Evaluation of the volume. The function supports derivative calculation
    // through order 2; that is, order <= 2 is required. If you want only the
    // position, pass in order of 0. If you want the position and first-order
    // derivatives, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order (use createJet()). The
    // values are ordered as: position X; first-order derivatives dX/du,
    // dX/dv, dX/dw; second-order derivatives d2X/du2, d2X/dv2, d2X/dw2,
    // d2X/dudv, d2X/dudw, d2X/dvdw.
    evaluate(u: number, v: number, w: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed || order >= BSplineVolume.SUP_ORDER) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i < BSplineVolume.SUP_ORDER; ++i) {
                jet[i].makeZero();
            }
            return;
        }

        const ru = this.mBasisFunction[0].evaluate(u, order);
        const rv = this.mBasisFunction[1].evaluate(v, order);
        const rw = this.mBasisFunction[2].evaluate(w, order);
        const a = ru.minIndex, b = ru.maxIndex;
        const c = rv.minIndex, d = rv.maxIndex;
        const e = rw.minIndex, f = rw.maxIndex;

        // Compute position.
        jet[0] = this.compute(0, 0, 0, a, b, c, d, e, f);
        if (order >= 1) {
            // Compute first-order derivatives.
            jet[1] = this.compute(1, 0, 0, a, b, c, d, e, f);
            jet[2] = this.compute(0, 1, 0, a, b, c, d, e, f);
            jet[3] = this.compute(0, 0, 1, a, b, c, d, e, f);
            if (order >= 2) {
                // Compute second-order derivatives.
                jet[4] = this.compute(2, 0, 0, a, b, c, d, e, f);
                jet[5] = this.compute(0, 2, 0, a, b, c, d, e, f);
                jet[6] = this.compute(0, 0, 2, a, b, c, d, e, f);
                jet[7] = this.compute(1, 1, 0, a, b, c, d, e, f);
                jet[8] = this.compute(1, 0, 1, a, b, c, d, e, f);
                jet[9] = this.compute(0, 1, 1, a, b, c, d, e, f);
            }
        }
    }

    // Support for evaluate(...).
    private compute(uOrder: number, vOrder: number, wOrder: number,
        iumin: number, iumax: number, ivmin: number, ivmax: number,
        iwmin: number, iwmax: number): Vector {
        // The j*-indices introduce a tiny amount of overhead in order to
        // handle both aperiodic and periodic splines. For aperiodic splines,
        // j* = i* always.
        const numControls0 = this.mNumControls[0];
        const numControls1 = this.mNumControls[1];
        const numControls2 = this.mNumControls[2];
        const result = new Vector(this.mDimension);
        for (let iw = iwmin; iw <= iwmax; ++iw) {
            const tmpw = this.mBasisFunction[2].getValue(wOrder, iw);
            const jw = (iw >= numControls2 ? iw - numControls2 : iw);
            for (let iv = ivmin; iv <= ivmax; ++iv) {
                const tmpv = this.mBasisFunction[1].getValue(vOrder, iv);
                const tmpvw = tmpv * tmpw;
                const jv = (iv >= numControls1 ? iv - numControls1 : iv);
                for (let iu = iumin; iu <= iumax; ++iu) {
                    const tmpu = this.mBasisFunction[0].getValue(uOrder, iu);
                    const ju = (iu >= numControls0 ? iu - numControls0 : iu);
                    const control = this.mControls[ju + numControls0 * (jv + numControls1 * jw)];
                    const scalar = tmpu * tmpvw;
                    for (let k = 0; k < this.mDimension; ++k) {
                        result.values[k] += scalar * control.values[k];
                    }
                }
            }
        }
        return result;
    }
}
