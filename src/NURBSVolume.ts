// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NURBSVolume.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A NURBS (nonuniform rational B-spline) volume
//   X(u,v,w) = [sum_{i0,i1,i2} N_{i0}(u)*N_{i1}(v)*N_{i2}(w)*h*C] / h(u,v,w)
//   h(u,v,w) = sum_{i0,i1,i2} N_{i0}(u)*N_{i1}(v)*N_{i2}(w)*h[i0][i1][i2]
// where the N are the B-spline basis functions of the three BasisFunction
// objects, C are the control points and h are the (positive) weights. When
// all the weights are equal, the volume reduces to the B-spline volume with
// the same control net.
//
// Port notes (matching the BSplineVolume precedent): upstream
// 'template <int32_t N, typename Real>' becomes a runtime dimension carried
// by the control points (Vector objects), passed to the constructor. The
// three BasisFunctionInput arguments become an array of three inputs. The C++
// 'controls' and 'weights' pointers (null to defer setting them) become
// optional arrays; when omitted the control points are zero vectors and the
// weights are zero, as upstream. C++ 'operator bool' has no TS equivalent; it
// is ported as isConstructed().

import { BasisFunction, BasisFunctionInput } from './BasisFunction.js';
import { Vector } from './Vector.js';

export class NURBSVolume {
    // The number of entries a 'jet' array must have: position X; first-order
    // derivatives dX/du, dX/dv, dX/dw; second-order derivatives d2X/du2,
    // d2X/dv2, d2X/dw2, d2X/dudv, d2X/dudw, d2X/dvdw.
    static readonly SUP_ORDER = 10;

    private mDimension: number;
    private mBasisFunction: BasisFunction[];
    private mNumControls: number[];
    private mControls: Vector[];
    private mWeights: number[];
    private mConstructed: boolean;

    // Construction. If 'controls' is provided, a copy is made of the control
    // points; similarly for 'weights'. To defer setting the control points or
    // weights, omit the arguments and later access them via getControls(),
    // getWeights(), setControl() or setWeight(). The 'controls' and 'weights'
    // must be stored in lexicographical order,
    //   attribute[i0 + numControls0 * (i1 + numControls1 * i2)]
    // As a 3D array, this corresponds to attribute3D[i2][i1][i0].
    constructor(dimension: number, input: BasisFunctionInput[],
        controls?: readonly Vector[], weights?: readonly number[]) {
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
        this.mWeights = new Array<number>(numControls);
        for (let i = 0; i < numControls; ++i) {
            this.mControls[i] = (controls !== undefined
                ? controls[i].clone() : new Vector(dimension));
            this.mWeights[i] = (weights !== undefined ? weights[i] : 0);
        }

        this.mConstructed = true;
    }

    // To validate construction, create an object as shown:
    //     const volume = new NURBSVolume(parameters);
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

    // The returned arrays alias internal storage (upstream returns pointers),
    // so writes through them modify the volume.
    getControls(): Vector[] {
        return this.mControls;
    }

    getWeights(): number[] {
        return this.mWeights;
    }

    // Port addition: upstream NURBSVolume exposes only the raw controls and
    // weights pointers, with no SetControl/GetControl/SetWeight/GetWeight (an
    // inconsistency with NURBSSurface and with BSplineVolume, both of which
    // have the indexed accessors). The accessors below mirror BSplineVolume's
    // and use the same lexicographical index as evaluate(...), so callers need
    // not redo the index arithmetic.
    setControl(i0: number, i1: number, i2: number, control: Vector): void {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            this.mControls[this.getIndex(i0, i1, i2)] = control.clone();
        }
    }

    getControl(i0: number, i1: number, i2: number): Vector {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            return this.mControls[this.getIndex(i0, i1, i2)];
        }
        return this.mControls[0];
    }

    setWeight(i0: number, i1: number, i2: number, weight: number): void {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            this.mWeights[this.getIndex(i0, i1, i2)] = weight;
        }
    }

    getWeight(i0: number, i1: number, i2: number): number {
        if (0 <= i0 && i0 < this.getNumControls(0)
            && 0 <= i1 && i1 < this.getNumControls(1)
            && 0 <= i2 && i2 < this.getNumControls(2)) {
            return this.mWeights[this.getIndex(i0, i1, i2)];
        }
        return this.mWeights[0];
    }

    private getIndex(i0: number, i1: number, i2: number): number {
        return i0 + this.mNumControls[0] * (i1 + this.mNumControls[1] * i2);
    }

    // Allocate storage for a jet of this volume: SUP_ORDER zero vectors of the
    // volume dimension. This is the port of the caller-side
    // 'std::array<Vector<N, Real>, SUP_ORDER> jet{}'.
    createJet(): Vector[] {
        const jet = new Array<Vector>(NURBSVolume.SUP_ORDER);
        for (let i = 0; i < NURBSVolume.SUP_ORDER; ++i) {
            jet[i] = new Vector(this.mDimension);
        }
        return jet;
    }

    // Evaluation of the volume. The function supports derivative calculation
    // through order 2; that is, order <= 2 is required. If you want only the
    // position, pass in order of 0. If you want the position and first-order
    // derivatives, pass in order of 1, and so on. The output array 'jet' must
    // have enough storage to support the maximum order (use createJet()). The
    // values are ordered as: position X; first-order derivatives dX/du, dX/dv,
    // dX/dw; second-order derivatives d2X/du2, d2X/dv2, d2X/dw2, d2X/dudv,
    // d2X/dudw, d2X/dvdw.
    evaluate(u: number, v: number, w: number, order: number, jet: Vector[]): void {
        if (!this.mConstructed || order >= NURBSVolume.SUP_ORDER) {
            // Return a zero-valued jet for invalid state.
            for (let i = 0; i < NURBSVolume.SUP_ORDER; ++i) {
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
        const p = this.compute(0, 0, 0, a, b, c, d, e, f);
        const invH = 1 / p.h;
        jet[0] = scale(invH, p.x);

        if (order >= 1) {
            // Compute first-order derivatives.
            const du = this.compute(1, 0, 0, a, b, c, d, e, f);
            jet[1] = scale(invH, combine(du.x, -du.h, jet[0]));

            const dv = this.compute(0, 1, 0, a, b, c, d, e, f);
            jet[2] = scale(invH, combine(dv.x, -dv.h, jet[0]));

            const dw = this.compute(0, 0, 1, a, b, c, d, e, f);
            jet[3] = scale(invH, combine(dw.x, -dw.h, jet[0]));

            if (order >= 2) {
                // Compute second-order derivatives.
                const duu = this.compute(2, 0, 0, a, b, c, d, e, f);
                jet[4] = scale(invH, combine(combine(duu.x, -2 * du.h, jet[1]),
                    -duu.h, jet[0]));

                const dvv = this.compute(0, 2, 0, a, b, c, d, e, f);
                jet[5] = scale(invH, combine(combine(dvv.x, -2 * dv.h, jet[2]),
                    -dvv.h, jet[0]));

                const dww = this.compute(0, 0, 2, a, b, c, d, e, f);
                jet[6] = scale(invH, combine(combine(dww.x, -2 * dw.h, jet[3]),
                    -dww.h, jet[0]));

                const duv = this.compute(1, 1, 0, a, b, c, d, e, f);
                jet[7] = scale(invH, combine(combine(combine(duv.x, -du.h, jet[2]),
                    -dv.h, jet[1]), -duv.h, jet[0]));

                const duw = this.compute(1, 0, 1, a, b, c, d, e, f);
                jet[8] = scale(invH, combine(combine(combine(duw.x, -du.h, jet[3]),
                    -dw.h, jet[1]), -duw.h, jet[0]));

                const dvw = this.compute(0, 1, 1, a, b, c, d, e, f);
                jet[9] = scale(invH, combine(combine(combine(dvw.x, -dv.h, jet[3]),
                    -dw.h, jet[2]), -dvw.h, jet[0]));
            }
        }
    }

    // Support for evaluate(...). Upstream writes the numerator vector X and
    // the weight h through reference parameters; the port returns them.
    private compute(uOrder: number, vOrder: number, wOrder: number,
        iumin: number, iumax: number, ivmin: number, ivmax: number,
        iwmin: number, iwmax: number): { x: Vector, h: number } {
        // The j*-indices introduce a tiny amount of overhead in order to
        // handle both aperiodic and periodic splines. For aperiodic splines,
        // j* = i* always.
        const numControls0 = this.mNumControls[0];
        const numControls1 = this.mNumControls[1];
        const numControls2 = this.mNumControls[2];
        const x = new Vector(this.mDimension);
        let h = 0;
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
                    const index = ju + numControls0 * (jv + numControls1 * jw);
                    const tmp = (tmpu * tmpvw) * this.mWeights[index];
                    const control = this.mControls[index];
                    for (let k = 0; k < this.mDimension; ++k) {
                        x.values[k] += tmp * control.values[k];
                    }
                    h += tmp;
                }
            }
        }
        return { x, h };
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
