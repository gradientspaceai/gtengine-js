// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NURBSSphere.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for representing a circle as a NURBS curve or a sphere as a
// NURBS surface is described in
//   https://www.geometrictools.com/Documentation/NURBSCircleSphere.pdf
// The implementations are related to the documents as shown next.
//   NURBSEighthSphereDegree4 implements Section 3.1.2 (triangular domain)
//   NURBSHalfSphereDegree3 implements Section 3.2 (rectangular domain)
//   NURBSFullSphereDegree3 implements Section 2.3 (rectangular domain)
// TODO (upstream): The class NURBSSurface currently assumes a rectangular
// domain. Once support is added for triangular domains, make that new class a
// base class of the sphere-representing NURBS. This will allow sharing of the
// NURBS basis functions and evaluation framework.
//
// Port notes: the header holds three class templates; all three live in this
// file per the one-file-per-header rule. NURBSEighthSphereDegree4 has a
// triangular domain and is not a NURBSSurface; its 'Evaluate(u, v, maxOrder,
// Vector<3,Real> values[6])' keeps the caller-supplied output array, which
// createEighthSphereValues() allocates. Its 5-by-5 C++ arrays keep upstream's
// [j1][j0] index order (nested arrays) so the formulas read as upstream; only
// entries with j0 + j1 <= 4 are used. The two NURBSSurface subclasses set the
// protected mControls/mWeights arrays that the port's NURBSSurface exposes,
// with the same row-major layout attribute[i0 + numControls0*i1].

import { BasisFunctionInput, UniqueKnot } from './BasisFunction';
import { NURBSSurface } from './NURBSSurface';
import { Vector } from './Vector';

// Allocate the 'values' output of NURBSEighthSphereDegree4.evaluate: six
// zero-valued 3D vectors (position; dX/du, dX/dv; d2X/du2, d2X/dudv,
// d2X/dv2).
export function createEighthSphereValues(): Vector[] {
    const values = new Array<Vector>(6);
    for (let i = 0; i < 6; ++i) {
        values[i] = new Vector(3);
    }
    return values;
}

function makeTable(): number[][] {
    const table = new Array<number[]>(5);
    for (let i = 0; i < 5; ++i) {
        table[i] = new Array<number>(5).fill(0);
    }
    return table;
}

export class NURBSEighthSphereDegree4 {
    // For simplicity of the implementation, 2-dimensional arrays of size
    // 5-by-5 are used. Only array[r][c] is used where 0 <= r <= 4 and
    // 0 <= c <= 4 - r.
    private mControls: Vector[][];
    private mWeights: number[][];

    // The eighth sphere is x^2 + y^2 + z^2 = 1 for x >= 0, y >= 0, z >= 0.
    constructor() {
        const sqrt2 = Math.sqrt(2);
        const sqrt3 = Math.sqrt(3);
        const a0 = (sqrt3 - 1) / sqrt3;
        const a1 = (sqrt3 + 1) / (2 * sqrt3);
        const a2 = 1 - (5 - sqrt2) * (7 - sqrt3) / 46;
        const b0 = 4 * sqrt3 * (sqrt3 - 1);
        const b1 = 3 * sqrt2;
        const b2 = 4;
        const b3 = sqrt2 * (3 + 2 * sqrt2 - sqrt3) / sqrt3;

        this.mControls = new Array<Vector[]>(5);
        for (let i = 0; i < 5; ++i) {
            this.mControls[i] = new Array<Vector>(5);
            for (let j = 0; j < 5; ++j) {
                this.mControls[i][j] = new Vector(3);
            }
        }

        const c = (j1: number, j0: number, x: number, y: number,
            z: number): void => {
            this.mControls[j1][j0] = Vector.fromArray([x, y, z]);
        };

        c(0, 0, 0, 0, 1);       // P004
        c(0, 1, 0, a0, 1);      // P013
        c(0, 2, 0, a1, a1);     // P022
        c(0, 3, 0, 1, a0);      // P031
        c(0, 4, 0, 1, 0);       // P040

        c(1, 0, a0, 0, 1);      // P103
        c(1, 1, a2, a2, 1);     // P112
        c(1, 2, a2, 1, a2);     // P121
        c(1, 3, a0, 1, 0);      // P130

        c(2, 0, a1, 0, a1);     // P202
        c(2, 1, 1, a2, a2);     // P211
        c(2, 2, a1, a1, 0);     // P220

        c(3, 0, 1, 0, a0);      // P301
        c(3, 1, 1, a0, 0);      // P310

        c(4, 0, 1, 0, 0);       // P400

        this.mWeights = makeTable();
        this.mWeights[0][0] = b0;   // w004
        this.mWeights[0][1] = b1;   // w013
        this.mWeights[0][2] = b2;   // w022
        this.mWeights[0][3] = b1;   // w031
        this.mWeights[0][4] = b0;   // w040

        this.mWeights[1][0] = b1;   // w103
        this.mWeights[1][1] = b3;   // w112
        this.mWeights[1][2] = b3;   // w121
        this.mWeights[1][3] = b1;   // w130

        this.mWeights[2][0] = b2;   // w202
        this.mWeights[2][1] = b3;   // w211
        this.mWeights[2][2] = b2;   // w220

        this.mWeights[3][0] = b1;   // w301
        this.mWeights[3][1] = b1;   // w310

        this.mWeights[4][0] = b0;   // w400
    }

    // Member access. The returned arrays alias internal storage.
    getControls(): readonly Vector[][] {
        return this.mControls;
    }

    getWeights(): readonly number[][] {
        return this.mWeights;
    }

    // Evaluation of the surface. The function supports derivative calculation
    // through order 2; that is, maxOrder <= 2 is required. If you want only
    // the position, pass in maxOrder of 0. If you want the position and
    // first-order derivatives, pass in maxOrder of 1, and so on. The output
    // 'values' are ordered as: position X; first-order derivatives dX/du,
    // dX/dv; second-order derivatives d2X/du2, d2X/dudv, d2X/dv2. Allocate
    // 'values' with createEighthSphereValues().
    evaluate(u: number, v: number, maxOrder: number, values: Vector[]): void {
        const w = 1 - u - v;
        const uu = u * u, uv = u * v, uw = u * w;
        const vv = v * v, vw = v * w, ww = w * w;

        // Compute the order-0 polynomials. Only the elements to be used are
        // filled in. The other terms are zero but never used.
        const B = makeTable();
        B[0][0] = ww * ww;
        B[0][1] = 4 * vw * ww;
        B[0][2] = 6 * vv * ww;
        B[0][3] = 4 * vv * vw;
        B[0][4] = vv * vv;
        B[1][0] = 4 * uw * ww;
        B[1][1] = 12 * uv * ww;
        B[1][2] = 12 * uv * vw;
        B[1][3] = 4 * uv * vv;
        B[2][0] = 6 * uu * ww;
        B[2][1] = 12 * uu * vw;
        B[2][2] = 6 * uu * vv;
        B[3][0] = 4 * uu * uw;
        B[3][1] = 4 * uu * uv;
        B[4][0] = uu * uu;

        // Compute the NURBS position.
        const N = new Vector(3);
        let D = 0;
        for (let j1 = 0; j1 <= 4; ++j1) {
            for (let j0 = 0; j0 <= 4 - j1; ++j0) {
                const product = this.mWeights[j1][j0] * B[j1][j0];
                const control = this.mControls[j1][j0].values;
                for (let k = 0; k < 3; ++k) {
                    N.values[k] += product * control[k];
                }
                D += product;
            }
        }
        values[0] = divideVec(N, D);

        if (maxOrder >= 1) {
            // Compute the order-1 polynomials. Only the elements to be used
            // are filled in. The other terms are zero but never used.
            const WmU = w - u;
            const WmTwoU = WmU - u;
            const WmThreeU = WmTwoU - u;
            const TwoWmU = w + WmU;
            const ThreeWmU = w + TwoWmU;
            const WmV = w - v;
            const WmTwoV = WmV - v;
            const WmThreeV = WmTwoV - v;
            const TwoWmV = w + WmV;
            const ThreeWmV = w + TwoWmV;
            const Dsqr = D * D;

            const Bu = makeTable();
            Bu[0][0] = -4 * ww * w;
            Bu[0][1] = -12 * v * ww;
            Bu[0][2] = -12 * vv * w;
            Bu[0][3] = -4 * v * vv;
            Bu[0][4] = 0;
            Bu[1][0] = 4 * ww * WmThreeU;
            Bu[1][1] = 12 * vw * WmTwoU;
            Bu[1][2] = 12 * vv * WmU;
            Bu[1][3] = 4 * v * vv;
            Bu[2][0] = 12 * uw * WmU;
            Bu[2][1] = 12 * uv * TwoWmU;
            Bu[2][2] = 12 * u * vv;
            Bu[3][0] = 4 * uu * ThreeWmU;
            Bu[3][1] = 12 * uu * v;
            Bu[4][0] = 4 * uu * u;

            const Bv = makeTable();
            Bv[0][0] = -4 * ww * w;
            Bv[0][1] = 4 * ww * WmThreeV;
            Bv[0][2] = 12 * vw * WmV;
            Bv[0][3] = 4 * vv * ThreeWmV;
            Bv[0][4] = 4 * vv * v;
            Bv[1][0] = -12 * u * ww;
            Bv[1][1] = 12 * uw * WmTwoV;
            Bv[1][2] = 12 * uv * TwoWmV;
            Bv[1][3] = 12 * u * vv;
            Bv[2][0] = -12 * uu * w;
            Bv[2][1] = 12 * uu * WmV;
            Bv[2][2] = 12 * uu * v;
            Bv[3][0] = -4 * uu * u;
            Bv[3][1] = 4 * uu * u;
            Bv[4][0] = 0;

            const Nu = new Vector(3);
            const Nv = new Vector(3);
            let Du = 0, Dv = 0;
            for (let j1 = 0; j1 <= 4; ++j1) {
                for (let j0 = 0; j0 <= 4 - j1; ++j0) {
                    const weight = this.mWeights[j1][j0];
                    const control = this.mControls[j1][j0].values;
                    let product = weight * Bu[j1][j0];
                    for (let k = 0; k < 3; ++k) {
                        Nu.values[k] += product * control[k];
                    }
                    Du += product;
                    product = weight * Bv[j1][j0];
                    for (let k = 0; k < 3; ++k) {
                        Nv.values[k] += product * control[k];
                    }
                    Dv += product;
                }
            }
            const numerDU = combine(D, Nu, -Du, N);
            const numerDV = combine(D, Nv, -Dv, N);
            values[1] = divideVec(numerDU, Dsqr);
            values[2] = divideVec(numerDV, Dsqr);

            if (maxOrder >= 2) {
                // Compute the order-2 polynomials. Only the elements to be
                // used are filled in. The other terms are zero but never
                // used.
                const Dcub = Dsqr * D;

                const Buu = makeTable();
                Buu[0][0] = 12 * ww;
                Buu[0][1] = 24 * vw;
                Buu[0][2] = 12 * vv;
                Buu[0][3] = 0;
                Buu[0][4] = 0;
                Buu[1][0] = -24 * w * WmU;
                Buu[1][1] = -24 * v * TwoWmU;
                Buu[1][2] = -24 * vv;
                Buu[1][3] = 0;
                Buu[2][0] = 12 * (ww - 4 * uw + uu);
                Buu[2][1] = 24 * v * WmTwoU;
                Buu[2][2] = 12 * vv;
                Buu[3][0] = 24 * u * WmU;
                Buu[3][1] = 24 * uv;
                Buu[4][0] = 12 * uu;

                const Buv = makeTable();
                Buv[0][0] = 12 * ww;
                Buv[0][1] = -12 * w * WmTwoV;
                Buv[0][2] = -12 * v * TwoWmV;
                Buv[0][3] = -12 * vv;
                Buv[0][4] = 0;
                Buv[1][0] = -12 * w * WmTwoU;
                Buv[1][1] = 12 * (ww + 2 * (uv - uw - vw));
                Buv[1][2] = 12 * v * (2 * WmU - v);
                Buv[1][3] = 12 * vv;
                Buv[2][0] = -12 * u * TwoWmU;
                Buv[2][1] = 12 * u * (2 * WmV - u);
                Buv[2][2] = 24 * uv;
                Buv[3][0] = -12 * uu;
                Buv[3][1] = 12 * uu;
                Buv[4][0] = 0;

                const Bvv = makeTable();
                Bvv[0][0] = 12 * ww;
                Bvv[0][1] = -24 * w * WmV;
                Bvv[0][2] = 12 * (ww - 4 * vw + vv);
                Bvv[0][3] = 24 * v * WmV;
                Bvv[0][4] = 12 * vv;
                Bvv[1][0] = 24 * uw;
                Bvv[1][1] = -24 * u * TwoWmV;
                Bvv[1][2] = 24 * u * WmTwoV;
                Bvv[1][3] = 24 * uv;
                Bvv[2][0] = 12 * uu;
                Bvv[2][1] = -24 * uu;
                Bvv[2][2] = 12 * uu;
                Bvv[3][0] = 0;
                Bvv[3][1] = 0;
                Bvv[4][0] = 0;

                const Nuu = new Vector(3);
                const Nuv = new Vector(3);
                const Nvv = new Vector(3);
                let Duu = 0, Duv = 0, Dvv = 0;
                for (let j1 = 0; j1 <= 4; ++j1) {
                    for (let j0 = 0; j0 <= 4 - j1; ++j0) {
                        const weight = this.mWeights[j1][j0];
                        const control = this.mControls[j1][j0].values;
                        let product = weight * Buu[j1][j0];
                        for (let k = 0; k < 3; ++k) {
                            Nuu.values[k] += product * control[k];
                        }
                        Duu += product;
                        product = weight * Buv[j1][j0];
                        for (let k = 0; k < 3; ++k) {
                            Nuv.values[k] += product * control[k];
                        }
                        Duv += product;
                        product = weight * Bvv[j1][j0];
                        for (let k = 0; k < 3; ++k) {
                            Nvv.values[k] += product * control[k];
                        }
                        Dvv += product;
                    }
                }
                const termDuu = combine(D, Nuu, -Duu, N);
                const termDuv = new Vector(3);
                for (let k = 0; k < 3; ++k) {
                    termDuv.values[k] = D * Nuv.values[k] - Dv * Nu.values[k]
                        - Du * Nv.values[k] - Duv * N.values[k];
                }
                const termDvv = combine(D, Nvv, -Dvv, N);
                values[3] = divideVec(combine(D, termDuu, -2 * Du, numerDU), Dcub);
                values[4] = divideVec(combine(D, termDuv, 2 * Du * Dv, N), Dcub);
                values[5] = divideVec(combine(D, termDvv, -2 * Dv, numerDV), Dcub);
            }
        }
    }
}

function divideVec(v: Vector, s: number): Vector {
    const result = new Vector(v.size);
    for (let k = 0; k < v.size; ++k) {
        result.values[k] = v.values[k] / s;
    }
    return result;
}

function combine(a0: number, v0: Vector, a1: number, v1: Vector): Vector {
    const result = new Vector(v0.size);
    for (let k = 0; k < v0.size; ++k) {
        result.values[k] = a0 * v0.values[k] + a1 * v1.values[k];
    }
    return result;
}

export class NURBSHalfSphereDegree3 extends NURBSSurface {
    constructor() {
        super(3, [new BasisFunctionInput(4, 3), new BasisFunctionInput(4, 3)]);

        // weight[j][i] is mWeights[i + 4 * j], 0 <= i < 4, 0 <= j < 4
        const oneThird = 1 / 3;
        const oneNinth = 1 / 9;
        this.mWeights[0] = 1;
        this.mWeights[1] = oneThird;
        this.mWeights[2] = oneThird;
        this.mWeights[3] = 1;
        this.mWeights[4] = oneThird;
        this.mWeights[5] = oneNinth;
        this.mWeights[6] = oneNinth;
        this.mWeights[7] = oneThird;
        this.mWeights[8] = oneThird;
        this.mWeights[9] = oneNinth;
        this.mWeights[10] = oneNinth;
        this.mWeights[11] = oneThird;
        this.mWeights[12] = 1;
        this.mWeights[13] = oneThird;
        this.mWeights[14] = oneThird;
        this.mWeights[15] = 1;

        // control[j][i] is mControls[i + 4 * j], 0 <= i < 4, 0 <= j < 4
        this.mControls[0] = Vector.fromArray([0, 0, 1]);
        this.mControls[1] = Vector.fromArray([0, 0, 1]);
        this.mControls[2] = Vector.fromArray([0, 0, 1]);
        this.mControls[3] = Vector.fromArray([0, 0, 1]);
        this.mControls[4] = Vector.fromArray([2, 0, 1]);
        this.mControls[5] = Vector.fromArray([2, 4, 1]);
        this.mControls[6] = Vector.fromArray([-2, 4, 1]);
        this.mControls[7] = Vector.fromArray([-2, 0, 1]);
        this.mControls[8] = Vector.fromArray([2, 0, -1]);
        this.mControls[9] = Vector.fromArray([2, 4, -1]);
        this.mControls[10] = Vector.fromArray([-2, 4, -1]);
        this.mControls[11] = Vector.fromArray([-2, 0, -1]);
        this.mControls[12] = Vector.fromArray([0, 0, -1]);
        this.mControls[13] = Vector.fromArray([0, 0, -1]);
        this.mControls[14] = Vector.fromArray([0, 0, -1]);
        this.mControls[15] = Vector.fromArray([0, 0, -1]);
    }
}

export class NURBSFullSphereDegree3 extends NURBSSurface {
    constructor() {
        super(3, [new BasisFunctionInput(4, 3),
            NURBSFullSphereDegree3.createBasisFunctionInputV()]);

        // weight[j][i] is mWeights[i + 4 * j], 0 <= i < 4, 0 <= j < 7
        const oneThird = 1 / 3;
        const oneNinth = 1 / 9;
        this.mWeights[0] = 1;
        this.mWeights[4] = oneThird;
        this.mWeights[8] = oneThird;
        this.mWeights[12] = 1;
        this.mWeights[16] = oneThird;
        this.mWeights[20] = oneThird;
        this.mWeights[24] = 1;
        this.mWeights[1] = oneThird;
        this.mWeights[5] = oneNinth;
        this.mWeights[9] = oneNinth;
        this.mWeights[13] = oneThird;
        this.mWeights[17] = oneNinth;
        this.mWeights[21] = oneNinth;
        this.mWeights[25] = oneThird;
        this.mWeights[2] = oneThird;
        this.mWeights[6] = oneNinth;
        this.mWeights[10] = oneNinth;
        this.mWeights[14] = oneThird;
        this.mWeights[18] = oneNinth;
        this.mWeights[22] = oneNinth;
        this.mWeights[26] = oneThird;
        this.mWeights[3] = 1;
        this.mWeights[7] = oneThird;
        this.mWeights[11] = oneThird;
        this.mWeights[15] = 1;
        this.mWeights[19] = oneThird;
        this.mWeights[23] = oneThird;
        this.mWeights[27] = 1;

        // control[j][i] is mControls[i + 4 * j], 0 <= i < 4, 0 <= j < 7
        this.mControls[0] = Vector.fromArray([0, 0, 1]);
        this.mControls[4] = Vector.fromArray([0, 0, 1]);
        this.mControls[8] = Vector.fromArray([0, 0, 1]);
        this.mControls[12] = Vector.fromArray([0, 0, 1]);
        this.mControls[16] = Vector.fromArray([0, 0, 1]);
        this.mControls[20] = Vector.fromArray([0, 0, 1]);
        this.mControls[24] = Vector.fromArray([0, 0, 1]);
        this.mControls[1] = Vector.fromArray([2, 0, 1]);
        this.mControls[5] = Vector.fromArray([2, 4, 1]);
        this.mControls[9] = Vector.fromArray([-2, 4, 1]);
        this.mControls[13] = Vector.fromArray([-2, 0, 1]);
        this.mControls[17] = Vector.fromArray([-2, -4, 1]);
        this.mControls[21] = Vector.fromArray([2, -4, 1]);
        this.mControls[25] = Vector.fromArray([2, 0, 1]);
        this.mControls[2] = Vector.fromArray([2, 0, -1]);
        this.mControls[6] = Vector.fromArray([2, 4, -1]);
        this.mControls[10] = Vector.fromArray([-2, 4, -1]);
        this.mControls[14] = Vector.fromArray([-2, 0, -1]);
        this.mControls[18] = Vector.fromArray([-2, -4, -1]);
        this.mControls[22] = Vector.fromArray([2, -4, -1]);
        this.mControls[26] = Vector.fromArray([2, 0, -1]);
        this.mControls[3] = Vector.fromArray([0, 0, -1]);
        this.mControls[7] = Vector.fromArray([0, 0, -1]);
        this.mControls[11] = Vector.fromArray([0, 0, -1]);
        this.mControls[15] = Vector.fromArray([0, 0, -1]);
        this.mControls[19] = Vector.fromArray([0, 0, -1]);
        this.mControls[23] = Vector.fromArray([0, 0, -1]);
        this.mControls[27] = Vector.fromArray([0, 0, -1]);
    }

    private static createBasisFunctionInputV(): BasisFunctionInput {
        const input = new BasisFunctionInput();
        input.numControls = 7;
        input.degree = 3;
        input.uniform = true;
        input.periodic = false;
        input.numUniqueKnots = 3;
        input.uniqueKnots = [
            new UniqueKnot(0, 4),
            new UniqueKnot(0.5, 3),
            new UniqueKnot(1, 4)
        ];
        return input;
    }
}
