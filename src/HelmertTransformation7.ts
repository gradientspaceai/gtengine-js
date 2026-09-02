// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HelmertTransformation7.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Implementation of the 7-parameter Helmert transformation. It is designed
// to rotate, translate, and uniformly scale one 3D point set to be as close
// as possible to another 3D point set. Details are provided in
//   https://www.geometrictools.com/Documentation/HelmertTransformation.pdf
//
// Port notes: upstream 'Vector3<T>' is the runtime-sized Vector with size 3
// and 'Matrix3x3<T>' is the runtime-sized Matrix with size 3-by-3 (see
// Matrix3x3.ts). The output reference parameters of Execute (rotate,
// translate, scale, function) become the fields of the returned
// HelmertTransformation7Result; the upstream 'function' output is named
// 'functionValue' here because 'function' is a reserved word. The
// UpdateEulerAngle* helpers take 'F' by reference upstream; the port returns
// the updated value alongside the boolean.

import { logAssert } from './Logger';
import { Matrix, mulMatrix, multiplyAB } from './Matrix';
import { Vector, add, div, dot, mul, sub } from './Vector';

export interface HelmertTransformation7Result {
    // The number of iterations actually performed.
    iterations: number;

    // The rotation matrix R of the transformation y = s * R * x + t.
    rotate: Matrix;

    // The translation t of the transformation.
    translate: Vector;

    // The uniform scale s of the transformation.
    scale: number;

    // The value of the objective function at the returned parameters. This
    // is the upstream 'function' output parameter.
    functionValue: number;
}

function identity3x3(): Matrix {
    const m = new Matrix(3, 3);
    m.makeIdentity();
    return m;
}

export class HelmertTransformation7 {
    private mNumPoints: number;
    private mU: Vector[];
    private mV: Vector[];
    private mLeft: Vector[];
    private mRight: Vector[];
    private mRotate: Matrix;
    private mRotate0: Matrix;
    private mRotate1: Matrix;
    private mRotate2: Matrix;
    private mTranslate: Vector;

    constructor() {
        this.mNumPoints = 0;
        this.mU = [];
        this.mV = [];
        this.mLeft = [];
        this.mRight = [];
        this.mRotate = new Matrix(3, 3);
        this.mRotate0 = new Matrix(3, 3);
        this.mRotate1 = new Matrix(3, 3);
        this.mRotate2 = new Matrix(3, 3);
        this.mTranslate = new Vector(3);
    }

    // The input points p[i] and q[i] must correspond for 0 <= i < n, where n
    // is the number of points (n = p.length = q.length). The outputs require
    // 7 parameters: 3 for rotation (Euler angles), 3 for translation, and 1
    // for uniform scale.
    execute(p: readonly Vector[], q: readonly Vector[],
        numIterations: number): HelmertTransformation7Result {
        this.mNumPoints = p.length;
        logAssert(this.mNumPoints >= 7 && q.length === this.mNumPoints,
            'Invalid input.');

        // Translate the centroid of the q-points to the origin. This
        // simplifies the function to be minimized, the only parameter being
        // the rotation matrix (function of 3 Euler angles). Also, compute the
        // centroid of the p-points.
        let pAverage = new Vector(3);
        let qAverage = new Vector(3);
        for (let i = 0; i < this.mNumPoints; ++i) {
            pAverage = add(pAverage, p[i]);
            qAverage = add(qAverage, q[i]);
        }
        pAverage = div(pAverage, this.mNumPoints);
        qAverage = div(qAverage, this.mNumPoints);

        // Translate by the centroid of q. The p-values are translated to
        // u-values and the q-values are translated to the v-values. The
        // average of the v-values is the zero vector.
        this.mU = new Array<Vector>(this.mNumPoints);
        this.mV = new Array<Vector>(this.mNumPoints);
        this.mLeft = new Array<Vector>(this.mNumPoints);
        this.mRight = new Array<Vector>(this.mNumPoints);
        for (let i = 0; i < this.mNumPoints; ++i) {
            this.mU[i] = sub(p[i], qAverage);
            this.mV[i] = sub(q[i], qAverage);
            this.mLeft[i] = new Vector(3);
            this.mRight[i] = new Vector(3);
        }

        // Upstream computes a 3-by-3 matrix A = sum_i u_i * v_i^T here but
        // never uses it. The dead code is omitted by the port; see the
        // "Upstream bug suspects" notes.

        // The initial rotation matrix is the identity.
        this.mRotate = identity3x3();
        this.mRotate0 = identity3x3();
        this.mRotate1 = identity3x3();
        this.mRotate2 = identity3x3();

        // The translation does not vary during the iterations.
        this.mTranslate = sub(pAverage, qAverage);

        let F = this.updateF(this.mRotate);
        let iteration = 0;
        for (iteration = 0; iteration < numIterations; ++iteration) {
            const result0 = this.updateEulerAngle0(F);
            F = result0.F;
            const result1 = this.updateEulerAngle1(F);
            F = result1.F;
            const result2 = this.updateEulerAngle2(F);
            F = result2.F;
            if (!result0.updated && !result1.updated && !result2.updated) {
                break;
            }
        }

        const rotate = this.mRotate;
        const scale = this.updateScale(this.mRotate);
        const translate = sub(add(this.mTranslate, qAverage),
            mul(mulMatrix(rotate, qAverage), scale));
        return { iterations: iteration, rotate, translate, scale, functionValue: F };
    }

    private updateScale(rotate: Matrix): number {
        let numer = 0;
        let denom = 0;
        for (let i = 0; i < this.mNumPoints; ++i) {
            numer += dot(this.mU[i], mulMatrix(rotate, this.mV[i]));
            denom += dot(this.mV[i], this.mV[i]);
        }
        return numer / denom;
    }

    private updateF(rotate: Matrix): number {
        const scale = this.updateScale(rotate);
        let functionValue = 0;
        for (let i = 0; i < this.mNumPoints; ++i) {
            const term = sub(add(mul(mulMatrix(rotate, this.mV[i]), scale),
                this.mTranslate), this.mU[i]);
            functionValue += dot(term, term);
        }
        functionValue /= this.mNumPoints;
        return functionValue;
    }

    private updateEulerAngle0(F: number): { updated: boolean, F: number } {
        const R1R2 = multiplyAB(this.mRotate1, this.mRotate2);
        const left = this.mU;
        const right = this.mRight;
        for (let i = 0; i < this.mNumPoints; ++i) {
            right[i] = mulMatrix(R1R2, this.mV[i]);
        }

        let sn = 0;
        let cs = 0;
        for (let i = 0; i < this.mNumPoints; ++i) {
            sn += left[i].values[1] * right[i].values[0]
                - left[i].values[0] * right[i].values[1];
            cs += left[i].values[0] * right[i].values[0]
                + left[i].values[1] * right[i].values[1];
        }

        const length = Math.sqrt(sn * sn + cs * cs);
        if (length > 0) {
            sn /= length;
            cs /= length;
        } else {
            sn = 0;
            cs = 1;
        }

        const rotate0 = new Matrix(3, 3);
        rotate0.set(0, 0, cs); rotate0.set(0, 1, -sn); rotate0.set(0, 2, 0);
        rotate0.set(1, 0, sn); rotate0.set(1, 1, cs); rotate0.set(1, 2, 0);
        rotate0.set(2, 0, 0); rotate0.set(2, 1, 0); rotate0.set(2, 2, 1);

        const updateRotate = multiplyAB(rotate0, R1R2);
        const updateF = this.updateF(updateRotate);
        if (updateF < F) {
            this.mRotate0 = rotate0;
            this.mRotate = updateRotate;
            return { updated: true, F: updateF };
        }
        return { updated: false, F };
    }

    private updateEulerAngle1(F: number): { updated: boolean, F: number } {
        const left = this.mLeft;
        const right = this.mRight;
        for (let i = 0; i < this.mNumPoints; ++i) {
            left[i] = mulMatrix(this.mU[i], this.mRotate0);
            right[i] = mulMatrix(this.mRotate2, this.mV[i]);
        }

        let sn = 0;
        let cs = 0;
        for (let i = 0; i < this.mNumPoints; ++i) {
            sn += left[i].values[0] * right[i].values[2]
                - left[i].values[2] * right[i].values[0];
            cs += left[i].values[0] * right[i].values[0]
                + left[i].values[2] * right[i].values[2];
        }

        const length = Math.sqrt(sn * sn + cs * cs);
        if (length > 0) {
            sn /= length;
            cs /= length;
        } else {
            sn = 0;
            cs = 1;
        }

        const rotate1 = new Matrix(3, 3);
        rotate1.set(0, 0, cs); rotate1.set(0, 1, 0); rotate1.set(0, 2, sn);
        rotate1.set(1, 0, 0); rotate1.set(1, 1, 1); rotate1.set(1, 2, 0);
        rotate1.set(2, 0, -sn); rotate1.set(2, 1, 0); rotate1.set(2, 2, cs);

        const updateRotate = multiplyAB(multiplyAB(this.mRotate0, rotate1),
            this.mRotate2);
        const updateF = this.updateF(updateRotate);
        if (updateF < F) {
            this.mRotate1 = rotate1;
            this.mRotate = updateRotate;
            return { updated: true, F: updateF };
        }
        return { updated: false, F };
    }

    private updateEulerAngle2(F: number): { updated: boolean, F: number } {
        const R0R1 = multiplyAB(this.mRotate0, this.mRotate1);
        const left = this.mLeft;
        const right = this.mV;
        for (let i = 0; i < this.mNumPoints; ++i) {
            left[i] = mulMatrix(this.mU[i], R0R1);
        }

        let sn = 0;
        let cs = 0;
        for (let i = 0; i < this.mNumPoints; ++i) {
            sn += left[i].values[2] * right[i].values[1]
                - left[i].values[1] * right[i].values[2];
            cs += left[i].values[1] * right[i].values[1]
                + left[i].values[2] * right[i].values[2];
        }

        const length = Math.sqrt(sn * sn + cs * cs);
        if (length > 0) {
            sn /= length;
            cs /= length;
        } else {
            sn = 0;
            cs = 1;
        }

        const rotate2 = new Matrix(3, 3);
        rotate2.set(0, 0, 1); rotate2.set(0, 1, 0); rotate2.set(0, 2, 0);
        rotate2.set(1, 0, 0); rotate2.set(1, 1, cs); rotate2.set(1, 2, -sn);
        rotate2.set(2, 0, 0); rotate2.set(2, 1, sn); rotate2.set(2, 2, cs);

        const updateRotate = multiplyAB(R0R1, rotate2);
        const updateF = this.updateF(updateRotate);
        if (updateF < F) {
            this.mRotate2 = rotate2;
            this.mRotate = updateRotate;
            return { updated: true, F: updateF };
        }
        return { updated: false, F };
    }
}
