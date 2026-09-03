// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ImplicitCurve2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The curve is defined by F(x,y) = 0. In all member functions it is the
// application's responsibility to ensure that (x,y) is a solution to F = 0.
// The class is abstract, so you must derive from it and implement the
// function and derivative evaluations.
//
// Port notes: the pure virtual C++ member functions F, FX, FY, FXX, FXY and
// FYY become the abstract methods f, fx, fy, fxx, fxy and fyy (camelCase per
// the API mapping). 'GetFrame', which writes the tangent and normal through
// reference parameters, becomes 'getFrame(position)' returning a named object
// literal. 'GetCurvature' returns '{ valid, curvature }' rather than writing
// through a reference and returning a bool. Matrix2x2<T> becomes a
// 2-by-2 'Matrix' (the B56 fixed-size matrix precedent).

import { Matrix } from './Matrix.js';
import { Vector } from './Vector.js';
import { computeOrthogonalComplement2 } from './Vector2.js';

export interface ImplicitCurve2Frame {
    tangent: Vector;
    normal: Vector;
}

export abstract class ImplicitCurve2 {
    protected constructor() {
    }

    // Evaluate the implicit function.
    abstract f(position: Vector): number;

    // Evaluate the first-order partial derivatives.
    abstract fx(position: Vector): number;
    abstract fy(position: Vector): number;

    // Evaluate the second-order partial derivatives.
    abstract fxx(position: Vector): number;
    abstract fxy(position: Vector): number;
    abstract fyy(position: Vector): number;

    // Verify the point is on the curve within the tolerance specified by
    // epsilon.
    isOnCurve(position: Vector, epsilon: number): boolean {
        return Math.abs(this.f(position)) <= epsilon;
    }

    // Compute all first-order derivatives.
    getGradient(position: Vector): Vector {
        const fx = this.fx(position);
        const fy = this.fy(position);
        const gradient = new Vector(2);
        gradient.values[0] = fx;
        gradient.values[1] = fy;
        return gradient;
    }

    // Compute all second-order derivatives.
    getHessian(position: Vector): Matrix {
        const fxx = this.fxx(position);
        const fxy = this.fxy(position);
        const fyy = this.fyy(position);
        return Matrix.fromArray(2, 2, [fxx, fxy, fxy, fyy]);
    }

    // Compute a coordinate frame. The set {T, N} is a right-handed
    // orthonormal basis.
    getFrame(position: Vector): ImplicitCurve2Frame {
        const basis: Vector[] = [this.getGradient(position), new Vector(2)];
        computeOrthogonalComplement2(1, basis);
        return { tangent: basis[1], normal: basis[0] };
    }

    // Compute the curvature at a point on the curve. The 'valid' field is
    // false when the gradient is the zero vector, in which case the
    // curvature is set to zero.
    getCurvature(position: Vector): { valid: boolean, curvature: number } {
        // The curvature is
        // (-Fy^2*Fxx + 2*Fx*Fy*Fxy - Fx^2*Fyy) / (Fx^2+Fy^2)^{3/2}

        // Evaluate the first derivatives.
        const fx = this.fx(position);
        const fy = this.fy(position);

        // Evaluate the denominator.
        const fxSqr = fx * fx;
        const fySqr = fy * fy;
        const denom = Math.pow(fxSqr + fySqr, 1.5);
        if (denom === 0) {
            return { valid: false, curvature: 0 };
        }

        // Evaluate the second derivatives.
        const fxx = this.fxx(position);
        const fxy = this.fxy(position);
        const fyy = this.fyy(position);

        // Evaluate the numerator.
        const numer = -fySqr * fxx + 2 * fx * fy * fxy - fxSqr * fyy;

        return { valid: true, curvature: numer / denom };
    }
}
