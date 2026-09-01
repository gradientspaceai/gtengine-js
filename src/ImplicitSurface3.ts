// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ImplicitSurface3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The surface is defined by F(x,y,z) = 0. In all member functions it is the
// application's responsibility to ensure that (x,y,z) is a solution to F = 0.
// The class is abstract, so you must derive from it and implement the
// function and derivative evaluations.
//
// The computation of principal curvature and principal directions is based
// on the document
// https://www.geometrictools.com/Documentation/PrincipalCurvature.pdf
//
// Port notes: the upstream member functions F, FX, ..., FZZ become the
// camelCase abstract methods f, fx, ..., fzz. The output reference
// parameters of GetFrame and GetPrincipalInformation become returned object
// literals per PORTING.md. Positions and directions are 3-dimensional
// Vector objects and the Hessian is a 3-by-3 Matrix.

import { Matrix, multiplyATB, mulMatrix, divMatrix } from './Matrix';
import { SymmetricEigensolver2x2 } from './SymmetricEigensolver2x2';
import { Vector, normalize } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

// The set {tangent0, tangent1, normal} is a right-handed orthonormal basis.
export interface ImplicitSurface3Frame {
    tangent0: Vector;
    tangent1: Vector;
    normal: Vector;
}

// The scalars are the principal curvatures and the vectors are the
// corresponding principal directions. When the gradient is the zero vector,
// the differential geometric quantities are undefined; 'valid' is false and
// the curvatures and directions are zero.
export interface ImplicitSurface3PrincipalInformation {
    valid: boolean;
    curvature0: number;
    curvature1: number;
    direction0: Vector;
    direction1: Vector;
}

// Abstract base class. Upstream declares a protected default constructor to
// prevent direct instantiation; in TypeScript the 'abstract' modifier already
// does that, and a protected constructor would keep derived instances from
// being used where an ImplicitSurface3 is expected.
export abstract class ImplicitSurface3 {
    // Evaluate the implicit function.
    abstract f(position: Vector): number;

    // Evaluate the first-order partial derivatives.
    abstract fx(position: Vector): number;
    abstract fy(position: Vector): number;
    abstract fz(position: Vector): number;

    // Evaluate the second-order partial derivatives.
    abstract fxx(position: Vector): number;
    abstract fxy(position: Vector): number;
    abstract fxz(position: Vector): number;
    abstract fyy(position: Vector): number;
    abstract fyz(position: Vector): number;
    abstract fzz(position: Vector): number;

    // Verify the point is on the surface within the tolerance specified by
    // epsilon.
    isOnSurface(position: Vector, epsilon: number): boolean {
        return Math.abs(this.f(position)) <= epsilon;
    }

    // Compute all first-order derivatives.
    getGradient(position: Vector): Vector {
        const fx = this.fx(position);
        const fy = this.fy(position);
        const fz = this.fz(position);
        return Vector.fromArray([fx, fy, fz]);
    }

    // Compute all second-order derivatives.
    getHessian(position: Vector): Matrix {
        const fxx = this.fxx(position);
        const fxy = this.fxy(position);
        const fxz = this.fxz(position);
        const fyy = this.fyy(position);
        const fyz = this.fyz(position);
        const fzz = this.fzz(position);
        return Matrix.fromArray(3, 3, [
            fxx, fxy, fxz,
            fxy, fyy, fyz,
            fxz, fyz, fzz
        ]);
    }

    // Compute a coordinate frame. The set {T0,T1,N} is a right-handed
    // orthonormal basis.
    getFrame(position: Vector): ImplicitSurface3Frame {
        const basis: Vector[] = [
            this.getGradient(position),
            new Vector(3),
            new Vector(3)
        ];
        computeOrthogonalComplement3(1, basis);
        return { tangent0: basis[1], tangent1: basis[2], normal: basis[0] };
    }

    // Differential geometric quantities. The returned scalars are the
    // principal curvatures and the returned vectors are the corresponding
    // principal directions.
    getPrincipalInformation(position: Vector): ImplicitSurface3PrincipalInformation {
        // Compute the normal N.
        const zero = 0;
        const normal = this.getGradient(position);
        const gradientLength = normalize(normal);
        if (gradientLength === zero) {
            return {
                valid: false,
                curvature0: zero,
                curvature1: zero,
                direction0: new Vector(3),
                direction1: new Vector(3)
            };
        }

        // Compute the matrix A.
        const A = divMatrix(this.getHessian(position), gradientLength);

        // Solve for the eigensystem of equation (8) of the PDF referenced at
        // the top of this file.
        const basis: Vector[] = [normal, new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        // basis[1] = tangent0
        // basis[2] = tangent1
        const J = new Matrix(3, 2);
        J.setCol(0, basis[1]);
        J.setCol(1, basis[2]);
        const barA = multiplyATB(J, mulMatrix(A, J) as Matrix);

        const eigensolver = new SymmetricEigensolver2x2();
        const { evals, evecs } = eigensolver.solve(barA.get(0, 0),
            barA.get(0, 1), barA.get(1, 1), +1);
        const v0 = Vector.fromArray([evecs[0][0], evecs[0][1]]);
        const v1 = Vector.fromArray([evecs[1][0], evecs[1][1]]);
        return {
            valid: true,
            curvature0: evals[0],
            curvature1: evals[1],
            direction0: mulMatrix(J, v0) as Vector,
            direction1: mulMatrix(J, v1) as Vector
        };
    }
}
