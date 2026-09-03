// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpSphere2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Interpolation of a scalar-valued function defined on a sphere. Although
// the sphere lives in 3D, the interpolation is a 2D method whose input
// points are angles (theta,phi) from spherical coordinates. The domains of
// the angles are -pi <= theta <= pi and 0 <= phi <= pi.
//
// Port notes:
// * Upstream declares a variadic class template with two specializations,
//   the deprecated IntpSphere2<InputType, ComputeType, RationalType> and the
//   replacement IntpSphere2<T>. Only the replacement is ported, matching the
//   Delaunay2 and Delaunay2Mesh ports.
// * The static GetSphericalCoordinates writes theta and phi through
//   reference parameters; the port returns { theta, phi } and keeps the
//   function as a static method so the name does not enter the flat export
//   namespace.
// * The C++ operator()(theta, phi, F&) returning bool becomes
//   evaluate(theta, phi), which returns { valid, F }; F is meaningless when
//   valid is false.

import { GTE_C_PI, GTE_C_TWO_PI } from './Constants.js';
import { Delaunay2 } from './Delaunay2.js';
import { Delaunay2Mesh } from './Delaunay2Mesh.js';
import { IntpQuadraticNonuniform2 } from './IntpQuadraticNonuniform2.js';
import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export interface IntpSphere2Result {
    // Valid is true if and only if the input point is in the convex hull of
    // the input (theta,phi) array, in which case the interpolation is valid.
    valid: boolean;

    // The interpolated function value; meaningful only when valid is true.
    F: number;
}

export interface SphericalCoordinates2 {
    theta: number;
    phi: number;
}

export class IntpSphere2 {
    private mWrapAngles: Vector[];
    private mWrapF: number[];
    private mDelaunay: Delaunay2;
    private mMesh: Delaunay2Mesh;
    private mInterp: IntpQuadraticNonuniform2;

    // Construction. For complete spherical coverage, include the two
    // antipodal (theta,phi) points (-pi,0,F(-pi,0)) and (-pi,pi,F(-pi,pi))
    // in the input data. These correspond to the sphere poles x = 0, y = 0,
    // and |z| = 1.
    constructor(theta: readonly number[], phi: readonly number[],
        F: readonly number[]) {
        const numPoints = theta.length;
        logAssert(numPoints > 0 && phi.length === numPoints
            && F.length === numPoints,
            'IntpSphere2: the theta, phi and F arrays must have the same '
            + 'positive length.');

        // Copy the input data. The larger arrays are used to support
        // wrap-around in the Delaunay triangulation for the interpolator.
        const totalPoints = 3 * numPoints;
        this.mWrapAngles = new Array<Vector>(totalPoints);
        this.mWrapF = new Array<number>(totalPoints).fill(0);
        for (let i = 0; i < numPoints; ++i) {
            this.mWrapAngles[i] = Vector.fromArray([theta[i], phi[i]]);
            this.mWrapF[i] = F[i];
        }

        // Use periodicity to get wrap-around in the Delaunay triangulation.
        for (let i0 = 0, i1 = numPoints, i2 = 2 * numPoints; i0 < numPoints;
            ++i0, ++i1, ++i2) {
            const a0 = this.mWrapAngles[i0].values;
            this.mWrapAngles[i1] =
                Vector.fromArray([a0[0] + GTE_C_TWO_PI, a0[1]]);
            this.mWrapAngles[i2] =
                Vector.fromArray([a0[0] - GTE_C_TWO_PI, a0[1]]);
            this.mWrapF[i1] = this.mWrapF[i0];
            this.mWrapF[i2] = this.mWrapF[i0];
        }

        this.mDelaunay = new Delaunay2();
        this.mDelaunay.compute(this.mWrapAngles);
        this.mMesh = new Delaunay2Mesh(this.mDelaunay);
        this.mInterp = IntpQuadraticNonuniform2.fromSpatialDelta(this.mMesh,
            this.mWrapF, 1);
    }

    // Spherical coordinates are
    //   x = cos(theta)*sin(phi)
    //   y = sin(theta)*sin(phi)
    //   z = cos(phi)
    // for -pi <= theta <= pi, 0 <= phi <= pi. The application can use this
    // function to convert unit-length vectors (x,y,z) to (theta,phi). The
    // input (x,y,z) is assumed to be unit length.
    static getSphericalCoordinates(x: number, y: number, z: number):
        SphericalCoordinates2 {
        if (z < 1) {
            if (z > -1) {
                return { theta: Math.atan2(y, x), phi: Math.acos(z) };
            }
            return { theta: -GTE_C_PI, phi: GTE_C_PI };
        }
        return { theta: -GTE_C_PI, phi: 0 };
    }

    // The 'valid' field is true if and only if the input point is in the
    // convex hull of the input (theta,phi) array, in which case the
    // interpolation is valid.
    evaluate(theta: number, phi: number): IntpSphere2Result {
        const angles = Vector.fromArray([theta, phi]);
        const result = this.mInterp.evaluate(angles);
        return { valid: result.valid, F: result.F };
    }
}
