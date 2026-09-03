// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DarbouxFrame.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Darboux frame of a parametric surface in 3D: the moving orthonormal
// frame {T0, T1, N} at a surface point, together with the principal
// curvatures and principal directions.
//
// Port notes (following the FrenetFrame precedent set by B64):
// - Upstream takes a 'std::shared_ptr<ParametricSurface<3, Real>>'; the port
//   takes the surface object directly (the caller keeps it alive).
// - The 'operator()' that writes position/tangent0/tangent1/normal through
//   reference parameters becomes 'compute(u, v)' returning a named object
//   literal; 'GetPrincipalInformation' likewise returns an object literal.
// - The surface's dimension is asserted at construction because the port's
//   ParametricSurface carries its dimension at runtime rather than as a
//   template parameter.
// - Matrix2x2<Real> becomes 'new Matrix(2, 2)' with determinant2x2 (see the
//   fixed-size matrix naming precedent from B56).

import { logAssert } from './Logger.js';
import { Matrix } from './Matrix.js';
import { determinant2x2 } from './Matrix2x2.js';
import { ParametricSurface } from './ParametricSurface.js';
import { Vector, add, dot, mul, normalize } from './Vector.js';
import { cross, unitCross } from './Vector3.js';

export interface DarbouxFrame3Result {
    position: Vector;
    tangent0: Vector;
    tangent1: Vector;
    normal: Vector;
}

export interface DarbouxFrame3PrincipalResult {
    curvature0: number;
    curvature1: number;
    direction0: Vector;
    direction1: Vector;
}

export class DarbouxFrame3 {
    private mSurface: ParametricSurface;

    // Construction. The surface must persist as long as the DarbouxFrame3
    // object does.
    constructor(surface: ParametricSurface) {
        logAssert(surface.getDimension() === 3,
            'DarbouxFrame3: the surface must be 3-dimensional.');
        this.mSurface = surface;
    }

    getSurface(): ParametricSurface {
        return this.mSurface;
    }

    // Get a coordinate frame, {T0, T1, N}. At a nondegenerate surface point,
    // dX/du and dX/dv are linearly independent tangent vectors. The frame is
    // constructed as
    //   T0 = (dX/du)/|dX/du|
    //   N  = Cross(dX/du,dX/dv)/|Cross(dX/du,dX/dv)|
    //   T1 = Cross(N, T0)
    // so that {T0, T1, N} is a right-handed orthonormal set.
    compute(u: number, v: number): DarbouxFrame3Result {
        const jet = this.mSurface.createJet();
        this.mSurface.evaluate(u, v, 1, jet);
        const position = jet[0].clone();
        const tangent0 = jet[1].clone();
        normalize(tangent0);
        let tangent1 = jet[2].clone();
        normalize(tangent1);
        const normal = unitCross(tangent0, tangent1);
        tangent1 = cross(normal, tangent0);
        return { position, tangent0, tangent1, normal };
    }

    // Compute the principal curvatures and principal directions.
    getPrincipalInformation(u: number, v: number): DarbouxFrame3PrincipalResult {
        // Tangents:  T0 = (x_u,y_u,z_u), T1 = (x_v,y_v,z_v)
        // Normal:    N = Cross(T0,T1)/Length(Cross(T0,T1))
        // Metric Tensor:    G = +-                      -+
        //                       | Dot(T0,T0)  Dot(T0,T1) |
        //                       | Dot(T1,T0)  Dot(T1,T1) |
        //                       +-                      -+
        //
        // Curvature Tensor:  B = +-                          -+
        //                        | -Dot(N,T0_u)  -Dot(N,T0_v) |
        //                        | -Dot(N,T1_u)  -Dot(N,T1_v) |
        //                        +-                          -+
        //
        // Principal curvatures k are the generalized eigenvalues of
        //
        //     Bw = kGw
        //
        // If k is a curvature and w=(a,b) is the corresponding solution to
        // Bw = kGw, then the principal direction as a 3D vector is
        // d = a*U+b*V.
        //
        // Let k1 and k2 be the principal curvatures. The mean curvature is
        // (k1+k2)/2 and the Gaussian curvature is k1*k2.

        // Compute derivatives.
        const jet = this.mSurface.createJet();
        this.mSurface.evaluate(u, v, 2, jet);
        const derU = jet[1];
        const derV = jet[2];
        const derUU = jet[3];
        const derUV = jet[4];
        const derVV = jet[5];

        // Compute the metric tensor.
        const metricTensor = new Matrix(2, 2);
        metricTensor.set(0, 0, dot(jet[1], jet[1]));
        metricTensor.set(0, 1, dot(jet[1], jet[2]));
        metricTensor.set(1, 0, metricTensor.get(0, 1));
        metricTensor.set(1, 1, dot(jet[2], jet[2]));

        // Compute the curvature tensor.
        const normal = unitCross(jet[1], jet[2]);
        const curvatureTensor = new Matrix(2, 2);
        curvatureTensor.set(0, 0, -dot(normal, derUU));
        curvatureTensor.set(0, 1, -dot(normal, derUV));
        curvatureTensor.set(1, 0, curvatureTensor.get(0, 1));
        curvatureTensor.set(1, 1, -dot(normal, derVV));

        // Characteristic polynomial is 0 = det(B-kG) = c2*k^2+c1*k+c0.
        const c0 = determinant2x2(curvatureTensor);
        const c1 = 2 * curvatureTensor.get(0, 1) * metricTensor.get(0, 1)
            - curvatureTensor.get(0, 0) * metricTensor.get(1, 1)
            - curvatureTensor.get(1, 1) * metricTensor.get(0, 0);
        const c2 = determinant2x2(metricTensor);

        // Principal curvatures are roots of the characteristic polynomial.
        const temp = Math.sqrt(Math.max(c1 * c1 - 4 * c0 * c2, 0));
        const mult = 0.5 / c2;
        const curvature0 = -mult * (c1 + temp);
        const curvature1 = -mult * (c1 - temp);

        // Principal directions are solutions to (B-kG)w = 0,
        // w1 = (b12-k1*g12,-(b11-k1*g11)) OR (b22-k1*g22,-(b12-k1*g12)).
        let direction0: Vector;
        let a0 = curvatureTensor.get(0, 1) - curvature0 * metricTensor.get(0, 1);
        let a1 = curvature0 * metricTensor.get(0, 0) - curvatureTensor.get(0, 0);
        let len = Math.sqrt(a0 * a0 + a1 * a1);
        if (len > 0) {
            direction0 = add(mul(a0, derU), mul(a1, derV));
        }
        else {
            a0 = curvatureTensor.get(1, 1) - curvature0 * metricTensor.get(1, 1);
            a1 = curvature0 * metricTensor.get(0, 1) - curvatureTensor.get(0, 1);
            len = Math.sqrt(a0 * a0 + a1 * a1);
            if (len > 0) {
                direction0 = add(mul(a0, derU), mul(a1, derV));
            }
            else {
                // Umbilic (surface is locally a sphere, any direction is
                // principal).
                direction0 = derU.clone();
            }
        }
        normalize(direction0);

        // The second tangent is the cross product of the first tangent and
        // the normal.
        const direction1 = cross(direction0, normal);

        return { curvature0, curvature1, direction0, direction1 };
    }
}
