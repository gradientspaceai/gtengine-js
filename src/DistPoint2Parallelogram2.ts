// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint2Parallelogram2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Implementation of a point-parallelogram distance and closest-point query.
// The details are described in
//   https://www.geometrictools.com/Documentation/DistancePointParallelogram.pdf
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector2<T>, Parallelogram2<T>>'
// becomes the class DistPoint2Parallelogram2 with the result type
// DistPoint2Parallelogram2Result. The public member function 'GetMinimizer'
// stays a public method; the private static 'Clamp' is the library's
// 'clamp' from Functions.ts, which has the same
// 'u <= umin ? umin : (u >= umax ? umax : u)' semantics.

import type { DCPQuery } from './DCPQuery';
import { clamp } from './Functions';
import { Matrix, multiplyATB, mulMatrix } from './Matrix';
import { inverse2x2 } from './Matrix2x2';
import type { Parallelogram2 } from './Parallelogram2';
import { Vector, add, mul, sub } from './Vector';

export interface DistPoint2Parallelogram2Result {
    // The point closest[0] is the query point. The point closest[1] is the
    // parallelogram point closest to the query point. The two points are the
    // same when the query point is contained by the parallelogram.
    distance: number;
    sqrDistance: number;
    closest: [Vector, Vector];
}

export class DistPoint2Parallelogram2
    implements DCPQuery<Vector, Parallelogram2, DistPoint2Parallelogram2Result> {
    compute(point: Vector, pgm: Parallelogram2):
        DistPoint2Parallelogram2Result {
        // For a parallelogram point X, let Y = {Dot(V0,X-C),Dot(V1,X-C)}.
        // Compute the quadratic function q(Y) = (Y-Z)^T * A * (Y-Z) / 2
        // where A = B^T * B is a symmetric matrix.
        const B = new Matrix(2, 2);
        B.setCol(0, pgm.axis[0]);
        B.setCol(1, pgm.axis[1]);
        const A = multiplyATB(B, B);

        // Transform the query point to parallelogram coordinates,
        // Z = Inverse(B) * (P - C).
        const delta = sub(point, pgm.center);
        const Z = mulMatrix(inverse2x2(B).inverse, delta) as Vector;

        // Get the minimizer for q(Y).
        const K = this.getMinimizer(A, Z);

        const closest0 = point.clone();
        const closest1 = add(pgm.center,
            add(mul(K.values[0], pgm.axis[0]), mul(K.values[1], pgm.axis[1])));
        const diff = sub(closest0, closest1);
        const sqrDistance = diff.values[0] * diff.values[0]
            + diff.values[1] * diff.values[1];
        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest: [closest0, closest1]
        };
    }

    // Compute the minimizer (in parallelogram coordinates) of the quadratic
    // function q. The domain is the square [-1,1]^2.
    getMinimizer(A: Matrix, Z: Vector): Vector {
        const a00 = A.get(0, 0);
        const a01 = A.get(0, 1);
        const a11 = A.get(1, 1);
        const z0 = Z.values[0];
        const z1 = Z.values[1];
        let root: number;
        const K = new Vector(2);

        if (z1 < -1) {
            // Examine the bottom edge.
            root = z0 - a01 * (-1 - z1) / a00;
            K.values[0] = clamp(root, -1, +1);
            K.values[1] = -1;

            if (z0 < -1) {
                if (K.values[0] === -1) {
                    // Examine the left edge.
                    root = z1 - a01 * (-1 - z0) / a11;
                    K.values[1] = clamp(root, -1, +1);
                    K.values[0] = -1;
                }
            }
            else if (+1 < z0) {
                if (K.values[0] === +1) {
                    // Examine the right edge.
                    root = z1 - a01 * (+1 - z0) / a11;
                    K.values[1] = clamp(root, -1, +1);
                    K.values[0] = +1;
                }
            }
        }
        else if (z1 <= +1) {
            if (z0 < -1) {
                // Examine the left edge.
                root = z1 - a01 * (-1 - z0) / a11;
                K.values[1] = clamp(root, -1, +1);
                K.values[0] = -1;
            }
            else if (z0 <= +1) {
                // The query point is inside the parallelogram.
                K.values[0] = z0;
                K.values[1] = z1;
            }
            else {
                // Examine the right edge.
                root = z1 - a01 * (+1 - z0) / a11;
                K.values[1] = clamp(root, -1, +1);
                K.values[0] = +1;
            }
        }
        else {
            // Examine the top edge.
            root = z0 - a01 * (+1 - z1) / a00;
            K.values[0] = clamp(root, -1, +1);
            K.values[1] = +1;

            if (z0 < -1) {
                if (K.values[0] === -1) {
                    // Examine the left edge.
                    root = z1 - a01 * (-1 - z0) / a11;
                    K.values[1] = clamp(root, -1, +1);
                    K.values[0] = -1;
                }
            }
            else if (+1 < z0) {
                if (K.values[0] === +1) {
                    // Examine the right edge.
                    root = z1 - a01 * (+1 - z0) / a11;
                    K.values[1] = clamp(root, -1, +1);
                    K.values[0] = +1;
                }
            }
        }

        return K;
    }
}
