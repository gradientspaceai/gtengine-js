// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Vector4.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'Vector4<Real>' is the alias 'Vector<4, Real>'; the
// port uses the runtime-sized Vector with size 4 (create with
// 'Vector.fromArray([x, y, z, w])'). The compile-time dimension checks become
// runtime asserts. Upstream ComputeOrthogonalComplement collides with the 2D
// and 3D versions in Vector2.h and Vector3.h under the library-wide flat
// export, so the port suffixes it with the dimension:
// computeOrthogonalComplement4.
//
// In Vector3.h, the Vector3 cross, unitCross and dotCross have a dimension
// that should be 3 or 4. The latter case supports affine vectors in 4D (last
// component w = 0) when you want to use 4-tuples and 4x4 matrices for affine
// algebra. Thus, you may use those functions for 4-tuples.

import { logAssert } from './Logger.js';
import { Vector, dot, normalize, orthonormalize } from './Vector.js';

function assertSize4(v: Vector): void {
    logAssert(v.size === 4, 'Vector4: vector must have size 4.');
}

// Compute the hypercross product using the formal determinant:
//   hcross = det{{e0,e1,e2,e3},{x0,x1,x2,x3},{y0,y1,y2,y3},{z0,z1,z2,z3}}
// where e0 = (1,0,0,0), e1 = (0,1,0,0), e2 = (0,0,1,0), e3 = (0,0,0,1),
// v0 = (x0,x1,x2,x3), v1 = (y0,y1,y2,y3), and v2 = (z0,z1,z2,z3).
export function hyperCross(v0: Vector, v1: Vector, v2: Vector): Vector {
    assertSize4(v0);
    assertSize4(v1);
    assertSize4(v2);

    const a0 = v0.values, a1 = v1.values, a2 = v2.values;
    const m01 = a0[0] * a1[1] - a0[1] * a1[0];  // x0*y1 - y0*x1
    const m02 = a0[0] * a1[2] - a0[2] * a1[0];  // x0*z1 - z0*x1
    const m03 = a0[0] * a1[3] - a0[3] * a1[0];  // x0*w1 - w0*x1
    const m12 = a0[1] * a1[2] - a0[2] * a1[1];  // y0*z1 - z0*y1
    const m13 = a0[1] * a1[3] - a0[3] * a1[1];  // y0*w1 - w0*y1
    const m23 = a0[2] * a1[3] - a0[3] * a1[2];  // z0*w1 - w0*z1
    return Vector.fromArray([
        +m23 * a2[1] - m13 * a2[2] + m12 * a2[3],  // +m23*y2 - m13*z2 + m12*w2
        -m23 * a2[0] + m03 * a2[2] - m02 * a2[3],  // -m23*x2 + m03*z2 - m02*w2
        +m13 * a2[0] - m03 * a2[1] + m01 * a2[3],  // +m13*x2 - m03*y2 + m01*w2
        -m12 * a2[0] + m02 * a2[1] - m01 * a2[2]   // -m12*x2 + m02*y2 - m01*z2
    ]);
}

// Compute the normalized hypercross product.
export function unitHyperCross(v0: Vector, v1: Vector, v2: Vector,
    robust: boolean = false): Vector {
    const result = hyperCross(v0, v1, v2);
    normalize(result, robust);
    return result;
}

// Compute Dot(HyperCross((x0,x1,x2,x3),(y0,y1,y2,y3),(z0,z1,z2,z3)),
// (w0,w1,w2,w3)), where v0 = (x0,x1,x2,x3), v1 = (y0,y1,y2,y3),
// v2 = (z0,z1,z2,z3), and v3 = (w0,w1,w2,w3).
export function dotHyperCross(v0: Vector, v1: Vector, v2: Vector,
    v3: Vector): number {
    return dot(hyperCross(v0, v1, v2), v3);
}

// Compute a right-handed orthonormal basis for the orthogonal complement of
// the input vectors. The function returns the smallest length of the
// unnormalized vectors computed during the process. If this value is nearly
// zero, it is possible that the inputs are linearly dependent (within
// numerical round-off errors). On input, numInputs must be 1, 2 or 3, and
// v[0] through v[numInputs-1] must be initialized. On output, the vectors
// v[0] through v[3] form an orthonormal set (v is mutated in place; the
// vectors from index numInputs on are assigned).
export function computeOrthogonalComplement4(numInputs: number, v: Vector[],
    robust: boolean = false): number {
    if (numInputs === 1) {
        assertSize4(v[0]);
        const a = v[0].values;
        let maxIndex = 0;
        let maxAbsValue = Math.abs(a[0]);
        for (let i = 1; i < 4; ++i) {
            const absValue = Math.abs(a[i]);
            if (absValue > maxAbsValue) {
                maxIndex = i;
                maxAbsValue = absValue;
            }
        }

        if (maxIndex < 2) {
            v[1] = Vector.fromArray([-a[1], +a[0], 0, 0]);
        } else if (maxIndex === 3) {
            // Generally, you can skip this clause and swap the last two
            // components. However, by swapping 2 and 3 in this case, we allow
            // the function to work properly when the inputs are 3D vectors
            // represented as 4D affine vectors (w = 0).
            v[1] = Vector.fromArray([0, +a[2], -a[1], 0]);
        } else {
            v[1] = Vector.fromArray([0, 0, -a[3], +a[2]]);
        }

        numInputs = 2;
    }

    if (numInputs === 2) {
        assertSize4(v[0]);
        assertSize4(v[1]);
        const a0 = v[0].values, a1 = v[1].values;
        const det = [
            a0[0] * a1[1] - a1[0] * a0[1],
            a0[0] * a1[2] - a1[0] * a0[2],
            a0[0] * a1[3] - a1[0] * a0[3],
            a0[1] * a1[2] - a1[1] * a0[2],
            a0[1] * a1[3] - a1[1] * a0[3],
            a0[2] * a1[3] - a1[2] * a0[3]
        ];

        let maxIndex = 0;
        let maxAbsValue = Math.abs(det[0]);
        for (let i = 1; i < 6; ++i) {
            const absValue = Math.abs(det[i]);
            if (absValue > maxAbsValue) {
                maxIndex = i;
                maxAbsValue = absValue;
            }
        }

        if (maxIndex === 0) {
            v[2] = Vector.fromArray([-det[4], +det[2], 0, -det[0]]);
        } else if (maxIndex <= 2) {
            v[2] = Vector.fromArray([+det[5], 0, -det[2], +det[1]]);
        } else {
            v[2] = Vector.fromArray([0, -det[5], +det[4], -det[3]]);
        }

        numInputs = 3;
    }

    if (numInputs === 3) {
        v[3] = hyperCross(v[0], v[1], v[2]);
        return orthonormalize(4, v, robust);
    }

    return 0;
}
