// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Vector3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'Vector3<Real>' is the alias 'Vector<3, Real>'; the
// port uses the runtime-sized Vector with size 3 (create with
// 'Vector.fromArray([x, y, z])'). The compile-time dimension checks become
// runtime asserts. Upstream ComputeOrthogonalComplement and
// ComputeBarycentrics collide with the 2D versions in Vector2.h under the
// library-wide flat export, so the port suffixes them with the dimension:
// computeOrthogonalComplement3 and computeBarycentrics3. The output
// parameters 'bary' of ComputeBarycentrics and 'v0'/'v1' of
// FastComputeOrthogonalComplement become fields of the returned objects.

import { logAssert } from './Logger.js';
import {
    Vector, dot, length, mul, normalize, orthonormalize, sub
} from './Vector.js';

function assertSize3(v: Vector): void {
    logAssert(v.size === 3, 'Vector3: vector must have size 3.');
}

// Cross, UnitCross, and DotCross accept vectors of size 3 or 4 (the port of
// the upstream template parameter N). The latter case supports affine
// vectors in 4D (last component w = 0) when you want to use 4-tuples and
// 4x4 matrices for affine algebra.

// Compute the cross product using the formal determinant:
//   cross = det{{e0,e1,e2},{x0,x1,x2},{y0,y1,y2}}
//         = (x1*y2-x2*y1, x2*y0-x0*y2, x0*y1-x1*y0)
// where e0 = (1,0,0), e1 = (0,1,0), e2 = (0,0,1), v0 = (x0,x1,x2), and
// v1 = (y0,y1,y2).
export function cross(v0: Vector, v1: Vector): Vector {
    logAssert(v0.size === 3 || v0.size === 4, 'Dimension must be 3 or 4.');
    logAssert(v0.size === v1.size, 'Vector3: mismatched sizes.');

    const result = new Vector(v0.size);
    result.values[0] = v0.values[1] * v1.values[2] - v0.values[2] * v1.values[1];
    result.values[1] = v0.values[2] * v1.values[0] - v0.values[0] * v1.values[2];
    result.values[2] = v0.values[0] * v1.values[1] - v0.values[1] * v1.values[0];
    return result;
}

// Compute the normalized cross product.
export function unitCross(v0: Vector, v1: Vector,
    robust: boolean = false): Vector {
    const result = cross(v0, v1);
    normalize(result, robust);
    return result;
}

// Compute Dot((x0,x1,x2),Cross((y0,y1,y2),(z0,z1,z2)), the triple scalar
// product of three vectors, where v0 = (x0,x1,x2), v1 = (y0,y1,y2), and
// v2 is (z0,z1,z2).
export function dotCross(v0: Vector, v1: Vector, v2: Vector): number {
    return dot(v0, cross(v1, v2));
}

// Compute a right-handed orthonormal basis for the orthogonal complement
// of the input vectors. The function returns the smallest length of the
// unnormalized vectors computed during the process. If this value is
// nearly zero, it is possible that the inputs are linearly dependent
// (within numerical round-off errors). On input, numInputs must be 1 or
// 2 and v[0] through v[numInputs-1] must be initialized. On output, the
// vectors v[0] through v[2] form an orthonormal set (v is mutated in
// place; v[1] and/or v[2] are assigned).
export function computeOrthogonalComplement3(numInputs: number, v: Vector[],
    robust: boolean = false): number {
    if (numInputs === 1) {
        assertSize3(v[0]);
        if (Math.abs(v[0].values[0]) > Math.abs(v[0].values[1])) {
            v[1] = Vector.fromArray([-v[0].values[2], 0, +v[0].values[0]]);
        } else {
            v[1] = Vector.fromArray([0, +v[0].values[2], -v[0].values[1]]);
        }
        numInputs = 2;
    }

    if (numInputs === 2) {
        v[2] = cross(v[0], v[1]);
        return orthonormalize(3, v, robust);
    }

    return 0;
}

// Compute a right-handed orthonormal basis {v0,v1,v2} for the orthogonal
// complement of a unit-length vector v2. See
// https://www.geometrictools.com/Documentation/FastOrthogonalComplement.pdf
// The upstream output parameters v0 and v1 are returned as an object.
export function fastComputeOrthogonalComplement(v2: Vector):
    { v0: Vector, v1: Vector } {
    assertSize3(v2);
    let temp0: number, temp1: number, temp2: number;
    let v0: Vector, v1: Vector;
    if (v2.values[2] >= 0) {
        temp0 = 1 + v2.values[2];
        temp1 = -v2.values[0] * v2.values[1] / temp0;
        temp2 = v2.values[1] * v2.values[1] / temp0;
        v0 = Vector.fromArray([v2.values[2] + temp2, temp1, -v2.values[0]]);
        v1 = Vector.fromArray([temp1, 1 - temp2, -v2.values[1]]);
    } else {
        temp0 = 1 - v2.values[2];
        temp1 = v2.values[0] * v2.values[1] / temp0;
        temp2 = v2.values[1] * v2.values[1] / temp0;
        v0 = Vector.fromArray([-v2.values[2] + temp2, -temp1, v2.values[0]]);
        v1 = Vector.fromArray([temp1, -1 + temp2, -v2.values[1]]);
    }
    return { v0, v1 };
}

// Compute the barycentric coordinates of the point P with respect to the
// tetrahedron <V0,V1,V2,V3>, P = b0*V0 + b1*V1 + b2*V2 + b3*V3, where
// b0 + b1 + b2 + b3 = 1. The 'valid' field is 'true' iff {V0,V1,V2,V3} is
// a linearly independent set. Numerically, this is measured by
// |det[V0 V1 V2 V3]| <= epsilon. The values bary[] are valid only when
// 'valid' is 'true' but set to zero when it is 'false'.
export function computeBarycentrics3(p: Vector, v0: Vector, v1: Vector,
    v2: Vector, v3: Vector, epsilon: number = 0):
    { valid: boolean, bary: [number, number, number, number] } {
    // Compute the vectors relative to V3 of the tetrahedron.
    const diff: [Vector, Vector, Vector, Vector] =
        [sub(v0, v3), sub(v1, v3), sub(v2, v3), sub(p, v3)];

    const det = dotCross(diff[0], diff[1], diff[2]);
    if (det < -epsilon || det > epsilon) {
        const bary: [number, number, number, number] = [0, 0, 0, 0];
        bary[0] = dotCross(diff[3], diff[1], diff[2]) / det;
        bary[1] = dotCross(diff[3], diff[2], diff[0]) / det;
        bary[2] = dotCross(diff[3], diff[0], diff[1]) / det;
        bary[3] = 1 - bary[0] - bary[1] - bary[2];
        return { valid: true, bary };
    }

    return { valid: false, bary: [0, 0, 0, 0] };
}

// Get intrinsic information about the input array of vectors. The members
// are valid iff the inputs are valid (vectors is nonempty and epsilon >= 0);
// otherwise they keep their zero defaults, as upstream. (The upstream
// numVectors/pointer pair becomes the array itself.)
export class IntrinsicsVector3 {
    // A nonnegative tolerance that is used to determine the intrinsic
    // dimension of the set.
    epsilon: number;

    // The intrinsic dimension of the input set, computed based on the
    // nonnegative tolerance epsilon.
    dimension: number;

    // Axis-aligned bounding box of the input set. The maximum range is
    // the larger of max[0]-min[0], max[1]-min[1], and max[2]-min[2].
    min: [number, number, number];
    max: [number, number, number];
    maxRange: number;

    // Coordinate system. The origin is valid for any dimension d. The
    // unit-length direction vector is valid only for 0 <= i < d. The
    // extreme index is relative to the array of input points, and is also
    // valid only for 0 <= i < d. If d = 0, all points are effectively
    // the same, but the use of an epsilon may lead to an extreme index
    // that is not zero. If d = 1, all points effectively lie on a line
    // segment. If d = 2, all points effectively lie on a plane. If
    // d = 3, the points are not coplanar.
    origin: Vector;
    direction: [Vector, Vector, Vector];

    // The indices that define the maximum dimensional extents. The
    // values extreme[0] and extreme[1] are the indices for the points
    // that define the largest extent in one of the coordinate axis
    // directions. If the dimension is 2, then extreme[2] is the index
    // for the point that generates the largest extent in the direction
    // perpendicular to the line through the points corresponding to
    // extreme[0] and extreme[1]. Furthermore, if the dimension is 3,
    // then extreme[3] is the index for the point that generates the
    // largest extent in the direction perpendicular to the triangle
    // defined by the other extreme points. The tetrahedron formed by the
    // points V[extreme[0]], V[extreme[1]], V[extreme[2]], and
    // V[extreme[3]] is clockwise or counterclockwise, the condition
    // stored in extremeCCW.
    extreme: [number, number, number, number];
    extremeCCW: boolean;

    // The constructor sets the class members based on the input set.
    constructor(v: readonly Vector[], inEpsilon: number) {
        this.epsilon = inEpsilon;
        this.dimension = 0;
        this.min = [0, 0, 0];
        this.max = [0, 0, 0];
        this.maxRange = 0;
        this.origin = new Vector(3);
        this.direction = [new Vector(3), new Vector(3), new Vector(3)];
        this.extreme = [0, 0, 0, 0];
        this.extremeCCW = false;

        const numVectors = v.length;
        if (numVectors > 0 && this.epsilon >= 0) {
            // Compute the axis-aligned bounding box for the input vectors.
            // Keep track of the indices into 'vectors' for the current
            // min and max.
            let j: number;
            const indexMin: [number, number, number] = [0, 0, 0];
            const indexMax: [number, number, number] = [0, 0, 0];
            for (j = 0; j < 3; ++j) {
                assertSize3(v[0]);
                this.min[j] = v[0].values[j];
                this.max[j] = this.min[j];
                indexMin[j] = 0;
                indexMax[j] = 0;
            }

            let i: number;
            for (i = 1; i < numVectors; ++i) {
                assertSize3(v[i]);
                for (j = 0; j < 3; ++j) {
                    if (v[i].values[j] < this.min[j]) {
                        this.min[j] = v[i].values[j];
                        indexMin[j] = i;
                    } else if (v[i].values[j] > this.max[j]) {
                        this.max[j] = v[i].values[j];
                        indexMax[j] = i;
                    }
                }
            }

            // Determine the maximum range for the bounding box.
            this.maxRange = this.max[0] - this.min[0];
            this.extreme[0] = indexMin[0];
            this.extreme[1] = indexMax[0];
            let range = this.max[1] - this.min[1];
            if (range > this.maxRange) {
                this.maxRange = range;
                this.extreme[0] = indexMin[1];
                this.extreme[1] = indexMax[1];
            }
            range = this.max[2] - this.min[2];
            if (range > this.maxRange) {
                this.maxRange = range;
                this.extreme[0] = indexMin[2];
                this.extreme[1] = indexMax[2];
            }

            // The origin is either the vector of minimum x0-value, vector
            // of minimum x1-value, or vector of minimum x2-value.
            this.origin = v[this.extreme[0]].clone();

            // Test whether the vector set is (nearly) a vector.
            if (this.maxRange <= this.epsilon) {
                this.dimension = 0;
                for (j = 0; j < 3; ++j) {
                    this.extreme[j + 1] = this.extreme[0];
                }
                return;
            }

            // Test whether the vector set is (nearly) a line segment. We
            // need {direction[1],direction[2]} to span the orthogonal
            // complement of direction[0].
            this.direction[0] = sub(v[this.extreme[1]], this.origin);
            normalize(this.direction[0], false);
            if (Math.abs(this.direction[0].values[0])
                > Math.abs(this.direction[0].values[1])) {
                this.direction[1] = Vector.fromArray([
                    -this.direction[0].values[2],
                    0,
                    +this.direction[0].values[0]]);
            } else {
                this.direction[1] = Vector.fromArray([
                    0,
                    +this.direction[0].values[2],
                    -this.direction[0].values[1]]);
            }
            normalize(this.direction[1], false);
            this.direction[2] = cross(this.direction[0], this.direction[1]);

            // Compute the maximum distance of the points from the line
            // origin + t * direction[0].
            let maxDistance = 0;
            let distance: number;
            let d: number;
            this.extreme[2] = this.extreme[0];
            for (i = 0; i < numVectors; ++i) {
                const diff = sub(v[i], this.origin);
                d = dot(this.direction[0], diff);
                const proj = sub(diff, mul(this.direction[0], d));
                distance = length(proj, false);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    this.extreme[2] = i;
                }
            }

            if (maxDistance <= this.epsilon * this.maxRange) {
                // The points are (nearly) on the line
                // origin + t * direction[0].
                this.dimension = 1;
                this.extreme[2] = this.extreme[1];
                this.extreme[3] = this.extreme[1];
                return;
            }

            // Test whether the vector set is (nearly) a planar polygon.
            // The point v[extreme[2]] is farthest from the line:
            // origin + t * direction[0]. The vector
            // v[extreme[2]] - origin is not necessarily perpendicular to
            // direction[0], so project out the direction[0] component so
            // that the result is perpendicular to direction[0].
            this.direction[1] = sub(v[this.extreme[2]], this.origin);
            d = dot(this.direction[0], this.direction[1]);
            this.direction[1] = sub(this.direction[1],
                mul(this.direction[0], d));
            normalize(this.direction[1], false);

            // We need direction[2] to span the orthogonal complement of
            // {direction[0],direction[1]}.
            this.direction[2] = cross(this.direction[0], this.direction[1]);

            // Compute the maximum distance of the points from the plane
            // origin + t0 * direction[0] + t1 * direction[1].
            maxDistance = 0;
            let maxSign = 0;
            this.extreme[3] = this.extreme[0];
            for (i = 0; i < numVectors; ++i) {
                const diff = sub(v[i], this.origin);
                distance = dot(this.direction[2], diff);
                const sign = (distance > 0 ? 1 : (distance < 0 ? -1 : 0));
                distance = Math.abs(distance);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    maxSign = sign;
                    this.extreme[3] = i;
                }
            }

            if (maxDistance <= this.epsilon * this.maxRange) {
                // The points are (nearly) on the plane
                // origin + t0 * direction[0] + t1 * direction[1].
                this.dimension = 2;
                this.extreme[3] = this.extreme[2];
                return;
            }

            this.dimension = 3;
            this.extremeCCW = (maxSign > 0);
            return;
        }
    }
}
