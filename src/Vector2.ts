// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Vector2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'Vector2<Real>' is the alias 'Vector<2, Real>'; the
// port uses the runtime-sized Vector with size 2 (create with
// 'Vector.fromArray([x, y])'). The compile-time dimension checks become
// runtime asserts. Upstream ComputeOrthogonalComplement and
// ComputeBarycentrics collide with the 3D versions in Vector3.h under the
// library-wide flat export, so the port suffixes them with the dimension:
// computeOrthogonalComplement2 and computeBarycentrics2. The output
// parameter 'bary' of ComputeBarycentrics becomes a field of the returned
// object.

import { logAssert } from './Logger';
import { Vector, dot, negate, normalize, orthonormalize, sub } from './Vector';

function assertSize2(v: Vector): void {
    logAssert(v.size === 2, 'Vector2: vector must have size 2.');
}

// Compute the perpendicular using the formal determinant,
//   perp = det{{e0,e1},{x0,x1}} = (x1,-x0)
// where e0 = (1,0), e1 = (0,1), and v = (x0,x1).
export function perp(v: Vector): Vector {
    assertSize2(v);
    return Vector.fromArray([v.values[1], -v.values[0]]);
}

// Compute the normalized perpendicular.
export function unitPerp(v: Vector, robust: boolean = false): Vector {
    const result = perp(v);
    normalize(result, robust);
    return result;
}

// Compute Dot((x0,x1),Perp(y0,y1)) = x0*y1 - x1*y0, where v0 = (x0,x1)
// and v1 = (y0,y1).
export function dotPerp(v0: Vector, v1: Vector): number {
    return dot(v0, perp(v1));
}

// Compute a right-handed orthonormal basis for the orthogonal complement
// of the input vectors. The function returns the smallest length of the
// unnormalized vectors computed during the process. If this value is
// nearly zero, it is possible that the inputs are linearly dependent
// (within numerical round-off errors). On input, numInputs must be 1 and
// v[0] must be initialized. On output, the vectors v[0] and v[1] form an
// orthonormal set (v is mutated in place; v[1] is assigned).
export function computeOrthogonalComplement2(numInputs: number, v: Vector[],
    robust: boolean = false): number {
    if (numInputs === 1) {
        v[1] = negate(perp(v[0]));
        return orthonormalize(2, v, robust);
    }

    return 0;
}

// Compute the barycentric coordinates of the point P with respect to the
// triangle <V0,V1,V2>, P = b0*V0 + b1*V1 + b2*V2, where b0 + b1 + b2 = 1.
// The 'valid' field is 'true' iff {V0,V1,V2} is a linearly independent set.
// Numerically, this is measured by |det[V0 V1 V2]| <= epsilon. The values
// bary[] are valid only when 'valid' is 'true' but set to zero when it is
// 'false'.
export function computeBarycentrics2(p: Vector, v0: Vector, v1: Vector,
    v2: Vector, epsilon: number = 0):
    { valid: boolean, bary: [number, number, number] } {
    // Compute the vectors relative to V2 of the triangle.
    const diff: [Vector, Vector, Vector] =
        [sub(v0, v2), sub(v1, v2), sub(p, v2)];

    const det = dotPerp(diff[0], diff[1]);
    if (det < -epsilon || det > epsilon) {
        const bary: [number, number, number] = [0, 0, 0];
        bary[0] = dotPerp(diff[2], diff[1]) / det;
        bary[1] = dotPerp(diff[0], diff[2]) / det;
        bary[2] = 1 - bary[0] - bary[1];
        return { valid: true, bary };
    }

    return { valid: false, bary: [0, 0, 0] };
}

// Get intrinsic information about the input array of vectors. The members
// are valid iff the inputs are valid (vectors is nonempty and epsilon >= 0);
// otherwise they keep their zero defaults, as upstream. (The upstream
// numVectors/pointer pair becomes the array itself.)
export class IntrinsicsVector2 {
    // A nonnegative tolerance that is used to determine the intrinsic
    // dimension of the set.
    epsilon: number;

    // The intrinsic dimension of the input set, computed based on the
    // nonnegative tolerance epsilon.
    dimension: number;

    // Axis-aligned bounding box of the input set. The maximum range is
    // the larger of max[0]-min[0] and max[1]-min[1].
    min: [number, number];
    max: [number, number];
    maxRange: number;

    // Coordinate system. The origin is valid for any dimension d. The
    // unit-length direction vector is valid only for 0 <= i < d. The
    // extreme index is relative to the array of input points, and is also
    // valid only for 0 <= i < d. If d = 0, all points are effectively
    // the same, but the use of an epsilon may lead to an extreme index
    // that is not zero. If d = 1, all points effectively lie on a line
    // segment. If d = 2, the points are not collinear.
    origin: Vector;
    direction: [Vector, Vector];

    // The indices that define the maximum dimensional extents. The
    // values extreme[0] and extreme[1] are the indices for the points
    // that define the largest extent in one of the coordinate axis
    // directions. If the dimension is 2, then extreme[2] is the index
    // for the point that generates the largest extent in the direction
    // perpendicular to the line through the points corresponding to
    // extreme[0] and extreme[1]. The triangle formed by the points
    // V[extreme[0]], V[extreme[1]], and V[extreme[2]] is clockwise or
    // counterclockwise, the condition stored in extremeCCW.
    extreme: [number, number, number];
    extremeCCW: boolean;

    // The constructor sets the class members based on the input set.
    constructor(v: readonly Vector[], inEpsilon: number) {
        this.epsilon = inEpsilon;
        this.dimension = 0;
        this.min = [0, 0];
        this.max = [0, 0];
        this.maxRange = 0;
        this.origin = new Vector(2);
        this.direction = [new Vector(2), new Vector(2)];
        this.extreme = [0, 0, 0];
        this.extremeCCW = false;

        const numVectors = v.length;
        if (numVectors > 0 && this.epsilon >= 0) {
            // Compute the axis-aligned bounding box for the input
            // vectors. Keep track of the indices into 'vectors' for the
            // current min and max.
            let j: number;
            const indexMin: [number, number] = [0, 0];
            const indexMax: [number, number] = [0, 0];
            for (j = 0; j < 2; ++j) {
                assertSize2(v[0]);
                this.min[j] = v[0].values[j];
                this.max[j] = this.min[j];
                indexMin[j] = 0;
                indexMax[j] = 0;
            }

            let i: number;
            for (i = 1; i < numVectors; ++i) {
                assertSize2(v[i]);
                for (j = 0; j < 2; ++j) {
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
            const range = this.max[1] - this.min[1];
            if (range > this.maxRange) {
                this.maxRange = range;
                this.extreme[0] = indexMin[1];
                this.extreme[1] = indexMax[1];
            }

            // The origin is either the vector of minimum x0-value or
            // vector of minimum x1-value.
            this.origin = v[this.extreme[0]].clone();

            // Test whether the vector set is (nearly) a vector.
            if (this.maxRange <= this.epsilon) {
                this.dimension = 0;
                for (j = 0; j < 2; ++j) {
                    this.extreme[j + 1] = this.extreme[0];
                }
                return;
            }

            // Test whether the vector set is (nearly) a line segment. We
            // need direction[1] to span the orthogonal complement of
            // direction[0].
            this.direction[0] = sub(v[this.extreme[1]], this.origin);
            normalize(this.direction[0], false);
            this.direction[1] = negate(perp(this.direction[0]));

            // Compute the maximum distance of the points from the line
            // origin + t * direction[0].
            let maxDistance = 0;
            let maxSign = 0;
            this.extreme[2] = this.extreme[0];
            for (i = 0; i < numVectors; ++i) {
                const diff = sub(v[i], this.origin);
                let distance = dot(this.direction[1], diff);
                const sign = (distance > 0 ? 1 : (distance < 0 ? -1 : 0));
                distance = Math.abs(distance);
                if (distance > maxDistance) {
                    maxDistance = distance;
                    maxSign = sign;
                    this.extreme[2] = i;
                }
            }

            if (maxDistance <= this.epsilon * this.maxRange) {
                // The points are (nearly) on the line
                // origin + t * direction[0].
                this.dimension = 1;
                this.extreme[2] = this.extreme[1];
                return;
            }

            this.dimension = 2;
            this.extremeCCW = (maxSign > 0);
            return;
        }
    }
}
