// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumAreaCircle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the minimum area circle containing the input set of points. The
// algorithm randomly permutes the input points so that the construction
// occurs in 'expected' O(N) time. All internal minimal circle calculations
// store the squared radius in the radius member of Circle2. Only at the end
// is a sqrt computed.
//
// Upstream templates the class on a ComputeType that is intended to be
// BSRational<T> for exact rational arithmetic. With exact arithmetic the
// result is guaranteed to be the minimum-area circle. This port instantiates
// the floating-point compute type (number), matching the upstream
// double-precision instantiation, so floating-point rounding errors can
// cause the updateSupport{2,3} functions to fail. The failure is trapped in
// those functions and a simple bounding circle is computed using
// getContainerCircle2 from ContCircle2.ts. That circle is generally not the
// minimum-area circle containing the points, and the minimum-area algorithm
// is terminated immediately. The circle is returned along with a boolean
// 'success' that is true when the circle is minimum area and false when the
// failure is trapped. When false is returned, another call to compute(...)
// can be tried: the random shuffle uses the persistent engine state, so it
// is highly likely to differ from the previous shuffle and there is a chance
// the algorithm succeeds just because of the different ordering of points.
//
// Port notes:
// * operator() becomes compute(points), returning { minimal, success }
//   instead of writing through an output reference and returning bool.
// * GetNumSupport/GetSupport become the numSupport/support accessors.
// * std::default_random_engine becomes a module-private minstd_rand0-style
//   Lehmer generator seeded with 1, so the shuffle is deterministic for a
//   given call sequence in the port (as it is upstream, though the exact
//   permutation is implementation-defined in C++ and therefore differs).
// * The std::sort comparator used to find unique points is augmented with an
//   index tie-break so the deduplication is fully deterministic (std::sort is
//   not stable, so upstream's choice of representative among equal points is
//   unspecified).

import { getContainerCircle2 } from './ContCircle2.js';
import { Hypersphere, type Circle2 } from './Hypersphere.js';
import { LinearSystem } from './LinearSystem.js';
import { logAssert, logError } from './Logger.js';
import { Matrix } from './Matrix.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// The result of the minimum-area-circle query. 'success' is upstream's
// boolean return value: true when 'minimal' is the minimum-area circle,
// false when floating-point round-off caused a trapped failure and 'minimal'
// is the (generally larger) circle from getContainerCircle2.
export interface MinimumAreaCircle2Result {
    minimal: Circle2;
    success: boolean;
}

// The multiplicative congruential generator x -> 16807 * x mod (2^31 - 1)
// with seed 1 (std::minstd_rand0, the common choice for
// std::default_random_engine). The products stay below 2^45, so the
// arithmetic is exact in IEEE doubles.
class DefaultRandomEngine {
    private mState = 1;

    // Returns a pseudo-random number in [1, 2^31 - 2].
    next(): number {
        this.mState = (16807 * this.mState) % 2147483647;
        return this.mState;
    }
}

// Randomly permute the array in place (the port of std::shuffle over the
// whole array): Fisher-Yates using the engine to draw an index in [0, i].
function shuffle(values: number[], engine: DefaultRandomEngine): void {
    for (let i = values.length - 1; i > 0; --i) {
        const j = engine.next() % (i + 1);
        const temp = values[i];
        values[i] = values[j];
        values[j] = temp;
    }
}

// The internal pairing of a candidate circle with the flag that reports
// whether the update succeeded (upstream's std::pair<Circle2,bool>).
interface UpdateResult {
    circle: Circle2;
    valid: boolean;
}

export class MinimumAreaCircle2 {
    private mNumSupport: number;
    private mSupport: number[];
    private mDRE: DefaultRandomEngine;
    private mComputePoints: Vector[];

    constructor() {
        this.mNumSupport = 0;
        this.mSupport = [0, 0, 0];
        this.mDRE = new DefaultRandomEngine();
        this.mComputePoints = [];
    }

    compute(points: readonly Vector[]): MinimumAreaCircle2Result {
        logAssert(points.length >= 1, 'Input must contain points.');
        for (const point of points) {
            logAssert(point.size === 2,
                'MinimumAreaCircle2: points must be 2D.');
        }

        // Function array to avoid a switch statement in the main loop.
        const update: Array<(i: number) => UpdateResult> = [
            () => logError('unreachable'),
            (i: number) => this.updateSupport1(i),
            (i: number) => this.updateSupport2(i),
            (i: number) => this.updateSupport3(i)
        ];

        // Process only the unique points.
        let permuted: number[] = [];
        for (let i = 0; i < points.length; ++i) {
            permuted.push(i);
        }
        permuted.sort((i0, i1) => {
            if (points[i0].lessThan(points[i1])) {
                return -1;
            }
            if (points[i1].lessThan(points[i0])) {
                return 1;
            }
            return i0 - i1;
        });
        const unique: number[] = [];
        for (const index of permuted) {
            if (unique.length === 0
                || !points[unique[unique.length - 1]].equals(points[index])) {
                unique.push(index);
            }
        }
        permuted = unique;
        const numPoints = permuted.length;

        // Create a random permutation of the points.
        shuffle(permuted, this.mDRE);

        // Convert to the compute type, which is a simple copy when the
        // compute type is the same as the input type.
        this.mComputePoints = new Array<Vector>(numPoints);
        for (let i = 0; i < numPoints; ++i) {
            this.mComputePoints[i] = points[permuted[i]].clone();
        }

        // Start with the first point.
        let ctMinimal = this.exactCircle1(0);
        this.mNumSupport = 1;
        this.mSupport[0] = 0;

        // The loop restarts from the beginning of the point list each time
        // the circle needs updating. Linus Kallberg (Computer Science at
        // Malardalen University in Sweden) discovered that performance is
        // better when the remaining points in the array are processed before
        // restarting. The points processed before the point that caused the
        // update are likely to be enclosed by the new circle (or near the
        // circle boundary) because they were enclosed by the previous circle.
        // The chances are better that points after the current one will cause
        // growth of the bounding circle.
        for (let i = 1 % numPoints, n = 0; i !== n; i = (i + 1) % numPoints) {
            if (!this.supportContains(i)) {
                if (!this.contains(i, ctMinimal)) {
                    const result = update[this.mNumSupport](i);
                    if (result.valid) {
                        if (result.circle.radius > ctMinimal.radius) {
                            ctMinimal = result.circle;
                            n = i;
                        }
                    } else {
                        // This case can happen because the compute type is
                        // floating point; see the comments at the beginning
                        // of this file. Return a non-minimal circle.
                        //
                        // Upstream bug (MinimumAreaCircle2.h, the trapped
                        // failure branch of operator()): the call is
                        // GetContainer(numPoints, points, minimal), but
                        // numPoints has been overwritten with the number of
                        // *unique* points while points[] is still the full
                        // input array. The bounding circle is therefore
                        // computed from only the first numPoints entries of
                        // the input and need not contain all of the points.
                        // The port fixes this by using all input points.
                        this.mNumSupport = 0;
                        this.mSupport = [0, 0, 0];
                        return {
                            minimal: getContainerCircle2(points),
                            success: false
                        };
                    }
                }
            }
        }

        const minimal = new Hypersphere(2);
        minimal.center = ctMinimal.center.clone();
        minimal.radius = Math.sqrt(ctMinimal.radius);

        for (let i = 0; i < this.mNumSupport; ++i) {
            this.mSupport[i] = permuted[this.mSupport[i]];
        }
        return { minimal, success: true };
    }

    // Member access. The support indices are lookups into the points[] array
    // passed to compute(...). Only the first numSupport entries are valid.
    get numSupport(): number {
        return this.mNumSupport;
    }

    get support(): readonly number[] {
        return this.mSupport;
    }

    // Test whether point P is inside circle C using squared distance and
    // squared radius. In this algorithm, circle.radius is the *squared
    // radius* until compute(...) returns, at which time a square root is
    // applied.
    private contains(i: number, circle: Circle2): boolean {
        const diff = sub(this.mComputePoints[i], circle.center);
        return dot(diff, diff) <= circle.radius;
    }

    private exactCircle1(i0: number): Circle2 {
        const minimal = new Hypersphere(2);
        minimal.center = this.mComputePoints[i0].clone();
        minimal.radius = 0;
        return minimal;
    }

    private exactCircle2(i0: number, i1: number): Circle2 {
        const P0 = this.mComputePoints[i0];
        const P1 = this.mComputePoints[i1];
        const diff = sub(P1, P0);
        const minimal = new Hypersphere(2);
        minimal.center = mul(0.5, add(P0, P1));
        minimal.radius = 0.25 * dot(diff, diff);
        return minimal;
    }

    private exactCircle3(i0: number, i1: number, i2: number): Circle2 {
        // Compute the 2D circle containing P0, P1, and P2. The center in
        // barycentric coordinates is C = x0*P0 + x1*P1 + x2*P2, where
        // x0 + x1 + x2 = 1. The center is equidistant from the three points,
        // so |C - P0| = |C - P1| = |C - P2| = R, where R is the radius of the
        // circle. From these conditions,
        //   C - P0 = x0*E0 + x1*E1 - E0
        //   C - P1 = x0*E0 + x1*E1 - E1
        //   C - P2 = x0*E0 + x1*E1
        // where E0 = P0 - P2 and E1 = P1 - P2, which leads to
        //   r^2 = |x0*E0 + x1*E1|^2 - 2*Dot(E0, x0*E0 + x1*E1) + |E0|^2
        //   r^2 = |x0*E0 + x1*E1|^2 - 2*Dot(E1, x0*E0 + x1*E1) + |E1|^2
        //   r^2 = |x0*E0 + x1*E1|^2
        // Subtracting the last equation from the first two and writing the
        // equations as a linear system,
        //
        // +-                     -++   -+       +-          -+
        // | Dot(E0,E0) Dot(E0,E1) || x0 | = 0.5 | Dot(E0,E0) |
        // | Dot(E1,E0) Dot(E1,E1) || x1 |       | Dot(E1,E1) |
        // +-                     -++   -+       +-          -+
        //
        // The following code solves this system for x0 and x1 and then
        // evaluates the third equation in r^2 to obtain r.
        const P0 = this.mComputePoints[i0];
        const P1 = this.mComputePoints[i1];
        const P2 = this.mComputePoints[i2];

        const E0 = sub(P0, P2);
        const E1 = sub(P1, P2);

        const A = new Matrix(2, 2);
        A.set(0, 0, dot(E0, E0));
        A.set(0, 1, dot(E0, E1));
        A.set(1, 0, A.get(0, 1));
        A.set(1, 1, dot(E1, E1));

        const half = 0.5;
        const B = Vector.fromArray([half * A.get(0, 0), half * A.get(1, 1)]);

        const minimal = new Hypersphere(2);
        const { X, invertible } = LinearSystem.solve2x2(A, B);
        if (invertible) {
            const x2 = 1 - X.get(0) - X.get(1);
            minimal.center = add(add(mul(X.get(0), P0), mul(X.get(1), P1)),
                mul(x2, P2));
            const tmp = add(mul(X.get(0), E0), mul(X.get(1), E1));
            minimal.radius = dot(tmp, tmp);
        } else {
            minimal.center = Vector.zero(2);
            minimal.radius = Number.MAX_VALUE;
        }

        return minimal;
    }

    private updateSupport1(i: number): UpdateResult {
        const minimal = this.exactCircle2(this.mSupport[0], i);
        this.mNumSupport = 2;
        this.mSupport[1] = i;
        return { circle: minimal, valid: true };
    }

    private updateSupport2(i: number): UpdateResult {
        // Permutations of type 2, used for calling exactCircle2(...).
        const type2: readonly (readonly number[])[] = [
            [0, /*2*/ 1],
            [1, /*2*/ 0]
        ];

        // Permutations of type 3, used for calling exactCircle3(...): the
        // single permutation {0, 1, 2}.

        const circle: Circle2[] = [];
        let minRSqr = Number.MAX_VALUE;
        let iCircle = 0;
        let iMinRSqr = -1;
        let k0: number;
        let k1: number;

        // Permutations of type 2.
        for (let j = 0; j < type2.length; ++j, ++iCircle) {
            k0 = this.mSupport[type2[j][0]];
            circle[iCircle] = this.exactCircle2(k0, i);
            if (circle[iCircle].radius < minRSqr) {
                k1 = this.mSupport[type2[j][1]];
                if (this.contains(k1, circle[iCircle])) {
                    minRSqr = circle[iCircle].radius;
                    iMinRSqr = iCircle;
                }
            }
        }

        // Permutations of type 3.
        k0 = this.mSupport[0];
        k1 = this.mSupport[1];
        circle[iCircle] = this.exactCircle3(k0, k1, i);
        if (circle[iCircle].radius < minRSqr) {
            minRSqr = circle[iCircle].radius;
            iMinRSqr = iCircle;
        }

        switch (iMinRSqr) {
            case 0:
                this.mSupport[1] = i;
                break;
            case 1:
                this.mSupport[0] = i;
                break;
            case 2:
                this.mNumSupport = 3;
                this.mSupport[2] = i;
                break;
            default:
                // For exact arithmetic, iMinRSqr >= 0, but for floating-point
                // arithmetic round-off errors can lead to iMinRSqr == -1.
                // When this happens, use a simple bounding circle for the
                // result and terminate the minimum-area algorithm.
                return { circle: new Hypersphere(2), valid: false };
        }

        return { circle: circle[iMinRSqr], valid: true };
    }

    private updateSupport3(i: number): UpdateResult {
        // Permutations of type 2, used for calling exactCircle2(...).
        const type2: readonly (readonly number[])[] = [
            [0, /*3*/ 1, 2],
            [1, /*3*/ 0, 2],
            [2, /*3*/ 0, 1]
        ];

        // Permutations of type 3, used for calling exactCircle3(...).
        const type3: readonly (readonly number[])[] = [
            [0, 1, /*3*/ 2],
            [0, 2, /*3*/ 1],
            [1, 2, /*3*/ 0]
        ];

        const circle: Circle2[] = [];
        let minRSqr = Number.MAX_VALUE;
        let iCircle = 0;
        let iMinRSqr = -1;
        let k0: number;
        let k1: number;
        let k2: number;

        // Permutations of type 2.
        for (let j = 0; j < type2.length; ++j, ++iCircle) {
            k0 = this.mSupport[type2[j][0]];
            circle[iCircle] = this.exactCircle2(k0, i);
            if (circle[iCircle].radius < minRSqr) {
                k1 = this.mSupport[type2[j][1]];
                k2 = this.mSupport[type2[j][2]];
                if (this.contains(k1, circle[iCircle])
                    && this.contains(k2, circle[iCircle])) {
                    minRSqr = circle[iCircle].radius;
                    iMinRSqr = iCircle;
                }
            }
        }

        // Permutations of type 3.
        for (let j = 0; j < type3.length; ++j, ++iCircle) {
            k0 = this.mSupport[type3[j][0]];
            k1 = this.mSupport[type3[j][1]];
            circle[iCircle] = this.exactCircle3(k0, k1, i);
            if (circle[iCircle].radius < minRSqr) {
                k2 = this.mSupport[type3[j][2]];
                if (this.contains(k2, circle[iCircle])) {
                    minRSqr = circle[iCircle].radius;
                    iMinRSqr = iCircle;
                }
            }
        }

        switch (iMinRSqr) {
            case 0:
                this.mNumSupport = 2;
                this.mSupport[1] = i;
                break;
            case 1:
                this.mNumSupport = 2;
                this.mSupport[0] = i;
                break;
            case 2:
                this.mNumSupport = 2;
                this.mSupport[0] = this.mSupport[2];
                this.mSupport[1] = i;
                break;
            case 3:
                this.mSupport[2] = i;
                break;
            case 4:
                this.mSupport[1] = i;
                break;
            case 5:
                this.mSupport[0] = i;
                break;
            default:
                // See the comment in updateSupport2 for iMinRSqr == -1.
                return { circle: new Hypersphere(2), valid: false };
        }

        return { circle: circle[iMinRSqr], valid: true };
    }

    // Indices of points that support the current minimum area circle.
    private supportContains(j: number): boolean {
        for (let i = 0; i < this.mNumSupport; ++i) {
            if (j === this.mSupport[i]) {
                return true;
            }
        }
        return false;
    }
}
