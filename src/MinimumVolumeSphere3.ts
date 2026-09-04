// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumVolumeSphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the minimum volume sphere containing the input set of points. The
// algorithm randomly permutes the input points so that the construction
// occurs in 'expected' O(N) time. All internal minimal sphere calculations
// store the squared radius in the radius member of Sphere3. Only at the end
// is a sqrt computed.
//
// Upstream templates the class on a ComputeType that is intended to be
// BSRational<T> for exact rational arithmetic. With exact arithmetic the
// result is guaranteed to be the minimum-volume sphere. This port
// instantiates the floating-point compute type (number), matching the
// upstream double-precision instantiation, so floating-point rounding errors
// can cause the updateSupport{2,3,4} functions to fail. The failure is
// trapped in those functions and a simple bounding sphere is computed using
// getContainerSphere3 from ContSphere3.ts. That sphere is generally not the
// minimum-volume sphere containing the points, and the minimum-volume
// algorithm is terminated immediately. The sphere is returned along with a
// boolean 'success' that is true when the sphere is minimum volume and false
// when the failure is trapped. When false is returned, another call to
// compute(...) can be tried: the random shuffle uses the persistent engine
// state, so it is highly likely to differ from the previous shuffle and
// there is a chance the algorithm succeeds just because of the different
// ordering of points.
//
// KNOWN UPSTREAM DEFECT (preserved), the 3D form of the one documented in
// MinimumAreaCircle2.ts. The claim above that a floating-point failure is
// always trapped is false: compute(...) can return success = true together
// with a sphere that does not contain every input point.
//
// updateSupport{2,3,4} always rewrite mSupport/mNumSupport, while the caller
// keeps the returned sphere only when its radius exceeds the current one.
// When round-off makes contains() reject a point that lies exactly on the
// current sphere (cospherical input) or makes exactSphere{3,4} return a
// garbage solve instead of the "radius = MAX_VALUE" sentinel (degenerate
// input, where LinearSystem's determinant != 0 test sees a cancellation
// residue), an update runs that would not run in exact arithmetic and
// returns a *smaller* sphere. The sphere is discarded but the support set is
// not, so a point can be left uncovered - either hidden inside the support
// set (supportContains skips it for the rest of the loop) or, because the
// last update index is not revisited either, simply never rechecked against
// the final sphere.
//
// Two reproductions are pinned in test/MinimumVolumeSphere3.test.ts: seven
// integer points where the missed point is a support point (outside by
// 0.084) and eleven where it is not (outside by 0.40). Surveys of 20000
// random integer point sets: no failure for coordinates in [-8,8]^3, about
// 0.05% for the much more degenerate [-2,2]^3. No failure was observed for
// random real-valued coordinates (80000 sets in 2D and 3D). A final
// containment check cannot repair this: the returned radius is the square
// root of the internal squared radius, so re-squaring it rejects a large
// fraction of the *correct* results. The real fix is upstream's exact
// ComputeType (BSRational), which this port does not instantiate. The
// behavior is preserved verbatim; see the "Upstream bug suspects" note of
// the V10 verification PR.
//
// Port notes: see MinimumAreaCircle2.ts, whose port decisions this file
// mirrors (compute(points) returning { minimal, success }, the
// numSupport/support accessors, the deterministic minstd_rand0-style
// shuffle, and the index tie-break in the uniqueness sort).

import { getContainerSphere3 } from './ContSphere3.js';
import { Hypersphere, type Sphere3 } from './Hypersphere.js';
import { LinearSystem } from './LinearSystem.js';
import { logAssert, logError } from './Logger.js';
import { Matrix } from './Matrix.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// The result of the minimum-volume-sphere query. 'success' is upstream's
// boolean return value: true when 'minimal' is the minimum-volume sphere,
// false when floating-point round-off caused a trapped failure and 'minimal'
// is the (generally larger) sphere from getContainerSphere3.
export interface MinimumVolumeSphere3Result {
    minimal: Sphere3;
    success: boolean;
}

// The multiplicative congruential generator x -> 16807 * x mod (2^31 - 1)
// with seed 1 (std::minstd_rand0, the common choice for
// std::default_random_engine).
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

// The internal pairing of a candidate sphere with the flag that reports
// whether the update succeeded (upstream's std::pair<Sphere3,bool>).
interface UpdateResult {
    sphere: Sphere3;
    valid: boolean;
}

export class MinimumVolumeSphere3 {
    private mNumSupport: number;
    private mSupport: number[];
    private mDRE: DefaultRandomEngine;
    private mComputePoints: Vector[];

    constructor() {
        this.mNumSupport = 0;
        this.mSupport = [0, 0, 0, 0];
        this.mDRE = new DefaultRandomEngine();
        this.mComputePoints = [];
    }

    compute(points: readonly Vector[]): MinimumVolumeSphere3Result {
        logAssert(points.length >= 1, 'Input must contain points.');
        for (const point of points) {
            logAssert(point.size === 3,
                'MinimumVolumeSphere3: points must be 3D.');
        }

        // Function array to avoid a switch statement in the main loop.
        const update: Array<(i: number) => UpdateResult> = [
            () => logError('unreachable'),
            (i: number) => this.updateSupport1(i),
            (i: number) => this.updateSupport2(i),
            (i: number) => this.updateSupport3(i),
            (i: number) => this.updateSupport4(i)
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
        let ctMinimal = this.exactSphere1(0);
        this.mNumSupport = 1;
        this.mSupport[0] = 0;

        // The loop restarts from the beginning of the point list each time
        // the sphere needs updating. Linus Kallberg (Computer Science at
        // Malardalen University in Sweden) discovered that performance is
        // better when the remaining points in the array are processed before
        // restarting. The points processed before the point that caused the
        // update are likely to be enclosed by the new sphere (or near the
        // sphere boundary) because they were enclosed by the previous sphere.
        // The chances are better that points after the current one will cause
        // growth of the bounding sphere.
        for (let i = 1 % numPoints, n = 0; i !== n; i = (i + 1) % numPoints) {
            if (!this.supportContains(i)) {
                if (!this.contains(i, ctMinimal)) {
                    const result = update[this.mNumSupport](i);
                    if (result.valid) {
                        if (result.sphere.radius > ctMinimal.radius) {
                            ctMinimal = result.sphere;
                            n = i;
                        }
                    } else {
                        // This case can happen because the compute type is
                        // floating point; see the comments at the beginning
                        // of this file. Return a non-minimal sphere.
                        //
                        // Upstream bug (MinimumVolumeSphere3.h, the trapped
                        // failure branch of operator()): the call is
                        // GetContainer(numPoints, points, minimal), but
                        // numPoints has been overwritten with the number of
                        // *unique* points while points[] is still the full
                        // input array. The bounding sphere is therefore
                        // computed from only the first numPoints entries of
                        // the input and need not contain all of the points.
                        // The port fixes this by using all input points.
                        this.mNumSupport = 0;
                        this.mSupport = [0, 0, 0, 0];
                        return {
                            minimal: getContainerSphere3(points),
                            success: false
                        };
                    }
                }
            }
        }

        const minimal = new Hypersphere(3);
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

    // Test whether point P is inside sphere S using squared distance and
    // squared radius. In this algorithm, sphere.radius is the *squared
    // radius* until compute(...) returns, at which time a square root is
    // applied.
    private contains(i: number, sphere: Sphere3): boolean {
        const diff = sub(this.mComputePoints[i], sphere.center);
        return dot(diff, diff) <= sphere.radius;
    }

    private exactSphere1(i0: number): Sphere3 {
        const minimal = new Hypersphere(3);
        minimal.center = this.mComputePoints[i0].clone();
        minimal.radius = 0;
        return minimal;
    }

    private exactSphere2(i0: number, i1: number): Sphere3 {
        const P0 = this.mComputePoints[i0];
        const P1 = this.mComputePoints[i1];
        const minimal = new Hypersphere(3);
        minimal.center = mul(0.5, add(P0, P1));
        const diff = sub(P1, P0);
        minimal.radius = 0.25 * dot(diff, diff);
        return minimal;
    }

    private exactSphere3(i0: number, i1: number, i2: number): Sphere3 {
        // Compute the circle containing P0, P1, and P2. The center in
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

        const minimal = new Hypersphere(3);
        const { X, invertible } = LinearSystem.solve2x2(A, B);
        if (invertible) {
            const x2 = 1 - X.get(0) - X.get(1);
            minimal.center = add(add(mul(X.get(0), P0), mul(X.get(1), P1)),
                mul(x2, P2));
            const tmp = add(mul(X.get(0), E0), mul(X.get(1), E1));
            minimal.radius = dot(tmp, tmp);
        } else {
            minimal.center = Vector.zero(3);
            minimal.radius = Number.MAX_VALUE;
        }
        return minimal;
    }

    private exactSphere4(i0: number, i1: number, i2: number, i3: number):
        Sphere3 {
        // Compute the sphere containing P0, P1, P2, and P3. The center in
        // barycentric coordinates is C = x0*P0 + x1*P1 + x2*P2 + x3*P3, where
        // x0 + x1 + x2 + x3 = 1. The center is equidistant from the four
        // points, so |C - P0| = |C - P1| = |C - P2| = |C - P3| = R, where R
        // is the radius of the sphere. From these conditions,
        //   C - P0 = x0*E0 + x1*E1 + x2*E2 - E0
        //   C - P1 = x0*E0 + x1*E1 + x2*E2 - E1
        //   C - P2 = x0*E0 + x1*E1 + x2*E2 - E2
        //   C - P3 = x0*E0 + x1*E1 + x2*E2
        // where E0 = P0 - P3, E1 = P1 - P3, and E2 = P2 - P3, which leads to
        //  r^2 = |x0*E0+x1*E1+x2*E2|^2-2*Dot(E0,x0*E0+x1*E1+x2*E2)+|E0|^2
        //  r^2 = |x0*E0+x1*E1+x2*E2|^2-2*Dot(E1,x0*E0+x1*E1+x2*E2)+|E1|^2
        //  r^2 = |x0*E0+x1*E1+x2*E2|^2-2*Dot(E2,x0*E0+x1*E1+x2*E2)+|E2|^2
        //  r^2 = |x0*E0+x1*E1+x2*E2|^2
        // Subtracting the last equation from the first three and writing the
        // equations as a linear system,
        //
        // +-                                -++   -+       +-          -+
        // | Dot(E0,E0) Dot(E0,E1) Dot(E0,E2) || x0 | = 0.5 | Dot(E0,E0) |
        // | Dot(E1,E0) Dot(E1,E1) Dot(E1,E2) || x1 |       | Dot(E1,E1) |
        // | Dot(E2,E0) Dot(E2,E1) Dot(E2,E2) || x2 |       | Dot(E2,E2) |
        // +-                                -++   -+       +-          -+
        //
        // The following code solves this system for x0, x1, and x2 and then
        // evaluates the fourth equation in r^2 to obtain r.
        const P0 = this.mComputePoints[i0];
        const P1 = this.mComputePoints[i1];
        const P2 = this.mComputePoints[i2];
        const P3 = this.mComputePoints[i3];

        const E0 = sub(P0, P3);
        const E1 = sub(P1, P3);
        const E2 = sub(P2, P3);

        const A = new Matrix(3, 3);
        A.set(0, 0, dot(E0, E0));
        A.set(0, 1, dot(E0, E1));
        A.set(0, 2, dot(E0, E2));
        A.set(1, 0, A.get(0, 1));
        A.set(1, 1, dot(E1, E1));
        A.set(1, 2, dot(E1, E2));
        A.set(2, 0, A.get(0, 2));
        A.set(2, 1, A.get(1, 2));
        A.set(2, 2, dot(E2, E2));

        const half = 0.5;
        const B = Vector.fromArray([
            half * A.get(0, 0), half * A.get(1, 1), half * A.get(2, 2)]);

        const minimal = new Hypersphere(3);
        const { X, invertible } = LinearSystem.solve3x3(A, B);
        if (invertible) {
            const x3 = 1 - X.get(0) - X.get(1) - X.get(2);
            minimal.center = add(
                add(mul(X.get(0), P0), mul(X.get(1), P1)),
                add(mul(X.get(2), P2), mul(x3, P3)));
            const tmp = add(
                add(mul(X.get(0), E0), mul(X.get(1), E1)),
                mul(X.get(2), E2));
            minimal.radius = dot(tmp, tmp);
        } else {
            minimal.center = Vector.zero(3);
            minimal.radius = Number.MAX_VALUE;
        }
        return minimal;
    }

    private updateSupport1(i: number): UpdateResult {
        const minimal = this.exactSphere2(this.mSupport[0], i);
        this.mNumSupport = 2;
        this.mSupport[1] = i;
        return { sphere: minimal, valid: true };
    }

    private updateSupport2(i: number): UpdateResult {
        // Permutations of type 2, used for calling exactSphere2(...).
        const type2: readonly (readonly number[])[] = [
            [0, /*2*/ 1],
            [1, /*2*/ 0]
        ];

        // Permutations of type 3, used for calling exactSphere3(...): the
        // single permutation {0, 1, 2}.

        const sphere: Sphere3[] = [];
        let minRSqr = Number.MAX_VALUE;
        let iSphere = 0;
        let iMinRSqr = -1;
        let k0: number;
        let k1: number;

        // Permutations of type 2.
        for (let j = 0; j < type2.length; ++j, ++iSphere) {
            k0 = this.mSupport[type2[j][0]];
            sphere[iSphere] = this.exactSphere2(k0, i);
            if (sphere[iSphere].radius < minRSqr) {
                k1 = this.mSupport[type2[j][1]];
                if (this.contains(k1, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
                }
            }
        }

        // Permutations of type 3.
        k0 = this.mSupport[0];
        k1 = this.mSupport[1];
        sphere[iSphere] = this.exactSphere3(k0, k1, i);
        if (sphere[iSphere].radius < minRSqr) {
            minRSqr = sphere[iSphere].radius;
            iMinRSqr = iSphere;
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
                // When this happens, use a simple bounding sphere for the
                // result and terminate the minimum-volume algorithm.
                return { sphere: new Hypersphere(3), valid: false };
        }

        return { sphere: sphere[iMinRSqr], valid: true };
    }

    private updateSupport3(i: number): UpdateResult {
        // Permutations of type 2, used for calling exactSphere2(...).
        const type2: readonly (readonly number[])[] = [
            [0, /*3*/ 1, 2],
            [1, /*3*/ 0, 2],
            [2, /*3*/ 0, 1]
        ];

        // Permutations of type 3, used for calling exactSphere3(...).
        const type3: readonly (readonly number[])[] = [
            [0, 1, /*3*/ 2],
            [0, 2, /*3*/ 1],
            [1, 2, /*3*/ 0]
        ];

        // Permutations of type 4, used for calling exactSphere4(...): the
        // single permutation {0, 1, 2, 3}.

        const sphere: Sphere3[] = [];
        let minRSqr = Number.MAX_VALUE;
        let iSphere = 0;
        let iMinRSqr = -1;
        let k0: number;
        let k1: number;
        let k2: number;

        // Permutations of type 2.
        for (let j = 0; j < type2.length; ++j, ++iSphere) {
            k0 = this.mSupport[type2[j][0]];
            sphere[iSphere] = this.exactSphere2(k0, i);
            if (sphere[iSphere].radius < minRSqr) {
                k1 = this.mSupport[type2[j][1]];
                k2 = this.mSupport[type2[j][2]];
                if (this.contains(k1, sphere[iSphere])
                    && this.contains(k2, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
                }
            }
        }

        // Permutations of type 3.
        for (let j = 0; j < type3.length; ++j, ++iSphere) {
            k0 = this.mSupport[type3[j][0]];
            k1 = this.mSupport[type3[j][1]];
            sphere[iSphere] = this.exactSphere3(k0, k1, i);
            if (sphere[iSphere].radius < minRSqr) {
                k2 = this.mSupport[type3[j][2]];
                if (this.contains(k2, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
                }
            }
        }

        // Permutations of type 4.
        k0 = this.mSupport[0];
        k1 = this.mSupport[1];
        k2 = this.mSupport[2];
        sphere[iSphere] = this.exactSphere4(k0, k1, k2, i);
        if (sphere[iSphere].radius < minRSqr) {
            minRSqr = sphere[iSphere].radius;
            iMinRSqr = iSphere;
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
            case 6:
                this.mNumSupport = 4;
                this.mSupport[3] = i;
                break;
            default:
                // See the comment in updateSupport2 for iMinRSqr == -1.
                return { sphere: new Hypersphere(3), valid: false };
        }

        return { sphere: sphere[iMinRSqr], valid: true };
    }

    private updateSupport4(i: number): UpdateResult {
        // Permutations of type 2, used for calling exactSphere2(...).
        const type2: readonly (readonly number[])[] = [
            [0, /*4*/ 1, 2, 3],
            [1, /*4*/ 0, 2, 3],
            [2, /*4*/ 0, 1, 3],
            [3, /*4*/ 0, 1, 2]
        ];

        // Permutations of type 3, used for calling exactSphere3(...).
        const type3: readonly (readonly number[])[] = [
            [0, 1, /*4*/ 2, 3],
            [0, 2, /*4*/ 1, 3],
            [0, 3, /*4*/ 1, 2],
            [1, 2, /*4*/ 0, 3],
            [1, 3, /*4*/ 0, 2],
            [2, 3, /*4*/ 0, 1]
        ];

        // Permutations of type 4, used for calling exactSphere4(...).
        const type4: readonly (readonly number[])[] = [
            [0, 1, 2, /*4*/ 3],
            [0, 1, 3, /*4*/ 2],
            [0, 2, 3, /*4*/ 1],
            [1, 2, 3, /*4*/ 0]
        ];

        const sphere: Sphere3[] = [];
        let minRSqr = Number.MAX_VALUE;
        let iSphere = 0;
        let iMinRSqr = -1;
        let k0: number;
        let k1: number;
        let k2: number;
        let k3: number;

        // Permutations of type 2.
        for (let j = 0; j < type2.length; ++j, ++iSphere) {
            k0 = this.mSupport[type2[j][0]];
            sphere[iSphere] = this.exactSphere2(k0, i);
            if (sphere[iSphere].radius < minRSqr) {
                k1 = this.mSupport[type2[j][1]];
                k2 = this.mSupport[type2[j][2]];
                k3 = this.mSupport[type2[j][3]];
                if (this.contains(k1, sphere[iSphere])
                    && this.contains(k2, sphere[iSphere])
                    && this.contains(k3, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
                }
            }
        }

        // Permutations of type 3.
        for (let j = 0; j < type3.length; ++j, ++iSphere) {
            k0 = this.mSupport[type3[j][0]];
            k1 = this.mSupport[type3[j][1]];
            sphere[iSphere] = this.exactSphere3(k0, k1, i);
            if (sphere[iSphere].radius < minRSqr) {
                k2 = this.mSupport[type3[j][2]];
                k3 = this.mSupport[type3[j][3]];
                if (this.contains(k2, sphere[iSphere])
                    && this.contains(k3, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
                }
            }
        }

        // Permutations of type 4.
        for (let j = 0; j < type4.length; ++j, ++iSphere) {
            k0 = this.mSupport[type4[j][0]];
            k1 = this.mSupport[type4[j][1]];
            k2 = this.mSupport[type4[j][2]];
            sphere[iSphere] = this.exactSphere4(k0, k1, k2, i);
            if (sphere[iSphere].radius < minRSqr) {
                k3 = this.mSupport[type4[j][3]];
                if (this.contains(k3, sphere[iSphere])) {
                    minRSqr = sphere[iSphere].radius;
                    iMinRSqr = iSphere;
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
                this.mNumSupport = 2;
                this.mSupport[0] = this.mSupport[3];
                this.mSupport[1] = i;
                break;
            case 4:
                this.mNumSupport = 3;
                this.mSupport[2] = i;
                break;
            case 5:
                this.mNumSupport = 3;
                this.mSupport[1] = i;
                break;
            case 6:
                this.mNumSupport = 3;
                this.mSupport[1] = this.mSupport[3];
                this.mSupport[2] = i;
                break;
            case 7:
                this.mNumSupport = 3;
                this.mSupport[0] = i;
                break;
            case 8:
                this.mNumSupport = 3;
                this.mSupport[0] = this.mSupport[3];
                this.mSupport[2] = i;
                break;
            case 9:
                this.mNumSupport = 3;
                this.mSupport[0] = this.mSupport[3];
                this.mSupport[1] = i;
                break;
            case 10:
                this.mSupport[3] = i;
                break;
            case 11:
                this.mSupport[2] = i;
                break;
            case 12:
                this.mSupport[1] = i;
                break;
            case 13:
                this.mSupport[0] = i;
                break;
            default:
                // See the comment in updateSupport2 for iMinRSqr == -1.
                return { sphere: new Hypersphere(3), valid: false };
        }

        return { sphere: sphere[iMinRSqr], valid: true };
    }

    // Indices of points that support the current minimum volume sphere.
    private supportContains(j: number): boolean {
        for (let i = 0; i < this.mNumSupport; ++i) {
            if (j === this.mSupport[i]) {
                return true;
            }
        }
        return false;
    }
}
