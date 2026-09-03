// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrCylinder3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test for intersection of two finite cylinders using the method of
// separating axes. The algorithm is described in the document
// https://www.geometrictools.com/Documentation/IntersectionOfCylinders.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrCylinder3Cylinder3TI.
// The result member is 'separated' (with the separating direction) rather
// than 'intersect'; a false value means no separating direction was found
// among those sampled, so the cylinders are (reported as) intersecting.
//
// The upstream query can distribute the hemisphere sampling over
// std::thread instances. The port keeps the 'numThreads' constructor
// parameter for API compatibility but always runs the single-threaded loop,
// which visits the same directions in the same order. The multithreaded
// upstream path can report a different separating direction than the
// single-threaded path (it takes the first thread that found one, not the
// first direction in scan order), so the port is deterministic.
//
// Upstream bug (preserved): the constructor comment documents
// phi[j] = pi * j / numPhi, but the code uses GTE_C_HALF_PI as the
// multiplier. Half pi is correct for sampling a hemisphere of directions
// (D and -D give the same separation test), so the comment is stale, not the
// code. The port keeps the code and drops the stale claim.
//
// Upstream bug (guarded here): the C++ uses cylinder.height unconditionally,
// but Cylinder3.h represents an infinite cylinder with the sentinel
// height = -1. Feeding that sentinel into the finite formulas silently
// produces wrong results. The port asserts that both cylinders are finite,
// using the same message as the upstream assertion in
// IntrCanonicalBox3Cylinder3.h. See also upstream issue #187.

import type { Cylinder3 } from './Cylinder3.js';
import { GTE_C_HALF_PI, GTE_C_TWO_PI } from './Constants.js';
import { logAssert } from './Logger.js';
import { Vector, add, sub, mul, dot, length, normalize } from './Vector.js';
import { cross, computeOrthogonalComplement3 } from './Vector3.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrCylinder3Cylinder3TI.test.
export interface IntrCylinder3Cylinder3TIResult {
    separated: boolean;
    separatingDirection: Vector;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrCylinder3Cylinder3TIResult {
    return { separated: false, separatingDirection: Vector.zero(3) };
}

export class IntrCylinder3Cylinder3TI implements
    TIQuery<Cylinder3, Cylinder3, IntrCylinder3Cylinder3TIResult> {

    private readonly mNumThreads: number;
    private readonly mNumTheta: number;
    private readonly mNumPhi: number;

    // Cylinder 0.
    private mW0: Vector = Vector.zero(3);  // W0
    private mR0: number = 0;               // r0
    private mHalfH0: number = 0;           // h0/2

    // Cylinder 1.
    private mW1: Vector = Vector.zero(3);  // W1
    private mR1: number = 0;               // r1
    private mHalfH1: number = 0;           // h1/2

    // Members dependent on both cylinders.
    private mDelta: Vector = Vector.zero(3);  // C1 - C0
    private mW0xW1: Vector = Vector.zero(3);  // Cross(W0, W1)

    // The potential separating directions are
    //   D(theta[i],phi[j]) = c0*s1 * U + s0*s1 * V + c1 * N
    // where {U,V,N} is a right-handed orthonormal basis with N the north pole
    // of a hemisphere. The parameters are theta[i] = 2*pi*i/numTheta with
    // 0 <= i < numTheta, phi[j] = (pi/2)*j/numPhi with 0 <= j < numPhi,
    // c0 = cos(theta[i]), s0 = sin(theta[i]), c1 = cos(phi[j]) and
    // s1 = sin(phi[j]).
    constructor(numThreads: number, numTheta: number, numPhi: number) {
        logAssert(numTheta > 0 && numPhi > 0, 'Invalid number of angles.');
        this.mNumThreads = numThreads;
        this.mNumTheta = numTheta;
        this.mNumPhi = numPhi;
    }

    test(cylinder0: Cylinder3, cylinder1: Cylinder3):
        IntrCylinder3Cylinder3TIResult {
        logAssert(cylinder0.isFinite() && cylinder1.isFinite(),
            'Infinite cylinders are not yet supported.');

        // The default result has separated = false and separatingDirection
        // = (0,0,0).
        const result = defaultTIResult();

        this.mDelta = sub(cylinder1.axis.origin, cylinder0.axis.origin);
        if (length(this.mDelta) === 0) {
            return result;
        }

        this.mW0 = cylinder0.axis.direction;
        this.mR0 = cylinder0.radius;
        this.mHalfH0 = 0.5 * cylinder0.height;
        this.mW1 = cylinder1.axis.direction;
        this.mR1 = cylinder1.radius;
        this.mHalfH1 = 0.5 * cylinder1.height;
        this.mW0xW1 = cross(this.mW0, this.mW1);
        const lengthW0xW1 = length(this.mW0xW1);
        if (lengthW0xW1 > 0) {
            // The cylinder directions are not parallel.

            // Test for separation by W0.
            const absDotW0W1 = Math.abs(dot(this.mW0, this.mW1));
            const absDotW0Delta = Math.abs(dot(this.mW0, this.mDelta));
            let test = this.mR1 * lengthW0xW1 + this.mHalfH0
                + this.mHalfH1 * absDotW0W1 - absDotW0Delta;
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = this.mW0.clone();
                return result;
            }

            // Test for separation by W1.
            const absDotW1Delta = Math.abs(dot(this.mW1, this.mDelta));
            test = this.mR0 * lengthW0xW1 + this.mHalfH0 * absDotW0W1
                + this.mHalfH1 - absDotW1Delta;
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = this.mW1.clone();
                return result;
            }

            // Test for separation by W0xW1.
            const absDotW0xW1Delta = Math.abs(dot(this.mW0xW1, this.mDelta));
            test = (this.mR0 + this.mR1) * lengthW0xW1 - absDotW0xW1Delta;
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = this.mW0xW1.clone();
                normalize(result.separatingDirection);
                return result;
            }

            // Test for separation by Delta.
            test = this.mR0 * length(cross(this.mDelta, this.mW0))
                + this.mR1 * length(cross(this.mDelta, this.mW1))
                + this.mHalfH0 * absDotW0Delta
                + this.mHalfH1 * absDotW1Delta
                - dot(this.mDelta, this.mDelta);
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = this.mDelta.clone();
                normalize(result.separatingDirection);
                return result;
            }

            // Test for separation by other directions.
            this.testForSeparation(result);
        } else {
            // The cylinder directions are parallel.

            // Test for separation by height.
            const dotDeltaW0 = dot(this.mDelta, this.mW0);
            let test = this.mHalfH0 + this.mHalfH1 - Math.abs(dotDeltaW0);
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = this.mW0.clone();
                return result;
            }

            // Test for separation radially.
            test = this.mR0 + this.mR1
                - length(cross(this.mDelta, this.mW0));
            if (test < 0) {
                result.separated = true;
                result.separatingDirection = sub(this.mDelta,
                    mul(dotDeltaW0, this.mW0));
                normalize(result.separatingDirection);
                return result;
            }
        }

        return result;
    }

    // The number of threads requested by the caller. The port always runs
    // single-threaded; see the port notes at the top of this file.
    get numThreads(): number {
        return this.mNumThreads;
    }

    private testForSeparation(result: IntrCylinder3Cylinder3TIResult): void {
        // Compute a right-handed orthonormal basis {U,V,N} so that N is the
        // north pole of a hemisphere.
        const basis: Vector[] = [this.mDelta.clone(), new Vector(3),
            new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const U = basis[1];
        const V = basis[2];
        const N = basis[0];

        // Sample the hemisphere for potential separating directions.
        const phiMultiplier = GTE_C_HALF_PI / this.mNumPhi;
        const thetaMultiplier = GTE_C_TWO_PI / this.mNumTheta;
        for (let j = 1; j < this.mNumPhi; ++j) {
            const phi = phiMultiplier * j;
            const c1 = Math.cos(phi);
            const s1 = Math.sin(phi);
            for (let i = 0; i < this.mNumTheta; ++i) {
                // Compute the potential separating direction.
                const theta = thetaMultiplier * i;
                const c0 = Math.cos(theta);
                const s0 = Math.sin(theta);
                const D = add(add(mul(c0 * s1, U), mul(s0 * s1, V)),
                    mul(c1, N));

                // Test for separation. If test is negative, the direction is
                // separating.
                const test =
                    this.mR0 * length(cross(this.mW0, D))
                    + this.mR1 * length(cross(this.mW1, D))
                    + this.mHalfH0 * Math.abs(dot(this.mW0, D))
                    + this.mHalfH1 * Math.abs(dot(this.mW1, D))
                    - Math.abs(dot(this.mDelta, D));
                if (test < 0) {
                    result.separated = true;
                    result.separatingDirection = D;
                    return;
                }
            }
        }
    }
}
