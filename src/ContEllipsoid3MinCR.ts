// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContEllipsoid3MinCR.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the minimum-volume ellipsoid, (X-C)^T R D R^T (X-C) = 1, given the
// center C and orientation matrix R. The columns of R are the axes of the
// ellipsoid. The algorithm computes the diagonal matrix D. The minimum
// volume is (4*pi/3)/sqrt(D[0]*D[1]*D[2]), where D = diag(D[0],D[1],D[2]).
// The problem is equivalent to maximizing the product D[0]*D[1]*D[2] for a
// given C and R, and subject to the constraints
//   (P[i]-C)^T R D R^T (P[i]-C) <= 1
// for all input points P[i] with 0 <= i < N. Each constraint has the form
//   A[0]*D[0] + A[1]*D[1] + A[2]*D[2] <= 1
// where A[0] >= 0, A[1] >= 0, and A[2] >= 0.
//
// Port notes: upstream's functor class ContEllipsoid3MinCR becomes the free
// function getContainerEllipsoid3MinCR, following the Cont* naming precedent
// (container type as the suffix so the exported name is globally unique; see
// ContAlignedBox.ts and ContOrientedBox2.ts). Upstream writes the result into
// a caller-provided 'Real D[3]'; the port returns the three diagonal entries
// as a number[] of length 3. The private FindEdgeMax/FindFacetMax/MaxProduct
// helpers become module-private functions. Upstream passes the plane indices
// by non-const reference; the writes through those references are never read
// by any caller (MaxProduct ignores its 'plane' variable afterwards and
// FindFacetMax returns immediately after calling FindEdgeMax), so the port
// passes them by value. An empty point set is rejected with logAssert,
// matching the Cont* empty-input precedent (upstream's LogAssert on
// 'plane != -1' fires in that case).
//
// Randomness: upstream jitters the constraint planes with std::mt19937 and
// std::uniform_real_distribution<Real>(0,1). The port reproduces the
// Mersenne twister MT19937 with the C++ default seed (5489) and the standard
// std::generate_canonical<double,53> two-draw construction that libstdc++
// uses for uniform_real_distribution<double>(0,1), so the jitter sequence is
// deterministic and matches upstream's up to floating-point rounding of the
// final division. The generator is created fresh per call, exactly as
// upstream creates a local std::mt19937 in MaxProduct.
//
// One upstream bug is fixed. Both plane-sorting loops compute the slack
// 'numer = 1 - A[i].D' of a candidate blocking plane and their comments say
// "some numerical error may make this a small negative number. In that case
// set tmax = 0 (no change in position)". The code does not do that; it calls
// LogAssert(numer >= 0) instead. A plane that is already active at the
// current point routinely has a slack of about -1e-16, so upstream either
// aborts on ordinary input or (with assertions removed) computes a slightly
// negative t, which fails the subsequent '0 <= t' test. The active plane is
// then not treated as a blocker at all and the walk steps straight through
// it, after which the point is infeasible by an arbitrarily large amount and
// the returned ellipsoid no longer contains the input points. The port
// implements the clamp the comments describe (numer is clamped to 0, so a
// just-violated plane blocks with t = 0), which restores the invariant.
//
// With that invariant restored the walk always returns a containing
// ellipsoid, but it remains a heuristic: FindEdgeMax returns with the point
// unchanged when tMax is 0 (upstream's "tmax == 0, so return" case), so the
// walk can stall at a vertex of the constraint polytope and report an
// ellipsoid larger than the minimum-volume one. Roughly 4% of random point
// clouds stall this way. This is upstream behavior and is preserved.

import { logAssert } from './Logger.js';
import { type Matrix3x3 } from './Matrix3x3.js';
import { mulMatrix } from './Matrix.js';
import { Vector, compMul, sub } from './Vector.js';

// Compute the diagonal matrix D = diag(D[0],D[1],D[2]) of the minimum-volume
// ellipsoid with the specified center and orientation that contains the input
// points. The returned array has length 3.
export function getContainerEllipsoid3MinCR(points: readonly Vector[],
    C: Vector, R: Matrix3x3): number[] {
    logAssert(points.length > 0, 'getContainerEllipsoid3MinCR: no points.');
    logAssert(C.size === 3, 'getContainerEllipsoid3MinCR: center must be 3D.');
    logAssert(R.numRows === 3 && R.numCols === 3,
        'getContainerEllipsoid3MinCR: rotation must be 3x3.');
    for (const point of points) {
        logAssert(point.size === 3,
            'getContainerEllipsoid3MinCR: points must be 3D.');
    }

    // Compute the constraint coefficients, of the form (A[0],A[1],A[2]) for
    // each i.
    const A: Vector[] = new Array<Vector>(points.length);
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], C);  // P[i] - C
        const prod = mulMatrix(diff, R);  // R^T*(P[i] - C) = (u,v,w)
        A[i] = compMul(prod, prod);  // (u^2, v^2, w^2)
    }

    // TODO (upstream): Sort the constraints to eliminate redundant ones. It
    // is clear how to do this in ContEllipse2MinCR. How to do this in 3D?

    const D = [0, 0, 0];
    maxProduct(A, D);
    return D;
}

function findEdgeMax(A: readonly Vector[], plane0: number, plane1: number,
    D: number[]): void {
    // Compute direction to local maximum point on line of intersection.
    const a0 = A[plane0].values, a1v = A[plane1].values;
    let xDir = a0[1] * a1v[2] - a1v[1] * a0[2];
    let yDir = a0[2] * a1v[0] - a1v[2] * a0[0];
    let zDir = a0[0] * a1v[1] - a1v[0] * a0[1];

    // Build quadratic Q'(t) = (d/dt)(x(t)y(t)z(t)) = a0+a1*t+a2*t^2.
    const q0 = D[0] * D[1] * zDir + D[0] * D[2] * yDir + D[1] * D[2] * xDir;
    const q1 = 2 * (D[2] * xDir * yDir + D[1] * xDir * zDir +
        D[0] * yDir * zDir);
    const q2 = 3 * (xDir * yDir * zDir);

    // Find root to Q'(t) = 0 corresponding to maximum.
    let tFinal: number;
    if (q2 !== 0) {
        const invQ2 = 1 / q2;
        let discr = q1 * q1 - 4 * q0 * q2;
        discr = Math.sqrt(Math.max(discr, 0));
        tFinal = -0.5 * (q1 + discr) * invQ2;
        if (q1 + 2 * q2 * tFinal > 0) {
            tFinal = 0.5 * (-q1 + discr) * invQ2;
        }
    } else if (q1 !== 0) {
        tFinal = -q0 / q1;
    } else if (q0 !== 0) {
        const fmax = Number.MAX_VALUE;
        tFinal = (q0 >= 0 ? fmax : -fmax);
    } else {
        return;
    }

    if (tFinal < 0) {
        // Make (xDir,yDir,zDir) point in direction of increase of Q.
        tFinal = -tFinal;
        xDir = -xDir;
        yDir = -yDir;
        zDir = -zDir;
    }

    // Sort remaining planes along line from current point to local maximum.
    let tMax = tFinal;
    let plane2 = -1;
    const numPoints = A.length;
    for (let i = 0; i < numPoints; ++i) {
        if (i === plane0 || i === plane1) {
            continue;
        }

        const ai = A[i].values;
        const norDotDir = ai[0] * xDir + ai[1] * yDir + ai[2] * zDir;
        if (norDotDir <= 0) {
            continue;
        }

        // Theoretically the numerator must be nonnegative since an invariant
        // in the algorithm is that (x0,y0,z0) is on the convex hull of the
        // constraints. However, some numerical error may make this a small
        // negative number. In that case set t = 0 (no change in position).
        // See the port note about upstream's LogAssert here.
        const numer = Math.max(1 - ai[0] * D[0] - ai[1] * D[1] - ai[2] * D[2],
            0);

        const t = numer / norDotDir;
        if (0 <= t && t < tMax) {
            plane2 = i;
            tMax = t;
        }
    }

    D[0] += tMax * xDir;
    D[1] += tMax * yDir;
    D[2] += tMax * zDir;

    if (tMax === tFinal) {
        return;
    }

    if (tMax > 0) {
        findFacetMax(A, plane2, D);
        return;
    }

    // tMax == 0, so return with D[0], D[1] and D[2] unchanged.
}

function findFacetMax(A: readonly Vector[], plane0: number,
    D: number[]): void {
    let tFinal: number, xDir: number, yDir: number, zDir: number;
    const a0 = A[plane0].values;

    if (a0[0] > 0 && a0[1] > 0 && a0[2] > 0) {
        // Compute local maximum point on plane.
        const oneThird = 1 / 3;
        const xMax = oneThird / a0[0];
        const yMax = oneThird / a0[1];
        const zMax = oneThird / a0[2];

        // Compute direction to local maximum point on plane.
        tFinal = 1;
        xDir = xMax - D[0];
        yDir = yMax - D[1];
        zDir = zMax - D[2];
    } else {
        tFinal = Number.MAX_VALUE;
        xDir = (a0[0] > 0 ? 0 : 1);
        yDir = (a0[1] > 0 ? 0 : 1);
        zDir = (a0[2] > 0 ? 0 : 1);
    }

    // Sort remaining planes along line from current point.
    let tMax = tFinal;
    let plane1 = -1;
    const numPoints = A.length;
    for (let i = 0; i < numPoints; ++i) {
        if (i === plane0) {
            continue;
        }

        const ai = A[i].values;
        const norDotDir = ai[0] * xDir + ai[1] * yDir + ai[2] * zDir;
        if (norDotDir <= 0) {
            continue;
        }

        // Theoretically the numerator must be nonnegative because an
        // invariant in the algorithm is that (x0,y0,z0) is on the convex hull
        // of the constraints. However, some numerical error may make this a
        // small negative number. In that case set t = 0 (no change in
        // position). See the port note about upstream's LogAssert here.
        const numer = Math.max(1 - ai[0] * D[0] - ai[1] * D[1] - ai[2] * D[2],
            0);

        const t = numer / norDotDir;
        if (0 <= t && t < tMax) {
            plane1 = i;
            tMax = t;
        }
    }

    D[0] += tMax * xDir;
    D[1] += tMax * yDir;
    D[2] += tMax * zDir;

    if (tMax === 1) {
        return;
    }

    if (tMax > 0) {
        findFacetMax(A, plane1, D);
        return;
    }

    findEdgeMax(A, plane0, plane1, D);
}

function maxProduct(A: readonly Vector[], D: number[]): void {
    // Maximize x*y*z subject to x >= 0, y >= 0, z >= 0, and
    // A[i]*x+B[i]*y+C[i]*z <= 1 for 0 <= i < N where A[i] >= 0, B[i] >= 0 and
    // C[i] >= 0.

    // Jitter the lines to avoid cases where more than three planes intersect
    // at the same point. Should also break parallelism and planes parallel to
    // the coordinate planes.
    const mte = new MT19937();
    const maxJitter = 1e-12;
    const numPoints = A.length;
    let i: number;
    for (i = 0; i < numPoints; ++i) {
        const ai = A[i].values;
        ai[0] += maxJitter * mte.nextCanonical();
        ai[1] += maxJitter * mte.nextCanonical();
        ai[2] += maxJitter * mte.nextCanonical();
    }

    // Sort lines along the z-axis (x = 0 and y = 0).
    let plane = -1;
    let zmax = 0;
    for (i = 0; i < numPoints; ++i) {
        if (A[i].values[2] > zmax) {
            zmax = A[i].values[2];
            plane = i;
        }
    }
    logAssert(plane !== -1, 'Unexpected condition.');

    // Walk along convex hull searching for maximum.
    D[0] = 0;
    D[1] = 0;
    D[2] = 1 / zmax;
    findFacetMax(A, plane, D);
}

// The port of std::mt19937 (the 32-bit Mersenne twister with the standard
// parameters and the C++ default seed) together with the
// std::generate_canonical<double,53> construction used by
// std::uniform_real_distribution<double>(0,1).
class MT19937 {
    private static readonly N = 624;
    private static readonly M = 397;
    private static readonly MATRIX_A = 0x9908b0df;
    private static readonly UPPER_MASK = 0x80000000;
    private static readonly LOWER_MASK = 0x7fffffff;

    private readonly mt: Uint32Array;
    private mti: number;

    constructor(seed: number = 5489) {
        this.mt = new Uint32Array(MT19937.N);
        this.mt[0] = seed >>> 0;
        for (let i = 1; i < MT19937.N; ++i) {
            const prev = this.mt[i - 1] ^ (this.mt[i - 1] >>> 30);
            this.mt[i] = (Math.imul(1812433253, prev) + i) >>> 0;
        }
        this.mti = MT19937.N;
    }

    // The next 32-bit unsigned random number.
    next(): number {
        if (this.mti >= MT19937.N) {
            for (let k = 0; k < MT19937.N; ++k) {
                const y = ((this.mt[k] & MT19937.UPPER_MASK) |
                    (this.mt[(k + 1) % MT19937.N] & MT19937.LOWER_MASK)) >>> 0;
                let next = (this.mt[(k + MT19937.M) % MT19937.N] ^
                    (y >>> 1)) >>> 0;
                if ((y & 1) !== 0) {
                    next = (next ^ MT19937.MATRIX_A) >>> 0;
                }
                this.mt[k] = next;
            }
            this.mti = 0;
        }

        let y = this.mt[this.mti++];
        y = (y ^ (y >>> 11)) >>> 0;
        y = (y ^ ((y << 7) & 0x9d2c5680)) >>> 0;
        y = (y ^ ((y << 15) & 0xefc60000)) >>> 0;
        y = (y ^ (y >>> 18)) >>> 0;
        return y >>> 0;
    }

    // std::generate_canonical<double, 53>: with R = 2^32 and b = 53 bits of
    // mantissa, k = ceil(b/log2(R)) = 2 draws are summed and scaled by R^k.
    nextCanonical(): number {
        const g0 = this.next();
        const g1 = this.next();
        return (g0 + g1 * 4294967296) / 18446744073709551616;
    }
}
