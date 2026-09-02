// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrEllipsoid3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query classifies the relationship of two solid ellipsoids: separated,
// intersecting, or one containing the other.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrEllipsoid3Ellipsoid3TI.
// The nested enum class 'Classification' becomes the file-qualified exported
// enum IntrEllipsoid3Ellipsoid3Classification. The three private 'GetRoots'
// overloads (distinguished in C++ by arity) become the module-private
// functions getRoots1, getRoots2 and getRoots3, each returning the array of
// roots rather than filling a caller-supplied buffer.

import type { Ellipsoid3 } from './Hyperellipsoid';
import { logAssert, logError } from './Logger';
import { Matrix, multiplyAB, multiplyATB, mulMatrix } from './Matrix';
import { RootsBisection } from './RootsBisection';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3';
import { Vector, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The relationship of the two solid ellipsoids.
export enum IntrEllipsoid3Ellipsoid3Classification {
    ELLIPSOIDS_SEPARATED,
    ELLIPSOIDS_INTERSECTING,
    ELLIPSOID0_CONTAINS_ELLIPSOID1,
    ELLIPSOID1_CONTAINS_ELLIPSOID0,
    INVALID
}

// The result of IntrEllipsoid3Ellipsoid3TI.test.
export interface IntrEllipsoid3Ellipsoid3TIResult {
    // As solids, the ellipsoids intersect as long as they are not separated.
    intersect: boolean;

    // This is one of the five enumerations listed above.
    classification: IntrEllipsoid3Ellipsoid3Classification;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrEllipsoid3Ellipsoid3TIResult {
    return {
        intersect: false,
        classification: IntrEllipsoid3Ellipsoid3Classification.INVALID
    };
}

const MAX_BISECTION_ITERATIONS = 1024;

// f(s) = d0*c0/(d0*s-1)^2 - 1
function getRoots1(d0: number, c0: number): number[] {
    const temp = Math.sqrt(d0 * c0);
    const inv = 1 / d0;
    return [(1 - temp) * inv, (1 + temp) * inv];
}

// f(s) = d0*c0/(d0*s-1)^2 + d1*c1/(d1*s-1)^2 - 1, with d0 > d1
function getRoots2(d0: number, d1: number, c0: number, c1: number): number[] {
    const d0c0 = d0 * c0;
    const d1c1 = d1 * c1;

    const F = (s: number): number => {
        const invN0 = 1 / (d0 * s - 1);
        const invN1 = 1 / (d1 * s - 1);
        const term0 = d0c0 * invN0 * invN0;
        const term1 = d1c1 * invN1 * invN1;
        return term0 + term1 - 1;
    };

    const DF = (s: number): number => {
        const invN0 = 1 / (d0 * s - 1);
        const invN1 = 1 / (d1 * s - 1);
        const term0 = d0 * d0c0 * invN0 * invN0 * invN0;
        const term1 = d1 * d1c1 * invN1 * invN1 * invN1;
        return -2 * (term0 + term1);
    };

    const roots: number[] = [];

    // Upstream TODO: what role does epsilon play?
    const epsilon = 0.001;
    const multiplier0 = Math.sqrt(2 / (1 - epsilon));
    const multiplier1 = Math.sqrt(1 / (1 + epsilon));
    const sqrtd0c0 = Math.sqrt(d0c0);
    const sqrtd1c1 = Math.sqrt(d1c1);
    const invD0 = 1 / d0;
    const invD1 = 1 / d1;

    // Compute the root in (-infinity,1/d0).
    let temp0 = (1 - multiplier0 * sqrtd0c0) * invD0;
    let temp1 = (1 - multiplier0 * sqrtd1c1) * invD1;
    let smin = Math.min(temp0, temp1);
    logAssert(F(smin) < 0, 'Unexpected condition.');
    let smax = (1 - multiplier1 * sqrtd0c0) * invD0;
    logAssert(F(smax) > 0, 'Unexpected condition.');
    let bisect = RootsBisection.find(F, smin, smax, MAX_BISECTION_ITERATIONS);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    // Compute roots (if any) in (1/d0,1/d1). It is the case that
    //   F(1/d0) = +infinity, F'(1/d0) = -infinity
    //   F(1/d1) = +infinity, F'(1/d1) = +infinity
    //   F"(s) > 0 for all s in the domain of F
    // Compute the unique root r of F'(s) on (1/d0,1/d1). The bisector needs
    // only the signs at the endpoints, so we pass -1 and +1 instead of the
    // infinite values. If F(r) < 0, F(s) has two roots in the interval. If
    // F(r) = 0, F(s) has only one root in the interval.
    const mid = RootsBisection.find(DF, invD0, invD1, -1, 1,
        MAX_BISECTION_ITERATIONS);
    logAssert(mid.iterations > 0, 'Unexpected condition.');
    if (F(mid.root) < 0) {
        // Pass in signs rather than infinities, because the bisector cares
        // only about the signs.
        bisect = RootsBisection.find(F, invD0, mid.root, 1, -1,
            MAX_BISECTION_ITERATIONS);
        logAssert(bisect.iterations > 0, 'Unexpected condition.');
        roots.push(bisect.root);
        bisect = RootsBisection.find(F, mid.root, invD1, -1, 1,
            MAX_BISECTION_ITERATIONS);
        logAssert(bisect.iterations > 0, 'Unexpected condition.');
        roots.push(bisect.root);
    }

    // Compute the root in (1/d1,+infinity).
    temp0 = (1 + multiplier0 * sqrtd0c0) * invD0;
    temp1 = (1 + multiplier0 * sqrtd1c1) * invD1;
    smax = Math.max(temp0, temp1);
    logAssert(F(smax) < 0, 'Unexpected condition.');
    smin = (1 + multiplier1 * sqrtd1c1) * invD1;
    logAssert(F(smin) > 0, 'Unexpected condition.');
    bisect = RootsBisection.find(F, smin, smax, MAX_BISECTION_ITERATIONS);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    return roots;
}

// f(s) = d0*c0/(d0*s-1)^2 + d1*c1/(d1*s-1)^2 + d2*c2/(d2*s-1)^2 - 1,
// with d0 > d1 > d2
function getRoots3(d0: number, d1: number, d2: number, c0: number, c1: number,
    c2: number): number[] {
    const d0c0 = d0 * c0;
    const d1c1 = d1 * c1;
    const d2c2 = d2 * c2;

    const F = (s: number): number => {
        const invN0 = 1 / (d0 * s - 1);
        const invN1 = 1 / (d1 * s - 1);
        const invN2 = 1 / (d2 * s - 1);
        const term0 = d0c0 * invN0 * invN0;
        const term1 = d1c1 * invN1 * invN1;
        const term2 = d2c2 * invN2 * invN2;
        return term0 + term1 + term2 - 1;
    };

    const DF = (s: number): number => {
        const invN0 = 1 / (d0 * s - 1);
        const invN1 = 1 / (d1 * s - 1);
        const invN2 = 1 / (d2 * s - 1);
        const term0 = d0 * d0c0 * invN0 * invN0 * invN0;
        const term1 = d1 * d1c1 * invN1 * invN1 * invN1;
        const term2 = d2 * d2c2 * invN2 * invN2 * invN2;
        return -2 * (term0 + term1 + term2);
    };

    const roots: number[] = [];

    // Upstream TODO: what role does epsilon play?
    const epsilon = 0.001;
    const multiplier0 = Math.sqrt(3 / (1 - epsilon));
    const multiplier1 = Math.sqrt(1 / (1 + epsilon));
    const sqrtd0c0 = Math.sqrt(d0c0);
    const sqrtd1c1 = Math.sqrt(d1c1);
    const sqrtd2c2 = Math.sqrt(d2c2);
    const invD0 = 1 / d0;
    const invD1 = 1 / d1;
    const invD2 = 1 / d2;

    // Compute the root in (-infinity,1/d0).
    let temp0 = (1 - multiplier0 * sqrtd0c0) * invD0;
    let temp1 = (1 - multiplier0 * sqrtd1c1) * invD1;
    let temp2 = (1 - multiplier0 * sqrtd2c2) * invD2;
    let smin = Math.min(Math.min(temp0, temp1), temp2);
    logAssert(F(smin) < 0, 'Unexpected condition.');
    let smax = (1 - multiplier1 * sqrtd0c0) * invD0;
    logAssert(F(smax) > 0, 'Unexpected condition.');
    let bisect = RootsBisection.find(F, smin, smax, MAX_BISECTION_ITERATIONS);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    // Compute roots (if any) in (1/d0,1/d1) and then in (1/d1,1/d2). In each
    // interval it is the case that F is +infinity at the endpoints, F' runs
    // from -infinity to +infinity, and F" > 0 on the domain of F. Compute the
    // unique root r of F'(s) in the interval. The bisector needs only the
    // signs at the endpoints, so we pass -1 and +1 instead of the infinite
    // values. If F(r) < 0, F(s) has two roots in the interval. If F(r) = 0,
    // F(s) has only one root in the interval.
    const intervals: [number, number][] = [[invD0, invD1], [invD1, invD2]];
    for (const [lo, hi] of intervals) {
        const mid = RootsBisection.find(DF, lo, hi, -1, 1,
            MAX_BISECTION_ITERATIONS);
        logAssert(mid.iterations > 0, 'Unexpected condition.');
        if (F(mid.root) < 0) {
            // Pass in signs rather than infinities, because the bisector
            // cares only about the signs.
            bisect = RootsBisection.find(F, lo, mid.root, 1, -1,
                MAX_BISECTION_ITERATIONS);
            logAssert(bisect.iterations > 0, 'Unexpected condition.');
            roots.push(bisect.root);
            bisect = RootsBisection.find(F, mid.root, hi, -1, 1,
                MAX_BISECTION_ITERATIONS);
            logAssert(bisect.iterations > 0, 'Unexpected condition.');
            roots.push(bisect.root);
        }
    }

    // Compute the root in (1/d2,+infinity).
    temp0 = (1 + multiplier0 * sqrtd0c0) * invD0;
    temp1 = (1 + multiplier0 * sqrtd1c1) * invD1;
    temp2 = (1 + multiplier0 * sqrtd2c2) * invD2;
    smax = Math.max(Math.max(temp0, temp1), temp2);
    logAssert(F(smax) < 0, 'Unexpected condition.');
    smin = (1 + multiplier1 * sqrtd2c2) * invD2;
    logAssert(F(smin) > 0, 'Unexpected condition.');
    bisect = RootsBisection.find(F, smin, smax, MAX_BISECTION_ITERATIONS);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    return roots;
}

// Build a 3x3 matrix whose columns are the specified vectors.
function matrixFromColumns(c0: Vector, c1: Vector, c2: Vector): Matrix {
    const M = Matrix.zero(3, 3);
    M.setCol(0, c0);
    M.setCol(1, c1);
    M.setCol(2, c2);
    return M;
}

// Build a 3x3 diagonal matrix.
function diagonal(d0: number, d1: number, d2: number): Matrix {
    const M = Matrix.zero(3, 3);
    M.set(0, 0, d0);
    M.set(1, 1, d1);
    M.set(2, 2, d2);
    return M;
}

// Test-intersection query for two solid ellipsoids in 3D.
export class IntrEllipsoid3Ellipsoid3TI implements
    TIQuery<Ellipsoid3, Ellipsoid3, IntrEllipsoid3Ellipsoid3TIResult> {

    test(ellipsoid0: Ellipsoid3, ellipsoid1: Ellipsoid3):
        IntrEllipsoid3Ellipsoid3TIResult {
        logAssert(ellipsoid0.dimension === 3 && ellipsoid1.dimension === 3,
            'IntrEllipsoid3Ellipsoid3TI: mismatched sizes.');

        const result = defaultTIResult();

        // Get the parameters of ellipsoid0.
        const K0 = ellipsoid0.center;
        const e0 = ellipsoid0.extent.values;
        const R0 = matrixFromColumns(ellipsoid0.axis[0], ellipsoid0.axis[1],
            ellipsoid0.axis[2]);

        // Get the parameters of ellipsoid1.
        const K1 = ellipsoid1.center;
        const e1 = ellipsoid1.extent.values;
        const R1 = matrixFromColumns(ellipsoid1.axis[0], ellipsoid1.axis[1],
            ellipsoid1.axis[2]);
        const D1 = diagonal(1 / (e1[0] * e1[0]), 1 / (e1[1] * e1[1]),
            1 / (e1[2] * e1[2]));

        // Compute K2.
        const D0NegHalf = diagonal(e0[0], e0[1], e0[2]);
        const D0Half = diagonal(1 / e0[0], 1 / e0[1], 1 / e0[2]);
        const K2 = mulMatrix(D0Half,
            mulMatrix(sub(K1, K0), R0) as Vector) as Vector;

        // Compute M2.
        const R1TR0D0NegHalf = multiplyATB(R1, multiplyAB(R0, D0NegHalf));
        const M2 = multiplyAB(multiplyATB(R1TR0D0NegHalf, D1), R1TR0D0NegHalf);

        // Factor M2 = R*D*R^T.
        const es = new SymmetricEigensolver3x3();
        const { evals: D, evecs: evec } = es.solve(M2.get(0, 0), M2.get(0, 1),
            M2.get(0, 2), M2.get(1, 1), M2.get(1, 2), M2.get(2, 2), false, +1);
        const R = matrixFromColumns(
            Vector.fromArray(evec[0]), Vector.fromArray(evec[1]),
            Vector.fromArray(evec[2]));

        // Compute K = R^T*K2.
        const K = mulMatrix(K2, R) as Vector;

        // Transformed ellipsoid0 is Z^T*Z = 1 and transformed ellipsoid1 is
        // (Z-K)^T*D*(Z-K) = 0.
        //
        // The minimum and maximum squared distances from the origin of points
        // on transformed ellipsoid1 are used to determine whether the
        // ellipsoids intersect, are separated, or one contains the other.
        let minSqrDistance = Number.MAX_VALUE;
        let maxSqrDistance = 0;

        if (K.values[0] === 0 && K.values[1] === 0 && K.values[2] === 0) {
            // The special case of common centers must be handled separately.
            // It is not possible for the ellipsoids to be separated.
            for (let i = 0; i < 3; ++i) {
                const invD = 1 / D[i];
                if (invD < minSqrDistance) {
                    minSqrDistance = invD;
                }
                if (invD > maxSqrDistance) {
                    maxSqrDistance = invD;
                }
            }

            if (maxSqrDistance < 1) {
                result.classification = IntrEllipsoid3Ellipsoid3Classification
                    .ELLIPSOID0_CONTAINS_ELLIPSOID1;
            }
            else if (minSqrDistance > 1) {
                result.classification = IntrEllipsoid3Ellipsoid3Classification
                    .ELLIPSOID1_CONTAINS_ELLIPSOID0;
            }
            else {
                result.classification = IntrEllipsoid3Ellipsoid3Classification
                    .ELLIPSOIDS_INTERSECTING;
            }
            result.intersect = true;
            return result;
        }

        // The closest point P0 and farthest point P1 are solutions to
        // s0*D*(P0 - K) = P0 and s1*D*(P1 - K) = P1 for some scalars s0 and
        // s1 that are roots to the function
        //   f(s) = d0*k0^2/(d0*s-1)^2 + d1*k1^2/(d1*s-1)^2
        //          + d2*k2^2/(d2*s-1)^2 - 1
        // where D = diagonal(d0,d1,d2) and K = (k0,k1,k2).
        const d0 = D[0], d1 = D[1], d2 = D[2];
        const c0 = K.values[0] * K.values[0];
        const c1 = K.values[1] * K.values[1];
        const c2 = K.values[2] * K.values[2];

        // Sort the values so that d0 >= d1 >= d2. This allows us to bound the
        // roots of f(s), of which there are at most 6.
        const param: [number, number][] = [[d0, c0], [d1, c1], [d2, c2]];
        param.sort((a, b) => {
            // The port of std::greater on std::pair (descending lexicographic
            // order).
            if (a[0] !== b[0]) {
                return a[0] > b[0] ? -1 : 1;
            }
            if (a[1] !== b[1]) {
                return a[1] > b[1] ? -1 : 1;
            }
            return 0;
        });

        const valid: [number, number][] = [];
        if (param[0][0] > param[1][0]) {
            if (param[1][0] > param[2][0]) {
                // d0 > d1 > d2
                for (let i = 0; i < 3; ++i) {
                    if (param[i][1] > 0) {
                        valid.push(param[i]);
                    }
                }
            }
            else {
                // d0 > d1 = d2
                if (param[0][1] > 0) {
                    valid.push(param[0]);
                }
                param[1][1] += param[0][1];
                if (param[1][1] > 0) {
                    valid.push(param[1]);
                }
            }
        }
        else {
            if (param[1][0] > param[2][0]) {
                // d0 = d1 > d2
                param[0][1] += param[1][1];
                if (param[0][1] > 0) {
                    valid.push(param[0]);
                }
                if (param[2][1] > 0) {
                    valid.push(param[2]);
                }
            }
            else {
                // d0 = d1 = d2
                param[0][1] += param[1][1] + param[2][1];
                if (param[0][1] > 0) {
                    valid.push(param[0]);
                }
            }
        }

        let roots: number[];
        if (valid.length === 3) {
            roots = getRoots3(valid[0][0], valid[1][0], valid[2][0],
                valid[0][1], valid[1][1], valid[2][1]);
        }
        else if (valid.length === 2) {
            roots = getRoots2(valid[0][0], valid[1][0], valid[0][1],
                valid[1][1]);
        }
        else if (valid.length === 1) {
            roots = getRoots1(valid[0][0], valid[0][1]);
        }
        else {
            // The number of valid pairs cannot be zero because the case K = 0
            // was already handled.
            logError('Unexpected condition.');
            roots = [];
        }

        for (const s of roots) {
            const p0 = d0 * K.values[0] * s / (d0 * s - 1);
            const p1 = d1 * K.values[1] * s / (d1 * s - 1);
            const p2 = d2 * K.values[2] * s / (d2 * s - 1);
            const sqrDistance = p0 * p0 + p1 * p1 + p2 * p2;
            if (sqrDistance < minSqrDistance) {
                minSqrDistance = sqrDistance;
            }
            if (sqrDistance > maxSqrDistance) {
                maxSqrDistance = sqrDistance;
            }
        }

        if (maxSqrDistance < 1) {
            result.intersect = true;
            result.classification = IntrEllipsoid3Ellipsoid3Classification
                .ELLIPSOID0_CONTAINS_ELLIPSOID1;
        }
        else if (minSqrDistance > 1) {
            if (d0 * c0 + d1 * c1 + d2 * c2 > 1) {
                result.intersect = false;
                result.classification = IntrEllipsoid3Ellipsoid3Classification
                    .ELLIPSOIDS_SEPARATED;
            }
            else {
                result.intersect = true;
                result.classification = IntrEllipsoid3Ellipsoid3Classification
                    .ELLIPSOID1_CONTAINS_ELLIPSOID0;
            }
        }
        else {
            result.intersect = true;
            result.classification = IntrEllipsoid3Ellipsoid3Classification
                .ELLIPSOIDS_INTERSECTING;
        }

        return result;
    }
}
