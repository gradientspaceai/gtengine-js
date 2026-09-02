// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrEllipse2Ellipse2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection and find-intersection queries implemented here are
// discussed in the document
// https://www.geometrictools.com/Documentation/IntersectionOfEllipses.pdf
// The number type should support exact rational arithmetic in order for the
// polynomial root construction to be robust. The classification of the
// intersections depends on various sign tests of computed values. When those
// values are computed with floating-point arithmetic (as in this port, which
// uses IEEE double), the sign tests can lead to misclassification.
//
// The find-intersection query had some robustness issues when computing with
// floating-point only. The current implementation fixes those. The algorithm
// is described in
// https://www.geometrictools.com/Documentation/RobustIntersectionOfEllipses.pdf
//
// Port notes:
// - Per the Intr* precedent the two upstream query classes become
//   IntrEllipse2Ellipse2TI and IntrEllipse2Ellipse2FI. The TI query returns
//   the classification enum directly (upstream has no Result struct for it).
// - The FI query has two 'operator()' overloads. The canonical
//   ellipse-versus-ellipse query keeps the name 'find'; the overload taking
//   the standard forms (center, matrix) is 'findStandardForm'. Likewise the
//   two 'ComputeAlignedBox' overloads become TypeScript overloads of
//   'computeAlignedBox' distinguished by arity.
// - Upstream reports "the ellipses are the same" with
//   numPoints = std::numeric_limits<size_t>::max(); the port uses the
//   exported constant intrEllipse2Ellipse2InfinitePoints
//   (Number.MAX_SAFE_INTEGER), following the SIZE_MAX precedent of BVTree.ts.
// - The private case handlers and root finders become module-private
//   functions.

import { AlignedBox } from './AlignedBox';
import { fma, robustDOP, robustSOP } from './Functions';
import type { Hyperellipsoid } from './Hyperellipsoid';
import { logAssert } from './Logger';
import { Matrix, multiplyAB, multiplyATB, mulMatrix, outerProduct } from './Matrix';
import { trace2x2 } from './Matrix2x2';
import { Polynomial1 } from './Polynomial1';
import { RootsBisection } from './RootsBisection';
import { RootsPolynomial } from './RootsPolynomial';
import type { RootMultiplicity } from './RootsPolynomial';
import { SymmetricEigensolver2x2 } from './SymmetricEigensolver2x2';
import { Vector, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The relationship of the two solid ellipses.
export enum IntrEllipse2Ellipse2Classification {
    ELLIPSES_SEPARATED,
    ELLIPSES_OVERLAP,
    ELLIPSE0_OUTSIDE_ELLIPSE1_BUT_TANGENT,
    ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1,
    ELLIPSE0_CONTAINS_ELLIPSE1_BUT_TANGENT,
    ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0,
    ELLIPSE1_CONTAINS_ELLIPSE0_BUT_TANGENT,
    ELLIPSES_EQUAL
}

// The maximum number of bisection iterations. Upstream computes
//   3 + std::numeric_limits<T>::digits - std::numeric_limits<T>::min_exponent
// which for 'double' is 3 + 53 - (-1021) = 1077.
const maxBisectionIterations = 1077;

// The query tests the relationship between the ellipses as solid objects.
export class IntrEllipse2Ellipse2TI implements
    TIQuery<Hyperellipsoid, Hyperellipsoid, IntrEllipse2Ellipse2Classification> {

    // The ellipse axes are already normalized, which most likely introduced
    // rounding errors.
    test(ellipse0: Hyperellipsoid, ellipse1: Hyperellipsoid):
        IntrEllipse2Ellipse2Classification {
        logAssert(ellipse0.dimension === 2 && ellipse1.dimension === 2,
            'IntrEllipse2Ellipse2TI: mismatched sizes.');

        const zero = 0, one = 1;

        // Get the parameters of ellipse0.
        const K0 = ellipse0.center;
        const R0 = Matrix.zero(2, 2);
        R0.setCol(0, ellipse0.axis[0]);
        R0.setCol(1, ellipse0.axis[1]);

        // Get the parameters of ellipse1.
        const K1 = ellipse1.center;
        const R1 = Matrix.zero(2, 2);
        R1.setCol(0, ellipse1.axis[0]);
        R1.setCol(1, ellipse1.axis[1]);
        const D1 = Matrix.fromArray(2, 2, [
            one / (ellipse1.extent.values[0] * ellipse1.extent.values[0]), zero,
            zero, one / (ellipse1.extent.values[1] * ellipse1.extent.values[1])
        ]);

        // Compute K2 = D0^{1/2}*R0^T*(K1-K0), where the product
        // R0^T*(K1-K0) is computed as the vector-on-the-left product
        // (K1-K0)*R0.
        const D0NegHalf = Matrix.fromArray(2, 2, [
            ellipse0.extent.values[0], zero,
            zero, ellipse0.extent.values[1]
        ]);

        const D0Half = Matrix.fromArray(2, 2, [
            one / ellipse0.extent.values[0], zero,
            zero, one / ellipse0.extent.values[1]
        ]);

        const K2 = mulMatrix(D0Half, mulMatrix(sub(K1, K0), R0));

        // Compute M2.
        const R1TR0D0NegHalf = multiplyATB(R1, multiplyAB(R0, D0NegHalf));
        const M2 = multiplyAB(multiplyATB(R1TR0D0NegHalf, D1), R1TR0D0NegHalf);

        // Factor M2 = R*D*R^T.
        const es = new SymmetricEigensolver2x2();
        const { evals: D, evecs: evec } =
            es.solve(M2.get(0, 0), M2.get(0, 1), M2.get(1, 1), +1);
        const R = Matrix.zero(2, 2);
        R.setCol(0, Vector.fromArray(evec[0]));
        R.setCol(1, Vector.fromArray(evec[1]));

        // Compute K = R^T*K2.
        const K = mulMatrix(K2, R);

        // Transformed ellipse0 is Z^T*Z = 1 and transformed ellipse1 is
        // (Z-K)^T*D*(Z-K) = 0.

        // The minimum and maximum squared distances from the origin of points
        // on transformed ellipse1 are used to determine whether the ellipses
        // intersect, are separated or one contains the other.
        let minSqrDistance = Number.MAX_VALUE;
        let maxSqrDistance = zero;

        if (K.values[0] === 0 && K.values[1] === 0) {
            // The special case of common centers must be handled separately.
            // It is not possible for the ellipses to be separated.
            for (let i = 0; i < 2; ++i) {
                const invD = one / D[i];
                if (invD < minSqrDistance) {
                    minSqrDistance = invD;
                }
                if (invD > maxSqrDistance) {
                    maxSqrDistance = invD;
                }
            }
            return classify(minSqrDistance, maxSqrDistance, zero);
        }

        // The closest point P0 and farthest point P1 are solutions to
        // s0*D*(P0 - K) = P0 and s1*D1*(P1 - K) = P1 for some scalars s0 and
        // s1 that are roots to the function
        //   f(s) = d0*k0^2/(d0*s-1)^2 + d1*k1^2/(d1*s-1)^2 - 1
        // where D = diagonal(d0,d1) and K = (k0,k1).
        const d0 = D[0], d1 = D[1];
        const c0 = K.values[0] * K.values[0], c1 = K.values[1] * K.values[1];

        // Sort the values so that d0 >= d1. This allows us to bound the roots
        // of f(s), of which there are at most 4.
        const param: Array<[number, number]> = (d0 >= d1
            ? [[d0, c0], [d1, c1]]
            : [[d1, c1], [d0, c0]]);

        const valid: Array<[number, number]> = [];
        if (param[0][0] > param[1][0]) {
            // d0 > d1
            for (let i = 0; i < 2; ++i) {
                if (param[i][1] > zero) {
                    valid.push(param[i]);
                }
            }
        }
        else {
            // d0 = d1
            param[0][1] += param[1][1];
            if (param[0][1] > zero) {
                valid.push(param[0]);
            }
        }

        let roots: number[] = [];
        if (valid.length === 2) {
            roots = getRootsTwoTerms(valid[0][0], valid[1][0], valid[0][1],
                valid[1][1]);
        }
        else if (valid.length === 1) {
            roots = getRootsOneTerm(valid[0][0], valid[0][1]);
        }
        // else: valid.length cannot be zero because the case K = 0 has
        // already been handled.

        for (const s of roots) {
            const p0 = d0 * K.values[0] * s / (d0 * s - 1);
            const p1 = d1 * K.values[1] * s / (d1 * s - 1);
            const sqrDistance = p0 * p0 + p1 * p1;
            if (sqrDistance < minSqrDistance) {
                minSqrDistance = sqrDistance;
            }
            if (sqrDistance > maxSqrDistance) {
                maxSqrDistance = sqrDistance;
            }
        }

        return classify(minSqrDistance, maxSqrDistance, d0 * c0 + d1 * c1);
    }
}

// The roots of f(s) = d0*c0/(d0*s-1)^2 - 1.
function getRootsOneTerm(d0: number, c0: number): number[] {
    const temp = Math.sqrt(d0 * c0);
    const inv = 1 / d0;
    return [(1 - temp) * inv, (1 + temp) * inv];
}

// The roots of f(s) = d0*c0/(d0*s-1)^2 + d1*c1/(d1*s-1)^2 - 1 with d0 > d1.
function getRootsTwoTerms(d0: number, d1: number, c0: number, c1: number):
    number[] {
    const zero = 0, one = 1;
    const d0c0 = d0 * c0;
    const d1c1 = d1 * c1;
    const sum = d0c0 + d1c1;
    const sqrtsum = Math.sqrt(sum);

    const F = (s: number): number => {
        const invN0 = one / (d0 * s - one);
        const invN1 = one / (d1 * s - one);
        const term0 = d0c0 * invN0 * invN0;
        const term1 = d1c1 * invN1 * invN1;
        return term0 + term1 - one;
    };

    const roots: number[] = [];

    const invD0 = one / d0;
    const invD1 = one / d1;
    let smin: number, smax: number, fval: number;

    // Compute the root in (-infinity,1/d0). Obtain a lower bound for the root
    // better than -MAX_VALUE.
    smax = invD0;
    fval = sum - one;
    if (fval > zero) {
        smin = (one - sqrtsum) * invD1;  // < 0
        fval = F(smin);
        logAssert(fval <= zero, 'Unexpected condition.');
    }
    else {
        smin = zero;
    }
    let bisect = RootsBisection.find(F, smin, smax, -one, one,
        maxBisectionIterations);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    // Compute the roots (if any) in (1/d0,1/d1). It is the case that
    //   F(1/d0) = +infinity, F'(1/d0) = -infinity
    //   F(1/d1) = +infinity, F'(1/d1) = +infinity
    //   F"(s) > 0 for all s in the domain of F
    // Compute the unique root r of F'(s) on (1/d0,1/d1). The bisector needs
    // only the signs at the endpoints, so we pass -1 and +1 instead of the
    // infinite values. If F(r) < 0, F(s) has two roots in the interval. If
    // F(r) = 0, F(s) has only one root in the interval.
    const oneThird = 1 / 3;
    const rho = Math.pow(d0 * d0c0 / (d1 * d1c1), oneThird);
    const smid = (one + rho) / (d0 + rho * d1);
    const fmid = F(smid);
    if (fmid < zero) {
        // Pass in signs rather than infinities, because the bisector cares
        // only about the signs.
        bisect = RootsBisection.find(F, invD0, smid, one, -one,
            maxBisectionIterations);
        logAssert(bisect.iterations > 0, 'Unexpected condition.');
        roots.push(bisect.root);
        bisect = RootsBisection.find(F, smid, invD1, -one, one,
            maxBisectionIterations);
        logAssert(bisect.iterations > 0, 'Unexpected condition.');
        roots.push(bisect.root);
    }
    else if (fmid === zero) {
        roots.push(smid);
    }

    // Compute the root in (1/d1,+infinity). Obtain an upper bound for the
    // root better than MAX_VALUE.
    smin = invD1;
    smax = (one + sqrtsum) * invD1;  // > 1/d1
    fval = F(smax);
    logAssert(fval <= zero, 'Unexpected condition.');
    bisect = RootsBisection.find(F, smin, smax, one, -one,
        maxBisectionIterations);
    logAssert(bisect.iterations > 0, 'Unexpected condition.');
    roots.push(bisect.root);

    return roots;
}

function classify(minSqrDistance: number, maxSqrDistance: number,
    d0c0pd1c1: number): IntrEllipse2Ellipse2Classification {
    const C = IntrEllipse2Ellipse2Classification;
    const one = 1;

    if (maxSqrDistance < one) {
        return C.ELLIPSE0_STRICTLY_CONTAINS_ELLIPSE1;
    }
    else if (maxSqrDistance > one) {
        if (minSqrDistance < one) {
            return C.ELLIPSES_OVERLAP;
        }
        else if (minSqrDistance > one) {
            if (d0c0pd1c1 > one) {
                return C.ELLIPSES_SEPARATED;
            }
            else {
                return C.ELLIPSE1_STRICTLY_CONTAINS_ELLIPSE0;
            }
        }
        else {  // minSqrDistance = 1
            if (d0c0pd1c1 > one) {
                return C.ELLIPSE0_OUTSIDE_ELLIPSE1_BUT_TANGENT;
            }
            else {
                return C.ELLIPSE1_CONTAINS_ELLIPSE0_BUT_TANGENT;
            }
        }
    }
    else {  // maxSqrDistance = 1
        if (minSqrDistance < one) {
            return C.ELLIPSE0_CONTAINS_ELLIPSE1_BUT_TANGENT;
        }
        else {  // minSqrDistance = 1
            return C.ELLIPSES_EQUAL;
        }
    }
}

// The value of 'numPoints' when the ellipses are the same, in which case the
// intersection has infinitely many points and 'points' is invalid. Upstream
// uses std::numeric_limits<size_t>::max().
export const intrEllipse2Ellipse2InfinitePoints = Number.MAX_SAFE_INTEGER;

// The result of the IntrEllipse2Ellipse2FI queries.
export interface IntrEllipse2Ellipse2FIResult {
    // This value is true when the ellipses intersect in at least one point.
    intersect: boolean;

    // If the ellipses are not the same, numPoints is 0 through 4 and that
    // number of elements of 'points' are valid. If the ellipses are the same,
    // numPoints is intrEllipse2Ellipse2InfinitePoints and 'points' is invalid
    // (set to zero-valued vectors).
    numPoints: number;
    points: [Vector, Vector, Vector, Vector];
    isTransverse: [boolean, boolean, boolean, boolean];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrEllipse2Ellipse2FIResult():
    IntrEllipse2Ellipse2FIResult {
    return {
        intersect: false,
        numPoints: 0,
        points: [Vector.zero(2), Vector.zero(2), Vector.zero(2), Vector.zero(2)],
        isTransverse: [false, false, false, false]
    };
}

// The queries find the intersections (if any) of the ellipses treated as
// hollow objects.
export class IntrEllipse2Ellipse2FI {
    // Compute the standard form (C,M) of the ellipse, where the ellipse
    // points X satisfy (X-C)^T*M*(X-C) = 1.
    getStandardForm(ellipse: Hyperellipsoid): { C: Vector, M: Matrix } {
        logAssert(ellipse.dimension === 2,
            'IntrEllipse2Ellipse2FI: mismatched sizes.');
        const UUTrn = outerProduct(ellipse.axis[0], ellipse.axis[0]);
        const VVTrn = outerProduct(ellipse.axis[1], ellipse.axis[1]);
        const USqrLen = trace2x2(UUTrn);
        const aSqr = ellipse.extent.values[0] * ellipse.extent.values[0];
        const bSqr = ellipse.extent.values[1] * ellipse.extent.values[1];
        const C = ellipse.center.clone();
        const M = Matrix.zero(2, 2);
        for (let r = 0; r < 2; ++r) {
            for (let c = 0; c < 2; ++c) {
                M.set(r, c,
                    (UUTrn.get(r, c) / aSqr + VVTrn.get(r, c) / bSqr) / USqrLen);
            }
        }
        return { C, M };
    }

    // The axis-aligned bounding box of the ellipse.
    computeAlignedBox(ellipse: Hyperellipsoid): AlignedBox;
    computeAlignedBox(C: Vector, M: Matrix): AlignedBox;
    computeAlignedBox(arg0: Hyperellipsoid | Vector, arg1?: Matrix): AlignedBox {
        if (arg1 === undefined) {
            const { C, M } = this.getStandardForm(arg0 as Hyperellipsoid);
            return this.computeAlignedBox(C, M);
        }

        const C = arg0 as Vector;
        const M = arg1;
        const determinant = M.get(0, 0) * M.get(1, 1) - M.get(0, 1) * M.get(0, 1);
        const distance: [number, number] = [
            Math.sqrt(M.get(1, 1) / determinant),
            Math.sqrt(M.get(0, 0) / determinant)
        ];

        const box = new AlignedBox(2);
        for (let i = 0; i < 2; ++i) {
            box.min.values[i] = C.values[i] - distance[i];
            box.max.values[i] = C.values[i] + distance[i];
        }
        return box;
    }

    // The query for the ellipses in standard form (center, matrix).
    findStandardForm(C0: Vector, M0: Matrix, C1: Vector, M1: Matrix,
        useEarlyExitNoIntersectionTest: boolean = true):
        IntrEllipse2Ellipse2FIResult {
        const result = defaultIntrEllipse2Ellipse2FIResult();

        // Test whether the ellipses are the same. If so, report that there
        // are infinitely many points of intersection.
        if (C0.equals(C1) && M0.equals(M1)) {
            result.numPoints = intrEllipse2Ellipse2InfinitePoints;
            return result;
        }

        if (useEarlyExitNoIntersectionTest) {
            // Test whether the axis-aligned bounding boxes are disjoint. If
            // so, the ellipses do not intersect.
            const box0 = this.computeAlignedBox(C0, M0);
            const box1 = this.computeAlignedBox(C1, M1);
            for (let i = 0; i < 2; ++i) {
                if (box0.max.values[i] < box1.min.values[i]
                    || box0.min.values[i] > box1.max.values[i]) {
                    // The member result.intersect is already false.
                    return result;
                }
            }
        }

        const zero = 0, one = 1, two = 2;

        const ell = M0.get(0, 1) / M0.get(0, 0);
        const d0 = M0.get(0, 0);
        const d1 = robustDOP(M0.get(0, 0), M0.get(1, 1), M0.get(0, 1),
            M0.get(0, 1)) / M0.get(0, 0);
        const k0 = C1.values[0] - C0.values[0];
        const k1 = C1.values[1] - C0.values[1];
        const term0 = robustSOP(k0, M1.get(0, 0), k1, M1.get(0, 1));
        const term1 = robustSOP(k0, M1.get(0, 1), k1, M1.get(1, 1));
        const g0 = robustSOP(k0, term0, k1, term1) - one;
        const g1 = -two * term0;
        const g2 = two * fma(term0, ell, -term1);
        const g3 = M1.get(0, 0);
        const g4 = -two * fma(M1.get(0, 0), ell, -M1.get(0, 1));
        const g5 = fma(-ell, robustDOP(two, M1.get(0, 1), ell, M1.get(0, 0)),
            M1.get(1, 1));
        const e0 = fma(d1, g0, g5);
        const e1 = d1 * g1;
        const e2 = d1 * g2;
        const e3 = robustDOP(d1, g3, d0, g5);
        const e4 = d1 * g4;

        if (e4 !== zero) {
            caseE4NotZero(C0, ell, d0, d1, e0, e1, e2, e3, e4, result);
        }
        else {
            if (e2 !== zero) {
                if (e3 !== zero) {
                    caseE4ZeroE2NotZeroE3NotZero(C0, ell, d0, d1, e0, e1, e2,
                        e3, result);
                }
                else {
                    caseE4ZeroE2NotZeroE3Zero(C0, ell, d0, d1, e0, e1, e2,
                        result);
                }
            }
            else {
                if (e3 !== zero) {
                    caseE4ZeroE2ZeroE3NotZero(C0, ell, d0, d1, e0, e1, e3,
                        result);
                }
                else if (e1 !== zero) {
                    caseE4ZeroE2ZeroE3Zero(C0, ell, d0, d1, e0, e1, result);
                }
                // else: The ellipses are axis-aligned and have the same
                // center. The extent vectors are parallel but not equal. One
                // ellipse is strictly inside the other, so there is no
                // intersection.
            }
        }

        return result;
    }

    // The canonical query for two ellipses.
    find(ellipse0: Hyperellipsoid, ellipse1: Hyperellipsoid,
        useEarlyExitNoIntersectionTest: boolean = true):
        IntrEllipse2Ellipse2FIResult {
        const sf0 = this.getStandardForm(ellipse0);
        const sf1 = this.getStandardForm(ellipse1);
        return this.findStandardForm(sf0.C, sf0.M, sf1.C, sf1.M,
            useEarlyExitNoIntersectionTest);
    }
}

// Read the coefficient of t^i of the polynomial, where a polynomial whose
// leading coefficients cancelled has a smaller degree than expected. Upstream
// indexes the coefficient array directly, which is out-of-range access when
// the degree is smaller; that cannot happen for the polynomials of this file,
// because the leading coefficients are positive by construction.
function coeff(p: Polynomial1, i: number): number {
    return (i <= p.getDegree() ? p.get(i) : 0);
}

function caseE4ZeroE2ZeroE3NotZero(C0: Vector, ell: number, d0: number,
    d1: number, e0: number, e1: number, e3: number,
    result: IntrEllipse2Ellipse2FIResult): void {
    const rmMap: RootMultiplicity[] = RootsPolynomial.solveQuadratic(e0, e1, e3);
    for (const rm of rmMap) {
        const y0 = rm.root;
        const lambda = fma(-d0, y0 * y0, 1);
        if (lambda < 0) {
            continue;
        }

        if (lambda > 0) {
            let y1 = -Math.sqrt(lambda / d1);
            result.points[result.numPoints] =
                Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                    y1 + C0.values[1]]);
            result.isTransverse[result.numPoints] = (rm.multiplicity === 1);
            ++result.numPoints;
            y1 = -y1;
            result.points[result.numPoints] =
                Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                    y1 + C0.values[1]]);
            result.isTransverse[result.numPoints] = (rm.multiplicity === 1);
            ++result.numPoints;
        }
        else {
            result.points[result.numPoints] =
                Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
            result.isTransverse[result.numPoints] = false;
            ++result.numPoints;
        }

        result.intersect = true;
    }
}

function caseE4ZeroE2ZeroE3Zero(C0: Vector, ell: number, d0: number,
    d1: number, e0: number, e1: number,
    result: IntrEllipse2Ellipse2FIResult): void {
    const y0 = -e0 / e1;
    const lambda = fma(-d0, y0 * y0, 1);
    if (lambda < 0) {
        return;
    }

    if (lambda > 0) {
        let y1 = -Math.sqrt(lambda / d1);
        result.points[result.numPoints] =
            Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                y1 + C0.values[1]]);
        result.isTransverse[result.numPoints] = true;
        ++result.numPoints;
        y1 = -y1;
        result.points[result.numPoints] =
            Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                y1 + C0.values[1]]);
        result.isTransverse[result.numPoints] = true;
        ++result.numPoints;
    }
    else {
        result.points[result.numPoints] =
            Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
        result.isTransverse[result.numPoints] = false;
        ++result.numPoints;
    }

    result.intersect = true;
}

function caseE4ZeroE2NotZeroE3Zero(C0: Vector, ell: number, d0: number,
    d1: number, e0: number, e1: number, e2: number,
    result: IntrEllipse2Ellipse2FIResult): void {
    const poly0 = Polynomial1.fromCoefficients([-1, 0, d0]);
    const poly1 = Polynomial1.fromCoefficients([e0, e1]);
    const H = poly0.mul(e2 * e2).add(poly1.mul(poly1).mul(d1));
    const rmMap: RootMultiplicity[] = RootsPolynomial.solveQuadratic(
        coeff(H, 0), coeff(H, 1), coeff(H, 2));
    for (const rm of rmMap) {
        const y0 = rm.root;
        const lambda = fma(-d0, y0 * y0, 1);
        if (lambda < 0) {
            continue;
        }

        if (lambda > 0) {
            // Choose the y1-root with smallest |(e0 + e1 * y0) + (e2) * y1|.
            const y1cand0 = -Math.sqrt(lambda / d1);
            const test0 = Math.abs(e0 + robustSOP(e1, y0, e2, y1cand0));
            const y1cand1 = -y1cand0;
            const test1 = Math.abs(e0 + robustSOP(e1, y0, e2, y1cand1));
            const y1 = (test0 <= test1 ? y1cand0 : y1cand1);
            result.points[result.numPoints] =
                Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                    y1 + C0.values[1]]);
        }
        else {
            result.points[result.numPoints] =
                Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
        }

        result.isTransverse[result.numPoints] = (rm.multiplicity === 1);
        ++result.numPoints;

        result.intersect = true;
    }
}

function caseE4ZeroE2NotZeroE3NotZero(C0: Vector, ell: number, d0: number,
    d1: number, e0: number, e1: number, e2: number, e3: number,
    result: IntrEllipse2Ellipse2FIResult): void {
    const poly0 = Polynomial1.fromCoefficients([-1, 0, d0]);
    const poly1 = Polynomial1.fromCoefficients([e0, e1, e3]);
    const H = poly0.mul(e2 * e2).add(poly1.mul(poly1).mul(d1));
    const rmMap: RootMultiplicity[] = RootsPolynomial.solveQuartic(
        coeff(H, 0), coeff(H, 1), coeff(H, 2), coeff(H, 3), coeff(H, 4));
    for (const rm of rmMap) {
        const y0 = rm.root;
        const lambda = fma(-d0, y0 * y0, 1);
        if (lambda < 0) {
            continue;
        }

        if (lambda > 0) {
            // Choose the y1-root with smallest
            // |(e0 + e1 * y0 + e3 * y0^2) + (e2) * y1|.
            const term0 = fma(e3, y0, e1);
            const term1 = fma(term0, y0, e0);
            const y1cand0 = -Math.sqrt(lambda / d1);
            const test0 = Math.abs(fma(e2, y1cand0, term1));
            const y1cand1 = -y1cand0;
            const test1 = Math.abs(fma(e2, y1cand1, term1));
            const y1 = (test0 < test1 ? y1cand0 : y1cand1);
            result.points[result.numPoints] =
                Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                    y1 + C0.values[1]]);
        }
        else {
            result.points[result.numPoints] =
                Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
        }

        result.isTransverse[result.numPoints] = (rm.multiplicity === 1);
        ++result.numPoints;

        result.intersect = true;
    }
}

function caseE4NotZero(C0: Vector, ell: number, d0: number, d1: number,
    e0: number, e1: number, e2: number, e3: number, e4: number,
    result: IntrEllipse2Ellipse2FIResult): void {
    const poly0 = Polynomial1.fromCoefficients([-1, 0, d0]);
    const poly1 = Polynomial1.fromCoefficients([e0, e1, e3]);
    const poly2 = Polynomial1.fromCoefficients([e2, e4]);
    const H = poly2.mul(poly2).mul(poly0).add(poly1.mul(poly1).mul(d1));
    const rmMap: RootMultiplicity[] = RootsPolynomial.solveQuartic(
        coeff(H, 0), coeff(H, 1), coeff(H, 2), coeff(H, 3), coeff(H, 4));
    for (const rm of rmMap) {
        const y0 = rm.root;
        const lambda = fma(-d0, y0 * y0, 1);
        if (lambda < 0) {
            continue;
        }

        const divisor = e2 + e4 * y0;
        if (divisor !== 0) {
            if (lambda > 0) {
                // Choose the y1-root with smallest
                // |(e0 + e1 * y0 + e3 * y0^2) + (e2 + e4 * y0) * y1|.
                const y1cand0 = -Math.sqrt(lambda / d1);
                const term0 = fma(e3, y0, e1);
                const term1 = fma(term0, y0, e0);
                const test0 = Math.abs(fma(divisor, y1cand0, term1));
                const y1cand1 = -y1cand0;
                const test1 = Math.abs(fma(divisor, y1cand1, term1));
                const y1 = (test0 < test1 ? y1cand0 : y1cand1);
                result.points[result.numPoints] =
                    Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                        y1 + C0.values[1]]);
            }
            else {
                result.points[result.numPoints] =
                    Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
            }

            result.isTransverse[result.numPoints] = (rm.multiplicity === 1);
            ++result.numPoints;

            result.intersect = true;
        }
        else {
            // Upstream bug (fixed here): in this branch upstream writes both
            // of the symmetric points to result.points[result.numPoints]
            // without incrementing numPoints between the two writes, so the
            // first point is overwritten by the second and only one of the
            // two intersection points is reported. The port stores both
            // points, each with its own isTransverse flag (which upstream
            // sets to lambda > 0).
            if (lambda > 0) {
                let y1 = -Math.sqrt(lambda / d1);
                result.points[result.numPoints] =
                    Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                        y1 + C0.values[1]]);
                result.isTransverse[result.numPoints] = true;
                ++result.numPoints;
                y1 = -y1;
                result.points[result.numPoints] =
                    Vector.fromArray([fma(-ell, y1, y0) + C0.values[0],
                        y1 + C0.values[1]]);
                result.isTransverse[result.numPoints] = true;
                ++result.numPoints;
            }
            else {
                result.points[result.numPoints] =
                    Vector.fromArray([y0 + C0.values[0], C0.values[1]]);
                result.isTransverse[result.numPoints] = false;
                ++result.numPoints;
            }

            result.intersect = true;
        }
    }
}
