// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsPolynomial.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The find functions return the roots, if any. If the polynomial is
// identically zero, find returns the single root 0.
//
// Some root-bounding algorithms for real-valued roots are mentioned next for
// the polynomial p(t) = c[0] + c[1]*t + ... + c[d-1]*t^{d-1} + c[d]*t^d.
//
// 1. The roots must be contained by the interval [-M,M] where
//   M = 1 + max{|c[0]|, ..., |c[d-1]|}/|c[d]| >= 1
// is called the Cauchy bound.
//
// 2. You may search for roots in the interval [-1,1]. Define
//   q(t) = t^d*p(1/t) = c[0]*t^d + c[1]*t^{d-1} + ... + c[d-1]*t + c[d]
// The roots of p(t) not in [-1,1] are the roots of q(t) in [-1,1].
//
// 3. Between two consecutive roots of the derivative p'(t), say, r0 < r1,
// the function p(t) is strictly monotonic on the open interval (r0,r1).
// If additionally, p(r0) * p(r1) <= 0, then p(x) has a unique root on the
// closed interval [r0,r1]. Thus, one can compute the derivatives through
// order d for p(t), find roots for the derivative of order k+1, then use
// these to bound roots for the derivative of order k.
//
// 4. Sturm sequences of polynomials may be used to determine bounds on the
// roots. This is a more sophisticated approach to root bounding than item 3.
// Moreover, a Sturm sequence allows you to compute the number of real-valued
// roots on a specified interval.
//
// 5. For the low-degree solve* functions, see
// https://www.geometrictools.com/Documentation/LowDegreePolynomialRoots.pdf
//
// Port notes:
// - Upstream templates the Solve*/GetRootInfo* functions on a Rational type
//   for exact arithmetic in the root classification (sign tests on the
//   discriminant-like quantities). This port instantiates the floating-point
//   path (Rational = number, IEEE double). The root classification is
//   therefore subject to rounding: polynomials whose classifier quantities
//   (delta, a0, a1, c0, c1, ...) are not exactly representable may be
//   misclassified near degenerate configurations (e.g. a double root may be
//   reported as two nearby simple roots or as a complex pair). For exact
//   classification, use rational arithmetic as upstream suggests
//   (BSRational); that path is not ported here.
// - Upstream accumulates roots in std::map<Real, int32_t>, whose iteration
//   is ordered by key. The port replicates this explicitly with an array of
//   { root, multiplicity } entries kept sorted ascending by root, and with
//   std::map::insert semantics (inserting an existing key is a no-op).
// - Upstream's GTE_ROOTS_LOW_DEGREE_BLOCK unit-testing hook is ported as the
//   optional static callback rootsLowDegreeBlock.
//
// NOTE: Upstream RootsPolynomial is deprecated in favor of RootsLinear,
// RootsQuadratic, RootsCubic, RootsQuartic and RootsGeneralPolynomial. It is
// ported for parity and for dependent files.

export interface RootMultiplicity {
    root: number;
    multiplicity: number;
}

// The port of std::map<Real, int32_t>: entries sorted ascending by root.
type RootMap = RootMultiplicity[];

// std::map::insert semantics: if the key already exists, do nothing.
function rmInsert(rmMap: RootMap, root: number, multiplicity: number): void {
    for (let i = 0; i < rmMap.length; ++i) {
        if (root === rmMap[i].root) {
            return;
        }
        if (root < rmMap[i].root) {
            rmMap.splice(i, 0, { root, multiplicity });
            return;
        }
    }
    rmMap.push({ root, multiplicity });
}

// std::map::find semantics: return the entry with the given key or null.
function rmFind(rmMap: RootMap, root: number): RootMultiplicity | null {
    for (const rm of rmMap) {
        if (rm.root === root) {
            return rm;
        }
    }
    return null;
}

export class RootsPolynomial {
    // FOR INTERNAL USE ONLY (unit testing). The port of upstream's
    // GTE_ROOTS_LOW_DEGREE_BLOCK macro: when non-null, the callback is
    // invoked with the block number of each classification branch taken by
    // the low-degree solvers.
    static rootsLowDegreeBlock: ((block: number) => void) | null = null;

    private static lowDegreeBlock(block: number): void {
        if (RootsPolynomial.rootsLowDegreeBlock !== null) {
            RootsPolynomial.rootsLowDegreeBlock(block);
        }
    }

    // Low-degree root finders. Upstream uses exact rational arithmetic for
    // theoretically correct root classification; this port uses IEEE double
    // (see the port notes above). The returned array contains the distinct
    // real roots with their multiplicities, sorted ascending by root (the
    // iteration order of upstream's std::map<Real, int32_t>). The
    // highest-order coefficients must be nonzero (p2 != 0 for quadratic,
    // p3 != 0 for cubic, and p4 != 0 for quartic).

    static solveQuadratic(p0: number, p1: number, p2: number): RootMultiplicity[] {
        const q0 = p0 / p2;
        const q1 = p1 / p2;
        const q1half = q1 / 2;
        const c0 = q0 - q1half * q1half;

        const rmLocalMap: RootMap = [];
        RootsPolynomial.solveDepressedQuadratic(c0, rmLocalMap);

        const rmMap: RootMap = [];
        for (const rm of rmLocalMap) {
            const root = rm.root - q1half;
            rmInsert(rmMap, root, rm.multiplicity);
        }
        return rmMap;
    }

    static solveCubic(p0: number, p1: number, p2: number, p3: number): RootMultiplicity[] {
        const q0 = p0 / p3;
        const q1 = p1 / p3;
        const q2 = p2 / p3;
        const q2third = q2 / 3;
        const c0 = q0 - q2third * (q1 - 2 * q2third * q2third);
        const c1 = q1 - q2 * q2third;

        const rmLocalMap: RootMap = [];
        RootsPolynomial.solveDepressedCubic(c0, c1, rmLocalMap);

        const rmMap: RootMap = [];
        for (const rm of rmLocalMap) {
            const root = rm.root - q2third;
            rmInsert(rmMap, root, rm.multiplicity);
        }
        return rmMap;
    }

    static solveQuartic(p0: number, p1: number, p2: number, p3: number,
        p4: number): RootMultiplicity[] {
        const q0 = p0 / p4;
        const q1 = p1 / p4;
        const q2 = p2 / p4;
        const q3 = p3 / p4;
        const q3fourth = q3 / 4;
        const q3fourthSqr = q3fourth * q3fourth;
        const c0 = q0 - q3fourth * (q1 - q3fourth * (q2 - q3fourthSqr * 3));
        const c1 = q1 - 2 * q3fourth * (q2 - 4 * q3fourthSqr);
        const c2 = q2 - 6 * q3fourthSqr;

        const rmLocalMap: RootMap = [];
        RootsPolynomial.solveDepressedQuartic(c0, c1, c2, rmLocalMap);

        const rmMap: RootMap = [];
        for (const rm of rmLocalMap) {
            const root = rm.root - q3fourth;
            rmInsert(rmMap, root, rm.multiplicity);
        }
        return rmMap;
    }

    // Return only the number of real-valued roots and their multiplicities.
    // info.length is the number of real-valued roots and info[i] is the
    // multiplicity of the root corresponding to index i.

    static getRootInfoQuadratic(p0: number, p1: number, p2: number): number[] {
        const q0 = p0 / p2;
        const q1 = p1 / p2;
        const q1half = q1 / 2;
        const c0 = q0 - q1half * q1half;

        const info: number[] = [];
        RootsPolynomial.getRootInfoDepressedQuadratic(c0, info);
        return info;
    }

    static getRootInfoCubic(p0: number, p1: number, p2: number, p3: number): number[] {
        const q0 = p0 / p3;
        const q1 = p1 / p3;
        const q2 = p2 / p3;
        const q2third = q2 / 3;
        const c0 = q0 - q2third * (q1 - 2 * q2third * q2third);
        const c1 = q1 - q2 * q2third;

        const info: number[] = [];
        RootsPolynomial.getRootInfoDepressedCubic(c0, c1, info);
        return info;
    }

    static getRootInfoQuartic(p0: number, p1: number, p2: number, p3: number,
        p4: number): number[] {
        const q0 = p0 / p4;
        const q1 = p1 / p4;
        const q2 = p2 / p4;
        const q3 = p3 / p4;
        const q3fourth = q3 / 4;
        const q3fourthSqr = q3fourth * q3fourth;
        const c0 = q0 - q3fourth * (q1 - q3fourth * (q2 - q3fourthSqr * 3));
        const c1 = q1 - 2 * q3fourth * (q2 - 4 * q3fourthSqr);
        const c2 = q2 - 6 * q3fourthSqr;

        const info: number[] = [];
        RootsPolynomial.getRootInfoDepressedQuartic(c0, c1, c2, info);
        return info;
    }

    // General equations: sum_{i=0}^{d} c(i)*t^i = 0. The input array 'c'
    // must have at least degree+1 elements.

    // Find the roots on (-infinity,+infinity). Returns the array of roots
    // found, in ascending order. NOTE: as in upstream, a root of even
    // multiplicity might be missed (no sign change) and a root of odd
    // multiplicity larger than 1 can be reported multiple times (once per
    // bounding subinterval whose endpoint function value is zero).
    static find(degree: number, c: readonly number[], maxIterations: number): number[];

    // If you know that p(tmin) * p(tmax) <= 0, then there must be at least
    // one root in [tmin, tmax]. Compute it using bisection. The result root
    // is valid only when found is true.
    static find(degree: number, c: readonly number[], tmin: number, tmax: number,
        maxIterations: number): { found: boolean; root: number };

    static find(degree: number, c: readonly number[], arg2: number, arg3?: number,
        arg4?: number): number[] | { found: boolean; root: number } {
        if (arg3 !== undefined && arg4 !== undefined) {
            return RootsPolynomial.findBounded(degree, c, arg2, arg3, arg4);
        }

        const maxIterations = arg2;
        if (degree >= 0 && c.length >= degree + 1) {
            while (degree >= 0 && c[degree] === 0) {
                --degree;
            }

            if (degree > 0) {
                // Compute the Cauchy bound.
                const invLeading = 1 / c[degree];
                let maxValue = 0;
                for (let i = 0; i < degree; ++i) {
                    const value = Math.abs(c[i] * invLeading);
                    if (value > maxValue) {
                        maxValue = value;
                    }
                }
                const bound = 1 + maxValue;

                const roots: number[] = [];
                RootsPolynomial.findRecursive(degree, c, -bound, bound,
                    maxIterations, roots);
                return roots;
            } else if (degree === 0) {
                // The polynomial is a nonzero constant.
                return [];
            } else {
                // The polynomial is identically zero.
                return [0];
            }
        } else {
            // Invalid degree or c.
            return [];
        }
    }

    // Support for the solve* functions.

    private static solveDepressedQuadratic(c0: number, rmMap: RootMap): void {
        if (c0 < 0) {
            // Two simple roots.
            const root1 = Math.sqrt(-c0);
            const root0 = -root1;
            rmInsert(rmMap, root0, 1);
            rmInsert(rmMap, root1, 1);
            RootsPolynomial.lowDegreeBlock(0);
        } else if (c0 === 0) {
            // One double root.
            rmInsert(rmMap, 0, 2);
            RootsPolynomial.lowDegreeBlock(1);
        } else {  // c0 > 0
            // A complex-conjugate pair of roots.
            // Complex z0 = -q1/2 - i*sqrt(c0);
            // Complex z0conj = -q1/2 + i*sqrt(c0);
            RootsPolynomial.lowDegreeBlock(2);
        }
    }

    private static solveDepressedCubic(c0: number, c1: number, rmMap: RootMap): void {
        // Handle the special case of c0 = 0, in which case the polynomial
        // reduces to a depressed quadratic.
        if (c0 === 0) {
            RootsPolynomial.solveDepressedQuadratic(c1, rmMap);
            const iter = rmFind(rmMap, 0);
            if (iter !== null) {
                // The quadratic has a root of zero, so the multiplicity
                // must be increased.
                ++iter.multiplicity;
                RootsPolynomial.lowDegreeBlock(3);
            } else {
                // The quadratic does not have a root of zero. Insert the
                // one for the cubic.
                rmInsert(rmMap, 0, 1);
                RootsPolynomial.lowDegreeBlock(4);
            }
            return;
        }

        // Handle the special case of c0 != 0 and c1 = 0.
        const oneThird = 1 / 3;
        if (c1 === 0) {
            // One simple real root.
            let root0: number;
            if (c0 > 0) {
                root0 = -Math.pow(c0, oneThird);
                RootsPolynomial.lowDegreeBlock(5);
            } else {
                root0 = Math.pow(-c0, oneThird);
                RootsPolynomial.lowDegreeBlock(6);
            }
            rmInsert(rmMap, root0, 1);

            // One complex conjugate pair.
            // Complex z0 = root0*(-1 - i*sqrt(3))/2;
            // Complex z0conj = root0*(-1 + i*sqrt(3))/2;
            return;
        }

        // At this time, c0 != 0 and c1 != 0.
        const delta = -(4 * c1 * c1 * c1 + 27 * c0 * c0);
        if (delta > 0) {
            // Three simple roots.
            const deltaDiv108 = delta / 108;
            const betaRe = -c0 / 2;
            const betaIm = Math.sqrt(deltaDiv108);
            const theta = Math.atan2(betaIm, betaRe);
            const thetaDiv3 = theta / 3;
            const cs = Math.cos(thetaDiv3);
            const sn = Math.sin(thetaDiv3);
            const rhoSqr = betaRe * betaRe + betaIm * betaIm;
            const rhoPowThird = Math.pow(rhoSqr, 1 / 6);
            const temp0 = rhoPowThird * cs;
            const temp1 = rhoPowThird * sn * Math.sqrt(3);
            const root0 = 2 * temp0;
            const root1 = -temp0 - temp1;
            const root2 = -temp0 + temp1;
            rmInsert(rmMap, root0, 1);
            rmInsert(rmMap, root1, 1);
            rmInsert(rmMap, root2, 1);
            RootsPolynomial.lowDegreeBlock(7);
        } else if (delta < 0) {
            // One simple root.
            const deltaDiv108 = delta / 108;
            const temp0 = -c0 / 2;
            const temp1 = Math.sqrt(-deltaDiv108);
            let temp2 = temp0 - temp1;
            let temp3 = temp0 + temp1;
            if (temp2 >= 0) {
                temp2 = Math.pow(temp2, oneThird);
                RootsPolynomial.lowDegreeBlock(8);
            } else {
                temp2 = -Math.pow(-temp2, oneThird);
                RootsPolynomial.lowDegreeBlock(9);
            }
            if (temp3 >= 0) {
                temp3 = Math.pow(temp3, oneThird);
                RootsPolynomial.lowDegreeBlock(10);
            } else {
                temp3 = -Math.pow(-temp3, oneThird);
                RootsPolynomial.lowDegreeBlock(11);
            }
            const root0 = temp2 + temp3;
            rmInsert(rmMap, root0, 1);

            // One complex conjugate pair.
            // Complex z0 = (-root0 - i*sqrt(3*root0*root0+4*c1))/2;
            // Complex z0conj = (-root0 + i*sqrt(3*root0*root0+4*c1))/2;
        } else {  // delta = 0
            // One simple root and one double root.
            const root0 = -3 * c0 / (2 * c1);
            const root1 = -2 * root0;
            rmInsert(rmMap, root0, 2);
            rmInsert(rmMap, root1, 1);
            RootsPolynomial.lowDegreeBlock(12);
        }
    }

    private static solveDepressedQuartic(c0: number, c1: number, c2: number,
        rmMap: RootMap): void {
        // Handle the special case of c0 = 0, in which case the polynomial
        // reduces to a depressed cubic.
        if (c0 === 0) {
            RootsPolynomial.solveDepressedCubic(c1, c2, rmMap);
            const iter = rmFind(rmMap, 0);
            if (iter !== null) {
                // The cubic has a root of zero, so the multiplicity must
                // be increased.
                ++iter.multiplicity;
                RootsPolynomial.lowDegreeBlock(13);
            } else {
                // The cubic does not have a root of zero. Insert the one
                // for the quartic.
                rmInsert(rmMap, 0, 1);
                RootsPolynomial.lowDegreeBlock(14);
            }
            return;
        }

        // Handle the special case of c1 = 0, in which case the quartic is
        // a biquadratic
        //   x^4 + c1*x^2 + c0 = (x^2 + c2/2)^2 + (c0 - c2^2/4)
        if (c1 === 0) {
            RootsPolynomial.solveBiquadratic(c0, c2, rmMap);
            return;
        }

        // At this time, c0 != 0 and c1 != 0, which is a requirement for the
        // general solver that must use a root of a special cubic polynomial.
        const c0sqr = c0 * c0, c1sqr = c1 * c1, c2sqr = c2 * c2;
        const delta = c1sqr * (-27 * c1sqr + 4 * c2 *
            (36 * c0 - c2sqr)) + 16 * c0 * (c2sqr * (c2sqr - 8 * c0) +
            16 * c0sqr);
        const a0 = 12 * c0 + c2sqr;
        const a1 = 4 * c0 - c2sqr;

        if (delta > 0) {
            if (c2 < 0 && a1 < 0) {
                // Four simple real roots.
                const rmCubicMap = RootsPolynomial.solveCubic(
                    c1sqr - 4 * c0 * c2, 8 * c0, 4 * c2, -8);
                const t = rmCubicMap[rmCubicMap.length - 1].root;
                const alphaSqr = 2 * t - c2;
                const alpha = Math.sqrt(Math.max(alphaSqr, 0));
                let sgnC1: number;
                if (c1 > 0) {
                    sgnC1 = 1;
                    RootsPolynomial.lowDegreeBlock(15);
                } else {
                    sgnC1 = -1;
                    RootsPolynomial.lowDegreeBlock(16);
                }
                const arg = t * t - c0;
                const beta = sgnC1 * Math.sqrt(Math.max(arg, 0));
                const D0 = alphaSqr - 4 * (t + beta);
                const sqrtD0 = Math.sqrt(Math.max(D0, 0));
                const D1 = alphaSqr - 4 * (t - beta);
                const sqrtD1 = Math.sqrt(Math.max(D1, 0));
                const root0 = (alpha - sqrtD0) / 2;
                const root1 = (alpha + sqrtD0) / 2;
                const root2 = (-alpha - sqrtD1) / 2;
                const root3 = (-alpha + sqrtD1) / 2;
                rmInsert(rmMap, root0, 1);
                rmInsert(rmMap, root1, 1);
                rmInsert(rmMap, root2, 1);
                rmInsert(rmMap, root3, 1);
            } else {  // c2 >= 0 or a1 >= 0
                // Two complex-conjugate pairs. The values alpha, D0 and D1
                // are those of the if-block.
                // Complex z0 = (alpha - i*sqrt(-D0))/2;
                // Complex z0conj = (alpha + i*sqrt(-D0))/2;
                // Complex z1 = (-alpha - i*sqrt(-D1))/2;
                // Complex z1conj = (-alpha + i*sqrt(-D1))/2;
                RootsPolynomial.lowDegreeBlock(17);
            }
        } else if (delta < 0) {
            // Two simple real roots, one complex-conjugate pair.
            const rmCubicMap = RootsPolynomial.solveCubic(
                c1sqr - 4 * c0 * c2, 8 * c0, 4 * c2, -8);
            const t = rmCubicMap[rmCubicMap.length - 1].root;
            const alphaSqr = 2 * t - c2;
            const alpha = Math.sqrt(Math.max(alphaSqr, 0));
            let sgnC1: number;
            if (c1 > 0) {
                sgnC1 = 1;  // Leads to block 18.
            } else {
                sgnC1 = -1;  // Leads to block 19.
            }
            const arg = t * t - c0;
            const beta = sgnC1 * Math.sqrt(Math.max(arg, 0));
            let root0: number, root1: number;
            if (sgnC1 > 0) {
                const D1 = alphaSqr - 4 * (t - beta);
                const sqrtD1 = Math.sqrt(Math.max(D1, 0));
                root0 = (-alpha - sqrtD1) / 2;
                root1 = (-alpha + sqrtD1) / 2;

                // One complex conjugate pair.
                // Complex z0 = (alpha - i*sqrt(-D0))/2;
                // Complex z0conj = (alpha + i*sqrt(-D0))/2;
                RootsPolynomial.lowDegreeBlock(18);
            } else {
                const D0 = alphaSqr - 4 * (t + beta);
                const sqrtD0 = Math.sqrt(Math.max(D0, 0));
                root0 = (alpha - sqrtD0) / 2;
                root1 = (alpha + sqrtD0) / 2;

                // One complex conjugate pair.
                // Complex z0 = (-alpha - i*sqrt(-D1))/2;
                // Complex z0conj = (-alpha + i*sqrt(-D1))/2;
                RootsPolynomial.lowDegreeBlock(19);
            }
            rmInsert(rmMap, root0, 1);
            rmInsert(rmMap, root1, 1);
        } else {  // delta = 0
            if (a1 > 0 || (c2 > 0 && (a1 !== 0 || c1 !== 0))) {
                // One double real root, one complex-conjugate pair.
                const root0 = -c1 * a0 / (9 * c1sqr - 2 * c2 * a1);
                rmInsert(rmMap, root0, 2);

                // One complex conjugate pair.
                // Complex z0 = -root0 - i*sqrt(c2 + root0^2);
                // Complex z0conj = -root0 + i*sqrt(c2 + root0^2);
                RootsPolynomial.lowDegreeBlock(20);
            } else {
                if (a0 !== 0) {
                    // One double real root, two simple real roots.
                    const root0 = -c1 * a0 / (9 * c1sqr - 2 * c2 * a1);
                    const alpha = 2 * root0;
                    const beta = c2 + 3 * root0 * root0;
                    const discr = alpha * alpha - 4 * beta;
                    const temp1 = Math.sqrt(Math.max(discr, 0));
                    const root1 = (-alpha - temp1) / 2;
                    const root2 = (-alpha + temp1) / 2;
                    rmInsert(rmMap, root0, 2);
                    rmInsert(rmMap, root1, 1);
                    rmInsert(rmMap, root2, 1);
                    RootsPolynomial.lowDegreeBlock(21);
                } else {
                    // One triple real root, one simple real root.
                    const root0 = -3 * c1 / (4 * c2);
                    const root1 = -3 * root0;
                    rmInsert(rmMap, root0, 3);
                    rmInsert(rmMap, root1, 1);
                    RootsPolynomial.lowDegreeBlock(22);
                }
            }
        }
    }

    private static solveBiquadratic(c0: number, c2: number, rmMap: RootMap): void {
        // Solve x^4 + c2*x^2 + c0 = 0. We know that c0 != 0 at the time of
        // the solveBiquadratic call, so x = 0 is not a root. Define
        // u = -c2/2 and v = c2^2/4 - c0 = u^2 - c0. Using the quadratic
        // formula,
        //   x^2 is in { u-sqrt(v), u+sqrt(v) }
        // Computing the square root,
        //   x is in { -sqrt(u-sqrt(v)), sqrt(u-sqrt(v)),
        //             -sqrt(u+sqrt(v)), sqrt(u+sqrt(v)) }
        // Because we know c0 != 0, which implies 0 is not a root, it must
        // be that u-sqrt(v) != 0 and u+sqrt(v) != 0.
        //
        // v > 0, u-sqrt(v) > 0 [implies u+sqrt(v) > 0]: (block 23)
        //   Four real roots: r0, -r0, r1, -r1
        //     r0 = sqrt(u-sqrt(v))
        //     r1 = sqrt(u+sqrt(v))
        //
        // v > 0, u+sqrt(v) < 0 [implies u-sqrt(v) < 0]: (block 24)
        //   Two complex conjugate pairs: z0, conj(z0), -z1, -conj(z1)
        //     z0 = sqrt(-u+sqrt(v)) * i
        //     z1 = sqrt(-u-sqrt(v)) * i
        //
        // v > 0, u-sqrt(v) < 0, u+sqrt(v) > 0: (block 25)
        //   Two real roots, one complex conjugate pair: r0, -r0, z0, conj(z0)
        //     r0 = sqrt(u+sqrt(v))
        //     z0 = sqrt(-u+sqrt(v)) * i
        //
        // v < 0: (block 26)
        //   Two complex conjugate pairs: z0, conj(z0), -z0, -conj(z0)
        //     z0 = sqrt((u+sqrt(u^2-v))/2) - sqrt((-u+sqrt(u^2-v))/2) * i
        //        = sqrt((-c2/2+sqrt(c0))/2) - sqrt((c2/2+sqrt(c0))/2) * i
        //
        // v = 0, u > 0: (block 27)
        //   Two real roots, each of multiplicity 2: r0, -r0
        //     r0 = sqrt(u) = sqrt(-c2/2)
        //
        // v = 0, u < 0: (block 28)
        //   Two complex conjugate pairs: z0, conj(z0), -z0, -conj(z0)
        //     z0 = sqrt(-u) * i = sqrt(c2/2) * i

        const u = c2 / -2;
        const v = u * u - c0;
        if (v > 0) {
            const sqrtv = Math.sqrt(v);
            const upsqrtv = u + sqrtv;
            // Compute u - sqrt(v) = c0 / (u + sqrt(v)) to avoid subtractive
            // cancellation.
            const umsqrtv = c0 / upsqrtv;
            if (umsqrtv > 0) {
                // Real roots: r0, -r0, r1, -r1
                // r0 = sqrt(u-sqrt(v))
                // r1 = sqrt(u+sqrt(v))
                const r0 = Math.sqrt(umsqrtv);
                const r1 = Math.sqrt(upsqrtv);
                rmInsert(rmMap, r0, 1);
                rmInsert(rmMap, -r0, 1);
                rmInsert(rmMap, r1, 1);
                rmInsert(rmMap, -r1, 1);
                RootsPolynomial.lowDegreeBlock(23);
            } else if (upsqrtv < 0) {
                // Complex roots: z0, conj(z0), -z1, -conj(z1)
                // z0 = sqrt(-u+sqrt(v)) * i
                // z1 = sqrt(-u-sqrt(v)) * i
                RootsPolynomial.lowDegreeBlock(24);
            } else {  // umsqrtv < 0 and upsqrtv > 0
                // Real roots: r0, -r0
                // Complex roots: z0, conj(z0)
                // r0 = sqrt(u+sqrt(v))
                // z0 = sqrt(-u+sqrt(v)) * i
                const r0 = Math.sqrt(upsqrtv);
                rmInsert(rmMap, r0, 1);
                rmInsert(rmMap, -r0, 1);
                RootsPolynomial.lowDegreeBlock(25);
            }
        } else if (v < 0) {
            // Complex roots: z0, conj(z0), -z0, -conj(z0)
            // z0 = sqrt((u+sqrt(u^2-v))/2)
            //      - sqrt((-u+sqrt(u^2-v))/2) * i
            RootsPolynomial.lowDegreeBlock(26);
        } else {  // v = 0
            if (u > 0) {
                // Real roots: r0, r0, -r0, -r0
                // r0 = sqrt(u)
                const r0 = Math.sqrt(u);
                rmInsert(rmMap, r0, 2);
                rmInsert(rmMap, -r0, 2);
                RootsPolynomial.lowDegreeBlock(27);
            } else {  // u < 0
                // Complex roots: z0, conj(z0), z0, conj(z0)
                // z0 = sqrt(-u) * i
                RootsPolynomial.lowDegreeBlock(28);
            }
        }
    }

    // Support for the getRootInfo* functions.

    private static getRootInfoDepressedQuadratic(c0: number, info: number[]): void {
        if (c0 < 0) {
            // Two simple roots.
            info.push(1);
            info.push(1);
        } else if (c0 === 0) {
            // One double root.
            info.push(2);  // root is zero
        } else {  // c0 > 0
            // A complex-conjugate pair of roots.
        }
    }

    private static getRootInfoDepressedCubic(c0: number, c1: number,
        info: number[]): void {
        // Handle the special case of c0 = 0, in which case the polynomial
        // reduces to a depressed quadratic.
        if (c0 === 0) {
            if (c1 === 0) {
                info.push(3);  // triple root of zero
            } else {
                info.push(1);  // simple root of zero
                RootsPolynomial.getRootInfoDepressedQuadratic(c1, info);
            }
            return;
        }

        const delta = -(4 * c1 * c1 * c1 + 27 * c0 * c0);
        if (delta > 0) {
            // Three simple real roots.
            info.push(1);
            info.push(1);
            info.push(1);
        } else if (delta < 0) {
            // One simple real root.
            info.push(1);
        } else {  // delta = 0
            // One simple real root and one double real root.
            info.push(1);
            info.push(2);
        }
    }

    private static getRootInfoDepressedQuartic(c0: number, c1: number, c2: number,
        info: number[]): void {
        // Handle the special case of c0 = 0, in which case the polynomial
        // reduces to a depressed cubic.
        if (c0 === 0) {
            if (c1 === 0) {
                if (c2 === 0) {
                    info.push(4);  // quadruple root of zero
                } else {
                    info.push(2);  // double root of zero
                    RootsPolynomial.getRootInfoDepressedQuadratic(c2, info);
                }
            } else {
                info.push(1);  // simple root of zero
                RootsPolynomial.getRootInfoDepressedCubic(c1, c2, info);
            }
            return;
        }

        // Handle the special case of c1 = 0, in which case the quartic is
        // a biquadratic
        //   x^4 + c1*x^2 + c0 = (x^2 + c2/2)^2 + (c0 - c2^2/4)
        if (c1 === 0) {
            RootsPolynomial.getRootInfoBiquadratic(c0, c2, info);
            return;
        }

        // At this time, c0 != 0 and c1 != 0, which is a requirement for the
        // general solver that must use a root of a special cubic polynomial.
        const c0sqr = c0 * c0, c1sqr = c1 * c1, c2sqr = c2 * c2;
        const delta = c1sqr * (-27 * c1sqr + 4 * c2 *
            (36 * c0 - c2sqr)) + 16 * c0 * (c2sqr * (c2sqr - 8 * c0) +
            16 * c0sqr);
        const a0 = 12 * c0 + c2sqr;
        const a1 = 4 * c0 - c2sqr;

        if (delta > 0) {
            if (c2 < 0 && a1 < 0) {
                // Four simple real roots.
                info.push(1);
                info.push(1);
                info.push(1);
                info.push(1);
            } else {  // c2 >= 0 or a1 >= 0
                // Two complex-conjugate pairs.
            }
        } else if (delta < 0) {
            // Two simple real roots, one complex-conjugate pair.
            info.push(1);
            info.push(1);
        } else {  // delta = 0
            if (a1 > 0 || (c2 > 0 && (a1 !== 0 || c1 !== 0))) {
                // One double real root, one complex-conjugate pair.
                info.push(2);
            } else {
                if (a0 !== 0) {
                    // One double real root, two simple real roots.
                    info.push(2);
                    info.push(1);
                    info.push(1);
                } else {
                    // One triple real root, one simple real root.
                    info.push(3);
                    info.push(1);
                }
            }
        }
    }

    private static getRootInfoBiquadratic(c0: number, c2: number,
        info: number[]): void {
        const u = c2 / -2;
        const v = u * u - c0;
        if (v > 0) {
            const sqrtv = Math.sqrt(v);
            const upsqrtv = u + sqrtv;
            const umsqrtv = c0 / upsqrtv;
            if (umsqrtv > 0) {
                // Four simple roots.
                info.push(1);
                info.push(1);
                info.push(1);
                info.push(1);
            } else if (upsqrtv < 0) {
                // Two simple complex conjugate pairs.
            } else {  // umsqrtv < 0 and upsqrtv > 0
                // Two simple real roots, one complex conjugate pair.
                info.push(1);
                info.push(1);
            }
        } else if (v < 0) {
            // Two simple complex conjugate pairs.
        } else {  // v = 0
            if (u > 0) {
                // Two double real roots.
                info.push(2);
                info.push(2);
            } else {  // u < 0
                // Double complex conjugate pairs.
            }
        }
    }

    // Support for the find functions. The port of the upstream bounded
    // bisection Find overload.
    private static findBounded(degree: number, c: readonly number[], tmin: number,
        tmax: number, maxIterations: number): { found: boolean; root: number } {
        let pmin = RootsPolynomial.evaluate(degree, c, tmin);
        if (pmin === 0) {
            return { found: true, root: tmin };
        }
        let pmax = RootsPolynomial.evaluate(degree, c, tmax);
        if (pmax === 0) {
            return { found: true, root: tmax };
        }

        if (pmin * pmax > 0) {
            // It is not known whether the interval bounds a root.
            return { found: false, root: 0 };
        }

        if (tmin >= tmax) {
            // Invalid ordering of interval endpoints.
            return { found: false, root: 0 };
        }

        if (maxIterations === 0) {
            // It is expected that the caller set maxIterations to a
            // positive number.
            return { found: false, root: 0 };
        }

        let root = 0;
        for (let i = 1; i <= maxIterations; ++i) {
            root = 0.5 * (tmin + tmax);

            // This test is designed for the case when tmin and tmax are
            // consecutive floating-point numbers.
            if (root === tmin || root === tmax) {
                break;
            }

            const p = RootsPolynomial.evaluate(degree, c, root);
            const product = p * pmin;
            if (product < 0) {
                tmax = root;
                pmax = p;
            } else if (product > 0) {
                tmin = root;
                pmin = p;
            } else {
                break;
            }
        }

        return { found: true, root };
    }

    private static findRecursive(degree: number, c: readonly number[], tmin: number,
        tmax: number, maxIterations: number, roots: number[]): number {
        // The base of the recursion.
        if (degree === 1) {
            let root = 0;
            let numRoots: number;
            if (c[1] !== 0) {
                root = -c[0] / c[1];
                numRoots = 1;
            } else if (c[0] === 0) {
                root = 0;
                numRoots = 1;
            } else {
                numRoots = 0;
            }

            if (numRoots > 0 && tmin <= root && root <= tmax) {
                roots.push(root);
                return 1;
            }
            return 0;
        }

        // Find the roots of the derivative polynomial scaled by 1/degree.
        // The scaling avoids the factorial growth in the coefficients; for
        // example, without the scaling, the high-order term x^d becomes
        // (d!)*x through multiple differentiations. With the scaling we
        // instead get x. This leads to better numerical behavior of the
        // root finder.
        const derivDegree = degree - 1;
        const derivCoeff = new Array<number>(derivDegree + 1);
        for (let i = 0, ip1 = 1; i <= derivDegree; ++i, ++ip1) {
            derivCoeff[i] = c[ip1] * ip1 / degree;
        }
        const derivRoots: number[] = [];
        const numDerivRoots = RootsPolynomial.findRecursive(degree - 1, derivCoeff,
            tmin, tmax, maxIterations, derivRoots);

        let numRoots = 0;
        if (numDerivRoots > 0) {
            // Find root on [tmin,derivRoots[0]].
            let result = RootsPolynomial.findBounded(degree, c, tmin, derivRoots[0],
                maxIterations);
            if (result.found) {
                roots.push(result.root);
                ++numRoots;
            }

            // Find root on [derivRoots[i],derivRoots[i+1]].
            for (let i = 0, ip1 = 1; i <= numDerivRoots - 2; ++i, ++ip1) {
                result = RootsPolynomial.findBounded(degree, c, derivRoots[i],
                    derivRoots[ip1], maxIterations);
                if (result.found) {
                    roots.push(result.root);
                    ++numRoots;
                }
            }

            // Find root on [derivRoots[numDerivRoots-1],tmax].
            result = RootsPolynomial.findBounded(degree, c,
                derivRoots[numDerivRoots - 1], tmax, maxIterations);
            if (result.found) {
                roots.push(result.root);
                ++numRoots;
            }
        } else {
            // The polynomial is monotone on [tmin,tmax], so has at most
            // one root.
            const result = RootsPolynomial.findBounded(degree, c, tmin, tmax,
                maxIterations);
            if (result.found) {
                roots.push(result.root);
                ++numRoots;
            }
        }
        return numRoots;
    }

    private static evaluate(degree: number, c: readonly number[], t: number): number {
        let i = degree;
        let result = c[i];
        while (--i >= 0) {
            result = t * result + c[i];
        }
        return result;
    }
}
