// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsBisection1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Estimate a root on an interval [tMin,tMax] for a continuous function F(t)
// defined on that interval. If a root is found, the function returns it via
// 'root'. Additionally, fAtRoot = F(root) is returned in case the caller
// wants to know how close to zero the function is at the root; numerical
// rounding errors can cause fAtRoot not to be exactly zero. The returned
// 'iterations' is the number of iterations used by the bisector. If that
// number is 0, F(tMin)*F(tMax) > 0 and it is unknown whether [tMin,tMax]
// contains a root. If that number is 1, either F(tMin) = 0 or F(tMax) = 0
// (exactly), and 'root' is the corresponding interval endpoint. If that
// number is 2 or larger, the bisection is applied until 'root' is found for
// which F(root) is exactly 0 or until the current root estimate is equal to
// tMin or tMax. The latter conditions can occur because of the fixed
// precision used in the computations (53-bit precision for the IEEE
// double-precision numbers used by this port).
//
// Port notes:
// - Upstream has two constructors selected by type traits, one for
//   floating-point Real and one for arbitrary-precision Real. Only the
//   floating-point instantiation is ported, so the constructor takes just
//   maxIterations and the private RoundInitial/RoundAverage helpers (which
//   are identity and the plain average for floating-point Real) are inlined.
// - The two operator() overloads become the method find, distinguished by
//   arity; upstream writes tRoot/fAtTRoot to reference parameters and
//   returns the iteration count, so the port returns
//   { iterations, root, fAtRoot }.

import { logAssert } from './Logger.js';

export interface RootsBisection1Result {
    // The status/iteration count described above.
    iterations: number;

    // The root estimate. When 'iterations' is 0 the root is not valid
    // (it is 0, upstream's zero-valued output for that case).
    root: number;

    // F(root). When 'iterations' is 0 or 1 this is 0, as upstream.
    fAtRoot: number;
}

export class RootsBisection1 {
    private mMaxIterations: number;

    constructor(maxIterations: number) {
        logAssert(maxIterations > 0, 'Invalid maximum iterations.');
        this.mMaxIterations = maxIterations;
    }

    // Use this overload when F(tMin) and F(tMax) are not already known.
    find(F: (t: number) => number, tMin: number, tMax: number): RootsBisection1Result;

    // Use this overload when fMin = F(tMin) and fMax = F(tMax) are already
    // known. This is useful when |fMin| or |fMax| is infinite, whereby you
    // can pass sign(fMin) or sign(fMax) rather than an infinity because the
    // bisector cares only about the signs of F(t).
    find(F: (t: number) => number, tMin: number, tMax: number,
        fMin: number, fMax: number): RootsBisection1Result;

    find(F: (t: number) => number, tMin: number, tMax: number,
        fMin?: number, fMax?: number): RootsBisection1Result {
        logAssert(tMin < tMax, 'Invalid ordering of t-interval endpoints.');

        if (fMin === undefined || fMax === undefined) {
            // Floating-point inputs are used as is (upstream rounds only
            // arbitrary-precision inputs to the specified precision).
            fMin = F(tMin);
            fMax = F(tMax);
        }

        const sign0 = (fMin > 0 ? +1 : (fMin < 0 ? -1 : 0));
        if (sign0 === 0) {
            return { iterations: 1, root: tMin, fAtRoot: 0 };
        }

        const sign1 = (fMax > 0 ? +1 : (fMax < 0 ? -1 : 0));
        if (sign1 === 0) {
            return { iterations: 1, root: tMax, fAtRoot: 0 };
        }

        if (sign0 === sign1) {
            // It is unknown whether the interval contains a root.
            return { iterations: 0, root: 0, fAtRoot: 0 };
        }

        // The bisection steps. Upstream leaves tRoot and fAtTRoot unassigned
        // when mMaxIterations is 1 (the loop body never executes); the port
        // initializes them to zero rather than returning garbage.
        let root = 0, fAtRoot = 0;
        let t0 = tMin, t1 = tMax;
        let iteration: number;
        for (iteration = 2; iteration <= this.mMaxIterations; ++iteration) {
            // Use the floating-point average as is.
            root = 0.5 * (t0 + t1);
            fAtRoot = F(root);

            // If the function is exactly zero, a root is found. For fixed
            // precision, the average of two consecutive numbers might be one
            // of the current interval endpoints.
            const signRoot = (fAtRoot > 0 ? +1 : (fAtRoot < 0 ? -1 : 0));
            if (signRoot === 0 || root === t0 || root === t1) {
                break;
            }

            // Update the correct endpoint to the midpoint.
            if (signRoot === sign0) {
                t0 = root;
            } else { // signRoot === sign1
                t1 = root;
            }
        }
        return { iterations: iteration, root, fAtRoot };
    }
}
