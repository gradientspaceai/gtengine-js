// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsBisection.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute a root of a function F(t) on an interval [t0, t1]. The caller
// specifies the maximum number of iterations, in case you want limited
// accuracy for the root. However, the function is designed for native
// floating-point types. If you specify a sufficiently large number of
// iterations, the root finder bisects until either F(t) is identically zero
// [a condition dependent on how you structure F(t) for evaluation] or the
// midpoint (t0 + t1)/2 rounds numerically to t0 or t1. Of course, it is
// required that t0 < t1. The 'iterations' member of the result is:
//   0: F(t0)*F(t1) > 0, we cannot determine a root
//   1: F(t0) = 0 or F(t1) = 0
//   2..maxIterations:  the number of bisections plus one
//   maxIterations+1:  the loop executed without a break (no convergence)
//
// Port notes: upstream returns the iteration count and writes the root to a
// reference parameter; the port returns { iterations, root }. The two C++
// overloads of Find become TypeScript overloads distinguished by arity.

export interface RootsBisectionResult {
    // The status/iteration count described above.
    iterations: number;

    // The root estimate. When iterations is 0 the root is not valid (it is
    // upstream's initialization value t0).
    root: number;
}

export class RootsBisection {
    // Use this overload when F(t0) and F(t1) are not already known.
    static find(F: (t: number) => number, t0: number, t1: number,
        maxIterations: number): RootsBisectionResult;

    // If f0 = F(t0) and f1 = F(t1) are already known, pass them to the
    // bisector. This is useful when |f0| or |f1| is infinite, and you can
    // pass sign(f0) or sign(f1) rather than the infinity because the
    // bisector cares only about the signs of f.
    static find(F: (t: number) => number, t0: number, t1: number,
        f0: number, f1: number, maxIterations: number): RootsBisectionResult;

    static find(F: (t: number) => number, t0: number, t1: number,
        arg3: number, arg4?: number, arg5?: number): RootsBisectionResult {
        if (!(t0 < t1)) {
            // The interval endpoints are invalid.
            return { iterations: 0, root: t0 };
        }

        let f0: number, f1: number, maxIterations: number;
        if (arg4 === undefined || arg5 === undefined) {
            maxIterations = arg3;

            // Test the endpoints to see whether F(t) is zero.
            f0 = F(t0);
            if (f0 === 0) {
                return { iterations: 1, root: t0 };
            }

            f1 = F(t1);
            if (f1 === 0) {
                return { iterations: 1, root: t1 };
            }
        } else {
            f0 = arg3;
            f1 = arg4;
            maxIterations = arg5;

            // Test the endpoints to see whether F(t) is zero.
            if (f0 === 0) {
                return { iterations: 1, root: t0 };
            }

            if (f1 === 0) {
                return { iterations: 1, root: t1 };
            }
        }

        if (f0 * f1 > 0) {
            // It is not known whether the interval bounds a root.
            return { iterations: 0, root: t0 };
        }

        let root = 0.5 * (t0 + t1);

        let i: number;
        for (i = 2; i <= maxIterations; ++i) {
            root = 0.5 * (t0 + t1);
            if (root === t0 || root === t1) {
                // The numbers t0 and t1 are consecutive floating-point
                // numbers.
                break;
            }

            const fm = F(root);
            const product = fm * f0;
            if (product < 0) {
                t1 = root;
                f1 = fm;
            } else if (product > 0) {
                t0 = root;
                f0 = fm;
            } else {
                break;
            }
        }
        return { iterations: i <= maxIterations ? i : maxIterations + 1, root };
    }
}
