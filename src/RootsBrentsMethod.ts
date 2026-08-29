// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RootsBrentsMethod.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This is an implementation of Brent's Method for computing a root of a
// function on an interval [t0,t1] for which F(t0)*F(t1) < 0. The method
// uses inverse quadratic interpolation to generate a root estimate but
// falls back to inverse linear interpolation (secant method) if necessary.
// Moreover, based on previous iterates, the method will fall back to
// bisection when it appears the interpolated estimate is not of sufficient
// quality.
//
//   maxIterations:
//       The maximum number of iterations used to locate a root. This
//       should be positive.
//   negFTolerance, posFTolerance:
//       The root estimate t is accepted when the function value F(t)
//       satisfies negFTolerance <= F(t) <= posFTolerance. The values must
//       satisfy: negFTolerance <= 0, posFTolerance >= 0.
//   stepTTolerance:
//       Brent's Method requires additional tests before an interpolated
//       t-value is accepted as the next root estimate. One of these tests
//       compares the difference of consecutive iterates and requires it to
//       be larger than a user-specified t-tolerance (to ensure progress is
//       made). This parameter is that tolerance and should be nonnegative.
//   convTTolerance:
//       The root search is allowed to terminate when the current
//       subinterval [tsub0,tsub1] is sufficiently small, say,
//       |tsub1 - tsub0| <= tolerance. This parameter is that tolerance and
//       should be nonnegative.
//
// Port notes: upstream returns bool and writes the root to a reference
// parameter; the port returns { found, root } where root is valid only when
// found is true (root is 0 otherwise).

export interface RootsBrentsMethodResult {
    found: boolean;

    // The root estimate; valid only when found is true.
    root: number;
}

export class RootsBrentsMethod {
    // It is necessary that F(t0)*F(t1) <= 0, in which case the function
    // returns found = true and the 'root' is valid; otherwise, the function
    // returns found = false and 'root' is invalid (do not use it). When
    // F(t0)*F(t1) > 0, the interval may very well contain a root but we
    // cannot know that. The function also returns found = false if t0 >= t1.
    static find(F: (t: number) => number, t0: number, t1: number,
        maxIterations: number, negFTolerance: number, posFTolerance: number,
        stepTTolerance: number, convTTolerance: number): RootsBrentsMethodResult {
        // Parameter validation.
        if (t1 <= t0
            || maxIterations === 0
            || negFTolerance > 0
            || posFTolerance < 0
            || stepTTolerance < 0
            || convTTolerance < 0) {
            // The input is invalid.
            return { found: false, root: 0 };
        }

        let f0 = F(t0);
        if (negFTolerance <= f0 && f0 <= posFTolerance) {
            // This endpoint is an approximate root that satisfies the
            // function tolerance.
            return { found: true, root: t0 };
        }

        let f1 = F(t1);
        if (negFTolerance <= f1 && f1 <= posFTolerance) {
            // This endpoint is an approximate root that satisfies the
            // function tolerance.
            return { found: true, root: t1 };
        }

        if (f0 * f1 > 0) {
            // The input interval must bound a root.
            return { found: false, root: 0 };
        }

        if (Math.abs(f0) < Math.abs(f1)) {
            // Swap t0 and t1 so that |F(t1)| <= |F(t0)|. The number t1 is
            // considered to be the best estimate of the root.
            let temp = t0; t0 = t1; t1 = temp;
            temp = f0; f0 = f1; f1 = temp;
        }

        // Initialize values for the root search.
        let t2 = t0, t3 = t0, f2 = f0;
        let prevBisected = true;

        // The root search.
        for (let i = 0; i < maxIterations; ++i) {
            const fDiff01 = f0 - f1, fDiff02 = f0 - f2, fDiff12 = f1 - f2;
            const invFDiff01 = 1 / fDiff01;
            let s: number;
            if (fDiff02 !== 0 && fDiff12 !== 0) {
                // Use inverse quadratic interpolation.
                const invFDiff02 = 1 / fDiff02;
                const invFDiff12 = 1 / fDiff12;
                s =
                    t0 * f1 * f2 * invFDiff01 * invFDiff02 -
                    t1 * f0 * f2 * invFDiff01 * invFDiff12 +
                    t2 * f0 * f1 * invFDiff02 * invFDiff12;
            } else {
                // Use inverse linear interpolation (secant method).
                s = (t1 * f0 - t0 * f1) * invFDiff01;
            }

            // Compute values needed in the accept-or-reject tests.
            const tDiffSAvr = s - 0.75 * t0 - 0.25 * t1;
            const tDiffS1 = s - t1;
            const absTDiffS1 = Math.abs(tDiffS1);
            const absTDiff12 = Math.abs(t1 - t2);
            const absTDiff23 = Math.abs(t2 - t3);

            let currBisected = false;
            if (tDiffSAvr * tDiffS1 > 0) {
                // The value s is not between 0.75*t0 + 0.25*t1 and t1.
                // NOTE: The algorithm sometimes has t0 < t1 but sometimes
                // t1 < t0, so the between-ness test does not use simple
                // comparisons.
                currBisected = true;
            } else if (prevBisected) {
                // The first of Brent's tests to determine whether to accept
                // the interpolated s-value.
                currBisected =
                    (absTDiffS1 >= 0.5 * absTDiff12) ||
                    (absTDiff12 <= stepTTolerance);
            } else {
                // The second of Brent's tests to determine whether to
                // accept the interpolated s-value.
                currBisected =
                    (absTDiffS1 >= 0.5 * absTDiff23) ||
                    (absTDiff23 <= stepTTolerance);
            }

            if (currBisected) {
                // One of the additional tests failed, so reject the
                // interpolated s-value and use bisection instead.
                s = 0.5 * (t0 + t1);
                if (s === t0 || s === t1) {
                    // The numbers t0 and t1 are consecutive floating-point
                    // numbers.
                    return { found: true, root: s };
                }
                prevBisected = true;
            } else {
                prevBisected = false;
            }

            // Evaluate the function at the new estimate and test for
            // convergence.
            const fs = F(s);
            if (negFTolerance <= fs && fs <= posFTolerance) {
                return { found: true, root: s };
            }

            // Update the subinterval to include the new estimate as an
            // endpoint.
            t3 = t2;
            t2 = t1;
            f2 = f1;
            if (f0 * fs < 0) {
                t1 = s;
                f1 = fs;
            } else {
                t0 = s;
                f0 = fs;
            }

            // Allow the algorithm to terminate when the subinterval is
            // sufficiently small.
            if (Math.abs(t1 - t0) <= convTTolerance) {
                return { found: true, root: t1 };
            }

            // A loop invariant is that t1 is the root estimate,
            // F(t0)*F(t1) < 0 and |F(t1)| <= |F(t0)|.
            if (Math.abs(f0) < Math.abs(f1)) {
                let temp = t0; t0 = t1; t1 = temp;
                temp = f0; f0 = f1; f1 = temp;
            }
        }

        // Failed to converge in the specified number of iterations.
        return { found: false, root: 0 };
    }
}
