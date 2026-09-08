// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Minimize1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Search for a minimum of F(t) on [t0,t1] using successive parabolic
// interpolation. The search is recursive based on the polyline associated
// with (t,F(t)) at the endpoints and the midpoint of an interval. Let
// f0 = F(t0), f1 = F(t1), tm is in (t0,t1) and fm = F(tm). The polyline is
// {(t0,f0),(tm,fm),(t1,f1)}.
//
// If the polyline is V-shaped, the interval [t0,t1] contains a minimum
// point. The polyline is fit with a parabola whose vertex tv is in (t0,t1).
// Let fv = F(tv). If {(t0,f0),(tv,fv),(tm,fm)} is a minimum bracket, the
// parabolic interpolation continues in [t0,tm]. If instead
// {(tm,fm),(tv,fv),(t1,f1)} is a minimum bracket, the parabolic
// interpolation continues in [tm,t1].
//
// If the polyline is not V-shaped, both subintervals [t0,tm] and [tm,t1]
// are searched for a minimum.
//
// Port notes.
//
// Upstream bug (fixed here, two parts). The bracket update in
// GetBracketedMinimum is unusable whenever F(t0) and F(t1) agree to within
// round-off, which is the normal situation for a function sampled over an
// interval whose two endpoints describe the same configuration. In exact
// arithmetic the vertex of the interpolating parabola is then the midpoint
// of the bracket, and upstream has a branch for exactly that case: it
// examines a neighborhood of the vertex, which keeps both halves of the
// bracket in play. In floating-point arithmetic the vertex misses the
// midpoint by a few ulps, so:
//
//  * The exact 'tv == tm' test fails, one of the asymmetric branches runs,
//    and it discards everything on one side of a point that is a few ulps
//    from the middle sample - a bracket collapse driven entirely by
//    round-off, which throws away the half that may contain the minimum.
//    The port compares tv with tm to within the floating-point resolution
//    of the bracket (8 ulps of its largest endpoint) instead.
//  * On the following iteration the collapsed bracket makes the parabola
//    degenerate, |denom| <= epsilon, and upstream returns, abandoning the
//    search after four function evaluations. The port instead continues
//    with the same neighborhood examination that the vertex-at-the-midpoint
//    case uses, and only the caller's maxBisections bounds the work.
//
// This is upstream issue #298, found through DistOrientedBox3Cone3, which
// sweeps a quad angle over [-pi/2, pi/2] whose endpoint slices are the same
// point set: for the configuration pinned in that file's tests the search
// reported a distance of 0.3628 instead of the true 0.0775 near angle 0.22.
// Both parts are needed; either alone leaves a range of near-symmetric
// inputs broken. Neither changes the iteration when the parabola vertex is
// resolvable, so a well-conditioned search is bit-for-bit as upstream.
//
// The output reference parameters and the two GetMinimum overloads are
// mapped as described at their declarations below.

import { logAssert } from './Logger.js';

// Upstream returns the location and value of the minimum through the output
// reference parameters tMin and fMin; per PORTING.md these become the named
// fields of a returned object.
export interface Minimize1Result {
    // The location of the estimated minimum.
    tMin: number;

    // The value of the function at the estimated minimum, fMin = F(tMin).
    fMin: number;
}

// The port of 'std::function<T(T)>'.
export type Minimize1Function = (t: number) => number;

// The three samples {(t0,f0),(tm,fm),(t1,f1)} that bracket a minimum. Used
// only to return the six updated values of the bracket refinement, which
// upstream performs by assigning to local variables in place.
interface Minimize1Bracket {
    t0: number; f0: number; tm: number; fm: number; t1: number; f1: number;
}

export class Minimize1 {
    private mFunction: Minimize1Function;
    private mMaxSubdivisions: number;
    private mMaxBisections: number;
    private mTMin: number;
    private mFMin: number;
    private mEpsilon: number;
    private mTolerance: number;

    // Construction.
    constructor(F: Minimize1Function, maxSubdivisions: number, maxBisections: number,
        epsilon: number = 1e-08, tolerance: number = 1e-04) {
        this.mFunction = F;
        this.mMaxSubdivisions = maxSubdivisions;
        this.mMaxBisections = maxBisections;
        this.mTMin = 0;
        this.mFMin = 0;
        this.mEpsilon = Math.max(epsilon, 0);
        this.mTolerance = Math.max(tolerance, 0);

        logAssert(
            this.mMaxSubdivisions > 0 && this.mMaxBisections > 0,
            'Invalid argument.');

        this.setEpsilon(epsilon);
        this.setTolerance(tolerance);
    }

    // Member access.
    setEpsilon(epsilon: number): void {
        this.mEpsilon = Math.max(epsilon, 0);
    }

    setTolerance(tolerance: number): void {
        this.mTolerance = Math.max(tolerance, 0);
    }

    getEpsilon(): number {
        return this.mEpsilon;
    }

    getTolerance(): number {
        return this.mTolerance;
    }

    // Search for a minimum of F(t) on the interval [t0,t1] using an initial
    // guess of tInitial. The location of the minimum is tMin and the value
    // of the minimum is fMin = F(tMin). Upstream has two GetMinimum overloads
    // where the 3-input one uses the initial guess (t0+t1)/2; because the
    // output parameters become a returned object, the port merges them into
    // a single function with a defaulted tInitial.
    getMinimum(t0: number, t1: number, tInitial: number = 0.5 * (t0 + t1)): Minimize1Result {
        logAssert(
            t0 <= tInitial && tInitial <= t1,
            'Invalid initial t value.');

        // Compute the minimum for the 3 initial points.
        this.mTMin = Number.MAX_VALUE;
        this.mFMin = Number.MAX_VALUE;

        const f0 = this.mFunction(t0);
        if (f0 < this.mFMin) {
            this.mTMin = t0;
            this.mFMin = f0;
        }

        const fInitial = this.mFunction(tInitial);
        if (fInitial < this.mFMin) {
            this.mTMin = tInitial;
            this.mFMin = fInitial;
        }

        const f1 = this.mFunction(t1);
        if (f1 < this.mFMin) {
            this.mTMin = t1;
            this.mFMin = f1;
        }

        // Search for the global minimum on [t0,t1] with tInitial chosen
        // hopefully to start with a minimum bracket.
        if (((fInitial < f0) && (f1 >= fInitial)) ||
            ((f1 > fInitial) && (f0 >= fInitial))) {
            // The polyline {(t0,f0), (tInitial,fInitial), (t1,f1)} is
            // V-shaped.
            this.getBracketedMinimum(t0, f0, tInitial, fInitial, t1, f1);
        }
        else {
            // The polyline {(t0,f0), (tInitial,fInitial), (t1,f1)} is not
            // V-shaped, so continue searching in subintervals
            // [t0,tInitial] and [tInitial,t1].
            this.subdivide(t0, f0, tInitial, fInitial, this.mMaxSubdivisions);
            this.subdivide(tInitial, fInitial, t1, f1, this.mMaxSubdivisions);
        }

        return { tMin: this.mTMin, fMin: this.mFMin };
    }

    // Search [t0,t1] recursively for a global minimum.
    private subdivide(t0: number, f0: number, t1: number, f1: number,
        subdivisionsRemaining: number): void {
        if (subdivisionsRemaining-- === 0) {
            // The maximum number of subdivisions has been reached.
            return;
        }

        // Compute the function at the midpoint of [t0,t1].
        const tm = 0.5 * (t0 + t1);
        const fm = this.mFunction(tm);
        if (fm < this.mFMin) {
            this.mTMin = tm;
            this.mFMin = fm;
        }

        if (((fm < f0) && (f1 >= fm)) || ((f1 > fm) && (f0 >= fm))) {
            // The polyline {(t0,f0), (tm,fm), (t1,f1)} is V-shaped.
            this.getBracketedMinimum(t0, f0, tm, fm, t1, f1);
        }
        else {
            // The polyline {(t0,f0), (tm,fm), (t1,f1)} is not V-shaped, so
            // continue searching in subintervals [t0,tm] and [tm,t1].
            this.subdivide(t0, f0, tm, fm, subdivisionsRemaining);
            this.subdivide(tm, fm, t1, f1, subdivisionsRemaining);
        }
    }

    // This is called when {f0,fm,f1} brackets a minimum.
    private getBracketedMinimum(t0: number, f0: number, tm: number, fm: number,
        t1: number, f1: number): void {
        const two = 2;
        const half = 0.5;

        for (let i = 0; i < this.mMaxBisections; ++i) {
            // Update the minimum location and value.
            if (fm < this.mFMin) {
                this.mTMin = tm;
                this.mFMin = fm;
            }

            // Test for convergence.
            const dt10 = t1 - t0;
            const dtBound = two * this.mTolerance * Math.abs(tm) + this.mEpsilon;
            if (dt10 <= dtBound) {
                break;
            }

            // Compute the vertex of the interpolating parabola.
            const dt0m = t0 - tm;
            const dt1m = t1 - tm;
            const df0m = f0 - fm;
            const df1m = f1 - fm;
            const tmp0 = dt0m * df1m;
            const tmp1 = dt1m * df0m;
            const denom = tmp1 - tmp0;
            if (Math.abs(denom) <= this.mEpsilon) {
                // The parabola through the three samples is degenerate, so
                // its vertex carries no information. Upstream returns here,
                // abandoning the search; see the port note at the top of the
                // file. The port instead continues with the same
                // neighborhood examination it uses when the vertex lands on
                // the middle sample.
                ({ t0, f0, tm, fm, t1, f1 } =
                    this.refineAroundMidpoint(t0, f0, tm, fm, t1, f1));
                continue;
            }

            // Compute tv and clamp to [t0,t1] to offset floating-point
            // rounding errors.
            let tv = tm + half * (dt1m * tmp1 - dt0m * tmp0) / denom;
            tv = Math.max(t0, Math.min(tv, t1));
            const fv = this.mFunction(tv);
            if (fv < this.mFMin) {
                this.mTMin = tv;
                this.mFMin = fv;
            }

            // The vertex is distinguishable from the middle sample only when
            // it is farther from it than the floating-point resolution of the
            // bracket. Upstream compares tv with tm exactly; see the port
            // note at the top of the file.
            const tvBound = 8 * Number.EPSILON
                * Math.max(Math.abs(t0), Math.abs(t1));

            if (tv < tm - tvBound) {
                if (fv < fm) {
                    t1 = tm;
                    f1 = fm;
                    tm = tv;
                    fm = fv;
                }
                else {
                    t0 = tv;
                    f0 = fv;
                }
            }
            else if (tv > tm + tvBound) {
                if (fv < fm) {
                    t0 = tm;
                    f0 = fm;
                    tm = tv;
                    fm = fv;
                }
                else {
                    t1 = tv;
                    f1 = fv;
                }
            }
            else {
                // The vertex of the parabola is located at the middle sample
                // point. A minimum could occur on either subinterval, but it
                // is also possible the minimum occurs at the vertex. In
                // either case, the search is continued by examining a
                // neighborhood of the vertex.
                ({ t0, f0, tm, fm, t1, f1 } =
                    this.refineAroundMidpoint(t0, f0, tm, fm, t1, f1));
            }
        }
    }

    // Examine a neighborhood of the middle sample and return the
    // subinterval that brackets a minimum. When two choices exist for a
    // bracket, the one with the smallest function value at the midpoint is
    // used. This is upstream's 'vertex of the parabola is at the middle
    // sample point' case, extracted so that the degenerate-parabola case can
    // use it as well.
    private refineAroundMidpoint(t0: number, f0: number, tm: number,
        fm: number, t1: number, f1: number): Minimize1Bracket {
        const half = 0.5;
        const tm0 = half * (t0 + tm);
        const fm0 = this.mFunction(tm0);
        const tm1 = half * (tm + t1);
        const fm1 = this.mFunction(tm1);

        if (fm0 < fm) {
            if (fm1 < fm) {
                if (fm0 < fm1) {
                    // {(t0,f0),(tm0,fm0),(tm,fm)}
                    return { t0, f0, tm: tm0, fm: fm0, t1: tm, f1: fm };
                }
                // {(tm,fm),(tm1,fm1),(t1,f1)}
                return { t0: tm, f0: fm, tm: tm1, fm: fm1, t1, f1 };
            }
            // fm1 >= fm: {(t0,f0),(tm0,fm0),(tm,fm)}
            return { t0, f0, tm: tm0, fm: fm0, t1: tm, f1: fm };
        }
        else if (fm0 > fm) {
            if (fm1 < fm) {
                // {(tm,fm),(tm1,fm1),(t1,f1)}
                return { t0: tm, f0: fm, tm: tm1, fm: fm1, t1, f1 };
            }
            // fm1 >= fm: {(tm0,fm0),(tm,fm),(tm1,fm1)}
            return { t0: tm0, f0: fm0, tm, fm, t1: tm1, f1: fm1 };
        }
        else {
            // fm0 = fm
            if (fm1 < fm) {
                // {(tm,fm),(tm1,fm1),(t1,f1)}
                return { t0: tm, f0: fm, tm: tm1, fm: fm1, t1, f1 };
            }
            // fm1 >= fm: {(tm0,fm0),(tm,fm),(tm1,fm1)}
            return { t0: tm0, f0: fm0, tm, fm, t1: tm1, f1: fm1 };
        }
    }
}
