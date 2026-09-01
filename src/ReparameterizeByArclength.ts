// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ReparameterizeByArclength.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The reparameterization by arclength of a curve can be used for moving along
// a curve at constant speed. The documentation for the algorithms is
// https://www.geometrictools.com/Documentation/MovingAlongCurveSpecifiedSpeed.pdf
//
// Port notes: upstream takes a 'std::shared_ptr<ParametricCurve<N, T>>'; the
// port takes the curve object directly (as ApprCurveByArcs does). The nested
// 'Output' struct becomes the exported interface
// ReparameterizeByArclengthOutput. The 'std::set<T>' of Newton iterates
// becomes a 'Set<number>'; both reject a duplicate insert, which is how a
// cycle is detected.
//
// Numerical accuracy: the arclength F(t,s) = Arclength(tMin,t) - s is
// evaluated by the curve's Romberg integration, so the residual f in the
// output is limited by the accuracy of that quadrature rather than by the
// root finder. As upstream notes, the iteration counts for 'double' are no
// larger than about 53 for bisection.

import { logAssert } from './Logger';
import { ParametricCurve } from './ParametricCurve';

// The output object stores the curve's t-parameter corresponding to a
// user-specified arclength s. The 'f' member is F(t, s). The
// 'numIterations' member is the number of iterations used to compute t.
export interface ReparameterizeByArclengthOutput {
    t: number;
    f: number;
    numIterations: number;
}

// Choose maxIterations sufficiently large for convergence. The value 4096 is
// sufficient. In practice, the number of iterations for type 'float' is no
// larger than approximately 24 and for type 'double' is no larger than
// approximately 53.
const maxIterations = 4096;

export class ReparameterizeByArclength {
    private mCurve: ParametricCurve;
    private mTMin: number;
    private mTMax: number;
    private mTotalArclength: number;

    // The ParametricCurve interface already contains the support for
    // computing arclength s from t. It has support for the inversion (compute
    // t from s), but the code here that uses the hybrid Newton's method and
    // bisection will eventually replace it.
    constructor(curve: ParametricCurve) {
        logAssert(curve !== undefined && curve !== null,
            'The input curve must exist.');
        this.mCurve = curve;
        this.mTMin = curve.getTMin();
        this.mTMax = curve.getTMax();
        this.mTotalArclength = curve.getTotalLength();
    }

    getCurve(): ParametricCurve {
        return this.mCurve;
    }

    getTMin(): number {
        return this.mTMin;
    }

    getTMax(): number {
        return this.mTMax;
    }

    getTotalArclength(): number {
        return this.mTotalArclength;
    }

    // Given an arclength s in [0,L] where the total arclength of the curve is
    // L = Arclength(tMin,tMax), the function returns the root t for
    // F(t,s) = Arclength(tMin,t) - s. Set 'useBisection' to true to use
    // bisection only. Set it to false to use the hybrid of Newton's method
    // and bisection.
    getT(s: number, useBisection: boolean): ReparameterizeByArclengthOutput {
        // Clamp the input to the valid interval.
        if (s <= 0) {
            return { t: this.mTMin, f: 0, numIterations: 0 };
        }

        if (s >= this.mTotalArclength) {
            return { t: this.mTMax, f: 0, numIterations: 0 };
        }

        // Compute a t-root of F(t, s) for the specified s-value. We know that
        // F(mTMin) < 0 and F(mTMax) > 0. Rather than use the initial interval
        // [mTMin,mTMax], choose a subinterval using an initial guess for the
        // t-root.
        let tMin = this.mTMin;
        let tMax = this.mTMax;
        const tMid = tMin + (tMax - tMin) * (s / this.mTotalArclength);
        const fMid = this.F(tMid, s);
        if (fMid > 0) {
            tMax = tMid;
        }
        else {
            tMin = tMid;
        }

        if (useBisection) {
            return this.doBisection(tMin, tMax, s);
        }
        else {
            return this.doNewtonsMethod(tMin, tMax, tMid, s);
        }
    }

    private F(t: number, s: number): number {
        return this.mCurve.getLength(this.mTMin, t) - s;
    }

    private DFDT(t: number): number {
        return this.mCurve.getSpeed(t);
    }

    // Upstream writes tMid and fMid through reference parameters; the port
    // returns them alongside the convergence flag.
    private bisectionConverged(tMin: number, tMax: number, s: number,
        tMid: number, fMid: number):
        { converged: boolean, tMid: number, fMid: number } {
        if (tMid === tMin || tMid === tMax) {
            // The precision of the number type is such that tMin and tMax are
            // consecutive floating-point numbers. Their average cannot be a
            // floating-point number strictly between them. This is the best
            // you can do. Return the t-endpoint whose f-value has smaller
            // magnitude.
            const fMinValue = this.F(tMin, s);
            const fMaxValue = this.F(tMax, s);
            if (Math.abs(fMinValue) <= Math.abs(fMaxValue)) {
                return { converged: true, tMid: tMin, fMid: fMinValue };
            }
            return { converged: true, tMid: tMax, fMid: fMaxValue };
        }
        return { converged: false, tMid, fMid };
    }

    private doBisection(tMin: number, tMax: number, s: number):
        ReparameterizeByArclengthOutput {
        let tMid = 0;
        let fMid = 0;
        let numIterations = 1;
        for (; numIterations <= maxIterations; ++numIterations) {
            // Compute the t-midpoint and the corresponding f-value. Exit
            // early if the f-value is zero.
            tMid = 0.5 * (tMin + tMax);
            fMid = this.F(tMid, s);
            if (fMid === 0) {
                break;
            }

            // Convergence occurs when tMid is tMin or tMax.
            const converged = this.bisectionConverged(tMin, tMax, s, tMid, fMid);
            tMid = converged.tMid;
            fMid = converged.fMid;
            if (converged.converged) {
                break;
            }

            // Update the correct t-endpoint using the t-midpoint.
            if (fMid > 0) {
                tMax = tMid;
            }
            else {
                tMin = tMid;
            }
        }

        return { t: tMid, f: fMid, numIterations };
    }

    private doNewtonsMethod(tMin: number, tMax: number, tMid: number,
        s: number): ReparameterizeByArclengthOutput {
        // Store the iterates from Newton's method in order to determine
        // whether a cycle has occurred. If it does, further iterates will
        // already be in the set, so the function returns when a cycle is
        // detected.
        const tIterates = new Set<number>();

        let fMid = 0;
        let numIterations = 1;
        for (; numIterations <= maxIterations; ++numIterations) {
            // Test whether tMid is an iterate visited previously. If so, a
            // cycle has occurred.
            if (tIterates.has(tMid)) {
                break;
            }
            tIterates.add(tMid);

            // Evaluate F(tMid). Exit early if it is zero.
            fMid = this.F(tMid, s);
            if (fMid === 0) {
                break;
            }

            // Update the bisection interval knowing the sign of F(tMid). The
            // current tMid becomes an endpoint of this interval.
            if (fMid > 0) {
                tMax = tMid;
            }
            else {
                tMin = tMid;
            }

            // Evaluate F'(tMid) >= 0. A bisection step must be taken when
            // F'(tMid) = 0 to avoid the division by zero.
            const dfdt = this.DFDT(tMid);
            if (dfdt === 0) {
                // Division by zero is not allowed. Try the bisection step.
                tMid = 0.5 * (tMin + tMax);
                const converged = this.bisectionConverged(tMin, tMax, s, tMid, fMid);
                tMid = converged.tMid;
                fMid = converged.fMid;
                if (converged.converged) {
                    break;
                }
            }

            // As upstream, this division is reached even when dfdt is zero
            // (fMid is nonzero at this point, so the quotient is an infinity
            // rather than a NaN). The subsequent out-of-interval test then
            // forces the bisection step that was just computed.
            const tNext = tMid - fMid / dfdt;
            if (tNext === tMid) {
                // The precision of the number type is not large enough to
                // disambiguate tMid and tNext. This is the best you can do.
                break;
            }

            // Determine whether to accept the Newton step or take the
            // bisection step.
            tMid = tNext;
            if (tMid < tMin || tMid > tMax) {
                // The iterate is outside the root-bounding interval. Try the
                // bisection step.
                tMid = 0.5 * (tMin + tMax);
                const converged = this.bisectionConverged(tMin, tMax, s, tMid, fMid);
                tMid = converged.tMid;
                fMid = converged.fMid;
                if (converged.converged) {
                    break;
                }
            }
        }

        return { t: tMid, f: fMid, numIterations };
    }
}
