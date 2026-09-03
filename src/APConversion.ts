// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) APConversion.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The conversion functions here are used to obtain arbitrary-precision
// approximations to rational numbers and to quadratic field numbers. The
// arbitrary-precision arithmetic is described in
// https://www.geometrictools.com/Documentation/ArbitraryPrecision.pdf
// The quadratic field numbers and conversions are described in
// https://www.geometrictools.com/Documentation/QuadraticFields.pdf
//
// Port notes:
//   * Upstream is the class template APConversion<Rational>. The estimates
//     rely on exact rational arithmetic (division, std::frexp, std::ldexp
//     and the Convert functions with directed rounding), so the only
//     meaningful instantiation is Rational = BSRational. Following the
//     precedent that a template layer with a single instantiation is
//     dropped (B34), this port is the concrete BSRational instantiation.
//   * The functions with Rational output parameters return object literals
//     with named fields, and the two overloads of EstimateSqrt and Estimate
//     (bounding interval versus single estimate) become distinct methods:
//       EstimateSqrt(aSqr, aMin, aMax) -> estimateSqrt(aSqr)
//         -> { numIterates, aMin, aMax }
//       EstimateSqrt(aSqr, a)          -> estimateSqrtValue(aSqr)
//         -> { numIterates, a }
//       EstimateApB(...)               -> estimateApB(aSqr, bSqr)
//         -> { numIterates, tMin, tMax }
//       EstimateAmB(...)               -> estimateAmB(aSqr, bSqr)
//         -> { numIterates, tMin, tMax }
//       Estimate(q, qMin, qMax)        -> estimate(q)
//         -> { numIterates, qMin, qMax }
//       Estimate(q, qEstimate)         -> estimateValue(q)
//         -> { numIterates, qEstimate }
//   * The input of Estimate is upstream's QFN1 = QFNumber<Rational, 1>. The
//     ported QFNumber is the 'number'-coefficient path (see QFNumber.ts), so
//     it cannot carry the BSRational coefficients required here; the input
//     is described instead by the structural interface APConversionQFN1,
//     which has the same members (x[0], x[1] and d) that upstream's Estimate
//     reads. The upstream alias QFN2 is declared but never used and is not
//     ported.
//   * std::nextafter(x, +-max) is not available in JavaScript; the port
//     reuses the IEEE binary64 helpers nextUp and nextDown from SWInterval.ts,
//     which reproduce it exactly.
//   * The C++ Convert(input, precision, mode, output) becomes
//     convertBSRational(input, precision, mode) and Convert(input, mode,
//     double& output) becomes convertBSRationalToNumber(input, mode).

import { BSNumberRoundingMode } from './BSNumber.js';
import { BSRational, convertBSRational, convertBSRationalToNumber } from './BSRational.js';
import { logAssert, logError } from './Logger.js';
import { nextDown, nextUp } from './SWInterval.js';

// The structural form of upstream's QFN1 = QFNumber<Rational, 1> with
// Rational = BSRational: the quadratic field number x[0] + x[1] * sqrt(d).
export interface APConversionQFN1 {
    x: readonly [BSRational, BSRational];
    d: BSRational;
}

export class APConversion {
    // Convenient constants.
    private readonly mZero: BSRational;
    private readonly mOne: BSRational;
    private readonly mThree: BSRational;
    private readonly mFive: BSRational;

    private mPrecision: number;
    private mMaxIterations: number;
    private mThreshold: BSRational;

    constructor(precision: number, maxIterations: number) {
        logAssert(precision > 0, 'Invalid precision.');
        logAssert(maxIterations > 0, 'Invalid maximum iterations.');

        this.mZero = new BSRational();
        this.mOne = BSRational.fromNumber(1);
        this.mThree = BSRational.fromNumber(3);
        this.mFive = BSRational.fromNumber(5);
        this.mPrecision = precision;
        this.mMaxIterations = maxIterations;
        this.mThreshold = BSRational.ldexp(this.mOne, -precision);
    }

    // Member access.
    setPrecision(precision: number): void {
        logAssert(precision > 0, 'Invalid precision.');
        this.mPrecision = precision;
        this.mThreshold = BSRational.ldexp(this.mOne, -precision);
    }

    setMaxIterations(maxIterations: number): void {
        logAssert(maxIterations > 0, 'Invalid maximum iterations.');
        this.mMaxIterations = maxIterations;
    }

    getPrecision(): number {
        return this.mPrecision;
    }

    getMaxIterations(): number {
        return this.mMaxIterations;
    }

    // The input a^2 is rational, but a itself is usually irrational,
    // although a rational value is allowed. Compute a bounding interval for
    // the root, aMin <= a <= aMax, where the endpoints are both within the
    // specified precision.
    estimateSqrt(aSqr: BSRational): { numIterates: number, aMin: BSRational, aMax: BSRational } {
        // Factor a^2 = r^2 * 2^e, where r^2 in [1/2,1). Compute s^2 and the
        // exponent used to generate the estimate of sqrt(a^2).
        const { rSqr: sSqr, exponent: exponentA } = APConversion.preprocessSqr(aSqr);

        // Use the FPU to estimate s = sqrt(sSqr) to 53-bit precision with
        // rounding up. Multiply by the appropriate exponent to obtain the
        // upper bound aMax > a.
        let aMax = APConversion.getMaxOfSqrt(sSqr, exponentA);

        // Compute a lower bound aMin < a.
        let aMin = aSqr.div(aMax);

        // Compute Newton iterates until convergence. The estimate closest to
        // a is aMin with aMin <= a <= aMax and a - aMin <= aMax - a.
        let iterate = 1;
        for (; iterate <= this.mMaxIterations; ++iterate) {
            if (aMax.sub(aMin).lessThan(this.mThreshold)) {
                break;
            }
            // Compute the average aMax = (aMin + aMax) / 2. Round up to twice
            // the precision to avoid quadratic growth in the number of bits
            // and to ensure that aMin can increase.
            aMax = BSRational.ldexp(aMin.add(aMax), -1);
            aMax = convertBSRational(aMax, 2 * this.mPrecision,
                BSNumberRoundingMode.FE_UPWARD);
            aMin = aSqr.div(aMax);
        }
        return { numIterates: iterate, aMin, aMax };
    }

    // Compute an estimate of the root when you do not need a bounding
    // interval.
    estimateSqrtValue(aSqr: BSRational): { numIterates: number, a: BSRational } {
        // Compute a bounding interval aMin <= a <= aMax.
        const { numIterates, aMin, aMax } = this.estimateSqrt(aSqr);
        // Use the average of the interval endpoints as the estimate.
        return { numIterates, a: BSRational.ldexp(aMin.add(aMax), -1) };
    }

    // Compute a bounding interval tMin <= a + b <= tMax for the sum of the
    // square roots a = sqrt(aSqr) and b = sqrt(bSqr).
    estimateApB(aSqr: BSRational, bSqr: BSRational):
        { numIterates: number, tMin: BSRational, tMax: BSRational } {
        // Factor a^2 = r^2 * 2^e, where r^2 in [1/2,1). Compute u^2 and the
        // exponent used to generate the estimate of sqrt(a^2).
        const { rSqr: uSqr, exponent: exponentA } = APConversion.preprocessSqr(aSqr);

        // Factor b^2 = s^2 * 2^e, where s^2 in [1/2,1). Compute v^2 and the
        // exponent used to generate the estimate of sqrt(b^2).
        const { rSqr: vSqr, exponent: exponentB } = APConversion.preprocessSqr(bSqr);

        // Use the FPU to estimate u = sqrt(u^2) and v = sqrt(v^2) to 53 bits
        // of precision with rounding up. Multiply by the appropriate
        // exponents to obtain upper bounds aMax > a and bMax > b. This
        // ensures tMax = aMax + bMax > a + b.
        const aMax = APConversion.getMaxOfSqrt(uSqr, exponentA);
        const bMax = APConversion.getMaxOfSqrt(vSqr, exponentB);
        let tMax = aMax.add(bMax);

        // Compute a lower bound tMin < a + b.
        const a2pb2 = aSqr.add(bSqr);
        const a2mb2 = aSqr.sub(bSqr);
        const a2mb2Sqr = a2mb2.mul(a2mb2);
        let tMaxSqr = tMax.mul(tMax);
        let tMin = a2pb2.mul(tMaxSqr).sub(a2mb2Sqr).div(tMax.mul(tMaxSqr.sub(a2pb2)));

        // Compute Newton iterates until convergence. The estimate closest to
        // a + b is tMin with tMin < a + b < tMax and
        // (a + b) - tMin < tMax - (a + b).
        let iterate = 1;
        for (; iterate <= this.mMaxIterations; ++iterate) {
            if (tMax.sub(tMin).lessThan(this.mThreshold)) {
                break;
            }
            // Compute the weighted average tMax = (3*tMin + tMax) / 4. Round
            // up to twice the precision to avoid quadratic growth in the
            // number of bits and to ensure that tMin can increase.
            tMax = BSRational.ldexp(this.mThree.mul(tMax).add(tMin), -2);
            tMax = convertBSRational(tMax, 2 * this.mPrecision,
                BSNumberRoundingMode.FE_UPWARD);
            tMaxSqr = tMax.mul(tMax);
            tMin = a2pb2.mul(tMaxSqr).sub(a2mb2Sqr).div(tMax.mul(tMaxSqr.sub(a2pb2)));
        }
        return { numIterates: iterate, tMin, tMax };
    }

    // Compute a bounding interval tMin <= a - b <= tMax for the difference
    // of the square roots a = sqrt(aSqr) and b = sqrt(bSqr).
    //
    // Upstream requires (but does not verify) that aSqr >= bSqr, so that
    // a - b >= 0: the f" > 0 branch clamps the initial lower bound to zero
    // to "stay on the nonnegative t-axis", which is invalid when a < b. The
    // port preserves the behavior and documents the precondition rather than
    // adding a guard upstream does not have.
    estimateAmB(aSqr: BSRational, bSqr: BSRational):
        { numIterates: number, tMin: BSRational, tMax: BSRational } {
        // The return value of the function.
        let iterate = 0;

        // Compute various quantities that are used later in the code.
        const a2tb2 = aSqr.mul(bSqr);                       // a^2 * b^2
        const a2pb2 = aSqr.add(bSqr);                       // a^2 + b^2
        const a2mb2 = aSqr.sub(bSqr);                       // a^2 - b^2
        const a2mb2Sqr = a2mb2.mul(a2mb2);                  // (a^2 - b^2)^2
        const twoa2pb2 = BSRational.ldexp(a2pb2, 1);        // 2 * (a^2 + b^2)

        // Factor a^2 = r^2 * 2^e, where r^2 in [1/2,1). Compute u^2 and the
        // exponent used to generate the estimate of sqrt(a^2).
        const { rSqr: uSqr, exponent: exponentA } = APConversion.preprocessSqr(aSqr);

        // Factor b^2 = s^2 * 2^e, where s^2 in [1/2,1). Compute v^2 and the
        // exponent used to generate the estimate of sqrt(b^2).
        const { rSqr: vSqr, exponent: exponentB } = APConversion.preprocessSqr(bSqr);

        // Compute the sign of f''(a-b)/8 = a^2 - 3*a*b + b^2. It can be shown
        // that Sign(a^2-3*a*b+b^2) = Sign(a^4-7*a^2*b^2+b^4) =
        // Sign((a^2-b^2)^2-5*a^2*b^2).
        let signSecDer = a2mb2Sqr.sub(this.mFive.mul(a2tb2));

        // Local variables shared by the two main blocks of code.
        let aMin: BSRational, aMax: BSRational, bMin: BSRational, bMax: BSRational;
        let tMin: BSRational, tMax: BSRational;
        let tMinSqr: BSRational, tMaxSqr: BSRational, tMid: BSRational, tMidSqr: BSRational;
        let f: BSRational;

        if (signSecDer.greaterThan(this.mZero)) {
            // Choose an initial guess tMin < a-b. Use the FPU to estimate
            // u = sqrt(u^2) and v = sqrt(v^2) to 53 bits of precision with
            // specified rounding. Multiply by the appropriate exponents to
            // obtain tMin = aMin - bMax < a-b.
            aMin = APConversion.getMinOfSqrt(uSqr, exponentA);
            bMax = APConversion.getMaxOfSqrt(vSqr, exponentB);
            tMin = aMin.sub(bMax);

            // When a-b is nearly zero, it is possible the lower bound is
            // negative. Clamp tMin to zero to stay on the nonnegative t-axis
            // where the f"-positive basin is.
            if (tMin.lessThan(this.mZero)) {
                tMin = this.mZero.clone();
            }

            // Test whether tMin is in the positive f"(t) basin containing
            // a-b. If it is not, compute a tMin that is in the basin. The
            // sign test is applied to f"(t)/4 = 3*t^2 - (a^2+b^2).
            tMinSqr = tMin.mul(tMin);
            signSecDer = this.mThree.mul(tMinSqr).sub(a2pb2);
            if (signSecDer.lessThan(this.mZero)) {
                // The initial guess satisfies f"(tMin) < 0. Compute an upper
                // bound tMax > a-b and bisect [tMin,tMax] until either the
                // t-value is an estimate to a-b within the specified
                // precision or until f"(t) >= 0 and f(t) >= 0. In the latter
                // case, continue on to Newton's method, which is then
                // guaranteed to converge.
                aMax = APConversion.getMaxOfSqrt(uSqr, exponentA);
                bMin = APConversion.getMinOfSqrt(vSqr, exponentB);
                tMax = aMax.sub(bMin);

                for (iterate = 1; iterate <= this.mMaxIterations; ++iterate) {
                    if (tMax.sub(tMin).lessThan(this.mThreshold)) {
                        return { numIterates: iterate, tMin, tMax };
                    }

                    tMid = BSRational.ldexp(tMin.add(tMax), -1);
                    tMidSqr = tMid.mul(tMid);
                    signSecDer = this.mThree.mul(tMidSqr).sub(a2pb2);
                    if (signSecDer.greaterThanOrEqual(this.mZero)) {
                        f = tMidSqr.mul(tMidSqr.sub(twoa2pb2)).add(a2mb2Sqr);
                        if (f.greaterThanOrEqual(this.mZero)) {
                            tMin = tMid;
                            tMinSqr = tMidSqr;
                            break;
                        } else {
                            // Round up to twice the precision to avoid
                            // quadratic growth in the number of bits.
                            tMax = convertBSRational(tMid, 2 * this.mPrecision,
                                BSNumberRoundingMode.FE_UPWARD);
                        }
                    } else {
                        // Round down to twice the precision to avoid
                        // quadratic growth in the number of bits.
                        tMin = convertBSRational(tMid, 2 * this.mPrecision,
                            BSNumberRoundingMode.FE_DOWNWARD);
                    }
                }

                // PORT FIX (upstream): when the bisection loop above exits by
                // exhausting mMaxIterations rather than by the break, tMin
                // was last updated in the rounding-down branch, which does
                // not recompute tMinSqr. Upstream then uses a tMinSqr that
                // does not correspond to tMin, which corrupts the bound
                // computed next. Recompute it here; this is a no-op on every
                // path where upstream is already consistent (the break sets
                // tMinSqr = tMidSqr = tMin*tMin exactly).
                tMinSqr = tMin.mul(tMin);
            }

            // Compute an upper bound tMax > a-b.
            tMax = a2pb2.mul(tMinSqr).sub(a2mb2Sqr).div(tMin.mul(tMinSqr.sub(a2pb2)));

            // Compute Newton iterates until convergence. The estimate closest
            // to a-b is tMax with tMin < a-b < tMax and
            // tMax - (a-b) < (a-b) - tMin.
            for (iterate = 1; iterate <= this.mMaxIterations; ++iterate) {
                if (tMax.sub(tMin).lessThan(this.mThreshold)) {
                    break;
                }
                // Compute the weighted average tMin = (3*tMin+tMax)/4. Round
                // down to twice the precision to avoid quadratic growth in
                // the number of bits and to ensure that tMax can decrease.
                tMin = BSRational.ldexp(this.mThree.mul(tMin).add(tMax), -2);
                tMin = convertBSRational(tMin, 2 * this.mPrecision,
                    BSNumberRoundingMode.FE_DOWNWARD);
                tMinSqr = tMin.mul(tMin);
                tMax = a2pb2.mul(tMinSqr).sub(a2mb2Sqr).div(tMin.mul(tMinSqr.sub(a2pb2)));
            }
            return { numIterates: iterate, tMin, tMax };
        }

        if (signSecDer.lessThan(this.mZero)) {
            // Choose an initial guess tMax > a-b. Use the FPU to estimate
            // u = sqrt(u^2) and v = sqrt(v^2) to 53 bits of precision with
            // specified rounding. Multiply by the appropriate exponents to
            // obtain tMax = aMax - bMin > a-b.
            aMax = APConversion.getMaxOfSqrt(uSqr, exponentA);
            bMin = APConversion.getMinOfSqrt(vSqr, exponentB);
            tMax = aMax.sub(bMin);

            // Test whether tMax is in the negative f"(t) basin containing
            // a-b. If it is not, compute a tMax that is in the basin. The
            // sign test is applied to f"(t)/4 = 3*t^2 - (a^2+b^2).
            tMaxSqr = tMax.mul(tMax);
            signSecDer = this.mThree.mul(tMaxSqr).sub(a2pb2);
            if (signSecDer.greaterThan(this.mZero)) {
                // The initial guess satisfies f"(tMax) > 0. Compute a lower
                // bound tMin < a-b and bisect [tMin,tMax] until either the
                // t-value is an estimate to a-b within the specified
                // precision or until f"(t) <= 0 and f(t) <= 0. In the latter
                // case, continue on to Newton's method, which is then
                // guaranteed to converge.
                aMin = APConversion.getMinOfSqrt(uSqr, exponentA);
                bMax = APConversion.getMaxOfSqrt(vSqr, exponentB);
                tMin = aMin.sub(bMax);

                for (iterate = 1; iterate <= this.mMaxIterations; ++iterate) {
                    if (tMax.sub(tMin).lessThan(this.mThreshold)) {
                        return { numIterates: iterate, tMin, tMax };
                    }

                    tMid = BSRational.ldexp(tMin.add(tMax), -1);
                    tMidSqr = tMid.mul(tMid);
                    signSecDer = this.mThree.mul(tMidSqr).sub(a2pb2);
                    if (signSecDer.lessThanOrEqual(this.mZero)) {
                        f = tMidSqr.mul(tMidSqr.sub(twoa2pb2)).add(a2mb2Sqr);
                        if (f.lessThanOrEqual(this.mZero)) {
                            tMax = tMid;
                            tMaxSqr = tMidSqr;
                            break;
                        } else {
                            // Round down to twice the precision to avoid
                            // quadratic growth in the number of bits.
                            tMin = convertBSRational(tMid, 2 * this.mPrecision,
                                BSNumberRoundingMode.FE_DOWNWARD);
                        }
                    } else {
                        // Round up to twice the precision to avoid quadratic
                        // growth in the number of bits.
                        tMax = convertBSRational(tMid, 2 * this.mPrecision,
                            BSNumberRoundingMode.FE_UPWARD);
                    }
                }

                // PORT FIX (upstream): see the corresponding comment in the
                // f" > 0 block. On exhaustion of the bisection iterations,
                // upstream's tMaxSqr does not correspond to its tMax.
                tMaxSqr = tMax.mul(tMax);
            }

            // Compute a lower bound tMin < a-b.
            tMin = a2pb2.mul(tMaxSqr).sub(a2mb2Sqr).div(tMax.mul(tMaxSqr.sub(a2pb2)));

            // Compute Newton iterates until convergence. The estimate closest
            // to a-b is tMin with tMin < a - b < tMax and
            // (a-b) - tMin < tMax - (a-b).
            for (iterate = 1; iterate <= this.mMaxIterations; ++iterate) {
                if (tMax.sub(tMin).lessThan(this.mThreshold)) {
                    break;
                }
                // Compute the weighted average tMax = (3*tMax+tMin)/4. Round
                // up to twice the precision to avoid quadratic growth in the
                // number of bits and to ensure that tMin can increase.
                tMax = BSRational.ldexp(this.mThree.mul(tMax).add(tMin), -2);
                tMax = convertBSRational(tMax, 2 * this.mPrecision,
                    BSNumberRoundingMode.FE_UPWARD);
                tMaxSqr = tMax.mul(tMax);
                tMin = a2pb2.mul(tMaxSqr).sub(a2mb2Sqr).div(tMax.mul(tMaxSqr.sub(a2pb2)));
            }
            return { numIterates: iterate, tMin, tMax };
        }

        // The sign of the second derivative is Sign(a^4-7*a^2*b^2+b^4) and
        // cannot be zero. Define the rational r = a^2/b^2 so that
        // a^4-7*a^2*b^2+b^4 = 0. This implies r^2 - 7*r + 1 = 0. The
        // irrational roots are r = (7 +- sqrt(45))/2, which is a
        // contradiction.
        logError('This second derivative cannot be zero at a-b.');
    }

    // Compute a bounding interval for the quadratic field number
    // q = x[0] + x[1] * sqrt(d), namely qMin <= q <= qMax, where the
    // endpoints are both within the specified precision.
    estimate(q: APConversionQFN1): { numIterates: number, qMin: BSRational, qMax: BSRational } {
        const x = q.x[0];
        const y = q.x[1];
        const d = q.d;

        let numIterates: number;
        let qMin: BSRational, qMax: BSRational;
        if (d.notEquals(this.mZero) && y.notEquals(this.mZero)) {
            const aSqr = y.mul(y).mul(d);
            const estimated = this.estimateSqrt(aSqr);
            numIterates = estimated.numIterates;
            if (y.greaterThan(this.mZero)) {
                qMin = x.add(estimated.aMin);
                qMax = x.add(estimated.aMax);
            } else {
                qMin = x.sub(estimated.aMax);
                qMax = x.sub(estimated.aMin);
            }
        } else {
            numIterates = 0;
            qMin = x.clone();
            qMax = x.clone();
        }

        return { numIterates, qMin, qMax };
    }

    // Compute an estimate of the quadratic field number when you do not need
    // a bounding interval.
    estimateValue(q: APConversionQFN1): { numIterates: number, qEstimate: BSRational } {
        // Compute a bounding interval qMin <= q <= qMax.
        const { numIterates, qMin, qMax } = this.estimate(q);
        // Use the average of the interval endpoints as the estimate.
        return { numIterates, qEstimate: BSRational.ldexp(qMin.add(qMax), -1) };
    }

    // Factor a^2 = r^2 * 2^e, where r^2 in [1/2,1), and return r^2 (possibly
    // doubled so that a = sqrt(rSqr) * 2^exponent).
    private static preprocessSqr(aSqr: BSRational): { rSqr: BSRational, exponent: number } {
        const { result, exponent: exponentASqr } = BSRational.frexp(aSqr);
        if ((exponentASqr & 1) !== 0) {  // odd exponent
            // a = sqrt(2*r^2) * 2^{(e-1)/2}, with 2*r^2 in [1,2). The
            // division is exact because e-1 is even.
            return {
                rSqr: BSRational.ldexp(result, 1),
                exponent: (exponentASqr - 1) / 2
            };
        } else {  // even exponent
            // a = sqrt(r^2) * 2^{e/2}, with r^2 in [1/2,1).
            return { rSqr: result, exponent: exponentASqr / 2 };
        }
    }

    // Compute a lower bound on the square root of r^2, scaled by 2^exponent.
    private static getMinOfSqrt(rSqr: BSRational, exponent: number): BSRational {
        const lowerRSqr = convertBSRationalToNumber(rSqr, BSNumberRoundingMode.FE_DOWNWARD);
        const sqrtLowerRSqr = Math.sqrt(lowerRSqr);
        const aMin = BSRational.fromNumber(nextDown(sqrtLowerRSqr));
        return BSRational.ldexp(aMin, exponent);
    }

    // Compute an upper bound on the square root of r^2, scaled by
    // 2^exponent.
    private static getMaxOfSqrt(rSqr: BSRational, exponent: number): BSRational {
        const upperRSqr = convertBSRationalToNumber(rSqr, BSNumberRoundingMode.FE_UPWARD);
        const sqrtUpperRSqr = Math.sqrt(upperRSqr);
        const aMax = BSRational.fromNumber(nextUp(sqrtUpperRSqr));
        return BSRational.ldexp(aMax, exponent);
    }
}
