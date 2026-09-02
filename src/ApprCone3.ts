// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The cone vertex is V, the unit-length axis direction is U and the cone
// angle is A in (0,pi/2). The cone is defined algebraically by those points
// X for which
//     Dot(U,X-V)/Length(X-V) = cos(A)
// This can be written as a quadratic equation
//     (V-X)^T * (cos(A)^2 - U * U^T) * (V-X) = 0
// with the implicit constraint that Dot(U, X-V) > 0 (X is on the "positive"
// cone). Define W = U/cos(A), so Length(W) > 1 and
//     F(X;V,W) = (V-X)^T * (I - W * W^T) * (V-X) = 0
// The nonlinear least squares fitting of points {X[i]}_{i=0}^{n-1} computes
// V and W to minimize the error function
//     E(V,W) = sum_{i=0}^{n-1} F(X[i];V,W)^2
// I recommend using the Gauss-Newton minimizer when your cone points are
// truly nearly a cone; otherwise, try the Levenberg-Marquardt minimizer.
//
// The mathematics used in this implementation are found in
//   https://www.geometrictools.com/Documentation/LeastSquaresFitting.pdf
//
// Port notes:
// * The two 'operator()' overloads become 'computeGaussNewton' and
//   'computeLevenbergMarquardt'; neither is more canonical than the other,
//   so both get descriptive names.
// * The in/out reference parameters (coneVertex, coneAxis, coneAngle) become
//   a single in/out ApprCone3Parameters object that the method reads (when
//   useConeInputAsInitialGuess is true) and overwrites with the fit, the
//   same convention used by ApprSphere3.ts. The minimizer Result is the
//   return value, as upstream.

import { ApprHeightLine2 } from './ApprHeightLine2';
import { GaussNewtonMinimizer } from './GaussNewtonMinimizer';
import type { GaussNewtonMinimizerResult } from './GaussNewtonMinimizer';
import { LevenbergMarquardtMinimizer } from './LevenbergMarquardtMinimizer';
import type {
    LevenbergMarquardtMinimizerResult
} from './LevenbergMarquardtMinimizer';
import { Matrix } from './Matrix';
import { Vector, dot, length, mul, normalize, sub } from './Vector';

// The cone parameters. On input to the fitting methods these are the initial
// guesses when useConeInputAsInitialGuess is true (in which case 'axis' is
// normalized in place). On output they are the fitted cone.
export interface ApprCone3Parameters {
    vertex: Vector;
    axis: Vector;
    angle: number;
}

export class ApprCone3 {
    private mPoints: readonly Vector[] = [];

    // F[i](V,W) = D^T * (I - W * W^T) * D, D = V - X[i], P = (V,W)
    private readonly mFFunction = (P: Vector, F: Vector): void => {
        const V = P.values;
        const w0 = P.values[3], w1 = P.values[4], w2 = P.values[5];
        for (let i = 0; i < this.mPoints.length; ++i) {
            const X = this.mPoints[i].values;
            const d0 = V[0] - X[0], d1 = V[1] - X[1], d2 = V[2] - X[2];
            const deltaDotW = d0 * w0 + d1 * w1 + d2 * w2;
            F.values[i] = d0 * d0 + d1 * d1 + d2 * d2 - deltaDotW * deltaDotW;
        }
    };

    // dF[i]/dV = 2 * (D - Dot(W, D) * W)
    // dF[i]/dW = -2 * Dot(W, D) * D
    private readonly mJFunction = (P: Vector, J: Matrix): void => {
        const two = 2;
        const V = P.values;
        const W = [P.values[3], P.values[4], P.values[5]];
        for (let row = 0; row < this.mPoints.length; ++row) {
            const X = this.mPoints[row].values;
            const delta = [V[0] - X[0], V[1] - X[1], V[2] - X[2]];
            const deltaDotW =
                delta[0] * W[0] + delta[1] * W[1] + delta[2] * W[2];
            for (let col = 0; col < 3; ++col) {
                const temp0 = delta[col] - deltaDotW * W[col];
                const temp1 = deltaDotW * delta[col];
                J.set(row, col, two * temp0);
                J.set(row, col + 3, -two * temp1);
            }
        }
    };

    // If you want to specify that cone.vertex, cone.axis and cone.angle are
    // the initial guesses for the minimizer, set the parameter
    // useConeInputAsInitialGuess to 'true'. If you want the function to
    // compute initial guesses, set that parameter to 'false'. A Gauss-Newton
    // minimizer is used to fit a cone using nonlinear least-squares. The
    // fitted cone is stored in 'cone'. See GaussNewtonMinimizer.ts for a
    // description of the least-squares algorithm and the parameters that it
    // requires.
    computeGaussNewton(points: readonly Vector[], maxIterations: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        useConeInputAsInitialGuess: boolean,
        cone: ApprCone3Parameters): GaussNewtonMinimizerResult {
        this.mPoints = points;
        const minimizer = GaussNewtonMinimizer.fromJFunction(6, points.length,
            this.mFFunction, this.mJFunction);

        const initial = this.computeInitialGuess(useConeInputAsInitialGuess,
            cone);

        const result = minimizer.minimize(initial, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance);

        // No test is made for result.converged so that we return some
        // estimates of the cone. The caller can decide how to respond when
        // result.converged is false.
        extractCone(result.minLocation, cone);

        this.mPoints = [];
        return result;
    }

    // The cone parameters are in/out variables; see computeGaussNewton. See
    // GaussNewtonMinimizer.ts for a description of the least-squares
    // algorithm and the parameters that it requires. (The file
    // LevenbergMarquardtMinimizer.ts directs you to the Gauss-Newton file to
    // read about the parameters.)
    computeLevenbergMarquardt(points: readonly Vector[], maxIterations: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        lambdaFactor: number, lambdaAdjust: number, maxAdjustments: number,
        useConeInputAsInitialGuess: boolean,
        cone: ApprCone3Parameters): LevenbergMarquardtMinimizerResult {
        this.mPoints = points;
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(6,
            points.length, this.mFFunction, this.mJFunction);

        const initial = this.computeInitialGuess(useConeInputAsInitialGuess,
            cone);

        const result = minimizer.minimize(initial, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance, lambdaFactor,
            lambdaAdjust, maxAdjustments);

        // No test is made for result.converged so that we return some
        // estimates of the cone. The caller can decide how to respond when
        // result.converged is false.
        extractCone(result.minLocation, cone);

        this.mPoints = [];
        return result;
    }

    // The initial guess (V,W) shared by the two minimizers.
    private computeInitialGuess(useConeInputAsInitialGuess: boolean,
        cone: ApprCone3Parameters): Vector {
        if (useConeInputAsInitialGuess) {
            normalize(cone.axis);
        } else {
            this.computeInitialCone(cone);
        }

        const initial = new Vector(6);

        // The initial guess for the cone vertex.
        initial.values[0] = cone.vertex.values[0];
        initial.values[1] = cone.vertex.values[1];
        initial.values[2] = cone.vertex.values[2];

        // The initial guess for the weighted cone axis.
        const coneCosAngle = Math.cos(cone.angle);
        initial.values[3] = cone.axis.values[0] / coneCosAngle;
        initial.values[4] = cone.axis.values[1] / coneCosAngle;
        initial.values[5] = cone.axis.values[2] / coneCosAngle;
        return initial;
    }

    private computeInitialCone(cone: ApprCone3Parameters): void {
        const points = this.mPoints;
        const numPoints = points.length;

        // Compute the average of the sample points.
        const center = new Vector(3);
        for (let i = 0; i < numPoints; ++i) {
            center.values[0] += points[i].values[0];
            center.values[1] += points[i].values[1];
            center.values[2] += points[i].values[2];
        }
        center.values[0] /= numPoints;
        center.values[1] /= numPoints;
        center.values[2] /= numPoints;

        // The cone axis is estimated from ZZTZ (see the PDF).
        let coneAxis = new Vector(3);
        for (let i = 0; i < numPoints; ++i) {
            const delta = sub(points[i], center);
            const dd = dot(delta, delta);
            coneAxis.values[0] += delta.values[0] * dd;
            coneAxis.values[1] += delta.values[1] * dd;
            coneAxis.values[2] += delta.values[2] * dd;
        }
        normalize(coneAxis);

        // Compute the signed heights of the points along the cone axis
        // relative to C. These are the projections of the points onto the
        // line C+t*U. Also compute the radial distances of the points from
        // the line C+t*U.
        const hrPairs: Vector[] = new Array<Vector>(numPoints);
        let hMin = Number.MAX_VALUE, hMax = -Number.MAX_VALUE;
        for (let i = 0; i < numPoints; ++i) {
            const delta = sub(points[i], center);
            const h = dot(coneAxis, delta);
            hMin = Math.min(hMin, h);
            hMax = Math.max(hMax, h);
            const projection = sub(delta, mul(coneAxis, h));
            const r = length(projection);
            hrPairs[i] = Vector.fromArray([h, r]);
        }

        // The radial distance is considered to be a function of height. Fit
        // the (h,r) pairs with a line:
        //   r - rAverage = hrSlope * (h - hAverage)
        const fitter = new ApprHeightLine2();
        fitter.fit(hrPairs);
        const parameters = fitter.getParameters();
        const hAverage = parameters.average.values[0];
        const rAverage = parameters.average.values[1];
        let hrSlope = parameters.coefficients.values[0];

        // If U is directed so that r increases as h increases, U is the
        // correct cone axis estimate. However, if r decreases as h increases,
        // -U is the correct cone axis estimate.
        if (hrSlope < 0) {
            coneAxis = mul(coneAxis, -1);
            hrSlope = -hrSlope;
            const temp = hMin;
            hMin = -hMax;
            hMax = -temp;
        }

        // Compute the extreme radial distance values for the points.
        const rMin = rAverage + hrSlope * (hMin - hAverage);
        const rMax = rAverage + hrSlope * (hMax - hAverage);
        const hRange = hMax - hMin;
        const rRange = rMax - rMin;

        // Using trigonometry and right triangles, compute the tangent
        // function of the cone angle.
        const tanAngle = rRange / hRange;
        cone.angle = Math.atan2(rRange, hRange);

        // Compute the cone vertex.
        const offset = rMax / tanAngle - hMax;
        cone.axis = coneAxis;
        cone.vertex = sub(center, mul(coneAxis, offset));
    }
}

// Extract (V,W) from the minimizer location and convert W to the unit-length
// cone axis and the cone angle.
function extractCone(minLocation: Vector, cone: ApprCone3Parameters): void {
    cone.vertex = Vector.fromArray([
        minLocation.values[0], minLocation.values[1], minLocation.values[2]
    ]);
    cone.axis = Vector.fromArray([
        minLocation.values[3], minLocation.values[4], minLocation.values[5]
    ]);

    // We know that coneCosAngle will be nonnegative. The Math.min call guards
    // against rounding errors leading to a number slightly larger than 1. The
    // clamping ensures Math.acos will not return a NaN.
    const coneCosAngle = Math.min(1 / normalize(cone.axis), 1);
    cone.angle = Math.acos(coneCosAngle);
}
