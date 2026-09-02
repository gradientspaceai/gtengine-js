// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprTorus3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Let the torus center be C with plane of symmetry containing C and having
// directions D0 and D1. The axis of symmetry is the line containing C and
// having direction N (the plane normal). The radius from the center of the
// torus is r0 and the radius of the tube of the torus is r1. A point P may
// be written as P = C + x*D0 + y*D1 + z*N, where matrix [D0 D1 N] is
// orthogonal and has determinant 1. Thus, x = Dot(D0,P-C), y = Dot(D1,P-C)
// and z = Dot(N,P-C). The implicit equation defining the torus is
//     (|P-C|^2 + r0^2 - r1^2)^2 - 4*r0^2*(|P-C|^2 - (Dot(N,P-C))^2) = 0
// Observe that D0 and D1 are not present in the equation, which is to be
// expected by the symmetry.
//
// Define u = r0^2 and v = r0^2 - r1^2. Define
//     F(X;C,N,u,v) = (|P-C|^2 + v)^2 - 4*u*(|P-C|^2 - (Dot(N,P-C))^2)
// The nonlinear least-squares fitting of points {X[i]}_{i=0}^{n-1} computes
// C, N, u and v to minimize the error function
//     E(C,N,u,v) = sum_{i=0}^{n-1} F(X[i];C,N,u,v)^2
// When the sample points are distributed so that there is large coverage
// by a purported fitted torus, a variation on fitting is the following.
// Compute the least-squares plane with origin C and normal N that fits the
// points. Define G(X;u,v) = F(X;C,N,u,v); the only variables now are u and
// v. Define L[i] = |X[i]-C|^2 and S[i] = 4 * (L[i] - (Dot(N,X[i]-C))^2).
// Define the error function
//     H(u,v) = sum_{i=0}^{n-1} G(X[i];u,v)^2
//            = sum_{i=0}^{n-1} ((v + L[i])^2 - S[i]*u)^2
// The first-order partial derivatives are
//     dH/du = -2 sum_{i=0}^{n-1} ((v + L[i])^2 - S[i]*u) * S[i]
//     dH/dv =  4 sum_{i=0}^{n-1} ((v + L[i])^2 - S[i]*u) * (v + L[i])
// Setting these to zero and expanding the terms, we have
//     0 = a2 * v^2 + a1 * v + a0 - b0 * u
//     0 = c3 * v^3 + c2 * v^2 + c1 * v + c0 - u * (d1 * v + d0)
// where a2 = sum(S[i]), a1 = 2*sum(S[i]*L[i]), a0 = sum(S[i]*L[i]^2),
// b0 = sum(S[i]^2), c3 = sum(1) = n, c2 = 3*sum(L[i]), c1 = 3*sum(L[i]^2),
// c0 = sum(L[i]^3), d1 = sum(S[i]) = a2 and d0 = sum(S[i]*L[i]) = a1/2.
// (Upstream's comment names the third coefficient a2 rather than a0; the
// code is correct.) The first equation is solved for
//     u = (a2 * v^2 + a1 * v + a0) / b0 = e2 * v^2 + e1 * v + e0
// and substituted into the second equation to obtain a cubic polynomial
// equation
//     0 = f3 * v^3 + f2 * v^2 + f1 * v + f0
// where f3 = c3 - d1 * e2, f2 = c2 - d1 * e1 - d0 * e2,
// f1 = c1 - d1 * e0 - d0 * e1 and f0 = c0 - d0 * e0. The positive v-roots
// are computed. For each root compute the corresponding u. For all pairs
// (u,v) with u > v > 0, evaluate H(u,v) and choose the pair that minimizes
// H(u,v). The torus radii are r0 = sqrt(u) and r1 = sqrt(u - v).
//
// Port notes:
// * Upstream converts the cubic coefficients f0..f3 to BSRational<UIntegerAP32>
//   so that RootsPolynomial::SolveCubic classifies the roots with exact
//   arithmetic. The port's RootsPolynomial is the double instantiation (the
//   precedent set by RootsPolynomial.ts and TriangulateEC.ts), so the
//   coefficients are passed as doubles. The structure is unchanged, so an
//   exact path can be added when BSRational is ported.
// * The three 'operator()' overloads become 'compute' (the plane-based
//   estimate), 'computeGaussNewton' and 'computeLevenbergMarquardt'.
// * The in/out reference parameters (C, N, r0, r1) become a single in/out
//   ApprTorus3Parameters object, the same convention used by ApprSphere3.ts.

import { ApprOrthogonalPlane3 } from './ApprOrthogonalPlane3';
import { GaussNewtonMinimizer } from './GaussNewtonMinimizer';
import type { GaussNewtonMinimizerResult } from './GaussNewtonMinimizer';
import { LevenbergMarquardtMinimizer } from './LevenbergMarquardtMinimizer';
import type {
    LevenbergMarquardtMinimizerResult
} from './LevenbergMarquardtMinimizer';
import { Matrix } from './Matrix';
import { RootsPolynomial } from './RootsPolynomial';
import { Vector, dot, sub } from './Vector';

// The torus parameters: the center C, the unit-length axis of symmetry N,
// the radius r0 from the center and the radius r1 of the tube.
export interface ApprTorus3Parameters {
    C: Vector;
    N: Vector;
    r0: number;
    r1: number;
}

// The result of ApprTorus3.compute. The 'success' member is 'true' when the
// estimate is valid, in which case 'error' is the least-squares error for
// that estimate. If any unexpected condition occurs that prevents computing
// an estimate, 'success' is false and 'error' is Number.MAX_VALUE.
export interface ApprTorus3ComputeResult {
    success: boolean;
    error: number;
}

export class ApprTorus3 {
    private mPoints: readonly Vector[] = [];

    // The unit-length normal is
    //   N = (cos(theta)*sin(phi), sin(theta)*sin(phi), cos(phi))
    // for theta in [0,2*pi) and phi in [0,pi). The radii are encoded as
    //   u = r0^2, v = r0^2 - r1^2
    // with 0 < v < u. Let D = C - X[i] where X[i] is a sample point. The
    // parameters are P = (C0,C1,C2,theta,phi,u,v).

    // F[i](C,theta,phi,u,v) = (|D|^2 + v)^2 - 4*u*(|D|^2 - Dot(N,D)^2)
    private readonly mFFunction = (P: Vector, F: Vector): void => {
        const p = P.values;
        const csTheta = Math.cos(p[3]);
        const snTheta = Math.sin(p[3]);
        const csPhi = Math.cos(p[4]);
        const snPhi = Math.sin(p[4]);
        const n0 = csTheta * snPhi, n1 = snTheta * snPhi, n2 = csPhi;
        const u = p[5];
        const v = p[6];
        for (let i = 0; i < this.mPoints.length; ++i) {
            const X = this.mPoints[i].values;
            const d0 = p[0] - X[0], d1 = p[1] - X[1], d2 = p[2] - X[2];
            const DdotD = d0 * d0 + d1 * d1 + d2 * d2;
            const NdotD = n0 * d0 + n1 * d1 + n2 * d2;
            const sum = DdotD + v;
            F.values[i] = sum * sum - 4 * u * (DdotD - NdotD * NdotD);
        }
    };

    // dF[i]/dC = 4 * (|D|^2 + v) * D - 8 * u * (I - N*N^T) * D
    // dF[i]/dTheta = 8 * u * Dot(N,D) * Dot(dN/dTheta, D)
    // dF[i]/dPhi = 8 * u * Dot(N,D) * Dot(dN/dPhi, D)
    // dF[i]/du = -4 * (|D|^2 - Dot(N,D)^2)
    // dF[i]/dv = 2 * (|D|^2 + v)
    private readonly mJFunction = (P: Vector, J: Matrix): void => {
        const p = P.values;
        const csTheta = Math.cos(p[3]);
        const snTheta = Math.sin(p[3]);
        const csPhi = Math.cos(p[4]);
        const snPhi = Math.sin(p[4]);
        const n0 = csTheta * snPhi, n1 = snTheta * snPhi, n2 = csPhi;
        const u = p[5];
        const v = p[6];
        const dNdTheta = [-snTheta * snPhi, csTheta * snPhi, 0];
        const dNdPhi = [csTheta * csPhi, snTheta * csPhi, -snPhi];
        for (let row = 0; row < this.mPoints.length; ++row) {
            const X = this.mPoints[row].values;
            const d0 = p[0] - X[0], d1 = p[1] - X[1], d2 = p[2] - X[2];
            const DdotD = d0 * d0 + d1 * d1 + d2 * d2;
            const NdotD = n0 * d0 + n1 * d1 + n2 * d2;
            const sum = DdotD + v;
            J.set(row, 0, 4 * sum * d0 - 8 * u * (d0 - NdotD * n0));
            J.set(row, 1, 4 * sum * d1 - 8 * u * (d1 - NdotD * n1));
            J.set(row, 2, 4 * sum * d2 - 8 * u * (d2 - NdotD * n2));
            J.set(row, 3, 8 * u * NdotD *
                (dNdTheta[0] * d0 + dNdTheta[1] * d1 + dNdTheta[2] * d2));
            J.set(row, 4, 8 * u * NdotD *
                (dNdPhi[0] * d0 + dNdPhi[1] * d1 + dNdPhi[2] * d2));
            J.set(row, 5, -4 * (DdotD - NdotD * NdotD));
            J.set(row, 6, 2 * sum);
        }
    };

    // When the samples are distributed approximately uniformly near a torus,
    // use this method. For example, if the purported torus has center (0,0,0)
    // and normal (0,0,1), you want the (x,y,z) samples to occur in all 8
    // octants. If the samples occur, say, only in one octant, this method
    // will estimate a C and N that are nowhere near (0,0,0) and (0,0,1). The
    // function sets the members of 'torus' as the fitted torus.
    compute(points: readonly Vector[],
        torus: ApprTorus3Parameters): ApprTorus3ComputeResult {
        const numPoints = points.length;
        const fitter = new ApprOrthogonalPlane3();
        if (!fitter.fit(points)) {
            return { success: false, error: Number.MAX_VALUE };
        }
        // Clone the fitter parameters; upstream copies by value.
        torus.C = fitter.getParameters().origin.clone();
        torus.N = fitter.getParameters().normal.clone();
        const C = torus.C, N = torus.N;

        let a0 = 0, a1 = 0, a2 = 0, b0 = 0;
        let c0 = 0, c1 = 0, c2 = 0;
        const c3 = numPoints;
        for (let i = 0; i < numPoints; ++i) {
            const delta = sub(points[i], C);
            const d = dot(N, delta);
            const L = dot(delta, delta), L2 = L * L, L3 = L * L2;
            const S = 4 * (L - d * d), S2 = S * S;
            a2 += S;
            a1 += S * L;
            a0 += S * L2;
            b0 += S2;
            c2 += L;
            c1 += L2;
            c0 += L3;
        }
        const d1 = a2;
        const d0 = a1;
        a1 *= 2;
        c2 *= 3;
        c1 *= 3;
        const invB0 = 1 / b0;
        const e0 = a0 * invB0;
        const e1 = a1 * invB0;
        const e2 = a2 * invB0;

        const f0 = c0 - d0 * e0;
        const f1 = c1 - d1 * e0 - d0 * e1;
        const f2 = c2 - d1 * e1 - d0 * e2;
        const f3 = c3 - d1 * e2;
        const rmMap = RootsPolynomial.solveCubic(f0, f1, f2, f3);

        let hmin = Number.MAX_VALUE;
        let umin = 0, vmin = 0;
        for (const element of rmMap) {
            const v = element.root;
            if (v > 0) {
                const u = e0 + v * (e1 + v * e2);
                if (u > v) {
                    let h = 0;
                    for (let i = 0; i < numPoints; ++i) {
                        const delta = sub(points[i], C);
                        const d = dot(N, delta);
                        const L = dot(delta, delta);
                        const S = 4 * (L - d * d);
                        const sum = v + L;
                        const term = sum * sum - S * u;
                        h += term * term;
                    }
                    if (h < hmin) {
                        hmin = h;
                        umin = u;
                        vmin = v;
                    }
                }
            }
        }

        if (hmin === Number.MAX_VALUE) {
            return { success: false, error: Number.MAX_VALUE };
        }

        torus.r0 = Math.sqrt(umin);
        torus.r1 = Math.sqrt(umin - vmin);
        return { success: true, error: hmin };
    }

    // If you want to specify that the members of 'torus' are the initial
    // guesses for the minimizer, set useTorusInputAsInitialGuess to 'true'.
    // If you want the function to compute initial guesses, set it to 'false'.
    // A Gauss-Newton minimizer is used to fit a torus using nonlinear
    // least-squares. The fitted torus is stored in 'torus'. See
    // GaussNewtonMinimizer.ts for a description of the least-squares
    // algorithm and the parameters that it requires.
    computeGaussNewton(points: readonly Vector[], maxIterations: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        useTorusInputAsInitialGuess: boolean,
        torus: ApprTorus3Parameters): GaussNewtonMinimizerResult {
        this.mPoints = points;
        const minimizer = GaussNewtonMinimizer.fromJFunction(7, points.length,
            this.mFFunction, this.mJFunction);

        if (!useTorusInputAsInitialGuess) {
            this.compute(points, torus);
        }

        const initial = computeInitialGuess(torus);

        const result = minimizer.minimize(initial, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance);

        // No test is made for result.converged so that we return some
        // estimates of the torus. The caller can decide how to respond when
        // result.converged is false.
        extractTorus(result.minLocation, torus);

        this.mPoints = [];
        return result;
    }

    // The torus parameters are in/out variables; see computeGaussNewton. A
    // Levenberg-Marquardt minimizer is used to fit a torus using nonlinear
    // least-squares.
    computeLevenbergMarquardt(points: readonly Vector[], maxIterations: number,
        updateLengthTolerance: number, errorDifferenceTolerance: number,
        lambdaFactor: number, lambdaAdjust: number, maxAdjustments: number,
        useTorusInputAsInitialGuess: boolean,
        torus: ApprTorus3Parameters): LevenbergMarquardtMinimizerResult {
        this.mPoints = points;
        const minimizer = LevenbergMarquardtMinimizer.fromJFunction(7,
            points.length, this.mFFunction, this.mJFunction);

        if (!useTorusInputAsInitialGuess) {
            this.compute(points, torus);
        }

        const initial = computeInitialGuess(torus);

        const result = minimizer.minimize(initial, maxIterations,
            updateLengthTolerance, errorDifferenceTolerance, lambdaFactor,
            lambdaAdjust, maxAdjustments);

        // No test is made for result.converged so that we return some
        // estimates of the torus. The caller can decide how to respond when
        // result.converged is false.
        extractTorus(result.minLocation, torus);

        this.mPoints = [];
        return result;
    }
}

// The initial guess P = (C0,C1,C2,theta,phi,u,v) shared by the minimizers.
function computeInitialGuess(torus: ApprTorus3Parameters): Vector {
    const initial = new Vector(7);

    // The initial guess for the plane origin.
    initial.values[0] = torus.C.values[0];
    initial.values[1] = torus.C.values[1];
    initial.values[2] = torus.C.values[2];

    // The initial guess for the plane normal. The angles must be extracted
    // for spherical coordinates.
    if (Math.abs(torus.N.values[2]) < 1) {
        initial.values[3] = Math.atan2(torus.N.values[1], torus.N.values[0]);
        initial.values[4] = Math.acos(torus.N.values[2]);
    } else {
        initial.values[3] = 0;
        initial.values[4] = 0;
    }

    // The initial guess for the radii-related parameters.
    initial.values[5] = torus.r0 * torus.r0;
    initial.values[6] = initial.values[5] - torus.r1 * torus.r1;
    return initial;
}

// Convert the minimizer location to the torus parameters.
function extractTorus(minLocation: Vector,
    torus: ApprTorus3Parameters): void {
    const p = minLocation.values;
    torus.C = Vector.fromArray([p[0], p[1], p[2]]);

    const theta = p[3];
    const phi = p[4];
    const csTheta = Math.cos(theta);
    const snTheta = Math.sin(theta);
    const csPhi = Math.cos(phi);
    const snPhi = Math.sin(phi);
    torus.N = Vector.fromArray([
        csTheta * snPhi, snTheta * snPhi, csPhi
    ]);

    const u = p[5];
    const v = p[6];
    torus.r0 = Math.sqrt(u);
    torus.r1 = Math.sqrt(u - v);
}
