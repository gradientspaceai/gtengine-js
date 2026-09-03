// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprCone3EllipseAndPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An infinite single-sided cone is fit to a 3D ellipse that is known to be
// the intersection of a plane with the cone. The ellipse itself is not
// enough information, producing the cone vertex and cone direction as a
// function of the cone angle. Additional points on the cone are required
// to determine the cone angle. The algorithm description is
// https://www.geometrictools.com/Documentation/FitConeToEllipseAndPoints.pdf
//
// Port notes:
// * The upstream header contains two classes. Only 'ApprCone3EllipseAndPoints'
//   is ported here. The companion 'ApprCone3ExtractEllipses' requires
//   OrientedBoxTreeOfPoints.h, which is not yet part of the port; it is
//   deferred until that dependency lands.
// * The nested 'Control' struct becomes the exported class
//   'ApprCone3EllipseAndPointsControl' (global export uniqueness), whose
//   constructor supplies the upstream default values.
// * The private static 'ComputeCone' is a module-private function.

import { Cone } from './Cone';
import type { Cone3 } from './Cone';
import { GTE_C_HALF_PI } from './Constants';
import { Ellipse3 } from './Ellipse3';
import { logAssert } from './Logger';
import { Minimize1 } from './Minimize1';
import { Vector, add, dot, mul, sub } from './Vector';

// The default control parameters appear to be reasonable for applications,
// but they are exposed to the caller for tuning.
export class ApprCone3EllipseAndPointsControl {
    // The least-squares error function is updated with the penalty value for
    // a points[i] that is below the plane supporting the cone; that is, when
    // the dot product Dot(coneDirection, points[i] - coneVertex) < 0.
    penalty: number;

    // Parameters for Minimize1.
    maxSubdivisions: number;
    maxBisections: number;
    epsilon: number;
    tolerance: number;

    // Search for the minimum on [0 + padding, pi/2 - padding] to avoid
    // divisions by zero of the least-squares error function at the endpoints
    // of [0,pi/2].
    padding: number;

    constructor() {
        this.penalty = 1;
        this.maxSubdivisions = 8;
        this.maxBisections = 64;
        this.epsilon = 1e-08;
        this.tolerance = 1e-04;
        this.padding = 1e-03;
    }

    validParameters(): boolean {
        return this.penalty > 0
            && this.maxSubdivisions > 0
            && this.maxBisections > 0
            && this.epsilon > 0
            && this.tolerance > 0
            && this.padding > 0;
    }
}

export class ApprCone3EllipseAndPoints {
    // The ellipse must be the intersection of a plane with the cone. In an
    // application, typically the ellipse is estimated from point samples of
    // the intersection which are then fitted with the ellipse.
    static fit(ellipse: Ellipse3, points: readonly Vector[],
        control: ApprCone3EllipseAndPointsControl =
            new ApprCone3EllipseAndPointsControl()): Cone3 {
        logAssert(control.validParameters(), 'Invalid control parameter.');

        // Upstream divides the accumulated error by points.size() and the
        // ellipse extents appear in denominators of ComputeCone, neither of
        // which is guarded upstream. The port asserts, following the
        // established Appr*/Cont* precedent of an explicit empty-input guard.
        logAssert(points.length > 0, 'ApprCone3EllipseAndPoints: no points.');
        logAssert(ellipse.extent.get(0) > 0 && ellipse.extent.get(1) > 0,
            'ApprCone3EllipseAndPoints: ellipse extents must be positive.');

        // The sign pair (sigma0,sigma1) is captured by the error function and
        // is varied between the four minimizer runs, as upstream does with
        // its lambda captures by reference.
        let sigma0 = 0;
        let sigma1 = 0;

        const F = (theta: number): number => {
            const cone = computeCone(theta, sigma0, sigma1, ellipse);

            let error = 0;
            for (const point of points) {
                const diff = sub(point, cone.ray.origin);
                const d = dot(cone.ray.direction, diff);
                if (d >= 0) {
                    const sqrLen = dot(diff, diff);
                    const quad = d * d - cone.cosAngleSqr * sqrLen;
                    error += quad * quad;
                }
                else {
                    error += control.penalty;
                }
            }

            return Math.sqrt(error) / points.length;
        };

        const minimizer = new Minimize1(F, control.maxSubdivisions,
            control.maxBisections, control.epsilon, control.tolerance);
        const t0 = control.padding;
        const t1 = GTE_C_HALF_PI - control.padding;
        let minError = -1;
        let minCone = new Cone(3);

        const signs: ReadonlyArray<readonly [number, number]> =
            [[1, 1], [1, -1], [-1, 1], [-1, -1]];
        for (const [s0, s1] of signs) {
            sigma0 = s0;
            sigma1 = s1;
            const { tMin, fMin } = minimizer.getMinimum(t0, t1);
            if (t0 < tMin && tMin < t1) {
                if (minError === -1 || fMin < minError) {
                    minError = fMin;
                    minCone = computeCone(tMin, sigma0, sigma1, ellipse);
                }
            }
        }

        logAssert(minError !== -1, 'Failed to find fitted cone.');

        return minCone;
    }
}

// For a cone angle theta, the ellipse determines the cone vertex and cone
// axis direction up to the sign choices (sigma0,sigma1) for sin(phi) and
// cos(phi), where phi is the angle between the cone axis and the ellipse
// plane normal. The relationship sin(phi) = e * cos(theta) holds, where e is
// the eccentricity of the ellipse.
function computeCone(theta: number, sigma0: number, sigma1: number,
    ellipse: Ellipse3): Cone3 {
    const C = ellipse.center;
    const N = ellipse.normal;
    const U = ellipse.axis[0];
    const a = ellipse.extent.get(0);
    const b = ellipse.extent.get(1);
    const bDivA = b / a;
    const eSqr = Math.max(0, 1 - bDivA * bDivA);
    const omesqr = 1 - eSqr;
    const e = Math.sqrt(eSqr);

    const snTheta = Math.sin(theta);
    const csTheta = Math.cos(theta);
    const snPhi = sigma0 * e * csTheta;
    const snPhiSqr = snPhi * snPhi;
    const csPhi = sigma1 * Math.sqrt(Math.max(0, 1 - snPhiSqr));
    const h = a * omesqr * csTheta / (snTheta * Math.abs(csPhi));
    const D = add(mul(csPhi, N), mul(snPhi, U));
    const snThetaSqr = snTheta * snTheta;
    const csThetaSqr = csTheta * csTheta;
    const Q = sub(C, mul((h * snPhi * snThetaSqr) / (csThetaSqr - snPhiSqr), U));
    const K = sub(Q, mul(h, D));

    const cone = new Cone(3);
    cone.makeInfiniteCone();
    cone.setAngle(theta);
    cone.ray.origin = K;
    cone.ray.direction = D;
    return cone;
}
