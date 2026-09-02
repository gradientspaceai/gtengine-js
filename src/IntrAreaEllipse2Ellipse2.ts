// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAreaEllipse2Ellipse2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the area of intersection for two ellipses in 2D. The algorithm is
// discussed in the document
//   https://www.geometrictools.com/Documentation/AreaIntersectingEllipses.pdf
//
// Port notes: the upstream class is AreaEllipse2Ellipse2 (it is not a
// TIQuery/FIQuery specialization), so the class keeps that name and
// operator() becomes compute(). The nested Result struct becomes the exported
// interface AreaEllipse2Ellipse2Result and the nested Configuration enum
// class becomes the exported enum AreaEllipse2Ellipse2Configuration. The
// ellipse type Ellipse2<T> is the port's Hyperellipsoid of dimension 2.

import { GTE_C_PI, GTE_C_TWO_PI } from './Constants';
import type { Hyperellipsoid } from './Hyperellipsoid';
import {
    IntrEllipse2Ellipse2FI, defaultIntrEllipse2Ellipse2FIResult
} from './IntrEllipse2Ellipse2';
import type { IntrEllipse2Ellipse2FIResult } from './IntrEllipse2Ellipse2';
import { logAssert } from './Logger';
import { Matrix, addMatrix, divMatrix, mulMatrix, outerProduct } from './Matrix';
import { Vector, dot, sub } from './Vector';
import { dotPerp } from './Vector2';

// The configuration of the two ellipses.
export enum AreaEllipse2Ellipse2Configuration {
    ELLIPSES_ARE_EQUAL,
    ELLIPSES_ARE_SEPARATED,
    E0_CONTAINS_E1,
    E1_CONTAINS_E0,
    ONE_CHORD_REGION,
    FOUR_CHORD_REGION,
    INVALID
}

// The result of the area-of-intersection query.
export interface AreaEllipse2Ellipse2Result {
    // One of the enumerates, determined in the call to areaDispatch.
    configuration: AreaEllipse2Ellipse2Configuration;

    // Information about the ellipse-ellipse intersection points.
    findResult: IntrEllipse2Ellipse2FIResult;

    // The area of intersection of the ellipses.
    area: number;
}

// The port of the upstream Result default constructor.
export function defaultAreaEllipse2Ellipse2Result(): AreaEllipse2Ellipse2Result {
    return {
        configuration: AreaEllipse2Ellipse2Configuration.INVALID,
        findResult: defaultIntrEllipse2Ellipse2FIResult(),
        area: 0
    };
}

// The precomputed information for one ellipse. This is the port of the
// upstream private struct EllipseInfo.
interface EllipseInfo {
    center: Vector;
    axis: [Vector, Vector];
    extent: Vector;
    M: Matrix;
    AB: number;      // extent[0] * extent[1]
    halfAB: number;  // extent[0] * extent[1] / 2
    BpA: number;     // extent[1] + extent[0]
    BmA: number;     // extent[1] - extent[0]
}

export class AreaEllipse2Ellipse2 {
    // Compute the area of intersection of the ellipses.
    compute(ellipse0: Hyperellipsoid, ellipse1: Hyperellipsoid):
        AreaEllipse2Ellipse2Result {
        logAssert(ellipse0.dimension === 2 && ellipse1.dimension === 2,
            'AreaEllipse2Ellipse2: mismatched sizes.');

        const E0 = makeEllipseInfo(ellipse0);
        const E1 = makeEllipseInfo(ellipse1);

        const ar = defaultAreaEllipse2Ellipse2Result();
        ar.configuration = AreaEllipse2Ellipse2Configuration.INVALID;
        ar.findResult = new IntrEllipse2Ellipse2FI().find(ellipse0, ellipse1);
        ar.area = 0;
        this.areaDispatch(E0, E1, ar);
        return ar;
    }

    private areaDispatch(E0: EllipseInfo, E1: EllipseInfo,
        ar: AreaEllipse2Ellipse2Result): void {
        if (ar.findResult.intersect) {
            if (ar.findResult.numPoints === 1) {
                // Containment or separation.
                this.areaCS(E0, E1, ar);
            }
            else if (ar.findResult.numPoints === 2) {
                if (ar.findResult.isTransverse[0]) {
                    // Both intersection points are transverse.
                    this.area2(E0, E1, 0, 1, ar);
                }
                else {
                    // Both intersection points are tangential, so one ellipse
                    // is contained in the other.
                    this.areaCS(E0, E1, ar);
                }
            }
            else if (ar.findResult.numPoints === 3) {
                // The tangential intersection is irrelevant in the area
                // computation.
                if (!ar.findResult.isTransverse[0]) {
                    this.area2(E0, E1, 1, 2, ar);
                }
                else if (!ar.findResult.isTransverse[1]) {
                    this.area2(E0, E1, 2, 0, ar);
                }
                else {  // ar.findResult.isTransverse[2] === false
                    this.area2(E0, E1, 0, 1, ar);
                }
            }
            else {  // ar.findResult.numPoints === 4
                this.area4(E0, E1, ar);
            }
        }
        else {
            // Containment, separation, or same ellipse.
            this.areaCS(E0, E1, ar);
        }
    }

    private areaCS(E0: EllipseInfo, E1: EllipseInfo,
        ar: AreaEllipse2Ellipse2Result): void {
        if (ar.findResult.numPoints <= 1) {
            const diff = sub(E0.center, E1.center);
            const qform0 = dot(diff, mulMatrix(E0.M, diff));
            const qform1 = dot(diff, mulMatrix(E1.M, diff));
            if (qform0 > 1 && qform1 > 1) {
                // Each ellipse center is outside the other ellipse, so the
                // ellipses are separated (numPoints == 0) or outside each
                // other and just touching (numPoints == 1).
                ar.configuration =
                    AreaEllipse2Ellipse2Configuration.ELLIPSES_ARE_SEPARATED;
                ar.area = 0;
            }
            else {
                // One ellipse is inside the other. Determine this indirectly
                // by comparing areas.
                if (E0.AB < E1.AB) {
                    ar.configuration =
                        AreaEllipse2Ellipse2Configuration.E1_CONTAINS_E0;
                    ar.area = GTE_C_PI * E0.AB;
                }
                else {
                    ar.configuration =
                        AreaEllipse2Ellipse2Configuration.E0_CONTAINS_E1;
                    ar.area = GTE_C_PI * E1.AB;
                }
            }
        }
        else {
            ar.configuration =
                AreaEllipse2Ellipse2Configuration.ELLIPSES_ARE_EQUAL;
            ar.area = GTE_C_PI * E0.AB;
        }
    }

    private area2(E0: EllipseInfo, E1: EllipseInfo, i0: number, i1: number,
        ar: AreaEllipse2Ellipse2Result): void {
        ar.configuration = AreaEllipse2Ellipse2Configuration.ONE_CHORD_REGION;

        // The endpoints of the chord.
        const P0 = ar.findResult.points[i0];
        const P1 = ar.findResult.points[i1];

        // Compute locations relative to the ellipses.
        const P0mC0 = sub(P0, E0.center);
        const P0mC1 = sub(P0, E1.center);
        const P1mC0 = sub(P1, E0.center);
        const P1mC1 = sub(P1, E1.center);

        // Compute the ellipse normal vectors at endpoint P0. This is
        // sufficient information to determine chord endpoint order.
        const N0 = mulMatrix(E0.M, P0mC0);
        const N1 = mulMatrix(E1.M, P0mC1);
        const dp = dotPerp(N1, N0);

        // Choose the endpoint order for the chord region associated with E0.
        if (dp > 0) {
            // The chord order for E0 is <P0,P1> and for E1 is <P1,P0>.
            ar.area =
                computeAreaChordRegion(E0, P0mC0, P1mC0) +
                computeAreaChordRegion(E1, P1mC1, P0mC1);
        }
        else {
            // The chord order for E0 is <P1,P0> and for E1 is <P0,P1>.
            ar.area =
                computeAreaChordRegion(E0, P1mC0, P0mC0) +
                computeAreaChordRegion(E1, P0mC1, P1mC1);
        }
    }

    private area4(E0: EllipseInfo, E1: EllipseInfo,
        ar: AreaEllipse2Ellipse2Result): void {
        ar.configuration = AreaEllipse2Ellipse2Configuration.FOUR_CHORD_REGION;

        // Select a counterclockwise ordering of the points of intersection.
        // Use the polar coordinates for E0 to do this. Upstream uses a
        // std::multimap in the event that computing the intersections
        // involved numerical rounding errors that lead to a duplicate
        // intersection, even though the intersections are all labeled as
        // transverse. The port sorts an array by angle; Array.prototype.sort
        // is stable, so equal angles keep their insertion order as the
        // multimap does.
        const ordering: { theta: number, index: number }[] = [];
        for (let i = 0; i < 4; ++i) {
            const PmC = sub(ar.findResult.points[i], E0.center);
            const x = dot(E0.axis[0], PmC);
            const y = dot(E0.axis[1], PmC);
            ordering.push({ theta: Math.atan2(y, x), index: i });
        }
        ordering.sort((a, b) => a.theta - b.theta);

        const permute: number[] = ordering.map(element => element.index);

        // Start with the area of the convex quadrilateral.
        const diag20 = sub(ar.findResult.points[permute[2]],
            ar.findResult.points[permute[0]]);
        const diag31 = sub(ar.findResult.points[permute[3]],
            ar.findResult.points[permute[1]]);
        ar.area = Math.abs(dotPerp(diag20, diag31)) / 2;

        // Visit each pair of consecutive points. The selection of ellipse for
        // the chord-region area calculation uses the "most counterclockwise"
        // tangent vector.
        for (let i0 = 3, i1 = 0; i1 < 4; i0 = i1++) {
            // Get a pair of consecutive points.
            const P0 = ar.findResult.points[permute[i0]];
            const P1 = ar.findResult.points[permute[i1]];

            // Compute locations relative to the ellipses.
            const P0mC0 = sub(P0, E0.center);
            const P0mC1 = sub(P0, E1.center);
            const P1mC0 = sub(P1, E0.center);
            const P1mC1 = sub(P1, E1.center);

            // Compute the ellipse normal vectors at endpoint P0.
            const N0 = mulMatrix(E0.M, P0mC0);
            const N1 = mulMatrix(E1.M, P0mC1);
            const dp = dotPerp(N1, N0);
            if (dp > 0) {
                // The chord goes with ellipse E0.
                ar.area += computeAreaChordRegion(E0, P0mC0, P1mC0);
            }
            else {
                // The chord goes with ellipse E1.
                ar.area += computeAreaChordRegion(E1, P0mC1, P1mC1);
            }
        }
    }
}

// The port of the upstream FinishEllipseInfo, which computes the derived
// members of an EllipseInfo.
//
// Port fix for an upstream documentation/robustness defect. The upstream
// comment for operator() states "The ellipse axes are not required to be
// normalized." The matrix M is indeed computed in an axis-length-invariant
// way, but the polar angles theta = atan2(Dot(axis[1],X), Dot(axis[0],X)) and
// the sector-area integral are correct only for unit-length axes: with axes
// of differing lengths the angles are distorted and the reported areas are
// wrong. The port normalizes copies of the axes here, which leaves the
// results unchanged for the documented (unit-length) input and makes the
// stated contract true.
function makeEllipseInfo(ellipse: Hyperellipsoid): EllipseInfo {
    const axis0 = ellipse.axis[0].clone();
    const axis1 = ellipse.axis[1].clone();
    const M = addMatrix(
        divMatrix(outerProduct(axis0, axis0),
            ellipse.extent.values[0] * ellipse.extent.values[0] *
            dot(axis0, axis0)),
        divMatrix(outerProduct(axis1, axis1),
            ellipse.extent.values[1] * ellipse.extent.values[1] *
            dot(axis1, axis1)));

    const length0 = Math.sqrt(dot(axis0, axis0));
    const length1 = Math.sqrt(dot(axis1, axis1));
    logAssert(length0 > 0 && length1 > 0,
        'AreaEllipse2Ellipse2: degenerate ellipse axes.');
    for (let d = 0; d < 2; ++d) {
        axis0.values[d] /= length0;
        axis1.values[d] /= length1;
    }

    const AB = ellipse.extent.values[0] * ellipse.extent.values[1];
    return {
        center: ellipse.center.clone(),
        axis: [axis0, axis1],
        extent: ellipse.extent.clone(),
        M,
        AB,
        halfAB: AB / 2,
        BpA: ellipse.extent.values[1] + ellipse.extent.values[0],
        BmA: ellipse.extent.values[1] - ellipse.extent.values[0]
    };
}

// Compute the area of the region of the ellipse E bounded by the chord from
// P0 to P1 (both relative to the ellipse center) and the elliptical arc
// traversed counterclockwise from P0 to P1.
function computeAreaChordRegion(E: EllipseInfo, P0mC: Vector,
    P1mC: Vector): number {
    // Compute polar coordinates for P0 and P1 on the ellipse.
    const x0 = dot(E.axis[0], P0mC);
    const y0 = dot(E.axis[1], P0mC);
    let theta0 = Math.atan2(y0, x0);
    const x1 = dot(E.axis[0], P1mC);
    const y1 = dot(E.axis[1], P1mC);
    let theta1 = Math.atan2(y1, x1);

    // The arc straddles the atan2 discontinuity on the negative x-axis. Wrap
    // the second angle to be larger than the first angle.
    if (theta1 < theta0) {
        theta1 += GTE_C_TWO_PI;
    }

    // Compute the area portion of the sector due to the triangle.
    const triArea = Math.abs(dotPerp(P0mC, P1mC)) / 2;

    // Compute the chord region area.
    const dtheta = theta1 - theta0;
    if (dtheta <= GTE_C_PI) {
        // Use the area formula directly.
        //   area(theta0,theta1) = F(theta1) - F(theta0) - area(triangle)
        const F0 = computeIntegral(E, theta0);
        const F1 = computeIntegral(E, theta1);
        const sectorArea = F1 - F0;
        return sectorArea - triArea;
    }
    else {
        // The angle of the elliptical sector is larger than pi radians. Use
        // the area formula
        //   area(theta0,theta1) = pi*a*b - area(theta1,theta0)
        theta0 += GTE_C_TWO_PI;  // ensure theta0 > theta1
        const F0 = computeIntegral(E, theta0);
        const F1 = computeIntegral(E, theta1);
        const sectorArea = F0 - F1;
        return GTE_C_PI * E.AB - (sectorArea - triArea);
    }
}

// The antiderivative used by computeAreaChordRegion.
function computeIntegral(E: EllipseInfo, theta: number): number {
    const twoTheta = 2 * theta;
    const sn = Math.sin(twoTheta);
    const cs = Math.cos(twoTheta);
    const arg = E.BmA * sn / (E.BpA + E.BmA * cs);
    return E.halfAB * (theta - Math.atan(arg));
}
