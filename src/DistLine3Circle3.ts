// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The 3D line-circle distance algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
// The notation used in the code matches that of the document. The circle has
// center C and the plane of the circle has unit-length normal N. The line has
// origin B and non-zero direction M. The parameterization is P(t) = t*M+B. It
// is not necessary that M be a unit-length vector. This allows for the
// line-circle query to be used in the segment-circle query for the two-point
// form of a segment where M is the difference of the endpoints, which avoids
// a normalization of M that has numerical rounding errors.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line3<T>, Circle3<T>>' becomes the
// class DistLine3Circle3 with the result type DistLine3Circle3Result. The
// private 'Critical' struct that upstream shares with its friend classes
// DCPQuery<T, Ray3, Circle3> and DCPQuery<T, Segment3, Circle3> is exported
// as DistLine3Circle3Critical, and the friend-accessible 'Execute' is
// exported as distLine3Circle3Execute (the Dist* precedent for
// friend-granted helpers). The remaining private members become
// module-private functions.
//
// The upstream bisection is always performed in double precision, which the
// port matches trivially since 'number' is IEEE f64.

import type { Circle3 } from './Circle3';
import type { DCPQuery } from './DCPQuery';
import type { Line3 } from './Line';
import { RootsBisection1 } from './RootsBisection1';
import {
    Vector, add, dot, getOrthogonal, length, mul, normalize, sub
} from './Vector';
import { cross } from './Vector3';

export interface DistLine3Circle3Result {
    // The possible number of closest line-circle pairs is 1, 2 or all circle
    // points. If 1 or 2, numClosestPairs is set to this number and
    // 'equidistant' is false; the number of valid elements in linearClosest[]
    // and circularClosest[] is numClosestPairs. If all circle points are
    // closest, the line must be C+s*N where C is the circle center and N is
    // the normal to the plane of the circle, and linearClosest[0] is set to
    // C. In this case, 'equidistant' is true and circularClosest[0] is set to
    // C+r*U, where r is the circle radius and U is a vector perpendicular
    // to N.
    //
    // This structure is also used by the ray-circle and segment-circle
    // distance queries. For line-circle, linearClosest[] refers to the
    // closest line points. For ray-circle, linearClosest[] refers to the
    // closest ray points. For segment-circle, linearClosest[] refers to the
    // closest segment points.
    numClosestPairs: number;
    linearClosest: [Vector, Vector];
    circularClosest: [Vector, Vector];
    distance: number;
    sqrDistance: number;
    equidistant: boolean;
}

// The critical points of the squared-distance function. This is upstream's
// private 'Critical' struct, shared with the ray-circle and segment-circle
// queries through friendship.
export interface DistLine3Circle3Critical {
    numPoints: number;
    linearPoint: [Vector, Vector];
    circularPoint: [Vector, Vector];
    parameter: [number, number];
    distance: [number, number];
}

export function distLine3Circle3DefaultResult(): DistLine3Circle3Result {
    return {
        numClosestPairs: 0,
        linearClosest: [new Vector(3), new Vector(3)],
        circularClosest: [new Vector(3), new Vector(3)],
        distance: 0,
        sqrDistance: 0,
        equidistant: false
    };
}

export function distLine3Circle3DefaultCritical(): DistLine3Circle3Critical {
    return {
        numPoints: 0,
        linearPoint: [new Vector(3), new Vector(3)],
        circularPoint: [new Vector(3), new Vector(3)],
        parameter: [0, 0],
        distance: [0, 0]
    };
}

function isZero3(v: Vector): boolean {
    return v.values[0] === 0 && v.values[1] === 0 && v.values[2] === 0;
}

// Bisect the function Phi(t) = t + a0 - a1 * t / sqrt(a2 * t^2 + a3) on the
// specified interval [tauMin,tauMax]. Bisection using double precision is
// much faster than using exact rational numbers.
function bisect(a0: number, a1: number, a2: number, a3: number,
    tauMin: number, tauMax: number): number {
    const maxIterations = 4096;

    const phi = (tau: number): number => {
        return tau + a0 - a1 * tau / Math.sqrt(a2 * tau * tau + a3);
    };

    // The function is known to be increasing, so we can specify -1 and +1 as
    // the function values at the bounding interval endpoints.
    const bisector = new RootsBisection1(maxIterations);
    const bisection = bisector.find(phi, tauMin, tauMax, -1, 1);
    return bisection.root;
}

function finalize(line: Line3, circle: Circle3, D: Vector,
    result: DistLine3Circle3Result,
    critical: DistLine3Circle3Critical): void {
    for (let i = 0; i < critical.numPoints; ++i) {
        // Get the closest pair of line and circle points.
        let linearPoint = add(mul(critical.parameter[i], line.direction), D);
        const d = dot(circle.normal, linearPoint);
        const project = sub(linearPoint, mul(d, circle.normal));
        normalize(project);
        linearPoint = add(linearPoint, circle.center);
        critical.linearPoint[i] = linearPoint;
        critical.circularPoint[i] = add(circle.center,
            mul(circle.radius, project));
        const diff = sub(critical.linearPoint[i], critical.circularPoint[i]);
        critical.distance[i] = length(diff);
    }

    if (critical.numPoints === 1) {
        result.numClosestPairs = 1;
        result.distance = critical.distance[0];
        result.linearClosest[0] = critical.linearPoint[0];
        result.circularClosest[0] = critical.circularPoint[0];
    }
    else {
        // critical.numPoints is 2.
        if (critical.distance[0] < critical.distance[1]) {
            result.numClosestPairs = 1;
            result.distance = critical.distance[0];
            result.linearClosest[0] = critical.linearPoint[0];
            result.circularClosest[0] = critical.circularPoint[0];
        }
        else if (critical.distance[0] > critical.distance[1]) {
            result.numClosestPairs = 1;
            result.distance = critical.distance[1];
            result.linearClosest[0] = critical.linearPoint[1];
            result.circularClosest[0] = critical.circularPoint[1];
        }
        else {
            // critical.distance[0] equals critical.distance[1].
            result.numClosestPairs = 2;
            result.distance = critical.distance[0];
            result.linearClosest = [critical.linearPoint[0],
                critical.linearPoint[1]];
            result.circularClosest = [critical.circularPoint[0],
                critical.circularPoint[1]];
        }
    }

    result.sqrDistance = result.distance * result.distance;
}

// The line is perpendicular to the plane of the circle and contains the
// circle center.
function pdfSection411(line: Line3, circle: Circle3, D: Vector,
    result: DistLine3Circle3Result,
    critical: DistLine3Circle3Critical): void {
    const M = line.direction;
    const C = circle.center;
    const N = circle.normal;
    const r = circle.radius;

    result.numClosestPairs = 1;
    result.linearClosest[0] = C.clone();
    const U = getOrthogonal(N, true);
    result.circularClosest[0] = add(C, mul(r, U));
    const diff = sub(result.linearClosest[0], result.circularClosest[0]);
    result.sqrDistance = dot(diff, diff);
    result.distance = Math.sqrt(result.sqrDistance);
    result.equidistant = true;
    critical.numPoints = 1;
    critical.linearPoint[0] = result.linearClosest[0];
    critical.circularPoint[0] = result.circularClosest[0];
    critical.parameter[0] = -dot(M, D) / dot(M, M);
    critical.distance[0] = result.distance;
}

// The line is perpendicular to the plane of the circle and does not contain
// the circle center.
function pdfSection412(line: Line3, circle: Circle3, D: Vector,
    result: DistLine3Circle3Result,
    critical: DistLine3Circle3Critical): void {
    const M = line.direction;

    critical.numPoints = 1;
    critical.parameter[0] = -dot(M, D) / dot(M, M);
    finalize(line, circle, D, result, critical);
}

// The line is not perpendicular to the plane of the circle but the line
// origin is on the normal line through the circle center.
function pdfSection421(line: Line3, circle: Circle3, D: Vector, NxM: Vector,
    result: DistLine3Circle3Result,
    critical: DistLine3Circle3Critical): void {
    const M = line.direction;
    const r = circle.radius;

    const MdD = dot(M, D);
    const MdM = dot(M, M);
    const rLenCrossMN = r * length(NxM);
    critical.numPoints = 2;
    critical.parameter[0] = (-MdD - rLenCrossMN) / MdM;
    critical.parameter[1] = (-MdD + rLenCrossMN) / MdM;
    finalize(line, circle, D, result, critical);
}

// The line is not perpendicular to the plane of the circle and the line
// origin is not on the normal line through the circle center.
function pdfSection422(line: Line3, circle: Circle3, D: Vector, NxM: Vector,
    NxD: Vector, result: DistLine3Circle3Result,
    critical: DistLine3Circle3Critical): void {
    const M = line.direction;
    const N = circle.normal;
    const r = circle.radius;

    // Choose a new line origin E so that P(t) = E + t * M and
    // Dot(NxM, NxE) = 0.
    const NxMdNxM = dot(NxM, NxM);
    const s = -dot(NxM, NxD) / NxMdNxM;
    const E = add(mul(s, M), D);

    // Phi(t) = (t + a0) - a1 * t / (a2 * t^2 + a3)^{1/2}
    // G(t) = a1 * t / (a2 * t^2 + a3)^{1/2}
    // G'(t) = a1 * a2 / (a2 * t^2 + a3)^{3/2}
    // G"(t) = -3 * a1 * a2^2 * t / (a2 * t^2 + a3)^{5/2}
    const MdM = dot(M, M);
    const NxE = cross(N, E);
    const a0 = dot(M, E) / MdM;
    const a1 = r * NxMdNxM / MdM;   // a1 > 0
    const a2 = NxMdNxM;             // a2 > 0
    const a3 = dot(NxE, NxE);       // a3 >= 0
    let tau: number;

    if (a1 > Math.sqrt(a3)) {
        // G'(0) > 1; Math.abs guards against numerical rounding errors
        // causing the argument of Math.sqrt to be negative.
        const twoThirds = 2 / 3;
        const tauHat = Math.sqrt(
            Math.abs(Math.pow(a1 * a3, twoThirds) - a3));
        const gTauHat = a1 * tauHat / Math.sqrt(a2 * tauHat * tauHat + a3);
        const intercept = gTauHat - tauHat;  // Theoretically positive.
        if (a0 <= -intercept) {
            tau = bisect(a0, a1, a2, a3, -a0, -a0 + a1 / Math.sqrt(a2));
            if (a0 < -intercept) {
                critical.numPoints = 1;
                critical.parameter[0] = tau + s;
            }
            else {
                critical.numPoints = 2;
                critical.parameter[0] = tau + s;
                critical.parameter[1] = -tauHat + s;
            }
        }
        else if (a0 >= intercept) {
            tau = bisect(a0, a1, a2, a3, -a0 - a1 / Math.sqrt(a2), -a0);
            if (a0 > intercept) {
                critical.numPoints = 1;
                critical.parameter[0] = tau + s;
            }
            else {
                critical.numPoints = 2;
                critical.parameter[0] = tauHat + s;
                critical.parameter[1] = tau + s;
            }
        }
        else {
            critical.numPoints = 2;
            if (a0 > 0) {
                tau = bisect(a0, a1, a2, a3, -a0 - a1 / Math.sqrt(a2), -a0);
                critical.parameter[0] = tau + s;
                tau = bisect(a0, a1, a2, a3, tauHat, -a0 + a1 / Math.sqrt(a2));
                critical.parameter[1] = tau + s;
            }
            else if (a0 < 0) {
                tau = bisect(a0, a1, a2, a3, -a0 - a1 / Math.sqrt(a2),
                    -tauHat);
                critical.parameter[0] = tau + s;
                tau = bisect(a0, a1, a2, a3, -a0, -a0 + a1 / Math.sqrt(a2));
                critical.parameter[1] = tau + s;
            }
            else {
                // a0 is 0.
                tau = Math.sqrt((a1 * a1 - a3) / a2);
                critical.parameter[0] = s - tau;
                critical.parameter[1] = s + tau;
            }
        }
    }
    else {
        // G'(0) <= 1
        if (a0 < 0) {
            tau = bisect(a0, a1, a2, a3, -a0, -a0 + a1 / Math.sqrt(a2));
        }
        else if (a0 > 0) {
            tau = bisect(a0, a1, a2, a3, -a0 - a1 / Math.sqrt(a2), -a0);
        }
        else {
            tau = 0;
        }
        critical.numPoints = 1;
        critical.parameter[0] = tau + s;
    }

    finalize(line, circle, D, result, critical);
}

// The port of the friend-accessible 'Execute'. It returns both the query
// result and the critical points, the latter needed by the ray-circle and
// segment-circle queries.
export function distLine3Circle3Execute(line: Line3, circle: Circle3):
    { result: DistLine3Circle3Result, critical: DistLine3Circle3Critical } {
    const result = distLine3Circle3DefaultResult();
    const critical = distLine3Circle3DefaultCritical();

    // Translate the line and circle so that the circle center is the origin
    // (0,0,0). D is the translated line origin.
    const N = circle.normal;
    const M = line.direction;
    const D = sub(line.origin, circle.center);
    const NxM = cross(N, M);
    const NxD = cross(N, D);

    if (!isZero3(NxM)) {
        // The line is not perpendicular to the plane of the circle.
        if (!isZero3(NxD)) {
            // The line origin is not on the normal line through the circle
            // center.
            pdfSection422(line, circle, D, NxM, NxD, result, critical);
        }
        else {
            // The line origin is on the normal line through the circle
            // center.
            pdfSection421(line, circle, D, NxM, result, critical);
        }
    }
    else {
        // The line is perpendicular to the plane of the circle.
        if (!isZero3(NxD)) {
            // The line does not contain the circle center.
            pdfSection412(line, circle, D, result, critical);
        }
        else {
            // The line contains the circle center.
            pdfSection411(line, circle, D, result, critical);
        }
    }

    return { result, critical };
}

export class DistLine3Circle3
    implements DCPQuery<Line3, Circle3, DistLine3Circle3Result> {
    compute(line: Line3, circle: Circle3): DistLine3Circle3Result {
        return distLine3Circle3Execute(line, circle).result;
    }
}
