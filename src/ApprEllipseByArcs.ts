// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ApprEllipseByArcs.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The ellipse is (x/a)^2 + (y/b)^2 = 1, but only the portion in the first
// quadrant (x >= 0 and y >= 0) is approximated. Generate numArcs >= 2 arcs by
// constructing points corresponding to the weighted averages of the
// curvatures at the ellipse points (a,0) and (0,b). The returned 'points'
// array has numArcs+1 elements and the returned 'centers' and 'radii' arrays
// each have numArcs elements. The arc associated with points[i] and
// points[i+1] has center centers[i] and radius radii[i]. The algorithm is
// described in
//   https://www.geometrictools.com/Documentation/ApproximateEllipse.pdf
//
// Port notes: upstream returns a bool and writes three std::vector output
// references that it clears on failure. The port returns the three containers
// in an object literal, or null on failure, per PORTING.md. Upstream's
// 'Circumscribe' becomes 'circumscribeCircle2' (the Cont* naming precedent),
// which returns null instead of writing an output reference.

import { circumscribeCircle2 } from './ContScribeCircle2';
import { Vector } from './Vector';

export interface ApproximateEllipseByArcsResult {
    // The numArcs+1 points on the ellipse arc in the first quadrant, ordered
    // counterclockwise from (a,0) to (0,b).
    points: Vector[];

    // The numArcs arc centers; centers[i] belongs to the arc with endpoints
    // points[i] and points[i+1].
    centers: Vector[];

    // The numArcs arc radii; radii[i] belongs to the arc with endpoints
    // points[i] and points[i+1].
    radii: number[];
}

// The function returns the arc data when the approximation succeeded, in
// which case the output arrays are nonempty. If 'numArcs' is smaller than 2
// or a == b or one of the calls to circumscribeCircle2 fails, the function
// returns null.
export function approximateEllipseByArcs(a: number, b: number,
    numArcs: number): ApproximateEllipseByArcsResult | null {
    if (numArcs < 2 || a === b) {
        // At least 2 arcs are required. The ellipse cannot already be a
        // circle.
        return null;
    }

    const points = new Array<Vector>(numArcs + 1);
    const centers = new Array<Vector>(numArcs);
    const radii = new Array<number>(numArcs);

    // Compute intermediate ellipse quantities.
    const a2 = a * a, b2 = b * b, ab = a * b;
    const invB2mA2 = 1 / (b2 - a2);

    // Compute the endpoints of the ellipse in the first quadrant. The points
    // are generated in counterclockwise order.
    points[0] = Vector.fromArray([a, 0]);
    points[numArcs] = Vector.fromArray([0, b]);

    // Compute the curvature at the endpoints. These are used when computing
    // the arcs.
    const curv0 = a / b2;
    const curv1 = b / a2;

    // Select the ellipse points based on curvature properties.
    const invNumArcs = 1 / numArcs;
    for (let i = 1; i < numArcs; ++i) {
        // The curvature at a new point is a weighted average of curvature at
        // the endpoints.
        const weight1 = i * invNumArcs;
        const weight0 = 1 - weight1;
        const curv = weight0 * curv0 + weight1 * curv1;

        // Compute the point having this curvature.
        const tmp = Math.pow(ab / curv, 2 / 3);
        points[i] = Vector.fromArray([
            a * Math.sqrt(Math.abs((tmp - a2) * invB2mA2)),
            b * Math.sqrt(Math.abs((tmp - b2) * invB2mA2))
        ]);
    }

    // Compute the arc at (a,0). The reflection of points[1] about the x-axis
    // is used so that the arc is tangent to the ellipse at (a,0).
    const p0 = points[0];
    const p1 = points[1];
    let circle = circumscribeCircle2(
        Vector.fromArray([p1.values[0], -p1.values[1]]), p0, p1);
    if (circle === null) {
        // This should not happen for the arc-fitting algorithm.
        return null;
    }
    centers[0] = circle.center;
    radii[0] = circle.radius;

    // Compute the arc at (0,b). The reflection of points[numArcs-1] about the
    // y-axis is used so that the arc is tangent to the ellipse at (0,b).
    const last = numArcs - 1;
    const pNm1 = points[last];
    const pN = points[numArcs];
    circle = circumscribeCircle2(
        Vector.fromArray([-pNm1.values[0], pNm1.values[1]]), pN, pNm1);
    if (circle === null) {
        // This should not happen for the arc-fitting algorithm.
        return null;
    }
    centers[last] = circle.center;
    radii[last] = circle.radius;

    // Compute the arcs at intermediate points between (a,0) and (0,b).
    for (let iM = 0, i = 1, iP = 2; i < last; ++iM, ++i, ++iP) {
        // Upstream ignores the boolean returned by this Circumscribe call,
        // which leaves 'circle' holding the previous arc's data when the
        // point triple is collinear. The port reports the failure instead of
        // silently storing a stale circle.
        circle = circumscribeCircle2(points[iM], points[i], points[iP]);
        if (circle === null) {
            return null;
        }
        centers[i] = circle.center;
        radii[i] = circle.radius;
    }
    return { points, centers, radii };
}
