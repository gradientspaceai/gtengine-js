// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SeparatePoints2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Separate two point sets, if possible, by computing a line for which the
// point sets lie on opposite sides. The algorithm computes the convex hull
// of the point sets, then uses the method of separating axes to determine
// whether the two convex polygons are disjoint.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
//
// Port notes:
// * Upstream templates the class on <Real, ComputeType>, but ComputeType is
//   never used: the current ConvexHull2 is templated only on the
//   floating-point type and performs its own exact predicates internally.
//   The port drops the vestigial parameter.
// * operator() returns a bool and writes the separating line through an
//   output reference. The port's compute(...) returns
//   { separated, separatingLine }. Upstream leaves 'separatingLine' holding
//   the last candidate line tested when the return value is false; the port
//   returns a line in the same state, but the field is meaningful only when
//   'separated' is true.
// * The code assumes each point set has at least 3 noncollinear points. When
//   either hull is 0- or 1-dimensional, the query reports "not separated",
//   as upstream does.
//
// Upstream bug (SeparatePoints2.h, OnSameSide and WhichSide): the side
// classification compares Dot(lineNormal,P) to lineConstant using exact
// floating-point comparisons, where lineNormal = Perp(Normalize(P1 - P0))
// and lineConstant = Dot(lineNormal,P0) are both rounded. The points of the
// hull that owns the candidate edge should classify as "on the line" or
// "negative side" (the normal of a counterclockwise hull edge points
// outward), but round-off makes the far endpoint of the edge classify as
// strictly positive. WhichSide then returns +1 instead of -1 for that hull
// and the test 'side0 * side1 <= 0' succeeds, so the query reports a
// separation for point sets that plainly overlap. For example, the regular
// hexagons with center (0,0), radius 1, phase 0.3 and center (1,0.5),
// radius 1.5, phase 1.1 overlap heavily, yet upstream reports them
// separated. The port fixes this by evaluating the side of a point relative
// to the edge <P0,P1> with an exact orientation predicate on the
// unnormalized edge normal Perp(P1 - P0), which has the same sign as the
// upstream expression in exact arithmetic. This follows the ConvexHull2
// (B100) precedent of computing the geometric predicates exactly. The
// separating line that is returned is still the normalized floating-point
// line of upstream.

import { BSNumber } from './BSNumber.js';
import { ConvexHull2 } from './ConvexHull2.js';
import { Line, type Line2 } from './Line.js';
import { Vector, normalize, sub } from './Vector.js';

// The result of the separation query. 'separated' is upstream's boolean
// return value: true if and only if 'separatingLine' separates the two point
// sets.
export interface SeparatePoints2Result {
    separated: boolean;
    separatingLine: Line2;
}

// The exact port of Vector2<Real> coordinates for the orientation predicate.
type ExactPoint2 = [BSNumber, BSNumber];

// The sign of Dot(Perp(p1 - p0), q - p0), which is negative when q is on the
// inner side of the counterclockwise hull edge <p0,p1>, positive when q is
// on the outer side, and zero when q is on the edge line.
function sideSign(p0: ExactPoint2, p1: ExactPoint2, q: ExactPoint2): number {
    const dx = p1[0].sub(p0[0]);
    const dy = p1[1].sub(p0[1]);
    return dy.mul(q[0].sub(p0[0])).sub(dx.mul(q[1].sub(p0[1]))).getSign();
}

// Test whether all hull points are on the same side of the edge line. The
// return value is 0 when the line splits the point set, +1 when the points
// are on the positive side and -1 otherwise. A hull whose points all lie
// exactly on the line reports -1, matching upstream.
function onSameSide(p0: ExactPoint2, p1: ExactPoint2,
    hull: readonly ExactPoint2[]): number {
    let posSide = 0;
    let negSide = 0;

    for (const q of hull) {
        const s = sideSign(p0, p1, q);
        if (s > 0) {
            ++posSide;
        }
        else if (s < 0) {
            ++negSide;
        }

        if (posSide !== 0 && negSide !== 0) {
            // The line splits the point set.
            return 0;
        }
    }

    return posSide !== 0 ? +1 : -1;
}

// Establish which side of the edge line the hull is on. The return value is
// 0 when the hull is effectively collinear with the line.
function whichSide(p0: ExactPoint2, p1: ExactPoint2,
    hull: readonly ExactPoint2[]): number {
    for (const q of hull) {
        const s = sideSign(p0, p1, q);
        if (s > 0) {
            // The hull is on the positive side.
            return +1;
        }
        if (s < 0) {
            // The hull is on the negative side.
            return -1;
        }
    }

    // The hull is effectively collinear.
    return 0;
}

// Convert the hull vertices of a point set to exact coordinates.
function exactHull(points: readonly Vector[],
    edges: readonly number[]): ExactPoint2[] {
    return edges.map(i => [
        BSNumber.fromNumber(points[i].get(0)),
        BSNumber.fromNumber(points[i].get(1))
    ] as ExactPoint2);
}

export class SeparatePoints2 {
    // The return value 'separated' is true if and only if there is a
    // separation. If true, the returned line is a separating line. The code
    // assumes that each point set has at least 3 noncollinear points.
    compute(points0: readonly Vector[], points1: readonly Vector[]):
        SeparatePoints2Result {
        const separatingLine = new Line(2);

        // Construct the convex hull of point set 0.
        const ch0 = new ConvexHull2();
        ch0.compute(points0);
        if (ch0.getDimension() !== 2) {
            return { separated: false, separatingLine };
        }

        // Construct the convex hull of point set 1.
        const ch1 = new ConvexHull2();
        ch1.compute(points1);
        if (ch1.getDimension() !== 2) {
            return { separated: false, separatingLine };
        }

        const edges0 = ch0.getHull();
        const edges1 = ch1.getHull();
        const hull0 = exactHull(points0, edges0);
        const hull1 = exactHull(points1, edges1);

        // Test the edges of hull 0 for possible separation of the points.
        const numEdges0 = edges0.length;
        for (let j1 = 0, j0 = numEdges0 - 1; j1 < numEdges0; j0 = j1++) {
            // Look up the edge (assert: i0 != i1).
            const i0 = edges0[j0];
            const i1 = edges0[j1];

            // Determine whether hull 1 is on the same side of the line
            // (assert: Perp(P1 - P0) != (0,0)).
            const side1 = onSameSide(hull0[j0], hull0[j1], hull1);

            if (side1 !== 0) {
                // Determine on which side of the line hull 0 lies.
                const side0 = whichSide(hull0[j0], hull0[j1], hull0);

                if (side0 * side1 <= 0) {
                    // The line separates the hulls.
                    setLine(separatingLine, points0[i0], points0[i1]);
                    return { separated: true, separatingLine };
                }
            }
        }

        // Test the edges of hull 1 for possible separation of the points.
        const numEdges1 = edges1.length;
        for (let j1 = 0, j0 = numEdges1 - 1; j1 < numEdges1; j0 = j1++) {
            // Look up the edge (assert: i0 != i1).
            const i0 = edges1[j0];
            const i1 = edges1[j1];

            // Determine whether hull 0 is on the same side of the line
            // (assert: Perp(P1 - P0) != (0,0)).
            const side0 = onSameSide(hull1[j0], hull1[j1], hull0);

            if (side0 !== 0) {
                // Determine on which side of the line hull 1 lies.
                const side1 = whichSide(hull1[j0], hull1[j1], hull1);

                if (side0 * side1 <= 0) {
                    // The line separates the hulls.
                    setLine(separatingLine, points1[i0], points1[i1]);
                    return { separated: true, separatingLine };
                }
            }
        }

        return { separated: false, separatingLine };
    }
}

// The separating line has the edge's first endpoint as origin and the
// unit-length edge direction, as upstream.
function setLine(line: Line2, p0: Vector, p1: Vector): void {
    line.origin = p0.clone();
    line.direction = sub(p1, p0);
    normalize(line.direction);
}
