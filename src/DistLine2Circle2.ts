// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a circle in 2D. The circle is
// considered to be a curve, not a solid disk.
//
// The line is P + t * D, where P is a point on the line and D is not required
// to be unit length. The t-value is any real number.
//
// The circle is C + r * U(s), where C is the center, r > 0 is the radius, and
// U(s) = (cos(s), sin(s)) for s in [0,2*pi).
//
// The number of pairs of closest points is result.numClosestPairs, which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the line
// t-value for its closest point result.closest[0][0]. The circle closest
// point is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the line t-values for its
// closest points result.closest[0][0] and result.closest[1][0]. The circle
// closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line2<T>, Circle2<T>>' becomes the
// class DistLine2Circle2 with the result type DistLine2Circle2Result. The
// protected static helper 'DoQuery' becomes the module-private function
// doQuery.

import type { DCPQuery } from './DCPQuery';
import type { Circle2 } from './Hypersphere';
import type { Line2 } from './Line';
import { Vector, add, dot, mul, normalize, sub } from './Vector';
import { dotPerp } from './Vector2';

export interface DistLine2Circle2Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, 1 or 2.
    numClosestPairs: number;

    // parameter[j] is the line t-value of closest[j][0].
    parameter: [number, number];

    // closest[j][0] is on the line, closest[j][1] is on the circle.
    closest: [[Vector, Vector], [Vector, Vector]];
}

function defaultResult(): DistLine2Circle2Result {
    return {
        distance: 0,
        sqrDistance: 0,
        numClosestPairs: 0,
        parameter: [0, 0],
        closest: [
            [new Vector(2), new Vector(2)],
            [new Vector(2), new Vector(2)]
        ]
    };
}

// Compute the distance and closest points between a line and a circle whose
// center is the origin.
function doQuery(delta: Vector, direction: Vector, radius: number,
    result: DistLine2Circle2Result): void {
    // Compute the distance from the line to the origin. The distance is
    // d = |Dot(D,Perp(D))|/|D|. The line does not intersect the circle when
    // d > r. The line is tangent to the circle when d = r. The line
    // intersects the circle in 2 points when d < r. Rather than normalize D
    // at this time, replace the comparisons by sign tests for
    // |Dot(D,Perp(D))|^2 - r^2 * |D|^2. This allows for theoretically correct
    // classification of line-circle tangency when using rational arithmetic.
    const dotDirDir = dot(direction, direction);
    const dotDirDel = dot(direction, delta);
    const dotPerpDirDel = dotPerp(direction, delta);
    const rSqr = radius * radius;
    const test = dotPerpDirDel * dotPerpDirDel - rSqr * dotDirDir;
    if (test >= 0) {
        // When the line-origin distance equals the radius, the line is
        // tangent to the circle; there is 1 point of intersection and the
        // line-circle distance is 0. When the line-origin distance is larger
        // than the radius, the line and circle do not intersect. The closest
        // circle point is the tangent point if the line were to be translated
        // in its normal direction to just touch the circle. In this case, the
        // distance between the circle and line is the difference between the
        // line-origin distance and the radius.

        // Compute the line point closest to the circle.
        result.numClosestPairs = 1;
        result.parameter[0] = -dotDirDel / dotDirDir;
        result.closest[0][0] = add(delta, mul(result.parameter[0], direction));
        result.closest[0][1] = result.closest[0][0].clone();

        // Compute the circle point closest to the line.
        if (test > 0) {
            normalize(result.closest[0][1]);
            result.closest[0][1] = mul(result.closest[0][1], radius);
        }
    }
    else {  // lineOriginDistance < radius
        // The line and circle intersect in 2 points. Solve the quadratic
        // equation a2*t^2 + 2*a1*t + a0 = 0. The solutions are
        // (-a1 +/- sqrt(a1 * a1 - a0 * a2)) / a2. Theoretically, discr > 0.
        // Guard against a negative floating-point result.
        const a0 = dot(delta, delta) - radius * radius;
        const a1 = dotDirDel;
        const a2 = dotDirDir;
        const discr = Math.max(a1 * a1 - a0 * a2, 0);
        const sqrtDiscr = Math.sqrt(discr);

        // Evaluate the line parameters but do so to avoid subtractive
        // cancellation.
        const temp = -dotDirDel + (dotDirDel > 0 ? -sqrtDiscr : sqrtDiscr);
        result.numClosestPairs = 2;
        result.parameter[0] = temp / dotDirDir;
        result.parameter[1] = a0 / temp;
        if (result.parameter[0] > result.parameter[1]) {
            const save = result.parameter[0];
            result.parameter[0] = result.parameter[1];
            result.parameter[1] = save;
        }

        // Compute the intersection points.
        result.closest[0][0] = add(delta, mul(result.parameter[0], direction));
        result.closest[0][1] = result.closest[0][0].clone();
        result.closest[1][0] = add(delta, mul(result.parameter[1], direction));
        result.closest[1][1] = result.closest[1][0].clone();
    }
}

export class DistLine2Circle2
    implements DCPQuery<Line2, Circle2, DistLine2Circle2Result> {
    compute(line: Line2, circle: Circle2): DistLine2Circle2Result {
        const result = defaultResult();

        // Translate the line and circle so that the circle has center at the
        // origin.
        const delta = sub(line.origin, circle.center);

        // The query computes 'result' relative to the circle with center at
        // the origin.
        doQuery(delta, line.direction, circle.radius, result);

        // Translate the closest points to the original coordinates and then
        // compute the distance and squared distance.
        for (let j = 0; j < result.numClosestPairs; ++j) {
            for (let i = 0; i < 2; ++i) {
                result.closest[j][i] =
                    add(result.closest[j][i], circle.center);
            }
        }

        const diff = sub(result.closest[0][0], result.closest[0][1]);
        result.sqrDistance = dot(diff, diff);
        result.distance = Math.sqrt(result.sqrDistance);
        return result;
    }
}
