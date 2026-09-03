// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and an arc in 2D.
//
// The line is P + t * D, where P is a point on the line and D is not required
// to be unit length. The t-value is any real number.
//
// The circle containing the arc has center C and radius r. The arc has two
// endpoints E0 and E1 on the circle so that E1 is obtained from E0 by
// traversing counterclockwise. The application is responsible for ensuring
// that E0 and E1 are on the circle and that they are properly ordered.
//
// The number of pairs of closest points is result.numClosestPairs which is
// 1 or 2. If result.numClosestPairs is 1, result.parameter[0] is the line
// t-value for its closest point result.closest[0][0]. The arc closest point
// is result.closest[0][1]. If result.numClosestPairs is 2,
// result.parameter[0] and result.parameter[1] are the line t-values for its
// closest points result.closest[0][0] and result.closest[1][0]. The arc
// closest points are result.closest[0][1] and result.closest[1][1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line2<T>, Arc2<T>>' becomes the class
// DistLine2Arc2 with the result type DistLine2Arc2Result, which is
// structurally identical to the line-circle result that upstream aliases. The
// upstream single-argument 'Arc2::Contains' (which assumes the point is on
// the circle) is the port's 'containsOnCircle'.

import { Arc2 } from './Arc2.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistLine2Circle2 } from './DistLine2Circle2.js';
import { DistPointLine } from './DistPointLine.js';
import { Hypersphere } from './Hypersphere.js';
import type { Line2 } from './Line.js';
import { Vector } from './Vector.js';

export interface DistLine2Arc2Result {
    distance: number;
    sqrDistance: number;

    // The number of pairs of closest points, 1 or 2.
    numClosestPairs: number;

    // parameter[j] is the line t-value of closest[j][0].
    parameter: [number, number];

    // closest[j][0] is on the line, closest[j][1] is on the arc.
    closest: [[Vector, Vector], [Vector, Vector]];
}

export class DistLine2Arc2
    implements DCPQuery<Line2, Arc2, DistLine2Arc2Result> {
    compute(line: Line2, arc: Arc2): DistLine2Arc2Result {
        const result: DistLine2Arc2Result = {
            distance: 0,
            sqrDistance: 0,
            numClosestPairs: 0,
            parameter: [0, 0],
            closest: [
                [new Vector(2), new Vector(2)],
                [new Vector(2), new Vector(2)]
            ]
        };

        // Execute the query for line-circle. Test whether the circle closest
        // points are on or off the arc. If any closest point is on the arc,
        // there is no need to test arc endpoints for closeness.
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const lcResult = new DistLine2Circle2().compute(line, circle);
        for (let i = 0; i < lcResult.numClosestPairs; ++i) {
            if (arc.containsOnCircle(lcResult.closest[i][1])) {
                const j = result.numClosestPairs++;
                result.distance = lcResult.distance;
                result.sqrDistance = lcResult.sqrDistance;
                result.parameter[j] = lcResult.parameter[i];
                result.closest[j][0] = lcResult.closest[i][0];
                result.closest[j][1] = lcResult.closest[i][1];
            }
        }

        if (result.numClosestPairs > 0) {
            // At least one circle closest point is on the arc. There is no
            // need to test arc endpoints.
            return result;
        }

        // No circle closest points are on the arc. Compute distances to the
        // arc endpoints and select the minima.
        const plQuery = new DistPointLine();
        const plResult0 = plQuery.compute(arc.end[0], line);
        const plResult1 = plQuery.compute(arc.end[1], line);
        if (plResult0.sqrDistance < plResult1.sqrDistance) {
            result.distance = Math.sqrt(plResult0.sqrDistance);
            result.sqrDistance = plResult0.sqrDistance;
            result.numClosestPairs = 1;
            result.parameter[0] = plResult0.parameter;
            result.closest[0][0] = plResult0.closest[1];
            result.closest[0][1] = arc.end[0].clone();
        }
        else if (plResult1.sqrDistance < plResult0.sqrDistance) {
            result.distance = Math.sqrt(plResult1.sqrDistance);
            result.sqrDistance = plResult1.sqrDistance;
            result.numClosestPairs = 1;
            result.parameter[0] = plResult1.parameter;
            result.closest[0][0] = plResult1.closest[1];
            result.closest[0][1] = arc.end[1].clone();
        }
        else {
            result.distance = Math.sqrt(plResult0.sqrDistance);
            result.sqrDistance = plResult0.sqrDistance;
            result.numClosestPairs = 2;
            result.parameter[0] = plResult0.parameter;
            result.parameter[1] = plResult1.parameter;
            result.closest[0][0] = plResult0.closest[1];
            result.closest[0][1] = arc.end[0].clone();
            result.closest[1][0] = plResult1.closest[1];
            result.closest[1][1] = arc.end[1].clone();
        }

        return result;
    }
}
