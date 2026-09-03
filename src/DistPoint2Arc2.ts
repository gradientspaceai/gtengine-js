// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and an arc in 2D.
//
// The input point is stored in the member closest[0]. If a single point on
// the arc is closest to the input point, the member closest[1] is set to that
// point and the equidistant member is set to false. If the entire arc is
// equidistant to the point, the member closest[1] is set to the endpoint E0
// of the arc and the equidistant member is set to true.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector2<T>, Arc2<T>>' becomes the
// class DistPoint2Arc2 with the result type DistPoint2Arc2Result. The
// upstream single-argument 'Arc2::Contains' (which assumes the point is on
// the circle) is the port's 'containsOnCircle'.

import { Arc2 } from './Arc2.js';
import type { DCPQuery } from './DCPQuery.js';
import { DistPoint2Circle2 } from './DistPoint2Circle2.js';
import { Hypersphere } from './Hypersphere.js';
import { Vector, dot, sub } from './Vector.js';

export interface DistPoint2Arc2Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest arc point.
    closest: [Vector, Vector];

    // True when the input point is the center of the circle containing the
    // arc, in which case every arc point is equidistant from it.
    equidistant: boolean;
}

export class DistPoint2Arc2
    implements DCPQuery<Vector, Arc2, DistPoint2Arc2Result> {
    compute(point: Vector, arc: Arc2): DistPoint2Arc2Result {
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const pcResult = new DistPoint2Circle2().compute(point, circle);
        if (!pcResult.equidistant) {
            // Test whether the closest circle point is on the arc. If it is,
            // that point is the closest arc point. If it is not, the closest
            // arc point is an arc endpoint. Determine which endpoint that is.
            if (arc.containsOnCircle(pcResult.closest[1])) {
                return {
                    distance: pcResult.distance,
                    sqrDistance: pcResult.sqrDistance,
                    closest: [pcResult.closest[0], pcResult.closest[1]],
                    equidistant: pcResult.equidistant
                };
            }

            const diff0 = sub(arc.end[0], point);
            const diff1 = sub(arc.end[1], point);
            const sqrLength0 = dot(diff0, diff0);
            const sqrLength1 = dot(diff1, diff1);
            if (sqrLength0 <= sqrLength1) {
                return {
                    distance: Math.sqrt(sqrLength0),
                    sqrDistance: sqrLength0,
                    closest: [point.clone(), arc.end[0].clone()],
                    equidistant: false
                };
            }
            return {
                distance: Math.sqrt(sqrLength1),
                sqrDistance: sqrLength1,
                closest: [point.clone(), arc.end[1].clone()],
                equidistant: false
            };
        }

        // The point is the center of the circle containing the arc.
        return {
            distance: arc.radius,
            sqrDistance: arc.radius * arc.radius,
            closest: [point.clone(), arc.end[0].clone()],
            equidistant: true
        };
    }
}
