// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrCircle2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Find-intersection query for a circle and an arc.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// an FIQuery specialization, so the port has only IntrCircle2Arc2FI. The
// single-argument upstream Arc2::Contains (which assumes the point is on the
// circle of the arc) is Arc2.containsOnCircle in the port.

import { Arc2 } from './Arc2';
import { Hypersphere } from './Hypersphere';
import { Vector } from './Vector';
import { IntrCircle2Circle2FI } from './IntrCircle2Circle2';
import type { FIQuery } from './FIQuery';

// The port of std::numeric_limits<int32_t>::max(), the 'numIntersections'
// value meaning "the arc is on the circle".
const INT32_MAX = 2147483647;

// The result of IntrCircle2Arc2FI.find.
export interface IntrCircle2Arc2FIResult {
    intersect: boolean;

    // The number of intersections is 0, 1, 2 or 2147483647 (the port of
    // std::numeric_limits<int32_t>::max()). When 1, the arc and circle
    // intersect in a single point. When 2, the arc is not on the circle and
    // they intersect in two points. When 2147483647, the arc is on the circle.
    numIntersections: number;

    // Valid only when numIntersections is 1 or 2.
    point: [Vector, Vector];

    // Valid only when numIntersections is 2147483647.
    arc: Arc2;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrCircle2Arc2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        point: [Vector.zero(2), Vector.zero(2)],
        arc: Arc2.fromCenterRadiusEnds(Vector.zero(2), 0, Vector.zero(2),
            Vector.zero(2))
    };
}

export class IntrCircle2Arc2FI implements
    FIQuery<Hypersphere, Arc2, IntrCircle2Arc2FIResult> {

    find(circle: Hypersphere, arc: Arc2): IntrCircle2Arc2FIResult {
        const result = defaultFIResult();

        const circleOfArc = Hypersphere.fromCenterRadius(arc.center,
            arc.radius);
        const ccQuery = new IntrCircle2Circle2FI();
        const ccResult = ccQuery.find(circle, circleOfArc);
        if (!ccResult.intersect) {
            result.intersect = false;
            result.numIntersections = 0;
            return result;
        }

        if (ccResult.numIntersections === INT32_MAX) {
            // The arc is on the circle.
            result.intersect = true;
            result.numIntersections = INT32_MAX;
            result.arc = arc.clone();
            return result;
        }

        // Test whether circle-circle intersection points are on the arc.
        result.numIntersections = 0;
        for (let i = 0; i < ccResult.numIntersections; ++i) {
            if (arc.containsOnCircle(ccResult.point[i])) {
                result.point[result.numIntersections++] = ccResult.point[i];
                result.intersect = true;
            }
        }

        return result;
    }
}
