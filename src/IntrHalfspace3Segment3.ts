// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Segment3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent.

import { Halfspace } from './Halfspace.js';
import { Segment } from './Segment.js';
import { Vector, add, sub, mul, dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The result of IntrHalfspace3Segment3TI.test.
export interface IntrHalfspace3Segment3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Segment3TIResult {
    return { intersect: false };
}

// The result of IntrHalfspace3Segment3FI.find.
export interface IntrHalfspace3Segment3FIResult {
    intersect: boolean;

    // The segment is clipped against the plane defining the halfspace. The
    // 'numPoints' is either 0 (no intersection), 1 (point), or 2 (segment).
    numPoints: number;
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrHalfspace3Segment3FIResult {
    return {
        intersect: false,
        numPoints: 0,
        point: [Vector.zero(3), Vector.zero(3)]
    };
}

export class IntrHalfspace3Segment3TI implements
    TIQuery<Halfspace, Segment, IntrHalfspace3Segment3TIResult> {

    test(halfspace: Halfspace, segment: Segment): IntrHalfspace3Segment3TIResult {
        const result = defaultTIResult();

        // Project the segment endpoints onto the normal line. The plane of
        // the halfspace occurs at the origin (zero) of the normal line.
        const s: number[] = [0, 0];
        for (let i = 0; i < 2; ++i) {
            s[i] = dot(halfspace.normal, segment.p[i]) - halfspace.constant;
        }

        // The segment and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (Math.max(s[0], s[1]) >= 0);
        return result;
    }
}

export class IntrHalfspace3Segment3FI implements
    FIQuery<Halfspace, Segment, IntrHalfspace3Segment3FIResult> {

    find(halfspace: Halfspace, segment: Segment): IntrHalfspace3Segment3FIResult {
        // Determine on which side of the plane the endpoints lie. The table
        // of possibilities is listed next with n = numNegative,
        // p = numPositive, and z = numZero.
        //
        //   n p z  intersection
        //   -------------------------
        //   0 2 0  segment (original)
        //   0 1 1  segment (original)
        //   0 0 2  segment (original)
        //   1 1 0  segment (clipped)
        //   1 0 1  point (endpoint)
        //   2 0 0  none
        //
        // UPSTREAM BUG (preserved): in the (n,p,z) = (1,1,0) case, the
        // intersection is the clipped segment from the positive endpoint to
        // the point where the segment crosses the plane, but upstream reports
        // numPoints = 1 and stores only the crossing point. The port keeps
        // this behavior so results match the C++ exactly; see the PR notes.

        const s: number[] = [0, 0];
        let numPositive = 0, numNegative = 0;
        for (let i = 0; i < 2; ++i) {
            s[i] = dot(halfspace.normal, segment.p[i]) - halfspace.constant;
            if (s[i] > 0) {
                ++numPositive;
            } else if (s[i] < 0) {
                ++numNegative;
            }
        }

        const result = defaultFIResult();

        if (numNegative === 0) {
            // The segment is in the halfspace.
            result.intersect = true;
            result.numPoints = 2;
            result.point[0] = segment.p[0].clone();
            result.point[1] = segment.p[1].clone();
        } else if (numNegative === 1) {
            result.intersect = true;
            result.numPoints = 1;
            if (numPositive === 1) {
                // The segment is intersected at an interior point.
                result.point[0] = add(segment.p[0],
                    mul(s[0] / (s[0] - s[1]), sub(segment.p[1], segment.p[0])));
            } else {
                // numZero = 1. One segment endpoint is on the plane.
                if (s[0] === 0) {
                    result.point[0] = segment.p[0].clone();
                } else {
                    result.point[0] = segment.p[1].clone();
                }
            }
        } else {
            // The segment is outside the halfspace (numNegative == 2).
            result.intersect = false;
            result.numPoints = 0;
        }

        return result;
    }
}
