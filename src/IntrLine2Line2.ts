// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Line2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection and find-intersection queries for two lines. The line
// directions are nonzero but not required to be unit length.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent.

import { Line } from './Line';
import { Vector, add, sub, mul } from './Vector';
import { dotPerp } from './Vector2';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';

// The port of std::numeric_limits<int32_t>::max(), the 'numIntersections'
// value meaning "the lines are the same".
const INT32_MAX = 2147483647;

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrLine2Line2TI.test.
//
// If the lines do not intersect,
//   intersect = false
//   numIntersections = 0
//
// If the lines intersect in a single point,
//   intersect = true
//   numIntersections = 1
//
// If the lines are the same,
//   intersect = true
//   numIntersections = 2147483647
export interface IntrLine2Line2TIResult {
    intersect: boolean;
    numIntersections: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine2Line2TIResult {
    return { intersect: false, numIntersections: 0 };
}

// The result of IntrLine2Line2FI.find.
//
// If the lines do not intersect,
//   intersect = false
//   numIntersections = 0
//   line0Parameter = [0, 0]  // invalid
//   line1Parameter = [0, 0]  // invalid
//   point = (0, 0)  // invalid
//
// If the lines intersect in a single point, the parameter for line0 is s0 and
// the parameter for line1 is s1,
//   intersect = true
//   numIntersections = 1
//   line0Parameter = [s0, s0]
//   line1Parameter = [s1, s1]
//   point = line0.origin + s0 * line0.direction
//         = line1.origin + s1 * line1.direction
//
// If the lines are the same, let maxT = Number.MAX_VALUE,
//   intersect = true
//   numIntersections = 2147483647
//   line0Parameter = [-maxT, +maxT]
//   line1Parameter = [-maxT, +maxT]
//   point = (0, 0)  // invalid
export interface IntrLine2Line2FIResult {
    intersect: boolean;
    numIntersections: number;
    line0Parameter: [number, number];
    line1Parameter: [number, number];
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine2Line2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        line0Parameter: [0, 0],
        line1Parameter: [0, 0],
        point: Vector.zero(2)
    };
}

export class IntrLine2Line2TI implements
    TIQuery<Line, Line, IntrLine2Line2TIResult> {

    test(line0: Line, line1: Line): IntrLine2Line2TIResult {
        const result = defaultTIResult();

        // The intersection of two lines is a solution to P0 + s0 * D0 =
        // P1 + s1 * D1. Rewrite this as s0*D0 - s1*D1 = P1 - P0 = Q. If
        // DotPerp(D0, D1) = 0, the lines are parallel. Additionally, if
        // DotPerp(Q, D1) = 0, the lines are the same. If DotPerp(D0, D1) is
        // not zero, then the lines intersect in a single point where
        //   s0 = DotPerp(Q, D1)/DotPerp(D0, D1)
        //   s1 = DotPerp(Q, D0)/DotPerp(D0, D1)

        const dotD0PerpD1 = dotPerp(line0.direction, line1.direction);
        if (dotD0PerpD1 !== 0) {
            // The lines are not parallel.
            result.intersect = true;
            result.numIntersections = 1;
        } else {
            // The lines are parallel.
            const Q = sub(line1.origin, line0.origin);
            const dotQDotPerpD1 = dotPerp(Q, line1.direction);
            if (dotQDotPerpD1 !== 0) {
                // The lines are parallel but distinct.
                result.intersect = false;
                result.numIntersections = 0;
            } else {
                // The lines are the same.
                result.intersect = true;
                result.numIntersections = INT32_MAX;
            }
        }

        return result;
    }
}

export class IntrLine2Line2FI implements
    FIQuery<Line, Line, IntrLine2Line2FIResult> {

    find(line0: Line, line1: Line): IntrLine2Line2FIResult {
        const result = defaultFIResult();

        // See the comments in IntrLine2Line2TI.test for the derivation.
        const Q = sub(line1.origin, line0.origin);
        const dotD0PerpD1 = dotPerp(line0.direction, line1.direction);
        if (dotD0PerpD1 !== 0) {
            // The lines are not parallel.
            result.intersect = true;
            result.numIntersections = 1;
            const dotQPerpD0 = dotPerp(Q, line0.direction);
            const dotQPerpD1 = dotPerp(Q, line1.direction);
            const s0 = dotQPerpD1 / dotD0PerpD1;
            const s1 = dotQPerpD0 / dotD0PerpD1;
            result.line0Parameter = [s0, s0];
            result.line1Parameter = [s1, s1];
            result.point = add(line0.origin, mul(s0, line0.direction));
        } else {
            // The lines are parallel.
            const dotQPerpD1 = dotPerp(Q, line1.direction);
            if (Math.abs(dotQPerpD1) !== 0) {
                // The lines are parallel but distinct.
                result.intersect = false;
                result.numIntersections = 0;
            } else {
                // The lines are the same.
                result.intersect = true;
                result.numIntersections = INT32_MAX;
                result.line0Parameter = [-MAX_T, MAX_T];
                result.line1Parameter = [-MAX_T, MAX_T];
            }
        }

        return result;
    }
}
