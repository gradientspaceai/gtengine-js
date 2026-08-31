// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLineLine.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two lines in nD.
//
// The lines are P[i] + s[i] * D[i], where D[i] is not required to be unit
// length.
//
// The closest point on line[i] is stored in closest[i] with parameter[i]
// storing s[i]. When there are infinitely many choices for the pair of
// closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line<N,T>, Line<N,T>>' becomes the
// class DistLineLine with the result type DistLineLineResult.

import type { DCPQuery } from './DCPQuery';
import type { Line } from './Line';
import { Vector, add, dot, mul, sub } from './Vector';

export interface DistLineLineResult {
    distance: number;
    sqrDistance: number;

    // parameter[i] is the s[i] of the closest point on line[i].
    parameter: [number, number];

    // closest[i] is the closest point on line[i].
    closest: [Vector, Vector];
}

export class DistLineLine implements DCPQuery<Line, Line, DistLineLineResult> {
    compute(line0: Line, line1: Line): DistLineLineResult {
        let diff = sub(line0.origin, line1.origin);
        const a00 = dot(line0.direction, line0.direction);
        const a01 = -dot(line0.direction, line1.direction);
        const a11 = dot(line1.direction, line1.direction);
        const b0 = dot(line0.direction, diff);
        const det = Math.max(a00 * a11 - a01 * a01, 0);
        let s0: number;
        let s1: number;

        if (det > 0) {
            // The lines are not parallel.
            const b1 = -dot(line1.direction, diff);
            s0 = (a01 * b1 - a11 * b0) / det;
            s1 = (a01 * b0 - a00 * b1) / det;
        }
        else {
            // The lines are parallel. Select any pair of closest points.
            s0 = -b0 / a00;
            s1 = 0;
        }

        const closest0 = add(line0.origin, mul(s0, line0.direction));
        const closest1 = add(line1.origin, mul(s1, line1.direction));
        diff = sub(closest0, closest1);
        const sqrDistance = dot(diff, diff);

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            parameter: [s0, s1],
            closest: [closest0, closest1]
        };
    }
}
