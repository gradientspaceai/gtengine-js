// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the circle to be a solid (disk).
//
// The line is P + t * D, where D is required to be unit length. The FIQuery
// intersection parameters are valid only when D is unit length.
//
// The circle has center C and radius R.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// 'protected void DoQuery(...)' helper (used by the Ray2/Segment2 versus
// Circle2 queries, which derive from this class) becomes the protected method
// 'doQuery' that mutates the passed-in result, as upstream does.

import type { Hypersphere } from './Hypersphere.js';
import type { Line } from './Line.js';
import { DistPointLine } from './DistPointLine.js';
import { Vector, add, sub, mul, dot } from './Vector.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The result of IntrLine2Circle2TI.test.
export interface IntrLine2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrLine2Circle2TIResult(): IntrLine2Circle2TIResult {
    return { intersect: false };
}

// The result of IntrLine2Circle2FI.find.
export interface IntrLine2Circle2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine2Circle2FIResult(): IntrLine2Circle2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

export class IntrLine2Circle2TI implements
    TIQuery<Line, Hypersphere, IntrLine2Circle2TIResult> {

    test(line: Line, circle: Hypersphere): IntrLine2Circle2TIResult {
        const result = defaultIntrLine2Circle2TIResult();
        const plQuery = new DistPointLine();
        const plResult = plQuery.compute(circle.center, line);
        result.intersect = (plResult.distance <= circle.radius);
        return result;
    }
}

export class IntrLine2Circle2FI implements
    FIQuery<Line, Hypersphere, IntrLine2Circle2FIResult> {

    find(line: Line, circle: Hypersphere): IntrLine2Circle2FIResult {
        const result = defaultIntrLine2Circle2FIResult();
        this.doQuery(line.origin, line.direction, circle, result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(line.origin,
                mul(result.parameter[i], line.direction));
        }
        return result;
    }

    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        circle: Hypersphere, result: IntrLine2Circle2FIResult): void {
        // Intersection of the line P+t*D and the circle |X-C| = R. The line
        // direction is unit length. The t-value is a real-valued root to the
        // quadratic equation
        //   0 = |t*D+P-C|^2 - R^2
        //     = t^2 + 2*Dot(D,P-C)*t + |P-C|^2-R^2
        //     = t^2 + 2*a1*t + a0
        // If there are two distinct roots, the order is t0 < t1.
        const diff = sub(lineOrigin, circle.center);
        const a0 = dot(diff, diff) - circle.radius * circle.radius;
        const a1 = dot(lineDirection, diff);
        const discr = a1 * a1 - a0;
        if (discr > 0) {
            const root = Math.sqrt(discr);
            result.intersect = true;
            result.numIntersections = 2;
            result.parameter[0] = -a1 - root;
            result.parameter[1] = -a1 + root;
        } else if (discr < 0) {
            result.intersect = false;
            result.numIntersections = 0;
        } else {
            // The line is tangent to the circle. Set the parameters to the
            // same number because other queries involving linear components
            // and circular components use interval-interval intersection
            // tests which consume both parameters.
            result.intersect = true;
            result.numIntersections = 1;
            result.parameter[0] = -a1;
            result.parameter[1] = -a1;
        }
    }
}
