// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the sphere to be a solid.
//
// The sphere is (X-C)^T*(X-C)-r^2 = 0 and the line is X = P+t*D. Substitute
// the line equation into the sphere equation to obtain a quadratic equation
// Q(t) = t^2 + 2*a1*t + a0 = 0, where a1 = D^T*(P-C) and
// a0 = (P-C)^T*(P-C)-r^2. The algorithm involves an analysis of the
// real-valued roots of Q(t) for all real t.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrLine3Sphere3TI (test) and
// IntrLine3Sphere3FI (find), with the result types IntrLine3Sphere3TIResult
// and IntrLine3Sphere3FIResult. The protected 'DoQuery' member, which
// upstream derived classes (ray/segment vs sphere) call, is exported as the
// module function 'intrLine3Sphere3DoQuery'.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Hypersphere } from './Hypersphere.js';
import type { Line } from './Line.js';
import { Vector, add, dot, mul, sub } from './Vector.js';

// The result of IntrLine3Sphere3TI queries.
export interface IntrLine3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine3Sphere3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Sphere3FI queries.
export interface IntrLine3Sphere3FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Sphere3FIResult(): IntrLine3Sphere3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(3), Vector.zero(3)]
    };
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrLine3Sphere3DoQuery(lineOrigin: Vector,
    lineDirection: Vector, sphere: Hypersphere,
    result: IntrLine3Sphere3FIResult): void {
    const diff = sub(lineOrigin, sphere.center);
    const a0 = dot(diff, diff) - sphere.radius * sphere.radius;
    const a1 = dot(lineDirection, diff);

    // Intersection occurs when Q(t) has real roots.
    const discr = a1 * a1 - a0;
    if (discr > 0) {
        // The line intersects the sphere in 2 distinct points.
        result.intersect = true;
        result.numIntersections = 2;
        const root = Math.sqrt(discr);
        result.parameter[0] = -a1 - root;
        result.parameter[1] = -a1 + root;
    }
    else if (discr === 0) {
        // The line is tangent to the sphere, so the intersection is a single
        // point. The parameter[1] value is set, because callers will access
        // the degenerate interval [-a1,-a1].
        result.intersect = true;
        result.numIntersections = 1;
        result.parameter[0] = -a1;
        result.parameter[1] = result.parameter[0];
    }
    // else: The line is outside the sphere, no intersection.
}

// Test-intersection query for a line and a solid sphere in 3D.
export class IntrLine3Sphere3TI implements TIQuery<Line, Hypersphere, IntrLine3Sphere3TIResult> {
    test(line: Line, sphere: Hypersphere): IntrLine3Sphere3TIResult {
        const result = defaultTIResult();

        const diff = sub(line.origin, sphere.center);
        const a0 = dot(diff, diff) - sphere.radius * sphere.radius;
        const a1 = dot(line.direction, diff);

        // An intersection occurs when Q(t) has real roots.
        const discr = a1 * a1 - a0;
        result.intersect = (discr >= 0);
        return result;
    }
}

// Find-intersection query for a line and a solid sphere in 3D.
export class IntrLine3Sphere3FI implements FIQuery<Line, Hypersphere, IntrLine3Sphere3FIResult> {
    find(line: Line, sphere: Hypersphere): IntrLine3Sphere3FIResult {
        const result = defaultIntrLine3Sphere3FIResult();
        intrLine3Sphere3DoQuery(line.origin, line.direction, sphere, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(line.origin,
                    mul(result.parameter[i], line.direction));
            }
        }
        return result;
    }
}
