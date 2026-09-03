// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the ellipsoid to be a solid.
//
// The ellipsoid is (X-C)^T*M*(X-C)-1 = 0 and the line is X = P+t*D.
// Substitute the line equation into the ellipsoid equation to obtain a
// quadratic equation Q(t) = a2*t^2 + 2*a1*t + a0 = 0, where a2 = D^T*M*D,
// a1 = D^T*M*(P-C) and a0 = (P-C)^T*M*(P-C)-1. The algorithm involves an
// analysis of the real-valued roots of Q(t) for all real t.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrLine3Ellipsoid3TI
// and IntrLine3Ellipsoid3FI. The protected FIQuery::DoQuery is also exported
// as the module function 'intrLine3Ellipsoid3FIDoQuery'.

import type { FIQuery } from './FIQuery.js';
import type { Ellipsoid3 } from './Hyperellipsoid.js';
import type { Line3 } from './Line.js';
import { logAssert } from './Logger.js';
import { mulMatrix } from './Matrix.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrLine3Ellipsoid3TI.test.
export interface IntrLine3Ellipsoid3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine3Ellipsoid3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Ellipsoid3FI.find.
export interface IntrLine3Ellipsoid3FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Ellipsoid3FIResult():
    IntrLine3Ellipsoid3FIResult {
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
export function intrLine3Ellipsoid3FIDoQuery(lineOrigin: Vector,
    lineDirection: Vector, ellipsoid: Ellipsoid3,
    result: IntrLine3Ellipsoid3FIResult): void {
    logAssert(ellipsoid.dimension === 3,
        'IntrLine3Ellipsoid3FI: mismatched sizes.');

    const M = ellipsoid.getM();
    const diff = sub(lineOrigin, ellipsoid.center);
    const matDir = mulMatrix(M, lineDirection) as Vector;
    const matDiff = mulMatrix(M, diff) as Vector;
    const a2 = dot(lineDirection, matDir);
    const a1 = dot(lineDirection, matDiff);
    const a0 = dot(diff, matDiff) - 1;

    // Intersection occurs when Q(t) has real roots.
    const discr = a1 * a1 - a0 * a2;
    if (discr > 0) {
        // The line intersects the ellipsoid in 2 distinct points.
        result.intersect = true;
        result.numIntersections = 2;
        const root = Math.sqrt(discr);
        result.parameter[0] = (-a1 - root) / a2;
        result.parameter[1] = (-a1 + root) / a2;
    }
    else if (discr === 0) {
        // The line is tangent to the ellipsoid, so the intersection is a
        // single point. The parameter[1] value is set, because callers will
        // access the degenerate interval [-a1/a2, -a1/a2].
        result.intersect = true;
        result.numIntersections = 1;
        result.parameter[0] = -a1 / a2;
        result.parameter[1] = result.parameter[0];
    }
    // else: The line is outside the ellipsoid, no intersection.
}

// Test-intersection query for a line and a solid ellipsoid in 3D.
export class IntrLine3Ellipsoid3TI implements
    TIQuery<Line3, Ellipsoid3, IntrLine3Ellipsoid3TIResult> {

    test(line: Line3, ellipsoid: Ellipsoid3): IntrLine3Ellipsoid3TIResult {
        logAssert(ellipsoid.dimension === 3,
            'IntrLine3Ellipsoid3TI: mismatched sizes.');

        const result = defaultTIResult();

        const M = ellipsoid.getM();
        const diff = sub(line.origin, ellipsoid.center);
        const matDir = mulMatrix(M, line.direction) as Vector;
        const matDiff = mulMatrix(M, diff) as Vector;
        const a2 = dot(line.direction, matDir);
        const a1 = dot(line.direction, matDiff);
        const a0 = dot(diff, matDiff) - 1;

        // An intersection occurs when Q(t) has real roots.
        const discr = a1 * a1 - a0 * a2;
        result.intersect = (discr >= 0);
        return result;
    }
}

// Find-intersection query for a line and a solid ellipsoid in 3D.
export class IntrLine3Ellipsoid3FI implements
    FIQuery<Line3, Ellipsoid3, IntrLine3Ellipsoid3FIResult> {

    find(line: Line3, ellipsoid: Ellipsoid3): IntrLine3Ellipsoid3FIResult {
        const result = defaultIntrLine3Ellipsoid3FIResult();
        this.doQuery(line.origin, line.direction, ellipsoid, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(line.origin,
                    mul(result.parameter[i], line.direction));
            }
        }
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        ellipsoid: Ellipsoid3, result: IntrLine3Ellipsoid3FIResult): void {
        intrLine3Ellipsoid3FIDoQuery(lineOrigin, lineDirection, ellipsoid,
            result);
    }
}
