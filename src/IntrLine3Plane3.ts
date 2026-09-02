// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Plane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the plane to be a 2-dimensional object embedded in 3D
// and the line to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrLine3Plane3TI and
// IntrLine3Plane3FI. The protected 'FIQuery::DoQuery' becomes the exported
// module function 'intrLine3Plane3FIDoQuery' (the precedent set by
// IntrLine3Sphere3.ts) so that derived ray/segment queries can reuse it.
// 'std::numeric_limits<int32_t>::max()' becomes the literal 2147483647 (the
// precedent set by IntrLine2Line2.ts).

import { DistPointHyperplane } from './DistPointHyperplane';
import type { FIQuery } from './FIQuery';
import type { Line3 } from './Line';
import type { Plane3 } from './Hyperplane';
import { Vector, add, dot, mul } from './Vector';
import type { TIQuery } from './TIQuery';

const INT32_MAX = 2147483647;

// The result of IntrLine3Plane3TI.test.
export interface IntrLine3Plane3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrLine3Plane3TIResult(): IntrLine3Plane3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Plane3FI.find.
export interface IntrLine3Plane3FIResult {
    intersect: boolean;

    // The number of intersections is 0 (no intersection), 1 (the linear
    // component and the plane intersect in a point), or 2147483647 (the
    // linear component is on the plane). When the linear component is on the
    // plane, 'point' is the component origin and 'parameter' is zero.
    numIntersections: number;
    parameter: number;
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Plane3FIResult(): IntrLine3Plane3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: 0,
        point: Vector.zero(3)
    };
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrLine3Plane3FIDoQuery(lineOrigin: Vector,
    lineDirection: Vector, plane: Plane3,
    result: IntrLine3Plane3FIResult): void {
    const DdN = dot(lineDirection, plane.normal);
    const vpQuery = new DistPointHyperplane();
    const vpResult = vpQuery.compute(lineOrigin, plane);

    if (DdN !== 0) {
        // The line is not parallel to the plane, so they must intersect.
        result.intersect = true;
        result.numIntersections = 1;
        result.parameter = -vpResult.signedDistance / DdN;
    }
    else {
        // The line and plane are parallel. Determine whether the line is on
        // the plane.
        if (vpResult.distance === 0) {
            // The line is coincident with the plane, so choose t = 0 for the
            // parameter.
            result.intersect = true;
            result.numIntersections = INT32_MAX;
            result.parameter = 0;
        }
        else {
            // The line is not on the plane.
            result.intersect = false;
            result.numIntersections = 0;
        }
    }
}

// Test-intersection query for a line and a plane in 3D.
export class IntrLine3Plane3TI implements
    TIQuery<Line3, Plane3, IntrLine3Plane3TIResult> {

    test(line: Line3, plane: Plane3): IntrLine3Plane3TIResult {
        const result = defaultIntrLine3Plane3TIResult();

        const DdN = dot(line.direction, plane.normal);
        if (DdN !== 0) {
            // The line is not parallel to the plane, so they must intersect.
            result.intersect = true;
        }
        else {
            // The line and plane are parallel.
            const vpQuery = new DistPointHyperplane();
            result.intersect = (vpQuery.compute(line.origin, plane).distance === 0);
        }

        return result;
    }
}

// Find-intersection query for a line and a plane in 3D.
export class IntrLine3Plane3FI implements
    FIQuery<Line3, Plane3, IntrLine3Plane3FIResult> {

    find(line: Line3, plane: Plane3): IntrLine3Plane3FIResult {
        const result = defaultIntrLine3Plane3FIResult();
        intrLine3Plane3FIDoQuery(line.origin, line.direction, plane, result);
        if (result.intersect) {
            result.point = add(line.origin, mul(result.parameter, line.direction));
        }
        return result;
    }
}
