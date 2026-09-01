// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a line and a solid triangle in 3D. The
// line is P + t * D. When the line is in the plane of the triangle, the
// queries state that there are no intersections.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrLine3Triangle3TI (test) and
// IntrLine3Triangle3FI (find), with the result types
// IntrLine3Triangle3TIResult and IntrLine3Triangle3FIResult.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { Line } from './Line';
import type { Triangle } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';
import { cross, dotCross } from './Vector3';

// The result of IntrLine3Triangle3TI queries.
export interface IntrLine3Triangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine3Triangle3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Triangle3FI queries.
export interface IntrLine3Triangle3FIResult {
    intersect: boolean;

    // The line parameter t at the intersection point.
    parameter: number;

    // The barycentric coordinates (b0,b1,b2) of the intersection point
    // relative to the triangle vertices, b0 + b1 + b2 = 1.
    triangleBary: [number, number, number];

    // The intersection point.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine3Triangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        triangleBary: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a line and a solid triangle in 3D.
export class IntrLine3Triangle3TI implements TIQuery<Line, Triangle, IntrLine3Triangle3TIResult> {
    test(line: Line, triangle: Triangle): IntrLine3Triangle3TIResult {
        const result = defaultTIResult();

        // Compute the offset origin, edges, and normal.
        const diff = sub(line.origin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(line.direction, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Line and triangle are parallel, call it a "no intersection"
            // even if the line and triangle are coplanar and intersecting.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(line.direction, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(line.direction, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle.
                    result.intersect = true;
                    return result;
                }
                // else: b1+b2 > 1, no intersection
            }
            // else: b2 < 0, no intersection
        }
        // else: b1 < 0, no intersection

        result.intersect = false;
        return result;
    }
}

// Find-intersection query for a line and a solid triangle in 3D.
export class IntrLine3Triangle3FI implements FIQuery<Line, Triangle, IntrLine3Triangle3FIResult> {
    find(line: Line, triangle: Triangle): IntrLine3Triangle3FIResult {
        const result = defaultFIResult();

        // Compute the offset origin, edges, and normal.
        const diff = sub(line.origin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = line direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(line.direction, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Line and triangle are parallel, call it a "no intersection"
            // even if the line and triangle are coplanar and intersecting.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(line.direction, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(line.direction, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle.
                    const QdN = -sign * dot(diff, normal);
                    const inv = 1 / DdN;

                    result.intersect = true;
                    result.parameter = QdN * inv;
                    result.triangleBary[1] = DdQxE2 * inv;
                    result.triangleBary[2] = DdE1xQ * inv;
                    result.triangleBary[0] =
                        1 - result.triangleBary[1] - result.triangleBary[2];
                    result.point = add(line.origin,
                        mul(result.parameter, line.direction));
                    return result;
                }
                // else: b1+b2 > 1, no intersection
            }
            // else: b2 < 0, no intersection
        }
        // else: b1 < 0, no intersection

        result.intersect = false;
        return result;
    }
}
