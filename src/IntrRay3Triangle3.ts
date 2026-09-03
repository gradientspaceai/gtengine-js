// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a ray and a solid triangle in 3D. The ray
// is P + t * D for t >= 0. When the ray is in the plane of the triangle, the
// queries state that there are no intersections.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrRay3Triangle3TI (test) and
// IntrRay3Triangle3FI (find), with the result types IntrRay3Triangle3TIResult
// and IntrRay3Triangle3FIResult.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Ray } from './Ray.js';
import type { Triangle } from './Triangle.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { cross, dotCross } from './Vector3.js';

// The result of IntrRay3Triangle3TI queries.
export interface IntrRay3Triangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay3Triangle3TIResult {
    return { intersect: false };
}

// The result of IntrRay3Triangle3FI queries.
export interface IntrRay3Triangle3FIResult {
    intersect: boolean;

    // The ray parameter t at the intersection point.
    parameter: number;

    // The barycentric coordinates (b0,b1,b2) of the intersection point
    // relative to the triangle vertices, b0 + b1 + b2 = 1.
    triangleBary: [number, number, number];

    // The intersection point.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrRay3Triangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        triangleBary: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a ray and a solid triangle in 3D.
export class IntrRay3Triangle3TI implements TIQuery<Ray, Triangle, IntrRay3Triangle3TIResult> {
    test(ray: Ray, triangle: Triangle): IntrRay3Triangle3TIResult {
        const result = defaultTIResult();

        // Compute the offset origin, edges, and normal.
        const diff = sub(ray.origin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = ray direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(ray.direction, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Ray and triangle are parallel, call it a "no intersection"
            // even if the ray does intersect.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(ray.direction, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(ray.direction, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle, check whether ray does.
                    const QdN = -sign * dot(diff, normal);
                    if (QdN >= 0) {
                        // Ray intersects triangle.
                        result.intersect = true;
                        return result;
                    }
                    // else: t < 0, no intersection
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

// Find-intersection query for a ray and a solid triangle in 3D.
export class IntrRay3Triangle3FI implements FIQuery<Ray, Triangle, IntrRay3Triangle3FIResult> {
    find(ray: Ray, triangle: Triangle): IntrRay3Triangle3FIResult {
        const result = defaultFIResult();

        // Compute the offset origin, edges, and normal.
        const diff = sub(ray.origin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = ray direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(ray.direction, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Ray and triangle are parallel, call it a "no intersection"
            // even if the ray does intersect.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(ray.direction, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(ray.direction, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle, check whether ray does.
                    const QdN = -sign * dot(diff, normal);
                    if (QdN >= 0) {
                        // Ray intersects triangle.
                        result.intersect = true;
                        result.parameter = QdN / DdN;
                        result.triangleBary[1] = DdQxE2 / DdN;
                        result.triangleBary[2] = DdE1xQ / DdN;
                        result.triangleBary[0] =
                            1 - result.triangleBary[1] - result.triangleBary[2];
                        result.point = add(ray.origin,
                            mul(result.parameter, ray.direction));
                        return result;
                    }
                    // else: t < 0, no intersection
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
