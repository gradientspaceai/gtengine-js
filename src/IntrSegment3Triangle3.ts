// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a segment and a solid triangle in 3D. The
// queries use the centered form of the segment, C + s * D with |s| <= e.
// When the segment is in the plane of the triangle, the queries state that
// there are no intersections.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrSegment3Triangle3TI (test)
// and IntrSegment3Triangle3FI (find), with the result types
// IntrSegment3Triangle3TIResult and IntrSegment3Triangle3FIResult. The
// reported 'parameter' is the parameter s of the centered segment form, as
// upstream.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { Segment } from './Segment';
import type { Triangle } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';
import { cross, dotCross } from './Vector3';

// The result of IntrSegment3Triangle3TI queries.
export interface IntrSegment3Triangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrSegment3Triangle3TIResult {
    return { intersect: false };
}

// The result of IntrSegment3Triangle3FI queries.
export interface IntrSegment3Triangle3FIResult {
    intersect: boolean;

    // The parameter s of the centered segment form C + s * D at the
    // intersection point, |s| <= extent.
    parameter: number;

    // The barycentric coordinates (b0,b1,b2) of the intersection point
    // relative to the triangle vertices, b0 + b1 + b2 = 1.
    triangleBary: [number, number, number];

    // The intersection point.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrSegment3Triangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        triangleBary: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a segment and a solid triangle in 3D.
export class IntrSegment3Triangle3TI implements TIQuery<Segment, Triangle, IntrSegment3Triangle3TIResult> {
    test(segment: Segment, triangle: Triangle): IntrSegment3Triangle3TIResult {
        const result = defaultTIResult();

        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        // Compute the offset origin, edges, and normal.
        const diff = sub(segOrigin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(segDirection, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Segment and triangle are parallel, call it a "no intersection"
            // even if the segment does intersect.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(segDirection, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(segDirection, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle, check whether segment does.
                    const QdN = -sign * dot(diff, normal);
                    const extDdN = segExtent * DdN;
                    if (-extDdN <= QdN && QdN <= extDdN) {
                        // Segment intersects triangle.
                        result.intersect = true;
                        return result;
                    }
                    // else: |t| > extent, no intersection
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

// Find-intersection query for a segment and a solid triangle in 3D.
export class IntrSegment3Triangle3FI implements FIQuery<Segment, Triangle, IntrSegment3Triangle3FIResult> {
    find(segment: Segment, triangle: Triangle): IntrSegment3Triangle3FIResult {
        const result = defaultFIResult();

        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        // Compute the offset origin, edges, and normal.
        const diff = sub(segOrigin, triangle.v[0]);
        const edge1 = sub(triangle.v[1], triangle.v[0]);
        const edge2 = sub(triangle.v[2], triangle.v[0]);
        const normal = cross(edge1, edge2);

        // Solve Q + t*D = b1*E1 + b2*E2 (Q = diff, D = segment direction,
        // E1 = edge1, E2 = edge2, N = Cross(E1,E2)) by
        //   |Dot(D,N)|*b1 = sign(Dot(D,N))*Dot(D,Cross(Q,E2))
        //   |Dot(D,N)|*b2 = sign(Dot(D,N))*Dot(D,Cross(E1,Q))
        //   |Dot(D,N)|*t = -sign(Dot(D,N))*Dot(Q,N)
        let DdN = dot(segDirection, normal);
        let sign: number;
        if (DdN > 0) {
            sign = 1;
        }
        else if (DdN < 0) {
            sign = -1;
            DdN = -DdN;
        }
        else {
            // Segment and triangle are parallel, call it a "no intersection"
            // even if the segment does intersect.
            result.intersect = false;
            return result;
        }

        const DdQxE2 = sign * dotCross(segDirection, diff, edge2);
        if (DdQxE2 >= 0) {
            const DdE1xQ = sign * dotCross(segDirection, edge1, diff);
            if (DdE1xQ >= 0) {
                if (DdQxE2 + DdE1xQ <= DdN) {
                    // Line intersects triangle, check whether segment does.
                    const QdN = -sign * dot(diff, normal);
                    const extDdN = segExtent * DdN;
                    if (-extDdN <= QdN && QdN <= extDdN) {
                        // Segment intersects triangle.
                        result.intersect = true;
                        result.parameter = QdN / DdN;
                        result.triangleBary[1] = DdQxE2 / DdN;
                        result.triangleBary[2] = DdE1xQ / DdN;
                        result.triangleBary[0] =
                            1 - result.triangleBary[1] - result.triangleBary[2];
                        result.point = add(segOrigin,
                            mul(result.parameter, segDirection));
                        return result;
                    }
                    // else: |t| > extent, no intersection
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
