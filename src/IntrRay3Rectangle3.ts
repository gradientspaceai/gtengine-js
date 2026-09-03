// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a ray and a solid rectangle in 3D.
//
// The ray is P + t * D for t >= 0, where D is not required to be unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The intersection point, if any, is stored in result.point. The
// corresponding ray parameter t is stored in result.parameter. The
// corresponding rectangle parameters s[] are stored in result.rectCoord.
// When the ray is in the plane of the rectangle and intersects the rectangle,
// the queries state that there are no intersections.
//
// Upstream TODO: modify to support non-unit-length W[], and return the point
// or segment of intersection when the ray is in the plane of the rectangle
// and intersects the rectangle.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// specializations become IntrRay3Rectangle3TI and IntrRay3Rectangle3FI. As in
// IntrLine3Rectangle3.ts, 'rectCoord' keeps the upstream length-3 array even
// though only the first two entries have meaning.

import { Line } from './Line.js';
import type { Ray } from './Ray.js';
import type { Rectangle } from './Rectangle.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector } from './Vector.js';
import { IntrLine3Rectangle3FI } from './IntrLine3Rectangle3.js';

// The result of IntrRay3Rectangle3TI queries.
export interface IntrRay3Rectangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrRay3Rectangle3TIResult {
    return { intersect: false };
}

// The result of IntrRay3Rectangle3FI queries.
export interface IntrRay3Rectangle3FIResult {
    intersect: boolean;

    // The ray parameter t at the intersection point.
    parameter: number;

    // The rectangle coordinates (s[0],s[1]) of the intersection point. See
    // the note in IntrLine3Rectangle3.ts about the length-3 array.
    rectCoord: [number, number, number];

    // The intersection point.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrRay3Rectangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        rectCoord: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a ray and a solid rectangle in 3D.
export class IntrRay3Rectangle3TI implements
    TIQuery<Ray, Rectangle, IntrRay3Rectangle3TIResult> {

    test(ray: Ray, rectangle: Rectangle): IntrRay3Rectangle3TIResult {
        const result = defaultTIResult();

        const lrQuery = new IntrLine3Rectangle3FI();
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = lrQuery.find(line, rectangle);
        if (lrResult.intersect) {
            if (lrResult.parameter >= 0) {
                // The line-rectangle intersection is on the ray.
                result.intersect = true;
                return result;
            }
        }

        result.intersect = false;
        return result;
    }
}

// Find-intersection query for a ray and a solid rectangle in 3D.
export class IntrRay3Rectangle3FI implements
    FIQuery<Ray, Rectangle, IntrRay3Rectangle3FIResult> {

    find(ray: Ray, rectangle: Rectangle): IntrRay3Rectangle3FIResult {
        const result = defaultFIResult();

        const lrQuery = new IntrLine3Rectangle3FI();
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const lrResult = lrQuery.find(line, rectangle);
        if (lrResult.intersect) {
            if (lrResult.parameter >= 0) {
                // The line-rectangle intersection is on the ray.
                result.intersect = true;
                result.parameter = lrResult.parameter;
                result.rectCoord = [
                    lrResult.rectCoord[0],
                    lrResult.rectCoord[1],
                    lrResult.rectCoord[2]
                ];
                result.point = lrResult.point;
                return result;
            }
        }

        result.intersect = false;
        return result;
    }
}
