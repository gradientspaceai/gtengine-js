// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the arc to be a 1-dimensional object.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrRay2Arc2TI and
// IntrRay2Arc2FI. Upstream calls the single-argument Arc2::Contains, which
// the port names 'containsOnCircle' (see IntrLine2Arc2.ts).
//
// Port fix for an upstream bug. Upstream runs FIQuery<Ray2,Circle2>, which
// treats the circle as a SOLID disk: it clips the line t-interval to
// [0,+infinity), so when the ray origin is inside the disk the first reported
// "intersection point" is the ray origin itself, which is not on the circle.
// Upstream then hands that interior point to Arc2::Contains, whose
// single-argument form assumes the point is on the circle, and the query can
// report an intersection at a point that is not on the arc. The port instead
// intersects the ray with the circular CURVE (the line-circle query filtered
// by t >= 0, the same technique upstream itself uses in
// IntrRay2SegmentMesh2.h), which yields the same answers whenever the ray
// origin is outside the disk and the correct answers when it is inside.

import type { Arc2 } from './Arc2';
import type { FIQuery } from './FIQuery';
import { Hypersphere } from './Hypersphere';
import { IntrLine2Circle2FI } from './IntrLine2Circle2';
import { Line } from './Line';
import type { Ray2 } from './Ray';
import { Vector } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrRay2Arc2TI.test.
export interface IntrRay2Arc2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrRay2Arc2TIResult(): IntrRay2Arc2TIResult {
    return { intersect: false };
}

// The result of IntrRay2Arc2FI.find.
export interface IntrRay2Arc2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay2Arc2FIResult(): IntrRay2Arc2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test-intersection query for a ray and an arc in 2D.
export class IntrRay2Arc2TI implements
    TIQuery<Ray2, Arc2, IntrRay2Arc2TIResult> {

    test(ray: Ray2, arc: Arc2): IntrRay2Arc2TIResult {
        const result = defaultIntrRay2Arc2TIResult();
        const raQuery = new IntrRay2Arc2FI();
        const raResult = raQuery.find(ray, arc);
        result.intersect = raResult.intersect;
        return result;
    }
}

// Find-intersection query for a ray and an arc in 2D.
export class IntrRay2Arc2FI implements
    FIQuery<Ray2, Arc2, IntrRay2Arc2FIResult> {

    find(ray: Ray2, arc: Arc2): IntrRay2Arc2FIResult {
        const result = defaultIntrRay2Arc2FIResult();

        const lcQuery = new IntrLine2Circle2FI();
        const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const rcResult = lcQuery.find(line, circle);
        if (rcResult.intersect) {
            // Test whether the line-circle intersections are on the ray and
            // on the arc.
            result.numIntersections = 0;
            for (let i = 0; i < rcResult.numIntersections; ++i) {
                if (rcResult.parameter[i] >= 0
                    && arc.containsOnCircle(rcResult.point[i])) {
                    result.intersect = true;
                    result.parameter[result.numIntersections] =
                        rcResult.parameter[i];
                    result.point[result.numIntersections] =
                        rcResult.point[i].clone();
                    ++result.numIntersections;
                }
            }
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }

        return result;
    }
}
