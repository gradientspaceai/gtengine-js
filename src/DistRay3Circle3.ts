// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistRay3Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The 3D ray-circle distance algorithm is described in
// https://www.geometrictools.com/Documentation/DistanceToCircle3.pdf
// The notation used in the code matches that of the document. The circle has
// center C and the plane of the circle has unit-length normal N. The ray has
// origin B and non-zero direction M. The parameterization is P(t) = t*M+B. It
// is not necessary that M be a unit-length vector.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Ray3<T>, Circle3<T>>' becomes the
// class DistRay3Circle3. As upstream does, the result type is the line-circle
// result type (with 'linearClosest' naming the ray points), re-exported here
// as the alias DistRay3Circle3Result. The upstream friend access to
// 'DCPQuery<T,Line3,Circle3>::Execute' and its private 'Critical' struct are
// the exported 'distLine3Circle3Execute' and 'DistLine3Circle3Critical' from
// DistLine3Circle3.ts. The private helpers Execute, HasOneCriticalPoint,
// HasTwoCriticalPoints, RayOriginClosest and SelectClosestPoint become
// module-private functions.

import type { Circle3 } from './Circle3';
import type { DCPQuery } from './DCPQuery';
import { distLine3Circle3Execute } from './DistLine3Circle3';
import type { DistLine3Circle3Critical, DistLine3Circle3Result }
    from './DistLine3Circle3';
import { DistPoint3Circle3 } from './DistPoint3Circle3';
import { Line } from './Line';
import type { Ray3 } from './Ray';
import { Vector } from './Vector';

// Upstream reuses the line-circle result type ('using Result = typename
// LCQuery::Result').
export type DistRay3Circle3Result = DistLine3Circle3Result;

function rayOriginClosest(rayOrigin: Vector, circle: Circle3,
    result: DistRay3Circle3Result): void {
    const pcResult = new DistPoint3Circle3().compute(rayOrigin, circle);
    result.numClosestPairs = 1;
    result.linearClosest[0] = rayOrigin.clone();
    result.linearClosest[1] = new Vector(3);
    result.circularClosest[0] = pcResult.closest[1];
    result.circularClosest[1] = new Vector(3);
    result.distance = pcResult.distance;
    result.sqrDistance = result.distance * result.distance;
}

function selectClosestPoint(point0: Vector, point1: Vector, circle: Circle3,
    result: DistRay3Circle3Result): void {
    const pcQuery = new DistPoint3Circle3();
    const pcResult0 = pcQuery.compute(point0, circle);
    const pcResult1 = pcQuery.compute(point1, circle);
    if (pcResult0.distance < pcResult1.distance) {
        result.numClosestPairs = 1;
        result.linearClosest[0] = point0.clone();
        result.linearClosest[1] = new Vector(3);
        result.circularClosest[0] = pcResult0.closest[1];
        result.circularClosest[1] = new Vector(3);
        result.distance = pcResult0.distance;
        result.sqrDistance = result.distance * result.distance;
    }
    else if (pcResult0.distance > pcResult1.distance) {
        result.numClosestPairs = 1;
        result.linearClosest[0] = point1.clone();
        result.linearClosest[1] = new Vector(3);
        result.circularClosest[0] = pcResult1.closest[1];
        result.circularClosest[1] = new Vector(3);
        result.distance = pcResult1.distance;
        result.sqrDistance = result.distance * result.distance;
    }
    else {
        // pcResult0.distance = pcResult1.distance
        result.numClosestPairs = 2;
        result.linearClosest[0] = point0.clone();
        result.linearClosest[1] = point1.clone();
        result.circularClosest[0] = pcResult0.closest[1];
        result.circularClosest[1] = pcResult1.closest[1];
        result.distance = pcResult0.distance;
        result.sqrDistance = result.distance * result.distance;
    }
}

function hasOneCriticalPoint(ray: Ray3, circle: Circle3,
    critical: DistLine3Circle3Critical, result: DistRay3Circle3Result): void {
    const t0 = critical.parameter[0];

    if (t0 <= 0) {
        // The critical point is not on the ray. The ray origin is the ray
        // point closest to the circle. See the red ray of the
        // one-critical-point graph of figure 7 in the PDF.
        rayOriginClosest(ray.origin, circle, result);
        return;
    }

    // At this time, t0 > 0. The closest line-circle pair is the closest
    // ray-circle pair. The output does not need to be modified. See the green
    // ray of the one-critical-point graph of figure 7 in the PDF.
}

function hasTwoCriticalPoints(ray: Ray3, circle: Circle3,
    critical: DistLine3Circle3Critical, result: DistRay3Circle3Result): void {
    const t0 = critical.parameter[0];
    const t1 = critical.parameter[1];

    if (t0 >= 0) {
        // The critical points are on the ray. The ray point closest to the
        // circle is the line point closest to the circle. The output remains
        // unchanged. See the green rays of the two-critical-point graphs of
        // figure 7 in the PDF.
        return;
    }

    if (t1 <= 0) {
        // The critical points are not on the ray. The ray origin is the ray
        // point closest to the circle. See the red rays of the
        // two-critical-point graphs of figure 7 in the PDF.
        rayOriginClosest(ray.origin, circle, result);
        return;
    }

    // The ray point closest to the circle is either the ray origin or the
    // second critical point, whichever has minimum distance. See the orange
    // and purple rays of the two-critical-point graphs of figure 7 in the
    // PDF.
    selectClosestPoint(ray.origin, critical.linearPoint[1], circle, result);
}

export class DistRay3Circle3
    implements DCPQuery<Ray3, Circle3, DistRay3Circle3Result> {
    compute(ray: Ray3, circle: Circle3): DistRay3Circle3Result {
        // Compute the line points closest to the circle. The line is
        // L(t) = P + t * D for any real-valued t. The ray restricts t >= 0
        // and has origin P = L(0).
        const line = Line.fromOriginDirection(ray.origin, ray.direction);
        const { result, critical } = distLine3Circle3Execute(line, circle);

        // Clamp the query output to the ray domain.
        if (critical.numPoints === 1) {
            hasOneCriticalPoint(ray, circle, critical, result);
        }
        else {
            hasTwoCriticalPoints(ray, circle, critical, result);
        }

        return result;
    }
}
