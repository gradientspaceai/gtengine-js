// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPoint3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to a cylinder that is finite or infinite.
// The queries consider the cylinder to be a solid.
//
// The input point is stored in the member closest[0]. The cylinder point
// closest to it is stored in the member closest[1].
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector3<T>, Cylinder3<T>>' becomes the
// class DistPoint3Cylinder3 with the result type DistPoint3Cylinder3Result.
// The private helpers DoQueryInfiniteCylinder and DoQueryFiniteCylinder
// become module-private functions.
//
// Upstream bug (fixed here): DistPoint3Cylinder3.h detects an infinite
// cylinder with 'cylinder.height == std::numeric_limits<T>::max()', but
// Cylinder3.h documents and implements the infinite-cylinder sentinel as
// height = -1 (MakeInfiniteCylinder). With the upstream test, an infinite
// cylinder built by MakeInfiniteCylinder falls into the finite branch and
// trips the 'positive height' assertion. The port uses the Cylinder3
// isInfinite() predicate instead.

import type { Cylinder3 } from './Cylinder3';
import type { DCPQuery } from './DCPQuery';
import { logAssert } from './Logger';
import { Vector, add, dot, mul, sub } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

export interface DistPoint3Cylinder3Result {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest cylinder
    // point.
    closest: [Vector, Vector];
}

// The point P and the returned closest point are in cylinder coordinates,
// where (0,0,0) is the cylinder axis origin and (0,0,1) is the cylinder axis
// direction.
function doQueryInfiniteCylinder(P: Vector, radius: number):
    { distance: number, sqrDistance: number, closest: Vector } {
    const sqrRadius = radius * radius;
    const sqrDistance = P.values[0] * P.values[0] + P.values[1] * P.values[1];
    if (sqrDistance >= sqrRadius) {
        // The point is outside the cylinder or on the cylinder wall.
        const dist = Math.sqrt(sqrDistance);
        const distance = dist - radius;
        const temp = radius / dist;
        return {
            distance,
            sqrDistance: distance * distance,
            closest: Vector.fromArray([
                P.values[0] * temp, P.values[1] * temp, P.values[2]])
        };
    }

    // The point is inside the cylinder.
    return { distance: 0, sqrDistance: 0, closest: P.clone() };
}

function doQueryFiniteCylinder(P: Vector, radius: number, height: number):
    { distance: number, sqrDistance: number, closest: Vector } {
    const query = doQueryInfiniteCylinder(P, radius);

    // Clamp the infinite cylinder's closest point to the finite cylinder.
    const halfHeight = 0.5 * height;
    if (query.closest.values[2] > halfHeight) {
        query.closest.values[2] = halfHeight;
    }
    else if (query.closest.values[2] < -halfHeight) {
        query.closest.values[2] = -halfHeight;
    }
    else {
        return query;
    }

    const diff = sub(query.closest, P);
    query.sqrDistance = dot(diff, diff);
    query.distance = Math.sqrt(query.sqrDistance);
    return query;
}

export class DistPoint3Cylinder3
    implements DCPQuery<Vector, Cylinder3, DistPoint3Cylinder3Result> {
    compute(point: Vector, cylinder: Cylinder3): DistPoint3Cylinder3Result {
        // Convert the point to the cylinder coordinate system. In this
        // system, the point believes (0,0,0) is the cylinder axis origin and
        // (0,0,1) is the cylinder axis direction.
        const basis: Vector[] = [cylinder.axis.direction.clone(),
            new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, basis);

        const delta = sub(point, cylinder.axis.origin);
        const P = Vector.fromArray([
            dot(basis[1], delta),
            dot(basis[2], delta),
            dot(basis[0], delta)
        ]);

        logAssert(cylinder.radius > 0,
            'The cylinder must have a positive radius.');

        let query: { distance: number, sqrDistance: number, closest: Vector };
        if (cylinder.isInfinite()) {
            query = doQueryInfiniteCylinder(P, cylinder.radius);
        }
        else {
            logAssert(cylinder.height > 0,
                'The cylinder must have a positive height.');

            query = doQueryFiniteCylinder(P, cylinder.radius, cylinder.height);
        }

        // Convert the closest point from the cylinder coordinate system to
        // the original coordinate system.
        const closest1 = add(cylinder.axis.origin,
            add(mul(query.closest.values[0], basis[1]),
                add(mul(query.closest.values[1], basis[2]),
                    mul(query.closest.values[2], basis[0]))));

        return {
            distance: query.distance,
            sqrDistance: query.sqrDistance,
            closest: [point.clone(), closest1]
        };
    }
}
