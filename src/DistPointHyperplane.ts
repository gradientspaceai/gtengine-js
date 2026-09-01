// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointHyperplane.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a point and a line (N = 2), between a point
// and a plane (N = 3) or generally between a point and a hyperplane (N >= 2).
//
// The plane is defined by Dot(N, X - P) = 0, where P is the plane origin and
// N is a unit-length normal for the plane.
//
// TODO (upstream): Modify to support non-unit-length N.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, Hyperplane<N,T>>' becomes
// the class DistPointHyperplane with the result type
// DistPointHyperplaneResult. The upstream 'static_assert(N >= 2)' becomes a
// runtime logAssert, and the dimension aliases are dropped because the
// runtime-dimension Vector serves every N.

import type { DCPQuery } from './DCPQuery';
import type { Hyperplane } from './Hyperplane';
import { logAssert } from './Logger';
import { Vector, dot, mul, sub } from './Vector';

export interface DistPointHyperplaneResult {
    distance: number;

    // The signed distance Dot(N, P) - c. It is positive when the point is on
    // the side of the hyperplane to which the normal points, negative when it
    // is on the other side, and zero when the point is on the hyperplane.
    signedDistance: number;

    // closest[0] is the input point, closest[1] is the closest hyperplane
    // point.
    closest: [Vector, Vector];
}

export class DistPointHyperplane
    implements DCPQuery<Vector, Hyperplane, DistPointHyperplaneResult> {
    compute(point: Vector, plane: Hyperplane): DistPointHyperplaneResult {
        logAssert(point.size >= 2 && point.size === plane.normal.size,
            'Invalid dimension.');

        const signedDistance = dot(plane.normal, point) - plane.constant;
        return {
            distance: Math.abs(signedDistance),
            signedDistance,
            closest: [
                point.clone(),
                sub(point, mul(signedDistance, plane.normal))
            ]
        };
    }
}
