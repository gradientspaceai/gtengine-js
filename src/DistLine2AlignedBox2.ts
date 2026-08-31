// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistLine2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between a line and a solid aligned box in 2D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The aligned box has minimum corner A and maximum corner B. A box point is X
// where A <= X <= B; the comparisons are componentwise.
//
// The closest point on the line is stored in closest[0] with parameter t. The
// closest point on the box is stored in closest[1]. When there are infinitely
// many choices for the pair of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Line2<T>, AlignedBox2<T>>' becomes the
// class DistLine2AlignedBox2 with the result type
// DistLine2AlignedBox2Result. The protected/private static helpers DoQuery,
// DoQuery2D, DoQuery1D and DoQuery0D become the exported function
// distLine2AlignedBox2DoQuery (the upstream 'friend class' grant to
// DCPQuery<T, Line2<T>, OrientedBox2<T>> becomes an export, so that
// DistLine2OrientedBox2 can call it without derivation) and the
// module-private functions doQuery2D, doQuery1D and doQuery0D. Upstream takes
// the origin and direction by non-const reference so the reflections can be
// applied in place; the port does the same on Vector copies made by the
// caller.

import type { AlignedBox2 } from './AlignedBox';
import type { DCPQuery } from './DCPQuery';
import { clamp } from './Functions';
import type { Line2 } from './Line';
import { Vector, add, dot, mul, sub } from './Vector';
import { dotPerp } from './Vector2';

export interface DistLine2AlignedBox2Result {
    distance: number;
    sqrDistance: number;

    // The line parameter t of the closest line point.
    parameter: number;

    // closest[0] is on the line, closest[1] is on the box.
    closest: [Vector, Vector];
}

function defaultResult(): DistLine2AlignedBox2Result {
    return {
        distance: 0,
        sqrDistance: 0,
        parameter: 0,
        closest: [new Vector(2), new Vector(2)]
    };
}

function doQuery2D(origin: Vector, direction: Vector, extent: Vector,
    result: DistLine2AlignedBox2Result): void {
    const K0 = Vector.fromArray([-extent.values[0], extent.values[1]]);
    let delta = sub(K0, origin);
    const K0dotPerpD = dotPerp(delta, direction);
    if (K0dotPerpD >= 0) {
        result.parameter = dot(delta, direction) / dot(direction, direction);
        result.closest[0] = add(origin, mul(result.parameter, direction));
        result.closest[1] = K0;
    }
    else {
        const K1 = Vector.fromArray([extent.values[0], -extent.values[1]]);
        delta = sub(K1, origin);
        const K1dotPerpD = dotPerp(delta, direction);
        if (K1dotPerpD <= 0) {
            result.parameter =
                dot(delta, direction) / dot(direction, direction);
            result.closest[0] = add(origin, mul(result.parameter, direction));
            result.closest[1] = K1;
        }
        else {
            const K2 = Vector.fromArray([extent.values[0], extent.values[1]]);
            delta = sub(K2, origin);
            const K2dotPerpD = dotPerp(delta, direction);
            if (K2dotPerpD >= 0) {
                result.parameter =
                    (extent.values[1] - origin.values[1]) / direction.values[1];
                result.closest[0] =
                    add(origin, mul(result.parameter, direction));
                result.closest[1].values[0] =
                    origin.values[0] + result.parameter * direction.values[0];
                result.closest[1].values[1] = extent.values[1];
            }
            else {
                result.parameter =
                    (extent.values[0] - origin.values[0]) / direction.values[0];
                result.closest[0] =
                    add(origin, mul(result.parameter, direction));
                result.closest[1].values[0] = extent.values[0];
                result.closest[1].values[1] =
                    origin.values[1] + result.parameter * direction.values[1];
            }
        }
    }
}

function doQuery1D(i0: number, i1: number, origin: Vector, direction: Vector,
    extent: Vector, result: DistLine2AlignedBox2Result): void {
    result.parameter =
        (extent.values[i0] - origin.values[i0]) / direction.values[i0];
    result.closest[0] = add(origin, mul(result.parameter, direction));
    result.closest[1].values[i0] = extent.values[i0];
    result.closest[1].values[i1] =
        clamp(origin.values[i1], -extent.values[i1], extent.values[i1]);
}

function doQuery0D(origin: Vector, extent: Vector,
    result: DistLine2AlignedBox2Result): void {
    result.parameter = 0;
    result.closest[0] = origin.clone();
    result.closest[1].values[0] =
        clamp(origin.values[0], -extent.values[0], extent.values[0]);
    result.closest[1].values[1] =
        clamp(origin.values[1], -extent.values[1], extent.values[1]);
}

// Compute the distance and closest point between a line and an aligned box
// whose center is the origin. The origin and direction are mutated (they are
// reflected in place to eliminate complicated sign logic), so the caller must
// pass copies it does not need afterwards.
export function distLine2AlignedBox2DoQuery(origin: Vector,
    direction: Vector, extent: Vector,
    result: DistLine2AlignedBox2Result): void {
    // Apply reflections so that the direction has nonnegative components.
    const reflect: boolean[] = [false, false];
    for (let i = 0; i < 2; ++i) {
        if (direction.values[i] < 0) {
            origin.values[i] = -origin.values[i];
            direction.values[i] = -direction.values[i];
            reflect[i] = true;
        }
    }

    // Compute the line-box distance and closest points. The doQueryND calls
    // compute result.parameter and result.closest[1]. The result.closest[0]
    // can be computed after these calls.
    if (direction.values[0] > 0) {
        if (direction.values[1] > 0) {
            // The direction signs are (+,+). If the line does not intersect
            // the box, the only possible closest box points are K[0] =
            // (-e0,e1) or K[1] = (e0,-e1). If the line intersects the box,
            // the closest points are the same and chosen to be the
            // intersection with box edge x0 = e0 or x1 = e1. For the
            // remaining discussion, define K[2] = (e0,e1).
            //
            // Test where the candidate corners are relative to the line. If
            // D = (d0,d1), then Perp(D) = (d1,-d0). The corner K[i] =
            // P + t[i] * D + s[i] * Perp(D), where
            // s[i] = Dot(K[i]-P,Perp(D))/|D|^2. K[0] is closest when
            // s[0] >= 0 or K[1] is closest when s[1] <= 0. Otherwise, the
            // line intersects the box. If s[2] >= 0, the common closest point
            // is chosen to be (p0+(e1-p1)*d0/d1,e1). If s[2] < 0, the common
            // closest point is chosen to be (e0,p1+(e0-p0)*d1/d0).
            //
            // It is sufficient to test the signs of Dot(K[i],Perp(D)) and
            // defer the division by |D|^2 until needed for computing the
            // closest point.
            doQuery2D(origin, direction, extent, result);
        }
        else {
            // The direction signs are (+,0). The parameter is the value of t
            // for which P + t * D = (e0, p1).
            doQuery1D(0, 1, origin, direction, extent, result);
        }
    }
    else {
        if (direction.values[1] > 0) {
            // The direction signs are (0,+). The parameter is the value of t
            // for which P + t * D = (p0, e1).
            doQuery1D(1, 0, origin, direction, extent, result);
        }
        else {
            // The direction signs are (0,0). The line is degenerate to a
            // point (its origin). Clamp the origin to the box to obtain the
            // closest point.
            doQuery0D(origin, extent, result);
        }
    }

    result.closest[0] = add(origin, mul(result.parameter, direction));

    // Undo the reflections. The origin and direction are not consumed by the
    // caller, so these do not need to be reflected. However, the closest
    // points are consumed.
    for (let i = 0; i < 2; ++i) {
        if (reflect[i]) {
            for (let j = 0; j < 2; ++j) {
                result.closest[j].values[i] = -result.closest[j].values[i];
            }
        }
    }
}

export class DistLine2AlignedBox2
    implements DCPQuery<Line2, AlignedBox2, DistLine2AlignedBox2Result> {
    compute(line: Line2, box: AlignedBox2): DistLine2AlignedBox2Result {
        const result = defaultResult();

        // Translate the line and box so that the box has center at the
        // origin.
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();
        const origin = sub(line.origin, boxCenter);
        const direction = line.direction.clone();

        // The query computes 'result' relative to the box with center at the
        // origin.
        distLine2AlignedBox2DoQuery(origin, direction, boxExtent, result);

        // Translate the closest points to the original coordinates.
        for (let i = 0; i < 2; ++i) {
            result.closest[i] = add(result.closest[i], boxCenter);
        }

        // Compute the distance and squared distance.
        const diff = sub(result.closest[0], result.closest[1]);
        result.sqrDistance = dot(diff, diff);
        result.distance = Math.sqrt(result.sqrDistance);
        return result;
    }
}
