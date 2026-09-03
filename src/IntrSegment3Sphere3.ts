// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the sphere to be a solid.
//
// The sphere is (X-C)^T*(X-C)-r^2 = 0. The segment has endpoints P0 and P1.
// The segment origin (center) is P = (P0+P1)/2, the segment direction is
// D = (P1-P0)/|P1-P0| and the segment extent (half the segment length) is
// e = |P1-P0|/2. The segment is X = P+t*D for t in [-e,e]. Substitute the
// segment equation into the sphere equation to obtain a quadratic equation
// Q(t) = t^2 + 2*a1*t + a0 = 0, where a1 = D^T*(P-C) and
// a0 = (P-C)^T*(P-C)-r^2. The algorithm involves an analysis of the
// real-valued roots of Q(t) for -e <= t <= e.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// FIQuery<Segment3,Sphere3> from FIQuery<Line3,Sphere3> only to reuse the
// protected DoQuery; the derived Result adds no members, so the result type
// is an alias of the line-sphere result type. The protected DoQuery is
// exported as the module function 'intrSegment3Sphere3DoQuery'. The reported
// parameters are relative to the centered form of the segment, C + t * D with
// |t| <= e, as upstream reports them.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Hypersphere } from './Hypersphere.js';
import type { Segment } from './Segment.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import {
    intrLine3Sphere3DoQuery,
    defaultIntrLine3Sphere3FIResult
} from './IntrLine3Sphere3.js';
import type { IntrLine3Sphere3FIResult } from './IntrLine3Sphere3.js';
import { IntrIntervalsFI } from './IntrIntervals.js';

// The result of IntrSegment3Sphere3TI queries.
export interface IntrSegment3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3Sphere3TIResult(): IntrSegment3Sphere3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3Sphere3FIResult = IntrLine3Sphere3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Sphere3FIResult(): IntrSegment3Sphere3FIResult {
    return defaultIntrLine3Sphere3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment3Sphere3DoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, sphere: Hypersphere,
    result: IntrSegment3Sphere3FIResult): void {
    intrLine3Sphere3DoQuery(segOrigin, segDirection, sphere, result);

    if (result.intersect) {
        // The line containing the segment intersects the sphere; the
        // t-interval is [t0,t1]. The segment intersects the sphere as long as
        // [t0,t1] overlaps the segment t-interval [-segExtent,+segExtent].
        const iiQuery = new IntrIntervalsFI();
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the segment does not intersect the sphere.
            const empty = defaultIntrSegment3Sphere3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a segment and a solid sphere in 3D.
export class IntrSegment3Sphere3TI implements
    TIQuery<Segment, Hypersphere, IntrSegment3Sphere3TIResult> {

    test(segment: Segment, sphere: Hypersphere): IntrSegment3Sphere3TIResult {
        const result = defaultIntrSegment3Sphere3TIResult();

        const { center: segOrigin, direction: segDirection, extent: segExtent } =
            segment.getCenteredForm();

        const diff = sub(segOrigin, sphere.center);
        const a0 = dot(diff, diff) - sphere.radius * sphere.radius;
        const a1 = dot(segDirection, diff);
        const discr = a1 * a1 - a0;
        if (discr < 0) {
            // Q(t) has no real-valued roots. The segment does not intersect
            // the sphere.
            result.intersect = false;
            return result;
        }

        // Q(-e) = e^2 - 2*a1*e + a0, Q(e) = e^2 + 2*a1*e + a0
        const tmp0 = segExtent * segExtent + a0;  // e^2 + a0
        const tmp1 = 2 * a1 * segExtent;  // 2*a1*e
        const qm = tmp0 - tmp1;  // Q(-e)
        const qp = tmp0 + tmp1;  // Q(e)
        // Q(-e) <= 0 or Q(e) <= 0 means that endpoint is inside or on the
        // solid sphere, so the segment intersects the sphere.
        //
        // Upstream bug (FIXED; see upstream-bug issue (B71)): upstream tests
        // 'qm * qp <= 0' here, which catches only the case where the two
        // endpoint values have opposite signs (one endpoint inside, one
        // outside). It then concludes, in the comment below, that "When Q at
        // the endpoints is negative, Q(t) < 0 for all t in [-e,e] and the
        // segment does not intersect the sphere", and its final test
        // 'qm > 0 && |a1| < e' returns false for that case. That is wrong for
        // a solid sphere: Q(t) < 0 means the segment point is strictly inside
        // the sphere, so a segment contained in the sphere does intersect it.
        // As written, the TI query disagreed with the FI query below, which
        // reports the whole segment. The fix follows IntrRay3Sphere3TI, which
        // explicitly returns true when the ray origin is inside the sphere:
        // test the endpoints for containment first, then fall through to the
        // upstream test for the remaining case.
        if (qm <= 0 || qp <= 0) {
            result.intersect = true;
            return result;
        }

        // Both Q(-e) > 0 and Q(e) > 0, so both endpoints are strictly outside
        // the sphere. The minimum of Q(t) occurs at t = -a1. We know that
        // discr >= 0, so Q(t) has a root on (-e,e) when -a1 is in (-e,e).
        // (This is upstream's 'qm > 0 && |a1| < e' with the qm > 0 conjunct
        // already established.)
        result.intersect = (Math.abs(a1) < segExtent);
        return result;
    }
}

// Find-intersection query for a segment and a solid sphere in 3D.
export class IntrSegment3Sphere3FI implements
    FIQuery<Segment, Hypersphere, IntrSegment3Sphere3FIResult> {

    find(segment: Segment, sphere: Hypersphere): IntrSegment3Sphere3FIResult {
        const { center: segOrigin, direction: segDirection, extent: segExtent } =
            segment.getCenteredForm();

        const result = defaultIntrSegment3Sphere3FIResult();
        intrSegment3Sphere3DoQuery(segOrigin, segDirection, segExtent, sphere,
            result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(segOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
