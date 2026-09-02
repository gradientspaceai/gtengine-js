// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Ellipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the ellipsoid to be a solid.
//
// The ellipsoid is (X-C)^T*M*(X-C)-1 = 0. The segment has endpoints P0 and
// P1. The segment origin (center) is P = (P0+P1)/2, the segment direction is
// D = (P1-P0)/|P1-P0| and the segment extent (half the segment length) is
// e = |P1-P0|/2. The segment is X = P+t*D for t in [-e,e]. Substitute the
// segment equation into the ellipsoid equation to obtain a quadratic equation
// Q(t) = a2*t^2 + 2*a1*t + a0 = 0, where a2 = D^T*M*D, a1 = D^T*M*(P-C) and
// a0 = (P-C)^T*M*(P-C)-1. The algorithm involves an analysis of the
// real-valued roots of Q(t) for -e <= t <= e.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// FIQuery<Segment3,Ellipsoid3> from FIQuery<Line3,Ellipsoid3> only to reuse
// the protected DoQuery; the derived Result adds no members, so the result
// type is an alias of the line-ellipsoid result type. The line-ellipsoid
// DoQuery is the exported module function 'intrLine3Ellipsoid3FIDoQuery', and
// the segment-ellipsoid one is exported here as
// 'intrSegment3Ellipsoid3FIDoQuery' for the same reason (the precedent set by
// IntrSegment3Sphere3.ts).
//
// Upstream's TIQuery comment for the coefficients describes a0 as
// '(P0-C)^T*M*(P0-C)-r^2' and a1 as '(P1-P0)^T*M*(P0-C)'; the code actually
// uses the centered form (P is the segment center and the constant term is 1,
// not r^2). The comment above is corrected to match the code. The comment in
// the final branch of the TIQuery mentions '-a1/2' and 'a3*e' where the code
// (correctly) uses the root -a1/a2 and the bound a2*e.

import type { Ellipsoid3 } from './Hyperellipsoid';
import type { FIQuery } from './FIQuery';
import { IntrIntervalsFI } from './IntrIntervals';
import {
    intrLine3Ellipsoid3FIDoQuery,
    defaultIntrLine3Ellipsoid3FIResult
} from './IntrLine3Ellipsoid3';
import type { IntrLine3Ellipsoid3FIResult } from './IntrLine3Ellipsoid3';
import { logAssert } from './Logger';
import { mulMatrix } from './Matrix';
import type { Segment3 } from './Segment';
import { Vector, add, dot, mul, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrSegment3Ellipsoid3TI.test.
export interface IntrSegment3Ellipsoid3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSegment3Ellipsoid3TIResult():
    IntrSegment3Ellipsoid3TIResult {
    return { intersect: false };
}

// The upstream derived FIQuery::Result adds no members.
export type IntrSegment3Ellipsoid3FIResult = IntrLine3Ellipsoid3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Ellipsoid3FIResult():
    IntrSegment3Ellipsoid3FIResult {
    return defaultIntrLine3Ellipsoid3FIResult();
}

// The port of the protected 'FIQuery::DoQuery'. The caller must ensure that
// on entry, 'result' is default constructed as if there is no intersection.
// If an intersection is found, the 'result' values are modified accordingly.
export function intrSegment3Ellipsoid3FIDoQuery(segOrigin: Vector,
    segDirection: Vector, segExtent: number, ellipsoid: Ellipsoid3,
    result: IntrSegment3Ellipsoid3FIResult): void {
    intrLine3Ellipsoid3FIDoQuery(segOrigin, segDirection, ellipsoid, result);

    if (result.intersect) {
        // The line containing the segment intersects the ellipsoid; the
        // t-interval is [t0,t1]. The segment intersects the ellipsoid as long
        // as [t0,t1] overlaps the segment t-interval [-segExtent,+segExtent].
        const segInterval: [number, number] = [-segExtent, segExtent];
        const iiQuery = new IntrIntervalsFI();
        const iiResult = iiQuery.find(result.parameter, segInterval);
        if (iiResult.intersect) {
            result.numIntersections = iiResult.numIntersections;
            result.parameter = [iiResult.overlap[0], iiResult.overlap[1]];
        }
        else {
            // The line containing the segment does not intersect the
            // ellipsoid.
            const empty = defaultIntrSegment3Ellipsoid3FIResult();
            result.intersect = empty.intersect;
            result.numIntersections = empty.numIntersections;
            result.parameter = empty.parameter;
            result.point = empty.point;
        }
    }
}

// Test-intersection query for a segment and a solid ellipsoid in 3D.
export class IntrSegment3Ellipsoid3TI implements
    TIQuery<Segment3, Ellipsoid3, IntrSegment3Ellipsoid3TIResult> {

    test(segment: Segment3, ellipsoid: Ellipsoid3):
        IntrSegment3Ellipsoid3TIResult {
        logAssert(ellipsoid.dimension === 3 && segment.p[0].size === 3,
            'IntrSegment3Ellipsoid3TI: mismatched sizes.');
        const result = defaultIntrSegment3Ellipsoid3TIResult();

        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const M = ellipsoid.getM();
        const diff = sub(segOrigin, ellipsoid.center);
        const matDir = mulMatrix(M, segDirection) as Vector;
        const matDiff = mulMatrix(M, diff) as Vector;
        const a0 = dot(diff, matDiff) - 1;
        const a1 = dot(segDirection, matDiff);
        const a2 = dot(segDirection, matDir);
        const discr = a1 * a1 - a0 * a2;
        if (discr < 0) {
            // Q(t) has no real-valued roots. The segment does not intersect
            // the ellipsoid.
            result.intersect = false;
            return result;
        }

        // Q(-e) = a2*e^2 - 2*a1*e + a0, Q(e) = a2*e^2 + 2*a1*e + a0
        const a2e = a2 * segExtent;
        const tmp0 = a2e * segExtent + a0;  // a2*e^2 + a0
        const tmp1 = 2 * a1 * segExtent;    // 2*a1*e
        const qm = tmp0 - tmp1;  // Q(-e)
        const qp = tmp0 + tmp1;  // Q(e)
        // Q(-e) <= 0 or Q(e) <= 0 means that endpoint is inside or on the
        // solid ellipsoid, so the segment intersects the ellipsoid.
        //
        // Upstream bug (FIXED, the same bug and the same fix as in
        // IntrSegment3Sphere3.ts, reported for B71): upstream tests
        // 'qm * qp <= 0' here, which catches only the case where the two
        // endpoint values have opposite signs (one endpoint inside, one
        // outside). It then concludes, in the comment below, that when Q at
        // the endpoints is negative the segment does not intersect the
        // ellipsoid, and its final test 'qm > 0 && |a1| < a2*e' returns false
        // for that case. That is wrong for a solid ellipsoid: Q(t) < 0 means
        // the segment point is strictly inside, so a segment contained in the
        // ellipsoid does intersect it. As written, the TI query disagreed
        // with the FI query below, which reports the whole segment.
        if (qm <= 0 || qp <= 0) {
            result.intersect = true;
            return result;
        }

        // Both Q(-e) > 0 and Q(e) > 0, so both endpoints are strictly outside
        // the ellipsoid. The minimum of Q(t) occurs at t = -a1/a2. We know
        // that discr >= 0, so Q(t) has a root on (-e,e) when -a1/a2 is in
        // (-e,e). (This is upstream's 'qm > 0 && |a1| < a2*e' with the qm > 0
        // conjunct already established.)
        result.intersect = (Math.abs(a1) < a2e);
        return result;
    }
}

// Find-intersection query for a segment and a solid ellipsoid in 3D.
export class IntrSegment3Ellipsoid3FI implements
    FIQuery<Segment3, Ellipsoid3, IntrSegment3Ellipsoid3FIResult> {

    find(segment: Segment3, ellipsoid: Ellipsoid3):
        IntrSegment3Ellipsoid3FIResult {
        const { center: segOrigin, direction: segDirection,
            extent: segExtent } = segment.getCenteredForm();

        const result = defaultIntrSegment3Ellipsoid3FIResult();
        intrSegment3Ellipsoid3FIDoQuery(segOrigin, segDirection, segExtent,
            ellipsoid, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(segOrigin,
                    mul(result.parameter[i], segDirection));
            }
        }
        return result;
    }
}
