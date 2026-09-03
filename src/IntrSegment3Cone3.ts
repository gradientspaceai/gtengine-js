// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSegment3Cone3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the cone to be single sided and solid. The cone height
// range is [hmin,hmax]. The cone can be infinite where hmin = 0 and
// hmax = +infinity, infinite truncated where hmin > 0 and hmax = +infinity,
// finite where hmin = 0 and hmax < +infinity, or a cone frustum where
// hmin > 0 and hmax < +infinity. The algorithm details are found in
// https://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
//
// Port notes: upstream derives FIQuery<T, Segment3, Cone3> from
// FIQuery<T, Line3, Cone3> and adds no result members, so the port aliases
// IntrLine3Cone3FIResult (the precedent set by IntrSegment2Triangle2.ts) and
// calls the exported 'intrLine3Cone3FIDoQuery' and
// 'intrLine3Cone3ComputePoints'. The inherited protected SetEmpty, SetPoint
// and SetSegment helpers are module-private in IntrLine3Cone3.ts, so they are
// replicated here verbatim as module-private functions (a port PR may only
// touch the files of its own batch).
//
// The query parameterizes the segment as p[0] + t * (p[1] - p[0]) for t in
// [0,1], so the reported result.t values are relative to that (non-unit)
// direction.

import type { Cone3 } from './Cone';
import type { FIQuery } from './FIQuery';
import {
    IntrLine3Cone3FIResultType, defaultIntrLine3Cone3FIResult,
    intrLine3Cone3ComputePoints, intrLine3Cone3FIDoQuery
} from './IntrLine3Cone3';
import type { IntrLine3Cone3FIResult } from './IntrLine3Cone3';
import { logAssert } from './Logger';
import { QFNumber } from './QFNumber';
import type { Segment3 } from './Segment';
import { sub } from './Vector';

// The result of IntrSegment3Cone3FI.find. Upstream derives its Result from
// the line-cone Result and adds no members.
export type IntrSegment3Cone3FIResult = IntrLine3Cone3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSegment3Cone3FIResult(): IntrSegment3Cone3FIResult {
    return defaultIntrLine3Cone3FIResult();
}

// Replicas of the module-private Set* helpers of IntrLine3Cone3.ts.
function setEmpty(result: IntrSegment3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isEmpty;
    result.t[0] = new QFNumber();
    result.t[1] = new QFNumber();
}

function setPoint(t: QFNumber, result: IntrSegment3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isPoint;
    result.t[0] = t;
    result.t[1] = result.t[0].clone();
}

function setSegment(t0: QFNumber, t1: QFNumber,
    result: IntrSegment3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isSegment;
    result.t[0] = t0;
    result.t[1] = t1;
}

// The ports of 'std::max' and 'std::min' for quadratic-field numbers, which
// return the first argument when the values compare equal (as the standard
// library functions do).
function qfMax(a: QFNumber, b: QFNumber): QFNumber {
    return a.lessThan(b) ? b : a;
}

function qfMin(a: QFNumber, b: QFNumber): QFNumber {
    return b.lessThan(a) ? b : a;
}

// Find-intersection query for a segment and a solid cone in 3D.
export class IntrSegment3Cone3FI implements
    FIQuery<Segment3, Cone3, IntrSegment3Cone3FIResult> {

    find(segment: Segment3, cone: Cone3): IntrSegment3Cone3FIResult {
        logAssert(cone.dimension === 3 && segment.dimension === 3,
            'IntrSegment3Cone3FI: mismatched sizes.');

        // Execute the line-cone query.
        const result = defaultIntrSegment3Cone3FIResult();
        const segOrigin = segment.p[0];
        const segDirection = sub(segment.p[1], segment.p[0]);
        intrLine3Cone3FIDoQuery(segOrigin, segDirection, cone, result);

        // Adjust the t-interval depending on whether the line-cone t-interval
        // overlaps the segment interval [0,1]. The block numbers are a
        // continuation of those in IntrRay3Cone3.ts, which themselves are a
        // continuation of those in IntrLine3Cone3.ts.
        if (result.type !== IntrLine3Cone3FIResultType.isEmpty) {
            const zero = new QFNumber(0, 0, result.t[0].d);
            const one = new QFNumber(1, 0, result.t[0].d);

            if (result.type === IntrLine3Cone3FIResultType.isPoint) {
                if (result.t[0].lessThan(zero) || result.t[0].greaterThan(one)) {
                    // Block 21.
                    setEmpty(result);
                }
                // else: Block 22.
            }
            else if (result.type === IntrLine3Cone3FIResultType.isSegment) {
                if (result.t[1].lessThan(zero) || result.t[0].greaterThan(one)) {
                    // Block 23.
                    setEmpty(result);
                }
                else {
                    const t0 = qfMax(zero, result.t[0]);
                    const t1 = qfMin(one, result.t[1]);
                    if (t0.lessThan(t1)) {
                        // Block 24.
                        setSegment(t0, t1, result);
                    }
                    else {
                        // Block 25.
                        setPoint(t0, result);
                    }
                }
            }
            else if (result.type === IntrLine3Cone3FIResultType.isRayPositive) {
                if (one.lessThan(result.t[0])) {
                    // Block 26.
                    setEmpty(result);
                }
                else if (one.greaterThan(result.t[0])) {
                    // Block 27.
                    setSegment(qfMax(zero, result.t[0]), one, result);
                }
                else {
                    // Block 28.
                    setPoint(one, result);
                }
            }
            else {  // result.type == isRayNegative
                if (zero.greaterThan(result.t[1])) {
                    // Block 29.
                    setEmpty(result);
                }
                else if (zero.lessThan(result.t[1])) {
                    // Block 30.
                    setSegment(zero, qfMin(one, result.t[1]), result);
                }
                else {
                    // Block 31.
                    setPoint(zero, result);
                }
            }
        }

        intrLine3Cone3ComputePoints(segment.p[0], segDirection, result);
        result.intersect = (result.type !== IntrLine3Cone3FIResultType.isEmpty);
        return result;
    }
}
