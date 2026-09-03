// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrRay3Cone3.h
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
// Port notes: upstream derives FIQuery<Real, Ray3, Cone3> from
// FIQuery<Real, Line3, Cone3> and adds no result members, so the port aliases
// IntrLine3Cone3FIResult (the precedent set by IntrSegment2Triangle2.ts) and
// calls the exported 'intrLine3Cone3FIDoQuery' and
// 'intrLine3Cone3ComputePoints'. The inherited protected SetEmpty, SetPoint,
// SetSegment and SetRayPositive helpers are module-private in
// IntrLine3Cone3.ts, so they are replicated here verbatim as module-private
// functions (a port PR may only touch the files of its own batch).

import type { Cone3 } from './Cone.js';
import type { FIQuery } from './FIQuery.js';
import {
    IntrLine3Cone3FIResultType, defaultIntrLine3Cone3FIResult,
    intrLine3Cone3ComputePoints, intrLine3Cone3FIDoQuery
} from './IntrLine3Cone3.js';
import type { IntrLine3Cone3FIResult } from './IntrLine3Cone3.js';
import { logAssert } from './Logger.js';
import { QFNumber } from './QFNumber.js';
import type { Ray3 } from './Ray.js';

// The result of IntrRay3Cone3FI.find. Upstream derives its Result from the
// line-cone Result and adds no members.
export type IntrRay3Cone3FIResult = IntrLine3Cone3FIResult;

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrRay3Cone3FIResult(): IntrRay3Cone3FIResult {
    return defaultIntrLine3Cone3FIResult();
}

// Replicas of the module-private Set* helpers of IntrLine3Cone3.ts.
function setEmpty(result: IntrRay3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isEmpty;
    result.t[0] = new QFNumber();
    result.t[1] = new QFNumber();
}

function setPoint(t: QFNumber, result: IntrRay3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isPoint;
    result.t[0] = t;
    result.t[1] = result.t[0].clone();
}

function setSegment(t0: QFNumber, t1: QFNumber,
    result: IntrRay3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isSegment;
    result.t[0] = t0;
    result.t[1] = t1;
}

function setRayPositive(t: QFNumber, result: IntrRay3Cone3FIResult): void {
    result.type = IntrLine3Cone3FIResultType.isRayPositive;
    result.t[0] = t;
    result.t[1] = new QFNumber(+1, 0, t.d);  // +infinity
}

// The port of 'std::max' for quadratic-field numbers, which returns the
// first argument when the values compare equal (as std::max does).
function qfMax(a: QFNumber, b: QFNumber): QFNumber {
    return a.lessThan(b) ? b : a;
}

// Find-intersection query for a ray and a solid cone in 3D.
export class IntrRay3Cone3FI implements
    FIQuery<Ray3, Cone3, IntrRay3Cone3FIResult> {

    find(ray: Ray3, cone: Cone3): IntrRay3Cone3FIResult {
        logAssert(cone.dimension === 3 && ray.origin.size === 3,
            'IntrRay3Cone3FI: mismatched sizes.');

        // Execute the line-cone query.
        const result = defaultIntrRay3Cone3FIResult();
        intrLine3Cone3FIDoQuery(ray.origin, ray.direction, cone, result);

        // Adjust the t-interval depending on whether the line-cone t-interval
        // overlaps the ray interval [0,+infinity). The block numbers are a
        // continuation of those in IntrLine3Cone3.ts.
        if (result.type !== IntrLine3Cone3FIResultType.isEmpty) {
            const zero = new QFNumber(0, 0, result.t[0].d);

            if (result.type === IntrLine3Cone3FIResultType.isPoint) {
                if (result.t[0].lessThan(zero)) {
                    // Block 12.
                    setEmpty(result);
                }
                // else: Block 13.
            }
            else if (result.type === IntrLine3Cone3FIResultType.isSegment) {
                if (result.t[1].greaterThan(zero)) {
                    // Block 14.
                    setSegment(qfMax(result.t[0], zero), result.t[1], result);
                }
                else if (result.t[1].lessThan(zero)) {
                    // Block 15.
                    setEmpty(result);
                }
                else {  // result.t[1] == zero
                    // Block 16.
                    setPoint(zero, result);
                }
            }
            else if (result.type === IntrLine3Cone3FIResultType.isRayPositive) {
                // Block 17.
                setRayPositive(qfMax(result.t[0], zero), result);
            }
            else {  // result.type == isRayNegative
                if (result.t[1].greaterThan(zero)) {
                    // Block 18.
                    setSegment(zero, result.t[1], result);
                }
                else if (result.t[1].lessThan(zero)) {
                    // Block 19.
                    setEmpty(result);
                }
                else {  // result.t[1] == zero
                    // Block 20.
                    setPoint(zero, result);
                }
            }
        }

        intrLine3Cone3ComputePoints(ray.origin, ray.direction, result);
        result.intersect = (result.type !== IntrLine3Cone3FIResultType.isEmpty);
        return result;
    }
}
