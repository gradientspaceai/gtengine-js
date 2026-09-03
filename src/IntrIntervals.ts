// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrIntervals.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The intervals are of the form [t0,t1], [t0,+infinity) or (-infinity,t1].
// Degenerate intervals are allowed (t0 = t1). The queries do not perform
// validation on the input intervals to test whether t0 <= t1.
//
// Port notes (this file is the first Intr* port and sets the precedent for
// the remaining intersection headers):
//
// * Upstream defines two template specializations,
//   'TIQuery<Real, std::array<Real,2>, std::array<Real,2>>' and
//   'FIQuery<Real, std::array<Real,2>, std::array<Real,2>>' (aliased as
//   TIIntervalInterval and FIIntervalInterval). TypeScript has no template
//   specialization, so each becomes a concrete class named after the file
//   with a 'TI'/'FI' suffix: IntrIntervalsTI and IntrIntervalsFI. The nested
//   Result structs become the exported types IntrIntervalsTIResult and
//   IntrIntervalsFIResult.
// * 'operator()' becomes 'test(...)' (TI) and 'find(...)' (FI), per
//   PORTING.md. Upstream overloads 'operator()' four times per class. Only
//   the two-finite-intervals overload matches the TIQuery/FIQuery interface
//   signature, so it keeps the name 'test'/'find'; the other three overloads
//   become explicitly named methods (TypeScript overload resolution on a
//   single name would require an untyped variadic implementation signature):
//     operator()(finite, a, isPositiveInfinite)
//         -> testFiniteSemiInfinite / findFiniteSemiInfinite
//     operator()(a0, isPositiveInfinite0, a1, isPositiveInfinite1)
//         -> testSemiInfiniteSemiInfinite / findSemiInfiniteSemiInfinite
//     operator()(maxTime, interval0, speed0, interval1, speed1)
//         -> testDynamic / findDynamic
// * 'std::array<Real,2>' is 'readonly number[]' for inputs (the queries only
//   read elements 0 and 1) and the tuple '[number, number]' for the result
//   'overlap' field.
// * The 'static int32_t const' members of FIQuery::Result (isEmpty, isPoint,
//   ...) become the exported enum IntrIntervalsFIResultType; a nested enum is
//   exported with a file-qualified name because src/index.ts star-exports
//   every file (see PdeFilterScaleType, BSPrecisionType).
// * Upstream's default Result constructor becomes the module-private factory
//   functions defaultTIResult()/defaultFIResult().

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The kind of intersection set reported by IntrIntervalsFI. Upstream stores
// these as 'static int32_t const' members of FIQuery::Result.
export enum IntrIntervalsFIResultType {
    // No intersection.
    isEmpty = 0,

    // The intervals touch at an endpoint, [t0,t0].
    isPoint = 1,

    // Finite-length interval of intersection, [t0,t1].
    isFinite = 2,

    // Semiinfinite interval of intersection, [t0,+infinity). The overlap[0]
    // is t0 and overlap[1] is +1 as a message that the right endpoint is
    // +infinity (you still need the type to know this interpretation).
    isPositiveInfinite = 3,

    // Semiinfinite interval of intersection, (-infinity,t1]. The overlap[0]
    // is -1 as a message that the left endpoint is -infinity (you still need
    // the type to know this interpretation). The overlap[1] is t1.
    isNegativeInfinite = 4,

    // The dynamic queries all set the type to isDynamicQuery because the
    // queries look for time of first and last contact.
    isDynamicQuery = 5
}

// The result of IntrIntervalsTI queries.
export interface IntrIntervalsTIResult {
    intersect: boolean;

    // Dynamic queries (intervals moving with constant speeds). If 'intersect'
    // is true, the contact times are valid and
    //     0 <= firstTime <= lastTime,  firstTime <= maxTime
    // If 'intersect' is false, there are two cases reported. If the intervals
    // will intersect at firstTime > maxTime, the contact times are reported
    // just as when 'intersect' is true. However, if the intervals will not
    // intersect, then firstTime and lastTime are both set to zero (invalid
    // because 'intersect' is false).
    firstTime: number;
    lastTime: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrIntervalsTIResult {
    return { intersect: false, firstTime: 0, lastTime: 0 };
}

// The result of IntrIntervalsFI queries.
export interface IntrIntervalsFIResult {
    intersect: boolean;

    // Static queries (no motion of intervals over time). The number of
    // intersections is 0 (no overlap), 1 (intervals are just touching), or 2
    // (intervals overlap in an interval). If 'intersect' is false,
    // numIntersections is 0 and 'overlap' is set to [0,0]. If 'intersect' is
    // true, numIntersections is 1 or 2. When 1, 'overlap' is set to [x,x],
    // which is degenerate and represents the single intersection point x.
    // When 2, 'overlap' is the interval of intersection.
    numIntersections: number;
    overlap: [number, number];

    // One of isEmpty, isPoint, isFinite, isPositiveInfinite,
    // isNegativeInfinite or isDynamicQuery.
    type: IntrIntervalsFIResultType;

    // Dynamic queries (intervals moving with constant speeds). If 'intersect'
    // is true, the contact times are valid and
    //     0 <= firstTime <= lastTime,  firstTime <= maxTime
    // If 'intersect' is false, there are two cases reported. If the intervals
    // will intersect at firstTime > maxTime, the contact times are reported
    // just as when 'intersect' is true. However, if the intervals will not
    // intersect, then firstTime and lastTime are both set to zero (invalid
    // because 'intersect' is false).
    firstTime: number;
    lastTime: number;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrIntervalsFIResult {
    return {
        intersect: false,
        numIntersections: 0,
        overlap: [0, 0],
        type: IntrIntervalsFIResultType.isEmpty,
        firstTime: 0,
        lastTime: 0
    };
}

// Test-intersection queries for intervals. The queries test overlap, whether
// a single point or an entire interval.
export class IntrIntervalsTI implements TIQuery<readonly number[], readonly number[], IntrIntervalsTIResult> {
    // Static query for two finite intervals. The firstTime and lastTime
    // values are set to zero, but they are invalid for the static query
    // regardless of the value of 'intersect'.
    test(interval0: readonly number[], interval1: readonly number[]): IntrIntervalsTIResult {
        const result = defaultTIResult();
        result.intersect = (interval0[0] <= interval1[1] && interval0[1] >= interval1[0]);
        return result;
    }

    // Static query where one interval is finite and the other is
    // semiinfinite. The two types of semiinfinite intervals are
    // [a,+infinity), called a positive-infinite interval, and (-infinity,a],
    // called a negative-infinite interval. The firstTime and lastTime values
    // are invalid for the static query regardless of the value of
    // 'intersect'.
    testFiniteSemiInfinite(finite: readonly number[], a: number, isPositiveInfinite: boolean): IntrIntervalsTIResult {
        const result = defaultTIResult();

        if (isPositiveInfinite) {
            result.intersect = (finite[1] >= a);
        }
        else {  // is negative-infinite
            result.intersect = (finite[0] <= a);
        }

        return result;
    }

    // Static query where both intervals are semiinfinite.
    testSemiInfiniteSemiInfinite(a0: number, isPositiveInfinite0: boolean,
        a1: number, isPositiveInfinite1: boolean): IntrIntervalsTIResult {
        const result = defaultTIResult();

        if (isPositiveInfinite0) {
            if (isPositiveInfinite1) {
                result.intersect = true;
            }
            else {  // interval1 is negative-infinite
                result.intersect = (a0 <= a1);
            }
        }
        else {  // interval0 is negative-infinite
            if (isPositiveInfinite1) {
                result.intersect = (a0 >= a1);
            }
            else {  // interval1 is negative-infinite
                result.intersect = true;
            }
        }

        return result;
    }

    // Dynamic query for two finite intervals moving with constant speeds.
    // The current time is 0 and maxTime > 0 is required.
    testDynamic(maxTime: number, interval0: readonly number[], speed0: number,
        interval1: readonly number[], speed1: number): IntrIntervalsTIResult {
        const zero = 0;
        const result = defaultTIResult();

        if (interval0[1] < interval1[0]) {
            // interval0 initially to the left of interval1.
            const diffSpeed = speed0 - speed1;
            if (diffSpeed > zero) {
                // The intervals must move towards each other. 'intersect' is
                // true when the intervals will intersect by maxTime.
                const diffPos = interval1[0] - interval0[1];
                result.intersect = (diffPos <= maxTime * diffSpeed);
                result.firstTime = diffPos / diffSpeed;
                result.lastTime = (interval1[1] - interval0[0]) / diffSpeed;
                return result;
            }
        }
        else if (interval0[0] > interval1[1]) {
            // interval0 initially to the right of interval1.
            const diffSpeed = speed1 - speed0;
            if (diffSpeed > zero) {
                // The intervals must move towards each other. 'intersect' is
                // true when the intervals will intersect by maxTime.
                const diffPos = interval0[0] - interval1[1];
                result.intersect = (diffPos <= maxTime * diffSpeed);
                result.firstTime = diffPos / diffSpeed;
                result.lastTime = (interval0[1] - interval1[0]) / diffSpeed;
                return result;
            }
        }
        else {
            // The intervals are initially intersecting.
            result.intersect = true;
            result.firstTime = zero;
            if (speed1 > speed0) {
                result.lastTime = (interval0[1] - interval1[0]) / (speed1 - speed0);
            }
            else if (speed1 < speed0) {
                result.lastTime = (interval1[1] - interval0[0]) / (speed0 - speed1);
            }
            else {
                result.lastTime = Number.MAX_VALUE;
            }
            return result;
        }

        // The default result has 'intersect' false and the 'firstTime' and
        // 'lastTime' zero.
        return result;
    }
}

// Find-intersection queries for intervals. The queries find the overlap,
// whether a single point or an entire interval.
export class IntrIntervalsFI implements FIQuery<readonly number[], readonly number[], IntrIntervalsFIResult> {
    // Static query for two finite intervals.
    find(interval0: readonly number[], interval1: readonly number[]): IntrIntervalsFIResult {
        const result = defaultFIResult();

        if (interval0[1] < interval1[0] || interval0[0] > interval1[1]) {
            result.numIntersections = 0;
            result.overlap[0] = 0;
            result.overlap[1] = 0;
            result.type = IntrIntervalsFIResultType.isEmpty;
        }
        else if (interval0[1] > interval1[0]) {
            if (interval0[0] < interval1[1]) {
                result.overlap[0] = (interval0[0] < interval1[0] ? interval1[0] : interval0[0]);
                result.overlap[1] = (interval0[1] > interval1[1] ? interval1[1] : interval0[1]);
                if (result.overlap[0] < result.overlap[1]) {
                    result.numIntersections = 2;
                    result.type = IntrIntervalsFIResultType.isFinite;
                }
                else {
                    result.numIntersections = 1;
                    result.type = IntrIntervalsFIResultType.isPoint;
                }
            }
            else {  // interval0[0] == interval1[1]
                result.numIntersections = 1;
                result.overlap[0] = interval0[0];
                result.overlap[1] = result.overlap[0];
                result.type = IntrIntervalsFIResultType.isPoint;
            }
        }
        else {  // interval0[1] == interval1[0]
            result.numIntersections = 1;
            result.overlap[0] = interval0[1];
            result.overlap[1] = result.overlap[0];
            result.type = IntrIntervalsFIResultType.isPoint;
        }

        result.intersect = (result.numIntersections > 0);
        return result;
    }

    // Static query where one interval is finite and the other is
    // semiinfinite. The two types of semiinfinite intervals are
    // [a,+infinity), called a positive-infinite interval, and (-infinity,a],
    // called a negative-infinite interval.
    findFiniteSemiInfinite(finite: readonly number[], a: number, isPositiveInfinite: boolean): IntrIntervalsFIResult {
        const result = defaultFIResult();

        if (isPositiveInfinite) {
            if (finite[1] > a) {
                result.overlap[0] = Math.max(finite[0], a);
                result.overlap[1] = finite[1];
                if (result.overlap[0] < result.overlap[1]) {
                    result.numIntersections = 2;
                    result.type = IntrIntervalsFIResultType.isFinite;
                }
                else {
                    result.numIntersections = 1;
                    result.type = IntrIntervalsFIResultType.isPoint;
                }
            }
            else if (finite[1] === a) {
                result.numIntersections = 1;
                result.overlap[0] = a;
                result.overlap[1] = result.overlap[0];
                result.type = IntrIntervalsFIResultType.isPoint;
            }
            else {
                result.numIntersections = 0;
                result.overlap[0] = 0;
                result.overlap[1] = 0;
                result.type = IntrIntervalsFIResultType.isEmpty;
            }
        }
        else {  // is negative-infinite
            if (finite[0] < a) {
                result.overlap[0] = finite[0];
                result.overlap[1] = Math.min(finite[1], a);
                if (result.overlap[0] < result.overlap[1]) {
                    result.numIntersections = 2;
                    result.type = IntrIntervalsFIResultType.isFinite;
                }
                else {
                    result.numIntersections = 1;
                    result.type = IntrIntervalsFIResultType.isPoint;
                }
            }
            else if (finite[0] === a) {
                result.numIntersections = 1;
                result.overlap[0] = a;
                result.overlap[1] = result.overlap[0];
                result.type = IntrIntervalsFIResultType.isPoint;
            }
            else {
                result.numIntersections = 0;
                result.overlap[0] = 0;
                result.overlap[1] = 0;
                result.type = IntrIntervalsFIResultType.isEmpty;
            }
        }

        result.intersect = (result.numIntersections > 0);
        return result;
    }

    // Static query where both intervals are semiinfinite.
    findSemiInfiniteSemiInfinite(a0: number, isPositiveInfinite0: boolean,
        a1: number, isPositiveInfinite1: boolean): IntrIntervalsFIResult {
        const result = defaultFIResult();

        if (isPositiveInfinite0) {
            if (isPositiveInfinite1) {
                // The overlap[1] is +infinity, but set it to +1 because the
                // upstream Real type might not have a representation for
                // +infinity. The type indicates the interval is
                // positive-infinite, so the +1 is a reminder that overlap[1]
                // is +infinity.
                result.numIntersections = 1;
                result.overlap[0] = Math.max(a0, a1);
                result.overlap[1] = +1;
                result.type = IntrIntervalsFIResultType.isPositiveInfinite;
            }
            else {  // interval1 is negative-infinite
                if (a0 > a1) {
                    result.numIntersections = 0;
                    result.overlap[0] = 0;
                    result.overlap[1] = 0;
                    result.type = IntrIntervalsFIResultType.isEmpty;
                }
                else if (a0 < a1) {
                    result.numIntersections = 2;
                    result.overlap[0] = a0;
                    result.overlap[1] = a1;
                    result.type = IntrIntervalsFIResultType.isFinite;
                }
                else {  // a0 == a1
                    result.numIntersections = 1;
                    result.overlap[0] = a0;
                    result.overlap[1] = result.overlap[0];
                    result.type = IntrIntervalsFIResultType.isPoint;
                }
            }
        }
        else {  // interval0 is negative-infinite
            if (isPositiveInfinite1) {
                if (a0 < a1) {
                    result.numIntersections = 0;
                    result.overlap[0] = 0;
                    result.overlap[1] = 0;
                    result.type = IntrIntervalsFIResultType.isEmpty;
                }
                else if (a0 > a1) {
                    result.numIntersections = 2;
                    result.overlap[0] = a1;
                    result.overlap[1] = a0;
                    result.type = IntrIntervalsFIResultType.isFinite;
                }
                else {
                    result.numIntersections = 1;
                    result.overlap[0] = a1;
                    result.overlap[1] = result.overlap[0];
                    result.type = IntrIntervalsFIResultType.isPoint;
                }
                // Upstream assigns 'result.intersect = (a0 >= a1)' here, but
                // the assignment at the end of the function overwrites it
                // with the equivalent '(numIntersections > 0)'. The dead
                // store is not ported.
            }
            else {  // interval1 is negative-infinite
                // The overlap[0] is -infinity, but set it to -1 because the
                // upstream Real type might not have a representation for
                // -infinity. The type indicates the interval is
                // negative-infinite, so the -1 is a reminder that overlap[0]
                // is -infinity.
                result.numIntersections = 1;
                result.overlap[0] = -1;
                result.overlap[1] = Math.min(a0, a1);
                result.type = IntrIntervalsFIResultType.isNegativeInfinite;
            }
        }

        result.intersect = (result.numIntersections > 0);
        return result;
    }

    // Dynamic query for two finite intervals moving with constant speeds.
    // The current time is 0 and maxTime > 0 is required.
    findDynamic(maxTime: number, interval0: readonly number[], speed0: number,
        interval1: readonly number[], speed1: number): IntrIntervalsFIResult {
        const result = defaultFIResult();
        result.type = IntrIntervalsFIResultType.isDynamicQuery;

        if (interval0[1] < interval1[0]) {
            // interval0 initially to the left of interval1.
            const diffSpeed = speed0 - speed1;
            if (diffSpeed > 0) {
                // The intervals must move towards each other. 'intersect' is
                // true when the intervals will intersect by maxTime.
                const diffPos = interval1[0] - interval0[1];
                result.intersect = (diffPos <= maxTime * diffSpeed);
                result.numIntersections = 1;
                result.firstTime = diffPos / diffSpeed;
                result.lastTime = (interval1[1] - interval0[0]) / diffSpeed;
                // The contact point at firstTime is where the right endpoint
                // of interval0 meets the left endpoint of interval1, that is
                //   interval0[1] + firstTime * speed0
                //     = interval1[0] + firstTime * speed1
                // Upstream computes 'interval0[0] + firstTime * speed0',
                // which is the moved *left* endpoint of interval0 and is not
                // the contact point unless interval0 is degenerate. This is
                // an upstream bug; the port uses the contact point, which is
                // also the mirror image of the interval0-on-the-right case
                // below.
                result.overlap[0] = interval0[1] + result.firstTime * speed0;
                result.overlap[1] = result.overlap[0];
                return result;
            }
        }
        else if (interval0[0] > interval1[1]) {
            // interval0 initially to the right of interval1.
            const diffSpeed = speed1 - speed0;
            if (diffSpeed > 0) {
                // The intervals must move towards each other. 'intersect' is
                // true when the intervals will intersect by maxTime.
                const diffPos = interval0[0] - interval1[1];
                result.intersect = (diffPos <= maxTime * diffSpeed);
                result.numIntersections = 1;
                result.firstTime = diffPos / diffSpeed;
                result.lastTime = (interval0[1] - interval1[0]) / diffSpeed;
                result.overlap[0] = interval1[1] + result.firstTime * speed1;
                result.overlap[1] = result.overlap[0];
                return result;
            }
        }
        else {
            // The intervals are initially intersecting.
            result.intersect = true;
            result.firstTime = 0;
            if (speed1 > speed0) {
                result.lastTime = (interval0[1] - interval1[0]) / (speed1 - speed0);
            }
            else if (speed1 < speed0) {
                result.lastTime = (interval1[1] - interval0[0]) / (speed0 - speed1);
            }
            else {
                result.lastTime = Number.MAX_VALUE;
            }

            if (interval0[1] > interval1[0]) {
                if (interval0[0] < interval1[1]) {
                    result.numIntersections = 2;
                    result.overlap[0] = (interval0[0] < interval1[0] ? interval1[0] : interval0[0]);
                    result.overlap[1] = (interval0[1] > interval1[1] ? interval1[1] : interval0[1]);
                }
                else {  // interval0[0] == interval1[1]
                    result.numIntersections = 1;
                    result.overlap[0] = interval0[0];
                    result.overlap[1] = result.overlap[0];
                }
            }
            else {  // interval0[1] == interval1[0]
                result.numIntersections = 1;
                result.overlap[0] = interval0[1];
                result.overlap[1] = result.overlap[0];
            }
            return result;
        }

        // The default result has the correct state for no-intersection.
        return result;
    }
}
