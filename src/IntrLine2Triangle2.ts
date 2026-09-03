// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2Triangle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the triangle to be a solid. The algorithms are based
// on determining on which side of the line the vertices lie. The test uses
// the sign of the projections of the vertices onto a normal line that is
// perpendicular to the specified line. The table of possibilities is listed
// next with n = numNegative, p = numPositive and z = numZero.
//
//   n p z  intersection
//   ------------------------------------
//   0 3 0  none
//   0 2 1  vertex
//   0 1 2  edge
//   0 0 3  none (degenerate triangle)
//   1 2 0  segment (2 edges clipped)
//   1 1 1  segment (1 edge clipped)
//   1 0 2  edge
//   2 1 0  segment (2 edges clipped)
//   2 0 1  vertex
//   3 0 0  none
//
// The case (n,p,z) = (0,0,3) is treated as a no-intersection because the
// triangle is degenerate.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// 'protected void DoQuery(...)' helper (used by the Ray2/Segment2 versus
// Triangle2 queries, which derive from this class) becomes the protected
// method 'doQuery' that mutates the passed-in result, as upstream does.

import { Line } from './Line.js';
import { Triangle } from './Triangle.js';
import { Vector, add, sub, mul, dot } from './Vector.js';
import { dotPerp } from './Vector2.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The result of IntrLine2Triangle2TI.test.
export interface IntrLine2Triangle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrLine2Triangle2TIResult(): IntrLine2Triangle2TIResult {
    return { intersect: false };
}

// The result of IntrLine2Triangle2FI.find.
export interface IntrLine2Triangle2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine2Triangle2FIResult(): IntrLine2Triangle2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

export class IntrLine2Triangle2TI implements
    TIQuery<Line, Triangle, IntrLine2Triangle2TIResult> {

    // The line is P + t * D, where P is a point on the line and D is a
    // direction vector that does not have to be unit length. This is useful
    // when using a 2-point representation P0 + t * (P1 - P0).
    test(line: Line, triangle: Triangle): IntrLine2Triangle2TIResult {
        const result = defaultIntrLine2Triangle2TIResult();

        let numPositive = 0, numNegative = 0, numZero = 0;
        for (let i = 0; i < 3; ++i) {
            const diff = sub(triangle.v[i], line.origin);
            const s = dotPerp(line.direction, diff);
            if (s > 0) {
                ++numPositive;
            } else if (s < 0) {
                ++numNegative;
            } else {
                ++numZero;
            }
        }

        result.intersect =
            (numZero === 0 && numPositive > 0 && numNegative > 0) ||
            (numZero === 1) ||
            (numZero === 2);

        return result;
    }
}

export class IntrLine2Triangle2FI implements
    FIQuery<Line, Triangle, IntrLine2Triangle2FIResult> {

    // The line is P + t * D, where P is a point on the line and D is a
    // direction vector that does not have to be unit length. This is useful
    // when using a 2-point representation P0 + t * (P1 - P0).
    find(line: Line, triangle: Triangle): IntrLine2Triangle2FIResult {
        const result = defaultIntrLine2Triangle2FIResult();
        this.doQuery(line.origin, line.direction, triangle, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(line.origin,
                    mul(result.parameter[i], line.direction));
            }
        }
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(origin: Vector, direction: Vector, triangle: Triangle,
        result: IntrLine2Triangle2FIResult): void {
        const s: number[] = [0, 0, 0];
        let numPositive = 0, numNegative = 0, numZero = 0;
        for (let i = 0; i < 3; ++i) {
            const diff = sub(triangle.v[i], origin);
            s[i] = dotPerp(direction, diff);
            if (s[i] > 0) {
                ++numPositive;
            } else if (s[i] < 0) {
                ++numNegative;
            } else {
                ++numZero;
            }
        }

        if (numZero === 0 && numPositive > 0 && numNegative > 0) {
            // (n,p,z) is (1,2,0) or (2,1,0).
            result.intersect = true;
            result.numIntersections = 2;

            // sign is +1 when (n,p) is (2,1) or -1 when (n,p) is (1,2).
            const sign = (3 > numPositive * 2 ? 1 : -1);
            for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                if (sign * s[i2] > 0) {
                    const diffVi0P0 = sub(triangle.v[i0], origin);
                    const diffVi2Vi0 = sub(triangle.v[i2], triangle.v[i0]);
                    const lambda0 = s[i0] / (s[i0] - s[i2]);
                    const q0 = add(diffVi0P0, mul(lambda0, diffVi2Vi0));
                    result.parameter[0] = dot(direction, q0);
                    const diffVi1P0 = sub(triangle.v[i1], origin);
                    const diffVi2Vi1 = sub(triangle.v[i2], triangle.v[i1]);
                    const lambda1 = s[i1] / (s[i1] - s[i2]);
                    const q1 = add(diffVi1P0, mul(lambda1, diffVi2Vi1));
                    result.parameter[1] = dot(direction, q1);
                    break;
                }
            }
        } else if (numZero === 1) {
            // (n,p,z) is (1,1,1), (2,0,1) or (0,2,1).
            result.intersect = true;
            for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                if (s[i2] === 0) {
                    const diffVi2P0 = sub(triangle.v[i2], origin);
                    result.parameter[0] = dot(direction, diffVi2P0);
                    if (numPositive === 2 || numNegative === 2) {
                        // (n,p,z) is (2,0,1) or (0,2,1).
                        result.numIntersections = 1;
                        result.parameter[1] = result.parameter[0];
                    } else {
                        // (n,p,z) is (1,1,1).
                        result.numIntersections = 2;
                        const diffVi0P0 = sub(triangle.v[i0], origin);
                        const diffVi1Vi0 = sub(triangle.v[i1], triangle.v[i0]);
                        const lambda0 = s[i0] / (s[i0] - s[i1]);
                        const q = add(diffVi0P0, mul(lambda0, diffVi1Vi0));
                        result.parameter[1] = dot(direction, q);
                    }
                    break;
                }
            }
        } else if (numZero === 2) {
            // (n,p,z) is (1,0,2) or (0,1,2).
            result.intersect = true;
            result.numIntersections = 2;
            for (let i0 = 1, i1 = 2, i2 = 0; i2 < 3; i0 = i1, i1 = i2++) {
                if (s[i2] !== 0) {
                    const diffVi0P0 = sub(triangle.v[i0], origin);
                    result.parameter[0] = dot(direction, diffVi0P0);
                    const diffVi1P0 = sub(triangle.v[i1], origin);
                    result.parameter[1] = dot(direction, diffVi1P0);
                    break;
                }
            }
        }
        // else: (n,p,z) is (3,0,0), (0,3,0) or (0,0,3). The default result
        // has all members zero, so no additional assignments are needed.

        if (result.intersect) {
            const directionSqrLength = dot(direction, direction);
            result.parameter[0] /= directionSqrLength;
            result.parameter[1] /= directionSqrLength;
            if (result.parameter[0] > result.parameter[1]) {
                const tmp = result.parameter[0];
                result.parameter[0] = result.parameter[1];
                result.parameter[1] = tmp;
            }
        }
    }
}
