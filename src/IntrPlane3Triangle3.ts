// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrPlane3Triangle3TI
// and IntrPlane3Triangle3FI.

import type { FIQuery } from './FIQuery';
import type { Plane3 } from './Hyperplane';
import { logAssert } from './Logger';
import type { Triangle3 } from './Triangle';
import { Vector, add, dot, mul, sub } from './Vector';
import type { TIQuery } from './TIQuery';

// The result of IntrPlane3Triangle3TI.test.
export interface IntrPlane3Triangle3TIResult {
    intersect: boolean;

    // The number is 0 (no intersection), 1 (plane and triangle intersect at a
    // single point [vertex]), 2 (plane and triangle intersect in a segment),
    // or 3 (triangle is in the plane). When the number is 2, the segment is
    // either interior to the triangle or is an edge of the triangle, the
    // distinction stored in 'isInterior'.
    numIntersections: number;
    isInterior: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrPlane3Triangle3TIResult {
    return { intersect: false, numIntersections: 0, isInterior: false };
}

// The result of IntrPlane3Triangle3FI.find.
export interface IntrPlane3Triangle3FIResult {
    intersect: boolean;

    // See the comment for IntrPlane3Triangle3TIResult.numIntersections.
    numIntersections: number;
    isInterior: boolean;
    point: [Vector, Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrPlane3Triangle3FIResult():
    IntrPlane3Triangle3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        isInterior: false,
        point: [Vector.zero(3), Vector.zero(3), Vector.zero(3)]
    };
}

// Determine on which side of the plane the vertices lie. The table of
// possibilities is listed next with n = numNegative, p = numPositive, and
// z = numZero.
//
//   n p z  intersection
//   ------------------------------------
//   0 3 0  none
//   0 2 1  vertex
//   0 1 2  edge
//   0 0 3  triangle in the plane
//   1 2 0  segment (2 edges clipped)
//   1 1 1  segment (1 edge clipped)
//   1 0 2  edge
//   2 1 0  segment (2 edges clipped)
//   2 0 1  vertex
//   3 0 0  none
function classifyVertices(plane: Plane3, triangle: Triangle3):
    { s: [number, number, number], numPositive: number, numNegative: number,
        numZero: number } {
    const s: [number, number, number] = [0, 0, 0];
    let numPositive = 0, numNegative = 0, numZero = 0;
    for (let i = 0; i < 3; ++i) {
        s[i] = dot(plane.normal, triangle.v[i]) - plane.constant;
        if (s[i] > 0) {
            ++numPositive;
        }
        else if (s[i] < 0) {
            ++numNegative;
        }
        else {
            ++numZero;
        }
    }
    return { s, numPositive, numNegative, numZero };
}

// Test-intersection query for a plane and a triangle in 3D.
export class IntrPlane3Triangle3TI implements
    TIQuery<Plane3, Triangle3, IntrPlane3Triangle3TIResult> {

    test(plane: Plane3, triangle: Triangle3): IntrPlane3Triangle3TIResult {
        logAssert(plane.dimension === 3 && triangle.v[0].size === 3,
            'IntrPlane3Triangle3TI: mismatched sizes.');

        const result = defaultTIResult();
        const { s, numPositive, numNegative, numZero } =
            classifyVertices(plane, triangle);

        if (numZero === 0 && numPositive > 0 && numNegative > 0) {
            result.intersect = true;
            result.numIntersections = 2;
            result.isInterior = true;
            return result;
        }

        if (numZero === 1) {
            result.intersect = true;
            for (let i = 0; i < 3; ++i) {
                if (s[i] === 0) {
                    if (numPositive === 2 || numNegative === 2) {
                        result.numIntersections = 1;
                    }
                    else {
                        result.numIntersections = 2;
                        result.isInterior = true;
                    }
                    break;
                }
            }
            return result;
        }

        if (numZero === 2) {
            result.intersect = true;
            result.numIntersections = 2;
            result.isInterior = false;
            return result;
        }

        if (numZero === 3) {
            result.intersect = true;
            result.numIntersections = 3;
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }
        return result;
    }
}

// Find-intersection query for a plane and a triangle in 3D.
export class IntrPlane3Triangle3FI implements
    FIQuery<Plane3, Triangle3, IntrPlane3Triangle3FIResult> {

    find(plane: Plane3, triangle: Triangle3): IntrPlane3Triangle3FIResult {
        logAssert(plane.dimension === 3 && triangle.v[0].size === 3,
            'IntrPlane3Triangle3FI: mismatched sizes.');

        const result = defaultIntrPlane3Triangle3FIResult();
        const { s, numPositive, numNegative, numZero } =
            classifyVertices(plane, triangle);

        if (numZero === 0 && numPositive > 0 && numNegative > 0) {
            result.intersect = true;
            result.numIntersections = 2;
            result.isInterior = true;
            const sign = 3 - numPositive * 2;
            for (let i0 = 0; i0 < 3; ++i0) {
                if (sign * s[i0] > 0) {
                    const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                    const t1 = s[i1] / (s[i1] - s[i0]);
                    const t2 = s[i2] / (s[i2] - s[i0]);
                    result.point[0] = add(triangle.v[i1],
                        mul(t1, sub(triangle.v[i0], triangle.v[i1])));
                    result.point[1] = add(triangle.v[i2],
                        mul(t2, sub(triangle.v[i0], triangle.v[i2])));
                    break;
                }
            }
            return result;
        }

        if (numZero === 1) {
            result.intersect = true;
            for (let i0 = 0; i0 < 3; ++i0) {
                if (s[i0] === 0) {
                    const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                    result.point[0] = triangle.v[i0].clone();
                    if (numPositive === 2 || numNegative === 2) {
                        result.numIntersections = 1;
                    }
                    else {
                        result.numIntersections = 2;
                        result.isInterior = true;
                        const t = s[i1] / (s[i1] - s[i2]);
                        result.point[1] = add(triangle.v[i1],
                            mul(t, sub(triangle.v[i2], triangle.v[i1])));
                    }
                    break;
                }
            }
            return result;
        }

        if (numZero === 2) {
            result.intersect = true;
            result.numIntersections = 2;
            result.isInterior = false;
            for (let i0 = 0; i0 < 3; ++i0) {
                if (s[i0] !== 0) {
                    const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                    result.point[0] = triangle.v[i1].clone();
                    result.point[1] = triangle.v[i2].clone();
                    break;
                }
            }
            return result;
        }

        if (numZero === 3) {
            result.intersect = true;
            result.numIntersections = 3;
            result.point[0] = triangle.v[0].clone();
            result.point[1] = triangle.v[1].clone();
            result.point[2] = triangle.v[2].clone();
        }
        else {
            result.intersect = false;
            result.numIntersections = 0;
        }
        return result;
    }
}
