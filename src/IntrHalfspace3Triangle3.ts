// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Queries for intersection of objects with halfspaces. These are useful for
// containment testing, object culling, and clipping.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent.

import { Halfspace } from './Halfspace';
import { Triangle } from './Triangle';
import { Vector, add, sub, mul, dot } from './Vector';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';

// The result of IntrHalfspace3Triangle3TI.test.
export interface IntrHalfspace3Triangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrHalfspace3Triangle3TIResult {
    return { intersect: false };
}

// The result of IntrHalfspace3Triangle3FI.find.
export interface IntrHalfspace3Triangle3FIResult {
    intersect: boolean;

    // The triangle is clipped against the plane defining the halfspace. The
    // 'numPoints' is either 0 (no intersection), 1 (point), 2 (segment),
    // 3 (triangle), or 4 (quadrilateral).
    numPoints: number;
    point: [Vector, Vector, Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrHalfspace3Triangle3FIResult {
    return {
        intersect: false,
        numPoints: 0,
        point: [Vector.zero(3), Vector.zero(3), Vector.zero(3), Vector.zero(3)]
    };
}

export class IntrHalfspace3Triangle3TI implements
    TIQuery<Halfspace, Triangle, IntrHalfspace3Triangle3TIResult> {

    test(halfspace: Halfspace, triangle: Triangle): IntrHalfspace3Triangle3TIResult {
        const result = defaultTIResult();

        // Project the triangle vertices onto the normal line. The plane of
        // the halfspace occurs at the origin (zero) of the normal line.
        const s: number[] = [0, 0, 0];
        for (let i = 0; i < 3; ++i) {
            s[i] = dot(halfspace.normal, triangle.v[i]) - halfspace.constant;
        }

        // The triangle and halfspace intersect when the projection interval
        // maximum is nonnegative.
        result.intersect = (Math.max(Math.max(s[0], s[1]), s[2]) >= 0);
        return result;
    }
}

export class IntrHalfspace3Triangle3FI implements
    FIQuery<Halfspace, Triangle, IntrHalfspace3Triangle3FIResult> {

    find(halfspace: Halfspace, triangle: Triangle): IntrHalfspace3Triangle3FIResult {
        const result = defaultFIResult();

        // Determine on which side of the plane the vertices lie. The table of
        // possibilities is listed next with n = numNegative, p = numPositive,
        // and z = numZero.
        //
        //   n p z  intersection
        //   ---------------------------------
        //   0 3 0  triangle (original)
        //   0 2 1  triangle (original)
        //   0 1 2  triangle (original)
        //   0 0 3  triangle (original)
        //   1 2 0  quad (2 edges clipped)
        //   1 1 1  triangle (1 edge clipped)
        //   1 0 2  edge
        //   2 1 0  triangle (2 edges clipped)
        //   2 0 1  vertex
        //   3 0 0  none

        const s: number[] = [0, 0, 0];
        let numPositive = 0, numNegative = 0;
        for (let i = 0; i < 3; ++i) {
            s[i] = dot(halfspace.normal, triangle.v[i]) - halfspace.constant;
            if (s[i] > 0) {
                ++numPositive;
            } else if (s[i] < 0) {
                ++numNegative;
            }
        }

        if (numNegative === 0) {
            // The triangle is in the halfspace.
            result.intersect = true;
            result.numPoints = 3;
            result.point[0] = triangle.v[0].clone();
            result.point[1] = triangle.v[1].clone();
            result.point[2] = triangle.v[2].clone();
        } else if (numNegative === 1) {
            result.intersect = true;
            if (numPositive === 2) {
                // The portion of the triangle in the halfspace is a
                // quadrilateral.
                result.numPoints = 4;
                for (let i0 = 0; i0 < 3; ++i0) {
                    if (s[i0] < 0) {
                        const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                        result.point[0] = triangle.v[i1].clone();
                        result.point[1] = triangle.v[i2].clone();
                        const t2 = s[i2] / (s[i2] - s[i0]);
                        const t0 = s[i0] / (s[i0] - s[i1]);
                        result.point[2] = add(triangle.v[i2],
                            mul(t2, sub(triangle.v[i0], triangle.v[i2])));
                        result.point[3] = add(triangle.v[i0],
                            mul(t0, sub(triangle.v[i1], triangle.v[i0])));
                        break;
                    }
                }
            } else if (numPositive === 1) {
                // The portion of the triangle in the halfspace is a triangle
                // with one vertex on the plane.
                result.numPoints = 3;
                for (let i0 = 0; i0 < 3; ++i0) {
                    if (s[i0] === 0) {
                        const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                        result.point[0] = triangle.v[i0].clone();
                        const t1 = s[i1] / (s[i1] - s[i2]);
                        const p = add(triangle.v[i1],
                            mul(t1, sub(triangle.v[i2], triangle.v[i1])));
                        if (s[i1] > 0) {
                            result.point[1] = triangle.v[i1].clone();
                            result.point[2] = p;
                        } else {
                            result.point[1] = p;
                            result.point[2] = triangle.v[i2].clone();
                        }
                        break;
                    }
                }
            } else {
                // Only an edge of the triangle is in the halfspace.
                result.numPoints = 0;
                for (let i = 0; i < 3; ++i) {
                    if (s[i] === 0) {
                        result.point[result.numPoints++] = triangle.v[i].clone();
                    }
                }
            }
        } else if (numNegative === 2) {
            result.intersect = true;
            if (numPositive === 1) {
                // The portion of the triangle in the halfspace is a triangle.
                result.numPoints = 3;
                for (let i0 = 0; i0 < 3; ++i0) {
                    if (s[i0] > 0) {
                        const i1 = (i0 + 1) % 3, i2 = (i0 + 2) % 3;
                        result.point[0] = triangle.v[i0].clone();
                        const t0 = s[i0] / (s[i0] - s[i1]);
                        const t2 = s[i2] / (s[i2] - s[i0]);
                        result.point[1] = add(triangle.v[i0],
                            mul(t0, sub(triangle.v[i1], triangle.v[i0])));
                        result.point[2] = add(triangle.v[i2],
                            mul(t2, sub(triangle.v[i0], triangle.v[i2])));
                        break;
                    }
                }
            } else {
                // Only a vertex of the triangle is in the halfspace.
                result.numPoints = 1;
                for (let i = 0; i < 3; ++i) {
                    if (s[i] === 0) {
                        result.point[0] = triangle.v[i].clone();
                        break;
                    }
                }
            }
        } else {
            // The triangle is outside the halfspace (numNegative == 3).
            result.intersect = false;
            result.numPoints = 0;
        }

        return result;
    }
}
