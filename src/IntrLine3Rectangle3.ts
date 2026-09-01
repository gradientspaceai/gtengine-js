// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Rectangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the intersection between a line and a solid rectangle in 3D.
//
// The line is P + t * D, where D is not required to be unit length.
//
// The rectangle has center C, unit-length axis directions W[0] and W[1], and
// extents e[0] and e[1]. A rectangle point is X = C + sum_{i=0}^1 s[i] * W[i]
// where |s[i]| <= e[i] for all i.
//
// The intersection point, if any, is stored in result.point. The
// corresponding line parameter t is stored in result.parameter. The
// corresponding rectangle parameters s[] are stored in result.rectCoord[].
// When the line is in the plane of the rectangle and intersects the
// rectangle, the queries state that there are no intersections.
//
// Upstream TODO: support non-unit-length W[]; return the point or segment of
// intersection when the line is in the plane of the rectangle and intersects
// the rectangle.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrLine3Rectangle3TI (test)
// and IntrLine3Rectangle3FI (find), with result types
// IntrLine3Rectangle3TIResult and IntrLine3Rectangle3FIResult.

import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import type { Line } from './Line';
import type { Rectangle } from './Rectangle';
import { Vector, add, dot, mul, sub } from './Vector';
import { cross } from './Vector3';

// The result of IntrLine3Rectangle3TI queries.
export interface IntrLine3Rectangle3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine3Rectangle3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Rectangle3FI queries.
export interface IntrLine3Rectangle3FIResult {
    intersect: boolean;

    // The line parameter t at the intersection point.
    parameter: number;

    // The rectangle coordinates (s[0],s[1]) of the intersection point.
    // Upstream declares this as 'std::array<T,3>' although only the first two
    // components have meaning for a rectangle; the third is always zero. The
    // port preserves the length-3 array so callers that mirror upstream code
    // behave identically.
    rectCoord: [number, number, number];

    // The intersection point.
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrLine3Rectangle3FIResult {
    return {
        intersect: false,
        parameter: 0,
        rectCoord: [0, 0, 0],
        point: Vector.zero(3)
    };
}

// Test-intersection query for a line and a solid rectangle in 3D.
export class IntrLine3Rectangle3TI implements TIQuery<Line, Rectangle, IntrLine3Rectangle3TIResult> {
    test(line: Line, rectangle: Rectangle): IntrLine3Rectangle3TIResult {
        const result = defaultTIResult();

        // Compute the offset origin and rectangle normal.
        const diff = sub(line.origin, rectangle.center);
        const normal = cross(rectangle.axis[0], rectangle.axis[1]);

        // Solve Q + t*D = s0*W0 + s1*W1 (Q = diff, D = line direction,
        // W0 = edge 0 direction, W1 = edge 1 direction, N = Cross(W0,W1)) by
        //   s0 = Dot(W1,Cross(D,Q)) / Dot(D,N)
        //   s1 = -Dot(W0,Cross(D,Q)) / Dot(D,N)
        //   t = -Dot(Q,N) / Dot(D,N)
        const DdN = dot(line.direction, normal);
        if (DdN === 0) {
            // Line and rectangle are parallel, call it a "no intersection"
            // even if the line and rectangle are coplanar and intersecting.
            result.intersect = false;
            return result;
        }

        const absDdN = Math.abs(DdN);
        const DxQ = cross(line.direction, diff);
        const W1dDxQ = dot(rectangle.axis[1], DxQ);
        if (Math.abs(W1dDxQ) <= rectangle.extent.values[0] * absDdN) {
            const W0dDxQ = dot(rectangle.axis[0], DxQ);
            if (Math.abs(W0dDxQ) <= rectangle.extent.values[1] * absDdN) {
                result.intersect = true;
                return result;
            }
        }
        result.intersect = false;
        return result;
    }
}

// Find-intersection query for a line and a solid rectangle in 3D.
export class IntrLine3Rectangle3FI implements FIQuery<Line, Rectangle, IntrLine3Rectangle3FIResult> {
    find(line: Line, rectangle: Rectangle): IntrLine3Rectangle3FIResult {
        const result = defaultFIResult();

        // Compute the offset origin and rectangle normal.
        const diff = sub(line.origin, rectangle.center);
        const normal = cross(rectangle.axis[0], rectangle.axis[1]);

        // Solve Q + t*D = s0*W0 + s1*W1 (Q = diff, D = line direction,
        // W0 = edge 0 direction, W1 = edge 1 direction, N = Cross(W0,W1)) by
        //   s0 = Dot(W1,Cross(D,Q)) / Dot(D,N)
        //   s1 = -Dot(W0,Cross(D,Q)) / Dot(D,N)
        //   t = -Dot(Q,N) / Dot(D,N)
        const DdN = dot(line.direction, normal);
        if (DdN === 0) {
            // Line and rectangle are parallel, call it a "no intersection"
            // even if the line and rectangle are coplanar and intersecting.
            result.intersect = false;
            return result;
        }

        const absDdN = Math.abs(DdN);
        const DxQ = cross(line.direction, diff);
        const W1dDxQ = dot(rectangle.axis[1], DxQ);
        if (Math.abs(W1dDxQ) <= rectangle.extent.values[0] * absDdN) {
            const W0dDxQ = dot(rectangle.axis[0], DxQ);
            if (Math.abs(W0dDxQ) <= rectangle.extent.values[1] * absDdN) {
                result.intersect = true;
                result.parameter = -dot(diff, normal) / DdN;
                result.rectCoord[0] = W1dDxQ / DdN;
                result.rectCoord[1] = -W0dDxQ / DdN;
                result.point = add(line.origin,
                    mul(result.parameter, line.direction));
                return result;
            }
        }
        result.intersect = false;
        return result;
    }
}
