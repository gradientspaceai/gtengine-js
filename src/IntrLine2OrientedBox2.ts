// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the four
// edges of the box.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream derives
// these queries from the Line2-vs-AlignedBox2 queries solely to reuse the
// protected DoQuery helpers; the derived Result adds no members. In TypeScript
// a subclass cannot narrow the second parameter of the inherited test/find
// (that would be an incompatible override), so the port instead calls the
// exported aligned-box module functions 'intrLine2AlignedBox2TIDoQuery' and
// 'intrLine2AlignedBox2FIDoQuery'. The result types are aliases of the
// aligned-box result types, matching the empty derived Result upstream.

import type { Line } from './Line.js';
import type { OrientedBox } from './OrientedBox.js';
import { Vector, add, sub, mul, dot } from './Vector.js';
import {
    intrLine2AlignedBox2TIDoQuery,
    intrLine2AlignedBox2FIDoQuery,
    defaultIntrLine2AlignedBox2TIResult,
    defaultIntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import type {
    IntrLine2AlignedBox2TIResult,
    IntrLine2AlignedBox2FIResult
} from './IntrLine2AlignedBox2.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The upstream derived Result structs add no members.
export type IntrLine2OrientedBox2TIResult = IntrLine2AlignedBox2TIResult;
export type IntrLine2OrientedBox2FIResult = IntrLine2AlignedBox2FIResult;

// Transform the line to the oriented-box coordinate system.
function transformLine(line: Line, box: OrientedBox):
    { lineOrigin: Vector, lineDirection: Vector } {
    const diff = sub(line.origin, box.center);
    return {
        lineOrigin: Vector.fromArray([
            dot(diff, box.axis[0]),
            dot(diff, box.axis[1])
        ]),
        lineDirection: Vector.fromArray([
            dot(line.direction, box.axis[0]),
            dot(line.direction, box.axis[1])
        ])
    };
}

export class IntrLine2OrientedBox2TI implements
    TIQuery<Line, OrientedBox, IntrLine2OrientedBox2TIResult> {

    test(line: Line, box: OrientedBox): IntrLine2OrientedBox2TIResult {
        const { lineOrigin, lineDirection } = transformLine(line, box);
        const result = defaultIntrLine2AlignedBox2TIResult();
        intrLine2AlignedBox2TIDoQuery(lineOrigin, lineDirection, box.extent,
            result);
        return result;
    }
}

export class IntrLine2OrientedBox2FI implements
    FIQuery<Line, OrientedBox, IntrLine2OrientedBox2FIResult> {

    find(line: Line, box: OrientedBox): IntrLine2OrientedBox2FIResult {
        const { lineOrigin, lineDirection } = transformLine(line, box);
        const result = defaultIntrLine2AlignedBox2FIResult();
        intrLine2AlignedBox2FIDoQuery(lineOrigin, lineDirection, box.extent,
            result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(line.origin,
                mul(result.parameter[i], line.direction));
        }
        return result;
    }
}
