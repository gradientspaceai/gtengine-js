// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3OrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the six
// faces of the box (Liang-Barsky clipping). The queries consider the box to
// be a solid. The algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): upstream derives
// the oriented-box queries from the aligned-box queries only to reuse the
// protected DoQuery members. In TypeScript the derived query cannot keep the
// canonical test()/find() names while changing the second parameter type, so
// the port calls the exported aligned-box module functions
// 'intrLine3AlignedBox3TIDoQuery' and 'intrLine3AlignedBox3FIDoQuery', and the
// public classes implement TIQuery/FIQuery for (Line, OrientedBox) directly.
// The upstream Result structs add no members to the aligned-box results, so
// the port exports type aliases.

import type { Line } from './Line.js';
import type { OrientedBox } from './OrientedBox.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import {
    intrLine3AlignedBox3TIDoQuery,
    intrLine3AlignedBox3FIDoQuery,
    defaultIntrLine3AlignedBox3TIResult,
    defaultIntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3.js';
import type {
    IntrLine3AlignedBox3TIResult,
    IntrLine3AlignedBox3FIResult
} from './IntrLine3AlignedBox3.js';

// The result of IntrLine3OrientedBox3TI.test. Upstream adds no members to
// the aligned-box result.
export type IntrLine3OrientedBox3TIResult = IntrLine3AlignedBox3TIResult;

// The result of IntrLine3OrientedBox3FI.find. Upstream adds no members to
// the aligned-box result.
export type IntrLine3OrientedBox3FIResult = IntrLine3AlignedBox3FIResult;

// Transform the line to the oriented-box coordinate system.
function toBoxCoordinates(line: Line, box: OrientedBox):
    { origin: Vector, direction: Vector } {
    const diff = sub(line.origin, box.center);
    return {
        origin: Vector.fromArray([
            dot(diff, box.axis[0]),
            dot(diff, box.axis[1]),
            dot(diff, box.axis[2])
        ]),
        direction: Vector.fromArray([
            dot(line.direction, box.axis[0]),
            dot(line.direction, box.axis[1]),
            dot(line.direction, box.axis[2])
        ])
    };
}

// Test-intersection query for a line and a solid oriented box in 3D.
export class IntrLine3OrientedBox3TI implements
    TIQuery<Line, OrientedBox, IntrLine3OrientedBox3TIResult> {

    test(line: Line, box: OrientedBox): IntrLine3OrientedBox3TIResult {
        const { origin, direction } = toBoxCoordinates(line, box);
        const result = defaultIntrLine3AlignedBox3TIResult();
        intrLine3AlignedBox3TIDoQuery(origin, direction, box.extent, result);
        return result;
    }
}

// Find-intersection query for a line and a solid oriented box in 3D.
export class IntrLine3OrientedBox3FI implements
    FIQuery<Line, OrientedBox, IntrLine3OrientedBox3FIResult> {

    find(line: Line, box: OrientedBox): IntrLine3OrientedBox3FIResult {
        const { origin, direction } = toBoxCoordinates(line, box);
        const result = defaultIntrLine3AlignedBox3FIResult();
        intrLine3AlignedBox3FIDoQuery(origin, direction, box.extent, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(line.origin,
                    mul(result.parameter[i], line.direction));
            }
        }
        return result;
    }
}
