// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine2AlignedBox2.h
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
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// 'protected void DoQuery(...)' helpers (used by the Ray2/Segment2 versus
// AlignedBox2 queries, which derive from these classes) become the protected
// methods 'doQuery' that mutate the passed-in result, as upstream does. The
// private static 'Clip' becomes the module-private function 'clip', which
// takes the mutable [t0,t1] interval as an object.

import { AlignedBox } from './AlignedBox.js';
import { Line } from './Line.js';
import { Vector, add, sub, mul } from './Vector.js';
import { dotPerp } from './Vector2.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrLine2AlignedBox2TI.test.
export interface IntrLine2AlignedBox2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrLine2AlignedBox2TIResult(): IntrLine2AlignedBox2TIResult {
    return { intersect: false };
}

// The result of IntrLine2AlignedBox2FI.find.
export interface IntrLine2AlignedBox2FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine2AlignedBox2FIResult(): IntrLine2AlignedBox2FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(2), Vector.zero(2)]
    };
}

// Test whether the current clipped segment intersects the current test plane.
// If the return value is true, the segment does intersect the plane and is
// clipped; otherwise, the segment is culled (no intersection with the box).
function clip(denom: number, numer: number, t: { t0: number, t1: number }): boolean {
    if (denom > 0) {
        if (numer > denom * t.t1) {
            return false;
        }
        if (numer > denom * t.t0) {
            t.t0 = numer / denom;
        }
        return true;
    } else if (denom < 0) {
        if (numer > denom * t.t0) {
            return false;
        }
        if (numer > denom * t.t1) {
            t.t1 = numer / denom;
        }
        return true;
    } else {
        return numer <= 0;
    }
}

export class IntrLine2AlignedBox2TI implements
    TIQuery<Line, AlignedBox, IntrLine2AlignedBox2TIResult> {

    test(line: Line, box: AlignedBox): IntrLine2AlignedBox2TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the line to the aligned-box coordinate system.
        const lineOrigin = sub(line.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2TIResult();
        this.doQuery(lineOrigin, line.direction, boxExtent, result);
        return result;
    }

    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        boxExtent: Vector, result: IntrLine2AlignedBox2TIResult): void {
        const LHS = Math.abs(dotPerp(lineDirection, lineOrigin));
        const RHS =
            boxExtent.values[0] * Math.abs(lineDirection.values[1]) +
            boxExtent.values[1] * Math.abs(lineDirection.values[0]);
        result.intersect = (LHS <= RHS);
    }
}

export class IntrLine2AlignedBox2FI implements
    FIQuery<Line, AlignedBox, IntrLine2AlignedBox2FIResult> {

    find(line: Line, box: AlignedBox): IntrLine2AlignedBox2FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector2::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the line to the aligned-box coordinate system.
        const lineOrigin = sub(line.origin, boxCenter);

        const result = defaultIntrLine2AlignedBox2FIResult();
        this.doQuery(lineOrigin, line.direction, boxExtent, result);
        for (let i = 0; i < result.numIntersections; ++i) {
            result.point[i] = add(line.origin,
                mul(result.parameter[i], line.direction));
        }
        return result;
    }

    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        boxExtent: Vector, result: IntrLine2AlignedBox2FIResult): void {
        // The line t-values are in the interval (-infinity,+infinity). Clip
        // the line against all four planes of an aligned box in centered
        // form. The result.numIntersections is
        //   0, no intersection
        //   1, intersect in a single point (t0 is line parameter of point)
        //   2, intersect in a segment (line parameter interval is [t0,t1])
        const t = { t0: -MAX_T, t1: MAX_T };
        const o = lineOrigin.values;
        const d = lineDirection.values;
        const e = boxExtent.values;
        if (clip(+d[0], -o[0] - e[0], t) &&
            clip(-d[0], +o[0] - e[0], t) &&
            clip(+d[1], -o[1] - e[1], t) &&
            clip(-d[1], +o[1] - e[1], t)) {
            result.intersect = true;
            if (t.t1 > t.t0) {
                result.numIntersections = 2;
                result.parameter[0] = t.t0;
                result.parameter[1] = t.t1;
            } else {
                result.numIntersections = 1;
                result.parameter[0] = t.t0;
                result.parameter[1] = t.t0;  // Used by derived classes.
            }
            return;
        }

        result.intersect = false;
        result.numIntersections = 0;
    }
}
