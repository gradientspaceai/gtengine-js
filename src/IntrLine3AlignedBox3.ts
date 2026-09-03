// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3AlignedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The test-intersection queries use the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The find-intersection queries use parametric clipping against the six faces
// of the box, via Liang-Barsky clipping. The queries consider the box to be a
// solid. The algorithms are described in
// https://www.geometrictools.com/Documentation/IntersectionLineBox.pdf
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// 'protected void DoQuery(...)' helpers (used by the Ray3/Segment3 versus
// AlignedBox3 queries, which derive from these classes) become the protected
// methods 'doQuery' that mutate the passed-in result, as upstream does. The
// private static 'Clip' becomes the module-private function 'clip', which
// takes the mutable [t0,t1] interval as an object.

import { AlignedBox } from './AlignedBox.js';
import { Line } from './Line.js';
import { Vector, add, sub, mul } from './Vector.js';
import { cross } from './Vector3.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The port of std::numeric_limits<T>::max() for T = double.
const MAX_T = Number.MAX_VALUE;

// The result of IntrLine3AlignedBox3TI.test.
export interface IntrLine3AlignedBox3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrLine3AlignedBox3TIResult(): IntrLine3AlignedBox3TIResult {
    return { intersect: false };
}

// The result of IntrLine3AlignedBox3FI.find.
export interface IntrLine3AlignedBox3FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3AlignedBox3FIResult(): IntrLine3AlignedBox3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(3), Vector.zero(3)]
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

export class IntrLine3AlignedBox3TI implements
    TIQuery<Line, AlignedBox, IntrLine3AlignedBox3TIResult> {

    test(line: Line, box: AlignedBox): IntrLine3AlignedBox3TIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the line to the aligned-box coordinate system.
        const lineOrigin = sub(line.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3TIResult();
        this.doQuery(lineOrigin, line.direction, boxExtent, result);
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        boxExtent: Vector, result: IntrLine3AlignedBox3TIResult): void {
        const WxD = cross(lineDirection, lineOrigin).values;
        const absWdU = [
            Math.abs(lineDirection.values[0]),
            Math.abs(lineDirection.values[1]),
            Math.abs(lineDirection.values[2])
        ];
        const e = boxExtent.values;

        if (Math.abs(WxD[0]) > e[1] * absWdU[2] + e[2] * absWdU[1]) {
            return;
        }

        if (Math.abs(WxD[1]) > e[0] * absWdU[2] + e[2] * absWdU[0]) {
            return;
        }

        if (Math.abs(WxD[2]) > e[0] * absWdU[1] + e[1] * absWdU[0]) {
            return;
        }

        result.intersect = true;
    }
}

export class IntrLine3AlignedBox3FI implements
    FIQuery<Line, AlignedBox, IntrLine3AlignedBox3FIResult> {

    find(line: Line, box: AlignedBox): IntrLine3AlignedBox3FIResult {
        // Get the centered form of the aligned box. The axes are implicitly
        // axis[d] = Vector3::Unit(d).
        const { center: boxCenter, extent: boxExtent } = box.getCenteredForm();

        // Transform the line to the aligned-box coordinate system.
        const lineOrigin = sub(line.origin, boxCenter);

        const result = defaultIntrLine3AlignedBox3FIResult();
        this.doQuery(lineOrigin, line.direction, boxExtent, result);
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
    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        boxExtent: Vector, result: IntrLine3AlignedBox3FIResult): void {
        // The line t-values are in the interval (-infinity,+infinity). Clip
        // the line against all six planes of an aligned box in centered form.
        // The result.numIntersections is
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
            clip(-d[1], +o[1] - e[1], t) &&
            clip(+d[2], -o[2] - e[2], t) &&
            clip(-d[2], +o[2] - e[2], t)) {
            result.intersect = true;
            if (t.t1 > t.t0) {
                result.numIntersections = 2;
                result.parameter[0] = t.t0;
                result.parameter[1] = t.t1;
            } else {
                result.numIntersections = 1;
                result.parameter[0] = t.t0;
                result.parameter[1] = t.t0;
            }
        }
    }
}
