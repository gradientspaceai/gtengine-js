// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PolygonWindingOrder.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Determine the winding order of a simple polygon. It is either
// counterclockwise (CCW) or clockwise (CW). If the polygon has one ordering
// but you need the opposite ordering for your application, reverse the array
// of vertices or iterate over it in reverse.
//
// The algorithm finds the lexicographically smallest vertex (the "lower
// left" vertex), which is necessarily a convex vertex of the polygon, and
// classifies the turn there. The sign of DotPerp(next - v, prev - v) is
// positive exactly when the polygon is counterclockwise.
//
// Port notes: the upstream operator() becomes isCounterClockwise(polygon).
// The polygon vertices are Vector objects of size 2; the upstream
// 'polygon[i] < polygon[lowerLeft]' is std::array's lexicographic
// comparison, which is Vector.lessThan in the port.

import { Vector, sub } from './Vector';
import { dotPerp } from './Vector2';

export class PolygonWindingOrder {
    // The polygon vertices must be ordered, either CCW or CW. The function
    // returns true when the ordering is CCW or false when the ordering
    // is CW.
    isCounterClockwise(polygon: Vector[]): boolean {
        const n = polygon.length;
        let lowerLeft = 0;
        for (let i = 1; i < n; ++i) {
            if (polygon[i].lessThan(polygon[lowerLeft])) {
                lowerLeft = i;
            }
        }

        const vLowerLeft = polygon[lowerLeft];
        const vNext = polygon[(lowerLeft + 1) % n];
        const vPrev = polygon[(lowerLeft + n - 1) % n];
        const diffNext = sub(vNext, vLowerLeft);
        const diffPrev = sub(vPrev, vLowerLeft);
        const dp = dotPerp(diffNext, diffPrev);
        return dp > 0;
    }
}
