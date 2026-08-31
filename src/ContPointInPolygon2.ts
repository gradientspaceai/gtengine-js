// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContPointInPolygon2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Given a polygon as an ordered list of vertices (x[i],y[i]) for
// 0 <= i < N and a test point (xt,yt), the queries return 'true' if
// (xt,yt) is in the polygon and 'false' if it is not. All queries require
// that the number of vertices satisfies N >= 3.
//
// Port notes:
// - The (numPoints, points) pair becomes a single array; the array length is
//   the vertex count.
// - Upstream's comment claims the class stores a copy of 'points', but the
//   member is 'Vector2<Real> const*', a pointer to the caller's data. The
//   port stores the array by reference, which has the same aliasing
//   behavior, and keeps the caveat about the caller's data staying alive
//   (and unmodified) for the lifetime of the object.

import { logAssert } from './Logger';
import { Vector } from './Vector';

export class PointInPolygon2 {
    // The port keeps a reference to the caller's array (see the port notes),
    // so be careful about the persistence of 'points'.
    private readonly mPoints: readonly Vector[];

    constructor(points: readonly Vector[]) {
        logAssert(points.length >= 3,
            'PointInPolygon2: at least 3 vertices are required.');
        for (const point of points) {
            logAssert(point.size === 2,
                'PointInPolygon2: the vertices must be 2D.');
        }
        this.mPoints = points;
    }

    // The number of polygon vertices.
    get numPoints(): number {
        return this.mPoints.length;
    }

    // Simple polygons (ray-intersection counting). The polygon may be convex
    // or concave, and the vertices may be in either orientation.
    contains(p: Vector): boolean {
        logAssert(p.size === 2, 'PointInPolygon2: the point must be 2D.');

        const numPoints = this.mPoints.length;
        let inside = false;
        for (let i = 0, j = numPoints - 1; i < numPoints; j = i++) {
            const u0 = this.mPoints[i].values;
            const u1 = this.mPoints[j].values;
            let lhs: number, rhs: number;

            if (p.values[1] < u1[1]) {
                // U1 is above the ray.
                if (u0[1] <= p.values[1]) {
                    // U0 is on or below the ray.
                    lhs = (p.values[1] - u0[1]) * (u1[0] - u0[0]);
                    rhs = (p.values[0] - u0[0]) * (u1[1] - u0[1]);
                    if (lhs > rhs) {
                        inside = !inside;
                    }
                }
            } else if (p.values[1] < u0[1]) {
                // U1 is on or below the ray, U0 is above the ray.
                lhs = (p.values[1] - u0[1]) * (u1[0] - u0[0]);
                rhs = (p.values[0] - u0[0]) * (u1[1] - u0[1]);
                if (lhs < rhs) {
                    inside = !inside;
                }
            }
        }
        return inside;
    }

    // Algorithms for convex polygons. The input polygons must have vertices
    // in counterclockwise order.

    // O(N) algorithm (which-side-of-edge tests).
    containsConvexOrderN(p: Vector): boolean {
        logAssert(p.size === 2, 'PointInPolygon2: the point must be 2D.');

        const numPoints = this.mPoints.length;
        for (let i1 = 0, i0 = numPoints - 1; i1 < numPoints; i0 = i1++) {
            const v0 = this.mPoints[i0].values;
            const v1 = this.mPoints[i1].values;
            const nx = v1[1] - v0[1];
            const ny = v0[0] - v1[0];
            const dx = p.values[0] - v0[0];
            const dy = p.values[1] - v0[1];
            if (nx * dx + ny * dy > 0) {
                return false;
            }
        }
        return true;
    }

    // O(log N) algorithm (bisection and recursion, like a BSP tree).
    containsConvexOrderLogN(p: Vector): boolean {
        logAssert(p.size === 2, 'PointInPolygon2: the point must be 2D.');
        return this.subContainsPoint(p, 0, 0);
    }

    // The polygon must have exactly four vertices. This method is like the
    // O(log N) one and uses three which-side-of-segment tests instead of
    // four which-side-of-edge tests. If the polygon does not have four
    // vertices, the function returns false.
    containsQuadrilateral(p: Vector): boolean {
        logAssert(p.size === 2, 'PointInPolygon2: the point must be 2D.');

        if (this.mPoints.length !== 4) {
            return false;
        }

        const v0 = this.mPoints[0].values;
        const v1 = this.mPoints[1].values;
        const v2 = this.mPoints[2].values;
        const v3 = this.mPoints[3].values;

        let nx = v2[1] - v0[1];
        let ny = v0[0] - v2[0];
        let dx = p.values[0] - v0[0];
        let dy = p.values[1] - v0[1];

        if (nx * dx + ny * dy > 0) {
            // P is potentially in <V0,V1,V2>.
            nx = v1[1] - v0[1];
            ny = v0[0] - v1[0];
            if (nx * dx + ny * dy > 0) {
                return false;
            }

            nx = v2[1] - v1[1];
            ny = v1[0] - v2[0];
            dx = p.values[0] - v1[0];
            dy = p.values[1] - v1[1];
            if (nx * dx + ny * dy > 0) {
                return false;
            }
        } else {
            // P is potentially in <V0,V2,V3>.
            nx = v0[1] - v3[1];
            ny = v3[0] - v0[0];
            if (nx * dx + ny * dy > 0) {
                return false;
            }

            nx = v3[1] - v2[1];
            ny = v2[0] - v3[0];
            dx = p.values[0] - v3[0];
            dy = p.values[1] - v3[1];
            if (nx * dx + ny * dy > 0) {
                return false;
            }
        }
        return true;
    }

    // For recursion in containsConvexOrderLogN.
    private subContainsPoint(p: Vector, i0: number, i1: number): boolean {
        const numPoints = this.mPoints.length;
        let nx: number, ny: number, dx: number, dy: number;

        const diff = i1 - i0;
        if (diff === 1 || (diff < 0 && diff + numPoints === 1)) {
            const v0 = this.mPoints[i0].values;
            const v1 = this.mPoints[i1].values;
            nx = v1[1] - v0[1];
            ny = v0[0] - v1[0];
            dx = p.values[0] - v0[0];
            dy = p.values[1] - v0[1];
            return nx * dx + ny * dy <= 0;
        }

        // Bisect the index range.
        let mid: number;
        if (i0 < i1) {
            mid = (i0 + i1) >> 1;
        } else {
            mid = (i0 + i1 + numPoints) >> 1;
            if (mid >= numPoints) {
                mid -= numPoints;
            }
        }

        // Determine which side of the splitting line contains the point.
        const vi0 = this.mPoints[i0].values;
        const vmid = this.mPoints[mid].values;
        nx = vmid[1] - vi0[1];
        ny = vi0[0] - vmid[0];
        dx = p.values[0] - vi0[0];
        dy = p.values[1] - vi0[1];
        if (nx * dx + ny * dy > 0) {
            // P is potentially in <V(i0),V(i0+1),...,V(mid-1),V(mid)>.
            return this.subContainsPoint(p, i0, mid);
        } else {
            // P is potentially in <V(mid),V(mid+1),...,V(i1-1),V(i1)>.
            return this.subContainsPoint(p, mid, i1);
        }
    }
}
