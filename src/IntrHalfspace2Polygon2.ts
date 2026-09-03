// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrHalfspace2Polygon2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the halfspace to be a solid and the polygon to be a
// convex solid.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only an FIQuery specialization for this pair of primitives. The C++
// 'std::vector<Vector2<T>>' polygon becomes a Vector[] (readonly on input).

import { Halfspace } from './Halfspace.js';
import { Vector, add, sub, mul, dot } from './Vector.js';
import type { FIQuery } from './FIQuery.js';

// The result of IntrHalfspace2Polygon2FI.find.
export interface IntrHalfspace2Polygon2FIResult {
    intersect: boolean;

    // If 'intersect' is true, the halfspace and polygon intersect in a convex
    // polygon. NOTE (upstream behavior, preserved): when the polygon lies
    // entirely in the closed halfspace, no clipping is necessary and upstream
    // leaves this array empty rather than copying the input polygon.
    polygon: Vector[];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrHalfspace2Polygon2FIResult {
    return { intersect: false, polygon: [] };
}

export class IntrHalfspace2Polygon2FI implements
    FIQuery<Halfspace, readonly Vector[], IntrHalfspace2Polygon2FIResult> {

    find(halfspace: Halfspace, polygon: readonly Vector[]):
        IntrHalfspace2Polygon2FIResult {
        const result = defaultFIResult();

        // Determine whether the polygon vertices are outside the halfspace,
        // inside the halfspace, or on the halfspace boundary.
        const numVertices = polygon.length;
        const distance = new Array<number>(numVertices).fill(0);
        let positive = 0, negative = 0, positiveIndex = -1;
        for (let i = 0; i < numVertices; ++i) {
            distance[i] = dot(halfspace.normal, polygon[i]) - halfspace.constant;
            if (distance[i] > 0) {
                ++positive;
                if (positiveIndex === -1) {
                    positiveIndex = i;
                }
            } else if (distance[i] < 0) {
                ++negative;
            }
        }

        if (positive === 0) {
            // The polygon is strictly outside the halfspace.
            result.intersect = false;
            return result;
        }

        if (negative === 0) {
            // The polygon is contained in the closed halfspace, so it is
            // fully visible and no clipping is necessary.
            result.intersect = true;
            return result;
        }

        // The line transversely intersects the polygon. Clip the polygon.
        let vertex: Vector;
        let curr: number, prev: number;
        let t: number;

        if (positiveIndex > 0) {
            // Compute the first clip vertex on the line.
            curr = positiveIndex;
            prev = curr - 1;
            t = distance[curr] / (distance[curr] - distance[prev]);
            vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
            result.polygon.push(vertex);

            // Include the vertices on the positive side of the line.
            while (curr < numVertices && distance[curr] > 0) {
                result.polygon.push(polygon[curr++].clone());
            }

            // Compute the last clip vertex on the line.
            if (curr < numVertices) {
                prev = curr - 1;
            } else {
                curr = 0;
                prev = numVertices - 1;
            }
            t = distance[curr] / (distance[curr] - distance[prev]);
            vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
            result.polygon.push(vertex);
        } else {
            // positiveIndex is 0.
            // Include the vertices on the positive side of the line.
            curr = 0;
            while (curr < numVertices && distance[curr] > 0) {
                result.polygon.push(polygon[curr++].clone());
            }

            // Compute the last clip vertex on the line.
            prev = curr - 1;
            t = distance[curr] / (distance[curr] - distance[prev]);
            vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
            result.polygon.push(vertex);

            // Skip the vertices on the negative side of the line.
            while (curr < numVertices && distance[curr] <= 0) {
                curr++;
            }

            // Compute the first clip vertex on the line.
            if (curr < numVertices) {
                prev = curr - 1;
                t = distance[curr] / (distance[curr] - distance[prev]);
                vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
                result.polygon.push(vertex);

                // Keep the vertices on the positive side of the line.
                while (curr < numVertices && distance[curr] > 0) {
                    result.polygon.push(polygon[curr++].clone());
                }
            } else {
                curr = 0;
                prev = numVertices - 1;
                t = distance[curr] / (distance[curr] - distance[prev]);
                vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
                result.polygon.push(vertex);
            }
        }

        result.intersect = true;
        return result;
    }
}
