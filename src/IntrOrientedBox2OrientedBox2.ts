// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrOrientedBox2OrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The test-intersection query uses the method of separating axes.
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
// The set of potential separating directions includes the 2 edge normals of
// box0 and the 2 edge normals of box1. The integer 'separating' identifies
// the axis that reported separation; there may be more than one but only one
// is reported. The value is 0 when box0.axis[0] separates, 1 when
// box0.axis[1] separates, 2 when box1.axis[0] separates or 3 when
// box1.axis[1] separates.
//
// Port notes (see IntrIntervals.ts for the Intr* precedent): the two upstream
// template specializations become the classes IntrOrientedBox2OrientedBox2TI
// (test) and IntrOrientedBox2OrientedBox2FI (find), with the result types
// IntrOrientedBox2OrientedBox2TIResult and
// IntrOrientedBox2OrientedBox2FIResult. The private 'Outside' clipping helper
// becomes a module-private function; because it both clips the polygon and
// reports the outside state, it returns the pair { outside, polygon } rather
// than mutating a reference argument.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { OrientedBox } from './OrientedBox.js';
import { Vector, add, dot, mul, negate, sub } from './Vector.js';

// The result of IntrOrientedBox2OrientedBox2TI queries.
export interface IntrOrientedBox2OrientedBox2TIResult {
    intersect: boolean;

    // The index of the separating axis, valid only when 'intersect' is false.
    separating: number;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrOrientedBox2OrientedBox2TIResult {
    return { intersect: false, separating: 0 };
}

// The result of IntrOrientedBox2OrientedBox2FI queries.
export interface IntrOrientedBox2OrientedBox2FIResult {
    intersect: boolean;

    // If 'intersect' is true, the boxes intersect in a convex 'polygon'.
    polygon: Vector[];
}

// The port of the upstream FIQuery::Result default constructor.
function defaultFIResult(): IntrOrientedBox2OrientedBox2FIResult {
    return { intersect: false, polygon: [] };
}

// The line normals are inner pointing. The function reports 'outside' true
// when the incoming polygon is outside the line, in which case the boxes do
// not intersect. When 'outside' is false, the returned polygon is the
// incoming polygon intersected with the closed halfspace defined by the line.
function outside(origin: Vector, normal: Vector, polygon: Vector[]):
    { outside: boolean, polygon: Vector[] } {
    // Determine whether the polygon vertices are outside the line, inside the
    // line, or on the line.
    const numVertices = polygon.length;
    const distance = new Array<number>(numVertices);
    let positive = 0, negative = 0, positiveIndex = -1;
    for (let i = 0; i < numVertices; ++i) {
        distance[i] = dot(normal, sub(polygon[i], origin));
        if (distance[i] > 0) {
            ++positive;
            if (positiveIndex === -1) {
                positiveIndex = i;
            }
        }
        else if (distance[i] < 0) {
            ++negative;
        }
    }

    if (positive === 0) {
        // The polygon is strictly outside the line.
        return { outside: true, polygon };
    }

    if (negative === 0) {
        // The polygon is contained in the closed halfspace whose boundary is
        // the line. It is fully visible and no clipping is necessary.
        return { outside: false, polygon };
    }

    // The line transversely intersects the polygon. Clip the polygon.
    const clipPolygon: Vector[] = [];
    let vertex: Vector;
    let curr: number, prev: number;
    let t: number;

    if (positiveIndex > 0) {
        // Compute the first clip vertex on the line.
        curr = positiveIndex;
        prev = curr - 1;
        t = distance[curr] / (distance[curr] - distance[prev]);
        vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
        clipPolygon.push(vertex);

        // Include the vertices on the positive side of the line.
        while (curr < numVertices && distance[curr] > 0) {
            clipPolygon.push(polygon[curr++]);
        }

        // Compute the last clip vertex on the line.
        if (curr < numVertices) {
            prev = curr - 1;
        }
        else {
            curr = 0;
            prev = numVertices - 1;
        }
        t = distance[curr] / (distance[curr] - distance[prev]);
        vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
        clipPolygon.push(vertex);
    }
    else {  // positiveIndex is 0
        // Include the vertices on the positive side of the line.
        curr = 0;
        while (curr < numVertices && distance[curr] > 0) {
            clipPolygon.push(polygon[curr++]);
        }

        // Compute the last clip vertex on the line.
        prev = curr - 1;
        t = distance[curr] / (distance[curr] - distance[prev]);
        vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
        clipPolygon.push(vertex);

        // Skip the vertices on the negative side of the line.
        while (curr < numVertices && distance[curr] <= 0) {
            curr++;
        }

        // Compute the first clip vertex on the line.
        if (curr < numVertices) {
            prev = curr - 1;
            t = distance[curr] / (distance[curr] - distance[prev]);
            vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
            clipPolygon.push(vertex);

            // Keep the vertices on the positive side of the line.
            while (curr < numVertices && distance[curr] > 0) {
                clipPolygon.push(polygon[curr++]);
            }
        }
        else {
            curr = 0;
            prev = numVertices - 1;
            t = distance[curr] / (distance[curr] - distance[prev]);
            vertex = add(polygon[curr], mul(t, sub(polygon[prev], polygon[curr])));
            clipPolygon.push(vertex);
        }
    }

    return { outside: false, polygon: clipPolygon };
}

// Test-intersection query for two solid oriented boxes in 2D.
export class IntrOrientedBox2OrientedBox2TI implements TIQuery<OrientedBox, OrientedBox, IntrOrientedBox2OrientedBox2TIResult> {
    test(box0: OrientedBox, box1: OrientedBox): IntrOrientedBox2OrientedBox2TIResult {
        const result = defaultTIResult();

        // Convenience variables.
        const A0 = box0.axis;
        const A1 = box1.axis;
        const E0 = box0.extent.values;
        const E1 = box1.extent.values;

        // Compute difference of box centers, D = C1-C0.
        const D = sub(box1.center, box0.center);

        const absA0dA1 = [[0, 0], [0, 0]];
        let rSum: number;

        // Test axis box0.axis[0].
        absA0dA1[0][0] = Math.abs(dot(A0[0], A1[0]));
        absA0dA1[0][1] = Math.abs(dot(A0[0], A1[1]));
        rSum = E0[0] + E1[0] * absA0dA1[0][0] + E1[1] * absA0dA1[0][1];
        if (Math.abs(dot(A0[0], D)) > rSum) {
            result.intersect = false;
            result.separating = 0;
            return result;
        }

        // Test axis box0.axis[1].
        absA0dA1[1][0] = Math.abs(dot(A0[1], A1[0]));
        absA0dA1[1][1] = Math.abs(dot(A0[1], A1[1]));
        rSum = E0[1] + E1[0] * absA0dA1[1][0] + E1[1] * absA0dA1[1][1];
        if (Math.abs(dot(A0[1], D)) > rSum) {
            result.intersect = false;
            result.separating = 1;
            return result;
        }

        // Test axis box1.axis[0].
        rSum = E1[0] + E0[0] * absA0dA1[0][0] + E0[1] * absA0dA1[1][0];
        if (Math.abs(dot(A1[0], D)) > rSum) {
            result.intersect = false;
            result.separating = 2;
            return result;
        }

        // Test axis box1.axis[1].
        rSum = E1[1] + E0[0] * absA0dA1[0][1] + E0[1] * absA0dA1[1][1];
        if (Math.abs(dot(A1[1], D)) > rSum) {
            result.intersect = false;
            result.separating = 3;
            return result;
        }

        result.intersect = true;
        return result;
    }
}

// Find-intersection query for two solid oriented boxes in 2D. The
// intersection is a convex polygon computed by clipping box0 against the
// four edge lines of box1.
export class IntrOrientedBox2OrientedBox2FI implements FIQuery<OrientedBox, OrientedBox, IntrOrientedBox2OrientedBox2FIResult> {
    find(box0: OrientedBox, box1: OrientedBox): IntrOrientedBox2OrientedBox2FIResult {
        const result = defaultFIResult();
        result.intersect = true;

        // Initialize the intersection polygon to box0, listing the vertices
        // in counterclockwise order.
        let vertex = box0.getVertices();
        result.polygon.push(vertex[0]);  // C - e0 * U0 - e1 * U1
        result.polygon.push(vertex[1]);  // C + e0 * U0 - e1 * U1
        result.polygon.push(vertex[3]);  // C + e0 * U0 + e1 * U1
        result.polygon.push(vertex[2]);  // C - e0 * U0 + e1 * U1

        // Clip the polygon using the lines defining the edges of box1. The
        // line normal points inside box1. The line origin is the first vertex
        // of the edge when traversing box1 counterclockwise.
        vertex = box1.getVertices();
        const normal = [
            box1.axis[1], negate(box1.axis[0]), box1.axis[0],
            negate(box1.axis[1])
        ];

        for (let i = 0; i < 4; ++i) {
            const clipped = outside(vertex[i], normal[i], result.polygon);
            if (clipped.outside) {
                // The boxes are separated.
                result.intersect = false;
                result.polygon = [];
                break;
            }
            result.polygon = clipped.polygon;
        }

        return result;
    }
}
