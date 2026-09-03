// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox2Circle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The find-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionMovingCircleRectangle.pdf
//
// Port notes: the FIQuery is a dynamic (moving-objects) query whose
// 'operator()' takes four arguments, so it does not implement the two-argument
// FIQuery interface; the method is named 'find' with the four arguments, as
// in IntrSphere3Triangle3.ts. Upstream's 'intersectionType' integer becomes
// the exported enum IntrAlignedBox2Circle2FIResultType. The protected DoQuery
// and the private case handlers become module-private functions that write
// into the result object (upstream passes 'Result&').

import type { AlignedBox } from './AlignedBox.js';
import { DistPointAlignedBox } from './DistPointAlignedBox.js';
import type { Hypersphere } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, sub } from './Vector.js';
import { perp } from './Vector2.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrAlignedBox2Circle2TI.test.
export interface IntrAlignedBox2Circle2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrAlignedBox2Circle2TIResult():
    IntrAlignedBox2Circle2TIResult {
    return { intersect: false };
}

// The intersection query considers the box and circle to be solids; that is,
// the circle object includes the region inside the circular boundary and the
// box object includes the region inside the rectangular boundary. If the
// circle object and box object overlap, the objects intersect.
export class IntrAlignedBox2Circle2TI implements
    TIQuery<AlignedBox, Hypersphere, IntrAlignedBox2Circle2TIResult> {

    test(box: AlignedBox, circle: Hypersphere): IntrAlignedBox2Circle2TIResult {
        logAssert(box.dimension === 2 && circle.center.size === 2,
            'IntrAlignedBox2Circle2TI: mismatched sizes.');
        const pbQuery = new DistPointAlignedBox();
        const pbResult = pbQuery.compute(circle.center, box);
        return {
            intersect: (pbResult.sqrDistance <= circle.radius * circle.radius)
        };
    }
}

// The kind of contact reported by IntrAlignedBox2Circle2FI.
export enum IntrAlignedBox2Circle2FIResultType {
    // The objects are initially overlapping. The contactPoint is only one of
    // infinitely many points in the overlap. The contactTime is 0.
    initiallyOverlapping = -1,

    // The objects are initially separated and do not intersect later. The
    // contactTime and contactPoint are invalid (both zero).
    noContact = 0,

    // The objects are initially separated but intersect later. The
    // contactTime is the first time T > 0 of contact and contactPoint is the
    // corresponding first point of contact.
    contact = 1
}

// The result of IntrAlignedBox2Circle2FI.find.
export interface IntrAlignedBox2Circle2FIResult {
    intersectionType: IntrAlignedBox2Circle2FIResultType;
    contactTime: number;
    contactPoint: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrAlignedBox2Circle2FIResult():
    IntrAlignedBox2Circle2FIResult {
    return {
        intersectionType: IntrAlignedBox2Circle2FIResultType.noContact,
        contactTime: 0,
        contactPoint: Vector.zero(2)
    };
}

// Currently, only a dynamic query is supported. A static query will need to
// compute the intersection set of (solid) box and circle.
export class IntrAlignedBox2Circle2FI {
    find(box: AlignedBox, boxVelocity: Vector, circle: Hypersphere,
        circleVelocity: Vector): IntrAlignedBox2Circle2FIResult {
        logAssert(box.dimension === 2 && boxVelocity.size === 2
            && circle.center.size === 2 && circleVelocity.size === 2,
            'IntrAlignedBox2Circle2FI: mismatched sizes.');
        const result = defaultIntrAlignedBox2Circle2FIResult();

        // Translate the circle and box so that the box center becomes the
        // origin. Compute the velocity of the circle relative to the box.
        const { center: boxCenter, extent } = box.getCenteredForm();
        const C = sub(circle.center, boxCenter);
        const V = sub(circleVelocity, boxVelocity);

        // Change signs on components, if necessary, to transform C to the
        // first quadrant. Adjust the velocity accordingly.
        const sign: [number, number] = [0, 0];
        for (let i = 0; i < 2; ++i) {
            if (C.values[i] >= 0) {
                sign[i] = 1;
            }
            else {
                C.values[i] = -C.values[i];
                V.values[i] = -V.values[i];
                sign[i] = -1;
            }
        }

        doQuery(extent, C, circle.radius, V, result);

        if (result.intersectionType !== IntrAlignedBox2Circle2FIResultType.noContact) {
            // Translate back to the original coordinate system.
            for (let i = 0; i < 2; ++i) {
                if (sign[i] < 0) {
                    result.contactPoint.values[i] = -result.contactPoint.values[i];
                }
            }

            result.contactPoint = add(result.contactPoint, boxCenter);
        }
        return result;
    }
}

// The query assumes the box is axis-aligned with center at the origin and
// with extent K, and that the circle center C is in the first quadrant. The
// Voronoi region of C relative to the rounded box determines which case
// handler applies.
function doQuery(K: Vector, C: Vector, radius: number, V: Vector,
    result: IntrAlignedBox2Circle2FIResult): void {
    const delta = sub(C, K);
    if (delta.values[1] <= radius) {
        if (delta.values[0] <= radius) {
            if (delta.values[1] <= 0) {
                if (delta.values[0] <= 0) {
                    interiorOverlap(C, result);
                }
                else {
                    edgeOverlap(0, 1, K, C, delta, radius, result);
                }
            }
            else {
                if (delta.values[0] <= 0) {
                    edgeOverlap(1, 0, K, C, delta, radius, result);
                }
                else {
                    if (dot(delta, delta) <= radius * radius) {
                        vertexOverlap(K, delta, radius, result);
                    }
                    else {
                        vertexSeparated(K, delta, V, radius, result);
                    }
                }
            }
        }
        else {
            edgeUnbounded(0, 1, K, C, radius, delta, V, result);
        }
    }
    else {
        if (delta.values[0] <= radius) {
            edgeUnbounded(1, 0, K, C, radius, delta, V, result);
        }
        else {
            vertexUnbounded(K, C, radius, delta, V, result);
        }
    }
}

function interiorOverlap(C: Vector,
    result: IntrAlignedBox2Circle2FIResult): void {
    result.intersectionType = IntrAlignedBox2Circle2FIResultType.initiallyOverlapping;
    result.contactTime = 0;
    result.contactPoint = C.clone();
}

function edgeOverlap(i0: number, i1: number, K: Vector, C: Vector,
    delta: Vector, radius: number,
    result: IntrAlignedBox2Circle2FIResult): void {
    result.intersectionType = (delta.values[i0] < radius
        ? IntrAlignedBox2Circle2FIResultType.initiallyOverlapping
        : IntrAlignedBox2Circle2FIResultType.contact);
    result.contactTime = 0;
    result.contactPoint.values[i0] = K.values[i0];
    result.contactPoint.values[i1] = C.values[i1];
}

function vertexOverlap(K0: Vector, delta: Vector, radius: number,
    result: IntrAlignedBox2Circle2FIResult): void {
    const sqrDistance = delta.values[0] * delta.values[0]
        + delta.values[1] * delta.values[1];
    const sqrRadius = radius * radius;
    result.intersectionType = (sqrDistance < sqrRadius
        ? IntrAlignedBox2Circle2FIResultType.initiallyOverlapping
        : IntrAlignedBox2Circle2FIResultType.contact);
    result.contactTime = 0;
    result.contactPoint = K0.clone();
}

function vertexSeparated(K0: Vector, delta0: Vector, V: Vector, radius: number,
    result: IntrAlignedBox2Circle2FIResult): void {
    const q0 = -dot(V, delta0);
    if (q0 > 0) {
        const dotVPerpD0 = dot(V, perp(delta0));
        const q2 = dot(V, V);
        const q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
        if (q1 >= 0) {
            intersectsVertex(0, 1, K0, q0, q1, q2, result);
        }
    }
}

function edgeUnbounded(i0: number, i1: number, K0: Vector, C: Vector,
    radius: number, delta0: Vector, V: Vector,
    result: IntrAlignedBox2Circle2FIResult): void {
    if (V.values[i0] < 0) {
        const dotVPerpD0 = V.values[i0] * delta0.values[i1]
            - V.values[i1] * delta0.values[i0];
        if (radius * V.values[i1] + dotVPerpD0 >= 0) {
            const K1 = Vector.zero(2);
            const delta1 = Vector.zero(2);
            K1.values[i0] = K0.values[i0];
            K1.values[i1] = -K0.values[i1];
            delta1.values[i0] = C.values[i0] - K1.values[i0];
            delta1.values[i1] = C.values[i1] - K1.values[i1];
            const dotVPerpD1 = V.values[i0] * delta1.values[i1]
                - V.values[i1] * delta1.values[i0];
            if (radius * V.values[i1] + dotVPerpD1 <= 0) {
                intersectsEdge(i0, i1, K0, C, radius, V, result);
            }
            else {
                const q2 = dot(V, V);
                const q1 = radius * radius * q2 - dotVPerpD1 * dotVPerpD1;
                if (q1 >= 0) {
                    const q0 = -(V.values[i0] * delta1.values[i0]
                        + V.values[i1] * delta1.values[i1]);
                    intersectsVertex(i0, i1, K1, q0, q1, q2, result);
                }
            }
        }
        else {
            const q2 = dot(V, V);
            const q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
            if (q1 >= 0) {
                const q0 = -(V.values[i0] * delta0.values[i0]
                    + V.values[i1] * delta0.values[i1]);
                intersectsVertex(i0, i1, K0, q0, q1, q2, result);
            }
        }
    }
}

function vertexUnbounded(K0: Vector, C: Vector, radius: number, delta0: Vector,
    V: Vector, result: IntrAlignedBox2Circle2FIResult): void {
    if (V.values[0] < 0 && V.values[1] < 0) {
        const dotVPerpD0 = dot(V, perp(delta0));
        if (radius * V.values[0] - dotVPerpD0 <= 0) {
            if (-radius * V.values[1] - dotVPerpD0 >= 0) {
                const q2 = dot(V, V);
                const q1 = radius * radius * q2 - dotVPerpD0 * dotVPerpD0;
                const q0 = -dot(V, delta0);
                intersectsVertex(0, 1, K0, q0, q1, q2, result);
            }
            else {
                const K1 = Vector.fromArray([K0.values[0], -K0.values[1]]);
                const delta1 = sub(C, K1);
                const dotVPerpD1 = dot(V, perp(delta1));
                if (-radius * V.values[1] - dotVPerpD1 >= 0) {
                    intersectsEdge(0, 1, K0, C, radius, V, result);
                }
                else {
                    const q2 = dot(V, V);
                    const q1 = radius * radius * q2 - dotVPerpD1 * dotVPerpD1;
                    if (q1 >= 0) {
                        const q0 = -dot(V, delta1);
                        intersectsVertex(0, 1, K1, q0, q1, q2, result);
                    }
                }
            }
        }
        else {
            const K2 = Vector.fromArray([-K0.values[0], K0.values[1]]);
            const delta2 = sub(C, K2);
            const dotVPerpD2 = dot(V, perp(delta2));
            if (radius * V.values[0] - dotVPerpD2 <= 0) {
                intersectsEdge(1, 0, K0, C, radius, V, result);
            }
            else {
                const q2 = dot(V, V);
                const q1 = radius * radius * q2 - dotVPerpD2 * dotVPerpD2;
                if (q1 >= 0) {
                    const q0 = -dot(V, delta2);
                    intersectsVertex(1, 0, K2, q0, q1, q2, result);
                }
            }
        }
    }
}

function intersectsVertex(i0: number, i1: number, K: Vector, q0: number,
    q1: number, q2: number,
    result: IntrAlignedBox2Circle2FIResult): void {
    result.intersectionType = IntrAlignedBox2Circle2FIResultType.contact;
    result.contactTime = (q0 - Math.sqrt(q1)) / q2;
    result.contactPoint.values[i0] = K.values[i0];
    result.contactPoint.values[i1] = K.values[i1];
}

function intersectsEdge(i0: number, i1: number, K0: Vector, C: Vector,
    radius: number, V: Vector,
    result: IntrAlignedBox2Circle2FIResult): void {
    result.intersectionType = IntrAlignedBox2Circle2FIResultType.contact;
    result.contactTime = (K0.values[i0] + radius - C.values[i0]) / V.values[i0];
    result.contactPoint.values[i0] = K0.values[i0];
    result.contactPoint.values[i1] = C.values[i1]
        + result.contactTime * V.values[i1];
}
