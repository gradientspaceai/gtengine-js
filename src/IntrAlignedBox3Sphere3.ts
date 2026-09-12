// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The find-intersection query is based on the document
// https://www.geometrictools.com/Documentation/IntersectionMovingSphereBox.pdf
// and also uses the method of separating axes,
// https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf
//
// Port notes: the FIQuery is a dynamic (moving-objects) query whose
// 'operator()' takes four arguments, so it does not implement the
// two-argument FIQuery interface; the method is named 'find' with the four
// arguments, as in IntrSphere3Triangle3.ts. Upstream's 'intersectionType'
// integer becomes the exported enum IntrAlignedBox3Sphere3FIResultType. The
// protected DoQuery and the private case handlers become module-private
// functions that write into the result object (upstream passes 'Result&').

import { AlignedBox } from './AlignedBox.js';
import { DistPointAlignedBox } from './DistPointAlignedBox.js';
import type { Hypersphere } from './Hypersphere.js';
import { IntrRay3AlignedBox3TI } from './IntrRay3AlignedBox3.js';
import { Ray } from './Ray.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrAlignedBox3Sphere3TI.test.
export interface IntrAlignedBox3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrAlignedBox3Sphere3TIResult():
    IntrAlignedBox3Sphere3TIResult {
    return { intersect: false };
}

// The intersection query considers the box and sphere to be solids; that is,
// the sphere object includes the region inside the spherical boundary and the
// box object includes the region inside the cuboid boundary. If the sphere
// object and box object overlap, the objects intersect.
export class IntrAlignedBox3Sphere3TI implements
    TIQuery<AlignedBox, Hypersphere, IntrAlignedBox3Sphere3TIResult> {

    test(box: AlignedBox, sphere: Hypersphere): IntrAlignedBox3Sphere3TIResult {
        logAssert(box.dimension === 3 && sphere.center.size === 3,
            'IntrAlignedBox3Sphere3TI: mismatched sizes.');
        const pbQuery = new DistPointAlignedBox();
        const pbResult = pbQuery.compute(sphere.center, box);
        return {
            intersect: (pbResult.sqrDistance <= sphere.radius * sphere.radius)
        };
    }
}

// The kind of contact reported by IntrAlignedBox3Sphere3FI.
export enum IntrAlignedBox3Sphere3FIResultType {
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

// The result of IntrAlignedBox3Sphere3FI.find.
export interface IntrAlignedBox3Sphere3FIResult {
    intersectionType: IntrAlignedBox3Sphere3FIResultType;
    contactTime: number;
    contactPoint: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrAlignedBox3Sphere3FIResult():
    IntrAlignedBox3Sphere3FIResult {
    return {
        intersectionType: IntrAlignedBox3Sphere3FIResultType.noContact,
        contactTime: 0,
        contactPoint: Vector.zero(3)
    };
}

// Currently, only a dynamic query is supported. A static query will need to
// compute the intersection set of (solid) box and sphere.
export class IntrAlignedBox3Sphere3FI {
    find(box: AlignedBox, boxVelocity: Vector, sphere: Hypersphere,
        sphereVelocity: Vector): IntrAlignedBox3Sphere3FIResult {
        logAssert(box.dimension === 3 && boxVelocity.size === 3
            && sphere.center.size === 3 && sphereVelocity.size === 3,
            'IntrAlignedBox3Sphere3FI: mismatched sizes.');
        const result = defaultIntrAlignedBox3Sphere3FIResult();

        // Translate the sphere and box so that the box center becomes the
        // origin. Compute the velocity of the sphere relative to the box.
        const { center: boxCenter, extent } = box.getCenteredForm();
        const C = sub(sphere.center, boxCenter);
        const V = sub(sphereVelocity, boxVelocity);

        // Test for no-intersection that leads to an early exit. The test is
        // fast, using the method of separating axes.
        const superBox = new AlignedBox(3);
        for (let i = 0; i < 3; ++i) {
            superBox.max.values[i] = extent.values[i] + sphere.radius;
            superBox.min.values[i] = -superBox.max.values[i];
        }
        const rbQuery = new IntrRay3AlignedBox3TI();
        const rbResult = rbQuery.test(Ray.fromOriginDirection(C, V), superBox);
        if (rbResult.intersect) {
            doQuery(extent, C, sphere.radius, V, result);

            // Translate the contact point back to the coordinate system of
            // the original sphere and box. Upstream performs this translation
            // unconditionally, so when doQuery reports 'noContact' the
            // contactPoint is boxCenter rather than the zero vector the
            // Result documentation claims. The point is documented as invalid
            // in that case, so the quirk is preserved here (see the PR's
            // upstream bug suspects); callers must test intersectionType.
            result.contactPoint = add(result.contactPoint, boxCenter);
        }
        return result;
    }
}

// The query assumes the box is axis-aligned with center at the origin and
// with extent K. Callers need to convert the results back to the original
// coordinate system of the query.
function doQuery(K: Vector, inC: Vector, radius: number, inV: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    // Change signs on components, if necessary, to transform C to the first
    // octant. Adjust the velocity accordingly.
    const C = inC.clone();
    const V = inV.clone();
    const sign: [number, number, number] = [0, 0, 0];
    for (let i = 0; i < 3; ++i) {
        if (C.values[i] >= 0) {
            sign[i] = 1;
        }
        else {
            C.values[i] = -C.values[i];
            V.values[i] = -V.values[i];
            sign[i] = -1;
        }
    }

    const delta = sub(C, K);
    const d = delta.values;
    if (d[2] <= radius) {
        if (d[1] <= radius) {
            if (d[0] <= radius) {
                if (d[2] <= 0) {
                    if (d[1] <= 0) {
                        if (d[0] <= 0) {
                            interiorOverlap(C, result);
                        }
                        else {
                            // x-face
                            faceOverlap(0, 1, 2, K, C, radius, delta, result);
                        }
                    }
                    else {
                        if (d[0] <= 0) {
                            // y-face
                            faceOverlap(1, 2, 0, K, C, radius, delta, result);
                        }
                        else {
                            // xy-edge
                            if (d[0] * d[0] + d[1] * d[1] <= radius * radius) {
                                edgeOverlap(0, 1, 2, K, C, radius, delta, result);
                            }
                            else {
                                edgeSeparated(0, 1, 2, K, C, radius, delta, V, result);
                            }
                        }
                    }
                }
                else {
                    if (d[1] <= 0) {
                        if (d[0] <= 0) {
                            // z-face
                            faceOverlap(2, 0, 1, K, C, radius, delta, result);
                        }
                        else {
                            // xz-edge
                            if (d[0] * d[0] + d[2] * d[2] <= radius * radius) {
                                edgeOverlap(2, 0, 1, K, C, radius, delta, result);
                            }
                            else {
                                edgeSeparated(2, 0, 1, K, C, radius, delta, V, result);
                            }
                        }
                    }
                    else {
                        if (d[0] <= 0) {
                            // yz-edge
                            if (d[1] * d[1] + d[2] * d[2] <= radius * radius) {
                                edgeOverlap(1, 2, 0, K, C, radius, delta, result);
                            }
                            else {
                                edgeSeparated(1, 2, 0, K, C, radius, delta, V, result);
                            }
                        }
                        else {
                            // xyz-vertex
                            if (dot(delta, delta) <= radius * radius) {
                                vertexOverlap(K, radius, delta, result);
                            }
                            else {
                                vertexSeparated(K, C, radius, delta, V,
                                    result);
                            }
                        }
                    }
                }
            }
            else {
                // x-face
                faceUnbounded(0, 1, 2, K, C, radius, delta, V, result);
            }
        }
        else {
            if (d[0] <= radius) {
                // y-face
                faceUnbounded(1, 2, 0, K, C, radius, delta, V, result);
            }
            else {
                // xy-edge
                edgeUnbounded(0, 1, 2, K, C, radius, delta, V, result);
            }
        }
    }
    else {
        if (d[1] <= radius) {
            if (d[0] <= radius) {
                // z-face
                faceUnbounded(2, 0, 1, K, C, radius, delta, V, result);
            }
            else {
                // xz-edge
                edgeUnbounded(2, 0, 1, K, C, radius, delta, V, result);
            }
        }
        else {
            if (d[0] <= radius) {
                // yz-edge
                edgeUnbounded(1, 2, 0, K, C, radius, delta, V, result);
            }
            else {
                // xyz-vertex
                vertexUnbounded(K, C, radius, delta, V, result);
            }
        }
    }

    if (result.intersectionType !== IntrAlignedBox3Sphere3FIResultType.noContact) {
        // Translate back to the coordinate system of the translated box and
        // sphere.
        for (let i = 0; i < 3; ++i) {
            if (sign[i] < 0) {
                result.contactPoint.values[i] = -result.contactPoint.values[i];
            }
        }
    }
}

function interiorOverlap(C: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    result.intersectionType = IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping;
    result.contactTime = 0;
    result.contactPoint = C.clone();
}

function vertexOverlap(K: Vector, radius: number, delta: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    result.intersectionType = (dot(delta, delta) < radius * radius
        ? IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping
        : IntrAlignedBox3Sphere3FIResultType.contact);
    result.contactTime = 0;
    result.contactPoint = K.clone();
}

function edgeOverlap(i0: number, i1: number, i2: number, K: Vector, C: Vector,
    radius: number, delta: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    const d = delta.values;
    result.intersectionType = (d[i0] * d[i0] + d[i1] * d[i1] < radius * radius
        ? IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping
        : IntrAlignedBox3Sphere3FIResultType.contact);
    result.contactTime = 0;
    result.contactPoint.values[i0] = K.values[i0];
    result.contactPoint.values[i1] = K.values[i1];
    result.contactPoint.values[i2] = C.values[i2];
}

function faceOverlap(i0: number, i1: number, i2: number, K: Vector, C: Vector,
    radius: number, delta: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    result.intersectionType = (delta.values[i0] < radius
        ? IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping
        : IntrAlignedBox3Sphere3FIResultType.contact);
    result.contactTime = 0;
    result.contactPoint.values[i0] = K.values[i0];
    result.contactPoint.values[i1] = C.values[i1];
    result.contactPoint.values[i2] = C.values[i2];
}

function vertexSeparated(K: Vector, C: Vector, radius: number, delta: Vector,
    V: Vector, result: IntrAlignedBox3Sphere3FIResult): void {
    if (V.values[0] < 0 || V.values[1] < 0 || V.values[2] < 0) {
        // Upstream probes only the rounded vertex nearest the sphere center.
        // That vertex is the nearest piece of the rounded box at time zero,
        // but it need not be the piece the ray enters first: the sphere can
        // slide along a box edge and touch the rounded edge (the cylinder)
        // earlier, or miss the vertex altogether and still hit the cylinder
        // (a third instance of the defect of #458 item 5 / #465 item 2, found
        // by the Minkowski-sum oracle of the test file). The three rounded
        // edges that meet the vertex are therefore probed as well and the
        // earliest contact is kept. The vertex is probed first, so a tie
        // keeps upstream's contact point.
        doQueryRayRoundedVertex(K, radius, delta, V, result);
        probeRoundedEdge(0, 1, 2, K, C, radius, V, 1, 1, result);
        probeRoundedEdge(1, 2, 0, K, C, radius, V, 1, 1, result);
        probeRoundedEdge(2, 0, 1, K, C, radius, V, 1, 1, result);
    }
}

function edgeSeparated(i0: number, i1: number, i2: number, K: Vector,
    C: Vector, radius: number, delta: Vector, V: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    if (V.values[i0] < 0 || V.values[i1] < 0) {
        doQueryRayRoundedEdge(i0, i1, i2, K, C, radius, delta, V, result);
    }
}

function vertexUnbounded(K: Vector, C: Vector, radius: number, delta: Vector,
    V: Vector, result: IntrAlignedBox3Sphere3FIResult): void {
    if (V.values[0] < 0 && V.values[1] < 0 && V.values[2] < 0) {
        // Determine the face of the rounded box that is intersected by the
        // ray C+T*V.
        let tmax = (radius - delta.values[0]) / V.values[0];
        let j0 = 0;
        let temp = (radius - delta.values[1]) / V.values[1];
        if (temp > tmax) {
            tmax = temp;
            j0 = 1;
        }
        temp = (radius - delta.values[2]) / V.values[2];
        if (temp > tmax) {
            tmax = temp;
            j0 = 2;
        }

        // The j0-rounded face is the candidate for intersection.
        const j1 = (j0 + 1) % 3;
        const j2 = (j1 + 1) % 3;
        doQueryRayRoundedFace(j0, j1, j2, K, C, radius, delta, V, result);
    }
}

function edgeUnbounded(i0: number, i1: number, _i2: number, K: Vector,
    C: Vector, radius: number, delta: Vector, V: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    if (V.values[i0] < 0 && V.values[i1] < 0) {
        // Determine the face of the rounded box that is intersected by the
        // ray C+T*V.
        let tmax = (radius - delta.values[i0]) / V.values[i0];
        let j0 = i0;
        const temp = (radius - delta.values[i1]) / V.values[i1];
        if (temp > tmax) {
            tmax = temp;
            j0 = i1;
        }

        // The j0-rounded face is the candidate for intersection.
        const j1 = (j0 + 1) % 3;
        const j2 = (j1 + 1) % 3;
        doQueryRayRoundedFace(j0, j1, j2, K, C, radius, delta, V, result);
    }
}

function faceUnbounded(i0: number, i1: number, i2: number, K: Vector,
    C: Vector, radius: number, delta: Vector, V: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    if (V.values[i0] < 0) {
        doQueryRayRoundedFace(i0, i1, i2, K, C, radius, delta, V, result);
    }
}

function doQueryRayRoundedVertex(K: Vector, radius: number, delta: Vector,
    V: Vector, result: IntrAlignedBox3Sphere3FIResult): void {
    const a1 = dot(V, delta);
    if (a1 < 0) {
        // The caller must ensure that a0 > 0 and a2 > 0.
        const a0 = dot(delta, delta) - radius * radius;
        const a2 = dot(V, V);
        const adiscr = a1 * a1 - a2 * a0;
        if (adiscr >= 0) {
            // The ray intersects the rounded vertex, so the sphere-box
            // contact point is the vertex.
            result.intersectionType = IntrAlignedBox3Sphere3FIResultType.contact;
            result.contactTime = -(a1 + Math.sqrt(adiscr)) / a2;
            result.contactPoint = K.clone();
        }
    }
}

function doQueryRayRoundedEdge(i0: number, i1: number, i2: number, K: Vector,
    C: Vector, radius: number, delta: Vector, V: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    const d = delta.values;
    const v = V.values;
    const b1 = v[i0] * d[i0] + v[i1] * d[i1];
    if (b1 < 0) {
        // The caller must ensure that b0 > 0 and b2 > 0.
        const b0 = d[i0] * d[i0] + d[i1] * d[i1] - radius * radius;
        const b2 = v[i0] * v[i0] + v[i1] * v[i1];
        const bdiscr = b1 * b1 - b2 * b0;
        if (bdiscr >= 0) {
            const tmax = -(b1 + Math.sqrt(bdiscr)) / b2;
            const p2 = C.values[i2] + tmax * v[i2];
            if (-K.values[i2] <= p2) {
                if (p2 <= K.values[i2]) {
                    // The ray intersects the finite cylinder of the rounded
                    // edge, so the sphere-box contact point is on the
                    // corresponding box edge.
                    result.intersectionType = IntrAlignedBox3Sphere3FIResultType.contact;
                    result.contactTime = tmax;
                    result.contactPoint.values[i0] = K.values[i0];
                    result.contactPoint.values[i1] = K.values[i1];
                    result.contactPoint.values[i2] = p2;
                }
                else {
                    // The ray intersects the infinite cylinder but not the
                    // finite cylinder of the rounded edge. It is possible the
                    // ray intersects the rounded vertex for K.
                    doQueryRayRoundedVertex(K, radius, delta, V, result);
                }
            }
            else {
                // The ray intersects the infinite cylinder but not the finite
                // cylinder of the rounded edge. It is possible the ray
                // intersects the rounded vertex for otherK.
                const otherK = Vector.zero(3);
                const otherDelta = Vector.zero(3);
                otherK.values[i0] = K.values[i0];
                otherK.values[i1] = K.values[i1];
                otherK.values[i2] = -K.values[i2];
                otherDelta.values[i0] = C.values[i0] - otherK.values[i0];
                otherDelta.values[i1] = C.values[i1] - otherK.values[i1];
                otherDelta.values[i2] = C.values[i2] - otherK.values[i2];
                doQueryRayRoundedVertex(otherK, radius, otherDelta, V, result);
            }
        }
    }
}

// Keep the earliest of the candidate contacts. Upstream accepts the first
// rounded-edge probe that reports a contact and never looks at the other
// candidates of the same face; see the comment in doQueryRayRoundedFace.
function keepEarliest(result: IntrAlignedBox3Sphere3FIResult,
    candidate: IntrAlignedBox3Sphere3FIResult): void {
    if (candidate.intersectionType === IntrAlignedBox3Sphere3FIResultType.noContact) {
        return;
    }
    if (result.intersectionType === IntrAlignedBox3Sphere3FIResultType.noContact
        || candidate.contactTime < result.contactTime) {
        result.intersectionType = candidate.intersectionType;
        result.contactTime = candidate.contactTime;
        result.contactPoint = candidate.contactPoint;
    }
}

// Probe the rounded edge parallel to axis i2 whose cylinder axis passes
// through (sign0 * K[i0], sign1 * K[i1]) and merge the candidate contact into
// 'result', keeping the earliest contact found so far. The rounded vertices
// at the two ends of that edge are probed by doQueryRayRoundedEdge itself.
function probeRoundedEdge(i0: number, i1: number, i2: number, K: Vector,
    C: Vector, radius: number, V: Vector, sign0: number, sign1: number,
    result: IntrAlignedBox3Sphere3FIResult): void {
    const edgeK = Vector.zero(3);
    const edgeDelta = Vector.zero(3);
    edgeK.values[i0] = sign0 * K.values[i0];
    edgeK.values[i1] = sign1 * K.values[i1];
    edgeK.values[i2] = K.values[i2];
    edgeDelta.values[i0] = C.values[i0] - edgeK.values[i0];
    edgeDelta.values[i1] = C.values[i1] - edgeK.values[i1];
    edgeDelta.values[i2] = C.values[i2] - edgeK.values[i2];
    // doQueryRayRoundedEdge requires b0 > 0, that is, the sphere center must
    // be outside the infinite cylinder of the rounded edge; otherwise the
    // quadratic has roots of opposite signs and the "first" root is negative.
    // A center inside that cylinder is beyond one of the cylinder's ends (it
    // is outside the rounded box), so the rounded vertex there - probed by
    // the other candidates - is the piece the ray can enter first.
    const b0 = edgeDelta.values[i0] * edgeDelta.values[i0]
        + edgeDelta.values[i1] * edgeDelta.values[i1] - radius * radius;
    if (b0 <= 0) {
        return;
    }

    const candidate = defaultIntrAlignedBox3Sphere3FIResult();
    doQueryRayRoundedEdge(i0, i1, i2, edgeK, C, radius, edgeDelta, V,
        candidate);
    keepEarliest(result, candidate);
}

function doQueryRayRoundedFace(i0: number, i1: number, i2: number, K: Vector,
    C: Vector, radius: number, delta: Vector, V: Vector,
    result: IntrAlignedBox3Sphere3FIResult): void {
    // The ray enters the box expanded by 'radius' (the bounding box of the
    // rounded box) through the i0-face at time tmax, so the first contact
    // cannot be earlier than tmax.
    const tmax = (radius - delta.values[i0]) / V.values[i0];
    const p1 = C.values[i1] + tmax * V.values[i1];
    const p2 = C.values[i2] + tmax * V.values[i2];

    // The side of the i0-face the ray leaves on, 0 when it stays within the
    // face. Upstream's three-way tests on p1 and p2 are reproduced exactly,
    // so the thresholds and their inclusive/exclusive ends are unchanged.
    const s1 = (p1 < -K.values[i1] ? -1 : (p1 > K.values[i1] ? 1 : 0));
    const s2 = (p2 < -K.values[i2] ? -1 : (p2 > K.values[i2] ? 1 : 0));

    if (s1 === 0 && s2 === 0) {
        // The ray intersects the i0-face of the rounded box, so the
        // sphere-box contact point is on the corresponding box face. This is
        // the first contact, because no contact can precede tmax.
        result.intersectionType = IntrAlignedBox3Sphere3FIResultType.contact;
        result.contactTime = tmax;
        result.contactPoint.values[i0] = K.values[i0];
        result.contactPoint.values[i1] = p1;
        result.contactPoint.values[i2] = p2;
        return;
    }

    // The ray leaves the i0-face, so the first contact - if any - is on a
    // rounded edge or rounded vertex of the face's overhang region.
    //
    // Upstream probes at most the two rounded edges that bound the i0-face
    // and stops at the first probe that reports a contact, which is wrong
    // twice over (see the PR for #458 item 5 and #465 item 2):
    //   * a probe that ends in its rounded-vertex fallback reports a contact
    //     with the vertex, which can be LATER than the contact with the other
    //     rounded edge (the accepted probe is only an upper bound for the
    //     first contact, since every piece is part of the rounded box);
    //   * when the ray leaves the face on both sides it can enter through the
    //     rounded edge PARALLEL to i0 at the overhang corner, which is not
    //     one of the two edges bounding the i0-face and which upstream never
    //     probes, so the contact is missed entirely.
    // All candidate pieces are therefore probed and the earliest contact is
    // kept. The pieces are subsets of the rounded box, so the earliest of
    // their contact times is the first time the ray enters the rounded box.
    if (s1 !== 0) {
        // The rounded (i0,i1)-edge on the p1 side of the face; its cylinder
        // is parallel to i2.
        probeRoundedEdge(i0, i1, i2, K, C, radius, V, 1, s1, result);
    }
    if (s2 !== 0) {
        // The rounded (i2,i0)-edge on the p2 side of the face; its cylinder
        // is parallel to i1.
        probeRoundedEdge(i2, i0, i1, K, C, radius, V, s2, 1, result);
    }
    if (s1 !== 0 && s2 !== 0) {
        // The rounded (i1,i2)-edge at the overhang corner; its cylinder is
        // parallel to i0 and it is not adjacent to the i0-face.
        probeRoundedEdge(i1, i2, i0, K, C, radius, V, s1, s2, result);
    }
}
