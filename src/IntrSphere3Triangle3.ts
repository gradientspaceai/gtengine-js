// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSphere3Triangle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Currently, only a dynamic query is supported. A static query will need to
// compute the intersection set of triangle and sphere.
//
// The query computes the first time of contact between a moving sphere and a
// moving triangle. It reduces to a ray-versus-sphere-swept-volume query,
// where the sphere-swept volume is the Minkowski sum of the triangle and the
// sphere: two triangular faces offset by the radius along the triangle
// normal, three half cylinders along the edges and three sphere wedges at the
// vertices.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// two 'operator()' overloads selected by 'is_arbitrary_precision', one for
// floating-point types and one for arbitrary-precision types (the latter
// reports the contact time and point as QFNumber quadratic-field values).
// Per PORTING.md ("port only the floating-point instantiation unless the
// rational/exact path is required by a dependent file"), only the
// floating-point query is ported; no file in the library uses the exact path.
// There is no TIQuery upstream, so this file exports only the FI class.

import type { Hypersphere } from './Hypersphere.js';
import type { Triangle } from './Triangle.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { cross, unitCross } from './Vector3.js';
import { DistPointTriangle } from './DistPointTriangle.js';

// The kind of contact reported by IntrSphere3Triangle3FI.
export enum IntrSphere3Triangle3FIResultType {
    // The objects are initially overlapping. The contactPoint is only one of
    // infinitely many points in the overlap: the triangle point closest to
    // the sphere center. The contactTime is 0.
    initiallyOverlapping = -1,

    // The objects are initially separated and do not intersect later. The
    // contactTime and contactPoint are invalid (both zero).
    noContact = 0,

    // The objects are initially separated but intersect later. The
    // contactTime is the first time T >= 0 of contact.
    //
    // Upstream quirk (preserved): for this case upstream sets
    // contactPoint = sphere.center + contactTime * sphereVelocity, which is
    // the position of the sphere CENTER at the contact time, not the point
    // where the two surfaces touch (the documentation says "corresponding
    // first contact"). The surface contact point is the closest point of the
    // moved triangle to that center. For the initiallyOverlapping case the
    // contactPoint IS a triangle point, so the two cases are inconsistent.
    contact = 1
}

// The result of IntrSphere3Triangle3FI queries.
export interface IntrSphere3Triangle3FIResult {
    intersectionType: IntrSphere3Triangle3FIResultType;
    contactTime: number;
    contactPoint: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSphere3Triangle3FIResult(): IntrSphere3Triangle3FIResult {
    return {
        intersectionType: IntrSphere3Triangle3FIResultType.noContact,
        contactTime: 0,
        contactPoint: Vector.zero(3)
    };
}

// Find-intersection (dynamic) query for a moving sphere and a moving triangle
// in 3D.
export class IntrSphere3Triangle3FI {
    find(sphere: Hypersphere, sphereVelocity: Vector, triangle: Triangle,
        triangleVelocity: Vector): IntrSphere3Triangle3FIResult {
        const result = defaultIntrSphere3Triangle3FIResult();

        // Test for initial overlap or contact.
        const ptQuery = new DistPointTriangle();
        const ptResult = ptQuery.compute(sphere.center, triangle);
        const rsqr = sphere.radius * sphere.radius;
        if (ptResult.sqrDistance <= rsqr) {
            result.intersectionType = (ptResult.sqrDistance < rsqr
                ? IntrSphere3Triangle3FIResultType.initiallyOverlapping
                : IntrSphere3Triangle3FIResultType.contact);
            result.contactTime = 0;
            result.contactPoint = ptResult.closest[1].clone();
            return result;
        }

        // To reach here, the sphere and triangle are initially separated.
        // Compute the velocity of the sphere relative to the triangle.
        const V = sub(sphereVelocity, triangleVelocity);
        const sqrLenV = dot(V, V);
        if (sqrLenV === 0) {
            // The sphere and triangle are separated and the sphere is not
            // moving relative to the triangle, so there is no contact. The
            // 'result' is already set to the correct state for this case.
            return result;
        }

        // Compute the triangle edge directions E[], the vector U normal to
        // the plane of the triangle, and the normals to the edges in the
        // plane of the triangle.
        const E = [
            sub(triangle.v[1], triangle.v[0]),
            sub(triangle.v[2], triangle.v[1]),
            sub(triangle.v[0], triangle.v[2])
        ];
        const sqrLenE = [dot(E[0], E[0]), dot(E[1], E[1]), dot(E[2], E[2])];
        const U = unitCross(E[0], E[1]);
        const ExU = [cross(E[0], U), cross(E[1], U), cross(E[2], U)];

        // Compute the vectors from the triangle vertices to the sphere
        // center.
        const Delta = [
            sub(sphere.center, triangle.v[0]),
            sub(sphere.center, triangle.v[1]),
            sub(sphere.center, triangle.v[2])
        ];

        // Determine where the sphere center is located relative to the planes
        // of the triangle offset faces of the sphere-swept volume.
        const dotUDelta0 = dot(U, Delta[0]);
        if (dotUDelta0 >= sphere.radius) {
            // The sphere is on the positive side of Dot(U,X-C) = r. If the
            // sphere will contact the sphere-swept volume at a triangular
            // face, it can do so only on the face of the aforementioned
            // plane.
            const dotUV = dot(U, V);
            if (dotUV >= 0) {
                // The sphere is moving away from, or parallel to, the plane
                // of the triangle. The 'result' is already set to the correct
                // state for this case.
                return result;
            }

            const tbar = (sphere.radius - dotUDelta0) / dotUV;
            let foundContact = true;
            for (let i = 0; i < 3; ++i) {
                const phi = dot(ExU[i], Delta[i]);
                const psi = dot(ExU[i], V);
                if (phi + psi * tbar > 0) {
                    foundContact = false;
                    break;
                }
            }
            if (foundContact) {
                result.intersectionType = IntrSphere3Triangle3FIResultType.contact;
                result.contactTime = tbar;
                result.contactPoint = add(sphere.center,
                    mul(tbar, sphereVelocity));
                return result;
            }
        }
        else if (dotUDelta0 <= -sphere.radius) {
            // The sphere is on the positive side of Dot(-U,X-C) = r. If the
            // sphere will contact the sphere-swept volume at a triangular
            // face, it can do so only on the face of the aforementioned
            // plane.
            const dotUV = dot(U, V);
            if (dotUV <= 0) {
                // The sphere is moving away from, or parallel to, the plane
                // of the triangle. The 'result' is already set to the correct
                // state for this case.
                return result;
            }

            const tbar = (-sphere.radius - dotUDelta0) / dotUV;
            let foundContact = true;
            for (let i = 0; i < 3; ++i) {
                const phi = dot(ExU[i], Delta[i]);
                const psi = dot(ExU[i], V);
                if (phi + psi * tbar > 0) {
                    foundContact = false;
                    break;
                }
            }
            if (foundContact) {
                result.intersectionType = IntrSphere3Triangle3FIResultType.contact;
                result.contactTime = tbar;
                result.contactPoint = add(sphere.center,
                    mul(tbar, sphereVelocity));
                return result;
            }
        }
        // else: The ray-sphere-swept-volume contact point (if any) cannot be
        // on a triangular face of the sphere-swept-volume.

        // The sphere is moving towards the slab between the two planes of the
        // sphere-swept volume triangular faces. Determine whether the ray
        // intersects the half cylinders or sphere wedges of the sphere-swept
        // volume.

        // Test for contact with half cylinders of the sphere-swept volume.
        // First, precompute some dot products required in the computations.
        const del = [0, 0, 0];
        const delp = [0, 0, 0];
        const nu = [0, 0, 0];
        for (let im1 = 2, i = 0; i < 3; im1 = i++) {
            del[i] = dot(E[i], Delta[i]);
            delp[im1] = dot(E[im1], Delta[i]);
            nu[i] = dot(E[i], V);
        }

        for (let i = 2, ip1 = 0; ip1 < 3; i = ip1++) {
            const hatV = sub(V, mul(nu[i] / sqrLenE[i], E[i]));
            const sqrLenHatV = dot(hatV, hatV);
            if (sqrLenHatV > 0) {
                const hatDelta = sub(Delta[i], mul(del[i] / sqrLenE[i], E[i]));
                const alpha = -dot(hatV, hatDelta);
                if (alpha >= 0) {
                    const sqrLenHatDelta = dot(hatDelta, hatDelta);
                    const beta = alpha * alpha -
                        sqrLenHatV * (sqrLenHatDelta - rsqr);
                    if (beta >= 0) {
                        const tbar = (alpha - Math.sqrt(beta)) / sqrLenHatV;

                        const mu = dot(ExU[i], Delta[i]);
                        const omega = dot(ExU[i], hatV);
                        if (mu + omega * tbar >= 0) {
                            if (del[i] + nu[i] * tbar >= 0) {
                                if (delp[i] + nu[i] * tbar <= 0) {
                                    // The constraints are satisfied, so tbar
                                    // is the first time of contact.
                                    result.intersectionType =
                                        IntrSphere3Triangle3FIResultType.contact;
                                    result.contactTime = tbar;
                                    result.contactPoint = add(sphere.center,
                                        mul(tbar, sphereVelocity));
                                    return result;
                                }
                            }
                        }
                    }
                }
            }
        }

        // Test for contact with sphere wedges of the sphere-swept volume. We
        // know that |V|^2 > 0 because of a previous early-exit test.
        for (let im1 = 2, i = 0; i < 3; im1 = i++) {
            const alpha = -dot(V, Delta[i]);
            if (alpha >= 0) {
                const sqrLenDelta = dot(Delta[i], Delta[i]);
                const beta = alpha * alpha - sqrLenV * (sqrLenDelta - rsqr);
                if (beta >= 0) {
                    const tbar = (alpha - Math.sqrt(beta)) / sqrLenV;
                    if (delp[im1] + nu[im1] * tbar >= 0) {
                        if (del[i] + nu[i] * tbar <= 0) {
                            // The constraints are satisfied, so tbar is the
                            // first time of contact.
                            result.intersectionType =
                                IntrSphere3Triangle3FIResultType.contact;
                            result.contactTime = tbar;
                            result.contactPoint = add(sphere.center,
                                mul(tbar, sphereVelocity));
                            return result;
                        }
                    }
                }
            }
        }

        // The ray and sphere-swept volume do not intersect, so the sphere and
        // triangle do not come into contact. The 'result' is already set to
        // the correct state for this case.
        return result;
    }
}
