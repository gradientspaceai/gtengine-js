// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrSphere3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the spheres to be solids.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The upstream
// 'int32_t type' documentation block becomes the exported enum
// IntrSphere3Sphere3FIResultType, whose name is file-qualified because
// src/index.ts star-exports every file.

import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';
import type { Hypersphere } from './Hypersphere.js';
import { Circle3 } from './Circle3.js';
import { Vector, add, dot, mul, normalize, sub } from './Vector.js';

// The type of intersection reported by IntrSphere3Sphere3FI. Upstream stores
// these as raw integers, with -1 for the default-constructed result.
export enum IntrSphere3Sphere3FIResultType {
    // The default-constructed value, meaning the query has not been run.
    invalid = -1,

    // The spheres are disjoint and separated.
    separated = 0,

    // The spheres touch at a point, each sphere outside the other.
    touchingOutside = 1,

    // The spheres intersect in a circle.
    circle = 2,

    // Sphere0 is strictly contained in sphere1.
    sphere0StrictlyInside = 3,

    // Sphere0 is contained in sphere1 and they share a common point.
    sphere0InsideTouching = 4,

    // Sphere1 is strictly contained in sphere0.
    sphere1StrictlyInside = 5,

    // Sphere1 is contained in sphere0 and they share a common point.
    sphere1InsideTouching = 6
}

// The result of IntrSphere3Sphere3TI queries.
export interface IntrSphere3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrSphere3Sphere3TIResult(): IntrSphere3Sphere3TIResult {
    return { intersect: false };
}

// The result of IntrSphere3Sphere3FI queries.
export interface IntrSphere3Sphere3FIResult {
    intersect: boolean;

    // The type of intersection.
    type: IntrSphere3Sphere3FIResultType;

    // Valid for the types touchingOutside, sphere0InsideTouching,
    // sphere1InsideTouching and (as a representative overlap point) for the
    // strictly-contained types.
    point: Vector;

    // Valid for the type 'circle'.
    circle: Circle3;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrSphere3Sphere3FIResult(): IntrSphere3Sphere3FIResult {
    return {
        intersect: false,
        type: IntrSphere3Sphere3FIResultType.invalid,
        point: Vector.zero(3),
        circle: Circle3.fromCenterNormalRadius(Vector.zero(3),
            Vector.zero(3), 0)
    };
}

// Test-intersection query for two solid spheres in 3D.
export class IntrSphere3Sphere3TI implements
    TIQuery<Hypersphere, Hypersphere, IntrSphere3Sphere3TIResult> {

    test(sphere0: Hypersphere, sphere1: Hypersphere): IntrSphere3Sphere3TIResult {
        const result = defaultIntrSphere3Sphere3TIResult();
        const diff = sub(sphere1.center, sphere0.center);
        const rSum = sphere0.radius + sphere1.radius;
        result.intersect = (dot(diff, diff) <= rSum * rSum);
        return result;
    }
}

// Find-intersection query for two solid spheres in 3D.
export class IntrSphere3Sphere3FI implements
    FIQuery<Hypersphere, Hypersphere, IntrSphere3Sphere3FIResult> {

    find(sphere0: Hypersphere, sphere1: Hypersphere): IntrSphere3Sphere3FIResult {
        const result = defaultIntrSphere3Sphere3FIResult();

        // The plane of intersection must have C1-C0 as its normal direction.
        const C1mC0 = sub(sphere1.center, sphere0.center);
        const sqrLen = dot(C1mC0, C1mC0);
        const r0 = sphere0.radius, r1 = sphere1.radius;
        const rSum = r0 + r1;
        const rSumSqr = rSum * rSum;

        if (sqrLen > rSumSqr) {
            // The spheres are disjoint/separated.
            result.intersect = false;
            result.type = IntrSphere3Sphere3FIResultType.separated;
            return result;
        }

        if (sqrLen === rSumSqr) {
            // The spheres are just touching with each sphere outside the
            // other.
            normalize(C1mC0);
            result.intersect = true;
            result.type = IntrSphere3Sphere3FIResultType.touchingOutside;
            result.point = add(sphere0.center, mul(r0, C1mC0));
            return result;
        }

        const rDif = r0 - r1;
        const rDifSqr = rDif * rDif;
        if (sqrLen < rDifSqr) {
            // One sphere is strictly contained in the other. Compute a point
            // in the intersection set.
            result.intersect = true;
            result.type = (rDif <= 0
                ? IntrSphere3Sphere3FIResultType.sphere0StrictlyInside
                : IntrSphere3Sphere3FIResultType.sphere1StrictlyInside);
            result.point = mul(0.5, add(sphere0.center, sphere1.center));
            return result;
        }
        if (sqrLen === rDifSqr) {
            // One sphere is contained in the other sphere but with a single
            // point of contact.
            normalize(C1mC0);
            result.intersect = true;
            if (rDif <= 0) {
                // Upstream bug (FIXED; see upstream-bug issue (B71)):
                // upstream computes 'sphere1.center + r1 * C1mC0', which is
                // the antipode of the true contact point on sphere1. With
                // sphere0 = ((2,0,0),1) and sphere1 = ((0,0,0),3), the spheres
                // touch at (3,0,0) but upstream reports (-3,0,0). The correct
                // point is C1 - r1 * C1mC0, equivalently
                // C0 + r0 * (C0-C1)/|C0-C1|; the 'else' branch below, which
                // handles sphere1 inside sphere0, already uses this sign.
                result.type = IntrSphere3Sphere3FIResultType.sphere0InsideTouching;
                result.point = sub(sphere1.center, mul(r1, C1mC0));
            }
            else {
                result.type = IntrSphere3Sphere3FIResultType.sphere1InsideTouching;
                result.point = add(sphere0.center, mul(r0, C1mC0));
            }
            return result;
        }

        // Compute t for which the circle of intersection has center
        // K = C0 + t*(C1 - C0).
        const t = 0.5 * (1 + rDif * rSum / sqrLen);

        // Compute the center and radius of the circle of intersection.
        result.circle.center = add(sphere0.center, mul(t, C1mC0));
        result.circle.radius =
            Math.sqrt(Math.max(r0 * r0 - t * t * sqrLen, 0));

        // Compute the normal for the plane of the circle.
        normalize(C1mC0);
        result.circle.normal = C1mC0;

        // The intersection is a circle.
        result.intersect = true;
        result.type = IntrSphere3Sphere3FIResultType.circle;
        return result;
    }
}
