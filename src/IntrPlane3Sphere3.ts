// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Sphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the sphere to be a solid. The comparisons of the
// point-plane distance to the sphere radius assume the plane normal is unit
// length.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrPlane3Sphere3TI
// and IntrPlane3Sphere3FI.

import { Circle3 } from './Circle3.js';
import { DistPointHyperplane } from './DistPointHyperplane.js';
import type { FIQuery } from './FIQuery.js';
import type { Hypersphere } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import type { Plane3 } from './Hyperplane.js';
import { Vector, mul, sub } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrPlane3Sphere3TI.test.
export interface IntrPlane3Sphere3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
export function defaultIntrPlane3Sphere3TIResult(): IntrPlane3Sphere3TIResult {
    return { intersect: false };
}

// The result of IntrPlane3Sphere3FI.find.
export interface IntrPlane3Sphere3FIResult {
    intersect: boolean;

    // If 'intersect' is true, the intersection is either a point or a circle.
    // When 'isCircle' is true, 'circle' is valid. When 'isCircle' is false,
    // 'point' is valid.
    isCircle: boolean;
    circle: Circle3;
    point: Vector;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrPlane3Sphere3FIResult(): IntrPlane3Sphere3FIResult {
    return {
        intersect: false,
        isCircle: false,
        // Upstream's Result constructor builds the circle from an explicit
        // (zero center, zero normal, zero radius); that is not the Circle3
        // default constructor, which has a unit normal and radius 1.
        circle: Circle3.fromCenterNormalRadius(Vector.zero(3), Vector.zero(3), 0),
        point: Vector.zero(3)
    };
}

// Test-intersection query for a plane and a solid sphere in 3D.
export class IntrPlane3Sphere3TI implements
    TIQuery<Plane3, Hypersphere, IntrPlane3Sphere3TIResult> {

    test(plane: Plane3, sphere: Hypersphere): IntrPlane3Sphere3TIResult {
        logAssert(plane.dimension === 3 && sphere.dimension === 3,
            'IntrPlane3Sphere3TI: mismatched sizes.');
        const result = defaultIntrPlane3Sphere3TIResult();
        const ppQuery = new DistPointHyperplane();
        const ppResult = ppQuery.compute(sphere.center, plane);
        result.intersect = (ppResult.distance <= sphere.radius);
        return result;
    }
}

// Find-intersection query for a plane and a solid sphere in 3D.
export class IntrPlane3Sphere3FI implements
    FIQuery<Plane3, Hypersphere, IntrPlane3Sphere3FIResult> {

    find(plane: Plane3, sphere: Hypersphere): IntrPlane3Sphere3FIResult {
        logAssert(plane.dimension === 3 && sphere.dimension === 3,
            'IntrPlane3Sphere3FI: mismatched sizes.');
        const result = defaultIntrPlane3Sphere3FIResult();
        const ppQuery = new DistPointHyperplane();
        const ppResult = ppQuery.compute(sphere.center, plane);
        if (ppResult.distance < sphere.radius) {
            result.intersect = true;
            result.isCircle = true;
            result.circle.center = sub(sphere.center,
                mul(ppResult.signedDistance, plane.normal));
            result.circle.normal = plane.normal.clone();

            // The sum and difference are both positive numbers.
            const sum = sphere.radius + ppResult.distance;
            const dif = sphere.radius - ppResult.distance;

            // arg = sqr(sphere.radius) - sqr(ppResult.distance)
            const arg = sum * dif;

            result.circle.radius = Math.sqrt(arg);
            return result;
        }
        else if (ppResult.distance === sphere.radius) {
            result.intersect = true;
            result.isCircle = false;
            result.point = sub(sphere.center,
                mul(ppResult.signedDistance, plane.normal));
            return result;
        }
        else {
            result.intersect = false;
            return result;
        }
    }
}
