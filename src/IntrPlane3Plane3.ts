// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrPlane3Plane3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has both
// a TIQuery and an FIQuery specialization, which become IntrPlane3Plane3TI
// and IntrPlane3Plane3FI.

import type { FIQuery } from './FIQuery';
import { Hyperplane } from './Hyperplane';
import type { Plane3 } from './Hyperplane';
import { Line } from './Line';
import type { Line3 } from './Line';
import { logAssert } from './Logger';
import { Vector, add, dot, mul } from './Vector';
import { unitCross } from './Vector3';
import type { TIQuery } from './TIQuery';

// The result of IntrPlane3Plane3TI.test.
export interface IntrPlane3Plane3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrPlane3Plane3TIResult {
    return { intersect: false };
}

// The result of IntrPlane3Plane3FI.find.
export interface IntrPlane3Plane3FIResult {
    intersect: boolean;

    // If 'intersect' is true, the intersection is either a line or the planes
    // are the same. When a line, 'line' is valid. When the same plane,
    // 'plane' is set to one of the planes.
    isLine: boolean;
    line: Line3;
    plane: Plane3;
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrPlane3Plane3FIResult(): IntrPlane3Plane3FIResult {
    return {
        intersect: false,
        isLine: false,
        line: Line.fromOriginDirection(Vector.zero(3), Vector.zero(3)),
        plane: Hyperplane.fromNormalConstant(Vector.zero(3), 0)
    };
}

// Test-intersection query for two planes in 3D.
export class IntrPlane3Plane3TI implements
    TIQuery<Plane3, Plane3, IntrPlane3Plane3TIResult> {

    test(plane0: Plane3, plane1: Plane3): IntrPlane3Plane3TIResult {
        logAssert(plane0.dimension === 3 && plane1.dimension === 3,
            'IntrPlane3Plane3TI: mismatched sizes.');

        // If Cross(N0,N1) is zero, then either the planes are parallel and
        // separated or they are the same plane. In both cases, 'false' is
        // returned. Otherwise, the planes intersect. To avoid subtle
        // differences in reporting between the test and find queries, the
        // same parallel test is used. Mathematically,
        //   |Cross(N0,N1)|^2 = Dot(N0,N0)*Dot(N1,N1)-Dot(N0,N1)^2
        //                    = 1 - Dot(N0,N1)^2
        // The last equality is true since planes are required to have
        // unit-length normal vectors. The test |Cross(N0,N1)| = 0 is the same
        // as |Dot(N0,N1)| = 1.
        const result = defaultTIResult();
        const d = dot(plane0.normal, plane1.normal);
        if (Math.abs(d) < 1) {
            result.intersect = true;
            return result;
        }

        // The planes are parallel. Check whether they are coplanar.
        let cDiff: number;
        if (d >= 0) {
            // Normals are in the same direction, need to look at c0-c1.
            cDiff = plane0.constant - plane1.constant;
        }
        else {
            // Normals are in opposite directions, need to look at c0+c1.
            cDiff = plane0.constant + plane1.constant;
        }

        result.intersect = (Math.abs(cDiff) === 0);
        return result;
    }
}

// Find-intersection query for two planes in 3D.
export class IntrPlane3Plane3FI implements
    FIQuery<Plane3, Plane3, IntrPlane3Plane3FIResult> {

    find(plane0: Plane3, plane1: Plane3): IntrPlane3Plane3FIResult {
        logAssert(plane0.dimension === 3 && plane1.dimension === 3,
            'IntrPlane3Plane3FI: mismatched sizes.');

        // If N0 and N1 are parallel, either the planes are parallel and
        // separated or they are the same plane. In both cases, 'false' is
        // returned. Otherwise, the intersection line is
        //   L(t) = t*Cross(N0,N1)/|Cross(N0,N1)| + c0*N0 + c1*N1
        // for some coefficients c0 and c1 and for t any real number (the line
        // parameter). Taking dot products with the normals,
        //   d0 = Dot(N0,L) = c0*Dot(N0,N0) + c1*Dot(N0,N1) = c0 + c1*d
        //   d1 = Dot(N1,L) = c0*Dot(N0,N1) + c1*Dot(N1,N1) = c0*d + c1
        // where d = Dot(N0,N1). These are two equations in two unknowns. The
        // solution is
        //   c0 = (d0 - d*d1)/det
        //   c1 = (d1 - d*d0)/det
        // where det = 1 - d^2.
        const result = defaultIntrPlane3Plane3FIResult();

        const d = dot(plane0.normal, plane1.normal);
        if (Math.abs(d) >= 1) {
            // The planes are parallel. Check if they are coplanar.
            let cDiff: number;
            if (d >= 0) {
                // Normals are in the same direction, need to look at c0-c1.
                cDiff = plane0.constant - plane1.constant;
            }
            else {
                // Normals are in opposite directions, need to look at c0+c1.
                cDiff = plane0.constant + plane1.constant;
            }

            if (Math.abs(cDiff) === 0) {
                // The planes are coplanar.
                result.intersect = true;
                result.isLine = false;
                result.plane = plane0.clone();
                return result;
            }

            // The planes are parallel but distinct.
            result.intersect = false;
            return result;
        }

        const invDet = 1 / (1 - d * d);
        const c0 = (plane0.constant - d * plane1.constant) * invDet;
        const c1 = (plane1.constant - d * plane0.constant) * invDet;
        result.intersect = true;
        result.isLine = true;
        result.line.origin = add(mul(c0, plane0.normal),
            mul(c1, plane1.normal));
        result.line.direction = unitCross(plane0.normal, plane1.normal);
        return result;
    }
}
