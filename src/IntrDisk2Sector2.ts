// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrDisk2Sector2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Circle2 object is considered to be a disk whose points X satisfy the
// constraint |X-C| <= R, where C is the disk center and R is the disk radius.
// The Sector2 object is also considered to be a solid. Also, the Sector2
// object is required to be convex, so the sector angle must be in (0,pi/2],
// even though the Sector2 definition allows for angles larger than pi/2
// (leading to nonconvex sectors). The sector vertex is V, the radius is L,
// the axis direction is D, and the angle is A. Sector points X satisfy
// |X-V| <= L and Dot(D,X-V) >= cos(A)|X-V| >= 0.
//
// A subproblem for the test-intersection query is to determine whether the
// disk intersects the cone of the sector. Although the query is in 2D, it is
// analogous to the 3D problem of determining whether a sphere and cone
// overlap. That algorithm is described in
//   https://www.geometrictools.com/Documentation/IntersectionSphereCone.pdf
// The algorithm leads to coordinate-free pseudocode that applies to 2D as
// well as 3D. That function is the first SphereIntersectsCone on page 4 of
// the PDF.
//
// If the disk is outside the cone, there is no intersection. If the disk
// overlaps the cone, we then need to test whether the disk overlaps the disk
// of the sector.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrDisk2Sector2TI. The
// upstream file name uses 'Disk2' although the query type is Circle2, which
// is the ported Hypersphere with dimension 2.

import type { Hypersphere } from './Hypersphere.js';
import type { Sector2 } from './Sector2.js';
import { Vector, sub, mul, dot, length } from './Vector.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrDisk2Sector2TI.test.
export interface IntrDisk2Sector2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrDisk2Sector2TIResult {
    return { intersect: false };
}

export class IntrDisk2Sector2TI implements
    TIQuery<Hypersphere, Sector2, IntrDisk2Sector2TIResult> {

    test(disk: Hypersphere, sector: Sector2): IntrDisk2Sector2TIResult {
        const result = defaultTIResult();

        // Test whether the disk and the disk of the sector overlap.
        const CmV = sub(disk.center, sector.vertex);
        const sqrLengthCmV = dot(CmV, CmV);
        const lengthCmV = Math.sqrt(sqrLengthCmV);
        if (lengthCmV > disk.radius + sector.radius) {
            // The disk is outside the disk of the sector.
            result.intersect = false;
            return result;
        }

        // Test whether the disk and cone of the sector overlap. The comments
        // about K, K' and K" refer to the PDF mentioned previously.
        const U = sub(sector.vertex,
            mul(disk.radius / sector.sinAngle, sector.direction));
        const CmU = sub(disk.center, U);
        const lengthCmU = length(CmU);
        if (dot(sector.direction, CmU) < lengthCmU * sector.cosAngle) {
            // The disk center is outside K" (in the white or gray regions).
            result.intersect = false;
            return result;
        }

        // The disk center is inside K" (in the red, orange, blue or green
        // regions).
        const dotDirCmV = dot(sector.direction, CmV);
        if (-dotDirCmV >= lengthCmV * sector.sinAngle) {
            // The disk center is inside K" and inside K' (in the blue or
            // green regions).
            if (lengthCmV <= disk.radius) {
                // The disk center is in the blue region, in which case the
                // disk contains the sector vertex.
                result.intersect = true;
            } else {
                // The disk center is in the green region.
                result.intersect = false;
            }
            return result;
        }

        // To reach here, we know that the disk overlaps the sector disk and
        // the sector cone. The disk center is in the orange region or in the
        // red region (not including the segments that separate the red and
        // blue regions).

        // Test whether the ray of the right boundary of the sector overlaps
        // the disk. The ray direction U0 is a clockwise rotation of the cone
        // axis by the cone angle.
        const d = sector.direction.values;
        const U0 = Vector.fromArray([
            +sector.cosAngle * d[0] + sector.sinAngle * d[1],
            -sector.sinAngle * d[0] + sector.cosAngle * d[1]
        ]);
        const dp0 = dot(U0, CmV);
        const discr0 = disk.radius * disk.radius + dp0 * dp0 - sqrLengthCmV;
        if (discr0 >= 0) {
            // The ray intersects the disk. Now test whether the sector
            // boundary segment contained by the ray overlaps the disk. The
            // quadratic root tmin generates the ray-disk point of
            // intersection closest to the sector vertex.
            const tmin = dp0 - Math.sqrt(discr0);
            if (sector.radius >= tmin) {
                // The segment overlaps the disk.
                result.intersect = true;
                return result;
            } else {
                // The segment does not overlap the disk. We know the disks
                // overlap, so if the disk center is outside the sector cone
                // or on the right-boundary ray, the overlap occurs outside
                // the cone, which implies the disk and sector do not
                // intersect.
                if (dotDirCmV <= lengthCmV * sector.cosAngle) {
                    // The disk center is not inside the sector cone.
                    result.intersect = false;
                    return result;
                }
            }
        }

        // Test whether the ray of the left boundary of the sector overlaps
        // the disk. The ray direction U1 is a counterclockwise rotation of
        // the cone axis by the cone angle.
        const U1 = Vector.fromArray([
            +sector.cosAngle * d[0] - sector.sinAngle * d[1],
            +sector.sinAngle * d[0] + sector.cosAngle * d[1]
        ]);
        const dp1 = dot(U1, CmV);
        const discr1 = disk.radius * disk.radius + dp1 * dp1 - sqrLengthCmV;
        if (discr1 >= 0) {
            // The ray intersects the disk. Now test whether the sector
            // boundary segment contained by the ray overlaps the disk. The
            // quadratic root tmin generates the ray-disk point of
            // intersection closest to the sector vertex.
            const tmin = dp1 - Math.sqrt(discr1);
            if (sector.radius >= tmin) {
                result.intersect = true;
                return result;
            } else {
                // The segment does not overlap the disk. We know the disks
                // overlap, so if the disk center is outside the sector cone
                // or on the left-boundary ray, the overlap occurs outside the
                // cone, which implies the disk and sector do not intersect.
                if (dotDirCmV <= lengthCmV * sector.cosAngle) {
                    // The disk center is not inside the sector cone.
                    result.intersect = false;
                    return result;
                }
            }
        }

        // To reach here, a strict subset of the sector arc boundary must
        // intersect the disk.
        result.intersect = true;
        return result;
    }
}
