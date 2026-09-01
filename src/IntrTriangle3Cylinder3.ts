// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrTriangle3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// An algorithm for the test-intersection query between a triangle and a
// finite cylinder is described in
// https://www.geometrictools.com/Documentation/IntersectionTriangleCylinder.pdf
// The code here is an implementation of that algorithm. The comments include
// references to Figure 1 of the PDF.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrTriangle3Cylinder3TI.
// The private member helpers (DiskOverlapsPoint, DiskOverlapsSegment,
// DiskOverlapsPolygon) become module-private functions; they are stateless.
// The 2D points of the projection onto the plane perpendicular to the
// cylinder axis are Vector objects of size 2, so Vector2's dotPerp applies.
//
// Upstream does not guard against an infinite cylinder (the height = -1
// sentinel of Cylinder3), which would silently be treated as a finite
// cylinder of height -1 (an empty slab, -h/2 = 0.5 > 0.5 = h/2). Following
// the IntrCanonicalBox3Cylinder3 precedent (B69) and upstream issues
// #187/#197, the port asserts that the cylinder is finite.

import type { Cylinder3 } from './Cylinder3';
import { logAssert } from './Logger';
import type { TIQuery } from './TIQuery';
import type { Triangle } from './Triangle';
import { Vector, add, div, dot, mul, sub } from './Vector';
import { dotPerp } from './Vector2';
import { computeOrthogonalComplement3 } from './Vector3';

// The result of IntrTriangle3Cylinder3TI.test.
export interface IntrTriangle3Cylinder3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrTriangle3Cylinder3TIResult {
    return { intersect: false };
}

// Support for the static test query.
function diskOverlapsPoint(Q: Vector, radius: number): boolean {
    return dot(Q, Q) <= radius * radius;
}

function diskOverlapsSegment(Q0: Vector, Q1: Vector, radius: number): boolean {
    const sqrRadius = radius * radius;
    const direction = sub(Q0, Q1);
    let d = dot(Q0, direction);
    if (d <= 0) {
        return dot(Q0, Q0) <= sqrRadius;
    }

    const sqrLength = dot(direction, direction);
    if (d >= sqrLength) {
        return dot(Q1, Q1) <= sqrRadius;
    }

    d = dotPerp(direction, Q0);
    return d * d <= sqrLength * sqrRadius;
}

function diskOverlapsPolygon(numVertices: number, Q: readonly Vector[],
    radius: number): boolean {
    // Test whether the polygon contains (0,0).
    let positive = 0, negative = 0;
    let i0: number, i1: number;
    for (i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++) {
        const d = dotPerp(Q[i0], sub(Q[i0], Q[i1]));
        if (d > 0) {
            ++positive;
        } else if (d < 0) {
            ++negative;
        }
    }
    if (positive === 0 || negative === 0) {
        // The polygon contains (0,0), so the disk and polygon overlap.
        //
        // NOTE: When the polygon degenerates to a single point, all of the
        // DotPerp values are zero, so positive == negative == 0 and the
        // upstream code reports containment of (0,0) even when the point is
        // far from the origin. This cannot occur for a nondegenerate
        // triangle: the projection along the cylinder axis has a
        // one-dimensional kernel, so it cannot collapse a triangle to a
        // point. It does occur for a triangle degenerated to a segment
        // parallel to the cylinder axis, for which the query then reports a
        // false positive. The port preserves the upstream behavior.
        return true;
    }

    // Test whether any edge is overlapped by the polygon.
    for (i0 = numVertices - 1, i1 = 0; i1 < numVertices; i0 = i1++) {
        if (diskOverlapsSegment(Q[i0], Q[i1], radius)) {
            return true;
        }
    }

    return false;
}

export class IntrTriangle3Cylinder3TI implements
    TIQuery<Triangle, Cylinder3, IntrTriangle3Cylinder3TIResult> {

    test(triangle: Triangle, cylinder: Cylinder3):
        IntrTriangle3Cylinder3TIResult {
        logAssert(triangle.dimension === 3,
            'IntrTriangle3Cylinder3TI: mismatched sizes.');
        logAssert(cylinder.isFinite(),
            'Infinite cylinders are not yet supported.');

        const result = defaultTIResult();

        // Get a right-handed orthonormal basis from the cylinder axis
        // direction. The basis is {U2,U0,U1}.
        const basis: Vector[] = [
            cylinder.axis.direction.clone(),
            Vector.zero(3),
            Vector.zero(3)
        ];
        computeOrthogonalComplement3(1, basis);

        // Compute coordinates of the triangle vertices in the coordinate
        // system {C;U0,U1,U2}, where C is the cylinder center and U2 is the
        // cylinder direction. The basis {U0,U1,U2} is orthonormal and
        // right-handed.
        const P: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        for (let i = 0; i < 3; ++i) {
            const delta = sub(triangle.v[i], cylinder.axis.origin);
            P[i][0] = dot(basis[1], delta);  // x[i]
            P[i][1] = dot(basis[2], delta);  // y[i]
            P[i][2] = dot(basis[0], delta);  // z[i]
        }

        // Sort the triangle vertices so that z[0] <= z[1] <= z[2].
        let j0 = 0, j1 = 0, j2 = 0;
        if (P[0][2] < P[1][2]) {
            if (P[2][2] < P[0][2]) {
                j0 = 2;
                j1 = 0;
                j2 = 1;
            } else if (P[2][2] < P[1][2]) {
                j0 = 0;
                j1 = 2;
                j2 = 1;
            } else {
                j0 = 0;
                j1 = 1;
                j2 = 2;
            }
        } else {
            if (P[2][2] < P[1][2]) {
                j0 = 2;
                j1 = 1;
                j2 = 0;
            } else if (P[2][2] < P[0][2]) {
                j0 = 1;
                j1 = 2;
                j2 = 0;
            } else {
                j0 = 1;
                j1 = 0;
                j2 = 2;
            }
        }

        const z: number[] = [P[j0][2], P[j1][2], P[j2][2]];

        // Maintain the xy-components and z-components separately. The
        // z-components are used for clipping against bottom and top planes
        // of the cylinder. The xy-components are used for disk-containment
        // tests x * x + y * y <= r * r.

        // Attempt an early exit by testing whether the triangle is strictly
        // outside the cylinder slab -h/2 < z < h/2.
        const hhalf = 0.5 * cylinder.height;
        if (z[2] < -hhalf) {
            // The triangle is strictly below the bottom-disk plane of the
            // cylinder. See case 0a of Figure 1 in the PDF.
            result.intersect = false;
            return result;
        }

        if (z[0] > hhalf) {
            // The triangle is strictly above the top-disk plane of the
            // cylinder. See case 0b of Figure 1 in the PDF.
            result.intersect = false;
            return result;
        }

        // Project the triangle vertices onto the xy-plane.
        const Q: Vector[] = [
            Vector.fromArray([P[j0][0], P[j0][1]]),
            Vector.fromArray([P[j1][0], P[j1][1]]),
            Vector.fromArray([P[j2][0], P[j2][1]])
        ];

        // Attempt an early exit when the triangle does not have to be
        // clipped.
        const radius = cylinder.radius;
        if (-hhalf <= z[0] && z[2] <= hhalf) {
            // The triangle is between the planes of the top-disk and the
            // bottom disk of the cylinder. Determine whether the projection
            // of the triangle onto a plane perpendicular to the cylinder
            // axis overlaps the disk of projection of the cylinder onto the
            // same plane. See case 3a of Figure 1 of the PDF.
            result.intersect = diskOverlapsPolygon(3, Q, radius);
            return result;
        }

        // Clip against |z| <= h/2. At this point we know that z2 >= -h/2 and
        // z0 <= h/2 with either z0 < -h/2 or z2 > h/2 or both. The
        // test-intersection query involves testing for overlap between the
        // xy-projection of the clipped triangle and the xy-projection of the
        // cylinder (a disk in the projection plane). The code below computes
        // the vertices of the projection of the clipped triangle. The
        // t-values of the triangle-edge parameterizations satisfy
        // 0 <= t <= 1.
        if (z[0] < -hhalf) {
            if (z[2] > hhalf) {
                if (z[1] >= hhalf) {
                    // Cases 4a and 4b of Figure 1 in the PDF.
                    //
                    // The edge <V0,V1> is parameterized by V0+t*(V1-V0). On
                    // the bottom of the slab,
                    //   -h/2 = z0 + t * (z1 - z0)
                    //   t = (-h/2 - z0) / (z1 - z0) = numerNeg0 / denom10
                    // and on the top of the slab,
                    //   +h/2 = z0 + t * (z1 - z0)
                    //   t = (+h/2 - z0) / (z1 - z0) = numerPos0 / denom10
                    //
                    // The edge <V0,V2> is parameterized by V0+t*(V2-V0). On
                    // the bottom of the slab,
                    //   -h/2 = z0 + t * (z2 - z0)
                    //   t = (-h/2 - z0) / (z2 - z0) = numerNeg0 / denom20
                    // and on the top of the slab,
                    //   +h/2 = z0 + t * (z2 - z0)
                    //   t = (+h/2 - z0) / (z2 - z0) = numerPos0 / denom20
                    const numerNeg0 = -hhalf - z[0];
                    const numerPos0 = +hhalf - z[0];
                    const denom10 = z[1] - z[0];
                    const denom20 = z[2] - z[0];
                    const dir20 = div(sub(Q[2], Q[0]), denom20);
                    const dir10 = div(sub(Q[1], Q[0]), denom10);
                    const polygon: Vector[] = [
                        add(Q[0], mul(numerNeg0, dir20)),
                        add(Q[0], mul(numerNeg0, dir10)),
                        add(Q[0], mul(numerPos0, dir10)),
                        add(Q[0], mul(numerPos0, dir20))
                    ];
                    result.intersect = diskOverlapsPolygon(4, polygon, radius);
                } else if (z[1] <= -hhalf) {
                    // Cases 4c and 4d of Figure 1 of the PDF.
                    //
                    // The edge <V2,V0> is parameterized by V2+t*(V0-V2). On
                    // the bottom of the slab,
                    //   -h/2 = z2 + t * (z0 - z2)
                    //   t = (-h/2 - z2) / (z0 - z2) = numerNeg2 / denom02
                    // and on the top of the slab,
                    //   +h/2 = z2 + t * (z0 - z2)
                    //   t = (+h/2 - z2) / (z0 - z2) = numerPos2 / denom02
                    //
                    // The edge <V2,V1> is parameterized by V2+t*(V1-V2). On
                    // the bottom of the slab,
                    //   -h/2 = z2 + t * (z1 - z2)
                    //   t = (-h/2 - z2) / (z1 - z2) = numerNeg2 / denom12
                    // and on the top of the slab,
                    //   +h/2 = z2 + t * (z1 - z2)
                    //   t = (+h/2 - z2) / (z1 - z2) = numerPos2 / denom12
                    const numerNeg2 = -hhalf - z[2];
                    const numerPos2 = +hhalf - z[2];
                    const denom02 = z[0] - z[2];
                    const denom12 = z[1] - z[2];
                    const dir02 = div(sub(Q[0], Q[2]), denom02);
                    const dir12 = div(sub(Q[1], Q[2]), denom12);
                    const polygon: Vector[] = [
                        add(Q[2], mul(numerNeg2, dir02)),
                        add(Q[2], mul(numerNeg2, dir12)),
                        add(Q[2], mul(numerPos2, dir12)),
                        add(Q[2], mul(numerPos2, dir02))
                    ];
                    result.intersect = diskOverlapsPolygon(4, polygon, radius);
                } else {
                    // -hhalf < z[1] < hhalf. Case 5 of Figure 1 of the PDF.
                    //
                    // The edge <V0,V2> is parameterized by V0+t*(V2-V0). On
                    // the bottom of the slab,
                    //   -h/2 = z0 + t * (z2 - z0)
                    //   t = (-h/2 - z0) / (z2 - z0) = numerNeg0 / denom20
                    // and on the top of the slab,
                    //   +h/2 = z0 + t * (z2 - z0)
                    //   t = (+h/2 - z0) / (z2 - z0) = numerPos0 / denom20
                    //
                    // The edge <V1,V0> is parameterized by V1+t*(V0-V1). On
                    // the bottom of the slab,
                    //   -h/2 = z1 + t * (z0 - z1)
                    //   t = (-h/2 - z1) / (z0 - z1) = numerNeg1 / denom01
                    //
                    // The edge <V1,V2> is parameterized by V1+t*(V2-V1). On
                    // the top of the slab,
                    //   +h/2 = z1 + t * (z2 - z1)
                    //   t = (+h/2 - z1) / (z2 - z1) = numerPos1 / denom21
                    const numerNeg0 = -hhalf - z[0];
                    const numerPos0 = +hhalf - z[0];
                    const numerNeg1 = -hhalf - z[1];
                    const numerPos1 = +hhalf - z[1];
                    const denom20 = z[2] - z[0];
                    const denom01 = z[0] - z[1];
                    const denom21 = z[2] - z[1];
                    const dir20 = div(sub(Q[2], Q[0]), denom20);
                    const dir01 = div(sub(Q[0], Q[1]), denom01);
                    const dir21 = div(sub(Q[2], Q[1]), denom21);
                    const polygon: Vector[] = [
                        add(Q[0], mul(numerNeg0, dir20)),
                        add(Q[1], mul(numerNeg1, dir01)),
                        Q[1],
                        add(Q[1], mul(numerPos1, dir21)),
                        add(Q[0], mul(numerPos0, dir20))
                    ];
                    result.intersect = diskOverlapsPolygon(5, polygon, radius);
                }
            } else if (z[2] > -hhalf) {
                if (z[1] <= -hhalf) {
                    // Cases 3b and 3c of Figure 1 of the PDF.
                    //
                    // The edge <V2,V0> is parameterized by V2+t*(V0-V2). On
                    // the bottom of the slab,
                    //   -h/2 = z2 + t * (z0 - z2)
                    //   t = (-h/2 - z2) / (z0 - z2) = numerNeg2 / denom02
                    //
                    // The edge <V2,V1> is parameterized by V2+t*(V1-V2). On
                    // the bottom of the slab,
                    //   -h/2 = z2 + t * (z1 - z2)
                    //   t = (-h/2 - z2) / (z1 - z2) = numerNeg2 / denom12
                    const numerNeg2 = -hhalf - z[2];
                    const denom02 = z[0] - z[2];
                    const denom12 = z[1] - z[2];
                    const dir02 = div(sub(Q[0], Q[2]), denom02);
                    const dir12 = div(sub(Q[1], Q[2]), denom12);
                    const polygon: Vector[] = [
                        Q[2],
                        add(Q[2], mul(numerNeg2, dir02)),
                        add(Q[2], mul(numerNeg2, dir12))
                    ];
                    result.intersect = diskOverlapsPolygon(3, polygon, radius);
                } else {
                    // z[1] > -hhalf. Case 4e of Figure 1 of the PDF.
                    //
                    // The edge <V0,V1> is parameterized by V0+t*(V1-V0). On
                    // the bottom of the slab,
                    //   -h/2 = z0 + t * (z1 - z0)
                    //   t = (-h/2 - z0) / (z1 - z0) = numerNeg0 / denom10
                    //
                    // The edge <V0,V2> is parameterized by V0+t*(V2-V0). On
                    // the bottom of the slab,
                    //   -h/2 = z0 + t * (z2 - z0)
                    //   t = (-h/2 - z0) / (z2 - z0) = numerNeg0 / denom20
                    const numerNeg0 = -hhalf - z[0];
                    const denom10 = z[1] - z[0];
                    const denom20 = z[2] - z[0];
                    const dir20 = div(sub(Q[2], Q[0]), denom20);
                    const dir10 = div(sub(Q[1], Q[0]), denom10);
                    const polygon: Vector[] = [
                        add(Q[0], mul(numerNeg0, dir20)),
                        add(Q[0], mul(numerNeg0, dir10)),
                        Q[1],
                        Q[2]
                    ];
                    result.intersect = diskOverlapsPolygon(4, polygon, radius);
                }
            } else {
                // z[2] == -hhalf
                if (z[1] < -hhalf) {
                    // Case 1a of Figure 1 of the PDF.
                    result.intersect = diskOverlapsPoint(Q[2], radius);
                } else {
                    // Case 2a of Figure 1 of the PDF.
                    result.intersect = diskOverlapsSegment(Q[1], Q[2], radius);
                }
            }
        } else if (z[0] < hhalf) {
            if (z[1] >= hhalf) {
                // Cases 3d and 3e of Figure 1 of the PDF.
                //
                // The edge <V0,V1> is parameterized by V0+t*(V1-V0). On the
                // top of the slab,
                //   +h/2 = z0 + t * (z1 - z0)
                //   t = (+h/2 - z0) / (z1 - z0) = numerPos0 / denom10
                //
                // The edge <V0,V2> is parameterized by V0+t*(V2-V0). On the
                // top of the slab,
                //   +h/2 = z0 + t * (z2 - z0)
                //   t = (+h/2 - z0) / (z2 - z0) = numerPos0 / denom20
                const numerPos0 = +hhalf - z[0];
                const denom10 = z[1] - z[0];
                const denom20 = z[2] - z[0];
                const dir10 = div(sub(Q[1], Q[0]), denom10);
                const dir20 = div(sub(Q[2], Q[0]), denom20);
                const polygon: Vector[] = [
                    Q[0],
                    add(Q[0], mul(numerPos0, dir10)),
                    add(Q[0], mul(numerPos0, dir20))
                ];
                result.intersect = diskOverlapsPolygon(3, polygon, radius);
            } else {
                // z[1] < hhalf. Case 4f of Figure 1 of the PDF.
                //
                // The edge <V2,V0> is parameterized by V2+t*(V0-V2). On the
                // top of the slab,
                //   +h/2 = z2 + t * (z0 - z2)
                //   t = (+h/2 - z2) / (z0 - z2) = numerPos2 / denom02
                //
                // The edge <V2,V1> is parameterized by V2+t*(V1-V2). On the
                // top of the slab,
                //   +h/2 = z2 + t * (z1 - z2)
                //   t = (+h/2 - z2) / (z1 - z2) = numerPos2 / denom12
                const numerPos2 = +hhalf - z[2];
                const denom02 = z[0] - z[2];
                const denom12 = z[1] - z[2];
                const dir02 = div(sub(Q[0], Q[2]), denom02);
                const dir12 = div(sub(Q[1], Q[2]), denom12);
                const polygon: Vector[] = [
                    Q[0],
                    Q[1],
                    add(Q[2], mul(numerPos2, dir12)),
                    add(Q[2], mul(numerPos2, dir02))
                ];
                result.intersect = diskOverlapsPolygon(4, polygon, radius);
            }
        } else {
            // z[0] == hhalf
            if (z[1] > hhalf) {
                // Case 1b of Figure 1 of the PDF.
                result.intersect = diskOverlapsPoint(Q[0], radius);
            } else {
                // Case 2b of Figure 1 of the PDF.
                result.intersect = diskOverlapsSegment(Q[0], Q[1], radius);
            }
        }

        return result;
    }
}
