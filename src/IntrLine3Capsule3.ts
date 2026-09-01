// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the capsule to be a solid.
//
// The test-intersection query is based on a distance computation.
//
// NOTE (upstream limitation, preserved): the find query builds a coordinate
// frame from the capsule segment's centered-form direction. A zero-length
// capsule segment (a sphere) has no such direction, the frame degenerates and
// the reported parameters are meaningless. The distance-based test query is
// unaffected. Callers of the find query must supply a nondegenerate capsule
// segment.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. The two upstream
// specializations become IntrLine3Capsule3TI and IntrLine3Capsule3FI. The
// protected DoQuery member becomes the protected doQuery() method.

import type { Capsule } from './Capsule';
import type { Line } from './Line';
import type { TIQuery } from './TIQuery';
import type { FIQuery } from './FIQuery';
import { Vector, add, dot, mul, sub } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';
import { DistLineSegment } from './DistLineSegment';

// The result of IntrLine3Capsule3TI queries.
export interface IntrLine3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrLine3Capsule3TIResult {
    return { intersect: false };
}

// The result of IntrLine3Capsule3FI queries. When 'intersect' is true, the
// intersection is the segment of line parameters [parameter[0],parameter[1]]
// (the two are equal for a single-point intersection) and 'point' holds the
// corresponding points on the line.
export interface IntrLine3Capsule3FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Capsule3FIResult(): IntrLine3Capsule3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(3), Vector.zero(3)]
    };
}

// Test-intersection query for a line and a solid capsule in 3D.
export class IntrLine3Capsule3TI implements
    TIQuery<Line, Capsule, IntrLine3Capsule3TIResult> {

    test(line: Line, capsule: Capsule): IntrLine3Capsule3TIResult {
        const result = defaultTIResult();
        const lsQuery = new DistLineSegment();
        const lsResult = lsQuery.compute(line, capsule.segment);
        result.intersect = (lsResult.distance <= capsule.radius);
        return result;
    }
}

// Find-intersection query for a line and a solid capsule in 3D.
export class IntrLine3Capsule3FI implements
    FIQuery<Line, Capsule, IntrLine3Capsule3FIResult> {

    find(line: Line, capsule: Capsule): IntrLine3Capsule3FIResult {
        const result = defaultIntrLine3Capsule3FIResult();
        this.doQuery(line.origin, line.direction, capsule, result);
        if (result.intersect) {
            for (let i = 0; i < 2; ++i) {
                result.point[i] = add(line.origin,
                    mul(result.parameter[i], line.direction));
            }
        }
        return result;
    }

    // The caller must ensure that on entry, 'result' is default constructed
    // as if there is no intersection. If an intersection is found, the
    // 'result' values are modified accordingly.
    protected doQuery(lineOrigin: Vector, lineDirection: Vector,
        capsule: Capsule, result: IntrLine3Capsule3FIResult): void {
        // Create a coordinate system for the capsule. In this system, the
        // capsule segment center C is the origin and the capsule axis
        // direction W is the z-axis. U and V are the other coordinate axis
        // directions. If P = x*U+y*V+z*W, the cylinder containing the capsule
        // wall is x^2 + y^2 = r^2, where r is the capsule radius. The finite
        // cylinder that makes up the capsule minus its hemispherical end caps
        // has z-values |z| <= e, where e is the extent of the capsule segment.
        // The top hemisphere cap is x^2+y^2+(z-e)^2 = r^2 for z >= e and the
        // bottom hemisphere cap is x^2+y^2+(z+e)^2 = r^2 for z <= -e.
        const {
            center: segOrigin,      // P
            direction: segDirection,  // D
            extent: segExtent       // e
        } = capsule.segment.getCenteredForm();

        // {W, U, V}
        const basis: Vector[] = [segDirection.clone(), new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const rSqr = capsule.radius * capsule.radius;
        const W = basis[0];
        const U = basis[1];
        const V = basis[2];

        // Convert incoming line origin to capsule coordinates.
        const diff = sub(lineOrigin, segOrigin);
        const P = [dot(U, diff), dot(V, diff), dot(W, diff)];

        // Get the z-value, in capsule coordinates, of the incoming line's
        // unit-length direction.
        const dz = dot(W, lineDirection);
        if (Math.abs(dz) === 1) {
            // The line is parallel to the capsule axis. Determine whether the
            // line intersects the capsule hemispheres.
            const radialSqrDist = rSqr - P[0] * P[0] - P[1] * P[1];
            if (radialSqrDist >= 0) {
                // The line intersects the hemispherical caps.
                result.intersect = true;
                result.numIntersections = 2;
                const zOffset = Math.sqrt(radialSqrDist) + segExtent;
                if (dz > 0) {
                    result.parameter[0] = -P[2] - zOffset;
                    result.parameter[1] = -P[2] + zOffset;
                }
                else {
                    result.parameter[0] = P[2] - zOffset;
                    result.parameter[1] = P[2] + zOffset;
                }
            }
            // else: The line is outside the capsule's cylinder, no
            // intersection.
            return;
        }

        // Convert the incoming line unit-length direction to capsule
        // coordinates.
        const D = [dot(U, lineDirection), dot(V, lineDirection), dz];

        // Test intersection of line P+t*D with infinite cylinder
        // x^2+y^2 = r^2. This reduces to computing the roots of a quadratic
        // equation. If P = (px,py,pz) and D = (dx,dy,dz), then the quadratic
        // equation is
        //   (dx^2+dy^2)*t^2 + 2*(px*dx+py*dy)*t + (px^2+py^2-r^2) = 0
        let a0 = P[0] * P[0] + P[1] * P[1] - rSqr;
        let a1 = P[0] * D[0] + P[1] * D[1];
        const a2 = D[0] * D[0] + D[1] * D[1];
        let discr = a1 * a1 - a0 * a2;
        if (discr < 0) {
            // The line does not intersect the infinite cylinder, so it cannot
            // intersect the capsule.
            return;
        }

        let root: number, tValue: number, zValue: number;
        if (discr > 0) {
            // The line intersects the infinite cylinder in two places.
            root = Math.sqrt(discr);
            tValue = (-a1 - root) / a2;
            zValue = P[2] + tValue * D[2];
            if (Math.abs(zValue) <= segExtent) {
                result.intersect = true;
                result.parameter[result.numIntersections++] = tValue;
            }

            tValue = (-a1 + root) / a2;
            zValue = P[2] + tValue * D[2];
            if (Math.abs(zValue) <= segExtent) {
                result.intersect = true;
                result.parameter[result.numIntersections++] = tValue;
            }

            if (result.numIntersections === 2) {
                // The line intersects the capsule wall in two places.
                return;
            }
        }
        else {
            // The line is tangent to the infinite cylinder but intersects the
            // cylinder in a single point.
            tValue = -a1 / a2;
            zValue = P[2] + tValue * D[2];
            if (Math.abs(zValue) <= segExtent) {
                result.intersect = true;
                result.numIntersections = 1;
                result.parameter[0] = tValue;
                result.parameter[1] = result.parameter[0];
                return;
            }
        }

        // Test intersection with the bottom hemisphere. The quadratic
        // equation is
        //   t^2 + 2*(px*dx+py*dy+(pz+e)*dz)*t + (px^2+py^2+(pz+e)^2-r^2) = 0
        // Use the fact that currently a1 = px*dx+py*dy and
        // a0 = px^2+py^2-r^2. The leading coefficient is a2 = 1, so there is
        // no need to include it in the construction.
        const PZpE = P[2] + segExtent;
        a1 += PZpE * D[2];
        a0 += PZpE * PZpE;
        discr = a1 * a1 - a0;
        if (discr > 0) {
            root = Math.sqrt(discr);
            tValue = -a1 - root;
            zValue = P[2] + tValue * D[2];
            if (zValue <= -segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }

            tValue = -a1 + root;
            zValue = P[2] + tValue * D[2];
            if (zValue <= -segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }
        }
        else if (discr === 0) {
            tValue = -a1;
            zValue = P[2] + tValue * D[2];
            if (zValue <= -segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }
        }

        // Test intersection with the top hemisphere. The quadratic equation
        // is
        //   t^2 + 2*(px*dx+py*dy+(pz-e)*dz)*t + (px^2+py^2+(pz-e)^2-r^2) = 0
        // Use the fact that currently a1 = px*dx+py*dy+(pz+e)*dz and
        // a0 = px^2+py^2+(pz+e)^2-r^2. The leading coefficient is a2 = 1, so
        // there is no need to include it in the construction.
        a1 -= 2 * segExtent * D[2];
        a0 -= 4 * segExtent * P[2];
        discr = a1 * a1 - a0;
        if (discr > 0) {
            root = Math.sqrt(discr);
            tValue = -a1 - root;
            zValue = P[2] + tValue * D[2];
            if (zValue >= segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }

            tValue = -a1 + root;
            zValue = P[2] + tValue * D[2];
            if (zValue >= segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }
        }
        else if (discr === 0) {
            tValue = -a1;
            zValue = P[2] + tValue * D[2];
            if (zValue >= segExtent) {
                result.parameter[result.numIntersections++] = tValue;
                if (result.numIntersections === 2) {
                    result.intersect = true;
                    sortParameters(result);
                    return;
                }
            }
        }

        if (result.numIntersections === 1) {
            result.parameter[1] = result.parameter[0];
        }
    }
}

// The port of the repeated upstream 'if (parameter[0] > parameter[1]) swap'.
function sortParameters(result: IntrLine3Capsule3FIResult): void {
    if (result.parameter[0] > result.parameter[1]) {
        const save = result.parameter[0];
        result.parameter[0] = result.parameter[1];
        result.parameter[1] = save;
    }
}
