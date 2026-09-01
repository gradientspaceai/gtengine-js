// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrLine3Cylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The query considers the cylinder to be a solid.
//
// NOTE (upstream behavior, preserved): the query supports only finite
// cylinders. It reads cylinder.height directly and has no infinite-cylinder
// branch, so passing a cylinder for which isInfinite() is true (Cylinder3
// stores height = -1 for that state) produces meaningless results. Callers
// must pass a finite cylinder.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream provides
// only an FIQuery specialization for this pair of primitives, which becomes
// IntrLine3Cylinder3FI. The protected DoQuery member becomes the protected
// doQuery() method.

import type { Cylinder3 } from './Cylinder3';
import type { Line } from './Line';
import type { FIQuery } from './FIQuery';
import { Vector, add, dot, mul, sub } from './Vector';
import { computeOrthogonalComplement3 } from './Vector3';

// The result of IntrLine3Cylinder3FI queries. When 'intersect' is true, the
// intersection is the segment of line parameters [parameter[0],parameter[1]]
// (the two are equal for a single-point intersection) and 'point' holds the
// corresponding points on the line.
export interface IntrLine3Cylinder3FIResult {
    intersect: boolean;
    numIntersections: number;
    parameter: [number, number];
    point: [Vector, Vector];
}

// The port of the upstream FIQuery::Result default constructor.
export function defaultIntrLine3Cylinder3FIResult(): IntrLine3Cylinder3FIResult {
    return {
        intersect: false,
        numIntersections: 0,
        parameter: [0, 0],
        point: [Vector.zero(3), Vector.zero(3)]
    };
}

// Find-intersection query for a line and a solid cylinder in 3D.
export class IntrLine3Cylinder3FI implements
    FIQuery<Line, Cylinder3, IntrLine3Cylinder3FIResult> {

    find(line: Line, cylinder: Cylinder3): IntrLine3Cylinder3FIResult {
        const result = defaultIntrLine3Cylinder3FIResult();
        this.doQuery(line.origin, line.direction, cylinder, result);
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
        cylinder: Cylinder3, result: IntrLine3Cylinder3FIResult): void {
        // Create a coordinate system for the cylinder. In this system, the
        // cylinder segment center C is the origin and the cylinder axis
        // direction W is the z-axis. U and V are the other coordinate axis
        // directions. If P = x*U+y*V+z*W, the cylinder is x^2 + y^2 = r^2,
        // where r is the cylinder radius. The end caps are |z| = h/2, where h
        // is the cylinder height.
        const basis: Vector[] = [
            cylinder.axis.direction.clone(), new Vector(3), new Vector(3)];
        computeOrthogonalComplement3(1, basis);
        const W = basis[0];
        const U = basis[1];
        const V = basis[2];
        const halfHeight = 0.5 * cylinder.height;
        const rSqr = cylinder.radius * cylinder.radius;

        // Convert the incoming line origin to cylinder coordinates.
        const diff = sub(lineOrigin, cylinder.axis.origin);
        const P = [dot(U, diff), dot(V, diff), dot(W, diff)];

        // Get the z-value, in cylinder coordinates, of the incoming line's
        // unit-length direction.
        const dz = dot(W, lineDirection);
        if (Math.abs(dz) === 1) {
            // The line is parallel to the cylinder axis. Determine whether
            // the line intersects the cylinder end disks.
            const radialSqrDist = rSqr - P[0] * P[0] - P[1] * P[1];
            if (radialSqrDist >= 0) {
                // The line intersects the cylinder end disks.
                result.intersect = true;
                result.numIntersections = 2;
                if (dz > 0) {
                    result.parameter[0] = -P[2] - halfHeight;
                    result.parameter[1] = -P[2] + halfHeight;
                }
                else {
                    result.parameter[0] = P[2] - halfHeight;
                    result.parameter[1] = P[2] + halfHeight;
                }
            }
            // else: The line is outside the cylinder, no intersection.
            return;
        }

        // Convert the incoming line unit-length direction to cylinder
        // coordinates.
        const D = [dot(U, lineDirection), dot(V, lineDirection), dz];
        if (D[2] === 0) {
            // The line is perpendicular to the cylinder axis.
            if (Math.abs(P[2]) <= halfHeight) {
                // Test intersection of line P+t*D with infinite cylinder
                // x^2+y^2 = r^2. This reduces to computing the roots of a
                // quadratic equation. If P = (px,py,pz) and D = (dx,dy,dz),
                // then the quadratic equation is
                //   (dx^2+dy^2)*t^2+2*(px*dx+py*dy)*t+(px^2+py^2-r^2) = 0
                const a0 = P[0] * P[0] + P[1] * P[1] - rSqr;
                const a1 = P[0] * D[0] + P[1] * D[1];
                const a2 = D[0] * D[0] + D[1] * D[1];
                const discr = a1 * a1 - a0 * a2;
                if (discr > 0) {
                    // The line intersects the cylinder in two places.
                    result.intersect = true;
                    result.numIntersections = 2;
                    const root = Math.sqrt(discr);
                    result.parameter[0] = (-a1 - root) / a2;
                    result.parameter[1] = (-a1 + root) / a2;
                }
                else if (discr === 0) {
                    // The line is tangent to the cylinder.
                    result.intersect = true;
                    result.numIntersections = 1;
                    result.parameter[0] = -a1 / a2;
                    result.parameter[1] = result.parameter[0];
                }
                // else: The line does not intersect the cylinder.
            }
            // else: The line is outside the planes of the cylinder end disks.
            return;
        }

        // At this time, the line direction is neither parallel nor
        // perpendicular to the cylinder axis. The line must intersect both
        // planes of the end disks, the intersection with the cylinder being a
        // segment. The t-interval of the segment is [t0,t1].

        // Test for intersections with the planes of the end disks.
        const t0 = (-halfHeight - P[2]) / D[2];
        let xTmp = P[0] + t0 * D[0];
        let yTmp = P[1] + t0 * D[1];
        if (xTmp * xTmp + yTmp * yTmp <= rSqr) {
            // Plane intersection inside the bottom cylinder end disk.
            result.parameter[result.numIntersections++] = t0;
        }

        const t1 = (+halfHeight - P[2]) / D[2];
        xTmp = P[0] + t1 * D[0];
        yTmp = P[1] + t1 * D[1];
        if (xTmp * xTmp + yTmp * yTmp <= rSqr) {
            // Plane intersection inside the top cylinder end disk.
            result.parameter[result.numIntersections++] = t1;
        }

        if (result.numIntersections < 2) {
            // Test for intersection with the cylinder wall.
            const a0 = P[0] * P[0] + P[1] * P[1] - rSqr;
            const a1 = P[0] * D[0] + P[1] * D[1];
            const a2 = D[0] * D[0] + D[1] * D[1];
            const discr = a1 * a1 - a0 * a2;
            if (discr > 0) {
                const root = Math.sqrt(discr);
                let tValue = (-a1 - root) / a2;
                if (t0 <= t1) {
                    if (t0 <= tValue && tValue <= t1) {
                        result.parameter[result.numIntersections++] = tValue;
                    }
                }
                else {
                    if (t1 <= tValue && tValue <= t0) {
                        result.parameter[result.numIntersections++] = tValue;
                    }
                }

                if (result.numIntersections < 2) {
                    tValue = (-a1 + root) / a2;
                    if (t0 <= t1) {
                        if (t0 <= tValue && tValue <= t1) {
                            result.parameter[result.numIntersections++] = tValue;
                        }
                    }
                    else {
                        if (t1 <= tValue && tValue <= t0) {
                            result.parameter[result.numIntersections++] = tValue;
                        }
                    }
                }
                // else: The line intersects an end disk and the cylinder
                // wall.
            }
            else if (discr === 0) {
                const tValue = -a1 / a2;
                if (t0 <= t1) {
                    if (t0 <= tValue && tValue <= t1) {
                        result.parameter[result.numIntersections++] = tValue;
                    }
                }
                else {
                    if (t1 <= tValue && tValue <= t0) {
                        result.parameter[result.numIntersections++] = tValue;
                    }
                }
            }
            // else: The line does not intersect the cylinder wall.
        }
        // else: The line intersects both the top and bottom cylinder end
        // disks.

        if (result.numIntersections === 2) {
            result.intersect = true;
            if (result.parameter[0] > result.parameter[1]) {
                const save = result.parameter[0];
                result.parameter[0] = result.parameter[1];
                result.parameter[1] = save;
            }
        }
        else if (result.numIntersections === 1) {
            result.intersect = true;
            result.parameter[1] = result.parameter[0];
        }
    }
}
