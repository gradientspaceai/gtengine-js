// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContCylinder3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for cylinders in 3D.
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent
// (getContainerCylinder3, inContainerCylinder3, with the container type as a
// suffix so the exported names are globally unique). Upstream returns a
// vestigial 'true' and fills an output reference; the port returns the
// cylinder. The number-of-points argument is dropped in favor of the array
// length.

import { ApprOrthogonalLine3 } from './ApprOrthogonalLine3.js';
import { Cylinder3 } from './Cylinder3.js';
import { DistPointLine } from './DistPointLine.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, length, mul, sub } from './Vector.js';

// Compute the cylinder axis segment using least-squares fit. The radius is
// the maximum distance from points to the axis. The height is determined by
// projection of points onto the axis and determining the containing
// interval.
//
// Upstream reads points[0] unconditionally, so at least one point is
// required.
export function getContainerCylinder3(points: readonly Vector[]): Cylinder3 {
    logAssert(points.length > 0, 'getContainerCylinder3: no points.');
    for (const point of points) {
        logAssert(point.size === 3, 'getContainerCylinder3: points must be 3D.');
    }

    const fitter = new ApprOrthogonalLine3();
    fitter.fit(points);
    const line = fitter.getParameters();

    const plQuery = new DistPointLine();
    let maxRadiusSqr = 0;
    for (let i = 0; i < points.length; ++i) {
        const result = plQuery.compute(points[i], line);
        if (result.sqrDistance > maxRadiusSqr) {
            maxRadiusSqr = result.sqrDistance;
        }
    }

    let diff = sub(points[0], line.origin);
    let wMin = dot(line.direction, diff);
    let wMax = wMin;
    for (let i = 1; i < points.length; ++i) {
        diff = sub(points[i], line.origin);
        const w = dot(line.direction, diff);
        if (w < wMin) {
            wMin = w;
        }
        else if (w > wMax) {
            wMax = w;
        }
    }

    const cylinder = new Cylinder3();
    cylinder.axis.origin = add(line.origin, mul(0.5 * (wMax + wMin), line.direction));
    cylinder.axis.direction = line.direction.clone();
    cylinder.radius = Math.sqrt(maxRadiusSqr);
    cylinder.height = wMax - wMin;
    return cylinder;
}

// Test for containment of a point by a cylinder.
export function inContainerCylinder3(point: Vector, cylinder: Cylinder3): boolean {
    logAssert(point.size === 3, 'inContainerCylinder3: point must be 3D.');

    const diff = sub(point, cylinder.axis.origin);
    const zProj = dot(diff, cylinder.axis.direction);
    if (Math.abs(zProj) * 2 > cylinder.height) {
        return false;
    }

    const xyProj = sub(diff, mul(zProj, cylinder.axis.direction));
    return length(xyProj) <= cylinder.radius;
}
