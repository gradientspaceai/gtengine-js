// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContAlignedBox2Arc2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the smallest-area axis-aligned box containing an arc. Let the arc
// have endpoints E[0] and E[1] and live on a circle with center C and radius
// r. The extreme circle points in the axis directions are P[0] = C+(r,0),
// P[1] = C-(r,0), P[2] = C+(0,r) and P[3] = C-(0,r). The box is supported by
// E0, E1 and the points P[i] that are on the arc.
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent (each
// GetContainer/InContainer/MergeContainers overload is suffixed so the
// exported names are globally unique). This file's function is named after
// the header, getContainerAlignedBox2Arc2, because getContainerAlignedBox is
// already taken by the point-set overload in ContAlignedBox.ts. Upstream
// returns 'bool' and fills an output reference; the port returns the box.

import { AlignedBox } from './AlignedBox';
import { Arc2 } from './Arc2';
import { Vector, computeExtremes } from './Vector';

// Compute the smallest-area axis-aligned box containing the arc.
//
// NOTE: Upstream tests '0 < numPoints && numPoints <= 6' and has an
// unreachable 'return false' after it. numPoints starts at 2 (the arc
// endpoints) and is incremented at most four times, so the test is always
// true. The port drops the dead branch and always returns a box.
export function getContainerAlignedBox2Arc2(arc: Arc2): AlignedBox {
    // Store the arc endpoints.
    const points: Vector[] = [arc.end[0].clone(), arc.end[1].clone()];

    // Store the circle points that are on the arc. The upstream call is the
    // single-argument Arc2::Contains, which assumes the point is on the
    // circle; the port names that overload containsOnCircle.
    const candidates: Vector[] = [
        Vector.fromArray([arc.center.values[0] + arc.radius, arc.center.values[1]]),
        Vector.fromArray([arc.center.values[0] - arc.radius, arc.center.values[1]]),
        Vector.fromArray([arc.center.values[0], arc.center.values[1] + arc.radius]),
        Vector.fromArray([arc.center.values[0], arc.center.values[1] - arc.radius])
    ];

    for (const candidate of candidates) {
        if (arc.containsOnCircle(candidate)) {
            points.push(candidate);
        }
    }

    // Compute the aligned bounding box. There are always at least the two
    // arc endpoints, so computeExtremes cannot fail.
    const extremes = computeExtremes(points);
    if (extremes === null) {
        throw new Error('getContainerAlignedBox2Arc2: unreachable.');
    }
    return AlignedBox.fromMinMax(extremes.vmin, extremes.vmax);
}
