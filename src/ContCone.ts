// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContCone.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment query for a point in a cone.
//
// Port notes: the Cont* naming precedent suffixes the query with the
// container type, so InContainer becomes inContainerCone. The upstream
// template is dimension-generic (Cone<N,Real>), and the port's Cone carries
// its dimension at runtime, so a single function serves every N.

import { Cone } from './Cone.js';
import { logAssert } from './Logger.js';
import { Vector, dot, sub } from './Vector.js';

// Test for containment of a point by a cone. The cone axis direction must be
// unit length and the cone angle must have been set (see Cone.setAngle),
// because the test uses the derived constant cosAngleSqr.
export function inContainerCone(point: Vector, cone: Cone): boolean {
    logAssert(point.size === cone.dimension,
        'inContainerCone: mismatched dimensions.');

    const diff = sub(point, cone.ray.origin);
    const h = dot(cone.ray.direction, diff);
    return cone.heightInRange(h) && h * h >= cone.cosAngleSqr * dot(diff, diff);
}
