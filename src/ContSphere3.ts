// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContSphere3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for spheres in 3D. These are the fast approximations
// that upstream provides; the bounding sphere from a point set is the
// average-center sphere, which is generally larger than the true
// minimum-volume sphere (see MinimumVolumeSphere3 for the exact
// algorithm).
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent
// (getContainerSphere3, inContainerSphere3, mergeContainersSphere3, with the
// container type as a suffix so the exported names are globally unique).
// Sphere3 is the 3D alias of Hypersphere, so the port asserts that the input
// vectors have dimension 3, which the C++ Vector3<Real> type guarantees.

import { Hypersphere, type Sphere3 } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import { Vector, add, div, mul, sub, dot, length } from './Vector.js';

// Compute the smallest bounding sphere whose center is the average of the
// input points. Upstream reads points[0] unconditionally, so at least one
// point is required.
export function getContainerSphere3(points: readonly Vector[]): Sphere3 {
    logAssert(points.length > 0, 'getContainerSphere3: no points.');

    let center = points[0].clone();
    logAssert(center.size === 3, 'getContainerSphere3: points must be 3D.');
    for (let i = 1; i < points.length; ++i) {
        logAssert(points[i].size === 3, 'getContainerSphere3: points must be 3D.');
        center = add(center, points[i]);
    }
    center = div(center, points.length);

    // Upstream reuses the radius field to accumulate the maximum squared
    // radius before taking the square root.
    let radiusSqr = 0;
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], center);
        const sqr = dot(diff, diff);
        if (sqr > radiusSqr) {
            radiusSqr = sqr;
        }
    }

    const sphere = new Hypersphere(3);
    sphere.center = center;
    sphere.radius = Math.sqrt(radiusSqr);
    return sphere;
}

// Test for containment of a point inside a sphere. The boundary is part of
// the sphere.
export function inContainerSphere3(point: Vector, sphere: Sphere3): boolean {
    logAssert(point.size === 3 && sphere.dimension === 3,
        'inContainerSphere3: inputs must be 3D.');

    const diff = sub(point, sphere.center);
    return length(diff) <= sphere.radius;
}

// Compute the smallest bounding sphere that contains the input spheres. When
// one sphere contains the other, that sphere is the result.
export function mergeContainersSphere3(sphere0: Sphere3, sphere1: Sphere3): Sphere3 {
    logAssert(sphere0.dimension === 3 && sphere1.dimension === 3,
        'mergeContainersSphere3: inputs must be 3D.');

    const cenDiff = sub(sphere1.center, sphere0.center);
    const lenSqr = dot(cenDiff, cenDiff);
    const rDiff = sphere1.radius - sphere0.radius;
    const rDiffSqr = rDiff * rDiff;

    if (rDiffSqr >= lenSqr) {
        // One sphere contains the other.
        return (rDiff >= 0 ? sphere1 : sphere0).clone();
    }

    const merge = new Hypersphere(3);
    const len = Math.sqrt(lenSqr);
    if (len > 0) {
        const coeff = (len + rDiff) / (2 * len);
        merge.center = add(sphere0.center, mul(coeff, cenDiff));
    } else {
        merge.center = sphere0.center.clone();
    }

    merge.radius = 0.5 * (len + sphere0.radius + sphere1.radius);
    return merge;
}

