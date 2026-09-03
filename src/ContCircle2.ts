// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContCircle2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for circles in 2D. These are the fast approximations
// that upstream provides; the bounding circle from a point set is the
// average-center circle, which is generally larger than the true
// minimum-area circle (see MinimumAreaCircle2 for the exact algorithm).
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent
// (getContainerCircle2, inContainerCircle2, mergeContainersCircle2, with the
// container type as a suffix so the exported names are globally unique).
// Circle2 is the 2D alias of Hypersphere, so the port asserts that the input
// vectors have dimension 2, which the C++ Vector2<Real> type guarantees.

import { Hypersphere, type Circle2 } from './Hypersphere.js';
import { logAssert } from './Logger.js';
import { Vector, add, div, mul, sub, dot, length } from './Vector.js';

// Compute the smallest bounding circle whose center is the average of the
// input points. Upstream reads points[0] unconditionally, so at least one
// point is required.
export function getContainerCircle2(points: readonly Vector[]): Circle2 {
    logAssert(points.length > 0, 'getContainerCircle2: no points.');

    let center = points[0].clone();
    logAssert(center.size === 2, 'getContainerCircle2: points must be 2D.');
    for (let i = 1; i < points.length; ++i) {
        logAssert(points[i].size === 2, 'getContainerCircle2: points must be 2D.');
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

    const circle = new Hypersphere(2);
    circle.center = center;
    circle.radius = Math.sqrt(radiusSqr);
    return circle;
}

// Test for containment of a point inside a circle. The boundary is part of
// the circle.
export function inContainerCircle2(point: Vector, circle: Circle2): boolean {
    logAssert(point.size === 2 && circle.dimension === 2,
        'inContainerCircle2: inputs must be 2D.');

    const diff = sub(point, circle.center);
    return length(diff) <= circle.radius;
}

// Compute the smallest bounding circle that contains the input circles. When
// one circle contains the other, that circle is the result.
export function mergeContainersCircle2(circle0: Circle2, circle1: Circle2): Circle2 {
    logAssert(circle0.dimension === 2 && circle1.dimension === 2,
        'mergeContainersCircle2: inputs must be 2D.');

    const cenDiff = sub(circle1.center, circle0.center);
    const lenSqr = dot(cenDiff, cenDiff);
    const rDiff = circle1.radius - circle0.radius;
    const rDiffSqr = rDiff * rDiff;

    if (rDiffSqr >= lenSqr) {
        // One circle contains the other.
        return (rDiff >= 0 ? circle1 : circle0).clone();
    }

    const merge = new Hypersphere(2);
    const len = Math.sqrt(lenSqr);
    if (len > 0) {
        const coeff = (len + rDiff) / (2 * len);
        merge.center = add(circle0.center, mul(coeff, cenDiff));
    } else {
        merge.center = circle0.center.clone();
    }

    merge.radius = 0.5 * (len + circle0.radius + circle1.radius);
    return merge;
}

