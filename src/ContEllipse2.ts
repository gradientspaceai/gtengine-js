// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContEllipse2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for ellipses in 2D.
//
// Port notes: see ContAlignedBox.ts and ContOrientedBox2.ts for the Cont*
// naming precedent (getContainerEllipse2, inContainerEllipse2,
// mergeContainersEllipse2, with the container type as a suffix so the
// exported names are globally unique). GetContainer returns 'false' when the
// Gaussian fit fails, so the port returns null in that case. MergeContainers
// returns a vestigial 'true', so the port returns the merged ellipse.

import { ApprGaussian2 } from './ApprGaussian2.js';
import { Hyperellipsoid, type Ellipse2 } from './Hyperellipsoid.js';
import { Line } from './Line.js';
import { logAssert } from './Logger.js';
import { projectEllipse2 } from './Projection.js';
import { Vector, add, dot, length, mul, negate, normalize, sub } from './Vector.js';
import { perp } from './Vector2.js';

// The input points are fit with a Gaussian distribution. The center C of the
// ellipse is chosen to be the mean of the distribution. The axes of the
// ellipse are chosen to be the eigenvectors of the covariance matrix M. The
// shape of the ellipse is determined by the absolute values of the
// eigenvalues. NOTE: The construction is ill-conditioned if the points are
// (nearly) collinear. In this case M has a (nearly) zero eigenvalue, so
// inverting M can be a problem numerically.
//
// NOTE (upstream issue #106, the Cont* empty-input family): with the default
// ApprQuery index validation disabled, ApprGaussian2::Fit(0, nullptr) reports
// success and divides by zero. The port requires at least one point, matching
// the guards in ContCircle2.ts, ContSphere3.ts and ContOrientedBox2.ts.
export function getContainerEllipse2(points: readonly Vector[]):
    Ellipse2 | null {
    logAssert(points.length > 0, 'getContainerEllipse2: no points.');
    for (const point of points) {
        logAssert(point.size === 2, 'getContainerEllipse2: points must be 2D.');
    }

    // Fit the points with a Gaussian distribution. The covariance matrix is
    // M = sum_j D[j]*U[j]*U[j]^T, where D[j] are the eigenvalues and U[j] are
    // corresponding unit-length eigenvectors.
    const fitter = new ApprGaussian2();
    if (!fitter.fit(points)) {
        return null;
    }

    // Upstream copies the fitted OrientedBox2 by value; getParameters returns
    // the fitter's own object, so the port clones it before modifying it.
    const box = fitter.getParameters().clone();

    // If either eigenvalue is nonpositive, adjust the D[] values so that we
    // actually build an ellipse.
    for (let j = 0; j < 2; ++j) {
        if (box.extent.values[j] < 0) {
            box.extent.values[j] = -box.extent.values[j];
        }
    }

    // Grow the ellipse, while retaining its shape determined by the
    // covariance matrix, to enclose all the input points. The quadratic form
    // that is used for the ellipse construction is
    //   Q(X) = (X-C)^T*M*(X-C)
    //        = (X-C)^T*(sum_j D[j]*U[j]*U[j]^T)*(X-C)
    //        = sum_j D[j]*Dot(U[j],X-C)^2
    // If the maximum value of Q(X[i]) for all input points is V^2, then a
    // bounding ellipse is Q(X) = V^2, because Q(X[i]) <= V^2 for all i.
    let maxValue = 0;
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], box.center);
        const dot0 = dot(box.axis[0], diff);
        const dot1 = dot(box.axis[1], diff);
        const value = box.extent.values[0] * dot0 * dot0
            + box.extent.values[1] * dot1 * dot1;
        if (value > maxValue) {
            maxValue = value;
        }
    }

    // Arrange for the quadratic to satisfy Q(X) <= 1.
    const ellipse = new Hyperellipsoid(2);
    ellipse.center = box.center.clone();
    for (let j = 0; j < 2; ++j) {
        ellipse.axis[j] = box.axis[j].clone();
        ellipse.extent.values[j] = Math.sqrt(maxValue / box.extent.values[j]);
    }
    return ellipse;
}

// Test for containment of a point inside an ellipse.
export function inContainerEllipse2(point: Vector, ellipse: Ellipse2): boolean {
    logAssert(point.size === 2 && ellipse.dimension === 2,
        'inContainerEllipse2: inputs must be 2D.');

    const diff = sub(point, ellipse.center);
    const standardized = Vector.fromArray([
        dot(diff, ellipse.axis[0]) / ellipse.extent.values[0],
        dot(diff, ellipse.axis[1]) / ellipse.extent.values[1]
    ]);
    return length(standardized) <= 1;
}

// Construct a bounding ellipse for the two input ellipses. The result is not
// necessarily the minimum-area ellipse containing the two ellipses.
export function mergeContainersEllipse2(ellipse0: Ellipse2,
    ellipse1: Ellipse2): Ellipse2 {
    logAssert(ellipse0.dimension === 2 && ellipse1.dimension === 2,
        'mergeContainersEllipse2: inputs must be 2D.');

    const merge = new Hyperellipsoid(2);

    // Compute the average of the input centers.
    merge.center = mul(0.5, add(ellipse0.center, ellipse1.center));

    // The bounding ellipse orientation is the average of the input
    // orientations.
    if (dot(ellipse0.axis[0], ellipse1.axis[0]) >= 0) {
        merge.axis[0] = mul(0.5, add(ellipse0.axis[0], ellipse1.axis[0]));
    }
    else {
        merge.axis[0] = mul(0.5, sub(ellipse0.axis[0], ellipse1.axis[0]));
    }
    normalize(merge.axis[0]);
    merge.axis[1] = negate(perp(merge.axis[0]));

    // Project the input ellipses onto the axes obtained by the average of the
    // orientations and that go through the center obtained by the average of
    // the centers.
    for (let j = 0; j < 2; ++j) {
        // Projection axis.
        const line = Line.fromOriginDirection(merge.center, merge.axis[j]);

        // Project the ellipses onto the axis.
        const p0 = projectEllipse2(ellipse0, line);
        const p1 = projectEllipse2(ellipse1, line);

        // Determine the smallest interval containing the projected intervals.
        const maxIntr = p0.smax >= p1.smax ? p0.smax : p1.smax;
        const minIntr = p0.smin <= p1.smin ? p0.smin : p1.smin;

        // Update the average center to be the center of the bounding box
        // defined by the projected intervals.
        merge.center = add(merge.center,
            mul(line.direction, 0.5 * (minIntr + maxIntr)));

        // Compute the extents of the box based on the new center.
        merge.extent.values[j] = 0.5 * (maxIntr - minIntr);
    }

    return merge;
}
