// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContEllipsoid3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for ellipsoids in 3D.
//
// Port notes: see ContAlignedBox.ts and ContOrientedBox2.ts for the Cont*
// naming precedent (getContainerEllipsoid3, inContainerEllipsoid3,
// mergeContainersEllipsoid3, with the container type as a suffix so the
// exported names are globally unique). GetContainer returns 'false' when the
// Gaussian fit fails, so the port returns null in that case. MergeContainers
// returns a vestigial 'true', so the port returns the merged ellipsoid.
// 'Rotation<3,Real>(matrix)' becomes 'Rotation.fromMatrix(m).toQuaternion()'
// and 'Rotation<3,Real>(quaternion)' becomes
// 'Rotation.fromQuaternion(q).toMatrix()'.

import { ApprGaussian3 } from './ApprGaussian3';
import { Hyperellipsoid, type Ellipsoid3 } from './Hyperellipsoid';
import { Line } from './Line';
import { logAssert } from './Logger';
import { Matrix } from './Matrix';
import { projectEllipsoid3 } from './Projection';
import { addQuaternion, negateQuaternion } from './Quaternion';
import { Rotation } from './Rotation';
import { Vector, add, dot, length, mul, normalize, sub } from './Vector';

// The input points are fit with a Gaussian distribution. The center C of the
// ellipsoid is chosen to be the mean of the distribution. The axes of the
// ellipsoid are chosen to be the eigenvectors of the covariance matrix M. The
// shape of the ellipsoid is determined by the absolute values of the
// eigenvalues. NOTE: The construction is ill-conditioned if the points are
// (nearly) collinear or (nearly) planar. In this case M has a (nearly) zero
// eigenvalue, so inverting M is problematic.
//
// NOTE (upstream issue #106, the Cont* empty-input family): the port requires
// at least one point, matching the guards in the other Cont* fitters.
export function getContainerEllipsoid3(points: readonly Vector[]):
    Ellipsoid3 | null {
    logAssert(points.length > 0, 'getContainerEllipsoid3: no points.');
    for (const point of points) {
        logAssert(point.size === 3,
            'getContainerEllipsoid3: points must be 3D.');
    }

    // Fit the points with a Gaussian distribution. The covariance matrix is
    // M = sum_j D[j]*U[j]*U[j]^T, where D[j] are the eigenvalues and U[j] are
    // corresponding unit-length eigenvectors.
    const fitter = new ApprGaussian3();
    if (!fitter.fit(points)) {
        return null;
    }

    // Upstream copies the fitted OrientedBox3 by value; getParameters returns
    // the fitter's own object, so the port clones it before modifying it.
    const box = fitter.getParameters().clone();

    // If any eigenvalue is nonpositive, adjust the D[] values so that we
    // actually build an ellipsoid.
    for (let j = 0; j < 3; ++j) {
        if (box.extent.values[j] < 0) {
            box.extent.values[j] = -box.extent.values[j];
        }
    }

    // Grow the ellipsoid, while retaining its shape determined by the
    // covariance matrix, to enclose all the input points. The quadratic form
    // that is used for the ellipsoid construction is
    //   Q(X) = (X-C)^T*M*(X-C)
    //        = (X-C)^T*(sum_j D[j]*U[j]*U[j]^T)*(X-C)
    //        = sum_j D[j]*Dot(U[j],X-C)^2
    // If the maximum value of Q(X[i]) for all input points is V^2, then a
    // bounding ellipsoid is Q(X) = V^2 since Q(X[i]) <= V^2 for all i.
    let maxValue = 0;
    for (let i = 0; i < points.length; ++i) {
        const diff = sub(points[i], box.center);
        let value = 0;
        for (let j = 0; j < 3; ++j) {
            const d = dot(box.axis[j], diff);
            value += box.extent.values[j] * d * d;
        }
        if (value > maxValue) {
            maxValue = value;
        }
    }

    // Arrange for the quadratic to satisfy Q(X) <= 1.
    const ellipsoid = new Hyperellipsoid(3);
    ellipsoid.center = box.center.clone();
    for (let j = 0; j < 3; ++j) {
        ellipsoid.axis[j] = box.axis[j].clone();
        ellipsoid.extent.values[j] = Math.sqrt(maxValue / box.extent.values[j]);
    }
    return ellipsoid;
}

// Test for containment of a point inside an ellipsoid.
export function inContainerEllipsoid3(point: Vector,
    ellipsoid: Ellipsoid3): boolean {
    logAssert(point.size === 3 && ellipsoid.dimension === 3,
        'inContainerEllipsoid3: inputs must be 3D.');

    const diff = sub(point, ellipsoid.center);
    const standardized = Vector.fromArray([
        dot(diff, ellipsoid.axis[0]) / ellipsoid.extent.values[0],
        dot(diff, ellipsoid.axis[1]) / ellipsoid.extent.values[1],
        dot(diff, ellipsoid.axis[2]) / ellipsoid.extent.values[2]
    ]);
    return length(standardized) <= 1;
}

// Construct a bounding ellipsoid for the two input ellipsoids. The result is
// not necessarily the minimum-volume ellipsoid containing the two ellipsoids.
export function mergeContainersEllipsoid3(ellipsoid0: Ellipsoid3,
    ellipsoid1: Ellipsoid3): Ellipsoid3 {
    logAssert(ellipsoid0.dimension === 3 && ellipsoid1.dimension === 3,
        'mergeContainersEllipsoid3: inputs must be 3D.');

    const merge = new Hyperellipsoid(3);

    // Compute the average of the input centers.
    merge.center = mul(0.5, add(ellipsoid0.center, ellipsoid1.center));

    // The bounding ellipsoid orientation is the average of the input
    // orientations. The axes, viewed as the columns of a matrix, form a
    // rotation matrix; the matrices are converted to quaternions, averaged,
    // normalized (a slerp with t = 1/2) and converted back.
    const rot0 = Matrix.zero(3, 3);
    const rot1 = Matrix.zero(3, 3);
    for (let j = 0; j < 3; ++j) {
        rot0.setCol(j, ellipsoid0.axis[j]);
        rot1.setCol(j, ellipsoid1.axis[j]);
    }
    const q0 = Rotation.fromMatrix(rot0).toQuaternion();
    let q1 = Rotation.fromMatrix(rot1).toQuaternion();
    if (dot(q0, q1) < 0) {
        q1 = negateQuaternion(q1);
    }

    const q = addQuaternion(q0, q1);
    normalize(q);
    const rot = Rotation.fromQuaternion(q).toMatrix();
    for (let j = 0; j < 3; ++j) {
        merge.axis[j] = rot.getCol(j);
    }

    // Project the input ellipsoids onto the axes obtained by the average of
    // the orientations and that go through the center obtained by the average
    // of the centers.
    for (let i = 0; i < 3; ++i) {
        // Projection axis.
        const line = Line.fromOriginDirection(merge.center, merge.axis[i]);

        // Project the ellipsoids onto the axis.
        const p0 = projectEllipsoid3(ellipsoid0, line);
        const p1 = projectEllipsoid3(ellipsoid1, line);

        // Determine the smallest interval containing the projected intervals.
        const maxIntr = p0.smax >= p1.smax ? p0.smax : p1.smax;
        const minIntr = p0.smin <= p1.smin ? p0.smin : p1.smin;

        // Update the average center to be the center of the bounding box
        // defined by the projected intervals.
        merge.center = add(merge.center,
            mul(line.direction, 0.5 * (minIntr + maxIntr)));

        // Compute the extents of the box based on the new center.
        merge.extent.values[i] = 0.5 * (maxIntr - minIntr);
    }

    return merge;
}
