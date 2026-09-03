// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContOrientedBox3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for oriented boxes in 3D.
//
// Port notes: see ContOrientedBox2.ts for the Cont* naming precedent
// (getContainerOrientedBox3, inContainerOrientedBox3,
// mergeContainersOrientedBox3, with the container type as a suffix so the
// exported names are globally unique). Upstream's second GetContainer
// overload takes a std::vector and forwards to the pointer overload; the port
// has a single function taking an array. GetContainer returns 'false' when
// the Gaussian fit fails, so the port returns null in that case.
// MergeContainers returns a vestigial 'true', so the port returns the merged
// box.

import { ApprGaussian3 } from './ApprGaussian3.js';
import { logAssert } from './Logger.js';
import { Matrix } from './Matrix.js';
import { OrientedBox, type OrientedBox3 } from './OrientedBox.js';
import { addQuaternion, negateQuaternion } from './Quaternion.js';
import { Rotation } from './Rotation.js';
import { Vector, add, dot, mul, normalize, sub } from './Vector.js';

// Compute an oriented bounding box of the points. The box center is the
// average of the points. The box axes are the eigenvectors of the covariance
// matrix.
//
// NOTE (upstream issue #106, the Cont* empty-input family): with the default
// ApprQuery index validation disabled, ApprGaussian3::Fit(0, nullptr) reports
// success and divides by zero, after which upstream reads points[0] out of
// bounds. The port requires at least one point, matching the guard in
// ContOrientedBox2.ts.
export function getContainerOrientedBox3(points: readonly Vector[]):
    OrientedBox3 | null {
    logAssert(points.length > 0, 'getContainerOrientedBox3: no points.');
    for (const point of points) {
        logAssert(point.size === 3,
            'getContainerOrientedBox3: points must be 3D.');
    }

    // Fit the points with a Gaussian distribution.
    const fitter = new ApprGaussian3();
    if (!fitter.fit(points)) {
        return null;
    }

    // Upstream copies the fitted OrientedBox3 by value; getParameters returns
    // the fitter's own object, so the port clones it.
    const box = fitter.getParameters().clone();

    // Let C be the box center and let U0, U1 and U2 be the box axes. Each
    // input point is of the form X = C + y0*U0 + y1*U1 + y2*U2. The following
    // code computes min(y_j) and max(y_j) for each j. The box center is then
    // adjusted to be
    //   C' = C + sum_j 0.5*(min(y_j)+max(y_j))*U_j
    let diff = sub(points[0], box.center);
    const pmin = [dot(diff, box.axis[0]), dot(diff, box.axis[1]),
        dot(diff, box.axis[2])];
    const pmax = [pmin[0], pmin[1], pmin[2]];
    for (let i = 1; i < points.length; ++i) {
        diff = sub(points[i], box.center);
        for (let j = 0; j < 3; ++j) {
            const d = dot(diff, box.axis[j]);
            if (d < pmin[j]) {
                pmin[j] = d;
            }
            else if (d > pmax[j]) {
                pmax[j] = d;
            }
        }
    }

    for (let j = 0; j < 3; ++j) {
        box.center = add(box.center, mul(0.5 * (pmin[j] + pmax[j]), box.axis[j]));
        box.extent.values[j] = 0.5 * (pmax[j] - pmin[j]);
    }
    return box;
}

// Test for containment. Let X = C + y0*U0 + y1*U1 + y2*U2 where C is the box
// center and U0, U1, U2 are the orthonormal axes of the box. X is in the box
// when |y_i| <= E_i for all i, where E_i are the extents of the box.
export function inContainerOrientedBox3(point: Vector,
    box: OrientedBox3): boolean {
    logAssert(point.size === 3 && box.dimension === 3,
        'inContainerOrientedBox3: inputs must be 3D.');

    const diff = sub(point, box.center);
    for (let i = 0; i < 3; ++i) {
        const coeff = dot(diff, box.axis[i]);
        if (Math.abs(coeff) > box.extent.values[i]) {
            return false;
        }
    }
    return true;
}

// Construct an oriented box that contains two other oriented boxes. The
// result is not guaranteed to be the minimum volume box containing the input
// boxes.
export function mergeContainersOrientedBox3(box0: OrientedBox3,
    box1: OrientedBox3): OrientedBox3 {
    logAssert(box0.dimension === 3 && box1.dimension === 3,
        'mergeContainersOrientedBox3: inputs must be 3D.');

    const merge = new OrientedBox(3);

    // The first guess at the box center. This value will be updated later
    // after the input box vertices are projected onto axes determined by an
    // average of box axes.
    merge.center = mul(0.5, add(box0.center, box1.center));

    // A box's axes, when viewed as the columns of a matrix, form a rotation
    // matrix. The input box axes are converted to quaternions. The average
    // quaternion is computed, then normalized to unit length. The result is
    // the slerp of the two input quaternions with t-value of 1/2. The result
    // is converted back to a rotation matrix and its columns are selected as
    // the merged box axes.
    const rot0 = Matrix.zero(3, 3);
    const rot1 = Matrix.zero(3, 3);
    for (let j = 0; j < 3; ++j) {
        rot0.setCol(j, box0.axis[j]);
        rot1.setCol(j, box1.axis[j]);
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

    // Project the input box vertices onto the merged-box axes. Each axis D[i]
    // containing the current center C has a minimum projected value min[i]
    // and a maximum projected value max[i]. The corresponding endpoints on
    // the axes are C+min[i]*D[i] and C+max[i]*D[i]. The point C is not
    // necessarily the midpoint for any of the intervals. The actual box
    // center will be adjusted from C to a point C' that is the midpoint of
    // each interval,
    //   C' = C + sum_{i=0}^2 0.5*(min[i]+max[i])*D[i]
    // The box extents are
    //   e[i] = 0.5*(max[i]-min[i])
    //
    // NOTE: Upstream seeds pmin and pmax with zero rather than with the
    // projections of the first vertex. This is valid only when the current
    // center C projects into each interval, which holds because C is the
    // midpoint of the two input box centers. The port preserves this (the
    // same choice is made in ContOrientedBox2.ts).
    const pmin = [0, 0, 0];
    const pmax = [0, 0, 0];

    for (const box of [box0, box1]) {
        const vertex = box.getVertices();
        for (let i = 0; i < 8; ++i) {
            const diff = sub(vertex[i], merge.center);
            for (let j = 0; j < 3; ++j) {
                const d = dot(diff, merge.axis[j]);
                if (d > pmax[j]) {
                    pmax[j] = d;
                }
                else if (d < pmin[j]) {
                    pmin[j] = d;
                }
            }
        }
    }

    // [min,max] is the axis-aligned box in the coordinate system of the
    // merged box axes. Update the current box center to be the center of the
    // new box. Compute the extents based on the new center.
    const half = 0.5;
    for (let j = 0; j < 3; ++j) {
        merge.center = add(merge.center,
            mul(half * (pmax[j] + pmin[j]), merge.axis[j]));
        merge.extent.values[j] = half * (pmax[j] - pmin[j]);
    }

    return merge;
}
