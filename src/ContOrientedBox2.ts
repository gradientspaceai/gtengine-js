// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContOrientedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for oriented boxes in 2D.
//
// Port notes: see ContAlignedBox.ts for the Cont* naming precedent
// (getContainerOrientedBox2, inContainerOrientedBox2,
// mergeContainersOrientedBox2, with the container type as a suffix so the
// exported names are globally unique). Upstream's second GetContainer
// overload takes a std::vector and forwards to the pointer overload; the port
// has a single function taking an array. GetContainer returns 'false' when
// the Gaussian fit fails, so the port returns null in that case (matching
// getContainerAlignedBox). MergeContainers returns a vestigial 'true', so the
// port returns the merged box. Upstream copies the fitted OrientedBox2 by
// value; ApprGaussian2::getParameters returns the fitter's own object, so the
// port clones it.

import { ApprGaussian2 } from './ApprGaussian2';
import { logAssert } from './Logger';
import { OrientedBox, type OrientedBox2 } from './OrientedBox';
import { Vector, add, dot, mul, negate, normalize, sub } from './Vector';
import { perp } from './Vector2';

// Compute an oriented bounding box of the points. The box center is the
// average of the points. The box axes are the eigenvectors of the covariance
// matrix.
// NOTE (upstream issue #106, the Cont* empty-input family): with the default
// ApprQuery index validation disabled, ApprGaussian2::Fit(0, nullptr) reports
// success and divides by zero, after which upstream reads points[0] out of
// bounds. The port requires at least one point, matching the guards in
// ContCircle2.ts and ContSphere3.ts.
export function getContainerOrientedBox2(points: readonly Vector[]): OrientedBox2 | null {
    logAssert(points.length > 0, 'getContainerOrientedBox2: no points.');
    for (const point of points) {
        logAssert(point.size === 2, 'getContainerOrientedBox2: points must be 2D.');
    }

    // Fit the points with a Gaussian distribution.
    const fitter = new ApprGaussian2();
    if (!fitter.fit(points)) {
        return null;
    }

    const box = fitter.getParameters().clone();

    // Let C be the box center and let U0 and U1 be the box axes. Each input
    // point is of the form X = C + y0*U0 + y1*U1. The following code computes
    // min(y0), max(y0), min(y1), and max(y1). The box center is then adjusted
    // to be
    //   C' = C + 0.5*(min(y0)+max(y0))*U0 + 0.5*(min(y1)+max(y1))*U1
    let diff = sub(points[0], box.center);
    const pmin = [dot(diff, box.axis[0]), dot(diff, box.axis[1])];
    const pmax = [pmin[0], pmin[1]];
    for (let i = 1; i < points.length; ++i) {
        diff = sub(points[i], box.center);
        for (let j = 0; j < 2; ++j) {
            const d = dot(diff, box.axis[j]);
            if (d < pmin[j]) {
                pmin[j] = d;
            }
            else if (d > pmax[j]) {
                pmax[j] = d;
            }
        }
    }

    for (let j = 0; j < 2; ++j) {
        box.center = add(box.center, mul(0.5 * (pmin[j] + pmax[j]), box.axis[j]));
        box.extent.values[j] = 0.5 * (pmax[j] - pmin[j]);
    }
    return box;
}

// Test for containment. Let X = C + y0*U0 + y1*U1 where C is the box center
// and U0 and U1 are the orthonormal axes of the box. X is in the box when
// |y_i| <= E_i for all i, where E_i are the extents of the box.
export function inContainerOrientedBox2(point: Vector, box: OrientedBox2): boolean {
    logAssert(point.size === 2 && box.dimension === 2,
        'inContainerOrientedBox2: inputs must be 2D.');

    const diff = sub(point, box.center);
    for (let i = 0; i < 2; ++i) {
        const coeff = dot(diff, box.axis[i]);
        if (Math.abs(coeff) > box.extent.values[i]) {
            return false;
        }
    }
    return true;
}

// Construct an oriented box that contains two other oriented boxes. The
// result is not guaranteed to be the minimum area box containing the input
// boxes.
export function mergeContainersOrientedBox2(box0: OrientedBox2,
    box1: OrientedBox2): OrientedBox2 {
    logAssert(box0.dimension === 2 && box1.dimension === 2,
        'mergeContainersOrientedBox2: inputs must be 2D.');

    const merge = new OrientedBox(2);

    // The first guess at the box center. This value will be updated later
    // after the input box vertices are projected onto axes determined by an
    // average of box axes.
    merge.center = mul(0.5, add(box0.center, box1.center));

    // The merged box axes are the averages of the input box axes. The axes of
    // the second box are negated, if necessary, so they form acute angles
    // with the axes of the first box.
    if (dot(box0.axis[0], box1.axis[0]) >= 0) {
        merge.axis[0] = mul(0.5, add(box0.axis[0], box1.axis[0]));
    }
    else {
        merge.axis[0] = mul(0.5, sub(box0.axis[0], box1.axis[0]));
    }
    normalize(merge.axis[0]);
    merge.axis[1] = negate(perp(merge.axis[0]));

    // Project the input box vertices onto the merged-box axes. Each axis D[i]
    // containing the current center C has a minimum projected value min[i]
    // and a maximum projected value max[i]. The corresponding endpoints on
    // the axes are C+min[i]*D[i] and C+max[i]*D[i]. The point C is not
    // necessarily the midpoint for any of the intervals. The actual box
    // center will be adjusted from C to a point C' that is the midpoint of
    // each interval,
    //   C' = C + sum_{i=0}^1 0.5*(min[i]+max[i])*D[i]
    // The box extents are
    //   e[i] = 0.5*(max[i]-min[i])
    //
    // NOTE: Upstream seeds pmin and pmax with zero rather than with the
    // projections of the first vertex. This is valid only when the current
    // center C projects into each interval, which holds because C is the
    // midpoint of the two input box centers. The port preserves this.
    const pmin = [0, 0];
    const pmax = [0, 0];

    for (const box of [box0, box1]) {
        const vertex = box.getVertices();
        for (let i = 0; i < 4; ++i) {
            const diff = sub(vertex[i], merge.center);
            for (let j = 0; j < 2; ++j) {
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
    for (let j = 0; j < 2; ++j) {
        merge.center = add(merge.center,
            mul(half * (pmax[j] + pmin[j]), merge.axis[j]));
        merge.extent.values[j] = half * (pmax[j] - pmin[j]);
    }

    return merge;
}
