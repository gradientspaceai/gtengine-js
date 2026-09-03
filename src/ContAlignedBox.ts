// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ContAlignedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Containment queries for axis-aligned bounding boxes.
//
// Port notes (the Cont* naming precedent):
// - Upstream overloads the free functions GetContainer, InContainer and
//   MergeContainers across many Cont*.h headers, one set per bounding
//   volume. The port requires globally unique exported symbols, so each
//   function is suffixed with the container type it applies to:
//   getContainerAlignedBox, inContainerAlignedBox, mergeContainersAlignedBox.
// - Upstream returns 'bool' and fills an output reference. The port returns
//   the constructed container instead (or null where upstream's 'true' is
//   not guaranteed, namely an empty point set).
// - The number-of-points argument is dropped in favor of the array length.

import { AlignedBox } from './AlignedBox.js';
import { logAssert } from './Logger.js';
import { Vector, computeExtremes } from './Vector.js';

// Compute the minimum size aligned bounding box of the points. The extreme
// values are the minima and maxima of the point coordinates. Upstream
// returns the value of ComputeExtremes, which is 'false' when there are no
// points; the port returns null in that case.
export function getContainerAlignedBox(points: readonly Vector[]): AlignedBox | null {
    const extremes = computeExtremes(points);
    if (extremes === null) {
        return null;
    }
    return AlignedBox.fromMinMax(extremes.vmin, extremes.vmax);
}

// Test for containment. The box boundary is part of the box.
export function inContainerAlignedBox(point: Vector, box: AlignedBox): boolean {
    logAssert(point.size === box.dimension,
        'inContainerAlignedBox: mismatched dimensions.');

    for (let i = 0; i < point.size; ++i) {
        const value = point.values[i];
        if (value < box.min.values[i] || value > box.max.values[i]) {
            return false;
        }
    }
    return true;
}

// Construct an aligned box that contains two other aligned boxes. The result
// is the minimum size box containing the input boxes.
export function mergeContainersAlignedBox(box0: AlignedBox, box1: AlignedBox): AlignedBox {
    logAssert(box0.dimension === box1.dimension,
        'mergeContainersAlignedBox: mismatched dimensions.');

    const merge = new AlignedBox(box0.dimension);
    for (let i = 0; i < box0.dimension; ++i) {
        merge.min.values[i] = Math.min(box0.min.values[i], box1.min.values[i]);
        merge.max.values[i] = Math.max(box0.max.values[i], box1.max.values[i]);
    }
    return merge;
}
