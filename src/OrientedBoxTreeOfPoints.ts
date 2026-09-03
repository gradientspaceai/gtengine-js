// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OrientedBoxTreeOfPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of points whose bounding volumes are oriented
// bounding boxes; see BVTree.ts and BVTreeOfPoints.ts for the tree
// construction and for the queries.
//
// Port notes:
//   - The upstream 'BoundingVolume' template parameter of BVTree is split by
//     the port into the instance interface 'BVTreeBoundingVolume' and the
//     static-function object 'BVTreeVolumeOps'; see BVTree.ts.
//     OrientedBoxBV.ts exports the ops object 'orientedBoxBVOps'.
//   - The upstream overrides are 'protected virtual'; they are 'protected
//     override' here.
//   - Upstream 'GetContainer(points, box)' writes through a reference and
//     returns a bool; the port's 'getContainerOrientedBox3' returns the box
//     or null when the Gaussian fit fails. As upstream, a failed fit leaves
//     the bounding volume unchanged.

import { getContainerOrientedBox3 } from './ContOrientedBox3.js';
import { BVTreeOfPoints } from './BVTreeOfPoints.js';
import { OrientedBoxBV, orientedBoxBVOps } from './OrientedBoxBV.js';
import { Vector } from './Vector.js';

export class OrientedBoxTreeOfPoints extends BVTreeOfPoints<OrientedBoxBV> {
    constructor() {
        super(orientedBoxBVOps);
    }

    // Let C be the box center and let U0, U1 and U2 be the box axes. Each
    // input point is of the form X = C + y0*U0 + y1*U1 + y2*U2. The container
    // computation determines min(yj) and max(yj) for each j and adjusts the
    // box center to be
    //   C' = C + sum_j 0.5 * (min(yj) + max(yj)) * Uj
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: OrientedBoxBV): void {
        const localCentroids: Vector[] = [];
        for (let i = i0; i <= i1; ++i) {
            localCentroids.push(this.mCentroids[this.mPartition[i]]);
        }

        const box = getContainerOrientedBox3(localCentroids);
        if (box !== null) {
            boundingVolume.box = box;
        }
    }

    // The bounding volume for a single point primitive is the degenerate box
    // centered at that point with zero extents and the standard axes.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: OrientedBoxBV): void {
        const box = boundingVolume.box;

        box.center = this.mCentroids[this.mPartition[i]].clone();
        for (let j = 0; j < 3; ++j) {
            box.axis[j] = Vector.unit(3, j);
            box.extent.set(j, 0);
        }
    }
}
