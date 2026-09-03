// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AlignedBoxTreeOfPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of points whose bounding volumes are axis-aligned
// bounding boxes; see BVTree.ts and BVTreeOfPoints.ts for the tree
// construction and for the queries.
//
// Port notes:
//   - The upstream 'BoundingVolume' template parameter of BVTree is split by
//     the port into the instance interface 'BVTreeBoundingVolume' and the
//     static-function object 'BVTreeVolumeOps'; see BVTree.ts. The ops object
//     for AlignedBoxBV is module-private here because AlignedBoxBV.ts (ported
//     earlier) does not export one.
//   - The upstream overrides are 'protected virtual'; they are 'protected
//     override' here.
//   - The bounding volumes are modified in place by
//     computeInteriorBoundingVolume/computeLeafBoundingVolume, as upstream.
//     C++ assignment from a vertex copies, so the port clones explicitly.

import { AlignedBoxBV } from './AlignedBoxBV.js';
import type { BVTreeVolumeOps } from './BVTree.js';
import { BVTreeOfPoints } from './BVTreeOfPoints.js';

const alignedBoxBVOps: BVTreeVolumeOps<AlignedBoxBV> = {
    create: () => new AlignedBoxBV(),
    intersectLine: (P, Q, bv) => AlignedBoxBV.intersectLine(P, Q, bv),
    intersectRay: (P, Q, bv) => AlignedBoxBV.intersectRay(P, Q, bv),
    intersectSegment: (P, Q, bv) => AlignedBoxBV.intersectSegment(P, Q, bv)
};

export class AlignedBoxTreeOfPoints extends BVTreeOfPoints<AlignedBoxBV> {
    constructor() {
        super(alignedBoxBVOps);
    }

    // The bounding volume for the primitives' vertices. The box is seeded
    // with the first vertex of the range and then grown over all vertices in
    // the range, including that first vertex.
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        box.min = this.mVertices[this.mPartition[i0]].clone();
        box.max = box.min.clone();

        for (let i = i0; i <= i1; ++i) {
            const vertex = this.mVertices[this.mPartition[i]];
            for (let k = 0; k < 3; ++k) {
                if (vertex.get(k) < box.min.get(k)) {
                    box.min.set(k, vertex.get(k));
                } else if (vertex.get(k) > box.max.get(k)) {
                    box.max.set(k, vertex.get(k));
                }
            }
        }
    }

    // The bounding volume for a single point primitive is the degenerate box
    // whose minimum and maximum are that point.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        box.min = this.mVertices[this.mPartition[i]].clone();
        box.max = box.min.clone();
    }
}
