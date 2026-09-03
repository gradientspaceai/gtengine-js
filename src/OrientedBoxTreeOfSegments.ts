// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OrientedBoxTreeOfSegments.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of segments whose bounding volumes are oriented
// bounding boxes; see BVTree.ts and BVTreeOfSegments.ts for the tree
// construction and for the queries.
//
// Port notes: see OrientedBoxTreeOfPoints.ts; the same remarks apply.

import { getContainerOrientedBox3 } from './ContOrientedBox3';
import { BVTreeOfSegments } from './BVTreeOfSegments';
import { OrientedBoxBV, orientedBoxBVOps } from './OrientedBoxBV';
import { Vector } from './Vector';

export class OrientedBoxTreeOfSegments extends BVTreeOfSegments<OrientedBoxBV> {
    constructor() {
        super(orientedBoxBVOps);
    }

    // The bounding volume for the endpoints of the segments in the range.
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: OrientedBoxBV): void {
        const localVertices: Vector[] = [];
        for (let i = i0; i <= i1; ++i) {
            const seg = this.mSegments[this.mPartition[i]];
            localVertices.push(this.mVertices[seg[0]]);
            localVertices.push(this.mVertices[seg[1]]);
        }

        const box = getContainerOrientedBox3(localVertices);
        if (box !== null) {
            boundingVolume.box = box;
        }
    }

    // The bounding volume for the endpoints of a single segment.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: OrientedBoxBV): void {
        const seg = this.mSegments[this.mPartition[i]];
        const localVertices: Vector[] = [
            this.mVertices[seg[0]],
            this.mVertices[seg[1]]
        ];

        const fitBox = getContainerOrientedBox3(localVertices);
        if (fitBox === null) {
            return;
        }
        boundingVolume.box = fitBox;
        const box = boundingVolume.box;

        // Numerical rounding errors in the Gaussian fit of the container
        // computation will lead to 2 extents nearly zero. Locate them and set
        // them to zero.
        let maxAbsExtent = Math.abs(box.extent.get(0));
        let maxIndex = 0;
        let absExtent = Math.abs(box.extent.get(1));
        if (absExtent > maxAbsExtent) {
            maxAbsExtent = absExtent;
            maxIndex = 1;
        }
        absExtent = Math.abs(box.extent.get(2));
        if (absExtent > maxAbsExtent) {
            maxAbsExtent = absExtent;
            maxIndex = 2;
        }

        box.extent.set((maxIndex + 1) % 3, 0);
        box.extent.set((maxIndex + 2) % 3, 0);
    }
}
