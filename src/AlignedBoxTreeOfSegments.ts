// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AlignedBoxTreeOfSegments.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of segments whose bounding volumes are axis-aligned
// bounding boxes; see BVTree.ts and BVTreeOfSegments.ts for the tree
// construction and for the queries.
//
// Port notes: see AlignedBoxTreeOfPoints.ts; the same remarks apply.

import { AlignedBoxBV } from './AlignedBoxBV.js';
import type { BVTreeVolumeOps } from './BVTree.js';
import { BVTreeOfSegments } from './BVTreeOfSegments.js';

const alignedBoxBVOps: BVTreeVolumeOps<AlignedBoxBV> = {
    create: () => new AlignedBoxBV(),
    intersectLine: (P, Q, bv) => AlignedBoxBV.intersectLine(P, Q, bv),
    intersectRay: (P, Q, bv) => AlignedBoxBV.intersectRay(P, Q, bv),
    intersectSegment: (P, Q, bv) => AlignedBoxBV.intersectSegment(P, Q, bv)
};

export class AlignedBoxTreeOfSegments extends BVTreeOfSegments<AlignedBoxBV> {
    constructor() {
        super(alignedBoxBVOps);
    }

    // The bounding volume for the primitives' vertices. The box is seeded
    // with the first endpoint of the first segment of the range and then
    // grown over the endpoints of all segments in the range, including that
    // first segment.
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        const initialSeg = this.mSegments[this.mPartition[i0]];
        box.min = this.mVertices[initialSeg[0]].clone();
        box.max = box.min.clone();

        for (let i = i0; i <= i1; ++i) {
            const seg = this.mSegments[this.mPartition[i]];
            for (let j = 0; j < 2; ++j) {
                const vertex = this.mVertices[seg[j]];
                for (let k = 0; k < 3; ++k) {
                    if (vertex.get(k) < box.min.get(k)) {
                        box.min.set(k, vertex.get(k));
                    } else if (vertex.get(k) > box.max.get(k)) {
                        box.max.set(k, vertex.get(k));
                    }
                }
            }
        }
    }

    // The bounding volume for a single segment's endpoints.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        const seg = this.mSegments[this.mPartition[i]];
        box.min = this.mVertices[seg[0]].clone();
        box.max = box.min.clone();

        const vertex = this.mVertices[seg[1]];
        for (let k = 0; k < 3; ++k) {
            if (vertex.get(k) < box.min.get(k)) {
                box.min.set(k, vertex.get(k));
            } else if (vertex.get(k) > box.max.get(k)) {
                box.max.set(k, vertex.get(k));
            }
        }
    }
}
