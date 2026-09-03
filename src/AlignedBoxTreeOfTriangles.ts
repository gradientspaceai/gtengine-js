// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AlignedBoxTreeOfTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of triangles whose bounding volumes are axis-aligned
// bounding boxes; see BVTree.ts and BVTreeOfTriangles.ts for the tree
// construction and for the queries. This class is the AlignedBoxBV-based
// counterpart of AABBBVTreeOfTriangles, which uses the older upstream
// AABBBoundingVolume type.
//
// Port notes: see AlignedBoxTreeOfPoints.ts; the same remarks apply. The
// linear-component/triangle intersections reported by execute() keep every
// hit, including coincident ones; see the note about upstream issue #167 in
// BVTreeOfTriangles.ts.

import { AlignedBoxBV } from './AlignedBoxBV.js';
import type { BVTreeVolumeOps } from './BVTree.js';
import { BVTreeOfTriangles } from './BVTreeOfTriangles.js';

const alignedBoxBVOps: BVTreeVolumeOps<AlignedBoxBV> = {
    create: () => new AlignedBoxBV(),
    intersectLine: (P, Q, bv) => AlignedBoxBV.intersectLine(P, Q, bv),
    intersectRay: (P, Q, bv) => AlignedBoxBV.intersectRay(P, Q, bv),
    intersectSegment: (P, Q, bv) => AlignedBoxBV.intersectSegment(P, Q, bv)
};

export class AlignedBoxTreeOfTriangles extends BVTreeOfTriangles<AlignedBoxBV> {
    constructor() {
        super(alignedBoxBVOps);
    }

    // The bounding volume for the primitives' vertices. The box is seeded
    // with the first vertex of the first triangle of the range and then grown
    // over the vertices of all triangles in the range, including that first
    // triangle.
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        const initialTri = this.mTriangles[this.mPartition[i0]];
        box.min = this.mVertices[initialTri[0]].clone();
        box.max = box.min.clone();

        for (let i = i0; i <= i1; ++i) {
            const tri = this.mTriangles[this.mPartition[i]];
            for (let j = 0; j < 3; ++j) {
                const vertex = this.mVertices[tri[j]];
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

    // The bounding volume for a single triangle's vertices.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: AlignedBoxBV): void {
        const box = boundingVolume.box;
        const tri = this.mTriangles[this.mPartition[i]];
        box.min = this.mVertices[tri[0]].clone();
        box.max = box.min.clone();
        for (let j = 1; j < 3; ++j) {
            const vertex = this.mVertices[tri[j]];
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
