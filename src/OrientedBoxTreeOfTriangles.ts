// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OrientedBoxTreeOfTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of triangles whose bounding volumes are oriented
// bounding boxes; see BVTree.ts and BVTreeOfTriangles.ts for the tree
// construction and for the queries.
//
// Port notes: see OrientedBoxTreeOfPoints.ts; the same remarks apply. The
// linear-component/triangle intersections reported by execute() keep every
// hit, including coincident ones; see the note about upstream issue #167 in
// BVTreeOfTriangles.ts.

import { getContainerOrientedBox3 } from './ContOrientedBox3';
import { BVTreeOfTriangles } from './BVTreeOfTriangles';
import { OrientedBoxBV, orientedBoxBVOps } from './OrientedBoxBV';
import { Vector } from './Vector';

export class OrientedBoxTreeOfTriangles extends BVTreeOfTriangles<OrientedBoxBV> {
    constructor() {
        super(orientedBoxBVOps);
    }

    // The bounding volume for the vertices of the triangles in the range.
    protected override computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: OrientedBoxBV): void {
        const localVertices: Vector[] = [];
        for (let i = i0; i <= i1; ++i) {
            const tri = this.mTriangles[this.mPartition[i]];
            localVertices.push(this.mVertices[tri[0]]);
            localVertices.push(this.mVertices[tri[1]]);
            localVertices.push(this.mVertices[tri[2]]);
        }

        const box = getContainerOrientedBox3(localVertices);
        if (box !== null) {
            boundingVolume.box = box;
        }
    }

    // The bounding volume for the vertices of a single triangle.
    protected override computeLeafBoundingVolume(i: number,
        boundingVolume: OrientedBoxBV): void {
        const tri = this.mTriangles[this.mPartition[i]];
        const localVertices: Vector[] = [
            this.mVertices[tri[0]],
            this.mVertices[tri[1]],
            this.mVertices[tri[2]]
        ];

        const fitBox = getContainerOrientedBox3(localVertices);
        if (fitBox === null) {
            return;
        }
        boundingVolume.box = fitBox;
        const box = boundingVolume.box;

        // Numerical rounding errors in the Gaussian fit of the container
        // computation will lead to 1 extent nearly zero. Locate it and set it
        // to zero.
        //
        // NOTE (upstream bug): the last comparison in
        // OrientedBoxTreeOfTriangles.h is 'absExtent > minAbsExtent', a
        // copy-and-paste slip from the "find the maximum" code in
        // OrientedBoxTreeOfSegments.h. With '>' the search selects the
        // *largest* extent whenever extent[2] is not the smallest, and the
        // leaf box is then collapsed along its longest axis, so it no longer
        // contains the triangle. The port uses '<' so the minimum is found.
        let minAbsExtent = Math.abs(box.extent.get(0));
        let minIndex = 0;
        let absExtent = Math.abs(box.extent.get(1));
        if (absExtent < minAbsExtent) {
            minAbsExtent = absExtent;
            minIndex = 1;
        }
        absExtent = Math.abs(box.extent.get(2));
        if (absExtent < minAbsExtent) {
            minAbsExtent = absExtent;
            minIndex = 2;
        }

        box.extent.set(minIndex, 0);
    }
}


