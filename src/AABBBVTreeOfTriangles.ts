// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AABBBVTreeOfTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A bounding volume tree of triangles whose bounding volumes are
// axis-aligned bounding boxes. The splitting axis at a node is the coordinate
// axis of largest box extent.
//
// Port notes:
//   - The upstream 'BoundingVolume' template parameter carries both instance
//     requirements (GetSplittingAxis) and static requirements (default
//     construction, IntersectLine/IntersectRay/IntersectSegment). The port of
//     BVTree splits these into the 'BVTreeBoundingVolume' interface and a
//     'BVTreeVolumeOps' object passed to the base constructor; see BVTree.ts.
//     AABBBoundingVolume implements the former and keeps the three upstream
//     static functions as static methods; 'aabbBoundingVolumeOps' bundles them
//     with default construction for the base class.
//   - 'GetSplittingAxis(origin, direction)' has output reference parameters
//     upstream; the port returns a 'BVTreeSplittingAxis' object, as the base
//     interface requires.
//   - The bounding volumes are modified in place by
//     computeInteriorBoundingVolume/computeLeafBoundingVolume, as upstream.
//     C++ assignments from a vertex copy, so the port clones explicitly.

import { AlignedBox } from './AlignedBox.js';
import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from './BVTree.js';
import { BVTreeOfTriangles } from './BVTreeOfTriangles.js';
import { IntrLine3AlignedBox3TI } from './IntrLine3AlignedBox3.js';
import { IntrRay3AlignedBox3TI } from './IntrRay3AlignedBox3.js';
import { IntrSegment3AlignedBox3TI } from './IntrSegment3AlignedBox3.js';
import { Line } from './Line.js';
import { Ray } from './Ray.js';
import { Segment } from './Segment.js';
import { Vector, add, mul, sub } from './Vector.js';

// The port of struct AABBBoundingVolume<T>. The default constructor is the
// upstream 'box{}', which value-initializes an AlignedBox3<T> and so runs its
// default constructor (min = (-1,-1,-1), max = (1,1,1)).
export class AABBBoundingVolume implements BVTreeBoundingVolume {
    box: AlignedBox;

    constructor() {
        this.box = new AlignedBox(3);
    }

    // The splitting axis passes through the box center in the direction of
    // the coordinate axis of largest extent. Ties are broken in favor of the
    // smaller axis index, because the comparisons are strict inequalities.
    getSplittingAxis(): BVTreeSplittingAxis {
        const zero = 0;
        const one = 1;
        const half = 0.5;

        const origin = mul(half, add(this.box.max, this.box.min));
        const extents = mul(half, sub(this.box.max, this.box.min));
        let projectionExtent = extents.get(0);
        let direction = Vector.fromArray([one, zero, zero]);
        if (extents.get(1) > projectionExtent) {
            projectionExtent = extents.get(1);
            direction = Vector.fromArray([zero, one, zero]);
        }
        if (extents.get(2) > projectionExtent) {
            projectionExtent = extents.get(2);
            direction = Vector.fromArray([zero, zero, one]);
        }

        return { origin: origin, direction: direction };
    }

    // The line is P + t * Q for all real t.
    static intersectLine(P: Vector, Q: Vector,
        boundingVolume: AABBBoundingVolume): boolean {
        const query = new IntrLine3AlignedBox3TI();
        const result = query.test(Line.fromOriginDirection(P, Q), boundingVolume.box);
        return result.intersect;
    }

    // The ray is P + t * Q for t >= 0.
    static intersectRay(P: Vector, Q: Vector,
        boundingVolume: AABBBoundingVolume): boolean {
        const query = new IntrRay3AlignedBox3TI();
        const result = query.test(Ray.fromOriginDirection(P, Q), boundingVolume.box);
        return result.intersect;
    }

    // The segment has endpoints P and Q.
    static intersectSegment(P: Vector, Q: Vector,
        boundingVolume: AABBBoundingVolume): boolean {
        const query = new IntrSegment3AlignedBox3TI();
        const result = query.test(Segment.fromEndpoints(P, Q), boundingVolume.box);
        return result.intersect;
    }
}

// The static-side requirements of the upstream BoundingVolume template
// parameter, bundled for the BVTree base class; see BVTree.ts.
export const aabbBoundingVolumeOps: BVTreeVolumeOps<AABBBoundingVolume> = {
    create: () => new AABBBoundingVolume(),
    intersectLine: (P, Q, bv) => AABBBoundingVolume.intersectLine(P, Q, bv),
    intersectRay: (P, Q, bv) => AABBBoundingVolume.intersectRay(P, Q, bv),
    intersectSegment: (P, Q, bv) => AABBBoundingVolume.intersectSegment(P, Q, bv)
};

export class AABBBVTreeOfTriangles extends BVTreeOfTriangles<AABBBoundingVolume> {
    constructor() {
        super(aabbBoundingVolumeOps);
    }

    // The bounding volume for the primitives' vertices. The box is seeded
    // with the first vertex of the first triangle of the range and then grown
    // over the vertices of all triangles in the range, including that first
    // triangle.
    protected computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: AABBBoundingVolume): void {
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
    protected computeLeafBoundingVolume(i: number,
        boundingVolume: AABBBoundingVolume): void {
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
