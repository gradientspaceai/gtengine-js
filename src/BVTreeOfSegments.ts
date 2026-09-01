// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BVTreeOfSegments.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in BVTree.ts regarding tree construction. Although this
// class appears to be non-abstract upstream, the BoundingVolume type has
// requirements for its interface. In this sense, BVTreeOfSegments is
// abstract.
//
// Port notes:
//   - See BVTreeOfPoints.ts for the BoundingVolume type-parameter split and
//     for why the class is declared 'abstract'.
//   - Upstream 'Create(vertices, segments, height)' hides the base class
//     'Create(centroids, height)' by C++ name hiding. TypeScript requires an
//     override to be assignable to the base member, so the segment-specific
//     creation function is named 'createFromSegments'. The inherited
//     'create(centroids, height)' remains callable but is not meaningful for
//     a segment tree.
//   - 'std::array<std::size_t, 2>' becomes the tuple type [number, number].
//   - 'Execute' becomes 'execute' and returns the node indices rather than
//     filling an output std::vector.

import { BVTree } from './BVTree';
import type { BVTreeBoundingVolume, BVTreeVolumeOps } from './BVTree';
import { logAssert } from './Logger';
import { Vector, add, mul } from './Vector';

export abstract class BVTreeOfSegments<BV extends BVTreeBoundingVolume> extends BVTree<BV> {
    protected mVertices: Vector[];
    protected mSegments: [number, number][];

    protected constructor(ops: BVTreeVolumeOps<BV>) {
        super(ops);
        this.mVertices = [];
        this.mSegments = [];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If BVTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from segments.length. If
    // larger than 31, the height is clamped to 31.
    createFromSegments(vertices: readonly Vector[],
        segments: readonly (readonly [number, number])[],
        height: number = BVTree.fullHeight): void {
        logAssert(vertices.length > 0,
            'Expecting vertices to create a bounding volume tree.');

        this.mVertices = new Array<Vector>(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            this.mVertices[i] = vertices[i].clone();
        }
        this.mSegments = new Array<[number, number]>(segments.length);
        for (let i = 0; i < segments.length; ++i) {
            this.mSegments[i] = [segments[i][0], segments[i][1]];
        }

        // Compute the segment centroids.
        const centroids = new Array<Vector>(this.mSegments.length);
        const half = 0.5;
        for (let i = 0; i < this.mSegments.length; ++i) {
            const seg = this.mSegments[i];
            centroids[i] = mul(half,
                add(this.mVertices[seg[0]], this.mVertices[seg[1]]));
        }

        // Create the bounding volume tree for centroids.
        super.create(centroids, height);
    }

    // Member access.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getSegments(): readonly (readonly [number, number])[] {
        return this.mSegments;
    }

    // Compute intersections of the linear component and leaf nodes. The
    // returned indices are lookups into the mNodes[] member of the base
    // class. They are ordered according to the depth-first traversal of the
    // tree.
    execute(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}
