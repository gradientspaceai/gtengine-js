// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BVTreeOfPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in BVTree.ts regarding tree construction. Although this
// class appears to be non-abstract upstream, the BoundingVolume type has
// requirements for its interface. In this sense, BVTreeOfPoints is abstract.
//
// Port notes:
//   - The upstream template parameter BoundingVolume becomes the type
//     parameter BV (constrained by BVTreeBoundingVolume) plus the
//     BVTreeVolumeOps<BV> object that BVTree requires; see BVTree.ts. The
//     class is declared 'abstract' in TypeScript because BVTree's
//     computeInteriorBoundingVolume and computeLeafBoundingVolume are
//     abstract and are supplied by the concrete bounding-volume-specific
//     derived classes, exactly as upstream.
//   - 'Execute' becomes 'execute' and returns the node indices rather than
//     filling an output std::vector.

import { BVTree } from './BVTree';
import type { BVTreeBoundingVolume, BVTreeVolumeOps } from './BVTree';
import { logAssert } from './Logger';
import { Vector } from './Vector';

export abstract class BVTreeOfPoints<BV extends BVTreeBoundingVolume> extends BVTree<BV> {
    protected mVertices: Vector[];

    protected constructor(ops: BVTreeVolumeOps<BV>) {
        super(ops);
        this.mVertices = [];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If BVTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from vertices.length. If
    // larger than 31, the height is clamped to 31.
    override create(vertices: readonly Vector[], height: number = BVTree.fullHeight): void {
        logAssert(vertices.length > 0,
            'Expecting vertices to create a bounding volume tree.');

        this.mVertices = new Array<Vector>(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            this.mVertices[i] = vertices[i].clone();
        }

        // The vertices are already the centroids. BVTree.create copies the
        // centroids it is given, so the upstream 'centroids' copy is implicit
        // here.
        super.create(this.mVertices, height);
    }

    // Member access.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    // Compute intersections of the linear component and leaf nodes. The
    // returned indices are lookups into the mNodes[] member of the base
    // class. They are ordered according to the depth-first traversal of the
    // tree.
    execute(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}
