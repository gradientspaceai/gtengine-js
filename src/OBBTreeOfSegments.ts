// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OBBTreeOfSegments.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in OBBTree.ts regarding tree construction.
//
// Port notes:
//   - Upstream 'Create(vertices, segments, height)' hides the base class
//     'Create(centroids, height)' by C++ name hiding. TypeScript requires an
//     override to be assignable to the base member, so the segment-specific
//     creation function is named 'createFromSegments'. The inherited
//     'create(centroids, height)' remains callable but is not meaningful for
//     a segment tree.
//   - 'std::array<size_t, 2>' becomes the tuple type [number, number].
//   - The upstream overrides ComputeInteriorBox and ComputeLeafBox are
//     'private virtual'; they are 'protected override' here because
//     TypeScript requires an override to be at least as visible as the base
//     member.
//   - Upstream initializes pmin and pmax to the zero vector rather than to
//     the projections of the first primitive vertex. Combined with the
//     'if (dot < pmin[j]) ... else if (dot > pmax[j]) ...' update, every
//     projection is still covered, so the box does contain all the segment
//     endpoints; it is merely conservative. The quirk is preserved; see
//     upstream issue #103.

import { logAssert } from './Logger.js';
import { OBBTree } from './OBBTree.js';
import { OrientedBox } from './OrientedBox.js';
import { Vector, add, dot, mul, normalize, sub } from './Vector.js';
import { computeOrthogonalComplement3 } from './Vector3.js';

export class OBBTreeOfSegments extends OBBTree {
    protected mVertices: Vector[];
    protected mSegments: [number, number][];

    constructor() {
        super();
        this.mVertices = [];
        this.mSegments = [];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If OBBTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from segments.length. If
    // larger than 31, the height is clamped to 31.
    createFromSegments(vertices: readonly Vector[],
        segments: readonly (readonly [number, number])[],
        height: number = OBBTree.fullHeight): void {
        logAssert(vertices.length >= 2 && segments.length > 0, 'Invalid input.');

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

        // Create the OBB tree for centroids.
        super.create(centroids, height);
    }

    // Member access.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getSegments(): readonly (readonly [number, number])[] {
        return this.mSegments;
    }

    // Let C be the box center and let U0, U1 and U2 be the box axes. Each
    // input point is of the form X = C + y0*U0 + y1*U1 + y2*U2. The following
    // code computes min(y0), max(y0), min(y1), max(y1), min(y2) and max(y2).
    // The box center is then adjusted to be
    //   C' = C + 0.5*(min(y0)+max(y0))*U0 + 0.5*(min(y1)+max(y1))*U1
    //        + 0.5*(min(y2)+max(y2))*U2
    protected override computeInteriorBox(i0: number, i1: number, box: OrientedBox): void {
        super.computeInteriorBox(i0, i1, box);

        const pmin = Vector.zero(3);
        const pmax = Vector.zero(3);
        for (let i = i0; i <= i1; ++i) {
            const seg = this.mSegments[this.mPartition[i]];
            for (let k = 0; k < 2; ++k) {
                const diff = sub(this.mVertices[seg[k]], box.center);
                for (let j = 0; j < 3; ++j) {
                    const d = dot(diff, box.axis[j]);
                    if (d < pmin.values[j]) {
                        pmin.values[j] = d;
                    } else if (d > pmax.values[j]) {
                        pmax.values[j] = d;
                    }
                }
            }
        }

        const half = 0.5;
        for (let j = 0; j < 3; ++j) {
            const scalar = half * (pmin.values[j] + pmax.values[j]);
            const axis = box.axis[j];
            for (let k = 0; k < 3; ++k) {
                box.center.values[k] += scalar * axis.values[k];
            }
            box.extent.values[j] = half * (pmax.values[j] - pmin.values[j]);
        }
    }

    protected override computeLeafBox(i: number, box: OrientedBox): void {
        // Create a degenerate box whose center is the midpoint of the segment
        // primitive, whose axis[0] is the segment direction and whose
        // extent[0] is half the length of the segment.
        const seg = this.mSegments[this.mPartition[i]];
        box.center = this.mCentroids[this.mPartition[i]].clone();
        box.axis[0] = sub(this.mVertices[seg[1]], this.mVertices[seg[0]]);
        box.extent.values[0] = 0.5 * normalize(box.axis[0]);
        computeOrthogonalComplement3(1, box.axis);
        box.extent.values[1] = 0;
        box.extent.values[2] = 0;
    }
}
