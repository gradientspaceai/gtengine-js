// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OBBTreeOfPoints.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in OBBTree.ts regarding tree construction.
//
// Port notes:
//   - The upstream IndexType template discussion does not apply; indices are
//     'number' here, as everywhere in the port.
//   - The upstream overrides ComputeInteriorBox and ComputeLeafBox are
//     'private virtual'; they are 'protected override' here because
//     TypeScript requires an override to be at least as visible as the base
//     member.
//   - Upstream initializes pmin and pmax to the zero vector rather than to
//     the projections of the first primitive vertex. Combined with the
//     'if (dot < pmin[j]) ... else if (dot > pmax[j]) ...' update, every
//     projection is still covered, so the box does contain all the points;
//     it is merely conservative (it is forced to contain the mean of the
//     centroids). The quirk is preserved; see upstream issue #103.

import { OBBTree } from './OBBTree';
import { OrientedBox } from './OrientedBox';
import { Vector, dot, sub } from './Vector';

export class OBBTreeOfPoints extends OBBTree {
    constructor() {
        super();
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If OBBTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from points.length. If
    // larger than 31, the height is clamped to 31.
    override create(points: readonly Vector[], height: number = OBBTree.fullHeight): void {
        // Create the OBB tree for centroids. The points are already the
        // centroids.
        super.create(points, height);
    }

    // Member access.
    getPoints(): readonly Vector[] {
        return this.mCentroids;
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
            const diff = sub(this.mCentroids[this.mPartition[i]], box.center);
            for (let j = 0; j < 3; ++j) {
                const d = dot(diff, box.axis[j]);
                if (d < pmin.values[j]) {
                    pmin.values[j] = d;
                } else if (d > pmax.values[j]) {
                    pmax.values[j] = d;
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
        // Create a degenerate box whose center is the point primitive.
        box.center = this.mCentroids[this.mPartition[i]].clone();
        box.axis[0] = Vector.unit(3, 0);
        box.axis[1] = Vector.unit(3, 1);
        box.axis[2] = Vector.unit(3, 2);
        box.extent = Vector.zero(3);
    }
}
