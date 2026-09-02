// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OrientedBoxBV.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Class OrientedBoxBV is a bounding volume that supports the queries based on
// BVTree and its derived classes.
//
// Port notes:
//   - The upstream 'BoundingVolume' template parameter of BVTree carries both
//     instance requirements (GetSplittingAxis) and static requirements
//     (default construction, IntersectLine/IntersectRay/IntersectSegment).
//     The port of BVTree splits these into the 'BVTreeBoundingVolume'
//     interface and a 'BVTreeVolumeOps' object passed to the base
//     constructor; see BVTree.ts. OrientedBoxBV implements the former and
//     keeps the three upstream static functions as static methods;
//     'orientedBoxBVOps' bundles them with default construction.
//   - 'GetSplittingAxis(origin, direction)' has output reference parameters
//     upstream; the port returns a 'BVTreeSplittingAxis' object.
//   - The static intersection predicates take the same (P, Q) pairs as
//     upstream: for a line and a ray, P is the origin and Q is the
//     unit-length direction; for a segment, P and Q are the endpoints.
//   - The upstream default constructor is 'box{}', which value-initializes an
//     OrientedBox3<T> and so runs its default constructor (center (0,0,0),
//     axis d = Unit(d), extent (1,1,1)); the port is 'new OrientedBox(3)'.

import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from './BVTree';
import { IntrLine3OrientedBox3TI } from './IntrLine3OrientedBox3';
import { IntrRay3OrientedBox3TI } from './IntrRay3OrientedBox3';
import { IntrSegment3OrientedBox3TI } from './IntrSegment3OrientedBox3';
import { Line } from './Line';
import { OrientedBox } from './OrientedBox';
import { Ray } from './Ray';
import { Segment } from './Segment';
import { Vector } from './Vector';

export class OrientedBoxBV implements BVTreeBoundingVolume {
    // Public member access.
    box: OrientedBox;

    constructor() {
        this.box = new OrientedBox(3);
    }

    // The port of 'OrientedBoxBV' constructed around an existing box. The box
    // is copied, matching C++ value semantics.
    static fromBox(inBox: OrientedBox): OrientedBoxBV {
        const bv = new OrientedBoxBV();
        bv.box = inBox.clone();
        return bv;
    }

    // The splitting axis passes through the box center in the direction of
    // the box axis of largest extent. Ties are broken in favor of the smaller
    // axis index, because the comparisons are strict inequalities.
    getSplittingAxis(): BVTreeSplittingAxis {
        const origin = this.box.center.clone();

        let maxExtent = this.box.extent.get(0);
        let maxIndex = 0;
        if (this.box.extent.get(1) > maxExtent) {
            maxExtent = this.box.extent.get(1);
            maxIndex = 1;
        }
        if (this.box.extent.get(2) > maxExtent) {
            maxIndex = 2;
        }
        const direction = this.box.axis[maxIndex].clone();
        return { origin: origin, direction: direction };
    }

    // The line is P + t * Q for all real t.
    static intersectLine(P: Vector, Q: Vector,
        boundingVolume: OrientedBoxBV): boolean {
        const query = new IntrLine3OrientedBox3TI();
        const output = query.test(Line.fromOriginDirection(P, Q),
            boundingVolume.box);
        return output.intersect;
    }

    // The ray is P + t * Q for t >= 0.
    static intersectRay(P: Vector, Q: Vector,
        boundingVolume: OrientedBoxBV): boolean {
        const query = new IntrRay3OrientedBox3TI();
        const output = query.test(Ray.fromOriginDirection(P, Q),
            boundingVolume.box);
        return output.intersect;
    }

    // The segment has endpoints P and Q.
    static intersectSegment(P: Vector, Q: Vector,
        boundingVolume: OrientedBoxBV): boolean {
        const query = new IntrSegment3OrientedBox3TI();
        const output = query.test(Segment.fromEndpoints(P, Q),
            boundingVolume.box);
        return output.intersect;
    }
}

// The static-side requirements of the upstream BoundingVolume template
// parameter, bundled for the BVTree base class; see BVTree.ts.
export const orientedBoxBVOps: BVTreeVolumeOps<OrientedBoxBV> = {
    create: () => new OrientedBoxBV(),
    intersectLine: (P, Q, bv) => OrientedBoxBV.intersectLine(P, Q, bv),
    intersectRay: (P, Q, bv) => OrientedBoxBV.intersectRay(P, Q, bv),
    intersectSegment: (P, Q, bv) => OrientedBoxBV.intersectSegment(P, Q, bv)
};
