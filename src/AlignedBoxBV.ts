// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AlignedBoxBV.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Class AlignedBoxBV is a bounding volume that supports the queries based on
// BVTree and its derived classes.
//
// Port notes: the upstream GetSplittingAxis has two output reference
// parameters; the port returns an object literal { origin, direction }. The
// static intersection predicates keep their upstream names in camelCase and
// take the same (P, Q) pairs: for a line and a ray, P is the origin and Q is
// the unit-length direction; for a segment, P and Q are the endpoints.

import { AlignedBox } from './AlignedBox';
import { IntrLine3AlignedBox3TI } from './IntrLine3AlignedBox3';
import { IntrRay3AlignedBox3TI } from './IntrRay3AlignedBox3';
import { IntrSegment3AlignedBox3TI } from './IntrSegment3AlignedBox3';
import { Line } from './Line';
import { Ray } from './Ray';
import { Segment } from './Segment';
import { Vector, add, mul, sub } from './Vector';

export class AlignedBoxBV {
    // Public member access.
    box: AlignedBox;

    // The port of the default constructor, which value-initializes the box.
    // The port's AlignedBox default constructor sets the minimum values to
    // -1 and the maximum values to +1, so the box is set explicitly to the
    // zero-initialized C++ state.
    constructor() {
        this.box = AlignedBox.fromMinMax(new Vector(3), new Vector(3));
    }

    // The port of 'AlignedBoxBV' constructed around an existing box. The box
    // is copied, matching C++ value semantics.
    static fromBox(inBox: AlignedBox): AlignedBoxBV {
        const bv = new AlignedBoxBV();
        bv.box = AlignedBox.fromMinMax(inBox.min, inBox.max);
        return bv;
    }

    getSplittingAxis(): { origin: Vector, direction: Vector } {
        const half = 0.5;

        const origin = mul(add(this.box.max, this.box.min), half);

        const extents = mul(sub(this.box.max, this.box.min), half);
        let maxExtent = extents.values[0];
        let maxIndex = 0;
        if (extents.values[1] > maxExtent) {
            maxExtent = extents.values[1];
            maxIndex = 1;
        }
        if (extents.values[2] > maxExtent) {
            maxIndex = 2;
        }
        const direction = Vector.unit(3, maxIndex);
        return { origin, direction };
    }

    static intersectLine(P: Vector, Q: Vector,
        boundingVolume: AlignedBoxBV): boolean {
        const query = new IntrLine3AlignedBox3TI();
        const output = query.test(Line.fromOriginDirection(P, Q),
            boundingVolume.box);
        return output.intersect;
    }

    static intersectRay(P: Vector, Q: Vector,
        boundingVolume: AlignedBoxBV): boolean {
        const query = new IntrRay3AlignedBox3TI();
        const output = query.test(Ray.fromOriginDirection(P, Q),
            boundingVolume.box);
        return output.intersect;
    }

    static intersectSegment(P: Vector, Q: Vector,
        boundingVolume: AlignedBoxBV): boolean {
        const query = new IntrSegment3AlignedBox3TI();
        const output = query.test(Segment.fromEndpoints(P, Q),
            boundingVolume.box);
        return output.intersect;
    }
}
