// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistPointCanonicalBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance from a point to a solid canonical box in nD.
//
// The canonical box has center at the origin and is aligned with the
// coordinate axes. The extents are E = (e[0],e[1],...,e[n-1]). A box point is
// Y = (y[0],y[1],...,y[n-1]) with |y[i]| <= e[i] for all i.
//
// The input point P is stored in closest[0]. The closest point on the box is
// stored in closest[1]. When there are infinitely many choices for the pair
// of closest points, only one of them is returned.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization 'DCPQuery<T, Vector<N,T>, CanonicalBox<N,T>>'
// becomes the class DistPointCanonicalBox with the result type
// DistPointCanonicalBoxResult.

import type { CanonicalBox } from './CanonicalBox';
import type { DCPQuery } from './DCPQuery';
import { Vector } from './Vector';

export interface DistPointCanonicalBoxResult {
    distance: number;
    sqrDistance: number;

    // closest[0] is the input point, closest[1] is the closest box point.
    closest: [Vector, Vector];
}

export class DistPointCanonicalBox
    implements DCPQuery<Vector, CanonicalBox, DistPointCanonicalBoxResult> {
    compute(point: Vector, box: CanonicalBox): DistPointCanonicalBoxResult {
        const closest0 = point.clone();
        const closest1 = point.clone();
        let sqrDistance = 0;
        const n = box.extent.size;
        for (let i = 0; i < n; ++i) {
            const e = box.extent.values[i];
            if (point.values[i] < -e) {
                const delta = closest1.values[i] + e;
                sqrDistance += delta * delta;
                closest1.values[i] = -e;
            }
            else if (point.values[i] > e) {
                const delta = closest1.values[i] - e;
                sqrDistance += delta * delta;
                closest1.values[i] = e;
            }
        }

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest: [closest0, closest1]
        };
    }
}
