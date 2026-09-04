// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) DistAlignedBoxAlignedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the distance between two solid aligned boxes in nD.
//
// Each aligned box has minimum corner A and maximum corner B. A box point is
// X where A <= X <= B; the comparisons are componentwise.
//
// The algorithm computes two aligned boxes of closest points, closest[0] for
// input box0 and closest[1] for input box1. One reasonable choice of closest
// pair is
//   P0 = (closest[0].min + closest[0].max)/2;
//   P1 = (closest[1].min + closest[1].max)/2;
//
// The upstream comment states more generally that "any choice of P0 in
// closest[0] and any choice of P1 in closest[1]" is a closest pair. That is
// an overstatement: on an axis where the two input boxes overlap, both
// closest boxes carry the same nondegenerate interval, and only pairs whose
// coordinates on that axis agree realize the distance. The precise statement
// is that P0 and P1 must use the same interpolation parameter per axis; the
// midpoint choice above is that statement with parameter 1/2 and is the one
// upstream recommends.
//
// Port notes: see DistPointLine.ts for the Dist* family conventions. The
// upstream specialization
// 'DCPQuery<T, AlignedBox<N,T>, AlignedBox<N,T>>' becomes the class
// DistAlignedBoxAlignedBox with the result type
// DistAlignedBoxAlignedBoxResult. Upstream default-constructs the two
// AlignedBox members of Result (which for N-dimensional boxes leaves the
// [-1,1]^N default box) and then overwrites every component; the port
// allocates zero-filled boxes of the correct runtime dimension for the same
// reason.

import { AlignedBox } from './AlignedBox.js';
import type { DCPQuery } from './DCPQuery.js';
import { IntrIntervalsFI } from './IntrIntervals.js';
import { Vector } from './Vector.js';

export interface DistAlignedBoxAlignedBoxResult {
    distance: number;
    sqrDistance: number;

    // closest[i] is the box of points of box[i] that realize the distance.
    closest: [AlignedBox, AlignedBox];
}

export class DistAlignedBoxAlignedBox
    implements DCPQuery<AlignedBox, AlignedBox, DistAlignedBoxAlignedBoxResult> {
    compute(box0: AlignedBox, box1: AlignedBox): DistAlignedBoxAlignedBoxResult {
        const n = box0.min.size;
        const closest: [AlignedBox, AlignedBox] = [
            AlignedBox.fromMinMax(new Vector(n), new Vector(n)),
            AlignedBox.fromMinMax(new Vector(n), new Vector(n))
        ];

        let sqrDistance = 0;
        const iiQuery = new IntrIntervalsFI();
        for (let i = 0; i < n; ++i) {
            const min0 = box0.min.values[i];
            const max0 = box0.max.values[i];
            const min1 = box1.min.values[i];
            const max1 = box1.max.values[i];

            if (min0 >= max1) {
                const delta = min0 - max1;
                sqrDistance += delta * delta;
                closest[0].min.values[i] = min0;
                closest[0].max.values[i] = min0;
                closest[1].min.values[i] = max1;
                closest[1].max.values[i] = max1;
            }
            else if (min1 >= max0) {
                const delta = min1 - max0;
                sqrDistance += delta * delta;
                closest[0].min.values[i] = max0;
                closest[0].max.values[i] = max0;
                closest[1].min.values[i] = min1;
                closest[1].max.values[i] = min1;
            }
            else {  // min0 <= max1 and min1 <= max0
                const iiResult = iiQuery.find([min0, max0], [min1, max1]);
                for (let j = 0; j < 2; ++j) {
                    closest[j].min.values[i] = iiResult.overlap[0];
                    closest[j].max.values[i] = iiResult.overlap[1];
                }
            }
        }

        return {
            distance: Math.sqrt(sqrDistance),
            sqrDistance,
            closest
        };
    }
}
