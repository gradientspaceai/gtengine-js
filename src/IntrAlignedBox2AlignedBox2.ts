// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrAlignedBox2AlignedBox2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The queries consider the box to be a solid.
//
// The aligned-aligned queries use simple min-max comparisons. The
// intersection of aligned boxes is an aligned box, possibly degenerate,
// where min[d] == max[d] for at least one dimension d.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent (one class per
// upstream TIQuery/FIQuery specialization, suffixed TI/FI, with test()/find()
// and exported <ClassName>Result types).

import { AlignedBox } from './AlignedBox.js';
import type { TIQuery } from './TIQuery.js';
import type { FIQuery } from './FIQuery.js';

// The result of IntrAlignedBox2AlignedBox2TI.test.
export interface IntrAlignedBox2AlignedBox2TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrAlignedBox2AlignedBox2TIResult {
    return { intersect: false };
}

// The result of IntrAlignedBox2AlignedBox2FI.find.
export interface IntrAlignedBox2AlignedBox2FIResult {
    intersect: boolean;

    // Valid only when 'intersect' is true. The intersection of two aligned
    // boxes is an aligned box, possibly degenerate.
    box: AlignedBox;
}

// The port of the upstream FIQuery::Result default constructor. Upstream
// value-initializes 'box', which calls the AlignedBox default constructor
// (min = -1 and max = +1 in each dimension).
function defaultFIResult(): IntrAlignedBox2AlignedBox2FIResult {
    return { intersect: false, box: new AlignedBox(2) };
}

export class IntrAlignedBox2AlignedBox2TI implements
    TIQuery<AlignedBox, AlignedBox, IntrAlignedBox2AlignedBox2TIResult> {

    test(box0: AlignedBox, box1: AlignedBox): IntrAlignedBox2AlignedBox2TIResult {
        const result = defaultTIResult();
        for (let i = 0; i < 2; ++i) {
            if (box0.max.values[i] < box1.min.values[i] ||
                box0.min.values[i] > box1.max.values[i]) {
                result.intersect = false;
                return result;
            }
        }
        result.intersect = true;
        return result;
    }
}

export class IntrAlignedBox2AlignedBox2FI implements
    FIQuery<AlignedBox, AlignedBox, IntrAlignedBox2AlignedBox2FIResult> {

    find(box0: AlignedBox, box1: AlignedBox): IntrAlignedBox2AlignedBox2FIResult {
        const result = defaultFIResult();
        for (let i = 0; i < 2; ++i) {
            if (box0.max.values[i] < box1.min.values[i] ||
                box0.min.values[i] > box1.max.values[i]) {
                result.intersect = false;
                return result;
            }
        }

        for (let i = 0; i < 2; ++i) {
            if (box0.max.values[i] <= box1.max.values[i]) {
                result.box.max.values[i] = box0.max.values[i];
            } else {
                result.box.max.values[i] = box1.max.values[i];
            }

            if (box0.min.values[i] <= box1.min.values[i]) {
                result.box.min.values[i] = box1.min.values[i];
            } else {
                result.box.min.values[i] = box0.min.values[i];
            }
        }
        result.intersect = true;
        return result;
    }
}
