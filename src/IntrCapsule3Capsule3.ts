// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntrCapsule3Capsule3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test-intersection query for two solid capsules. Two capsules intersect when
// the distance between their medial segments is no larger than the sum of the
// capsule radii.
//
// Port notes: see IntrIntervals.ts for the Intr* precedent. Upstream has only
// a TIQuery specialization, so the port has only IntrCapsule3Capsule3TI.

import type { Capsule } from './Capsule.js';
import { DistSegmentSegment } from './DistSegmentSegment.js';
import { logAssert } from './Logger.js';
import type { TIQuery } from './TIQuery.js';

// The result of IntrCapsule3Capsule3TI.test.
export interface IntrCapsule3Capsule3TIResult {
    intersect: boolean;
}

// The port of the upstream TIQuery::Result default constructor.
function defaultTIResult(): IntrCapsule3Capsule3TIResult {
    return { intersect: false };
}

export class IntrCapsule3Capsule3TI implements
    TIQuery<Capsule, Capsule, IntrCapsule3Capsule3TIResult> {

    test(capsule0: Capsule, capsule1: Capsule): IntrCapsule3Capsule3TIResult {
        logAssert(capsule0.dimension === 3 && capsule1.dimension === 3,
            'IntrCapsule3Capsule3TI: mismatched sizes.');

        const result = defaultTIResult();
        const ssQuery = new DistSegmentSegment();
        const ssResult = ssQuery.compute(capsule0.segment, capsule1.segment);
        const rSum = capsule0.radius + capsule1.radius;
        result.intersect = (ssResult.distance <= rSum);
        return result;
    }
}
