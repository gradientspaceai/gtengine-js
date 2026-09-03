// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Capsule.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A capsule is the set of points that are equidistant from a segment, the
// common distance called the radius.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Capsule(n)' for the default
// constructor, named static factories that copy their arguments, comparison
// methods).

import { Segment } from './Segment.js';

export class Capsule {
    // Public member access.
    segment: Segment;
    radius: number;

    // The port of the default constructor, which sets the segment to have
    // endpoints p0 = (-1,0,...,0) and p1 = (1,0,...,0), and the radius to 1.
    // The dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        this.segment = new Segment(n);
        this.radius = 1;
    }

    // The port of 'Capsule(inSegment, inRadius)'. The segment is copied,
    // matching C++ value semantics.
    static fromSegmentRadius(inSegment: Segment, inRadius: number): Capsule {
        const capsule = new Capsule(inSegment.dimension);
        capsule.segment = inSegment.clone();
        capsule.radius = inRadius;
        return capsule;
    }

    // The dimension N of the capsule.
    get dimension(): number {
        return this.segment.dimension;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Capsule {
        return Capsule.fromSegmentRadius(this.segment, this.radius);
    }

    // Comparisons to support sorted containers.
    equals(capsule: Capsule): boolean {
        return this.segment.equals(capsule.segment)
            && this.radius === capsule.radius;
    }

    notEquals(capsule: Capsule): boolean {
        return !this.equals(capsule);
    }

    lessThan(capsule: Capsule): boolean {
        if (this.segment.lessThan(capsule.segment)) {
            return true;
        }

        if (this.segment.greaterThan(capsule.segment)) {
            return false;
        }

        return this.radius < capsule.radius;
    }

    lessThanOrEqual(capsule: Capsule): boolean {
        return !capsule.lessThan(this);
    }

    greaterThan(capsule: Capsule): boolean {
        return capsule.lessThan(this);
    }

    greaterThanOrEqual(capsule: Capsule): boolean {
        return !this.lessThan(capsule);
    }
}

// Alias for convenience (the port of the upstream template alias).
export type Capsule3 = Capsule;
