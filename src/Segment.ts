// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Segment.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The segment is represented by (1-t)*P0 + t*P1, where P0 and P1 are the
// endpoints of the segment and 0 <= t <= 1. Some algorithms prefer a centered
// representation that is similar to how oriented bounding boxes are defined.
// This representation is C + s*D, where C = (P0 + P1)/2 is the center of the
// segment, D = (P1 - P0)/|P1 - P0| is a unit-length direction vector for the
// segment, and |s| <= e. The value e = |P1 - P0|/2 is the extent (or radius
// or half-length) of the segment.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Segment(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods). The C++ 'std::array<Vector<N, Real>, 2> p' becomes a
// Vector[] of length 2; its comparisons are lexicographic over the elements,
// as std::array's are. 'GetCenteredForm' returns { center, direction, extent }
// instead of writing to output references.

import { logAssert } from './Logger';
import { Vector, add, sub, mul, negate, normalize } from './Vector';

// Lexicographic comparison of the two-element endpoint arrays (the port of
// std::array's relational operators). Returns -1, 0 or +1.
function comparePoints(p0: readonly Vector[], p1: readonly Vector[]): number {
    for (let i = 0; i < 2; ++i) {
        if (p0[i].lessThan(p1[i])) {
            return -1;
        }
        if (p1[i].lessThan(p0[i])) {
            return +1;
        }
    }
    return 0;
}

export class Segment {
    // Public member access: the two endpoints.
    p: Vector[];

    // The port of the default constructor, which sets p0 to (-1,0,...,0) and
    // p1 to (1,0,...,0). The dimension N of the C++ template is a constructor
    // argument here. NOTE: If you set p0 and p1; compute C, D, and e; and then
    // recompute q0 = C-e*D and q1 = C+e*D, numerical round-off errors can lead
    // to q0 not exactly equal to p0 and q1 not exactly equal to p1.
    constructor(n: number) {
        const p1 = Vector.unit(n, 0);
        this.p = [negate(p1), p1];
    }

    // The port of 'Segment(p0, p1)'. The vectors are copied, matching C++
    // value semantics.
    static fromEndpoints(p0: Vector, p1: Vector): Segment {
        logAssert(p0.size === p1.size, 'Segment: mismatched sizes.');
        const segment = new Segment(p0.size);
        segment.p = [p0.clone(), p1.clone()];
        return segment;
    }

    // The port of 'Segment(std::array<Vector<N, Real>, 2> const& inP)'.
    static fromPointArray(inP: readonly Vector[]): Segment {
        logAssert(inP.length === 2, 'Segment: invalid number of endpoints.');
        return Segment.fromEndpoints(inP[0], inP[1]);
    }

    // The port of 'Segment(center, direction, extent)'.
    static fromCenteredForm(center: Vector, direction: Vector,
        extent: number): Segment {
        logAssert(center.size === direction.size,
            'Segment: mismatched sizes.');
        const segment = new Segment(center.size);
        segment.setCenteredForm(center, direction, extent);
        return segment;
    }

    // The dimension N of the segment.
    get dimension(): number {
        return this.p[0].size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Segment {
        return Segment.fromEndpoints(this.p[0], this.p[1]);
    }

    // Manipulation via the centered form.
    setCenteredForm(center: Vector, direction: Vector, extent: number): void {
        logAssert(center.size === direction.size,
            'Segment: mismatched sizes.');
        this.p = [
            sub(center, mul(extent, direction)),
            add(center, mul(extent, direction))
        ];
    }

    // The port of 'GetCenteredForm(center, direction, extent)'; the C++
    // output references become the fields of the returned object.
    getCenteredForm(): { center: Vector, direction: Vector, extent: number } {
        const center = mul(0.5, add(this.p[0], this.p[1]));
        const direction = sub(this.p[1], this.p[0]);
        const extent = 0.5 * normalize(direction);
        return { center, direction, extent };
    }

    // Comparisons to support sorted containers.
    equals(segment: Segment): boolean {
        return comparePoints(this.p, segment.p) === 0;
    }

    notEquals(segment: Segment): boolean {
        return !this.equals(segment);
    }

    lessThan(segment: Segment): boolean {
        return comparePoints(this.p, segment.p) < 0;
    }

    lessThanOrEqual(segment: Segment): boolean {
        return comparePoints(this.p, segment.p) <= 0;
    }

    greaterThan(segment: Segment): boolean {
        return comparePoints(this.p, segment.p) > 0;
    }

    greaterThanOrEqual(segment: Segment): boolean {
        return comparePoints(this.p, segment.p) >= 0;
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Segment2 = Segment;
export type Segment3 = Segment;
