// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Parallelogram2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A parallelogram has center C and axes A0 and A1. Its points are
// C + s0*A0 + s1*A1 with |s0| <= 1 and |s1| <= 1. The axes form a
// right-handed basis but need be neither orthogonal nor unit length.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 2D vectors. The
// C++ 'std::array<Vector2<T>, 2> axis' becomes a Vector[] of length 2 and
// its comparisons are lexicographic over the elements, as std::array's are.
// 'GetVertices' returns a new array of 4 vertices instead of filling a
// caller-supplied std::array.

import { logAssert } from './Logger.js';
import { Vector, add, sub } from './Vector.js';
import { dotPerp } from './Vector2.js';

// Lexicographic comparison of the two-element axis arrays (the port of
// std::array's relational operators). Returns -1, 0 or +1.
function compareAxes(a0: readonly Vector[], a1: readonly Vector[]): number {
    for (let i = 0; i < 2; ++i) {
        if (a0[i].lessThan(a1[i])) {
            return -1;
        }
        if (a1[i].lessThan(a0[i])) {
            return +1;
        }
    }
    return 0;
}

export class Parallelogram2 {
    // Public member access.
    center: Vector;
    axis: Vector[];

    // The port of the default constructor, which sets the center to (0,0),
    // axis[0] to (1,0), and axis[1] to (0,1).
    constructor() {
        this.center = new Vector(2);
        this.axis = [Vector.unit(2, 0), Vector.unit(2, 1)];
    }

    // The port of 'Parallelogram2(inCenter, inAxis)'. The axes must form a
    // right-handed basis. The axes do not have to be orthogonal. The axis
    // lengths do not have to be unit length. The vectors are copied, matching
    // C++ value semantics.
    static fromCenterAxis(inCenter: Vector,
        inAxis: readonly Vector[]): Parallelogram2 {
        logAssert(inCenter.size === 2 && inAxis.length === 2
            && inAxis[0].size === 2 && inAxis[1].size === 2,
            'Parallelogram2: mismatched sizes.');
        logAssert(dotPerp(inAxis[0], inAxis[1]) > 0,
            'The axes must form a right-handed basis.');
        const parallelogram = new Parallelogram2();
        parallelogram.center = inCenter.clone();
        parallelogram.axis = [inAxis[0].clone(), inAxis[1].clone()];
        return parallelogram;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Parallelogram2 {
        const parallelogram = new Parallelogram2();
        parallelogram.center = this.center.clone();
        parallelogram.axis = [this.axis[0].clone(), this.axis[1].clone()];
        return parallelogram;
    }

    // If index i has the bit pattern i = b[1]b[0], then
    //   vertices[i] = center + sum_{d=0}^{1} sign[d] * axis[d]
    // where sign[d] = 2*b[d] - 1. (The upstream comment claims the vertices
    // are listed in counterclockwise order, but the bit-pattern order visits
    // vertices 2 and 3 in the opposite of that order. The port preserves the
    // upstream ordering.)
    getVertices(): Vector[] {
        const a0 = this.axis[0];
        const a1 = this.axis[1];
        const c = this.center;
        return [
            sub(sub(c, a0), a1),
            sub(add(c, a0), a1),
            add(sub(c, a0), a1),
            add(add(c, a0), a1)
        ];
    }

    // Comparisons to support sorted containers.
    equals(other: Parallelogram2): boolean {
        return this.center.equals(other.center)
            && compareAxes(this.axis, other.axis) === 0;
    }

    notEquals(other: Parallelogram2): boolean {
        return !this.equals(other);
    }

    lessThan(other: Parallelogram2): boolean {
        if (this.center.lessThan(other.center)) {
            return true;
        }

        if (this.center.greaterThan(other.center)) {
            return false;
        }

        return compareAxes(this.axis, other.axis) < 0;
    }

    lessThanOrEqual(other: Parallelogram2): boolean {
        return !other.lessThan(this);
    }

    greaterThan(other: Parallelogram2): boolean {
        return other.lessThan(this);
    }

    greaterThanOrEqual(other: Parallelogram2): boolean {
        return !this.lessThan(other);
    }
}
