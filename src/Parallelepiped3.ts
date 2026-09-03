// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Parallelepiped3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A parallelepiped has center C and axes A0, A1 and A2. Its points are
// C + s0*A0 + s1*A1 + s2*A2 with |si| <= 1. The axes form a right-handed
// basis but need be neither orthogonal nor unit length.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 3D vectors. The
// C++ 'std::array<Vector3<T>, 3> axis' becomes a Vector[] of length 3 and
// its comparisons are lexicographic over the elements, as std::array's are.
// 'GetVertices' returns a new array of 8 vertices instead of filling a
// caller-supplied std::array.

import { logAssert } from './Logger.js';
import { Vector, add, sub } from './Vector.js';
import { dotCross } from './Vector3.js';

// Lexicographic comparison of the three-element axis arrays (the port of
// std::array's relational operators). Returns -1, 0 or +1.
function compareAxes(a0: readonly Vector[], a1: readonly Vector[]): number {
    for (let i = 0; i < 3; ++i) {
        if (a0[i].lessThan(a1[i])) {
            return -1;
        }
        if (a1[i].lessThan(a0[i])) {
            return +1;
        }
    }
    return 0;
}

export class Parallelepiped3 {
    // Public member access.
    center: Vector;
    axis: Vector[];

    // The port of the default constructor, which sets the center to (0,0,0),
    // axis[0] to (1,0,0), axis[1] to (0,1,0), and axis[2] to (0,0,1).
    constructor() {
        this.center = new Vector(3);
        this.axis = [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)];
    }

    // The port of 'Parallelepiped3(inCenter, inAxis)'. The axes must form a
    // right-handed basis. The axes do not have to be orthogonal. The axis
    // lengths do not have to be unit length. The vectors are copied, matching
    // C++ value semantics.
    static fromCenterAxis(inCenter: Vector,
        inAxis: readonly Vector[]): Parallelepiped3 {
        logAssert(inCenter.size === 3 && inAxis.length === 3
            && inAxis[0].size === 3 && inAxis[1].size === 3
            && inAxis[2].size === 3, 'Parallelepiped3: mismatched sizes.');
        logAssert(dotCross(inAxis[0], inAxis[1], inAxis[2]) > 0,
            'The axes must form a right-handed basis.');
        const parallelepiped = new Parallelepiped3();
        parallelepiped.center = inCenter.clone();
        parallelepiped.axis = [inAxis[0].clone(), inAxis[1].clone(),
            inAxis[2].clone()];
        return parallelepiped;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Parallelepiped3 {
        const parallelepiped = new Parallelepiped3();
        parallelepiped.center = this.center.clone();
        parallelepiped.axis = [this.axis[0].clone(), this.axis[1].clone(),
            this.axis[2].clone()];
        return parallelepiped;
    }

    // If index i has the bit pattern i = b[2]b[1]b[0], then
    //   vertices[i] = center + sum_{d=0}^{2} sign[d] * axis[d]
    // where sign[d] = 2*b[d] - 1. (The upstream comment claims the vertices
    // are listed in counterclockwise order, but the ordering is the
    // bit-pattern order shown here. The port preserves the upstream
    // ordering.)
    getVertices(): Vector[] {
        const a0 = this.axis[0];
        const a1 = this.axis[1];
        const a2 = this.axis[2];
        const c = this.center;
        return [
            sub(sub(sub(c, a0), a1), a2),
            sub(sub(add(c, a0), a1), a2),
            sub(add(sub(c, a0), a1), a2),
            sub(add(add(c, a0), a1), a2),
            add(sub(sub(c, a0), a1), a2),
            add(sub(add(c, a0), a1), a2),
            add(add(sub(c, a0), a1), a2),
            add(add(add(c, a0), a1), a2)
        ];
    }

    // Comparisons to support sorted containers.
    equals(other: Parallelepiped3): boolean {
        return this.center.equals(other.center)
            && compareAxes(this.axis, other.axis) === 0;
    }

    notEquals(other: Parallelepiped3): boolean {
        return !this.equals(other);
    }

    lessThan(other: Parallelepiped3): boolean {
        if (this.center.lessThan(other.center)) {
            return true;
        }

        if (this.center.greaterThan(other.center)) {
            return false;
        }

        return compareAxes(this.axis, other.axis) < 0;
    }

    lessThanOrEqual(other: Parallelepiped3): boolean {
        return !other.lessThan(this);
    }

    greaterThan(other: Parallelepiped3): boolean {
        return other.lessThan(this);
    }

    greaterThanOrEqual(other: Parallelepiped3): boolean {
        return !this.lessThan(other);
    }
}
