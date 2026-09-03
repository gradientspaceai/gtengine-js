// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Halfspace.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The halfspace is represented as Dot(N,X) >= c where N is a unit-length
// normal vector, c is the plane constant, and X is any point in space. The
// user must ensure that the normal vector is unit length.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Halfspace(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods).

import { Vector } from './Vector.js';

export class Halfspace {
    // Public member access.
    normal: Vector;
    constant: number;

    // The port of the default constructor, which sets the normal to
    // (0,...,0,1) and the constant to zero (halfspace x[N-1] >= 0). The
    // dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        this.normal = Vector.unit(n, n - 1);
        this.constant = 0;
    }

    // Specify N and c directly. The vector is copied, matching C++ value
    // semantics.
    static fromNormalConstant(inNormal: Vector, inConstant: number): Halfspace {
        const halfspace = new Halfspace(inNormal.size);
        halfspace.normal = inNormal.clone();
        halfspace.constant = inConstant;
        return halfspace;
    }

    // The dimension N of the halfspace.
    get dimension(): number {
        return this.normal.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Halfspace {
        return Halfspace.fromNormalConstant(this.normal, this.constant);
    }

    // Comparisons to support sorted containers.
    equals(halfspace: Halfspace): boolean {
        return this.normal.equals(halfspace.normal)
            && this.constant === halfspace.constant;
    }

    notEquals(halfspace: Halfspace): boolean {
        return !this.equals(halfspace);
    }

    lessThan(halfspace: Halfspace): boolean {
        if (this.normal.lessThan(halfspace.normal)) {
            return true;
        }

        if (this.normal.greaterThan(halfspace.normal)) {
            return false;
        }

        return this.constant < halfspace.constant;
    }

    lessThanOrEqual(halfspace: Halfspace): boolean {
        return !halfspace.lessThan(this);
    }

    greaterThan(halfspace: Halfspace): boolean {
        return halfspace.lessThan(this);
    }

    greaterThanOrEqual(halfspace: Halfspace): boolean {
        return !this.lessThan(halfspace);
    }
}

// Alias for convenience (the port of the upstream template alias).
export type Halfspace3 = Halfspace;
