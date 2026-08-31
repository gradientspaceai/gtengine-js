// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Hypersphere.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The hypersphere is represented as |X-C| = R where C is the center and R is
// the radius. The hypersphere is a circle for dimension 2 or a sphere for
// dimension 3.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Hypersphere(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods).

import { Vector } from './Vector';

export class Hypersphere {
    // Public member access.
    center: Vector;
    radius: number;

    // The port of the default constructor, which sets the center to
    // (0,...,0) and the radius to 1. The dimension N of the C++ template is
    // a constructor argument here.
    constructor(n: number) {
        this.center = new Vector(n);
        this.radius = 1;
    }

    // The port of 'Hypersphere(inCenter, inRadius)'. The vector is copied,
    // matching C++ value semantics.
    static fromCenterRadius(inCenter: Vector, inRadius: number): Hypersphere {
        const hypersphere = new Hypersphere(inCenter.size);
        hypersphere.center = inCenter.clone();
        hypersphere.radius = inRadius;
        return hypersphere;
    }

    // The dimension N of the hypersphere.
    get dimension(): number {
        return this.center.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Hypersphere {
        return Hypersphere.fromCenterRadius(this.center, this.radius);
    }

    // Comparisons to support sorted containers.
    equals(hypersphere: Hypersphere): boolean {
        return this.center.equals(hypersphere.center)
            && this.radius === hypersphere.radius;
    }

    notEquals(hypersphere: Hypersphere): boolean {
        return !this.equals(hypersphere);
    }

    lessThan(hypersphere: Hypersphere): boolean {
        if (this.center.lessThan(hypersphere.center)) {
            return true;
        }

        if (this.center.greaterThan(hypersphere.center)) {
            return false;
        }

        return this.radius < hypersphere.radius;
    }

    lessThanOrEqual(hypersphere: Hypersphere): boolean {
        return !hypersphere.lessThan(this);
    }

    greaterThan(hypersphere: Hypersphere): boolean {
        return hypersphere.lessThan(this);
    }

    greaterThanOrEqual(hypersphere: Hypersphere): boolean {
        return !this.lessThan(hypersphere);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Circle2 = Hypersphere;
export type Sphere3 = Hypersphere;
