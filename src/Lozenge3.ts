// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Lozenge3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A lozenge is the set of points that are equidistant from a rectangle, the
// common distance called the radius.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their arguments, comparison
// methods). The class is not templated on the dimension upstream, so the
// default constructor takes no arguments and builds a 3D rectangle.

import { logAssert } from './Logger';
import { Rectangle } from './Rectangle';

export class Lozenge3 {
    // Public member access.
    rectangle: Rectangle;
    radius: number;

    // The port of the default constructor, which sets the rectangle to have
    // origin (0,0,0), axes (1,0,0) and (0,1,0), and both extents 1. The
    // default radius is 1.
    constructor() {
        this.rectangle = new Rectangle(3);
        this.radius = 1;
    }

    // The port of 'Lozenge3(inRectangle, inRadius)'. The rectangle is copied,
    // matching C++ value semantics.
    static fromRectangleRadius(inRectangle: Rectangle,
        inRadius: number): Lozenge3 {
        logAssert(inRectangle.dimension === 3, 'Lozenge3: mismatched sizes.');
        const lozenge = new Lozenge3();
        lozenge.rectangle = inRectangle.clone();
        lozenge.radius = inRadius;
        return lozenge;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Lozenge3 {
        return Lozenge3.fromRectangleRadius(this.rectangle, this.radius);
    }

    // Comparisons to support sorted containers.
    equals(other: Lozenge3): boolean {
        return this.rectangle.equals(other.rectangle)
            && this.radius === other.radius;
    }

    notEquals(other: Lozenge3): boolean {
        return !this.equals(other);
    }

    lessThan(other: Lozenge3): boolean {
        if (this.rectangle.lessThan(other.rectangle)) {
            return true;
        }

        if (this.rectangle.greaterThan(other.rectangle)) {
            return false;
        }

        return this.radius < other.radius;
    }

    lessThanOrEqual(other: Lozenge3): boolean {
        return !other.lessThan(this);
    }

    greaterThan(other: Lozenge3): boolean {
        return other.lessThan(this);
    }

    greaterThanOrEqual(other: Lozenge3): boolean {
        return !this.lessThan(other);
    }
}
