// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Circle3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The circle is the intersection of the sphere |X-C|^2 = r^2 and the plane
// Dot(N,X-C) = 0, where C is the circle center, r is the radius, and N is a
// unit-length plane normal.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 3D vectors.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export class Circle3 {
    // Public member access.
    center: Vector;
    normal: Vector;
    radius: number;

    // The port of the default constructor, which sets the center to (0,0,0),
    // the normal to (0,0,1) and the radius to 1.
    constructor() {
        this.center = new Vector(3);
        this.normal = Vector.unit(3, 2);
        this.radius = 1;
    }

    // The port of 'Circle3(inCenter, inNormal, inRadius)'. The vectors are
    // copied, matching C++ value semantics.
    static fromCenterNormalRadius(inCenter: Vector, inNormal: Vector,
        inRadius: number): Circle3 {
        logAssert(inCenter.size === 3 && inNormal.size === 3,
            'Circle3: mismatched sizes.');
        const circle = new Circle3();
        circle.center = inCenter.clone();
        circle.normal = inNormal.clone();
        circle.radius = inRadius;
        return circle;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Circle3 {
        return Circle3.fromCenterNormalRadius(this.center, this.normal,
            this.radius);
    }

    // Comparisons to support sorted containers.
    equals(circle: Circle3): boolean {
        return this.center.equals(circle.center)
            && this.normal.equals(circle.normal)
            && this.radius === circle.radius;
    }

    notEquals(circle: Circle3): boolean {
        return !this.equals(circle);
    }

    lessThan(circle: Circle3): boolean {
        if (this.center.lessThan(circle.center)) {
            return true;
        }

        if (this.center.greaterThan(circle.center)) {
            return false;
        }

        if (this.normal.lessThan(circle.normal)) {
            return true;
        }

        if (this.normal.greaterThan(circle.normal)) {
            return false;
        }

        return this.radius < circle.radius;
    }

    lessThanOrEqual(circle: Circle3): boolean {
        return !circle.lessThan(this);
    }

    greaterThan(circle: Circle3): boolean {
        return circle.lessThan(this);
    }

    greaterThanOrEqual(circle: Circle3): boolean {
        return !this.lessThan(circle);
    }
}
