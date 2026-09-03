// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Ray.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The ray is represented as P+t*D, where P is the ray origin, D is a
// unit-length direction vector, and t >= 0. The user must ensure that D is
// unit length.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Ray(n)' for the default constructor,
// named static factories that copy their Vector arguments, comparison
// methods).

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export class Ray {
    // Public member access. The direction must be unit length.
    origin: Vector;
    direction: Vector;

    // The port of the default constructor, which sets the origin to
    // (0,...,0) and the ray direction to (1,0,...,0). The dimension N of the
    // C++ template is a constructor argument here.
    constructor(n: number) {
        this.origin = new Vector(n);
        this.direction = Vector.unit(n, 0);
    }

    // The port of 'Ray(inOrigin, inDirection)'. The vectors are copied,
    // matching C++ value semantics.
    static fromOriginDirection(inOrigin: Vector, inDirection: Vector): Ray {
        logAssert(inOrigin.size === inDirection.size,
            'Ray: mismatched sizes.');
        const ray = new Ray(inOrigin.size);
        ray.origin = inOrigin.clone();
        ray.direction = inDirection.clone();
        return ray;
    }

    // The dimension N of the ray.
    get dimension(): number {
        return this.origin.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Ray {
        return Ray.fromOriginDirection(this.origin, this.direction);
    }

    // Comparisons to support sorted containers.
    equals(ray: Ray): boolean {
        return this.origin.equals(ray.origin)
            && this.direction.equals(ray.direction);
    }

    notEquals(ray: Ray): boolean {
        return !this.equals(ray);
    }

    lessThan(ray: Ray): boolean {
        if (this.origin.lessThan(ray.origin)) {
            return true;
        }

        if (this.origin.greaterThan(ray.origin)) {
            return false;
        }

        return this.direction.lessThan(ray.direction);
    }

    lessThanOrEqual(ray: Ray): boolean {
        return !ray.lessThan(this);
    }

    greaterThan(ray: Ray): boolean {
        return ray.lessThan(this);
    }

    greaterThanOrEqual(ray: Ray): boolean {
        return !this.lessThan(ray);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Ray2 = Ray;
export type Ray3 = Ray;
