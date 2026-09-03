// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Line.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The line is represented by P+t*D, where P is an origin point, D is a
// unit-length direction vector, and t is any real number. The user must
// ensure that D is unit length.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Line(n)' for the default constructor,
// named static factories that copy their Vector arguments, comparison
// methods).

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export class Line {
    // Public member access. The direction must be unit length.
    origin: Vector;
    direction: Vector;

    // The port of the default constructor, which sets the origin to
    // (0,...,0) and the line direction to (1,0,...,0). The dimension N of the
    // C++ template is a constructor argument here.
    constructor(n: number) {
        this.origin = new Vector(n);
        this.direction = Vector.unit(n, 0);
    }

    // The port of 'Line(inOrigin, inDirection)'. The vectors are copied,
    // matching C++ value semantics.
    static fromOriginDirection(inOrigin: Vector, inDirection: Vector): Line {
        logAssert(inOrigin.size === inDirection.size,
            'Line: mismatched sizes.');
        const line = new Line(inOrigin.size);
        line.origin = inOrigin.clone();
        line.direction = inDirection.clone();
        return line;
    }

    // The dimension N of the line.
    get dimension(): number {
        return this.origin.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Line {
        return Line.fromOriginDirection(this.origin, this.direction);
    }

    // Comparisons to support sorted containers.
    equals(line: Line): boolean {
        return this.origin.equals(line.origin)
            && this.direction.equals(line.direction);
    }

    notEquals(line: Line): boolean {
        return !this.equals(line);
    }

    lessThan(line: Line): boolean {
        if (this.origin.lessThan(line.origin)) {
            return true;
        }

        if (this.origin.greaterThan(line.origin)) {
            return false;
        }

        return this.direction.lessThan(line.direction);
    }

    lessThanOrEqual(line: Line): boolean {
        return !line.lessThan(this);
    }

    greaterThan(line: Line): boolean {
        return line.lessThan(this);
    }

    greaterThanOrEqual(line: Line): boolean {
        return !this.lessThan(line);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Line2 = Line;
export type Line3 = Line;
