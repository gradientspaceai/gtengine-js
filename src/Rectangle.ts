// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Rectangle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Points are R(s0,s1) = C + s0*A0 + s1*A1, where C is the center of the
// rectangle and A0 and A1 are unit-length and perpendicular axes. The
// parameters s0 and s1 are constrained by |s0| <= e0 and |s1| <= e1, where
// e0 > 0 and e1 > 0 are the extents of the rectangle.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Rectangle(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods). The rectangle lives in N-space but always has two
// axes and a 2-dimensional extent vector, matching upstream's
// 'std::array<Vector<N, Real>, 2> axis' and 'Vector<2, Real> extent'.

import { logAssert } from './Logger.js';
import { Vector, add, sub, mul } from './Vector.js';

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

export class Rectangle {
    // Public member access. The axes are vectors in N-space; the extent has
    // two components.
    center: Vector;
    axis: Vector[];
    extent: Vector;

    // The port of the default constructor, which sets the center to
    // (0,...,0), axis A0 to (1,0,...,0), axis A1 to (0,1,0,...,0) and both
    // extents to 1. The dimension N of the C++ template is a constructor
    // argument here.
    constructor(n: number) {
        this.center = new Vector(n);
        this.axis = [Vector.unit(n, 0), Vector.unit(n, 1)];
        this.extent = Vector.filled(2, 1);
    }

    // The port of 'Rectangle(inCenter, inAxis, inExtent)'. The vectors are
    // copied, matching C++ value semantics.
    static fromCenterAxisExtent(inCenter: Vector, inAxis: readonly Vector[],
        inExtent: Vector): Rectangle {
        const n = inCenter.size;
        logAssert(inAxis.length === 2 && inExtent.size === 2,
            'Rectangle: mismatched sizes.');
        logAssert(inAxis[0].size === n && inAxis[1].size === n,
            'Rectangle: mismatched sizes.');
        const rectangle = new Rectangle(n);
        rectangle.center = inCenter.clone();
        rectangle.axis = [inAxis[0].clone(), inAxis[1].clone()];
        rectangle.extent = inExtent.clone();
        return rectangle;
    }

    // The dimension N of the space containing the rectangle.
    get dimension(): number {
        return this.center.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Rectangle {
        return Rectangle.fromCenterAxisExtent(this.center, this.axis,
            this.extent);
    }

    // Compute the vertices of the rectangle. If index i has the bit pattern
    // i = b[1]b[0], then
    //   vertex[i] = center + sum_{d=0}^{1} sign[d] * extent[d] * axis[d]
    // where sign[d] = 2*b[d] - 1.
    getVertices(): Vector[] {
        const product0 = mul(this.extent.values[0], this.axis[0]);
        const product1 = mul(this.extent.values[1], this.axis[1]);
        const sum = add(product0, product1);
        const dif = sub(product0, product1);

        return [
            sub(this.center, sum),
            add(this.center, dif),
            sub(this.center, dif),
            add(this.center, sum)
        ];
    }

    // Comparisons to support sorted containers.
    equals(rectangle: Rectangle): boolean {
        if (this.center.notEquals(rectangle.center)) {
            return false;
        }

        for (let i = 0; i < 2; ++i) {
            if (this.axis[i].notEquals(rectangle.axis[i])) {
                return false;
            }
        }

        for (let i = 0; i < 2; ++i) {
            if (this.extent.values[i] !== rectangle.extent.values[i]) {
                return false;
            }
        }

        return true;
    }

    notEquals(rectangle: Rectangle): boolean {
        return !this.equals(rectangle);
    }

    lessThan(rectangle: Rectangle): boolean {
        if (this.center.lessThan(rectangle.center)) {
            return true;
        }

        if (this.center.greaterThan(rectangle.center)) {
            return false;
        }

        const axisOrder = compareAxes(this.axis, rectangle.axis);
        if (axisOrder < 0) {
            return true;
        }

        if (axisOrder > 0) {
            return false;
        }

        return this.extent.lessThan(rectangle.extent);
    }

    lessThanOrEqual(rectangle: Rectangle): boolean {
        return !rectangle.lessThan(this);
    }

    greaterThan(rectangle: Rectangle): boolean {
        return rectangle.lessThan(this);
    }

    greaterThanOrEqual(rectangle: Rectangle): boolean {
        return !this.lessThan(rectangle);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Rectangle2 = Rectangle;
export type Rectangle3 = Rectangle;
