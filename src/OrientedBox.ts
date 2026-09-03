// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OrientedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A box has center C, axis directions U[i], and extents e[i]. The set
// {U[0],...,U[N-1]} is orthonormal, which means the vectors are unit-length
// and mutually perpendicular. The extents are nonnegative; zero is allowed,
// meaning the box is degenerate in the corresponding direction. A point X is
// represented in box coordinates by X = C + y[0]*U[0] + y[1]*U[1]. This point
// is inside or on the box whenever |y[i]| <= e[i] for all i.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new OrientedBox(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods). The C++ 'std::array<Vector<N, T>, N> axis' becomes a
// Vector[] of length N; its comparisons are lexicographic over the elements,
// as std::array's are.

import { logAssert } from './Logger.js';
import { Vector, add, sub, mul } from './Vector.js';

// Lexicographic comparison of equal-length vector arrays (the port of
// std::array's relational operators). Returns -1, 0 or +1.
function compareAxes(a0: readonly Vector[], a1: readonly Vector[]): number {
    logAssert(a0.length === a1.length, 'OrientedBox: mismatched sizes.');
    for (let i = 0; i < a0.length; ++i) {
        if (a0[i].lessThan(a1[i])) {
            return -1;
        }
        if (a1[i].lessThan(a0[i])) {
            return +1;
        }
    }
    return 0;
}

export class OrientedBox {
    // Public member access. It is required that extent[i] >= 0.
    center: Vector;
    axis: Vector[];
    extent: Vector;

    // The port of the default constructor, which sets the center to
    // (0,...,0), axis d to Vector<N,T>::Unit(d) and extent d to +1. The
    // dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        this.center = new Vector(n);
        this.axis = new Array<Vector>(n);
        this.extent = Vector.filled(n, 1);
        for (let i = 0; i < n; ++i) {
            this.axis[i] = Vector.unit(n, i);
        }
    }

    // The port of 'OrientedBox(inCenter, inAxis, inExtent)'. The vectors are
    // copied, matching C++ value semantics.
    static fromCenterAxisExtent(inCenter: Vector, inAxis: readonly Vector[],
        inExtent: Vector): OrientedBox {
        const n = inCenter.size;
        logAssert(inAxis.length === n && inExtent.size === n,
            'OrientedBox: mismatched sizes.');
        const box = new OrientedBox(n);
        box.center = inCenter.clone();
        box.extent = inExtent.clone();
        for (let i = 0; i < n; ++i) {
            logAssert(inAxis[i].size === n, 'OrientedBox: mismatched sizes.');
            box.axis[i] = inAxis[i].clone();
        }
        return box;
    }

    // The dimension N of the box.
    get dimension(): number {
        return this.center.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): OrientedBox {
        return OrientedBox.fromCenterAxisExtent(this.center, this.axis,
            this.extent);
    }

    // Compute the vertices of the box. If index i has the bit pattern
    // i = b[N-1]...b[0], then
    // vertex[i] = center + sum_{d=0}^{N-1} sign[d] * extent[d] * axis[d]
    // where sign[d] = 2*b[d] - 1.
    getVertices(): Vector[] {
        const n = this.dimension;
        const product: Vector[] = new Array<Vector>(n);
        for (let d = 0; d < n; ++d) {
            product[d] = mul(this.extent.values[d], this.axis[d]);
        }

        const imax = 1 << n;
        const vertex: Vector[] = new Array<Vector>(imax);
        for (let i = 0; i < imax; ++i) {
            let v = this.center.clone();
            for (let d = 0, mask = 1; d < n; ++d, mask <<= 1) {
                v = (i & mask) > 0 ? add(v, product[d]) : sub(v, product[d]);
            }
            vertex[i] = v;
        }
        return vertex;
    }

    // Comparisons to support sorted containers.
    equals(box: OrientedBox): boolean {
        return this.center.equals(box.center)
            && compareAxes(this.axis, box.axis) === 0
            && this.extent.equals(box.extent);
    }

    notEquals(box: OrientedBox): boolean {
        return !this.equals(box);
    }

    lessThan(box: OrientedBox): boolean {
        if (this.center.lessThan(box.center)) {
            return true;
        }

        if (this.center.greaterThan(box.center)) {
            return false;
        }

        const axisOrder = compareAxes(this.axis, box.axis);
        if (axisOrder < 0) {
            return true;
        }

        if (axisOrder > 0) {
            return false;
        }

        return this.extent.lessThan(box.extent);
    }

    lessThanOrEqual(box: OrientedBox): boolean {
        return !box.lessThan(this);
    }

    greaterThan(box: OrientedBox): boolean {
        return box.lessThan(this);
    }

    greaterThanOrEqual(box: OrientedBox): boolean {
        return !this.lessThan(box);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type OrientedBox2 = OrientedBox;
export type OrientedBox3 = OrientedBox;
