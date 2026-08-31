// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) AlignedBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The box is aligned with the standard coordinate axes, which allows us to
// represent it using minimum and maximum values along each axis. Some
// algorithms prefer the centered representation that is used for oriented
// boxes. The center is C and the extents are the half-lengths in each
// coordinate-axis direction.
//
// Port notes (the geometric-primitive precedent, shared by the B29 batch):
// - Upstream 'template <int32_t N, typename T> class AlignedBox' becomes a
//   class with a runtime dimension, over the ported Vector.
// - The default constructor takes the dimension: 'new AlignedBox(n)' is the
//   port of 'AlignedBox<N, T>()'.
// - The value constructors become named static factories that copy their
//   Vector arguments, preserving C++ value semantics; 'clone()' is the port
//   of C++ copy construction.
// - Comparison operators become the methods equals, notEquals, lessThan,
//   lessThanOrEqual, greaterThan, greaterThanOrEqual.
// - Out-parameters become returned object literals ('GetCenteredForm') and
//   filled output arrays become returned arrays ('GetVertices').

import { logAssert } from './Logger';
import { Vector, add, sub, mul } from './Vector';

export class AlignedBox {
    // Public member access. It is required that min[i] <= max[i].
    min: Vector;
    max: Vector;

    // The port of the default constructor, which sets the minimum values to
    // -1 and the maximum values to +1. The dimension N of the C++ template
    // is a constructor argument here.
    constructor(n: number) {
        this.min = Vector.filled(n, -1);
        this.max = Vector.filled(n, 1);
    }

    // Please ensure that inMin[i] <= inMax[i] for all i. The vectors are
    // copied, matching the C++ member-by-member assignment.
    static fromMinMax(inMin: Vector, inMax: Vector): AlignedBox {
        logAssert(inMin.size === inMax.size, 'AlignedBox: mismatched sizes.');
        const box = new AlignedBox(inMin.size);
        box.min = inMin.clone();
        box.max = inMax.clone();
        return box;
    }

    // The dimension N of the box.
    get dimension(): number {
        return this.min.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): AlignedBox {
        return AlignedBox.fromMinMax(this.min, this.max);
    }

    // Compute the centered representation. NOTE: If you set the minimum and
    // maximum values, compute C and extents, and then recompute the minimum
    // and maximum values, the numerical round-off errors can lead to results
    // different from what you started with.
    getCenteredForm(): { center: Vector, extent: Vector } {
        const half = 0.5;
        return {
            center: mul(add(this.max, this.min), half),
            extent: mul(sub(this.max, this.min), half)
        };
    }

    // Compute the vertices of the box. If index i has the bit pattern
    // i = b[N-1]...b[0], then the corner at index i is vertex[i], where
    // vertex[i][d] = min[d] when b[d] = 0 or vertex[i][d] = max[d] when
    // b[d] = 1.
    getVertices(): Vector[] {
        const n = this.dimension;
        const imax = 1 << n;
        const vertex: Vector[] = new Array<Vector>(imax);
        for (let i = 0; i < imax; ++i) {
            const v = new Vector(n);
            for (let d = 0, mask = 1; d < n; ++d, mask <<= 1) {
                v.values[d] = ((i & mask) > 0 ? this.max : this.min).values[d];
            }
            vertex[i] = v;
        }
        return vertex;
    }

    // Comparisons to support sorted containers.
    equals(box: AlignedBox): boolean {
        return this.min.equals(box.min) && this.max.equals(box.max);
    }

    notEquals(box: AlignedBox): boolean {
        return !this.equals(box);
    }

    lessThan(box: AlignedBox): boolean {
        if (this.min.lessThan(box.min)) {
            return true;
        }

        if (this.min.greaterThan(box.min)) {
            return false;
        }

        return this.max.lessThan(box.max);
    }

    lessThanOrEqual(box: AlignedBox): boolean {
        return !box.lessThan(this);
    }

    greaterThan(box: AlignedBox): boolean {
        return box.lessThan(this);
    }

    greaterThanOrEqual(box: AlignedBox): boolean {
        return !this.lessThan(box);
    }
}

// Aliases for convenience (the ports of the upstream template aliases). The
// dimension is a runtime value, so these are documentation aliases for the
// single AlignedBox type rather than distinct types.
export type AlignedBox2 = AlignedBox;
export type AlignedBox3 = AlignedBox;
