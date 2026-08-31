// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) CanonicalBox.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A canonical box has center at the origin and is aligned with the standard
// Euclidean basis vectors. It has E = (e[0],e[1],...,e[N-1]) with e[i] >= 0
// for all i. A zero extent is allowed, meaning the box is degenerate in the
// corresponding direction. A box point is X = (x[0],x[1],...,x[N-1]) with
// |x[i]| <= e[i] for all i.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new CanonicalBox(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods).

import { Vector } from './Vector';

export class CanonicalBox {
    // It is required that extent[i] >= 0.
    extent: Vector;

    // The port of the default constructor, which sets all members to zero.
    // The dimension N of the C++ template is a constructor argument here.
    constructor(n: number) {
        this.extent = new Vector(n);
    }

    // The port of 'CanonicalBox(Vector<N, T> const& inExtent)'. The vector is
    // copied, matching C++ value semantics.
    static fromExtent(inExtent: Vector): CanonicalBox {
        const box = new CanonicalBox(inExtent.size);
        box.extent = inExtent.clone();
        return box;
    }

    // The dimension N of the box.
    get dimension(): number {
        return this.extent.size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): CanonicalBox {
        return CanonicalBox.fromExtent(this.extent);
    }

    // Compute the vertices of the box. If index i has the bit pattern
    // i = b[N-1]...b[0], then the corner at index i is
    //   vertex[i] = center + sum_{d=0}^{N-1} sign[d]*extent[d]*axis[d]
    // where sign[d] = 2*b[d] - 1.
    getVertices(): Vector[] {
        const n = this.dimension;
        const imax = 1 << n;
        const vertex: Vector[] = new Array<Vector>(imax);
        for (let i = 0; i < imax; ++i) {
            const v = new Vector(n);
            for (let d = 0, mask = 1; d < n; ++d, mask <<= 1) {
                if ((i & mask) > 0) {
                    v.values[d] += this.extent.values[d];
                } else {
                    v.values[d] -= this.extent.values[d];
                }
            }
            vertex[i] = v;
        }
        return vertex;
    }

    // Comparisons to support sorted containers.
    equals(box: CanonicalBox): boolean {
        return this.extent.equals(box.extent);
    }

    notEquals(box: CanonicalBox): boolean {
        return !this.equals(box);
    }

    lessThan(box: CanonicalBox): boolean {
        return this.extent.lessThan(box.extent);
    }

    lessThanOrEqual(box: CanonicalBox): boolean {
        return !box.lessThan(this);
    }

    greaterThan(box: CanonicalBox): boolean {
        return box.lessThan(this);
    }

    greaterThanOrEqual(box: CanonicalBox): boolean {
        return !this.lessThan(box);
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type CanonicalBox2 = CanonicalBox;
export type CanonicalBox3 = CanonicalBox;
