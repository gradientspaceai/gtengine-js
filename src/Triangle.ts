// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Triangle.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The triangle is represented as an array of three vertices. The dimension N
// must be 2 or larger.
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (runtime dimension, 'new Triangle(n)' for the default
// constructor, named static factories that copy their Vector arguments,
// comparison methods). The C++ 'std::array<Vector<N, Real>, 3> v' becomes a
// Vector[] of length 3; its comparisons are lexicographic over the elements,
// as std::array's are.

import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

// Lexicographic comparison of the three-element vertex arrays (the port of
// std::array's relational operators). Returns -1, 0 or +1.
function compareVertices(v0: readonly Vector[], v1: readonly Vector[]): number {
    for (let i = 0; i < 3; ++i) {
        if (v0[i].lessThan(v1[i])) {
            return -1;
        }
        if (v1[i].lessThan(v0[i])) {
            return +1;
        }
    }
    return 0;
}

export class Triangle {
    // Public member access: the three vertices.
    v: Vector[];

    // The port of the default constructor, which sets the vertices to
    // (0,...,0), (1,0,...,0) and (0,1,0,...,0). The dimension N of the C++
    // template is a constructor argument here.
    constructor(n: number) {
        this.v = [Vector.zero(n), Vector.unit(n, 0), Vector.unit(n, 1)];
    }

    // The port of 'Triangle(v0, v1, v2)'. The vectors are copied, matching
    // C++ value semantics.
    static fromVertices(v0: Vector, v1: Vector, v2: Vector): Triangle {
        logAssert(v0.size === v1.size && v0.size === v2.size,
            'Triangle: mismatched sizes.');
        const triangle = new Triangle(v0.size);
        triangle.v = [v0.clone(), v1.clone(), v2.clone()];
        return triangle;
    }

    // The port of 'Triangle(std::array<Vector<N, Real>, 3> const& inV)'.
    static fromVertexArray(inV: readonly Vector[]): Triangle {
        logAssert(inV.length === 3, 'Triangle: invalid number of vertices.');
        return Triangle.fromVertices(inV[0], inV[1], inV[2]);
    }

    // The dimension N of the space containing the triangle.
    get dimension(): number {
        return this.v[0].size;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Triangle {
        return Triangle.fromVertices(this.v[0], this.v[1], this.v[2]);
    }

    // Comparisons to support sorted containers.
    equals(triangle: Triangle): boolean {
        return compareVertices(this.v, triangle.v) === 0;
    }

    notEquals(triangle: Triangle): boolean {
        return !this.equals(triangle);
    }

    lessThan(triangle: Triangle): boolean {
        return compareVertices(this.v, triangle.v) < 0;
    }

    lessThanOrEqual(triangle: Triangle): boolean {
        return compareVertices(this.v, triangle.v) <= 0;
    }

    greaterThan(triangle: Triangle): boolean {
        return compareVertices(this.v, triangle.v) > 0;
    }

    greaterThanOrEqual(triangle: Triangle): boolean {
        return compareVertices(this.v, triangle.v) >= 0;
    }
}

// Aliases for convenience (the ports of the upstream template aliases).
export type Triangle2 = Triangle;
export type Triangle3 = Triangle;
