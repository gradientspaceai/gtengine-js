// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GVector.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream GVector<Real> duplicates the entire Vector<N, Real>
// API with a std::vector-backed, resizable tuple. The port's Vector class
// already has a runtime dimension, so GVector is a subclass that adds only
// the resizing behavior (setSize) and the std::vector comparison semantics.
// All of Vector's module free functions (negate, add, sub, mul, div, compMul,
// compDiv, dot, length, normalize, orthonormalize, computeExtremes, hlift,
// hproject, lift, project) apply to GVector unchanged and are NOT duplicated
// here; import them from './Vector' (upstream GVector.h re-implements them
// with identical algorithms).
//
// Behavioral differences from Vector, matching upstream GVector.h:
// - The default constructor produces a size-0 tuple (upstream GVector()).
// - setSize(size) resizes the tuple; new elements are zero (the port of
//   std::vector::resize, which value-initializes).
// - Comparisons use std::vector semantics: vectors of different sizes may be
//   compared (lexicographic with the size as tie-breaker) instead of
//   throwing as the fixed-size Vector comparisons do.
//
// Known deviations (documented rather than duplicated, since the free
// function names are owned by Vector.ts):
// - 'v / scalar' with scalar 0 throws upstream ("Division by zero.") but the
//   shared div() from Vector.ts returns the zero vector (Vector.h behavior).
// - HProject/Project of a size<=1 tuple return an empty GVector upstream;
//   the shared hproject()/project() throw as Vector.h's static_assert does.

import { logAssert } from './Logger';
import { Vector } from './Vector';

export class GVector extends Vector {
    // Create a tuple of the given size (default 0, matching the upstream
    // default constructor) with all components 0. Upstream leaves the
    // elements uninitialized; the port zero-fills for safety. When d is
    // provided, component d is set to 1 (the port of GVector(size, d)); if d
    // is invalid, the zero vector is created.
    constructor(size: number = 0, d?: number) {
        logAssert(size >= 0, 'Invalid size.');
        super(size);
        if (d !== undefined) {
            this.makeUnit(d);
        }
    }

    // Create a GVector that copies the input values.
    static override fromArray(values: readonly number[]): GVector {
        const result = new GVector(values.length);
        for (let i = 0; i < values.length; ++i) {
            result.values[i] = values[i];
        }
        return result;
    }

    // All components are 0.
    static override zero(size: number): GVector {
        return new GVector(size);
    }

    // Component d is 1, all others are 0.
    static override unit(size: number, d: number): GVector {
        return new GVector(size, d);
    }

    // Resize the tuple. Upstream SetSize does not initialize new elements
    // (std::vector::resize value-initializes them to 0); the port sets them
    // to 0. Existing elements are preserved.
    setSize(size: number): void {
        logAssert(size >= 0, 'Invalid size.');
        const oldSize = this.values.length;
        this.values.length = size;
        for (let i = oldSize; i < size; ++i) {
            this.values[i] = 0;
        }
    }

    // A deep copy (the port of C++ copy construction/assignment).
    override clone(): GVector {
        return GVector.fromArray(this.values);
    }

    // Comparisons (for use by sorted containers). Unlike the fixed-size
    // Vector comparisons, which throw on mismatched sizes, these follow
    // std::vector: lexicographic over the common prefix, with the shorter
    // vector ordered first when the prefix ties; equality requires equal
    // sizes.
    private vectorCompare(vec: Vector): number {
        const n0 = this.values.length;
        const n1 = vec.values.length;
        const n = Math.min(n0, n1);
        for (let i = 0; i < n; ++i) {
            if (this.values[i] < vec.values[i]) {
                return -1;
            }
            if (this.values[i] > vec.values[i]) {
                return +1;
            }
        }
        return n0 < n1 ? -1 : (n0 > n1 ? +1 : 0);
    }

    override equals(vec: Vector): boolean {
        return this.vectorCompare(vec) === 0;
    }

    override notEquals(vec: Vector): boolean {
        return this.vectorCompare(vec) !== 0;
    }

    override lessThan(vec: Vector): boolean {
        return this.vectorCompare(vec) < 0;
    }

    override lessThanOrEqual(vec: Vector): boolean {
        return this.vectorCompare(vec) <= 0;
    }

    override greaterThan(vec: Vector): boolean {
        return this.vectorCompare(vec) > 0;
    }

    override greaterThanOrEqual(vec: Vector): boolean {
        return this.vectorCompare(vec) >= 0;
    }
}
