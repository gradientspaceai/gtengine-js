// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Array2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Array2 class represents a 2-dimensional array stored in a contiguous
// flat array in lexicographical order (dimension 0 varies fastest).
//
// Port notes: upstream returns row pointers from operator[] so that elements
// are accessed as 'myArray[i1][i0]' with i0 the index in dimension 0 (bound0
// columns) and i1 the index in dimension 1 (bound1 rows). Pointer indirection
// does not port, so element access is 'myArray.get(i0, i1)' and
// 'myArray.set(i0, i1, value)' -- note the index order is (i0, i1), matching
// the constructor's bound order, so C++ 'a[i1][i0]' becomes 'a.get(i0, i1)'.
// As upstream, the storage may be owned (allocated by the constructor) or
// caller-owned (a flat array passed in, aliased rather than copied). Copy and
// move machinery is not needed in TypeScript and is not ported. Accesses are
// not bounds-checked, matching upstream.

import { logAssert } from './Logger.js';

export class Array2<T> {
    private mBound0: number;
    private mBound1: number;
    private mObjects: T[];

    // Construction. When 'objects' is provided, the flat array (of length
    // bound0 * bound1, in lexicographical order) is owned by the caller and
    // aliased. Otherwise storage is allocated as a sparse array (holes read
    // as undefined). NOTE: upstream's owning constructor value-initializes a
    // std::vector<T>, which zero-fills numeric T; the generic port cannot
    // choose a zero for arbitrary T, so callers must fill() or set() before
    // get() (accumulating with += on a fresh array yields NaN, see
    // IntpAkimaUniform3).
    constructor(bound0: number = 0, bound1: number = 0, objects?: T[]) {
        this.mBound0 = bound0;
        this.mBound1 = bound1;
        if (objects !== undefined) {
            logAssert(objects.length === bound0 * bound1,
                'Array2: objects.length must equal bound0 * bound1.');
            this.mObjects = objects;
        } else {
            this.mObjects = new Array<T>(bound0 * bound1);
        }
    }

    getBound0(): number {
        return this.mBound0;
    }

    getBound1(): number {
        return this.mBound1;
    }

    // Element access; C++ 'a[i1][i0]' is 'a.get(i0, i1)'.
    get(i0: number, i1: number): T {
        return this.mObjects[i0 + this.mBound0 * i1];
    }

    set(i0: number, i1: number, value: T): void {
        this.mObjects[i0 + this.mBound0 * i1] = value;
    }

    // Access to the flat storage (aliased, not a copy).
    data(): T[] {
        return this.mObjects;
    }

    // Set all elements to the given value.
    fill(value: T): void {
        this.mObjects.fill(value, 0, this.mBound0 * this.mBound1);
    }
}
