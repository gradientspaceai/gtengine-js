// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Array4.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Array4 class represents a 4-dimensional array stored in a contiguous
// flat array in lexicographical order (dimension 0 varies fastest).
//
// Port notes: see Array2.ts. Upstream element access 'a[i3][i2][i1][i0]'
// becomes 'a.get(i0, i1, i2, i3)' -- the index order matches the
// constructor's bound order. The flat index is
// i0 + bound0 * (i1 + bound1 * (i2 + bound2 * i3)). Accesses are not
// bounds-checked, matching upstream.

import { logAssert } from './Logger.js';

export class Array4<T> {
    private mBound0: number;
    private mBound1: number;
    private mBound2: number;
    private mBound3: number;
    private mObjects: T[];

    // Construction. When 'objects' is provided, the flat array (of length
    // bound0 * bound1 * bound2 * bound3, in lexicographical order) is owned
    // by the caller and aliased. Otherwise storage is allocated and left
    // uninitialized as upstream does for native types; use fill() or set()
    // before get().
    constructor(bound0: number = 0, bound1: number = 0, bound2: number = 0,
        bound3: number = 0, objects?: T[]) {
        this.mBound0 = bound0;
        this.mBound1 = bound1;
        this.mBound2 = bound2;
        this.mBound3 = bound3;
        if (objects !== undefined) {
            logAssert(objects.length === bound0 * bound1 * bound2 * bound3,
                'Array4: objects.length must equal bound0 * bound1 * bound2 * bound3.');
            this.mObjects = objects;
        } else {
            this.mObjects = new Array<T>(bound0 * bound1 * bound2 * bound3);
        }
    }

    getBound0(): number {
        return this.mBound0;
    }

    getBound1(): number {
        return this.mBound1;
    }

    getBound2(): number {
        return this.mBound2;
    }

    getBound3(): number {
        return this.mBound3;
    }

    // Element access; C++ 'a[i3][i2][i1][i0]' is 'a.get(i0, i1, i2, i3)'.
    get(i0: number, i1: number, i2: number, i3: number): T {
        return this.mObjects[i0 + this.mBound0 * (i1 + this.mBound1 * (i2 + this.mBound2 * i3))];
    }

    set(i0: number, i1: number, i2: number, i3: number, value: T): void {
        this.mObjects[i0 + this.mBound0 * (i1 + this.mBound1 * (i2 + this.mBound2 * i3))] = value;
    }

    // Access to the flat storage (aliased, not a copy).
    data(): T[] {
        return this.mObjects;
    }

    // Set all elements to the given value.
    fill(value: T): void {
        this.mObjects.fill(value, 0, this.mBound0 * this.mBound1 * this.mBound2 * this.mBound3);
    }
}
