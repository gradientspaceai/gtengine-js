// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Frustum3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Orthogonal frustum. Let E be the origin, D be the direction vector, U be
// the up vector, and R be the right vector. Let u > 0 and r > 0 be the
// extents in the U and R directions, respectively. Let n and f be the extents
// in the D direction with 0 < n < f. The four corners of the frustum in the
// near plane are E + n*D + s0*u*U + s1*r*R where |s0| = |s1| = 1 (four
// choices). The four corners of the frustum in the far plane are
// E + f*D + (f/n)*(s0*u*U + s1*r*R) where |s0| = |s1| = 1 (four choices).
//
// Port notes: see AlignedBox.ts for the shared geometric-primitive
// conventions (named static factories that copy their Vector arguments,
// comparison methods). The class is not templated on the dimension upstream,
// so the default constructor takes no arguments and builds 3D vectors.
// 'ComputeVertices' returns a new array of 8 vertices instead of filling a
// caller-supplied std::array.

import { logAssert } from './Logger.js';
import { Vector, add, sub, mul } from './Vector.js';

export class Frustum3 {
    // Public member access.
    origin: Vector;
    dVector: Vector;
    uVector: Vector;
    rVector: Vector;
    dMin: number;
    dMax: number;
    uBound: number;
    rBound: number;

    // Quantities derived from the constructor inputs.
    private mDRatio: number;
    private mMTwoUF: number;
    private mMTwoRF: number;

    // The port of the default constructor, which sets the origin (E) to
    // (0,0,0), dVector (D) to (0,0,1), uVector (U) to (0,1,0), rVector (R) to
    // (1,0,0), dMin (n) to 1, dMax (f) to 2, uBound (u) to 1, and rBound (r)
    // to 1.
    constructor() {
        this.origin = new Vector(3);
        this.dVector = Vector.unit(3, 2);
        this.uVector = Vector.unit(3, 1);
        this.rVector = Vector.unit(3, 0);
        this.dMin = 1;
        this.dMax = 2;
        this.uBound = 1;
        this.rBound = 1;
        this.mDRatio = 0;
        this.mMTwoUF = 0;
        this.mMTwoRF = 0;
        this.update();
    }

    // The port of the eight-argument constructor. The vectors are copied,
    // matching C++ value semantics.
    static fromParameters(inOrigin: Vector, inDVector: Vector,
        inUVector: Vector, inRVector: Vector, inDMin: number, inDMax: number,
        inUBound: number, inRBound: number): Frustum3 {
        logAssert(inOrigin.size === 3 && inDVector.size === 3
            && inUVector.size === 3 && inRVector.size === 3,
            'Frustum3: mismatched sizes.');
        const frustum = new Frustum3();
        frustum.origin = inOrigin.clone();
        frustum.dVector = inDVector.clone();
        frustum.uVector = inUVector.clone();
        frustum.rVector = inRVector.clone();
        frustum.dMin = inDMin;
        frustum.dMax = inDMax;
        frustum.uBound = inUBound;
        frustum.rBound = inRBound;
        frustum.update();
        return frustum;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Frustum3 {
        return Frustum3.fromParameters(this.origin, this.dVector,
            this.uVector, this.rVector, this.dMin, this.dMax, this.uBound,
            this.rBound);
    }

    // The update() function must be called whenever changes are made to
    // dMin, dMax, uBound or rBound. The values mDRatio, mMTwoUF and mMTwoRF
    // are dependent on the changes, so call the get*() accessors only after
    // the update() call.
    update(): void {
        this.mDRatio = this.dMax / this.dMin;
        this.mMTwoUF = -2 * this.uBound * this.dMax;
        this.mMTwoRF = -2 * this.rBound * this.dMax;
    }

    getDRatio(): number {
        return this.mDRatio;
    }

    getMTwoUF(): number {
        return this.mMTwoUF;
    }

    getMTwoRF(): number {
        return this.mMTwoRF;
    }

    // The four near-plane vertices are indices 0..3 and the four far-plane
    // vertices are indices 4..7.
    computeVertices(): Vector[] {
        const dScaled = mul(this.dMin, this.dVector);
        const uScaled = mul(this.uBound, this.uVector);
        const rScaled = mul(this.rBound, this.rVector);

        const vertex: Vector[] = new Array<Vector>(8);
        vertex[0] = sub(sub(dScaled, uScaled), rScaled);
        vertex[1] = add(sub(dScaled, uScaled), rScaled);
        vertex[2] = add(add(dScaled, uScaled), rScaled);
        vertex[3] = sub(add(dScaled, uScaled), rScaled);

        for (let i = 0, ip = 4; i < 4; ++i, ++ip) {
            vertex[ip] = add(this.origin, mul(this.mDRatio, vertex[i]));
            vertex[i] = add(vertex[i], this.origin);
        }

        return vertex;
    }

    // Comparisons to support sorted containers.
    equals(frustum: Frustum3): boolean {
        return this.origin.equals(frustum.origin)
            && this.dVector.equals(frustum.dVector)
            && this.uVector.equals(frustum.uVector)
            && this.rVector.equals(frustum.rVector)
            && this.dMin === frustum.dMin
            && this.dMax === frustum.dMax
            && this.uBound === frustum.uBound
            && this.rBound === frustum.rBound;
    }

    notEquals(frustum: Frustum3): boolean {
        return !this.equals(frustum);
    }

    lessThan(frustum: Frustum3): boolean {
        if (this.origin.lessThan(frustum.origin)) {
            return true;
        }

        if (this.origin.greaterThan(frustum.origin)) {
            return false;
        }

        if (this.dVector.lessThan(frustum.dVector)) {
            return true;
        }

        if (this.dVector.greaterThan(frustum.dVector)) {
            return false;
        }

        if (this.uVector.lessThan(frustum.uVector)) {
            return true;
        }

        if (this.uVector.greaterThan(frustum.uVector)) {
            return false;
        }

        if (this.rVector.lessThan(frustum.rVector)) {
            return true;
        }

        if (this.rVector.greaterThan(frustum.rVector)) {
            return false;
        }

        if (this.dMin < frustum.dMin) {
            return true;
        }

        if (this.dMin > frustum.dMin) {
            return false;
        }

        if (this.dMax < frustum.dMax) {
            return true;
        }

        if (this.dMax > frustum.dMax) {
            return false;
        }

        if (this.uBound < frustum.uBound) {
            return true;
        }

        if (this.uBound > frustum.uBound) {
            return false;
        }

        return this.rBound < frustum.rBound;
    }

    lessThanOrEqual(frustum: Frustum3): boolean {
        return !frustum.lessThan(this);
    }

    greaterThan(frustum: Frustum3): boolean {
        return frustum.lessThan(this);
    }

    greaterThanOrEqual(frustum: Frustum3): boolean {
        return !this.lessThan(frustum);
    }
}
