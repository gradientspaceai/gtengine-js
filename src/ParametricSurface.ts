// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ParametricSurface.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Abstract base class for a parameterized surface X(u,v). The parametric
// domain is either rectangular or triangular. Valid (u,v) values for a
// rectangular domain satisfy
//   umin <= u <= umax,  vmin <= v <= vmax
// and valid (u,v) values for a triangular domain satisfy
//   umin <= u <= umax,  vmin <= v <= vmax,
//   (vmax-vmin)*(u-umin)+(umax-umin)*(v-vmax) <= 0
//
// Port notes: upstream 'template <int32_t N, typename Real>' becomes a
// runtime dimension (the N-tuple positions are 'Vector' of that dimension),
// so the protected constructor takes the dimension as its first argument.
// The 'jet' output of Evaluate is a 'Vector[]' of at least SUP_ORDER entries
// (see createJet). C++ 'operator bool' has no TS equivalent; it is ported as
// the isConstructed() accessor.

import { Vector, normalize } from './Vector.js';

export abstract class ParametricSurface {
    // The number of entries a 'jet' array must have: position X;
    // first-order derivatives dX/du, dX/dv; second-order derivatives
    // d2X/du2, d2X/dudv, d2X/dv2.
    static readonly SUP_ORDER = 6;

    protected mDimension: number;
    protected mUMin: number;
    protected mUMax: number;
    protected mVMin: number;
    protected mVMax: number;
    protected mRectangular: boolean;
    protected mConstructed: boolean;

    protected constructor(dimension: number, umin: number, umax: number,
        vmin: number, vmax: number, rectangular: boolean) {
        this.mDimension = dimension;
        this.mUMin = umin;
        this.mUMax = umax;
        this.mVMin = vmin;
        this.mVMax = vmax;
        this.mRectangular = rectangular;
        this.mConstructed = false;
    }

    // To validate construction, create an object as shown:
    //     const surface = new DerivedClassSurface(parameters);
    //     if (!surface.isConstructed()) { /* constructor failed */ }
    isConstructed(): boolean {
        return this.mConstructed;
    }

    // Member access.
    getDimension(): number {
        return this.mDimension;
    }

    getUMin(): number {
        return this.mUMin;
    }

    getUMax(): number {
        return this.mUMax;
    }

    getVMin(): number {
        return this.mVMin;
    }

    getVMax(): number {
        return this.mVMax;
    }

    isRectangular(): boolean {
        return this.mRectangular;
    }

    // Allocate storage for a jet of this surface: SUP_ORDER zero vectors of
    // the surface dimension. This is the port of upstream's
    // 'std::array<Vector<N, Real>, SUP_ORDER> jet{}'; TS callers of
    // evaluate() need it because the array cannot be sized from a template
    // parameter.
    createJet(): Vector[] {
        const jet = new Array<Vector>(ParametricSurface.SUP_ORDER);
        for (let i = 0; i < ParametricSurface.SUP_ORDER; ++i) {
            jet[i] = new Vector(this.mDimension);
        }
        return jet;
    }

    // Evaluation of the surface. The function supports derivative
    // calculation through order 2; that is, order <= 2 is required. If you
    // want only the position, pass in order of 0. If you want the position
    // and first-order derivatives, pass in order of 1, and so on. The output
    // array 'jet' must have enough storage to support the maximum order. The
    // values are ordered as: position X; first-order derivatives dX/du,
    // dX/dv; second-order derivatives d2X/du2, d2X/dudv, d2X/dv2.
    abstract evaluate(u: number, v: number, order: number, jet: Vector[]): void;

    // Differential geometric quantities. Upstream's jet is an array of
    // value objects, so an implementation of Evaluate that writes a stored
    // vector into a jet slot writes a copy. TS arrays hold references, so
    // these accessors clone the jet entry before returning (and before
    // normalizing) to keep the C++ value semantics.
    getPosition(u: number, v: number): Vector {
        const jet = this.createJet();
        this.evaluate(u, v, 0, jet);
        return jet[0].clone();
    }

    getUTangent(u: number, v: number): Vector {
        const jet = this.createJet();
        this.evaluate(u, v, 1, jet);
        const tangent = jet[1].clone();
        normalize(tangent);
        return tangent;
    }

    getVTangent(u: number, v: number): Vector {
        const jet = this.createJet();
        this.evaluate(u, v, 1, jet);
        const tangent = jet[2].clone();
        normalize(tangent);
        return tangent;
    }
}
