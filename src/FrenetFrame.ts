// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FrenetFrame.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Frenet frame of a parametric curve: the moving orthonormal frame
// {tangent, normal} in 2D and {tangent, normal, binormal} in 3D, together
// with the curvature (2D and 3D) and torsion (3D).
//
// Port notes: upstream has two class templates in this header, FrenetFrame2
// and FrenetFrame3; both live in this file per the one-file-per-header rule.
// Upstream takes a 'std::shared_ptr<ParametricCurve<N, Real>>'; the port
// takes the curve object directly (the caller keeps it alive). The
// 'operator()' that writes position/tangent/normal(/binormal) through
// reference parameters becomes 'compute(t)' returning a named object literal.
// The curve's dimension is asserted at construction because the port's
// ParametricCurve carries its dimension at runtime rather than as a template
// parameter.

import { logAssert } from './Logger';
import { ParametricCurve } from './ParametricCurve';
import { Vector, dot, length, normalize, mul, sub } from './Vector';
import { perp, dotPerp } from './Vector2';
import { cross } from './Vector3';

export interface FrenetFrame2Result {
    position: Vector;
    tangent: Vector;
    normal: Vector;
}

export interface FrenetFrame3Result {
    position: Vector;
    tangent: Vector;
    normal: Vector;
    binormal: Vector;
}

export class FrenetFrame2 {
    private mCurve: ParametricCurve;

    // Construction. The curve must persist as long as the FrenetFrame2
    // object does.
    constructor(curve: ParametricCurve) {
        logAssert(curve.getDimension() === 2,
            'FrenetFrame2: the curve must be 2-dimensional.');
        this.mCurve = curve;
    }

    getCurve(): ParametricCurve {
        return this.mCurve;
    }

    // The normal is perpendicular to the tangent, rotated clockwise by
    // pi/2 radians.
    compute(t: number): FrenetFrame2Result {
        const jet = this.mCurve.createJet();
        this.mCurve.evaluate(t, 1, jet);
        const position = jet[0].clone();
        const tangent = jet[1].clone();
        normalize(tangent);
        const normal = perp(tangent);
        return { position, tangent, normal };
    }

    getCurvature(t: number): number {
        const jet = this.mCurve.createJet();
        this.mCurve.evaluate(t, 2, jet);
        const speedSqr = dot(jet[1], jet[1]);
        if (speedSqr > 0) {
            const numer = dotPerp(jet[1], jet[2]);
            const denom = Math.pow(speedSqr, 1.5);
            return numer / denom;
        }
        else {
            // Curvature is indeterminate, just return 0.
            return 0;
        }
    }
}

export class FrenetFrame3 {
    private mCurve: ParametricCurve;

    // Construction. The curve must persist as long as the FrenetFrame3
    // object does.
    constructor(curve: ParametricCurve) {
        logAssert(curve.getDimension() === 3,
            'FrenetFrame3: the curve must be 3-dimensional.');
        this.mCurve = curve;
    }

    getCurve(): ParametricCurve {
        return this.mCurve;
    }

    // The binormal is Cross(tangent, normal).
    compute(t: number): FrenetFrame3Result {
        const jet = this.mCurve.createJet();
        this.mCurve.evaluate(t, 2, jet);
        const position = jet[0].clone();
        const VDotV = dot(jet[1], jet[1]);
        const VDotA = dot(jet[1], jet[2]);
        const normal = sub(mul(VDotV, jet[2]), mul(VDotA, jet[1]));
        normalize(normal);
        const tangent = jet[1].clone();
        normalize(tangent);
        const binormal = cross(tangent, normal);
        return { position, tangent, normal, binormal };
    }

    getCurvature(t: number): number {
        const jet = this.mCurve.createJet();
        this.mCurve.evaluate(t, 2, jet);
        const speedSqr = dot(jet[1], jet[1]);
        if (speedSqr > 0) {
            const numer = length(cross(jet[1], jet[2]));
            const denom = Math.pow(speedSqr, 1.5);
            return numer / denom;
        }
        else {
            // Curvature is indeterminate, just return 0.
            return 0;
        }
    }

    getTorsion(t: number): number {
        const jet = this.mCurve.createJet();
        this.mCurve.evaluate(t, 3, jet);
        const c = cross(jet[1], jet[2]);
        const denom = dot(c, c);
        if (denom > 0) {
            const numer = dot(c, jet[3]);
            return numer / denom;
        }
        else {
            // Torsion is indeterminate, just return 0.
            return 0;
        }
    }
}
