// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IntpVectorField2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Given points (x0[i],y0[i]) which are mapped to (x1[i],y1[i]) for
// 0 <= i < N, interpolate positions (xIn,yIn) to (xOut,yOut).
//
// Port notes:
// * Upstream declares a variadic class template with two specializations,
//   the deprecated IntpVectorField2<InputType, ComputeType, RationalType>
//   and the replacement IntpVectorField2<T>. Only the replacement is ported,
//   matching the Delaunay2 and Delaunay2Mesh ports.
// * The C++ operator()(input, output&) returning bool becomes
//   evaluate(input), which returns { valid, output }; the output is
//   meaningless when valid is false.

import { Delaunay2 } from './Delaunay2.js';
import { Delaunay2Mesh } from './Delaunay2Mesh.js';
import { IntpQuadraticNonuniform2 } from './IntpQuadraticNonuniform2.js';
import { logAssert } from './Logger.js';
import { Vector } from './Vector.js';

export interface IntpVectorField2Result {
    // Valid is true if and only if the input point is in the convex hull of
    // the input domain points, in which case the interpolation is valid.
    valid: boolean;

    // The interpolated range point; meaningful only when valid is true.
    output: Vector;
}

export class IntpVectorField2 {
    private mDelaunay: Delaunay2;
    private mMesh: Delaunay2Mesh;
    private mXRange: number[];
    private mYRange: number[];
    private mXInterp: IntpQuadraticNonuniform2;
    private mYInterp: IntpQuadraticNonuniform2;

    // Construction. The domain and range arrays must have the same number of
    // 2D points.
    constructor(domain: readonly Vector[], range: readonly Vector[]) {
        logAssert(domain.length === range.length && domain.length > 0,
            'IntpVectorField2: the domain and range must have the same '
            + 'positive number of points.');

        // Repackage the output vectors into individual components. This is
        // required because of the format that the quadratic interpolator
        // expects for its input data.
        this.mXRange = [];
        this.mYRange = [];
        for (let i = 0; i < range.length; ++i) {
            logAssert(range[i].size === 2 && domain[i].size === 2,
                'IntpVectorField2: the points must be 2D.');
            this.mXRange.push(range[i].values[0]);
            this.mYRange.push(range[i].values[1]);
        }

        // Common triangulator for the interpolators.
        this.mDelaunay = new Delaunay2();
        this.mDelaunay.compute(domain);
        this.mMesh = new Delaunay2Mesh(this.mDelaunay);

        // Create the interpolator for the x-coordinate of the vector field.
        this.mXInterp = IntpQuadraticNonuniform2.fromSpatialDelta(this.mMesh,
            this.mXRange, 1);

        // Create the interpolator for the y-coordinate of the vector field,
        // but share the triangulation already created for the x-interpolator.
        this.mYInterp = IntpQuadraticNonuniform2.fromSpatialDelta(this.mMesh,
            this.mYRange, 1);
    }

    // The 'valid' field is true if and only if the input point is in the
    // convex hull of the input domain points, in which case the
    // interpolation is valid.
    evaluate(input: Vector): IntpVectorField2Result {
        const xResult = this.mXInterp.evaluate(input);
        if (!xResult.valid) {
            return { valid: false, output: Vector.zero(2) };
        }
        const yResult = this.mYInterp.evaluate(input);
        if (!yResult.valid) {
            return { valid: false, output: Vector.zero(2) };
        }
        return {
            valid: true,
            output: Vector.fromArray([xResult.F, yResult.F])
        };
    }
}
