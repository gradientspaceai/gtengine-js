// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSplineGeodesic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Computing geodesics on a surface is a differential geometric topic that
// involves Riemannian geometry. The algorithm for constructing geodesics that
// is implemented here uses a multiresolution approach. A description of the
// algorithm is in the document
// https://www.geometrictools.com/Documentation/RiemannianGeodesics.pdf
//
// Port notes: upstream's 'BSplineSurface<3,Real> const*' member becomes a
// reference to the port's runtime-dimension BSplineSurface, whose dimension
// must be 3. The 'std::array<Vector<3,Real>,6> mJet' cache becomes the jet
// array allocated by ParametricSurface.createJet(); note that the port's
// BSplineSurface.evaluate replaces the jet entries rather than assigning into
// them, so the cached vectors are re-read from the array in
// computeChristoffel1.

import { BSplineSurface } from './BSplineSurface.js';
import { GVector } from './GVector.js';
import { logAssert } from './Logger.js';
import { RiemannianGeodesic } from './RiemannianGeodesic.js';
import { Vector, dot } from './Vector.js';

export class BSplineGeodesic extends RiemannianGeodesic {
    private mSpline: BSplineSurface;

    // We are guaranteed that RiemannianGeodesic calls computeMetric before
    // computeChristoffel1. Thus, we can compute the B-spline first- and
    // second-order derivatives in computeMetric and cache the results for use
    // in computeChristoffel1.
    private mJet: Vector[];

    constructor(spline: BSplineSurface) {
        super(2);
        logAssert(spline.getDimension() === 3,
            'BSplineGeodesic: the spline surface must be 3-dimensional.');
        this.mSpline = spline;
        this.mJet = spline.createJet();
    }

    protected override computeMetric(point: GVector): void {
        this.mSpline.evaluate(point.values[0], point.values[1], 2, this.mJet);
        const der0 = this.mJet[1];
        const der1 = this.mJet[2];

        this.mMetric.set(0, 0, dot(der0, der0));
        this.mMetric.set(0, 1, dot(der0, der1));
        this.mMetric.set(1, 0, this.mMetric.get(0, 1));
        this.mMetric.set(1, 1, dot(der1, der1));
    }

    protected override computeChristoffel1(_point: GVector): void {
        const der0 = this.mJet[1];
        const der1 = this.mJet[2];
        const der00 = this.mJet[3];
        const der01 = this.mJet[4];
        const der11 = this.mJet[5];

        this.mChristoffel1[0].set(0, 0, dot(der00, der0));
        this.mChristoffel1[0].set(0, 1, dot(der01, der0));
        this.mChristoffel1[0].set(1, 0, this.mChristoffel1[0].get(0, 1));
        this.mChristoffel1[0].set(1, 1, dot(der11, der0));

        this.mChristoffel1[1].set(0, 0, dot(der00, der1));
        this.mChristoffel1[1].set(0, 1, dot(der01, der1));
        this.mChristoffel1[1].set(1, 0, this.mChristoffel1[1].get(0, 1));
        this.mChristoffel1[1].set(1, 1, dot(der11, der1));
    }
}
