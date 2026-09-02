// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) EllipsoidGeodesic.h
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
// Port notes: the Christoffel symbols of the first kind use the convention
// Gamma_{k,ij} = Dot(P_ij, P_k) of RiemannianGeodesic.ts (see issue #295),
// which is exactly what upstream computes here.

import { GVector } from './GVector';
import { RiemannianGeodesic } from './RiemannianGeodesic';
import { Vector, dot } from './Vector';

export class EllipsoidGeodesic extends RiemannianGeodesic {
    // The ellipsoid axis half-lengths.
    private mXExtent: number;
    private mYExtent: number;
    private mZExtent: number;

    // We are guaranteed that RiemannianGeodesic calls computeMetric before
    // computeChristoffel1. Thus, we can compute the surface first- and
    // second-order derivatives in computeMetric and cache the results for use
    // in computeChristoffel1.
    private mCos0: number;
    private mSin0: number;
    private mCos1: number;
    private mSin1: number;
    private mDer0: Vector;
    private mDer1: Vector;

    // The ellipsoid is (x/a)^2 + (y/b)^2 + (z/c)^2 = 1, where xExtent is 'a',
    // yExtent is 'b' and zExtent is 'c'. The surface is represented
    // parametrically by angles u and v, say
    //   P(u,v) = (x(u,v),y(u,v),z(u,v)),
    //   P(u,v) = (a*cos(u)*sin(v), b*sin(u)*sin(v), c*cos(v))
    // with 0 <= u < 2*pi and 0 <= v <= pi. The first-order derivatives are
    //   dP/du = (-a*sin(u)*sin(v), b*cos(u)*sin(v), 0)
    //   dP/dv = (a*cos(u)*cos(v), b*sin(u)*cos(v), -c*sin(v))
    // The metric tensor elements are
    //   g_{00} = Dot(dP/du,dP/du)
    //   g_{01} = Dot(dP/du,dP/dv)
    //   g_{10} = g_{01}
    //   g_{11} = Dot(dP/dv,dP/dv)
    constructor(xExtent: number, yExtent: number, zExtent: number) {
        super(2);
        this.mXExtent = xExtent;
        this.mYExtent = yExtent;
        this.mZExtent = zExtent;
        this.mCos0 = 0;
        this.mSin0 = 0;
        this.mCos1 = 0;
        this.mSin1 = 0;
        this.mDer0 = new Vector(3);
        this.mDer1 = new Vector(3);
    }

    // The surface point P(u,v) for point = (u,v).
    computePosition(point: GVector): Vector {
        const cos0 = Math.cos(point.values[0]);
        const sin0 = Math.sin(point.values[0]);
        const cos1 = Math.cos(point.values[1]);
        const sin1 = Math.sin(point.values[1]);

        return Vector.fromArray([
            this.mXExtent * cos0 * sin1,
            this.mYExtent * sin0 * sin1,
            this.mZExtent * cos1
        ]);
    }

    // To compute the geodesic path connecting two parameter points (u0,v0)
    // and (u1,v1):
    //
    //   const eg = new EllipsoidGeodesic(a, b, c);
    //   const param0 = GVector.fromArray([u0, v0]);
    //   const param1 = GVector.fromArray([u1, v1]);
    //   const { quantity, path } = eg.computeGeodesic(param0, param1);

    protected override computeMetric(point: GVector): void {
        this.mCos0 = Math.cos(point.values[0]);
        this.mSin0 = Math.sin(point.values[0]);
        this.mCos1 = Math.cos(point.values[1]);
        this.mSin1 = Math.sin(point.values[1]);

        this.mDer0 = Vector.fromArray([
            -this.mXExtent * this.mSin0 * this.mSin1,
            this.mYExtent * this.mCos0 * this.mSin1,
            0
        ]);
        this.mDer1 = Vector.fromArray([
            this.mXExtent * this.mCos0 * this.mCos1,
            this.mYExtent * this.mSin0 * this.mCos1,
            -this.mZExtent * this.mSin1
        ]);

        this.mMetric.set(0, 0, dot(this.mDer0, this.mDer0));
        this.mMetric.set(0, 1, dot(this.mDer0, this.mDer1));
        this.mMetric.set(1, 0, this.mMetric.get(0, 1));
        this.mMetric.set(1, 1, dot(this.mDer1, this.mDer1));
    }

    protected override computeChristoffel1(_point: GVector): void {
        const der00 = Vector.fromArray([
            -this.mXExtent * this.mCos0 * this.mSin1,
            -this.mYExtent * this.mSin0 * this.mSin1,
            0
        ]);

        const der01 = Vector.fromArray([
            -this.mXExtent * this.mSin0 * this.mCos1,
            this.mYExtent * this.mCos0 * this.mCos1,
            0
        ]);

        const der11 = Vector.fromArray([
            -this.mXExtent * this.mCos0 * this.mSin1,
            -this.mYExtent * this.mSin0 * this.mSin1,
            -this.mZExtent * this.mCos1
        ]);

        this.mChristoffel1[0].set(0, 0, dot(der00, this.mDer0));
        this.mChristoffel1[0].set(0, 1, dot(der01, this.mDer0));
        this.mChristoffel1[0].set(1, 0, this.mChristoffel1[0].get(0, 1));
        this.mChristoffel1[0].set(1, 1, dot(der11, this.mDer0));

        this.mChristoffel1[1].set(0, 0, dot(der00, this.mDer1));
        this.mChristoffel1[1].set(0, 1, dot(der01, this.mDer1));
        this.mChristoffel1[1].set(1, 0, this.mChristoffel1[1].get(0, 1));
        this.mChristoffel1[1].set(1, 1, dot(der11, this.mDer1));
    }
}
