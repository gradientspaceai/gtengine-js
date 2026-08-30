// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) EulerAngles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Euler angle data structure for representing rotations. See the
// document
//   https://www.geometrictools.com/Documentation/EulerAngles.pdf

// Factorization into Euler angles is not necessarily unique. Let the
// integer indices for the axes be (N0,N1,N2), which must be in the set
//   {(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0),
//    (0,1,0),(0,2,0),(1,0,1),(1,2,1),(2,0,2),(2,1,2)}
// Let the corresponding angles be (angleN0,angleN1,angleN2). If the
// result is NOT_UNIQUE_SUM, then the multiple solutions occur because
// angleN2+angleN0 is constant. If the result is NOT_UNIQUE_DIF, then
// the multiple solutions occur because angleN2-angleN0 is constant.
// In either type of nonuniqueness, the function returns angleN0=0.
export enum EulerResult {
    // The solution is invalid (incorrect axis indices).
    INVALID,

    // The solution is unique.
    UNIQUE,

    // The solution is not unique. A sum of angles is constant.
    NOT_UNIQUE_SUM,

    // The solution is not unique. A difference of angles is constant.
    NOT_UNIQUE_DIF
}

export class EulerAngles {
    // The axis indices (each 0, 1 or 2) of the factorization.
    axis: [number, number, number];

    // The angles (in radians) about the corresponding axes.
    angle: [number, number, number];

    // This member is set during conversions from rotation matrices,
    // quaternions, or axis-angles.
    result: EulerResult;

    // The default constructor produces the INVALID state with zeroed axes
    // and angles; the 6-argument constructor produces the UNIQUE state, as
    // upstream.
    constructor();
    constructor(i0: number, i1: number, i2: number,
        a0: number, a1: number, a2: number);
    constructor(i0?: number, i1?: number, i2?: number,
        a0?: number, a1?: number, a2?: number) {
        if (i0 === undefined) {
            this.axis = [0, 0, 0];
            this.angle = [0, 0, 0];
            this.result = EulerResult.INVALID;
        } else {
            this.axis = [i0, i1 as number, i2 as number];
            this.angle = [a0 as number, a1 as number, a2 as number];
            this.result = EulerResult.UNIQUE;
        }
    }
}
