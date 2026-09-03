// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteTriquintic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite triquintic polynomial is
//   G(x,y,z) = sum_{i=0}^5 sum_{j=0}^5 sum_{k=0}^5
//              c[i][j][k] * P(i,x) * P(j,y) * P(k,z)
// where P(i,t) = (1-t)^{5-i} * t^i. The domain is (x,y,z) in [0,1]^3.
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteTriquintic::Sample is exported
// as HermiteTriquinticSample (global export uniqueness). Fields F, Fx, ...,
// Fxxyyzz are ported camelCase as f, fx, ..., fxxyyzz.
// 'operator()(xOrder, yOrder, zOrder, x, y, z)' is ported as
// 'evaluate(xOrder, yOrder, zOrder, x, y, z)'. The private GenerateSingle
// with 27 output references returns an object with fields v000..v222. The
// basis polynomial evaluation reuses HermiteQuintic.p, exactly as upstream
// reuses HermiteQuintic<T>::P.

import { HermiteQuintic } from './HermiteQuintic.js';

// Sample of a function at a lattice voxel: the function value f and the
// derivatives through order (2,2,2). The default sample is all-zero.
export class HermiteTriquinticSample {
    f: number;
    fx: number;
    fy: number;
    fz: number;
    fxx: number;
    fxy: number;
    fxz: number;
    fyy: number;
    fyz: number;
    fzz: number;
    fxxy: number;
    fxxz: number;
    fxyy: number;
    fxyz: number;
    fxzz: number;
    fyyz: number;
    fyzz: number;
    fxxyy: number;
    fxxyz: number;
    fxxzz: number;
    fxyyz: number;
    fxyzz: number;
    fyyzz: number;
    fxxyyz: number;
    fxxyzz: number;
    fxyyzz: number;
    fxxyyzz: number;

    constructor(
        f: number = 0,
        fx: number = 0, fy: number = 0, fz: number = 0,
        fxx: number = 0, fxy: number = 0, fxz: number = 0,
        fyy: number = 0, fyz: number = 0, fzz: number = 0,
        fxxy: number = 0, fxxz: number = 0, fxyy: number = 0,
        fxyz: number = 0, fxzz: number = 0, fyyz: number = 0, fyzz: number = 0,
        fxxyy: number = 0, fxxyz: number = 0, fxxzz: number = 0,
        fxyyz: number = 0, fxyzz: number = 0, fyyzz: number = 0,
        fxxyyz: number = 0, fxxyzz: number = 0, fxyyzz: number = 0,
        fxxyyzz: number = 0) {
        this.f = f;
        this.fx = fx;
        this.fy = fy;
        this.fz = fz;
        this.fxx = fxx;
        this.fxy = fxy;
        this.fxz = fxz;
        this.fyy = fyy;
        this.fyz = fyz;
        this.fzz = fzz;
        this.fxxy = fxxy;
        this.fxxz = fxxz;
        this.fxyy = fxyy;
        this.fxyz = fxyz;
        this.fxzz = fxzz;
        this.fyyz = fyyz;
        this.fyzz = fyzz;
        this.fxxyy = fxxyy;
        this.fxxyz = fxxyz;
        this.fxxzz = fxxzz;
        this.fxyyz = fxyyz;
        this.fxyzz = fxyzz;
        this.fyyzz = fyyzz;
        this.fxxyyz = fxxyyz;
        this.fxxyzz = fxxyzz;
        this.fxyyzz = fxyyzz;
        this.fxxyyzz = fxxyyzz;
    }
}

type HermiteTriquinticPair = readonly [HermiteTriquinticSample, HermiteTriquinticSample];
type HermiteTriquinticBlocks = readonly [
    readonly [HermiteTriquinticPair, HermiteTriquinticPair],
    readonly [HermiteTriquinticPair, HermiteTriquinticPair]
];

interface HermiteTriquinticCorner {
    v000: number; v100: number; v010: number; v001: number;
    v200: number; v110: number; v101: number; v020: number;
    v011: number; v002: number; v210: number; v201: number;
    v120: number; v111: number; v102: number; v021: number;
    v012: number; v220: number; v211: number; v202: number;
    v121: number; v112: number; v022: number; v221: number;
    v212: number; v122: number; v222: number;
}

export class HermiteTriquintic {
    // Set the coefficients manually as desired. For Hermite triquintic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C2-continuous.
    private c: number[][][];

    // The default polynomial is identically zero. When 'blocks' is provided,
    // generate the 6x6x6 coefficients c[][][] for a cell of the lattice with
    // voxels at (x,y,z), (x+1,y,z), (x,y+1,z), (x+1,y+1,z), (x,y,z+1),
    // (x+1,y,z+1), (x,y+1,z+1) and (x+1,y+1,z+1); blocks[b0][b1][b2] is the
    // sample at (x+b0, y+b1, z+b2). The caller is responsible for tracking
    // the voxel (x,y,z) that is associated with the coefficients.
    constructor(blocks?: HermiteTriquinticBlocks) {
        this.c = [];
        for (let i = 0; i < 6; ++i) {
            const ci: number[][] = [];
            for (let j = 0; j < 6; ++j) {
                ci.push([0, 0, 0, 0, 0, 0]);
            }
            this.c.push(ci);
        }
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified orders. The returned value
    // is zero if xOrder >= 6 or yOrder >= 6 or zOrder >= 6. Otherwise, some
    // examples are the following where 'hermite' is of type
    // HermiteTriquintic:
    //   G(x, y, z) = hermite.evaluate(0, 0, 0, x, y, z)
    //   Gx(x, y, z) = hermite.evaluate(1, 0, 0, x, y, z)
    //   Gy(x, y, z) = hermite.evaluate(0, 1, 0, x, y, z)
    //   Gz(x, y, z) = hermite.evaluate(0, 0, 1, x, y, z)
    //   Gxx(x, y, z) = hermite.evaluate(2, 0, 0, x, y, z)
    //   Gxy(x, y, z) = hermite.evaluate(1, 1, 0, x, y, z)
    //   Gxz(x, y, z) = hermite.evaluate(1, 0, 1, x, y, z)
    //   Gyy(x, y, z) = hermite.evaluate(0, 2, 0, x, y, z)
    //   Gyz(x, y, z) = hermite.evaluate(0, 1, 1, x, y, z)
    //   Gzz(x, y, z) = hermite.evaluate(0, 0, 2, x, y, z)
    evaluate(xOrder: number, yOrder: number, zOrder: number,
        x: number, y: number, z: number): number {
        if (xOrder <= 5 && yOrder <= 5 && zOrder <= 5) {
            let result = 0;
            for (let i = 0; i < 6; ++i) {
                const xValue = HermiteQuintic.p(i, xOrder, x);
                for (let j = 0; j < 6; ++j) {
                    const yValue = HermiteQuintic.p(j, yOrder, y);
                    const xyValue = xValue * yValue;
                    for (let k = 0; k < 6; ++k) {
                        const zValue = HermiteQuintic.p(k, zOrder, z);
                        result += this.c[i][j][k] * xyValue * zValue;
                    }
                }
            }
            return result;
        }
        else {
            return 0;
        }
    }

    generate(blocks: HermiteTriquinticBlocks): void {
        for (let b0 = 0; b0 <= 1; ++b0) {
            const z0 = 5 * b0 + 0;
            const p0 = 3 * b0 + 1;
            const q0 = 1 * b0 + 2;
            const s0 = 1 - 2 * b0;

            for (let b1 = 0; b1 <= 1; ++b1) {
                const z1 = 5 * b1 + 0;
                const p1 = 3 * b1 + 1;
                const q1 = 1 * b1 + 2;
                const s1 = 1 - 2 * b1;
                const s0s1 = s0 * s1;

                for (let b2 = 0; b2 <= 1; ++b2) {
                    const z2 = 5 * b2 + 0;
                    const p2 = 3 * b2 + 1;
                    const q2 = 1 * b2 + 2;
                    const s2 = 1 - 2 * b2;
                    const s0s2 = s0 * s2;
                    const s1s2 = s1 * s2;
                    const s0s1s2 = s0 * s1s2;

                    const b = blocks[b0][b1][b2];
                    const input = new HermiteTriquinticSample(
                        b.f,
                        s0 * b.fx,
                        s1 * b.fy,
                        s2 * b.fz,
                        b.fxx,
                        s0s1 * b.fxy,
                        s0s2 * b.fxz,
                        b.fyy,
                        s1s2 * b.fyz,
                        b.fzz,
                        s1 * b.fxxy,
                        s2 * b.fxxz,
                        s0 * b.fxyy,
                        s0s1s2 * b.fxyz,
                        s0 * b.fxzz,
                        s2 * b.fyyz,
                        s1 * b.fyzz,
                        b.fxxyy,
                        s1s2 * b.fxxyz,
                        b.fxxzz,
                        s0s2 * b.fxyyz,
                        s0s1 * b.fxyzz,
                        b.fyyzz,
                        s2 * b.fxxyyz,
                        s1 * b.fxxyzz,
                        s0 * b.fxyyzz,
                        b.fxxyyzz);

                    const v = HermiteTriquintic.generateSingle(input);
                    this.c[z0][z1][z2] = v.v000;
                    this.c[p0][z1][z2] = v.v100;
                    this.c[z0][p1][z2] = v.v010;
                    this.c[z0][z1][p2] = v.v001;
                    this.c[q0][z1][z2] = v.v200;
                    this.c[p0][p1][z2] = v.v110;
                    this.c[p0][z1][p2] = v.v101;
                    this.c[z0][q1][z2] = v.v020;
                    this.c[z0][p1][p2] = v.v011;
                    this.c[z0][z1][q2] = v.v002;
                    this.c[q0][p1][z2] = v.v210;
                    this.c[q0][z1][p2] = v.v201;
                    this.c[p0][q1][z2] = v.v120;
                    this.c[p0][p1][p2] = v.v111;
                    this.c[p0][z1][q2] = v.v102;
                    this.c[z0][q1][p2] = v.v021;
                    this.c[z0][p1][q2] = v.v012;
                    this.c[q0][q1][z2] = v.v220;
                    this.c[q0][p1][p2] = v.v211;
                    this.c[q0][z1][q2] = v.v202;
                    this.c[p0][q1][p2] = v.v121;
                    this.c[p0][p1][q2] = v.v112;
                    this.c[z0][q1][q2] = v.v022;
                    this.c[q0][q1][p2] = v.v221;
                    this.c[q0][p1][q2] = v.v212;
                    this.c[p0][q1][q2] = v.v122;
                    this.c[q0][q1][q2] = v.v222;
                }
            }
        }
    }

    private static generateSingle(
        input: HermiteTriquinticSample): HermiteTriquinticCorner {
        const f = input.f;
        const fx = input.fx;
        const fy = input.fy;
        const fz = input.fz;
        const fxx = input.fxx;
        const fxy = input.fxy;
        const fxz = input.fxz;
        const fyy = input.fyy;
        const fyz = input.fyz;
        const fzz = input.fzz;
        const fxxy = input.fxxy;
        const fxxz = input.fxxz;
        const fxyy = input.fxyy;
        const fxyz = input.fxyz;
        const fxzz = input.fxzz;
        const fyyz = input.fyyz;
        const fyzz = input.fyzz;
        const fxxyy = input.fxxyy;
        const fxxyz = input.fxxyz;
        const fxxzz = input.fxxzz;
        const fxyyz = input.fxyyz;
        const fxyzz = input.fxyzz;
        const fyyzz = input.fyyzz;
        const fxxyyz = input.fxxyyz;
        const fxxyzz = input.fxxyzz;
        const fxyyzz = input.fxyyzz;
        const fxxyyzz = input.fxxyyzz;

        const v000 = f;

        const v100 = 5 * f + fx;
        const v010 = 5 * f + fy;
        const v001 = 5 * f + fz;

        const v200 = 10 * f + 4 * fx + 0.5 * fxx;
        const v110 = 25 * f + 5 * fx + 5 * fy + fxy;
        const v101 = 25 * f + 5 * fx + 5 * fz + fxz;
        const v020 = 10 * f + 4 * fy + 0.5 * fyy;
        const v011 = 25 * f + 5 * fy + 5 * fz + fyz;
        const v002 = 10 * f + 4 * fz + 0.5 * fzz;

        const v210 = 50 * f + 20 * fx + 10 * fy + 2.5 * fxx + 4 * fxy + 0.5 * fxxy;
        const v201 = 50 * f + 20 * fx + 10 * fz + 2.5 * fxx + 4 * fxz + 0.5 * fxxz;
        const v120 = 50 * f + 10 * fx + 20 * fy + 4 * fxy + 2.5 * fyy + 0.5 * fxyy;
        const v111 = 125 * f + 25 * fx + 25 * fy + 25 * fz + 5 * fxy + 5 * fxz
            + 5 * fyz + fxyz;
        const v102 = 50 * f + 10 * fx + 20 * fz + 4 * fxz + 2.5 * fzz + 0.5 * fxzz;
        const v021 = 50 * f + 20 * fy + 10 * fz + 2.5 * fyy + 4 * fyz + 0.5 * fyyz;
        const v012 = 50 * f + 10 * fy + 20 * fz + 4 * fyz + 2.5 * fzz + 0.5 * fyzz;

        const v220 = 100 * f + 40 * fx + 40 * fy + 5 * fxx + 16 * fxy + 5 * fyy
            + 2 * fxxy + 2 * fxyy + 0.25 * fxxyy;
        const v211 = 250 * f + 100 * fx + 50 * fy + 50 * fz + 12.5 * fxx
            + 20 * fxy + 20 * fxz + 10 * fyz + 2.5 * fxxy + 2.5 * fxxz
            + 4 * fxyz + 0.5 * fxxyz;
        const v202 = 100 * f + 40 * fx + 40 * fz + 5 * fxx + 16 * fxz + 5 * fzz
            + 2 * fxxz + 2 * fxzz + 0.25 * fxxzz;
        const v121 = 250 * f + 50 * fx + 100 * fy + 50 * fz + 20 * fxy
            + 10 * fxz + 12.5 * fyy + 20 * fyz + 2.5 * fxyy + 4 * fxyz
            + 2.5 * fyyz + 0.5 * fxyyz;
        const v112 = 250 * f + 50 * fx + 50 * fy + 100 * fz + 10 * fxy
            + 20 * fxz + 20 * fyz + 12.5 * fzz + 4 * fxyz + 2.5 * fxzz
            + 2.5 * fyzz + 0.5 * fxyzz;
        const v022 = 100 * f + 40 * fy + 40 * fz + 5 * fyy + 16 * fyz + 5 * fzz
            + 2 * fyyz + 2 * fyzz + 0.25 * fyyzz;

        const v221 = 500 * f + 200 * fx + 200 * fy + 100 * fz + 25 * fxx
            + 80 * fxy + 40 * fxz + 25 * fyy + 40 * fyz + 10 * fxxy + 5 * fxxz
            + 10 * fxyy + 16 * fxyz + 5 * fyyz + 1.25 * fxxyy + 2 * fxxyz
            + 2 * fxyyz + 0.25 * fxxyyz;
        const v212 = 500 * f + 200 * fx + 100 * fy + 200 * fz + 25 * fxx
            + 40 * fxy + 80 * fxz + 40 * fyz + 25 * fzz + 5 * fxxy + 10 * fxxz
            + 16 * fxyz + 10 * fxzz + 5 * fyzz + 2 * fxxyz + 1.25 * fxxzz
            + 2 * fxyzz + 0.25 * fxxyzz;
        const v122 = 500 * f + 100 * fx + 200 * fy + 200 * fz + 40 * fxy
            + 40 * fxz + 25 * fyy + 80 * fyz + 25 * fzz + 5 * fxyy + 16 * fxyz
            + 5 * fxzz + 10 * fyyz + 10 * fyzz + 2 * fxyyz + 2 * fxyzz
            + 1.25 * fyyzz + 0.25 * fxyyzz;

        const v222 = 1000 * f + 400 * fx + 400 * fy + 400 * fz + 50 * fxx
            + 160 * fxy + 160 * fxz + 50 * fyy + 160 * fyz + 50 * fzz
            + 20 * fxxy + 20 * fxxz + 20 * fxyy + 64 * fxyz + 20 * fxzz
            + 20 * fyyz + 20 * fyzz + 2.5 * fxxyy + 8 * fxxyz + 2.5 * fxxzz
            + 8 * fxyyz + 8 * fxyzz + 2.5 * fyyzz + fxxyyz + fxxyzz + fxyyzz
            + 0.125 * fxxyyzz;

        return {
            v000, v100, v010, v001, v200, v110, v101, v020, v011, v002,
            v210, v201, v120, v111, v102, v021, v012, v220, v211, v202,
            v121, v112, v022, v221, v212, v122, v222
        };
    }
}
