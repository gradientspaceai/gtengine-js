// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteTricubic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite tricubic polynomial is
//   H(x,y,z) = sum_{i=0}^3 sum_{j=0}^3 sum_{k=0}^3
//              c[i][j][k] * P(i,x) * P(j,y) * P(k,z)
// where P(i,t) = (1-t)^{3-i} * t^i. The domain is (x,y,z) in [0,1]^3.
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteTricubic::Sample is exported as
// HermiteTricubicSample (global export uniqueness). Fields F, Fx, ..., Fxyz
// are ported camelCase as f, fx, ..., fxyz.
// 'operator()(xOrder, yOrder, zOrder, x, y, z)' is ported as
// 'evaluate(xOrder, yOrder, zOrder, x, y, z)'. The basis polynomial
// evaluation reuses HermiteCubic.p, exactly as upstream reuses
// HermiteCubic<T>::P.

import { HermiteCubic } from './HermiteCubic.js';

// Sample of a function at a lattice voxel: the function value f and the
// derivatives fx, fy, fz, fxy, fxz, fyz, fxyz. The default sample is
// all-zero.
export class HermiteTricubicSample {
    f: number;
    fx: number;
    fy: number;
    fz: number;
    fxy: number;
    fxz: number;
    fyz: number;
    fxyz: number;

    constructor(f: number = 0, fx: number = 0, fy: number = 0, fz: number = 0,
        fxy: number = 0, fxz: number = 0, fyz: number = 0, fxyz: number = 0) {
        this.f = f;
        this.fx = fx;
        this.fy = fy;
        this.fz = fz;
        this.fxy = fxy;
        this.fxz = fxz;
        this.fyz = fyz;
        this.fxyz = fxyz;
    }
}

type HermiteTricubicPair = readonly [HermiteTricubicSample, HermiteTricubicSample];
type HermiteTricubicBlocks = readonly [
    readonly [HermiteTricubicPair, HermiteTricubicPair],
    readonly [HermiteTricubicPair, HermiteTricubicPair]
];

export class HermiteTricubic {
    // Set the coefficients manually as desired (public, as upstream). For Hermite tricubic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C1-continuous.
    c: number[][][];

    // The default polynomial is identically zero. When 'blocks' is provided,
    // generate the 4x4x4 coefficients c[][][] for a cell of the lattice with
    // voxels at (x,y,z), (x+1,y,z), (x,y+1,z), (x+1,y+1,z), (x,y,z+1),
    // (x+1,y,z+1), (x,y+1,z+1) and (x+1,y+1,z+1); blocks[b0][b1][b2] is the
    // sample at (x+b0, y+b1, z+b2). The caller is responsible for tracking
    // the voxel (x,y,z) that is associated with the coefficients.
    constructor(blocks?: HermiteTricubicBlocks) {
        this.c = [];
        for (let i = 0; i < 4; ++i) {
            const ci: number[][] = [];
            for (let j = 0; j < 4; ++j) {
                ci.push([0, 0, 0, 0]);
            }
            this.c.push(ci);
        }
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified orders. The returned value
    // is zero if xOrder >= 4 or yOrder >= 4 or zOrder >= 4. Otherwise, some
    // examples are the following where 'hermite' is of type HermiteTricubic:
    //   H(x, y, z) = hermite.evaluate(0, 0, 0, x, y, z)
    //   Hx(x, y, z) = hermite.evaluate(1, 0, 0, x, y, z)
    //   Hy(x, y, z) = hermite.evaluate(0, 1, 0, x, y, z)
    //   Hz(x, y, z) = hermite.evaluate(0, 0, 1, x, y, z)
    //   Hxx(x, y, z) = hermite.evaluate(2, 0, 0, x, y, z)
    //   Hxy(x, y, z) = hermite.evaluate(1, 1, 0, x, y, z)
    //   Hxz(x, y, z) = hermite.evaluate(1, 0, 1, x, y, z)
    //   Hyy(x, y, z) = hermite.evaluate(0, 2, 0, x, y, z)
    //   Hyz(x, y, z) = hermite.evaluate(0, 1, 1, x, y, z)
    //   Hzz(x, y, z) = hermite.evaluate(0, 0, 2, x, y, z)
    evaluate(xOrder: number, yOrder: number, zOrder: number,
        x: number, y: number, z: number): number {
        if (xOrder <= 3 && yOrder <= 3 && zOrder <= 3) {
            let result = 0;
            for (let i = 0; i < 4; ++i) {
                const xValue = HermiteCubic.p(i, xOrder, x);
                for (let j = 0; j < 4; ++j) {
                    const yValue = HermiteCubic.p(j, yOrder, y);
                    const xyValue = xValue * yValue;
                    for (let k = 0; k < 4; ++k) {
                        const zValue = HermiteCubic.p(k, zOrder, z);
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

    generate(blocks: HermiteTricubicBlocks): void {
        for (let b0 = 0; b0 <= 1; ++b0) {
            const z0 = 3 * b0 + 0;
            const p0 = 1 * b0 + 1;
            const s0 = 1 - 2 * b0;

            for (let b1 = 0; b1 <= 1; ++b1) {
                const z1 = 3 * b1 + 0;
                const p1 = 1 * b1 + 1;
                const s1 = 1 - 2 * b1;
                const s0s1 = s0 * s1;

                for (let b2 = 0; b2 <= 1; ++b2) {
                    const z2 = 3 * b2 + 0;
                    const p2 = 1 * b2 + 1;
                    const s2 = 1 - 2 * b2;
                    const s0s2 = s0 * s2;
                    const s1s2 = s1 * s2;
                    const s0s1s2 = s0 * s1s2;

                    const b = blocks[b0][b1][b2];
                    const input = new HermiteTricubicSample(b.f, s0 * b.fx,
                        s1 * b.fy, s2 * b.fz, s0s1 * b.fxy, s0s2 * b.fxz,
                        s1s2 * b.fyz, s0s1s2 * b.fxyz);

                    const [v000, v100, v010, v001, v110, v101, v011, v111] =
                        HermiteTricubic.generateSingle(input);
                    this.c[z0][z1][z2] = v000;
                    this.c[p0][z1][z2] = v100;
                    this.c[z0][p1][z2] = v010;
                    this.c[z0][z1][p2] = v001;
                    this.c[p0][p1][z2] = v110;
                    this.c[p0][z1][p2] = v101;
                    this.c[z0][p1][p2] = v011;
                    this.c[p0][p1][p2] = v111;
                }
            }
        }
    }

    private static generateSingle(input: HermiteTricubicSample):
        [number, number, number, number, number, number, number, number] {
        const f = input.f;
        const fx = input.fx;
        const fy = input.fy;
        const fz = input.fz;
        const fxy = input.fxy;
        const fxz = input.fxz;
        const fyz = input.fyz;
        const fxyz = input.fxyz;

        const v000 = f;
        const v100 = 3 * f + fx;
        const v010 = 3 * f + fy;
        const v001 = 3 * f + fz;
        const v110 = 9 * f + 3 * (fx + fy) + fxy;
        const v101 = 9 * f + 3 * (fx + fz) + fxz;
        const v011 = 9 * f + 3 * (fy + fz) + fyz;
        const v111 = 27 * f + 9 * (fx + fy + fz) + 3 * (fxy + fxz + fyz) + fxyz;
        return [v000, v100, v010, v001, v110, v101, v011, v111];
    }
}
