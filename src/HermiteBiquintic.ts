// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteBiquintic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite biquintic polynomial is
//   G(x,y) = sum_{i=0}^5 sum_{j=0}^5 c[i][j] * P(i,x) * P(j,y)
// where P(i,t) = (1-t)^{5-i} * t^i. The domain is (x,y) in [0,1]^2.
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteBiquintic::Sample is exported as
// HermiteBiquinticSample (global export uniqueness). Fields F, Fx, ...,
// Fxxyy are ported camelCase as f, fx, ..., fxxyy. The overloaded private
// 'Generate(input, v00, ...)' is ported as 'generateSingle' (output
// references become a returned tuple). 'operator()(xOrder, yOrder, x, y)'
// is ported as 'evaluate(xOrder, yOrder, x, y)'. The basis polynomial
// evaluation reuses HermiteQuintic.p, exactly as upstream reuses
// HermiteQuintic<T>::P.

import { HermiteQuintic } from './HermiteQuintic';

// Sample of a function at a lattice pixel: the function value f and the
// derivatives fx, fy, fxx, fxy, fyy, fxxy, fxyy, fxxyy. The default sample
// is all-zero.
export class HermiteBiquinticSample {
    f: number;
    fx: number;
    fy: number;
    fxx: number;
    fxy: number;
    fyy: number;
    fxxy: number;
    fxyy: number;
    fxxyy: number;

    constructor(f: number = 0, fx: number = 0, fy: number = 0,
        fxx: number = 0, fxy: number = 0, fyy: number = 0,
        fxxy: number = 0, fxyy: number = 0, fxxyy: number = 0) {
        this.f = f;
        this.fx = fx;
        this.fy = fy;
        this.fxx = fxx;
        this.fxy = fxy;
        this.fyy = fyy;
        this.fxxy = fxxy;
        this.fxyy = fxyy;
        this.fxxyy = fxxyy;
    }
}

type HermiteBiquinticBlocks = readonly [
    readonly [HermiteBiquinticSample, HermiteBiquinticSample],
    readonly [HermiteBiquinticSample, HermiteBiquinticSample]
];

export class HermiteBiquintic {
    // Set the coefficients manually as desired. For Hermite biquintic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C2-continuous.
    private c: number[][];

    // The default polynomial is identically zero. When 'blocks' is provided,
    // generate the 6x6 coefficients c[][] for a cell of the lattice with
    // pixels at (x,y), (x+1,y), (x,y+1) and (x+1,y+1); blocks[b0][b1] is the
    // sample at (x+b0, y+b1). The caller is responsible for tracking the
    // pixel (x,y) that is associated with the coefficients.
    constructor(blocks?: HermiteBiquinticBlocks) {
        this.c = [];
        for (let i = 0; i < 6; ++i) {
            this.c.push([0, 0, 0, 0, 0, 0]);
        }
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified orders. The returned value
    // is zero if xOrder >= 6 or yOrder >= 6. Otherwise, some examples are
    // the following where 'hermite' is of type HermiteBiquintic:
    //   G(x, y) = hermite.evaluate(0, 0, x, y)
    //   Gx(x, y) = hermite.evaluate(1, 0, x, y)
    //   Gy(x, y) = hermite.evaluate(0, 1, x, y)
    //   Gxx(x, y) = hermite.evaluate(2, 0, x, y)
    //   Gxy(x, y) = hermite.evaluate(1, 1, x, y)
    //   Gyy(x, y) = hermite.evaluate(0, 2, x, y)
    evaluate(xOrder: number, yOrder: number, x: number, y: number): number {
        if (xOrder <= 5 && yOrder <= 5) {
            let result = 0;
            for (let i = 0; i < 6; ++i) {
                const xValue = HermiteQuintic.p(i, xOrder, x);
                for (let j = 0; j < 6; ++j) {
                    const yValue = HermiteQuintic.p(j, yOrder, y);
                    result += this.c[i][j] * xValue * yValue;
                }
            }
            return result;
        }
        else {
            return 0;
        }
    }

    generate(blocks: HermiteBiquinticBlocks): void {
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

                const b = blocks[b0][b1];
                const input = new HermiteBiquinticSample(b.f, s0 * b.fx,
                    s1 * b.fy, b.fxx, s0s1 * b.fxy, b.fyy, s1 * b.fxxy,
                    s0 * b.fxyy, b.fxxyy);

                const [v00, v10, v01, v20, v11, v02, v21, v12, v22] =
                    HermiteBiquintic.generateSingle(input);
                this.c[z0][z1] = v00;
                this.c[p0][z1] = v10;
                this.c[z0][p1] = v01;
                this.c[q0][z1] = v20;
                this.c[p0][p1] = v11;
                this.c[z0][q1] = v02;
                this.c[q0][p1] = v21;
                this.c[p0][q1] = v12;
                this.c[q0][q1] = v22;
            }
        }
    }

    private static generateSingle(input: HermiteBiquinticSample):
        [number, number, number, number, number, number, number, number, number] {
        const f = input.f;
        const fx = input.fx;
        const fy = input.fy;
        const fxx = input.fxx;
        const fxy = input.fxy;
        const fyy = input.fyy;
        const fxxy = input.fxxy;
        const fxyy = input.fxyy;
        const fxxyy = input.fxxyy;

        const v00 = f;
        const v10 = 5 * f + fx;
        const v01 = 5 * f + fy;
        const v20 = 10 * f + 4 * fx + 0.5 * fxx;
        const v11 = 25 * f + 5 * (fx + fy) + fxy;
        const v02 = 10 * f + 4 * fy + 0.5 * fyy;
        const v21 = 50 * f + 20 * fx + 10 * fy + 2.5 * fxx + 4 * fxy + 0.5 * fxxy;
        const v12 = 50 * f + 10 * fx + 20 * fy + 4 * fxy + 2.5 * fyy + 0.5 * fxyy;
        const v22 = 100 * f + 40 * (fx + fy) + 5 * (fxx + fyy) + 16 * fxy
            + 2 * (fxxy + fxyy) + 0.25 * fxxyy;
        return [v00, v10, v01, v20, v11, v02, v21, v12, v22];
    }
}
