// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteBicubic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite bicubic polynomial is
//   H(x,y) = sum_{i=0}^3 sum_{j=0}^3 c[i][j] * P(i,x) * P(j,y)
// where P(i,t) = (1-t)^{3-i} * t^i. The domain is (x,y) in [0,1]^2.
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteBicubic::Sample is exported as
// HermiteBicubicSample (global export uniqueness). Fields F, Fx, Fy, Fxy are
// ported camelCase as f, fx, fy, fxy. 'operator()(xOrder, yOrder, x, y)' is
// ported as 'evaluate(xOrder, yOrder, x, y)'. The basis polynomial evaluation
// reuses HermiteCubic.p, exactly as upstream reuses HermiteCubic<T>::P.

import { HermiteCubic } from './HermiteCubic.js';

// Sample of a function at a lattice pixel: the function value f, the first
// derivatives fx, fy and the mixed second derivative fxy. The default sample
// is all-zero.
export class HermiteBicubicSample {
    f: number;
    fx: number;
    fy: number;
    fxy: number;

    constructor(f: number = 0, fx: number = 0, fy: number = 0, fxy: number = 0) {
        this.f = f;
        this.fx = fx;
        this.fy = fy;
        this.fxy = fxy;
    }
}

type HermiteBicubicBlocks = readonly [
    readonly [HermiteBicubicSample, HermiteBicubicSample],
    readonly [HermiteBicubicSample, HermiteBicubicSample]
];

export class HermiteBicubic {
    // Set the coefficients manually as desired (public, as upstream). For Hermite bicubic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C1-continuous.
    c: number[][];

    // The default polynomial is identically zero. When 'blocks' is provided,
    // generate the 4x4 coefficients c[][] for a cell of the lattice with
    // pixels at (x,y), (x+1,y), (x,y+1) and (x+1,y+1); blocks[b0][b1] is the
    // sample at (x+b0, y+b1). The caller is responsible for tracking the
    // pixel (x,y) that is associated with the coefficients.
    constructor(blocks?: HermiteBicubicBlocks) {
        this.c = [];
        for (let i = 0; i < 4; ++i) {
            this.c.push([0, 0, 0, 0]);
        }
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified orders. The returned value
    // is zero if xOrder >= 4 or yOrder >= 4. Otherwise, some examples are
    // the following where 'hermite' is of type HermiteBicubic:
    //   H(x, y) = hermite.evaluate(0, 0, x, y)
    //   Hx(x, y) = hermite.evaluate(1, 0, x, y)
    //   Hy(x, y) = hermite.evaluate(0, 1, x, y)
    //   Hxx(x, y) = hermite.evaluate(2, 0, x, y)
    //   Hxy(x, y) = hermite.evaluate(1, 1, x, y)
    //   Hyy(x, y) = hermite.evaluate(0, 2, x, y)
    evaluate(xOrder: number, yOrder: number, x: number, y: number): number {
        if (xOrder <= 3 && yOrder <= 3) {
            let result = 0;
            for (let i = 0; i < 4; ++i) {
                const xValue = HermiteCubic.p(i, xOrder, x);
                for (let j = 0; j < 4; ++j) {
                    const yValue = HermiteCubic.p(j, yOrder, y);
                    result += this.c[i][j] * xValue * yValue;
                }
            }
            return result;
        }
        else {
            return 0;
        }
    }

    generate(blocks: HermiteBicubicBlocks): void {
        for (let b0 = 0; b0 <= 1; ++b0) {
            const z0 = 3 * b0 + 0;
            const p0 = 1 * b0 + 1;
            const s0 = 1 - 2 * b0;

            for (let b1 = 0; b1 <= 1; ++b1) {
                const z1 = 3 * b1 + 0;
                const p1 = 1 * b1 + 1;
                const s1 = 1 - 2 * b1;
                const s0s1 = s0 * s1;

                const b = blocks[b0][b1];
                const input = new HermiteBicubicSample(
                    b.f, s0 * b.fx, s1 * b.fy, s0s1 * b.fxy);

                const [v00, v10, v01, v11] = HermiteBicubic.generateSingle(input);
                this.c[z0][z1] = v00;
                this.c[p0][z1] = v10;
                this.c[z0][p1] = v01;
                this.c[p0][p1] = v11;
            }
        }
    }

    private static generateSingle(
        input: HermiteBicubicSample): [number, number, number, number] {
        const f = input.f;
        const fx = input.fx;
        const fy = input.fy;
        const fxy = input.fxy;

        const v00 = f;
        const v10 = 3 * f + fx;
        const v01 = 3 * f + fy;
        const v11 = 9 * f + 3 * (fx + fy) + fxy;
        return [v00, v10, v01, v11];
    }
}
