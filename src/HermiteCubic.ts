// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteCubic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite cubic polynomial is
//   H(x) = sum_{i=0}^3 c[i] * P(i,x)
// where P(i,x) = (1-x)^{3-i} * x^i. The domain is x in [0,1].
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteCubic::Sample is exported as
// HermiteCubicSample (global export uniqueness; the Bi/Tri variants have
// their own Sample structs). Fields F, Fx are ported camelCase as f, fx.
// 'operator()(xOrder, x)' is ported as 'evaluate(xOrder, x)'.

// Sample of a function at a lattice pixel: the function value f and the
// first derivative fx. The default sample is all-zero.
export class HermiteCubicSample {
    f: number;
    fx: number;

    constructor(f: number = 0, fx: number = 0) {
        this.f = f;
        this.fx = fx;
    }
}

// The basis polynomials P(i,t) = (1-t)^{3-i} * t^i and their derivatives.
// pIDJ is the J-th derivative of P(I,t).
function p0d0(t: number): number {
    return +p3d0(1 - t);
}

function p0d1(t: number): number {
    return -p3d1(1 - t);
}

function p0d2(t: number): number {
    return +p3d2(1 - t);
}

function p0d3(t: number): number {
    return -p3d3(1 - t);
}

function p1d0(t: number): number {
    return +p2d0(1 - t);
}

function p1d1(t: number): number {
    return -p2d1(1 - t);
}

function p1d2(t: number): number {
    return +p2d2(1 - t);
}

function p1d3(t: number): number {
    return -p2d3(1 - t);
}

function p2d0(t: number): number {
    return (1 - t) * t * t;
}

function p2d1(t: number): number {
    return t * (2 - 3 * t);
}

function p2d2(t: number): number {
    return 2 * (1 - 3 * t);
}

function p2d3(_t: number): number {
    return -6;
}

function p3d0(t: number): number {
    return t * t * t;
}

function p3d1(t: number): number {
    return 3 * t * t;
}

function p3d2(t: number): number {
    return 6 * t;
}

function p3d3(_t: number): number {
    return 6;
}

const table: ReadonlyArray<ReadonlyArray<(t: number) => number>> = [
    [p0d0, p0d1, p0d2, p0d3],
    [p1d0, p1d1, p1d2, p1d3],
    [p2d0, p2d1, p2d2, p2d3],
    [p3d0, p3d1, p3d2, p3d3]
];

export class HermiteCubic {
    // Set the coefficients manually as desired. For Hermite cubic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C1-continuous.
    private c: number[];

    // The default polynomial is identically zero. When 'blocks' is
    // provided, generate the 4x1 coefficients c[] for a cell of the
    // lattice with pixels at (x) and (x+1). The caller is responsible for
    // tracking the pixel (x) that is associated with the coefficients.
    constructor(blocks?: readonly [HermiteCubicSample, HermiteCubicSample]) {
        this.c = [0, 0, 0, 0];
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified order. The returned
    // value is zero if xOrder >= 4. Otherwise, some examples are the
    // following where 'hermite' is of type HermiteCubic:
    //   H(x) = hermite.evaluate(0, x)
    //   Hx(x) = hermite.evaluate(1, x)
    //   Hxx(x) = hermite.evaluate(2, x)
    //   Hxxx(x) = hermite.evaluate(3, x)
    evaluate(xOrder: number, x: number): number {
        if (xOrder <= 3) {
            let result = 0;
            for (let i = 0; i < 4; ++i) {
                result += this.c[i] * HermiteCubic.p(i, xOrder, x);
            }
            return result;
        }
        else {
            return 0;
        }
    }

    generate(blocks: readonly [HermiteCubicSample, HermiteCubicSample]): void {
        for (let b0 = 0; b0 <= 1; ++b0) {
            const z0 = 3 * b0 + 0;
            const p0 = 1 * b0 + 1;
            const s0 = 1 - 2 * b0;
            const b = blocks[b0];
            const input = new HermiteCubicSample(b.f, s0 * b.fx);
            const [v0, v1] = HermiteCubic.generateSingle(input);
            this.c[z0] = v0;
            this.c[p0] = v1;
        }
    }

    private static generateSingle(input: HermiteCubicSample): [number, number] {
        const f = input.f;
        const fx = input.fx;
        const v0 = f;
        const v1 = 3 * f + fx;
        return [v0, v1];
    }

    // For internal use in HermiteBicubic and HermiteTricubic. The 'select'
    // parameter is the i argument for the polynomial and the 'order'
    // parameter is the order of the derivative.
    static p(select: number, order: number, t: number): number {
        return order <= 3 ? table[select][order](t) : 0;
    }
}
