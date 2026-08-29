// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) HermiteQuintic.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Hermite quintic polynomial is
//   H(x) = sum_{i=0}^5 c[i] * P(i,x)
// where P(i,x) = (1-x)^{5-i} * x^i. The domain is x in [0,1].
// Interpolation using these polynomials is described in
// https://www.geometrictools.com/Documentation/SmoothLatticeInterpolation.pdf
//
// Port notes: upstream nested struct HermiteQuintic::Sample is exported as
// HermiteQuinticSample (global export uniqueness; the Bi/Tri variants have
// their own Sample structs). Fields F, Fx, Fxx are ported camelCase as
// f, fx, fxx. 'operator()(xOrder, x)' is ported as 'evaluate(xOrder, x)'.

// Sample of a function at a lattice pixel: the function value f, the first
// derivative fx and the second derivative fxx. The default sample is
// all-zero.
export class HermiteQuinticSample {
    f: number;
    fx: number;
    fxx: number;

    constructor(f: number = 0, fx: number = 0, fxx: number = 0) {
        this.f = f;
        this.fx = fx;
        this.fxx = fxx;
    }
}

// The basis polynomials P(i,t) = (1-t)^{5-i} * t^i and their derivatives.
// pIDJ is the J-th derivative of P(I,t).
function p0d0(t: number): number {
    return +p5d0(1 - t);
}

function p0d1(t: number): number {
    return -p5d1(1 - t);
}

function p0d2(t: number): number {
    return +p5d2(1 - t);
}

function p0d3(t: number): number {
    return -p5d3(1 - t);
}

function p0d4(t: number): number {
    return +p5d4(1 - t);
}

function p0d5(t: number): number {
    return -p5d5(1 - t);
}

function p1d0(t: number): number {
    return +p4d0(1 - t);
}

function p1d1(t: number): number {
    return -p4d1(1 - t);
}

function p1d2(t: number): number {
    return +p4d2(1 - t);
}

function p1d3(t: number): number {
    return -p4d3(1 - t);
}

function p1d4(t: number): number {
    return +p4d4(1 - t);
}

function p1d5(t: number): number {
    return -p4d5(1 - t);
}

function p2d0(t: number): number {
    return +p3d0(1 - t);
}

function p2d1(t: number): number {
    return -p3d1(1 - t);
}

function p2d2(t: number): number {
    return +p3d2(1 - t);
}

function p2d3(t: number): number {
    return -p3d3(1 - t);
}

function p2d4(t: number): number {
    return +p3d4(1 - t);
}

function p2d5(t: number): number {
    return -p3d5(1 - t);
}

function p3d0(t: number): number {
    return ((1 - t) * t) * ((1 - t) * t) * t;
}

function p3d1(t: number): number {
    return (1 - t) * t * t * (3 - 5 * t);
}

function p3d2(t: number): number {
    return t * (6 + t * (-24 + 20 * t));
}

function p3d3(t: number): number {
    return 6 + t * (-48 + 60 * t);
}

function p3d4(t: number): number {
    return -48 + 120 * t;
}

function p3d5(_t: number): number {
    return 120;
}

function p4d0(t: number): number {
    return (1 - t) * (t * t) * (t * t);
}

function p4d1(t: number): number {
    return t * t * t * (4 - 5 * t);
}

function p4d2(t: number): number {
    return t * t * (12 - 20 * t);
}

function p4d3(t: number): number {
    return t * (24 - 60 * t);
}

function p4d4(t: number): number {
    return 24 - 120 * t;
}

function p4d5(_t: number): number {
    return -120;
}

function p5d0(t: number): number {
    return t * (t * t) * (t * t);
}

function p5d1(t: number): number {
    return 5 * (t * t) * (t * t);
}

function p5d2(t: number): number {
    return 20 * t * (t * t);
}

function p5d3(t: number): number {
    return 60 * t * t;
}

function p5d4(t: number): number {
    return 120 * t;
}

function p5d5(_t: number): number {
    return 120;
}

const table: ReadonlyArray<ReadonlyArray<(t: number) => number>> = [
    [p0d0, p0d1, p0d2, p0d3, p0d4, p0d5],
    [p1d0, p1d1, p1d2, p1d3, p1d4, p1d5],
    [p2d0, p2d1, p2d2, p2d3, p2d4, p2d5],
    [p3d0, p3d1, p3d2, p3d3, p3d4, p3d5],
    [p4d0, p4d1, p4d2, p4d3, p4d4, p4d5],
    [p5d0, p5d1, p5d2, p5d3, p5d4, p5d5]
];

export class HermiteQuintic {
    // Set the coefficients manually as desired. For Hermite quintic
    // interpolation on a lattice, use generate(...). The lattice
    // interpolator is globally C1-continuous.
    private c: number[];

    // The default polynomial is identically zero. When 'blocks' is
    // provided, generate the 6x1 coefficients c[] for a cell of the
    // lattice with pixels at (x) and (x+1). The caller is responsible for
    // tracking the pixel (x) that is associated with the coefficients.
    constructor(blocks?: readonly [HermiteQuinticSample, HermiteQuinticSample]) {
        this.c = [0, 0, 0, 0, 0, 0];
        if (blocks !== undefined) {
            this.generate(blocks);
        }
    }

    // Evaluate the polynomial with the specified order. The returned
    // value is zero if xOrder >= 6. Otherwise, some examples are the
    // following where 'hermite' is of type HermiteQuintic:
    //   H(x) = hermite.evaluate(0, x)
    //   Hx(x) = hermite.evaluate(1, x)
    //   Hxx(x) = hermite.evaluate(2, x)
    //   Hxxx(x) = hermite.evaluate(3, x)
    //   Hxxxx(x) = hermite.evaluate(4, x)
    //   Hxxxxx(x) = hermite.evaluate(5, x)
    evaluate(xOrder: number, x: number): number {
        if (xOrder <= 5) {
            let result = 0;
            for (let i = 0; i < 6; ++i) {
                result += this.c[i] * HermiteQuintic.p(i, xOrder, x);
            }
            return result;
        }
        else {
            return 0;
        }
    }

    generate(blocks: readonly [HermiteQuinticSample, HermiteQuinticSample]): void {
        for (let b0 = 0; b0 <= 1; ++b0) {
            const z0 = 5 * b0 + 0;
            const p0 = 3 * b0 + 1;
            const q0 = 1 * b0 + 2;
            const s0 = 1 - 2 * b0;
            const b = blocks[b0];
            const input = new HermiteQuinticSample(b.f, s0 * b.fx, b.fxx);
            const [v0, v1, v2] = HermiteQuintic.generateSingle(input);
            this.c[z0] = v0;
            this.c[p0] = v1;
            this.c[q0] = v2;
        }
    }

    private static generateSingle(input: HermiteQuinticSample): [number, number, number] {
        const f = input.f;
        const fx = input.fx;
        const fxx = input.fxx;
        const v0 = f;
        const v1 = 5 * f + fx;
        const v2 = 10 * f + 4 * fx + 0.5 * fxx;
        return [v0, v1, v2];
    }

    // For internal use in HermiteBiquintic and HermiteTriquintic. The
    // 'select' parameter is the i argument for the polynomial and the
    // 'order' parameter is the order of the derivative.
    static p(select: number, order: number, t: number): number {
        return order <= 5 ? table[select][order](t) : 0;
    }
}
