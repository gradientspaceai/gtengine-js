// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FastGaussianBlur2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithms here are based on solving the linear heat equation using
// finite differences in scale, not in time. See FastGaussianBlur1.ts for the
// full description of the method (shared by dimensions 1, 2 and 3); the
// original document is
//   https://www.geometrictools.com/Documentation/FastGaussianBlur.pdf
//
// For iterative application of these functions, the caller is responsible
// for constructing a geometric sequence of scales,
//   s0, s1 = s0*b, s2 = s1*b = s0*b^2, ...
// where the base b satisfies 1 < b < exp(0.5*d) where d is the dimension of
// the image.  The upper bound on b guarantees stability of the finite
// difference method used to approximate the partial differential equation.
// The method assumes a pixel size of h = 1.
//
// Port notes: the upstream image type T is one of int16_t, int32_t, float or
// double with computations performed in double. Here the images are numeric
// arrays (number[] or typed arrays); computations are IEEE double. When the
// output is an Int16Array or Int32Array, the store performs the integer
// conversion that upstream's static_cast<T> performs.

// The input and output images must both have xBound*yBound elements and be
// stored in lexicographical order. The indexing is i = x + xBound * y.
export class FastGaussianBlur2 {
    private mXBound: number = 0;
    private mYBound: number = 0;
    private mInput: ArrayLike<number> | null = null;
    private mOutput: number[] | Float64Array | Float32Array | Int32Array | Int16Array | null = null;

    execute(xBound: number, yBound: number, input: ArrayLike<number>,
        output: number[] | Float64Array | Float32Array | Int32Array | Int16Array,
        scale: number, logBase: number): void {
        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mInput = input;
        this.mOutput = output;

        const xBoundM1 = xBound - 1, yBoundM1 = yBound - 1;
        for (let y = 0; y < yBound; ++y) {
            const ryps = y + scale;
            const ryms = y - scale;
            const yp1 = Math.floor(ryps);
            const ym1 = Math.ceil(ryms);

            for (let x = 0; x < xBound; ++x) {
                const rxps = x + scale;
                const rxms = x - scale;
                const xp1 = Math.floor(rxps);
                const xm1 = Math.ceil(rxms);

                const center = this.input(x, y);
                let xsum = -2.0 * center, ysum = xsum;

                // x portion of second central difference
                if (xp1 >= xBoundM1) {  // use boundary value
                    xsum += this.input(xBoundM1, y);
                } else {  // linearly interpolate
                    const imgXp1 = this.input(xp1, y);
                    const imgXp2 = this.input(xp1 + 1, y);
                    const delta = rxps - xp1;
                    xsum += imgXp1 + delta * (imgXp2 - imgXp1);
                }

                if (xm1 <= 0) {  // use boundary value
                    xsum += this.input(0, y);
                } else {  // linearly interpolate
                    const imgXm1 = this.input(xm1, y);
                    const imgXm2 = this.input(xm1 - 1, y);
                    const delta = rxms - xm1;
                    xsum += imgXm1 + delta * (imgXm1 - imgXm2);
                }

                // y portion of second central difference
                if (yp1 >= yBoundM1) {  // use boundary value
                    ysum += this.input(x, yBoundM1);
                } else {  // linearly interpolate
                    const imgYp1 = this.input(x, yp1);
                    const imgYp2 = this.input(x, yp1 + 1);
                    const delta = ryps - yp1;
                    ysum += imgYp1 + delta * (imgYp2 - imgYp1);
                }

                if (ym1 <= 0) {  // use boundary value
                    ysum += this.input(x, 0);
                } else {  // linearly interpolate
                    const imgYm1 = this.input(x, ym1);
                    const imgYm2 = this.input(x, ym1 - 1);
                    const delta = ryms - ym1;
                    ysum += imgYm1 + delta * (imgYm1 - imgYm2);
                }

                this.setOutput(x, y, center + logBase * (xsum + ysum));
            }
        }

        this.mXBound = 0;
        this.mYBound = 0;
        this.mInput = null;
        this.mOutput = null;
    }

    private input(x: number, y: number): number {
        return (this.mInput as ArrayLike<number>)[x + this.mXBound * y];
    }

    private setOutput(x: number, y: number, value: number): void {
        (this.mOutput as number[])[x + this.mXBound * y] = value;
    }
}
