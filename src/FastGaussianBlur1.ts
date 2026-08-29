// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FastGaussianBlur1.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithms here are based on solving the linear heat equation using
// finite differences in scale, not in time.  The following document has
// a brief summary of the concept,
//   https://www.geometrictools.com/Documentation/FastGaussianBlur.pdf
// The idea is to represent the blurred image as f(x,s) in terms of position
// x and scale s.  Gaussian blurring is accomplished by using the input image
// I(x,s0) as the initial image (of scale s0 > 0) for the partial differential
// equation
//   s*df/ds = s^2*Laplacian(f)
// where the Laplacian operator is
//   Laplacian = (d/dx)^2, dimension 1
//   Laplacian = (d/dx)^2+(d/dy)^2, dimension 2
//   Laplacian = (d/dx)^2+(d/dy)^2+(d/dz)^2, dimension 3
//
// The term s*df/ds is approximated by
//   s*df(x,s)/ds = (f(x,b*s)-f(x,s))/ln(b)
// for b > 1, but close to 1, where ln(b) is the natural logarithm of b.  If
// you take the limit of the right-hand side as b approaches 1, you get the
// left-hand side.
//
// The term s^2*((d/dx)^2)f is approximated by
//   s^2*((d/dx)^2)f = (f(x+h*s,s)-2*f(x,s)+f(x-h*s,s))/h^2
// for h > 0, but close to zero.
//
// Equating the approximations for the left-hand side and the right-hand side
// of the partial differential equation leads to the numerical method used in
// this code.
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

// The input and output images must both have xBound elements.
export class FastGaussianBlur1 {
    execute(xBound: number, input: ArrayLike<number>,
        output: number[] | Float64Array | Float32Array | Int32Array | Int16Array,
        scale: number, logBase: number): void {
        const xBoundM1 = xBound - 1;
        for (let x = 0; x < xBound; ++x) {
            const rxps = x + scale;
            const rxms = x - scale;
            const xp1 = Math.floor(rxps);
            const xm1 = Math.ceil(rxms);

            const center = input[x];
            let xsum = -2.0 * center;

            if (xp1 >= xBoundM1) {  // use boundary value
                xsum += input[xBoundM1];
            } else {  // linearly interpolate
                const imgXp1 = input[xp1];
                const imgXp2 = input[xp1 + 1];
                const delta = rxps - xp1;
                xsum += imgXp1 + delta * (imgXp2 - imgXp1);
            }

            if (xm1 <= 0) {  // use boundary value
                xsum += input[0];
            } else {  // linearly interpolate
                const imgXm1 = input[xm1];
                const imgXm2 = input[xm1 - 1];
                const delta = rxms - xm1;
                xsum += imgXm1 + delta * (imgXm1 - imgXm2);
            }

            output[x] = center + logBase * xsum;
        }
    }
}
