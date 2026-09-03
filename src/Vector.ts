// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Vector.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream 'template <int32_t N, typename Real> class Vector'
// becomes a class with a runtime dimension backed by number[]. Canonical
// translations used throughout the library:
//
// - 'Vector<N, Real> v' -> 'new Vector(n)' (zero-filled; upstream leaves the
//   tuple uninitialized) or 'Vector.fromArray([...])'.
// - 'v[i]' -> 'v.get(i)' / 'v.set(i, x)'; the backing array is exposed as
//   'v.values' for hot inner loops (do not change its length).
// - Comparison operators -> equals, notEquals, lessThan, lessThanOrEqual,
//   greaterThan, greaterThanOrEqual (lexicographic, as std::array's).
// - Value semantics: C++ assignment copies; use 'v.clone()' where upstream
//   copies a vector.
// - Free operators -> module functions returning new vectors: unary
//   'operator-' -> negate; 'operator+/-' -> add, sub; scalar 'operator*/' ->
//   mul, div; componentwise 'operator*(v0,v1)' and 'operator/(v0,v1)' ->
//   compMul, compDiv. Compound assignments (+=, -=, ...) have no in-place
//   ports; write 'v0 = add(v0, v1)'. Unary 'operator+' is the identity and
//   has no port.
// - The geometric free functions keep their names in camelCase: dot, length,
//   normalize, orthonormalize, getOrthogonal, computeExtremes, hlift,
//   hproject, lift, project. normalize and orthonormalize mutate their
//   vector arguments in place and return the length, exactly as upstream.
// - Mismatched dimensions (a compile error upstream) throw via logAssert.

import { logAssert } from './Logger.js';

export class Vector {
    // The components of the vector. The contents may be read and written
    // directly, but the length must never be changed.
    readonly values: number[];

    // Create a vector of the given size with all components 0. (Upstream's
    // default constructor leaves the tuple uninitialized; the port
    // zero-fills for safety.)
    constructor(size: number) {
        this.values = new Array<number>(size).fill(0);
    }

    // Create a vector that copies the input values (the port of the
    // std::array and initializer-list constructors).
    static fromArray(values: readonly number[]): Vector {
        const result = new Vector(values.length);
        for (let i = 0; i < values.length; ++i) {
            result.values[i] = values[i];
        }
        return result;
    }

    // Create a vector with all components set to 'value' (the port of the
    // Vector(Real value) constructor).
    static filled(size: number, value: number): Vector {
        const result = new Vector(size);
        result.values.fill(value);
        return result;
    }

    // Special vectors. All components are 0.
    static zero(size: number): Vector {
        return new Vector(size);
    }

    // All components are 1.
    static ones(size: number): Vector {
        return Vector.filled(size, 1);
    }

    // Component d is 1, all others are 0. If d is invalid, the zero vector
    // is created. This is a convenience for creating the standard Euclidean
    // basis vectors.
    static unit(size: number, d: number): Vector {
        const result = new Vector(size);
        result.makeUnit(d);
        return result;
    }

    // Member access.
    get size(): number {
        return this.values.length;
    }

    get(i: number): number {
        return this.values[i];
    }

    set(i: number, value: number): void {
        this.values[i] = value;
    }

    // A deep copy (the port of C++ copy construction/assignment).
    clone(): Vector {
        return Vector.fromArray(this.values);
    }

    // All components are 0.
    makeZero(): void {
        this.values.fill(0);
    }

    // All components are 1.
    makeOnes(): void {
        this.values.fill(1);
    }

    // Component d is 1, all others are zero.
    makeUnit(d: number): void {
        this.values.fill(0);
        if (0 <= d && d < this.values.length) {
            this.values[d] = 1;
        }
    }

    // Comparisons for sorted containers and geometric ordering
    // (lexicographic, matching std::array's operators). Comparing vectors of
    // different sizes is a compile error upstream and throws here.
    private compare(vec: Vector): number {
        logAssert(this.values.length === vec.values.length,
            'Vector: mismatched sizes.');
        for (let i = 0; i < this.values.length; ++i) {
            if (this.values[i] < vec.values[i]) {
                return -1;
            }
            if (this.values[i] > vec.values[i]) {
                return +1;
            }
        }
        return 0;
    }

    equals(vec: Vector): boolean {
        return this.compare(vec) === 0;
    }

    notEquals(vec: Vector): boolean {
        return this.compare(vec) !== 0;
    }

    lessThan(vec: Vector): boolean {
        return this.compare(vec) < 0;
    }

    lessThanOrEqual(vec: Vector): boolean {
        return this.compare(vec) <= 0;
    }

    greaterThan(vec: Vector): boolean {
        return this.compare(vec) > 0;
    }

    greaterThanOrEqual(vec: Vector): boolean {
        return this.compare(vec) >= 0;
    }
}

function assertSameSize(v0: Vector, v1: Vector): void {
    logAssert(v0.size === v1.size, 'Vector: mismatched sizes.');
}

// Unary operations. The port of unary 'operator-'.
export function negate(v: Vector): Vector {
    const result = new Vector(v.size);
    for (let i = 0; i < v.size; ++i) {
        result.values[i] = -v.values[i];
    }
    return result;
}

// Linear-algebraic operations.
export function add(v0: Vector, v1: Vector): Vector {
    assertSameSize(v0, v1);
    const result = new Vector(v0.size);
    for (let i = 0; i < v0.size; ++i) {
        result.values[i] = v0.values[i] + v1.values[i];
    }
    return result;
}

export function sub(v0: Vector, v1: Vector): Vector {
    assertSameSize(v0, v1);
    const result = new Vector(v0.size);
    for (let i = 0; i < v0.size; ++i) {
        result.values[i] = v0.values[i] - v1.values[i];
    }
    return result;
}

// The port of 'v * scalar' and 'scalar * v'; both argument orders are
// accepted.
export function mul(v: Vector, scalar: number): Vector;
export function mul(scalar: number, v: Vector): Vector;
export function mul(arg0: Vector | number, arg1: Vector | number): Vector {
    const v = (typeof arg0 === 'number' ? arg1 : arg0) as Vector;
    const scalar = (typeof arg0 === 'number' ? arg0 : arg1) as number;
    const result = new Vector(v.size);
    for (let i = 0; i < v.size; ++i) {
        result.values[i] = v.values[i] * scalar;
    }
    return result;
}

// The port of 'v / scalar'. As upstream, division by zero produces the zero
// vector, and the division is performed as multiplication by 1/scalar.
export function div(v: Vector, scalar: number): Vector {
    const result = new Vector(v.size);
    if (scalar !== 0) {
        const invScalar = 1 / scalar;
        for (let i = 0; i < v.size; ++i) {
            result.values[i] = v.values[i] * invScalar;
        }
    }
    return result;
}

// Componentwise algebraic operations (the ports of 'operator*(v0, v1)' and
// 'operator/(v0, v1)').
export function compMul(v0: Vector, v1: Vector): Vector {
    assertSameSize(v0, v1);
    const result = new Vector(v0.size);
    for (let i = 0; i < v0.size; ++i) {
        result.values[i] = v0.values[i] * v1.values[i];
    }
    return result;
}

export function compDiv(v0: Vector, v1: Vector): Vector {
    assertSameSize(v0, v1);
    const result = new Vector(v0.size);
    for (let i = 0; i < v0.size; ++i) {
        result.values[i] = v0.values[i] / v1.values[i];
    }
    return result;
}

// Geometric operations. The functions with 'robust' set to 'false' use the
// standard algorithm for normalizing a vector by computing the length as a
// square root of the squared length and dividing by it. The results can be
// infinite (or NaN) if the length is zero. When 'robust' is set to 'true',
// the algorithm is designed to avoid floating-point overflow and sets the
// normalized vector to zero when the length is zero.
export function dot(v0: Vector, v1: Vector): number {
    assertSameSize(v0, v1);
    let result = v0.values[0] * v1.values[0];
    for (let i = 1; i < v0.size; ++i) {
        result += v0.values[i] * v1.values[i];
    }
    return result;
}

export function length(v: Vector, robust: boolean = false): number {
    if (robust) {
        let maxAbsComp = Math.abs(v.values[0]);
        for (let i = 1; i < v.size; ++i) {
            const absComp = Math.abs(v.values[i]);
            if (absComp > maxAbsComp) {
                maxAbsComp = absComp;
            }
        }

        if (maxAbsComp > 0) {
            const scaled = div(v, maxAbsComp);
            return maxAbsComp * Math.sqrt(dot(scaled, scaled));
        }
        return 0;
    }
    return Math.sqrt(dot(v, v));
}

// Normalizes v in place and returns the length that v had on input (zero
// length leaves/sets v to the zero vector).
export function normalize(v: Vector, robust: boolean = false): number {
    if (robust) {
        let maxAbsComp = Math.abs(v.values[0]);
        for (let i = 1; i < v.size; ++i) {
            const absComp = Math.abs(v.values[i]);
            if (absComp > maxAbsComp) {
                maxAbsComp = absComp;
            }
        }

        if (maxAbsComp > 0) {
            const invMax = 1 / maxAbsComp;
            for (let i = 0; i < v.size; ++i) {
                v.values[i] *= invMax;
            }
            let len = Math.sqrt(dot(v, v));
            const invLen = 1 / len;
            for (let i = 0; i < v.size; ++i) {
                v.values[i] *= invLen;
            }
            len *= maxAbsComp;
            return len;
        }
        v.makeZero();
        return 0;
    }

    const len = Math.sqrt(dot(v, v));
    if (len > 0) {
        const invLen = 1 / len;
        for (let i = 0; i < v.size; ++i) {
            v.values[i] *= invLen;
        }
    } else {
        v.makeZero();
    }
    return len;
}

// Gram-Schmidt orthonormalization to generate orthonormal vectors from the
// linearly independent inputs. The function returns the smallest length of
// the unnormalized vectors computed during the process. If this value is
// nearly zero, it is possible that the inputs are linearly dependent (within
// numerical round-off errors). On input, 1 <= numInputs <= v[0].size and
// v[0] through v[numInputs-1] must be initialized. On output, the vectors
// v[0] through v[numInputs-1] form an orthonormal set (mutated in place).
export function orthonormalize(numInputs: number, v: Vector[],
    robust: boolean = false): number {
    if (v.length >= numInputs && 1 <= numInputs && numInputs <= v[0].size) {
        let minLength = normalize(v[0], robust);
        for (let i = 1; i < numInputs; ++i) {
            for (let j = 0; j < i; ++j) {
                const d = dot(v[i], v[j]);
                for (let k = 0; k < v[i].size; ++k) {
                    v[i].values[k] -= v[j].values[k] * d;
                }
            }
            const len = normalize(v[i], robust);
            if (len < minLength) {
                minLength = len;
            }
        }
        return minLength;
    }

    return 0;
}

// Construct a single vector orthogonal to the nonzero input vector. If the
// maximum absolute component occurs at index i, then the orthogonal vector
// U has u[i] = v[i+1], u[i+1] = -v[i], and all other components zero. The
// index addition i+1 is computed modulo N.
export function getOrthogonal(v: Vector, unitLength: boolean): Vector {
    const n = v.size;
    let cmax = Math.abs(v.values[0]);
    let imax = 0;
    for (let i = 1; i < n; ++i) {
        const c = Math.abs(v.values[i]);
        if (c > cmax) {
            cmax = c;
            imax = i;
        }
    }

    const result = new Vector(n);
    let inext = imax + 1;
    if (inext === n) {
        inext = 0;
    }
    result.values[imax] = v.values[inext];
    result.values[inext] = -v.values[imax];
    if (unitLength) {
        const sqrDistance = result.values[imax] * result.values[imax]
            + result.values[inext] * result.values[inext];
        const invLength = 1 / Math.sqrt(sqrDistance);
        result.values[imax] *= invLength;
        result.values[inext] *= invLength;
    }
    return result;
}

// Compute the axis-aligned bounding box of the vectors. Upstream returns
// 'true' and fills vmin/vmax when the inputs are valid; the port returns
// { vmin, vmax } or null when the input array is empty.
export function computeExtremes(vectors: readonly Vector[]):
    { vmin: Vector, vmax: Vector } | null {
    if (vectors.length === 0) {
        return null;
    }

    const vmin = vectors[0].clone();
    const vmax = vmin.clone();
    for (let j = 1; j < vectors.length; ++j) {
        const vec = vectors[j];
        assertSameSize(vmin, vec);
        for (let i = 0; i < vmin.size; ++i) {
            if (vec.values[i] < vmin.values[i]) {
                vmin.values[i] = vec.values[i];
            } else if (vec.values[i] > vmax.values[i]) {
                vmax.values[i] = vec.values[i];
            }
        }
    }
    return { vmin, vmax };
}

// Lift n-tuple v to homogeneous (n+1)-tuple (v, last).
export function hlift(v: Vector, last: number): Vector {
    const n = v.size;
    const result = new Vector(n + 1);
    for (let i = 0; i < n; ++i) {
        result.values[i] = v.values[i];
    }
    result.values[n] = last;
    return result;
}

// Project homogeneous n-tuple v = (u, v[n-1]) to (n-1)-tuple u.
export function hproject(v: Vector): Vector {
    const n = v.size;
    logAssert(n >= 2, 'Invalid dimension.');
    const result = new Vector(n - 1);
    for (let i = 0; i < n - 1; ++i) {
        result.values[i] = v.values[i];
    }
    return result;
}

// Lift n-tuple v = (w0, w1) to (n+1)-tuple u = (w0, u[inject], w1). By
// inference, w0 is a (inject)-tuple [nonexistent when inject=0] and w1 is a
// (n-inject)-tuple [nonexistent when inject=n].
export function lift(v: Vector, inject: number, value: number): Vector {
    const n = v.size;
    const result = new Vector(n + 1);
    let i: number;
    for (i = 0; i < inject; ++i) {
        result.values[i] = v.values[i];
    }
    result.values[i] = value;
    let j = i + 1;
    for (; i < n; ++i, ++j) {
        result.values[j] = v.values[i];
    }
    return result;
}

// Project n-tuple v = (w0, v[reject], w1) to (n-1)-tuple u = (w0, w1). By
// inference, w0 is a (reject)-tuple [nonexistent when reject=0] and w1 is a
// (n-1-reject)-tuple [nonexistent when reject=n-1].
export function project(v: Vector, reject: number): Vector {
    const n = v.size;
    logAssert(n >= 2, 'Invalid dimension.');
    const result = new Vector(n - 1);
    for (let i = 0, j = 0; i < n - 1; ++i, ++j) {
        if (j === reject) {
            ++j;
        }
        result.values[i] = v.values[j];
    }
    return result;
}
