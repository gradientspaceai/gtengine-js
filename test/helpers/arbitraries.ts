// Shared fast-check arbitraries and helpers for gtengine-js verification tests.
// Every generator produces finite, moderately scaled values so that properties
// stay numerically meaningful (no NaN, no overflow, no 1e300 coordinates).
import fc from 'fast-check';
import { expect } from 'vitest';
import { Vector, dot, length, normalize, orthonormalize } from '../../src/Vector.js';
import { Matrix } from '../../src/Matrix.js';
import { AlignedBox } from '../../src/AlignedBox.js';
import { OrientedBox } from '../../src/OrientedBox.js';
import { Line } from '../../src/Line.js';
import { Ray } from '../../src/Ray.js';
import { Segment } from '../../src/Segment.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Hyperplane } from '../../src/Hyperplane.js';
import { Triangle } from '../../src/Triangle.js';

// ---- scalars ---------------------------------------------------------------

/** Finite double in [min, max]. */
export const finite = (min = -10, max = 10): fc.Arbitrary<number> =>
    fc.double({ min, max, noNaN: true, noDefaultInfinity: true });

/** Finite double in [min, max] excluding a band around zero. */
export const nonzero = (min = -10, max = 10, eps = 1e-3): fc.Arbitrary<number> =>
    finite(min, max).filter(x => Math.abs(x) > eps);

/** Strictly positive double in (eps, max]. */
export const positive = (max = 10, eps = 1e-3): fc.Arbitrary<number> =>
    finite(eps, max).filter(x => x > eps);

/**
 * Finite double in [min, max] with magnitudes below eps snapped to exactly 0.
 * fast-check's double() samples the bit patterns of the range uniformly, so it
 * produces subnormals (1e-320 and smaller) very often. Products of such values
 * underflow, which makes relative error tolerances meaningless. Use this where
 * a property compares two floating-point computations of the same quantity.
 */
export const wellScaled = (min = -10, max = 10, eps = 1e-3):
    fc.Arbitrary<number> =>
    finite(min, max).map(x => (Math.abs(x) < eps ? 0 : x));

// ---- vectors ---------------------------------------------------------------

/** Vector of dimension n with components in [min, max]. */
export const vector = (n: number, min = -10, max = 10): fc.Arbitrary<Vector> =>
    fc.array(finite(min, max), { minLength: n, maxLength: n })
        .map(a => Vector.fromArray(a));

/** Vector of dimension n whose components are wellScaled(min, max). */
export const wellScaledVector = (n: number, min = -10, max = 10):
    fc.Arbitrary<Vector> =>
    fc.array(wellScaled(min, max), { minLength: n, maxLength: n })
        .map(a => Vector.fromArray(a));

// ---- matrices --------------------------------------------------------------

/** numRows-by-numCols matrix with row-major elements in [min, max]. */
export const matrix = (numRows: number, numCols: number, min = -10, max = 10):
    fc.Arbitrary<Matrix> =>
    fc.array(finite(min, max),
        { minLength: numRows * numCols, maxLength: numRows * numCols })
        .map(a => Matrix.fromArray(numRows, numCols, a));

/** numRows-by-numCols matrix whose elements are wellScaled(min, max). */
export const wellScaledMatrix = (numRows: number, numCols: number,
    min = -10, max = 10): fc.Arbitrary<Matrix> =>
    fc.array(wellScaled(min, max),
        { minLength: numRows * numCols, maxLength: numRows * numCols })
        .map(a => Matrix.fromArray(numRows, numCols, a));

/** Unit-length vector of dimension n. */
export const unitVector = (n: number): fc.Arbitrary<Vector> =>
    vector(n, -1, 1).filter(v => length(v) > 1e-2).map(v => {
        const u = v.clone();
        normalize(u);
        return u;
    });

/** Orthonormal frame: n mutually orthogonal unit vectors of dimension n. */
export const orthonormalFrame = (n: number): fc.Arbitrary<Vector[]> =>
    fc.array(vector(n, -1, 1), { minLength: n, maxLength: n })
        .map(vs => vs.map(v => v.clone()))
        .filter(vs => orthonormalize(n, vs) > 1e-2)
        .map(vs => {
            // orthonormalize mutated the filtered copy; redo on a fresh copy so
            // shrinking stays deterministic.
            const f = vs.map(v => v.clone());
            orthonormalize(n, f);
            return f;
        });

// ---- primitives ------------------------------------------------------------

export const alignedBox = (n: number, min = -10, max = 10): fc.Arbitrary<AlignedBox> =>
    fc.tuple(vector(n, min, max), vector(n, min, max)).map(([a, b]) => {
        const lo = new Vector(n);
        const hi = new Vector(n);
        for (let i = 0; i < n; ++i) {
            lo.set(i, Math.min(a.get(i), b.get(i)));
            hi.set(i, Math.max(a.get(i), b.get(i)));
        }
        return AlignedBox.fromMinMax(lo, hi);
    });

export const orientedBox = (n: number): fc.Arbitrary<OrientedBox> =>
    fc.tuple(vector(n), orthonormalFrame(n),
        fc.array(positive(5), { minLength: n, maxLength: n }))
        .map(([c, axis, ext]) =>
            OrientedBox.fromCenterAxisExtent(c, axis, Vector.fromArray(ext)));

export const line = (n: number): fc.Arbitrary<Line> =>
    fc.tuple(vector(n), unitVector(n))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

export const ray = (n: number): fc.Arbitrary<Ray> =>
    fc.tuple(vector(n), unitVector(n))
        .map(([o, d]) => Ray.fromOriginDirection(o, d));

export const segment = (n: number): fc.Arbitrary<Segment> =>
    fc.tuple(vector(n), vector(n))
        .filter(([p0, p1]) => {
            let s = 0;
            for (let i = 0; i < n; ++i) { s += (p1.get(i) - p0.get(i)) ** 2; }
            return s > 1e-4;
        })
        .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

export const sphere = (n: number): fc.Arbitrary<Hypersphere> =>
    fc.tuple(vector(n), positive(5))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

export const plane = (n: number): fc.Arbitrary<Hyperplane> =>
    fc.tuple(unitVector(n), vector(n))
        .map(([nrm, o]) => Hyperplane.fromNormalOrigin(nrm, o));

/** Non-degenerate triangle (area bounded away from zero). */
export const triangle = (n: number): fc.Arbitrary<Triangle> =>
    fc.tuple(vector(n), vector(n), vector(n))
        .filter(([a, b, c]) => {
            const e0 = new Vector(n);
            const e1 = new Vector(n);
            for (let i = 0; i < n; ++i) {
                e0.set(i, b.get(i) - a.get(i));
                e1.set(i, c.get(i) - a.get(i));
            }
            const d00 = dot(e0, e0), d11 = dot(e1, e1), d01 = dot(e0, e1);
            return d00 * d11 - d01 * d01 > 1e-3;   // squared (2*area)
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

// ---- deterministic pseudo-random for brute-force cross-checks -------------

/** mulberry32: small, fast, seedable PRNG returning doubles in [0, 1). */
export function seededRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// ---- assertions ------------------------------------------------------------

/** |a - b| <= abs + rel * max(|a|, |b|). */
export function expectClose(a: number, b: number, abs = 1e-9, rel = 1e-9): void {
    const tol = abs + rel * Math.max(Math.abs(a), Math.abs(b));
    expect(Math.abs(a - b), `expected ${a} ≈ ${b} (tol ${tol})`).toBeLessThanOrEqual(tol);
}

export function expectVectorClose(a: Vector, b: Vector, abs = 1e-9, rel = 1e-9): void {
    expect(a.size).toBe(b.size);
    for (let i = 0; i < a.size; ++i) { expectClose(a.get(i), b.get(i), abs, rel); }
}

/** Run a fast-check property with the project defaults (deterministic count). */
export function check<T>(arb: fc.Arbitrary<T>, predicate: (t: T) => void | boolean,
    numRuns = 200): void {
    fc.assert(fc.property(arb, predicate), { numRuns });
}

export { fc };
