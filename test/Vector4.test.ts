import { describe, it, expect } from 'vitest';
import {
    check, vector, wellScaledVector, expectClose, expectVectorClose, fc
} from './helpers/arbitraries.js';
import {
    hyperCross, unitHyperCross, dotHyperCross, computeOrthogonalComplement4
} from '../src/Vector4.js';
import { Vector, dot, length, mul, negate } from '../src/Vector.js';
import { Matrix, determinant } from '../src/Matrix.js';

const e0 = () => Vector.fromArray([1, 0, 0, 0]);
const e1 = () => Vector.fromArray([0, 1, 0, 0]);
const e2 = () => Vector.fromArray([0, 0, 1, 0]);
const e3 = () => Vector.fromArray([0, 0, 0, 1]);

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296 - 0.5;
    };
}

function randomVector4(rand: () => number): Vector {
    return Vector.fromArray([4 * rand(), 4 * rand(), 4 * rand(), 4 * rand()]);
}

// The hypercross product is the formal determinant
//   det{{e0,e1,e2,e3},{v0},{v1},{v2}},
// so Dot(HyperCross(v0,v1,v2), v3) is the determinant of the matrix whose
// rows are v3, v0, v1 and v2. This is the independent cross-check used
// below.
function rowDeterminant(v3: Vector, v0: Vector, v1: Vector,
    v2: Vector): number {
    const M = new Matrix(4, 4);
    const rows = [v3, v0, v1, v2];
    for (let r = 0; r < 4; ++r) {
        for (let c = 0; c < 4; ++c) {
            M.set(r, c, rows[r].values[c]);
        }
    }
    return determinant(M);
}

describe('hyperCross', () => {
    it('has its known value on the standard basis', () => {
        // det{{e0,e1,e2,e3},{e0},{e1},{e2}} = -e3.
        expect(hyperCross(e0(), e1(), e2()).values).toEqual([0, 0, 0, -1]);
        expect(hyperCross(e1(), e2(), e3()).values).toEqual([1, 0, 0, 0]);
        expect(hyperCross(e0(), e1(), e3()).values).toEqual([0, 0, 1, 0]);
    });

    it('is orthogonal to each of its inputs', () => {
        const rand = makeRandom(2718);
        for (let trial = 0; trial < 10; ++trial) {
            const v0 = randomVector4(rand);
            const v1 = randomVector4(rand);
            const v2 = randomVector4(rand);
            const h = hyperCross(v0, v1, v2);
            expect(dot(h, v0)).toBeCloseTo(0, 10);
            expect(dot(h, v1)).toBeCloseTo(0, 10);
            expect(dot(h, v2)).toBeCloseTo(0, 10);
        }
    });

    it('is alternating: swapping two inputs negates the result', () => {
        const rand = makeRandom(1414);
        const v0 = randomVector4(rand);
        const v1 = randomVector4(rand);
        const v2 = randomVector4(rand);
        const h = negate(hyperCross(v0, v1, v2));
        for (let i = 0; i < 4; ++i) {
            expect(hyperCross(v1, v0, v2).values[i]).toBeCloseTo(h.values[i], 12);
            expect(hyperCross(v0, v2, v1).values[i]).toBeCloseTo(h.values[i], 12);
            expect(hyperCross(v2, v1, v0).values[i]).toBeCloseTo(h.values[i], 12);
        }
    });

    it('vanishes for linearly dependent inputs', () => {
        const rand = makeRandom(31415);
        const v0 = randomVector4(rand);
        const v1 = randomVector4(rand);
        for (const value of hyperCross(v0, v1, v0).values) {
            expect(Math.abs(value)).toBeLessThan(1e-12);
        }
        const v2 = Vector.fromArray(v0.values.map(
            (x, i) => 2 * x + 3 * v1.values[i]));
        const h = hyperCross(v0, v1, v2);
        for (const value of h.values) {
            expect(Math.abs(value)).toBeLessThan(1e-12);
        }
    });

    it('is multilinear in its first argument', () => {
        const rand = makeRandom(6180);
        const a = randomVector4(rand);
        const b = randomVector4(rand);
        const v1 = randomVector4(rand);
        const v2 = randomVector4(rand);
        const sum = Vector.fromArray(a.values.map(
            (x, i) => x + 2 * b.values[i]));
        const h = hyperCross(sum, v1, v2);
        const ha = hyperCross(a, v1, v2);
        const hb = hyperCross(b, v1, v2);
        for (let i = 0; i < 4; ++i) {
            expect(h.values[i]).toBeCloseTo(
                ha.values[i] + 2 * hb.values[i], 10);
        }
    });

    it('requires 4-tuples', () => {
        expect(() => hyperCross(Vector.fromArray([1, 2, 3]), e1(), e2()))
            .toThrow('Vector4: vector must have size 4.');
        expect(() => hyperCross(e0(), Vector.fromArray([1, 2, 3]), e2()))
            .toThrow('Vector4: vector must have size 4.');
        expect(() => hyperCross(e0(), e1(), Vector.fromArray([1, 2, 3])))
            .toThrow('Vector4: vector must have size 4.');
    });
});

describe('unitHyperCross', () => {
    it('is the normalized hypercross product', () => {
        const rand = makeRandom(9001);
        const v0 = randomVector4(rand);
        const v1 = randomVector4(rand);
        const v2 = randomVector4(rand);
        const u = unitHyperCross(v0, v1, v2);
        expect(length(u)).toBeCloseTo(1, 12);
        const h = hyperCross(v0, v1, v2);
        const len = length(h);
        for (let i = 0; i < 4; ++i) {
            expect(u.values[i]).toBeCloseTo(h.values[i] / len, 10);
        }
        // The robust path produces the same unit-length vector.
        const ur = unitHyperCross(v0, v1, v2, true);
        for (let i = 0; i < 4; ++i) {
            expect(ur.values[i]).toBeCloseTo(u.values[i], 10);
        }
    });

    it('is zero when the inputs are linearly dependent and robust is set',
        () => {
            expect(unitHyperCross(e0(), e1(), e0(), true).values)
                .toEqual([0, 0, 0, 0]);
        });
});

describe('dotHyperCross', () => {
    it('has its known value on the standard basis', () => {
        expect(dotHyperCross(e0(), e1(), e2(), e3())).toBe(-1);
        expect(dotHyperCross(e1(), e2(), e3(), e0())).toBe(1);
    });

    it('agrees with the 4x4 determinant of the rows v3, v0, v1, v2', () => {
        const rand = makeRandom(20260831);
        for (let trial = 0; trial < 20; ++trial) {
            const v0 = randomVector4(rand);
            const v1 = randomVector4(rand);
            const v2 = randomVector4(rand);
            const v3 = randomVector4(rand);
            expect(dotHyperCross(v0, v1, v2, v3))
                .toBeCloseTo(rowDeterminant(v3, v0, v1, v2), 8);
        }
    });

    it('vanishes when two arguments are equal', () => {
        const rand = makeRandom(4711);
        const v0 = randomVector4(rand);
        const v1 = randomVector4(rand);
        const v2 = randomVector4(rand);
        expect(dotHyperCross(v0, v1, v2, v0)).toBeCloseTo(0, 10);
        expect(dotHyperCross(v0, v1, v2, v2)).toBeCloseTo(0, 10);
    });
});

describe('computeOrthogonalComplement4', () => {
    function expectOrthonormalBasis(v: Vector[]): void {
        for (let i = 0; i < 4; ++i) {
            expect(length(v[i])).toBeCloseTo(1, 10);
            for (let j = i + 1; j < 4; ++j) {
                expect(dot(v[i], v[j])).toBeCloseTo(0, 10);
            }
        }
        // The basis produced by the algorithm satisfies
        // Dot(HyperCross(v0,v1,v2), v3) = +1 because v3 is assigned the
        // hypercross of the first three vectors.
        expect(dotHyperCross(v[0], v[1], v[2], v[3])).toBeCloseTo(1, 10);
    }

    it('extends one input to an orthonormal basis', () => {
        const rand = makeRandom(101);
        for (let trial = 0; trial < 10; ++trial) {
            const v: Vector[] = [randomVector4(rand), new Vector(4),
                new Vector(4), new Vector(4)];
            const minLength = computeOrthogonalComplement4(1, v);
            expect(minLength).toBeGreaterThan(0);
            expectOrthonormalBasis(v);
        }
    });

    it('extends two inputs to an orthonormal basis', () => {
        const rand = makeRandom(202);
        for (let trial = 0; trial < 10; ++trial) {
            const v: Vector[] = [randomVector4(rand), randomVector4(rand),
                new Vector(4), new Vector(4)];
            const minLength = computeOrthogonalComplement4(2, v, true);
            expect(minLength).toBeGreaterThan(0);
            expectOrthonormalBasis(v);
        }
    });

    it('extends three inputs to an orthonormal basis', () => {
        const rand = makeRandom(303);
        for (let trial = 0; trial < 10; ++trial) {
            const v: Vector[] = [randomVector4(rand), randomVector4(rand),
                randomVector4(rand), new Vector(4)];
            const minLength = computeOrthogonalComplement4(3, v);
            expect(minLength).toBeGreaterThan(0);
            expectOrthonormalBasis(v);
        }
    });

    it('produces a known basis for the first standard basis vector', () => {
        const v: Vector[] = [e0(), new Vector(4), new Vector(4),
            new Vector(4)];
        computeOrthogonalComplement4(1, v);
        const expected = [
            [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1], [0, 0, -1, 0]
        ];
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 4; ++j) {
                expect(v[i].values[j]).toBeCloseTo(expected[i][j], 15);
            }
        }
    });

    it('handles a 3D vector represented as a 4D affine vector (w = 0)', () => {
        const v: Vector[] = [Vector.fromArray([0, 0, 5, 0]), new Vector(4),
            new Vector(4), new Vector(4)];
        const minLength = computeOrthogonalComplement4(1, v);
        expect(minLength).toBeGreaterThan(0);
        expectOrthonormalBasis(v);
    });

    it('reports linear dependence with a zero return value', () => {
        const v: Vector[] = [e0(), e0(), new Vector(4), new Vector(4)];
        expect(computeOrthogonalComplement4(2, v)).toBe(0);
    });

    it('returns 0 for an unsupported number of inputs', () => {
        const v: Vector[] = [e0(), e1(), e2(), e3()];
        expect(computeOrthogonalComplement4(0, v)).toBe(0);
        expect(computeOrthogonalComplement4(4, v)).toBe(0);
    });

    it('preserves the upstream degeneracy when the largest component is the '
        + 'last and components 1 and 2 are zero', () => {
        // Upstream takes a special branch when the largest-magnitude
        // component is at index 3, setting v[1] = (0, x2, -x1, 0). That
        // vector is zero whenever x1 = x2 = 0, so no basis can be built.
        // The port preserves this behavior; the zero return value signals
        // the failure to the caller. See the PR "Upstream bug suspects".
        const v: Vector[] = [Vector.fromArray([1, 0, 0, 5]), new Vector(4),
            new Vector(4), new Vector(4)];
        expect(computeOrthogonalComplement4(1, v)).toBe(0);
        expect(v[1].values).toEqual([0, 0, 0, 0]);
    });

    it('requires 4-tuples', () => {
        const v: Vector[] = [Vector.fromArray([1, 2, 3]), new Vector(4),
            new Vector(4), new Vector(4)];
        expect(() => computeOrthogonalComplement4(1, v))
            .toThrow('Vector4: vector must have size 4.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream Vector4.h.
// ---------------------------------------------------------------------------

describe('Vector4 verification', () => {
    it('hyperCross negates exactly when its first two arguments swap', () => {
        check(fc.tuple(vector(4), vector(4), vector(4)), ([a, b, c]) => {
            // Swapping v0 and v1 negates every 2x2 minor exactly (x*y is
            // commutative and u - v is the exact negation of v - u), and
            // round-to-nearest is symmetric under negation, so the whole
            // result negates bit for bit. Swaps involving v2 rearrange which
            // products are formed and are only algebraically antisymmetric;
            // they are checked with a tolerance below.
            const h = hyperCross(a, b, c);
            const swapped = hyperCross(b, a, c);
            for (let i = 0; i < 4; ++i) {
                // '+ 0' normalizes the -0 that -0 + -0 produces.
                expect(h.get(i) + swapped.get(i) + 0).toBe(0);
            }
            // Repeating the first two arguments makes every 2x2 minor
            // exactly zero, so the result is exactly the zero vector.
            expect(hyperCross(a, a, c).values.map(x => x + 0))
                .toEqual([0, 0, 0, 0]);
        });
    });

    it('hyperCross is orthogonal to each of its three arguments', () => {
        check(fc.tuple(wellScaledVector(4), wellScaledVector(4),
            wellScaledVector(4)), ([a, b, c]) => {
            const h = hyperCross(a, b, c);
            // Dot(h, a) is the determinant of a matrix with a repeated row:
            // exactly zero in exact arithmetic, so the residual is pure
            // round-off, bounded relative to the product of the row norms.
            const scale = length(a) * length(b) * length(c);
            expectClose(dot(h, a), 0, 1e-12 * scale * length(a), 0);
            expectClose(dot(h, b), 0, 1e-12 * scale * length(b), 0);
            expectClose(dot(h, c), 0, 1e-12 * scale * length(c), 0);
            // Repeating an argument other than the first two only cancels up
            // to rounding, since the 2x2 minors are not exactly zero there.
            expectVectorClose(hyperCross(a, b, b), new Vector(4),
                1e-12 * length(a) * length(b) * length(b), 0);
            expectVectorClose(hyperCross(a, b, a), new Vector(4),
                1e-12 * length(a) * length(a) * length(b), 0);
            // Swapping an argument into or out of the v2 slot rearranges the
            // products, so this antisymmetry holds only up to rounding.
            for (const swapped of [hyperCross(a, c, b), hyperCross(c, b, a)]) {
                for (let i = 0; i < 4; ++i) {
                    expectClose(h.get(i) + swapped.get(i), 0, 1e-12 * scale, 0);
                }
            }
        });
    });

    it('dotHyperCross is the 4x4 determinant of the rows v3, v0, v1, v2',
        () => {
            check(fc.tuple(wellScaledVector(4), wellScaledVector(4),
                wellScaledVector(4), wellScaledVector(4)),
                ([v0, v1, v2, v3]) => {
                    const M = Matrix.fromArray(4, 4, [...v3.values,
                        ...v0.values, ...v1.values, ...v2.values]);
                    // Gaussian elimination with full pivoting is an
                    // independent algorithm; the tolerance is relative to the
                    // Hadamard bound (the product of the row norms).
                    const bound = length(v0) * length(v1) * length(v2)
                        * length(v3);
                    expectClose(dotHyperCross(v0, v1, v2, v3), determinant(M),
                        1e-11 * bound, 0);
                    // Repeating the first two arguments zeroes every 2x2
                    // minor exactly, so the triple product is exactly zero.
                    expect(dotHyperCross(v0, v0, v2, v3) + 0).toBe(0);
                    // Repeating the last argument cancels only up to rounding.
                    expectClose(dotHyperCross(v0, v1, v2, v2), 0,
                        1e-12 * length(v0) * length(v1) * length(v2) ** 2, 0);
                });
        });

    it('unitHyperCross is the normalized hypercross product', () => {
        check(fc.tuple(vector(4), vector(4), vector(4), fc.boolean()),
            ([a, b, c, robust]) => {
                const h = hyperCross(a, b, c);
                const len = length(h);
                if (len < 1e-2) {
                    return;   // near-dependent inputs
                }
                const u = unitHyperCross(a, b, c, robust);
                expectClose(length(u), 1, 1e-12, 1e-12);
                expectVectorClose(mul(u, len), h, 1e-9, 1e-9);
                // Upstream does not modify its inputs.
                expect(hyperCross(a, b, c).values).toEqual(h.values);
            });
        // Linearly dependent inputs give the zero hypercross, which Normalize
        // maps to the zero vector.
        const e0 = Vector.fromArray([1, 0, 0, 0]);
        const e1 = Vector.fromArray([0, 1, 0, 0]);
        expect(unitHyperCross(e0, e1, e0).values.map(x => x + 0))
            .toEqual([0, 0, 0, 0]);
    });

    it('computeOrthogonalComplement4 extends 1, 2 or 3 inputs to a '
        + 'right-handed orthonormal basis', () => {
        // Well-scaled components: with robust = true and a subnormal largest
        // component, upstream Normalize returns NaN (see the upstream-quirk
        // test in test/Vector.test.ts), which makes the frame meaningless.
        check(fc.tuple(wellScaledVector(4), wellScaledVector(4),
            wellScaledVector(4), fc.integer({ min: 1, max: 3 }), fc.boolean()),
            ([a, b, c, numInputs, robust]) => {
                const v = [a.clone(), b.clone(), c.clone(),
                    new Vector(4)].slice(0, numInputs);
                const minLength = computeOrthogonalComplement4(numInputs, v,
                    robust);
                // '!(x >= eps)' rather than 'x < eps' so a NaN return, which
                // the degenerate branches can produce, is also skipped.
                if (!(minLength >= 1e-2)) {
                    return;   // near-dependent inputs; upstream returns ~0
                }
                expect(v.length).toBe(4);
                for (let i = 0; i < 4; ++i) {
                    expect(v[i].size).toBe(4);
                    expectClose(length(v[i]), 1, 1e-9, 1e-9);
                    for (let j = 0; j < i; ++j) {
                        expectClose(dot(v[i], v[j]), 0, 1e-9, 1e-9);
                    }
                }
                // v[3] is assigned HyperCross(v[0], v[1], v[2]) before the
                // Gram-Schmidt pass, so the frame is right handed.
                expectClose(dotHyperCross(v[0], v[1], v[2], v[3]), 1,
                    1e-9, 1e-9);
            });
    });

    it('computeOrthogonalComplement4 picks v[1] by the largest component',
        () => {
            // Upstream branches on the index of the largest |component| of
            // v[0] before any normalization: index 0 or 1 gives
            // (-x1, +x0, 0, 0), index 3 gives (0, +x2, -x1, 0) [so that
            // affine 3D vectors with w = 0 work], index 2 gives
            // (0, 0, -x3, +x2).
            const run = (a: number[]): number[] => {
                const v = [Vector.fromArray(a)];
                computeOrthogonalComplement4(1, v);
                return v[1].values;
            };
            // Largest at index 1: (-x1, +x0, 0, 0) before orthonormalizing,
            // which leaves v[1] in the e0-e1 plane.
            const v1a = run([1, 5, 2, 3]);
            expect(v1a[2]).toBe(0);
            expect(v1a[3]).toBe(0);
            // Largest at index 2: (0, 0, -x3, +x2), in the e2-e3 plane.
            const v1b = run([1, 2, 9, 3]);
            expect(v1b[0]).toBe(0);
            expect(v1b[1]).toBe(0);
            // Largest at index 3: (0, +x2, -x1, 0), in the e1-e2 plane.
            const v1c = run([1, 2, 3, 9]);
            expect(v1c[0]).toBe(0);
            expect(v1c[3]).toBe(0);
        });

    it('computeOrthogonalComplement4 handles affine 3D vectors (w = 0)',
        () => {
            // The maxIndex == 3 clause exists so that a 3-tuple written as a
            // 4-tuple with w = 0 still produces a usable basis.
            check(vector(3).filter(u => length(u) > 1e-2), u => {
                const v = [Vector.fromArray([...u.values, 0])];
                const minLength = computeOrthogonalComplement4(1, v);
                if (minLength < 1e-2) {
                    return;
                }
                for (let i = 0; i < 4; ++i) {
                    expectClose(length(v[i]), 1, 1e-9, 1e-9);
                    for (let j = 0; j < i; ++j) {
                        expectClose(dot(v[i], v[j]), 0, 1e-9, 1e-9);
                    }
                }
            });
        });

    it('computeOrthogonalComplement4 returns 0 for unsupported numInputs',
        () => {
            check(fc.integer({ min: -2, max: 6 }).filter(n => n < 1 || n > 3),
                n => {
                    const v = [Vector.fromArray([1, 0, 0, 0]),
                        Vector.fromArray([0, 1, 0, 0]),
                        Vector.fromArray([0, 0, 1, 0]),
                        Vector.fromArray([0, 0, 0, 1])];
                    expect(computeOrthogonalComplement4(n, v)).toBe(0);
                    // The inputs are untouched.
                    expect(v[0].values).toEqual([1, 0, 0, 0]);
                });
        });

    it('every entry point requires 4-tuples', () => {
        const bad = Vector.fromArray([1, 2, 3]);
        const ok = Vector.fromArray([1, 0, 0, 0]);
        const msg = 'Vector4: vector must have size 4.';
        expect(() => hyperCross(bad, ok, ok)).toThrow(msg);
        expect(() => hyperCross(ok, bad, ok)).toThrow(msg);
        expect(() => hyperCross(ok, ok, bad)).toThrow(msg);
        expect(() => unitHyperCross(bad, ok, ok)).toThrow(msg);
        expect(() => dotHyperCross(bad, ok, ok, ok)).toThrow(msg);
        expect(() => computeOrthogonalComplement4(1, [bad])).toThrow(msg);
        expect(() => computeOrthogonalComplement4(2, [ok, bad])).toThrow(msg);
    });
});
