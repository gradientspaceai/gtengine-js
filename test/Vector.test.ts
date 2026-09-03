import { describe, it, expect } from 'vitest';
import {
    check, finite, nonzero, vector, expectClose, expectVectorClose, fc
} from './helpers/arbitraries.js';
import {
    Vector, negate, add, sub, mul, div, compMul, compDiv,
    dot, length, normalize, orthonormalize, getOrthogonal, computeExtremes,
    hlift, hproject, lift, project
} from '../src/Vector.js';

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

function randomVector(rng: () => number, n: number): Vector {
    const v = new Vector(n);
    for (let i = 0; i < n; ++i) {
        v.values[i] = (rng() - 0.5) * 20;
    }
    return v;
}

describe('Vector construction and member access', () => {
    it('new Vector(n) is the zero vector of size n', () => {
        const v = new Vector(3);
        expect(v.size).toBe(3);
        expect(v.size).toBe(3);
        expect(v.values).toEqual([0, 0, 0]);
    });

    it('fromArray copies the input', () => {
        const source = [1, 2, 3];
        const v = Vector.fromArray(source);
        source[0] = 99;
        expect(v.values).toEqual([1, 2, 3]);
    });

    it('filled, zero, ones and unit produce the special vectors', () => {
        expect(Vector.filled(2, 7).values).toEqual([7, 7]);
        expect(Vector.zero(4).values).toEqual([0, 0, 0, 0]);
        expect(Vector.ones(3).values).toEqual([1, 1, 1]);
        expect(Vector.unit(3, 1).values).toEqual([0, 1, 0]);
        // Invalid d yields the zero vector, matching upstream MakeUnit.
        expect(Vector.unit(3, -1).values).toEqual([0, 0, 0]);
        expect(Vector.unit(3, 3).values).toEqual([0, 0, 0]);
    });

    it('get/set access components', () => {
        const v = new Vector(2);
        v.set(0, 5);
        v.set(1, -3);
        expect(v.get(0)).toBe(5);
        expect(v.get(1)).toBe(-3);
    });

    it('makeZero, makeOnes and makeUnit mutate in place', () => {
        const v = Vector.fromArray([4, 5, 6]);
        v.makeOnes();
        expect(v.values).toEqual([1, 1, 1]);
        v.makeUnit(2);
        expect(v.values).toEqual([0, 0, 1]);
        v.makeZero();
        expect(v.values).toEqual([0, 0, 0]);
    });

    it('clone is a deep copy', () => {
        const v = Vector.fromArray([1, 2]);
        const w = v.clone();
        w.set(0, 42);
        expect(v.get(0)).toBe(1);
    });
});

describe('Vector comparisons', () => {
    it('equality', () => {
        const a = Vector.fromArray([1, 2, 3]);
        expect(a.equals(Vector.fromArray([1, 2, 3]))).toBe(true);
        expect(a.notEquals(Vector.fromArray([1, 2, 3]))).toBe(false);
        expect(a.equals(Vector.fromArray([1, 2, 4]))).toBe(false);
        expect(a.notEquals(Vector.fromArray([1, 2, 4]))).toBe(true);
    });

    it('lexicographic ordering matches std::array semantics', () => {
        const a = Vector.fromArray([1, 2, 3]);
        const b = Vector.fromArray([1, 3, 0]);
        expect(a.lessThan(b)).toBe(true);
        expect(a.lessThanOrEqual(b)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(b.greaterThanOrEqual(a)).toBe(true);
        expect(a.lessThan(a.clone())).toBe(false);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);
    });

    it('comparing mismatched sizes throws', () => {
        const a = Vector.fromArray([1, 2]);
        const b = Vector.fromArray([1, 2, 3]);
        expect(() => a.equals(b)).toThrow('Vector: mismatched sizes.');
        expect(() => a.lessThan(b)).toThrow('Vector: mismatched sizes.');
    });
});

describe('Vector algebraic operations', () => {
    const v0 = Vector.fromArray([1, -2, 3]);
    const v1 = Vector.fromArray([4, 5, -6]);

    it('negate', () => {
        expect(negate(v0).values).toEqual([-1, 2, -3]);
    });

    it('add and sub', () => {
        expect(add(v0, v1).values).toEqual([5, 3, -3]);
        expect(sub(v0, v1).values).toEqual([-3, -7, 9]);
        // Inputs are not mutated.
        expect(v0.values).toEqual([1, -2, 3]);
    });

    it('mul accepts both argument orders', () => {
        expect(mul(v0, 2).values).toEqual([2, -4, 6]);
        expect(mul(2, v0).values).toEqual([2, -4, 6]);
    });

    it('div multiplies by the reciprocal; zero divisor yields the zero vector', () => {
        expect(div(Vector.fromArray([2, 4, -8]), 2).values).toEqual([1, 2, -4]);
        expect(div(v0, 0).values).toEqual([0, 0, 0]);
    });

    it('componentwise mul and div', () => {
        expect(compMul(v0, v1).values).toEqual([4, -10, -18]);
        expect(compDiv(Vector.fromArray([4, 9, -8]), Vector.fromArray([2, 3, 4])).values)
            .toEqual([2, 3, -2]);
    });

    it('mismatched sizes throw', () => {
        const small = Vector.fromArray([1, 2]);
        expect(() => add(v0, small)).toThrow('Vector: mismatched sizes.');
        expect(() => dot(v0, small)).toThrow('Vector: mismatched sizes.');
    });

    it('algebraic identities on random vectors', () => {
        const rng = makeRng(0xC0FFEE);
        for (let trial = 0; trial < 100; ++trial) {
            const n = 2 + Math.floor(rng() * 5);
            const a = randomVector(rng, n);
            const b = randomVector(rng, n);
            // a + b == b + a
            expect(add(a, b).equals(add(b, a))).toBe(true);
            // a - b == a + (-b)
            expect(sub(a, b).equals(add(a, negate(b)))).toBe(true);
            // (a + b) . c linearity through dot symmetry
            expect(dot(a, b)).toBeCloseTo(dot(b, a), 12);
        }
    });
});

describe('Vector geometric operations', () => {
    it('dot of known vectors', () => {
        expect(dot(Vector.fromArray([1, 2, 3]), Vector.fromArray([4, -5, 6]))).toBe(12);
        expect(dot(Vector.fromArray([1, 0]), Vector.fromArray([0, 1]))).toBe(0);
    });

    it('length of a 3-4-5 triangle', () => {
        expect(length(Vector.fromArray([3, 4]))).toBe(5);
        expect(length(Vector.fromArray([3, 4]), true)).toBeCloseTo(5, 14);
    });

    it('robust length avoids overflow where the standard path overflows', () => {
        const huge = Vector.fromArray([3e300, 4e300]);
        expect(length(huge)).toBe(Infinity);
        expect(length(huge, true) / 5e300).toBeCloseTo(1, 14);
        expect(length(Vector.zero(3), true)).toBe(0);
    });

    it('normalize mutates in place and returns the original length', () => {
        const v = Vector.fromArray([3, 4]);
        const len = normalize(v);
        expect(len).toBe(5);
        expect(v.values[0]).toBeCloseTo(0.6, 15);
        expect(v.values[1]).toBeCloseTo(0.8, 15);
        expect(length(v)).toBeCloseTo(1, 15);
    });

    it('normalize of the zero vector produces the zero vector and length 0', () => {
        const v = Vector.zero(3);
        expect(normalize(v)).toBe(0);
        expect(v.values).toEqual([0, 0, 0]);
        const w = Vector.zero(3);
        expect(normalize(w, true)).toBe(0);
        expect(w.values).toEqual([0, 0, 0]);
    });

    it('robust normalize handles vectors whose squared length underflows', () => {
        const tiny = Vector.fromArray([3e-200, 4e-200]);
        const standard = tiny.clone();
        expect(normalize(standard)).toBe(0); // dot underflows to 0
        const robust = tiny.clone();
        const len = normalize(robust, true);
        expect(len / 5e-200).toBeCloseTo(1, 14);
        expect(length(robust)).toBeCloseTo(1, 14);
    });

    it('orthonormalize produces an orthonormal set', () => {
        const rng = makeRng(0xBEEF);
        for (let trial = 0; trial < 20; ++trial) {
            const v = [randomVector(rng, 3), randomVector(rng, 3), randomVector(rng, 3)];
            const minLength = orthonormalize(3, v);
            expect(minLength).toBeGreaterThan(0);
            for (let i = 0; i < 3; ++i) {
                for (let j = 0; j < 3; ++j) {
                    expect(dot(v[i], v[j])).toBeCloseTo(i === j ? 1 : 0, 10);
                }
            }
        }
    });

    it('orthonormalize returns 0 for invalid inputs', () => {
        const v = [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])];
        expect(orthonormalize(0, v)).toBe(0);
        expect(orthonormalize(3, v)).toBe(0); // numInputs > dimension
    });

    it('getOrthogonal returns a vector orthogonal to the input', () => {
        const rng = makeRng(0xACE);
        for (let trial = 0; trial < 50; ++trial) {
            const n = 2 + Math.floor(rng() * 4);
            const v = randomVector(rng, n);
            const u = getOrthogonal(v, false);
            expect(dot(u, v)).toBeCloseTo(0, 10);
            const w = getOrthogonal(v, true);
            expect(dot(w, v)).toBeCloseTo(0, 10);
            expect(length(w)).toBeCloseTo(1, 12);
        }
    });

    it('getOrthogonal picks the maximum-magnitude component', () => {
        const u = getOrthogonal(Vector.fromArray([1, 5, 2]), false);
        // imax = 1, inext = 2: u[1] = v[2], u[2] = -v[1].
        expect(u.values).toEqual([0, 2, -5]);
    });

    it('computeExtremes finds the axis-aligned bounding box', () => {
        const result = computeExtremes([
            Vector.fromArray([1, 5]),
            Vector.fromArray([-2, 7]),
            Vector.fromArray([3, 6])
        ]);
        expect(result).not.toBeNull();
        expect(result!.vmin.values).toEqual([-2, 5]);
        expect(result!.vmax.values).toEqual([3, 7]);
    });

    it('computeExtremes of an empty set is null', () => {
        expect(computeExtremes([])).toBeNull();
    });
});

describe('Vector lift and project', () => {
    const v = Vector.fromArray([1, 2, 3]);

    it('hlift appends the last component', () => {
        expect(hlift(v, 9).values).toEqual([1, 2, 3, 9]);
    });

    it('hproject drops the last component', () => {
        expect(hproject(Vector.fromArray([1, 2, 3, 9])).values).toEqual([1, 2, 3]);
        expect(() => hproject(Vector.fromArray([1]))).toThrow('Invalid dimension.');
    });

    it('hproject inverts hlift', () => {
        expect(hproject(hlift(v, 5)).equals(v)).toBe(true);
    });

    it('lift injects a value at the given index', () => {
        expect(lift(v, 0, 9).values).toEqual([9, 1, 2, 3]);
        expect(lift(v, 1, 9).values).toEqual([1, 9, 2, 3]);
        expect(lift(v, 3, 9).values).toEqual([1, 2, 3, 9]);
    });

    it('project rejects the component at the given index', () => {
        const u = Vector.fromArray([1, 2, 3, 4]);
        expect(project(u, 0).values).toEqual([2, 3, 4]);
        expect(project(u, 2).values).toEqual([1, 2, 4]);
        expect(project(u, 3).values).toEqual([1, 2, 3]);
        expect(() => project(Vector.fromArray([1]), 0)).toThrow('Invalid dimension.');
    });

    it('project inverts lift for every inject index', () => {
        for (let inject = 0; inject <= v.size; ++inject) {
            expect(project(lift(v, inject, 7), inject).equals(v)).toBe(true);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream Vector.h.
// ---------------------------------------------------------------------------

describe('Vector verification', () => {
    it('dot is symmetric and bilinear', () => {
        check(fc.tuple(vector(4), vector(4), vector(4), finite()),
            ([u, v, w, s]) => {
                expect(dot(u, v)).toBe(dot(v, u));   // exact: same products
                // dot(u, v + s*w) = dot(u, v) + s*dot(u, w)
                expectClose(dot(u, add(v, mul(w, s))),
                    dot(u, v) + s * dot(u, w), 1e-9, 1e-9);
            });
    });

    it('length equals sqrt(dot(v, v)) and the robust path agrees', () => {
        check(vector(5), v => {
            expect(length(v)).toBe(Math.sqrt(dot(v, v)));
            // The robust path rescales by the largest |component| before
            // squaring, so it differs only by rounding of that scaling.
            expectClose(length(v, true), length(v, false), 1e-12, 1e-12);
        });
    });

    it('normalize returns the input length and leaves a unit vector', () => {
        check(fc.tuple(vector(4).filter(v => length(v) > 1e-3), fc.boolean()),
            ([v, robust]) => {
                const original = v.clone();
                const len = normalize(v, robust);
                expectClose(len, length(original), 1e-12, 1e-12);
                expectClose(length(v), 1, 1e-12, 1e-12);
                // v * len reproduces the input.
                expectVectorClose(mul(v, len), original, 1e-9, 1e-9);
            });
    });

    it('normalize of the zero vector yields the zero vector and length 0', () => {
        for (const robust of [false, true]) {
            const v = new Vector(3);
            expect(normalize(v, robust)).toBe(0);
            expect(v.values).toEqual([0, 0, 0]);
        }
    });

    it('div by a scalar equals multiplication by its reciprocal', () => {
        check(fc.tuple(vector(4), nonzero()), ([v, s]) => {
            // Upstream operator/= multiplies by 1/scalar, so this is exact.
            expect(div(v, s).values).toEqual(mul(v, 1 / s).values);
        });
        // Upstream operator/= zeroes the vector when the scalar is zero.
        expect(div(Vector.fromArray([1, 2, 3]), 0).values).toEqual([0, 0, 0]);
    });

    it('add/sub/negate satisfy the vector-space identities', () => {
        check(fc.tuple(vector(4), vector(4)), ([u, v]) => {
            expectVectorClose(sub(add(u, v), v), u, 1e-9, 1e-9);
            expect(add(u, negate(u)).values).toEqual(new Vector(4).values);
            expect(sub(u, v).values).toEqual(add(u, negate(v)).values);
        });
    });

    it('compMul/compDiv are componentwise inverses', () => {
        check(fc.tuple(vector(4),
            fc.array(nonzero(), { minLength: 4, maxLength: 4 })), ([u, ds]) => {
            const d = Vector.fromArray(ds);
            expectVectorClose(compDiv(compMul(u, d), d), u, 1e-9, 1e-9);
        });
    });

    it('orthonormalize produces an orthonormal set and the minimum length', () => {
        check(fc.tuple(fc.array(vector(4), { minLength: 4, maxLength: 4 }),
            fc.boolean()), ([raw, robust]) => {
            const v = raw.map(x => x.clone());
            const minLength = orthonormalize(4, v, robust);
            if (minLength < 1e-3) {
                return;   // near-dependent inputs; the result is ill-conditioned
            }
            for (let i = 0; i < 4; ++i) {
                expectClose(length(v[i]), 1, 1e-9, 1e-9);
                for (let j = 0; j < i; ++j) {
                    expectClose(dot(v[i], v[j]), 0, 1e-9, 1e-9);
                }
            }
            // v[0] keeps the direction of the first input.
            const u0 = raw[0].clone();
            normalize(u0, robust);
            expectVectorClose(v[0], u0, 1e-9, 1e-9);
        });
    });

    it('orthonormalize returns 0 for invalid numInputs', () => {
        const v = [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0])];
        expect(orthonormalize(0, v)).toBe(0);
        expect(orthonormalize(4, v)).toBe(0);
    });

    it('getOrthogonal is exactly orthogonal to its input', () => {
        check(vector(5), v => {
            // u has exactly two nonzero components, u[imax] = v[inext] and
            // u[inext] = -v[imax], so dot(u, v) = v[imax]*v[inext] -
            // v[inext]*v[imax] cancels bit-exactly.
            expect(dot(getOrthogonal(v, false), v)).toBe(0);
        });
    });

    it('getOrthogonal with unitLength returns a unit orthogonal vector', () => {
        check(vector(5).filter(v => length(v) > 1e-2), v => {
            const u = getOrthogonal(v, true);
            expectClose(length(u), 1, 1e-12, 1e-12);
            // Scaling the two nonzero components by 1/sqrt(sqrDistance)
            // rounds them independently, so the cancellation is no longer
            // exact; the residual is bounded by the rounding of |v|.
            expectClose(dot(u, v), 0, 1e-12 * length(v), 0);
        });
    });

    it('getOrthogonal picks the index of the largest absolute component', () => {
        // Upstream: u[imax] = v[inext], u[inext] = -v[imax] with inext =
        // (imax + 1) mod N; ties keep the first (strict > comparison).
        // v = (1,2,3): imax = 2, inext = 0, so u[2] = v[0] and u[0] = -v[2].
        expect(getOrthogonal(Vector.fromArray([1, 2, 3]), false).values)
            .toEqual([-3, 0, 1]);
        expect(getOrthogonal(Vector.fromArray([3, 2, 1]), false).values)
            .toEqual([2, -3, 0]);
        // imax = N-1 wraps to inext = 0.
        expect(getOrthogonal(Vector.fromArray([1, 2, 5]), false).values)
            .toEqual([-5, 0, 1]);
        // A tie keeps imax = 0.
        expect(getOrthogonal(Vector.fromArray([2, 2, 0]), false).values)
            .toEqual([2, -2, 0]);
    });

    it('computeExtremes agrees with a brute-force componentwise min/max', () => {
        check(fc.array(vector(3), { minLength: 1, maxLength: 12 }), vs => {
            const res = computeExtremes(vs)!;
            for (let i = 0; i < 3; ++i) {
                let lo = vs[0].get(i), hi = vs[0].get(i);
                for (const v of vs) {
                    lo = Math.min(lo, v.get(i));
                    hi = Math.max(hi, v.get(i));
                }
                expect(res.vmin.get(i)).toBe(lo);
                expect(res.vmax.get(i)).toBe(hi);
            }
        });
        expect(computeExtremes([])).toBeNull();
    });

    it('computeExtremes does not alias its inputs', () => {
        const a = Vector.fromArray([1, 2]);
        const b = Vector.fromArray([3, 4]);
        const res = computeExtremes([a, b])!;
        res.vmin.set(0, 99);
        res.vmax.set(0, -99);
        expect(a.values).toEqual([1, 2]);
        expect(res.vmin.get(0)).toBe(99);
        expect(res.vmax.get(0)).toBe(-99);
    });

    it('hproject inverts hlift and lift/project round trip at every index', () => {
        check(fc.tuple(vector(4), finite()), ([v, last]) => {
            const h = hlift(v, last);
            expect(h.size).toBe(5);
            expect(h.get(4)).toBe(last);
            expect(hproject(h).values).toEqual(v.values);
            for (let k = 0; k <= v.size; ++k) {
                const l = lift(v, k, last);
                expect(l.size).toBe(5);
                expect(l.get(k)).toBe(last);
                expect(project(l, k).values).toEqual(v.values);
            }
        });
    });

    it('comparisons are the lexicographic std::array ordering', () => {
        check(fc.tuple(vector(3), vector(3)), ([u, v]) => {
            const lt = u.lessThan(v), gt = u.greaterThan(v), eq = u.equals(v);
            // Exactly one of <, >, == holds for finite components.
            expect([lt, gt, eq].filter(Boolean).length).toBe(1);
            expect(u.lessThanOrEqual(v)).toBe(lt || eq);
            expect(u.greaterThanOrEqual(v)).toBe(gt || eq);
            expect(u.notEquals(v)).toBe(!eq);
            // Reference lexicographic comparison.
            let expected = 0;
            for (let i = 0; i < 3 && expected === 0; ++i) {
                expected = u.get(i) < v.get(i) ? -1
                    : (u.get(i) > v.get(i) ? 1 : 0);
            }
            expect(lt ? -1 : (gt ? 1 : 0)).toBe(expected);
        });
    });

    it('equals follows std::array: a NaN component is not equal to itself', () => {
        // Regression: the port compared lexicographically, which reports NaN
        // components as "equivalent"; C++ 'NaN == NaN' is false.
        const v = Vector.fromArray([1, NaN, 3]);
        expect(v.equals(v)).toBe(false);
        expect(v.notEquals(v)).toBe(true);
        expect(v.equals(v.clone())).toBe(false);
        // The ordering operators keep the lexicographic behavior, where NaN
        // is neither less than nor greater than anything.
        expect(v.lessThan(v)).toBe(false);
        expect(v.greaterThan(v)).toBe(false);
        // Ordinary vectors are unaffected.
        const w = Vector.fromArray([1, 2, 3]);
        expect(w.equals(w.clone())).toBe(true);
    });

    it('mismatched sizes throw for the operations that require agreement', () => {
        const a = Vector.fromArray([1, 2]);
        const b = Vector.fromArray([1, 2, 3]);
        expect(() => dot(a, b)).toThrow('mismatched sizes');
        expect(() => add(a, b)).toThrow('mismatched sizes');
        expect(() => sub(a, b)).toThrow('mismatched sizes');
        expect(() => compMul(a, b)).toThrow('mismatched sizes');
        expect(() => compDiv(a, b)).toThrow('mismatched sizes');
        expect(() => a.equals(b)).toThrow('mismatched sizes');
        expect(() => a.lessThan(b)).toThrow('mismatched sizes');
    });

    it('the special vectors match MakeZero/MakeOnes/MakeUnit', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 6 }),
            fc.integer({ min: -2, max: 8 })), ([n, d]) => {
            expect(Vector.zero(n).values).toEqual(new Array<number>(n).fill(0));
            expect(Vector.ones(n).values).toEqual(new Array<number>(n).fill(1));
            const u = Vector.unit(n, d);
            for (let i = 0; i < n; ++i) {
                expect(u.get(i)).toBe(i === d ? 1 : 0);
            }
        });
    });
});
