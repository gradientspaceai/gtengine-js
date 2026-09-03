import { describe, it, expect } from 'vitest';
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
