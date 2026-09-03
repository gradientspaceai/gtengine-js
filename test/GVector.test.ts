import { describe, it, expect } from 'vitest';
import { check, finite, fc } from './helpers/arbitraries.js';
import { GVector } from '../src/GVector.js';
import {
    Vector, add, sub, mul, div, dot, length, normalize, orthonormalize, hlift
} from '../src/Vector.js';

describe('GVector construction and sizing', () => {
    it('default constructor produces the size-0 tuple', () => {
        const v = new GVector();
        expect(v.size).toBe(0);
        expect(v.values).toEqual([]);
    });

    it('GVector(size) zero-fills', () => {
        const v = new GVector(3);
        expect(v.size).toBe(3);
        expect(v.values).toEqual([0, 0, 0]);
    });

    it('GVector(size, d) makes a Euclidean basis vector', () => {
        expect(new GVector(3, 1).values).toEqual([0, 1, 0]);
        // Invalid d yields the zero vector, matching upstream MakeUnit.
        expect(new GVector(3, -1).values).toEqual([0, 0, 0]);
        expect(new GVector(3, 3).values).toEqual([0, 0, 0]);
    });

    it('negative size throws', () => {
        expect(() => new GVector(-1)).toThrow('Invalid size.');
    });

    it('is a Vector (upstream API parity via subclassing)', () => {
        const v = new GVector(2);
        expect(v).toBeInstanceOf(Vector);
    });

    it('static factories return GVector instances', () => {
        const z = GVector.zero(4);
        const u = GVector.unit(4, 2);
        const f = GVector.fromArray([1, 2, 3]);
        expect(z).toBeInstanceOf(GVector);
        expect(u).toBeInstanceOf(GVector);
        expect(f).toBeInstanceOf(GVector);
        expect(z.values).toEqual([0, 0, 0, 0]);
        expect(u.values).toEqual([0, 0, 1, 0]);
        expect(f.values).toEqual([1, 2, 3]);
    });

    it('clone returns a deep GVector copy', () => {
        const v = GVector.fromArray([5, 6]);
        const c = v.clone();
        expect(c).toBeInstanceOf(GVector);
        c.set(0, -1);
        expect(v.values).toEqual([5, 6]);
    });

    it('setSize grows with zero-initialized elements, preserving data', () => {
        const v = GVector.fromArray([1, 2]);
        v.setSize(5);
        expect(v.values).toEqual([1, 2, 0, 0, 0]);
    });

    it('setSize shrinks, preserving the prefix', () => {
        const v = GVector.fromArray([1, 2, 3, 4]);
        v.setSize(2);
        expect(v.values).toEqual([1, 2]);
        // Regrown elements are zero, not stale.
        v.setSize(3);
        expect(v.values).toEqual([1, 2, 0]);
    });

    it('setSize rejects negative sizes', () => {
        const v = new GVector(2);
        expect(() => v.setSize(-3)).toThrow('Invalid size.');
    });

    it('makeZero and makeUnit are inherited', () => {
        const v = GVector.fromArray([1, 2, 3]);
        v.makeUnit(0);
        expect(v.values).toEqual([1, 0, 0]);
        v.makeZero();
        expect(v.values).toEqual([0, 0, 0]);
    });
});

describe('GVector comparisons (std::vector semantics)', () => {
    it('equality requires equal sizes and elements', () => {
        const a = GVector.fromArray([1, 2]);
        const b = GVector.fromArray([1, 2]);
        const c = GVector.fromArray([1, 3]);
        expect(a.equals(b)).toBe(true);
        expect(a.notEquals(b)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(c)).toBe(true);
    });

    it('vectors of different sizes compare without throwing', () => {
        const short = GVector.fromArray([1, 2]);
        const long = GVector.fromArray([1, 2, 0]);
        expect(short.equals(long)).toBe(false);
        expect(short.notEquals(long)).toBe(true);
        // Lexicographic: equal prefix, shorter orders first.
        expect(short.lessThan(long)).toBe(true);
        expect(short.lessThanOrEqual(long)).toBe(true);
        expect(long.greaterThan(short)).toBe(true);
        expect(long.greaterThanOrEqual(short)).toBe(true);
    });

    it('lexicographic ordering on the common prefix dominates size', () => {
        const a = GVector.fromArray([1, 5]);
        const b = GVector.fromArray([2, 0, 0]);
        expect(a.lessThan(b)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
    });

    it('same-size ordering matches Vector semantics', () => {
        const a = GVector.fromArray([1, 2, 3]);
        const b = GVector.fromArray([1, 2, 4]);
        expect(a.lessThan(b)).toBe(true);
        expect(a.lessThanOrEqual(a)).toBe(true);
        expect(a.greaterThanOrEqual(a)).toBe(true);
        expect(a.greaterThan(b)).toBe(false);
    });
});

describe('GVector reuses the Vector module free functions', () => {
    it('algebraic operators apply to GVector arguments', () => {
        const v0 = GVector.fromArray([1, 2, 3]);
        const v1 = GVector.fromArray([4, 5, 6]);
        expect(add(v0, v1).values).toEqual([5, 7, 9]);
        expect(sub(v0, v1).values).toEqual([-3, -3, -3]);
        expect(mul(v0, 2).values).toEqual([2, 4, 6]);
        expect(div(v1, 2).values).toEqual([2, 2.5, 3]);
        expect(dot(v0, v1)).toBe(32);
    });

    it('length and normalize apply to GVector', () => {
        const v = GVector.fromArray([3, 4]);
        expect(length(v)).toBe(5);
        const len = normalize(v);
        expect(len).toBe(5);
        expect(v.values[0]).toBeCloseTo(0.6, 15);
        expect(v.values[1]).toBeCloseTo(0.8, 15);
    });

    it('orthonormalize applies to arrays of GVector', () => {
        const v = [
            GVector.fromArray([1, 1, 0]),
            GVector.fromArray([1, 0, 1]),
            GVector.fromArray([0, 1, 1])
        ];
        const minLength = orthonormalize(3, v);
        expect(minLength).toBeGreaterThan(0);
        for (let i = 0; i < 3; ++i) {
            expect(length(v[i])).toBeCloseTo(1, 12);
            for (let j = i + 1; j < 3; ++j) {
                expect(dot(v[i], v[j])).toBeCloseTo(0, 12);
            }
        }
    });

    it('hlift produces the homogeneous tuple', () => {
        const v = GVector.fromArray([1, 2]);
        expect(hlift(v, 7).values).toEqual([1, 2, 7]);
    });

    it('mismatched sizes throw in the shared operators', () => {
        const v0 = GVector.fromArray([1, 2]);
        const v1 = GVector.fromArray([1, 2, 3]);
        expect(() => add(v0, v1)).toThrow('mismatched sizes');
        expect(() => dot(v0, v1)).toThrow('mismatched sizes');
    });
});

// ---------------------------------------------------------------------------
// Verification wave: property-based checks against upstream GVector.h.
// ---------------------------------------------------------------------------

describe('GVector verification', () => {
    it('setSize grows with zeros, keeps the prefix and truncates', () => {
        // Upstream SetSize calls std::vector::resize, which preserves the
        // leading elements and value-initializes any new ones.
        check(fc.tuple(fc.array(finite(), { minLength: 0, maxLength: 8 }),
            fc.integer({ min: 0, max: 10 })), ([vals, n]) => {
            const v = GVector.fromArray(vals);
            v.setSize(n);
            expect(v.size).toBe(n);
            for (let i = 0; i < n; ++i) {
                expect(v.get(i)).toBe(i < vals.length ? vals[i] : 0);
            }
        });
        expect(() => new GVector(0).setSize(-1)).toThrow('Invalid size.');
        expect(() => new GVector(-1)).toThrow('Invalid size.');
    });

    it('comparisons follow std::vector, not the fixed-size Vector', () => {
        // std::vector's relational operators are lexicographic over the
        // common prefix with the shorter range ordered first; unlike
        // Vector's, they accept mismatched sizes instead of throwing.
        check(fc.tuple(fc.array(finite(), { minLength: 0, maxLength: 4 }),
            fc.array(finite(), { minLength: 0, maxLength: 4 })),
            ([a, b]) => {
                const u = GVector.fromArray(a);
                const w = GVector.fromArray(b);
                let expected = 0;
                const n = Math.min(a.length, b.length);
                for (let i = 0; i < n && expected === 0; ++i) {
                    expected = a[i] < b[i] ? -1 : (a[i] > b[i] ? 1 : 0);
                }
                if (expected === 0) {
                    expected = a.length < b.length ? -1
                        : (a.length > b.length ? 1 : 0);
                }
                expect(u.lessThan(w)).toBe(expected < 0);
                expect(u.greaterThan(w)).toBe(expected > 0);
                expect(u.equals(w)).toBe(expected === 0);
                expect(u.notEquals(w)).toBe(expected !== 0);
                expect(u.lessThanOrEqual(w)).toBe(expected <= 0);
                expect(u.greaterThanOrEqual(w)).toBe(expected >= 0);
            });
    });

    it('equals follows std::vector: a NaN component breaks self-equality', () => {
        // Regression: the port compared lexicographically, which treats NaN
        // as equivalent; C++ 'NaN == NaN' is false.
        const v = GVector.fromArray([1, NaN]);
        expect(v.equals(v)).toBe(false);
        expect(v.notEquals(v)).toBe(true);
        expect(v.equals(v.clone())).toBe(false);
        expect(GVector.fromArray([1, 2]).equals(GVector.fromArray([1, 2])))
            .toBe(true);
    });

    it('GVector(size, d) is the d-th standard basis vector', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 6 }),
            fc.integer({ min: -2, max: 8 })), ([n, d]) => {
            for (const v of [new GVector(n, d), GVector.unit(n, d)]) {
                expect(v.size).toBe(n);
                for (let i = 0; i < n; ++i) {
                    expect(v.get(i)).toBe(i === d ? 1 : 0);
                }
            }
        });
    });

    it('clone is deep and preserves the GVector type', () => {
        check(fc.array(finite(), { minLength: 1, maxLength: 5 }), vals => {
            const v = GVector.fromArray(vals);
            const c = v.clone();
            expect(c).toBeInstanceOf(GVector);
            expect(c.values).toEqual(v.values);
            c.set(0, c.get(0) + 1);
            expect(v.get(0)).toBe(vals[0]);
            // clone() of a resizable tuple must itself be resizable.
            c.setSize(vals.length + 2);
            expect(v.size).toBe(vals.length);
        });
    });

    it('Vector.h free functions apply unchanged to GVector tuples', () => {
        // Upstream GVector.h re-implements Dot/Length/Normalize/... with the
        // same algorithms, so the shared implementations must agree.
        check(fc.array(finite(), { minLength: 1, maxLength: 6 }), vals => {
            const g = GVector.fromArray(vals);
            const v = Vector.fromArray(vals);
            expect(dot(g, g)).toBe(dot(v, v));
            expect(length(g)).toBe(length(v));
            expect(add(g, g).values).toEqual(add(v, v).values);
            expect(sub(g, g).values).toEqual(sub(v, v).values);
            expect(mul(g, 3).values).toEqual(mul(v, 3).values);
            expect(hlift(g, 5).values).toEqual(hlift(v, 5).values);
            const gn = g.clone(), vn = v.clone();
            expect(normalize(gn)).toBe(normalize(vn));
            expect(gn.values).toEqual(vn.values);
        });
    });

    it('orthonormalize over GVector tuples yields an orthonormal set', () => {
        check(fc.array(fc.array(finite(), { minLength: 3, maxLength: 3 }),
            { minLength: 3, maxLength: 3 }), rows => {
            const v = rows.map(r => GVector.fromArray(r));
            const minLength = orthonormalize(3, v);
            if (minLength < 1e-3) {
                return;   // near-dependent inputs
            }
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(length(v[i]) - 1)).toBeLessThan(1e-9);
                for (let j = 0; j < i; ++j) {
                    expect(Math.abs(dot(v[i], v[j]))).toBeLessThan(1e-9);
                }
            }
        });
    });

    it('div of a GVector by zero yields zero (Vector.h, not GVector.h)', () => {
        // Documented deviation: upstream GVector operator/= calls
        // LogError("Division by zero."); the shared div() keeps Vector.h's
        // behavior of producing the zero vector.
        expect(div(GVector.fromArray([1, 2]), 0).values).toEqual([0, 0]);
    });
});
