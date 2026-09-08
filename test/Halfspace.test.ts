import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Vector, dot, mul } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    unitVector, vector } from './helpers/arbitraries.js';

describe('Halfspace construction', () => {
    it('the default constructor is x[N-1] >= 0', () => {
        const halfspace = new Halfspace(3);
        expect(halfspace.dimension).toBe(3);
        expect(halfspace.normal.values).toEqual([0, 0, 1]);
        expect(halfspace.constant).toBe(0);
    });

    it('the default constructor works in 2D as well', () => {
        const halfspace = new Halfspace(2);
        expect(halfspace.normal.values).toEqual([0, 1]);
        expect(halfspace.constant).toBe(0);
    });

    it('fromNormalConstant copies the normal', () => {
        const normal = Vector.fromArray([1, 0, 0]);
        const halfspace = Halfspace.fromNormalConstant(normal, 5);
        normal.set(0, 99);
        expect(halfspace.normal.values).toEqual([1, 0, 0]);
        expect(halfspace.constant).toBe(5);
    });

    it('clone is a deep copy', () => {
        const halfspace = Halfspace.fromNormalConstant(
            Vector.fromArray([0, 1]), 2);
        const copy = halfspace.clone();
        copy.normal.set(0, 7);
        copy.constant = -1;
        expect(halfspace.normal.values).toEqual([0, 1]);
        expect(halfspace.constant).toBe(2);
    });
});

describe('Halfspace membership predicate', () => {
    // The halfspace is Dot(N,X) >= c with unit-length N.
    it('classifies points by Dot(N,X) - c', () => {
        const invSqrt2 = 1 / Math.sqrt(2);
        const halfspace = Halfspace.fromNormalConstant(
            Vector.fromArray([invSqrt2, invSqrt2]), invSqrt2);
        // The boundary line is x + y = 1.
        const onBoundary = Vector.fromArray([1, 0]);
        const inside = Vector.fromArray([2, 2]);
        const outside = Vector.fromArray([0, 0]);
        expect(dot(halfspace.normal, onBoundary) - halfspace.constant)
            .toBeCloseTo(0, 15);
        expect(dot(halfspace.normal, inside) - halfspace.constant)
            .toBeGreaterThan(0);
        expect(dot(halfspace.normal, outside) - halfspace.constant)
            .toBeLessThan(0);
    });
});

describe('Halfspace comparisons', () => {
    const a = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), 1);
    const sameAsA = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), 1);
    // Larger normal.
    const b = Halfspace.fromNormalConstant(Vector.fromArray([1, 0]), 1);
    // Same normal as a, larger constant.
    const c = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), 2);

    it('equals/notEquals compare normal and constant', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders by normal first, then by constant', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(c.greaterThan(a)).toBe(true);
    });
});

describe('Halfspace verification', () => {
    const halfspace = (n: number) =>
        fc.tuple(unitVector(n), finite(-5, 5))
            .map(([nrm, c]) => Halfspace.fromNormalConstant(nrm, c));
    const key = (h: Halfspace) => [...h.normal.values, h.constant];

    it('the sign convention is Dot(N,X) >= c', () => {
        // X = c*N + t*N lies in the halfspace exactly when t >= 0, because
        // Dot(N,X) = c + t for a unit-length N.
        check(fc.tuple(halfspace(3), finite(-5, 5)), ([h, t]) => {
            const x = mul(h.constant + t, h.normal);
            const value = dot(h.normal, x) - h.constant;
            // |N| = 1 up to rounding, so the value is t to within rounding of
            // the dot product of a unit vector with a vector of length |c+t|.
            expectClose(value, t, 1e-12, 1e-12);
        });
    });

    it('the comparisons follow the (normal, constant) member order', () => {
        check(fc.tuple(halfspace(3), halfspace(3)), ([a, b]) => {
            const c = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            expect(a.equals(b)).toBe(c === 0);
        });
    });

    it('a differing constant alone orders by the constant', () => {
        check(fc.tuple(unitVector(3), finite(-5, 5), finite(-5, 5)),
            ([nrm, c0, c1]) => {
                const a = Halfspace.fromNormalConstant(nrm, c0);
                const b = Halfspace.fromNormalConstant(nrm, c1);
                expect(a.lessThan(b)).toBe(c0 < c1);
            });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(halfspace(2), { minLength: 4, maxLength: 6 }), hs => {
            expectStrictWeakOrder(hs, (x, y) => x.lessThan(y));
        }, 50);
    });

    it('equals is element equality, so NaN members break self-equality', () => {
        const h = Halfspace.fromNormalConstant(Vector.fromArray([NaN, 1]), 0);
        expect(h.equals(h)).toBe(false);
        expect(h.lessThan(h)).toBe(false);
        const g = Halfspace.fromNormalConstant(Vector.fromArray([0, 1]), NaN);
        expect(g.equals(g)).toBe(false);
        // The upstream ordering chain ends with 'constant < constant', which
        // is false for NaN, so the lexicographic operators still order it as
        // equivalent to itself.
        expect(g.lessThan(g)).toBe(false);
        expect(g.lessThanOrEqual(g)).toBe(true);
    });

    it('fromNormalConstant and clone are independent of their inputs', () => {
        check(vector(3), nrm => {
            const h = Halfspace.fromNormalConstant(nrm, 1);
            const copy = h.clone();
            nrm.set(0, 999);
            h.normal.set(1, 777);
            expect(h.normal.get(0)).not.toBe(999);
            expect(copy.normal.get(1)).not.toBe(777);
        });
    });
});
