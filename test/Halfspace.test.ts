import { describe, it, expect } from 'vitest';
import { Halfspace } from '../src/Halfspace.js';
import { Vector, dot } from '../src/Vector.js';

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
