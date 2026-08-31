import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere';
import { Vector, sub, length } from '../src/Vector';

describe('Hypersphere construction', () => {
    it('the default constructor is the unit hypersphere at the origin', () => {
        const sphere = new Hypersphere(3);
        expect(sphere.dimension).toBe(3);
        expect(sphere.center.values).toEqual([0, 0, 0]);
        expect(sphere.radius).toBe(1);
    });

    it('fromCenterRadius copies the center', () => {
        const center = Vector.fromArray([1, 2]);
        const circle = Hypersphere.fromCenterRadius(center, 3);
        center.set(0, 99);
        expect(circle.center.values).toEqual([1, 2]);
        expect(circle.radius).toBe(3);
        expect(circle.dimension).toBe(2);
    });

    it('clone is a deep copy', () => {
        const sphere = Hypersphere.fromCenterRadius(
            Vector.fromArray([1, 1, 1]), 2);
        const copy = sphere.clone();
        copy.center.set(2, -5);
        copy.radius = 10;
        expect(sphere.center.values).toEqual([1, 1, 1]);
        expect(sphere.radius).toBe(2);
    });
});

describe('Hypersphere membership predicate', () => {
    it('|X - C| = R holds for points on the circle', () => {
        const circle = Hypersphere.fromCenterRadius(
            Vector.fromArray([1, -2]), 5);
        // (1,-2) + 5*(3/5, 4/5) = (4, 2) is on the circle.
        const onCircle = Vector.fromArray([4, 2]);
        expect(length(sub(onCircle, circle.center))).toBeCloseTo(5, 12);
        const inside = Vector.fromArray([1, -2]);
        expect(length(sub(inside, circle.center))).toBeLessThan(circle.radius);
    });
});

describe('Hypersphere comparisons', () => {
    const a = Hypersphere.fromCenterRadius(Vector.fromArray([0, 0]), 1);
    const sameAsA = Hypersphere.fromCenterRadius(Vector.fromArray([0, 0]), 1);
    // Larger center.
    const b = Hypersphere.fromCenterRadius(Vector.fromArray([0, 1]), 1);
    // Same center as a, larger radius.
    const c = Hypersphere.fromCenterRadius(Vector.fromArray([0, 0]), 2);

    it('equals/notEquals compare center and radius', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders by center first, then by radius', () => {
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThanOrEqual(sameAsA)).toBe(true);
        expect(a.greaterThanOrEqual(sameAsA)).toBe(true);
        expect(b.greaterThan(a)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
    });
});
