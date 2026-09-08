import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, add, sub, mul, length } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, positive,
    sphere, unitVector, vector } from './helpers/arbitraries.js';

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

describe('Hypersphere verification', () => {
    const key = (h: Hypersphere) => [...h.center.values, h.radius];

    it('points at parameter distance r are on the hypersphere', () => {
        check(fc.tuple(sphere(3), unitVector(3)), ([s, u]) => {
            const x = add(s.center, mul(s.radius, u));
            expectClose(length(sub(x, s.center)), s.radius, 1e-12, 1e-12);
        });
    });

    it('the comparisons follow the (center, radius) member order', () => {
        check(fc.tuple(sphere(3), sphere(3)), ([a, b]) => {
            const c = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            expect(a.equals(b)).toBe(c === 0);
        });
    });

    it('a differing radius alone orders by the radius', () => {
        check(fc.tuple(vector(3), positive(5), positive(5)),
            ([c, r0, r1]) => {
                const a = Hypersphere.fromCenterRadius(c, r0);
                const b = Hypersphere.fromCenterRadius(c, r1);
                expect(a.lessThan(b)).toBe(r0 < r1);
            });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(sphere(2), { minLength: 4, maxLength: 6 }), ss => {
            expectStrictWeakOrder(ss, (x, y) => x.lessThan(y));
        }, 50);
    });

    it('equals is element equality, so NaN members break self-equality', () => {
        const s = Hypersphere.fromCenterRadius(Vector.fromArray([NaN, 0]), 1);
        expect(s.equals(s)).toBe(false);
        expect(s.lessThan(s)).toBe(false);
        expect(s.lessThanOrEqual(s)).toBe(true);
    });

    it('fromCenterRadius and clone are independent of their inputs', () => {
        check(vector(3), c => {
            const s = Hypersphere.fromCenterRadius(c, 2);
            const copy = s.clone();
            c.set(0, 999);
            s.center.set(1, 777);
            expect(s.center.get(0)).not.toBe(999);
            expect(copy.center.get(1)).not.toBe(777);
        });
    });
});
