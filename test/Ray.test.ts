import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Vector, add, sub, mul, length } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, positive,
    ray, unitVector, vector } from './helpers/arbitraries.js';

describe('Ray construction', () => {
    it('the default constructor is the positive x-axis from the origin', () => {
        const ray = new Ray(3);
        expect(ray.dimension).toBe(3);
        expect(ray.origin.values).toEqual([0, 0, 0]);
        expect(ray.direction.values).toEqual([1, 0, 0]);
    });

    it('fromOriginDirection copies the input vectors', () => {
        const origin = Vector.fromArray([1, 2]);
        const direction = Vector.fromArray([0, 1]);
        const ray = Ray.fromOriginDirection(origin, direction);
        origin.set(0, 99);
        direction.set(1, 99);
        expect(ray.origin.values).toEqual([1, 2]);
        expect(ray.direction.values).toEqual([0, 1]);
    });

    it('fromOriginDirection throws on mismatched sizes', () => {
        expect(() => Ray.fromOriginDirection(new Vector(3), new Vector(2)))
            .toThrow('Ray: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const ray = Ray.fromOriginDirection(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([0, 0, 1]));
        const copy = ray.clone();
        copy.direction.set(2, -1);
        expect(ray.direction.values).toEqual([0, 0, 1]);
        expect(copy.direction.values).toEqual([0, 0, -1]);
    });
});

describe('Ray parameterization', () => {
    it('P + t*D with t >= 0 reaches hand-computed points', () => {
        const ray = Ray.fromOriginDirection(Vector.fromArray([2, 0]),
            Vector.fromArray([0.6, 0.8]));
        const at5 = add(ray.origin, mul(ray.direction, 5));
        expect(at5.get(0)).toBeCloseTo(5, 12);
        expect(at5.get(1)).toBeCloseTo(4, 12);
        const at0 = add(ray.origin, mul(ray.direction, 0));
        expect(at0.equals(ray.origin)).toBe(true);
    });
});

describe('Ray comparisons', () => {
    const a = Ray.fromOriginDirection(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    const sameAsA = Ray.fromOriginDirection(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    // Larger origin.
    const b = Ray.fromOriginDirection(Vector.fromArray([1, 0]),
        Vector.fromArray([1, 0]));
    // Same origin as a, larger direction.
    const c = Ray.fromOriginDirection(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 1]));

    it('equals/notEquals compare origin and direction', () => {
        expect(a.equals(sameAsA)).toBe(true);
        expect(a.notEquals(sameAsA)).toBe(false);
        expect(a.equals(c)).toBe(false);
        expect(a.notEquals(b)).toBe(true);
    });

    it('orders by origin first, then by direction', () => {
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

describe('Ray verification', () => {
    const key = (r: Ray) => [...r.origin.values, ...r.direction.values];

    it('the parameterisation P + t*D has unit speed for t >= 0', () => {
        check(fc.tuple(ray(3), positive(5), positive(5)), ([r, s, t]) => {
            const p0 = add(r.origin, mul(s, r.direction));
            const p1 = add(r.origin, mul(t, r.direction));
            expectClose(length(sub(p0, r.origin)), s, 1e-12, 1e-12);
            expectClose(length(sub(p0, p1)), Math.abs(s - t), 1e-12, 1e-12);
        });
    });

    it('the comparisons follow the (origin, direction) member order', () => {
        check(fc.tuple(ray(3), ray(3)), ([a, b]) => {
            const c = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            expect(a.equals(b)).toBe(c === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(ray(2), { minLength: 4, maxLength: 6 }), rs => {
            expectStrictWeakOrder(rs, (x, y) => x.lessThan(y));
        }, 50);
    });

    it('equals is element equality, so a NaN ray does not equal itself', () => {
        const r = Ray.fromOriginDirection(Vector.fromArray([0, 0]),
            Vector.fromArray([NaN, 1]));
        expect(r.equals(r)).toBe(false);
        expect(r.lessThan(r)).toBe(false);
        expect(r.lessThanOrEqual(r)).toBe(true);
    });

    it('fromOriginDirection and clone are independent of their inputs', () => {
        check(fc.tuple(vector(3), unitVector(3)), ([o, d]) => {
            const r = Ray.fromOriginDirection(o, d);
            const copy = r.clone();
            o.set(0, 999);
            d.set(1, 888);
            r.origin.set(2, 777);
            expect(r.origin.get(0)).not.toBe(999);
            expect(r.direction.get(1)).not.toBe(888);
            expect(copy.origin.get(2)).not.toBe(777);
        });
    });
});
