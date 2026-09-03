import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Vector, add, mul } from '../src/Vector.js';

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
