import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Vector, add, mul } from '../src/Vector.js';

describe('Line construction', () => {
    it('the default constructor is the x-axis through the origin', () => {
        const line = new Line(3);
        expect(line.dimension).toBe(3);
        expect(line.origin.values).toEqual([0, 0, 0]);
        expect(line.direction.values).toEqual([1, 0, 0]);
    });

    it('fromOriginDirection copies the input vectors', () => {
        const origin = Vector.fromArray([1, 2]);
        const direction = Vector.fromArray([0, 1]);
        const line = Line.fromOriginDirection(origin, direction);
        origin.set(0, 99);
        direction.set(1, 99);
        expect(line.origin.values).toEqual([1, 2]);
        expect(line.direction.values).toEqual([0, 1]);
    });

    it('fromOriginDirection throws on mismatched sizes', () => {
        expect(() => Line.fromOriginDirection(new Vector(2), new Vector(3)))
            .toThrow('Line: mismatched sizes.');
    });

    it('clone is a deep copy', () => {
        const line = Line.fromOriginDirection(Vector.fromArray([1, 1, 1]),
            Vector.fromArray([0, 0, 1]));
        const copy = line.clone();
        copy.origin.set(0, -3);
        expect(line.origin.values).toEqual([1, 1, 1]);
        expect(copy.origin.values).toEqual([-3, 1, 1]);
    });
});

describe('Line parameterization', () => {
    it('P + t*D reaches hand-computed points', () => {
        const line = Line.fromOriginDirection(Vector.fromArray([1, 2, 3]),
            Vector.fromArray([0, 1, 0]));
        const point = add(line.origin, mul(line.direction, 4));
        expect(point.values).toEqual([1, 6, 3]);
        const behind = add(line.origin, mul(line.direction, -2));
        expect(behind.values).toEqual([1, 0, 3]);
    });
});

describe('Line comparisons', () => {
    const a = Line.fromOriginDirection(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    const sameAsA = Line.fromOriginDirection(Vector.fromArray([0, 0]),
        Vector.fromArray([1, 0]));
    // Larger origin.
    const b = Line.fromOriginDirection(Vector.fromArray([0, 1]),
        Vector.fromArray([1, 0]));
    // Same origin as a, larger direction.
    const c = Line.fromOriginDirection(Vector.fromArray([0, 0]),
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
        expect(c.greaterThan(a)).toBe(true);
        expect(a.greaterThan(sameAsA)).toBe(false);
    });
});
