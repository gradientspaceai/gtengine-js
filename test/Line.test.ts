import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Vector, add, sub, mul, length } from '../src/Vector.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    line, unitVector, vector } from './helpers/arbitraries.js';

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

describe('Line verification', () => {
    const key = (l: Line) => [...l.origin.values, ...l.direction.values];

    it('the parameterisation P + t*D reproduces the origin at t = 0 and is '
        + 'affine in t', () => {
            check(fc.tuple(line(3), finite(-5, 5), finite(-5, 5)),
                ([l, s, t]) => {
                    const p0 = add(l.origin, mul(s, l.direction));
                    const p1 = add(l.origin, mul(t, l.direction));
                    // |P(s) - P(t)| = |s - t| for a unit-length direction.
                    expectClose(length(sub(p0, p1)), Math.abs(s - t),
                        1e-12, 1e-12);
                    // The point at t = 0 is the origin, exactly.
                    expect(add(l.origin, mul(0, l.direction)).values)
                        .toEqual(l.origin.values);
                });
        });

    it('the comparisons follow the (origin, direction) member order', () => {
        check(fc.tuple(line(3), line(3)), ([a, b]) => {
            const c = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(c < 0);
            expect(a.greaterThan(b)).toBe(c > 0);
            expect(a.lessThanOrEqual(b)).toBe(c <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(c >= 0);
            expect(a.equals(b)).toBe(c === 0);
        });
    });

    it('a differing direction alone orders by the direction', () => {
        check(fc.tuple(vector(3), unitVector(3), unitVector(3)),
            ([o, d0, d1]) => {
                const a = Line.fromOriginDirection(o, d0);
                const b = Line.fromOriginDirection(o, d1);
                expect(a.lessThan(b)).toBe(d0.lessThan(d1));
            });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(line(2), { minLength: 4, maxLength: 6 }), ls => {
            expectStrictWeakOrder(ls, (x, y) => x.lessThan(y));
        }, 50);
    });

    it('equals is element equality, so a NaN line does not equal itself',
        () => {
            const l = Line.fromOriginDirection(Vector.fromArray([0, 0]),
                Vector.fromArray([NaN, 1]));
            expect(l.equals(l)).toBe(false);
            expect(l.lessThan(l)).toBe(false);
            expect(l.lessThanOrEqual(l)).toBe(true);
        });

    it('fromOriginDirection and clone are independent of their inputs', () => {
        check(fc.tuple(vector(3), unitVector(3)), ([o, d]) => {
            const l = Line.fromOriginDirection(o, d);
            const copy = l.clone();
            o.set(0, 999);
            d.set(1, 888);
            l.origin.set(2, 777);
            expect(l.origin.get(0)).not.toBe(999);
            expect(l.direction.get(1)).not.toBe(888);
            expect(copy.origin.get(2)).not.toBe(777);
        });
    });
});
