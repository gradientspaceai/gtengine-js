import { describe, it, expect } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { Vector, dot, sub, add, mul, length } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import { check, compareKeys, expectClose, expectStrictWeakOrder, fc, finite,
    positive, unitVector, vector } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

describe('Circle3 construction', () => {
    it('the default constructor is the unit circle in the x-y plane', () => {
        const circle = new Circle3();
        expect(circle.center.values).toEqual([0, 0, 0]);
        expect(circle.normal.values).toEqual([0, 0, 1]);
        expect(circle.radius).toBe(1);
    });

    it('fromCenterNormalRadius copies the inputs', () => {
        const center = v3(1, 2, 3);
        const normal = v3(0, 1, 0);
        const circle = Circle3.fromCenterNormalRadius(center, normal, 4);
        center.set(0, 99);
        normal.set(1, 99);
        expect(circle.center.values).toEqual([1, 2, 3]);
        expect(circle.normal.values).toEqual([0, 1, 0]);
        expect(circle.radius).toBe(4);
    });

    it('rejects vectors that are not 3D', () => {
        expect(() => Circle3.fromCenterNormalRadius(Vector.fromArray([0, 0]),
            v3(0, 0, 1), 1)).toThrow();
    });

    it('clone is a deep copy', () => {
        const circle = new Circle3();
        const copy = circle.clone();
        copy.center.set(2, 5);
        copy.radius = 9;
        expect(circle.center.values).toEqual([0, 0, 0]);
        expect(circle.radius).toBe(1);
    });
});

describe('Circle3 geometry (the defining constraints)', () => {
    it('points of a circle satisfy |X-C| = r and Dot(N, X-C) = 0', () => {
        // Circle with center (1,2,3), normal (0,1,0), radius 5. Its plane is
        // spanned by (1,0,0) and (0,0,1).
        const circle = Circle3.fromCenterNormalRadius(v3(1, 2, 3),
            v3(0, 1, 0), 5);
        for (let k = 0; k < 16; ++k) {
            const t = (k / 16) * 2 * Math.PI;
            const x = v3(1 + 5 * Math.cos(t), 2, 3 + 5 * Math.sin(t));
            const diff = sub(x, circle.center);
            expect(length(diff)).toBeCloseTo(circle.radius, 12);
            expect(dot(circle.normal, diff)).toBeCloseTo(0, 12);
        }
    });
});

describe('Circle3 comparisons', () => {
    const base = new Circle3();

    it('equals compares center, normal and radius', () => {
        expect(base.equals(new Circle3())).toBe(true);
        expect(base.notEquals(new Circle3())).toBe(false);

        const other = base.clone();
        other.normal = v3(1, 0, 0);
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);
    });

    it('lessThan orders by center, then normal, then radius', () => {
        const smallCenter = base.clone();
        smallCenter.center = v3(-1, 0, 0);
        smallCenter.radius = 100;
        expect(smallCenter.lessThan(base)).toBe(true);

        const smallNormal = base.clone();
        smallNormal.normal = v3(0, 0, -1);
        smallNormal.radius = 100;
        expect(smallNormal.lessThan(base)).toBe(true);

        const smallRadius = base.clone();
        smallRadius.radius = 0.5;
        expect(smallRadius.lessThan(base)).toBe(true);
        expect(base.lessThan(smallRadius)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const bigger = base.clone();
        bigger.radius = 2;
        expect(base.lessThanOrEqual(bigger)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(bigger.greaterThan(base)).toBe(true);
        expect(bigger.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
    });
});

describe('Circle3 verification', () => {
    const circle = () => fc.tuple(vector(3, -5, 5), unitVector(3), positive(5))
        .map(([c, n, r]) => Circle3.fromCenterNormalRadius(c, n, r));
    const key = (c: Circle3) =>
        [...c.center.values, ...c.normal.values, c.radius];

    it('the members describe the intersection of a sphere and a plane', () => {
        // X = C + r*(cos(t)*U + sin(t)*V) for an orthonormal {U,V,N} lies on
        // the sphere |X-C| = r and in the plane Dot(N, X-C) = 0.
        check(fc.tuple(circle(), finite(-Math.PI, Math.PI)), ([c, t]) => {
            const basis = [c.normal.clone(), new Vector(3), new Vector(3)];
            computeOrthogonalComplement3(1, basis);
            const x = add(c.center, add(mul(c.radius * Math.cos(t), basis[1]),
                mul(c.radius * Math.sin(t), basis[2])));
            const delta = sub(x, c.center);
            expectClose(length(delta), c.radius, 1e-12, 1e-12);
            expectClose(dot(c.normal, delta), 0, 1e-12, 1e-12);
        });
    });

    it('the comparisons follow the (center, normal, radius) member order',
        () => {
            check(fc.tuple(circle(), circle()), ([a, b]) => {
                const cmp = compareKeys(key(a), key(b));
                expect(a.lessThan(b)).toBe(cmp < 0);
                expect(a.greaterThan(b)).toBe(cmp > 0);
                expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
                expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
                expect(a.equals(b)).toBe(cmp === 0);
            });
        });

    it('the normal is compared before the radius', () => {
        const c = Vector.fromArray([0, 0, 0]);
        const a = Circle3.fromCenterNormalRadius(c,
            Vector.fromArray([0, 0, 1]), 9);
        const b = Circle3.fromCenterNormalRadius(c,
            Vector.fromArray([1, 0, 0]), 1);
        expect(a.lessThan(b)).toBe(true);
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(circle(), { minLength: 4, maxLength: 5 }), cs => {
            expectStrictWeakOrder(cs, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('equals is element equality, so NaN members break self-equality', () => {
        const c = Circle3.fromCenterNormalRadius(Vector.fromArray([0, 0, 0]),
            Vector.fromArray([NaN, 0, 1]), 1);
        expect(c.equals(c)).toBe(false);
        expect(c.lessThan(c)).toBe(false);
        expect(c.lessThanOrEqual(c)).toBe(true);
    });

    it('the factory and clone are independent of their inputs', () => {
        check(fc.tuple(vector(3), unitVector(3)), ([c, n]) => {
            const circle3 = Circle3.fromCenterNormalRadius(c, n, 1);
            const cloned = circle3.clone();
            c.set(0, 999);
            n.set(1, 888);
            circle3.center.set(2, 777);
            expect(circle3.center.get(0)).not.toBe(999);
            expect(circle3.normal.get(1)).not.toBe(888);
            expect(cloned.center.get(2)).not.toBe(777);
        });
    });
});
