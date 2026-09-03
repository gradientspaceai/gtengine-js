import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { IntrLine3Cylinder3FI } from '../src/IntrLine3Cylinder3.js';
import {
    intrLine3Cylinder3FIDoQuery,
    defaultIntrLine3Cylinder3FIResult
} from '../src/IntrLine3Cylinder3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const dir = vec(direction);
    normalize(dir);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(vec(origin), dir), radius, height);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function insideCylinder(c: Cylinder3, x: Vector): boolean {
    const diff = sub(x, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    if (Math.abs(z) > 0.5 * c.height) {
        return false;
    }
    const radial = sub(diff, mul(z, c.axis.direction));
    return dot(radial, radial) <= c.radius * c.radius;
}

describe('IntrLine3Cylinder3', () => {
    const fi = new IntrLine3Cylinder3FI();

    // Cylinder about the z-axis, radius 1, height 4 (z in [-2,2]).
    const cyl = cylinder([0, 0, 0], [0, 0, 1], 1, 4);

    it('finds the chord of a line perpendicular to the axis', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('handles a line parallel to the cylinder axis', () => {
        const l = line([0.5, 0, -10], [0, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(8, 12);
        expect(result.parameter[1]).toBeCloseTo(12, 12);
        expect(result.point[0].values[2]).toBeCloseTo(-2, 12);
        expect(result.point[1].values[2]).toBeCloseTo(2, 12);
    });

    it('reports no intersection for a parallel line outside the radius', () => {
        const l = line([2, 0, -10], [0, 0, 1]);
        expect(fi.find(l, cyl).intersect).toBe(false);
    });

    it('reports no intersection for a perpendicular line beyond the end disk', () => {
        const l = line([0, 0, 3], [1, 0, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds the chord of an oblique line that exits through the wall', () => {
        const l = line([0, 0, 0], [1, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-Math.SQRT2, 10);
        expect(result.parameter[1]).toBeCloseTo(Math.SQRT2, 10);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 10);
        expect(result.point[1].values[0]).toBeCloseTo(1, 10);
    });

    it('finds an oblique line entering and leaving through the end disks', () => {
        const l = line([0, 0, -10], [0.01, 0, 1]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[2]).toBeCloseTo(-2, 10);
        expect(result.point[1].values[2]).toBeCloseTo(2, 10);
    });

    it('reports a single point for a line tangent to the wall', () => {
        const l = line([1, -10, 0], [0, 1, 0]);
        const result = fi.find(l, cyl);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(10, 10);
        expect(result.parameter[1]).toBeCloseTo(result.parameter[0], 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 10);
    });

    it('handles a zero-height cylinder as a disk', () => {
        const flat = cylinder([0, 0, 0], [0, 0, 1], 1, 0);
        const l = line([0.5, 0, -3], [0, 0, 1]);
        const result = fi.find(l, flat);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
    });

    it('does not support an infinite cylinder', () => {
        // Cylinder3 stores height = -1 for the infinite state and this
        // upstream query has no infinite branch, so it computes with
        // halfHeight = -0.5 and finds nothing. The behavior is preserved.
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, 1);
        infinite.makeInfiniteCylinder();
        expect(infinite.isInfinite()).toBe(true);
        const l = line([0, 0, 0], [1, 0, 0]);
        expect(fi.find(l, infinite).intersect).toBe(false);
    });

    it('agrees with dense sampling along the line', () => {
        const rand = makeRandom(112358);
        const cyl2 = cylinder([0.5, -1, 0.25], [1, 2, -1], 0.9, 3);
        for (let trial = 0; trial < 120; ++trial) {
            const l = line(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, cyl2);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = -10 + (20 * k) / n;
                const x = add(l.origin, mul(t, l.direction));
                if (insideCylinder(cyl2, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }

            if (result.intersect) {
                // Every reported point lies in the closed cylinder.
                for (let i = 0; i < 2; ++i) {
                    const p = add(l.origin, mul(result.parameter[i], l.direction));
                    const diff = sub(p, cyl2.axis.origin);
                    const z = dot(diff, cyl2.axis.direction);
                    const radial = sub(diff, mul(z, cyl2.axis.direction));
                    expect(Math.abs(z)).toBeLessThanOrEqual(1.5 + 1e-9);
                    expect(Math.sqrt(dot(radial, radial)))
                        .toBeLessThanOrEqual(0.9 + 1e-9);
                }
            }
        }
    });
});

describe('intrLine3Cylinder3FIDoQuery', () => {
    const c = cylinder([0, 0, 0], [0, 0, 1], 1, 2);

    it('matches the class query but does not compute points', () => {
        const l = line([-5, 0, 0], [1, 0, 0]);
        const result = defaultIntrLine3Cylinder3FIResult();
        intrLine3Cylinder3FIDoQuery(l.origin, l.direction, c, result);
        const expected = new IntrLine3Cylinder3FI().find(l, c);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter[0]).toBeCloseTo(expected.parameter[0], 12);
        expect(result.parameter[1]).toBeCloseTo(expected.parameter[1], 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('reports no intersection for a line missing the cylinder', () => {
        const result = defaultIntrLine3Cylinder3FIResult();
        intrLine3Cylinder3FIDoQuery(vec([-5, 3, 0]), vec([1, 0, 0]), c,
            result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});
