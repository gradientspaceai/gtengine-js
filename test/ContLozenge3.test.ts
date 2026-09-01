import { describe, it, expect } from 'vitest';
import {
    getContainerLozenge3,
    inContainerLozenge3
} from '../src/ContLozenge3';
import { Lozenge3 } from '../src/Lozenge3';
import { Rectangle } from '../src/Rectangle';
import { DistPointRectangle } from '../src/DistPointRectangle';
import { Vector, dot } from '../src/Vector';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function distanceToRectangle(point: Vector, rectangle: Rectangle): number {
    return new DistPointRectangle().compute(point, rectangle).distance;
}

function makeRandom(seed: number): () => number {
    let s = seed;
    return () => {
        s = (s * 1103515245 + 12345) % 2147483648;
        return s / 2147483648;
    };
}

describe('getContainerLozenge3', () => {
    it('contains every input point (random slab-like clouds)', () => {
        const rand = makeRandom(31337);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 50; ++i) {
                points.push(v(
                    6 * (rand() - 0.5),
                    4 * (rand() - 0.5),
                    1 * (rand() - 0.5)));
            }
            const lozenge = getContainerLozenge3(points);
            expect(lozenge.radius).toBeGreaterThan(0);
            for (const p of points) {
                expect(distanceToRectangle(p, lozenge.rectangle))
                    .toBeLessThanOrEqual(lozenge.radius + 1e-9);
            }
        }
    });

    it('contains every input point of a random isotropic cloud', () => {
        const rand = makeRandom(90210);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                points.push(v(
                    4 * (rand() - 0.5),
                    4 * (rand() - 0.5),
                    4 * (rand() - 0.5)));
            }
            const lozenge = getContainerLozenge3(points);
            for (const p of points) {
                // Extreme points land on the lozenge boundary, where the
                // strict upstream test can fail by an ulp; allow rounding.
                expect(distanceToRectangle(p, lozenge.rectangle))
                    .toBeLessThanOrEqual(lozenge.radius + 1e-12);
            }
            // Points pulled toward the center are strictly inside.
            for (const p of points) {
                const shrunk = v(
                    lozenge.rectangle.center.values[0]
                        + 0.5 * (p.values[0] - lozenge.rectangle.center.values[0]),
                    lozenge.rectangle.center.values[1]
                        + 0.5 * (p.values[1] - lozenge.rectangle.center.values[1]),
                    lozenge.rectangle.center.values[2]
                        + 0.5 * (p.values[2] - lozenge.rectangle.center.values[2]));
                expect(inContainerLozenge3(shrunk, lozenge)).toBe(true);
            }
        }
    });

    it('fits a planar point set with a near-zero radius', () => {
        // All points in the z = 0 plane, so the thin direction has zero
        // spread and the lozenge radius is 0.
        const points: Vector[] = [];
        for (let i = -3; i <= 3; ++i) {
            for (let j = -2; j <= 2; ++j) {
                points.push(v(i, j, 0));
            }
        }
        const lozenge = getContainerLozenge3(points);
        expect(lozenge.radius).toBeCloseTo(0, 10);
        // The rectangle plane is z = 0.
        expect(lozenge.rectangle.center.values[2]).toBeCloseTo(0, 10);
        expect(Math.abs(lozenge.rectangle.axis[0].values[2])).toBeCloseTo(0, 10);
        expect(Math.abs(lozenge.rectangle.axis[1].values[2])).toBeCloseTo(0, 10);
        // The rectangle covers the 6 x 4 extent of the data.
        const extents = [
            lozenge.rectangle.extent.values[0],
            lozenge.rectangle.extent.values[1]
        ].sort((a, b) => a - b);
        expect(extents[0]).toBeCloseTo(2, 8);
        expect(extents[1]).toBeCloseTo(3, 8);
        for (const p of points) {
            expect(distanceToRectangle(p, lozenge.rectangle))
                .toBeLessThanOrEqual(lozenge.radius + 1e-9);
        }
    });

    it('degenerates to a sphere-like lozenge for coincident points', () => {
        const p = v(-2, 3, 1);
        const lozenge = getContainerLozenge3([p, p.clone(), p.clone(), p.clone()]);
        expect(lozenge.radius).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.extent.values[0]).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.extent.values[1]).toBeCloseTo(0, 12);
        expect(lozenge.rectangle.center.values[0]).toBeCloseTo(-2, 12);
        expect(lozenge.rectangle.center.values[1]).toBeCloseTo(3, 12);
        expect(lozenge.rectangle.center.values[2]).toBeCloseTo(1, 12);
    });

    it('handles collinear points (a capsule-shaped lozenge)', () => {
        const points = [v(0, 0, 0), v(1, 0, 0), v(2, 0, 0), v(3, 0, 0)];
        const lozenge = getContainerLozenge3(points);
        expect(lozenge.radius).toBeCloseTo(0, 10);
        const extents = [
            lozenge.rectangle.extent.values[0],
            lozenge.rectangle.extent.values[1]
        ].sort((a, b) => a - b);
        expect(extents[0]).toBeCloseTo(0, 10);
        expect(extents[1]).toBeCloseTo(1.5, 10);
        for (const p of points) {
            expect(distanceToRectangle(p, lozenge.rectangle))
                .toBeLessThanOrEqual(1e-9);
        }
    });

    it('keeps the rectangle axes orthonormal', () => {
        const rand = makeRandom(1717);
        const points: Vector[] = [];
        for (let i = 0; i < 30; ++i) {
            points.push(v(5 * (rand() - 0.5), 3 * (rand() - 0.5),
                1.5 * (rand() - 0.5)));
        }
        const lozenge = getContainerLozenge3(points);
        const a0 = lozenge.rectangle.axis[0];
        const a1 = lozenge.rectangle.axis[1];
        expect(dot(a0, a0)).toBeCloseTo(1, 10);
        expect(dot(a1, a1)).toBeCloseTo(1, 10);
        expect(dot(a0, a1)).toBeCloseTo(0, 10);
    });

    it('is translation covariant', () => {
        const rand = makeRandom(606);
        const base: Vector[] = [];
        for (let i = 0; i < 25; ++i) {
            base.push(v(4 * (rand() - 0.5), 3 * (rand() - 0.5),
                (rand() - 0.5)));
        }
        const shifted = base.map(p => v(p.values[0] + 10, p.values[1] - 5,
            p.values[2] + 2));
        const l0 = getContainerLozenge3(base);
        const l1 = getContainerLozenge3(shifted);
        expect(l1.radius).toBeCloseTo(l0.radius, 8);
        expect(l1.rectangle.center.values[0])
            .toBeCloseTo(l0.rectangle.center.values[0] + 10, 8);
        expect(l1.rectangle.center.values[1])
            .toBeCloseTo(l0.rectangle.center.values[1] - 5, 8);
        expect(l1.rectangle.center.values[2])
            .toBeCloseTo(l0.rectangle.center.values[2] + 2, 8);
    });

    it('throws for an empty point set and for non-3D points', () => {
        expect(() => getContainerLozenge3([])).toThrow();
        expect(() => getContainerLozenge3([Vector.fromArray([1, 2])])).toThrow();
    });
});

describe('inContainerLozenge3', () => {
    const lozenge = Lozenge3.fromRectangleRadius(
        Rectangle.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0)], Vector.fromArray([2, 1])), 1);

    it('accepts points inside the slab over the rectangle', () => {
        expect(inContainerLozenge3(v(0, 0, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(2, 1, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 0, 1), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 0, -1), lozenge)).toBe(true);
    });

    it('accepts points in the rounded side and corner regions', () => {
        expect(inContainerLozenge3(v(3, 0, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(v(0, 2, 0), lozenge)).toBe(true);
        expect(inContainerLozenge3(
            v(2 + Math.SQRT1_2, 1 + Math.SQRT1_2, 0), lozenge)).toBe(true);
    });

    it('rejects points beyond the radius', () => {
        expect(inContainerLozenge3(v(0, 0, 1.0001), lozenge)).toBe(false);
        expect(inContainerLozenge3(v(3.0001, 0, 0), lozenge)).toBe(false);
        expect(inContainerLozenge3(v(3, 2, 0), lozenge)).toBe(false);
    });
});
