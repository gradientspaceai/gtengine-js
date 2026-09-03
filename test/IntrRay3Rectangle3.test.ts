import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { IntrLine3Rectangle3FI } from '../src/IntrLine3Rectangle3.js';
import {
    IntrRay3Rectangle3TI,
    IntrRay3Rectangle3FI
} from '../src/IntrRay3Rectangle3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The rectangle in the plane z = 0, centered at the origin, with extents
// 2 (along x) and 1 (along y).
const rect = Rectangle.fromCenterAxisExtent(
    vec([0, 0, 0]), [vec([1, 0, 0]), vec([0, 1, 0])], vec([2, 1]));

describe('IntrRay3Rectangle3', () => {
    const ti = new IntrRay3Rectangle3TI();
    const fi = new IntrRay3Rectangle3FI();

    it('finds the crossing of a ray through the rectangle interior', () => {
        const r = ray([0.5, -0.25, 3], [0, 0, -1]);
        expect(ti.test(r, rect).intersect).toBe(true);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.5, 12);
        expect(result.point.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports no intersection when the rectangle is behind the ray', () => {
        const r = ray([0.5, -0.25, 3], [0, 0, 1]);
        expect(ti.test(r, rect).intersect).toBe(false);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.rectCoord).toEqual([0, 0, 0]);
    });

    it('accepts a ray whose origin is on the rectangle', () => {
        const r = ray([1, 0.5, 0], [0, 0, 1]);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.rectCoord[0]).toBeCloseTo(1, 12);
        expect(result.rectCoord[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no intersection just outside an edge', () => {
        const r = ray([2.0001, 0, 3], [0, 0, -1]);
        expect(ti.test(r, rect).intersect).toBe(false);
        expect(fi.find(r, rect).intersect).toBe(false);
    });

    it('accepts a ray hitting a corner exactly', () => {
        const r = ray([2, 1, 3], [0, 0, -1]);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.rectCoord[0]).toBeCloseTo(2, 12);
        expect(result.rectCoord[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the ray is parallel to the plane', () => {
        // Upstream calls a coplanar ray a "no intersection" case.
        const r = ray([-5, 0, 0], [1, 0, 0]);
        expect(ti.test(r, rect).intersect).toBe(false);
        expect(fi.find(r, rect).intersect).toBe(false);
    });

    it('is the line query restricted to t >= 0', () => {
        const rand = makeRandom(5150);
        const lineFI = new IntrLine3Rectangle3FI();
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, rect);
            const rayResult = fi.find(r, rect);
            expect(ti.test(r, rect).intersect).toBe(rayResult.intersect);

            if (lineResult.intersect && lineResult.parameter >= 0) {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter).toBeCloseTo(lineResult.parameter, 12);
                expect(rayResult.rectCoord[0]).toBeCloseTo(lineResult.rectCoord[0], 12);
                expect(rayResult.rectCoord[1]).toBeCloseTo(lineResult.rectCoord[1], 12);
            }
            else {
                expect(rayResult.intersect).toBe(false);
            }
        }
    });

    it('produces points that lie on the ray and inside the rectangle', () => {
        const rand = makeRandom(60606);
        // A tilted rectangle to exercise the general case.
        const a0 = vec([1, 1, 0]);
        normalize(a0);
        const a1 = vec([-1, 1, 1]);
        normalize(a1);
        // Make a1 orthogonal to a0 (it already is, up to rounding).
        const tilted = Rectangle.fromCenterAxisExtent(
            vec([0.5, -1, 2]), [a0, a1], vec([1.5, 0.75]));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, tilted);
            if (!result.intersect) {
                continue;
            }
            ++hits;
            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let k = 0; k < 3; ++k) {
                expect(result.point.values[k]).toBeCloseTo(onRay.values[k], 9);
            }
            const diff = sub(result.point, tilted.center);
            expect(Math.abs(dot(diff, tilted.axis[0])))
                .toBeLessThanOrEqual(tilted.extent.values[0] + 1e-9);
            expect(Math.abs(dot(diff, tilted.axis[1])))
                .toBeLessThanOrEqual(tilted.extent.values[1] + 1e-9);
            expect(dot(diff, tilted.axis[0]))
                .toBeCloseTo(result.rectCoord[0], 9);
            expect(dot(diff, tilted.axis[1]))
                .toBeCloseTo(result.rectCoord[1], 9);
        }
        expect(hits).toBeGreaterThan(5);
    });
});
