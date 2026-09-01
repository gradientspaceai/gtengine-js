import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line';
import { Torus3 } from '../src/Torus3';
import { Vector, add, dot, mul, normalize, sub, length } from '../src/Vector';
import { IntrLine3Torus3FI } from '../src/IntrLine3Torus3';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function line(p: number[], d: number[]): Line {
    const dir = vec(d);
    normalize(dir);
    return Line.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The implicit torus function; it is zero on the surface.
function implicit(torus: Torus3, x: Vector): number {
    const delta = sub(x, torus.center);
    const sqrLen = dot(delta, delta);
    const dotN = dot(torus.normal, delta);
    const r0Sqr = torus.radius0 * torus.radius0;
    const r1Sqr = torus.radius1 * torus.radius1;
    const a = sqrLen + r0Sqr - r1Sqr;
    return a * a - 4 * r0Sqr * (sqrLen - dotN * dotN);
}

describe('IntrLine3Torus3', () => {
    const fi = new IntrLine3Torus3FI();

    // The standard torus: center at the origin, axis (0,0,1), r0 = 2, r1 = 1.
    const torus = new Torus3();

    it('finds four intersections along the plane of symmetry', () => {
        const l = line([0, 0, 0], [1, 0, 0]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(4);
        const expected = [-3, -1, 1, 3];
        for (let i = 0; i < 4; ++i) {
            expect(result.lineParameter[i]).toBeCloseTo(expected[i], 8);
            expect(result.point[i].values[0]).toBeCloseTo(expected[i], 8);
            expect(Math.abs(implicit(torus, result.point[i])))
                .toBeLessThan(1e-6);
        }
    });

    it('finds two intersections for a line through the tube', () => {
        // The line x = 2, y = 0 pierces the tube at z = -1 and z = 1.
        const l = line([2, 0, -5], [0, 0, 1]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.lineParameter[0]).toBeCloseTo(4, 8);
        expect(result.lineParameter[1]).toBeCloseTo(6, 8);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 8);
        expect(result.point[1].values[2]).toBeCloseTo(1, 8);
    });

    it('reports no intersection for a line through the hole', () => {
        const l = line([0, 0, -5], [0, 0, 1]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection for a line above the torus', () => {
        const l = line([0, 0, 5], [1, 0, 0]);
        expect(fi.find(l, torus).intersect).toBe(false);
    });

    it('reports the two tangential contacts in the plane z = 1', () => {
        // At z = 1 the torus meets the plane in the double circle
        // x^2 + y^2 = 4, so each contact is a double root.
        const l = line([0, 0, 1], [1, 0, 0]);
        const result = fi.find(l, torus);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.lineParameter[0]).toBeCloseTo(-2, 6);
        expect(result.lineParameter[1]).toBeCloseTo(2, 6);
    });

    it('reports torus parameters consistent with the surface evaluation', () => {
        const l = line([0, 0, 0], [1, 0.35, 0.2]);
        const result = fi.find(l, torus);
        expect(result.numIntersections).toBeGreaterThan(0);
        for (let i = 0; i < result.numIntersections; ++i) {
            const [u, v] = result.torusParameter[i];
            const jet = torus.evaluate(u, v, 0);
            expect(length(sub(jet[0], result.point[i]))).toBeLessThan(1e-6);
        }
    });

    it('handles a translated and reoriented torus', () => {
        // The torus has axis (0,1,0) and center (1,-2,3). A line through the
        // center along the axis passes through the hole and misses.
        const t = Torus3.fromCenterFrameRadii(
            vec([1, -2, 3]), vec([0, 0, 1]), vec([1, 0, 0]), vec([0, 1, 0]),
            2, 0.5);
        const miss = line([1, -10, 3], [0, 1, 0]);
        expect(fi.find(miss, t).intersect).toBe(false);

        // A line in the plane of symmetry crosses all four times.
        const hit = line([1, -2, 3], [0, 0, 1]);
        const result = fi.find(hit, t);
        expect(result.numIntersections).toBe(4);
        const expected = [-2.5, -1.5, 1.5, 2.5];
        for (let i = 0; i < 4; ++i) {
            expect(result.lineParameter[i]).toBeCloseTo(expected[i], 8);
            expect(Math.abs(implicit(t, result.point[i]))).toBeLessThan(1e-6);
        }
    });

    it('produces points on the torus for random lines', () => {
        const rand = makeRandom(271828);
        let total = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const l = line(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(l, torus);
            expect(result.numIntersections).toBeLessThanOrEqual(4);
            total += result.numIntersections;
            for (let i = 0; i < result.numIntersections; ++i) {
                // The point is on the line.
                const onLine = add(l.origin,
                    mul(result.lineParameter[i], l.direction));
                expect(length(sub(onLine, result.point[i]))).toBeLessThan(1e-9);
                // The point is on the torus.
                expect(Math.abs(implicit(torus, result.point[i])))
                    .toBeLessThan(1e-5);
            }
            // The roots are reported in increasing order.
            for (let i = 1; i < result.numIntersections; ++i) {
                expect(result.lineParameter[i])
                    .toBeGreaterThan(result.lineParameter[i - 1]);
            }
        }
        expect(total).toBeGreaterThan(50);
    });
});
