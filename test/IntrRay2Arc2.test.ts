import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2.js';
import { Line } from '../src/Line.js';
import {
    IntrRay2Arc2TI,
    IntrRay2Arc2FI,
    defaultIntrRay2Arc2FIResult
} from '../src/IntrRay2Arc2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, length, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function ray(origin: number[], direction: number[]): Ray {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Ray.fromOriginDirection(Vector.fromArray(origin), d);
}

// The arc of the unit circle in the first quadrant, from (1,0) to (0,1).
function firstQuadrantArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

const ti = new IntrRay2Arc2TI();
const fi = new IntrRay2Arc2FI();

describe('IntrRay2Arc2', () => {
    it('has an empty default result', () => {
        const result = defaultIntrRay2Arc2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds the single crossing of an arc at a known point', () => {
        // A ray from the origin along (1,1) meets the arc at
        // (sqrt(2)/2, sqrt(2)/2) at parameter 1.
        const result = fi.find(ray([0, 0], [1, 1]), firstQuadrantArc());
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(ti.test(ray([0, 0], [1, 1]), firstQuadrantArc()).intersect)
            .toBe(true);
    });

    it('discards circle crossings that are off the arc', () => {
        // A ray along -x from (2,0) crosses the circle at (1,0) and (-1,0).
        // Only (1,0) is on the first-quadrant arc.
        const R = ray([2, 0], [-1, 0]);
        const A = firstQuadrantArc();
        const result = fi.find(R, A);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('finds two crossings when both are on the arc', () => {
        // The arc from (0,-1) to (0,1) through (1,0) is the right half of the
        // circle. A ray along +x from (0.5,-2) does not hit it; a ray along
        // +y from (0.5,-2) hits it twice.
        const A = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(0, -1),
            vec(0, 1));
        const result = fi.find(ray([0.5, -2], [0, 1]), A);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        expect(result.point[0].values[1]).toBeCloseTo(-y, 12);
        expect(result.point[1].values[1]).toBeCloseTo(y, 12);
    });

    it('reports no intersection for a ray pointing away from the arc', () => {
        const R = ray([2, 0], [1, 0]);
        const A = firstQuadrantArc();
        expect(fi.find(R, A).intersect).toBe(false);
        expect(ti.test(R, A).intersect).toBe(false);
        // A ray that misses the circle entirely.
        const M = ray([0, 5], [1, 0]);
        expect(fi.find(M, A).intersect).toBe(false);
        expect(ti.test(M, A).intersect).toBe(false);
    });

    it('agrees with the circle query filtered by the arc on random inputs', () => {
        let state = 8080808;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        const lcQuery = new IntrLine2Circle2FI();
        let numOne = 0;
        let numTwo = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const center = vec(rand(), rand());
            const radius = 0.5 + Math.abs(rand()) * 2;
            const a0 = rand() * Math.PI;
            const a1 = a0 + 0.3 + Math.abs(rand()) * 2;
            const A = Arc2.fromCenterRadiusEnds(center, radius,
                add(center, vec(radius * Math.cos(a0), radius * Math.sin(a0))),
                add(center, vec(radius * Math.cos(a1), radius * Math.sin(a1))));
            const R = ray([rand() * 4, rand() * 4], [rand(), rand() + 0.001]);

            const result = fi.find(R, A);
            expect(ti.test(R, A).intersect).toBe(result.intersect);

            // The reference computes the ray-vs-circular-curve intersections
            // directly: the line-circle crossings with nonnegative parameter.
            const circleResult = lcQuery.find(
                Line.fromOriginDirection(R.origin, R.direction),
                Hypersphere.fromCenterRadius(center, radius));
            let expected = 0;
            for (let i = 0; i < circleResult.numIntersections; ++i) {
                if (circleResult.parameter[i] >= 0
                    && A.containsOnCircle(circleResult.point[i])) {
                    ++expected;
                }
            }
            expect(result.numIntersections).toBe(expected);
            if (expected === 1) {
                ++numOne;
            }
            else if (expected === 2) {
                ++numTwo;
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                // The point is on the circle, on the ray and on the arc.
                expect(length(sub(result.point[i], center)))
                    .toBeCloseTo(radius, 8);
                expect(result.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                const onRay = add(R.origin,
                    mul(result.parameter[i], R.direction));
                expect(length(sub(result.point[i], onRay)))
                    .toBeCloseTo(0, 10);
                expect(A.containsOnCircle(result.point[i])).toBe(true);
            }
        }
        expect(numOne).toBeGreaterThan(20);
        expect(numTwo).toBeGreaterThan(1);
    });
});
