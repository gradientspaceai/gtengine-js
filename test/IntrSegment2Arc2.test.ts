import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { Hypersphere } from '../src/Hypersphere';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2';
import {
    IntrSegment2Arc2TI,
    IntrSegment2Arc2FI,
    defaultIntrSegment2Arc2FIResult
} from '../src/IntrSegment2Arc2';
import { Line } from '../src/Line';
import { Segment } from '../src/Segment';
import { Vector, add, length, mul, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// The arc of the unit circle in the first quadrant, from (1,0) to (0,1).
function firstQuadrantArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

const ti = new IntrSegment2Arc2TI();
const fi = new IntrSegment2Arc2FI();

describe('IntrSegment2Arc2', () => {
    it('has an empty default result', () => {
        const result = defaultIntrSegment2Arc2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('finds a crossing at a known point', () => {
        // The segment from (0,0) to (2,2) crosses the arc at
        // (sqrt(2)/2, sqrt(2)/2). The centered form has origin (1,1), unit
        // direction (1,1)/sqrt(2) and extent sqrt(2), so the parameter is
        // 1 - sqrt(2).
        const result = fi.find(seg([0, 0], [2, 2]), firstQuadrantArc());
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(Math.SQRT1_2, 12);
        expect(result.parameter[0]).toBeCloseTo(1 - Math.SQRT2, 12);
        expect(ti.test(seg([0, 0], [2, 2]), firstQuadrantArc()).intersect)
            .toBe(true);
    });

    it('discards crossings outside the segment', () => {
        // The full line through (2,0) and (-2,0) crosses the circle at
        // (1,0) and (-1,0), but the segment from (2,0) to (1.5,0) reaches
        // neither.
        const S = seg([2, 0], [1.5, 0]);
        const A = firstQuadrantArc();
        expect(fi.find(S, A).intersect).toBe(false);
        expect(ti.test(S, A).intersect).toBe(false);
        // Extending the segment past (1,0) produces one crossing.
        expect(fi.find(seg([2, 0], [0.5, 0]), A).numIntersections).toBe(1);
    });

    it('discards crossings that are off the arc', () => {
        // The segment from (-2,0) to (2,0) crosses the circle at (-1,0) and
        // (1,0); only (1,0) is on the first-quadrant arc.
        const result = fi.find(seg([-2, 0], [2, 0]), firstQuadrantArc());
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('finds two crossings when both are on the arc and the segment', () => {
        const A = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(0, -1),
            vec(0, 1));
        const result = fi.find(seg([0.5, -2], [0.5, 2]), A);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        expect(result.point[0].values[1]).toBeCloseTo(-y, 12);
        expect(result.point[1].values[1]).toBeCloseTo(y, 12);
    });

    it('agrees with the line-circle curve query on random inputs', () => {
        let state = 5150;
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
            const S = seg([rand() * 4, rand() * 4], [rand() * 4, rand() * 4]);

            const result = fi.find(S, A);
            expect(ti.test(S, A).intersect).toBe(result.intersect);

            const { center: segOrigin, direction: segDirection,
                extent: segExtent } = S.getCenteredForm();
            const circleResult = lcQuery.find(
                Line.fromOriginDirection(segOrigin, segDirection),
                Hypersphere.fromCenterRadius(center, radius));
            let expected = 0;
            for (let i = 0; i < circleResult.numIntersections; ++i) {
                if (Math.abs(circleResult.parameter[i]) <= segExtent
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
                // The point is on the circle, on the segment and on the arc.
                expect(length(sub(result.point[i], center)))
                    .toBeCloseTo(radius, 8);
                expect(Math.abs(result.parameter[i]))
                    .toBeLessThanOrEqual(segExtent + 1e-12);
                const onSegment = add(segOrigin,
                    mul(result.parameter[i], segDirection));
                expect(length(sub(result.point[i], onSegment)))
                    .toBeCloseTo(0, 10);
                expect(A.containsOnCircle(result.point[i])).toBe(true);
            }
        }
        expect(numOne).toBeGreaterThan(20);
        expect(numTwo).toBeGreaterThanOrEqual(1);
    });
});
