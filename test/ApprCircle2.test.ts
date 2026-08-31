import { describe, it, expect } from 'vitest';
import { ApprCircle2 } from '../src/ApprCircle2';
import { Hypersphere } from '../src/Hypersphere';
import { Vector } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Points exactly on the circle with the given center and radius.
function circlePoints(cx: number, cy: number, r: number, count: number,
    startAngle: number = 0): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        const angle = startAngle + 2 * Math.PI * i / count;
        points.push(v2(cx + r * Math.cos(angle), cy + r * Math.sin(angle)));
    }
    return points;
}

// A deterministic pseudo-random generator for the noisy-data tests.
function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprCircle2.fitUsingSquaredLengths', () => {
    it('recovers a circle from points that lie exactly on it', () => {
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        const points = circlePoints(3, -2, 5, 16);
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(true);
        expect(circle.center.values[0]).toBeCloseTo(3, 10);
        expect(circle.center.values[1]).toBeCloseTo(-2, 10);
        expect(circle.radius).toBeCloseTo(5, 10);
    });

    it('recovers a unit circle centered at the origin', () => {
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        const points = circlePoints(0, 0, 1, 7, 0.3);
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(true);
        expect(circle.center.values[0]).toBeCloseTo(0, 12);
        expect(circle.center.values[1]).toBeCloseTo(0, 12);
        expect(circle.radius).toBeCloseTo(1, 12);
    });

    it('fits three noncollinear points to their circumscribed circle', () => {
        // The unique circle through the three points is the least-squares
        // solution, so the fitted radius is the distance to each point.
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        const points = [v2(1, 0), v2(0, 1), v2(-1, 0)];
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(true);
        expect(circle.center.values[0]).toBeCloseTo(0, 12);
        expect(circle.center.values[1]).toBeCloseTo(0, 12);
        expect(circle.radius).toBeCloseTo(1, 12);
    });

    it('fails for collinear points and zeroes the circle', () => {
        const fitter = new ApprCircle2();
        const circle = Hypersphere.fromCenterRadius(v2(9, 9), 9);
        const points = [v2(-2, -4), v2(0, 0), v2(1, 2), v2(3, 6)];
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(false);
        expect(circle.center.values).toEqual([0, 0]);
        expect(circle.radius).toBe(0);
    });

    it('fails for coincident points', () => {
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        const points = [v2(2, 3), v2(2, 3), v2(2, 3)];
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(false);
        expect(circle.center.values).toEqual([0, 0]);
        expect(circle.radius).toBe(0);
    });

    it('is close to the true circle for noisy samples', () => {
        const random = makeRandom(12345);
        const points: Vector[] = [];
        for (let i = 0; i < 200; ++i) {
            const angle = 2 * Math.PI * i / 200;
            const r = 4 + 0.01 * (2 * random() - 1);
            points.push(v2(-1 + r * Math.cos(angle), 6 + r * Math.sin(angle)));
        }
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        expect(fitter.fitUsingSquaredLengths(points, circle)).toBe(true);
        expect(circle.center.values[0]).toBeCloseTo(-1, 2);
        expect(circle.center.values[1]).toBeCloseTo(6, 2);
        expect(circle.radius).toBeCloseTo(4, 2);
    });
});

describe('ApprCircle2.fitUsingLengths', () => {
    it('converges to the exact circle when starting from the average', () => {
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        const points = circlePoints(3, -2, 5, 32);
        const iterations = fitter.fitUsingLengths(points, 1024, true, circle, 1e-14);
        expect(iterations).toBeLessThan(1024);
        expect(circle.center.values[0]).toBeCloseTo(3, 8);
        expect(circle.center.values[1]).toBeCloseTo(-2, 8);
        expect(circle.radius).toBeCloseTo(5, 8);
    });

    it('agrees with the squared-length fit for exact data', () => {
        const points = circlePoints(-4, 7, 2.5, 24, 0.17);
        const fitter = new ApprCircle2();
        const circleA = new Hypersphere(2);
        fitter.fitUsingSquaredLengths(points, circleA);
        const circleB = new Hypersphere(2);
        fitter.fitUsingLengths(points, 4096, true, circleB, 0);
        expect(circleB.center.values[0]).toBeCloseTo(circleA.center.values[0], 6);
        expect(circleB.center.values[1]).toBeCloseTo(circleA.center.values[1], 6);
        expect(circleB.radius).toBeCloseTo(circleA.radius, 6);
    });

    it('uses the incoming center as the initial guess when the flag is false', () => {
        // With maxIterations = 0 the loop body never executes, so the
        // incoming circle is left untouched and zero iterations are used.
        const fitter = new ApprCircle2();
        const circle = Hypersphere.fromCenterRadius(v2(1, 1), 3);
        const points = circlePoints(3, -2, 5, 8);
        expect(fitter.fitUsingLengths(points, 0, false, circle)).toBe(0);
        expect(circle.center.values).toEqual([1, 1]);
        expect(circle.radius).toBe(3);
    });

    it('overwrites the incoming center with the average when the flag is true', () => {
        // One iteration from the average of the samples: the average of a
        // symmetric point set is the true center, which is a fixed point of
        // the iteration, so the result is the exact circle after one step.
        const fitter = new ApprCircle2();
        const circle = Hypersphere.fromCenterRadius(v2(100, 100), 1000);
        const points = circlePoints(3, -2, 5, 8);
        expect(fitter.fitUsingLengths(points, 1, true, circle)).toBe(1);
        expect(circle.center.values[0]).toBeCloseTo(3, 12);
        expect(circle.center.values[1]).toBeCloseTo(-2, 12);
        expect(circle.radius).toBeCloseTo(5, 12);
    });

    it('polishes a fit when restarted from the previous center', () => {
        const points = circlePoints(0.5, 0.25, 1.75, 40, 0.05);
        const fitter = new ApprCircle2();
        const circle = new Hypersphere(2);
        fitter.fitUsingLengths(points, 4, true, circle, 0);
        const iterations = fitter.fitUsingLengths(points, 1024, false, circle, 1e-15);
        expect(iterations).toBeLessThan(1024);
        expect(circle.center.values[0]).toBeCloseTo(0.5, 8);
        expect(circle.center.values[1]).toBeCloseTo(0.25, 8);
        expect(circle.radius).toBeCloseTo(1.75, 8);
    });

    it('skips samples that coincide with the current center', () => {
        // The sample at the center contributes no length or derivative term.
        const fitter = new ApprCircle2();
        const circle = Hypersphere.fromCenterRadius(v2(0, 0), 1);
        const points = [v2(0, 0), v2(1, 0), v2(-1, 0), v2(0, 1), v2(0, -1)];
        fitter.fitUsingLengths(points, 1, false, circle);
        // Four unit-length samples and one zero-length sample.
        expect(circle.radius).toBeCloseTo(4 / 5, 12);
        expect(circle.center.values[0]).toBeCloseTo(0, 12);
        expect(circle.center.values[1]).toBeCloseTo(0, 12);
    });
});
