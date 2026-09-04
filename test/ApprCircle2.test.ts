import { describe, it, expect } from 'vitest';
import { ApprCircle2 } from '../src/ApprCircle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, fc, finite, positive, vector } from './helpers/arbitraries.js';

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

describe('ApprCircle2 verification', () => {
    // k points equally spaced on the circle with the given center/radius,
    // starting at the given angular offset.
    const circlePoints = (cx: number, cy: number, r: number, k: number,
        offset: number): Vector[] => {
        const points: Vector[] = [];
        for (let i = 0; i < k; ++i) {
            const t = offset + (2 * Math.PI * i) / k;
            points.push(Vector.fromArray(
                [cx + r * Math.cos(t), cy + r * Math.sin(t)]));
        }
        return points;
    };

    const circleArb = fc.tuple(finite(-5, 5), finite(-5, 5), positive(5, 0.5),
        fc.integer({ min: 3, max: 12 }), finite(0, 2 * Math.PI));

    it('fitUsingSquaredLengths recovers a circle its samples lie on', () => {
        // The algebraic error sum_i (|X_i-C|^2 - r^2)^2 is zero at the true
        // circle, which is therefore the global minimizer the normal
        // equations solve for. Equally spaced samples keep the 2x2 system
        // well conditioned, so the tolerance is a few ulps of the data.
        check(circleArb, ([cx, cy, r, k, offset]) => {
            const points = circlePoints(cx, cy, r, k, offset);
            const circle = new Hypersphere(2);
            expect(new ApprCircle2().fitUsingSquaredLengths(points, circle))
                .toBe(true);
            expectClose(circle.center.get(0), cx, 1e-9, 1e-9);
            expectClose(circle.center.get(1), cy, 1e-9, 1e-9);
            expectClose(circle.radius, r, 1e-9, 1e-9);
        });
    });

    it('fitUsingSquaredLengths is equivariant under translation', () => {
        check(fc.tuple(circleArb, vector(2, -20, 20)),
            ([[cx, cy, r, k, offset], t]) => {
                const points = circlePoints(cx, cy, r, k, offset);
                const shifted = points.map(p => Vector.fromArray(
                    [p.get(0) + t.get(0), p.get(1) + t.get(1)]));

                const a = new Hypersphere(2);
                const b = new Hypersphere(2);
                const fitter = new ApprCircle2();
                expect(fitter.fitUsingSquaredLengths(points, a)).toBe(true);
                expect(fitter.fitUsingSquaredLengths(shifted, b)).toBe(true);

                // Translation moves the fitted center and leaves the radius.
                expectClose(b.center.get(0), a.center.get(0) + t.get(0),
                    1e-8, 1e-8);
                expectClose(b.center.get(1), a.center.get(1) + t.get(1),
                    1e-8, 1e-8);
                expectClose(b.radius, a.radius, 1e-8, 1e-8);
            });
    });

    it('fitUsingSquaredLengths fails and zeroes the circle for degenerate '
        + 'samples', () => {
            // Samples on the x-axis give M01 = M11 = 0 exactly, so the 2x2
            // determinant is exactly zero and upstream takes the failure
            // branch. Coincident samples zero the whole covariance matrix.
            check(fc.tuple(fc.array(finite(-10, 10),
                { minLength: 1, maxLength: 8 }), vector(2, -10, 10)),
                ([xs, p]) => {
                    const onAxis = xs.map(x => Vector.fromArray([x, 0]));
                    const coincident = xs.map(() => p.clone());
                    for (const points of [onAxis, coincident]) {
                        const circle = Hypersphere.fromCenterRadius(
                            Vector.fromArray([7, 8]), 9);
                        expect(new ApprCircle2()
                            .fitUsingSquaredLengths(points, circle))
                            .toBe(false);
                        expect(circle.center.values).toEqual([0, 0]);
                        expect(circle.radius).toBe(0);
                    }
                });
        });

    it('fitUsingLengths keeps a symmetric exact fit at its fixed point',
        () => {
            // For equally spaced samples the average of the unit vectors
            // from the true center vanishes, so C = average + L*dL is the
            // true center and the radius is the true radius.
            check(circleArb, ([cx, cy, r, k, offset]) => {
                const points = circlePoints(cx, cy, r, k, offset);
                const circle = new Hypersphere(2);
                const iterations = new ApprCircle2().fitUsingLengths(
                    points, 8, true, circle);
                expect(iterations).toBeLessThanOrEqual(8);
                expectClose(circle.center.get(0), cx, 1e-8, 1e-8);
                expectClose(circle.center.get(1), cy, 1e-8, 1e-8);
                expectClose(circle.radius, r, 1e-8, 1e-8);
            });
        });

    it('fitUsingLengths honors maxIterations and epsilon', () => {
        check(circleArb, ([cx, cy, r, k, offset]) => {
            const points = circlePoints(cx, cy, r, k, offset);

            // Zero iterations leaves the incoming circle untouched.
            const untouched = Hypersphere.fromCenterRadius(
                Vector.fromArray([1, 2]), 3);
            expect(new ApprCircle2().fitUsingLengths(points, 0, false,
                untouched)).toBe(0);
            expect(untouched.center.values).toEqual([1, 2]);
            expect(untouched.radius).toBe(3);

            // A huge epsilon stops after the very first update, which
            // upstream reports as iteration 0.
            const early = new Hypersphere(2);
            expect(new ApprCircle2().fitUsingLengths(points, 25, true, early,
                1e6)).toBe(0);
        });
    });

    it('neither fit mutates its input samples', () => {
        check(circleArb, ([cx, cy, r, k, offset]) => {
            const points = circlePoints(cx, cy, r, k, offset);
            const before = points.map(p => [...p.values]);
            const circle = new Hypersphere(2);
            const fitter = new ApprCircle2();
            fitter.fitUsingSquaredLengths(points, circle);
            fitter.fitUsingLengths(points, 4, true, circle);
            expect(points.map(p => [...p.values])).toEqual(before);
        });
    });
});
