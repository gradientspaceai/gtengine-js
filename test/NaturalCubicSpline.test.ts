import { describe, it, expect } from 'vitest';
import { NaturalCubicSpline } from '../src/NaturalCubicSpline';
import { Vector, length, sub } from '../src/Vector';

function vec(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// The points of a smooth 3D test curve sampled at the given times.
function samplePoints(times: readonly number[],
    f: (t: number) => number[]): Vector[] {
    return times.map(t => vec(...f(t)));
}

const TIMES = [0, 1, 2.5, 4, 5];

describe('NaturalCubicSpline', () => {
    it('rejects invalid input', () => {
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1)], [0, 1])).toThrow();
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1), vec(2, 0)], [0, 1])).toThrow();
        expect(() => NaturalCubicSpline.createFree(
            [vec(0, 0), vec(1, 1, 1), vec(2, 0)], [0, 1, 2])).toThrow();
    });

    it('reports its construction state, domain and segment count', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        expect(spline.isConstructed()).toBe(true);
        expect(spline.getDimension()).toBe(2);
        expect(spline.getTMin()).toBe(0);
        expect(spline.getTMax()).toBe(5);
        expect(spline.getNumSegments()).toBe(TIMES.length - 1);
        expect(spline.getPolynomials().length).toBe(TIMES.length - 1);
        for (const poly of spline.getPolynomials()) {
            expect(poly.length).toBe(4);
        }
    });

    it('interpolates the knots for all three boundary conditions', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(t), Math.cos(t)]);
        const d0 = vec(1, 1, 0);
        const d1 = vec(1, Math.cos(5), -Math.sin(5));
        const splines = [
            NaturalCubicSpline.createFree(points, TIMES),
            NaturalCubicSpline.createClosed(points, TIMES),
            NaturalCubicSpline.createClamped(points, TIMES, d0, d1)
        ];
        for (const spline of splines) {
            const jet = spline.createJet();
            for (let i = 0; i < TIMES.length; ++i) {
                spline.evaluate(TIMES[i], 0, jet);
                expect(length(sub(jet[0], points[i]))).toBeLessThan(1e-10);
            }
        }
    });

    it('is C2 continuous at the interior knots', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jetL = spline.createJet();
        const jetR = spline.createJet();
        const h = 1e-7;
        for (let i = 1; i < TIMES.length - 1; ++i) {
            const t = TIMES[i];
            spline.evaluate(t - h, 2, jetL);
            spline.evaluate(t + h, 2, jetR);
            for (let order = 0; order <= 2; ++order) {
                expect(length(sub(jetL[order], jetR[order]))).toBeLessThan(1e-5);
            }
        }
    });

    it('has zero second derivative at the ends for a free spline', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t), t * t]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 2, jet);
        expect(length(jet[2])).toBeLessThan(1e-10);
        spline.evaluate(TIMES[TIMES.length - 1], 2, jet);
        expect(length(jet[2])).toBeLessThan(1e-10);
    });

    it('matches the requested end derivatives for a clamped spline', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t)]);
        const d0 = vec(1, 2);
        const d1 = vec(1, 2 * Math.cos(10));
        const spline = NaturalCubicSpline.createClamped(points, TIMES, d0, d1);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 1, jet);
        expect(length(sub(jet[1], d0))).toBeLessThan(1e-10);
        spline.evaluate(TIMES[TIMES.length - 1], 1, jet);
        expect(length(sub(jet[1], d1))).toBeLessThan(1e-10);
    });

    it('matches the first and second derivatives at the ends when closed', () => {
        const points = samplePoints(TIMES, t => [Math.cos(t), Math.sin(t)]);
        const spline = NaturalCubicSpline.createClosed(points, TIMES);
        const jet0 = spline.createJet();
        const jet1 = spline.createJet();
        spline.evaluate(TIMES[0], 2, jet0);
        spline.evaluate(TIMES[TIMES.length - 1], 2, jet1);
        expect(length(sub(jet0[1], jet1[1]))).toBeLessThan(1e-10);
        expect(length(sub(jet0[2], jet1[2]))).toBeLessThan(1e-10);
    });

    it('reproduces a cubic polynomial exactly when clamped', () => {
        const P = (t: number) => [
            1 + 2 * t - 0.5 * t * t + 0.25 * t * t * t,
            -3 + t * t * t,
            7 - t
        ];
        const D = (t: number) => [2 - t + 0.75 * t * t, 3 * t * t, -1];
        const points = samplePoints(TIMES, P);
        const spline = NaturalCubicSpline.createClamped(points, TIMES,
            vec(...D(0)), vec(...D(5)));
        const jet = spline.createJet();
        for (const t of [0, 0.3, 1, 1.7, 2.5, 3.3, 4, 4.9, 5]) {
            spline.evaluate(t, 3, jet);
            expect(length(sub(jet[0], vec(...P(t))))).toBeLessThan(1e-10);
            expect(length(sub(jet[1], vec(...D(t))))).toBeLessThan(1e-10);
            // The third derivative of the reproduced cubic is constant.
            expect(jet[3].get(0)).toBeCloseTo(1.5, 8);
            expect(jet[3].get(1)).toBeCloseTo(6, 8);
            expect(jet[3].get(2)).toBeCloseTo(0, 8);
        }
    });

    it('derivatives agree with finite differences of the position', () => {
        const points = samplePoints(TIMES, t => [t, Math.sin(2 * t), Math.exp(-t)]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        const jetP = spline.createJet();
        const jetM = spline.createJet();
        const h = 1e-6;
        for (const t of [0.4, 1.6, 2.9, 3.5, 4.6]) {
            spline.evaluate(t, 2, jet);
            spline.evaluate(t + h, 1, jetP);
            spline.evaluate(t - h, 1, jetM);
            for (let i = 0; i < 3; ++i) {
                const d1 = (jetP[0].get(i) - jetM[0].get(i)) / (2 * h);
                const d2 = (jetP[1].get(i) - jetM[1].get(i)) / (2 * h);
                expect(jet[1].get(i)).toBeCloseTo(d1, 6);
                expect(jet[2].get(i)).toBeCloseTo(d2, 5);
            }
        }
    });

    it('clamps evaluation outside the time interval to the end segments', () => {
        const points = samplePoints(TIMES, t => [t, t * t]);
        const spline = NaturalCubicSpline.createFree(points, TIMES);
        const jet = spline.createJet();
        spline.evaluate(-10, 0, jet);
        const before = jet[0].clone();
        spline.evaluate(TIMES[0], 0, jet);
        expect(length(sub(before, jet[0]))).toBeLessThan(1e-12);
        spline.evaluate(100, 0, jet);
        const after = jet[0].clone();
        spline.evaluate(TIMES[TIMES.length - 1], 0, jet);
        expect(length(sub(after, jet[0]))).toBeLessThan(1e-12);
    });

    it('interpolates random data with C2 continuity (randomized)', () => {
        let seed = 24680;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const n = 3 + Math.floor(6 * rand());
            const times: number[] = [0];
            for (let i = 1; i < n; ++i) {
                times.push(times[i - 1] + 0.25 + rand());
            }
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                points.push(vec(4 * rand() - 2, 4 * rand() - 2));
            }
            const spline = NaturalCubicSpline.createFree(points, times);
            const jet = spline.createJet();
            for (let i = 0; i < n; ++i) {
                spline.evaluate(times[i], 0, jet);
                expect(length(sub(jet[0], points[i]))).toBeLessThan(1e-8);
            }
            // Natural boundary conditions.
            spline.evaluate(times[0], 2, jet);
            expect(length(jet[2])).toBeLessThan(1e-8);
            spline.evaluate(times[n - 1], 2, jet);
            expect(length(jet[2])).toBeLessThan(1e-8);
            // C2 across interior knots.
            const h = 1e-7;
            const jetL = spline.createJet();
            const jetR = spline.createJet();
            for (let i = 1; i < n - 1; ++i) {
                spline.evaluate(times[i] - h, 2, jetL);
                spline.evaluate(times[i] + h, 2, jetR);
                expect(length(sub(jetL[2], jetR[2]))).toBeLessThan(1e-4);
            }
        }
    });
});
