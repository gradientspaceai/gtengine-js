import { describe, it, expect } from 'vitest';
import { NaturalQuinticSpline } from '../src/NaturalQuinticSpline';
import { Vector, length, sub } from '../src/Vector';

function vec(...values: number[]): Vector {
    return Vector.fromArray(values);
}

const TIMES = [0, 1, 2.5, 4, 5];

// A quintic polynomial and its first three derivatives, per component.
const Q = (t: number) => [
    1 + t - 0.3 * t * t + 0.1 * t ** 3 + 0.02 * t ** 4 - 0.003 * t ** 5,
    -2 + 0.5 * t ** 5,
    3 - t + 2 * t * t
];
const Q1 = (t: number) => [
    1 - 0.6 * t + 0.3 * t * t + 0.08 * t ** 3 - 0.015 * t ** 4,
    2.5 * t ** 4,
    -1 + 4 * t
];
const Q2 = (t: number) => [
    -0.6 + 0.6 * t + 0.24 * t * t - 0.06 * t ** 3,
    10 * t ** 3,
    4
];

function sample(times: readonly number[],
    f: (t: number) => number[]): Vector[] {
    return times.map(t => vec(...f(t)));
}

describe('NaturalQuinticSpline', () => {
    it('rejects invalid input', () => {
        expect(() => NaturalQuinticSpline.createFree(
            [vec(0), vec(1)], [vec(0), vec(0)], [0, 1])).toThrow();
        expect(() => NaturalQuinticSpline.createFree(
            [vec(0), vec(1), vec(2)], [vec(0), vec(0)], [0, 1, 2])).toThrow();
        expect(() => NaturalQuinticSpline.createFree(
            [vec(0), vec(1), vec(2)], [vec(0), vec(0), vec(0)],
            [0, 1])).toThrow();
    });

    it('reports its construction state, domain and segments', () => {
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        expect(spline.isConstructed()).toBe(true);
        expect(spline.getDimension()).toBe(3);
        expect(spline.getTMin()).toBe(0);
        expect(spline.getTMax()).toBe(5);
        expect(spline.getNumSegments()).toBe(TIMES.length - 1);
        expect(spline.getPolynomials().length).toBe(TIMES.length - 1);
        for (const poly of spline.getPolynomials()) {
            expect(poly.length).toBe(6);
        }
        // createJet allocates enough entries for the fifth derivative.
        expect(spline.createJet().length)
            .toBe(NaturalQuinticSpline.SUP_ORDER_QUINTIC);
    });

    it('interpolates the values and first derivatives at the knots', () => {
        const f0 = sample(TIMES, Q);
        const f1 = sample(TIMES, Q1);
        const splines = [
            NaturalQuinticSpline.createFree(f0, f1, TIMES),
            NaturalQuinticSpline.createClosed(f0, f1, TIMES),
            NaturalQuinticSpline.createClamped(f0, f1, TIMES,
                vec(...Q2(0)), vec(...Q2(5)))
        ];
        for (const spline of splines) {
            const jet = spline.createJet();
            for (let i = 0; i < TIMES.length; ++i) {
                spline.evaluate(TIMES[i], 1, jet);
                expect(length(sub(jet[0], f0[i]))).toBeLessThan(1e-9);
                expect(length(sub(jet[1], f1[i]))).toBeLessThan(1e-9);
            }
        }
    });

    it('is C3 continuous at the interior knots', () => {
        // The quintic spline has 6 coefficients per segment. Interpolating
        // the value and the first derivative at both ends of each segment
        // uses 4 of them, leaving 2(n-1) degrees of freedom against
        // 2(n-2) interior continuity conditions plus 2 boundary conditions.
        // Continuity therefore holds through order 3, and the fourth
        // derivative jumps in general.
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jetL = spline.createJet();
        const jetR = spline.createJet();
        const h = 1e-8;
        for (let i = 1; i < TIMES.length - 1; ++i) {
            spline.evaluate(TIMES[i] - h, 3, jetL);
            spline.evaluate(TIMES[i] + h, 3, jetR);
            for (let order = 0; order <= 3; ++order) {
                // The tolerance accommodates the round-off of dividing by
                // delta^order in the jet evaluation.
                expect(length(sub(jetL[order], jetR[order]))).toBeLessThan(1e-4);
            }
        }
    });

    it('has zero third derivative at the ends for a free spline', () => {
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 3, jet);
        expect(length(jet[3])).toBeLessThan(1e-9);
        spline.evaluate(TIMES[TIMES.length - 1], 3, jet);
        expect(length(jet[3])).toBeLessThan(1e-9);
    });

    it('matches the requested end second derivatives when clamped', () => {
        const d0 = vec(...Q2(0));
        const d1 = vec(...Q2(5));
        const spline = NaturalQuinticSpline.createClamped(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES, d0, d1);
        const jet = spline.createJet();
        spline.evaluate(TIMES[0], 2, jet);
        expect(length(sub(jet[2], d0))).toBeLessThan(1e-9);
        spline.evaluate(TIMES[TIMES.length - 1], 2, jet);
        expect(length(sub(jet[2], d1))).toBeLessThan(1e-9);
    });

    it('matches the second and third derivatives at the ends when closed', () => {
        const spline = NaturalQuinticSpline.createClosed(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jet0 = spline.createJet();
        const jet1 = spline.createJet();
        spline.evaluate(TIMES[0], 3, jet0);
        spline.evaluate(TIMES[TIMES.length - 1], 3, jet1);
        expect(length(sub(jet0[2], jet1[2]))).toBeLessThan(1e-8);
        expect(length(sub(jet0[3], jet1[3]))).toBeLessThan(1e-8);
    });

    it('reproduces a quintic polynomial exactly when clamped', () => {
        const spline = NaturalQuinticSpline.createClamped(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES, vec(...Q2(0)), vec(...Q2(5)));
        const jet = spline.createJet();
        for (const t of [0, 0.3, 1, 1.7, 2.5, 3.3, 4, 4.9, 5]) {
            spline.evaluate(t, 2, jet);
            expect(length(sub(jet[0], vec(...Q(t))))).toBeLessThan(1e-9);
            expect(length(sub(jet[1], vec(...Q1(t))))).toBeLessThan(1e-9);
            expect(length(sub(jet[2], vec(...Q2(t))))).toBeLessThan(1e-8);
        }
    });

    it('has the constant fifth derivative 120*c5/delta^5 on each segment', () => {
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jet = spline.createJet();
        const polys = spline.getPolynomials();
        for (let i = 0; i < polys.length; ++i) {
            const delta = TIMES[i + 1] - TIMES[i];
            const t = 0.5 * (TIMES[i] + TIMES[i + 1]);
            spline.evaluate(t, 5, jet);
            for (let k = 0; k < 3; ++k) {
                expect(jet[5].get(k)).toBeCloseTo(
                    120 * polys[i][5].get(k) / delta ** 5, 8);
            }
        }
    });

    it('derivatives of orders 0-3 agree with finite differences', () => {
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jet = spline.createJet();
        const jetP = spline.createJet();
        const jetM = spline.createJet();
        const h = 1e-6;
        for (const t of [0.4, 1.6, 2.9, 3.5, 4.6]) {
            spline.evaluate(t, 3, jet);
            spline.evaluate(t + h, 3, jetP);
            spline.evaluate(t - h, 3, jetM);
            for (let i = 0; i < 3; ++i) {
                for (let order = 1; order <= 3; ++order) {
                    const fd = (jetP[order - 1].get(i)
                        - jetM[order - 1].get(i)) / (2 * h);
                    expect(Math.abs(jet[order].get(i) - fd)).toBeLessThan(1e-3);
                }
            }
        }
    });

    it('clamps evaluation outside the time interval', () => {
        const spline = NaturalQuinticSpline.createFree(sample(TIMES, Q),
            sample(TIMES, Q1), TIMES);
        const jet = spline.createJet();
        spline.evaluate(-20, 0, jet);
        const before = jet[0].clone();
        spline.evaluate(TIMES[0], 0, jet);
        expect(length(sub(before, jet[0]))).toBeLessThan(1e-12);
        spline.evaluate(50, 0, jet);
        const after = jet[0].clone();
        spline.evaluate(TIMES[TIMES.length - 1], 0, jet);
        expect(length(sub(after, jet[0]))).toBeLessThan(1e-12);
    });

    it('interpolates random Hermite data with C3 continuity (randomized)', () => {
        let seed = 1357911;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 50; ++trial) {
            const n = 3 + Math.floor(5 * rand());
            const times: number[] = [0];
            for (let i = 1; i < n; ++i) {
                times.push(times[i - 1] + 0.5 + rand());
            }
            const f0: Vector[] = [];
            const f1: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                f0.push(vec(2 * rand() - 1, 2 * rand() - 1));
                f1.push(vec(2 * rand() - 1, 2 * rand() - 1));
            }
            const spline = NaturalQuinticSpline.createFree(f0, f1, times);
            const jet = spline.createJet();
            for (let i = 0; i < n; ++i) {
                spline.evaluate(times[i], 1, jet);
                expect(length(sub(jet[0], f0[i]))).toBeLessThan(1e-7);
                expect(length(sub(jet[1], f1[i]))).toBeLessThan(1e-7);
            }
            // Natural boundary conditions: third derivative vanishes.
            spline.evaluate(times[0], 3, jet);
            expect(length(jet[3])).toBeLessThan(1e-6);
            spline.evaluate(times[n - 1], 3, jet);
            expect(length(jet[3])).toBeLessThan(1e-6);
            // C3 across the interior knots.
            const h = 1e-8;
            const jetL = spline.createJet();
            const jetR = spline.createJet();
            for (let i = 1; i < n - 1; ++i) {
                spline.evaluate(times[i] - h, 3, jetL);
                spline.evaluate(times[i] + h, 3, jetR);
                for (let order = 0; order <= 3; ++order) {
                    expect(length(sub(jetL[order], jetR[order])))
                        .toBeLessThan(1e-4);
                }
            }
        }
    });
});
