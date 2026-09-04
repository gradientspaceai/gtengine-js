import { describe, it, expect } from 'vitest';
import { NaturalQuinticSpline } from '../src/NaturalQuinticSpline.js';
import { Vector, length, sub } from '../src/Vector.js';

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


// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). NaturalQuinticSpline.h was
// read line by line against src/NaturalQuinticSpline.ts. As for the cubic
// spline, the block solver is checked against the conditions that define the
// spline: on segment i the curve is the quintic
// p_i(u) = c0 + u*c1 + ... + u^5*c5 with u = (t - times[i])/delta[i], so
// interpolation of both f0 and f1, the C2 and C3 joins and each factory's
// boundary condition are exact linear identities on the coefficients returned
// by GetPolynomials.
import {
    check, fc, expectClose, expectVectorClose, wellScaledVector
} from './helpers/arbitraries.js';

interface QuinticInput {
    f0: Vector[];
    f1: Vector[];
    times: number[];
    deltas: number[];
}

const quinticInput = (dim: number, minPoints = 3, maxPoints = 6):
    fc.Arbitrary<QuinticInput> =>
    fc.integer({ min: minPoints, max: maxPoints }).chain(n => fc.tuple(
        fc.array(wellScaledVector(dim, -4, 4), { minLength: n, maxLength: n }),
        fc.array(wellScaledVector(dim, -3, 3), { minLength: n, maxLength: n }),
        fc.array(fc.double({ min: 0.4, max: 2, noNaN: true,
            noDefaultInfinity: true }),
        { minLength: n - 1, maxLength: n - 1 })))
        .map(([f0, f1, gaps]) => {
            const times = [0];
            for (const g of gaps) { times.push(times[times.length - 1] + g); }
            return { f0, f1, times, deltas: gaps };
        });

// The value and the first three derivatives with respect to t of segment i at
// its local parameter u.
function quinticJet(poly: readonly Vector[], delta: number, u: number):
    Vector[] {
    const n = poly[0].size;
    const out = [new Vector(n), new Vector(n), new Vector(n), new Vector(n)];
    for (let k = 0; k < n; ++k) {
        const c = [0, 1, 2, 3, 4, 5].map(j => poly[j].values[k]);
        out[0].values[k] = c[0] + u * (c[1] + u * (c[2] + u * (c[3] +
            u * (c[4] + u * c[5]))));
        out[1].values[k] = (c[1] + u * (2 * c[2] + u * (3 * c[3] +
            u * (4 * c[4] + u * 5 * c[5])))) / delta;
        out[2].values[k] = (2 * c[2] + u * (6 * c[3] + u * (12 * c[4] +
            u * 20 * c[5]))) / (delta * delta);
        out[3].values[k] = (6 * c[3] + u * (24 * c[4] + u * 60 * c[5])) /
            (delta * delta * delta);
    }
    return out;
}

function expectQuinticIdentities(curve: NaturalQuinticSpline,
    input: QuinticInput): void {
    const polys = curve.getPolynomials();
    const { f0, f1, deltas } = input;
    expect(polys.length).toBe(f0.length - 1);
    for (let i = 0; i < polys.length; ++i) {
        const begin = quinticJet(polys[i], deltas[i], 0);
        expectVectorClose(begin[0], f0[i], 1e-9, 1e-9);
        expectVectorClose(begin[1], f1[i], 1e-9, 1e-9);
    }
    const last = polys.length - 1;
    const end = quinticJet(polys[last], deltas[last], 1);
    expectVectorClose(end[0], f0[f0.length - 1], 1e-7, 1e-7);
    expectVectorClose(end[1], f1[f1.length - 1], 1e-7, 1e-7);

    // C0 through C3 joins at the interior knots. The spline is only C3: the
    // fourth derivative generally jumps, which is what distinguishes this
    // construction from a C4 quintic spline.
    for (let i = 0; i + 1 < polys.length; ++i) {
        const left = quinticJet(polys[i], deltas[i], 1);
        const right = quinticJet(polys[i + 1], deltas[i + 1], 0);
        for (let order = 0; order <= 3; ++order) {
            expectVectorClose(left[order], right[order], 1e-7, 1e-7);
        }
    }
}

describe('NaturalQuinticSpline verification', () => {
    it('free splines interpolate f0 and f1, are C3 and have zero end jerk',
        () => {
            check(quinticInput(3), input => {
                const curve = NaturalQuinticSpline.createFree(input.f0,
                    input.f1, input.times);
                expectQuinticIdentities(curve, input);
                const polys = curve.getPolynomials();
                const last = polys.length - 1;
                const zero = new Vector(3);
                expectVectorClose(
                    quinticJet(polys[0], input.deltas[0], 0)[3], zero,
                    1e-8, 1e-8);
                expectVectorClose(
                    quinticJet(polys[last], input.deltas[last], 1)[3], zero,
                    1e-7, 1e-7);
            });
        });

    it('clamped splines match the requested end second derivatives', () => {
        check(fc.tuple(quinticInput(2), wellScaledVector(2, -3, 3),
            wellScaledVector(2, -3, 3)), ([input, d0, d1]) => {
            const curve = NaturalQuinticSpline.createClamped(input.f0,
                input.f1, input.times, d0, d1);
            expectQuinticIdentities(curve, input);
            const polys = curve.getPolynomials();
            const last = polys.length - 1;
            expectVectorClose(quinticJet(polys[0], input.deltas[0], 0)[2], d0,
                1e-8, 1e-8);
            expectVectorClose(
                quinticJet(polys[last], input.deltas[last], 1)[2], d1,
                1e-7, 1e-7);
        });
    });

    it('closed splines match the end second and third derivatives', () => {
        check(quinticInput(2), input => {
            const curve = NaturalQuinticSpline.createClosed(input.f0,
                input.f1, input.times);
            expectQuinticIdentities(curve, input);
            const polys = curve.getPolynomials();
            const last = polys.length - 1;
            const begin = quinticJet(polys[0], input.deltas[0], 0);
            const end = quinticJet(polys[last], input.deltas[last], 1);
            expectVectorClose(begin[2], end[2], 1e-7, 1e-7);
            expectVectorClose(begin[3], end[3], 1e-7, 1e-7);
        });
    });

    it('evaluate agrees with the stored polynomial coefficients', () => {
        check(fc.tuple(quinticInput(3), fc.double({ min: 0, max: 0.999,
            noNaN: true, noDefaultInfinity: true })), ([input, u]) => {
            const curve = NaturalQuinticSpline.createFree(input.f0, input.f1,
                input.times);
            const polys = curve.getPolynomials();
            for (let i = 0; i < polys.length; ++i) {
                const t = input.times[i] + u * input.deltas[i];
                const expected = quinticJet(polys[i], input.deltas[i], u);
                const jet = curve.createJet();
                curve.evaluate(t, 3, jet);
                for (let order = 0; order <= 3; ++order) {
                    expectVectorClose(jet[order], expected[order],
                        1e-8, 1e-8);
                }
            }
        });
    });

    it('clamped interpolation reproduces quintic polynomials exactly', () => {
        // A quintic sampled with its own values and first derivatives, plus
        // its own end second derivatives, satisfies every condition of the
        // clamped construction, so the spline must reproduce it.
        check(fc.tuple(quinticInput(2, 3, 5),
            fc.array(wellScaledVector(2, -1.5, 1.5),
                { minLength: 6, maxLength: 6 }),
            fc.double({ min: 0, max: 1, noNaN: true,
                noDefaultInfinity: true })), ([input, c, s]) => {
            const poly = (t: number, r: number): Vector =>
                Vector.fromArray([0, 1].map(k => {
                    let sum = 0;
                    for (let j = r; j < 6; ++j) {
                        let factor = 1;
                        for (let q = 0; q < r; ++q) { factor *= j - q; }
                        sum += factor * c[j].values[k] * Math.pow(t, j - r);
                    }
                    return sum;
                }));
            const times = input.times;
            const tmax = times[times.length - 1];
            const curve = NaturalQuinticSpline.createClamped(
                times.map(t => poly(t, 0)), times.map(t => poly(t, 1)), times,
                poly(times[0], 2), poly(tmax, 2));
            const t = s * tmax;
            const jet = curve.createJet();
            curve.evaluate(t, 2, jet);
            expectVectorClose(jet[0], poly(t, 0), 1e-7, 1e-7);
            expectVectorClose(jet[1], poly(t, 1), 1e-7, 1e-7);
            expectVectorClose(jet[2], poly(t, 2), 1e-6, 1e-6);
        });
    });

    it('has a piecewise constant fifth derivative and clamps the domain', () => {
        check(fc.tuple(quinticInput(2), fc.double({ min: 0.001, max: 50,
            noNaN: true, noDefaultInfinity: true })), ([input, d]) => {
            const curve = NaturalQuinticSpline.createFree(input.f0, input.f1,
                input.times);
            const polys = curve.getPolynomials();
            const jet = curve.createJet();
            expect(jet.length).toBeGreaterThanOrEqual(6);
            for (let i = 0; i < polys.length; ++i) {
                const delta = input.deltas[i];
                const expected = new Vector(2);
                for (let k = 0; k < 2; ++k) {
                    expected.values[k] = 120 * polys[i][5].values[k] /
                        Math.pow(delta, 5);
                }
                for (const u of [0.1, 0.5, 0.9]) {
                    curve.evaluate(input.times[i] + u * delta, 5, jet);
                    expectVectorClose(jet[5], expected, 1e-8, 1e-8);
                }
            }
            curve.evaluate(input.times[0] - d, 0, jet);
            expectVectorClose(jet[0], input.f0[0], 1e-9, 1e-9);
            curve.evaluate(input.times[input.times.length - 1] + d, 0, jet);
            expectVectorClose(jet[0], input.f0[input.f0.length - 1],
                1e-7, 1e-7);
        }, 50);
    });
});
