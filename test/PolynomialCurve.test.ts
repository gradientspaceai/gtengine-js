import { describe, it, expect } from 'vitest';
import { PolynomialCurve } from '../src/PolynomialCurve.js';
import { Polynomial1 } from '../src/Polynomial1.js';
import { Vector, length as vectorLength, sub, dot } from '../src/Vector.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function poly(...coefficients: number[]): Polynomial1 {
    return Polynomial1.fromCoefficients(coefficients);
}

function jetOf(curve: PolynomialCurve, t: number, order: number): Vector[] {
    const jet = curve.createJet();
    curve.evaluate(t, order, jet);
    return jet;
}

// Direct evaluation of sum_k c[k] t^k and its derivatives, independent of
// the Polynomial1 derivative machinery.
function evalDerivative(coefficients: number[], t: number, d: number): number {
    let sum = 0;
    for (let k = d; k < coefficients.length; ++k) {
        let factor = 1;
        for (let j = 0; j < d; ++j) {
            factor *= k - j;
        }
        sum += coefficients[k] * factor * Math.pow(t, k - d);
    }
    return sum;
}

describe('PolynomialCurve', () => {
    it('default construction gives the constant zero curve', () => {
        const curve = new PolynomialCurve(3, 0, 1);
        expect(curve.isConstructed()).toBe(true);
        expect(curve.getDimension()).toBe(3);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);

        for (let i = 0; i < 3; ++i) {
            expect(curve.getPolynomial(i).getDegree()).toBe(0);
            expect(curve.getPolynomial(i).get(0)).toBe(0);
            expect(curve.getDer1Polynomial(i).get(0)).toBe(0);
            expect(curve.getDer2Polynomial(i).get(0)).toBe(0);
            expect(curve.getDer3Polynomial(i).get(0)).toBe(0);
        }

        const jet = jetOf(curve, 0.375, 3);
        for (let k = 0; k < 4; ++k) {
            for (let i = 0; i < 3; ++i) {
                expect(jet[k].values[i]).toBe(0);
            }
        }
    });

    it('evaluates a 2D curve against hand-computed values', () => {
        // x(t) = 1 + 2t + 3t^2, y(t) = -1 + t^3.
        const cx = [1, 2, 3];
        const cy = [-1, 0, 0, 1];
        const curve = new PolynomialCurve(2, -1, 2, [poly(...cx), poly(...cy)]);

        const t = 0.5;
        const jet = jetOf(curve, t, 3);

        // x: 1 + 1 + 0.75 = 2.75; x' = 2 + 6t = 5; x'' = 6; x''' = 0.
        expect(jet[0].values[0]).toBeCloseTo(2.75, 14);
        expect(jet[1].values[0]).toBeCloseTo(5, 14);
        expect(jet[2].values[0]).toBeCloseTo(6, 14);
        expect(jet[3].values[0]).toBeCloseTo(0, 14);

        // y = -1 + 0.125 = -0.875; y' = 3t^2 = 0.75; y'' = 6t = 3; y''' = 6.
        expect(jet[0].values[1]).toBeCloseTo(-0.875, 14);
        expect(jet[1].values[1]).toBeCloseTo(0.75, 14);
        expect(jet[2].values[1]).toBeCloseTo(3, 14);
        expect(jet[3].values[1]).toBeCloseTo(6, 14);
    });

    it('exposes the derivative polynomials with the expected coefficients', () => {
        // x(t) = 5 - 4t + 3t^2 - 2t^3 + t^4.
        const curve = new PolynomialCurve(1, 0, 1, [poly(5, -4, 3, -2, 1)]);
        expect(curve.getPolynomial(0).getCoefficients()).toEqual([5, -4, 3, -2, 1]);
        expect(curve.getDer1Polynomial(0).getCoefficients()).toEqual([-4, 6, -6, 4]);
        expect(curve.getDer2Polynomial(0).getCoefficients()).toEqual([6, -12, 12]);
        expect(curve.getDer3Polynomial(0).getCoefficients()).toEqual([-12, 24]);
    });

    it('setPolynomial copies the input and refreshes the derivatives', () => {
        const curve = new PolynomialCurve(2, 0, 1);
        const p = poly(0, 0, 4);  // 4 t^2
        curve.setPolynomial(1, p);
        expect(curve.getDer1Polynomial(1).getCoefficients()).toEqual([0, 8]);

        // Mutating the source polynomial must not affect the curve (C++ copy
        // semantics).
        p.set(2, 100);
        expect(curve.getPolynomial(1).get(2)).toBe(4);
        expect(curve.getDer1Polynomial(1).get(1)).toBe(8);

        const jet = jetOf(curve, 2, 2);
        expect(jet[0].values[1]).toBeCloseTo(16, 12);
        expect(jet[1].values[1]).toBeCloseTo(16, 12);
        expect(jet[2].values[1]).toBeCloseTo(8, 12);
        expect(jet[0].values[0]).toBe(0);
    });

    it('honors the requested evaluation order', () => {
        const curve = new PolynomialCurve(1, 0, 1, [poly(1, 1, 1, 1)]);

        const jet0 = jetOf(curve, 1, 0);
        expect(jet0[0].values[0]).toBeCloseTo(4, 14);
        expect(jet0[1].values[0]).toBe(0);
        expect(jet0[2].values[0]).toBe(0);
        expect(jet0[3].values[0]).toBe(0);

        const jet1 = jetOf(curve, 1, 1);
        expect(jet1[1].values[0]).toBeCloseTo(6, 14);
        expect(jet1[2].values[0]).toBe(0);
        expect(jet1[3].values[0]).toBe(0);

        const jet2 = jetOf(curve, 1, 2);
        expect(jet2[2].values[0]).toBeCloseTo(8, 14);
        expect(jet2[3].values[0]).toBe(0);

        const jet3 = jetOf(curve, 1, 3);
        expect(jet3[3].values[0]).toBeCloseTo(6, 14);
    });

    it('preserves the upstream order == 3 test for the third derivative', () => {
        // Upstream uses 'order == 3' rather than 'order >= 3', so an order
        // larger than the documented maximum leaves jet[3] untouched.
        const curve = new PolynomialCurve(1, 0, 1, [poly(0, 0, 0, 1)]);
        const jet = curve.createJet();
        jet[3].values[0] = 12345;
        curve.evaluate(1, 4, jet);
        expect(jet[0].values[0]).toBeCloseTo(1, 14);
        expect(jet[2].values[0]).toBeCloseTo(6, 14);
        expect(jet[3].values[0]).toBe(12345);
    });

    it('matches finite differences of the position (3D)', () => {
        const cx = [0.5, -1.25, 2, 0.75];
        const cy = [3, 0, -0.5, 0, 0.25];
        const cz = [-2, 1.5];
        const curve = new PolynomialCurve(3, -2, 2,
            [poly(...cx), poly(...cy), poly(...cz)]);
        const all = [cx, cy, cz];

        const h = 1e-4;
        for (const t of [-1.5, -0.5, 0, 0.25, 1.75]) {
            const jet = jetOf(curve, t, 3);
            for (let i = 0; i < 3; ++i) {
                const f = (u: number) => evalDerivative(all[i], u, 0);
                const d1 = (f(t + h) - f(t - h)) / (2 * h);
                const d2 = (f(t + h) - 2 * f(t) + f(t - h)) / (h * h);
                const d3 = (f(t + 2 * h) - 2 * f(t + h) + 2 * f(t - h)
                    - f(t - 2 * h)) / (2 * h * h * h);
                expect(jet[0].values[i]).toBeCloseTo(f(t), 12);
                expect(jet[1].values[i]).toBeCloseTo(d1, 6);
                expect(jet[2].values[i]).toBeCloseTo(d2, 4);
                expect(jet[3].values[i]).toBeCloseTo(d3, 2);
            }
        }
    });

    it('agrees with direct differentiation over random polynomials', () => {
        const rand = makeRandom(20260901);
        for (let trial = 0; trial < 200; ++trial) {
            const dim = 1 + Math.floor(rand() * 4);
            const coefficients: number[][] = [];
            const components: Polynomial1[] = [];
            for (let i = 0; i < dim; ++i) {
                const degree = Math.floor(rand() * 6);
                const c = new Array<number>(degree + 1);
                for (let k = 0; k <= degree; ++k) {
                    c[k] = 2 * rand() - 1;
                }
                // Keep the leading coefficient away from zero so the degree
                // is what we asked for.
                c[degree] += (c[degree] >= 0 ? 0.5 : -0.5);
                coefficients.push(c);
                components.push(poly(...c));
            }

            const curve = new PolynomialCurve(dim, 0, 1, components);
            const t = 2 * rand() - 1;
            const jet = jetOf(curve, t, 3);
            for (let i = 0; i < dim; ++i) {
                for (let d = 0; d <= 3; ++d) {
                    expect(jet[d].values[i]).toBeCloseTo(
                        evalDerivative(coefficients[i], t, d), 10);
                }
            }
        }
    });

    it('has zero curvature and constant tangent for a linear curve', () => {
        // X(t) = (1, 2, 3) + t * (2, -1, 2), speed 3.
        const curve = new PolynomialCurve(3, 0, 4,
            [poly(1, 2), poly(2, -1), poly(3, 2)]);

        for (const t of [0, 1.25, 4]) {
            const jet = jetOf(curve, t, 3);
            expect(vectorLength(jet[2])).toBe(0);
            expect(vectorLength(jet[3])).toBe(0);
            expect(curve.getSpeed(t)).toBeCloseTo(3, 14);

            const tangent = curve.getTangent(t);
            expect(tangent.values[0]).toBeCloseTo(2 / 3, 14);
            expect(tangent.values[1]).toBeCloseTo(-1 / 3, 14);
            expect(tangent.values[2]).toBeCloseTo(2 / 3, 14);
        }

        // Curvature |X' x X''| / |X'|^3 is zero because X'' is zero; the
        // tangential/normal decomposition of X'' is trivially zero as well.
        const jet = jetOf(curve, 2, 2);
        const speedSqr = dot(jet[1], jet[1]);
        const curvatureNumerator = dot(jet[2], jet[2]) * speedSqr
            - dot(jet[1], jet[2]) * dot(jet[1], jet[2]);
        expect(curvatureNumerator).toBe(0);
    });

    it('arclength of a straight-line curve equals the chord length', () => {
        const curve = new PolynomialCurve(3, 0, 4,
            [poly(1, 2), poly(2, -1), poly(3, 2)]);

        const p0 = curve.getPosition(0);
        const p1 = curve.getPosition(4);
        expect(vectorLength(sub(p1, p0))).toBeCloseTo(12, 12);
        expect(curve.getLength(0, 4)).toBeCloseTo(12, 10);
        expect(curve.getTotalLength()).toBeCloseTo(12, 10);
        expect(curve.getLength(1, 3)).toBeCloseTo(6, 10);

        // Arc-length reparameterization: t = Length^{-1}(s) = s / 3.
        expect(curve.getTime(6)).toBeCloseTo(2, 8);
        expect(curve.getTime(0)).toBe(0);
        expect(curve.getTime(20)).toBe(4);
    });

    it('respects the time-domain bounds', () => {
        const curve = new PolynomialCurve(2, -1, 3, [poly(0, 1), poly(0, 0, 1)]);
        expect(curve.getTMin()).toBe(-1);
        expect(curve.getTMax()).toBe(3);
        expect(curve.getNumSegments()).toBe(1);

        // getLength clamps to [tmin, tmax].
        expect(curve.getLength(-10, 10)).toBeCloseTo(curve.getLength(-1, 3), 10);

        // getTime returns values in [tmin, tmax].
        const total = curve.getTotalLength();
        expect(curve.getTime(-1)).toBe(-1);
        expect(curve.getTime(2 * total)).toBe(3);
        const tmid = curve.getTime(0.5 * total);
        expect(tmid).toBeGreaterThanOrEqual(-1);
        expect(tmid).toBeLessThanOrEqual(3);
        expect(curve.getLength(-1, tmid)).toBeCloseTo(0.5 * total, 6);

        // The domain can be reset (single-segment form).
        curve.setTimeInterval(0, 1);
        expect(curve.getTMin()).toBe(0);
        expect(curve.getTMax()).toBe(1);
    });

    it('handles degenerate constant polynomials', () => {
        const curve = new PolynomialCurve(2, 0, 1,
            [poly(7), Polynomial1.fromCoefficients([-3])]);

        for (const t of [0, 0.5, 1]) {
            const jet = jetOf(curve, t, 3);
            expect(jet[0].values[0]).toBe(7);
            expect(jet[0].values[1]).toBe(-3);
            for (let k = 1; k <= 3; ++k) {
                expect(jet[k].values[0]).toBe(0);
                expect(jet[k].values[1]).toBe(0);
            }
        }

        // A constant curve has zero speed and zero length. normalize of the
        // zero tangent yields the zero vector (the Vector convention).
        expect(curve.getSpeed(0.25)).toBe(0);
        expect(curve.getLength(0, 1)).toBe(0);
        const tangent = curve.getTangent(0.5);
        expect(tangent.values[0]).toBe(0);
        expect(tangent.values[1]).toBe(0);
    });

    it('a degree-0 polynomial keeps degree-0 zero derivatives', () => {
        const curve = new PolynomialCurve(1, 0, 1, [poly(9)]);
        expect(curve.getDer1Polynomial(0).getDegree()).toBe(0);
        expect(curve.getDer2Polynomial(0).getDegree()).toBe(0);
        expect(curve.getDer3Polynomial(0).getDegree()).toBe(0);
        expect(curve.getDer3Polynomial(0).get(0)).toBe(0);
    });

    it('rejects a mismatched number of components and bad indices', () => {
        // Upstream enforces components.size() == N at compile time via
        // std::array; the port asserts at runtime.
        expect(() => new PolynomialCurve(3, 0, 1, [poly(1), poly(2)]))
            .toThrow(/Invalid number of components/);
        expect(() => new PolynomialCurve(2, 0, 1,
            [poly(1), poly(2), poly(3)]))
            .toThrow(/Invalid number of components/);
        expect(() => new PolynomialCurve(0, 0, 1)).toThrow(/Invalid dimension/);

        const curve = new PolynomialCurve(2, 0, 1);
        expect(() => curve.setPolynomial(2, poly(1))).toThrow(/Invalid index/);
        expect(() => curve.setPolynomial(-1, poly(1))).toThrow(/Invalid index/);
    });

    it('subdivides by time and by length', () => {
        const curve = new PolynomialCurve(2, 0, 2, [poly(0, 3), poly(0, 4)]);

        const byTime = curve.subdivideByTime(5);
        expect(byTime.length).toBe(5);
        for (let i = 0; i < 5; ++i) {
            const t = (2 * i) / 4;
            expect(byTime[i].values[0]).toBeCloseTo(3 * t, 12);
            expect(byTime[i].values[1]).toBeCloseTo(4 * t, 12);
        }

        // Total length is 5 * 2 = 10; equal-length samples are equally
        // spaced in t for a straight line.
        expect(curve.getTotalLength()).toBeCloseTo(10, 10);
        const byLength = curve.subdivideByLength(5);
        expect(byLength.length).toBe(5);
        for (let i = 0; i < 5; ++i) {
            const t = (2 * i) / 4;
            expect(byLength[i].values[0]).toBeCloseTo(3 * t, 6);
            expect(byLength[i].values[1]).toBeCloseTo(4 * t, 6);
        }
    });

    it('getPosition returns a copy of the jet entry', () => {
        const curve = new PolynomialCurve(2, 0, 1, [poly(1, 1), poly(2, -1)]);
        const p = curve.getPosition(0.5);
        p.values[0] = 999;
        const q = curve.getPosition(0.5);
        expect(q.values[0]).toBeCloseTo(1.5, 14);
        expect(q.values[1]).toBeCloseTo(1.5, 14);
    });
});
