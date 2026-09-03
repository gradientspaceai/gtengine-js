import { describe, it, expect } from 'vitest';
import {
    slerpEstimate, slerpEstimateUsingCosAngle, slerpEstimateUsingMidpoint
} from '../src/SlerpEstimate.js';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../src/Slerp.js';
import { getChebyshevRatioEstimateMaxError } from '../src/ChebyshevRatioEstimate.js';

const DEGREES = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16] as const;

function dot(a: number[], b: number[]): number {
    let s = 0;
    for (let i = 0; i < a.length; ++i) {
        s += a[i] * b[i];
    }
    return s;
}

function normalize(v: number[]): number[] {
    const len = Math.sqrt(dot(v, v));
    return v.map((x) => x / len);
}

// Exact slerp reference on the great arc through q0 and q1.
function referenceSlerp(t: number, q0: number[], q1: number[]): number[] {
    const angle = Math.acos(Math.min(1, Math.max(-1, dot(q0, q1))));
    const sinA = Math.sin(angle);
    if (sinA === 0) {
        return q0.map((v, i) => (1 - t) * v + t * q1[i]);
    }
    const f0 = Math.sin((1 - t) * angle) / sinA;
    const f1 = Math.sin(t * angle) / sinA;
    return q0.map((v, i) => f0 * v + f1 * q1[i]);
}

function maxAbsDiff(a: number[], b: number[]): number {
    let m = 0;
    for (let i = 0; i < a.length; ++i) {
        m = Math.max(m, Math.abs(a[i] - b[i]));
    }
    return m;
}

// The estimate error of each Chebyshev ratio is bounded by the documented
// value; each component of the result is f[0]*q0[i] + f[1]*q1[i] with
// |q0[i]|, |q1[i]| <= 1, so the componentwise bound is twice that.
function componentBound(degree: number): number {
    return 2 * getChebyshevRatioEstimateMaxError(degree) + 1e-15;
}

// Deterministic pseudo-random generator so failures are reproducible.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// A random unit 4-vector with a non-negative dot product against q0, so the
// angle between them is in [0,pi/2] as the estimates require.
function randomPair(rand: () => number): { q0: number[]; q1: number[] } {
    const q0 = normalize([
        2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
    ]);
    let q1 = normalize([
        2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1
    ]);
    if (dot(q0, q1) < 0) {
        q1 = q1.map((v) => -v);
    }
    return { q0, q1 };
}

describe('slerpEstimate', () => {
    it('stays within the documented bound over angles in [0,pi/2]', () => {
        for (const degree of DEGREES) {
            const bound = componentBound(degree);
            for (let a = 0; a <= 32; ++a) {
                const angle = (a / 32) * (Math.PI / 2);
                const q0 = [1, 0, 0, 0];
                const q1 = [Math.cos(angle), Math.sin(angle), 0, 0];
                for (let i = 0; i <= 16; ++i) {
                    const t = i / 16;
                    const estimate = slerpEstimate(t, q0, q1, degree);
                    const exact = referenceSlerp(t, q0, q1);
                    expect(maxAbsDiff(estimate, exact)).toBeLessThanOrEqual(bound);
                }
            }
        }
    });

    it('agrees with exact slerp for random pairs at high degree', () => {
        const rand = makeRandom(31415);
        for (let k = 0; k < 150; ++k) {
            const { q0, q1 } = randomPair(rand);
            if (dot(q0, q1) >= 1) {
                continue;
            }
            for (const t of [0.05, 0.25, 0.5, 0.75, 0.95]) {
                expect(maxAbsDiff(slerpEstimate(t, q0, q1, 16), slerp(t, q0, q1)))
                    .toBeLessThanOrEqual(componentBound(16));
            }
        }
    });

    it('reproduces the endpoints exactly for every degree', () => {
        // At t = 0 the ratios are {1,0} and at t = 1 they are {0,1}, exactly,
        // because every series term carries a factor of t (respectively 1-t).
        const q0 = normalize([1, 2, 3, 4]);
        const q1 = normalize([4, 3, 2, 1]);
        for (const degree of DEGREES) {
            expect(slerpEstimate(0, q0, q1, degree)).toEqual(q0);
            expect(slerpEstimate(1, q0, q1, degree)).toEqual(q1);
        }
    });

    it('produces nearly unit-length results for unit-length inputs', () => {
        const rand = makeRandom(2718);
        for (let k = 0; k < 100; ++k) {
            const { q0, q1 } = randomPair(rand);
            for (const t of [0.2, 0.5, 0.8]) {
                const r = slerpEstimate(t, q0, q1, 12);
                expect(Math.abs(Math.sqrt(dot(r, r)) - 1))
                    .toBeLessThanOrEqual(componentBound(12));
            }
        }
    });

    it('improves with increasing degree', () => {
        const angle = Math.PI / 2;
        const q0 = [1, 0, 0, 0];
        const q1 = [Math.cos(angle), Math.sin(angle), 0, 0];
        let previous = Number.MAX_VALUE;
        for (const degree of [1, 2, 3, 4, 5, 6, 7, 8]) {
            let worst = 0;
            for (let i = 0; i <= 32; ++i) {
                const t = i / 32;
                worst = Math.max(worst, maxAbsDiff(
                    slerpEstimate(t, q0, q1, degree), referenceSlerp(t, q0, q1)));
            }
            expect(worst).toBeLessThan(previous);
            previous = worst;
        }
    });

    it('handles identical inputs (zero angle) as linear interpolation', () => {
        const q0 = [0, 0, 0, 1];
        for (const degree of [1, 4, 16]) {
            for (const t of [0, 0.3, 0.5, 1]) {
                const r = slerpEstimate(t, q0, q0, degree);
                expect(maxAbsDiff(r, q0)).toBeLessThan(1e-15);
            }
        }
    });

    it('works for dimensions other than 4', () => {
        const q0 = [1, 0];
        const q1 = [0, 1];
        expect(maxAbsDiff(slerpEstimate(0.5, q0, q1, 16),
            [Math.SQRT1_2, Math.SQRT1_2])).toBeLessThanOrEqual(componentBound(16));
    });

    it('validates dimension and degree', () => {
        expect(() => slerpEstimate(0.5, [1], [1], 4)).toThrow(/Invalid dimension/);
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1, 0], 4))
            .toThrow(/Mismatched dimensions/);
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 0)).toThrow(/Invalid degree/);
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 17)).toThrow(/Invalid degree/);
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 2.5)).toThrow(/Invalid degree/);
    });
});

describe('slerpEstimateUsingCosAngle', () => {
    it('matches slerpEstimate when given the exact dot product', () => {
        const rand = makeRandom(777);
        for (let k = 0; k < 100; ++k) {
            const { q0, q1 } = randomPair(rand);
            const cosA = dot(q0, q1);
            for (const t of [0, 0.4, 0.5, 1]) {
                expect(slerpEstimateUsingCosAngle(t, q0, q1, cosA, 10))
                    .toEqual(slerpEstimate(t, q0, q1, 10));
            }
        }
    });

    it('stays within the documented bound of exact slerp', () => {
        const q0 = [1, 0, 0, 0];
        for (const degree of DEGREES) {
            const bound = componentBound(degree);
            for (let a = 1; a <= 16; ++a) {
                const angle = (a / 16) * (Math.PI / 2);
                const q1 = [Math.cos(angle), Math.sin(angle), 0, 0];
                const cosA = dot(q0, q1);
                for (let i = 0; i <= 8; ++i) {
                    const t = i / 8;
                    const estimate = slerpEstimateUsingCosAngle(t, q0, q1, cosA, degree);
                    const exact = slerpUsingCosAngle(t, q0, q1, cosA);
                    expect(maxAbsDiff(estimate, exact)).toBeLessThanOrEqual(bound);
                }
            }
        }
    });

    it('validates dimension and degree', () => {
        expect(() => slerpEstimateUsingCosAngle(0.5, [1], [1], 0.5, 4))
            .toThrow(/Invalid dimension/);
        expect(() => slerpEstimateUsingCosAngle(0.5, [1, 0], [0, 1], 0, 17))
            .toThrow(/Invalid degree/);
    });
});

describe('slerpEstimateUsingMidpoint', () => {
    function midpointData(q0: number[], q1: number[]): { qh: number[]; cosAH: number } {
        const cosA = dot(q0, q1);
        const cosAH = Math.sqrt((1 + cosA) / 2);
        const qh = q0.map((v, i) => (v + q1[i]) / (2 * cosAH));
        return { qh, cosAH };
    }

    it('stays within the documented bound for angles in [0,pi)', () => {
        const q0 = [1, 0, 0, 0];
        for (const degree of DEGREES) {
            // Each branch interpolates over a half arc using two unit-length
            // endpoints, so the componentwise bound is the same.
            const bound = componentBound(degree);
            for (let a = 1; a <= 24; ++a) {
                const angle = (a / 25) * Math.PI;
                const q1 = [Math.cos(angle), Math.sin(angle), 0, 0];
                const { qh, cosAH } = midpointData(q0, q1);
                for (let i = 0; i <= 10; ++i) {
                    const t = i / 10;
                    const estimate = slerpEstimateUsingMidpoint(t, q0, q1, qh, cosAH, degree);
                    const exact = [Math.cos(t * angle), Math.sin(t * angle), 0, 0];
                    expect(maxAbsDiff(estimate, exact)).toBeLessThanOrEqual(bound);
                }
            }
        }
    });

    it('agrees with the exact midpoint form at high degree', () => {
        const rand = makeRandom(1618);
        for (let k = 0; k < 100; ++k) {
            const { q0, q1 } = randomPair(rand);
            if (dot(q0, q1) >= 1) {
                continue;
            }
            const { qh, cosAH } = midpointData(q0, q1);
            for (const t of [0.1, 0.5, 0.6, 0.9]) {
                expect(maxAbsDiff(
                    slerpEstimateUsingMidpoint(t, q0, q1, qh, cosAH, 16),
                    slerpUsingMidpoint(t, q0, q1, qh, cosAH)))
                    .toBeLessThanOrEqual(componentBound(16));
            }
        }
    });

    it('reproduces the endpoints and the midpoint exactly', () => {
        const q0 = normalize([1, 0.4, -0.2, 0.7]);
        const q1 = normalize([-0.3, 1, 0.5, 0.1]);
        const { qh, cosAH } = midpointData(q0, q1);
        for (const degree of DEGREES) {
            expect(slerpEstimateUsingMidpoint(0, q0, q1, qh, cosAH, degree)).toEqual(q0);
            // t = 1/2 takes the twoT <= 1 branch with argument 1, giving qh.
            expect(slerpEstimateUsingMidpoint(0.5, q0, q1, qh, cosAH, degree)).toEqual(qh);
            expect(slerpEstimateUsingMidpoint(1, q0, q1, qh, cosAH, degree)).toEqual(q1);
        }
    });

    it('is continuous across the t = 1/2 branch switch', () => {
        const q0 = normalize([1, 0.4, -0.2, 0.7]);
        const q1 = normalize([-0.3, 1, 0.5, 0.1]);
        const { qh, cosAH } = midpointData(q0, q1);
        const below = slerpEstimateUsingMidpoint(0.5 - 1e-9, q0, q1, qh, cosAH, 12);
        const above = slerpEstimateUsingMidpoint(0.5 + 1e-9, q0, q1, qh, cosAH, 12);
        expect(maxAbsDiff(below, above)).toBeLessThan(1e-8);
    });

    it('validates dimensions and degree', () => {
        expect(() => slerpEstimateUsingMidpoint(0.5, [1], [1], [1], 0.9, 4))
            .toThrow(/Invalid dimension/);
        expect(() => slerpEstimateUsingMidpoint(0.5, [1, 0], [0, 1], [1, 0, 0], 0.9, 4))
            .toThrow(/Mismatched dimensions/);
        expect(() => slerpEstimateUsingMidpoint(0.5, [1, 0], [0, 1], [1, 0], 0.9, 0))
            .toThrow(/Invalid degree/);
    });
});
