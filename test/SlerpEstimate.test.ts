import { describe, it, expect } from 'vitest';
import {
    slerpEstimate, slerpEstimateUsingCosAngle, slerpEstimateUsingMidpoint
} from '../src/SlerpEstimate.js';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../src/Slerp.js';
import {
    chebyshevRatioEstimate, getChebyshevRatioEstimateMaxError
} from '../src/ChebyshevRatioEstimate.js';
import { check, fc, scaled, unitVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification (V23): independent review against upstream SlerpEstimate.h.
// ---------------------------------------------------------------------------

describe('SlerpEstimate verification', () => {
    const unitArray = (n: number): fc.Arbitrary<number[]> =>
        unitVector(n).map(v => [...v.values]);
    // The first two estimates require an angle in [0,pi/2], i.e. a
    // non-negative dot product; upstream's documented preprocessing flips the
    // sign of q1 to guarantee it.
    const pairArb = fc.tuple(unitArray(4), unitArray(4))
        .map(([a, b]) => (dot(a, b) < 0 ? [a, b.map(v => -v)] : [a, b]))
        .filter(([a, b]) => dot(a, b) < 0.999);
    const tArb = scaled(0, 1);

    it('is exactly the Chebyshev estimate weighting of the two inputs', () => {
        // Pins the port's accumulation loop against the file it delegates to,
        // including the choice of q0/q1 order in each branch.
        for (const degree of [1, 5, 16]) {
            check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
                const f = chebyshevRatioEstimate(t, dot(q0, q1), degree);
                const r = slerpEstimate(t, q0, q1, degree);
                for (let i = 0; i < 4; ++i) {
                    if (r[i] + 0 !== f[0] * q0[i] + f[1] * q1[i] + 0) {
                        return false;
                    }
                }
                return true;
            });
        }
    });

    it('tracks exact slerp within twice the documented ratio error', () => {
        // Each component is f[0]*q0[i] + f[1]*q1[i] with |q| <= 1, so the two
        // ratio errors simply add.
        for (const degree of [1, 2, 4, 8, 12, 16]) {
            check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
                const exact = slerp(t, q0, q1);
                const est = slerpEstimate(t, q0, q1, degree);
                return maxAbsDiff(exact, est) <= componentBound(degree);
            });
        }
    });

    it('improves monotonically with the degree', () => {
        // Averaged over a deterministic sample: a swapped u-table or an
        // off-by-one degree index would break the ordering.
        const rand = makeRandom(20230423);
        const samples: { q0: number[]; q1: number[]; t: number }[] = [];
        for (let i = 0; i < 400; ++i) {
            const { q0, q1 } = randomPair(rand);
            samples.push({ q0, q1, t: rand() });
        }
        let previous = Number.POSITIVE_INFINITY;
        for (let degree = 1; degree <= 16; ++degree) {
            let worst = 0;
            for (const s of samples) {
                worst = Math.max(worst, maxAbsDiff(slerp(s.t, s.q0, s.q1),
                    slerpEstimate(s.t, s.q0, s.q1, degree)));
            }
            expect(worst).toBeLessThanOrEqual(componentBound(degree));
            expect(worst).toBeLessThan(previous);
            previous = worst;
        }
    }, 30000);

    it('stays close to the unit hypersphere', () => {
        // The exact slerp is unit length, so the estimate cannot deviate by
        // more than the componentwise bound times sqrt(N).
        for (const degree of [4, 10, 16]) {
            check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
                const r = slerpEstimate(t, q0, q1, degree);
                return Math.abs(Math.sqrt(dot(r, r)) - 1)
                    <= 2 * componentBound(degree);
            });
        }
    });

    it('is exactly the precomputed-cosine overload for the exact dot', () => {
        for (const degree of [3, 16]) {
            check(fc.tuple(pairArb, tArb), ([[q0, q1], t]) => {
                const a = slerpEstimate(t, q0, q1, degree);
                const b = slerpEstimateUsingCosAngle(t, q0, q1,
                    dot(q0, q1), degree);
                for (let i = 0; i < 4; ++i) {
                    if (a[i] + 0 !== b[i] + 0) { return false; }
                }
                return true;
            });
        }
    });

    it('handles angles beyond pi/2 through the midpoint overload', () => {
        // The plain estimate is only documented for [0,pi/2]; the midpoint
        // form halves the angle, so it covers [0,pi). This checks the branch
        // switch at t = 1/2 uses q0/qh below and qh/q1 above.
        const widePair = fc.tuple(unitArray(4), unitArray(4))
            .filter(([a, b]) => {
                const c = dot(a, b);
                return c > -0.98 && c < 0.999;
            });
        for (const degree of [8, 16]) {
            check(fc.tuple(widePair, tArb), ([[q0, q1], t]) => {
                const cosA = dot(q0, q1);
                const cosAH = Math.sqrt((1 + cosA) / 2);
                const qh = q0.map((v, i) => (v + q1[i]) / (2 * cosAH));
                const est = slerpEstimateUsingMidpoint(t, q0, q1, qh,
                    cosAH, degree);
                const exact = slerp(t, q0, q1);
                return maxAbsDiff(exact, est) <= componentBound(degree) + 1e-9;
            });
        }
    });

    it('reproduces the endpoints to the accuracy of the ratio estimate', () => {
        // The zero weight is exact (t = 0 kills term1 identically), while the
        // other weight is the estimate of f(1,x) = 1.
        for (const degree of [1, 8, 16]) {
            check(pairArb, ([q0, q1]) => {
                const at0 = slerpEstimate(0, q0, q1, degree);
                const at1 = slerpEstimate(1, q0, q1, degree);
                return maxAbsDiff(at0, q0) <= componentBound(degree)
                    && maxAbsDiff(at1, q1) <= componentBound(degree);
            });
        }
    });

    it('rejects mismatched dimensions and invalid degrees', () => {
        expect(() => slerpEstimate(0.5, [1], [1], 4))
            .toThrow('Invalid dimension.');
        expect(() => slerpEstimate(0.5, [1, 0], [1, 0, 0], 4))
            .toThrow('Mismatched dimensions.');
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 0))
            .toThrow('Invalid degree.');
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 17))
            .toThrow('Invalid degree.');
        expect(() => slerpEstimate(0.5, [1, 0], [0, 1], 2.5))
            .toThrow('Invalid degree.');
        expect(() => slerpEstimateUsingMidpoint(0.5, [1, 0], [0, 1],
            [1, 0, 0], 1, 4)).toThrow('Mismatched dimensions.');
    });
});
