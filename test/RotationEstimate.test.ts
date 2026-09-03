import { describe, it, expect } from 'vitest';
import { GTE_C_PI } from '../src/Constants.js';
import { Matrix, addMatrix, mulMatrix } from '../src/Matrix.js';
import {
    rotC0Estimate, rotC1Estimate, rotC2Estimate, rotC3Estimate,
    rotC4Estimate, getRotC0EstimateMaxError, getRotC1EstimateMaxError,
    getRotC2EstimateMaxError, getRotC3EstimateMaxError,
    getRotC4EstimateMaxError, rotationEstimate, rotationDerivativeEstimate,
    rotationAndDerivativeEstimate
} from '../src/RotationEstimate.js';
import { Vector, length } from '../src/Vector.js';

const DEGREES = [4, 6, 8, 10, 12, 14, 16] as const;

// ---------------------------------------------------------------------------
// Numerically stable references for the five rotc functions. The closed forms
// (1-cos(t))/t^2, (2*(1-cos(t))-t*sin(t))/t^4, ... suffer catastrophic
// cancellation as t -> 0 (the t^4 denominator loses ~16 digits by t = 1e-4),
// so the Maclaurin series is used on [0,1) where it converges in a handful of
// terms and the closed form is used on [1,pi].
// ---------------------------------------------------------------------------

function factorial(n: number): number {
    let f = 1;
    for (let i = 2; i <= n; ++i) {
        f *= i;
    }
    return f;
}

function evenSeries(term: (k: number) => number, kStart: number,
    kEnd: number): number {
    let sum = 0;
    for (let k = kStart; k <= kEnd; ++k) {
        sum += term(k);
    }
    return sum;
}

// rotc0(t) = sin(t)/t = sum_{k>=0} (-1)^k t^(2k)/(2k+1)!
function exactRotC0(t: number): number {
    if (t >= 1) {
        return Math.sin(t) / t;
    }
    return evenSeries(k => (k % 2 === 0 ? 1 : -1) * Math.pow(t, 2 * k)
        / factorial(2 * k + 1), 0, 20);
}

// rotc1(t) = (1-cos(t))/t^2 = sum_{k>=0} (-1)^k t^(2k)/(2k+2)!
function exactRotC1(t: number): number {
    if (t >= 1) {
        return (1 - Math.cos(t)) / (t * t);
    }
    return evenSeries(k => (k % 2 === 0 ? 1 : -1) * Math.pow(t, 2 * k)
        / factorial(2 * k + 2), 0, 20);
}

// rotc2(t) = (sin(t) - t*cos(t))/t^3
//          = sum_{k>=1} (-1)^(k+1) 2k t^(2k-2)/(2k+1)!
function exactRotC2(t: number): number {
    if (t >= 1) {
        return (Math.sin(t) - t * Math.cos(t)) / (t * t * t);
    }
    return evenSeries(k => ((k + 1) % 2 === 0 ? 1 : -1) * 2 * k
        * Math.pow(t, 2 * k - 2) / factorial(2 * k + 1), 1, 20);
}

// rotc3(t) = (2*(1-cos(t)) - t*sin(t))/t^4
//          = sum_{m>=2} (-1)^(m+1) (2-2m) t^(2m-4)/(2m)!
function exactRotC3(t: number): number {
    if (t >= 1) {
        return (2 * (1 - Math.cos(t)) - t * Math.sin(t)) / (t * t * t * t);
    }
    return evenSeries(m => ((m + 1) % 2 === 0 ? 1 : -1) * (2 - 2 * m)
        * Math.pow(t, 2 * m - 4) / factorial(2 * m), 2, 20);
}

// rotc4(t) = (t - sin(t))/t^3 = sum_{k>=1} (-1)^(k+1) t^(2k-2)/(2k+1)!
function exactRotC4(t: number): number {
    if (t >= 1) {
        return (t - Math.sin(t)) / (t * t * t);
    }
    return evenSeries(k => ((k + 1) % 2 === 0 ? 1 : -1)
        * Math.pow(t, 2 * k - 2) / factorial(2 * k + 1), 1, 20);
}

// The upstream-documented maximum errors on [0,pi] (C_ROTC*_EST_MAX_ERROR).
const TABULATED_MAX_ERROR: Record<string, Record<number, number>> = {
    rotc0: {
        4: 6.9656371186750e-03, 6: 2.2379506089580e-04,
        8: 4.8670096434722e-06, 10: 7.5654711606532e-08,
        12: 8.7939167753293e-10, 14: 1.8030021919913e-12,
        16: 6.8001160258291e-16
    },
    rotc1: {
        4: 9.2119010150538e-04, 6: 2.3251261806301e-05,
        8: 4.1693160884870e-07, 10: 5.5177887814395e-09,
        12: 5.5865700954172e-11, 14: 7.1609385088323e-15,
        16: 7.2164496600635e-16
    },
    rotc2: {
        4: 8.1461508460229e-04, 6: 2.1075025784856e-05,
        8: 3.8414838612888e-07, 10: 5.1435966597069e-09,
        12: 5.2533449812486e-11, 14: 7.7715611723761e-15,
        16: 2.2759572004816e-15
    },
    rotc3: {
        4: 8.4612036888886e-05, 6: 1.8051973185995e-06,
        8: 2.8016103950645e-08, 10: 3.2675391559156e-10,
        12: 1.3714029911682e-13, 14: 3.2078506517763e-14,
        16: 4.7774284528401e-14
    },
    rotc4: {
        4: 2.9118443863947e-06, 6: 1.1775194543557e-08,
        8: 3.3551378409470e-11, 10: 7.0975185732742e-14,
        12: 1.1586733901401e-16, 14: 1.5038700542518e-19,
        16: 1.5889884416916e-22
    }
};

// The true maxima of |estimate(t) - exact(t)| over [0,pi], computed offline
// with 60-decimal-digit arithmetic on a 4001-point grid. For most entries
// these agree with the upstream table to all printed digits; the entries that
// do not are the upstream defect reported in the PR (see the
// 'documents where the upstream max-error table is wrong' test below).
const TRUE_MAX_ERROR: Record<string, Record<number, number>> = {
    rotc0: {
        4: 6.96564e-03, 6: 2.23795e-04, 8: 4.86701e-06, 10: 7.56547e-08,
        12: 8.79395e-10, 14: 1.62275e-11, 16: 6.82690e-13
    },
    rotc1: {
        4: 9.21190e-04, 6: 2.32513e-05, 8: 4.16931e-07, 10: 5.51783e-09,
        12: 5.58687e-11, 14: 3.61721e-12, 16: 2.11557e-14
    },
    rotc2: {
        4: 8.14615e-04, 6: 2.10750e-05, 8: 3.84148e-07, 10: 5.14360e-09,
        12: 5.25398e-11, 14: 3.43280e-12, 16: 1.02106e-14
    },
    rotc3: {
        4: 8.46120e-05, 6: 1.80520e-06, 8: 2.80161e-08, 10: 3.26762e-10,
        12: 1.64047e-11, 14: 1.08672e-13, 16: 8.44349e-13
    },
    rotc4: {
        4: 1.76814e-03, 6: 1.48047e-04, 8: 8.65572e-06, 10: 3.73868e-07,
        12: 1.24258e-08, 14: 2.26397e-07, 16: 7.02836e-12
    }
};

const FUNCTIONS: ReadonlyArray<{
    name: string,
    estimate: (t: number, degree: number) => number,
    exact: (t: number) => number,
    maxError: (degree: number) => number,
    valueAtZero: number
}> = [
    {
        name: 'rotc0', estimate: rotC0Estimate, exact: exactRotC0,
        maxError: getRotC0EstimateMaxError, valueAtZero: 1
    },
    {
        name: 'rotc1', estimate: rotC1Estimate, exact: exactRotC1,
        maxError: getRotC1EstimateMaxError, valueAtZero: 0.5
    },
    {
        name: 'rotc2', estimate: rotC2Estimate, exact: exactRotC2,
        maxError: getRotC2EstimateMaxError, valueAtZero: 1 / 3
    },
    {
        name: 'rotc3', estimate: rotC3Estimate, exact: exactRotC3,
        maxError: getRotC3EstimateMaxError, valueAtZero: 1 / 12
    },
    {
        name: 'rotc4', estimate: rotC4Estimate, exact: exactRotC4,
        maxError: getRotC4EstimateMaxError, valueAtZero: 1 / 6
    }
];

function observedMaxError(f: (t: number) => number,
    exact: (t: number) => number, samples: number = 4000): number {
    let maxObserved = 0;
    for (let i = 0; i <= samples; ++i) {
        const t = (GTE_C_PI * i) / samples;
        const error = Math.abs(f(t) - exact(t));
        if (error > maxObserved) {
            maxObserved = error;
        }
    }
    return maxObserved;
}

describe('rotC*Estimate', () => {
    it('matches the exact functions to the verified minimax error on [0,pi]',
        () => {
            for (const f of FUNCTIONS) {
                for (const degree of DEGREES) {
                    const observed = observedMaxError(
                        t => f.estimate(t, degree), f.exact);
                    const trueBound = TRUE_MAX_ERROR[f.name][degree];
                    // The 2e-15 slack covers double-precision evaluation of
                    // the polynomial and of the reference function.
                    expect(observed,
                        `${f.name} degree ${degree}`)
                        .toBeLessThanOrEqual(1.05 * trueBound + 2e-15);
                    // The bound is attained (it is a minimax error), so a
                    // dense grid must come close to it.
                    expect(observed, `${f.name} degree ${degree} tightness`)
                        .toBeGreaterThan(0.4 * trueBound - 2e-15);
                }
            }
        });

    it('respects the upstream-tabulated bound wherever that bound is correct',
        () => {
            // Degrees 4..10 (and 12 except for rotc3) of rotc0..rotc3 are the
            // entries of C_ROTC*_EST_MAX_ERROR that really do bound the
            // approximation error.
            const correct: Record<string, number[]> = {
                rotc0: [4, 6, 8, 10, 12],
                rotc1: [4, 6, 8, 10, 12],
                rotc2: [4, 6, 8, 10, 12],
                rotc3: [4, 6, 8, 10]
            };
            for (const f of FUNCTIONS) {
                const degrees = correct[f.name];
                if (degrees === undefined) {
                    continue;
                }
                for (const degree of degrees) {
                    const observed = observedMaxError(
                        t => f.estimate(t, degree), f.exact);
                    expect(observed, `${f.name} degree ${degree}`)
                        .toBeLessThanOrEqual(
                            TABULATED_MAX_ERROR[f.name][degree] * (1 + 1e-3));
                }
            }
        });

    it('documents where the upstream max-error table is wrong', () => {
        // These entries of C_ROTC*_EST_MAX_ERROR understate the true error of
        // the tabulated coefficients on [0,pi]; the whole rotc4 table is
        // wrong. See the PR's "Upstream bug suspects" section. The port keeps
        // the upstream numbers verbatim (they are documentation, not used in
        // any computation), so this test pins the discrepancy.
        const understated: Array<[string, number]> = [
            ['rotc0', 14], ['rotc0', 16],
            ['rotc1', 14], ['rotc1', 16],
            ['rotc2', 14], ['rotc2', 16],
            ['rotc3', 12], ['rotc3', 14], ['rotc3', 16],
            ['rotc4', 4], ['rotc4', 6], ['rotc4', 8], ['rotc4', 10],
            ['rotc4', 12], ['rotc4', 14], ['rotc4', 16]
        ];
        for (const [name, degree] of understated) {
            const f = FUNCTIONS.find(g => g.name === name)!;
            const observed = observedMaxError(t => f.estimate(t, degree),
                f.exact);
            expect(observed, `${name} degree ${degree}`)
                .toBeGreaterThan(2 * TABULATED_MAX_ERROR[name][degree]);
        }
    });

    it('reports the upstream-tabulated max errors verbatim', () => {
        for (const f of FUNCTIONS) {
            for (const degree of DEGREES) {
                expect(f.maxError(degree))
                    .toBe(TABULATED_MAX_ERROR[f.name][degree]);
            }
        }
    });

    it('reproduces the limit values at t = 0', () => {
        for (const f of FUNCTIONS) {
            for (const degree of DEGREES) {
                // The estimate at t = 0 is coeff[0] exactly.
                expect(Math.abs(f.estimate(0, degree) - f.valueAtZero))
                    .toBeLessThanOrEqual(
                        Math.max(TRUE_MAX_ERROR[f.name][degree], 1e-15));
            }
        }
    });

    it('is an even function (only even powers of t appear)', () => {
        for (const f of FUNCTIONS) {
            for (const degree of DEGREES) {
                for (const t of [0.25, 1, 2, GTE_C_PI]) {
                    expect(f.estimate(-t, degree)).toBe(f.estimate(t, degree));
                }
            }
        }
    });

    it('throws for invalid degrees', () => {
        for (const f of FUNCTIONS) {
            expect(() => f.estimate(0.5, 3)).toThrow('Invalid degree.');
            expect(() => f.estimate(0.5, 2)).toThrow('Invalid degree.');
            expect(() => f.estimate(0.5, 18)).toThrow('Invalid degree.');
            expect(() => f.estimate(0.5, 6.5)).toThrow('Invalid degree.');
            expect(() => f.maxError(5)).toThrow('Invalid degree.');
            expect(() => f.maxError(20)).toThrow('Invalid degree.');
        }
    });
});

// ---------------------------------------------------------------------------
// Rotation matrix and derivative estimates.
// ---------------------------------------------------------------------------

function skew(p: Vector): Matrix {
    return Matrix.fromArray(3, 3, [
        0, -p.values[2], p.values[1],
        p.values[2], 0, -p.values[0],
        -p.values[1], p.values[0], 0
    ]);
}

// The exact R = exp(S) = I + rotc0(t)*S + rotc1(t)*S^2, using the stable
// references above.
function exactRotation(p: Vector): Matrix {
    const S = skew(p);
    const Ssqr = mulMatrix(S, S) as Matrix;
    const t = length(p);
    return addMatrix(addMatrix(Matrix.identity(3, 3),
        mulMatrix(exactRotC0(t), S) as Matrix),
        mulMatrix(exactRotC1(t), Ssqr) as Matrix);
}

function maxAbsDifference(A: Matrix, B: Matrix): number {
    let m = 0;
    for (let i = 0; i < A.numElements; ++i) {
        m = Math.max(m, Math.abs(A.values[i] - B.values[i]));
    }
    return m;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// A uniformly-distributed rotation vector with |p| <= pi.
function randomRotationVector(rand: () => number): Vector {
    const z = 2 * rand() - 1;
    const phi = 2 * GTE_C_PI * rand();
    const r = Math.sqrt(Math.max(0, 1 - z * z));
    const t = GTE_C_PI * rand();
    return Vector.fromArray([t * r * Math.cos(phi), t * r * Math.sin(phi),
        t * z]);
}

describe('rotationEstimate', () => {
    it('is the identity for the zero vector at every degree', () => {
        for (const degree of DEGREES) {
            const R = rotationEstimate(new Vector(3), degree);
            expect(maxAbsDifference(R, Matrix.identity(3, 3))).toBe(0);
        }
    });

    it('reproduces a known 90-degree rotation about z', () => {
        // p = (0,0,pi/2) gives R = {{0,-1,0},{1,0,0},{0,0,1}} for MAT_VEC.
        const p = Vector.fromArray([0, 0, GTE_C_PI / 2]);
        const R = rotationEstimate(p, 16);
        const expected = Matrix.fromArray(3, 3, [0, -1, 0, 1, 0, 0, 0, 0, 1]);
        expect(maxAbsDifference(R, expected)).toBeLessThanOrEqual(1e-12);
    });

    it('converges to the exact Rodrigues rotation as the degree grows', () => {
        const rand = makeRandom(20260901);
        const worst: Record<number, number> = {};
        for (const degree of DEGREES) {
            worst[degree] = 0;
        }
        for (let trial = 0; trial < 200; ++trial) {
            const p = randomRotationVector(rand);
            const exact = exactRotation(p);
            for (const degree of DEGREES) {
                const error = maxAbsDifference(rotationEstimate(p, degree),
                    exact);
                worst[degree] = Math.max(worst[degree], error);
            }
        }
        // |R_est - R| <= t^2 * (|rotc0 error| + |rotc1 error| * t) is bounded
        // by a small multiple of pi^2 times the rotc0/rotc1 errors.
        for (const degree of DEGREES) {
            const bound = 20 * (TRUE_MAX_ERROR.rotc0[degree]
                + TRUE_MAX_ERROR.rotc1[degree]) + 1e-14;
            expect(worst[degree], `degree ${degree}`)
                .toBeLessThanOrEqual(bound);
        }
        // The estimates really do get better with the degree.
        expect(worst[16]).toBeLessThan(worst[4]);
        expect(worst[10]).toBeLessThan(worst[6]);
    });

    it('produces nearly orthonormal matrices with determinant 1', () => {
        const rand = makeRandom(777);
        for (let trial = 0; trial < 100; ++trial) {
            const p = randomRotationVector(rand);
            const R = rotationEstimate(p, 16);
            const RtR = mulMatrix(
                Matrix.fromArray(3, 3, [
                    R.get(0, 0), R.get(1, 0), R.get(2, 0),
                    R.get(0, 1), R.get(1, 1), R.get(2, 1),
                    R.get(0, 2), R.get(1, 2), R.get(2, 2)
                ]), R) as Matrix;
            expect(maxAbsDifference(RtR, Matrix.identity(3, 3)))
                .toBeLessThanOrEqual(1e-11);
        }
    });

    it('throws for non-3D input', () => {
        expect(() => rotationEstimate(new Vector(2), 8))
            .toThrow('RotationEstimate: expecting a 3-tuple.');
        expect(() => rotationDerivativeEstimate(new Vector(4), 8))
            .toThrow('RotationEstimate: expecting a 3-tuple.');
        expect(() => rotationAndDerivativeEstimate(new Vector(4), 8))
            .toThrow('RotationEstimate: expecting a 3-tuple.');
    });
});

describe('rotationDerivativeEstimate', () => {
    it('matches central finite differences of the exact rotation', () => {
        const rand = makeRandom(4242);
        const h = 1e-5;
        let worst = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const p = randomRotationVector(rand);
            if (length(p) < 0.05) {
                continue;  // finite differences are ill-conditioned here
            }
            const Rder = rotationDerivativeEstimate(p, 16);
            for (let i = 0; i < 3; ++i) {
                const pPlus = p.clone();
                pPlus.values[i] += h;
                const pMinus = p.clone();
                pMinus.values[i] -= h;
                const Rp = exactRotation(pPlus);
                const Rm = exactRotation(pMinus);
                for (let e = 0; e < 9; ++e) {
                    const fd = (Rp.values[e] - Rm.values[e]) / (2 * h);
                    worst = Math.max(worst,
                        Math.abs(fd - Rder[i].values[e]));
                }
            }
        }
        // Central differences with h = 1e-5 carry O(h^2) ~ 1e-10 truncation
        // error plus O(eps/h) ~ 1e-11 roundoff.
        expect(worst).toBeLessThanOrEqual(1e-8);
    });

    it('gives skew-symmetric generators at p = 0', () => {
        // At p = 0 the derivative reduces to a*E[i] with a = rotc0(0) = 1.
        const Rder = rotationDerivativeEstimate(new Vector(3), 16);
        expect(maxAbsDifference(Rder[0],
            Matrix.fromArray(3, 3, [0, 0, 0, 0, 0, -1, 0, 1, 0])))
            .toBeLessThanOrEqual(1e-15);
        expect(maxAbsDifference(Rder[1],
            Matrix.fromArray(3, 3, [0, 0, 1, 0, 0, 0, -1, 0, 0])))
            .toBeLessThanOrEqual(1e-15);
        expect(maxAbsDifference(Rder[2],
            Matrix.fromArray(3, 3, [0, -1, 0, 1, 0, 0, 0, 0, 0])))
            .toBeLessThanOrEqual(1e-15);
    });

    it('is consistent with rotationAndDerivativeEstimate', () => {
        const rand = makeRandom(31337);
        for (let trial = 0; trial < 40; ++trial) {
            const p = randomRotationVector(rand);
            for (const degree of [4, 10, 16]) {
                const { R, Rder } = rotationAndDerivativeEstimate(p, degree);
                expect(R.values).toEqual(rotationEstimate(p, degree).values);
                const separate = rotationDerivativeEstimate(p, degree);
                for (let i = 0; i < 3; ++i) {
                    expect(Rder[i].values).toEqual(separate[i].values);
                }
            }
        }
    });
});
