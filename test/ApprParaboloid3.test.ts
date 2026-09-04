import { describe, it, expect } from 'vitest';
import {
    ApprParaboloid3, apprParaboloid3SolveGaussian
} from '../src/ApprParaboloid3.js';
import { Matrix } from '../src/Matrix.js';
import { Vector } from '../src/Vector.js';
import {
    check, expectClose, fc, finite, seededRandom
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function evaluate(u: readonly number[], x: number, y: number): number {
    return u[0] * x * x + u[1] * x * y + u[2] * y * y + u[3] * x + u[4] * y
        + u[5];
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

// Samples of the paraboloid on a grid in [-1,1]^2 scaled by 'extent'.
function gridSamples(u: readonly number[], extent: number, n: number,
    noise?: () => number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        const x = extent * (-1 + 2 * i / (n - 1));
        for (let j = 0; j < n; ++j) {
            const y = extent * (-1 + 2 * j / (n - 1));
            points.push(v3(x, y, evaluate(u, x, y) + (noise ? noise() : 0)));
        }
    }
    return points;
}

describe('ApprParaboloid3.fit', () => {
    it('recovers the coefficients of an exact paraboloid', () => {
        const u = [1.5, -0.5, 2, 3, -1, 0.25];
        const result = ApprParaboloid3.fit(gridSamples(u, 2, 7));
        expect(result.success).toBe(true);
        for (let i = 0; i < 6; ++i) {
            expect(result.u[i]).toBeCloseTo(u[i], 9);
        }
        expect(result.meanSquareError).toBeLessThan(1e-9);
    });

    it('recovers a plane (all quadratic coefficients zero)', () => {
        const u = [0, 0, 0, 2, -3, 7];
        const result = ApprParaboloid3.fit(gridSamples(u, 1, 5));
        expect(result.success).toBe(true);
        for (let i = 0; i < 6; ++i) {
            expect(result.u[i]).toBeCloseTo(u[i], 9);
        }
    });

    it('is the least-squares minimizer (normal equations hold)', () => {
        const rand = makeRandom(777);
        const points: Vector[] = [];
        for (let i = 0; i < 80; ++i) {
            const x = -2 + 4 * rand();
            const y = -2 + 4 * rand();
            points.push(v3(x, y,
                evaluate([1, 0.5, -1, 0.25, 2, -0.5], x, y)
                + 0.2 * (rand() - 0.5)));
        }
        const { success, u } = ApprParaboloid3.fit(points);
        expect(success).toBe(true);

        const gradient = [0, 0, 0, 0, 0, 0];
        for (const p of points) {
            const x = p.values[0], y = p.values[1], z = p.values[2];
            const e = evaluate(u, x, y) - z;
            const basis = [x * x, x * y, y * y, x, y, 1];
            for (let k = 0; k < 6; ++k) {
                gradient[k] += e * basis[k];
            }
        }
        for (let k = 0; k < 6; ++k) {
            expect(Math.abs(gradient[k]) / points.length).toBeLessThan(1e-10);
        }
    });

    it('recovers coefficients from noisy samples', () => {
        const rand = makeRandom(31415);
        const u = [0.75, 0.2, -1.25, 0.5, -0.5, 1];
        const points = gridSamples(u, 2, 15, () => 0.01 * (rand() - 0.5));
        const result = ApprParaboloid3.fit(points);
        expect(result.success).toBe(true);
        for (let i = 0; i < 6; ++i) {
            expect(result.u[i]).toBeCloseTo(u[i], 2);
        }
        expect(result.meanSquareError).toBeLessThan(0.01);
    });

    it('reports failure for a degenerate (collinear in xy) point set', () => {
        // All samples have y = x, so the six basis functions are linearly
        // dependent and the 6x6 system is singular.
        const points: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            const t = i;
            points.push(v3(t, t, t * t));
        }
        const result = ApprParaboloid3.fit(points);
        expect(result.success).toBe(false);
        expect(result.u).toEqual([0, 0, 0, 0, 0, 0]);
    });

    it('throws when there are too few points', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 5; ++i) {
            points.push(v3(i, i * i, 0));
        }
        expect(() => ApprParaboloid3.fit(points)).toThrow(
            /Insufficient points/);
    });

    it('recovers random paraboloids from exact samples', () => {
        const rand = makeRandom(20260901);
        for (let trial = 0; trial < 10; ++trial) {
            const u: number[] = [];
            for (let i = 0; i < 6; ++i) {
                u.push(-2 + 4 * rand());
            }
            const result = ApprParaboloid3.fit(gridSamples(u, 1.5, 6));
            expect(result.success).toBe(true);
            for (let i = 0; i < 6; ++i) {
                expect(result.u[i]).toBeCloseTo(u[i], 7);
            }
        }
    });
});

describe('ApprParaboloid3.fitRobust', () => {
    it('recovers the shifted coefficients and the average', () => {
        const u = [1.5, -0.5, 2, 3, -1, 0.25];
        const points = gridSamples(u, 2, 7);
        const result = ApprParaboloid3.fitRobust(points);
        expect(result.success).toBe(true);

        const avg = [0, 0, 0];
        for (const p of points) {
            for (let k = 0; k < 3; ++k) {
                avg[k] += p.values[k];
            }
        }
        for (let k = 0; k < 3; ++k) {
            avg[k] /= points.length;
            expect(result.average.values[k]).toBeCloseTo(avg[k], 10);
        }

        // Convert the v-polynomial back to the u-polynomial.
        const a = avg[0], b = avg[1], c = avg[2];
        const v = result.v;
        const back = [
            v[0],
            v[1],
            v[2],
            v[3] - v[0] * 2 * a - v[1] * b,
            v[4] - v[1] * a - v[2] * 2 * b,
            v[0] * a * a + v[1] * a * b + v[2] * b * b - v[3] * a - v[4] * b
                + v[5] + c
        ];
        for (let i = 0; i < 6; ++i) {
            expect(back[i]).toBeCloseTo(u[i], 8);
        }
        expect(result.meanSquareError).toBeLessThan(1e-9);
    });

    it('agrees with fit on well-conditioned data', () => {
        const rand = makeRandom(555);
        const u = [0.4, -0.3, 0.9, 1.2, -0.8, 2.5];
        const points: Vector[] = [];
        for (let i = 0; i < 100; ++i) {
            const x = -1.5 + 3 * rand();
            const y = -1.5 + 3 * rand();
            points.push(v3(x, y,
                evaluate(u, x, y) + 0.05 * (rand() - 0.5)));
        }
        const direct = ApprParaboloid3.fit(points);
        const robust = ApprParaboloid3.fitRobust(points);
        expect(direct.success && robust.success).toBe(true);

        const a = robust.average.values[0];
        const b = robust.average.values[1];
        const c = robust.average.values[2];
        const v = robust.v;
        const back = [
            v[0],
            v[1],
            v[2],
            v[3] - v[0] * 2 * a - v[1] * b,
            v[4] - v[1] * a - v[2] * 2 * b,
            v[0] * a * a + v[1] * a * b + v[2] * b * b - v[3] * a - v[4] * b
                + v[5] + c
        ];
        for (let i = 0; i < 6; ++i) {
            expect(back[i]).toBeCloseTo(direct.u[i], 6);
        }
    });

    it('throws when there are too few points', () => {
        expect(() => ApprParaboloid3.fitRobust([v3(0, 0, 0)])).toThrow(
            /Insufficient points/);
    });
});

describe('apprParaboloid3SolveGaussian', () => {
    it('solves the same system as the LDLT-based solver', () => {
        // A symmetric positive definite 6x6 system with a known solution.
        const rand = makeRandom(8080);
        const M = new Matrix(6, 6);
        for (let r = 0; r < 6; ++r) {
            for (let c = 0; c < 6; ++c) {
                M.set(r, c, -1 + 2 * rand());
            }
        }
        const A = new Matrix(6, 6);
        for (let r = 0; r < 6; ++r) {
            for (let c = 0; c < 6; ++c) {
                let sum = (r === c ? 6 : 0);
                for (let k = 0; k < 6; ++k) {
                    sum += M.get(k, r) * M.get(k, c);
                }
                A.set(r, c, sum);
            }
        }
        const expected = [1, -2, 3, -4, 5, -6];
        const B = new Vector(6);
        for (let r = 0; r < 6; ++r) {
            let sum = 0;
            for (let c = 0; c < 6; ++c) {
                sum += A.get(r, c) * expected[c];
            }
            B.values[r] = sum;
        }

        const { success, X } = apprParaboloid3SolveGaussian(A, B);
        expect(success).toBe(true);
        for (let i = 0; i < 6; ++i) {
            expect(X.values[i]).toBeCloseTo(expected[i], 8);
        }
    });
});

// ---------------------------------------------------------------------------
// Property-based verification (VERIFYING.md).

describe('ApprParaboloid3 verification', () => {
    // Random coefficients over a fixed 4x4 grid: the grid keeps the 6x6
    // moment matrix well conditioned, so the property tests the algebra of
    // the port rather than the conditioning of a random design matrix.
    const model = fc.tuple(
        fc.array(finite(-3, 3), { minLength: 6, maxLength: 6 }),
        fc.integer({ min: 3, max: 5 }), finite(0.5, 3));

    it('recovers the coefficients of an exact paraboloid', () => {
        check(model, ([u, n, extent]) => {
            const points = gridSamples(u, extent, n);
            const result = ApprParaboloid3.fit(points);
            expect(result.success).toBe(true);
            for (let k = 0; k < 6; ++k) {
                expectClose(result.u[k], u[k], 1e-7, 1e-7);
            }
            expect(result.meanSquareError).toBeLessThan(1e-7);
        }, 100);
    });

    it('reports sqrt(sum of squared residuals)/n, the upstream metric', () => {
        // Upstream computes sqrt(totalSqrError)/numPoints, which is neither
        // the mean square error nor the RMS error (upstream bug suspect).
        const rnd = seededRandom(112358);
        check(model, ([u, n, extent]) => {
            const points = gridSamples(u, extent, n,
                () => 0.1 * (2 * rnd() - 1));
            const result = ApprParaboloid3.fit(points);
            expect(result.success).toBe(true);
            let total = 0;
            for (const p of points) {
                const e = evaluate(result.u, p.values[0], p.values[1])
                    - p.values[2];
                total += e * e;
            }
            expectClose(result.meanSquareError,
                Math.sqrt(total) / points.length, 1e-14, 1e-9);
        }, 100);
    });

    it('agrees with fitRobust through the documented change of variables', () => {
        // u0 = v0, u1 = v1, u2 = v2, u3 = v3 - 2*v0*a - v1*b,
        // u4 = v4 - v1*a - 2*v2*b, u5 = v0*a^2 + v1*a*b + v2*b^2
        //      - v3*a - v4*b + v5 + c.
        const rnd = seededRandom(271828);
        check(model, ([u, n, extent]) => {
            const points = gridSamples(u, extent, n,
                () => 0.05 * (2 * rnd() - 1));
            const direct = ApprParaboloid3.fit(points);
            const robust = ApprParaboloid3.fitRobust(points);
            expect(direct.success).toBe(true);
            expect(robust.success).toBe(true);

            const a = robust.average.values[0];
            const b = robust.average.values[1];
            const c = robust.average.values[2];
            const v = robust.v;
            const converted = [
                v[0], v[1], v[2],
                v[3] - 2 * v[0] * a - v[1] * b,
                v[4] - v[1] * a - 2 * v[2] * b,
                v[0] * a * a + v[1] * a * b + v[2] * b * b
                    - v[3] * a - v[4] * b + v[5] + c
            ];
            // The centered system is better conditioned than the uncentered
            // one, so the two solutions agree only to the conditioning of the
            // uncentered normal equations.
            for (let k = 0; k < 6; ++k) {
                expectClose(direct.u[k], converted[k], 1e-6, 1e-6);
            }
        }, 100);
    });

    it('agrees with the Gaussian-elimination solve on the same system', () => {
        // Upstream selects LDLTDecomposition for floating point and
        // LinearSystem (Gaussian elimination) for the rational instantiation.
        // Both must solve the same 6x6 system, so the fitted coefficients
        // agree; this exercises apprParaboloid3SolveGaussian, the port of the
        // rational branch's solver.
        check(model, ([u, n, extent]) => {
            const points = gridSamples(u, extent, n);
            const A = new Matrix(6, 6);
            const B = new Vector(6);
            for (const p of points) {
                const x = p.values[0], y = p.values[1], z = p.values[2];
                const basis = [x * x, x * y, y * y, x, y, 1];
                for (let r = 0; r < 6; ++r) {
                    for (let cc = 0; cc < 6; ++cc) {
                        A.set(r, cc, A.get(r, cc) + basis[r] * basis[cc]);
                    }
                    B.values[r] += basis[r] * z;
                }
            }
            const gaussian = apprParaboloid3SolveGaussian(A, B);
            expect(gaussian.success).toBe(true);
            const fitted = ApprParaboloid3.fit(points);
            for (let k = 0; k < 6; ++k) {
                expectClose(gaussian.X.values[k], fitted.u[k], 1e-6, 1e-6);
            }
        }, 60);
    });

    it('behaves as documented on degenerate inputs', () => {
        const five = gridSamples([1, 0, 1, 0, 0, 0], 1, 3).slice(0, 5);
        expect(() => ApprParaboloid3.fit(five))
            .toThrow(/Insufficient points/);
        expect(() => ApprParaboloid3.fitRobust(five))
            .toThrow(/Insufficient points/);

        // Samples on a line in the xy-plane leave the 6x6 system singular.
        const collinear: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            collinear.push(v3(i, i, i * i));
        }
        expect(ApprParaboloid3.fit(collinear).success).toBe(false);
        expect(ApprParaboloid3.fit(collinear).u).toEqual([0, 0, 0, 0, 0, 0]);
        expect(ApprParaboloid3.fit(collinear).meanSquareError).toBe(0);
    });
});
