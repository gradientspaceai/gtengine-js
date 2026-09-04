import { describe, expect, it } from 'vitest';
import { ApprPolynomial2 } from '../src/ApprPolynomial2.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

// Evaluate sum_i c[i]*x^i.
function poly(c: readonly number[], x: number): number {
    let w = 0;
    for (let i = c.length - 1; i >= 0; --i) {
        w = c[i] + w * x;
    }
    return w;
}

// Solve the (small, well-conditioned) system A*X = B by Gaussian
// elimination with partial pivoting; an independent check of the fit.
function solve(A: number[][], B: number[]): number[] {
    const n = B.length;
    const M = A.map((row, i) => row.concat([B[i]]));
    for (let c = 0; c < n; ++c) {
        let piv = c;
        for (let r = c + 1; r < n; ++r) {
            if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) {
                piv = r;
            }
        }
        const t = M[c]; M[c] = M[piv]; M[piv] = t;
        for (let r = 0; r < n; ++r) {
            if (r !== c) {
                const f = M[r][c] / M[c][c];
                for (let k = c; k <= n; ++k) {
                    M[r][k] -= f * M[c][k];
                }
            }
        }
    }
    return M.map((row, i) => row[n] / row[i]);
}

describe('ApprPolynomial2', () => {
    it('initializes the parameters to zero and reports the domain sentinel', () => {
        const fitter = new ApprPolynomial2(3);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(4);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('recovers a cubic exactly from data sampled on it', () => {
        const c = [-2, 0.5, 3, -1.25];
        const observations: number[][] = [];
        for (let i = 0; i <= 10; ++i) {
            const x = -1 + 0.2 * i;
            observations.push([x, poly(c, x)]);
        }

        const fitter = new ApprPolynomial2(3);
        expect(fitter.fit(observations)).toBe(true);

        const p = fitter.getParameters();
        expect(p.length).toBe(4);
        for (let i = 0; i < 4; ++i) {
            expect(p[i]).toBeCloseTo(c[i], 10);
        }

        // The model error at every sample is (essentially) zero.
        for (const obs of observations) {
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }

        // Interpolation and extrapolation both agree with the cubic.
        for (const x of [-0.9, 0.13, 0.77, 3.5, -4]) {
            expect(fitter.evaluate(x)).toBeCloseTo(poly(c, x), 8);
        }
    });

    it('records the x-domain of the fitted samples', () => {
        const fitter = new ApprPolynomial2(1);
        fitter.fit([[-3, 1], [2, 4], [7.5, -1]]);
        expect(fitter.getXDomain()).toEqual([-3, 7.5]);
    });

    it('accumulates the x-domain across fits, as upstream does', () => {
        // Upstream never resets mXDomain between calls to FitIndexed; the
        // port preserves that behavior.
        const fitter = new ApprPolynomial2(1);
        fitter.fit([[0, 0], [1, 1]]);
        expect(fitter.getXDomain()).toEqual([0, 1]);
        fitter.fit([[5, 0], [6, 1]]);
        expect(fitter.getXDomain()).toEqual([0, 6]);
    });

    it('matches the normal equations for an overdetermined line fit', () => {
        const observations = [[0, 1], [1, 1.9], [2, 3.2], [3, 3.9], [4, 5.4]];

        // The least-squares system for w = c0 + c1*x is
        //   [ S1 Sx ] [c0]   [ Sw   ]
        //   [ Sx Sx2] [c1] = [ Sxw  ]
        let s1 = 0, sx = 0, sx2 = 0, sw = 0, sxw = 0;
        for (const [x, w] of observations) {
            s1 += 1; sx += x; sx2 += x * x; sw += w; sxw += x * w;
        }
        const expected = solve([[s1, sx], [sx, sx2]], [sw, sxw]);

        const fitter = new ApprPolynomial2(1);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(expected[0], 10);
        expect(p[1]).toBeCloseTo(expected[1], 10);

        // The residuals are orthogonal to the model basis (1, x), the
        // defining property of a least-squares fit.
        let r0 = 0, r1 = 0;
        for (const [x, w] of observations) {
            const r = fitter.evaluate(x) - w;
            r0 += r;
            r1 += r * x;
        }
        expect(r0).toBeCloseTo(0, 10);
        expect(r1).toBeCloseTo(0, 10);
    });

    it('fits a degree-0 polynomial as the sample mean', () => {
        const observations = [[0, 2], [1, 4], [2, 9]];
        const fitter = new ApprPolynomial2(0);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getParameters()[0]).toBeCloseTo(5, 12);
        expect(fitter.evaluate(1000)).toBeCloseTo(5, 12);
    });

    it('fails when the Vandermonde system is singular', () => {
        // Two samples cannot determine a cubic; the normal-equations matrix
        // is singular, GaussianElimination returns the zero inverse, and all
        // coefficients are zero, so the fit reports failure.
        const fitter = new ApprPolynomial2(3);
        expect(fitter.fit([[0, 1], [1, 2]])).toBe(false);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0]);
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomial2(1);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0], [1, 1]], [0, 5])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const source = new ApprPolynomial2(2);
        source.fit([[0, 1], [1, 3], [2, 9], [3, 19]]);
        const target = new ApprPolynomial2(2);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.evaluate(1.5)).toBe(source.evaluate(1.5));

        // The copy is deep: refitting the source leaves the target alone.
        const copied = target.getParameters().slice();
        source.fit([[0, 0], [1, 0], [2, 0], [3, 1]]);
        expect(target.getParameters()).toEqual(copied);
    });

    it('fits contiguous and indexed subsets of the observations', () => {
        // The first two observations are outliers; fit only the last three,
        // which lie exactly on w = 1 + 2*x.
        const observations = [[0, 100], [1, -100], [2, 5], [3, 7], [4, 9]];

        const contiguous = new ApprPolynomial2(1);
        expect(contiguous.fit(observations, 2, 4)).toBe(true);
        expect(contiguous.getParameters()[0]).toBeCloseTo(1, 10);
        expect(contiguous.getParameters()[1]).toBeCloseTo(2, 10);

        const indexed = new ApprPolynomial2(1);
        expect(indexed.fit(observations, [2, 3, 4])).toBe(true);
        expect(indexed.getParameters()[0]).toBeCloseTo(1, 10);
        expect(indexed.getParameters()[1]).toBeCloseTo(2, 10);
    });

    it('cross-checks randomized data against the normal equations', () => {
        let seed = 12345;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 20; ++trial) {
            const observations: number[][] = [];
            for (let i = 0; i < 12; ++i) {
                observations.push([-1 + 2 * rand(), -1 + 2 * rand()]);
            }

            const degree = 2;
            const size = degree + 1;
            const A: number[][] = [];
            const B: number[] = [];
            for (let r = 0; r < size; ++r) {
                A.push(new Array<number>(size).fill(0));
                B.push(0);
                for (const [x, w] of observations) {
                    B[r] += w * Math.pow(x, r);
                    for (let c = 0; c < size; ++c) {
                        A[r][c] += Math.pow(x, r + c);
                    }
                }
            }
            const expected = solve(A, B);

            const fitter = new ApprPolynomial2(degree);
            expect(fitter.fit(observations)).toBe(true);
            const p = fitter.getParameters();
            for (let i = 0; i < size; ++i) {
                expect(p[i]).toBeCloseTo(expected[i], 8);
            }
        }
    });
});

describe('ApprPolynomial2 verification', () => {
    // Strictly increasing abscissas with gaps bounded away from zero, so the
    // Vandermonde normal equations stay well conditioned. Upstream's own
    // header warns that the monomial basis is nonrobust for large degrees
    // and large data, so the degrees stay small and |x| <= 3.
    const abscissas = (count: number) =>
        fc.array(finite(0.4, 1.4), { minLength: count, maxLength: count })
            .map(gaps => {
                let acc = -3;
                return gaps.map(g => (acc += g));
            });

    const configArb = fc.integer({ min: 0, max: 3 })
        .chain(degree => fc.tuple(fc.constant(degree),
            // At least one coefficient must be nonzero: upstream returns
            // 'hasNonzero', so an exact fit of the zero polynomial is
            // reported as a failure (pinned separately below).
            fc.array(finite(-4, 4),
                { minLength: degree + 1, maxLength: degree + 1 })
                .filter(c => c.some(v => Math.abs(v) > 0.1)),
            fc.integer({ min: degree + 1, max: degree + 4 })
                .chain(n => abscissas(n))));

    const evaluatePoly = (c: readonly number[], x: number): number => {
        let w = 0;
        for (let i = c.length - 1; i >= 0; --i) { w = c[i] + w * x; }
        return w;
    };

    it('recovers a polynomial its samples lie on', () => {
        check(configArb, ([degree, coefficients, xs]) => {
            const observations = xs.map(
                x => [x, evaluatePoly(coefficients, x)]);
            const fitter = new ApprPolynomial2(degree);
            expect(fitter.fit(observations)).toBe(true);

            const fitted = fitter.getParameters();
            expect(fitted.length).toBe(degree + 1);
            const scale = coefficients.reduce(
                (u, c) => u + Math.abs(c), 0) + 1;
            for (let i = 0; i <= degree; ++i) {
                expectClose(fitted[i], coefficients[i], 1e-8 * scale, 1e-8);
            }

            // evaluate() reproduces the polynomial, and the model error
            // vanishes on the samples.
            for (const [x, w] of observations) {
                expectClose(fitter.evaluate(x), w, 1e-8 * scale, 1e-8);
                expectClose(fitter.error([x, w]), 0, 1e-7 * scale, 0);
            }
        });
    });

    it('satisfies the least-squares normal equations', () => {
        check(fc.tuple(configArb, fc.array(finite(-5, 5),
            { minLength: 8, maxLength: 8 })), ([[degree, , xs], ws]) => {
                const observations = xs.map(
                    (x, i) => [x, ws[i % ws.length]]);
                const fitter = new ApprPolynomial2(degree);
                // An all-zero solution is reported as a failure by upstream.
                if (!fitter.fit(observations)) { return; }

                // sum_s x_s^i * (p(x_s) - w_s) = 0 for i = 0..degree.
                for (let i = 0; i <= degree; ++i) {
                    let residual = 0, scale = 0;
                    for (const [x, w] of observations) {
                        const xp = Math.pow(x, i);
                        residual += xp * (fitter.evaluate(x) - w);
                        scale += Math.abs(xp) * (Math.abs(w) + 1);
                    }
                    expect(Math.abs(residual))
                        .toBeLessThanOrEqual(1e-7 + 1e-7 * scale);
                }
            });
    });

    it('evaluate is Horner over the reported parameters', () => {
        check(fc.tuple(configArb, finite(-6, 6)),
            ([[degree, coefficients, xs], x]) => {
                const fitter = new ApprPolynomial2(degree);
                fitter.fit(xs.map(u => [u, evaluatePoly(coefficients, u)]));
                const c = fitter.getParameters();
                let expected = c[degree];
                for (let i = degree - 1; i >= 0; --i) {
                    expected = c[i] + expected * x;
                }
                expect(fitter.evaluate(x)).toBe(expected);
            });
    });

    it('reports the minimum required and accumulates the x-domain', () => {
        check(fc.tuple(configArb, configArb),
            ([[degree, c0, xs0], [, c1, xs1]]) => {
                const fitter = new ApprPolynomial2(degree);
                expect(fitter.getMinimumRequired()).toBe(degree + 1);

                // The initial domain is the empty-interval sentinel.
                expect(fitter.getXDomain()[0]).toBe(Number.MAX_VALUE);
                expect(fitter.getXDomain()[1]).toBe(-Number.MAX_VALUE);

                fitter.fit(xs0.map(x => [x, evaluatePoly(c0, x)]));
                expect(fitter.getXDomain()[0]).toBe(Math.min(...xs0));
                expect(fitter.getXDomain()[1]).toBe(Math.max(...xs0));

                // Upstream never resets the domain between fits, so it grows
                // to the union of every fitted sample set.
                fitter.fit(xs1.map(x => [x + 20, evaluatePoly(c1, x)]));
                expect(fitter.getXDomain()[0]).toBe(Math.min(...xs0));
                expect(fitter.getXDomain()[1])
                    .toBe(Math.max(...xs1.map(x => x + 20)));
            });
    });

    it('zeroes the parameters when the indices are rejected', () => {
        ApprQuery.validateIndices = true;
        try {
            check(configArb, ([degree, coefficients, xs]) => {
                const fitter = new ApprPolynomial2(degree);
                expect(fitter.fit(xs.map(
                    x => [x, evaluatePoly(coefficients, x)]))).toBe(true);
                expect(fitter.fit([[0, 0]], [0, 99])).toBe(false);
                expect([...fitter.getParameters()])
                    .toEqual(new Array<number>(degree + 1).fill(0));
            }, 60);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copyParameters deep-copies the model', () => {
        check(configArb, ([degree, coefficients, xs]) => {
            const source = new ApprPolynomial2(degree);
            expect(source.fit(xs.map(
                x => [x, evaluatePoly(coefficients, x)]))).toBe(true);
            const target = new ApprPolynomial2(degree);
            target.copyParameters(source);

            expect([...target.getParameters()])
                .toEqual([...source.getParameters()]);
            expect(target.getXDomain()).toEqual(source.getXDomain());
            expect(target.getXDomain()).not.toBe(source.getXDomain());
            expect(target.getParameters()).not.toBe(source.getParameters());

            // Refitting the source must not disturb the copy.
            const copied = [...target.getParameters()];
            source.fit([[0, 1], [1, 2], [2, 4], [3, 8], [4, 16]]);
            expect([...target.getParameters()]).toEqual(copied);
        });
    });

    it('reports failure for an exact fit of the zero polynomial', () => {
        // Upstream returns 'hasNonzero': every coefficient being zero is
        // reported as a failure even though the fit is exact and the stored
        // parameters are the correct ones.
        check(fc.integer({ min: 0, max: 3 }).chain(degree =>
            fc.tuple(fc.constant(degree), abscissas(degree + 2))),
            ([degree, xs]) => {
                const fitter = new ApprPolynomial2(degree);
                expect(fitter.fit(xs.map(x => [x, 0]))).toBe(false);
                expect([...fitter.getParameters()])
                    .toEqual(new Array<number>(degree + 1).fill(0));
                expect(fitter.evaluate(1.25)).toBe(0);
            });
    });
});
