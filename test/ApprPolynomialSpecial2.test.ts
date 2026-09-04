import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial2 } from '../src/ApprPolynomialSpecial2.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

// Solve the 2x2 system A*X = B by Cramer's rule; an independent check.
function solve2(A: number[][], B: number[]): number[] {
    const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    return [
        (B[0] * A[1][1] - A[0][1] * B[1]) / det,
        (A[0][0] * B[1] - B[0] * A[1][0]) / det
    ];
}

describe('ApprPolynomialSpecial2', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial2([]))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial2([0, 2, 2]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([2, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([-1, 3]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial2([0, 1, 5])).not.toThrow();
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial2([0, 2, 5]);
        expect(fitter.getParameters()).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(3);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 3 + 5*x^2 with the term set {1, x^2}', () => {
        // The internal transform maps x in [-1,1] to itself and w affinely,
        // so the model {1, x^2} represents the data exactly.
        const observations = [-1, -0.5, 0, 0.5, 1].map(
            (x) => [x, 3 + 5 * x * x]);

        const fitter = new ApprPolynomialSpecial2([0, 2]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);

        // In the transformed space, w' = -1 + 2*(x')^2.
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(-1, 10);
        expect(p[1]).toBeCloseTo(2, 10);

        // evaluate() undoes the transform, so it reproduces the original
        // data and interpolates/extrapolates the underlying parabola.
        for (const [x, w] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(w, 10);
            expect(fitter.error([x, w])).toBeLessThan(1e-10);
        }
        for (const x of [0.25, -0.8, 2, -3]) {
            expect(fitter.evaluate(x)).toBeCloseTo(3 + 5 * x * x, 8);
        }
    });

    it('restricts the fit to the requested powers', () => {
        // With the term set {1, x^3} an odd cubic through the origin is
        // recovered but a quadratic component cannot be.
        const observations = [-1, -0.5, 0, 0.5, 1].map(
            (x) => [x, 4 * x * x * x]);
        const fitter = new ApprPolynomialSpecial2([0, 3]);
        expect(fitter.fit(observations)).toBe(true);
        for (const [x, w] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(w, 10);
        }
        expect(fitter.evaluate(0.75)).toBeCloseTo(4 * 0.75 ** 3, 8);
    });

    it('matches the transformed normal equations for the term set {1, x}', () => {
        const observations = [[0, 1], [1, 1.9], [2, 3.2], [3, 3.9], [4, 5.4]];

        // Reproduce the upstream transform to [-1,1]^2.
        const xs = observations.map((o) => o[0]);
        const ws = observations.map((o) => o[1]);
        const xmin = Math.min(...xs), xmax = Math.max(...xs);
        const wmin = Math.min(...ws), wmax = Math.max(...ws);
        const sx = 1 / (xmax - xmin), sw = 1 / (wmax - wmin);
        const tx = xs.map((x) => -1 + 2 * sx * (x - xmin));
        const tw = ws.map((w) => -1 + 2 * sw * (w - wmin));

        // A and B are the normalized sums over the terms (x^0, x^1).
        const n = observations.length;
        const A = [[0, 0], [0, 0]];
        const B = [0, 0];
        for (let i = 0; i < n; ++i) {
            for (let r = 0; r < 2; ++r) {
                for (let c = 0; c < 2; ++c) {
                    A[r][c] += Math.pow(tx[i], r + c);
                }
                B[r] += Math.pow(tx[i], r) * tw[i];
            }
        }
        for (let r = 0; r < 2; ++r) {
            B[r] /= n;
            for (let c = 0; c < 2; ++c) {
                A[r][c] /= n;
            }
        }
        const expected = solve2(A, B);

        const fitter = new ApprPolynomialSpecial2([0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(expected[0], 10);
        expect(p[1]).toBeCloseTo(expected[1], 10);

        // Because the transform is affine in both x and w, the {1, x} fit
        // agrees with the ordinary least-squares line in the original space.
        let s1 = 0, sxs = 0, sx2 = 0, sws = 0, sxw = 0;
        for (const [x, w] of observations) {
            s1 += 1; sxs += x; sx2 += x * x; sws += w; sxw += x * w;
        }
        const det = s1 * sx2 - sxs * sxs;
        const a0 = (sws * sx2 - sxs * sxw) / det;
        const a1 = (s1 * sxw - sws * sxs) / det;
        for (const [x] of observations) {
            expect(fitter.evaluate(x)).toBeCloseTo(a0 + a1 * x, 8);
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial2([0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0], [1, 1]], [0, 9])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const source = new ApprPolynomialSpecial2([0, 2]);
        source.fit([-1, -0.5, 0, 0.5, 1].map((x) => [x, 3 + 5 * x * x]));
        const target = new ApprPolynomialSpecial2([0, 2]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.evaluate(0.3)).toBe(source.evaluate(0.3));

        // The copy is deep, including the transform state.
        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3);
        source.fit([-2, -1, 0, 1, 2].map((x) => [x, 1 - x * x]));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3)).toBe(value);
    });
});

describe('ApprPolynomialSpecial2 verification', () => {
    // Strictly increasing nonnegative powers, as the constructor demands.
    const degreesArb = fc.array(fc.integer({ min: 0, max: 2 }),
        { minLength: 1, maxLength: 3 })
        .map(gaps => {
            let acc = -1;
            return gaps.map(g => (acc += g + 1));
        });

    // Well-separated abscissas so the transformed normal equations are
    // conditioned; the transform needs xmax > xmin to be finite at all.
    const abscissas = (count: number) =>
        fc.array(finite(0.4, 1.4), { minLength: count, maxLength: count })
            .map(gaps => {
                let acc = -2;
                return gaps.map(g => (acc += g));
            });

    const configArb = degreesArb.chain(degrees => fc.tuple(
        fc.constant(degrees),
        abscissas(degrees.length + 3),
        fc.array(finite(-4, 4),
            { minLength: degrees.length + 3, maxLength: degrees.length + 3 })
            .filter(ws => Math.max(...ws) - Math.min(...ws) > 0.5)));

    it('is a least-squares fit in the transformed coordinates', () => {
        // Upstream maps the samples to [-1,1]^2 before fitting, so the
        // residual of the transformed model is orthogonal to every basis
        // function (x')^{p_i}.
        check(configArb, ([degrees, xs, ws]) => {
            const observations = xs.map((x, i) => [x, ws[i]]);
            const fitter = new ApprPolynomialSpecial2(degrees);
            if (!fitter.fit(observations)) { return; }
            const c = fitter.getParameters();

            const xmin = Math.min(...xs), xmax = Math.max(...xs);
            const wmin = Math.min(...ws), wmax = Math.max(...ws);
            const tx = (x: number): number =>
                -1 + (2 * (x - xmin)) / (xmax - xmin);
            const tw = (w: number): number =>
                -1 + (2 * (w - wmin)) / (wmax - wmin);
            const model = (xp: number): number => degrees.reduce(
                (u, p, i) => u + c[i] * Math.pow(xp, p), 0);

            for (let i = 0; i < degrees.length; ++i) {
                let residual = 0, scale = 0;
                for (const [x, w] of observations) {
                    const xp = tx(x);
                    const basis = Math.pow(xp, degrees[i]);
                    residual += basis * (model(xp) - tw(w));
                    scale += Math.abs(basis) * 2;
                }
                expect(Math.abs(residual))
                    .toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
            }
        });
    });

    it('evaluate applies the stored transform in both directions', () => {
        // Evaluate maps x into [-1,1], sums c[i]*(x')^{p_i} and maps the
        // result back with (w'+1)*mInvTwoWScale + wmin.
        check(fc.tuple(configArb, finite(-6, 6)),
            ([[degrees, xs, ws], x]) => {
                const observations = xs.map((u, i) => [u, ws[i]]);
                const fitter = new ApprPolynomialSpecial2(degrees);
                if (!fitter.fit(observations)) { return; }
                const c = fitter.getParameters();

                const xmin = Math.min(...xs), xmax = Math.max(...xs);
                const wmin = Math.min(...ws), wmax = Math.max(...ws);
                const xp = -1 + (2 * (x - xmin)) / (xmax - xmin);
                const wp = degrees.reduce(
                    (u, p, i) => u + c[i] * Math.pow(xp, p), 0);
                const expected = ((wp + 1) * (wmax - wmin)) / 2 + wmin;

                const scale = 1 + Math.abs(expected);
                expectClose(fitter.evaluate(x), expected, 1e-7 * scale, 1e-7);
            });
    });

    it('reproduces samples of a full monomial basis exactly', () => {
        // The affine transform of the samples preserves the span only when
        // the term set is the full basis {1, x, ..., x^d}; then the exact fit
        // survives the change of coordinates.
        check(fc.tuple(fc.integer({ min: 0, max: 3 }), fc.array(finite(-3, 3),
            { minLength: 4, maxLength: 4 })).chain(([d, c]) =>
                fc.tuple(fc.constant(d), fc.constant(c.slice(0, d + 1)),
                    abscissas(d + 3))),
            ([d, c, xs]) => {
                const degrees = Array.from({ length: d + 1 }, (_, i) => i);
                const poly = (x: number): number =>
                    c.reduce((u, v, i) => u + v * Math.pow(x, i), 0);
                const observations = xs.map(x => [x, poly(x)]);
                // Transform divides by (wmax - wmin), so constant height
                // data has no fit at all (pinned separately below).
                const hs = observations.map(o => o[1]);
                if (Math.max(...hs) - Math.min(...hs) < 0.5) { return; }
                const fitter = new ApprPolynomialSpecial2(degrees);
                if (!fitter.fit(observations)) { return; }

                const scale = c.reduce((u, v) => u + Math.abs(v), 0) + 1;
                for (const [x, w] of observations) {
                    expectClose(fitter.evaluate(x), w, 1e-6 * scale, 1e-6);
                    expectClose(fitter.error([x, w]), 0, 1e-6 * scale, 0);
                }
            });
    });

    it('reports the term count and the domain of the most recent fit', () => {
        // Unlike ApprPolynomial2, Transform assigns (rather than accumulates)
        // the domain, so it tracks only the latest sample set.
        check(fc.tuple(configArb, configArb),
            ([[degrees, xs0, ws0], [, xs1, ws1]]) => {
                const fitter = new ApprPolynomialSpecial2(degrees);
                expect(fitter.getMinimumRequired()).toBe(degrees.length);
                expect(fitter.getXDomain()[0]).toBe(Number.MAX_VALUE);

                fitter.fit(xs0.map((x, i) => [x, ws0[i]]));
                expect(fitter.getXDomain())
                    .toEqual([Math.min(...xs0), Math.max(...xs0)]);

                const shifted = xs1.map(x => x + 50);
                fitter.fit(shifted.map((x, i) => [x, ws1[i % ws1.length]]));
                expect(fitter.getXDomain())
                    .toEqual([Math.min(...shifted), Math.max(...shifted)]);
            });
    });

    it('rejects degree lists that are not strictly increasing', () => {
        check(fc.array(fc.integer({ min: -2, max: 4 }),
            { minLength: 0, maxLength: 4 }), degrees => {
                const strictlyIncreasing = degrees.length > 0
                    && degrees[0] >= 0
                    && degrees.every((d, i) => i === 0 || d > degrees[i - 1]);
                if (strictlyIncreasing) {
                    expect(() => new ApprPolynomialSpecial2(degrees))
                        .not.toThrow();
                }
                else if (degrees.length === 0) {
                    expect(() => new ApprPolynomialSpecial2(degrees))
                        .toThrowError('The input array must have elements.');
                }
                else {
                    expect(() => new ApprPolynomialSpecial2(degrees))
                        .toThrowError('Degrees must be increasing.');
                }
            });
    });

    it('copyParameters deep-copies the model', () => {
        check(configArb, ([degrees, xs, ws]) => {
            const source = new ApprPolynomialSpecial2(degrees);
            if (!source.fit(xs.map((x, i) => [x, ws[i]]))) { return; }
            const target = new ApprPolynomialSpecial2(degrees);
            target.copyParameters(source);

            expect([...target.getParameters()])
                .toEqual([...source.getParameters()]);
            expect(target.getParameters()).not.toBe(source.getParameters());
            expect(target.getXDomain()).toEqual(source.getXDomain());
            expect(target.evaluate(0.5)).toBe(source.evaluate(0.5));

            const copied = [...target.getParameters()];
            source.fit(xs.map((x, i) => [x, ws[i] + 10]));
            expect([...target.getParameters()]).toEqual(copied);
        });
    });

    it('produces NaN for constant height data', () => {
        // Transform computes mScale[i] = 1 / (omax[i] - omin[i]) with no
        // guard, so samples whose heights are all equal divide by zero and
        // every transformed value is NaN. Preserved from upstream.
        check(fc.tuple(fc.array(finite(0.4, 1.4),
            { minLength: 4, maxLength: 4 }), finite(-5, 5)),
            ([gaps, w]) => {
                let acc = -2;
                const xs = gaps.map(g => (acc += g));
                const fitter = new ApprPolynomialSpecial2([0, 1]);
                fitter.fit(xs.map(x => [x, w]));
                expect(Number.isNaN(fitter.evaluate(xs[0]))).toBe(true);
            });
    });
});
