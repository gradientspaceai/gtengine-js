import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial3 } from '../src/ApprPolynomialSpecial3.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

const grid = [-1, -0.5, 0, 0.5, 1];

describe('ApprPolynomialSpecial3', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial3([0, 1], [0]))
            .toThrowError('The input arrays must have the same size.');
        expect(() => new ApprPolynomialSpecial3([], []))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial3([0, 0], [0, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial3([0, 1], [1, 0]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial3([0, 1, 2], [0, 2, 3]))
            .not.toThrow();
    });

    it('rejects the natural term set {(0,0),(1,0),(0,1)}, as upstream does', () => {
        // Upstream requires each degree list to be strictly increasing,
        // which is stronger than "distinct <p,q> pairs"; the port preserves
        // the restriction.
        expect(() => new ApprPolynomialSpecial3([0, 1, 0], [0, 0, 1]))
            .toThrowError('Degrees must be increasing.');
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 2]);
        expect(fitter.getParameters()).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 5 + 2*x*y with the term set {1, x*y}', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);
        expect(fitter.getYDomain()).toEqual([-1, 1]);

        // The transform is the identity on x and y here, and w' = x'*y', so
        // the transformed coefficients are (0, 1).
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(0, 10);
        expect(p[1]).toBeCloseTo(1, 10);

        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1])).toBeCloseTo(obs[2], 10);
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }
        // Interpolation and extrapolation follow the same surface.
        for (const [x, y] of [[0.3, -0.4], [0.75, 0.75], [2, 3]]) {
            expect(fitter.evaluate(x, y)).toBeCloseTo(5 + 2 * x * y, 8);
        }
    });

    it('reproduces data from a three-term set {1, x*y, x^2*y^2}', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                const t = x * y;
                observations.push([x, y, 1 + 3 * t + 4 * t * t]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1, 2], [0, 1, 2]);
        expect(fitter.fit(observations)).toBe(true);
        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1])).toBeCloseTo(obs[2], 8);
        }
        for (const [x, y] of [[0.3, -0.4], [0.9, 0.6]]) {
            const t = x * y;
            expect(fitter.evaluate(x, y)).toBeCloseTo(1 + 3 * t + 4 * t * t, 6);
        }
    });

    it('produces a genuine least-squares fit for data off the model', () => {
        // w = 5 + 2*x*y plus a symmetric perturbation at each sample, so the
        // fit returns the unperturbed surface.
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y + 0.4]);
                observations.push([x, y, 5 + 2 * x * y - 0.4]);
            }
        }

        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);

        // The residuals are orthogonal to both basis terms in the
        // transformed space; equivalently, they cancel in pairs here.
        for (const x of grid) {
            for (const y of grid) {
                expect(fitter.evaluate(x, y)).toBeCloseTo(5 + 2 * x * y, 8);
            }
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 0], [1, 1, 1]], [0, 4])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const observations: number[][] = [];
        for (const x of grid) {
            for (const y of grid) {
                observations.push([x, y, 5 + 2 * x * y]);
            }
        }

        const source = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        source.fit(observations);
        const target = new ApprPolynomialSpecial3([0, 1], [0, 1]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.getYDomain()).toEqual(source.getYDomain());
        expect(target.evaluate(0.3, 0.4)).toBe(source.evaluate(0.3, 0.4));

        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3, 0.4);
        source.fit(observations.map(([x, y]) => [x, y, 1 - x * y]));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3, 0.4)).toBe(value);
    });
});

describe('ApprPolynomialSpecial3 verification', () => {
    // Upstream asserts that each degree list is separately strictly
    // increasing, so only "diagonal" term sets are admissible.
    const degreesArb = fc.array(fc.integer({ min: 0, max: 2 }),
        { minLength: 1, maxLength: 3 })
        .map(gaps => {
            let acc = -1;
            return gaps.map(g => (acc += g + 1));
        });

    const spread = (count: number, start: number) =>
        fc.array(finite(0.4, 1.4), { minLength: count, maxLength: count })
            .map(gaps => {
                let acc = start;
                return gaps.map(g => (acc += g));
            });

    const configArb = fc.tuple(degreesArb, degreesArb)
        .chain(([xDegrees, yDegrees]) => {
            const n = Math.min(xDegrees.length, yDegrees.length);
            const count = n + 3;
            return fc.tuple(fc.constant(xDegrees.slice(0, n)),
                fc.constant(yDegrees.slice(0, n)),
                spread(count, -2), spread(count, -1.5),
                fc.array(finite(-4, 4), { minLength: count, maxLength: count })
                    .filter(ws => Math.max(...ws) - Math.min(...ws) > 0.5));
        });

    it('is a least-squares fit in the transformed coordinates', () => {
        // The samples are mapped to [-1,1]^3 before fitting, so the residual
        // of the transformed model is orthogonal to every basis function
        // (x')^{p_i} * (y')^{q_i}.
        check(configArb, ([xd, yd, xs, ys, ws]) => {
            const observations = xs.map((x, i) => [x, ys[i], ws[i]]);
            const fitter = new ApprPolynomialSpecial3(xd, yd);
            if (!fitter.fit(observations)) { return; }
            const c = fitter.getParameters();

            const map = (v: number, lo: number, hi: number): number =>
                -1 + (2 * (v - lo)) / (hi - lo);
            const xlo = Math.min(...xs), xhi = Math.max(...xs);
            const ylo = Math.min(...ys), yhi = Math.max(...ys);
            const wlo = Math.min(...ws), whi = Math.max(...ws);
            const model = (xp: number, yp: number): number => xd.reduce(
                (u, p, i) => u + c[i] * Math.pow(xp, p) * Math.pow(yp, yd[i]),
                0);

            for (let i = 0; i < xd.length; ++i) {
                let residual = 0, scale = 0;
                for (const [x, y, w] of observations) {
                    const xp = map(x, xlo, xhi), yp = map(y, ylo, yhi);
                    const basis = Math.pow(xp, xd[i]) * Math.pow(yp, yd[i]);
                    residual += basis * (model(xp, yp) - map(w, wlo, whi));
                    scale += Math.abs(basis) * 2;
                }
                expect(Math.abs(residual))
                    .toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
            }
        });
    });

    it('evaluate applies the stored transform in both directions', () => {
        check(fc.tuple(configArb, finite(-5, 5), finite(-5, 5)),
            ([[xd, yd, xs, ys, ws], x, y]) => {
                const observations = xs.map((u, i) => [u, ys[i], ws[i]]);
                const fitter = new ApprPolynomialSpecial3(xd, yd);
                if (!fitter.fit(observations)) { return; }
                const c = fitter.getParameters();

                const map = (v: number, lo: number, hi: number): number =>
                    -1 + (2 * (v - lo)) / (hi - lo);
                const wlo = Math.min(...ws), whi = Math.max(...ws);
                const xp = map(x, Math.min(...xs), Math.max(...xs));
                const yp = map(y, Math.min(...ys), Math.max(...ys));
                const wp = xd.reduce((u, p, i) =>
                    u + c[i] * Math.pow(xp, p) * Math.pow(yp, yd[i]), 0);
                const expected = ((wp + 1) * (whi - wlo)) / 2 + wlo;

                expectClose(fitter.evaluate(x, y), expected,
                    1e-7 * (1 + Math.abs(expected)), 1e-7);
            });
    });

    it('the model error is |evaluate(x,y) - w|', () => {
        check(fc.tuple(configArb, finite(-5, 5), finite(-5, 5), finite(-5, 5)),
            ([[xd, yd, xs, ys, ws], x, y, w]) => {
                const fitter = new ApprPolynomialSpecial3(xd, yd);
                if (!fitter.fit(xs.map((u, i) => [u, ys[i], ws[i]]))) {
                    return;
                }
                expect(fitter.error([x, y, w]))
                    .toBe(Math.abs(fitter.evaluate(x, y) - w));
            });
    });

    it('reports the term count and the domains of the most recent fit',
        () => {
            // Transform assigns rather than accumulates the domains, so they
            // track only the latest sample set.
            check(configArb, ([xd, yd, xs, ys, ws]) => {
                const fitter = new ApprPolynomialSpecial3(xd, yd);
                expect(fitter.getMinimumRequired()).toBe(xd.length);
                expect(fitter.getXDomain()[0]).toBe(Number.MAX_VALUE);
                expect(fitter.getYDomain()[1]).toBe(-Number.MAX_VALUE);

                fitter.fit(xs.map((x, i) => [x, ys[i], ws[i]]));
                expect(fitter.getXDomain())
                    .toEqual([Math.min(...xs), Math.max(...xs)]);
                expect(fitter.getYDomain())
                    .toEqual([Math.min(...ys), Math.max(...ys)]);

                const moved = xs.map(x => x + 60);
                fitter.fit(moved.map((x, i) => [x, ys[i], ws[i]]));
                expect(fitter.getXDomain())
                    .toEqual([Math.min(...moved), Math.max(...moved)]);
            });
        });

    it('accepts only term sets whose degree lists are separately strictly '
        + 'increasing', () => {
            // The header promises "distinct pairs <p[i],q[i]>", but the
            // constructor asserts each list is strictly increasing on its
            // own, which rejects ordinary term sets such as the affine model
            // {(0,0),(1,0),(0,1)}. Preserved from upstream, pinned here.
            check(fc.tuple(
                fc.array(fc.integer({ min: -1, max: 3 }), { maxLength: 3 }),
                fc.array(fc.integer({ min: -1, max: 3 }), { maxLength: 3 })),
                ([xd, yd]) => {
                    const rising = (d: readonly number[]): boolean =>
                        d.length > 0 && d[0] >= 0
                        && d.every((v, i) => i === 0 || v > d[i - 1]);
                    if (xd.length !== yd.length) {
                        expect(() => new ApprPolynomialSpecial3(xd, yd))
                            .toThrowError(
                                'The input arrays must have the same size.');
                    }
                    else if (xd.length === 0) {
                        expect(() => new ApprPolynomialSpecial3(xd, yd))
                            .toThrowError(
                                'The input array must have elements.');
                    }
                    else if (rising(xd) && rising(yd)) {
                        expect(() => new ApprPolynomialSpecial3(xd, yd))
                            .not.toThrow();
                    }
                    else {
                        expect(() => new ApprPolynomialSpecial3(xd, yd))
                            .toThrowError('Degrees must be increasing.');
                    }
                });

            // The documented distinct-pair term set {(0,0),(1,0),(0,1)}.
            expect(() => new ApprPolynomialSpecial3([0, 1, 0], [0, 0, 1]))
                .toThrowError('Degrees must be increasing.');
        });

    it('copyParameters deep-copies the model', () => {
        check(configArb, ([xd, yd, xs, ys, ws]) => {
            const source = new ApprPolynomialSpecial3(xd, yd);
            if (!source.fit(xs.map((x, i) => [x, ys[i], ws[i]]))) { return; }
            const target = new ApprPolynomialSpecial3(xd, yd);
            target.copyParameters(source);

            expect([...target.getParameters()])
                .toEqual([...source.getParameters()]);
            expect(target.getParameters()).not.toBe(source.getParameters());
            expect(target.evaluate(0.25, -0.5))
                .toBe(source.evaluate(0.25, -0.5));

            const copied = [...target.getParameters()];
            source.fit(xs.map((x, i) => [x, ys[i], ws[i] + 10]));
            expect([...target.getParameters()]).toEqual(copied);
        });
    });
});
