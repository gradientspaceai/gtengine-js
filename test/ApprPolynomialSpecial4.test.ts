import { describe, expect, it } from 'vitest';
import { ApprPolynomialSpecial4 } from '../src/ApprPolynomialSpecial4.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

const grid = [-1, -0.5, 0, 0.5, 1];

function cube(f: (x: number, y: number, z: number) => number): number[][] {
    const observations: number[][] = [];
    for (const x of grid) {
        for (const y of grid) {
            for (const z of grid) {
                observations.push([x, y, z, f(x, y, z)]);
            }
        }
    }
    return observations;
}

describe('ApprPolynomialSpecial4', () => {
    it('validates the degrees in the constructor', () => {
        expect(() => new ApprPolynomialSpecial4([0, 1], [0, 1], [0]))
            .toThrowError('The input arrays must have the same size.');
        expect(() => new ApprPolynomialSpecial4([], [], []))
            .toThrowError('The input array must have elements.');
        expect(() => new ApprPolynomialSpecial4([0, 1], [0, 1], [1, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial4([0, 2], [0, 1], [3, 1]))
            .toThrowError('Degrees must be increasing.');
        expect(() => new ApprPolynomialSpecial4([0, 1, 2], [0, 1, 3], [0, 2, 4]))
            .not.toThrow();
    });

    it('initializes the parameters to zero', () => {
        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 2], [0, 3]);
        expect(fitter.getParameters()).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getZDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('reproduces data from w = 7 + 3*x*y*z with the term set {1, x*y*z}', () => {
        const f = (x: number, y: number, z: number): number =>
            7 + 3 * x * y * z;
        const observations = cube(f);

        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        expect(fitter.getXDomain()).toEqual([-1, 1]);
        expect(fitter.getYDomain()).toEqual([-1, 1]);
        expect(fitter.getZDomain()).toEqual([-1, 1]);

        // With the identity transform on x, y and z, w' = x'*y'*z', so the
        // transformed coefficients are (0, 1).
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(0, 10);
        expect(p[1]).toBeCloseTo(1, 10);

        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1], obs[2]))
                .toBeCloseTo(obs[3], 10);
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }
        for (const [x, y, z] of [[0.3, -0.4, 0.7], [2, 1, -3]]) {
            expect(fitter.evaluate(x, y, z)).toBeCloseTo(f(x, y, z), 8);
        }
    });

    it('reproduces data from the three-term set {1, x*y*z, x^2*y^2*z^2}', () => {
        const f = (x: number, y: number, z: number): number => {
            const t = x * y * z;
            return 2 - t + 5 * t * t;
        };
        const observations = cube(f);

        const fitter = new ApprPolynomialSpecial4([0, 1, 2], [0, 1, 2], [0, 1, 2]);
        expect(fitter.fit(observations)).toBe(true);
        for (const obs of observations) {
            expect(fitter.evaluate(obs[0], obs[1], obs[2]))
                .toBeCloseTo(obs[3], 8);
        }
        for (const [x, y, z] of [[0.3, -0.4, 0.9], [0.8, 0.8, 0.8]]) {
            expect(fitter.evaluate(x, y, z)).toBeCloseTo(f(x, y, z), 6);
        }
    });

    it('produces a genuine least-squares fit for data off the model', () => {
        const observations: number[][] = [];
        for (const [x, y, z] of cube(() => 0).map((o) => [o[0], o[1], o[2]])) {
            const w = 7 + 3 * x * y * z;
            observations.push([x, y, z, w + 0.6]);
            observations.push([x, y, z, w - 0.6]);
        }

        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        expect(fitter.fit(observations)).toBe(true);
        for (const x of grid) {
            for (const y of grid) {
                for (const z of grid) {
                    expect(fitter.evaluate(x, y, z))
                        .toBeCloseTo(7 + 3 * x * y * z, 8);
                }
            }
        }
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 0, 0], [1, 1, 1, 1]], [0, 6])).toBe(false);
            expect(fitter.getParameters()).toEqual([0, 0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const observations = cube((x, y, z) => 7 + 3 * x * y * z);

        const source = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        source.fit(observations);
        const target = new ApprPolynomialSpecial4([0, 1], [0, 1], [0, 1]);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.getYDomain()).toEqual(source.getYDomain());
        expect(target.getZDomain()).toEqual(source.getZDomain());
        expect(target.evaluate(0.3, 0.4, 0.5))
            .toBe(source.evaluate(0.3, 0.4, 0.5));

        const copied = target.getParameters().slice();
        const value = target.evaluate(0.3, 0.4, 0.5);
        source.fit(cube((x, y, z) => 1 - x * y * z));
        expect(target.getParameters()).toEqual(copied);
        expect(target.evaluate(0.3, 0.4, 0.5)).toBe(value);
    });
});

describe('ApprPolynomialSpecial4 verification', () => {
    // Upstream asserts that each of the three degree lists is separately
    // strictly increasing, so only "diagonal" term sets are admissible.
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

    const configArb = fc.tuple(degreesArb, degreesArb, degreesArb)
        .chain(([xd, yd, zd]) => {
            const n = Math.min(xd.length, yd.length, zd.length);
            const count = n + 3;
            return fc.tuple(fc.constant(xd.slice(0, n)),
                fc.constant(yd.slice(0, n)), fc.constant(zd.slice(0, n)),
                spread(count, -2), spread(count, -1.5), spread(count, -1),
                fc.array(finite(-4, 4), { minLength: count, maxLength: count })
                    .filter(ws => Math.max(...ws) - Math.min(...ws) > 0.5));
        });

    const mapTo = (v: number, lo: number, hi: number): number =>
        -1 + (2 * (v - lo)) / (hi - lo);

    it('is a least-squares fit in the transformed coordinates', () => {
        // The samples are mapped to [-1,1]^4 before fitting, so the residual
        // of the transformed model is orthogonal to every basis function
        // (x')^{p_i} (y')^{q_i} (z')^{r_i}.
        check(configArb, ([xd, yd, zd, xs, ys, zs, ws]) => {
            const observations = xs.map(
                (x, i) => [x, ys[i], zs[i], ws[i]]);
            const fitter = new ApprPolynomialSpecial4(xd, yd, zd);
            if (!fitter.fit(observations)) { return; }
            const c = fitter.getParameters();

            const xlo = Math.min(...xs), xhi = Math.max(...xs);
            const ylo = Math.min(...ys), yhi = Math.max(...ys);
            const zlo = Math.min(...zs), zhi = Math.max(...zs);
            const wlo = Math.min(...ws), whi = Math.max(...ws);
            const basisAt = (i: number, xp: number, yp: number,
                zp: number): number => Math.pow(xp, xd[i])
                    * Math.pow(yp, yd[i]) * Math.pow(zp, zd[i]);
            const model = (xp: number, yp: number, zp: number): number =>
                c.reduce((u, v, i) => u + v * basisAt(i, xp, yp, zp), 0);

            for (let i = 0; i < xd.length; ++i) {
                let residual = 0, scale = 0;
                for (const [x, y, z, w] of observations) {
                    const xp = mapTo(x, xlo, xhi);
                    const yp = mapTo(y, ylo, yhi);
                    const zp = mapTo(z, zlo, zhi);
                    const basis = basisAt(i, xp, yp, zp);
                    residual += basis
                        * (model(xp, yp, zp) - mapTo(w, wlo, whi));
                    scale += Math.abs(basis) * 2;
                }
                expect(Math.abs(residual))
                    .toBeLessThanOrEqual(1e-8 + 1e-8 * scale);
            }
        });
    });

    it('evaluate applies the stored transform in both directions', () => {
        check(fc.tuple(configArb, finite(-4, 4), finite(-4, 4), finite(-4, 4)),
            ([[xd, yd, zd, xs, ys, zs, ws], x, y, z]) => {
                const fitter = new ApprPolynomialSpecial4(xd, yd, zd);
                if (!fitter.fit(xs.map(
                    (u, i) => [u, ys[i], zs[i], ws[i]]))) { return; }
                const c = fitter.getParameters();

                const wlo = Math.min(...ws), whi = Math.max(...ws);
                const xp = mapTo(x, Math.min(...xs), Math.max(...xs));
                const yp = mapTo(y, Math.min(...ys), Math.max(...ys));
                const zp = mapTo(z, Math.min(...zs), Math.max(...zs));
                const wp = c.reduce((u, v, i) => u + v * Math.pow(xp, xd[i])
                    * Math.pow(yp, yd[i]) * Math.pow(zp, zd[i]), 0);
                const expected = ((wp + 1) * (whi - wlo)) / 2 + wlo;

                expectClose(fitter.evaluate(x, y, z), expected,
                    1e-7 * (1 + Math.abs(expected)), 1e-7);
            });
    });

    it('the model error is |evaluate(x,y,z) - w|', () => {
        check(fc.tuple(configArb, finite(-4, 4), finite(-4, 4), finite(-4, 4),
            finite(-5, 5)), ([[xd, yd, zd, xs, ys, zs, ws], x, y, z, w]) => {
                const fitter = new ApprPolynomialSpecial4(xd, yd, zd);
                if (!fitter.fit(xs.map(
                    (u, i) => [u, ys[i], zs[i], ws[i]]))) { return; }
                expect(fitter.error([x, y, z, w]))
                    .toBe(Math.abs(fitter.evaluate(x, y, z) - w));
            });
    });

    it('reports the term count and the domains of the most recent fit',
        () => {
            check(configArb, ([xd, yd, zd, xs, ys, zs, ws]) => {
                const fitter = new ApprPolynomialSpecial4(xd, yd, zd);
                expect(fitter.getMinimumRequired()).toBe(xd.length);
                expect(fitter.getXDomain()[0]).toBe(Number.MAX_VALUE);
                expect(fitter.getZDomain()[1]).toBe(-Number.MAX_VALUE);

                fitter.fit(xs.map((x, i) => [x, ys[i], zs[i], ws[i]]));
                expect(fitter.getXDomain())
                    .toEqual([Math.min(...xs), Math.max(...xs)]);
                expect(fitter.getYDomain())
                    .toEqual([Math.min(...ys), Math.max(...ys)]);
                expect(fitter.getZDomain())
                    .toEqual([Math.min(...zs), Math.max(...zs)]);

                // The domains are assigned, not accumulated.
                const moved = zs.map(z => z + 70);
                fitter.fit(xs.map((x, i) => [x, ys[i], moved[i], ws[i]]));
                expect(fitter.getZDomain())
                    .toEqual([Math.min(...moved), Math.max(...moved)]);
            });
        });

    it('accepts only term sets whose degree lists are separately strictly '
        + 'increasing', () => {
            // As in ApprPolynomialSpecial3, the header promises "distinct
            // triples" but the constructor demands each list rise on its own.
            check(fc.tuple(
                fc.array(fc.integer({ min: -1, max: 3 }), { maxLength: 3 }),
                fc.array(fc.integer({ min: -1, max: 3 }), { maxLength: 3 }),
                fc.array(fc.integer({ min: -1, max: 3 }), { maxLength: 3 })),
                ([xd, yd, zd]) => {
                    const rising = (d: readonly number[]): boolean =>
                        d.length > 0 && d[0] >= 0
                        && d.every((v, i) => i === 0 || v > d[i - 1]);
                    if (xd.length !== yd.length || xd.length !== zd.length) {
                        expect(() => new ApprPolynomialSpecial4(xd, yd, zd))
                            .toThrowError(
                                'The input arrays must have the same size.');
                    }
                    else if (xd.length === 0) {
                        expect(() => new ApprPolynomialSpecial4(xd, yd, zd))
                            .toThrowError(
                                'The input array must have elements.');
                    }
                    else if (rising(xd) && rising(yd) && rising(zd)) {
                        expect(() => new ApprPolynomialSpecial4(xd, yd, zd))
                            .not.toThrow();
                    }
                    else {
                        expect(() => new ApprPolynomialSpecial4(xd, yd, zd))
                            .toThrowError('Degrees must be increasing.');
                    }
                });

            // The documented distinct-triple set {(0,0,0),(1,0,0),(0,1,0)}.
            expect(() => new ApprPolynomialSpecial4([0, 1, 0], [0, 0, 1],
                [0, 0, 0])).toThrowError('Degrees must be increasing.');
        });

    it('copyParameters deep-copies the model', () => {
        check(configArb, ([xd, yd, zd, xs, ys, zs, ws]) => {
            const source = new ApprPolynomialSpecial4(xd, yd, zd);
            if (!source.fit(xs.map(
                (x, i) => [x, ys[i], zs[i], ws[i]]))) { return; }
            const target = new ApprPolynomialSpecial4(xd, yd, zd);
            target.copyParameters(source);

            expect([...target.getParameters()])
                .toEqual([...source.getParameters()]);
            expect(target.getParameters()).not.toBe(source.getParameters());
            expect(target.evaluate(0.25, -0.5, 0.75))
                .toBe(source.evaluate(0.25, -0.5, 0.75));

            const copied = [...target.getParameters()];
            source.fit(xs.map((x, i) => [x, ys[i], zs[i], ws[i] + 10]));
            expect([...target.getParameters()]).toEqual(copied);
        });
    });
});
