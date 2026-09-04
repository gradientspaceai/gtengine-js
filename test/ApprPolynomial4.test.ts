import { describe, expect, it } from 'vitest';
import { ApprPolynomial4 } from '../src/ApprPolynomial4.js';
import { ApprQuery } from '../src/ApprQuery.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

// w = sum_{i,j,k} c[i + (dx+1)*(j + (dy+1)*k)] * x^i * y^j * z^k
function poly(c: readonly number[], dx: number, dy: number, dz: number,
    x: number, y: number, z: number): number {
    let w = 0;
    for (let k = 0; k <= dz; ++k) {
        for (let j = 0; j <= dy; ++j) {
            for (let i = 0; i <= dx; ++i) {
                w += c[i + (dx + 1) * (j + (dy + 1) * k)] *
                    Math.pow(x, i) * Math.pow(y, j) * Math.pow(z, k);
            }
        }
    }
    return w;
}

describe('ApprPolynomial4', () => {
    it('initializes the parameters to zero and reports the sizes', () => {
        const fitter = new ApprPolynomial4(1, 1, 1);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0, 0, 0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(8);
        expect(fitter.getXDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getYDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
        expect(fitter.getZDomain()).toEqual([Number.MAX_VALUE, -Number.MAX_VALUE]);
    });

    it('recovers a trilinear polynomial exactly from data sampled on it', () => {
        const c = [1, -2, 3, 0.5, -1.5, 2.25, 0.75, -0.25];
        const observations: number[][] = [];
        for (let a = 0; a <= 2; ++a) {
            for (let b = 0; b <= 2; ++b) {
                for (let d = 0; d <= 2; ++d) {
                    const x = -1 + a;
                    const y = -1 + b;
                    const z = -1 + d;
                    observations.push([x, y, z, poly(c, 1, 1, 1, x, y, z)]);
                }
            }
        }

        const fitter = new ApprPolynomial4(1, 1, 1);
        expect(fitter.fit(observations)).toBe(true);

        const p = fitter.getParameters();
        expect(p.length).toBe(8);
        for (let i = 0; i < 8; ++i) {
            expect(p[i]).toBeCloseTo(c[i], 10);
        }

        for (const obs of observations) {
            expect(fitter.error(obs)).toBeLessThan(1e-10);
        }

        // The nested-Horner evaluate() agrees with the direct sum.
        for (const [x, y, z] of [[0.3, -0.7, 0.2], [2, 3, -4], [-5, 1, 6]]) {
            expect(fitter.evaluate(x, y, z))
                .toBeCloseTo(poly(c, 1, 1, 1, x, y, z), 8);
        }
    });

    it('recovers a degree (2,1,0) polynomial exactly', () => {
        const dx = 2, dy = 1, dz = 0;
        const c = [0.5, 1.5, -2, 3, -0.5, 0.25];
        const observations: number[][] = [];
        for (let a = 0; a <= 3; ++a) {
            for (let b = 0; b <= 2; ++b) {
                for (let d = 0; d <= 1; ++d) {
                    const x = -1.5 + a;
                    const y = -1 + b;
                    const z = 4 + 2 * d;  // z is irrelevant when dz = 0
                    observations.push([x, y, z,
                        poly(c, dx, dy, dz, x, y, z)]);
                }
            }
        }

        const fitter = new ApprPolynomial4(dx, dy, dz);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        for (let i = 0; i < 6; ++i) {
            expect(p[i]).toBeCloseTo(c[i], 8);
        }
        expect(fitter.getZDomain()).toEqual([4, 6]);
    });

    it('records the domains of the fitted samples', () => {
        const fitter = new ApprPolynomial4(0, 0, 0);
        fitter.fit([[-1, 4, 9, 0], [2, -3, 11, 1], [0.5, 7, -2, 2]]);
        expect(fitter.getXDomain()).toEqual([-1, 2]);
        expect(fitter.getYDomain()).toEqual([-3, 7]);
        expect(fitter.getZDomain()).toEqual([-2, 11]);
        // Degree (0,0,0) fits the constant that is the mean of the heights.
        expect(fitter.getParameters()[0]).toBeCloseTo(1, 12);
    });

    it('finds the least-squares plane in three variables', () => {
        // Heights perturbed off w = 1 + 2x + 3y + 4z symmetrically about the
        // exact plane, so the least-squares fit returns the exact plane.
        const base: number[][] = [];
        for (let a = 0; a <= 1; ++a) {
            for (let b = 0; b <= 1; ++b) {
                for (let d = 0; d <= 1; ++d) {
                    base.push([a, b, d]);
                }
            }
        }
        const observations: number[][] = [];
        for (const [x, y, z] of base) {
            const w = 1 + 2 * x + 3 * y + 4 * z;
            observations.push([x, y, z, w + 0.7]);
            observations.push([x, y, z, w - 0.7]);
        }

        // The model with degrees (1,1,1) contains the plane's terms, so the
        // cross terms must vanish and the linear terms must be exact.
        const fitter = new ApprPolynomial4(1, 1, 1);
        expect(fitter.fit(observations)).toBe(true);
        const p = fitter.getParameters();
        expect(p[0]).toBeCloseTo(1, 8);   // 1
        expect(p[1]).toBeCloseTo(2, 8);   // x
        expect(p[2]).toBeCloseTo(3, 8);   // y
        expect(p[3]).toBeCloseTo(0, 8);   // x*y
        expect(p[4]).toBeCloseTo(4, 8);   // z
        expect(p[5]).toBeCloseTo(0, 8);   // x*z
        expect(p[6]).toBeCloseTo(0, 8);   // y*z
        expect(p[7]).toBeCloseTo(0, 8);   // x*y*z
    });

    it('fails when the system is singular', () => {
        const fitter = new ApprPolynomial4(1, 1, 1);
        expect(fitter.fit([[0, 0, 0, 1], [1, 1, 1, 2]])).toBe(false);
        expect(fitter.getParameters()).toEqual([0, 0, 0, 0, 0, 0, 0, 0]);
    });

    it('reports failure for invalid indices when validation is enabled', () => {
        const fitter = new ApprPolynomial4(0, 0, 0);
        ApprQuery.validateIndices = true;
        try {
            expect(fitter.fit([[0, 0, 0, 1]], [3])).toBe(false);
            expect(fitter.getParameters()).toEqual([0]);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });

    it('copies the parameters between models', () => {
        const source = new ApprPolynomial4(1, 0, 0);
        source.fit([[0, 0, 0, 1], [1, 0, 0, 3], [2, 0, 0, 5]]);
        const target = new ApprPolynomial4(1, 0, 0);
        target.copyParameters(source);
        expect(target.getParameters()).toEqual(source.getParameters());
        expect(target.getXDomain()).toEqual(source.getXDomain());
        expect(target.evaluate(2.5, 0, 0)).toBe(source.evaluate(2.5, 0, 0));

        const copied = target.getParameters().slice();
        source.fit([[0, 0, 0, 0], [1, 0, 0, 0], [2, 0, 0, 1]]);
        expect(target.getParameters()).toEqual(copied);
    });
});

describe('ApprPolynomial4 verification', () => {
    // Strictly increasing, well-separated coordinates keep the Vandermonde
    // normal equations conditioned. The degrees stay small: the system is
    // (dx+1)*(dy+1)*(dz+1) square and the monomial basis is nonrobust.
    const coords = (count: number) =>
        fc.array(finite(0.5, 1.3), { minLength: count, maxLength: count })
            .map(gaps => {
                let acc = -1.5;
                return gaps.map(g => (acc += g));
            });

    const configArb = fc.tuple(fc.integer({ min: 0, max: 2 }),
        fc.integer({ min: 0, max: 1 }), fc.integer({ min: 0, max: 1 }),
        fc.integer({ min: 0, max: 1 }))
        .chain(([xd, yd, zd, extra]) => {
            const size = (xd + 1) * (yd + 1) * (zd + 1);
            return fc.tuple(fc.constant(xd), fc.constant(yd),
                fc.constant(zd),
                fc.array(finite(-3, 3), { minLength: size, maxLength: size })
                    .filter(c => c.some(v => Math.abs(v) > 0.1)),
                coords(xd + 1 + extra), coords(yd + 1), coords(zd + 1));
        });

    const evaluatePoly = (c: readonly number[], xd: number, yd: number,
        zd: number, x: number, y: number, z: number): number => {
        let w = 0;
        for (let k = 0; k <= zd; ++k) {
            for (let j = 0; j <= yd; ++j) {
                for (let i = 0; i <= xd; ++i) {
                    w += c[i + (xd + 1) * (j + (yd + 1) * k)]
                        * Math.pow(x, i) * Math.pow(y, j) * Math.pow(z, k);
                }
            }
        }
        return w;
    };

    const grid = (c: readonly number[], xd: number, yd: number, zd: number,
        xs: readonly number[], ys: readonly number[],
        zs: readonly number[]): number[][] => {
        const observations: number[][] = [];
        for (const x of xs) {
            for (const y of ys) {
                for (const z of zs) {
                    observations.push([x, y, z,
                        evaluatePoly(c, xd, yd, zd, x, y, z)]);
                }
            }
        }
        return observations;
    };

    it('recovers a polynomial its samples lie on', () => {
        check(configArb, ([xd, yd, zd, c, xs, ys, zs]) => {
            const observations = grid(c, xd, yd, zd, xs, ys, zs);
            const fitter = new ApprPolynomial4(xd, yd, zd);
            expect(fitter.fit(observations)).toBe(true);

            const fitted = fitter.getParameters();
            expect(fitted.length).toBe((xd + 1) * (yd + 1) * (zd + 1));
            const scale = c.reduce((u, v) => u + Math.abs(v), 0) + 1;
            for (let i = 0; i < fitted.length; ++i) {
                expectClose(fitted[i], c[i], 1e-6 * scale, 1e-6);
            }

            // The index map is c[i][j][k] = parameters[i+(d0+1)*(j+(d1+1)*k)]
            // and evaluate()/error() agree with it on the samples.
            for (const [x, y, z, w] of observations) {
                expectClose(fitter.evaluate(x, y, z), w, 1e-6 * scale, 1e-6);
                expectClose(fitter.error([x, y, z, w]), 0, 1e-5 * scale, 0);
            }
        }, 40);
    });

    it('evaluate matches a direct monomial sum over the parameters', () => {
        check(fc.tuple(configArb, finite(-3, 3), finite(-3, 3), finite(-3, 3)),
            ([[xd, yd, zd, c, xs, ys, zs], x, y, z]) => {
                const fitter = new ApprPolynomial4(xd, yd, zd);
                expect(fitter.fit(grid(c, xd, yd, zd, xs, ys, zs))).toBe(true);
                const fitted = fitter.getParameters();

                const direct = evaluatePoly(fitted, xd, yd, zd, x, y, z);
                const scale = fitted.reduce((u, v) => u + Math.abs(v), 0) + 1;
                const power = Math.pow(
                    Math.max(Math.abs(x), Math.abs(y), Math.abs(z), 1),
                    xd + yd + zd);
                expectClose(fitter.evaluate(x, y, z), direct,
                    1e-8 * scale * power, 1e-8);
            }, 40);
    });

    it('satisfies the least-squares normal equations', () => {
        check(fc.tuple(configArb, fc.array(finite(-5, 5),
            { minLength: 13, maxLength: 13 })),
            ([[xd, yd, zd, c, xs, ys, zs], ws]) => {
                const observations = grid(c, xd, yd, zd, xs, ys, zs)
                    .map((o, i) => [o[0], o[1], o[2], ws[i % ws.length]]);
                const fitter = new ApprPolynomial4(xd, yd, zd);
                if (!fitter.fit(observations)) { return; }

                for (let k = 0; k <= zd; ++k) {
                    for (let j = 0; j <= yd; ++j) {
                        for (let i = 0; i <= xd; ++i) {
                            let residual = 0, scale = 0;
                            for (const [x, y, z, w] of observations) {
                                const basis = Math.pow(x, i) * Math.pow(y, j)
                                    * Math.pow(z, k);
                                residual += basis
                                    * (fitter.evaluate(x, y, z) - w);
                                scale += Math.abs(basis) * (Math.abs(w) + 1);
                            }
                            expect(Math.abs(residual))
                                .toBeLessThanOrEqual(1e-6 + 1e-6 * scale);
                        }
                    }
                }
            }, 30);
    });

    it('reports the minimum required and accumulates all three domains',
        () => {
            check(configArb, ([xd, yd, zd, c, xs, ys, zs]) => {
                const fitter = new ApprPolynomial4(xd, yd, zd);
                expect(fitter.getMinimumRequired())
                    .toBe((xd + 1) * (yd + 1) * (zd + 1));
                expect(fitter.getXDomain()[0]).toBe(Number.MAX_VALUE);
                expect(fitter.getZDomain()[1]).toBe(-Number.MAX_VALUE);

                fitter.fit(grid(c, xd, yd, zd, xs, ys, zs));
                expect(fitter.getXDomain()).toEqual(
                    [Math.min(...xs), Math.max(...xs)]);
                expect(fitter.getYDomain()).toEqual(
                    [Math.min(...ys), Math.max(...ys)]);
                expect(fitter.getZDomain()).toEqual(
                    [Math.min(...zs), Math.max(...zs)]);

                // Upstream never resets the domains between fits.
                fitter.fit(grid(c, xd, yd, zd, xs, ys, zs.map(z => z + 40)));
                expect(fitter.getZDomain()).toEqual(
                    [Math.min(...zs), Math.max(...zs) + 40]);
            }, 40);
        });

    it('copyParameters deep-copies the model', () => {
        check(configArb, ([xd, yd, zd, c, xs, ys, zs]) => {
            const source = new ApprPolynomial4(xd, yd, zd);
            expect(source.fit(grid(c, xd, yd, zd, xs, ys, zs))).toBe(true);
            const target = new ApprPolynomial4(xd, yd, zd);
            target.copyParameters(source);

            expect([...target.getParameters()])
                .toEqual([...source.getParameters()]);
            expect(target.getParameters()).not.toBe(source.getParameters());
            expect(target.getZDomain()).toEqual(source.getZDomain());

            const copied = [...target.getParameters()];
            source.fit(grid(c.map(v => v + 1), xd, yd, zd, xs, ys, zs));
            expect([...target.getParameters()]).toEqual(copied);
        }, 40);
    });

    it('zeroes the parameters when the indices are rejected', () => {
        ApprQuery.validateIndices = true;
        try {
            check(configArb, ([xd, yd, zd, c, xs, ys, zs]) => {
                const fitter = new ApprPolynomial4(xd, yd, zd);
                expect(fitter.fit(grid(c, xd, yd, zd, xs, ys, zs))).toBe(true);
                expect(fitter.fit([[0, 0, 0, 0]], [0, 99])).toBe(false);
                expect([...fitter.getParameters()]).toEqual(
                    new Array<number>((xd + 1) * (yd + 1) * (zd + 1)).fill(0));
            }, 30);
        }
        finally {
            ApprQuery.validateIndices = false;
        }
    });
});
