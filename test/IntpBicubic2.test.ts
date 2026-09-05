import { describe, it, expect } from 'vitest';
import { IntpBicubic2 } from '../src/IntpBicubic2.js';
import { check, expectClose, fc, finite, positive, scaled, wellScaled } from './helpers/arbitraries.js';

// Build the row-major sample array F[c + xBound*r] = f(xMin + c*dx, yMin + r*dy).
function makeSamples(xBound: number, yBound: number, xMin: number, dx: number,
    yMin: number, dy: number, f: (x: number, y: number) => number): number[] {
    const F: number[] = [];
    for (let r = 0; r < yBound; ++r) {
        for (let c = 0; c < xBound; ++c) {
            F.push(f(xMin + c * dx, yMin + r * dy));
        }
    }
    return F;
}

// Independent 1D evaluators written in the classic polynomial form rather
// than as the blending-matrix product used by the port.
function catmullRom1(f0: number, f1: number, f2: number, f3: number, t: number): number {
    return 0.5 * (2 * f1 + (-f0 + f2) * t
        + (2 * f0 - 5 * f1 + 4 * f2 - f3) * t * t
        + (-f0 + 3 * f1 - 3 * f2 + f3) * t * t * t);
}

function bspline1(f0: number, f1: number, f2: number, f3: number, t: number): number {
    const s = 1 - t;
    return (s * s * s * f0
        + (3 * t * t * t - 6 * t * t + 4) * f1
        + (-3 * t * t * t + 3 * t * t + 3 * t + 1) * f2
        + t * t * t * f3) / 6;
}

// Independent tensor-product evaluation with upstream's index clamping.
function reference(F: number[], xBound: number, yBound: number, xMin: number,
    dx: number, yMin: number, dy: number, catmullRom: boolean,
    x: number, y: number): number {
    const xIndex = (x - xMin) / dx;
    const yIndex = (y - yMin) / dy;
    let ix = Math.trunc(xIndex);
    ix = ix < 0 ? 0 : (ix >= xBound ? xBound - 1 : ix);
    let iy = Math.trunc(yIndex);
    iy = iy < 0 ? 0 : (iy >= yBound ? yBound - 1 : iy);
    const s = xIndex - ix;
    const t = yIndex - iy;
    const clamp = (i: number, n: number): number => (i < 0 ? 0 : (i > n - 1 ? n - 1 : i));
    const oneD = catmullRom ? catmullRom1 : bspline1;

    const rows: number[] = [];
    for (let row = 0; row < 4; ++row) {
        const r = clamp(iy - 1 + row, yBound);
        const v: number[] = [];
        for (let col = 0; col < 4; ++col) {
            const c = clamp(ix - 1 + col, xBound);
            v.push(F[c + xBound * r]);
        }
        rows.push(oneD(v[0], v[1], v[2], v[3], s));
    }
    return oneD(rows[0], rows[1], rows[2], rows[3], t);
}

describe('IntpBicubic2', () => {
    it('throws for invalid inputs', () => {
        const F9 = new Array<number>(9).fill(0);
        expect(() => new IntpBicubic2(2, 3, 0, 1, 0, 1, F9, true)).toThrow('Invalid input.');
        expect(() => new IntpBicubic2(3, 2, 0, 1, 0, 1, F9, true)).toThrow('Invalid input.');
        expect(() => new IntpBicubic2(3, 3, 0, 0, 0, 1, F9, true)).toThrow('Invalid input.');
        expect(() => new IntpBicubic2(3, 3, 0, 1, 0, -1, F9, false)).toThrow('Invalid input.');
        expect(() => new IntpBicubic2(4, 3, 0, 1, 0, 1, F9, false)).toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const F = makeSamples(4, 5, -1, 0.5, 2, 0.25, (x, y) => x + y);
        const interp = new IntpBicubic2(4, 5, -1, 0.5, 2, 0.25, F, true);
        expect(interp.getXBound()).toBe(4);
        expect(interp.getYBound()).toBe(5);
        expect(interp.getQuantity()).toBe(20);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXMax()).toBeCloseTo(0.5, 14);
        expect(interp.getXSpacing()).toBe(0.5);
        expect(interp.getYMin()).toBe(2);
        expect(interp.getYMax()).toBeCloseTo(3, 14);
        expect(interp.getYSpacing()).toBe(0.25);
    });

    it('interpolates the samples exactly with Catmull-Rom blending', () => {
        const xBound = 6, yBound = 5;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25;
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.3 * y);
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpBicubic2(xBound, yBound, xMin, dx, yMin, dy, F, true);

        for (let r = 0; r < yBound; ++r) {
            for (let c = 0; c < xBound; ++c) {
                const value = interp.evaluate(xMin + c * dx, yMin + r * dy);
                expect(value).toBeCloseTo(F[c + xBound * r], 12);
            }
        }
    });

    it('smooths rather than interpolates with B-spline blending', () => {
        // At a grid point the B-spline blend applies the (1,4,1)/6 stencil in
        // each direction. For the samples f(x,y) = x^2 (independent of y) at
        // the grid point x = 3 this is (4 + 4*9 + 16)/6 = 28/3.
        const xBound = 8, yBound = 8;
        const F = makeSamples(xBound, yBound, 0, 1, 0, 1, (x) => x * x);
        const interp = new IntpBicubic2(xBound, yBound, 0, 1, 0, 1, F, false);
        expect(interp.evaluate(3, 3)).toBeCloseTo(28 / 3, 12);
        // First and second derivatives of a quadratic are still exact.
        expect(interp.evaluate(1, 0, 3, 3)).toBeCloseTo(6, 12);
        expect(interp.evaluate(2, 0, 3, 3)).toBeCloseTo(2, 12);
        expect(interp.evaluate(3, 0, 3, 3)).toBeCloseTo(0, 12);
    });

    it('is exact for a tensor-product quadratic with Catmull-Rom blending', () => {
        // f(x,y) = g(x)*h(y) with g and h quadratic; the Catmull-Rom blend
        // reproduces quadratics, so both the values and the derivatives are
        // exact away from the boundary rows/columns where indices clamp.
        const g = (x: number): number => 1 + 2 * x - 0.5 * x * x;
        const dg = (x: number): number => 2 - x;
        const h = (y: number): number => 2 - y + 0.75 * y * y;
        const dh = (y: number): number => -1 + 1.5 * y;

        const xBound = 8, yBound = 8;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, (x, y) => g(x) * h(y));
        const interp = new IntpBicubic2(xBound, yBound, xMin, dx, yMin, dy, F, true);

        for (const [x, y] of [[-0.4, 2.4], [0.3, 2.9], [1.1, 3.1]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(g(x) * h(y), 11);
            expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(g(x) * h(y), 11);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(dg(x) * h(y), 11);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(g(x) * dh(y), 11);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(dg(x) * dh(y), 11);
            expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(-1 * h(y), 11);
            expect(interp.evaluate(0, 2, x, y)).toBeCloseTo(g(x) * 1.5, 11);
            expect(interp.evaluate(3, 0, x, y)).toBeCloseTo(0, 10);
            expect(interp.evaluate(0, 3, x, y)).toBeCloseTo(0, 10);
        }
    });

    it('reproduces a smoothed quadratic with B-spline blending', () => {
        // The cubic B-spline blend applied to samples of a quadratic q gives
        // q + (h^2/6) q'' in each direction (a known shift), so the tensor
        // product of the two shifted quadratics is reproduced exactly.
        const g = (x: number): number => 1 + 2 * x - 0.5 * x * x;
        const h = (y: number): number => 2 - y + 0.75 * y * y;
        const xBound = 8, yBound = 8;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, (x, y) => g(x) * h(y));
        const interp = new IntpBicubic2(xBound, yBound, xMin, dx, yMin, dy, F, false);

        const gs = (x: number): number => g(x) + (dx * dx / 6) * (-1);
        const hs = (y: number): number => h(y) + (dy * dy / 6) * 1.5;
        for (const [x, y] of [[-0.4, 2.4], [0.3, 2.9], [1.1, 3.1]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(gs(x) * hs(y), 11);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo((2 - x) * hs(y), 11);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(gs(x) * (-1 + 1.5 * y), 11);
        }
    });

    it('is C1 across cell boundaries for both blends', () => {
        const xBound = 8, yBound = 8;
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.2 * y);
        const F = makeSamples(xBound, yBound, 0, 0.5, 0, 0.5, f);
        const eps = 1e-7;
        for (const catmullRom of [true, false]) {
            const interp = new IntpBicubic2(xBound, yBound, 0, 0.5, 0, 0.5, F, catmullRom);
            for (const y of [0.7, 1.3, 2.1]) {
                for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                    const left = interp.evaluate(xo, yo, 1 - eps, y);
                    const right = interp.evaluate(xo, yo, 1 + eps, y);
                    expect(right - left).toBeCloseTo(0, 5);
                }
            }
            for (const x of [0.7, 1.3, 2.1]) {
                for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                    const below = interp.evaluate(xo, yo, x, 1.5 - eps);
                    const above = interp.evaluate(xo, yo, x, 1.5 + eps);
                    expect(above - below).toBeCloseTo(0, 5);
                }
            }
        }
    });

    it('agrees with an independent tensor-product evaluation, in and out of the domain', () => {
        const xBound = 7, yBound = 6;
        const xMin = -2, dx = 0.75, yMin = 1, dy = 0.4;
        let seed = 987654321;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const F: number[] = [];
        for (let i = 0; i < xBound * yBound; ++i) {
            F.push(2 * rand() - 1);
        }

        for (const catmullRom of [true, false]) {
            const interp = new IntpBicubic2(xBound, yBound, xMin, dx, yMin, dy, F, catmullRom);
            for (let trial = 0; trial < 200; ++trial) {
                // Sample a range wider than the domain to exercise the index
                // clamping on all four sides.
                const x = xMin + (1.6 * rand() - 0.3) * dx * (xBound - 1);
                const y = yMin + (1.6 * rand() - 0.3) * dy * (yBound - 1);
                const expected = reference(F, xBound, yBound, xMin, dx, yMin, dy,
                    catmullRom, x, y);
                expect(interp.evaluate(x, y)).toBeCloseTo(expected, 10);
                expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(expected, 10);
            }
        }
    });

    it('matches finite differences of its own values', () => {
        const xBound = 8, yBound = 8;
        const xMin = 0, dx = 0.5, yMin = 0, dy = 0.25;
        const f = (x: number, y: number): number => Math.cos(x) * (1 + y * y);
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        for (const catmullRom of [true, false]) {
            const interp = new IntpBicubic2(xBound, yBound, xMin, dx, yMin, dy, F, catmullRom);
            const eps = 1e-5;
            for (const [x, y] of [[0.7, 0.6], [1.3, 1.1], [2.2, 0.9]]) {
                const dfx = (interp.evaluate(x + eps, y) - interp.evaluate(x - eps, y)) / (2 * eps);
                expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(dfx, 6);
                const dfy = (interp.evaluate(x, y + eps) - interp.evaluate(x, y - eps)) / (2 * eps);
                expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(dfy, 6);
                const d2fx = (interp.evaluate(x + eps, y) - 2 * interp.evaluate(x, y)
                    + interp.evaluate(x - eps, y)) / (eps * eps);
                expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(d2fx, 3);
            }
        }
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const F = makeSamples(4, 4, 0, 1, 0, 1, (x, y) => x * y);
        const interp = new IntpBicubic2(4, 4, 0, 1, 0, 1, F, true);
        expect(interp.evaluate(4, 0, 1.5, 1.5)).toBe(0);
        expect(interp.evaluate(0, 4, 1.5, 1.5)).toBe(0);
        expect(interp.evaluate(-1, 0, 1.5, 1.5)).toBe(0);
    });

    it('scales derivatives by the inverse sample spacing', () => {
        // The same samples on grids with spacing 1 and 4 give derivatives that
        // differ exactly by the ratio of the spacings.
        const F = makeSamples(6, 6, 0, 1, 0, 1, (x, y) => x * x + x * y);
        const unit = new IntpBicubic2(6, 6, 0, 1, 0, 1, F, true);
        const wide = new IntpBicubic2(6, 6, 0, 4, 0, 4, F, true);
        for (const [c, r] of [[2.3, 3.1], [1.7, 2.5]]) {
            expect(wide.evaluate(4 * c, 4 * r)).toBeCloseTo(unit.evaluate(c, r), 11);
            expect(wide.evaluate(1, 0, 4 * c, 4 * r))
                .toBeCloseTo(unit.evaluate(1, 0, c, r) / 4, 11);
            expect(wide.evaluate(1, 1, 4 * c, 4 * r))
                .toBeCloseTo(unit.evaluate(1, 1, c, r) / 16, 11);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against IntpBicubic2.h.
// ---------------------------------------------------------------------------
describe('IntpBicubic2 verification', () => {
    const XB = 7, YB = 6;
    const geometry = () => fc.tuple(finite(-3, 3), positive(3, 0.25),
        finite(-3, 3), positive(3, 0.25));
    type Geom = [number, number, number, number];
    const samples = () => fc.array(wellScaled(-6, 6),
        { minLength: XB * YB, maxLength: XB * YB });
    const build = (F: number[], g: Geom, catmullRom: boolean) =>
        new IntpBicubic2(XB, YB, g[0], g[1], g[2], g[3], F, catmullRom);

    // A cell index and offset that keep the 4x4 stencil strictly inside the
    // grid, so no index clamping happens and the scheme is the plain uniform
    // one. Requires 1 <= ix <= xBound - 3.
    const interiorParam = (bound: number) =>
        fc.tuple(fc.integer({ min: 1, max: bound - 3 }), scaled(0.02, 0.98));

    it('reproduces biquadratics with the Catmull-Rom blend in the interior', () => {
        // The Catmull-Rom tangent at a knot is the centered difference, which
        // is exact for polynomials of degree <= 2 but not for cubics, so the
        // scheme has quadratic precision in each variable.
        check(fc.tuple(fc.array(wellScaled(-2, 2), { minLength: 9, maxLength: 9 }),
            geometry(), interiorParam(XB), interiorParam(YB)),
            ([a, gg, [ix, u], [iy, v]]) => {
                const g = gg as Geom;
                const f = (x: number, y: number) => {
                    let s = 0;
                    for (let i = 0; i <= 2; ++i) {
                        for (let j = 0; j <= 2; ++j) {
                            s += a[3 * i + j] * Math.pow(x, i) * Math.pow(y, j);
                        }
                    }
                    return s;
                };
                const F: number[] = [];
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        F.push(f(g[0] + g[1] * c, g[2] + g[3] * r));
                    }
                }
                const x = g[0] + g[1] * (ix + u);
                const y = g[2] + g[3] * (iy + v);
                // The samples reach |f| ~ 2*(1 + 21 + 441) ~ 900 near the far
                // corner, and the blend recombines 16 of them, so a few
                // hundred ulps of that scale is the achievable accuracy.
                expectClose(build(F, g, true).evaluate(x, y), f(x, y), 1e-9, 1e-11);
            });
    });

    it('reproduces linear data with the B-spline blend in the interior', () => {
        // The uniform cubic B-spline with the samples used directly as control
        // points is the Schoenberg operator, which has linear precision.
        check(fc.tuple(fc.array(wellScaled(-3, 3), { minLength: 3, maxLength: 3 }),
            geometry(), interiorParam(XB), interiorParam(YB)),
            ([k, gg, [ix, u], [iy, v]]) => {
                const g = gg as Geom;
                const f = (x: number, y: number) => k[0] + k[1] * x + k[2] * y;
                const F: number[] = [];
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        F.push(f(g[0] + g[1] * c, g[2] + g[3] * r));
                    }
                }
                const x = g[0] + g[1] * (ix + u);
                const y = g[2] + g[3] * (iy + v);
                expectClose(build(F, g, false).evaluate(x, y), f(x, y), 1e-10, 1e-12);
            });
    });

    it('interpolates the samples at interior grid points (Catmull-Rom only)', () => {
        // At an integer parameter the Catmull-Rom blend column is (0,1,0,0),
        // so the stencil collapses onto the sample itself. The B-spline blend
        // column is (1,4,1,0)/6 and only approximates.
        check(fc.tuple(samples(), geometry()), ([F, gg]) => {
            const g = gg as Geom;
            const cr = build(F, g, true);
            for (let iy = 1; iy < YB - 1; ++iy) {
                for (let ix = 1; ix < XB - 1; ++ix) {
                    expectClose(cr.evaluate(g[0] + g[1] * ix, g[2] + g[3] * iy),
                        F[ix + XB * iy], 1e-10, 1e-12);
                }
            }
        });
    });

    it('reproduces constants everywhere with both blends', () => {
        // Every blend column but the first sums to zero and the first sums to
        // one, so (M*U) is a partition of unity for any parameter, including
        // the out-of-range parameters produced outside the domain.
        check(fc.tuple(wellScaled(-5, 5), geometry(), fc.boolean(),
            scaled(-1.5, 2.5), scaled(-1.5, 2.5)), ([c, gg, cm, s, t]) => {
            const g = gg as Geom;
            const F = new Array<number>(XB * YB).fill(c);
            const x = g[0] + s * g[1] * (XB - 1);
            const y = g[2] + t * g[3] * (YB - 1);
            // Upstream clamps the cell index but not the cell parameter, so
            // outside the domain the monomials (1, p, p^2, p^3) grow without
            // bound and the partition-of-unity cancellation is only as good as
            // the largest weight. The tolerance is therefore |c| times the
            // largest product of monomials times a few hundred ulps.
            const pu = Math.abs((x - g[0]) / g[1]) + 1;
            const pv = Math.abs((y - g[2]) / g[3]) + 1;
            const tol = 1e-12 + Math.abs(c) * 1e-13 * pu * pu * pu * pv * pv * pv;
            expectClose(build(F, g, cm).evaluate(x, y), c, tol, 0);
        });
    });

    it('reports derivatives consistently with central differences', () => {
        const h = 1e-5;
        check(fc.tuple(samples(), geometry(), fc.boolean(),
            interiorParam(XB), interiorParam(YB)),
            ([F, gg, cm, [ix, u], [iy, v]]) => {
                const g = gg as Geom;
                const intp = build(F, gg as Geom, cm);
                const x = g[0] + g[1] * (ix + u);
                const y = g[2] + g[3] * (iy + v);
                const dfx = (intp.evaluate(x + h, y) - intp.evaluate(x - h, y)) / (2 * h);
                const dfy = (intp.evaluate(x, y + h) - intp.evaluate(x, y - h)) / (2 * h);
                // The offsets stay inside the same cell (u is at least 0.02
                // cell widths from a boundary and the spacing is at least
                // 0.25). Truncation is h^2/6 * P''' with P''' up to ~1e3;
                // roundoff is O(eps * |P| / h) ~ 1e-11.
                expectClose(intp.evaluate(1, 0, x, y), dfx, 1e-4, 1e-5);
                expectClose(intp.evaluate(0, 1, x, y), dfy, 1e-4, 1e-5);
            });
    });

    it('is affine in the sample values', () => {
        // Unlike the Akima schemes, the blend weights do not depend on the
        // data at all, so the interpolant is exactly linear in the samples.
        check(fc.tuple(samples(), samples(), wellScaled(-4, 4), geometry(),
            fc.boolean(), scaled(0, 1), scaled(0, 1)),
            ([P, Q, k, gg, cm, u, v]) => {
                const g = gg as Geom;
                const x = g[0] + u * g[1] * (XB - 1);
                const y = g[2] + v * g[3] * (YB - 1);
                const lhs = build(P.map((p, i) => p + k * Q[i]), g, cm).evaluate(x, y);
                const rhs = build(P, g, cm).evaluate(x, y)
                    + k * build(Q, g, cm).evaluate(x, y);
                expectClose(lhs, rhs, 1e-10, 1e-12);
            });
    });

    it('is equivariant under translation of the grid origin', () => {
        check(fc.tuple(samples(), geometry(), finite(-6, 6), finite(-6, 6),
            fc.boolean(), scaled(0, 1), scaled(0, 1)),
            ([F, gg, tx, ty, cm, u, v]) => {
                const g = gg as Geom;
                const a = build(F, g, cm);
                const b = build(F, [g[0] + tx, g[1], g[2] + ty, g[3]], cm);
                const x = g[0] + u * g[1] * (XB - 1);
                const y = g[2] + v * g[3] * (YB - 1);
                // (x + tx) - (xMin + tx) reassociates the shift, so the cell
                // parameter differs by a few ulps; the blend then amplifies
                // that by the derivative of the interpolant.
                expectClose(a.evaluate(x, y), b.evaluate(x + tx, y + ty),
                    1e-9, 1e-10);
            });
    });

    it('returns zero for derivative orders above three', () => {
        check(fc.tuple(samples(), geometry(), fc.boolean(),
            fc.integer({ min: 4, max: 20 }), scaled(0, 1), scaled(0, 1)),
            ([F, gg, cm, big, u, v]) => {
                const g = gg as Geom;
                const intp = build(F, g, cm);
                const x = g[0] + u * g[1] * (XB - 1);
                const y = g[2] + v * g[3] * (YB - 1);
                expect(intp.evaluate(big, 0, x, y)).toBe(0);
                expect(intp.evaluate(0, big, x, y)).toBe(0);
            });
    });

    it('rejects grids smaller than 3x3 and non-positive spacings', () => {
        const F = new Array<number>(9).fill(1);
        expect(() => new IntpBicubic2(2, 3, 0, 1, 0, 1, F, true)).toThrow();
        expect(() => new IntpBicubic2(3, 2, 0, 1, 0, 1, F, true)).toThrow();
        expect(() => new IntpBicubic2(3, 3, 0, 0, 0, 1, F, true)).toThrow();
        expect(() => new IntpBicubic2(3, 3, 0, 1, 0, -1, F, false)).toThrow();
        expect(() => new IntpBicubic2(3, 3, 0, 1, 0, 1, [1, 2, 3], true)).toThrow();
    });
});
