import { describe, it, expect } from 'vitest';
import { IntpAkimaUniform2 } from '../src/IntpAkimaUniform2.js';
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

describe('IntpAkimaUniform2', () => {
    it('throws for invalid inputs', () => {
        const F9 = new Array<number>(9).fill(0);
        expect(() => new IntpAkimaUniform2(2, 3, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 2, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 3, 0, 0, 0, 1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(3, 3, 0, 1, 0, -1, F9)).toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform2(4, 3, 0, 1, 0, 1, F9)).toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const F = makeSamples(4, 5, -1, 0.5, 2, 0.25, (x, y) => x + y);
        const interp = new IntpAkimaUniform2(4, 5, -1, 0.5, 2, 0.25, F);
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

    it('interpolates a hand-computed 3x3 grid of f(x,y) = x*y', () => {
        // For f = x*y all x-slopes in a row equal y (and y-slopes in a
        // column equal x), so the Akima estimates fx = y and fy = x are
        // exact everywhere, and the interior/min-boundary fxy estimates are
        // exactly 1. In the cell [0,1]^2 the polynomial is exactly x*y.
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x * y);
        const interp = new IntpAkimaUniform2(3, 3, 0, 1, 0, 1, F);

        for (const [x, y] of [[0.5, 0.5], [0.25, 0.75], [0, 0.3], [1, 1]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(x * y, 12);
            expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(x * y, 12);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(y, 12);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(x, 12);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(1, 12);
            expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(0, 12);
            expect(interp.evaluate(0, 2, x, y)).toBeCloseTo(0, 12);
        }

        // Regression pin for the corrected max-boundary sign in getFXY (see
        // the note there and upstream issue #58). Upstream reuses the
        // min-boundary one-sided difference coefficients at the max
        // boundaries without negating them for the reversed direction, which
        // makes the fxy estimates on the top row and right column come out as
        // -1 instead of +1 for f = x*y and leaves the interpolant equal to
        // 1.1015625 at (1.5, 0.75) instead of the exact 1.125. With the sign
        // corrected the bilinear function is reproduced over the whole
        // rectangle, including the cells that touch the max boundaries.
        expect(interp.evaluate(1.5, 0.75)).toBeCloseTo(1.5 * 0.75, 12);
        for (const [x, y] of [[1.5, 0.75], [1.75, 1.5], [2, 2], [0.5, 1.9]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(x * y, 12);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(1, 12);
        }
        // The samples themselves are still interpolated exactly.
        for (let r = 0; r < 3; ++r) {
            for (let c = 0; c < 3; ++c) {
                expect(interp.evaluate(c, r)).toBeCloseTo(c * r, 12);
            }
        }
    });

    it('passes through the samples of a non-polynomial function', () => {
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.5 * y);
        const xBound = 6, yBound = 5;
        const xMin = -1, dx = 0.5, yMin = 0, dy = 0.4;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (let r = 0; r < yBound; ++r) {
            for (let c = 0; c < xBound; ++c) {
                const x = xMin + c * dx;
                const y = yMin + r * dy;
                expect(interp.evaluate(x, y)).toBeCloseTo(F[c + xBound * r], 10);
                expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(F[c + xBound * r], 10);
            }
        }
    });

    it('exactly reproduces a per-variable quadratic away from the max boundaries', () => {
        // f(x,y) = (1 + 2x - x^2) * (3 - y + 0.5 y^2). Uniform slopes of a
        // quadratic are linear in the index, so the Akima weighted average
        // reduces to the exact central difference and the boundary
        // extrapolations are exact as well; the interpolant reproduces f in
        // every cell whose corners avoid the max-boundary fxy estimates
        // (see the upstream GetFXY sign note in the previous test).
        const g = (x: number): number => 1 + 2 * x - x * x;
        const dg = (x: number): number => 2 - 2 * x;
        const ddg = (): number => -2;
        const h = (y: number): number => 3 - y + 0.5 * y * y;
        const dh = (y: number): number => -1 + y;
        const ddh = (): number => 1;

        const xBound = 5, yBound = 6;
        const xMin = -0.5, dx = 0.5, yMin = 1, dy = 0.25;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy,
            (x, y) => g(x) * h(y));
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (const [x, y] of [[-0.3, 1.1], [0.42, 1.87], [0.9, 1.6], [0.75, 1.95]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(g(x) * h(y), 9);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(dg(x) * h(y), 9);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(g(x) * dh(y), 9);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(dg(x) * dh(y), 9);
            expect(interp.evaluate(2, 0, x, y)).toBeCloseTo(ddg() * h(y), 9);
            expect(interp.evaluate(0, 2, x, y)).toBeCloseTo(g(x) * ddh(), 9);
            expect(interp.evaluate(2, 2, x, y)).toBeCloseTo(ddg() * ddh(), 8);
        }
    });

    it('is C1 across interior cell boundaries', () => {
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.5 * y);
        const xBound = 6, yBound = 6;
        const xMin = 0, dx = 0.5, yMin = 0, dy = 0.5;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, F);

        const eps = 1e-7;
        // Cross the vertical grid line x = 1.0 and the horizontal grid line
        // y = 1.5 at several positions; value and first derivatives must
        // agree in the limit from both sides.
        for (const y of [0.3, 1.2, 2.1]) {
            for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                const left = interp.evaluate(xo, yo, 1 - eps, y);
                const right = interp.evaluate(xo, yo, 1 + eps, y);
                expect(right - left).toBeCloseTo(0, 5);
            }
        }
        for (const x of [0.4, 1.3, 2.2]) {
            for (const [xo, yo] of [[0, 0], [1, 0], [0, 1]]) {
                const below = interp.evaluate(xo, yo, x, 1.5 - eps);
                const above = interp.evaluate(xo, yo, x, 1.5 + eps);
                expect(above - below).toBeCloseTo(0, 5);
            }
        }
    });

    it('clamps evaluations to the domain', () => {
        const F = makeSamples(4, 4, 0, 1, 0, 1, (x, y) => x * x + y);
        const interp = new IntpAkimaUniform2(4, 4, 0, 1, 0, 1, F);
        expect(interp.evaluate(-10, 1.5)).toBeCloseTo(interp.evaluate(0, 1.5), 12);
        expect(interp.evaluate(10, 1.5)).toBeCloseTo(interp.evaluate(3, 1.5), 12);
        expect(interp.evaluate(1.5, -10)).toBeCloseTo(interp.evaluate(1.5, 0), 12);
        expect(interp.evaluate(1.5, 10)).toBeCloseTo(interp.evaluate(1.5, 3), 12);
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x + y);
        const interp = new IntpAkimaUniform2(3, 3, 0, 1, 0, 1, F);
        expect(interp.evaluate(4, 0, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(0, 4, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(3, 0, 0.5, 0.5)).toBeCloseTo(0, 12);
    });
});

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against IntpAkimaUniform2.h.
// ---------------------------------------------------------------------------
describe('IntpAkimaUniform2 verification', () => {
    const XB = 5;
    const YB = 4;
    // Grid samples in row-major order: F[ix + XB*iy] is f(x_ix, y_iy).
    // wellScaled snaps |x| < 1e-3 to zero so no sample is subnormal; the
    // construction divides sample differences by the spacing twice, and a
    // subnormal there would make relative tolerances meaningless.
    const samples = () => fc.array(wellScaled(-5, 5),
        { minLength: XB * YB, maxLength: XB * YB });
    // The spacing is bounded away from zero because construct() divides the
    // corner data by dx^2 dy^2; a 1e-3 spacing would amplify roundoff by 1e12.
    const latticeSamples = () => fc.array(fc.integer({ min: -20, max: 20 }),
        { minLength: XB * YB, maxLength: XB * YB });
    const signedPowerOfTwo = () =>
        fc.tuple(fc.integer({ min: -4, max: 4 }), fc.boolean())
            .map(([m, neg]) => (neg ? -1 : 1) * Math.pow(2, m));
    const geometry = () => fc.tuple(finite(-3, 3), positive(3, 0.5),
        finite(-3, 3), positive(3, 0.5));

    const build = (F: number[], g: [number, number, number, number]) =>
        new IntpAkimaUniform2(XB, YB, g[0], g[1], g[2], g[3], F);

    it('interpolates the samples at the grid points', () => {
        check(fc.tuple(samples(), geometry()), ([F, g]) => {
            const intp = build(F, g);
            for (let iy = 0; iy < YB; ++iy) {
                for (let ix = 0; ix < XB; ++ix) {
                    const x = g[0] + g[1] * ix;
                    const y = g[2] + g[3] * iy;
                    // The polynomial constant term is exactly F at the cell
                    // corner, so only the lookup's dx = x - (xMin + i*dx)
                    // roundoff enters, and it is multiplied by O(1) slopes.
                    expectClose(intp.evaluate(x, y), F[ix + XB * iy], 1e-10, 1e-11);
                }
            }
        });
    });

    it('reproduces a + b x + c y + d x y exactly', () => {
        // Every x-slope in a row is then b + d*y, so all four window slopes
        // agree and computeDerivative returns them unchanged; the mixed
        // difference for FXY is exact for the xy term. The bicubic through
        // exact F, FX, FY, FXY at the corners is therefore the data itself.
        check(fc.tuple(fc.array(finite(-2, 2), { minLength: 4, maxLength: 4 }),
            geometry(), scaled(0, 1), scaled(0, 1)), ([abcd, g, u, v]) => {
            const [a, b, c, d] = abcd;
            const f = (x: number, y: number) => a + b * x + c * y + d * x * y;
            const F: number[] = [];
            for (let iy = 0; iy < YB; ++iy) {
                for (let ix = 0; ix < XB; ++ix) {
                    F.push(f(g[0] + g[1] * ix, g[2] + g[3] * iy));
                }
            }
            const intp = build(F, g);
            const x = g[0] + u * g[1] * (XB - 1);
            const y = g[2] + v * g[3] * (YB - 1);
            // The cubic correction terms are differences of quantities that
            // cancel exactly in exact arithmetic, divided by dx^2 dy^2, so the
            // absolute error is about |f| * eps / spacing^2 -- with |f| <= 20
            // and spacing >= 0.5 that is a few times 1e-14 per term, and the
            // mixed second difference accumulates it over 16 terms and two
            // more divisions by the spacing.
            expectClose(intp.evaluate(x, y), f(x, y), 1e-8, 1e-10);
            expectClose(intp.evaluate(1, 0, x, y), b + d * y, 1e-7, 1e-9);
            expectClose(intp.evaluate(0, 1, x, y), c + d * x, 1e-7, 1e-9);
            expectClose(intp.evaluate(1, 1, x, y), d, 1e-7, 1e-9);
        });
    });

    it('clamps queries to the sample rectangle', () => {
        check(fc.tuple(samples(), geometry(), finite(-30, 30), finite(-30, 30)),
            ([F, g, x, y]) => {
                const intp = build(F, g);
                const cx = Math.min(Math.max(x, intp.getXMin()), intp.getXMax());
                const cy = Math.min(Math.max(y, intp.getYMin()), intp.getYMax());
                expect(intp.evaluate(x, y)).toBe(intp.evaluate(cx, cy));
                expect(intp.evaluate(1, 1, x, y)).toBe(intp.evaluate(1, 1, cx, cy));
            });
    });

    it('is homogeneous of degree one in the sample values', () => {
        // computeDerivative is homogeneous but not additive (the branch and
        // the |.| weights depend nonlinearly on the slopes), so scaling -- not
        // superposition -- is the invariant the Akima estimator satisfies.
    // Integer samples (and an integer offset, or a power-of-two scale) keep
    // every slope difference exact, so the branch that computeDerivative takes
    // is unchanged by the transformation. With general doubles the Akima
    // weights |s3-s2| and |s0-s1| are catastrophically cancelled differences
    // and the convex weight they define is not stable under either map.
        check(fc.tuple(latticeSamples(), signedPowerOfTwo(), geometry(),
            scaled(0, 1), scaled(0, 1)), ([F, k, g, u, v]) => {
            const x = g[0] + u * g[1] * (XB - 1);
            const y = g[2] + v * g[3] * (YB - 1);
            const lhs = build(F.map(f => k * f), g).evaluate(x, y);
            const rhs = k * build(F, g).evaluate(x, y);
            expectClose(lhs, rhs, 1e-9, 1e-11);
        });
    });

    it('adds a constant to the result when a constant is added to the samples', () => {
        // The derivative estimates are differences and so are unchanged; only
        // the polynomial's constant term moves.
        check(fc.tuple(latticeSamples(), fc.integer({ min: -30, max: 30 }),
            geometry(), scaled(0, 1), scaled(0, 1)), ([F, c, g, u, v]) => {
            const x = g[0] + u * g[1] * (XB - 1);
            const y = g[2] + v * g[3] * (YB - 1);
            const lhs = build(F.map(f => f + c), g).evaluate(x, y);
            const rhs = build(F, g).evaluate(x, y) + c;
            expectClose(lhs, rhs, 1e-9, 1e-11);
            // Derivatives are untouched by the shift.
            expectClose(build(F.map(f => f + c), g).evaluate(1, 0, x, y),
                build(F, g).evaluate(1, 0, x, y), 1e-9, 1e-11);
        });
    });

    it('is equivariant under translation of the grid origin', () => {
        check(fc.tuple(samples(), geometry(), finite(-6, 6), finite(-6, 6),
            scaled(0, 1), scaled(0, 1)), ([F, g, tx, ty, u, v]) => {
            const a = build(F, g);
            const b = build(F, [g[0] + tx, g[1], g[2] + ty, g[3]]);
            const x = g[0] + u * g[1] * (XB - 1);
            const y = g[2] + v * g[3] * (YB - 1);
            // The lookup subtracts the origin before evaluating, so the only
            // difference is the roundoff of (x + tx) - (xMin + tx).
            expectClose(a.evaluate(x, y), b.evaluate(x + tx, y + ty), 1e-9, 1e-10);
        });
    });

    it('ignores y when the samples depend on x only', () => {
        check(fc.tuple(fc.array(finite(-5, 5), { minLength: XB, maxLength: XB }),
            geometry(), scaled(0, 1), scaled(0, 1), scaled(0, 1)),
            ([row, g, u, v, w]) => {
                const F: number[] = [];
                for (let iy = 0; iy < YB; ++iy) { F.push(...row); }
                const intp = build(F, g);
                const x = g[0] + u * g[1] * (XB - 1);
                const y0 = g[2] + v * g[3] * (YB - 1);
                const y1 = g[2] + w * g[3] * (YB - 1);
                expectClose(intp.evaluate(x, y0), intp.evaluate(x, y1), 1e-10, 1e-11);
                expectClose(intp.evaluate(0, 1, x, y0), 0, 1e-9, 0);
            });
    });

    it('reports first derivatives consistently with central differences', () => {
        const h = 1e-5;
        check(fc.tuple(samples(), geometry(), scaled(0.2, 0.8), scaled(0.2, 0.8)),
            ([F, g, u, v]) => {
                const intp = build(F, g);
                // Stay inside one cell so both offsets use the same polynomial.
                const ix = 1 + Math.floor(u * (XB - 3));
                const iy = 1 + Math.floor(v * (YB - 3));
                const x = g[0] + g[1] * (ix + 0.5);
                const y = g[2] + g[3] * (iy + 0.5);
                const dfx = (intp.evaluate(x + h, y) - intp.evaluate(x - h, y)) / (2 * h);
                const dfy = (intp.evaluate(x, y + h) - intp.evaluate(x, y - h)) / (2 * h);
                // Truncation h^2/6 * P''' with P''' up to ~1e3 for spacing
                // 0.25 and |F| <= 5; roundoff is O(eps * |P| / h) ~ 1e-10.
                expectClose(intp.evaluate(1, 0, x, y), dfx, 1e-4, 1e-5);
                expectClose(intp.evaluate(0, 1, x, y), dfy, 1e-4, 1e-5);
            });
    });

    it('returns zero for derivative orders above three', () => {
        check(fc.tuple(samples(), geometry(), fc.integer({ min: 4, max: 20 }),
            scaled(0, 1)), ([F, g, big, u]) => {
            const intp = build(F, g);
            const x = g[0] + u * g[1] * (XB - 1);
            expect(intp.evaluate(big, 0, x, g[2])).toBe(0);
            expect(intp.evaluate(0, big, x, g[2])).toBe(0);
        });
    });

    it('rejects grids smaller than 3x3 and non-positive spacings', () => {
        const F = new Array<number>(9).fill(1);
        expect(() => new IntpAkimaUniform2(2, 3, 0, 1, 0, 1, F)).toThrow();
        expect(() => new IntpAkimaUniform2(3, 2, 0, 1, 0, 1, F)).toThrow();
        expect(() => new IntpAkimaUniform2(3, 3, 0, 0, 0, 1, F)).toThrow();
        expect(() => new IntpAkimaUniform2(3, 3, 0, 1, 0, -1, F)).toThrow();
        expect(() => new IntpAkimaUniform2(3, 3, 0, 1, 0, 1, [1, 2, 3])).toThrow();
    });
});
