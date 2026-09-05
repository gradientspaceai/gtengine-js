import { describe, it, expect } from 'vitest';
import { IntpBilinear2 } from '../src/IntpBilinear2.js';
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

describe('IntpBilinear2', () => {
    it('throws for invalid inputs', () => {
        const F4 = new Array<number>(4).fill(0);
        expect(() => new IntpBilinear2(1, 2, 0, 1, 0, 1, F4)).toThrow('Invalid input.');
        expect(() => new IntpBilinear2(2, 1, 0, 1, 0, 1, F4)).toThrow('Invalid input.');
        expect(() => new IntpBilinear2(2, 2, 0, 0, 0, 1, F4)).toThrow('Invalid input.');
        expect(() => new IntpBilinear2(2, 2, 0, 1, 0, -1, F4)).toThrow('Invalid input.');
        expect(() => new IntpBilinear2(3, 2, 0, 1, 0, 1, F4)).toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const F = makeSamples(4, 5, -1, 0.5, 2, 0.25, (x, y) => x + y);
        const interp = new IntpBilinear2(4, 5, -1, 0.5, 2, 0.25, F);
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

    it('reproduces the samples at the grid points', () => {
        const xBound = 5, yBound = 4;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25;
        const f = (x: number, y: number): number => Math.sin(x) * Math.exp(0.3 * y);
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, f);
        const interp = new IntpBilinear2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (let r = 0; r < yBound; ++r) {
            for (let c = 0; c < xBound; ++c) {
                const value = interp.evaluate(xMin + c * dx, yMin + r * dy);
                expect(value).toBeCloseTo(F[c + xBound * r], 12);
            }
        }
    });

    it('is exact for a bilinear function and its derivatives', () => {
        // p(x,y) = 3 - 2x + 0.5y + 4xy is bilinear, so the interpolant
        // reproduces it exactly inside the domain.
        const p = (x: number, y: number): number => 3 - 2 * x + 0.5 * y + 4 * x * y;
        const px = (_x: number, y: number): number => -2 + 4 * y;
        const py = (x: number, _y: number): number => 0.5 + 4 * x;

        const xBound = 4, yBound = 5;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25;
        const F = makeSamples(xBound, yBound, xMin, dx, yMin, dy, p);
        const interp = new IntpBilinear2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (const [x, y] of [[-0.9, 2.1], [-0.25, 2.4], [0.1, 2.75], [0.35, 2.9]]) {
            expect(interp.evaluate(x, y)).toBeCloseTo(p(x, y), 12);
            expect(interp.evaluate(0, 0, x, y)).toBeCloseTo(p(x, y), 12);
            expect(interp.evaluate(1, 0, x, y)).toBeCloseTo(px(x, y), 12);
            expect(interp.evaluate(0, 1, x, y)).toBeCloseTo(py(x, y), 12);
            expect(interp.evaluate(1, 1, x, y)).toBeCloseTo(4, 12);
        }
    });

    it('matches a hand-computed bilinear blend inside a cell', () => {
        // The 2x2 grid values are f00 = 1, f10 = 2, f01 = 4, f11 = 8 on the
        // unit square; at (0.25, 0.5) the blend is
        // (1-u)(1-v)f00 + u(1-v)f10 + (1-u)v f01 + u v f11 with u = 0.25,
        // v = 0.5, which equals 0.375*1 + 0.125*2 + 0.375*4 + 0.125*8 = 3.125.
        const F = [1, 2, 4, 8];
        const interp = new IntpBilinear2(2, 2, 0, 1, 0, 1, F);
        expect(interp.evaluate(0.25, 0.5)).toBeCloseTo(3.125, 12);
        // d/dx = (1-v)(f10-f00) + v(f11-f01) = 0.5*1 + 0.5*4 = 2.5.
        expect(interp.evaluate(1, 0, 0.25, 0.5)).toBeCloseTo(2.5, 12);
        // d/dy = (1-u)(f01-f00) + u(f11-f10) = 0.75*3 + 0.25*6 = 3.75.
        expect(interp.evaluate(0, 1, 0.25, 0.5)).toBeCloseTo(3.75, 12);
        // d2/dxdy = f00 - f10 - f01 + f11 = 1 - 2 - 4 + 8 = 3.
        expect(interp.evaluate(1, 1, 0.25, 0.5)).toBeCloseTo(3, 12);
    });

    it('scales the derivatives by the sample spacing', () => {
        // The same unit-square data on a grid with spacing (2, 4) gives
        // derivatives divided by the spacings.
        const F = [1, 2, 4, 8];
        const interp = new IntpBilinear2(2, 2, 0, 2, 0, 4, F);
        expect(interp.evaluate(0.5, 2)).toBeCloseTo(3.125, 12);
        expect(interp.evaluate(1, 0, 0.5, 2)).toBeCloseTo(2.5 / 2, 12);
        expect(interp.evaluate(0, 1, 0.5, 2)).toBeCloseTo(3.75 / 4, 12);
        expect(interp.evaluate(1, 1, 0.5, 2)).toBeCloseTo(3 / (2 * 4), 12);
    });

    it('clamps the sample indices, not the inputs, outside the domain', () => {
        // Upstream's comment says the inputs are clamped to the domain, but
        // the code clamps only the sample indices; the fractional parts
        // 'xIndex - ix' and 'yIndex - iy' are left unclamped. Above the max
        // boundary the clamped 2x2 block has two identical columns/rows and
        // the blend weights sum to one, so the value is constant. Below the
        // min boundary the fractional part is negative and the interpolant
        // extrapolates linearly off the first cell. The port preserves this.
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x + 2 * y);
        const interp = new IntpBilinear2(3, 3, 0, 1, 0, 1, F);
        expect(interp.evaluate(2, 1)).toBeCloseTo(4, 12);
        expect(interp.evaluate(5, 1)).toBeCloseTo(4, 12);
        expect(interp.evaluate(1, 7)).toBeCloseTo(5, 12);
        // Linear extrapolation below the min boundary.
        expect(interp.evaluate(-3, 1)).toBeCloseTo(-1, 12);
        expect(interp.evaluate(1, -7)).toBeCloseTo(-13, 12);
        // The x-derivative vanishes on the max boundary because the clamped
        // 2x2 block has two identical columns.
        expect(interp.evaluate(1, 0, 2, 1)).toBeCloseTo(0, 12);
        // It is still exact below the min boundary.
        expect(interp.evaluate(1, 0, -3, 1)).toBeCloseTo(1, 12);
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const F = makeSamples(3, 3, 0, 1, 0, 1, (x, y) => x * y);
        const interp = new IntpBilinear2(3, 3, 0, 1, 0, 1, F);
        expect(interp.evaluate(2, 0, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(0, 2, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(-1, 0, 0.5, 0.5)).toBe(0);
    });

    it('agrees with an independent bilinear evaluation on random inputs', () => {
        const xBound = 7, yBound = 6;
        const xMin = -2, dx = 0.75, yMin = 1, dy = 0.4;
        const F: number[] = [];
        let seed = 12345;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let i = 0; i < xBound * yBound; ++i) {
            F.push(2 * rand() - 1);
        }
        const interp = new IntpBilinear2(xBound, yBound, xMin, dx, yMin, dy, F);

        for (let trial = 0; trial < 200; ++trial) {
            const x = xMin + rand() * dx * (xBound - 1);
            const y = yMin + rand() * dy * (yBound - 1);
            const u = (x - xMin) / dx;
            const v = (y - yMin) / dy;
            const c = Math.min(Math.floor(u), xBound - 2);
            const r = Math.min(Math.floor(v), yBound - 2);
            const s = u - c;
            const t = v - r;
            const expected =
                (1 - s) * (1 - t) * F[c + xBound * r] +
                s * (1 - t) * F[c + 1 + xBound * r] +
                (1 - s) * t * F[c + xBound * (r + 1)] +
                s * t * F[c + 1 + xBound * (r + 1)];
            expect(interp.evaluate(x, y)).toBeCloseTo(expected, 12);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V28): property-based cross-checks against IntpBilinear2.h.
// ---------------------------------------------------------------------------
describe('IntpBilinear2 verification', () => {
    const XB = 5, YB = 4;
    const geometry = () => fc.tuple(finite(-3, 3), positive(3, 0.25),
        finite(-3, 3), positive(3, 0.25));
    type Geom = [number, number, number, number];
    const samples = () => fc.array(wellScaled(-6, 6),
        { minLength: XB * YB, maxLength: XB * YB });
    const build = (F: number[], g: Geom) =>
        new IntpBilinear2(XB, YB, g[0], g[1], g[2], g[3], F);

    // The closed-form bilinear blend with upstream's index clamping, written
    // independently of the blend-matrix formulation the port uses.
    const reference = (F: number[], g: Geom, x: number, y: number): number => {
        const s = (x - g[0]) / g[1];
        const t = (y - g[2]) / g[3];
        let ix = Math.trunc(s);
        ix = ix < 0 ? 0 : (ix >= XB ? XB - 1 : ix);
        let iy = Math.trunc(t);
        iy = iy < 0 ? 0 : (iy >= YB ? YB - 1 : iy);
        const u = s - ix, v = t - iy;
        const cx = (i: number) => (i >= XB ? XB - 1 : i);
        const cy = (j: number) => (j >= YB ? YB - 1 : j);
        const f = (i: number, j: number) => F[cx(i) + XB * cy(j)];
        return (1 - u) * (1 - v) * f(ix, iy) + u * (1 - v) * f(ix + 1, iy)
            + (1 - u) * v * f(ix, iy + 1) + u * v * f(ix + 1, iy + 1);
    };

    it('agrees with the closed-form bilinear blend', () => {
        check(fc.tuple(samples(), geometry(), scaled(0, 1), scaled(0, 1)),
            ([F, gg, u, v]) => {
                const g = gg as Geom;
                const x = g[0] + u * g[1] * (XB - 1);
                const y = g[2] + v * g[3] * (YB - 1);
                // The port factors the blend as (M*U)^t D (M*V) while the
                // reference multiplies the four weights out directly, so the
                // two differ only by the order of a handful of roundings.
                expectClose(build(F, g).evaluate(x, y), reference(F, g, x, y),
                    1e-12, 1e-13);
            });
    });

    it('interpolates the samples at the grid points', () => {
        check(fc.tuple(samples(), geometry()), ([F, gg]) => {
            const g = gg as Geom;
            const intp = build(F, g);
            for (let iy = 0; iy < YB; ++iy) {
                for (let ix = 0; ix < XB; ++ix) {
                    // At a node the fractional parts are exactly 0 (or, at the
                    // max boundary, the stencil collapses onto the last
                    // sample), so the blend returns the sample itself.
                    expectClose(intp.evaluate(g[0] + g[1] * ix, g[2] + g[3] * iy),
                        F[ix + XB * iy], 1e-11, 1e-12);
                }
            }
        });
    });

    it('reproduces a + b x + c y + d x y inside the domain', () => {
        check(fc.tuple(fc.array(wellScaled(-2, 2), { minLength: 4, maxLength: 4 }),
            geometry(), scaled(0, 1), scaled(0, 1), scaled(0, 0.9),
            scaled(0, 0.9)), ([k, gg, u, v, du, dv]) => {
            const g = gg as Geom;
            const f = (x: number, y: number) =>
                k[0] + k[1] * x + k[2] * y + k[3] * x * y;
            const F: number[] = [];
            for (let iy = 0; iy < YB; ++iy) {
                for (let ix = 0; ix < XB; ++ix) {
                    F.push(f(g[0] + g[1] * ix, g[2] + g[3] * iy));
                }
            }
            const intp = build(F, g);
            expectClose(intp.evaluate(g[0] + u * g[1] * (XB - 1),
                g[2] + v * g[3] * (YB - 1)),
                f(g[0] + u * g[1] * (XB - 1), g[2] + v * g[3] * (YB - 1)),
                1e-11, 1e-12);
            // The derivatives are sampled strictly inside the domain: at the
            // max boundary the 2x2 stencil collapses onto the last row and
            // column (upstream clamps ix + 1 and iy + 1 but not the fractional
            // coordinate), so the reported slope there is zero, not the true
            // slope. That degeneracy is checked separately below.
            const x = g[0] + du * g[1] * (XB - 1);
            const y = g[2] + dv * g[3] * (YB - 1);
            expectClose(intp.evaluate(1, 0, x, y), k[1] + k[3] * y, 1e-10, 1e-11);
            expectClose(intp.evaluate(0, 1, x, y), k[2] + k[3] * x, 1e-10, 1e-11);
        });
    });

    it('reports a zero slope at the max boundary, where the stencil collapses', () => {
        // With an exactly representable geometry, xIndex at x = xMax is
        // exactly xBound - 1, so ix stays at xBound - 1 with a zero fractional
        // part and the ix + 1 sample is clamped back onto ix. The difference
        // that forms the x-derivative therefore vanishes. This is upstream
        // behaviour, a direct consequence of clamping the index rather than
        // the input; the geometry is fixed here because (xMax - xMin) *
        // (1/xSpacing) lands just below xBound - 1 for most spacings, in which
        // case the last real cell is used and the slope is the true one.
        const F: number[] = [];
        for (let iy = 0; iy < YB; ++iy) {
            for (let ix = 0; ix < XB; ++ix) { F.push(ix + 10 * iy); }
        }
        const intp = new IntpBilinear2(XB, YB, 0, 1, 0, 1, F);
        expect(intp.evaluate(1, 0, XB - 1, 0) + 0).toBe(0);
        expect(intp.evaluate(0, 1, 0, YB - 1) + 0).toBe(0);
        // Just inside, the slopes are the true ones (1 in x, 10 in y).
        expect(intp.evaluate(1, 0, XB - 1.5, 0)).toBeCloseTo(1, 12);
        expect(intp.evaluate(0, 1, 0, YB - 1.5)).toBeCloseTo(10, 12);
        // The value at the far corner is still the corner sample.
        expect(intp.evaluate(XB - 1, YB - 1)).toBeCloseTo(F[XB * YB - 1], 12);
    });

    it('reproduces constants everywhere, including outside the domain', () => {
        // The blend rows sum to one for every parameter value, so a constant
        // field is reproduced even where the clamped stencil degenerates.
        check(fc.tuple(finite(-5, 5), geometry(), finite(-30, 30), finite(-30, 30)),
            ([c, gg, x, y]) => {
                const g = gg as Geom;
                const F = new Array<number>(XB * YB).fill(c);
                expectClose(build(F, g).evaluate(x, y), c, 1e-9, 1e-10);
            });
    });

    it('gives the cell slope as the first derivative', () => {
        check(fc.tuple(samples(), geometry(), fc.integer({ min: 0, max: XB - 2 }),
            fc.integer({ min: 0, max: YB - 2 }), scaled(0.05, 0.95),
            scaled(0.05, 0.95)), ([F, gg, ix, iy, u, v]) => {
            const g = gg as Geom;
            const x = g[0] + g[1] * (ix + u);
            const y = g[2] + g[3] * (iy + v);
            const at = (i: number, j: number) => F[i + XB * j];
            // d/dx of the bilinear blend is the y-interpolated x-difference.
            const dfx = ((1 - v) * (at(ix + 1, iy) - at(ix, iy))
                + v * (at(ix + 1, iy + 1) - at(ix, iy + 1))) / g[1];
            const dfy = ((1 - u) * (at(ix, iy + 1) - at(ix, iy))
                + u * (at(ix + 1, iy + 1) - at(ix + 1, iy))) / g[3];
            const intp = build(F, g);
            // The reconstruction of u, v from x, y costs a couple of ulps,
            // amplified by 1/spacing (bounded by 4).
            expectClose(intp.evaluate(1, 0, x, y), dfx, 1e-10, 1e-11);
            expectClose(intp.evaluate(0, 1, x, y), dfy, 1e-10, 1e-11);
        });
    });

    it('returns zero for derivative orders above one', () => {
        check(fc.tuple(samples(), geometry(), fc.integer({ min: 2, max: 20 }),
            scaled(0, 1), scaled(0, 1)), ([F, gg, big, u, v]) => {
            const g = gg as Geom;
            const intp = build(F, g);
            const x = g[0] + u * g[1] * (XB - 1);
            const y = g[2] + v * g[3] * (YB - 1);
            expect(intp.evaluate(big, 0, x, y)).toBe(0);
            expect(intp.evaluate(0, big, x, y)).toBe(0);
        });
    });

    it('extrapolates rather than clamping the query, as upstream does', () => {
        // Upstream's comment promises clamping of x and y, but only the cell
        // index is clamped; the fractional coordinate keeps its out-of-range
        // value. Below the minimum that extrapolates the first cell linearly.
        check(fc.tuple(samples(), geometry(), scaled(0.05, 0.95)),
            ([F, gg, s]) => {
                const g = gg as Geom;
                const intp = build(F, g);
                const y = g[2];
                const x = g[0] - s * g[1];
                const f0 = F[0], f1 = F[1];
                // xIndex in (-1, 0) truncates to ix = 0 and u = xIndex < 0, so
                // the blend evaluates the first cell's line at a negative
                // parameter.
                expectClose(intp.evaluate(x, y), f0 + (f1 - f0) * (-s),
                    1e-10, 1e-11);

            });
    });

    it('does not clamp the query: a fixed case below the minimum', () => {
        // A genuine clamp would return F[0] = 1; upstream extrapolates the
        // first cell's line to the parameter -0.5, giving 1 - 0.5*(3 - 1) = 0.
        const F = [1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7, 1, 3, 5, 7];
        const intp = new IntpBilinear2(XB, YB, 0, 1, 0, 1, F);
        expect(intp.evaluate(-0.5, 0)).toBeCloseTo(0, 12);
        expect(intp.evaluate(0, 0)).toBeCloseTo(1, 12);
    });

    it('rejects grids smaller than 2x2 and non-positive spacings', () => {
        const F = new Array<number>(4).fill(1);
        expect(() => new IntpBilinear2(1, 2, 0, 1, 0, 1, F)).toThrow();
        expect(() => new IntpBilinear2(2, 1, 0, 1, 0, 1, F)).toThrow();
        expect(() => new IntpBilinear2(2, 2, 0, 0, 0, 1, F)).toThrow();
        expect(() => new IntpBilinear2(2, 2, 0, 1, 0, -1, F)).toThrow();
        expect(() => new IntpBilinear2(2, 2, 0, 1, 0, 1, [1, 2])).toThrow();
    });
});
