import { describe, it, expect } from 'vitest';
import { IntpBilinear2 } from '../src/IntpBilinear2.js';

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
