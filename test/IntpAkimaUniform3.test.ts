import { describe, it, expect } from 'vitest';
import { IntpAkimaUniform3 } from '../src/IntpAkimaUniform3.js';
import { IntpAkimaUniform2 } from '../src/IntpAkimaUniform2.js';

// Build the lexicographic sample array F[c + xBound*(r + yBound*s)] =
// f(xMin + c*dx, yMin + r*dy, zMin + s*dz).
function makeSamples(xBound: number, yBound: number, zBound: number,
    xMin: number, dx: number, yMin: number, dy: number, zMin: number, dz: number,
    f: (x: number, y: number, z: number) => number): number[] {
    const F: number[] = [];
    for (let s = 0; s < zBound; ++s) {
        for (let r = 0; r < yBound; ++r) {
            for (let c = 0; c < xBound; ++c) {
                F.push(f(xMin + c * dx, yMin + r * dy, zMin + s * dz));
            }
        }
    }
    return F;
}

describe('IntpAkimaUniform3', () => {
    it('throws for invalid inputs', () => {
        const F27 = new Array<number>(27).fill(0);
        expect(() => new IntpAkimaUniform3(2, 3, 3, 0, 1, 0, 1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(3, 2, 3, 0, 1, 0, 1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(3, 3, 2, 0, 1, 0, 1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(4, 3, 3, 0, 1, 0, 1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(3, 3, 3, 0, 0, 0, 1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(3, 3, 3, 0, 1, 0, -1, 0, 1, F27))
            .toThrow('Invalid input.');
        expect(() => new IntpAkimaUniform3(3, 3, 3, 0, 1, 0, 1, 0, 0, F27))
            .toThrow('Invalid input.');
    });

    it('provides member access', () => {
        const F = makeSamples(3, 4, 5, -1, 0.5, 2, 0.25, 0, 2, (x, y, z) => x + y + z);
        const interp = new IntpAkimaUniform3(3, 4, 5, -1, 0.5, 2, 0.25, 0, 2, F);
        expect(interp.getXBound()).toBe(3);
        expect(interp.getYBound()).toBe(4);
        expect(interp.getZBound()).toBe(5);
        expect(interp.getQuantity()).toBe(60);
        expect(interp.getF()).toBe(F);
        expect(interp.getXMin()).toBe(-1);
        expect(interp.getXMax()).toBeCloseTo(0, 14);
        expect(interp.getXSpacing()).toBe(0.5);
        expect(interp.getYMin()).toBe(2);
        expect(interp.getYMax()).toBeCloseTo(2.75, 14);
        expect(interp.getYSpacing()).toBe(0.25);
        expect(interp.getZMin()).toBe(0);
        expect(interp.getZMax()).toBeCloseTo(8, 14);
        expect(interp.getZSpacing()).toBe(2);
    });

    it('passes through the samples of a non-polynomial function', () => {
        const f = (x: number, y: number, z: number): number =>
            Math.sin(x) * Math.exp(0.4 * y) + Math.cos(0.7 * z);
        const xBound = 5, yBound = 4, zBound = 6;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25, zMin = 0, dz = 0.75;
        const F = makeSamples(xBound, yBound, zBound, xMin, dx, yMin, dy, zMin, dz, f);
        const interp = new IntpAkimaUniform3(xBound, yBound, zBound,
            xMin, dx, yMin, dy, zMin, dz, F);

        for (let s = 0; s < zBound; ++s) {
            for (let r = 0; r < yBound; ++r) {
                for (let c = 0; c < xBound; ++c) {
                    const x = xMin + c * dx, y = yMin + r * dy, z = zMin + s * dz;
                    const expected = F[c + xBound * (r + yBound * s)];
                    expect(interp.evaluate(x, y, z)).toBeCloseTo(expected, 10);
                    expect(interp.evaluate(0, 0, 0, x, y, z)).toBeCloseTo(expected, 10);
                }
            }
        }
    });

    it('is exact for an affine function everywhere', () => {
        // All slopes are constant, so every Akima first-derivative estimate is
        // exact and every mixed-partial stencil (one-sided or centered)
        // annihilates the affine function, including at the boundaries where
        // the upstream sign quirk would otherwise show up.
        const f = (x: number, y: number, z: number): number => 2 - 3 * x + 0.5 * y + 4 * z;
        const xBound = 4, yBound = 3, zBound = 5;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25, zMin = 1, dz = 2;
        const F = makeSamples(xBound, yBound, zBound, xMin, dx, yMin, dy, zMin, dz, f);
        const interp = new IntpAkimaUniform3(xBound, yBound, zBound,
            xMin, dx, yMin, dy, zMin, dz, F);

        for (const [x, y, z] of [[-0.9, 2.1, 1.3], [-0.2, 2.4, 5.7], [0.4, 2.45, 8.9]]) {
            expect(interp.evaluate(x, y, z)).toBeCloseTo(f(x, y, z), 10);
            expect(interp.evaluate(1, 0, 0, x, y, z)).toBeCloseTo(-3, 10);
            expect(interp.evaluate(0, 1, 0, x, y, z)).toBeCloseTo(0.5, 10);
            expect(interp.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(4, 10);
            expect(interp.evaluate(1, 1, 0, x, y, z)).toBeCloseTo(0, 10);
            expect(interp.evaluate(1, 0, 1, x, y, z)).toBeCloseTo(0, 10);
            expect(interp.evaluate(0, 1, 1, x, y, z)).toBeCloseTo(0, 10);
            expect(interp.evaluate(1, 1, 1, x, y, z)).toBeCloseTo(0, 10);
            expect(interp.evaluate(2, 0, 0, x, y, z)).toBeCloseTo(0, 10);
        }
    });

    it('exactly reproduces a tensor-product quadratic away from the max boundaries', () => {
        // f(x,y,z) = g(x)*h(y)*k(z) with quadratic factors. The uniform slopes
        // of a quadratic are linear in the index, so the Akima weighted
        // average reduces to the exact central difference (and the boundary
        // slope extrapolations stay exact); the one-sided and centered masks
        // used for the mixed partials are exact for quadratics as well. Cells
        // that avoid the max boundaries therefore reproduce f exactly (see the
        // sign-quirk test below for the max boundaries).
        const g = (x: number): number => 1 + 2 * x - 0.5 * x * x;
        const dg = (x: number): number => 2 - x;
        const h = (y: number): number => 2 - y + 0.75 * y * y;
        const dh = (y: number): number => -1 + 1.5 * y;
        const k = (z: number): number => 0.5 + 0.25 * z + 0.125 * z * z;
        const dk = (z: number): number => 0.25 + 0.25 * z;

        const xBound = 6, yBound = 6, zBound = 6;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25, zMin = 0, dz = 1;
        const F = makeSamples(xBound, yBound, zBound, xMin, dx, yMin, dy, zMin, dz,
            (x, y, z) => g(x) * h(y) * k(z));
        const interp = new IntpAkimaUniform3(xBound, yBound, zBound,
            xMin, dx, yMin, dy, zMin, dz, F);

        // Points inside cells whose corners avoid index xBound-1, yBound-1 and
        // zBound-1.
        for (const [x, y, z] of [[-0.9, 2.1, 0.4], [-0.1, 2.6, 1.7], [0.6, 2.8, 2.9]]) {
            expect(interp.evaluate(x, y, z)).toBeCloseTo(g(x) * h(y) * k(z), 9);
            expect(interp.evaluate(1, 0, 0, x, y, z))
                .toBeCloseTo(dg(x) * h(y) * k(z), 9);
            expect(interp.evaluate(0, 1, 0, x, y, z))
                .toBeCloseTo(g(x) * dh(y) * k(z), 9);
            expect(interp.evaluate(0, 0, 1, x, y, z))
                .toBeCloseTo(g(x) * h(y) * dk(z), 9);
            expect(interp.evaluate(1, 1, 0, x, y, z))
                .toBeCloseTo(dg(x) * dh(y) * k(z), 9);
            expect(interp.evaluate(1, 0, 1, x, y, z))
                .toBeCloseTo(dg(x) * h(y) * dk(z), 9);
            expect(interp.evaluate(0, 1, 1, x, y, z))
                .toBeCloseTo(g(x) * dh(y) * dk(z), 9);
            expect(interp.evaluate(1, 1, 1, x, y, z))
                .toBeCloseTo(dg(x) * dh(y) * dk(z), 9);
            expect(interp.evaluate(2, 0, 0, x, y, z))
                .toBeCloseTo(-1 * h(y) * k(z), 9);
            expect(interp.evaluate(0, 2, 0, x, y, z))
                .toBeCloseTo(g(x) * 1.5 * k(z), 9);
            expect(interp.evaluate(0, 0, 2, x, y, z))
                .toBeCloseTo(g(x) * h(y) * 0.25, 9);
        }
    });

    it('preserves the upstream max-boundary sign quirk in the mixed partials', () => {
        // GetFXY/GetFXZ/GetFYZ/GetFXYZ reuse the min-boundary one-sided
        // difference coefficients at the max boundaries with reflected sample
        // indices and without negating the mask, so each reflected direction
        // flips the sign of the estimated mixed partial there (the 3D analogue
        // of the IntpAkimaUniform2 issue). With f = x*y (constant in z) the
        // exact f_xy is +1, yet the estimates on the x-max column and the
        // y-max row come out as -1, and the pinned value below is the
        // resulting upstream value (it matches the IntpAkimaUniform2 port,
        // whose 2D value 1.1015625 is hand-computed from the algorithm).
        const F = makeSamples(3, 3, 3, 0, 1, 0, 1, 0, 1, (x, y) => x * y);
        const interp = new IntpAkimaUniform3(3, 3, 3, 0, 1, 0, 1, 0, 1, F);

        // Cells away from the x- and y-max boundaries are exact.
        for (const z of [0, 0.5, 2]) {
            expect(interp.evaluate(0.5, 0.5, z)).toBeCloseTo(0.25, 11);
            expect(interp.evaluate(1, 1, 0, 0.5, 0.5, z)).toBeCloseTo(1, 11);
            expect(interp.evaluate(0, 0, 1, 0.5, 0.5, z)).toBeCloseTo(0, 11);
        }
        // The cell touching the x-max boundary deviates.
        for (const z of [0, 0.5, 2]) {
            expect(interp.evaluate(1.5, 0.75, z)).toBeCloseTo(1.1015625, 11);
        }
        // The samples themselves are still interpolated exactly.
        for (let s = 0; s < 3; ++s) {
            for (let r = 0; r < 3; ++r) {
                for (let c = 0; c < 3; ++c) {
                    expect(interp.evaluate(c, r, s)).toBeCloseTo(c * r, 11);
                }
            }
        }
    });

    it('agrees with IntpAkimaUniform2 for z-independent samples', () => {
        // For samples that do not vary with z, every z-derivative estimate
        // vanishes and the tricubic cell polynomial collapses to the bicubic
        // one, so the 3D interpolator must agree with the independent 2D port.
        const xBound = 6, yBound = 5, zBound = 4;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.4, zMin = 3, dz = 1.25;
        let seed = 24680;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const slice: number[] = [];
        for (let i = 0; i < xBound * yBound; ++i) {
            slice.push(2 * rand() - 1);
        }
        const F: number[] = [];
        for (let s = 0; s < zBound; ++s) {
            F.push(...slice);
        }

        const interp3 = new IntpAkimaUniform3(xBound, yBound, zBound,
            xMin, dx, yMin, dy, zMin, dz, F);
        const interp2 = new IntpAkimaUniform2(xBound, yBound, xMin, dx, yMin, dy, slice);

        for (let trial = 0; trial < 100; ++trial) {
            const x = xMin + rand() * dx * (xBound - 1);
            const y = yMin + rand() * dy * (yBound - 1);
            const z = zMin + rand() * dz * (zBound - 1);
            expect(interp3.evaluate(x, y, z)).toBeCloseTo(interp2.evaluate(x, y), 10);
            expect(interp3.evaluate(1, 0, 0, x, y, z))
                .toBeCloseTo(interp2.evaluate(1, 0, x, y), 9);
            expect(interp3.evaluate(0, 1, 0, x, y, z))
                .toBeCloseTo(interp2.evaluate(0, 1, x, y), 9);
            expect(interp3.evaluate(1, 1, 0, x, y, z))
                .toBeCloseTo(interp2.evaluate(1, 1, x, y), 9);
            expect(interp3.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(0, 9);
        }
    });

    it('is symmetric under permuting the axes', () => {
        // The algorithm treats the three axes the same way, so transposing the
        // samples (and the corresponding grid parameters) transposes the
        // interpolant.
        const xBound = 5, yBound = 4, zBound = 6;
        const xMin = -1, dx = 0.5, yMin = 2, dy = 0.25, zMin = 0, dz = 0.75;
        const f = (x: number, y: number, z: number): number =>
            Math.sin(x + 0.3 * y) * (1 + 0.2 * z * z);
        const F = makeSamples(xBound, yBound, zBound, xMin, dx, yMin, dy, zMin, dz, f);
        // Swap the roles of x and z: G(c,r,s) = F(s,r,c).
        const G = makeSamples(zBound, yBound, xBound, zMin, dz, yMin, dy, xMin, dx,
            (z, y, x) => f(x, y, z));

        const interpF = new IntpAkimaUniform3(xBound, yBound, zBound,
            xMin, dx, yMin, dy, zMin, dz, F);
        const interpG = new IntpAkimaUniform3(zBound, yBound, xBound,
            zMin, dz, yMin, dy, xMin, dx, G);

        for (const [x, y, z] of [[-0.7, 2.1, 0.4], [0.2, 2.6, 1.9], [0.9, 2.7, 3.2]]) {
            expect(interpG.evaluate(z, y, x)).toBeCloseTo(interpF.evaluate(x, y, z), 10);
            expect(interpG.evaluate(0, 0, 1, z, y, x))
                .toBeCloseTo(interpF.evaluate(1, 0, 0, x, y, z), 9);
            expect(interpG.evaluate(1, 0, 0, z, y, x))
                .toBeCloseTo(interpF.evaluate(0, 0, 1, x, y, z), 9);
            expect(interpG.evaluate(1, 1, 1, z, y, x))
                .toBeCloseTo(interpF.evaluate(1, 1, 1, x, y, z), 9);
        }
    });

    it('is C1 across interior cell boundaries', () => {
        const f = (x: number, y: number, z: number): number =>
            Math.sin(x) * Math.exp(0.3 * y) + z * z * x;
        const bound = 6, min = 0, spacing = 0.5;
        const F = makeSamples(bound, bound, bound, min, spacing, min, spacing,
            min, spacing, f);
        const interp = new IntpAkimaUniform3(bound, bound, bound,
            min, spacing, min, spacing, min, spacing, F);

        const eps = 1e-7;
        const orders = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
        for (const [a, b] of [[0.7, 1.2], [1.6, 0.9], [2.1, 1.8]]) {
            for (const [xo, yo, zo] of orders) {
                // Cross the grid plane x = 1.5.
                expect(interp.evaluate(xo, yo, zo, 1.5 + eps, a, b)
                    - interp.evaluate(xo, yo, zo, 1.5 - eps, a, b)).toBeCloseTo(0, 5);
                // Cross the grid plane y = 1.0.
                expect(interp.evaluate(xo, yo, zo, a, 1 + eps, b)
                    - interp.evaluate(xo, yo, zo, a, 1 - eps, b)).toBeCloseTo(0, 5);
                // Cross the grid plane z = 2.0.
                expect(interp.evaluate(xo, yo, zo, a, b, 2 + eps)
                    - interp.evaluate(xo, yo, zo, a, b, 2 - eps)).toBeCloseTo(0, 5);
            }
        }
    });

    it('clamps evaluations to the domain', () => {
        const F = makeSamples(4, 4, 4, 0, 1, 0, 1, 0, 1,
            (x, y, z) => x * x + y - 0.5 * z);
        const interp = new IntpAkimaUniform3(4, 4, 4, 0, 1, 0, 1, 0, 1, F);
        expect(interp.evaluate(-10, 1.5, 2.5)).toBeCloseTo(interp.evaluate(0, 1.5, 2.5), 12);
        expect(interp.evaluate(10, 1.5, 2.5)).toBeCloseTo(interp.evaluate(3, 1.5, 2.5), 12);
        expect(interp.evaluate(1.5, -10, 2.5)).toBeCloseTo(interp.evaluate(1.5, 0, 2.5), 12);
        expect(interp.evaluate(1.5, 10, 2.5)).toBeCloseTo(interp.evaluate(1.5, 3, 2.5), 12);
        expect(interp.evaluate(1.5, 2.5, -10)).toBeCloseTo(interp.evaluate(1.5, 2.5, 0), 12);
        expect(interp.evaluate(1.5, 2.5, 10)).toBeCloseTo(interp.evaluate(1.5, 2.5, 3), 12);
    });

    it('returns zero for derivative orders beyond the degree', () => {
        const F = makeSamples(3, 3, 3, 0, 1, 0, 1, 0, 1, (x, y, z) => x + y + z);
        const interp = new IntpAkimaUniform3(3, 3, 3, 0, 1, 0, 1, 0, 1, F);
        expect(interp.evaluate(4, 0, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(0, 4, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(0, 0, 4, 0.5, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(-1, 0, 0, 0.5, 0.5, 0.5)).toBe(0);
        expect(interp.evaluate(3, 0, 0, 0.5, 0.5, 0.5)).toBeCloseTo(0, 10);
    });
});
