import { describe, it, expect } from 'vitest';
import { IntpTrilinear3 } from '../src/IntpTrilinear3';

// Grid used by most of the tests. The bounds, origins and spacings are all
// distinct so that a transposed index or a swapped spacing is caught.
const XB = 5, YB = 4, ZB = 6;
const XMIN = -1, XSP = 0.5;
const YMIN = 2, YSP = 0.25;
const ZMIN = -3, ZSP = 1.5;

type Fn3 = (x: number, y: number, z: number) => number;

// Sample f on the lattice in the lexicographical order the class expects:
// F[c + xBound*(r + yBound*s)].
function sampleGrid(f: Fn3, xb = XB, yb = YB, zb = ZB): number[] {
    const F = new Array<number>(xb * yb * zb);
    for (let s = 0; s < zb; ++s) {
        for (let r = 0; r < yb; ++r) {
            for (let c = 0; c < xb; ++c) {
                F[c + xb * (r + yb * s)] = f(XMIN + XSP * c, YMIN + YSP * r, ZMIN + ZSP * s);
            }
        }
    }
    return F;
}

function makeInterp(F: ArrayLike<number>): IntpTrilinear3 {
    return new IntpTrilinear3(XB, YB, ZB, XMIN, XSP, YMIN, YSP, ZMIN, ZSP, F);
}

// An independent trilinear evaluation written in the (1-t, t) corner-weight
// form rather than with the blending matrix of the port.
function referenceTrilinear(F: ArrayLike<number>, x: number, y: number, z: number): number {
    const tx = (x - XMIN) / XSP, ty = (y - YMIN) / YSP, tz = (z - ZMIN) / ZSP;
    const ix = Math.min(Math.max(Math.trunc(tx), 0), XB - 1);
    const iy = Math.min(Math.max(Math.trunc(ty), 0), YB - 1);
    const iz = Math.min(Math.max(Math.trunc(tz), 0), ZB - 1);
    const u = tx - ix, v = ty - iy, w = tz - iz;

    let result = 0;
    for (let k = 0; k < 2; ++k) {
        const bz = (k === 0 ? 1 - w : w);
        const zc = Math.min(iz + k, ZB - 1);
        for (let j = 0; j < 2; ++j) {
            const by = (j === 0 ? 1 - v : v);
            const yc = Math.min(iy + j, YB - 1);
            for (let i = 0; i < 2; ++i) {
                const bx = (i === 0 ? 1 - u : u);
                const xc = Math.min(ix + i, XB - 1);
                result += bx * by * bz * F[xc + XB * (yc + YB * zc)];
            }
        }
    }
    return result;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// A general (non-polynomial) sample function.
const smooth: Fn3 = (x, y, z) => Math.sin(1.3 * x) * Math.exp(-0.2 * y) + Math.cos(0.7 * z) - 0.3 * x * z;

// A trilinear polynomial, which the interpolator must reproduce exactly.
const trilinearPoly: Fn3 = (x, y, z) =>
    2 - 3 * x + 0.5 * y + 4 * z + 1.5 * x * y - 2 * x * z + 0.75 * y * z + 1.25 * x * y * z;

describe('IntpTrilinear3', () => {
    describe('construction', () => {
        it('exposes the grid parameters and the derived maxima', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            expect(interp.getXBound()).toBe(XB);
            expect(interp.getYBound()).toBe(YB);
            expect(interp.getZBound()).toBe(ZB);
            expect(interp.getQuantity()).toBe(XB * YB * ZB);
            expect(interp.getF()).toBe(F);
            expect(interp.getXMin()).toBe(XMIN);
            expect(interp.getYMin()).toBe(YMIN);
            expect(interp.getZMin()).toBe(ZMIN);
            expect(interp.getXSpacing()).toBe(XSP);
            expect(interp.getYSpacing()).toBe(YSP);
            expect(interp.getZSpacing()).toBe(ZSP);
            expect(interp.getXMax()).toBeCloseTo(XMIN + XSP * (XB - 1), 12);
            expect(interp.getYMax()).toBeCloseTo(YMIN + YSP * (YB - 1), 12);
            expect(interp.getZMax()).toBeCloseTo(ZMIN + ZSP * (ZB - 1), 12);
        });

        it('requires at least a 2x2x2 block of samples', () => {
            const F = new Array<number>(2 * 2 * 2).fill(0);
            expect(() => new IntpTrilinear3(1, 2, 2, 0, 1, 0, 1, 0, 1, F)).toThrow('Invalid input.');
            expect(() => new IntpTrilinear3(2, 1, 2, 0, 1, 0, 1, 0, 1, F)).toThrow('Invalid input.');
            expect(() => new IntpTrilinear3(2, 2, 1, 0, 1, 0, 1, 0, 1, F)).toThrow('Invalid input.');
            expect(() => new IntpTrilinear3(2, 2, 2, 0, 1, 0, 1, 0, 1, F)).not.toThrow();
        });

        it('requires positive spacings', () => {
            const F = new Array<number>(8).fill(0);
            expect(() => new IntpTrilinear3(2, 2, 2, 0, 0, 0, 1, 0, 1, F)).toThrow('Invalid input.');
            expect(() => new IntpTrilinear3(2, 2, 2, 0, 1, 0, -1, 0, 1, F)).toThrow('Invalid input.');
            expect(() => new IntpTrilinear3(2, 2, 2, 0, 1, 0, 1, 0, 0, F)).toThrow('Invalid input.');
        });

        it('requires enough samples for the bounds', () => {
            const F = new Array<number>(7).fill(0);
            expect(() => new IntpTrilinear3(2, 2, 2, 0, 1, 0, 1, 0, 1, F)).toThrow('Invalid input.');
        });

        it('accepts a typed array of samples', () => {
            const F = Float64Array.from(sampleGrid(trilinearPoly));
            const interp = new IntpTrilinear3(XB, YB, ZB, XMIN, XSP, YMIN, YSP, ZMIN, ZSP, F);
            expect(interp.evaluate(0.3, 2.4, -1.1)).toBeCloseTo(trilinearPoly(0.3, 2.4, -1.1), 10);
        });

        it('references the samples rather than copying them', () => {
            const F = sampleGrid(() => 0);
            const interp = makeInterp(F);
            expect(interp.evaluate(0, 2.5, -1)).toBe(0);
            F.fill(5);
            expect(interp.evaluate(0, 2.5, -1)).toBeCloseTo(5, 12);
        });
    });

    describe('function evaluation', () => {
        it('passes through every sample point', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            for (let s = 0; s < ZB; ++s) {
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        const value = interp.evaluate(XMIN + XSP * c, YMIN + YSP * r, ZMIN + ZSP * s);
                        expect(value).toBeCloseTo(F[c + XB * (r + YB * s)], 12);
                    }
                }
            }
        });

        it('reproduces a trilinear polynomial exactly on the domain', () => {
            const F = sampleGrid(trilinearPoly);
            const interp = makeInterp(F);
            const rand = makeRandom(12345);
            for (let n = 0; n < 200; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(trilinearPoly(x, y, z), 10);
            }
        });

        it('agrees with an independent corner-weight evaluation', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const rand = makeRandom(999);
            for (let n = 0; n < 200; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(referenceTrilinear(F, x, y, z), 12);
            }
        });

        it('averages the eight corners at a cell center', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const c = 1, r = 2, s = 3;
            let average = 0;
            for (let k = 0; k <= 1; ++k) {
                for (let j = 0; j <= 1; ++j) {
                    for (let i = 0; i <= 1; ++i) {
                        average += F[(c + i) + XB * ((r + j) + YB * (s + k))];
                    }
                }
            }
            average /= 8;
            const value = interp.evaluate(XMIN + XSP * (c + 0.5), YMIN + YSP * (r + 0.5),
                ZMIN + ZSP * (s + 0.5));
            expect(value).toBeCloseTo(average, 12);
        });

        it('handles the minimal 2x2x2 grid', () => {
            // A single cell: the interpolant is the corner-weight formula.
            const F = [0, 1, 2, 3, 4, 5, 6, 7];
            const interp = new IntpTrilinear3(2, 2, 2, 0, 1, 0, 1, 0, 1, F);
            for (let k = 0; k < 2; ++k) {
                for (let j = 0; j < 2; ++j) {
                    for (let i = 0; i < 2; ++i) {
                        expect(interp.evaluate(i, j, k)).toBeCloseTo(F[i + 2 * (j + 2 * k)], 12);
                    }
                }
            }
            // The cell center is the average of the eight corners.
            expect(interp.evaluate(0.5, 0.5, 0.5)).toBeCloseTo(3.5, 12);
            // F is 1*i + 2*j + 4*k, so the partials are the coefficients.
            expect(interp.evaluate(0.25, 0.5, 0.75)).toBeCloseTo(0.25 + 2 * 0.5 + 4 * 0.75, 12);
            expect(interp.evaluate(1, 0, 0, 0.25, 0.5, 0.75)).toBeCloseTo(1, 12);
            expect(interp.evaluate(0, 1, 0, 0.25, 0.5, 0.75)).toBeCloseTo(2, 12);
            expect(interp.evaluate(0, 0, 1, 0.25, 0.5, 0.75)).toBeCloseTo(4, 12);
        });

        it('reproduces a constant field everywhere, including outside the domain', () => {
            const F = sampleGrid(() => -2.75);
            const interp = makeInterp(F);
            expect(interp.evaluate(XMIN, YMIN, ZMIN)).toBeCloseTo(-2.75, 12);
            expect(interp.evaluate(0.1, 2.4, -1.7)).toBeCloseTo(-2.75, 12);
            expect(interp.evaluate(-10, 100, 50)).toBeCloseTo(-2.75, 12);
        });
    });

    describe('derivative evaluation', () => {
        it('returns the function value when all orders are zero', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const rand = makeRandom(7);
            for (let n = 0; n < 50; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(0, 0, 0, x, y, z)).toBe(interp.evaluate(x, y, z));
            }
        });

        it('matches the analytic gradient of a trilinear polynomial', () => {
            const F = sampleGrid(trilinearPoly);
            const interp = makeInterp(F);
            const dfdx = (x: number, y: number, z: number) => -3 + 1.5 * y - 2 * z + 1.25 * y * z;
            const dfdy = (x: number, y: number, z: number) => 0.5 + 1.5 * x + 0.75 * z + 1.25 * x * z;
            const dfdz = (x: number, y: number, z: number) => 4 - 2 * x + 0.75 * y + 1.25 * x * y;
            const rand = makeRandom(4242);
            for (let n = 0; n < 100; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(1, 0, 0, x, y, z)).toBeCloseTo(dfdx(x, y, z), 10);
                expect(interp.evaluate(0, 1, 0, x, y, z)).toBeCloseTo(dfdy(x, y, z), 10);
                expect(interp.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(dfdz(x, y, z), 10);
            }
        });

        it('matches central differences of the interpolant inside a cell', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const h = 1e-4;
            const points: Array<[number, number, number]> = [
                [XMIN + 0.3 * XSP, YMIN + 0.7 * YSP, ZMIN + 0.4 * ZSP],
                [XMIN + 2.4 * XSP, YMIN + 1.5 * YSP, ZMIN + 3.6 * ZSP],
                [XMIN + 3.5 * XSP, YMIN + 2.2 * YSP, ZMIN + 4.5 * ZSP]
            ];
            for (const [x, y, z] of points) {
                const fdx = (interp.evaluate(x + h, y, z) - interp.evaluate(x - h, y, z)) / (2 * h);
                const fdy = (interp.evaluate(x, y + h, z) - interp.evaluate(x, y - h, z)) / (2 * h);
                const fdz = (interp.evaluate(x, y, z + h) - interp.evaluate(x, y, z - h)) / (2 * h);
                expect(interp.evaluate(1, 0, 0, x, y, z)).toBeCloseTo(fdx, 7);
                expect(interp.evaluate(0, 1, 0, x, y, z)).toBeCloseTo(fdy, 7);
                expect(interp.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(fdz, 7);
            }
        });

        it('matches central differences for the mixed partials', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const h = 1e-3;
            const x = XMIN + 1.4 * XSP, y = YMIN + 2.3 * YSP, z = ZMIN + 3.6 * ZSP;

            const fdxy = (interp.evaluate(1, 0, 0, x, y + h, z)
                - interp.evaluate(1, 0, 0, x, y - h, z)) / (2 * h);
            expect(interp.evaluate(1, 1, 0, x, y, z)).toBeCloseTo(fdxy, 7);

            const fdxz = (interp.evaluate(1, 0, 0, x, y, z + h)
                - interp.evaluate(1, 0, 0, x, y, z - h)) / (2 * h);
            expect(interp.evaluate(1, 0, 1, x, y, z)).toBeCloseTo(fdxz, 7);

            const fdxyz = (interp.evaluate(1, 1, 0, x, y, z + h)
                - interp.evaluate(1, 1, 0, x, y, z - h)) / (2 * h);
            expect(interp.evaluate(1, 1, 1, x, y, z)).toBeCloseTo(fdxyz, 7);
        });

        it('reproduces the analytic mixed partial of a trilinear polynomial', () => {
            const F = sampleGrid(trilinearPoly);
            const interp = makeInterp(F);
            const x = 0.13, y = 2.4, z = -1.1;
            // d2f/dxdy = 1.5 + 1.25*z, d3f/dxdydz = 1.25.
            expect(interp.evaluate(1, 1, 0, x, y, z)).toBeCloseTo(1.5 + 1.25 * z, 10);
            expect(interp.evaluate(1, 1, 1, x, y, z)).toBeCloseTo(1.25, 10);
        });

        it('returns zero for orders outside [0,1]', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const x = 0.13, y = 2.4, z = -1.1;
            expect(interp.evaluate(2, 0, 0, x, y, z)).toBe(0);
            expect(interp.evaluate(0, 2, 0, x, y, z)).toBe(0);
            expect(interp.evaluate(0, 0, 2, x, y, z)).toBe(0);
            expect(interp.evaluate(1, 1, 3, x, y, z)).toBe(0);
            expect(interp.evaluate(-1, 0, 0, x, y, z)).toBe(0);
        });

        it('scales the derivative with the inverse spacing', () => {
            // f(x) = x on a grid whose x-spacing is halved must still have
            // df/dx = 1.
            const F = sampleGrid((x) => x);
            const coarse = makeInterp(F);
            expect(coarse.evaluate(1, 0, 0, 0.1, 2.4, -1.1)).toBeCloseTo(1, 10);

            const G = new Array<number>(XB * YB * ZB);
            for (let s = 0; s < ZB; ++s) {
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        G[c + XB * (r + YB * s)] = XMIN + 0.25 * c;
                    }
                }
            }
            const fine = new IntpTrilinear3(XB, YB, ZB, XMIN, 0.25, YMIN, YSP, ZMIN, ZSP, G);
            expect(fine.evaluate(1, 0, 0, XMIN + 0.1, 2.4, -1.1)).toBeCloseTo(1, 10);
        });
    });

    describe('boundary behavior', () => {
        it('evaluates the extreme corners as the corner samples', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const xMax = interp.getXMax(), yMax = interp.getYMax(), zMax = interp.getZMax();
            expect(interp.evaluate(xMax, yMax, zMax)).toBeCloseTo(F[XB * YB * ZB - 1], 12);
            expect(interp.evaluate(XMIN, YMIN, ZMIN)).toBeCloseTo(F[0], 12);
            expect(interp.evaluate(xMax, YMIN, ZMIN)).toBeCloseTo(F[XB - 1], 12);
        });

        it('collapses the stencil at the upper boundary, giving a zero derivative there', () => {
            // At exactly x = xMax the cell index is the last sample index, so
            // both columns of the 2x2x2 stencil clamp to the same x sample.
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            expect(interp.evaluate(1, 0, 0, interp.getXMax(), 2.4, -1.1)).toBeCloseTo(0, 10);
            expect(interp.evaluate(0, 1, 0, 0.1, interp.getYMax(), -1.1)).toBeCloseTo(0, 10);
            expect(interp.evaluate(0, 0, 1, 0.1, 2.4, interp.getZMax())).toBeCloseTo(0, 10);
        });

        it('extrapolates below the domain but holds the value above it', () => {
            // The evaluators clamp the cell index but not the fractional cell
            // coordinate. Below the minimum the 2x2x2 stencil still spans two
            // distinct samples per axis, so a linear field is extrapolated
            // exactly. At or above the maximum the cell index is the last
            // sample index, so both entries of that axis clamp to the same
            // sample and the two blend weights sum to one: the value is held
            // at the boundary sample no matter how far outside the input is.
            const linear = (x: number, y: number, z: number) => 2 * x - 3 * y + 0.5 * z;
            const F = sampleGrid(linear);
            const interp = makeInterp(F);

            const below: Array<[number, number, number]> = [
                [XMIN - 0.4 * XSP, YMIN + 0.2, ZMIN + 0.3],
                [XMIN - 3 * XSP, YMIN - 2 * YSP, ZMIN - 4 * ZSP]
            ];
            for (const [x, y, z] of below) {
                expect(interp.evaluate(x, y, z)).toBeCloseTo(linear(x, y, z), 10);
            }

            const xMax = interp.getXMax(), yMax = interp.getYMax(), zMax = interp.getZMax();
            const held = linear(xMax, yMax, zMax);
            expect(interp.evaluate(xMax + 2 * XSP, yMax + 1.5 * YSP, zMax + 0.5 * ZSP))
                .toBeCloseTo(held, 10);
            expect(interp.evaluate(xMax + 100, yMax + 100, zMax + 100)).toBeCloseTo(held, 10);
            // Mixing an inside axis with an outside one clamps only that axis.
            const yInside = YMIN + 1.25 * YSP;
            expect(interp.evaluate(xMax + 7, yInside, zMax + 7))
                .toBeCloseTo(linear(xMax, yInside, zMax), 10);
        });

        it('keeps the derivative of a linear field constant outside the domain', () => {
            const F = sampleGrid((x, y, z) => 2 * x - 3 * y + 0.5 * z);
            const interp = makeInterp(F);
            expect(interp.evaluate(1, 0, 0, XMIN - 5, YMIN - 5, ZMIN - 5)).toBeCloseTo(2, 10);
            expect(interp.evaluate(0, 1, 0, XMIN - 5, YMIN - 5, ZMIN - 5)).toBeCloseTo(-3, 10);
            expect(interp.evaluate(0, 0, 1, XMIN - 5, YMIN - 5, ZMIN - 5)).toBeCloseTo(0.5, 10);
        });

        it('is continuous across a cell boundary', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F);
            const xEdge = XMIN + 2 * XSP;
            const y = YMIN + 1.3 * YSP, z = ZMIN + 2.7 * ZSP;
            const eps = 1e-10;
            const left = interp.evaluate(xEdge - eps, y, z);
            const right = interp.evaluate(xEdge + eps, y, z);
            expect(left).toBeCloseTo(right, 8);
            expect(left).toBeCloseTo(interp.evaluate(xEdge, y, z), 8);
        });
    });
});
