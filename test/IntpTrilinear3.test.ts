import { describe, it, expect } from 'vitest';
import { IntpTrilinear3 } from '../src/IntpTrilinear3.js';
import { check, expectClose, fc, finite, seededRandom }
    from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpTrilinear3.h.
// ---------------------------------------------------------------------------

interface TrilinearCase {
    xBound: number; yBound: number; zBound: number;
    xMin: number; xSpacing: number;
    yMin: number; ySpacing: number;
    zMin: number; zSpacing: number;
    F: number[];
    intp: IntpTrilinear3;
}

// Grid origins are integers and the spacings are exact binary fractions, so
// (x - min) * (1/spacing) is exact at the sample abscissae and the properties
// below need no allowance for index round-off.
const trilinearCase = fc.tuple(
    fc.integer({ min: 2, max: 4 }),
    fc.integer({ min: 2, max: 4 }),
    fc.integer({ min: 2, max: 4 }),
    fc.integer({ min: -3, max: 3 }),
    fc.constantFrom(0.5, 1, 2),
    fc.integer({ min: -3, max: 3 }),
    fc.constantFrom(0.5, 1, 2),
    fc.integer({ min: -3, max: 3 }),
    fc.constantFrom(0.5, 1, 2),
    fc.integer({ min: 0, max: 0xffff })
).map(([xBound, yBound, zBound, xMin, xSpacing, yMin, ySpacing, zMin,
    zSpacing, seed]): TrilinearCase => {
    const rand = seededRandom(seed + 1);
    const F: number[] = [];
    for (let i = 0; i < xBound * yBound * zBound; ++i) {
        F.push(Math.round(20 * (2 * rand() - 1)));
    }
    return {
        xBound, yBound, zBound, xMin, xSpacing, yMin, ySpacing, zMin, zSpacing,
        F,
        intp: new IntpTrilinear3(xBound, yBound, zBound, xMin, xSpacing, yMin,
            ySpacing, zMin, zSpacing, F)
    };
});

// An independent evaluation: clamp the cell index the way upstream does, then
// apply three successive one-dimensional linear blends. This exercises the
// blending matrix and the flat sample indexing without reusing either.
function trilinearReference(c: TrilinearCase, x: number, y: number,
    z: number): number {
    const clamp = (i: number, bound: number) =>
        i < 0 ? 0 : (i >= bound ? bound - 1 : i);
    const xIndex = (x - c.xMin) / c.xSpacing;
    const yIndex = (y - c.yMin) / c.ySpacing;
    const zIndex = (z - c.zMin) / c.zSpacing;
    const ix = clamp(Math.trunc(xIndex), c.xBound);
    const iy = clamp(Math.trunc(yIndex), c.yBound);
    const iz = clamp(Math.trunc(zIndex), c.zBound);
    const u = xIndex - ix, v = yIndex - iy, w = zIndex - iz;

    const at = (dx: number, dy: number, dz: number) =>
        c.F[clamp(ix + dx, c.xBound)
            + c.xBound * (clamp(iy + dy, c.yBound)
                + c.yBound * clamp(iz + dz, c.zBound))];
    const lerpX = (dy: number, dz: number) =>
        (1 - u) * at(0, dy, dz) + u * at(1, dy, dz);
    const lerpY = (dz: number) => (1 - v) * lerpX(0, dz) + v * lerpX(1, dz);
    return (1 - w) * lerpY(0) + w * lerpY(1);
}

describe('IntpTrilinear3 verification', () => {
    it('agrees with an independent successive-lerp evaluation, inside and '
        + 'outside the domain', () => {
            check(fc.tuple(trilinearCase, fc.integer({ min: 0, max: 0xffff })),
                ([c, seed]) => {
                    const rand = seededRandom(seed + 3);
                    for (let n = 0; n < 12; ++n) {
                        // Sample a box one cell wider than the domain in each
                        // direction so the clamping paths are exercised.
                        const x = c.xMin + c.xSpacing * (-1 + (c.xBound + 1) * rand());
                        const y = c.yMin + c.ySpacing * (-1 + (c.yBound + 1) * rand());
                        const z = c.zMin + c.zSpacing * (-1 + (c.zBound + 1) * rand());
                        expectClose(c.intp.evaluate(x, y, z),
                            trilinearReference(c, x, y, z), 1e-9, 1e-12);
                    }
                    return true;
                });
        });

    it('interpolates every sample of the grid', () => {
        check(trilinearCase, c => {
            for (let k = 0; k < c.zBound; ++k) {
                for (let j = 0; j < c.yBound; ++j) {
                    for (let i = 0; i < c.xBound; ++i) {
                        const value = c.intp.evaluate(c.xMin + c.xSpacing * i,
                            c.yMin + c.ySpacing * j, c.zMin + c.zSpacing * k);
                        expectClose(value,
                            c.F[i + c.xBound * (j + c.yBound * k)], 1e-9, 1e-12);
                    }
                }
            }
            return true;
        });
    });

    // Trilinear interpolation reproduces the eight-term trilinear polynomial
    // exactly, and its mixed partial d^3/dxdydz is the coefficient of xyz.
    it('reproduces a trilinear polynomial and its derivatives', () => {
        check(fc.tuple(trilinearCase,
            fc.array(fc.integer({ min: -5, max: 5 }),
                { minLength: 8, maxLength: 8 }),
            fc.integer({ min: 0, max: 0xffff })), ([c, k, seed]) => {
                const f = (x: number, y: number, z: number) =>
                    k[0] + k[1] * x + k[2] * y + k[3] * z + k[4] * x * y
                    + k[5] * x * z + k[6] * y * z + k[7] * x * y * z;
                const F: number[] = [];
                for (let s = 0; s < c.zBound; ++s) {
                    for (let r = 0; r < c.yBound; ++r) {
                        for (let q = 0; q < c.xBound; ++q) {
                            F.push(f(c.xMin + c.xSpacing * q,
                                c.yMin + c.ySpacing * r,
                                c.zMin + c.zSpacing * s));
                        }
                    }
                }
                const intp = new IntpTrilinear3(c.xBound, c.yBound, c.zBound,
                    c.xMin, c.xSpacing, c.yMin, c.ySpacing, c.zMin, c.zSpacing, F);
                const rand = seededRandom(seed + 5);
                for (let n = 0; n < 8; ++n) {
                    const x = c.xMin + c.xSpacing * (c.xBound - 1) * rand();
                    const y = c.yMin + c.ySpacing * (c.yBound - 1) * rand();
                    const z = c.zMin + c.zSpacing * (c.zBound - 1) * rand();
                    expectClose(intp.evaluate(x, y, z), f(x, y, z), 1e-8, 1e-11);
                    expectClose(intp.evaluate(1, 0, 0, x, y, z),
                        k[1] + k[4] * y + k[5] * z + k[7] * y * z, 1e-8, 1e-11);
                    expectClose(intp.evaluate(1, 1, 1, x, y, z), k[7],
                        1e-8, 1e-11);
                }
                return true;
            });
    });

    // The samples all carry weight one, so a constant field is reproduced
    // everywhere, including outside the domain where the stencil is clamped.
    it('reproduces a constant field everywhere', () => {
        check(fc.tuple(trilinearCase, fc.integer({ min: -9, max: 9 }),
            fc.integer({ min: 0, max: 0xffff })), ([c, value, seed]) => {
                const F = new Array<number>(c.xBound * c.yBound * c.zBound)
                    .fill(value);
                const intp = new IntpTrilinear3(c.xBound, c.yBound, c.zBound,
                    c.xMin, c.xSpacing, c.yMin, c.ySpacing, c.zMin, c.zSpacing, F);
                const rand = seededRandom(seed + 9);
                for (let n = 0; n < 8; ++n) {
                    const x = c.xMin + c.xSpacing * (-2 + (c.xBound + 3) * rand());
                    const y = c.yMin + c.ySpacing * (-2 + (c.yBound + 3) * rand());
                    const z = c.zMin + c.zSpacing * (-2 + (c.zBound + 3) * rand());
                    expectClose(intp.evaluate(x, y, z), value, 1e-9, 1e-12);
                }
                return true;
            });
    });

    // Upstream (issue #69): the evaluators claim to clamp the inputs to the
    // domain but only clamp the cell index. Above the maximum both entries of
    // the stencil collapse onto the last sample and the two blend weights sum
    // to one, so the boundary value is held; below the minimum the stencil
    // still spans two distinct samples and the boundary cell is extrapolated.
    it('holds the value above the domain and extrapolates below it', () => {
        check(fc.tuple(trilinearCase, finite(0.01, 5)), ([c, delta]) => {
            const yMid = c.yMin + c.ySpacing * (c.yBound - 1) / 2;
            const zMid = c.zMin + c.zSpacing * (c.zBound - 1) / 2;
            const xMax = c.intp.getXMax();
            expectClose(c.intp.evaluate(xMax + delta, yMid, zMid),
                c.intp.evaluate(xMax, yMid, zMid), 1e-9, 1e-12);

            // Below the minimum the interpolant continues the first cell's
            // linear polynomial: f(xMin - d) = f(xMin) - d * f'(xMin).
            const xMin = c.intp.getXMin();
            const slope = c.intp.evaluate(1, 0, 0, xMin, yMid, zMid);
            if (delta < c.xSpacing) {
                expectClose(c.intp.evaluate(xMin - delta, yMid, zMid),
                    c.intp.evaluate(xMin, yMid, zMid) - delta * slope,
                    1e-9, 1e-11);
            }
            return true;
        });
    });

    // The trilinear polynomial is degree one in each variable, so the second
    // derivative in any variable vanishes; upstream returns zero for any order
    // above one without even evaluating the stencil.
    it('returns zero for derivative orders outside [0,1]', () => {
        check(fc.tuple(trilinearCase, fc.integer({ min: 2, max: 6 }),
            fc.integer({ min: 0, max: 2 })), ([c, order, axis]) => {
                const orders = [0, 0, 0];
                orders[axis] = order;
                const x = c.xMin + c.xSpacing * 0.5;
                const y = c.yMin + c.ySpacing * 0.5;
                const z = c.zMin + c.zSpacing * 0.5;
                expect(c.intp.evaluate(orders[0], orders[1], orders[2], x, y, z))
                    .toBe(0);
                return true;
            });
    });
});
