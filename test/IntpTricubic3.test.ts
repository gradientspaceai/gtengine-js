import { describe, it, expect } from 'vitest';
import { IntpTricubic3 } from '../src/IntpTricubic3.js';

// Grid used by most of the tests. The bounds, origins and spacings are all
// distinct so that a transposed index or a swapped spacing is caught.
const XB = 6, YB = 7, ZB = 6;
const XMIN = -1, XSP = 0.5;
const YMIN = 2, YSP = 0.25;
const ZMIN = -3, ZSP = 1.5;

type Fn3 = (x: number, y: number, z: number) => number;

// Sample f on the lattice in the lexicographical order the class expects:
// F[c + xBound*(r + yBound*s)].
function sampleGrid(f: Fn3): number[] {
    const F = new Array<number>(XB * YB * ZB);
    for (let s = 0; s < ZB; ++s) {
        for (let r = 0; r < YB; ++r) {
            for (let c = 0; c < XB; ++c) {
                F[c + XB * (r + YB * s)] = f(XMIN + XSP * c, YMIN + YSP * r, ZMIN + ZSP * s);
            }
        }
    }
    return F;
}

function makeInterp(F: ArrayLike<number>, catmullRom: boolean): IntpTricubic3 {
    return new IntpTricubic3(XB, YB, ZB, XMIN, XSP, YMIN, YSP, ZMIN, ZSP, F, catmullRom);
}

// Catmull-Rom weights derived independently of the port's blending matrix,
// from the cubic Hermite form p(t) = h00*P1 + h10*m0 + h01*P2 + h11*m1 with
// tangents m0 = (P2 - P0)/2 and m1 = (P3 - P1)/2.
function catmullRomWeights(t: number): number[] {
    const t2 = t * t, t3 = t2 * t;
    const h00 = 2 * t3 - 3 * t2 + 1;
    const h10 = t3 - 2 * t2 + t;
    const h01 = -2 * t3 + 3 * t2;
    const h11 = t3 - t2;
    return [-0.5 * h10, h00 - 0.5 * h11, h01 + 0.5 * h10, 0.5 * h11];
}

// The uniform cubic B-spline basis in closed form.
function bsplineWeights(t: number): number[] {
    const t2 = t * t, t3 = t2 * t;
    const omt = 1 - t;
    return [
        omt * omt * omt / 6,
        (3 * t3 - 6 * t2 + 4) / 6,
        (-3 * t3 + 3 * t2 + 3 * t + 1) / 6,
        t3 / 6
    ];
}

function clampIndex(index: number, bound: number): number {
    return Math.min(Math.max(index, 0), bound - 1);
}

// An independent tricubic evaluation using the closed-form bases above.
function referenceTricubic(F: ArrayLike<number>, x: number, y: number, z: number,
    catmullRom: boolean): number {
    const weights = catmullRom ? catmullRomWeights : bsplineWeights;
    const tx = (x - XMIN) / XSP, ty = (y - YMIN) / YSP, tz = (z - ZMIN) / ZSP;
    const ix = clampIndex(Math.trunc(tx), XB);
    const iy = clampIndex(Math.trunc(ty), YB);
    const iz = clampIndex(Math.trunc(tz), ZB);
    const wx = weights(tx - ix), wy = weights(ty - iy), wz = weights(tz - iz);

    let result = 0;
    for (let k = 0; k < 4; ++k) {
        const zc = clampIndex(iz - 1 + k, ZB);
        for (let j = 0; j < 4; ++j) {
            const yc = clampIndex(iy - 1 + j, YB);
            for (let i = 0; i < 4; ++i) {
                const xc = clampIndex(ix - 1 + i, XB);
                result += wx[i] * wy[j] * wz[k] * F[xc + XB * (yc + YB * zc)];
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

// Points whose 4x4x4 stencil lies entirely inside the grid, so no clamping
// alters the reproduction properties of the scheme. The fractional cell
// coordinate stays away from the cell walls so that the finite-difference
// checks below do not straddle a cell boundary.
function interiorPoint(rand: () => number): [number, number, number] {
    const pick = (bound: number, min: number, spacing: number): number => {
        const cellCount = bound - 3;  // cells with index 1 .. bound-3
        const cell = 1 + Math.min(Math.floor(rand() * cellCount), cellCount - 1);
        return min + spacing * (cell + 0.15 + 0.7 * rand());
    };
    return [pick(XB, XMIN, XSP), pick(YB, YMIN, YSP), pick(ZB, ZMIN, ZSP)];
}

// A general (non-polynomial) sample function.
const smooth: Fn3 = (x, y, z) => Math.sin(1.3 * x) * Math.exp(-0.2 * y) + Math.cos(0.7 * z) - 0.3 * x * z;

// Degree at most two in each variable: reproduced exactly by Catmull-Rom.
const triquadratic: Fn3 = (x, y, z) =>
    1 + 0.5 * x - 0.75 * x * x + 0.25 * y + 0.4 * y * y - 0.6 * z + 0.2 * z * z
    + 1.1 * x * y - 0.3 * x * x * z * z + 0.7 * x * y * z + 0.9 * x * x * y * y * z * z;

// Degree at most one in each variable: reproduced exactly by both blendings.
const trilinearPoly: Fn3 = (x, y, z) =>
    2 - 3 * x + 0.5 * y + 4 * z + 1.5 * x * y - 2 * x * z + 0.75 * y * z + 1.25 * x * y * z;

describe('IntpTricubic3', () => {
    describe('construction', () => {
        it('exposes the grid parameters and the derived maxima', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
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

        it('requires at least a 4x4x4 block of samples', () => {
            const F = new Array<number>(4 * 4 * 4).fill(0);
            expect(() => new IntpTricubic3(3, 4, 4, 0, 1, 0, 1, 0, 1, F, true)).toThrow('Invalid input.');
            expect(() => new IntpTricubic3(4, 3, 4, 0, 1, 0, 1, 0, 1, F, true)).toThrow('Invalid input.');
            expect(() => new IntpTricubic3(4, 4, 3, 0, 1, 0, 1, 0, 1, F, true)).toThrow('Invalid input.');
            expect(() => new IntpTricubic3(4, 4, 4, 0, 1, 0, 1, 0, 1, F, true)).not.toThrow();
        });

        it('requires positive spacings and enough samples', () => {
            const F = new Array<number>(64).fill(0);
            expect(() => new IntpTricubic3(4, 4, 4, 0, 0, 0, 1, 0, 1, F, false)).toThrow('Invalid input.');
            expect(() => new IntpTricubic3(4, 4, 4, 0, 1, 0, -2, 0, 1, F, false)).toThrow('Invalid input.');
            expect(() => new IntpTricubic3(4, 4, 4, 0, 1, 0, 1, 0, 0, F, false)).toThrow('Invalid input.');
            const short = new Array<number>(63).fill(0);
            expect(() => new IntpTricubic3(4, 4, 4, 0, 1, 0, 1, 0, 1, short, true)).toThrow('Invalid input.');
        });

        it('accepts a typed array of samples', () => {
            const F = Float64Array.from(sampleGrid(triquadratic));
            const interp = new IntpTricubic3(XB, YB, ZB, XMIN, XSP, YMIN, YSP, ZMIN, ZSP, F, true);
            expect(interp.evaluate(0.3, 2.6, 0.4)).toBeCloseTo(triquadratic(0.3, 2.6, 0.4), 9);
        });

        it('references the samples rather than copying them', () => {
            const F = sampleGrid(() => 0);
            const interp = makeInterp(F, true);
            expect(interp.evaluate(0.1, 2.6, 0.4)).toBe(0);
            F.fill(5);
            expect(interp.evaluate(0.1, 2.6, 0.4)).toBeCloseTo(5, 12);
        });

        it('produces different results for the two blending matrices', () => {
            const F = sampleGrid(smooth);
            const catmull = makeInterp(F, true);
            const bspline = makeInterp(F, false);
            const [x, y, z] = interiorPoint(makeRandom(3));
            expect(Math.abs(catmull.evaluate(x, y, z) - bspline.evaluate(x, y, z)))
                .toBeGreaterThan(1e-6);
        });
    });

    describe('Catmull-Rom blending', () => {
        it('passes through every sample point', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            for (let s = 0; s < ZB; ++s) {
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        const value = interp.evaluate(XMIN + XSP * c, YMIN + YSP * r, ZMIN + ZSP * s);
                        expect(value).toBeCloseTo(F[c + XB * (r + YB * s)], 12);
                    }
                }
            }
        });

        it('reproduces a polynomial of degree two per variable in the interior', () => {
            const F = sampleGrid(triquadratic);
            const interp = makeInterp(F, true);
            const rand = makeRandom(2024);
            for (let n = 0; n < 200; ++n) {
                const [x, y, z] = interiorPoint(rand);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(triquadratic(x, y, z), 9);
            }
        });

        it('agrees with an independent Hermite-form evaluation over the whole domain', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const rand = makeRandom(555);
            for (let n = 0; n < 200; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(referenceTricubic(F, x, y, z, true), 11);
            }
        });

        it('has a continuous first derivative across a cell boundary', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const xEdge = XMIN + 3 * XSP;
            const y = YMIN + 2.3 * YSP, z = ZMIN + 2.7 * ZSP;
            const eps = 1e-7;
            expect(interp.evaluate(1, 0, 0, xEdge - eps, y, z))
                .toBeCloseTo(interp.evaluate(1, 0, 0, xEdge + eps, y, z), 5);
        });
    });

    describe('B-spline blending', () => {
        it('reproduces a polynomial of degree one per variable in the interior', () => {
            const F = sampleGrid(trilinearPoly);
            const interp = makeInterp(F, false);
            const rand = makeRandom(31337);
            for (let n = 0; n < 200; ++n) {
                const [x, y, z] = interiorPoint(rand);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(trilinearPoly(x, y, z), 9);
            }
        });

        it('smooths rather than interpolates the samples', () => {
            // At an interior lattice point the tensor weights are
            // (1/6, 4/6, 1/6) in each direction.
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, false);
            const c = 2, r = 3, s = 2;
            const w = [1 / 6, 4 / 6, 1 / 6];
            let expected = 0;
            for (let k = 0; k < 3; ++k) {
                for (let j = 0; j < 3; ++j) {
                    for (let i = 0; i < 3; ++i) {
                        expected += w[i] * w[j] * w[k]
                            * F[(c - 1 + i) + XB * ((r - 1 + j) + YB * (s - 1 + k))];
                    }
                }
            }
            const value = interp.evaluate(XMIN + XSP * c, YMIN + YSP * r, ZMIN + ZSP * s);
            expect(value).toBeCloseTo(expected, 12);
            expect(value).not.toBeCloseTo(F[c + XB * (r + YB * s)], 6);
        });

        it('agrees with an independent basis-function evaluation over the whole domain', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, false);
            const rand = makeRandom(4711);
            for (let n = 0; n < 200; ++n) {
                const x = XMIN + rand() * XSP * (XB - 1);
                const y = YMIN + rand() * YSP * (YB - 1);
                const z = ZMIN + rand() * ZSP * (ZB - 1);
                expect(interp.evaluate(x, y, z)).toBeCloseTo(referenceTricubic(F, x, y, z, false), 11);
            }
        });

        it('has a continuous second derivative across a cell boundary', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, false);
            const xEdge = XMIN + 3 * XSP;
            const y = YMIN + 2.3 * YSP, z = ZMIN + 2.7 * ZSP;
            const eps = 1e-7;
            expect(interp.evaluate(2, 0, 0, xEdge - eps, y, z))
                .toBeCloseTo(interp.evaluate(2, 0, 0, xEdge + eps, y, z), 5);
        });
    });

    describe('partition of unity', () => {
        for (const catmullRom of [true, false]) {
            it(`reproduces a constant field everywhere (catmullRom=${catmullRom})`, () => {
                const F = sampleGrid(() => 3.25);
                const interp = makeInterp(F, catmullRom);
                const rand = makeRandom(77);
                for (let n = 0; n < 50; ++n) {
                    const x = XMIN + rand() * XSP * (XB - 1);
                    const y = YMIN + rand() * YSP * (YB - 1);
                    const z = ZMIN + rand() * ZSP * (ZB - 1);
                    expect(interp.evaluate(x, y, z)).toBeCloseTo(3.25, 11);
                }
                expect(interp.evaluate(XMIN, YMIN, ZMIN)).toBeCloseTo(3.25, 11);
                // Modestly outside the domain the blend weights still sum to
                // one per axis, so the constant survives. Far outside, the
                // cubic weights grow like the cube of the (unclamped)
                // fractional coordinate and the exact cancellation to one is
                // lost to round-off, so this check stays near the domain.
                expect(interp.evaluate(XMIN - 2 * XSP, interp.getYMax() + 1.5 * YSP,
                    ZMIN - 0.5 * ZSP)).toBeCloseTo(3.25, 9);
            });
        }
    });

    describe('derivative evaluation', () => {
        it('returns the function value when all orders are zero', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const rand = makeRandom(11);
            for (let n = 0; n < 50; ++n) {
                const [x, y, z] = interiorPoint(rand);
                expect(interp.evaluate(0, 0, 0, x, y, z)).toBe(interp.evaluate(x, y, z));
            }
        });

        it('matches the analytic gradient of a reproduced polynomial', () => {
            const F = sampleGrid(triquadratic);
            const interp = makeInterp(F, true);
            const x = 0.13, y = 2.6, z = 0.4;
            const dfdx = 0.5 - 1.5 * x + 1.1 * y - 0.6 * x * z * z + 0.7 * y * z
                + 1.8 * x * y * y * z * z;
            const dfdy = 0.25 + 0.8 * y + 1.1 * x + 0.7 * x * z + 1.8 * x * x * y * z * z;
            const dfdz = -0.6 + 0.4 * z - 0.6 * x * x * z + 0.7 * x * y
                + 1.8 * x * x * y * y * z;
            expect(interp.evaluate(1, 0, 0, x, y, z)).toBeCloseTo(dfdx, 8);
            expect(interp.evaluate(0, 1, 0, x, y, z)).toBeCloseTo(dfdy, 8);
            expect(interp.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(dfdz, 8);
        });

        it('matches the analytic second derivative of a reproduced polynomial', () => {
            const F = sampleGrid(triquadratic);
            const interp = makeInterp(F, true);
            const x = 0.13, y = 2.6, z = 0.4;
            // d2f/dx2 = -1.5 - 0.6*z^2 + 1.8*y^2*z^2
            expect(interp.evaluate(2, 0, 0, x, y, z))
                .toBeCloseTo(-1.5 - 0.6 * z * z + 1.8 * y * y * z * z, 8);
            // d2f/dy2 = 0.8 + 1.8*x^2*z^2
            expect(interp.evaluate(0, 2, 0, x, y, z))
                .toBeCloseTo(0.8 + 1.8 * x * x * z * z, 8);
            // d3f/dxdydz = 0.7 + 7.2*x*y*z, from 0.7*x*y*z and 0.9*x^2*y^2*z^2.
            expect(interp.evaluate(1, 1, 1, x, y, z))
                .toBeCloseTo(0.7 + 7.2 * x * y * z, 8);
        });

        it('matches central differences of the interpolant for the first derivatives', () => {
            const F = sampleGrid(smooth);
            for (const catmullRom of [true, false]) {
                const interp = makeInterp(F, catmullRom);
                const h = 1e-5;
                const rand = makeRandom(8);
                for (let n = 0; n < 10; ++n) {
                    const [x, y, z] = interiorPoint(rand);
                    const fdx = (interp.evaluate(x + h, y, z) - interp.evaluate(x - h, y, z)) / (2 * h);
                    const fdy = (interp.evaluate(x, y + h, z) - interp.evaluate(x, y - h, z)) / (2 * h);
                    const fdz = (interp.evaluate(x, y, z + h) - interp.evaluate(x, y, z - h)) / (2 * h);
                    expect(interp.evaluate(1, 0, 0, x, y, z)).toBeCloseTo(fdx, 7);
                    expect(interp.evaluate(0, 1, 0, x, y, z)).toBeCloseTo(fdy, 7);
                    expect(interp.evaluate(0, 0, 1, x, y, z)).toBeCloseTo(fdz, 7);
                }
            }
        });

        it('matches central differences of the first derivatives for the second', () => {
            const F = sampleGrid(smooth);
            for (const catmullRom of [true, false]) {
                const interp = makeInterp(F, catmullRom);
                const h = 1e-3;
                const rand = makeRandom(21);
                for (let n = 0; n < 10; ++n) {
                    const [x, y, z] = interiorPoint(rand);
                    const fdxx = (interp.evaluate(1, 0, 0, x + h, y, z)
                        - interp.evaluate(1, 0, 0, x - h, y, z)) / (2 * h);
                    const fdzz = (interp.evaluate(0, 0, 1, x, y, z + h)
                        - interp.evaluate(0, 0, 1, x, y, z - h)) / (2 * h);
                    expect(interp.evaluate(2, 0, 0, x, y, z)).toBeCloseTo(fdxx, 7);
                    expect(interp.evaluate(0, 0, 2, x, y, z)).toBeCloseTo(fdzz, 7);
                }
            }
        });

        it('matches central differences of the second derivatives for the third', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, false);
            const h = 1e-3;
            const [x, y, z] = interiorPoint(makeRandom(99));
            const fdxxx = (interp.evaluate(2, 0, 0, x + h, y, z)
                - interp.evaluate(2, 0, 0, x - h, y, z)) / (2 * h);
            expect(interp.evaluate(3, 0, 0, x, y, z)).toBeCloseTo(fdxxx, 6);
        });

        it('makes the third derivative constant within a cell', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const y = YMIN + 2.3 * YSP, z = ZMIN + 2.7 * ZSP;
            const base = interp.evaluate(3, 0, 0, XMIN + 2.1 * XSP, y, z);
            expect(interp.evaluate(3, 0, 0, XMIN + 2.5 * XSP, y, z)).toBeCloseTo(base, 10);
            expect(interp.evaluate(3, 0, 0, XMIN + 2.9 * XSP, y, z)).toBeCloseTo(base, 10);
        });

        it('matches central differences for the mixed partials', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const h = 1e-3;
            const [x, y, z] = interiorPoint(makeRandom(64));
            const fdxy = (interp.evaluate(1, 0, 0, x, y + h, z)
                - interp.evaluate(1, 0, 0, x, y - h, z)) / (2 * h);
            expect(interp.evaluate(1, 1, 0, x, y, z)).toBeCloseTo(fdxy, 7);

            const fdxyz = (interp.evaluate(1, 1, 0, x, y, z + h)
                - interp.evaluate(1, 1, 0, x, y, z - h)) / (2 * h);
            expect(interp.evaluate(1, 1, 1, x, y, z)).toBeCloseTo(fdxyz, 7);

            const fdxxz = (interp.evaluate(2, 0, 0, x, y, z + h)
                - interp.evaluate(2, 0, 0, x, y, z - h)) / (2 * h);
            expect(interp.evaluate(2, 0, 1, x, y, z)).toBeCloseTo(fdxxz, 7);
        });

        it('scales the derivative with the inverse spacing', () => {
            // f(x,y,z) = x on grids that differ only in the x-spacing must
            // both report df/dx = 1.
            const F = sampleGrid((x) => x);
            const coarse = makeInterp(F, true);
            expect(coarse.evaluate(1, 0, 0, 0.1, 2.6, 0.4)).toBeCloseTo(1, 10);

            const G = new Array<number>(XB * YB * ZB);
            for (let s = 0; s < ZB; ++s) {
                for (let r = 0; r < YB; ++r) {
                    for (let c = 0; c < XB; ++c) {
                        G[c + XB * (r + YB * s)] = XMIN + 0.125 * c;
                    }
                }
            }
            const fine = new IntpTricubic3(XB, YB, ZB, XMIN, 0.125, YMIN, YSP, ZMIN, ZSP, G, true);
            expect(fine.evaluate(1, 0, 0, XMIN + 0.2, 2.6, 0.4)).toBeCloseTo(1, 10);
            // The second derivative of a linear field vanishes.
            expect(fine.evaluate(2, 0, 0, XMIN + 0.2, 2.6, 0.4)).toBeCloseTo(0, 10);
        });

        it('returns zero for orders outside [0,3]', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const x = 0.13, y = 2.6, z = 0.4;
            expect(interp.evaluate(4, 0, 0, x, y, z)).toBe(0);
            expect(interp.evaluate(0, 4, 0, x, y, z)).toBe(0);
            expect(interp.evaluate(0, 0, 4, x, y, z)).toBe(0);
            expect(interp.evaluate(3, 3, 5, x, y, z)).toBe(0);
            expect(interp.evaluate(-1, 0, 0, x, y, z)).toBe(0);
        });
    });

    describe('boundary behavior', () => {
        it('evaluates the extreme corners as the corner samples (Catmull-Rom)', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            expect(interp.evaluate(XMIN, YMIN, ZMIN)).toBeCloseTo(F[0], 12);
            expect(interp.evaluate(interp.getXMax(), interp.getYMax(), interp.getZMax()))
                .toBeCloseTo(F[XB * YB * ZB - 1], 12);
            expect(interp.evaluate(interp.getXMax(), YMIN, ZMIN)).toBeCloseTo(F[XB - 1], 12);
        });

        it('agrees with the clamped reference outside the domain', () => {
            const F = sampleGrid(smooth);
            for (const catmullRom of [true, false]) {
                const interp = makeInterp(F, catmullRom);
                const outside: Array<[number, number, number]> = [
                    [XMIN - 0.4 * XSP, YMIN + 0.2, ZMIN + 0.3],
                    [XMIN - 2 * XSP, YMIN - 1.5 * YSP, ZMIN - 3 * ZSP],
                    [XMIN + XSP * (XB - 1) + 1.5 * XSP, YMIN + YSP * (YB - 1) + 0.5 * YSP,
                        ZMIN + ZSP * (ZB - 1) + 2 * ZSP]
                ];
                for (const [x, y, z] of outside) {
                    expect(interp.evaluate(x, y, z))
                        .toBeCloseTo(referenceTricubic(F, x, y, z, catmullRom), 8);
                }
            }
        });

        it('handles the minimal 4x4x4 grid', () => {
            // Every stencil is clamped on at least one side here. Catmull-Rom
            // still interpolates the samples, and both blendings reproduce a
            // constant field.
            const G = new Array<number>(64);
            for (let k = 0; k < 4; ++k) {
                for (let j = 0; j < 4; ++j) {
                    for (let i = 0; i < 4; ++i) {
                        G[i + 4 * (j + 4 * k)] = i + 10 * j + 100 * k;
                    }
                }
            }
            const catmull = new IntpTricubic3(4, 4, 4, 0, 1, 0, 1, 0, 1, G, true);
            for (let k = 0; k < 4; ++k) {
                for (let j = 0; j < 4; ++j) {
                    for (let i = 0; i < 4; ++i) {
                        expect(catmull.evaluate(i, j, k)).toBeCloseTo(G[i + 4 * (j + 4 * k)], 10);
                    }
                }
            }
            const constant = new Array<number>(64).fill(-1.5);
            for (const catmullRom of [true, false]) {
                const interp = new IntpTricubic3(4, 4, 4, 0, 1, 0, 1, 0, 1, constant, catmullRom);
                expect(interp.evaluate(1.3, 2.7, 0.4)).toBeCloseTo(-1.5, 11);
                expect(interp.evaluate(0, 0, 0)).toBeCloseTo(-1.5, 11);
            }
        });

        it('is continuous across a cell boundary near the grid edge', () => {
            const F = sampleGrid(smooth);
            const interp = makeInterp(F, true);
            const xEdge = XMIN + XSP;  // first interior lattice plane
            const y = YMIN + 0.3 * YSP, z = ZMIN + 0.7 * ZSP;
            const eps = 1e-10;
            const left = interp.evaluate(xEdge - eps, y, z);
            const right = interp.evaluate(xEdge + eps, y, z);
            expect(left).toBeCloseTo(right, 8);
            expect(left).toBeCloseTo(interp.evaluate(xEdge, y, z), 8);
        });
    });
});
