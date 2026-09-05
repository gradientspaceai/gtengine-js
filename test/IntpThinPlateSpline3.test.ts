import { describe, it, expect } from 'vitest';
import {
    IntpThinPlateSpline3, intpThinPlateSpline3Kernel
} from '../src/IntpThinPlateSpline3.js';
import { check, expectClose, fc, latticeVector, rotationFrame, seededRandom }
    from './helpers/arbitraries.js';

// ---------------------------------------------------------------------------
// An independent implementation of the classical 3D thin-plate spline, used to
// cross-check the port. The saddle-point system
//     [ A  B ] [ a ]   [ z ]
//     [ B^T 0 ] [ b ] = [ 0 ]
// with A = M + lambda*I and B(i,:) = (1, x_i, y_i, z_i) is solved by a
// hand-rolled Gaussian elimination with partial pivoting (independent of
// GMatrix and of the port's GaussianElimination).
// ---------------------------------------------------------------------------

function solveLinear(n: number, A: number[][], rhs: number[]): number[] | null {
    const M: number[][] = [];
    for (let r = 0; r < n; ++r) {
        M.push([...A[r], rhs[r]]);
    }
    for (let c = 0; c < n; ++c) {
        let piv = c;
        for (let r = c + 1; r < n; ++r) {
            if (Math.abs(M[r][c]) > Math.abs(M[piv][c])) {
                piv = r;
            }
        }
        if (Math.abs(M[piv][c]) < 1e-14) {
            return null;
        }
        const t = M[c];
        M[c] = M[piv];
        M[piv] = t;
        for (let r = 0; r < n; ++r) {
            if (r === c) {
                continue;
            }
            const f = M[r][c] / M[c][c];
            for (let k = c; k <= n; ++k) {
                M[r][k] -= f * M[c][k];
            }
        }
    }
    const x: number[] = [];
    for (let r = 0; r < n; ++r) {
        x.push(M[r][n] / M[r][r]);
    }
    return x;
}

function kernel(t: number): number {
    return -Math.abs(t);
}

function dist(ax: number, ay: number, az: number, bx: number, by: number,
    bz: number): number {
    const dx = ax - bx, dy = ay - by, dz = az - bz;
    return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

interface RefSpline {
    a: number[];
    b: number[];
    x: number[];
    y: number[];
    z: number[];
    min: number[];
    invRange: number[];
    smooth: number;
}

function refBuild(X: number[], Y: number[], Z: number[], F: number[],
    smooth: number, transform: boolean): RefSpline | null {
    const n = X.length;
    const min = [0, 0, 0];
    const invRange = [1, 1, 1];
    if (transform) {
        const src = [X, Y, Z];
        for (let d = 0; d < 3; ++d) {
            min[d] = Math.min(...src[d]);
            invRange[d] = 1 / (Math.max(...src[d]) - min[d]);
        }
    }
    const x = X.map(v => (v - min[0]) * invRange[0]);
    const y = Y.map(v => (v - min[1]) * invRange[1]);
    const z = Z.map(v => (v - min[2]) * invRange[2]);

    const N = n + 4;
    const A: number[][] = [];
    for (let r = 0; r < N; ++r) {
        A.push(new Array<number>(N).fill(0));
    }
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            A[r][c] = r === c ? smooth :
                kernel(dist(x[r], y[r], z[r], x[c], y[c], z[c]));
        }
        const B = [1, x[r], y[r], z[r]];
        for (let k = 0; k < 4; ++k) {
            A[r][n + k] = B[k];
            A[n + k][r] = B[k];
        }
    }
    const sol = solveLinear(N, A, [...F, 0, 0, 0, 0]);
    if (sol === null) {
        return null;
    }
    return {
        a: sol.slice(0, n), b: sol.slice(n), x, y, z, min, invRange, smooth
    };
}

function refEval(s: RefSpline, px: number, py: number, pz: number): number {
    const u = (px - s.min[0]) * s.invRange[0];
    const v = (py - s.min[1]) * s.invRange[1];
    const w = (pz - s.min[2]) * s.invRange[2];
    let result = s.b[0] + s.b[1] * u + s.b[2] * v + s.b[3] * w;
    for (let i = 0; i < s.a.length; ++i) {
        result += s.a[i] * kernel(dist(u, v, w, s.x[i], s.y[i], s.z[i]));
    }
    return result;
}

function refFunctional(s: RefSpline): number {
    const n = s.a.length;
    let f = 0;
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            const m = r === c ? s.smooth :
                kernel(dist(s.x[r], s.y[r], s.z[r], s.x[c], s.y[c], s.z[c]));
            f += m * s.a[r] * s.a[c];
        }
    }
    return s.smooth > 0 ? f * s.smooth : f;
}

// A small deterministic pseudorandom generator for the randomized checks.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// A scattered sample set whose bounding box is exactly the unit cube.
const X = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.35, 0.8, 0.15, 0.62];
const Y = [0.0, 0.0, 1.0, 1.0, 0.25, 0.9, 0.7, 0.45, 0.25, 0.9];
const Z = [0.0, 1.0, 0.4, 0.6, 1.0, 0.0, 0.55, 0.2, 0.85, 0.33];
const F = X.map((x, i) =>
    Math.sin(3 * x) + Math.cos(2 * Y[i]) + 0.25 * x * Z[i] - 0.5 * Z[i] * Y[i]);

describe('intpThinPlateSpline3Kernel', () => {
    it('is -|t|', () => {
        expect(intpThinPlateSpline3Kernel(0)).toBe(-0);
        expect(intpThinPlateSpline3Kernel(1)).toBe(-1);
        expect(intpThinPlateSpline3Kernel(2.5)).toBe(-2.5);
        // Unlike the 2D kernel, the 3D kernel is not clamped at t = 0, so a
        // negative argument returns -|t| (upstream never passes one).
        expect(intpThinPlateSpline3Kernel(-3)).toBe(-3);
        // It is homogeneous of degree one, which is what makes the classical
        // spline invariant to uniform scaling.
        expect(intpThinPlateSpline3Kernel(7 * 1.3)).toBeCloseTo(
            7 * intpThinPlateSpline3Kernel(1.3), 12);
    });
});

describe('IntpThinPlateSpline3 construction validation', () => {
    it('requires at least four points', () => {
        expect(() => new IntpThinPlateSpline3(3, [0, 1, 2], [0, 1, 2],
            [0, 1, 2], [0, 1, 2], 0, false)).toThrow('Invalid input.');
    });

    it('requires a nonnegative smoothing parameter', () => {
        expect(() => new IntpThinPlateSpline3(X.length, X, Y, Z, F, -1e-8,
            false)).toThrow('Invalid input.');
    });

    it('requires arrays with at least numPoints entries', () => {
        expect(() => new IntpThinPlateSpline3(X.length + 1, X, Y, Z, F, 0,
            false)).toThrow('Invalid input.');
    });

    it('accepts smooth = 0 and reports initialization', () => {
        const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, false);
        expect(tps.isInitialized()).toBe(true);
    });
});

describe('IntpThinPlateSpline3 interpolation', () => {
    for (const transform of [false, true]) {
        it(`reproduces the samples exactly when smooth = 0 (transform=${transform})`,
            () => {
                const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0,
                    transform);
                expect(tps.isInitialized()).toBe(true);
                for (let i = 0; i < X.length; ++i) {
                    expect(tps.evaluate(X[i], Y[i], Z[i])).toBeCloseTo(F[i], 10);
                }
            });
    }

    it('does not reproduce the samples exactly when smooth > 0', () => {
        const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0.05, false);
        expect(tps.isInitialized()).toBe(true);
        let maxDiff = 0;
        for (let i = 0; i < X.length; ++i) {
            maxDiff = Math.max(maxDiff,
                Math.abs(tps.evaluate(X[i], Y[i], Z[i]) - F[i]));
        }
        expect(maxDiff).toBeGreaterThan(1e-6);
    });

    it('reproduces an affine function everywhere, with and without smoothing',
        () => {
            const affine = (x: number, y: number, z: number) =>
                2.5 - 1.25 * x + 0.75 * y + 0.5 * z;
            const G = X.map((x, i) => affine(x, Y[i], Z[i]));
            for (const transform of [false, true]) {
                for (const smooth of [0, 0.1, 5]) {
                    const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, G,
                        smooth, transform);
                    expect(tps.isInitialized()).toBe(true);
                    for (const [px, py, pz] of [[0.3, 0.4, 0.5], [-2, 3.5, 7],
                        [10, -7, -1], [0, 0, 0], [1, 1, 1]]) {
                        expect(tps.evaluate(px, py, pz)).toBeCloseTo(
                            affine(px, py, pz), 8);
                    }
                    // The Green's-function coefficients vanish, so the
                    // functional is zero.
                    expect(tps.computeFunctional()).toBeCloseTo(0, 8);
                }
            }
        });
});

describe('IntpThinPlateSpline3 versus an independent direct solve', () => {
    for (const transform of [false, true]) {
        for (const smooth of [0, 0.02, 0.5]) {
            it(`matches the saddle-point solution (transform=${transform}, smooth=${smooth})`,
                () => {
                    const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F,
                        smooth, transform);
                    const ref = refBuild(X, Y, Z, F, smooth, transform);
                    expect(tps.isInitialized()).toBe(true);
                    expect(ref).not.toBeNull();
                    const rnd = makeRandom(12345);
                    for (let k = 0; k < 40; ++k) {
                        const px = -0.5 + 2 * rnd();
                        const py = -0.5 + 2 * rnd();
                        const pz = -0.5 + 2 * rnd();
                        expect(tps.evaluate(px, py, pz)).toBeCloseTo(
                            refEval(ref!, px, py, pz), 8);
                    }
                    expect(tps.computeFunctional()).toBeCloseTo(
                        refFunctional(ref!), 8);
                });
        }
    }

    it('matches the direct solve on randomized sample sets', () => {
        const rnd = makeRandom(98765);
        let checked = 0;
        for (let trial = 0; trial < 8; ++trial) {
            const n = 6 + trial;
            const RX: number[] = [], RY: number[] = [], RZ: number[] = [],
                RF: number[] = [];
            for (let i = 0; i < n; ++i) {
                RX.push(rnd() * 4 - 2);
                RY.push(rnd() * 4 - 2);
                RZ.push(rnd() * 4 - 2);
                RF.push(rnd() * 10 - 5);
            }
            const tps = new IntpThinPlateSpline3(n, RX, RY, RZ, RF, 0, true);
            const ref = refBuild(RX, RY, RZ, RF, 0, true);
            if (!tps.isInitialized() || ref === null) {
                continue;
            }
            ++checked;
            for (let k = 0; k < 10; ++k) {
                const px = rnd() * 4 - 2;
                const py = rnd() * 4 - 2;
                const pz = rnd() * 4 - 2;
                expect(tps.evaluate(px, py, pz)).toBeCloseTo(
                    refEval(ref, px, py, pz), 6);
            }
        }
        expect(checked).toBeGreaterThan(0);
    });
});

describe('IntpThinPlateSpline3 known values', () => {
    it('matches a hand-computed four-point spline', () => {
        // Four non-coplanar points make B square and invertible, so the
        // side conditions B^T*a = 0 force a = 0 and b to be the unique affine
        // function through the four samples.
        const PX = [0, 2, 0, 0];
        const PY = [0, 0, 3, 0];
        const PZ = [0, 0, 0, 4];
        // f(x,y,z) = 1 + x - y + 0.5*z
        const PF = [1, 3, -2, 3];
        const tps = new IntpThinPlateSpline3(4, PX, PY, PZ, PF, 0, false);
        expect(tps.isInitialized()).toBe(true);
        expect(tps.evaluate(0, 0, 0)).toBeCloseTo(1, 10);
        expect(tps.evaluate(2, 0, 0)).toBeCloseTo(3, 10);
        expect(tps.evaluate(0, 3, 0)).toBeCloseTo(-2, 10);
        expect(tps.evaluate(0, 0, 4)).toBeCloseTo(3, 10);
        expect(tps.evaluate(4, -1, 2)).toBeCloseTo(1 + 4 + 1 + 1, 10);
        expect(tps.evaluate(-0.5, 0.25, -2)).toBeCloseTo(1 - 0.5 - 0.25 - 1, 10);
        expect(tps.computeFunctional()).toBeCloseTo(0, 10);
    });

    it('matches a hand-computed five-point spline', () => {
        // Five samples on the axes plus the origin. The coefficients are
        // recovered from the 9x9 saddle-point system by the reference solver
        // and the spline is evaluated by the explicit thin-plate formula
        //   s(p) = b0 + b1*x + b2*y + b3*z - sum_i a_i*|p - p_i|.
        const PX = [0, 1, -1, 0, 0];
        const PY = [0, 0, 0, 1, -1];
        const PZ = [0, 0, 0, 0, 1];
        const PF = [0, 1, 2, 3, 4];
        const tps = new IntpThinPlateSpline3(5, PX, PY, PZ, PF, 0, false);
        const ref = refBuild(PX, PY, PZ, PF, 0, false);
        expect(tps.isInitialized()).toBe(true);
        expect(ref).not.toBeNull();
        for (let i = 0; i < 5; ++i) {
            expect(tps.evaluate(PX[i], PY[i], PZ[i])).toBeCloseTo(PF[i], 9);
        }
        for (const [px, py, pz] of [[0.5, 0.5, 0.5], [-1, 2, 3], [0.1, 0, 0]]) {
            expect(tps.evaluate(px, py, pz)).toBeCloseTo(
                refEval(ref!, px, py, pz), 9);
        }
    });

    it('produces a nonzero functional for non-affine data', () => {
        const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, false);
        expect(Math.abs(tps.computeFunctional())).toBeGreaterThan(1e-6);
    });
});

describe('IntpThinPlateSpline3 smoothing behavior', () => {
    it('pulls the interpolant toward the affine least-squares fit', () => {
        const probes: Array<[number, number, number]> = [
            [0.2, 0.3, 0.4], [0.5, 0.5, 0.5], [0.9, 0.1, 0.7],
            [0.4, 0.85, 0.2], [0.7, 0.7, 0.9]
        ];
        const limit = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 1e8,
            false);
        let previous = Number.POSITIVE_INFINITY;
        for (const smooth of [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e4, 1e6]) {
            const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, smooth,
                false);
            expect(tps.isInitialized()).toBe(true);
            let deviation = 0;
            for (const [px, py, pz] of probes) {
                deviation = Math.max(deviation, Math.abs(
                    tps.evaluate(px, py, pz) - limit.evaluate(px, py, pz)));
            }
            expect(deviation).toBeLessThan(previous + 1e-12);
            previous = deviation;
        }
        expect(previous).toBeLessThan(1e-4);
    });

    it('reduces the bending energy a^T*M*a as smoothing grows', () => {
        // computeFunctional() returns lambda * a^T*A*a when lambda > 0, so it
        // is not directly comparable across lambda values. Dividing out
        // lambda recovers a^T*A*a, which decreases monotonically to zero.
        let previous = Number.POSITIVE_INFINITY;
        for (const smooth of [1e-6, 1e-4, 1e-2, 1, 100, 1e4, 1e6]) {
            const tps = new IntpThinPlateSpline3(X.length, X, Y, Z, F, smooth,
                false);
            const energy = Math.abs(tps.computeFunctional() / smooth);
            expect(energy).toBeLessThan(previous);
            previous = energy;
        }
        expect(previous).toBeLessThan(1e-4);
    });
});

describe('IntpThinPlateSpline3 invariance', () => {
    it('the classical spline is invariant to translation and rotation', () => {
        const base = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, false);
        // A rotation about the axis (1,1,1)/sqrt(3) followed by a translation.
        const angle = 0.9;
        const c = Math.cos(angle), s = Math.sin(angle), t = 1 - c;
        const u = 1 / Math.sqrt(3);
        const R = [
            [t * u * u + c, t * u * u - s * u, t * u * u + s * u],
            [t * u * u + s * u, t * u * u + c, t * u * u - s * u],
            [t * u * u - s * u, t * u * u + s * u, t * u * u + c]
        ];
        const T = [3.5, -1.25, 0.75];
        const map = (x: number, y: number, z: number): [number, number, number] =>
            [R[0][0] * x + R[0][1] * y + R[0][2] * z + T[0],
            R[1][0] * x + R[1][1] * y + R[1][2] * z + T[1],
            R[2][0] * x + R[2][1] * y + R[2][2] * z + T[2]];

        const TX: number[] = [], TY: number[] = [], TZ: number[] = [];
        for (let i = 0; i < X.length; ++i) {
            const p = map(X[i], Y[i], Z[i]);
            TX.push(p[0]);
            TY.push(p[1]);
            TZ.push(p[2]);
        }
        const moved = new IntpThinPlateSpline3(X.length, TX, TY, TZ, F, 0,
            false);
        const rnd = makeRandom(24680);
        for (let k = 0; k < 30; ++k) {
            const px = rnd() * 2 - 0.5;
            const py = rnd() * 2 - 0.5;
            const pz = rnd() * 2 - 0.5;
            const q = map(px, py, pz);
            expect(moved.evaluate(q[0], q[1], q[2])).toBeCloseTo(
                base.evaluate(px, py, pz), 6);
        }
    });

    it('the classical spline IS invariant to uniform scaling (contrary to the upstream comment)',
        () => {
            // Kernel(t) = -|t| is homogeneous of degree one, so scaling the
            // samples by s scales A by s and the coefficients a by 1/s; the
            // interpolated values are unchanged. The upstream WARNING comment
            // claims the opposite.
            const base = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0,
                false);
            for (const s of [4, 0.05, 1000]) {
                const scaled = new IntpThinPlateSpline3(X.length,
                    X.map(v => s * v), Y.map(v => s * v), Z.map(v => s * v),
                    F, 0, false);
                for (const [px, py, pz] of [[0.3, 0.4, 0.5], [0.6, 0.2, 0.9],
                    [-1, 2, 0.25]]) {
                    expect(scaled.evaluate(s * px, s * py, s * pz))
                        .toBeCloseTo(base.evaluate(px, py, pz), 6);
                }
            }
        });

    it('the classical spline is not invariant to anisotropic scaling', () => {
        const base = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, false);
        const sx = 4, sy = 1, sz = 1;
        const scaled = new IntpThinPlateSpline3(X.length, X.map(v => sx * v),
            Y.map(v => sy * v), Z.map(v => sz * v), F, 0, false);
        let maxDiff = 0;
        for (const [px, py, pz] of [[0.3, 0.4, 0.5], [0.6, 0.2, 0.9],
            [0.15, 0.85, 0.35]]) {
            maxDiff = Math.max(maxDiff, Math.abs(
                scaled.evaluate(sx * px, sy * py, sz * pz) -
                base.evaluate(px, py, pz)));
        }
        expect(maxDiff).toBeGreaterThan(1e-6);
    });

    it('the unit-cube transform makes the spline invariant to translation and to per-axis scaling',
        () => {
            const base = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0,
                true);
            const s = [7, 0.125, 2.5];
            const t = [-3, 11, 0.5];
            const TX = X.map(v => s[0] * v + t[0]);
            const TY = Y.map(v => s[1] * v + t[1]);
            const TZ = Z.map(v => s[2] * v + t[2]);
            const moved = new IntpThinPlateSpline3(X.length, TX, TY, TZ, F, 0,
                true);
            const rnd = makeRandom(13579);
            for (let k = 0; k < 30; ++k) {
                const px = rnd() * 1.5 - 0.25;
                const py = rnd() * 1.5 - 0.25;
                const pz = rnd() * 1.5 - 0.25;
                expect(moved.evaluate(s[0] * px + t[0], s[1] * py + t[1],
                    s[2] * pz + t[2])).toBeCloseTo(base.evaluate(px, py, pz), 8);
            }
        });

    it('the transform is the identity when the data already fills the unit cube',
        () => {
            const a = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, false);
            const b = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0, true);
            // min/max of X, Y and Z are 0 and 1 for the sample set, so the
            // transformed and untransformed splines agree.
            const rnd = makeRandom(11111);
            for (let k = 0; k < 20; ++k) {
                const px = rnd(), py = rnd(), pz = rnd();
                expect(b.evaluate(px, py, pz)).toBeCloseTo(
                    a.evaluate(px, py, pz), 9);
            }
            expect(b.computeFunctional()).toBeCloseTo(a.computeFunctional(), 9);
        });
});

describe('IntpThinPlateSpline3 degenerate inputs', () => {
    it('fails to initialize when two samples coincide', () => {
        // Two identical points make two rows of A identical, so A is exactly
        // singular and the constructor returns before setting mInitialized.
        const DX = [...X, X[2]];
        const DY = [...Y, Y[2]];
        const DZ = [...Z, Z[2]];
        const DF = [...F, 9];
        const tps = new IntpThinPlateSpline3(DX.length, DX, DY, DZ, DF, 0,
            false);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(0.5, 0.5, 0.5)).toBe(Number.MAX_VALUE);
        // The coefficients are all zero, so the functional is zero.
        expect(tps.computeFunctional()).toBe(0);
    });

    it('is singular in exact arithmetic for a coplanar sample set', () => {
        // Coplanar samples make B rank deficient, so Q = B^T*A^{-1}*B is
        // singular in exact arithmetic. The floating-point inverse does not
        // detect the rank deficiency for this configuration, so the spline
        // reports itself as initialized; the affine coefficients are then
        // ill-determined but the samples are still reproduced (the residual
        // of the linear solve is what the spline actually fits).
        const n = X.length;
        const CZ = X.map((x, i) => 0.5 * x + 0.25 * Y[i] + 1);
        const tps = new IntpThinPlateSpline3(n, X, Y, CZ, F, 0, false);
        expect(tps.isInitialized()).toBe(true);
        for (let i = 0; i < n; ++i) {
            expect(Number.isFinite(tps.evaluate(X[i], Y[i], CZ[i]))).toBe(true);
            expect(tps.evaluate(X[i], Y[i], CZ[i])).toBeCloseTo(F[i], 4);
        }
    });

    it('fails to initialize for a collinear sample set', () => {
        const CX = [0, 1, 2, 3, 4];
        const CY = [0, 1, 2, 3, 4];
        const CZ = [0, 1, 2, 3, 4];
        const CF = [0, 1, 4, 9, 16];
        const tps = new IntpThinPlateSpline3(5, CX, CY, CZ, CF, 0, false);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(1, 1, 1)).toBe(Number.MAX_VALUE);
    });

    it('fails to initialize when a transformed axis is degenerate', () => {
        // All samples share the same z, so mZInvRange is infinite and the
        // transform maps every z to NaN (0 * Infinity). The NaN propagates
        // into A and the inverse is reported as not invertible. Upstream has
        // the same behavior; the port preserves it.
        const CZ = X.map(() => 5);
        const tps = new IntpThinPlateSpline3(X.length, X, Y, CZ, F, 0, true);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(0.5, 0.5, 5)).toBe(Number.MAX_VALUE);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpThinPlateSpline3.h, cross-checked with the independent saddle-point
// solve above.
// ---------------------------------------------------------------------------

// Scattered sample sets on a lattice: distinct points that are not coplanar,
// so that B has full rank and Q is invertible. Integer coordinates and
// integer samples keep the input exactly representable.
const tps3Samples = fc.tuple(
    fc.array(latticeVector(3, -5, 5), { minLength: 6, maxLength: 9 }),
    fc.array(fc.integer({ min: -9, max: 9 }), { minLength: 9, maxLength: 9 })
).filter(([pts]) => {
    for (let i = 0; i < pts.length; ++i) {
        for (let j = i + 1; j < pts.length; ++j) {
            if (pts[i].values.every((v, k) => v === pts[j].values[k])) {
                return false;
            }
        }
    }
    // Reject coplanar sets: pick three points spanning a plane through pts[0]
    // and require a fourth off it.
    const p0 = pts[0].values;
    for (let i = 1; i < pts.length; ++i) {
        for (let j = i + 1; j < pts.length; ++j) {
            const u = pts[i].values.map((v, k) => v - p0[k]);
            const v = pts[j].values.map((w, k) => w - p0[k]);
            const nrm = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0]];
            if (nrm.every(t => t === 0)) {
                continue;
            }
            for (let k = j + 1; k < pts.length; ++k) {
                const w = pts[k].values.map((t, m) => t - p0[m]);
                if (Math.abs(nrm[0] * w[0] + nrm[1] * w[1] + nrm[2] * w[2])
                    > 0.5) {
                    return true;
                }
            }
        }
    }
    return false;
}).map(([pts, samples]) => ({
    X: pts.map(p => p.values[0]),
    Y: pts.map(p => p.values[1]),
    Z: pts.map(p => p.values[2]),
    F: pts.map((_, i) => samples[i % samples.length])
}));

describe('IntpThinPlateSpline3 verification', () => {
    // The port and the independent saddle-point solve must produce the same
    // interpolant; the tolerance is relative because the two solvers round
    // differently on a moderately ill-conditioned kernel matrix.
    it('agrees with the independent saddle-point solve', () => {
        let numChecked = 0;
        check(fc.tuple(tps3Samples, fc.boolean(),
            fc.constantFrom(0, 0.05, 0.5, 5),
            fc.integer({ min: 0, max: 0xffff })),
            ([{ X, Y, Z, F }, transform, smooth, seed]) => {
                const spline = new IntpThinPlateSpline3(X.length, X, Y, Z, F,
                    smooth, transform);
                const ref = refBuild(X, Y, Z, F, smooth, transform);
                if (!spline.isInitialized() || ref === null) {
                    return true;
                }
                ++numChecked;
                const rand = seededRandom(seed + 1);
                for (let n = 0; n < 6; ++n) {
                    const px = -5 + 10 * rand();
                    const py = -5 + 10 * rand();
                    const pz = -5 + 10 * rand();
                    expectClose(spline.evaluate(px, py, pz),
                        refEval(ref, px, py, pz), 1e-6, 1e-6);
                }
                expectClose(spline.computeFunctional(), refFunctional(ref),
                    1e-6, 1e-6);
                return true;
            }, 50);
        expect(numChecked).toBeGreaterThan(10);
    }, 30000);

    // With lambda = 0 the diagonal of A is zero and the system enforces the
    // interpolation conditions exactly.
    it('interpolates the samples exactly when smooth = 0', () => {
        check(fc.tuple(tps3Samples, fc.boolean()),
            ([{ X, Y, Z, F }, transform]) => {
                const spline = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0,
                    transform);
                if (!spline.isInitialized()) {
                    return true;
                }
                for (let i = 0; i < X.length; ++i) {
                    expectClose(spline.evaluate(X[i], Y[i], Z[i]), F[i],
                        1e-7, 1e-8);
                }
                return true;
            }, 50);
    }, 30000);

    // Affine data is reproduced for every lambda: the affine part of the
    // spline fits it with a = 0.
    it('reproduces affine data for every smoothing parameter', () => {
        check(fc.tuple(tps3Samples, fc.boolean(),
            fc.constantFrom(0, 0.25, 2),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 0, max: 0xffff })),
            ([{ X, Y, Z }, transform, smooth, a, b, c, d, seed]) => {
                const F = X.map((x, i) => a + b * x + c * Y[i] + d * Z[i]);
                const spline = new IntpThinPlateSpline3(X.length, X, Y, Z, F,
                    smooth, transform);
                if (!spline.isInitialized()) {
                    return true;
                }
                const rand = seededRandom(seed + 5);
                for (let n = 0; n < 4; ++n) {
                    const px = -5 + 10 * rand();
                    const py = -5 + 10 * rand();
                    const pz = -5 + 10 * rand();
                    expectClose(spline.evaluate(px, py, pz),
                        a + b * px + c * py + d * pz, 1e-6, 1e-7);
                }
                expectClose(spline.computeFunctional(), 0, 1e-6, 1e-6);
                return true;
            }, 50);
    }, 30000);

    // The classical (untransformed) spline is built only from the pairwise
    // distances and the affine basis, so it commutes with a rigid motion of
    // space. This is the half of the upstream WARNING comment that is correct
    // (issue #191).
    it('the classical spline is invariant to rigid motions', () => {
        check(fc.tuple(tps3Samples, rotationFrame(3),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: 0, max: 0xffff })),
            ([{ X, Y, Z, F }, frame, tx, ty, tz, seed]) => {
                const base = new IntpThinPlateSpline3(X.length, X, Y, Z, F, 0,
                    false);
                if (!base.isInitialized()) {
                    return true;
                }
                const move = (x: number, y: number, z: number) => [
                    frame[0].values[0] * x + frame[1].values[0] * y
                        + frame[2].values[0] * z + tx,
                    frame[0].values[1] * x + frame[1].values[1] * y
                        + frame[2].values[1] * z + ty,
                    frame[0].values[2] * x + frame[1].values[2] * y
                        + frame[2].values[2] * z + tz];
                const moved3 = X.map((x, i) => move(x, Y[i], Z[i]));
                const moved = new IntpThinPlateSpline3(X.length,
                    moved3.map(p => p[0]), moved3.map(p => p[1]),
                    moved3.map(p => p[2]), F, 0, false);
                if (!moved.isInitialized()) {
                    return true;
                }
                const rand = seededRandom(seed + 7);
                for (let n = 0; n < 4; ++n) {
                    const px = -5 + 10 * rand();
                    const py = -5 + 10 * rand();
                    const pz = -5 + 10 * rand();
                    const q = move(px, py, pz);
                    expectClose(moved.evaluate(q[0], q[1], q[2]),
                        base.evaluate(px, py, pz), 1e-5, 1e-6);
                }
                return true;
            }, 50);
    }, 30000);

    // Smoothing trades fidelity for bending energy: as lambda grows the
    // interpolant departs further from the samples.
    it('increases the sample residual monotonically with the smoothing parameter',
        () => {
            check(tps3Samples, ({ X, Y, Z, F }) => {
                let previous = -1;
                for (const lambda of [0, 0.01, 0.1, 1, 10]) {
                    const spline = new IntpThinPlateSpline3(X.length, X, Y, Z, F,
                        lambda, false);
                    if (!spline.isInitialized()) {
                        return true;
                    }
                    let residual = 0;
                    for (let i = 0; i < X.length; ++i) {
                        residual += (spline.evaluate(X[i], Y[i], Z[i]) - F[i]) ** 2;
                    }
                    expect(residual).toBeGreaterThan(previous - 1e-6);
                    previous = residual;
                }
                return true;
            }, 40);
        }, 30000);

    // Upstream (issue #191): when the data has zero range along an axis, the
    // unit-cube transform divides by zero. The port preserves that.
    it('produces a non-finite transform for a zero-range axis', () => {
        const Xc = [1, 1, 1, 1, 1];
        const Yc = [0, 1, 2, 3, 1];
        const Zc = [0, 1, 0, 1, 2];
        const Fc = [0, 1, 4, 9, 3];
        const spline = new IntpThinPlateSpline3(5, Xc, Yc, Zc, Fc, 0, true);
        expect(spline.isInitialized()).toBe(false);
        expect(spline.evaluate(1, 1, 1)).toBe(Number.MAX_VALUE);
    });
});
