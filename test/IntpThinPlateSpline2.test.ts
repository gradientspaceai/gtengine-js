import { describe, it, expect } from 'vitest';
import {
    IntpThinPlateSpline2, intpThinPlateSpline2Kernel
} from '../src/IntpThinPlateSpline2';

// ---------------------------------------------------------------------------
// An independent implementation of the classical thin-plate spline, used to
// cross-check the port. The saddle-point system
//     [ A  B ] [ a ]   [ z ]
//     [ B^T 0 ] [ b ] = [ 0 ]
// with A = M + lambda*I and B(i,:) = (1, x_i, y_i) is solved by a hand-rolled
// Gaussian elimination with partial pivoting (independent of GMatrix and of
// the port's GaussianElimination).
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
    return t > 0 ? t * t * Math.log(t * t) : 0;
}

interface RefSpline {
    a: number[];
    b: number[];
    x: number[];
    y: number[];
    xMin: number;
    xInvRange: number;
    yMin: number;
    yInvRange: number;
    smooth: number;
}

function refBuild(X: number[], Y: number[], F: number[], smooth: number,
    transform: boolean): RefSpline | null {
    const n = X.length;
    let xMin = 0, xInvRange = 1, yMin = 0, yInvRange = 1;
    if (transform) {
        xMin = Math.min(...X);
        yMin = Math.min(...Y);
        xInvRange = 1 / (Math.max(...X) - xMin);
        yInvRange = 1 / (Math.max(...Y) - yMin);
    }
    const x = X.map(v => (v - xMin) * xInvRange);
    const y = Y.map(v => (v - yMin) * yInvRange);

    const N = n + 3;
    const A: number[][] = [];
    for (let r = 0; r < N; ++r) {
        A.push(new Array<number>(N).fill(0));
    }
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            A[r][c] = r === c ? smooth :
                kernel(Math.hypot(x[r] - x[c], y[r] - y[c]));
        }
        const B = [1, x[r], y[r]];
        for (let k = 0; k < 3; ++k) {
            A[r][n + k] = B[k];
            A[n + k][r] = B[k];
        }
    }
    const rhs = [...F, 0, 0, 0];
    const sol = solveLinear(N, A, rhs);
    if (sol === null) {
        return null;
    }
    return {
        a: sol.slice(0, n), b: sol.slice(n), x, y,
        xMin, xInvRange, yMin, yInvRange, smooth
    };
}

function refEval(s: RefSpline, px: number, py: number): number {
    const u = (px - s.xMin) * s.xInvRange;
    const v = (py - s.yMin) * s.yInvRange;
    let result = s.b[0] + s.b[1] * u + s.b[2] * v;
    for (let i = 0; i < s.a.length; ++i) {
        result += s.a[i] * kernel(Math.hypot(u - s.x[i], v - s.y[i]));
    }
    return result;
}

function refFunctional(s: RefSpline): number {
    const n = s.a.length;
    let f = 0;
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            const m = r === c ? s.smooth :
                kernel(Math.hypot(s.x[r] - s.x[c], s.y[r] - s.y[c]));
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

// A well-conditioned scattered sample set used by most tests.
const X = [0.0, 1.0, 0.0, 1.0, 0.35, 0.8, 0.15, 0.62];
const Y = [0.0, 0.0, 1.0, 1.0, 0.7, 0.45, 0.25, 0.9];
const F = X.map((x, i) => Math.sin(3 * x) + Math.cos(2 * Y[i]) + 0.25 * x * Y[i]);

describe('intpThinPlateSpline2Kernel', () => {
    it('is t^2 * log(t^2) for positive t and 0 at t = 0', () => {
        expect(intpThinPlateSpline2Kernel(0)).toBe(0);
        // Kernel(1) = 1 * log(1) = 0.
        expect(intpThinPlateSpline2Kernel(1)).toBe(0);
        // Kernel(e) = e^2 * log(e^2) = 2 e^2.
        expect(intpThinPlateSpline2Kernel(Math.E)).toBeCloseTo(
            2 * Math.E * Math.E, 12);
        // Kernel(t) < 0 for 0 < t < 1.
        expect(intpThinPlateSpline2Kernel(0.5)).toBeCloseTo(
            0.25 * Math.log(0.25), 15);
        expect(intpThinPlateSpline2Kernel(0.5)).toBeLessThan(0);
        // Negative t is treated as t <= 0 by upstream (distances are never
        // negative, so the branch is unreachable in practice).
        expect(intpThinPlateSpline2Kernel(-2)).toBe(0);
    });
});

describe('IntpThinPlateSpline2 construction validation', () => {
    it('requires at least three points', () => {
        expect(() => new IntpThinPlateSpline2(2, [0, 1], [0, 1], [0, 1],
            0, false)).toThrow('Invalid input.');
    });

    it('requires a nonnegative smoothing parameter', () => {
        expect(() => new IntpThinPlateSpline2(3, X, Y, F, -1e-8, false))
            .toThrow('Invalid input.');
    });

    it('requires arrays with at least numPoints entries', () => {
        expect(() => new IntpThinPlateSpline2(9, X, Y, F, 0, false))
            .toThrow('Invalid input.');
    });

    it('accepts smooth = 0 and reports initialization', () => {
        const tps = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
        expect(tps.isInitialized()).toBe(true);
    });
});

describe('IntpThinPlateSpline2 interpolation', () => {
    for (const transform of [false, true]) {
        it(`reproduces the samples exactly when smooth = 0 (transform=${transform})`,
            () => {
                const tps = new IntpThinPlateSpline2(X.length, X, Y, F, 0,
                    transform);
                expect(tps.isInitialized()).toBe(true);
                for (let i = 0; i < X.length; ++i) {
                    expect(tps.evaluate(X[i], Y[i])).toBeCloseTo(F[i], 10);
                }
            });
    }

    it('does not reproduce the samples exactly when smooth > 0', () => {
        const tps = new IntpThinPlateSpline2(X.length, X, Y, F, 0.05, false);
        expect(tps.isInitialized()).toBe(true);
        let maxDiff = 0;
        for (let i = 0; i < X.length; ++i) {
            maxDiff = Math.max(maxDiff,
                Math.abs(tps.evaluate(X[i], Y[i]) - F[i]));
        }
        expect(maxDiff).toBeGreaterThan(1e-6);
    });

    it('reproduces an affine function everywhere, with and without smoothing',
        () => {
            const affine = (x: number, y: number) => 2.5 - 1.25 * x + 0.75 * y;
            const G = X.map((x, i) => affine(x, Y[i]));
            for (const transform of [false, true]) {
                for (const smooth of [0, 0.1, 5]) {
                    const tps = new IntpThinPlateSpline2(X.length, X, Y, G,
                        smooth, transform);
                    expect(tps.isInitialized()).toBe(true);
                    for (const [px, py] of [[0.3, 0.4], [-2, 3.5], [10, -7],
                        [0, 0], [1, 1]]) {
                        expect(tps.evaluate(px, py)).toBeCloseTo(
                            affine(px, py), 8);
                    }
                    // The Green's-function coefficients vanish, so the
                    // functional is zero.
                    expect(tps.computeFunctional()).toBeCloseTo(0, 8);
                }
            }
        });
});

describe('IntpThinPlateSpline2 versus an independent direct solve', () => {
    for (const transform of [false, true]) {
        for (const smooth of [0, 0.02, 0.5]) {
            it(`matches the saddle-point solution (transform=${transform}, smooth=${smooth})`,
                () => {
                    const tps = new IntpThinPlateSpline2(X.length, X, Y, F,
                        smooth, transform);
                    const ref = refBuild(X, Y, F, smooth, transform);
                    expect(tps.isInitialized()).toBe(true);
                    expect(ref).not.toBeNull();
                    const rnd = makeRandom(12345);
                    for (let k = 0; k < 40; ++k) {
                        const px = -0.5 + 2 * rnd();
                        const py = -0.5 + 2 * rnd();
                        expect(tps.evaluate(px, py)).toBeCloseTo(
                            refEval(ref!, px, py), 8);
                    }
                    expect(tps.computeFunctional()).toBeCloseTo(
                        refFunctional(ref!), 8);
                });
        }
    }

    it('matches the direct solve on randomized sample sets', () => {
        const rnd = makeRandom(98765);
        for (let trial = 0; trial < 8; ++trial) {
            const n = 5 + trial;
            const RX: number[] = [], RY: number[] = [], RF: number[] = [];
            for (let i = 0; i < n; ++i) {
                RX.push(rnd() * 4 - 2);
                RY.push(rnd() * 4 - 2);
                RF.push(rnd() * 10 - 5);
            }
            const tps = new IntpThinPlateSpline2(n, RX, RY, RF, 0, true);
            const ref = refBuild(RX, RY, RF, 0, true);
            if (!tps.isInitialized() || ref === null) {
                continue;
            }
            for (let k = 0; k < 10; ++k) {
                const px = rnd() * 4 - 2;
                const py = rnd() * 4 - 2;
                const expected = refEval(ref, px, py);
                expect(tps.evaluate(px, py)).toBeCloseTo(expected, 6);
            }
        }
    });
});

describe('IntpThinPlateSpline2 known values', () => {
    it('matches a hand-computed three-point spline', () => {
        // Three non-collinear points make B square and invertible, so the
        // saddle system B^T*a = 0 forces a = 0 and b to be the unique plane
        // through the three samples. Hence the spline is exactly that plane.
        const PX = [0, 2, 0];
        const PY = [0, 0, 3];
        const PF = [1, 3, -2];
        const tps = new IntpThinPlateSpline2(3, PX, PY, PF, 0, false);
        expect(tps.isInitialized()).toBe(true);
        // The plane is f(x,y) = 1 + x - y.
        expect(tps.evaluate(0, 0)).toBeCloseTo(1, 10);
        expect(tps.evaluate(2, 0)).toBeCloseTo(3, 10);
        expect(tps.evaluate(0, 3)).toBeCloseTo(-2, 10);
        expect(tps.evaluate(4, -1)).toBeCloseTo(6, 10);
        expect(tps.evaluate(-0.5, 0.25)).toBeCloseTo(0.25, 10);
        expect(tps.computeFunctional()).toBeCloseTo(0, 10);
    });

    it('produces a nonzero functional for non-affine data', () => {
        const tps = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
        expect(Math.abs(tps.computeFunctional())).toBeGreaterThan(1e-6);
    });
});

describe('IntpThinPlateSpline2 smoothing behavior', () => {
    it('pulls the interpolant toward the affine least-squares fit', () => {
        // Deviation from the (fixed) affine trend of the data shrinks as the
        // smoothing parameter grows.
        const probes: Array<[number, number]> = [
            [0.2, 0.3], [0.5, 0.5], [0.9, 0.1], [0.4, 0.85], [0.7, 0.7]
        ];
        let previous = Number.POSITIVE_INFINITY;
        for (const smooth of [1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e4,
            1e6]) {
            const tps = new IntpThinPlateSpline2(X.length, X, Y, F, smooth,
                false);
            expect(tps.isInitialized()).toBe(true);
            // Measure the non-affine content by the second difference of the
            // interpolant along a probe set, comparing to the large-smoothing
            // limit (an affine function).
            const limit = new IntpThinPlateSpline2(X.length, X, Y, F, 1e8,
                false);
            let deviation = 0;
            for (const [px, py] of probes) {
                deviation = Math.max(deviation,
                    Math.abs(tps.evaluate(px, py) - limit.evaluate(px, py)));
            }
            expect(deviation).toBeLessThan(previous + 1e-12);
            previous = deviation;
        }
        expect(previous).toBeLessThan(1e-4);
    });

    it('reduces the bending energy a^T*M*a as smoothing grows', () => {
        // computeFunctional() returns lambda * a^T*A*a when lambda > 0, so it
        // is not directly comparable across lambda values. The underlying
        // a^T*A*a is recovered by dividing out lambda, and that quantity
        // decreases monotonically toward zero.
        let previous = Number.POSITIVE_INFINITY;
        for (const smooth of [1e-6, 1e-4, 1e-2, 1, 100, 1e4, 1e6]) {
            const tps = new IntpThinPlateSpline2(X.length, X, Y, F, smooth,
                false);
            const energy = Math.abs(tps.computeFunctional() / smooth);
            expect(energy).toBeLessThan(previous);
            previous = energy;
        }
        expect(previous).toBeLessThan(1e-4);
    });
});

describe('IntpThinPlateSpline2 invariance', () => {
    it('the classical spline is invariant to translation and rotation', () => {
        const base = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
        const angle = 0.7;
        const cs = Math.cos(angle), sn = Math.sin(angle);
        const tx = 3.5, ty = -1.25;
        const map = (x: number, y: number): [number, number] =>
            [cs * x - sn * y + tx, sn * x + cs * y + ty];
        const TX: number[] = [], TY: number[] = [];
        for (let i = 0; i < X.length; ++i) {
            const [u, v] = map(X[i], Y[i]);
            TX.push(u);
            TY.push(v);
        }
        const moved = new IntpThinPlateSpline2(X.length, TX, TY, F, 0, false);
        const rnd = makeRandom(24680);
        for (let k = 0; k < 30; ++k) {
            const px = rnd() * 2 - 0.5;
            const py = rnd() * 2 - 0.5;
            const [u, v] = map(px, py);
            expect(moved.evaluate(u, v)).toBeCloseTo(base.evaluate(px, py), 6);
        }
    });

    it('the classical spline IS invariant to uniform scaling (contrary to the upstream comment)',
        () => {
            // Kernel(s*r) = s^2*Kernel(r) + 2*s^2*log(s)*r^2. Under the side
            // conditions B^T*a = 0 the extra r^2 term contributes only a
            // constant, which the affine part absorbs, so the interpolant is
            // unchanged. The upstream WARNING comment claims the opposite.
            const base = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
            for (const s of [4, 0.05, 1000]) {
                const scaled = new IntpThinPlateSpline2(X.length,
                    X.map(v => s * v), Y.map(v => s * v), F, 0, false);
                for (const [px, py] of [[0.3, 0.4], [0.6, 0.2], [0.15, 0.85],
                    [-1, 2]]) {
                    expect(scaled.evaluate(s * px, s * py)).toBeCloseTo(
                        base.evaluate(px, py), 6);
                }
            }
        });

    it('the classical spline is not invariant to anisotropic scaling', () => {
        const base = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
        const sx = 4, sy = 1;
        const scaled = new IntpThinPlateSpline2(X.length, X.map(v => sx * v),
            Y.map(v => sy * v), F, 0, false);
        let maxDiff = 0;
        for (const [px, py] of [[0.3, 0.4], [0.6, 0.2], [0.15, 0.85]]) {
            maxDiff = Math.max(maxDiff, Math.abs(
                scaled.evaluate(sx * px, sy * py) - base.evaluate(px, py)));
        }
        expect(maxDiff).toBeGreaterThan(1e-6);
    });

    it('the unit-square transform is not invariant to rotation (contrary to the upstream comment)',
        () => {
            // The transform rescales x and y independently, which is not a
            // similarity, so rotating the samples changes the interpolant.
            // The sample set must not have a square bounding box, otherwise
            // the two axis scales coincide and the transform happens to be a
            // similarity.
            const AX = X;
            const AY = Y.map(v => 3 * v);
            const base = new IntpThinPlateSpline2(AX.length, AX, AY, F, 0,
                true);
            const angle = 0.7;
            const cs = Math.cos(angle), sn = Math.sin(angle);
            const RX: number[] = [], RY: number[] = [];
            for (let i = 0; i < AX.length; ++i) {
                RX.push(cs * AX[i] - sn * AY[i]);
                RY.push(sn * AX[i] + cs * AY[i]);
            }
            const rotated = new IntpThinPlateSpline2(AX.length, RX, RY, F, 0,
                true);
            let maxDiff = 0;
            for (const [px, py] of [[0.3, 0.4], [0.6, 0.2], [0.15, 0.85]]) {
                maxDiff = Math.max(maxDiff, Math.abs(
                    rotated.evaluate(cs * px - sn * py, sn * px + cs * py) -
                    base.evaluate(px, py)));
            }
            expect(maxDiff).toBeGreaterThan(1e-6);
    });

    it('the unit-square transform makes the spline invariant to translation and to per-axis scaling',
        () => {
            const base = new IntpThinPlateSpline2(X.length, X, Y, F, 0, true);
            const sx = 7, sy = 0.125, tx = -3, ty = 11;
            const TX = X.map(v => sx * v + tx);
            const TY = Y.map(v => sy * v + ty);
            const moved = new IntpThinPlateSpline2(X.length, TX, TY, F, 0,
                true);
            const rnd = makeRandom(13579);
            for (let k = 0; k < 30; ++k) {
                const px = rnd() * 1.5 - 0.25;
                const py = rnd() * 1.5 - 0.25;
                expect(moved.evaluate(sx * px + tx, sy * py + ty))
                    .toBeCloseTo(base.evaluate(px, py), 8);
            }
        });

    it('the transform is the identity when the data already fills the unit square',
        () => {
            // min/max of X and of Y are 0 and 1 for the sample set, so the
            // transformed and untransformed splines agree.
            const a = new IntpThinPlateSpline2(X.length, X, Y, F, 0, false);
            const b = new IntpThinPlateSpline2(X.length, X, Y, F, 0, true);
            const rnd = makeRandom(11111);
            for (let k = 0; k < 20; ++k) {
                const px = rnd(), py = rnd();
                expect(b.evaluate(px, py)).toBeCloseTo(a.evaluate(px, py), 9);
            }
            expect(b.computeFunctional()).toBeCloseTo(a.computeFunctional(), 9);
        });
});

describe('IntpThinPlateSpline2 degenerate inputs', () => {
    it('fails to initialize for the unit right triangle', () => {
        // The legs have length 1 and Kernel(1) = 1*log(1) = 0, so the row of
        // A for the right-angle vertex is entirely zero and A is singular.
        const tps = new IntpThinPlateSpline2(3, [0, 1, 0], [0, 0, 1],
            [1, 3, -2], 0, false);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(0, 0)).toBe(Number.MAX_VALUE);
    });

    it('fails to initialize when two samples coincide', () => {
        // Two identical points make two rows of A identical, so A is exactly
        // singular and the constructor returns before setting mInitialized.
        const DX = [0, 1, 0, 1, 0.5, 0];
        const DY = [0, 0, 1, 1, 0.5, 1];
        const DF = [1, 2, 3, 4, 5, 9];
        const tps = new IntpThinPlateSpline2(DX.length, DX, DY, DF, 0, false);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(0.5, 0.5)).toBe(Number.MAX_VALUE);
        // The coefficients are all zero, so the functional is zero.
        expect(tps.computeFunctional()).toBe(0);
    });

    it('fails to initialize for a collinear sample set', () => {
        // Collinear samples make B rank deficient, so Q = B^T*A^{-1}*B is
        // singular and the constructor returns early. (Whether the
        // floating-point inverse detects the rank deficiency depends on the
        // configuration; this one is detected.)
        const CX = [0, 1, 2, 3, 4];
        const CY = [0, 1, 2, 3, 4];
        const CF = [0, 1, 4, 9, 16];
        const tps = new IntpThinPlateSpline2(5, CX, CY, CF, 0, false);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(1, 1)).toBe(Number.MAX_VALUE);
        expect(tps.evaluate(-100, 250)).toBe(Number.MAX_VALUE);
    });

    it('fails to initialize when a transformed axis is degenerate', () => {
        // All samples share the same y, so mYInvRange is infinite and the
        // transform maps every y to NaN (0 * Infinity). The NaN propagates
        // into A and the inverse is reported as not invertible. Upstream has
        // the same behavior; the port preserves it.
        const CX = [0, 1, 2, 3];
        const CY = [5, 5, 5, 5];
        const CF = [1, 2, 3, 4];
        const tps = new IntpThinPlateSpline2(4, CX, CY, CF, 0, true);
        expect(tps.isInitialized()).toBe(false);
        expect(tps.evaluate(1, 5)).toBe(Number.MAX_VALUE);
    });
});
