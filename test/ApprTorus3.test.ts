import { describe, it, expect } from 'vitest';
import { ApprTorus3 } from '../src/ApprTorus3.js';
import type { ApprTorus3Parameters } from '../src/ApprTorus3.js';
import { Vector, dot, normalize } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, orthonormalFrame,
    seededRandom, unitVector, vector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

// Points on the torus with center C, unit-length axis N, center radius r0 and
// tube radius r1:
//   X = C + (r0 + r1*cos(phi)) * (cos(theta)*D0 + sin(theta)*D1)
//       + r1*sin(phi) * N
function torusPoints(C: Vector, N: Vector, r0: number, r1: number,
    numTheta: number, numPhi: number, noise?: () => number): Vector[] {
    const basis = [N.clone(), new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    const D0 = basis[1], D1 = basis[2];
    const points: Vector[] = [];
    for (let i = 0; i < numTheta; ++i) {
        const theta = 2 * Math.PI * i / numTheta;
        const csT = Math.cos(theta), snT = Math.sin(theta);
        for (let j = 0; j < numPhi; ++j) {
            const phi = 2 * Math.PI * j / numPhi;
            const radial = r0 + r1 * Math.cos(phi);
            const height = r1 * Math.sin(phi);
            const p = new Vector(3);
            for (let k = 0; k < 3; ++k) {
                p.values[k] = C.values[k]
                    + radial * (csT * D0.values[k] + snT * D1.values[k])
                    + height * N.values[k]
                    + (noise ? noise() : 0);
            }
            points.push(p);
        }
    }
    return points;
}

// The implicit torus function; it vanishes for points on the torus.
function torusResidual(torus: ApprTorus3Parameters, X: Vector): number {
    const D = new Vector(3);
    for (let k = 0; k < 3; ++k) {
        D.values[k] = X.values[k] - torus.C.values[k];
    }
    const L = dot(D, D);
    const h = dot(torus.N, D);
    const u = torus.r0 * torus.r0;
    const v = u - torus.r1 * torus.r1;
    const sum = L + v;
    return sum * sum - 4 * u * (L - h * h);
}

function newParameters(): ApprTorus3Parameters {
    return { C: new Vector(3), N: new Vector(3), r0: 0, r1: 0 };
}

const truthC = v3(1, -2, 0.5);
const truthN = (() => {
    const n = v3(1, 2, 3);
    normalize(n);
    return n;
})();
const truthR0 = 3;
const truthR1 = 1;

describe('ApprTorus3.compute', () => {
    it('recovers a torus from points that lie exactly on it', () => {
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16);
        const torus = newParameters();
        const fitter = new ApprTorus3();
        const result = fitter.compute(points, torus);

        expect(result.success).toBe(true);
        expect(result.error).toBeLessThan(1e-12);
        for (let k = 0; k < 3; ++k) {
            expect(torus.C.values[k]).toBeCloseTo(truthC.values[k], 8);
        }
        expect(Math.abs(dot(torus.N, truthN))).toBeCloseTo(1, 8);
        expect(torus.r0).toBeCloseTo(truthR0, 8);
        expect(torus.r1).toBeCloseTo(truthR1, 8);
        for (const p of points) {
            expect(Math.abs(torusResidual(torus, p))).toBeLessThan(1e-6);
        }
    });

    it('recovers random tori (axis-aligned frames excluded)', () => {
        const rand = makeRandom(60606);
        for (let trial = 0; trial < 5; ++trial) {
            const C = v3(-1 + 2 * rand(), -1 + 2 * rand(), -1 + 2 * rand());
            const N = v3(-1 + 2 * rand(), -1 + 2 * rand(), 0.5 + rand());
            normalize(N);
            const r0 = 2 + 2 * rand();
            const r1 = 0.25 + 0.75 * rand();
            const points = torusPoints(C, N, r0, r1, 20, 12);

            const torus = newParameters();
            const fitter = new ApprTorus3();
            const result = fitter.compute(points, torus);
            expect(result.success).toBe(true);
            expect(Math.abs(dot(torus.N, N))).toBeCloseTo(1, 6);
            expect(torus.r0).toBeCloseTo(r0, 6);
            expect(torus.r1).toBeCloseTo(r1, 6);
            for (let k = 0; k < 3; ++k) {
                expect(torus.C.values[k]).toBeCloseTo(C.values[k], 6);
            }
        }
    });

    it('tolerates noise in the samples', () => {
        const rand = makeRandom(1357);
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16,
            () => 0.005 * (rand() - 0.5));
        const torus = newParameters();
        const result = new ApprTorus3().compute(points, torus);

        expect(result.success).toBe(true);
        expect(torus.r0).toBeCloseTo(truthR0, 2);
        expect(torus.r1).toBeCloseTo(truthR1, 2);
        expect(Math.abs(dot(torus.N, truthN))).toBeGreaterThan(0.999);
    });
});

describe('ApprTorus3.computeGaussNewton', () => {
    it('refines a perturbed initial guess of an exact torus', () => {
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16);
        const torus: ApprTorus3Parameters = {
            C: v3(truthC.values[0] + 0.05, truthC.values[1] - 0.05,
                truthC.values[2] + 0.05),
            N: (() => {
                const n = v3(truthN.values[0] + 0.02, truthN.values[1],
                    truthN.values[2] - 0.02);
                normalize(n);
                return n;
            })(),
            r0: truthR0 + 0.1,
            r1: truthR1 - 0.05
        };
        const fitter = new ApprTorus3();
        const result = fitter.computeGaussNewton(points, 16, 1e-12, 1e-12,
            true, torus);

        expect(result.minError).toBeLessThan(1e-8);
        expect(torus.r0).toBeCloseTo(truthR0, 4);
        expect(torus.r1).toBeCloseTo(truthR1, 4);
        expect(Math.abs(dot(torus.N, truthN))).toBeCloseTo(1, 5);
        for (let k = 0; k < 3; ++k) {
            expect(torus.C.values[k]).toBeCloseTo(truthC.values[k], 4);
        }
    });

    it('computes its own initial guess when asked', () => {
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16);
        const torus = newParameters();
        const fitter = new ApprTorus3();
        const result = fitter.computeGaussNewton(points, 8, 1e-12, 1e-12,
            false, torus);

        expect(result.minError).toBeLessThan(1e-8);
        expect(torus.r0).toBeCloseTo(truthR0, 5);
        expect(torus.r1).toBeCloseTo(truthR1, 5);
        expect(Math.abs(dot(torus.N, truthN))).toBeCloseTo(1, 6);
    });

    it('returns the initial guess when maxIterations is zero', () => {
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 20, 12);
        const torus = newParameters();
        const fitter = new ApprTorus3();
        const result = fitter.computeGaussNewton(points, 0, 1e-12, 1e-12,
            false, torus);

        expect(result.numIterations).toBe(1);
        expect(result.converged).toBe(false);
        expect(torus.r0).toBeCloseTo(truthR0, 6);
        expect(torus.r1).toBeCloseTo(truthR1, 6);
    });
});

describe('ApprTorus3.computeLevenbergMarquardt', () => {
    it('refines a perturbed initial guess of an exact torus', () => {
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16);
        const torus: ApprTorus3Parameters = {
            C: v3(truthC.values[0] + 0.05, truthC.values[1] - 0.05,
                truthC.values[2] + 0.05),
            N: truthN.clone(),
            r0: truthR0 + 0.1,
            r1: truthR1 - 0.05
        };
        const fitter = new ApprTorus3();
        const result = fitter.computeLevenbergMarquardt(points, 16, 1e-12,
            1e-12, 0.001, 10, 4, true, torus);

        expect(result.minError).toBeLessThan(1e-8);
        expect(torus.r0).toBeCloseTo(truthR0, 4);
        expect(torus.r1).toBeCloseTo(truthR1, 4);
        expect(Math.abs(dot(torus.N, truthN))).toBeCloseTo(1, 5);
    });

    it('fits noisy samples', () => {
        const rand = makeRandom(90210);
        const points = torusPoints(truthC, truthN, truthR0, truthR1, 24, 16,
            () => 0.004 * (rand() - 0.5));
        const torus = newParameters();
        const fitter = new ApprTorus3();
        fitter.computeLevenbergMarquardt(points, 16, 1e-12, 1e-12,
            0.001, 10, 4, false, torus);

        expect(torus.r0).toBeCloseTo(truthR0, 2);
        expect(torus.r1).toBeCloseTo(truthR1, 2);
        expect(Math.abs(dot(torus.N, truthN))).toBeGreaterThan(0.999);
        for (let k = 0; k < 3; ++k) {
            expect(torus.C.values[k]).toBeCloseTo(truthC.values[k], 2);
        }
    });
});

// ---------------------------------------------------------------------------
// Property-based verification (VERIFYING.md).

describe('ApprTorus3 verification', () => {
    // Tori with r1 well below r0 (a ring torus) sampled over the whole
    // surface, which is the regime the header documents for compute().
    const config = fc.tuple(vector(3, -4, 4), unitVector(3), finite(2, 5),
        finite(0.4, 1.2), fc.integer({ min: 12, max: 20 }),
        fc.integer({ min: 8, max: 14 }));

    function build(t: [Vector, Vector, number, number, number, number],
        noise?: () => number): Vector[] {
        const [C, N, r0, r1, numTheta, numPhi] = t;
        return torusPoints(C, N, r0, r1, numTheta, numPhi, noise);
    }

    it('recovers the torus from exact samples', () => {
        check(config, (t) => {
            const [C, N, r0, r1] = t;
            const torus = newParameters();
            const result = new ApprTorus3().compute(build(t), torus);
            expect(result.success).toBe(true);
            for (let k = 0; k < 3; ++k) {
                expectClose(torus.C.values[k], C.values[k], 1e-8, 1e-8);
            }
            expect(Math.abs(dot(torus.N, N))).toBeGreaterThan(1 - 1e-10);
            expectClose(torus.r0, r0, 1e-7, 1e-8);
            expectClose(torus.r1, r1, 1e-7, 1e-8);
            // Exact data means H(u,v) is zero at the recovered pair.
            expect(result.error).toBeLessThan(1e-8 * Math.pow(r0, 8) + 1e-8);
        }, 60);
    });

    it('reports H(u,v) for the (u,v) pair it selects', () => {
        // The returned error is the value of
        //   H(u,v) = sum_i ((v + L[i])^2 - S[i]*u)^2
        // at the chosen root. Recomputing it from the returned radii is an
        // independent check of the a/b/c/d/e/f coefficient chain that reduces
        // the two normal equations to a cubic in v.
        const rnd = seededRandom(161803);
        check(config, (t) => {
            const points = build(t, () => 0.02 * (2 * rnd() - 1));
            const torus = newParameters();
            const result = new ApprTorus3().compute(points, torus);
            expect(result.success).toBe(true);
            let h = 0;
            for (const p of points) {
                h += torusResidual(torus, p) ** 2;
            }
            // The residual is a difference of degree-4 quantities, so the
            // comparison is relative to the accumulated magnitude.
            expectClose(result.error, h, 1e-9, 1e-6);
        }, 40);
    });

    it('is equivariant under rigid motions of the samples', () => {
        // The plane fit, the L/S moments and the cubic solve all commute with
        // a rigid motion, so the recovered radii are invariant and the center
        // moves with the samples.
        check(fc.tuple(config, orthonormalFrame(3), vector(3, -6, 6)),
            ([t, R, shift]) => {
                const points = build(t);
                const moved = points.map(p => v3(
                    dot(R[0], p) + shift.values[0],
                    dot(R[1], p) + shift.values[1],
                    dot(R[2], p) + shift.values[2]));

                const a = newParameters();
                const ra = new ApprTorus3().compute(points, a);
                const b = newParameters();
                const rb = new ApprTorus3().compute(moved, b);
                expect(ra.success && rb.success).toBe(true);

                expectClose(a.r0, b.r0, 1e-7, 1e-7);
                expectClose(a.r1, b.r1, 1e-7, 1e-7);
                const movedCenter = v3(
                    dot(R[0], a.C) + shift.values[0],
                    dot(R[1], a.C) + shift.values[1],
                    dot(R[2], a.C) + shift.values[2]);
                expectVectorClose(movedCenter, b.C, 1e-6, 1e-7);
            }, 40);
    });

    it('keeps an exact torus stationary under the minimizers', () => {
        // The residuals F[i] vanish at the true parameters, so the
        // Gauss-Newton and Levenberg-Marquardt steps solve a homogeneous
        // system and cannot move away from the solution.
        check(config, (t) => {
            const [C, N, r0, r1] = t;
            for (const useLM of [false, true]) {
                const torus = newParameters();
                torus.C = C.clone();
                torus.N = N.clone();
                torus.r0 = r0;
                torus.r1 = r1;
                const points = build(t);
                const fitter = new ApprTorus3();
                if (useLM) {
                    fitter.computeLevenbergMarquardt(points, 4, 1e-12, 1e-12,
                        2, 0.5, 4, true, torus);
                } else {
                    fitter.computeGaussNewton(points, 4, 1e-12, 1e-12, true,
                        torus);
                }
                expectClose(torus.r0, r0, 1e-6, 1e-7);
                expectClose(torus.r1, r1, 1e-6, 1e-7);
                expect(Math.abs(dot(torus.N, N))).toBeGreaterThan(1 - 1e-8);
                expectVectorClose(torus.C, C, 1e-5, 1e-6);
            }
        }, 30);
    });

    it('reports failure on the documented degenerate inputs', () => {
        // Coincident samples: ApprOrthogonalPlane3 cannot fit a plane.
        const repeated = new Array<Vector>(12).fill(v3(1, 2, 3));
        const a = newParameters();
        const ra = new ApprTorus3().compute(repeated, a);
        expect(ra.success).toBe(false);
        expect(ra.error).toBe(Number.MAX_VALUE);

        // Coplanar samples on a circle: the plane fit succeeds but no cubic
        // root yields a valid u > v > 0 pair, so 'success' is false rather
        // than a garbage torus (upstream's unguarded 1/b0 and SolveCubic
        // leading coefficient produce NaNs that fail the same validity
        // tests).
        const ring: Vector[] = [];
        for (let j = 0; j < 24; ++j) {
            const t = 2 * Math.PI * j / 24;
            ring.push(v3(2 * Math.cos(t), 2 * Math.sin(t), 0));
        }
        const b = newParameters();
        const rb = new ApprTorus3().compute(ring, b);
        expect(rb.success).toBe(false);
        expect(rb.error).toBe(Number.MAX_VALUE);
    });
});
