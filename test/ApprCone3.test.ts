import { describe, it, expect } from 'vitest';
import { ApprCone3 } from '../src/ApprCone3.js';
import type { ApprCone3Parameters } from '../src/ApprCone3.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import { computeOrthogonalComplement3 } from '../src/Vector3.js';
import {
    check, expectClose, fc, finite, unitVector, vector
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

// Points on the cone with vertex V, unit-length axis U and angle A. The
// point at height h and polar angle t is
//   X = V + h * U + h * tan(A) * (cos(t) * E0 + sin(t) * E1)
function conePoints(V: Vector, U: Vector, angle: number, hMin: number,
    hMax: number, numHeights: number, numAngles: number,
    noise?: () => number): Vector[] {
    const basis = [U.clone(), new Vector(3), new Vector(3)];
    computeOrthogonalComplement3(1, basis);
    const E0 = basis[1], E1 = basis[2];
    const tanA = Math.tan(angle);
    const points: Vector[] = [];
    for (let i = 0; i < numHeights; ++i) {
        const h = hMin + (hMax - hMin) * i / (numHeights - 1);
        const r = h * tanA;
        for (let j = 0; j < numAngles; ++j) {
            const t = 2 * Math.PI * j / numAngles;
            const cs = Math.cos(t), sn = Math.sin(t);
            const p = new Vector(3);
            for (let k = 0; k < 3; ++k) {
                p.values[k] = V.values[k] + h * U.values[k]
                    + r * (cs * E0.values[k] + sn * E1.values[k])
                    + (noise ? noise() : 0);
            }
            points.push(p);
        }
    }
    return points;
}

// The implicit cone function F(X) = |D|^2 - Dot(D,W)^2 with D = V - X and
// W = U / cos(A). It vanishes for points on the (double) cone.
function coneResidual(cone: ApprCone3Parameters, X: Vector): number {
    const D = sub(cone.vertex, X);
    const W = cone.axis.clone();
    const invCos = 1 / Math.cos(cone.angle);
    for (let k = 0; k < 3; ++k) {
        W.values[k] *= invCos;
    }
    const dw = dot(D, W);
    return dot(D, D) - dw * dw;
}

const truthVertex = v3(1, 2, 3);
const truthAxis = (() => {
    const u = v3(1, 2, 2);
    normalize(u);
    return u;
})();
const truthAngle = 0.4;

describe('ApprCone3 initial estimate', () => {
    it('estimates the vertex, axis and angle before any iteration', () => {
        // maxIterations = 0 returns the initial guess, so this exercises the
        // ComputeInitialCone code path in isolation.
        const points = conePoints(truthVertex, truthAxis, truthAngle,
            1, 3, 9, 16);
        const cone: ApprCone3Parameters = {
            vertex: new Vector(3), axis: new Vector(3), angle: 0
        };
        const fitter = new ApprCone3();
        fitter.computeGaussNewton(points, 0, 1e-8, 1e-8, false, cone);

        expect(dot(cone.axis, truthAxis)).toBeGreaterThan(0.99);
        expect(cone.angle).toBeCloseTo(truthAngle, 2);
        for (let k = 0; k < 3; ++k) {
            expect(cone.vertex.values[k]).toBeCloseTo(truthVertex.values[k], 2);
        }
    });

    it('flips the axis estimate when the radius decreases with height', () => {
        // Sampling the same cone with the axis reversed (heights measured
        // along -U) must recover the same cone.
        const reversed = v3(-truthAxis.values[0], -truthAxis.values[1],
            -truthAxis.values[2]);
        const points = conePoints(truthVertex, reversed, truthAngle,
            1, 3, 9, 16);
        const cone: ApprCone3Parameters = {
            vertex: new Vector(3), axis: new Vector(3), angle: 0
        };
        const fitter = new ApprCone3();
        fitter.computeGaussNewton(points, 0, 1e-8, 1e-8, false, cone);

        expect(dot(cone.axis, reversed)).toBeGreaterThan(0.99);
        expect(cone.angle).toBeCloseTo(truthAngle, 2);
    });
});

describe('ApprCone3.computeGaussNewton', () => {
    it('recovers a cone from points that lie exactly on it', () => {
        const points = conePoints(truthVertex, truthAxis, truthAngle,
            1, 3, 9, 16);
        const cone: ApprCone3Parameters = {
            vertex: new Vector(3), axis: new Vector(3), angle: 0
        };
        const fitter = new ApprCone3();
        const result = fitter.computeGaussNewton(points, 32, 1e-12, 1e-12,
            false, cone);

        expect(result.minError).toBeLessThan(1e-12);
        expect(dot(cone.axis, truthAxis)).toBeCloseTo(1, 6);
        expect(cone.angle).toBeCloseTo(truthAngle, 6);
        for (let k = 0; k < 3; ++k) {
            expect(cone.vertex.values[k]).toBeCloseTo(truthVertex.values[k], 5);
        }
        for (const p of points) {
            expect(Math.abs(coneResidual(cone, p))).toBeLessThan(1e-6);
        }
    });

    it('accepts a caller-supplied initial guess and normalizes the axis', () => {
        const points = conePoints(truthVertex, truthAxis, truthAngle,
            1, 3, 9, 16);
        const cone: ApprCone3Parameters = {
            // A deliberately non-unit-length axis near the truth.
            vertex: v3(1.2, 1.8, 3.3),
            axis: v3(5 * truthAxis.values[0] + 0.1,
                5 * truthAxis.values[1], 5 * truthAxis.values[2]),
            angle: 0.45
        };
        const fitter = new ApprCone3();
        const result = fitter.computeGaussNewton(points, 32, 1e-12, 1e-12,
            true, cone);

        expect(result.minError).toBeLessThan(1e-10);
        expect(dot(cone.axis, truthAxis)).toBeCloseTo(1, 5);
        expect(cone.angle).toBeCloseTo(truthAngle, 5);
    });

    it('recovers random cones from exact samples', () => {
        const rand = makeRandom(1234567);
        for (let trial = 0; trial < 5; ++trial) {
            const V = v3(-2 + 4 * rand(), -2 + 4 * rand(), -2 + 4 * rand());
            const U = v3(-1 + 2 * rand(), -1 + 2 * rand(), -1 + 2 * rand());
            if (normalize(U) === 0) {
                continue;
            }
            const angle = 0.25 + 0.5 * rand();
            const points = conePoints(V, U, angle, 1, 3, 7, 12);

            const cone: ApprCone3Parameters = {
                vertex: new Vector(3), axis: new Vector(3), angle: 0
            };
            const fitter = new ApprCone3();
            fitter.computeGaussNewton(points, 32, 1e-12, 1e-12, false, cone);

            expect(dot(cone.axis, U)).toBeCloseTo(1, 4);
            expect(cone.angle).toBeCloseTo(angle, 4);
            for (let k = 0; k < 3; ++k) {
                expect(cone.vertex.values[k]).toBeCloseTo(V.values[k], 3);
            }
        }
    });
});

describe('ApprCone3.computeLevenbergMarquardt', () => {
    it('recovers a cone from points that lie exactly on it', () => {
        const points = conePoints(truthVertex, truthAxis, truthAngle,
            1, 3, 9, 16);
        const cone: ApprCone3Parameters = {
            vertex: new Vector(3), axis: new Vector(3), angle: 0
        };
        const fitter = new ApprCone3();
        const result = fitter.computeLevenbergMarquardt(points, 32, 1e-12,
            1e-12, 0.001, 10, 8, false, cone);

        expect(result.minError).toBeLessThan(1e-10);
        expect(dot(cone.axis, truthAxis)).toBeCloseTo(1, 5);
        expect(cone.angle).toBeCloseTo(truthAngle, 5);
        for (let k = 0; k < 3; ++k) {
            expect(cone.vertex.values[k]).toBeCloseTo(truthVertex.values[k], 4);
        }
    });

    it('recovers a cone from noisy samples', () => {
        const rand = makeRandom(24680);
        const points = conePoints(truthVertex, truthAxis, truthAngle,
            1, 3, 9, 16, () => 0.002 * (rand() - 0.5));
        const cone: ApprCone3Parameters = {
            vertex: new Vector(3), axis: new Vector(3), angle: 0
        };
        const fitter = new ApprCone3();
        fitter.computeLevenbergMarquardt(points, 32, 1e-12, 1e-12,
            0.001, 10, 8, false, cone);

        expect(dot(cone.axis, truthAxis)).toBeGreaterThan(0.999);
        expect(cone.angle).toBeCloseTo(truthAngle, 2);
        for (let k = 0; k < 3; ++k) {
            expect(cone.vertex.values[k]).toBeCloseTo(truthVertex.values[k], 1);
        }
    });
});

// ---------------------------------------------------------------------------
// Property-based verification (VERIFYING.md).

describe('ApprCone3 verification', () => {
    function newCone(): ApprCone3Parameters {
        return { vertex: new Vector(3), axis: new Vector(3), angle: 0 };
    }

    // Cones whose half-angle stays well inside (0, pi/2) and whose sampled
    // height band is bounded away from the vertex, so the fit is conditioned.
    const config = fc.tuple(vector(3, -4, 4), unitVector(3),
        finite(0.2, 1.1), finite(0.8, 2), finite(2.5, 5));

    function build(t: [Vector, Vector, number, number, number],
        noise?: () => number): Vector[] {
        const [V, U, angle, hMin, hMax] = t;
        return conePoints(V, U, angle, hMin, hMax, 7, 12, noise);
    }

    it('keeps an exact cone stationary under Gauss-Newton', () => {
        // The residuals F[i] vanish at the true (V,W), so every Gauss-Newton
        // step solves J^T J * delta = -J^T F = 0. The minimizer must stay at
        // the true cone and report a zero error.
        check(config, (t) => {
            const [V, U, angle] = t;
            const points = build(t);
            const cone = newCone();
            cone.vertex = V.clone();
            cone.axis = U.clone();
            cone.angle = angle;
            const result = new ApprCone3().computeGaussNewton(points, 8,
                1e-12, 1e-12, true, cone);
            expect(result.minError).toBeLessThan(1e-12);
            expect(dot(cone.axis, U)).toBeGreaterThan(1 - 1e-9);
            expectClose(cone.angle, angle, 1e-6, 1e-6);
            for (let k = 0; k < 3; ++k) {
                expectClose(cone.vertex.values[k], V.values[k], 1e-5, 1e-6);
            }
        }, 60);
    });

    it('gives the same fixed point for Levenberg-Marquardt', () => {
        check(config, (t) => {
            const [V, U, angle] = t;
            const points = build(t);
            const cone = newCone();
            cone.vertex = V.clone();
            cone.axis = U.clone();
            cone.angle = angle;
            const result = new ApprCone3().computeLevenbergMarquardt(points,
                8, 1e-12, 1e-12, 2, 0.5, 4, true, cone);
            expect(result.minError).toBeLessThan(1e-12);
            expect(dot(cone.axis, U)).toBeGreaterThan(1 - 1e-9);
            expectClose(cone.angle, angle, 1e-6, 1e-6);
        }, 40);
    });

    it('produces a residual-free cone for exact samples', () => {
        // Independent check of the fitted cone: the implicit cone function
        // F(X) = |D|^2 - Dot(D, U/cos A)^2 must vanish at every sample.
        check(config, (t) => {
            const points = build(t);
            const cone = newCone();
            new ApprCone3().computeGaussNewton(points, 32, 1e-14, 1e-14,
                false, cone);
            let scale = 0;
            for (const p of points) {
                const d = sub(cone.vertex, p);
                scale = Math.max(scale, dot(d, d));
            }
            for (const p of points) {
                // Relative to |D|^2: the residual is a difference of two
                // quantities of that size, so round-off scales with it.
                expect(Math.abs(coneResidual(cone, p)))
                    .toBeLessThan(1e-6 * scale + 1e-9);
            }
        }, 40);
    });

    it('is equivariant under translation of the samples', () => {
        check(fc.tuple(config, vector(3, -6, 6)), ([t, shift]) => {
            const points = build(t);
            const moved = points.map(p => v3(p.values[0] + shift.values[0],
                p.values[1] + shift.values[1], p.values[2] + shift.values[2]));

            const a = newCone();
            new ApprCone3().computeGaussNewton(points, 0, 1e-12, 1e-12,
                false, a);
            const b = newCone();
            new ApprCone3().computeGaussNewton(moved, 0, 1e-12, 1e-12,
                false, b);

            // maxIterations = 0 isolates ComputeInitialCone, whose steps -
            // the sample average, the ZZTZ axis estimate, the (h,r) height
            // line and the vertex offset - are all translation equivariant.
            expectClose(a.angle, b.angle, 1e-8, 1e-8);
            for (let k = 0; k < 3; ++k) {
                expectClose(a.axis.values[k], b.axis.values[k], 1e-8, 1e-8);
                expectClose(a.vertex.values[k] + shift.values[k],
                    b.vertex.values[k], 1e-6, 1e-7);
            }
        }, 60);
    });

    it('estimates the same cone when the sampling axis is reversed', () => {
        // ComputeInitialCone flips the ZZTZ axis estimate when the fitted
        // (h,r) slope is negative, mirroring hMin and hMax. Sampling the same
        // cone along -U must therefore recover the same cone.
        check(config, (t) => {
            const [V, U, angle, hMin, hMax] = t;
            const forward = conePoints(V, U, angle, hMin, hMax, 7, 12);
            const reversed = conePoints(V,
                v3(-U.values[0], -U.values[1], -U.values[2]), angle,
                hMin, hMax, 7, 12);

            const a = newCone();
            new ApprCone3().computeGaussNewton(forward, 0, 1e-12, 1e-12,
                false, a);
            const b = newCone();
            new ApprCone3().computeGaussNewton(reversed, 0, 1e-12, 1e-12,
                false, b);

            expect(dot(a.axis, U)).toBeGreaterThan(0.99);
            expect(dot(b.axis, U)).toBeLessThan(-0.99);
            expectClose(a.angle, b.angle, 1e-6, 1e-6);
            for (let k = 0; k < 3; ++k) {
                expectClose(a.vertex.values[k], b.vertex.values[k],
                    1e-5, 1e-6);
            }
        }, 60);
    });

    it('preserves the unguarded degenerate paths of ComputeInitialCone', () => {
        // Upstream bug suspect (preserved): ComputeInitialCone divides by
        // hRange and by tanAngle without guards. Coincident samples make both
        // ranges zero, so tanAngle is 0/0 and the initial cone is NaN. The
        // port reproduces that rather than adding a guard.
        const repeated = new Array<Vector>(12).fill(v3(1, 2, 3));
        const cone = newCone();
        new ApprCone3().computeGaussNewton(repeated, 0, 1e-8, 1e-8, false,
            cone);
        expect(Number.isNaN(cone.vertex.values[0])).toBe(true);
        // The ZZTZ axis estimate is the zero vector, so the reported angle is
        // acos(min(1/0, 1)) = 0 rather than NaN; only the vertex carries the
        // 0/0 from tanAngle.
        expect(cone.angle).toBe(0);
        expect(length(cone.axis)).toBe(0);

        // Note that rRange = hrSlope * hRange, so a zero hRange always makes
        // tanAngle a 0/0 NaN; the "all points at one height" case in the
        // upstream report can never merely produce a large quotient.
    });
});
