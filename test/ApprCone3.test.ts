import { describe, it, expect } from 'vitest';
import { ApprCone3 } from '../src/ApprCone3';
import type { ApprCone3Parameters } from '../src/ApprCone3';
import { Vector, dot, normalize, sub } from '../src/Vector';
import { computeOrthogonalComplement3 } from '../src/Vector3';

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
