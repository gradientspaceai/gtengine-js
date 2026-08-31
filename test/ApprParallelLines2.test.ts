import { describe, it, expect } from 'vitest';
import { ApprParallelLines2 } from '../src/ApprParallelLines2';
import { Vector, dot } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Samples that lie exactly on the two lines C -+ r*U + s*V, where
// V = (cos(angle),sin(angle)) and U = (-sin(angle),cos(angle)).
function twoLineSamples(angle: number, cx: number, cy: number, r: number,
    parameters: readonly number[]): Vector[] {
    const V = v2(Math.cos(angle), Math.sin(angle));
    const U = v2(-Math.sin(angle), Math.cos(angle));
    const points: Vector[] = [];
    for (const sign of [-1, 1]) {
        for (const s of parameters) {
            points.push(v2(
                cx + sign * r * U.values[0] + s * V.values[0],
                cy + sign * r * U.values[1] + s * V.values[1]));
        }
    }
    return points;
}

// The quantity the algorithm minimizes: the average of (u^2 - r^2)^2, where
// u is the component of the sample along the normal U of the two lines. It
// is zero exactly when every sample lies on one of the two lines.
function fitError(points: readonly Vector[], center: Vector,
    direction: Vector, radius: number): number {
    const U = v2(-direction.values[1], direction.values[0]);
    let sum = 0;
    for (const p of points) {
        const diff = v2(p.values[0] - center.values[0],
            p.values[1] - center.values[1]);
        const u = dot(diff, U);
        const d = u * u - radius * radius;
        sum += d * d;
    }
    return sum / points.length;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprParallelLines2', () => {
    it('recovers two axis-aligned parallel lines', () => {
        // The lines y = 3 and y = -3, so V = (1,0), C = (0,0), radius 3.
        const points = twoLineSamples(0, 0, 0, 3, [-4, -2, -1, 0, 1, 2, 4]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 1024);

        expect(Math.abs(result.direction.values[0])).toBeCloseTo(1, 8);
        expect(result.direction.values[1]).toBeCloseTo(0, 8);
        expect(result.radius).toBeCloseTo(3, 8);
        expect(result.center.values[1]).toBeCloseTo(0, 8);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 10);
    });

    it('recovers two rotated parallel lines', () => {
        const angle = 0.7;
        const points = twoLineSamples(angle, 2, -1, 1.5,
            [-3, -2, -1, 0, 1, 2, 3]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 1024);

        const V = v2(Math.cos(angle), Math.sin(angle));
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 6);
        expect(result.radius).toBeCloseTo(1.5, 6);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 8);
    });

    it('returns a unit-length direction with the center on the normal axis', () => {
        const angle = -0.35;
        const points = twoLineSamples(angle, -4, 5, 0.75, [-2, -1, 0, 1, 2, 5]);
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 1024);

        // V is unit length and C has no V-component (only a U-component).
        expect(dot(result.direction, result.direction)).toBeCloseTo(1, 10);
        expect(dot(result.center, result.direction)).toBeCloseTo(0, 8);
        expect(result.radius).toBeGreaterThan(0);
    });

    it('is invariant under translation of the samples along the lines', () => {
        const angle = 0.9;
        const base = twoLineSamples(angle, 0, 0, 2, [-2, -1, 0, 1, 2, 3]);
        const V = v2(Math.cos(angle), Math.sin(angle));
        const shifted = base.map(p => v2(p.values[0] + 7 * V.values[0],
            p.values[1] + 7 * V.values[1]));

        const fitter = new ApprParallelLines2();
        const a = fitter.fit(base, 1024);
        const b = fitter.fit(shifted, 1024);
        expect(b.radius).toBeCloseTo(a.radius, 6);
        expect(Math.abs(dot(a.direction, b.direction))).toBeCloseTo(1, 6);
    });

    it('recovers the separation of two clusters of noisy samples', () => {
        const random = makeRandom(20260830);
        const angle = 0.25;
        const V = v2(Math.cos(angle), Math.sin(angle));
        const U = v2(-Math.sin(angle), Math.cos(angle));
        const points: Vector[] = [];
        for (let i = 0; i < 40; ++i) {
            const sign = (i % 2 === 0) ? 1 : -1;
            const s = 10 * random() - 5;
            const noise = 0.01 * (2 * random() - 1);
            const u = sign * 2 + noise;
            points.push(v2(3 + u * U.values[0] + s * V.values[0],
                1 + u * U.values[1] + s * V.values[1]));
        }
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 2048);
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 2);
        expect(result.radius).toBeCloseTo(2, 2);
    });

    it('does not do worse than the trivial (0,1) starting estimate', () => {
        // The algorithm seeds the minimization with sigma = 0, gamma = 1 and
        // only replaces it when a candidate has a smaller error, so the
        // reported fit is at least as good as the seed for every input.
        const random = makeRandom(13579);
        for (let trial = 0; trial < 10; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 24; ++i) {
                points.push(v2(10 * random() - 5, 10 * random() - 5));
            }
            const fitter = new ApprParallelLines2();
            const result = fitter.fit(points, 1024);

            // The seed fit: direction (1,0) with the center and radius
            // computed from the moments of the mean-centered samples.
            let ax = 0, ay = 0;
            for (const p of points) {
                ax += p.values[0];
                ay += p.values[1];
            }
            ax /= points.length;
            ay /= points.length;
            let Z02 = 0, Z03 = 0;
            for (const p of points) {
                const y = p.values[1] - ay;
                Z02 += y * y;
                Z03 += y * y * y;
            }
            Z02 /= points.length;
            Z03 /= points.length;
            const k = Z03 / (2 * Z02);
            const seedCenter = v2(0, ay + k);
            const seedRadius = Math.sqrt(k * k + Z02);
            const seedError = fitError(points, seedCenter, v2(1, 0), seedRadius);

            const error = fitError(points, result.center, result.direction,
                result.radius);
            expect(error).toBeLessThanOrEqual(seedError + 1e-8);
            expect(dot(result.direction, result.direction)).toBeCloseTo(1, 8);
        }
    });

    it('recovers unevenly sampled lines with nonzero odd moments', () => {
        // Five samples on one line and three on the other, at an angle: the
        // mean-centered data is not centrally symmetric, so the odd moments
        // Z12, Z30, Z21 and Z03 that drive the polynomial F are all nonzero.
        const angle = 0.45;
        const V = v2(Math.cos(angle), Math.sin(angle));
        const U = v2(-Math.sin(angle), Math.cos(angle));
        const points: Vector[] = [];
        const push = (sign: number, s: number) => {
            points.push(v2(
                -2 + sign * 1.25 * U.values[0] + s * V.values[0],
                3 + sign * 1.25 * U.values[1] + s * V.values[1]));
        };
        for (const s of [-2, -1, 0, 1.5, 4]) {
            push(1, s);
        }
        for (const s of [-3, 0.5, 2]) {
            push(-1, s);
        }

        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 2048);
        expect(Math.abs(dot(result.direction, V))).toBeCloseTo(1, 6);
        expect(result.radius).toBeCloseTo(1.25, 6);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 6);
    });

    it('handles isotropic samples, for which no direction is preferred', () => {
        // Samples spread uniformly on a circle have Z11 = 0, Z20 = Z02,
        // vanishing odd moments and Z40 = Z04 = 3*Z22, so f0 and f1 are the
        // zero polynomial in exact arithmetic and the error is the same for
        // every direction. Whichever direction is reported, the geometry
        // must still be consistent: r^2 = average(u^2) = 1/2.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v2(Math.cos(angle), Math.sin(angle)));
        }
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 1024);
        expect(dot(result.direction, result.direction)).toBeCloseTo(1, 8);
        expect(result.radius).toBeCloseTo(Math.SQRT1_2, 8);
        expect(Number.isFinite(result.center.values[0])).toBe(true);
        expect(Number.isFinite(result.center.values[1])).toBe(true);
    });

    it('handles four points at the corners of a rectangle', () => {
        // The two long sides of the rectangle are the parallel lines.
        const points = [v2(-3, 1), v2(3, 1), v2(-3, -1), v2(3, -1)];
        const fitter = new ApprParallelLines2();
        const result = fitter.fit(points, 1024);
        expect(result.radius).toBeCloseTo(1, 8);
        expect(Math.abs(result.direction.values[0])).toBeCloseTo(1, 8);
        expect(fitError(points, result.center, result.direction,
            result.radius)).toBeCloseTo(0, 10);
    });
});
