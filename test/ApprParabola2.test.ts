import { describe, it, expect } from 'vitest';
import { ApprParabola2 } from '../src/ApprParabola2';
import { Vector } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Samples of y = u0*x^2 + u1*x + u2 at x = xMin + i * (xMax - xMin) / (n - 1).
function parabolaSamples(u: readonly number[], xMin: number, xMax: number,
    count: number, noise?: (i: number) => number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        const x = xMin + (xMax - xMin) * i / (count - 1);
        const y = u[0] * x * x + u[1] * x + u[2] + (noise ? noise(i) : 0);
        points.push(v2(x, y));
    }
    return points;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprParabola2.fit', () => {
    it('recovers the coefficients of an exact parabola', () => {
        const u = [2, -3, 1];
        const result = ApprParabola2.fit(parabolaSamples(u, -2, 3, 25));
        expect(result.success).toBe(true);
        expect(result.u[0]).toBeCloseTo(u[0], 10);
        expect(result.u[1]).toBeCloseTo(u[1], 10);
        expect(result.u[2]).toBeCloseTo(u[2], 10);
        expect(result.meanSquareError).toBeLessThan(1e-10);
    });

    it('fits a line exactly (the quadratic coefficient is zero)', () => {
        const result = ApprParabola2.fit(parabolaSamples([0, 4, -5], 0, 10, 11));
        expect(result.success).toBe(true);
        expect(result.u[0]).toBeCloseTo(0, 9);
        expect(result.u[1]).toBeCloseTo(4, 9);
        expect(result.u[2]).toBeCloseTo(-5, 9);
    });

    it('is the least-squares minimizer (normal equations hold)', () => {
        // The gradient of sum_i (u0*x^2 + u1*x + u2 - y)^2 with respect to
        // each of u0, u1 and u2 must vanish at the fit.
        const rand = makeRandom(12345);
        const points: Vector[] = [];
        for (let i = 0; i < 40; ++i) {
            const x = -3 + 6 * rand();
            points.push(v2(x, 1.5 * x * x - 0.5 * x + 2 + (rand() - 0.5)));
        }
        const { success, u } = ApprParabola2.fit(points);
        expect(success).toBe(true);

        let g0 = 0, g1 = 0, g2 = 0;
        for (const p of points) {
            const x = p.values[0], y = p.values[1];
            const e = u[0] * x * x + u[1] * x + u[2] - y;
            g0 += e * x * x;
            g1 += e * x;
            g2 += e;
        }
        expect(Math.abs(g0) / points.length).toBeLessThan(1e-12);
        expect(Math.abs(g1) / points.length).toBeLessThan(1e-12);
        expect(Math.abs(g2) / points.length).toBeLessThan(1e-12);
    });

    it('recovers coefficients from noisy samples', () => {
        const rand = makeRandom(999);
        const u = [-1.25, 0.75, 3];
        const points = parabolaSamples(u, -4, 4, 200,
            () => 0.02 * (rand() - 0.5));
        const result = ApprParabola2.fit(points);
        expect(result.success).toBe(true);
        expect(result.u[0]).toBeCloseTo(u[0], 3);
        expect(result.u[1]).toBeCloseTo(u[1], 3);
        expect(result.u[2]).toBeCloseTo(u[2], 3);
        expect(result.meanSquareError).toBeLessThan(0.01);
    });

    it('reports failure for a degenerate (single-abscissa) point set', () => {
        const points = [v2(2, 0), v2(2, 1), v2(2, -1), v2(2, 4)];
        const result = ApprParabola2.fit(points);
        expect(result.success).toBe(false);
        expect(result.u).toEqual([0, 0, 0]);
    });

    it('throws when there are too few points', () => {
        expect(() => ApprParabola2.fit([v2(0, 0), v2(1, 1)])).toThrow(
            /Insufficient points/);
    });

    it('recovers random parabolas from exact samples', () => {
        const rand = makeRandom(2024);
        for (let trial = 0; trial < 20; ++trial) {
            const u = [
                -2 + 4 * rand(),
                -2 + 4 * rand(),
                -2 + 4 * rand()
            ];
            if (Math.abs(u[0]) < 0.1) {
                u[0] = 0.5;
            }
            const points = parabolaSamples(u, -1.5, 2.5, 17);
            const result = ApprParabola2.fit(points);
            expect(result.success).toBe(true);
            for (let i = 0; i < 3; ++i) {
                expect(result.u[i]).toBeCloseTo(u[i], 8);
            }
        }
    });
});

describe('ApprParabola2.fitRobust', () => {
    it('recovers the shifted coefficients and the average', () => {
        const u = [2, -3, 1];
        const points = parabolaSamples(u, -2, 3, 25);
        const result = ApprParabola2.fitRobust(points);
        expect(result.success).toBe(true);

        let ax = 0, ay = 0;
        for (const p of points) {
            ax += p.values[0];
            ay += p.values[1];
        }
        ax /= points.length;
        ay /= points.length;
        expect(result.average.values[0]).toBeCloseTo(ax, 12);
        expect(result.average.values[1]).toBeCloseTo(ay, 12);

        // Convert the v-polynomial back to the u-polynomial:
        //   u0 = v0, u1 = v1 - 2*v0*a, u2 = v0*a^2 - v1*a + v2 + b.
        const a = result.average.values[0], b = result.average.values[1];
        const v = result.v;
        expect(v[0]).toBeCloseTo(u[0], 10);
        expect(v[1] - 2 * v[0] * a).toBeCloseTo(u[1], 9);
        expect(v[0] * a * a - v[1] * a + v[2] + b).toBeCloseTo(u[2], 9);
        expect(result.meanSquareError).toBeLessThan(1e-10);
    });

    it('agrees with fit on well-conditioned data', () => {
        const rand = makeRandom(4242);
        const points: Vector[] = [];
        for (let i = 0; i < 60; ++i) {
            const x = -2 + 4 * rand();
            points.push(v2(x, 0.3 * x * x + 1.1 * x - 0.4 + 0.1 * (rand() - 0.5)));
        }
        const direct = ApprParabola2.fit(points);
        const robust = ApprParabola2.fitRobust(points);
        expect(direct.success && robust.success).toBe(true);

        const a = robust.average.values[0], b = robust.average.values[1];
        const v = robust.v;
        expect(v[0]).toBeCloseTo(direct.u[0], 8);
        expect(v[1] - 2 * v[0] * a).toBeCloseTo(direct.u[1], 8);
        expect(v[0] * a * a - v[1] * a + v[2] + b).toBeCloseTo(direct.u[2], 8);
    });

    it('throws when there are too few points', () => {
        expect(() => ApprParabola2.fitRobust([v2(0, 0)])).toThrow(
            /Insufficient points/);
    });
});
