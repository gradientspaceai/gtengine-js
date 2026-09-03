import { describe, expect, it } from 'vitest';
import { MinimumAreaCircle2 } from '../src/MinimumAreaCircle2.js';
import { Vector, sub, dot } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function distance(p: Vector, q: Vector): number {
    const d = sub(p, q);
    return Math.sqrt(dot(d, d));
}

// An independent brute-force minimum-area circle: enumerate the circles
// defined by every pair (diametral circle) and every triple (circumcircle,
// computed with the determinant formula, which is a different derivation
// from the barycentric solve used by the port) and keep the smallest one
// that contains every point.
function bruteForceMinimumCircle(points: readonly Vector[]):
    { center: Vector, radius: number } {
    const n = points.length;
    if (n === 1) {
        return { center: points[0].clone(), radius: 0 };
    }

    let best: { center: Vector, radius: number } | null = null;
    const consider = (center: Vector, radius: number): void => {
        // Allow a small tolerance so that the defining points, which sit
        // exactly on the boundary in exact arithmetic, are not rejected.
        const tol = 1e-9 * Math.max(1, radius);
        for (const p of points) {
            if (distance(p, center) > radius + tol) {
                return;
            }
        }
        if (best === null || radius < best.radius) {
            best = { center, radius };
        }
    };

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const center = v2(
                0.5 * (points[i].get(0) + points[j].get(0)),
                0.5 * (points[i].get(1) + points[j].get(1)));
            consider(center, 0.5 * distance(points[i], points[j]));
        }
    }

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            for (let k = j + 1; k < n; ++k) {
                const ax = points[i].get(0), ay = points[i].get(1);
                const bx = points[j].get(0), by = points[j].get(1);
                const cx = points[k].get(0), cy = points[k].get(1);
                const d = 2 * (ax * (by - cy) + bx * (cy - ay)
                    + cx * (ay - by));
                if (Math.abs(d) < 1e-14) {
                    continue;
                }
                const a2 = ax * ax + ay * ay;
                const b2 = bx * bx + by * by;
                const c2 = cx * cx + cy * cy;
                const ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by))
                    / d;
                const uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax))
                    / d;
                const center = v2(ux, uy);
                consider(center, distance(center, points[i]));
            }
        }
    }

    if (best === null) {
        throw new Error('brute force failed');
    }
    return best;
}

// A deterministic pseudo-random generator for the randomized tests.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (1664525 * state + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('MinimumAreaCircle2', () => {
    it('handles a single point', () => {
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute([v2(3, -4)]);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(minimal.center.get(0)).toBe(3);
        expect(minimal.center.get(1)).toBe(-4);
        expect(mac.numSupport).toBe(1);
        expect(mac.support[0]).toBe(0);
    });

    it('computes the diametral circle of two points', () => {
        const mac = new MinimumAreaCircle2();
        const points = [v2(-1, 0), v2(1, 0)];
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0, 12);
        expect(mac.numSupport).toBe(2);
        expect(new Set(mac.support.slice(0, 2))).toEqual(new Set([0, 1]));
    });

    it('computes the circumcircle of an acute triangle', () => {
        // Equilateral triangle with circumradius 1.
        const points = [
            v2(1, 0),
            v2(Math.cos(2 * Math.PI / 3), Math.sin(2 * Math.PI / 3)),
            v2(Math.cos(4 * Math.PI / 3), Math.sin(4 * Math.PI / 3))
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 10);
        expect(minimal.center.get(0)).toBeCloseTo(0, 10);
        expect(minimal.center.get(1)).toBeCloseTo(0, 10);
        expect(mac.numSupport).toBe(3);
    });

    it('uses the longest side for an obtuse triangle', () => {
        // The circumcircle of this triangle is larger than the circle whose
        // diameter is the longest side, so the minimum circle has 2 support
        // points.
        const points = [v2(-2, 0), v2(2, 0), v2(0, 0.25)];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0, 12);
        expect(mac.numSupport).toBe(2);
        expect(new Set(mac.support.slice(0, 2))).toEqual(new Set([0, 1]));
    });

    it('computes the circle of a unit square', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const { minimal, success } = new MinimumAreaCircle2().compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(Math.SQRT1_2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0.5, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0.5, 12);
    });

    it('ignores duplicate points', () => {
        const points = [
            v2(0, 0), v2(1, 0), v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1),
            v2(0, 1)
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(Math.SQRT1_2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0.5, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0.5, 12);
        // The support indices index into the input array.
        for (let i = 0; i < mac.numSupport; ++i) {
            const index = mac.support[i];
            expect(index).toBeGreaterThanOrEqual(0);
            expect(index).toBeLessThan(points.length);
            expect(distance(points[index], minimal.center))
                .toBeCloseTo(minimal.radius, 10);
        }
    });

    it('handles a set of identical points', () => {
        const points = [v2(5, 5), v2(5, 5), v2(5, 5)];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(minimal.center.get(0)).toBe(5);
        expect(mac.numSupport).toBe(1);
    });

    it('handles collinear points', () => {
        const points = [
            v2(0, 0), v2(1, 1), v2(2, 2), v2(3, 3), v2(-1, -1)
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        // The extremes are (-1,-1) and (3,3).
        expect(minimal.center.get(0)).toBeCloseTo(1, 10);
        expect(minimal.center.get(1)).toBeCloseTo(1, 10);
        expect(minimal.radius).toBeCloseTo(2 * Math.SQRT2, 10);
        expect(mac.numSupport).toBe(2);
    });

    it('handles cocircular points', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v2(2 * Math.cos(angle), 2 * Math.sin(angle)));
        }
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(2, 8);
        expect(minimal.center.get(0)).toBeCloseTo(0, 8);
        expect(minimal.center.get(1)).toBeCloseTo(0, 8);
    });

    it('throws for an empty point set', () => {
        expect(() => new MinimumAreaCircle2().compute([]))
            .toThrow('Input must contain points.');
    });

    it('throws for non-2D points', () => {
        expect(() => new MinimumAreaCircle2()
            .compute([Vector.fromArray([1, 2, 3])])).toThrow();
    });

    it('matches brute force on random point sets', () => {
        const random = makeRandom(20260902);
        for (let trial = 0; trial < 40; ++trial) {
            const numPoints = 3 + (trial % 10);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                points.push(v2(10 * random() - 5, 10 * random() - 5));
            }

            const mac = new MinimumAreaCircle2();
            const { minimal, success } = mac.compute(points);
            expect(success).toBe(true);

            const expected = bruteForceMinimumCircle(points);
            expect(minimal.radius).toBeCloseTo(expected.radius, 8);
            expect(minimal.center.get(0))
                .toBeCloseTo(expected.center.get(0), 7);
            expect(minimal.center.get(1))
                .toBeCloseTo(expected.center.get(1), 7);

            // Every point is contained.
            for (const p of points) {
                expect(distance(p, minimal.center))
                    .toBeLessThanOrEqual(minimal.radius + 1e-9);
            }

            // Every support point is on the boundary and the center is a
            // convex combination of the support points (the optimality
            // certificate for the minimum enclosing ball).
            expect(mac.numSupport).toBeGreaterThanOrEqual(2);
            for (let i = 0; i < mac.numSupport; ++i) {
                expect(distance(points[mac.support[i]], minimal.center))
                    .toBeCloseTo(minimal.radius, 7);
            }
            if (mac.numSupport === 2) {
                const p0 = points[mac.support[0]];
                const p1 = points[mac.support[1]];
                expect(minimal.center.get(0))
                    .toBeCloseTo(0.5 * (p0.get(0) + p1.get(0)), 7);
                expect(minimal.center.get(1))
                    .toBeCloseTo(0.5 * (p0.get(1) + p1.get(1)), 7);
            } else {
                // Solve center = w0*P0 + w1*P1 + (1-w0-w1)*P2 and require
                // nonnegative barycentric weights.
                const P0 = points[mac.support[0]];
                const P1 = points[mac.support[1]];
                const P2 = points[mac.support[2]];
                const a00 = P0.get(0) - P2.get(0);
                const a01 = P1.get(0) - P2.get(0);
                const a10 = P0.get(1) - P2.get(1);
                const a11 = P1.get(1) - P2.get(1);
                const b0 = minimal.center.get(0) - P2.get(0);
                const b1 = minimal.center.get(1) - P2.get(1);
                const det = a00 * a11 - a01 * a10;
                expect(Math.abs(det)).toBeGreaterThan(1e-12);
                const w0 = (b0 * a11 - a01 * b1) / det;
                const w1 = (a00 * b1 - b0 * a10) / det;
                const w2 = 1 - w0 - w1;
                expect(w0).toBeGreaterThan(-1e-8);
                expect(w1).toBeGreaterThan(-1e-8);
                expect(w2).toBeGreaterThan(-1e-8);
            }
        }
    });

    it('is minimal: shrinking the radius excludes a point', () => {
        const random = makeRandom(777);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v2(random(), random()));
            }
            const { minimal, success } =
                new MinimumAreaCircle2().compute(points);
            expect(success).toBe(true);

            // Perturbing the center in any of four directions while keeping
            // the radius must leave some point outside.
            const eps = 1e-4;
            const directions = [[1, 0], [-1, 0], [0, 1], [0, -1],
                [0.7, 0.7], [-0.7, 0.7], [0.7, -0.7], [-0.7, -0.7]];
            for (const d of directions) {
                const c = v2(minimal.center.get(0) + eps * d[0],
                    minimal.center.get(1) + eps * d[1]);
                let maxDist = 0;
                for (const p of points) {
                    maxDist = Math.max(maxDist, distance(p, c));
                }
                expect(maxDist).toBeGreaterThan(minimal.radius - 1e-12);
            }
        }
    });
});
