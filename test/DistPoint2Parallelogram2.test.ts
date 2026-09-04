import { describe, expect, it } from 'vitest';
import { DistPoint2Parallelogram2 } from '../src/DistPoint2Parallelogram2.js';
import { Parallelogram2 } from '../src/Parallelogram2.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { dotPerp } from '../src/Vector2.js';
import { Matrix, multiplyATB, mulMatrix } from '../src/Matrix.js';
import { inverse2x2 } from '../src/Matrix2x2.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function pgm(center: number[], a0: number[], a1: number[]): Parallelogram2 {
    return Parallelogram2.fromCenterAxis(v(...center),
        [v(...a0), v(...a1)]);
}

// The exact minimum distance from a point to the solid parallelogram,
// obtained by a dense sampling of the (s0,s1) domain [-1,1]^2 followed by a
// local refinement.
function bruteForce(p: Vector, g: Parallelogram2): number {
    const at = (s0: number, s1: number): number => length(sub(p,
        add(g.center, add(mul(s0, g.axis[0]), mul(s1, g.axis[1])))));
    let best = Number.MAX_VALUE;
    let bs0 = 0;
    let bs1 = 0;
    const n = 200;
    for (let i = 0; i <= n; ++i) {
        const s0 = -1 + 2 * i / n;
        for (let j = 0; j <= n; ++j) {
            const s1 = -1 + 2 * j / n;
            const d = at(s0, s1);
            if (d < best) {
                best = d;
                bs0 = s0;
                bs1 = s1;
            }
        }
    }
    // Refine around the best grid point.
    let h = 2 / n;
    for (let pass = 0; pass < 60; ++pass) {
        for (const [d0, d1] of [[1, 0], [-1, 0], [0, 1], [0, -1],
            [1, 1], [1, -1], [-1, 1], [-1, -1]]) {
            const s0 = Math.max(-1, Math.min(1, bs0 + d0 * h));
            const s1 = Math.max(-1, Math.min(1, bs1 + d1 * h));
            const d = at(s0, s1);
            if (d < best) {
                best = d;
                bs0 = s0;
                bs1 = s1;
            }
        }
        h *= 0.6;
    }
    return best;
}

describe('DistPoint2Parallelogram2', () => {
    const query = new DistPoint2Parallelogram2();
    // The unit square [-1,1]^2.
    const square = pgm([0, 0], [1, 0], [0, 1]);

    it('returns zero distance for a point inside', () => {
        const result = query.compute(v(0.25, -0.5), square);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.sqrDistance).toBeCloseTo(0, 12);
        expect(result.closest[0].values).toEqual([0.25, -0.5]);
        expect(result.closest[1].values[0]).toBeCloseTo(0.25, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(-0.5, 12);
    });

    it('finds an edge point for a point beyond one edge', () => {
        const result = query.compute(v(3, 0.5), square);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 12);
    });

    it('finds a corner for a point beyond two edges', () => {
        const result = query.compute(v(4, 5), square);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
    });

    it('is symmetric under reflection of the query point', () => {
        const r0 = query.compute(v(2.5, 1.75), square);
        const r1 = query.compute(v(-2.5, -1.75), square);
        expect(r1.distance).toBeCloseTo(r0.distance, 12);
        expect(r1.closest[1].values[0]).toBeCloseTo(-r0.closest[1].values[0],
            12);
        expect(r1.closest[1].values[1]).toBeCloseTo(-r0.closest[1].values[1],
            12);
    });

    it('handles a sheared parallelogram', () => {
        // Axes (1,0) and (1,1); the parallelogram has vertices at
        // (2,1), (0,-1), (-2,-1) and (0,1).
        const g = pgm([0, 0], [1, 0], [1, 1]);
        const result = query.compute(v(0, 4), g);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
    });

    it('agrees with a dense sampling of the parallelogram', () => {
        let seed = 76543210;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 40; ++trial) {
            // Build a right-handed axis pair: DotPerp(a0,a1) > 0.
            const angle = 2 * Math.PI * rand();
            const len0 = 0.4 + 2 * rand();
            const len1 = 0.4 + 2 * rand();
            const spread = 0.3 + 2.4 * rand();
            const a0 = [len0 * Math.cos(angle), len0 * Math.sin(angle)];
            const a1 = [len1 * Math.cos(angle + spread),
                len1 * Math.sin(angle + spread)];
            const g = pgm([2 * rand() - 1, 2 * rand() - 1], a0, a1);
            const p = v(10 * rand() - 5, 10 * rand() - 5);
            const result = query.compute(p, g);
            const brute = bruteForce(p, g);
            expect(result.distance).toBeCloseTo(brute, 6);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 10);
            expect(result.closest[0].values).toEqual(p.values);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistPoint2Parallelogram2.h.
// ---------------------------------------------------------------------------

// A well-conditioned parallelogram: axis[0] = l0*U0 and axis[1] = l1*U1 + m*U0
// for a right-handed orthonormal frame {U0,U1}. Then DotPerp(axis0, axis1) =
// l0*l1 > 0, which is the right-handedness precondition, and the shear m
// keeps the axes non-orthogonal. Both l0 and l1 are bounded away from zero so
// that A(0,0) and A(1,1) (the divisors in GetMinimizer) are well scaled.
const pgmArb = fc.tuple(wellScaledVector(2, -6, 6), rotationFrame(2),
    positive(4, 0.25), positive(4, 0.25), wellScaled(-3, 3))
    .map(([c, U, l0, l1, m]) => Parallelogram2.fromCenterAxis(c,
        [mul(l0, U[0]), add(mul(l1, U[1]), mul(m, U[0]))]));

const pointArb = wellScaledVector(2, -10, 10);

// The (s0,s1) coordinates of X relative to the parallelogram: X = C + s0*A0 +
// s1*A1. Solved with Cramer's rule on the 2x2 system.
function coordinates(g: Parallelogram2, x: Vector): [number, number] {
    const d = sub(x, g.center);
    const det = dotPerp(g.axis[0], g.axis[1]);
    return [dotPerp(d, g.axis[1]) / det, dotPerp(g.axis[0], d) / det];
}

describe('DistPoint2Parallelogram2 verification', () => {
    const query = new DistPoint2Parallelogram2();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(pointArb, pgmArb), ([p, g]) => {
                const r = query.compute(p, g);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                expectVectorClose(r.closest[0], p, 0, 0);
                expectClose(r.distance,
                    length(sub(r.closest[0], r.closest[1])), 1e-9, 1e-9);
                const [s0, s1] = coordinates(g, r.closest[1]);
                expect(Math.abs(s0)).toBeLessThanOrEqual(1 + 1e-9);
                expect(Math.abs(s1)).toBeLessThanOrEqual(1 + 1e-9);
            });
        });

    it('matches a brute-force minimization over the parallelogram', () => {
        check(fc.tuple(pointArb, pgmArb), ([p, g]) => {
            expectClose(query.compute(p, g).distance, bruteForce(p, g), 1e-7,
                1e-7);
        }, 60);
    }, 30000);

    it('returns the query point itself when it is inside', () => {
        check(fc.tuple(pgmArb, wellScaled(-1, 1), wellScaled(-1, 1)),
            ([g, s0, s1]) => {
                const p = add(g.center,
                    add(mul(s0, g.axis[0]), mul(s1, g.axis[1])));
                const r = query.compute(p, g);
                expect(r.distance).toBeLessThanOrEqual(1e-9);
                expectVectorClose(r.closest[1], p, 1e-9, 1e-9);
            });
    });

    it('getMinimizer stays in the domain [-1,1]^2 and reproduces closest[1]',
        () => {
            check(fc.tuple(pointArb, pgmArb), ([p, g]) => {
                const B = Matrix.fromArray(2, 2, [
                    g.axis[0].values[0], g.axis[1].values[0],
                    g.axis[0].values[1], g.axis[1].values[1]]);
                const A = multiplyATB(B, B);
                const Z = mulMatrix(inverse2x2(B).inverse,
                    sub(p, g.center)) as Vector;
                const K = query.getMinimizer(A, Z);
                expect(Math.abs(K.values[0])).toBeLessThanOrEqual(1);
                expect(Math.abs(K.values[1])).toBeLessThanOrEqual(1);
                const x = add(g.center, add(mul(K.values[0], g.axis[0]),
                    mul(K.values[1], g.axis[1])));
                expectVectorClose(query.compute(p, g).closest[1], x, 1e-12,
                    1e-12);
            });
        });

    it('the residual is orthogonal to the free directions at the minimum',
        () => {
            // First-order optimality: the projection of the residual onto any
            // feasible direction at the minimizer cannot decrease the
            // distance. Test with the four axis directions, pushed inward
            // when the minimizer sits on a boundary of [-1,1]^2.
            check(fc.tuple(pointArb, pgmArb), ([p, g]) => {
                const r = query.compute(p, g);
                const [s0, s1] = coordinates(g, r.closest[1]);
                const residual = sub(p, r.closest[1]);
                const step = 1e-6;
                for (const [d0, d1] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
                    const t0 = s0 + d0 * step;
                    const t1 = s1 + d1 * step;
                    if (Math.abs(t0) > 1 || Math.abs(t1) > 1) {
                        continue;   // leaves the domain
                    }
                    const dir = add(mul(d0 * step, g.axis[0]),
                        mul(d1 * step, g.axis[1]));
                    // Moving by 'dir' must not reduce the squared distance
                    // beyond second order: -2*Dot(residual,dir) + |dir|^2 >= 0.
                    expect(-2 * dot(residual, dir) + dot(dir, dir))
                        .toBeGreaterThan(-1e-9);
                }
            });
        });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(pointArb, pgmArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([p, g, R, tr]) => {
            const xf = (x: Vector): Vector => add(tr, Vector.fromArray([
                R[0].values[0] * x.values[0] + R[1].values[0] * x.values[1],
                R[0].values[1] * x.values[0] + R[1].values[1] * x.values[1]]));
            const rot = (x: Vector): Vector => Vector.fromArray([
                R[0].values[0] * x.values[0] + R[1].values[0] * x.values[1],
                R[0].values[1] * x.values[0] + R[1].values[1] * x.values[1]]);
            const moved = Parallelogram2.fromCenterAxis(xf(g.center),
                [rot(g.axis[0]), rot(g.axis[1])]);
            const r0 = query.compute(p, g);
            const r1 = query.compute(xf(p), moved);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(xf(r0.closest[1]), r1.closest[1], 1e-8, 1e-8);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(pointArb, pgmArb), ([p, g]) => {
            const p0 = p.clone();
            const c = g.center.clone();
            const a0 = g.axis[0].clone();
            const a1 = g.axis[1].clone();
            const r = query.compute(p, g);
            expect(p.values).toEqual(p0.values);
            expect(g.center.values).toEqual(c.values);
            expect(g.axis[0].values).toEqual(a0.values);
            expect(g.axis[1].values).toEqual(a1.values);
            // closest[0] must be a copy of the query point, not an alias.
            r.closest[0].values[0] = 999;
            expect(p.values).toEqual(p0.values);
        });
    });
});
