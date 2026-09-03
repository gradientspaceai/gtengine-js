import { describe, expect, it } from 'vitest';
import { DistPoint2Parallelogram2 } from '../src/DistPoint2Parallelogram2.js';
import { Parallelogram2 } from '../src/Parallelogram2.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';

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
