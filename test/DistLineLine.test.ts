import { describe, expect, it } from 'vitest';
import { DistLineLine } from '../src/DistLineLine.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistLineLine', () => {
    const query = new DistLineLine();

    it('computes the distance between perpendicular skew lines', () => {
        // The x-axis and the line (0,0,5)+s*(0,1,0) are skew with distance 5.
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            line([0, 0, 5], [0, 1, 0]));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
    });

    it('computes zero distance for intersecting lines', () => {
        const result = query.compute(line([0, 0, 0], [1, 1, 0]),
            line([2, 0, 0], [-1, 1, 0]));
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[0].values[1]).toBeCloseTo(1, 10);
    });

    it('handles parallel lines by choosing s1 = 0', () => {
        const result = query.compute(line([0, 0, 0], [1, 0, 0]),
            line([3, 4, 0], [2, 0, 0]));
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.parameter[1]).toBe(0);
        expect(result.closest[1].values).toEqual([3, 4, 0]);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
    });

    it('handles coincident lines (zero distance, arbitrary pair)', () => {
        const result = query.compute(line([1, 2, 3], [1, 1, 1]),
            line([2, 3, 4], [2, 2, 2]));
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('is symmetric in the distance', () => {
        const l0 = line([1, -2, 3], [0.5, 1.5, -2]);
        const l1 = line([-4, 0, 1], [2, -1, 0.25]);
        const a = query.compute(l0, l1);
        const b = query.compute(l1, l0);
        expect(b.distance).toBeCloseTo(a.distance, 10);
        expect(b.parameter[0]).toBeCloseTo(a.parameter[1], 10);
        expect(b.parameter[1]).toBeCloseTo(a.parameter[0], 10);
    });

    it('produces a connecting segment orthogonal to both directions', () => {
        let seed = 5150;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const l0 = line([rand(), rand(), rand()],
                [rand() + 4, rand(), rand()]);
            const l1 = line([rand(), rand(), rand()],
                [rand(), rand() + 4, rand()]);
            const result = query.compute(l0, l1);

            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, l0.direction)).toBeCloseTo(0, 8);
            expect(dot(diff, l1.direction)).toBeCloseTo(0, 8);
            expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 9);

            // The closest points lie on their lines.
            const c0 = add(l0.origin, mul(result.parameter[0], l0.direction));
            const c1 = add(l1.origin, mul(result.parameter[1], l1.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[0].values[i]).toBeCloseTo(c0.values[i],
                    9);
                expect(result.closest[1].values[i]).toBeCloseTo(c1.values[i],
                    9);
            }

            // Perturbations do not reduce the squared distance.
            for (const ds of [-0.05, 0.05]) {
                const a = add(l0.origin,
                    mul(result.parameter[0] + ds, l0.direction));
                const d = sub(a, result.closest[1]);
                expect(dot(d, d)).toBeGreaterThanOrEqual(
                    result.sqrDistance - 1e-9);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistLineLine.ts
// against the upstream header DistLineLine.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

// A 3D line with a well-scaled, deliberately non-unit direction.
const nonUnitLine3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

// Two lines that are well away from parallel, so that the determinant
// a00*a11 - a01^2 carries significant digits and the closest-point
// parameters are well conditioned.
const skewPair = fc.tuple(nonUnitLine3, nonUnitLine3)
    .filter(([l0, l1]) =>
        length(cross(l0.direction, l1.direction))
        > 0.2 * length(l0.direction) * length(l1.direction));

describe('DistLineLine verification', () => {
    const query = new DistLineLine();
    const pointLine = new DistPointLine();

    it('result is self consistent', () => {
        check(fc.tuple(nonUnitLine3, nonUnitLine3), ([l0, l1]) => {
            const r = query.compute(l0, l1);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
            expectVectorClose(r.closest[0],
                add(l0.origin, mul(r.parameter[0], l0.direction)), 1e-12, 1e-12);
            expectVectorClose(r.closest[1],
                add(l1.origin, mul(r.parameter[1], l1.direction)), 1e-12, 1e-12);
        });
    });

    it('matches the closed-form skew-line distance |Dot(P0-P1,N)|/|N|', () => {
        check(skewPair, ([l0, l1]) => {
            const N = cross(l0.direction, l1.direction);
            const expected =
                Math.abs(dot(sub(l0.origin, l1.origin), N)) / length(N);
            expectClose(query.compute(l0, l1).distance, expected, 1e-9, 1e-9);
        });
    });

    it('produces a segment perpendicular to both directions', () => {
        check(skewPair, ([l0, l1]) => {
            const r = query.compute(l0, l1);
            const w = sub(r.closest[0], r.closest[1]);
            const scale = 1 + length(w);
            expect(Math.abs(dot(l0.direction, w)))
                .toBeLessThanOrEqual(1e-8 * scale * length(l0.direction));
            expect(Math.abs(dot(l1.direction, w)))
                .toBeLessThanOrEqual(1e-8 * scale * length(l1.direction));
        });
    });

    it('is symmetric under argument swap', () => {
        check(fc.tuple(nonUnitLine3, nonUnitLine3), ([l0, l1]) => {
            const r0 = query.compute(l0, l1);
            const r1 = query.compute(l1, l0);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    // The pair is restricted to skewPair: for nearly parallel lines the
    // determinant and both numerators cancel to rounding noise, so upstream's
    // s0 = (a01*b1 - a11*b0)/det is garbage of order one. The closest points
    // still lie on their lines, but they need not be near the minimum. That
    // conditioning limit is inherent in the upstream formula, not a port
    // defect, and is why DistSegmentSegment offers computeRobust.
    it('is minimal over sampled point pairs', () => {
        const rand = seededRandom(0x51d5);
        check(skewPair, ([l0, l1]) => {
            const r = query.compute(l0, l1);
            for (let k = 0; k < 20; ++k) {
                const p = add(l0.origin, mul(20 * (rand() - 0.5), l0.direction));
                const q = add(l1.origin, mul(20 * (rand() - 0.5), l1.direction));
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(skewPair, rotationFrame(3), wellScaledVector(3, -5, 5)),
            ([[l0, l1], R, tr]) => {
                const move = (l: Line): Line => Line.fromOriginDirection(
                    add(rot(R, l.origin), tr), rot(R, l.direction));
                expectClose(query.compute(l0, l1).distance,
                    query.compute(move(l0), move(l1)).distance, 1e-8, 1e-8);
            });
    });

    it('falls back to a point-line distance for parallel lines', () => {
        // The scale factor is a power of two so that a01 = -k*a00 and
        // a11 = k^2*a00 hold exactly in binary64 and the determinant is
        // exactly zero. With a generic k, rounding leaves a tiny positive
        // determinant and upstream takes the (ill-conditioned) nonparallel
        // branch instead.
        check(fc.tuple(nonUnitLine3, wellScaledVector(3, -8, 8),
            fc.constantFrom(1, -1, 2, -2, 0.5, -0.5, 4, -4, 0.25, -0.25)),
        ([l0, o1, k]) => {
            const l1 = Line.fromOriginDirection(o1, mul(k, l0.direction));
            const r = query.compute(l0, l1);
            // Upstream selects s1 = 0, so the closest pair is the projection
            // of the second origin onto the first line.
            expect(r.parameter[1]).toBe(0);
            expectVectorClose(r.closest[1], o1, 0, 0);
            expectClose(r.distance, pointLine.compute(o1, l0).distance, 1e-9,
                1e-9);
        });
    });
});
