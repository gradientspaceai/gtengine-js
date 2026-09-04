import { describe, expect, it } from 'vitest';
import { DistPointLine } from '../src/DistPointLine.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, nonzero, rotationFrame, seededRandom, vector, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistPointLine', () => {
    const query = new DistPointLine();

    it('computes the perpendicular distance to an axis line', () => {
        const result = query.compute(v(0, 2), line([0, 0], [1, 0]));
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.sqrDistance).toBeCloseTo(4, 12);
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.closest[1].values).toEqual([0, 0]);
    });

    it('reports the input point in closest[0]', () => {
        const point = v(3, -1, 7);
        const result = query.compute(point, line([0, 0, 0], [0, 0, 1]));
        expect(result.closest[0].values).toEqual([3, -1, 7]);
        expect(result.parameter).toBeCloseTo(7, 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(10), 12);
    });

    it('handles a non-unit-length direction (parameter scales)', () => {
        // The line is (1,1) + t*(2,0). The closest point to (5,4) is (5,1),
        // reached at t = 2.
        const result = query.compute(v(5, 4), line([1, 1], [2, 0]));
        expect(result.parameter).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(5, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the line', () => {
        const l = line([1, 2, 3], [1, -1, 2]);
        const point = add(l.origin, mul(-0.75, l.direction));
        const result = query.compute(point, l);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(-0.75, 12);
    });

    it('is translation invariant', () => {
        const point = v(2, -3, 5);
        const l = line([1, 0, -1], [0.3, 0.4, -0.5]);
        const shift = v(10, -20, 30);
        const shifted = Line.fromOriginDirection(add(l.origin, shift),
            l.direction);
        const r0 = query.compute(point, l);
        const r1 = query.compute(add(point, shift), shifted);
        expect(r1.distance).toBeCloseTo(r0.distance, 12);
        expect(r1.parameter).toBeCloseTo(r0.parameter, 12);
    });

    it('produces a closest point whose offset is orthogonal to the line',
        () => {
            let seed = 12345;
            const rand = () => {
                seed = (seed * 1103515245 + 12345) % 2147483648;
                return seed / 2147483648 * 4 - 2;
            };
            for (let trial = 0; trial < 50; ++trial) {
                const point = v(rand(), rand(), rand());
                const l = line([rand(), rand(), rand()],
                    [rand() + 3, rand(), rand()]);
                const result = query.compute(point, l);
                const diff = sub(result.closest[0], result.closest[1]);
                expect(dot(diff, l.direction)).toBeCloseTo(0, 9);
                expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 9);
                expect(result.distance).toBeCloseTo(
                    Math.sqrt(result.sqrDistance), 12);

                // The closest point lies on the line.
                const onLine = add(l.origin, mul(result.parameter,
                    l.direction));
                expect(result.closest[1].values[0]).toBeCloseTo(
                    onLine.values[0], 9);
                expect(result.closest[1].values[1]).toBeCloseTo(
                    onLine.values[1], 9);
                expect(result.closest[1].values[2]).toBeCloseTo(
                    onLine.values[2], 9);

                // No nearby line point is closer.
                for (const dt of [-0.1, -0.01, 0.01, 0.1]) {
                    const other = add(l.origin,
                        mul(result.parameter + dt, l.direction));
                    const od = sub(point, other);
                    expect(dot(od, od)).toBeGreaterThanOrEqual(
                        result.sqrDistance - 1e-12);
                }
            }
        });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistPointLine.ts
// against the upstream header DistPointLine.h.
// ---------------------------------------------------------------------------

/** Apply the orthonormal frame R (as columns) to p. */
function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

/**
 * A 3D line whose direction is well scaled but not unit length (upstream
 * explicitly allows non-unit directions). The length filter keeps
 * Dot(D,D) away from underflow so the division in the port is meaningful.
 */
const nonUnitLine3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

describe('DistPointLine verification', () => {
    const query = new DistPointLine();

    it('result is self consistent (distance, sqrDistance, closest)', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitLine3), ([p, l]) => {
            const r = query.compute(p, l);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.distance, Math.sqrt(dot(diff, diff)), 1e-12, 1e-12);
            // closest[0] is the input point and is a copy, not an alias.
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            // closest[1] is on the line at the reported parameter.
            expectVectorClose(r.closest[1],
                add(l.origin, mul(r.parameter, l.direction)), 1e-12, 1e-12);
        });
    });

    it('the closest line point is the orthogonal projection', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitLine3), ([p, l]) => {
            const resid = sub(p, query.compute(p, l).closest[1]);
            // The residual is formed by cancellation, so the tolerance is
            // relative to the magnitudes entering the dot product.
            const scale = 1 + length(l.direction) * length(resid);
            expect(Math.abs(dot(l.direction, resid)))
                .toBeLessThanOrEqual(1e-9 * scale);
        });
    });

    it('is minimal over sampled line points', () => {
        const rand = seededRandom(0x51d0);
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitLine3), ([p, l]) => {
            const r = query.compute(p, l);
            for (let k = 0; k < 20; ++k) {
                const t = 20 * (rand() - 0.5);
                const diff = sub(p, add(l.origin, mul(t, l.direction)));
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitLine3,
            rotationFrame(3), wellScaledVector(3, -5, 5)), ([p, l, R, tr]) => {
            const r0 = query.compute(p, l);
            const r1 = query.compute(add(rot(R, p), tr),
                Line.fromOriginDirection(add(rot(R, l.origin), tr),
                    rot(R, l.direction)));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.parameter, r1.parameter, 1e-9, 1e-8);
        });
    });

    it('scales the parameter by 1/k when the direction is scaled by k', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitLine3,
            nonzero(-4, 4, 0.25)), ([p, l, k]) => {
            const r0 = query.compute(p, l);
            const r1 = query.compute(p,
                Line.fromOriginDirection(l.origin, mul(k, l.direction)));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.parameter, k * r1.parameter, 1e-9, 1e-8);
        });
    });

    it('works for every dimension the runtime Vector supports', () => {
        check(fc.tuple(fc.integer({ min: 1, max: 5 })).chain(([n]) =>
            fc.tuple(vector(n, -5, 5), vector(n, -5, 5), vector(n, -5, 5))
                .filter(([, , d]) => length(d) > 0.25)),
        ([p, o, d]) => {
            const r = query.compute(p, Line.fromOriginDirection(o, d));
            expect(r.closest[0].size).toBe(p.size);
            expect(r.closest[1].size).toBe(p.size);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
        });
    });
});
