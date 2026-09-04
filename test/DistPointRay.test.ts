import { describe, expect, it } from 'vitest';
import { DistPointRay } from '../src/DistPointRay.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { Line } from '../src/Line.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

describe('DistPointRay', () => {
    const query = new DistPointRay();

    it('uses an interior ray point when the projection is positive', () => {
        const result = query.compute(v(3, 4), ray([0, 0], [1, 0]));
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([3, 0]);
        expect(result.distance).toBeCloseTo(4, 12);
    });

    it('clamps to the ray origin when the projection is nonpositive', () => {
        const result = query.compute(v(-3, 4), ray([0, 0], [1, 0]));
        expect(result.parameter).toBe(0);
        expect(result.closest[1].values).toEqual([0, 0]);
        expect(result.distance).toBeCloseTo(5, 12);
    });

    it('clamps at exactly the origin projection', () => {
        const result = query.compute(v(0, 2), ray([0, 0], [1, 0]));
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(2, 12);
    });

    it('handles a non-unit direction', () => {
        // Ray (1,1) + t*(0,3). The point (4,7) projects to t = 2.
        const result = query.compute(v(4, 7), ray([1, 1], [0, 3]));
        expect(result.parameter).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(7, 12);
        expect(result.distance).toBeCloseTo(3, 12);
    });

    it('reports zero distance for a point on the ray', () => {
        const r = ray([1, 2, 3], [1, -1, 2]);
        const point = add(r.origin, mul(1.5, r.direction));
        const result = query.compute(point, r);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeCloseTo(1.5, 12);
    });

    it('matches a sampled minimum over the ray', () => {
        let seed = 987;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 6 - 3;
        };
        for (let trial = 0; trial < 60; ++trial) {
            const point = v(rand(), rand(), rand());
            const r = ray([rand(), rand(), rand()],
                [rand(), rand() + 2, rand()]);
            const result = query.compute(point, r);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    onRay.values[i], 9);
            }

            let best = Number.MAX_VALUE;
            for (let k = 0; k <= 4000; ++k) {
                const t = k * 0.005;
                const q = add(r.origin, mul(t, r.direction));
                const d = sub(point, q);
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistPointRay.ts
// against the upstream header DistPointRay.h.
// ---------------------------------------------------------------------------

function rot(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < R.length; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

// A 3D ray whose direction is well scaled but not unit length; upstream
// explicitly allows non-unit directions.
const nonUnitRay3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Ray.fromOriginDirection(o, d));

describe('DistPointRay verification', () => {
    const query = new DistPointRay();
    const lineQuery = new DistPointLine();

    it('result is self consistent and the parameter is nonnegative', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitRay3), ([p, ray]) => {
            const r = query.compute(p, ray);
            expect(r.parameter).toBeGreaterThanOrEqual(0);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.distance, Math.sqrt(dot(diff, diff)), 1e-12, 1e-12);
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            expectVectorClose(r.closest[1],
                add(ray.origin, mul(r.parameter, ray.direction)), 1e-12, 1e-12);
        });
    });

    it('agrees with the point-line query exactly when the line parameter is '
        + 'positive, and clamps to the origin otherwise', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitRay3), ([p, ray]) => {
            const line = Line.fromOriginDirection(ray.origin, ray.direction);
            const rl = lineQuery.compute(p, line);
            const rr = query.compute(p, ray);
            if (rl.parameter > 0) {
                expectClose(rr.parameter, rl.parameter, 0, 0);
                expectClose(rr.distance, rl.distance, 0, 0);
            }
            else {
                expect(rr.parameter).toBe(0);
                expectVectorClose(rr.closest[1], ray.origin, 0, 0);
                // Clamping can only increase the distance.
                expect(rr.distance).toBeGreaterThanOrEqual(
                    rl.distance - 1e-9 * (1 + rl.distance));
            }
        });
    });

    it('is minimal over sampled ray points', () => {
        const rand = seededRandom(0x51d1);
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitRay3), ([p, ray]) => {
            const r = query.compute(p, ray);
            for (let k = 0; k < 20; ++k) {
                const t = 12 * rand();
                const diff = sub(p, add(ray.origin, mul(t, ray.direction)));
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is invariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), nonUnitRay3,
            rotationFrame(3), wellScaledVector(3, -5, 5)),
        ([p, ray, R, tr]) => {
            const r0 = query.compute(p, ray);
            const r1 = query.compute(add(rot(R, p), tr),
                Ray.fromOriginDirection(add(rot(R, ray.origin), tr),
                    rot(R, ray.direction)));
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(r0.parameter, r1.parameter, 1e-9, 1e-8);
        });
    });

    it('reports the origin for a degenerate (zero-direction) ray', () => {
        check(wellScaledVector(3, -8, 8), p => {
            const ray = Ray.fromOriginDirection(Vector.fromArray([1, 2, 3]),
                new Vector(3));
            const r = query.compute(p, ray);
            // Dot(D,D) is zero, so the parameter is NaN upstream too; the
            // comparison 'parameter > 0' is false and the origin is used.
            expect(r.parameter).toBe(0);
            expectVectorClose(r.closest[1], ray.origin, 0, 0);
            const diff = sub(p, ray.origin);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
        });
    });
});
