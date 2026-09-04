import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { DistRay3CanonicalBox3 } from '../src/DistRay3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistRay3CanonicalBox3', () => {
    const query = new DistRay3CanonicalBox3();
    const unitBox = CanonicalBox.fromExtent(v(1, 1, 1));

    it('reports zero distance for a ray that enters the box', () => {
        const result = query.compute(ray([-5, 0, 0], [1, 0, 0]), unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.parameter).toBeGreaterThanOrEqual(0);
    });

    it('reports zero distance for a ray originating inside the box', () => {
        const result = query.compute(ray([0, 0, 0], [1, 2, 3]), unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('uses the ray origin when the ray points away from the box', () => {
        // The origin is at (5,0,0) and the ray points to +x.
        const result = query.compute(ray([5, 0, 0], [1, 0, 0]), unitBox);
        expect(result.parameter).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0].values).toEqual([5, 0, 0]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('matches the point-box query when the ray points away', () => {
        const pbQuery = new DistPointCanonicalBox();
        const origin = v(3, -4, 5);
        const result = query.compute(
            Ray.fromOriginDirection(origin, v(1, -1, 1)), unitBox);
        const pbResult = pbQuery.compute(origin, unitBox);
        expect(result.distance).toBeCloseTo(pbResult.distance, 12);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(
                pbResult.closest[1].values[i], 12);
        }
    });

    it('matches the line query when the closest line point is on the ray',
        () => {
            const lbQuery = new DistLine3CanonicalBox3();
            const r = ray([-5, 3, 0], [1, 0, 0]);
            const result = query.compute(r, unitBox);
            const lbResult = lbQuery.compute(
                Line.fromOriginDirection(r.origin, r.direction), unitBox);
            expect(lbResult.parameter).toBeGreaterThanOrEqual(0);
            expect(result.distance).toBeCloseTo(lbResult.distance, 12);
        });

    it('handles a degenerate box with zero extents', () => {
        const box = CanonicalBox.fromExtent(v(0, 0, 0));
        const result = query.compute(ray([0, 3, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the ray and box', () => {
        const rnd = makeRandom(51515);
        const extent = v(1, 0.75, 1.5);
        const box = CanonicalBox.fromExtent(extent);

        for (let trial = 0; trial < 25; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r = Ray.fromOriginDirection(origin, dir);
            const result = query.compute(r, box);

            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onRay.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
                expect(Math.abs(result.closest[1].values[i]))
                    .toBeLessThanOrEqual(extent.values[i] + 1e-8);
            }
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled (ray point, box point) pair is closer.
            const nb = 10;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= nb; ++i) {
                for (let j = 0; j <= nb; ++j) {
                    for (let k = 0; k <= nb; ++k) {
                        const q = v(
                            (2 * i / nb - 1) * extent.values[0],
                            (2 * j / nb - 1) * extent.values[1],
                            (2 * k / nb - 1) * extent.values[2]);
                        // The closest ray point to q, clamped to t >= 0.
                        const w = sub(q, r.origin);
                        let t = dot(w, r.direction)
                            / dot(r.direction, r.direction);
                        if (t < 0) {
                            t = 0;
                        }
                        const f = sub(w, mul(t, r.direction));
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistRay3CanonicalBox3.h.
// ---------------------------------------------------------------------------

describe('DistRay3CanonicalBox3 verification', () => {
    const query = new DistRay3CanonicalBox3();
    const lineQuery = new DistLine3CanonicalBox3();
    const pointQuery = new DistPointCanonicalBox();

    const boxArb = fc.array(finite(0, 4), { minLength: 3, maxLength: 3 })
        .map(e => CanonicalBox.fromExtent(Vector.fromArray(e)));

    const rayArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Ray.fromOriginDirection(o, d));

    // Independent closed-form distance from a point to a solid canonical box.
    function pointBoxDistance(p: Vector, b: CanonicalBox): number {
        let sum = 0;
        for (let i = 0; i < 3; ++i) {
            const over = Math.abs(p.values[i]) - b.extent.values[i];
            if (over > 0) { sum += over * over; }
        }
        return Math.sqrt(sum);
    }

    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rayArb, boxArb), ([ray, b]) => {
            const r = query.compute(ray, b);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expect(r.parameter).toBeGreaterThanOrEqual(0);
            expectVectorClose(r.closest[0],
                add(ray.origin, mul(r.parameter, ray.direction)), 1e-9, 1e-9);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-8, 1e-8);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(r.closest[1].values[i]))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
            }
        });
    });

    it('matches an independent convex minimization along the ray', () => {
        check(fc.tuple(rayArb, boxArb), ([ray, b]) => {
            const r = query.compute(ray, b);
            const best = ternaryMin(t => pointBoxDistance(
                add(ray.origin, mul(t, ray.direction)), b), 0, 200);
            // See the DistLine3AlignedBox3 note: the incremental
            // canonical-box accumulation loses about half the mantissa for a
            // grazing line, so the tolerance is absolute rather than tight.
            expectClose(r.distance, best, 2e-6, 1e-9);
        }, 100);
    });

    it('clamps the line parameter to the ray', () => {
        check(fc.tuple(rayArb, boxArb), ([ray, b]) => {
            const line = Line.fromOriginDirection(ray.origin, ray.direction);
            const rl = lineQuery.compute(line, b);
            const rr = query.compute(ray, b);
            if (rl.parameter >= 0) {
                expectClose(rr.distance, rl.distance, 1e-12, 1e-12);
                expectClose(rr.parameter, rl.parameter, 1e-12, 1e-12);
                expectVectorClose(rr.closest[1], rl.closest[1], 1e-12, 1e-12);
            } else {
                expect(rr.parameter).toBe(0);
                expectVectorClose(rr.closest[0], ray.origin, 1e-12, 1e-12);
                const rp = pointQuery.compute(ray.origin, b);
                expectClose(rr.distance, rp.distance, 1e-12, 1e-12);
                expectVectorClose(rr.closest[1], rp.closest[1], 1e-12, 1e-12);
            }
            expect(rr.distance).toBeGreaterThanOrEqual(rl.distance - 1e-12);
        });
    });

    it('reports zero distance for a ray starting inside the box', () => {
        check(fc.tuple(boxArb, fc.array(finite(-1, 1),
            { minLength: 3, maxLength: 3 }), unitVector(3)),
            ([b, u, dir]) => {
                const p = Vector.fromArray([u[0] * b.extent.values[0],
                    u[1] * b.extent.values[1], u[2] * b.extent.values[2]]);
                const r = query.compute(Ray.fromOriginDirection(p, dir), b);
                expectClose(r.distance, 0, 1e-9, 1e-9);
            });
    });

    it('is invariant under reflection of any coordinate axis', () => {
        // The canonical box is symmetric about every coordinate plane.
        check(fc.tuple(rayArb, boxArb, fc.nat(2)), ([ray, b, k]) => {
            const flip = (p: Vector): Vector => {
                const q = p.clone();
                q.values[k] = -q.values[k];
                return q;
            };
            const r0 = query.compute(ray, b);
            const r1 = query.compute(Ray.fromOriginDirection(flip(ray.origin),
                flip(ray.direction)), b);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
        });
    });
});
