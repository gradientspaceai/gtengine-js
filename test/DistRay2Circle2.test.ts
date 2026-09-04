import { describe, expect, it } from 'vitest';
import { DistPoint2Circle2 } from '../src/DistPoint2Circle2.js';
import { DistRay2Circle2 } from '../src/DistRay2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Ray } from '../src/Ray.js';
import { DistLine2Circle2 } from '../src/DistLine2Circle2.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistRay2Circle2', () => {
    const query = new DistRay2Circle2();
    const unit = circle([0, 0], 1);

    it('reports two closest pairs when the ray crosses the circle', () => {
        const result = query.compute(ray([-5, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[0][1].values[0]).toBeCloseTo(-1, 10);
        expect(result.closest[1][1].values[0]).toBeCloseTo(1, 10);
    });

    it('drops the t0-point when the ray origin is inside the circle', () => {
        const result = query.compute(ray([0, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(1, 10);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('uses the ray origin when the ray points away from the circle', () => {
        // The ray origin at (5,0) pointing to +x; the whole line crosses the
        // circle only at negative parameters.
        const result = query.compute(ray([5, 0], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.parameter[0]).toBe(0);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0][0].values).toEqual([5, 0]);
        expect(result.closest[0][1].values[0]).toBeCloseTo(1, 10);
    });

    it('matches the point-circle query when the ray points away', () => {
        const pcQuery = new DistPoint2Circle2();
        const origin = v(3, 4);
        const result = query.compute(
            Ray.fromOriginDirection(origin, v(1, 1)), unit);
        const pcResult = pcQuery.compute(origin, unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(pcResult.distance, 12);
        for (let i = 0; i < 2; ++i) {
            expect(result.closest[0][1].values[i]).toBeCloseTo(
                pcResult.closest[1].values[i], 12);
        }
    });

    it('handles a line that misses the circle with a positive parameter',
        () => {
            const result = query.compute(ray([-5, 3], [1, 0]), unit);
            expect(result.numClosestPairs).toBe(1);
            expect(result.parameter[0]).toBeGreaterThan(0);
            expect(result.distance).toBeCloseTo(2, 10);
        });

    it('handles a tangent ray', () => {
        const result = query.compute(ray([-5, 1], [1, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('agrees with a dense sampling of the ray and circle', () => {
        const rnd = makeRandom(70001);
        const center = v(0.5, -0.75);
        const radius = 1.25;
        const c = Hypersphere.fromCenterRadius(center, radius);

        for (let trial = 0; trial < 60; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const r = Ray.fromOriginDirection(origin, dir);
            const result = query.compute(r, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                expect(result.parameter[j]).toBeGreaterThanOrEqual(-1e-12);
                const p = add(r.origin, mul(result.parameter[j], r.direction));
                for (let i = 0; i < 2; ++i) {
                    expect(p.values[i]).toBeCloseTo(
                        result.closest[j][0].values[i], 7);
                }
                const d = sub(result.closest[j][1], center);
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 7);
                const e = sub(result.closest[j][0], result.closest[j][1]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 6);
            }

            // No sampled (ray point, circle point) pair is closer.
            const n = 3000;
            const dd = dot(r.direction, r.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, v(radius * Math.cos(t),
                    radius * Math.sin(t)));
                const w = sub(q, r.origin);
                let s = dot(w, r.direction) / dd;
                if (s < 0) {
                    s = 0;
                }
                const f = sub(w, mul(s, r.direction));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistRay2Circle2.h.
// ---------------------------------------------------------------------------

describe('DistRay2Circle2 verification', () => {
    const query = new DistRay2Circle2();
    const lineQuery = new DistLine2Circle2();

    const circleArb = fc.tuple(wellScaledVector(2, -5, 5), finite(0.25, 4))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    const rayArb = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2))
        .map(([o, d]) => Ray.fromOriginDirection(o, d));

    // Distance from a point to the circle (a curve, not a disk).
    function pointCircleDistance(p: Vector, c: Hypersphere): number {
        return Math.abs(length(sub(p, c.center)) - c.radius);
    }

    // The distance along the ray is |P(t)-C| - r in absolute value; |P(t)-C|
    // is convex so the composite has at most two local minima. A dense scan
    // followed by a golden-section refinement therefore finds the global
    // minimum. Beyond tMax the ray is monotonically leaving the circle.
    function bruteForce(ray: Ray, c: Hypersphere): number {
        const f = (t: number): number => pointCircleDistance(
            add(ray.origin, mul(t, ray.direction)), c);
        const tMax = 2 * (length(sub(ray.origin, c.center)) + c.radius) + 1;
        const n = 20000;
        let bestI = 0;
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i <= n; ++i) {
            const y = f((i / n) * tMax);
            if (y < best) { best = y; bestI = i; }
        }
        const h = tMax / n;
        let lo = Math.max((bestI / n) * tMax - h, 0);
        let hi = lo + 2 * h;
        const phi = (Math.sqrt(5) - 1) / 2;
        for (let i = 0; i < 200; ++i) {
            const m0 = hi - phi * (hi - lo);
            const m1 = lo + phi * (hi - lo);
            if (f(m0) <= f(m1)) { hi = m1; } else { lo = m0; }
        }
        return Math.min(best, f(0.5 * (lo + hi)));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rayArb, circleArb), ([ray, circle]) => {
            const r = query.compute(ray, circle);
            expect(r.numClosestPairs === 1 || r.numClosestPairs === 2)
                .toBe(true);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            for (let j = 0; j < r.numClosestPairs; ++j) {
                expect(r.parameter[j]).toBeGreaterThanOrEqual(0);
                expectVectorClose(r.closest[j][0],
                    add(ray.origin, mul(r.parameter[j], ray.direction)),
                    1e-8, 1e-8);
                expectClose(length(sub(r.closest[j][1], circle.center)),
                    circle.radius, 1e-8, 1e-8);
                expectClose(length(sub(r.closest[j][0], r.closest[j][1])),
                    r.distance, 1e-8, 1e-8);
            }
        });
    });

    it('matches an independent minimization along the ray', () => {
        check(fc.tuple(rayArb, circleArb), ([ray, circle]) => {
            const r = query.compute(ray, circle);
            expectClose(r.distance, bruteForce(ray, circle), 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('is at least the line-circle distance', () => {
        check(fc.tuple(rayArb, circleArb), ([ray, circle]) => {
            const line = Line.fromOriginDirection(ray.origin, ray.direction);
            const rl = lineQuery.compute(line, circle);
            const rr = query.compute(ray, circle);
            expect(rr.distance).toBeGreaterThanOrEqual(rl.distance - 1e-12);
            // When every line closest point is on the ray the two agree.
            if (rl.numClosestPairs === 2 && rl.parameter[0] >= 0) {
                expectClose(rr.distance, rl.distance, 1e-12, 1e-12);
                expect(rr.numClosestPairs).toBe(2);
            }
            if (rl.numClosestPairs === 1 && rl.parameter[0] >= 0) {
                expectClose(rr.distance, rl.distance, 1e-12, 1e-12);
            }
        });
    });

    it('falls back to the origin when the line closest point is behind',
        () => {
            check(fc.tuple(rayArb, circleArb), ([ray, circle]) => {
                const line = Line.fromOriginDirection(ray.origin,
                    ray.direction);
                const rl = lineQuery.compute(line, circle);
                if (rl.numClosestPairs !== 1 || rl.parameter[0] >= 0) {
                    return;
                }
                const rr = query.compute(ray, circle);
                expect(rr.numClosestPairs).toBe(1);
                expect(rr.parameter[0]).toBe(0);
                expectClose(rr.distance,
                    pointCircleDistance(ray.origin, circle), 1e-12, 1e-12);
            });
        });

    it('drops only the negative intersection when the origin is inside',
        () => {
            check(fc.tuple(circleArb, unitVector(2), finite(0, 0.95),
                finite(-Math.PI, Math.PI)),
                ([circle, dir, frac, angle]) => {
                    const inside = add(circle.center,
                        v(frac * circle.radius * Math.cos(angle),
                            frac * circle.radius * Math.sin(angle)));
                    const ray = Ray.fromOriginDirection(inside, dir);
                    const r = query.compute(ray, circle);
                    // The ray leaves the disk through exactly one point.
                    expect(r.numClosestPairs).toBe(1);
                    expectClose(r.distance, 0, 1e-9, 1e-9);
                    expect(r.parameter[0]).toBeGreaterThanOrEqual(0);
                    expectClose(length(sub(r.closest[0][1], circle.center)),
                        circle.radius, 1e-9, 1e-9);
                });
        });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(rayArb, circleArb, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([ray, circle, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                v(frame[0].values[0] * p.values[0]
                    + frame[1].values[0] * p.values[1],
                    frame[0].values[1] * p.values[0]
                    + frame[1].values[1] * p.values[1]);
            const movedRay = Ray.fromOriginDirection(
                add(shift, rot(ray.origin)), rot(ray.direction));
            const movedCircle = Hypersphere.fromCenterRadius(
                add(shift, rot(circle.center)), circle.radius);
            const r0 = query.compute(ray, circle);
            const r1 = query.compute(movedRay, movedCircle);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expect(r0.numClosestPairs).toBe(r1.numClosestPairs);
        });
    });
});
