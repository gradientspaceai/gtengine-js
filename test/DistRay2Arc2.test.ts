import { describe, expect, it } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { DistRay2Arc2 } from '../src/DistRay2Arc2.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

// Build an arc on the circle with the given center and radius, traversing
// counterclockwise from angle a0 to angle a1.
function arc(center: number[], radius: number, a0: number,
    a1: number): Arc2 {
    const c = v(...center);
    const end0 = add(c, v(radius * Math.cos(a0), radius * Math.sin(a0)));
    const end1 = add(c, v(radius * Math.cos(a1), radius * Math.sin(a1)));
    return Arc2.fromCenterRadiusEnds(c, radius, end0, end1);
}

// The exact distance from a point to a ray, computed independently.
function pointRayDistance(p: Vector, r: Ray): number {
    const diff = sub(p, r.origin);
    let t = dot(r.direction, diff) / dot(r.direction, r.direction);
    if (t < 0) {
        t = 0;
    }
    return length(sub(p, add(r.origin, mul(t, r.direction))));
}

// Brute-force minimum distance between the ray and a dense sampling of the
// arc.
function bruteForce(r: Ray, a: Arc2, a0: number, a1: number): number {
    let hi = a1;
    while (hi < a0) {
        hi += 2 * Math.PI;
    }
    const at = (t: number): number => pointRayDistance(add(a.center,
        v(a.radius * Math.cos(t), a.radius * Math.sin(t))), r);
    const n = 4000;
    let best = Number.MAX_VALUE;
    let bt = a0;
    for (let i = 0; i <= n; ++i) {
        const t = a0 + (hi - a0) * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    // Refine locally around the best sample; the distance function is smooth
    // and locally convex near its minimum.
    let h = (hi - a0) / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.min(hi, Math.max(a0, bt + sign * h));
            const d = at(t);
            if (d < best) {
                best = d;
                bt = t;
            }
        }
        h *= 0.75;
    }
    return best;
}

describe('DistRay2Arc2', () => {
    const query = new DistRay2Arc2();
    // The quarter arc of the unit circle in the first quadrant.
    const quarter = arc([0, 0], 1, 0, Math.PI / 2);

    it('computes the distance when the closest circle point is on the arc',
        () => {
            const r = ray([5, 5], [1, 1]);
            const result = query.compute(r, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(5 * Math.SQRT2 - 1, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(Math.SQRT1_2,
                10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(Math.SQRT1_2,
                10);
        });

    it('reports zero distance when the ray meets the arc', () => {
        const r = ray([0, 0], [1, 1]);
        const result = query.compute(r, quarter);
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
    });

    it('falls back to an arc endpoint when the circle point is off the arc',
        () => {
            // The ray is far along the negative x-axis pointing away from
            // the arc; the closest arc point is the endpoint (0,1).
            const r = ray([-4, 1], [-1, 0]);
            const result = query.compute(r, quarter);
            expect(result.numClosestPairs).toBe(1);
            expect(result.distance).toBeCloseTo(4, 10);
            expect(result.closest[0][1].values[0]).toBeCloseTo(0, 10);
            expect(result.closest[0][1].values[1]).toBeCloseTo(1, 10);
        });

    it('reports two closest pairs when the arc endpoints are equidistant',
        () => {
            // The upper half of the unit circle and a ray below it pointing
            // down, so both endpoints (1,0) and (-1,0) are equidistant.
            const upper = arc([0, 0], 1, 0, Math.PI);
            const r = ray([0, -5], [0, -1]);
            const result = query.compute(r, upper);
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(Math.sqrt(26), 10);
            const xs = [result.closest[0][1].values[0],
                result.closest[1][1].values[0]].sort((p, q) => p - q);
            expect(xs[0]).toBeCloseTo(-1, 10);
            expect(xs[1]).toBeCloseTo(1, 10);
        });

    it('handles a degenerate zero-radius arc', () => {
        const point = Arc2.fromCenterRadiusEnds(v(2, 0), 0, v(2, 0), v(2, 0));
        const r = ray([0, 0], [0, 1]);
        const result = query.compute(r, point);
        expect(result.distance).toBeCloseTo(2, 8);
    });

    it('agrees with a dense sampling of the arc on random inputs', () => {
        let seed = 60708090;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 30; ++trial) {
            const radius = 0.5 + 2 * rand();
            const a0 = 2 * Math.PI * rand();
            const a1 = a0 + 0.3 + 5.5 * rand();
            const center = [4 * rand() - 2, 4 * rand() - 2];
            const a = arc(center, radius, a0, a1);
            const r = ray([10 * rand() - 5, 10 * rand() - 5],
                [2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-6) {
                continue;
            }
            const result = query.compute(r, a);
            const brute = bruteForce(r, a, a0, a1);
            expect(result.distance).toBeCloseTo(brute, 7);

            // Every reported pair realizes the reported distance, the ray
            // point is on the ray and the arc point is on the arc.
            for (let j = 0; j < result.numClosestPairs; ++j) {
                const onRay = add(r.origin,
                    mul(result.parameter[j], r.direction));
                expect(result.parameter[j]).toBeGreaterThanOrEqual(0);
                expect(length(sub(result.closest[j][0], onRay)))
                    .toBeLessThan(1e-8);
                expect(length(sub(result.closest[j][1], a.center)))
                    .toBeCloseTo(a.radius, 8);
                expect(length(sub(result.closest[j][0], result.closest[j][1])))
                    .toBeCloseTo(result.distance, 8);
            }
        }
    });
});
