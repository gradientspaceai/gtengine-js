import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3, distLine3Circle3Execute } from '../src/DistLine3Circle3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3Circle3', () => {
    const query = new DistLine3Circle3();
    const unit = circle([0, 0, 0], [0, 0, 1], 1);

    it('reports equidistance for the normal line through the center', () => {
        const result = query.compute(line([0, 0, -5], [0, 0, 1]), unit);
        expect(result.equidistant).toBe(true);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.linearClosest[0].values).toEqual([0, 0, 0]);
        const d = sub(result.circularClosest[0], unit.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(1, 10);
    });

    it('handles a line perpendicular to the plane but off the center', () => {
        // The vertical line through (3,0,z); the closest circle point is
        // (1,0,0) and the closest line point is (3,0,0).
        const result = query.compute(line([3, 0, -7], [0, 0, 1]), unit);
        expect(result.equidistant).toBe(false);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
        expect(result.circularClosest[0].values[0]).toBeCloseTo(1, 10);
    });

    it('reports two closest pairs for a line on the circle axis line',
        () => {
            // A line through the circle center, in the plane of the circle:
            // it meets the circle at two diametrically opposite points.
            const result = query.compute(line([0, 0, 0], [1, 0, 0]), unit);
            expect(result.numClosestPairs).toBe(2);
            expect(result.distance).toBeCloseTo(0, 10);
            expect(result.circularClosest[0].values[0]).toBeCloseTo(-1, 10);
            expect(result.circularClosest[1].values[0]).toBeCloseTo(1, 10);
        });

    it('reports two closest pairs for a symmetric skew line', () => {
        // A line parallel to the x-axis at height z = 2 through the axis of
        // the circle. It is symmetric about the yz-plane, so both (+-1,0,0)
        // are closest.
        const result = query.compute(line([0, 0, 2], [1, 0, 0]), unit);
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('measures a line in the plane of the circle that misses it', () => {
        const result = query.compute(line([0, 3, 0], [1, 0, 0]), unit);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 8);
        expect(result.circularClosest[0].values[1]).toBeCloseTo(1, 8);
    });

    it('reports zero distance for a line through a circle point', () => {
        const result = query.compute(line([1, 0, 0], [0, 0, 1]), unit);
        expect(result.distance).toBeCloseTo(0, 8);
    });

    it('exposes the critical points through the friend-style export', () => {
        const { result, critical } =
            distLine3Circle3Execute(line([3, 0, -7], [0, 0, 1]), unit);
        expect(critical.numPoints).toBe(1);
        expect(result.numClosestPairs).toBe(1);
        expect(critical.distance[0]).toBeCloseTo(result.distance, 12);
        // The critical line point matches the line parameter.
        const ln = line([3, 0, -7], [0, 0, 1]);
        const p = add(ln.origin, mul(critical.parameter[0], ln.direction));
        for (let i = 0; i < 3; ++i) {
            expect(p.values[i]).toBeCloseTo(critical.linearPoint[0].values[i],
                10);
        }
    });

    it('does not require a unit-length line direction', () => {
        const r0 = query.compute(line([3, 0, -7], [0, 0, 1]), unit);
        const r1 = query.compute(line([3, 0, -7], [0, 0, 5]), unit);
        expect(r1.distance).toBeCloseTo(r0.distance, 10);
        for (let i = 0; i < 3; ++i) {
            expect(r1.linearClosest[0].values[i]).toBeCloseTo(
                r0.linearClosest[0].values[i], 8);
        }
    });

    it('agrees with a dense sampling of the circle', () => {
        const rnd = makeRandom(160901);
        const center = v(0.5, -1, 0.25);
        const normal = v(1, 2, -0.5);
        normalize(normal);
        const radius = 1.5;
        const c = Circle3.fromCenterNormalRadius(center, normal, radius);

        const U = v(-normal.values[1], normal.values[0], 0);
        normalize(U);
        const W = cross(normal, U);

        for (let trial = 0; trial < 40; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                // The circle point is on the circle.
                const d = sub(result.circularClosest[j], center);
                expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 6);
                expect(Math.abs(dot(d, normal))).toBeLessThan(1e-6);

                // The line point is on the line.
                const w = sub(result.linearClosest[j], ln.origin);
                const f = sub(w, mul(dot(w, ln.direction)
                    / dot(ln.direction, ln.direction), ln.direction));
                expect(Math.sqrt(dot(f, f))).toBeLessThan(1e-6);

                // The pair realizes the reported distance.
                const e = sub(result.linearClosest[j],
                    result.circularClosest[j]);
                expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 5);
            }

            // No sampled circle point is closer to the line.
            const n = 4000;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, add(mul(radius * Math.cos(t), U),
                    mul(radius * Math.sin(t), W)));
                const w = sub(q, ln.origin);
                const s = dot(w, ln.direction) / dd;
                const f = sub(w, mul(s, ln.direction));
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-5);
        }
    });
});
