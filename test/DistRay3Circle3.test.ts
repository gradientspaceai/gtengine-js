import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3 } from '../src/DistLine3Circle3.js';
import { DistRay3Circle3 } from '../src/DistRay3Circle3.js';
import { Line, type Line3 } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ray(origin: number[], direction: number[]): Ray {
    return Ray.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

// The exact distance from a point to the circle (a curve, not a disk).
function pointCircleDistance(p: Vector, c: Circle3): number {
    const delta = sub(p, c.center);
    const h = dot(c.normal, delta);
    const inPlane = sub(delta, mul(h, c.normal));
    const radial = length(inPlane);
    const dr = radial - c.radius;
    return Math.sqrt(h * h + dr * dr);
}

// Verify that a reported closest pair is consistent: the circular point is on
// the circle and the pair realizes the reported distance.
function verifyPair(c: Circle3, linear: Vector, circular: Vector,
    distance: number): void {
    const delta = sub(circular, c.center);
    expect(dot(c.normal, delta)).toBeCloseTo(0, 8);
    expect(length(delta)).toBeCloseTo(c.radius, 8);
    expect(length(sub(linear, circular))).toBeCloseTo(distance, 8);
}

// The ray/segment queries delegate to DistLine3Circle3 for the critical
// points of the line-circle distance. Upstream's PDFSection422 computes
// tauHat = sqrt(|(a1*a3)^(2/3) - a3|) where the derivation requires
// tauHat = sqrt(|(a1*a3)^(2/3) - a3| / a2). The port fixes this (issue filed
// upstream), so every trial checks the line solver against a brute-force
// minimum over the whole line instead of skipping unreliable configurations.
function lineSolverIsReliable(line: Line3, c: Circle3): boolean {
    const lineDistance = new DistLine3Circle3().compute(line, c).distance;
    const at = (t: number): number =>
        pointCircleDistance(add(line.origin, mul(t, line.direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = -10 + 20 * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 20 / n;
    for (let pass = 0; pass < 80; ++pass) {
        for (const sign of [1, -1]) {
            const d = at(bt + sign * h);
            if (d < best) {
                best = d;
                bt = bt + sign * h;
            }
        }
        h *= 0.75;
    }
    return Math.abs(lineDistance - best) <= 1e-4;
}

// Brute-force minimum over the ray parameter: a coarse sampling followed by a
// local refinement (the distance function is smooth away from the axis).
function bruteForce(r: Ray, c: Circle3, tmax: number): number {
    const at = (t: number): number =>
        pointCircleDistance(add(r.origin, mul(t, r.direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = tmax * i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = tmax / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.max(0, bt + sign * h);
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

describe('DistRay3Circle3', () => {
    const query = new DistRay3Circle3();
    const unitCircle = circle([0, 0, 0], [0, 0, 1], 2);

    it('uses the ray origin when the ray points away from the circle', () => {
        const r = ray([5, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
        expect(result.circularClosest[0].values[0]).toBeCloseTo(2, 10);
    });

    it('reports zero distance when the ray meets the circle', () => {
        const r = ray([-10, 0, 0], [1, 0, 0]);
        const result = query.compute(r, unitCircle);
        expect(result.distance).toBeCloseTo(0, 10);
        for (let j = 0; j < result.numClosestPairs; ++j) {
            verifyPair(unitCircle, result.linearClosest[j],
                result.circularClosest[j], result.distance);
        }
    });

    it('handles a ray on the circle axis', () => {
        const r = ray([0, 0, 5], [0, 0, 1]);
        const result = query.compute(r, unitCircle);
        expect(result.distance).toBeCloseTo(Math.sqrt(25 + 4), 10);
        expect(result.linearClosest[0].values).toEqual([0, 0, 5]);
        verifyPair(unitCircle, result.linearClosest[0],
            result.circularClosest[0], result.distance);
    });

    it('keeps the line solution when the critical point is on the ray', () => {
        const r = ray([0, 0, -5], [0, 0, 1]);
        const result = query.compute(r, unitCircle);
        // The line meets the plane of the circle at the center, so the ray
        // point closest to the circle is the circle center.
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
    });

    it('handles a degenerate zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const r = ray([3, 0, 0], [1, 0, 0]);
        const result = query.compute(r, c);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 24681012;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        let compared = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const c = circle([2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 0.2 + rand()],
                0.5 + 2 * rand());
            const r = ray([8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            if (length(r.direction) < 1e-3) {
                continue;
            }
            const result = query.compute(r, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
            for (let j = 0; j < result.numClosestPairs; ++j) {
                verifyPair(c, result.linearClosest[j],
                    result.circularClosest[j], result.distance);
            }
            const line = Line.fromOriginDirection(r.origin, r.direction);
            expect(lineSolverIsReliable(line, c)).toBe(true);
            const brute = bruteForce(r, c, 60);
            expect(result.distance).toBeCloseTo(brute, 6);
            ++compared;
        }
        expect(compared).toBeGreaterThan(50);
    });
});
