import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistLine3Circle3 } from '../src/DistLine3Circle3.js';
import { DistSegment3Circle3 } from '../src/DistSegment3Circle3.js';
import { Line, type Line3 } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
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

// Brute-force minimum over the segment parameter: a dense sampling followed
// by a local refinement.
function bruteForce(s: Segment, c: Circle3): number {
    const direction = sub(s.p[1], s.p[0]);
    const at = (t: number): number =>
        pointCircleDistance(add(s.p[0], mul(t, direction)), c);
    const n = 8000;
    let best = Number.MAX_VALUE;
    let bt = 0;
    for (let i = 0; i <= n; ++i) {
        const t = i / n;
        const d = at(t);
        if (d < best) {
            best = d;
            bt = t;
        }
    }
    let h = 1 / n;
    for (let pass = 0; pass < 120; ++pass) {
        for (const sign of [1, -1]) {
            const t = Math.min(1, Math.max(0, bt + sign * h));
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

describe('DistSegment3Circle3', () => {
    const query = new DistSegment3Circle3();
    const unitCircle = circle([0, 0, 0], [0, 0, 1], 2);

    it('clamps to the first endpoint', () => {
        const s = segment([5, 0, 0], [10, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
    });

    it('clamps to the second endpoint', () => {
        const s = segment([10, 0, 0], [5, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.linearClosest[0].values).toEqual([5, 0, 0]);
    });

    it('reports zero distance when the segment meets the circle', () => {
        const s = segment([-5, 0, 0], [5, 0, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.distance).toBeCloseTo(0, 10);
        for (let j = 0; j < result.numClosestPairs; ++j) {
            verifyPair(unitCircle, result.linearClosest[j],
                result.circularClosest[j], result.distance);
        }
    });

    it('keeps the interior solution when the critical point is on the segment',
        () => {
            const s = segment([0, 0, -5], [0, 0, 5]);
            const result = query.compute(s, unitCircle);
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.linearClosest[0].values[2]).toBeCloseTo(0, 10);
        });

    it('handles a nearly degenerate segment', () => {
        // A zero-length segment is not a valid input (the algorithm divides
        // by Dot(M,M) with M = P1 - P0), but a very short segment behaves
        // like a point query.
        const s = segment([5, 0, 0], [5, 1e-9, 0]);
        const result = query.compute(s, unitCircle);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('handles a degenerate zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const s = segment([3, 0, 0], [6, 0, 0]);
        const result = query.compute(s, c);
        expect(result.distance).toBeCloseTo(3, 8);
    });

    it('agrees with a brute-force sampling on random inputs', () => {
        let seed = 97531864;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        let compared = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const c = circle([2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1],
                [2 * rand() - 1, 2 * rand() - 1, 0.2 + rand()],
                0.5 + 2 * rand());
            const s = segment(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4]);
            const result = query.compute(s, c);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
            for (let j = 0; j < result.numClosestPairs; ++j) {
                verifyPair(c, result.linearClosest[j],
                    result.circularClosest[j], result.distance);
            }
            const line = Line.fromOriginDirection(s.p[0],
                sub(s.p[1], s.p[0]));
            expect(lineSolverIsReliable(line, c)).toBe(true);
            const brute = bruteForce(s, c);
            expect(result.distance).toBeCloseTo(brute, 6);
            ++compared;
        }
        expect(compared).toBeGreaterThan(50);
    });
});
