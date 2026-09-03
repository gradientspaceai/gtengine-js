import { describe, it, expect } from 'vitest';
import {
    IntrLine3Sphere3TI,
    IntrLine3Sphere3FI,
    intrLine3Sphere3DoQuery,
    defaultIntrLine3Sphere3FIResult
} from '../src/IntrLine3Sphere3.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

const ti = new IntrLine3Sphere3TI();
const fi = new IntrLine3Sphere3FI();

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeLine(origin: Vector, direction: Vector): Line {
    const d = direction.clone();
    normalize(d);
    return Line.fromOriginDirection(origin, d);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The signed distance from the sphere center to the line (nonnegative).
function distanceToLine(line: Line, sphere: Hypersphere): number {
    const diff = sub(sphere.center, line.origin);
    const t = dot(diff, line.direction);
    const closest = add(line.origin, mul(t, line.direction));
    const delta = sub(sphere.center, closest);
    return Math.sqrt(dot(delta, delta));
}

describe('IntrLine3Sphere3TI', () => {
    it('reports a line through the sphere center as intersecting', () => {
        const line = makeLine(v3(-5, 0, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        expect(ti.test(line, sphere).intersect).toBe(true);
    });

    it('reports a tangent line as intersecting', () => {
        const line = makeLine(v3(-5, 2, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        expect(ti.test(line, sphere).intersect).toBe(true);
    });

    it('reports a missing line as not intersecting', () => {
        const line = makeLine(v3(-5, 2.5, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        expect(ti.test(line, sphere).intersect).toBe(false);
    });

    it('ignores the line origin position along the line (infinite extent)', () => {
        // The origin is far past the sphere; the line still intersects.
        const line = makeLine(v3(100, 0, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        expect(ti.test(line, sphere).intersect).toBe(true);
    });

    it('handles a zero-radius sphere on the line', () => {
        const line = makeLine(v3(0, 0, 0), v3(0, 0, 1));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 3), 0);
        expect(ti.test(line, sphere).intersect).toBe(true);
    });
});

describe('IntrLine3Sphere3FI', () => {
    it('computes the two intersection points of a central chord', () => {
        const line = makeLine(v3(-5, 0, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        const result = fi.find(line, sphere);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // P + t*D on the sphere: t = 3 and t = 7.
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(7, 12);
        expect(result.point[0].values).toEqual([-2, 0, 0]);
        expect(result.point[1].values).toEqual([2, 0, 0]);
    });

    it('computes a single point for a tangent line, with a degenerate interval', () => {
        const line = makeLine(v3(-5, 2, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        const result = fi.find(line, sphere);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(5);
        expect(result.parameter[1]).toBe(result.parameter[0]);
        expect(result.point[0].values).toEqual([0, 2, 0]);
        expect(result.point[1].values).toEqual([0, 2, 0]);
    });

    it('leaves the default result for a missing line', () => {
        const line = makeLine(v3(-5, 2.5, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        const result = fi.find(line, sphere);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([0, 0, 0]);
    });

    it('reports a chord of the expected half-length', () => {
        // Offset the line by 1 from the center of a radius-2 sphere. The
        // half-chord is sqrt(4-1) = sqrt(3).
        const line = makeLine(v3(-5, 1, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        const result = fi.find(line, sphere);
        const half = 0.5 * (result.parameter[1] - result.parameter[0]);
        expect(half).toBeCloseTo(Math.sqrt(3), 12);
    });

    it('places the intersection points on the sphere for an oblique line', () => {
        const line = makeLine(v3(1, -3, 2), v3(0.3, 0.9, -0.2));
        const sphere = Hypersphere.fromCenterRadius(v3(0.5, 0.25, 1), 1.75);
        const result = fi.find(line, sphere);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (const p of result.point) {
            const delta = sub(p, sphere.center);
            expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(sphere.radius, 12);
        }
    });

    it('exposes intrLine3Sphere3DoQuery, which does not compute points', () => {
        const line = makeLine(v3(-5, 0, 0), v3(1, 0, 0));
        const sphere = Hypersphere.fromCenterRadius(v3(0, 0, 0), 2);
        const result = defaultIntrLine3Sphere3FIResult();
        intrLine3Sphere3DoQuery(line.origin, line.direction, sphere, result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        // DoQuery leaves the points at their default values.
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });
});

describe('IntrLine3Sphere3 consistency', () => {
    it('agrees between TI and FI and with the distance oracle', () => {
        const rand = makeRandom(20250831);
        for (let trial = 0; trial < 500; ++trial) {
            const origin = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const direction = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(direction, direction) < 1e-8) {
                continue;
            }
            const line = makeLine(origin, direction);
            const sphere = Hypersphere.fromCenterRadius(
                v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2),
                0.25 + 1.5 * rand());

            const tiResult = ti.test(line, sphere);
            const fiResult = fi.find(line, sphere);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            const distance = distanceToLine(line, sphere);
            expect(tiResult.intersect).toBe(distance <= sphere.radius);

            if (fiResult.intersect) {
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const delta = sub(fiResult.point[i], sphere.center);
                    expect(Math.sqrt(dot(delta, delta)))
                        .toBeCloseTo(sphere.radius, 10);
                }
                expect(fiResult.parameter[0])
                    .toBeLessThanOrEqual(fiResult.parameter[1]);
            }
        }
    });
});
