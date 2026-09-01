import { describe, it, expect } from 'vitest';
import { Hypersphere } from '../src/Hypersphere';
import { Line } from '../src/Line';
import { Vector, add, sub, mul, dot, length, normalize } from '../src/Vector';
import {
    IntrLine2Circle2TI,
    IntrLine2Circle2FI
} from '../src/IntrLine2Circle2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(px: number, py: number, dx: number, dy: number): Line {
    const d = v2(dx, dy);
    normalize(d);
    return Line.fromOriginDirection(v2(px, py), d);
}

function circle(cx: number, cy: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v2(cx, cy), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrLine2Circle2', () => {
    const ti = new IntrLine2Circle2TI();
    const fi = new IntrLine2Circle2FI();
    const unit = circle(0, 0, 1);

    it('finds two intersections for a line through the center', () => {
        const result = fi.find(line(0, 0, 1, 0), unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('finds the tangent contact as a single point', () => {
        // The line y = 1 is tangent to the unit circle at (0,1).
        const result = fi.find(line(-5, 1, 1, 0), unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(result.parameter[1]);
        expect(result.parameter[0]).toBeCloseTo(5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a line that misses the circle', () => {
        const result = fi.find(line(-5, 2, 1, 0), unit);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(line(-5, 2, 1, 0), unit).intersect).toBe(false);
    });

    it('places the intersection points on both the line and the circle', () => {
        const rnd = makeRandom(2024);
        let numTwo = 0;
        for (let k = 0; k < 300; ++k) {
            const l = line(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            if (!Number.isFinite(l.direction.values[0])) {
                continue;
            }
            const c = circle(rnd() * 4 - 2, rnd() * 4 - 2, 0.2 + rnd() * 2);
            const result = fi.find(l, c);
            for (let i = 0; i < result.numIntersections; ++i) {
                const p = result.point[i];
                expect(length(sub(p, c.center))).toBeCloseTo(c.radius, 9);
                const q = add(l.origin, mul(result.parameter[i], l.direction));
                expect(p.values[0]).toBeCloseTo(q.values[0], 9);
                expect(p.values[1]).toBeCloseTo(q.values[1], 9);
            }
            if (result.numIntersections === 2) {
                expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
                ++numTwo;
            }
        }
        expect(numTwo).toBeGreaterThan(0);
    });

    it('agrees with the TI query on random configurations', () => {
        const rnd = makeRandom(99);
        let numHit = 0, numMiss = 0;
        for (let k = 0; k < 500; ++k) {
            const l = line(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const c = circle(rnd() * 4 - 2, rnd() * 4 - 2, 0.2 + rnd() * 2);
            const tiResult = ti.test(l, c).intersect;
            const fiResult = fi.find(l, c).intersect;
            // Both are decided by the sign of radius^2 - distance^2, so they
            // agree except possibly at exact tangency.
            expect(tiResult).toBe(fiResult);
            if (tiResult) {
                ++numHit;
            } else {
                ++numMiss;
            }
        }
        expect(numHit).toBeGreaterThan(0);
        expect(numMiss).toBeGreaterThan(0);
    });

    it('cross-checks the TI query against the point-line distance', () => {
        const rnd = makeRandom(31337);
        for (let k = 0; k < 300; ++k) {
            const l = line(rnd() * 6 - 3, rnd() * 6 - 3,
                rnd() * 2 - 1, rnd() * 2 - 1);
            const c = circle(rnd() * 4 - 2, rnd() * 4 - 2, 0.2 + rnd() * 2);
            const diff = sub(c.center, l.origin);
            const t = dot(diff, l.direction);
            const closest = add(l.origin, mul(t, l.direction));
            const distance = length(sub(c.center, closest));
            expect(ti.test(l, c).intersect).toBe(distance <= c.radius);
        }
    });

    it('handles a zero-radius circle (a point on the line)', () => {
        // The origin is on the x axis, so the degenerate circle is hit.
        const result = fi.find(line(-3, 0, 1, 0), circle(0, 0, 0));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(3, 12);

        // A point off the line is missed.
        expect(fi.find(line(-3, 0, 1, 0), circle(0, 1, 0)).intersect)
            .toBe(false);
    });

    it('handles a circle containing the line origin', () => {
        const result = fi.find(line(0.25, 0.25, 1, 0), unit);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeLessThan(0);
        expect(result.parameter[1]).toBeGreaterThan(0);
    });
});
