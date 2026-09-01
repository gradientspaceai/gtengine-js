import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { Hypersphere } from '../src/Hypersphere';
import { Vector, sub, length } from '../src/Vector';
import { IntrCircle2Arc2FI } from '../src/IntrCircle2Arc2';

const INT32_MAX = 2147483647;

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A point on the unit circle at the given angle in degrees. The four axis
// directions are exact so that endpoint equality tests are exact.
function pt(deg: number): Vector {
    switch (((deg % 360) + 360) % 360) {
        case 0: return v2(1, 0);
        case 90: return v2(0, 1);
        case 180: return v2(-1, 0);
        case 270: return v2(0, -1);
        default: {
            const a = (deg * Math.PI) / 180;
            return v2(Math.cos(a), Math.sin(a));
        }
    }
}

// The arc of the unit circle from deg0 counterclockwise to deg1.
function arc(deg0: number, deg1: number): Arc2 {
    return Arc2.fromCenterRadiusEnds(v2(0, 0), 1, pt(deg0), pt(deg1));
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

describe('IntrCircle2Arc2', () => {
    const fi = new IntrCircle2Arc2FI();

    it('reports the arc when the circle contains it (cocircular)', () => {
        const a = arc(0, 90);
        const result = fi.find(circle(0, 0, 1), a);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.arc.equals(a)).toBe(true);
    });

    it('reports no intersection for disjoint circles', () => {
        const result = fi.find(circle(5, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection when a nested circle misses', () => {
        const result = fi.find(circle(0, 0, 0.5), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds two points when both circle crossings are on the arc', () => {
        // The circle of radius 1 centered at (1,1) meets the unit circle at
        // (1,0) and (0,1), both endpoints of the first-quadrant arc.
        const result = fi.find(circle(1, 1, 1), arc(0, 90));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(length(result.point[i])).toBeCloseTo(1, 12);
        }
    });

    it('finds one point when only one crossing is on the arc', () => {
        // The unit circle centered at (1,1) meets the unit circle at (1,0)
        // and (0,1). The arc from 90 to 270 degrees contains only (0,1).
        const result = fi.find(circle(1, 1, 1), arc(90, 270));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('finds the tangent contact when it lies on the arc', () => {
        // The circle of radius 1 centered at (2,0) is tangent to the unit
        // circle at (1,0), which is the start of the first-quadrant arc.
        const result = fi.find(circle(2, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a tangent contact that is not on the arc', () => {
        // Tangency at (-1,0), which is not on the first-quadrant arc.
        const result = fi.find(circle(-2, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('returns points on both the circle and the arc', () => {
        const rnd = makeRandom(60607);
        let numFound = 0;
        for (let k = 0; k < 400; ++k) {
            const a = arc(rnd() * 360, rnd() * 360);
            const c = circle(rnd() * 4 - 2, rnd() * 4 - 2, 0.2 + rnd() * 2);
            const result = fi.find(c, a);
            if (result.numIntersections === INT32_MAX) {
                continue;
            }
            expect(result.intersect).toBe(result.numIntersections > 0);
            for (let i = 0; i < result.numIntersections; ++i) {
                const p = result.point[i];
                expect(length(sub(p, c.center))).toBeCloseTo(c.radius, 8);
                expect(length(sub(p, a.center))).toBeCloseTo(a.radius, 8);
                expect(a.containsOnCircle(p)).toBe(true);
                ++numFound;
            }
        }
        expect(numFound).toBeGreaterThan(0);
    });

    it('handles a zero-radius circle', () => {
        // A degenerate circle at a point of the arc.
        const result = fi.find(circle(1, 0, 0), arc(0, 90));
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);

        // A degenerate circle away from the unit circle.
        expect(fi.find(circle(0.5, 0, 0), arc(0, 90)).intersect).toBe(false);
    });
});
