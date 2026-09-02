import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { IntrLine2Arc2TI, IntrLine2Arc2FI } from '../src/IntrLine2Arc2';
import { Line } from '../src/Line';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

// The quarter arc of the unit circle from (1,0) to (0,1), which contains the
// points (cos t, sin t) for t in [0, pi/2].
function quarterArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

describe('IntrLine2Arc2', () => {
    const ti = new IntrLine2Arc2TI();
    const fi = new IntrLine2Arc2FI();

    it('finds the single intersection of a chord that crosses the arc once', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 0.5], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);

        const expectedX = Math.sqrt(0.75);
        expect(result.point[0].values[0]).toBeCloseTo(expectedX, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0.5, 12);
        expect(result.parameter[0]).toBeCloseTo(expectedX, 12);
        expect(ti.test(line([0, 0.5], [1, 0]), arc).intersect).toBe(true);
    });

    it('finds two intersections when both circle points are on the arc', () => {
        // The arc from (1,0) to (-1,0) going through (0,1) is the upper half
        // circle; the horizontal line y = 0.5 meets it twice.
        const arc = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0),
            vec(-1, 0));
        const result = fi.find(line([0, 0.5], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(result.point[i].values[1]).toBeCloseTo(0.5, 12);
            expect(Math.abs(result.point[i].values[0]))
                .toBeCloseTo(Math.sqrt(0.75), 12);
        }
        // Parameters are ordered as the underlying line-circle query orders
        // them, smallest first.
        expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
    });

    it('accepts an arc endpoint (the on-boundary test is inclusive)', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 0], [0, 1]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the circle points miss the arc', () => {
        const arc = quarterArc();
        // The line y = -0.5 meets the circle only in the lower half plane.
        const result = fi.find(line([0, -0.5], [1, 0]), arc);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(line([0, -0.5], [1, 0]), arc).intersect).toBe(false);
    });

    it('reports no intersection when the line misses the circle entirely', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 3], [1, 0]), arc);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(line([0, 3], [1, 0]), arc).intersect).toBe(false);
    });

    it('reports the tangent point when it lies on the arc', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 1], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 9);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('agrees with a direct arc-membership check on random lines', () => {
        let seed = 24680;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const arc = Arc2.fromCenterRadiusEnds(vec(0.5, -0.25), 2,
            add(vec(0.5, -0.25), mul(2, vec(Math.cos(0.3), Math.sin(0.3)))),
            add(vec(0.5, -0.25), mul(2, vec(Math.cos(2.1), Math.sin(2.1)))));

        for (let trial = 0; trial < 200; ++trial) {
            const l = line([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            if (dot(l.direction, l.direction) < 0.5) {
                continue;
            }
            const result = fi.find(l, arc);
            expect(ti.test(l, arc).intersect).toBe(result.intersect);

            for (let i = 0; i < result.numIntersections; ++i) {
                const p = result.point[i];
                // The point is on the line at the reported parameter.
                const q = add(l.origin, mul(result.parameter[i], l.direction));
                expect(sub(p, q).values[0]).toBeCloseTo(0, 9);
                expect(sub(p, q).values[1]).toBeCloseTo(0, 9);

                // The point is on the circle of the arc.
                const r = sub(p, arc.center);
                expect(Math.sqrt(dot(r, r))).toBeCloseTo(arc.radius, 9);

                // The point is on the arc.
                expect(arc.containsOnCircle(p)).toBe(true);
            }
        }
    });
});
