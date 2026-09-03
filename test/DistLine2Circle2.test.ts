import { describe, expect, it } from 'vitest';
import { DistLine2Circle2 } from '../src/DistLine2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v(...center), radius);
}

describe('DistLine2Circle2', () => {
    const query = new DistLine2Circle2();

    it('measures a line that misses the circle', () => {
        // The horizontal line y = 5 and the unit circle at the origin.
        const result = query.compute(line([0, 5], [1, 0]), circle([0, 0], 1));
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.closest[0][0].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[0][0].values[1]).toBeCloseTo(5, 10);
        expect(result.closest[0][1].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[0][1].values[1]).toBeCloseTo(1, 10);
    });

    it('reports one coincident pair for a tangent line', () => {
        const result = query.compute(line([0, 2], [1, 0]), circle([0, 0], 2));
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.closest[0][0].values[1]).toBeCloseTo(2, 10);
        expect(result.closest[0][1].values[1]).toBeCloseTo(2, 10);
    });

    it('reports two intersection points for a secant line', () => {
        // y = 0 meets the unit circle at (-1,0) and (1,0).
        const result = query.compute(line([0, 0], [1, 0]), circle([0, 0], 1));
        expect(result.numClosestPairs).toBe(2);
        expect(result.distance).toBe(0);
        expect(result.parameter[0]).toBeCloseTo(-1, 10);
        expect(result.parameter[1]).toBeCloseTo(1, 10);
        expect(result.closest[0][0].values[0]).toBeCloseTo(-1, 10);
        expect(result.closest[1][0].values[0]).toBeCloseTo(1, 10);
    });

    it('orders the two parameters increasingly', () => {
        const result = query.compute(line([3, 0], [-1, 0]),
            circle([0, 0], 1));
        expect(result.numClosestPairs).toBe(2);
        expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
        expect(result.parameter[0]).toBeCloseTo(2, 10);
        expect(result.parameter[1]).toBeCloseTo(4, 10);
    });

    it('translates correctly for an off-origin circle', () => {
        const c = circle([4, -3], 2);
        const result = query.compute(line([4, 7], [3, 0]), c);
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(8, 10);
        expect(result.closest[0][1].values[0]).toBeCloseTo(4, 10);
        expect(result.closest[0][1].values[1]).toBeCloseTo(-1, 10);
    });

    it('handles a non-unit direction', () => {
        const result = query.compute(line([0, 3], [5, 0]), circle([0, 0], 1));
        expect(result.numClosestPairs).toBe(1);
        expect(result.distance).toBeCloseTo(2, 10);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
    });

    it('reports valid closest points and matches a sampled minimum', () => {
        let seed = 24680;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 80; ++trial) {
            const c = circle([rand(), rand()], Math.abs(rand()) + 0.25);
            const l = line([rand(), rand()], [rand() + 5, rand()]);
            const result = query.compute(l, c);

            for (let j = 0; j < result.numClosestPairs; ++j) {
                // The line point matches the reported parameter.
                const onLine = add(l.origin,
                    mul(result.parameter[j], l.direction));
                expect(result.closest[j][0].values[0]).toBeCloseTo(
                    onLine.values[0], 8);
                expect(result.closest[j][0].values[1]).toBeCloseTo(
                    onLine.values[1], 8);

                // The circle point is on the circle.
                const radial = sub(result.closest[j][1], c.center);
                expect(Math.sqrt(dot(radial, radial))).toBeCloseTo(c.radius,
                    8);
            }

            // Compare against a sampled minimum over circle angles.
            let best = Number.MAX_VALUE;
            const a00 = dot(l.direction, l.direction);
            for (let k = 0; k < 3600; ++k) {
                const s = k * Math.PI / 1800;
                const q = v(c.center.values[0] + c.radius * Math.cos(s),
                    c.center.values[1] + c.radius * Math.sin(s));
                const w = sub(q, l.origin);
                const t = dot(l.direction, w) / a00;
                const d = sub(w, mul(t, l.direction));
                best = Math.min(best, dot(d, d));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
            expect(result.distance).toBeCloseTo(
                Math.sqrt(result.sqrDistance), 10);
        }
    });
});
