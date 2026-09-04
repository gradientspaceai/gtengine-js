import { describe, expect, it } from 'vitest';
import { DistLine2Circle2 } from '../src/DistLine2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistPointLine } from '../src/DistPointLine.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of DistLine2Circle2.ts
// against the upstream header DistLine2Circle2.h.
// ---------------------------------------------------------------------------

function rot2(R: readonly Vector[], p: Vector): Vector {
    return add(mul(p.values[0], R[0]), mul(p.values[1], R[1]));
}

// A 2D line with a well-scaled, deliberately non-unit direction.
const nonUnitLine2 = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

const circle2 = fc.tuple(wellScaledVector(2, -8, 8), positive(6))
    .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

describe('DistLine2Circle2 verification', () => {
    const query = new DistLine2Circle2();
    const pointLine = new DistPointLine();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(nonUnitLine2, circle2), ([line, circle]) => {
                const r = query.compute(line, circle);
                expect(r.numClosestPairs === 1 || r.numClosestPairs === 2)
                    .toBe(true);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                const diff = sub(r.closest[0][0], r.closest[0][1]);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
                for (let j = 0; j < r.numClosestPairs; ++j) {
                    // closest[j][0] is on the line at parameter[j].
                    expectVectorClose(r.closest[j][0],
                        add(line.origin, mul(r.parameter[j], line.direction)),
                        1e-9, 1e-9);
                    // closest[j][1] is on the circle.
                    expectClose(length(sub(r.closest[j][1], circle.center)),
                        circle.radius, 1e-7, 1e-7);
                }
            });
        });

    it('matches | distance(center, line) - radius | when they do not meet',
        () => {
            check(fc.tuple(nonUnitLine2, circle2), ([line, circle]) => {
                const r = query.compute(line, circle);
                const h = pointLine.compute(circle.center, line).distance;
                if (r.numClosestPairs === 1) {
                    expectClose(r.distance, Math.abs(h - circle.radius), 1e-7,
                        1e-7);
                    // The single closest line point is the projection of the
                    // circle center onto the line.
                    expectClose(r.parameter[0],
                        pointLine.compute(circle.center, line).parameter, 1e-9,
                        1e-9);
                }
                else {
                    expect(r.distance).toBe(0);
                    expect(h).toBeLessThanOrEqual(
                        circle.radius * (1 + 1e-9) + 1e-9);
                }
            });
        });

    it('returns two intersection points in increasing parameter order', () => {
        check(fc.tuple(nonUnitLine2, circle2)
            .filter(([line, circle]) =>
                pointLine.compute(circle.center, line).distance
                < 0.9 * circle.radius),
        ([line, circle]) => {
            const r = query.compute(line, circle);
            expect(r.numClosestPairs).toBe(2);
            expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            for (let j = 0; j < 2; ++j) {
                // Each point is on the line and on the circle.
                expectVectorClose(r.closest[j][0], r.closest[j][1], 0, 0);
                expectClose(length(sub(r.closest[j][0], circle.center)),
                    circle.radius, 1e-7, 1e-7);
            }
        });
    });

    it('is minimal over sampled line/circle point pairs', () => {
        const rand = seededRandom(0x51e1);
        check(fc.tuple(nonUnitLine2, circle2), ([line, circle]) => {
            const r = query.compute(line, circle);
            for (let k = 0; k < 24; ++k) {
                const t = 40 * (rand() - 0.5);
                const a = 2 * Math.PI * rand();
                const p = add(line.origin, mul(t, line.direction));
                const q = add(circle.center, Vector.fromArray(
                    [circle.radius * Math.cos(a),
                        circle.radius * Math.sin(a)]));
                const gap = length(sub(p, q));
                expect(r.distance).toBeLessThanOrEqual(gap + 1e-9 * (1 + gap));
            }
        }, 60);
    });

    it('is equivariant under rigid motions of the plane', () => {
        check(fc.tuple(nonUnitLine2, circle2, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([line, circle, R, tr]) => {
            const movedLine = Line.fromOriginDirection(
                add(rot2(R, line.origin), tr), rot2(R, line.direction));
            const movedCircle = Hypersphere.fromCenterRadius(
                add(rot2(R, circle.center), tr), circle.radius);
            const r0 = query.compute(line, circle);
            const r1 = query.compute(movedLine, movedCircle);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expect(r1.numClosestPairs).toBe(r0.numClosestPairs);
        });
    });

    it('reports a single tangency point for a tangent line', () => {
        // The tangency test is exact only when DotPerp(D,delta)^2 and
        // r^2*Dot(D,D) round identically, so the tangent line is built
        // axis-aligned: delta = (+-r, 0) or (0, +-r) with D the other axis
        // gives DotPerp(D,delta)^2 = r^2 and Dot(D,D) = 1 exactly. A rotated
        // frame would leave Dot(D,D) one ulp off 1 and report two nearly
        // coincident intersections instead.
        // The center and radius are integers so that C + r*N and the
        // subtraction (C + r*N) - C are both exact; with generic doubles the
        // cancellation leaves delta one ulp off r*N.
        check(fc.tuple(fc.integer({ min: -8, max: 8 }),
            fc.integer({ min: -8, max: 8 }), fc.integer({ min: 1, max: 6 }),
            fc.integer({ min: 0, max: 3 })),
        ([cx, cy, radius, which]) => {
            const circle = Hypersphere.fromCenterRadius(v(cx, cy), radius);
            {
                const sign = which < 2 ? 1 : -1;
                const N = which % 2 === 0
                    ? Vector.fromArray([sign, 0])
                    : Vector.fromArray([0, sign]);
                const D = which % 2 === 0
                    ? Vector.fromArray([0, 1])
                    : Vector.fromArray([1, 0]);
                const origin = add(circle.center, mul(circle.radius, N));
                const r = query.compute(Line.fromOriginDirection(origin, D),
                    circle);
                expect(r.numClosestPairs).toBe(1);
                expect(r.distance).toBe(0);
                // At tangency upstream sets closest[0][1] = closest[0][0].
                expectVectorClose(r.closest[0][0], r.closest[0][1], 0, 0);
                expectClose(length(sub(r.closest[0][0], circle.center)),
                    circle.radius, 1e-9, 1e-9);
            }
        });
    });
});
