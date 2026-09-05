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
import { length } from '../src/Vector.js';
import { check, expectVectorClose, fc, positive, unitVector, wellScaled } from './helpers/arbitraries.js';

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

describe('IntrLine3Sphere3 verification', () => {
    const sphereArb = fc.tuple(
        fc.array(wellScaled(-5, 5), { minLength: 3, maxLength: 3 }),
        positive(4)).map(([c, r]) =>
            Hypersphere.fromCenterRadius(Vector.fromArray(c), r));
    const lineArb = fc.tuple(
        fc.array(wellScaled(-6, 6), { minLength: 3, maxLength: 3 }),
        unitVector(3))
        .map(([p, d]) => Line.fromOriginDirection(Vector.fromArray(p), d));

    // Distance of the line from tangency, relative to the radius: the two
    // roots are -a1 -+ sqrt(a1^2 - a0) and the square root loses half its
    // digits as the discriminant goes to zero.
    function transversality(l: Line, s: Hypersphere): number {
        const diff = sub(l.origin, s.center);
        const a1 = dot(l.direction, diff);
        const a0 = dot(diff, diff) - s.radius * s.radius;
        return Math.abs(a1 * a1 - a0) / (s.radius * s.radius);
    }

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(lineArb, sphereArb), ([l, s]) => {
            expect(fi.find(l, s).intersect).toBe(ti.test(l, s).intersect);
        });
    });

    it('reported points lie on the line and on the sphere', () => {
        check(fc.tuple(lineArb, sphereArb), ([l, s]) => {
            const r = fi.find(l, s);
            if (!r.intersect) {
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            // Upstream fills both points whenever 'intersect' is true, and
            // sets parameter[1] = parameter[0] in the tangent case.
            for (let k = 0; k < 2; ++k) {
                const p = r.point[k];
                for (let i = 0; i < 3; ++i) {
                    expect(Number.isFinite(p.values[i])).toBe(true);
                }
                expectVectorClose(p,
                    add(l.origin, mul(r.parameter[k], l.direction)), 0, 0);
            }
            if (transversality(l, s) < 1e-6) {
                return;    // near tangency: the root loses half its digits
            }
            const scale = 1 + length(l.origin) + length(s.center) + s.radius;
            for (let k = 0; k < r.numIntersections; ++k) {
                expect(Math.abs(length(sub(r.point[k], s.center)) - s.radius))
                    .toBeLessThanOrEqual(1e-9 * scale);
            }
        });
    });

    it('the roots satisfy the quadratic they solve', () => {
        check(fc.tuple(lineArb, sphereArb), ([l, s]) => {
            const r = fi.find(l, s);
            if (!r.intersect) {
                return;
            }
            const diff = sub(l.origin, s.center);
            const a1 = dot(l.direction, diff);
            const a0 = dot(diff, diff) - s.radius * s.radius;
            // Q(t) = t^2 + 2*a1*t + a0 has roots t0, t1 with
            // t0 + t1 = -2*a1 and t0*t1 = a0. The sum is formed by exact
            // cancellation of the two square roots, so it is tight; the
            // product inherits the conditioning of the roots.
            const t0 = r.parameter[0], t1 = r.parameter[1];
            expect(t0 + t1).toBeCloseTo(-2 * a1, 9);
            const scale = 1 + Math.abs(a0) + Math.abs(a1) * Math.abs(a1);
            expect(Math.abs(t0 * t1 - a0)).toBeLessThanOrEqual(1e-7 * scale);
        });
    });

    it('the reported interval matches a fine sweep of the line', () => {
        const rnd = makeRandom(0x63b1ea);
        check(fc.tuple(lineArb, sphereArb), ([l, s]) => {
            const r = fi.find(l, s);
            for (let k = 0; k < 400; ++k) {
                const t = 24 * rnd() - 12;
                const p = add(l.origin, mul(t, l.direction));
                if (length(sub(p, s.center)) < s.radius * (1 - 1e-9)) {
                    expect(r.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(r.parameter[0] - 1e-9);
                    expect(t).toBeLessThanOrEqual(r.parameter[1] + 1e-9);
                }
            }
        }, 60);
    }, 30000);

    it('a zero-radius sphere is met only by a line through its center', () => {
        const s = Hypersphere.fromCenterRadius(v3(1, 2, 3), 0);
        const through = makeLine(v3(0, 2, 3), v3(1, 0, 0));
        const r = fi.find(through, s);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(1);
        expectVectorClose(r.point[0], v3(1, 2, 3), 1e-12, 1e-12);
        expect(fi.find(makeLine(v3(0, 5, 3), v3(1, 0, 0)), s).intersect)
            .toBe(false);
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(fc.tuple(lineArb, sphereArb), ([l, s]) => {
            const r = defaultIntrLine3Sphere3FIResult();
            intrLine3Sphere3DoQuery(l.origin, l.direction, s, r);
            const expected = fi.find(l, s);
            expect(r.intersect).toBe(expected.intersect);
            expect(r.numIntersections).toBe(expected.numIntersections);
            expect(r.parameter[0]).toBe(expected.parameter[0]);
            expect(r.parameter[1]).toBe(expected.parameter[1]);
            // DoQuery leaves the points untouched; only find() fills them.
            if (!expected.intersect) {
                expect(r.point[0].values).toEqual([0, 0, 0]);
            }
        });
    });
});
