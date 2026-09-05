import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { IntrLine3Rectangle3FI } from '../src/IntrLine3Rectangle3.js';
import {
    IntrRay3Rectangle3TI,
    IntrRay3Rectangle3FI
} from '../src/IntrRay3Rectangle3.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The rectangle in the plane z = 0, centered at the origin, with extents
// 2 (along x) and 1 (along y).
const rect = Rectangle.fromCenterAxisExtent(
    vec([0, 0, 0]), [vec([1, 0, 0]), vec([0, 1, 0])], vec([2, 1]));

describe('IntrRay3Rectangle3', () => {
    const ti = new IntrRay3Rectangle3TI();
    const fi = new IntrRay3Rectangle3FI();

    it('finds the crossing of a ray through the rectangle interior', () => {
        const r = ray([0.5, -0.25, 3], [0, 0, -1]);
        expect(ti.test(r, rect).intersect).toBe(true);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.5, 12);
        expect(result.point.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports no intersection when the rectangle is behind the ray', () => {
        const r = ray([0.5, -0.25, 3], [0, 0, 1]);
        expect(ti.test(r, rect).intersect).toBe(false);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.rectCoord).toEqual([0, 0, 0]);
    });

    it('accepts a ray whose origin is on the rectangle', () => {
        const r = ray([1, 0.5, 0], [0, 0, 1]);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.rectCoord[0]).toBeCloseTo(1, 12);
        expect(result.rectCoord[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no intersection just outside an edge', () => {
        const r = ray([2.0001, 0, 3], [0, 0, -1]);
        expect(ti.test(r, rect).intersect).toBe(false);
        expect(fi.find(r, rect).intersect).toBe(false);
    });

    it('accepts a ray hitting a corner exactly', () => {
        const r = ray([2, 1, 3], [0, 0, -1]);
        const result = fi.find(r, rect);
        expect(result.intersect).toBe(true);
        expect(result.rectCoord[0]).toBeCloseTo(2, 12);
        expect(result.rectCoord[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the ray is parallel to the plane', () => {
        // Upstream calls a coplanar ray a "no intersection" case.
        const r = ray([-5, 0, 0], [1, 0, 0]);
        expect(ti.test(r, rect).intersect).toBe(false);
        expect(fi.find(r, rect).intersect).toBe(false);
    });

    it('is the line query restricted to t >= 0', () => {
        const rand = makeRandom(5150);
        const lineFI = new IntrLine3Rectangle3FI();
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray(
                [6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, rect);
            const rayResult = fi.find(r, rect);
            expect(ti.test(r, rect).intersect).toBe(rayResult.intersect);

            if (lineResult.intersect && lineResult.parameter >= 0) {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter).toBeCloseTo(lineResult.parameter, 12);
                expect(rayResult.rectCoord[0]).toBeCloseTo(lineResult.rectCoord[0], 12);
                expect(rayResult.rectCoord[1]).toBeCloseTo(lineResult.rectCoord[1], 12);
            }
            else {
                expect(rayResult.intersect).toBe(false);
            }
        }
    });

    it('produces points that lie on the ray and inside the rectangle', () => {
        const rand = makeRandom(60606);
        // A tilted rectangle to exercise the general case.
        const a0 = vec([1, 1, 0]);
        normalize(a0);
        const a1 = vec([-1, 1, 1]);
        normalize(a1);
        // Make a1 orthogonal to a0 (it already is, up to rounding).
        const tilted = Rectangle.fromCenterAxisExtent(
            vec([0.5, -1, 2]), [a0, a1], vec([1.5, 0.75]));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray(
                [8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, tilted);
            if (!result.intersect) {
                continue;
            }
            ++hits;
            expect(result.parameter).toBeGreaterThanOrEqual(0);
            const onRay = add(r.origin, mul(result.parameter, r.direction));
            for (let k = 0; k < 3; ++k) {
                expect(result.point.values[k]).toBeCloseTo(onRay.values[k], 9);
            }
            const diff = sub(result.point, tilted.center);
            expect(Math.abs(dot(diff, tilted.axis[0])))
                .toBeLessThanOrEqual(tilted.extent.values[0] + 1e-9);
            expect(Math.abs(dot(diff, tilted.axis[1])))
                .toBeLessThanOrEqual(tilted.extent.values[1] + 1e-9);
            expect(dot(diff, tilted.axis[0]))
                .toBeCloseTo(result.rectCoord[0], 9);
            expect(dot(diff, tilted.axis[1]))
                .toBeCloseTo(result.rectCoord[1], 9);
        }
        expect(hits).toBeGreaterThan(5);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, rotationFrame, unitVector, wellScaledVector, expectClose,
    expectVectorClose
} from './helpers/arbitraries.js';
import { length as vlength } from '../src/Vector.js';

const rayRect3 = fc.tuple(wellScaledVector(3, -6, 6), unitVector(3),
    wellScaledVector(3, -4, 4), rotationFrame(3),
    fc.double({ min: 0.2, max: 4, noNaN: true }),
    fc.double({ min: 0.2, max: 4, noNaN: true }))
    .map(([o, d, c, R, e0, e1]) => ({
        ray: Ray.fromOriginDirection(o, d),
        rectangle: Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
            Vector.fromArray([e0, e1]))
    }));

describe('IntrRay3Rectangle3 verification', () => {
    const tiq = new IntrRay3Rectangle3TI();
    const fiq = new IntrRay3Rectangle3FI();
    const lfi = new IntrLine3Rectangle3FI();

    it('TI and FI agree on intersect', () => {
        check(rayRect3, ({ ray: r, rectangle: q }) => {
            expect(tiq.test(r, q).intersect).toBe(fiq.find(r, q).intersect);
        });
    });

    it('a ray hit is the line hit with a nonnegative parameter', () => {
        check(rayRect3, ({ ray: r, rectangle: q }) => {
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lf = lfi.find(l, q);
            const f = fiq.find(r, q);
            const hit = lf.intersect && lf.parameter >= 0;
            expect(f.intersect).toBe(hit);
            if (!hit) {
                // The default (unset) fields are kept on a miss.
                expect(f.parameter).toBe(0);
                expect(f.rectCoord).toEqual([0, 0, 0]);
                expect(f.point.values).toEqual([0, 0, 0]);
                return;
            }
            expect(f.parameter).toBe(lf.parameter);
            expectVectorClose(f.point, lf.point, 0, 0);
            expect(f.rectCoord[0]).toBe(lf.rectCoord[0]);
            expect(f.rectCoord[1]).toBe(lf.rectCoord[1]);
            // Upstream #141(2), preserved: rectCoord is a 3-array whose last
            // entry is always zero.
            expect(f.rectCoord[2]).toBe(0);
        });
    });

    it('the hit point lies on the ray and inside the rectangle', () => {
        check(rayRect3, ({ ray: r, rectangle: q }) => {
            const f = fiq.find(r, q);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter).toBeGreaterThanOrEqual(0);
            expectVectorClose(f.point,
                add(r.origin, mul(f.parameter, r.direction)), 1e-12, 1e-12);
            expect(Math.abs(f.rectCoord[0]))
                .toBeLessThanOrEqual(q.extent.values[0] + 1e-9);
            expect(Math.abs(f.rectCoord[1]))
                .toBeLessThanOrEqual(q.extent.values[1] + 1e-9);
            // The rectangle coordinates reconstruct the point.
            expectVectorClose(f.point,
                add(q.center, add(mul(f.rectCoord[0], q.axis[0]),
                    mul(f.rectCoord[1], q.axis[1]))), 1e-8, 1e-9);
            for (const x of f.point.values) {
                expect(Number.isNaN(x)).toBe(false);
            }
        });
    });

    it('a ray fired at a sampled rectangle point always hits', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }),
            fc.double({ min: -0.9, max: 0.9, noNaN: true }), unitVector(3),
            fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([c, R, e0, e1, s0, s1, d, back]) => {
                const q = Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
                    Vector.fromArray([e0, e1]));
                const target = add(c, add(mul(s0 * e0, R[0]),
                    mul(s1 * e1, R[1])));
                const r = Ray.fromOriginDirection(sub(target, mul(back, d)), d);
                // Skip rays nearly parallel to the rectangle plane; upstream
                // reports no intersection for a coplanar ray.
                const normal = R[2];
                if (Math.abs(dot(normal, d)) < 0.05) {
                    return;
                }
                const f = fiq.find(r, q);
                expect(f.intersect).toBe(true);
                expect(tiq.test(r, q).intersect).toBe(true);
                expectVectorClose(f.point, target, 1e-8, 1e-9);
                expectClose(f.rectCoord[0], s0 * e0, 1e-8, 1e-9);
                expectClose(f.rectCoord[1], s1 * e1, 1e-8, 1e-9);
                // The reversed ray points away from the rectangle.
                const rev = Ray.fromOriginDirection(r.origin, mul(-1, d));
                expect(fiq.find(rev, q).intersect).toBe(false);
            });
    });

    it('a ray in the plane of the rectangle reports no intersection', () => {
        // Upstream returns "no intersection" only when Dot(D,N) is exactly
        // zero, so the configuration is built with an axis-aligned rectangle
        // in the z = 0 plane and a direction whose z component is exactly 0.
        // (A direction merely *near* the plane does intersect: the line still
        // crosses the plane, just at a large parameter.)
        check(fc.tuple(wellScaledVector(3, -4, 4),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: 0.2, max: 4, noNaN: true }),
            fc.double({ min: -1, max: 1, noNaN: true }),
            fc.double({ min: -1, max: 1, noNaN: true }),
            fc.double({ min: -5, max: 5, noNaN: true })),
            ([c, e0, e1, a, b, shift]) => {
                if (Math.hypot(a, b) < 0.1) {
                    return;
                }
                const q = Rectangle.fromCenterAxisExtent(c,
                    [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0])],
                    Vector.fromArray([e0, e1]));
                const dir = Vector.fromArray([a, b, 0]);
                normalize(dir);
                expect(dir.values[2]).toBe(0);
                const origin = add(c, mul(shift, dir));
                const r = Ray.fromOriginDirection(origin, dir);
                expect(fiq.find(r, q).intersect).toBe(false);
                expect(tiq.test(r, q).intersect).toBe(false);
            });
    });

    it('a zero-extent rectangle is a single point', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
            unitVector(3), fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([c, R, d, back]) => {
                const q = Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
                    Vector.fromArray([0, 0]));
                const r = Ray.fromOriginDirection(sub(c, mul(back, d)), d);
                if (Math.abs(dot(R[2], d)) < 0.05) {
                    return;
                }
                const f = fiq.find(r, q);
                if (f.intersect) {
                    expect(vlength(sub(f.point, c))).toBeLessThan(1e-6);
                }
            });
    });
});
