import { describe, it, expect } from 'vitest';
import {
    IntrLine3Rectangle3TI,
    IntrLine3Rectangle3FI
} from '../src/IntrLine3Rectangle3.js';
import { Line } from '../src/Line.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { length } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { check, expectVectorClose, fc, positive, rotationFrame, unitVector, wellScaled } from './helpers/arbitraries.js';

const ti = new IntrLine3Rectangle3TI();
const fi = new IntrLine3Rectangle3FI();

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

// The unit square in the z = 0 plane, extents (2,1) about the origin.
function xyRectangle(center: Vector, e0: number, e1: number): Rectangle {
    return Rectangle.fromCenterAxisExtent(center,
        [v3(1, 0, 0), v3(0, 1, 0)], Vector.fromArray([e0, e1]));
}

describe('IntrLine3Rectangle3TI', () => {
    it('reports a line through the rectangle center', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(0, 0, -5), v3(0, 0, 1));
        expect(ti.test(line, rect).intersect).toBe(true);
    });

    it('reports a line piercing near a corner as intersecting', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(1.9, 0.9, -1), v3(0, 0, 1));
        expect(ti.test(line, rect).intersect).toBe(true);
    });

    it('accepts a line piercing exactly at a corner (closed rectangle)', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(2, 1, -1), v3(0, 0, 1));
        expect(ti.test(line, rect).intersect).toBe(true);
    });

    it('rejects a line piercing the plane outside the rectangle', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(2.01, 0, -1), v3(0, 0, 1));
        expect(ti.test(line, rect).intersect).toBe(false);

        const line2 = makeLine(v3(0, 1.01, -1), v3(0, 0, 1));
        expect(ti.test(line2, rect).intersect).toBe(false);
    });

    it('reports no intersection for a line parallel to the plane, even when coplanar', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        // Coplanar line through the rectangle center.
        const coplanar = makeLine(v3(-5, 0, 0), v3(1, 0, 0));
        expect(ti.test(coplanar, rect).intersect).toBe(false);
        // Parallel but offset.
        const offset = makeLine(v3(-5, 0, 3), v3(1, 0, 0));
        expect(ti.test(offset, rect).intersect).toBe(false);
    });

    it('handles a degenerate rectangle with zero extents', () => {
        const rect = xyRectangle(v3(0, 0, 0), 0, 0);
        const hit = makeLine(v3(0, 0, -1), v3(0, 0, 1));
        expect(ti.test(hit, rect).intersect).toBe(true);
        const miss = makeLine(v3(0.1, 0, -1), v3(0, 0, 1));
        expect(ti.test(miss, rect).intersect).toBe(false);
    });
});

describe('IntrLine3Rectangle3FI', () => {
    it('computes the parameter, rectangle coordinates and point', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(0.5, -0.25, -3), v3(0, 0, 1));
        const result = fi.find(line, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[1]).toBeCloseTo(-0.25, 12);
        // The third component is unused by the rectangle query.
        expect(result.rectCoord[2]).toBe(0);
        expect(result.point.values[0]).toBeCloseTo(0.5, 12);
        expect(result.point.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('gives a negative parameter when the rectangle is behind the origin', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(0, 0, 4), v3(0, 0, 1));
        const result = fi.find(line, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(-4, 12);
    });

    it('leaves the default result when there is no intersection', () => {
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const line = makeLine(v3(5, 5, -1), v3(0, 0, 1));
        const result = fi.find(line, rect);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.rectCoord).toEqual([0, 0, 0]);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('reconstructs the point from the rectangle coordinates', () => {
        const axis0 = v3(1, 1, 0);
        const axis1 = v3(-1, 1, 1);
        normalize(axis0);
        normalize(axis1);
        const rect = Rectangle.fromCenterAxisExtent(v3(1, -2, 0.5),
            [axis0, axis1], Vector.fromArray([1.5, 0.75]));
        const line = makeLine(v3(3, 1, 4), v3(-0.4, -0.7, -1));
        const result = fi.find(line, rect);
        expect(result.intersect).toBe(true);

        const fromCoords = add(rect.center,
            add(mul(result.rectCoord[0], rect.axis[0]),
                mul(result.rectCoord[1], rect.axis[1])));
        for (let i = 0; i < 3; ++i) {
            expect(fromCoords.values[i]).toBeCloseTo(result.point.values[i], 10);
        }
        expect(Math.abs(result.rectCoord[0]))
            .toBeLessThanOrEqual(rect.extent.values[0] + 1e-12);
        expect(Math.abs(result.rectCoord[1]))
            .toBeLessThanOrEqual(rect.extent.values[1] + 1e-12);
    });
});

describe('IntrLine3Rectangle3 consistency', () => {
    it('agrees between TI and FI on random configurations', () => {
        const rand = makeRandom(987654321);
        const rect = xyRectangle(v3(0.25, -0.5, 0), 1.25, 0.75);
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const origin = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const direction = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(direction, direction) < 1e-8) {
                continue;
            }
            const line = makeLine(origin, direction);

            const tiResult = ti.test(line, rect);
            const fiResult = fi.find(line, rect);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            if (fiResult.intersect) {
                ++hits;
                // The reported point lies on the line and in the rectangle.
                const onLine = add(line.origin,
                    mul(fiResult.parameter, line.direction));
                for (let i = 0; i < 3; ++i) {
                    expect(onLine.values[i])
                        .toBeCloseTo(fiResult.point.values[i], 10);
                }
                const delta = sub(fiResult.point, rect.center);
                expect(Math.abs(dot(delta, rect.axis[0])))
                    .toBeLessThanOrEqual(rect.extent.values[0] + 1e-9);
                expect(Math.abs(dot(delta, rect.axis[1])))
                    .toBeLessThanOrEqual(rect.extent.values[1] + 1e-9);
                expect(Math.abs(dot(delta, v3(0, 0, 1)))).toBeLessThan(1e-9);
            }
        }
        expect(hits).toBeGreaterThan(10);
    });
});

describe('IntrLine3Rectangle3 verification', () => {
    const rectArb = fc.tuple(
        fc.array(wellScaled(-4, 4), { minLength: 3, maxLength: 3 }),
        rotationFrame(3), positive(3), positive(3))
        .map(([c, axes, e0, e1]) => Rectangle.fromCenterAxisExtent(
            Vector.fromArray(c), [axes[0], axes[1]],
            Vector.fromArray([e0, e1])));
    const lineArb = fc.tuple(
        fc.array(wellScaled(-6, 6), { minLength: 3, maxLength: 3 }),
        unitVector(3))
        .map(([p, d]) => Line.fromOriginDirection(Vector.fromArray(p), d));

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(lineArb, rectArb), ([l, rect]) => {
            expect(fi.find(l, rect).intersect).toBe(ti.test(l, rect).intersect);
        });
    });

    it('the reported point is on the line and inside the rectangle', () => {
        check(fc.tuple(lineArb, rectArb), ([l, rect]) => {
            const r = fi.find(l, rect);
            if (!r.intersect) {
                return;
            }
            expect(Number.isFinite(r.parameter)).toBe(true);
            expectVectorClose(r.point,
                add(l.origin, mul(r.parameter, l.direction)), 0, 0);

            // The rectangle coordinates must be inside the extents and must
            // reconstruct the point. Both come from one shared denominator
            // Dot(D,N), so the residual scales with the geometry, not with
            // the conditioning of the query.
            expect(Math.abs(r.rectCoord[0]))
                .toBeLessThanOrEqual(rect.extent.values[0] * (1 + 1e-9) + 1e-9);
            expect(Math.abs(r.rectCoord[1]))
                .toBeLessThanOrEqual(rect.extent.values[1] * (1 + 1e-9) + 1e-9);
            const rebuilt = add(add(rect.center,
                mul(r.rectCoord[0], rect.axis[0])),
                mul(r.rectCoord[1], rect.axis[1]));
            const scale = 1 + length(l.origin) + length(rect.center)
                + Math.abs(r.parameter);
            expect(length(sub(rebuilt, r.point)))
                .toBeLessThanOrEqual(1e-9 * scale);

            // Upstream declares rectCoord as a 3-array copied from the
            // triangle query; the third entry is always zero (issue #141).
            expect(r.rectCoord[2]).toBe(0);
        });
    });

    it('a line through a sampled rectangle point always intersects', () => {
        const rnd = makeRandom(0x2f81cd);
        check(fc.tuple(rectArb, unitVector(3)), ([rect, d]) => {
            const s0 = (2 * rnd() - 1) * rect.extent.values[0] * 0.9;
            const s1 = (2 * rnd() - 1) * rect.extent.values[1] * 0.9;
            const target = add(add(rect.center, mul(s0, rect.axis[0])),
                mul(s1, rect.axis[1]));
            const normal = cross(rect.axis[0], rect.axis[1]);
            if (Math.abs(dot(d, normal)) < 1e-3) {
                return;    // nearly in the plane: upstream reports no hit
            }
            const l = Line.fromOriginDirection(
                sub(target, mul(3, d)), d);
            const r = fi.find(l, rect);
            expect(r.intersect).toBe(true);
            expect(ti.test(l, rect).intersect).toBe(true);
            const scale = 1 + length(target);
            expect(length(sub(r.point, target)))
                .toBeLessThanOrEqual(1e-9 * scale);
            expect(r.rectCoord[0]).toBeCloseTo(s0, 8);
            expect(r.rectCoord[1]).toBeCloseTo(s1, 8);
        });
    });

    it('a line in the plane of the rectangle reports no intersection', () => {
        // Upstream calls the parallel case "no intersection" even when the
        // line lies in the plane and crosses the rectangle.
        const rect = xyRectangle(v3(0, 0, 0), 2, 1);
        const inPlane = makeLine(v3(-5, 0, 0), v3(1, 0, 0));
        expect(ti.test(inPlane, rect).intersect).toBe(false);
        const r = fi.find(inPlane, rect);
        expect(r.intersect).toBe(false);
        expect(r.parameter).toBe(0);
        expect(r.point.values).toEqual([0, 0, 0]);
        // A line parallel to but off the plane is also reported as a miss.
        expect(ti.test(makeLine(v3(-5, 0, 1), v3(1, 0, 0)), rect).intersect)
            .toBe(false);
    });

    it('parameters scale with a non-unit direction while the point does not', () => {
        check(fc.tuple(lineArb, rectArb, fc.integer({ min: 2, max: 4 })),
            ([l, rect, k]) => {
                const a = fi.find(l, rect);
                if (!a.intersect) {
                    return;
                }
                const b = fi.find(Line.fromOriginDirection(l.origin,
                    mul(k, l.direction)), rect);
                expect(b.intersect).toBe(true);
                expect(b.parameter * k).toBeCloseTo(a.parameter, 9);
                expectVectorClose(b.point, a.point, 1e-9, 1e-9);
                expect(b.rectCoord[0]).toBeCloseTo(a.rectCoord[0], 9);
                expect(b.rectCoord[1]).toBeCloseTo(a.rectCoord[1], 9);
            });
    });

    it('a degenerate (zero-extent) rectangle is the point test', () => {
        const rect = xyRectangle(v3(1, 2, 0), 0, 0);
        expect(fi.find(makeLine(v3(1, 2, -3), v3(0, 0, 1)), rect).intersect)
            .toBe(true);
        expect(fi.find(makeLine(v3(1.5, 2, -3), v3(0, 0, 1)), rect).intersect)
            .toBe(false);
    });
});
