import { describe, it, expect } from 'vitest';
import { Rectangle } from '../src/Rectangle.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub, normalize } from '../src/Vector.js';
import {
    IntrSegment3Rectangle3TI,
    IntrSegment3Rectangle3FI
} from '../src/IntrSegment3Rectangle3.js';
import { Line } from '../src/Line.js';
import { IntrLine3Rectangle3FI } from '../src/IntrLine3Rectangle3.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The unit square in the z = 0 plane, extents 1 along x and y.
const square = Rectangle.fromCenterAxisExtent(vec(0, 0, 0),
    [vec(1, 0, 0), vec(0, 1, 0)], Vector.fromArray([1, 1]));

describe('IntrSegment3Rectangle3', () => {
    const ti = new IntrSegment3Rectangle3TI();
    const fi = new IntrSegment3Rectangle3FI();

    it('finds a crossing at the rectangle center', () => {
        const s = segment([0, 0, -2], [0, 0, 2]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0, 12);
        expect(result.rectCoord[1]).toBeCloseTo(0, 12);
        // The third component is unused by IntrLine3Rectangle3 (issue #141).
        expect(result.rectCoord[2]).toBe(0);
        expect(result.point.values[0]).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('finds an off-center crossing and reports rectangle coordinates', () => {
        const s = segment([0.5, -0.25, -1], [0.5, -0.25, 3]);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.25, 12);
        expect(result.rectCoord[0]).toBeCloseTo(0.5, 12);
        expect(result.rectCoord[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.5, 12);
        expect(result.point.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('accepts an intersection exactly at a rectangle corner', () => {
        const s = segment([1, 1, -1], [1, 1, 1]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.rectCoord[0]).toBeCloseTo(1, 12);
        expect(result.rectCoord[1]).toBeCloseTo(1, 12);
    });

    it('accepts an intersection exactly at a segment endpoint', () => {
        // The endpoint (0,0,0) is on the rectangle, so t = 0.
        const s = segment([0, 0, 0], [0, 0, 4]);
        expect(ti.test(s, square).intersect).toBe(true);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0, 12);
    });

    it('rejects an intersection beyond the segment endpoint', () => {
        // The supporting line hits the rectangle at t = 2.
        const s = segment([0, 0, -2], [0, 0, -1]);
        expect(ti.test(s, square).intersect).toBe(false);
        const result = fi.find(s, square);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('rejects a crossing outside the rectangle extents', () => {
        const s = segment([3, 0, -1], [3, 0, 1]);
        expect(ti.test(s, square).intersect).toBe(false);
        expect(fi.find(s, square).intersect).toBe(false);
    });

    it('reports no intersection when the segment is in the plane', () => {
        // The upstream queries deliberately report no intersection for a
        // segment lying in the plane of the rectangle.
        const s = segment([-2, 0, 0], [2, 0, 0]);
        expect(ti.test(s, square).intersect).toBe(false);
        expect(fi.find(s, square).intersect).toBe(false);
    });

    it('handles a degenerate (zero-length) segment', () => {
        // The direction is (0,0,0), so DdN = 0 and the line query reports no
        // intersection even for a point on the rectangle.
        const onIt = segment([0, 0, 0], [0, 0, 0]);
        expect(ti.test(onIt, square).intersect).toBe(false);
        expect(fi.find(onIt, square).intersect).toBe(false);
    });

    it('handles a rotated, translated rectangle', () => {
        const a0 = vec(1, 1, 0);
        normalize(a0);
        const a1 = vec(-1, 1, 1);
        normalize(a1);
        const rect = Rectangle.fromCenterAxisExtent(vec(2, -1, 3), [a0, a1],
            Vector.fromArray([1.5, 0.75]));
        // A point on the rectangle, then a segment straddling it.
        const target = add(rect.center,
            add(mul(0.6, a0), mul(-0.3, a1)));
        const n = vec(
            a0.values[1] * a1.values[2] - a0.values[2] * a1.values[1],
            a0.values[2] * a1.values[0] - a0.values[0] * a1.values[2],
            a0.values[0] * a1.values[1] - a0.values[1] * a1.values[0]);
        const s = Segment.fromEndpoints(sub(target, mul(2, n)),
            add(target, mul(2, n)));
        expect(ti.test(s, rect).intersect).toBe(true);
        const result = fi.find(s, rect);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0.5, 10);
        expect(result.rectCoord[0]).toBeCloseTo(0.6, 10);
        expect(result.rectCoord[1]).toBeCloseTo(-0.3, 10);
        const diff = sub(result.point, target);
        expect(Math.sqrt(dot(diff, diff))).toBeLessThan(1e-10);
    });

    it('agrees with an independent plane solve on random configurations', () => {
        const rnd = makeRandom(60221408);
        let tiFiMismatch = 0;
        let referenceMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;

        for (let trial = 0; trial < 400; ++trial) {
            const p0 = vec(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const p1 = vec(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-4 || Math.abs(d.values[2]) < 1e-3) {
                continue;
            }
            const s = Segment.fromEndpoints(p0, p1);

            const tiResult = ti.test(s, square);
            const fiResult = fi.find(s, square);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            // Reference: solve p0.z + t * d.z = 0 and check the bounds.
            const t = -p0.values[2] / d.values[2];
            const q = add(p0, mul(t, d));
            const reference = (t >= 0 && t <= 1 &&
                Math.abs(q.values[0]) <= 1 && Math.abs(q.values[1]) <= 1);
            if (reference !== fiResult.intersect) {
                ++referenceMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                if (Math.abs(fiResult.parameter - t) > 1e-9) {
                    ++pointMismatch;
                }
                const diff = sub(fiResult.point, q);
                if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                    ++pointMismatch;
                }
                if (Math.abs(fiResult.rectCoord[0] - q.values[0]) > 1e-9 ||
                    Math.abs(fiResult.rectCoord[1] - q.values[1]) > 1e-9) {
                    ++pointMismatch;
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, referenceMismatch, pointMismatch])
            .toEqual([0, 0, 0]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrSegment3Rectangle3.h.
// ---------------------------------------------------------------------------

describe('IntrSegment3Rectangle3 verification', () => {
    const ti = new IntrSegment3Rectangle3TI();
    const fi = new IntrSegment3Rectangle3FI();
    const lineFi = new IntrLine3Rectangle3FI();

    const arbRectangle = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
        fc.array(positive(3, 0.25), { minLength: 2, maxLength: 2 }))
        .map(([c, frame, e]) => Rectangle.fromCenterAxisExtent(c,
            [frame[0], frame[1]], Vector.fromArray(e)));
    const arbSegment = fc.tuple(wellScaledVector(3, -8, 8),
        wellScaledVector(3, -8, 8))
        .filter(([a, b]) => {
            const d = sub(b, a);
            return dot(d, d) > 1e-2;
        })
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbPair = fc.tuple(arbSegment, arbRectangle);

    it('TI and FI agree on intersect', () => {
        check(arbPair, ([s, r]) => {
            expect(ti.test(s, r).intersect).toBe(fi.find(s, r).intersect);
        });
    });

    it('the FI point lies on the segment and on the rectangle', () => {
        check(arbPair, ([s, r]) => {
            const res = fi.find(s, r);
            if (!res.intersect) {
                expect(res.parameter).toBe(0);
                expect(res.rectCoord).toEqual([0, 0, 0]);
                expect(res.point.equals(Vector.zero(3))).toBe(true);
                return;
            }
            const t = res.parameter;
            expect(t).toBeGreaterThanOrEqual(0);
            expect(t).toBeLessThanOrEqual(1);
            const d = sub(s.p[1], s.p[0]);
            expectVectorClose(res.point, add(s.p[0], mul(t, d)), 1e-9, 1e-9);
            // The point is inside the rectangle, and its rectangle
            // coordinates are the reported ones.
            const delta = sub(res.point, r.center);
            for (let i = 0; i < 2; ++i) {
                const c = dot(r.axis[i], delta);
                expectClose(c, res.rectCoord[i], 1e-8, 1e-8);
                expect(Math.abs(c)).toBeLessThanOrEqual(
                    r.extent.get(i) * (1 + 1e-8) + 1e-8);
            }
            // The upstream rectCoord has a vestigial third component that the
            // line-rectangle query never assigns (upstream issue #141).
            expect(res.rectCoord[2]).toBe(0);
            // The point is in the plane of the rectangle.
            const n = cross(r.axis[0], r.axis[1]);
            expectClose(dot(n, delta), 0, 1e-8, 1e-8);
        });
    });

    it('the segment result is the line result restricted to t in [0,1]', () => {
        check(arbPair, ([s, r]) => {
            const d = sub(s.p[1], s.p[0]);
            const lineRes = lineFi.find(Line.fromOriginDirection(s.p[0], d), r);
            const segRes = fi.find(s, r);
            const onSegment = lineRes.intersect
                && lineRes.parameter >= 0 && lineRes.parameter <= 1;
            expect(segRes.intersect).toBe(onSegment);
            if (onSegment) {
                expect(segRes.parameter).toBe(lineRes.parameter);
                expect(segRes.point.equals(lineRes.point)).toBe(true);
                expect(segRes.rectCoord).toEqual(lineRes.rectCoord);
            }
        });
    });

    it('intersect is true when the two endpoints are strictly on opposite'
        + ' sides of the plane and the crossing is well inside', () => {
        check(arbPair, ([s, r]) => {
            const n = cross(r.axis[0], r.axis[1]);
            const h0 = dot(n, sub(s.p[0], r.center));
            const h1 = dot(n, sub(s.p[1], r.center));
            if (!(h0 * h1 < 0)) {
                return;
            }
            // The crossing point of the segment with the plane of the
            // rectangle.
            const t = h0 / (h0 - h1);
            const p = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
            const delta = sub(p, r.center);
            let strictlyInside = true;
            for (let i = 0; i < 2 && strictlyInside; ++i) {
                strictlyInside = Math.abs(dot(r.axis[i], delta))
                    < r.extent.get(i) * (1 - 1e-6);
            }
            if (strictlyInside) {
                const res = fi.find(s, r);
                expect(res.intersect).toBe(true);
                expectVectorClose(res.point, p, 1e-6, 1e-6);
            }
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(arbPair, rotationFrame(3), wellScaledVector(3, -4, 4)),
            ([[s, r], frame, shift]) => {
                const rot = (v: Vector): Vector => Vector.fromArray([
                    dot(frame[0], v), dot(frame[1], v), dot(frame[2], v)]);
                const map = (v: Vector): Vector => add(rot(v), shift);
                const s2 = Segment.fromEndpoints(map(s.p[0]), map(s.p[1]));
                const r2 = Rectangle.fromCenterAxisExtent(map(r.center),
                    [rot(r.axis[0]), rot(r.axis[1])], r.extent.clone());
                const a = fi.find(s, r), b = fi.find(s2, r2);
                // Only assert for configurations where the crossing is well
                // inside the rectangle; on the boundary the exact comparisons
                // can flip under an arbitrarily small rotation.
                if (a.intersect) {
                    let wellInside = true;
                    for (let i = 0; i < 2; ++i) {
                        wellInside = wellInside
                            && Math.abs(a.rectCoord[i])
                                < r.extent.get(i) * (1 - 1e-4);
                    }
                    wellInside = wellInside
                        && a.parameter > 1e-4 && a.parameter < 1 - 1e-4;
                    if (wellInside) {
                        expect(b.intersect).toBe(true);
                        expectClose(b.parameter, a.parameter, 1e-8, 1e-8);
                        expectVectorClose(b.point, map(a.point), 1e-7, 1e-7);
                    }
                }
            });
    });

    it('reports no intersection for a segment in the plane of the rectangle',
        () => {
            // Upstream documents that a segment lying in the plane of the
            // rectangle is reported as no intersection.
            const inPlane = segment([-3, 0, 0], [3, 0, 0]);
            expect(ti.test(inPlane, square).intersect).toBe(false);
            expect(fi.find(inPlane, square).intersect).toBe(false);
        });
});
