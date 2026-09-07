import { describe, it, expect } from 'vitest';
import { Segment } from '../src/Segment.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    IntrSegment2Triangle2TI,
    IntrSegment2Triangle2FI,
    defaultIntrSegment2Triangle2FIResult,
    intrSegment2Triangle2DoQuery
} from '../src/IntrSegment2Triangle2.js';
import { Line } from '../src/Line.js';
import { IntrLine2Triangle2FI } from '../src/IntrLine2Triangle2.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function triangle(a: number[], b: number[], c: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(a), Vector.fromArray(b),
        Vector.fromArray(c));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Independent containment test using barycentric coordinates.
function inTriangle(p: Vector, t: Triangle, tol: number): boolean {
    const [v0, v1, v2] = t.v;
    const d = (v1.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1]);
    if (Math.abs(d) < 1e-12) {
        return false;
    }
    const b1 = ((p.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (p.values[1] - v0.values[1])) / d;
    const b2 = ((v1.values[0] - v0.values[0]) * (p.values[1] - v0.values[1]) -
        (p.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1])) / d;
    const b0 = 1 - b1 - b2;
    return b0 >= -tol && b1 >= -tol && b2 >= -tol;
}

describe('IntrSegment2Triangle2', () => {
    const ti = new IntrSegment2Triangle2TI();
    const fi = new IntrSegment2Triangle2FI();
    const tri = triangle([0, 0], [4, 0], [0, 4]);

    it('clips a segment crossing the triangle interior', () => {
        // The line y = 1 crosses the triangle from (0,1) to (3,1).
        const s = segment([-2, 1], [6, 1]);
        expect(ti.test(s, tri).intersect).toBe(true);
        const result = fi.find(s, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are for the endpoint form, direction (8,0).
        expect(result.parameter[0]).toBeCloseTo(0.25, 12);
        expect(result.parameter[1]).toBeCloseTo(0.625, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('reports the segment itself when it lies inside the triangle', () => {
        const s = segment([0.5, 0.5], [1, 1]);
        expect(ti.test(s, tri).intersect).toBe(true);
        const result = fi.find(s, tri);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips a segment with one endpoint inside the triangle', () => {
        const s = segment([1, 1], [6, 1]);
        const result = fi.find(s, tri);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('reports a segment endpoint that lies on an edge', () => {
        // The endpoint (2,0) is on the edge from (0,0) to (4,0).
        const s = segment([2, 0], [2, -5]);
        expect(ti.test(s, tri).intersect).toBe(true);
        const result = fi.find(s, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('reports an edge overlap when the segment lies along an edge', () => {
        const s = segment([-1, 0], [5, 0]);
        const result = fi.find(s, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(4, 12);
    });

    it('misses when the segment stops short of the triangle', () => {
        const s = segment([-5, 1], [-1, 1]);
        expect(ti.test(s, tri).intersect).toBe(false);
        const result = fi.find(s, tri);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('misses when the supporting line misses the triangle', () => {
        const s = segment([-5, 6], [5, 6]);
        expect(ti.test(s, tri).intersect).toBe(false);
        expect(fi.find(s, tri).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate triangle', () => {
        // The three vertices are collinear. A segment along the same line
        // gives (n,p,z) = (0,0,3), which upstream treats as no intersection.
        const degenerate = triangle([0, 0], [2, 0], [4, 0]);
        const along = segment([1, 0], [3, 0]);
        expect(ti.test(along, degenerate).intersect).toBe(false);
        expect(fi.find(along, degenerate).intersect).toBe(false);

        // A segment crossing the degenerate triangle transversally gives
        // (n,p,z) = (2,1,0), so upstream reports a (zero-length) crossing.
        const across = segment([1, -1], [1, 1]);
        expect(ti.test(across, degenerate).intersect).toBe(true);
        const result = fi.find(across, degenerate);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(result.point[1].values[1]).toBeCloseTo(0, 12);
    });

    it('the exported DoQuery matches the class query', () => {
        const s = segment([-1, 0.7], [5, 2.1]);
        const direct = defaultIntrSegment2Triangle2FIResult();
        intrSegment2Triangle2DoQuery(s.p[0], sub(s.p[1], s.p[0]), tri, direct);
        const viaClass = fi.find(s, tri);
        expect(direct.intersect).toBe(viaClass.intersect);
        expect(direct.numIntersections).toBe(viaClass.numIntersections);
        expect(direct.parameter[0]).toBeCloseTo(viaClass.parameter[0], 12);
        expect(direct.parameter[1]).toBeCloseTo(viaClass.parameter[1], 12);
    });

    it('agrees with brute-force sampling on random configurations', () => {
        const rnd = makeRandom(1618033);
        let tiFiMismatch = 0;
        let sampleMismatch = 0;
        let pointMismatch = 0;
        let hits = 0;
        const samples = 3000;

        for (let trial = 0; trial < 250; ++trial) {
            const t = Triangle.fromVertices(
                vec(4 * rnd() - 2, 4 * rnd() - 2),
                vec(4 * rnd() - 2, 4 * rnd() - 2),
                vec(4 * rnd() - 2, 4 * rnd() - 2));
            // Skip near-degenerate triangles, where the exact sign tests and
            // the sampled containment test disagree.
            const e0 = sub(t.v[1], t.v[0]);
            const e1 = sub(t.v[2], t.v[0]);
            const area2 = e0.values[0] * e1.values[1] -
                e1.values[0] * e0.values[1];
            if (Math.abs(area2) < 0.5) {
                continue;
            }

            const p0 = vec(6 * rnd() - 3, 6 * rnd() - 3);
            const centroid = mul(1 / 3, add(t.v[0], add(t.v[1], t.v[2])));
            const p1 = add(p0, add(mul(2, sub(centroid, p0)),
                vec(3 * rnd() - 1.5, 3 * rnd() - 1.5)));
            const d = sub(p1, p0);
            if (dot(d, d) < 1e-6) {
                continue;
            }
            const s = Segment.fromEndpoints(p0, p1);

            const tiResult = ti.test(s, t);
            const fiResult = fi.find(s, t);
            if (tiResult.intersect !== fiResult.intersect) {
                ++tiFiMismatch;
            }

            let sampled = false;
            for (let k = 0; k <= samples; ++k) {
                if (inTriangle(add(p0, mul(k / samples, d)), t, 0)) {
                    sampled = true;
                    break;
                }
            }
            if (sampled && !fiResult.intersect) {
                ++sampleMismatch;
            }

            if (fiResult.intersect) {
                ++hits;
                for (let i = 0; i < fiResult.numIntersections; ++i) {
                    const u = fiResult.parameter[i];
                    if (u < -1e-9 || u > 1 + 1e-9) {
                        ++pointMismatch;
                    }
                    const expected = add(p0, mul(u, d));
                    const diff = sub(fiResult.point[i], expected);
                    if (Math.sqrt(dot(diff, diff)) > 1e-9) {
                        ++pointMismatch;
                    }
                    if (!inTriangle(expected, t, 1e-7)) {
                        ++pointMismatch;
                    }
                }
            }
        }

        expect(hits).toBeGreaterThan(20);
        expect([tiFiMismatch, sampleMismatch, pointMismatch]).toEqual([0, 0, 0]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrSegment2Triangle2.h.
// ---------------------------------------------------------------------------

describe('IntrSegment2Triangle2 verification', () => {
    const ti = new IntrSegment2Triangle2TI();
    const fi = new IntrSegment2Triangle2FI();
    const lineFi = new IntrLine2Triangle2FI();

    const arbSegment = fc.tuple(wellScaledVector(2, -6, 6),
        wellScaledVector(2, -6, 6))
        .filter(([a, b]) => {
            const d = sub(b, a);
            return dot(d, d) > 1e-2;
        })
        .map(([a, b]) => Segment.fromEndpoints(a, b));
    const arbTriangle = fc.tuple(wellScaledVector(2, -5, 5),
        wellScaledVector(2, -5, 5), wellScaledVector(2, -5, 5))
        .filter(([a, b, c]) => {
            const e0 = sub(b, a), e1 = sub(c, a);
            return Math.abs(e0.get(0) * e1.get(1) - e0.get(1) * e1.get(0)) > 0.5;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));
    const arbPair = fc.tuple(arbSegment, arbTriangle);

    // Barycentric coordinates of p with respect to the triangle.
    function bary(t: Triangle, p: Vector): [number, number, number] {
        const v0 = sub(t.v[1], t.v[0]), v1 = sub(t.v[2], t.v[0]);
        const v2 = sub(p, t.v[0]);
        const den = v0.get(0) * v1.get(1) - v0.get(1) * v1.get(0);
        const b1 = (v2.get(0) * v1.get(1) - v2.get(1) * v1.get(0)) / den;
        const b2 = (v0.get(0) * v2.get(1) - v0.get(1) * v2.get(0)) / den;
        return [1 - b1 - b2, b1, b2];
    }

    it('TI and FI agree on intersect (TI delegates to FI upstream)', () => {
        check(arbPair, ([s, t]) => {
            expect(ti.test(s, t).intersect).toBe(fi.find(s, t).intersect);
        });
    });

    it('FI parameters are in [0,1] and the points are on the segment and'
        + ' in the triangle', () => {
        check(arbPair, ([s, t]) => {
            const res = fi.find(s, t);
            if (!res.intersect) {
                expect(res.numIntersections).toBe(0);
                return;
            }
            expect(res.numIntersections === 1 || res.numIntersections === 2)
                .toBe(true);
            const d = sub(s.p[1], s.p[0]);
            for (let i = 0; i < res.numIntersections; ++i) {
                const p = res.parameter[i];
                expect(Number.isFinite(p)).toBe(true);
                expect(p).toBeGreaterThanOrEqual(-1e-12);
                expect(p).toBeLessThanOrEqual(1 + 1e-12);
                expectVectorClose(res.point[i], add(s.p[0], mul(p, d)),
                    1e-12, 1e-12);
                const b = bary(t, res.point[i]);
                for (const bi of b) {
                    expect(bi).toBeGreaterThanOrEqual(-1e-6);
                    expect(bi).toBeLessThanOrEqual(1 + 1e-6);
                }
            }
            expect(res.parameter[0]).toBeLessThanOrEqual(res.parameter[1]);
        });
    });

    it('the segment result is the line result clipped to [0,1]', () => {
        check(arbPair, ([s, t]) => {
            const d = sub(s.p[1], s.p[0]);
            const lineRes = lineFi.find(
                Line.fromOriginDirection(s.p[0], d), t);
            const segRes = fi.find(s, t);
            if (!lineRes.intersect) {
                expect(segRes.intersect).toBe(false);
                return;
            }
            const c0 = Math.max(lineRes.parameter[0], 0);
            const c1 = Math.min(lineRes.parameter[1], 1);
            if (c0 > c1) {
                expect(segRes.intersect).toBe(false);
                return;
            }
            expect(segRes.intersect).toBe(true);
            expect(segRes.parameter[0] + 0).toBe(c0 + 0);
            expect(segRes.parameter[1] + 0).toBe(c1 + 0);
            expect(segRes.numIntersections).toBe(c0 < c1 ? 2 : 1);
        });
    });

    it('intersect is true whenever a finely sampled segment point is well'
        + ' inside the triangle', () => {
        check(arbPair, ([s, t]) => {
            const d = sub(s.p[1], s.p[0]);
            let inside = false;
            for (let k = 0; k <= 512 && !inside; ++k) {
                const b = bary(t, add(s.p[0], mul(k / 512, d)));
                inside = b[0] > 1e-6 && b[1] > 1e-6 && b[2] > 1e-6;
            }
            if (inside) {
                expect(ti.test(s, t).intersect).toBe(true);
            }
        });
    });

    // A configuration is generic when no segment endpoint is close to a
    // triangle edge line, no triangle vertex is close to the segment line and
    // the segment is not nearly parallel to an edge. Rigid-motion
    // equivariance is only meaningful for generic configurations: the query
    // classifies with exact sign tests, so a segment lying along an edge
    // changes its answer under an arbitrarily small rotation.
    function generic(s: Segment, t: Triangle): boolean {
        const cross2 = (a: Vector, b: Vector): number =>
            a.get(0) * b.get(1) - a.get(1) * b.get(0);
        const ds = sub(s.p[1], s.p[0]);
        const lenS = Math.sqrt(dot(ds, ds));
        for (let i = 0; i < 3; ++i) {
            const a = t.v[i], b = t.v[(i + 1) % 3];
            const de = sub(b, a);
            const lenE = Math.sqrt(dot(de, de));
            if (Math.abs(cross2(ds, de)) / (lenS * lenE) < 1e-2) {
                return false;
            }
            for (const p of [s.p[0], s.p[1]]) {
                if (Math.abs(cross2(de, sub(p, a))) / lenE < 1e-2) {
                    return false;
                }
            }
            if (Math.abs(cross2(ds, sub(a, s.p[0]))) / lenS < 1e-2) {
                return false;
            }
        }
        return true;
    }

    it('is equivariant under rigid motions of the plane', () => {
        check(fc.tuple(arbPair, rotationFrame(2), wellScaledVector(2, -4, 4)),
            ([[s, t], frame, shift]) => {
                if (!generic(s, t)) {
                    return;
                }
                const map = (v: Vector): Vector => Vector.fromArray([
                    dot(frame[0], v) + shift.get(0),
                    dot(frame[1], v) + shift.get(1)]);
                const s2 = Segment.fromEndpoints(map(s.p[0]), map(s.p[1]));
                const t2 = Triangle.fromVertices(map(t.v[0]), map(t.v[1]),
                    map(t.v[2]));
                const a = fi.find(s, t), b = fi.find(s2, t2);
                expect(b.intersect).toBe(a.intersect);
                if (a.intersect) {
                    expect(b.numIntersections).toBe(a.numIntersections);
                    for (let i = 0; i < a.numIntersections; ++i) {
                        expectClose(b.parameter[i], a.parameter[i], 1e-9, 1e-9);
                        expectVectorClose(b.point[i], map(a.point[i]),
                            1e-8, 1e-8);
                    }
                }
            });
    });

    it('reports a contained segment as its own endpoints', () => {
        const t = triangle([0, 0], [6, 0], [0, 6]);
        const s = segment([1, 1], [2, 2]);
        const res = fi.find(s, t);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(2);
        expect(res.parameter[0]).toBe(0);
        expect(res.parameter[1]).toBe(1);
        expectVectorClose(res.point[0], vec(1, 1));
        expectVectorClose(res.point[1], vec(2, 2));
    });

    it('reports a single point for a segment touching one vertex', () => {
        const t = triangle([0, 0], [6, 0], [0, 6]);
        const s = segment([-2, 2], [0, 0]);
        const res = fi.find(s, t);
        expect(res.intersect).toBe(true);
        expect(res.numIntersections).toBe(1);
        expectVectorClose(res.point[0], vec(0, 0));
    });

    it('the DoQuery helper resets a result that misses the segment domain',
        () => {
            const res = defaultIntrSegment2Triangle2FIResult();
            // The line hits the triangle only for t > 1.
            intrSegment2Triangle2DoQuery(vec(-10, 1), vec(1, 0),
                triangle([0, 0], [6, 0], [0, 6]), res);
            expect(res.intersect).toBe(false);
            expect(res.numIntersections).toBe(0);
            expect(res.parameter).toEqual([0, 0]);
        });
});
