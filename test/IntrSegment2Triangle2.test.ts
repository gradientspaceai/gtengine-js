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
