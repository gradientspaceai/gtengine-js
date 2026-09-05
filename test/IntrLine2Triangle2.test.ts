import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { dotPerp } from '../src/Vector2.js';
import {
    IntrLine2Triangle2TI,
    IntrLine2Triangle2FI
} from '../src/IntrLine2Triangle2.js';
import {
    intrLine2Triangle2DoQuery,
    defaultIntrLine2Triangle2FIResult
} from '../src/IntrLine2Triangle2.js';
import { sub } from '../src/Vector.js';
import { check, expectVectorClose, fc, seededRandom } from './helpers/arbitraries.js';

function line(px: number, py: number, dx: number, dy: number): Line {
    return Line.fromOriginDirection(Vector.fromArray([px, py]),
        Vector.fromArray([dx, dy]));
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

// Barycentric coordinates of p with respect to the triangle.
function barycentric(p: Vector, t: Triangle): number[] {
    const [v0, v1, v2] = t.v;
    const d = (v1.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1]);
    const b1 = ((p.values[0] - v0.values[0]) * (v2.values[1] - v0.values[1]) -
        (v2.values[0] - v0.values[0]) * (p.values[1] - v0.values[1])) / d;
    const b2 = ((v1.values[0] - v0.values[0]) * (p.values[1] - v0.values[1]) -
        (p.values[0] - v0.values[0]) * (v1.values[1] - v0.values[1])) / d;
    return [1 - b1 - b2, b1, b2];
}

describe('IntrLine2Triangle2', () => {
    const ti = new IntrLine2Triangle2TI();
    const fi = new IntrLine2Triangle2FI();
    const tri = triangle([0, 0], [4, 0], [0, 4]);

    it('clips a line crossing the interior, (n,p,z) = (2,1,0)', () => {
        // The line y = 1 crosses the triangle from (0,1) to (3,1).
        const l = line(0, 1, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('clips a line through a vertex and the opposite edge, (1,1,1)', () => {
        // The line through (0,0) with direction (1,1) exits at (2,2).
        const l = line(0, 0, 1, 1);
        const result = fi.find(l, tri);
        expect(result.numIntersections).toBe(2);
        // The direction is not unit length, so the parameters are divided by
        // |D|^2 = 2.
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[1]).toBeCloseTo(2, 12);
    });

    it('reports a single point when the line touches a vertex, (2,0,1)', () => {
        const l = line(0, 4, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBe(result.parameter[1]);
        expect(result.point[0].values).toEqual([0, 4]);
    });

    it('reports the whole edge when the line contains one, (0,1,2)', () => {
        const l = line(0, 0, 1, 0);
        expect(ti.test(l, tri).intersect).toBe(true);
        const result = fi.find(l, tri);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(4, 12);
    });

    it('reports no intersection when all vertices are on one side', () => {
        for (const l of [line(0, -1, 1, 0), line(0, 5, 1, 0)]) {
            expect(ti.test(l, tri).intersect).toBe(false);
            const result = fi.find(l, tri);
            expect(result.intersect).toBe(false);
            expect(result.numIntersections).toBe(0);
        }
    });

    it('treats a degenerate triangle on the line as no intersection, (0,0,3)', () => {
        const collinear = triangle([0, 0], [1, 0], [2, 0]);
        const l = line(0, 0, 1, 0);
        expect(ti.test(l, collinear).intersect).toBe(false);
        expect(fi.find(l, collinear).intersect).toBe(false);
    });

    it('orders the parameters and is direction-reversal consistent', () => {
        const forward = fi.find(line(0, 1, 1, 0), tri);
        const backward = fi.find(line(0, 1, -1, 0), tri);
        expect(forward.parameter[0]).toBeLessThanOrEqual(forward.parameter[1]);
        expect(backward.parameter[0]).toBeLessThanOrEqual(backward.parameter[1]);
        // The intersection segment is the same set of points.
        const fPoints = [forward.point[0].values, forward.point[1].values]
            .sort((a, b) => a[0] - b[0]);
        const bPoints = [backward.point[0].values, backward.point[1].values]
            .sort((a, b) => a[0] - b[0]);
        expect(fPoints[0][0]).toBeCloseTo(bPoints[0][0], 12);
        expect(fPoints[1][0]).toBeCloseTo(bPoints[1][0], 12);
    });

    it('keeps TI and FI consistent and clips inside the triangle', () => {
        const rand = makeRandom(13571113);
        let numHit = 0, numMiss = 0;
        for (let trial = 0; trial < 800; ++trial) {
            const t = triangle(
                [4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2],
                [4 * rand() - 2, 4 * rand() - 2]);
            const d = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1]);
            if (normalize(d) < 1e-6) {
                continue;
            }
            const l = Line.fromOriginDirection(
                Vector.fromArray([4 * rand() - 2, 4 * rand() - 2]), d);

            // Independent oracle: the line meets the solid triangle when the
            // vertices are not all strictly on one side.
            const s = [0, 1, 2].map(i =>
                dotPerp(l.direction, Vector.fromArray([
                    t.v[i].values[0] - l.origin.values[0],
                    t.v[i].values[1] - l.origin.values[1]])));
            const oracle = !(s.every(x => x > 0) || s.every(x => x < 0) ||
                s.every(x => x === 0));

            const tiResult = ti.test(l, t).intersect;
            const f = fi.find(l, t);
            expect(tiResult).toBe(oracle);
            expect(f.intersect).toBe(tiResult);

            if (f.intersect) {
                ++numHit;
                expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
                // Sampled points of the reported interval are in the triangle.
                for (let k = 0; k <= 4; ++k) {
                    const u = f.parameter[0] +
                        (k / 4) * (f.parameter[1] - f.parameter[0]);
                    const p = add(l.origin, mul(u, l.direction));
                    for (const b of barycentric(p, t)) {
                        expect(b).toBeGreaterThan(-1e-8);
                    }
                }
            } else {
                ++numMiss;
                expect(f.numIntersections).toBe(0);
            }
        }
        expect(numHit).toBeGreaterThan(50);
        expect(numMiss).toBeGreaterThan(50);
    });
});

describe('intrLine2Triangle2DoQuery', () => {
    const tri = triangle([0, 0], [4, 0], [0, 4]);

    it('matches the class query but does not compute points', () => {
        const l = line(-1, 1, 1, 0);
        const result = defaultIntrLine2Triangle2FIResult();
        intrLine2Triangle2DoQuery(l.origin, l.direction, tri, result);
        const expected = new IntrLine2Triangle2FI().find(l, tri);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter[0]).toBeCloseTo(expected.parameter[0], 12);
        expect(result.parameter[1]).toBeCloseTo(expected.parameter[1], 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('reports no intersection for a line missing the triangle', () => {
        const result = defaultIntrLine2Triangle2FIResult();
        intrLine2Triangle2DoQuery(Vector.fromArray([-1, 5]),
            Vector.fromArray([1, 0]), tri, result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

describe('IntrLine2Triangle2 verification', () => {
    const ti = new IntrLine2Triangle2TI();
    const fi = new IntrLine2Triangle2FI();

    // Integer origins, directions and vertices make every DotPerp exact, so
    // the (n,p,z) classification the queries branch on is decided without
    // round-off and the on-line vertex cases are common draws.
    const latticeLine = fc.tuple(fc.integer({ min: -4, max: 4 }),
        fc.integer({ min: -4, max: 4 }), fc.integer({ min: -3, max: 3 }),
        fc.integer({ min: -3, max: 3 }))
        .filter(([, , dx, dy]) => dx !== 0 || dy !== 0)
        .map(([px, py, dx, dy]) => line(px, py, dx, dy));
    const latticeTri = fc.array(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 2, maxLength: 2 }),
        { minLength: 3, maxLength: 3 })
        .map(vs => triangle(vs[0], vs[1], vs[2]));

    function signs(l: Line, t: Triangle): number[] {
        return t.v.map(v => {
            const s = dotPerp(l.direction, Vector.fromArray(
                [v.values[0] - l.origin.values[0],
                    v.values[1] - l.origin.values[1]]));
            return s > 0 ? 1 : (s < 0 ? -1 : 0);
        });
    }

    it('TI matches the exact (n,p,z) table', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const s = signs(l, t);
            const numZero = s.filter(x => x === 0).length;
            const numPos = s.filter(x => x > 0).length;
            const numNeg = s.filter(x => x < 0).length;
            const expected = (numZero === 0 && numPos > 0 && numNeg > 0)
                || numZero === 1 || numZero === 2;
            expect(ti.test(l, t).intersect).toBe(expected);
            expect(fi.find(l, t).intersect).toBe(expected);
        });
    });

    it('FI parameters are ordered and their points are on the line and in the triangle', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const r = fi.find(l, t);
            if (!r.intersect) {
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            if (r.numIntersections === 1) {
                expect(r.parameter[1]).toBe(r.parameter[0]);
            }
            const e1 = sub(t.v[1], t.v[0]), e2 = sub(t.v[2], t.v[0]);
            const det = e1.values[0] * e2.values[1] - e1.values[1] * e2.values[0];
            for (let k = 0; k < r.numIntersections; ++k) {
                const p = r.point[k];
                expect(Number.isFinite(p.values[0])).toBe(true);
                expect(Number.isFinite(p.values[1])).toBe(true);
                expectVectorClose(p,
                    add(l.origin, mul(r.parameter[k], l.direction)), 0, 0);
                if (det === 0) {
                    continue;    // degenerate triangle: no barycentrics
                }
                const w = sub(p, t.v[0]);
                const b1 = (w.values[0] * e2.values[1]
                    - w.values[1] * e2.values[0]) / det;
                const b2 = (e1.values[0] * w.values[1]
                    - e1.values[1] * w.values[0]) / det;
                // The clip parameters are single quotients of exact integer
                // combinations, so the barycentrics are accurate to ~1e-12
                // relative to the triangle size.
                expect(b1).toBeGreaterThanOrEqual(-1e-9);
                expect(b2).toBeGreaterThanOrEqual(-1e-9);
                expect(b1 + b2).toBeLessThanOrEqual(1 + 1e-9);
            }
        });
    });

    it('the reported interval contains every line parameter inside the triangle', () => {
        const rnd = seededRandom(0x4de210b);
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const e1 = sub(t.v[1], t.v[0]), e2 = sub(t.v[2], t.v[0]);
            const det = e1.values[0] * e2.values[1] - e1.values[1] * e2.values[0];
            if (det === 0) {
                return;
            }
            const r = fi.find(l, t);
            for (let k = 0; k < 400; ++k) {
                const tt = 24 * rnd() - 12;
                const p = add(l.origin, mul(tt, l.direction));
                const w = sub(p, t.v[0]);
                const b1 = (w.values[0] * e2.values[1]
                    - w.values[1] * e2.values[0]) / det;
                const b2 = (e1.values[0] * w.values[1]
                    - e1.values[1] * w.values[0]) / det;
                if (b1 > 1e-9 && b2 > 1e-9 && b1 + b2 < 1 - 1e-9) {
                    expect(r.intersect).toBe(true);
                    expect(tt).toBeGreaterThanOrEqual(r.parameter[0] - 1e-9);
                    expect(tt).toBeLessThanOrEqual(r.parameter[1] + 1e-9);
                }
            }
        }, 60);
    }, 30000);

    it('the parameters are those of the given direction, not of a unit direction', () => {
        check(fc.tuple(latticeLine, latticeTri, fc.integer({ min: 2, max: 4 })),
            ([l, t, k]) => {
                const a = fi.find(l, t);
                if (!a.intersect) {
                    return;
                }
                const scaled = Line.fromOriginDirection(l.origin,
                    mul(k, l.direction));
                const b = fi.find(scaled, t);
                expect(b.intersect).toBe(true);
                expect(b.numIntersections).toBe(a.numIntersections);
                for (let i = 0; i < a.numIntersections; ++i) {
                    expect(b.parameter[i] * k).toBeCloseTo(a.parameter[i], 9);
                    expectVectorClose(b.point[i], a.point[i], 1e-12, 1e-12);
                }
            });
    });

    it('a degenerate triangle on the line reports no intersection', () => {
        // (n,p,z) = (0,0,3): all three vertices project to zero.
        const l = line(0, 0, 1, 0);
        const t = triangle([0, 0], [2, 0], [5, 0]);
        expect(ti.test(l, t).intersect).toBe(false);
        const r = fi.find(l, t);
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
    });

    it('the exported DoQuery reproduces the class result', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const r = defaultIntrLine2Triangle2FIResult();
            intrLine2Triangle2DoQuery(l.origin, l.direction, t, r);
            const expected = fi.find(l, t);
            expect(r.intersect).toBe(expected.intersect);
            expect(r.numIntersections).toBe(expected.numIntersections);
            expect(r.parameter[0]).toBe(expected.parameter[0]);
            expect(r.parameter[1]).toBe(expected.parameter[1]);
        });
    });
});
