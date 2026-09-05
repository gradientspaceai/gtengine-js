import { describe, it, expect } from 'vitest';
import {
    IntrSegment3Triangle3TI,
    IntrSegment3Triangle3FI
} from '../src/IntrSegment3Triangle3.js';
import {
    IntrLine3Triangle3FI
} from '../src/IntrLine3Triangle3.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { length } from '../src/Vector.js';
import { check, expectVectorClose, fc } from './helpers/arbitraries.js';

const ti = new IntrSegment3Triangle3TI();
const fi = new IntrSegment3Triangle3FI();
const lineFI = new IntrLine3Triangle3FI();

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

const unitTriangle = Triangle.fromVertices(
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));

// An independent oracle: intersect the segment with the triangle plane and
// compute the barycentric coordinates from signed sub-triangle areas. The
// parameter u is in [0,1] along p0 -> p1.
function oracle(segment: Segment, triangle: Triangle):
    { hit: boolean, u: number, bary: number[] } | null {
    const e1 = sub(triangle.v[1], triangle.v[0]);
    const e2 = sub(triangle.v[2], triangle.v[0]);
    const n = cross(e1, e2);
    const dir = sub(segment.p[1], segment.p[0]);
    const den = dot(dir, n);
    if (den === 0) {
        return null;
    }
    const u = dot(sub(triangle.v[0], segment.p[0]), n) / den;
    const p = add(segment.p[0], mul(u, dir));
    const nn = dot(n, n);
    const b1 = dot(cross(sub(p, triangle.v[0]), e2), n) / nn;
    const b2 = dot(cross(e1, sub(p, triangle.v[0])), n) / nn;
    const b0 = 1 - b1 - b2;
    return {
        hit: u >= 0 && u <= 1 && b0 >= 0 && b1 >= 0 && b2 >= 0,
        u,
        bary: [b0, b1, b2]
    };
}

describe('IntrSegment3Triangle3TI', () => {
    it('reports a segment crossing the triangle interior', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -1), v3(0.25, 0.25, 1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(true);
    });

    it('rejects a segment that stops short of the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -3), v3(0.25, 0.25, -1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('accepts a segment whose endpoint touches the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -2), v3(0.25, 0.25, 0));
        expect(ti.test(segment, unitTriangle).intersect).toBe(true);
    });

    it('rejects a segment that starts past the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, 1), v3(0.25, 0.25, 3));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('rejects a segment crossing the plane outside the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.75, 0.75, -1), v3(0.75, 0.75, 1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection when the segment is parallel to the plane', () => {
        const coplanar = Segment.fromEndpoints(v3(-1, 0.25, 0), v3(1, 0.25, 0));
        expect(ti.test(coplanar, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate (zero-length) segment', () => {
        // The centered direction is the zero vector, so Dot(D,N) is zero and
        // the query takes the parallel branch.
        const degenerate = Segment.fromEndpoints(
            v3(0.25, 0.25, 0), v3(0.25, 0.25, 0));
        expect(ti.test(degenerate, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate triangle', () => {
        const degenerate = Triangle.fromVertices(
            v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0));
        const segment = Segment.fromEndpoints(
            v3(0.5, 0, -1), v3(0.5, 0, 1));
        expect(ti.test(segment, degenerate).intersect).toBe(false);
    });
});

describe('IntrSegment3Triangle3FI', () => {
    it('computes the centered parameter, barycentrics and point', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -1), v3(0.25, 0.25, 3));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(true);
        // The centered form has center (0.25,0.25,1), extent 2 and direction
        // (0,0,1), so the intersection is at s = -1.
        expect(result.parameter).toBeCloseTo(-1, 12);
        expect(result.triangleBary[0]).toBeCloseTo(0.5, 12);
        expect(result.triangleBary[1]).toBeCloseTo(0.25, 12);
        expect(result.triangleBary[2]).toBeCloseTo(0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.25, 12);
        expect(result.point.values[1]).toBeCloseTo(0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports the endpoint parameter for a touching endpoint', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -2), v3(0.25, 0.25, 0));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(true);
        const { extent } = segment.getCenteredForm();
        expect(result.parameter).toBeCloseTo(extent, 12);
    });

    it('leaves the default result when there is no intersection', () => {
        const segment = Segment.fromEndpoints(
            v3(5, 5, -1), v3(5, 5, 1));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.triangleBary).toEqual([0, 0, 0]);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('is watertight along an edge shared by two triangles', () => {
        const t0 = unitTriangle;
        const t1 = Triangle.fromVertices(
            v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0));
        for (let i = 0; i <= 10; ++i) {
            const s = i / 10;
            const p = add(mul(1 - s, v3(1, 0, 0)), mul(s, v3(0, 1, 0)));
            const segment = Segment.fromEndpoints(
                add(p, v3(0, 0, -1)), add(p, v3(0, 0, 1)));
            expect(fi.find(segment, t0).intersect
                || fi.find(segment, t1).intersect).toBe(true);
        }
    });
});

describe('IntrSegment3Triangle3 consistency', () => {
    it('agrees between TI, FI, the line query and the oracle', () => {
        const rand = makeRandom(97531);
        const triangle = Triangle.fromVertices(
            v3(0.5, -1, 0.25), v3(2, 0.5, -0.5), v3(-0.75, 1.5, 1));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const p0 = v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5);
            const p1 = v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5);
            const segment = Segment.fromEndpoints(p0, p1);
            const { center, direction, extent } = segment.getCenteredForm();
            if (extent < 1e-6) {
                continue;
            }

            const tiResult = ti.test(segment, triangle);
            const fiResult = fi.find(segment, triangle);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            const expected = oracle(segment, triangle);
            expect(expected).not.toBeNull();
            expect(tiResult.intersect).toBe(expected!.hit);

            // The segment hits exactly when the line hits with |s| <= extent.
            const line = Line.fromOriginDirection(center, direction);
            const lineResult = lineFI.find(line, triangle);
            expect(tiResult.intersect).toBe(lineResult.intersect
                && Math.abs(lineResult.parameter) <= extent);

            if (fiResult.intersect) {
                ++hits;
                expect(Math.abs(fiResult.parameter))
                    .toBeLessThanOrEqual(extent);
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.triangleBary[i])
                        .toBeCloseTo(expected!.bary[i], 8);
                }
                const combo = add(add(
                    mul(fiResult.triangleBary[0], triangle.v[0]),
                    mul(fiResult.triangleBary[1], triangle.v[1])),
                    mul(fiResult.triangleBary[2], triangle.v[2]));
                for (let i = 0; i < 3; ++i) {
                    expect(combo.values[i])
                        .toBeCloseTo(fiResult.point.values[i], 8);
                }
                // The centered parameter maps to the [0,1] parameter of the
                // endpoint form.
                const u = (fiResult.parameter + extent) / (2 * extent);
                expect(u).toBeCloseTo(expected!.u, 8);
            }
        }
        expect(hits).toBeGreaterThan(10);
    });
});

describe('IntrSegment3Triangle3 verification', () => {
    // Integer vertices and endpoints. The query works on the centered form
    // C + s*D with |s| <= e, where D is the normalized edge vector, so only
    // the normalization is inexact; the triangle-side sign tests still come
    // from small integer determinants.
    const latticeTri = fc.array(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        { minLength: 3, maxLength: 3 })
        .map(vs => Triangle.fromVertices(Vector.fromArray(vs[0]),
            Vector.fromArray(vs[1]), Vector.fromArray(vs[2])));
    const latticeSeg = fc.tuple(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }))
        .filter(([a, b]) => a[0] !== b[0] || a[1] !== b[1] || a[2] !== b[2])
        .map(([a, b]) => Segment.fromEndpoints(Vector.fromArray(a),
            Vector.fromArray(b)));

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(latticeSeg, latticeTri), ([s, t]) => {
            expect(fi.find(s, t).intersect).toBe(ti.test(s, t).intersect);
        });
    });

    it('the reported parameter is the centered-form s with |s| <= extent', () => {
        check(fc.tuple(latticeSeg, latticeTri), ([s, t]) => {
            const r = fi.find(s, t);
            if (!r.intersect) {
                return;
            }
            const { center, direction, extent } = s.getCenteredForm();
            // Upstream reports the parameter of C + s*D, not the [0,1]
            // parameter of (1-u)*P0 + u*P1 (see the port's file comments).
            expect(Number.isFinite(r.parameter)).toBe(true);
            expect(Math.abs(r.parameter))
                .toBeLessThanOrEqual(extent * (1 + 1e-12) + 1e-12);
            expectVectorClose(r.point,
                add(center, mul(r.parameter, direction)), 0, 0);
            // The equivalent [0,1] parameter, for callers that need it.
            const u = 0.5 + r.parameter / (2 * extent);
            expect(u).toBeGreaterThanOrEqual(-1e-12);
            expect(u).toBeLessThanOrEqual(1 + 1e-12);
            const scale = 1 + length(s.p[0]) + length(s.p[1]);
            expect(length(sub(r.point,
                add(mul(1 - u, s.p[0]), mul(u, s.p[1])))))
                .toBeLessThanOrEqual(1e-11 * scale);
        });
    });

    it('the barycentric coordinates reconstruct the reported point', () => {
        check(fc.tuple(latticeSeg, latticeTri), ([s, t]) => {
            const r = fi.find(s, t);
            if (!r.intersect) {
                return;
            }
            const [b0, b1, b2] = r.triangleBary;
            expect(b0 + b1 + b2).toBeCloseTo(1, 12);
            expect(b1).toBeGreaterThanOrEqual(0);
            expect(b2).toBeGreaterThanOrEqual(0);
            expect(b0).toBeGreaterThanOrEqual(-1e-12);
            const fromBary = add(add(mul(b0, t.v[0]), mul(b1, t.v[1])),
                mul(b2, t.v[2]));
            const scale = 1 + length(t.v[0]) + length(t.v[1]) + length(t.v[2]);
            expect(length(sub(r.point, fromBary)))
                .toBeLessThanOrEqual(1e-10 * scale);
        });
    });

    it('a segment hit is a line hit whose parameter lies inside the segment', () => {
        check(fc.tuple(latticeSeg, latticeTri), ([s, t]) => {
            const { center, direction, extent } = s.getCenteredForm();
            const lr = lineFI.find(
                Line.fromOriginDirection(center, direction), t);
            const sr = fi.find(s, t);
            if (sr.intersect) {
                expect(lr.intersect).toBe(true);
                expect(sr.parameter).toBeCloseTo(lr.parameter, 9);
            } else if (lr.intersect) {
                expect(Math.abs(lr.parameter))
                    .toBeGreaterThan(extent * (1 - 1e-9) - 1e-9);
            }
        });
    });

    it('reversing the endpoints negates the parameter and keeps the point', () => {
        check(fc.tuple(latticeSeg, latticeTri), ([s, t]) => {
            const a = fi.find(s, t);
            const b = fi.find(Segment.fromEndpoints(s.p[1], s.p[0]), t);
            expect(b.intersect).toBe(a.intersect);
            if (!a.intersect) {
                return;
            }
            expect(b.parameter).toBeCloseTo(-a.parameter, 9);
            expectVectorClose(b.point, a.point, 1e-10, 1e-10);
            for (let i = 0; i < 3; ++i) {
                expect(b.triangleBary[i]).toBeCloseTo(a.triangleBary[i], 9);
            }
        });
    });

    it('a segment crossing a sampled interior point hits at the right place', () => {
        const rnd = makeRandom(0x19fe32);
        check(fc.tuple(latticeTri), ([t]) => {
            const n = cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
            const nlen = Math.sqrt(dot(n, n));
            if (nlen < 1) {
                return;    // degenerate lattice triangle
            }
            const u = mul(1 / nlen, n);
            let a = rnd(), b = rnd();
            if (a + b > 1) {
                a = 1 - a; b = 1 - b;
            }
            a = 0.05 + 0.9 * a;
            b = 0.05 + 0.9 * b;
            if (a + b > 0.95) {
                return;
            }
            const target = add(add(mul(1 - a - b, t.v[0]), mul(a, t.v[1])),
                mul(b, t.v[2]));
            const s = Segment.fromEndpoints(sub(target, mul(2, u)),
                add(target, mul(3, u)));
            const r = fi.find(s, t);
            expect(r.intersect).toBe(true);
            expect(ti.test(s, t).intersect).toBe(true);
            const scale = 1 + length(target);
            expect(length(sub(r.point, target)))
                .toBeLessThanOrEqual(1e-9 * scale);
            // Shrinking the segment so it stops short of the plane must miss.
            const short = Segment.fromEndpoints(sub(target, mul(2, u)),
                sub(target, mul(0.5, u)));
            expect(fi.find(short, t).intersect).toBe(false);
        });
    });

    it('a segment in the plane of the triangle reports no intersection', () => {
        const s = Segment.fromEndpoints(v3(-5, 0.25, 0), v3(5, 0.25, 0));
        expect(ti.test(s, unitTriangle).intersect).toBe(false);
        const r = fi.find(s, unitTriangle);
        expect(r.intersect).toBe(false);
        expect(r.parameter).toBe(0);
        expect(r.triangleBary).toEqual([0, 0, 0]);
    });

    it('an endpoint exactly on the triangle is reported at |s| = extent', () => {
        const s = Segment.fromEndpoints(v3(0.25, 0.25, 0), v3(0.25, 0.25, 4));
        const r = fi.find(s, unitTriangle);
        expect(r.intersect).toBe(true);
        expect(r.parameter).toBeCloseTo(-2, 12);
        expectVectorClose(r.point, v3(0.25, 0.25, 0), 1e-12, 1e-12);
    });
});
