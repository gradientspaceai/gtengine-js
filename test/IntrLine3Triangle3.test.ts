import { describe, it, expect } from 'vitest';
import {
    IntrLine3Triangle3TI,
    IntrLine3Triangle3FI
} from '../src/IntrLine3Triangle3.js';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { length } from '../src/Vector.js';
import { dotCross } from '../src/Vector3.js';
import { check, expectVectorClose, fc, unitVector } from './helpers/arbitraries.js';

const ti = new IntrLine3Triangle3TI();
const fi = new IntrLine3Triangle3FI();

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

// The reference triangle in the z = 0 plane.
const unitTriangle = Triangle.fromVertices(
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));

// An independent oracle: intersect the line with the triangle plane and
// compute the barycentric coordinates from signed sub-triangle areas.
function oracle(line: Line, triangle: Triangle):
    { hit: boolean, t: number, bary: number[] } | null {
    const e1 = sub(triangle.v[1], triangle.v[0]);
    const e2 = sub(triangle.v[2], triangle.v[0]);
    const n = cross(e1, e2);
    const den = dot(line.direction, n);
    if (den === 0) {
        return null;
    }
    const t = dot(sub(triangle.v[0], line.origin), n) / den;
    const p = add(line.origin, mul(t, line.direction));
    const nn = dot(n, n);
    const b1 = dot(cross(sub(p, triangle.v[0]), e2), n) / nn;
    const b2 = dot(cross(e1, sub(p, triangle.v[0])), n) / nn;
    const b0 = 1 - b1 - b2;
    return { hit: b0 >= 0 && b1 >= 0 && b2 >= 0, t, bary: [b0, b1, b2] };
}

describe('IntrLine3Triangle3TI', () => {
    it('reports a line through the triangle interior', () => {
        const line = makeLine(v3(0.25, 0.25, -1), v3(0, 0, 1));
        expect(ti.test(line, unitTriangle).intersect).toBe(true);
    });

    it('reports a line through a vertex', () => {
        for (const v of unitTriangle.v) {
            const line = makeLine(add(v, v3(0, 0, -1)), v3(0, 0, 1));
            expect(ti.test(line, unitTriangle).intersect).toBe(true);
        }
    });

    it('reports a line through an edge midpoint', () => {
        const mid = mul(0.5, add(unitTriangle.v[1], unitTriangle.v[2]));
        const line = makeLine(add(mid, v3(0, 0, -1)), v3(0, 0, 1));
        expect(ti.test(line, unitTriangle).intersect).toBe(true);
    });

    it('rejects a line just outside the hypotenuse', () => {
        const line = makeLine(v3(0.51, 0.51, -1), v3(0, 0, 1));
        expect(ti.test(line, unitTriangle).intersect).toBe(false);
    });

    it('rejects a line with a negative barycentric coordinate', () => {
        expect(ti.test(makeLine(v3(-0.01, 0.5, -1), v3(0, 0, 1)),
            unitTriangle).intersect).toBe(false);
        expect(ti.test(makeLine(v3(0.5, -0.01, -1), v3(0, 0, 1)),
            unitTriangle).intersect).toBe(false);
    });

    it('reports the same result for the two line orientations', () => {
        const up = makeLine(v3(0.25, 0.25, -1), v3(0, 0, 1));
        const down = makeLine(v3(0.25, 0.25, -1), v3(0, 0, -1));
        expect(ti.test(up, unitTriangle).intersect).toBe(true);
        expect(ti.test(down, unitTriangle).intersect).toBe(true);
    });

    it('ignores the line origin position (infinite extent)', () => {
        const behind = makeLine(v3(0.25, 0.25, 100), v3(0, 0, 1));
        expect(ti.test(behind, unitTriangle).intersect).toBe(true);
    });

    it('reports no intersection when the line is parallel to the triangle plane', () => {
        // Coplanar and passing through the interior.
        const coplanar = makeLine(v3(-1, 0.25, 0), v3(1, 0, 0));
        expect(ti.test(coplanar, unitTriangle).intersect).toBe(false);
        // Parallel and offset.
        const offset = makeLine(v3(-1, 0.25, 1), v3(1, 0, 0));
        expect(ti.test(offset, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate triangle (zero normal)', () => {
        const degenerate = Triangle.fromVertices(
            v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0));
        const line = makeLine(v3(0.5, 0, -1), v3(0, 0, 1));
        expect(ti.test(line, degenerate).intersect).toBe(false);
    });
});

describe('IntrLine3Triangle3FI', () => {
    it('computes the parameter, barycentric coordinates and point', () => {
        const line = makeLine(v3(0.25, 0.25, -1), v3(0, 0, 1));
        const result = fi.find(line, unitTriangle);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(1, 12);
        expect(result.triangleBary[0]).toBeCloseTo(0.5, 12);
        expect(result.triangleBary[1]).toBeCloseTo(0.25, 12);
        expect(result.triangleBary[2]).toBeCloseTo(0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.25, 12);
        expect(result.point.values[1]).toBeCloseTo(0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports barycentric (1,0,0) at vertex 0', () => {
        const line = makeLine(v3(0, 0, -1), v3(0, 0, 1));
        const result = fi.find(line, unitTriangle);
        expect(result.intersect).toBe(true);
        expect(result.triangleBary).toEqual([1, 0, 0]);
    });

    it('reports a negative parameter when the triangle is behind the origin', () => {
        const line = makeLine(v3(0.25, 0.25, 2), v3(0, 0, 1));
        const result = fi.find(line, unitTriangle);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(-2, 12);
    });

    it('leaves the default result when there is no intersection', () => {
        const line = makeLine(v3(2, 2, -1), v3(0, 0, 1));
        const result = fi.find(line, unitTriangle);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.triangleBary).toEqual([0, 0, 0]);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('is watertight along an edge shared by two triangles', () => {
        // Two triangles sharing the edge from (1,0,0) to (0,1,0).
        const t0 = unitTriangle;
        const t1 = Triangle.fromVertices(
            v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0));
        for (let i = 0; i <= 10; ++i) {
            const s = i / 10;
            const p = add(mul(1 - s, v3(1, 0, 0)), mul(s, v3(0, 1, 0)));
            const line = makeLine(add(p, v3(0, 0, -1)), v3(0, 0, 1));
            const hit0 = fi.find(line, t0).intersect;
            const hit1 = fi.find(line, t1).intersect;
            expect(hit0 || hit1).toBe(true);
        }
    });
});

describe('IntrLine3Triangle3 consistency', () => {
    it('agrees between TI, FI and the barycentric oracle', () => {
        const rand = makeRandom(13579);
        const triangle = Triangle.fromVertices(
            v3(0.5, -1, 0.25), v3(2, 0.5, -0.5), v3(-0.75, 1.5, 1));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const origin = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const direction = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(direction, direction) < 1e-8) {
                continue;
            }
            const line = makeLine(origin, direction);

            const tiResult = ti.test(line, triangle);
            const fiResult = fi.find(line, triangle);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            const expected = oracle(line, triangle);
            expect(expected).not.toBeNull();
            expect(tiResult.intersect).toBe(expected!.hit);

            if (fiResult.intersect) {
                ++hits;
                expect(fiResult.parameter).toBeCloseTo(expected!.t, 8);
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.triangleBary[i])
                        .toBeCloseTo(expected!.bary[i], 8);
                }
                // The point equals the barycentric combination of vertices.
                const combo = add(add(
                    mul(fiResult.triangleBary[0], triangle.v[0]),
                    mul(fiResult.triangleBary[1], triangle.v[1])),
                    mul(fiResult.triangleBary[2], triangle.v[2]));
                for (let i = 0; i < 3; ++i) {
                    expect(combo.values[i])
                        .toBeCloseTo(fiResult.point.values[i], 8);
                }
            }
        }
        expect(hits).toBeGreaterThan(10);
    });
});

describe('IntrLine3Triangle3 verification', () => {
    // Integer vertices, origins and directions: DotCross is a small integer
    // determinant, so the three sign tests the query branches on are exact
    // and the edge/vertex hits are ordinary draws rather than accidents.
    const latticeTri = fc.array(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        { minLength: 3, maxLength: 3 })
        .map(vs => Triangle.fromVertices(Vector.fromArray(vs[0]),
            Vector.fromArray(vs[1]), Vector.fromArray(vs[2])));
    const latticeLine = fc.tuple(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        fc.array(fc.integer({ min: -3, max: 3 }), { minLength: 3, maxLength: 3 }))
        .filter(([, d]) => d[0] !== 0 || d[1] !== 0 || d[2] !== 0)
        .map(([p, d]) => Line.fromOriginDirection(Vector.fromArray(p),
            Vector.fromArray(d)));

    it('TI and FI agree, and match the exact sign classification', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const diff = sub(l.origin, t.v[0]);
            const e1 = sub(t.v[1], t.v[0]);
            const e2 = sub(t.v[2], t.v[0]);
            const DdN = dot(l.direction, cross(e1, e2));
            const sign = DdN > 0 ? 1 : (DdN < 0 ? -1 : 0);
            let expected: boolean;
            if (sign === 0) {
                expected = false;
            } else {
                const b1 = sign * dotCross(l.direction, diff, e2);
                const b2 = sign * dotCross(l.direction, e1, diff);
                expected = b1 >= 0 && b2 >= 0 && b1 + b2 <= sign * DdN;
            }
            expect(ti.test(l, t).intersect).toBe(expected);
            expect(fi.find(l, t).intersect).toBe(expected);
        });
    });

    it('the barycentric coordinates reconstruct the reported point', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const r = fi.find(l, t);
            if (!r.intersect) {
                return;
            }
            const [b0, b1, b2] = r.triangleBary;
            expect(Number.isFinite(r.parameter)).toBe(true);
            expect(b0 + b1 + b2).toBeCloseTo(1, 12);
            expect(b0).toBeGreaterThanOrEqual(-1e-12);
            expect(b1).toBeGreaterThanOrEqual(0);
            expect(b2).toBeGreaterThanOrEqual(0);
            expectVectorClose(r.point,
                add(l.origin, mul(r.parameter, l.direction)), 0, 0);
            const fromBary = add(add(mul(b0, t.v[0]), mul(b1, t.v[1])),
                mul(b2, t.v[2]));
            // Both forms divide by the same exactly computed Dot(D,N), so
            // the residual is a few ulps of the vertex magnitudes.
            const scale = 1 + length(t.v[0]) + length(t.v[1]) + length(t.v[2])
                + Math.abs(r.parameter);
            expect(length(sub(r.point, fromBary)))
                .toBeLessThanOrEqual(1e-11 * scale);
        });
    });

    it('a line through a sampled interior point always hits', () => {
        const rnd = makeRandom(0x51c07d);
        check(fc.tuple(latticeTri, unitVector(3)), ([t, d]) => {
            const n = cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
            if (length(n) < 1) {
                return;    // degenerate lattice triangle
            }
            normalize(n);
            if (Math.abs(dot(d, n)) < 1e-2) {
                return;    // nearly parallel: upstream reports no hit
            }
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
            const l = Line.fromOriginDirection(sub(target, mul(4, d)), d);
            const r = fi.find(l, t);
            expect(r.intersect).toBe(true);
            expect(ti.test(l, t).intersect).toBe(true);
            const scale = 1 + length(target);
            expect(length(sub(r.point, target)))
                .toBeLessThanOrEqual(1e-9 * scale);
            expect(r.triangleBary[1]).toBeCloseTo(a, 8);
            expect(r.triangleBary[2]).toBeCloseTo(b, 8);
        });
    });

    it('reversing the line direction keeps the hit and negates the parameter', () => {
        check(fc.tuple(latticeLine, latticeTri), ([l, t]) => {
            const a = fi.find(l, t);
            const b = fi.find(Line.fromOriginDirection(l.origin,
                mul(-1, l.direction)), t);
            expect(b.intersect).toBe(a.intersect);
            if (!a.intersect) {
                return;
            }
            expect(b.parameter).toBeCloseTo(-a.parameter, 9);
            expectVectorClose(b.point, a.point, 1e-11, 1e-11);
            for (let i = 0; i < 3; ++i) {
                expect(b.triangleBary[i]).toBeCloseTo(a.triangleBary[i], 9);
            }
        });
    });

    it('a line in the plane of the triangle reports no intersection', () => {
        // Upstream calls Dot(D,N) = 0 a "no intersection" even when the line
        // is coplanar and crosses the triangle.
        const l = makeLine(v3(-5, 0.25, 0), v3(1, 0, 0));
        expect(ti.test(l, unitTriangle).intersect).toBe(false);
        const r = fi.find(l, unitTriangle);
        expect(r.intersect).toBe(false);
        expect(r.parameter).toBe(0);
        expect(r.triangleBary).toEqual([0, 0, 0]);
        expect(r.point.values).toEqual([0, 0, 0]);
    });

    it('a degenerate triangle has a zero normal and never intersects', () => {
        const degenerate = Triangle.fromVertices(v3(0, 0, 0), v3(1, 0, 0),
            v3(2, 0, 0));
        const l = makeLine(v3(0.5, 0, -1), v3(0, 0, 1));
        expect(ti.test(l, degenerate).intersect).toBe(false);
        expect(fi.find(l, degenerate).intersect).toBe(false);
    });
});
