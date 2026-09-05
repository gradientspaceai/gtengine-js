import { describe, it, expect } from 'vitest';
import {
    IntrRay3Triangle3TI,
    IntrRay3Triangle3FI
} from '../src/IntrRay3Triangle3.js';
import {
    IntrLine3Triangle3TI,
    IntrLine3Triangle3FI
} from '../src/IntrLine3Triangle3.js';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import { length } from '../src/Vector.js';
import { dotCross } from '../src/Vector3.js';
import { check, expectVectorClose, fc, unitVector } from './helpers/arbitraries.js';

const ti = new IntrRay3Triangle3TI();
const fi = new IntrRay3Triangle3FI();
const lineTI = new IntrLine3Triangle3TI();
const lineFI = new IntrLine3Triangle3FI();

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRay(origin: Vector, direction: Vector): Ray {
    const d = direction.clone();
    normalize(d);
    return Ray.fromOriginDirection(origin, d);
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

// An independent oracle: intersect the ray with the triangle plane and
// compute the barycentric coordinates from signed sub-triangle areas.
function oracle(ray: Ray, triangle: Triangle):
    { hit: boolean, t: number, bary: number[] } | null {
    const e1 = sub(triangle.v[1], triangle.v[0]);
    const e2 = sub(triangle.v[2], triangle.v[0]);
    const n = cross(e1, e2);
    const den = dot(ray.direction, n);
    if (den === 0) {
        return null;
    }
    const t = dot(sub(triangle.v[0], ray.origin), n) / den;
    const p = add(ray.origin, mul(t, ray.direction));
    const nn = dot(n, n);
    const b1 = dot(cross(sub(p, triangle.v[0]), e2), n) / nn;
    const b2 = dot(cross(e1, sub(p, triangle.v[0])), n) / nn;
    const b0 = 1 - b1 - b2;
    return {
        hit: t >= 0 && b0 >= 0 && b1 >= 0 && b2 >= 0,
        t,
        bary: [b0, b1, b2]
    };
}

describe('IntrRay3Triangle3TI', () => {
    it('reports a ray pointing at the triangle interior', () => {
        const ray = makeRay(v3(0.25, 0.25, -1), v3(0, 0, 1));
        expect(ti.test(ray, unitTriangle).intersect).toBe(true);
    });

    it('rejects a ray pointing away from the triangle', () => {
        const ray = makeRay(v3(0.25, 0.25, -1), v3(0, 0, -1));
        expect(ti.test(ray, unitTriangle).intersect).toBe(false);
    });

    it('accepts a ray whose origin lies on the triangle (t = 0)', () => {
        const ray = makeRay(v3(0.25, 0.25, 0), v3(0, 0, 1));
        expect(ti.test(ray, unitTriangle).intersect).toBe(true);
    });

    it('accepts hits at vertices and edge midpoints', () => {
        for (const v of unitTriangle.v) {
            const ray = makeRay(add(v, v3(0, 0, -1)), v3(0, 0, 1));
            expect(ti.test(ray, unitTriangle).intersect).toBe(true);
        }
        const mid = mul(0.5, add(unitTriangle.v[1], unitTriangle.v[2]));
        const ray = makeRay(add(mid, v3(0, 0, -1)), v3(0, 0, 1));
        expect(ti.test(ray, unitTriangle).intersect).toBe(true);
    });

    it('rejects a ray missing the triangle laterally', () => {
        const ray = makeRay(v3(0.6, 0.6, -1), v3(0, 0, 1));
        expect(ti.test(ray, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection when the ray is parallel to the plane', () => {
        const coplanar = makeRay(v3(-1, 0.25, 0), v3(1, 0, 0));
        expect(ti.test(coplanar, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate triangle', () => {
        const degenerate = Triangle.fromVertices(
            v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2));
        const ray = makeRay(v3(1, 1, 0), v3(0, 0, 1));
        expect(ti.test(ray, degenerate).intersect).toBe(false);
    });
});

describe('IntrRay3Triangle3FI', () => {
    it('computes the parameter, barycentric coordinates and point', () => {
        const ray = makeRay(v3(0.25, 0.25, -2), v3(0, 0, 1));
        const result = fi.find(ray, unitTriangle);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(2, 12);
        expect(result.triangleBary[0]).toBeCloseTo(0.5, 12);
        expect(result.triangleBary[1]).toBeCloseTo(0.25, 12);
        expect(result.triangleBary[2]).toBeCloseTo(0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.25, 12);
        expect(result.point.values[1]).toBeCloseTo(0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports a zero parameter for an origin on the triangle', () => {
        const ray = makeRay(v3(0.25, 0.25, 0), v3(0, 0, 1));
        const result = fi.find(ray, unitTriangle);
        expect(result.intersect).toBe(true);
        expect(result.parameter).toBeCloseTo(0, 15);
    });

    it('leaves the default result when there is no intersection', () => {
        const ray = makeRay(v3(0.25, 0.25, -1), v3(0, 0, -1));
        const result = fi.find(ray, unitTriangle);
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
            const ray = makeRay(add(p, v3(0, 0, -1)), v3(0, 0, 1));
            expect(fi.find(ray, t0).intersect || fi.find(ray, t1).intersect)
                .toBe(true);
        }
    });
});

describe('IntrRay3Triangle3 consistency', () => {
    it('agrees with the line query for nonnegative parameters', () => {
        const rand = makeRandom(24680);
        const triangle = Triangle.fromVertices(
            v3(0.5, -1, 0.25), v3(2, 0.5, -0.5), v3(-0.75, 1.5, 1));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const origin = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const direction = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(direction, direction) < 1e-8) {
                continue;
            }
            const ray = makeRay(origin, direction);
            const line = Line.fromOriginDirection(ray.origin, ray.direction);

            const tiResult = ti.test(ray, triangle);
            const fiResult = fi.find(ray, triangle);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            const lineFIResult = lineFI.find(line, triangle);
            expect(lineTI.test(line, triangle).intersect)
                .toBe(lineFIResult.intersect);
            // The ray hits exactly when the line hits with t >= 0.
            expect(tiResult.intersect).toBe(
                lineFIResult.intersect && lineFIResult.parameter >= 0);

            const expected = oracle(ray, triangle);
            expect(expected).not.toBeNull();
            expect(tiResult.intersect).toBe(expected!.hit);

            if (fiResult.intersect) {
                ++hits;
                expect(fiResult.parameter).toBeGreaterThanOrEqual(0);
                expect(fiResult.parameter).toBeCloseTo(expected!.t, 8);
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.triangleBary[i])
                        .toBeCloseTo(expected!.bary[i], 8);
                }
                const onRay = add(ray.origin,
                    mul(fiResult.parameter, ray.direction));
                for (let i = 0; i < 3; ++i) {
                    expect(onRay.values[i])
                        .toBeCloseTo(fiResult.point.values[i], 10);
                }
            }
        }
        expect(hits).toBeGreaterThan(10);
    });
});

describe('IntrRay3Triangle3 verification', () => {
    // Integer data keeps every DotCross an exact small determinant, so the
    // sign tests the query branches on are decided without round-off.
    const latticeTri = fc.array(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        { minLength: 3, maxLength: 3 })
        .map(vs => Triangle.fromVertices(Vector.fromArray(vs[0]),
            Vector.fromArray(vs[1]), Vector.fromArray(vs[2])));
    const latticeRay = fc.tuple(
        fc.array(fc.integer({ min: -4, max: 4 }), { minLength: 3, maxLength: 3 }),
        fc.array(fc.integer({ min: -3, max: 3 }), { minLength: 3, maxLength: 3 }))
        .filter(([, d]) => d[0] !== 0 || d[1] !== 0 || d[2] !== 0)
        .map(([p, d]) => Ray.fromOriginDirection(Vector.fromArray(p),
            Vector.fromArray(d)));

    it('TI and FI agree and match the exact sign classification', () => {
        check(fc.tuple(latticeRay, latticeTri), ([r, t]) => {
            const diff = sub(r.origin, t.v[0]);
            const e1 = sub(t.v[1], t.v[0]);
            const e2 = sub(t.v[2], t.v[0]);
            const normal = cross(e1, e2);
            const DdN = dot(r.direction, normal);
            const sign = DdN > 0 ? 1 : (DdN < 0 ? -1 : 0);
            let expected: boolean;
            if (sign === 0) {
                expected = false;
            } else {
                const b1 = sign * dotCross(r.direction, diff, e2);
                const b2 = sign * dotCross(r.direction, e1, diff);
                expected = b1 >= 0 && b2 >= 0 && b1 + b2 <= sign * DdN
                    && -sign * dot(diff, normal) >= 0;
            }
            expect(ti.test(r, t).intersect).toBe(expected);
            expect(fi.find(r, t).intersect).toBe(expected);
        });
    });

    it('a ray hit is a line hit with a nonnegative parameter', () => {
        check(fc.tuple(latticeRay, latticeTri), ([r, t]) => {
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lr = lineFI.find(l, t);
            const rr = fi.find(r, t);
            expect(lineTI.test(l, t).intersect).toBe(lr.intersect);
            if (rr.intersect) {
                expect(lr.intersect).toBe(true);
                expect(rr.parameter).toBeGreaterThanOrEqual(0);
                // The two queries divide by the same exact Dot(D,N).
                expect(rr.parameter).toBeCloseTo(lr.parameter, 9);
                expectVectorClose(rr.point, lr.point, 1e-11, 1e-11);
                for (let i = 0; i < 3; ++i) {
                    expect(rr.triangleBary[i])
                        .toBeCloseTo(lr.triangleBary[i], 9);
                }
            } else if (lr.intersect) {
                expect(lr.parameter).toBeLessThan(0);
            }
        });
    });

    it('the barycentric coordinates reconstruct the reported point', () => {
        check(fc.tuple(latticeRay, latticeTri), ([r, t]) => {
            const res = fi.find(r, t);
            if (!res.intersect) {
                return;
            }
            const [b0, b1, b2] = res.triangleBary;
            expect(Number.isFinite(res.parameter)).toBe(true);
            expect(b0 + b1 + b2).toBeCloseTo(1, 12);
            expect(b1).toBeGreaterThanOrEqual(0);
            expect(b2).toBeGreaterThanOrEqual(0);
            expect(b0).toBeGreaterThanOrEqual(-1e-12);
            expectVectorClose(res.point,
                add(r.origin, mul(res.parameter, r.direction)), 0, 0);
            const fromBary = add(add(mul(b0, t.v[0]), mul(b1, t.v[1])),
                mul(b2, t.v[2]));
            const scale = 1 + length(t.v[0]) + length(t.v[1]) + length(t.v[2])
                + Math.abs(res.parameter);
            expect(length(sub(res.point, fromBary)))
                .toBeLessThanOrEqual(1e-11 * scale);
        });
    });

    it('a ray fired at a sampled interior point hits, and its reverse misses', () => {
        const rnd = makeRandom(0x71ab04);
        check(fc.tuple(latticeTri, unitVector(3)), ([t, d]) => {
            const n = cross(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
            if (length(n) < 1) {
                return;    // degenerate lattice triangle
            }
            normalize(n);
            if (Math.abs(dot(d, n)) < 1e-2) {
                return;
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
            const origin = sub(target, mul(4, d));
            const hit = fi.find(Ray.fromOriginDirection(origin, d), t);
            expect(hit.intersect).toBe(true);
            expect(hit.parameter).toBeGreaterThan(0);
            const scale = 1 + length(target);
            expect(length(sub(hit.point, target)))
                .toBeLessThanOrEqual(1e-9 * scale);
            // Firing away from the triangle must miss.
            expect(fi.find(Ray.fromOriginDirection(origin, mul(-1, d)), t)
                .intersect).toBe(false);
        });
    });

    it('a ray whose origin is on the triangle hits at parameter zero', () => {
        const origin = v3(0.25, 0.25, 0);
        for (const d of [v3(0, 0, 1), v3(0, 0, -1), v3(1, 1, 1)]) {
            const r = fi.find(makeRay(origin, d), unitTriangle);
            expect(r.intersect).toBe(true);
            expect(r.parameter).toBeCloseTo(0, 12);
            expectVectorClose(r.point, origin, 1e-12, 1e-12);
        }
    });

    it('a ray in the plane of the triangle reports no intersection', () => {
        const r = makeRay(v3(-5, 0.25, 0), v3(1, 0, 0));
        expect(ti.test(r, unitTriangle).intersect).toBe(false);
        const res = fi.find(r, unitTriangle);
        expect(res.intersect).toBe(false);
        expect(res.parameter).toBe(0);
        expect(res.triangleBary).toEqual([0, 0, 0]);
    });
});
