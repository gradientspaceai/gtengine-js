import { describe, it, expect } from 'vitest';
import { Line } from '../src/Line.js';
import { Ray } from '../src/Ray.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, mul, normalize } from '../src/Vector.js';
import { IntrLine2Triangle2FI } from '../src/IntrLine2Triangle2.js';
import {
    IntrRay2Triangle2TI,
    IntrRay2Triangle2FI
} from '../src/IntrRay2Triangle2.js';
import { intrRay2Triangle2DoQuery } from '../src/IntrRay2Triangle2.js';
import { defaultIntrLine2Triangle2FIResult }
    from '../src/IntrLine2Triangle2.js';

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

// The triangle (0,0), (4,0), (0,4).
const tri = Triangle.fromVertices(vec([0, 0]), vec([4, 0]), vec([0, 4]));

function insideTriangle(t: Triangle, x: Vector): boolean {
    // Barycentric sign test with a small tolerance.
    const [v0, v1, v2] = t.v;
    const d = (v1.values[1] - v2.values[1]) * (v0.values[0] - v2.values[0])
        + (v2.values[0] - v1.values[0]) * (v0.values[1] - v2.values[1]);
    const b0 = ((v1.values[1] - v2.values[1]) * (x.values[0] - v2.values[0])
        + (v2.values[0] - v1.values[0]) * (x.values[1] - v2.values[1])) / d;
    const b1 = ((v2.values[1] - v0.values[1]) * (x.values[0] - v2.values[0])
        + (v0.values[0] - v2.values[0]) * (x.values[1] - v2.values[1])) / d;
    const b2 = 1 - b0 - b1;
    const eps = 1e-12;
    return b0 >= -eps && b1 >= -eps && b2 >= -eps;
}

describe('IntrRay2Triangle2', () => {
    const ti = new IntrRay2Triangle2TI();
    const fi = new IntrRay2Triangle2FI();

    it('finds the chord of a ray that crosses the triangle', () => {
        // The ray y = 1, x increasing from -3, crosses the triangle from
        // x = 0 to x = 3.
        const r = ray([-3, 1], [1, 0]);
        expect(ti.test(r, tri).intersect).toBe(true);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(6, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
    });

    it('clips the near end when the ray origin is inside the triangle', () => {
        const r = ray([1, 1], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(2, 12);
    });

    it('reports no intersection when the ray points away', () => {
        const r = ray([-3, 1], [-1, 0]);
        expect(ti.test(r, tri).intersect).toBe(false);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('reports a single point when the ray touches only a vertex', () => {
        // The line y = 4 supports the triangle at the vertex (0,4); the ray
        // reaches it at t = 3.
        const r = ray([-3, 4], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(4, 12);
    });

    it('handles a ray collinear with an edge, clipped at the origin', () => {
        // The ray starts at the vertex (4,0) and runs along the supporting
        // line of the edge from (0,0) to (4,0), leaving the triangle at once.
        const r = ray([4, 0], [1, 0]);
        const result = fi.find(r, tri);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.parameter[0]).toBeCloseTo(0, 12);
        expect(result.parameter[1]).toBeCloseTo(0, 12);
    });

    it('misses a triangle that lies behind and to the side', () => {
        const r = ray([-1, -1], [0, -1]);
        expect(ti.test(r, tri).intersect).toBe(false);
        expect(fi.find(r, tri).intersect).toBe(false);
    });

    it('is the line query clipped to t >= 0', () => {
        const rand = makeRandom(2718281);
        const lineFI = new IntrLine2Triangle2FI();
        for (let trial = 0; trial < 400; ++trial) {
            const r = ray([10 * rand() - 5, 10 * rand() - 5],
                [2 * rand() - 1, 2 * rand() - 1]);
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lineResult = lineFI.find(l, tri);
            const rayResult = fi.find(r, tri);
            expect(ti.test(r, tri).intersect).toBe(rayResult.intersect);

            if (!lineResult.intersect || lineResult.parameter[1] < 0) {
                expect(rayResult.intersect).toBe(false);
            }
            else {
                expect(rayResult.intersect).toBe(true);
                expect(rayResult.parameter[0]).toBeCloseTo(
                    Math.max(lineResult.parameter[0], 0), 12);
                expect(rayResult.parameter[1]).toBeCloseTo(
                    lineResult.parameter[1], 12);
            }
        }
    });

    it('agrees with dense sampling along the ray', () => {
        const rand = makeRandom(161803);
        for (let trial = 0; trial < 80; ++trial) {
            const r = ray([8 * rand() - 4, 8 * rand() - 4],
                [2 * rand() - 1, 2 * rand() - 1]);
            const result = fi.find(r, tri);

            let tLo = Number.POSITIVE_INFINITY;
            let tHi = Number.NEGATIVE_INFINITY;
            const n = 20000;
            for (let k = 0; k <= n; ++k) {
                const t = (16 * k) / n;
                const x = add(r.origin, mul(t, r.direction));
                if (insideTriangle(tri, x)) {
                    if (t < tLo) { tLo = t; }
                    if (t > tHi) { tHi = t; }
                }
            }

            if (tLo <= tHi) {
                expect(result.intersect).toBe(true);
                expect(result.parameter[0]).toBeLessThanOrEqual(tLo + 1e-9);
                expect(result.parameter[1]).toBeGreaterThanOrEqual(tHi - 1e-9);
                expect(tLo - result.parameter[0]).toBeLessThan(3e-3);
                expect(result.parameter[1] - tHi).toBeLessThan(3e-3);
            }
        }
    });
});

describe('intrRay2Triangle2DoQuery', () => {
    const tri = Triangle.fromVertices(vec([0, 0]), vec([4, 0]), vec([0, 4]));

    it('matches the class query but does not compute points', () => {
        const r = ray([-1, 1], [1, 0]);
        const result = defaultIntrLine2Triangle2FIResult();
        intrRay2Triangle2DoQuery(r.origin, r.direction, tri, result);
        const expected = new IntrRay2Triangle2FI().find(r, tri);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter[0]).toBeCloseTo(expected.parameter[0], 12);
        expect(result.parameter[1]).toBeCloseTo(expected.parameter[1], 12);
        // DoQuery leaves 'point' at its default value.
        expect(result.point[0].values).toEqual([0, 0]);
        expect(result.point[1].values).toEqual([0, 0]);
    });

    it('rejects a ray whose supporting line hits behind the origin', () => {
        const r = ray([-1, 1], [-1, 0]);
        const result = defaultIntrLine2Triangle2FIResult();
        intrRay2Triangle2DoQuery(r.origin, r.direction, tri, result);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, unitVector, wellScaledVector, expectClose, expectVectorClose,
    seededRandom
} from './helpers/arbitraries.js';
import { sub, dot } from '../src/Vector.js';

// A non-degenerate triangle built from well-scaled vertices (the shared
// 'triangle' generator draws from fc.double, which emits subnormal
// coordinates whose edge products underflow).
const triangle2Arb = fc.tuple(wellScaledVector(2, -5, 5),
    wellScaledVector(2, -5, 5), wellScaledVector(2, -5, 5))
    .filter(([a, b, c]) => {
        const e0 = sub(b, a), e1 = sub(c, a);
        const d00 = dot(e0, e0), d11 = dot(e1, e1), d01 = dot(e0, e1);
        return d00 * d11 - d01 * d01 > 1;
    })
    .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

const rayTriangle2 = fc.tuple(wellScaledVector(2, -8, 8), unitVector(2),
    triangle2Arb)
    .map(([o, d, t]) => ({ ray: Ray.fromOriginDirection(o, d), triangle: t }));

// Barycentric coordinates of p with respect to the triangle.
function barycentric(t: Triangle, p: Vector): number[] {
    const e0 = sub(t.v[1], t.v[0]);
    const e1 = sub(t.v[2], t.v[0]);
    const q = sub(p, t.v[0]);
    const det = e0.values[0] * e1.values[1] - e0.values[1] * e1.values[0];
    const b1 = (q.values[0] * e1.values[1] - q.values[1] * e1.values[0]) / det;
    const b2 = (e0.values[0] * q.values[1] - e0.values[1] * q.values[0]) / det;
    return [1 - b1 - b2, b1, b2];
}

describe('IntrRay2Triangle2 verification', () => {
    const tiq = new IntrRay2Triangle2TI();
    const fiq = new IntrRay2Triangle2FI();
    const lfi = new IntrLine2Triangle2FI();

    it('TI and FI agree on intersect', () => {
        check(rayTriangle2, ({ ray: r, triangle: t }) => {
            expect(tiq.test(r, t).intersect).toBe(fiq.find(r, t).intersect);
        });
    });

    it('a ray hit is the line hit clipped to nonnegative parameters', () => {
        check(rayTriangle2, ({ ray: r, triangle: t }) => {
            const l = Line.fromOriginDirection(r.origin, r.direction);
            const lf = lfi.find(l, t);
            const f = fiq.find(r, t);
            if (!lf.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const lo = Math.max(lf.parameter[0], 0);
            const hi = lf.parameter[1];
            expect(f.intersect).toBe(hi >= 0);
            if (!f.intersect) {
                // The reset restores every default field.
                expect(f.numIntersections).toBe(0);
                expect(f.parameter).toEqual([0, 0]);
                return;
            }
            expectClose(f.parameter[0], lo, 0, 0);
            expectClose(f.parameter[1], hi, 0, 0);
            expect(f.numIntersections).toBe(lo < hi ? 2 : 1);
        });
    });

    it('the reported points are on the ray and inside the triangle', () => {
        check(rayTriangle2, ({ ray: r, triangle: t }) => {
            const f = fiq.find(r, t);
            if (!f.intersect) {
                return;
            }
            expect(f.parameter[0]).toBeGreaterThanOrEqual(0);
            expect(f.parameter[0]).toBeLessThanOrEqual(f.parameter[1]);
            for (let i = 0; i < f.numIntersections; ++i) {
                expectVectorClose(f.point[i],
                    add(r.origin, mul(f.parameter[i], r.direction)), 0, 0);
                const b = barycentric(t, f.point[i]);
                expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
                for (const bi of b) {
                    expect(bi).toBeGreaterThanOrEqual(-1e-9);
                }
            }
        });
    });

    it('a ray fired at an interior point hits and its reverse misses', () => {
        check(fc.tuple(triangle2Arb, unitVector(2),
            fc.double({ min: 0.1, max: 0.8, noNaN: true }),
            fc.double({ min: 0.1, max: 0.8, noNaN: true }),
            fc.double({ min: 0.5, max: 6, noNaN: true })),
            ([t, d, b1, b2, back]) => {
                if (b1 + b2 > 0.9) {
                    return;
                }
                const target = add(t.v[0],
                    add(mul(b1, sub(t.v[1], t.v[0])),
                        mul(b2, sub(t.v[2], t.v[0]))));
                const r = Ray.fromOriginDirection(sub(target, mul(back, d)), d);
                expect(fiq.find(r, t).intersect).toBe(true);
                expect(tiq.test(r, t).intersect).toBe(true);
                // The reverse ray misses, but only when its origin is
                // outside the triangle (a long enough step back can land
                // inside on the far side).
                const ob = barycentric(t, r.origin);
                if (ob[0] < -1e-6 || ob[1] < -1e-6 || ob[2] < -1e-6) {
                    const rev = Ray.fromOriginDirection(r.origin, mul(-1, d));
                    expect(fiq.find(rev, t).intersect).toBe(false);
                }
            });
    });

    it('a ray whose origin is inside the triangle hits at parameter 0', () => {
        check(fc.tuple(triangle2Arb, unitVector(2),
            fc.double({ min: 0.1, max: 0.7, noNaN: true }),
            fc.double({ min: 0.1, max: 0.7, noNaN: true })),
            ([t, d, b1, b2]) => {
                if (b1 + b2 > 0.9) {
                    return;
                }
                const o = add(t.v[0], add(mul(b1, sub(t.v[1], t.v[0])),
                    mul(b2, sub(t.v[2], t.v[0]))));
                const r = Ray.fromOriginDirection(o, d);
                const f = fiq.find(r, t);
                expect(f.intersect).toBe(true);
                expectClose(f.parameter[0], 0, 1e-12, 0);
            });
    });

    it('a fine sweep of the ray agrees with the reported interval', () => {
        const rnd = seededRandom(0x7f3a91c5);
        for (let trial = 0; trial < 150; ++trial) {
            const t = Triangle.fromVertices(
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]),
                Vector.fromArray([rnd() * 6 - 3, rnd() * 6 - 3]));
            const e0 = sub(t.v[1], t.v[0]), e1 = sub(t.v[2], t.v[0]);
            const area2 = e0.values[0] * e1.values[1]
                - e0.values[1] * e1.values[0];
            if (Math.abs(area2) < 1) {
                continue;
            }
            const a = rnd() * 2 * Math.PI;
            const d = Vector.fromArray([Math.cos(a), Math.sin(a)]);
            const r = Ray.fromOriginDirection(
                Vector.fromArray([rnd() * 10 - 5, rnd() * 10 - 5]), d);
            const f = fiq.find(r, t);
            for (let k = 0; k <= 600; ++k) {
                const s = (14 * k) / 600;
                const p = add(r.origin, mul(s, d));
                const b = barycentric(t, p);
                if (b[0] > 1e-6 && b[1] > 1e-6 && b[2] > 1e-6) {
                    expect(f.intersect).toBe(true);
                    expect(s).toBeGreaterThanOrEqual(f.parameter[0] - 1e-9);
                    expect(s).toBeLessThanOrEqual(f.parameter[1] + 1e-9);
                }
            }
        }
    }, 30000);

    it('the exported DoQuery reproduces the class result', () => {
        check(rayTriangle2, ({ ray: r, triangle: t }) => {
            const res = defaultIntrLine2Triangle2FIResult();
            intrRay2Triangle2DoQuery(r.origin, r.direction, t, res);
            const f = fiq.find(r, t);
            expect(res.intersect).toBe(f.intersect);
            expect(res.numIntersections).toBe(f.numIntersections);
            for (let i = 0; i < res.numIntersections; ++i) {
                expect(res.parameter[i]).toBe(f.parameter[i]);
            }
            // The helper does not fill point[]; the class does.
            if (!res.intersect) {
                expect(res.point[0].values).toEqual([0, 0]);
            }
        });
    });
});
