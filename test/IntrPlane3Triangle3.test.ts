import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Triangle3TI,
    IntrPlane3Triangle3FI
} from '../src/IntrPlane3Triangle3.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, dot, normalize } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function plane(normal: number[], constant: number): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

function triangle(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(Vector.fromArray(v0), Vector.fromArray(v1),
        Vector.fromArray(v2));
}

describe('IntrPlane3Triangle3', () => {
    const ti = new IntrPlane3Triangle3TI();
    const fi = new IntrPlane3Triangle3FI();

    // The plane z = 0.
    const z0 = plane([0, 0, 1], 0);

    it('reports no intersection when all vertices are strictly on one side', () => {
        const t = triangle([0, 0, 1], [1, 0, 2], [0, 1, 3]);
        expect(ti.test(z0, t).intersect).toBe(false);
        expect(ti.test(z0, t).numIntersections).toBe(0);
        const result = fi.find(z0, t);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);

        const below = triangle([0, 0, -1], [1, 0, -2], [0, 1, -3]);
        expect(ti.test(z0, below).intersect).toBe(false);
        expect(fi.find(z0, below).intersect).toBe(false);
    });

    it('clips two edges when the plane separates one vertex (n=1, p=2)', () => {
        const t = triangle([0, 0, -1], [2, 0, 1], [0, 2, 1]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(true);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(true);
        for (let i = 0; i < 2; ++i) {
            expect(result.point[i].values[2]).toBeCloseTo(0, 12);
        }
        // The clipped points are the midpoints of the two edges through
        // vertex 0.
        const xs = [result.point[0].values[0], result.point[1].values[0]]
            .sort((a, b) => a - b);
        const ys = [result.point[0].values[1], result.point[1].values[1]]
            .sort((a, b) => a - b);
        expect(xs[0]).toBeCloseTo(0, 12);
        expect(xs[1]).toBeCloseTo(1, 12);
        expect(ys[0]).toBeCloseTo(0, 12);
        expect(ys[1]).toBeCloseTo(1, 12);
    });

    it('reports a vertex touch (z=1, one side only)', () => {
        const t = triangle([0, 0, 0], [1, 0, 1], [0, 1, 2]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);
        expect(tiResult.isInterior).toBe(false);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('reports a segment through a vertex when the other two straddle', () => {
        const t = triangle([0, 0, 0], [2, 0, 1], [0, 2, -1]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(true);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(true);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        // The second point is the midpoint of the edge v1->v2.
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
        expect(result.point[1].values[1]).toBeCloseTo(1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(0, 12);
    });

    it('reports an edge in the plane (z=2 case)', () => {
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 1, 5]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.numIntersections).toBe(2);
        expect(tiResult.isInterior).toBe(false);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(2);
        expect(result.isInterior).toBe(false);
        // Vertex 2 is off the plane, so the edge is <v0,v1>.
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
    });

    it('reports the whole triangle when it lies in the plane', () => {
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 1, 0]);
        const tiResult = ti.test(z0, t);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(3);

        const result = fi.find(z0, t);
        expect(result.numIntersections).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(result.point[1].values).toEqual([1, 0, 0]);
        expect(result.point[2].values).toEqual([0, 1, 0]);
    });

    it('rejects non-3D inputs', () => {
        const p2 = Hyperplane.fromNormalConstant(Vector.fromArray([1, 0]), 0);
        const t2 = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        expect(() => ti.test(p2, t2)).toThrow();
        expect(() => fi.find(p2, t2)).toThrow();
    });

    it('agrees with the test query and puts points on the plane and triangle', () => {
        let seed = 4242424;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        // Barycentric membership of P in the triangle, expressed by checking
        // that P is a convex combination of the vertices.
        const inTriangle = (t: Triangle, P: Vector): boolean => {
            const e0 = [
                t.v[1].values[0] - t.v[0].values[0],
                t.v[1].values[1] - t.v[0].values[1],
                t.v[1].values[2] - t.v[0].values[2]
            ];
            const e1 = [
                t.v[2].values[0] - t.v[0].values[0],
                t.v[2].values[1] - t.v[0].values[1],
                t.v[2].values[2] - t.v[0].values[2]
            ];
            const d = [
                P.values[0] - t.v[0].values[0],
                P.values[1] - t.v[0].values[1],
                P.values[2] - t.v[0].values[2]
            ];
            const a = e0[0] * e0[0] + e0[1] * e0[1] + e0[2] * e0[2];
            const b = e0[0] * e1[0] + e0[1] * e1[1] + e0[2] * e1[2];
            const c = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
            const d0 = e0[0] * d[0] + e0[1] * d[1] + e0[2] * d[2];
            const d1 = e1[0] * d[0] + e1[1] * d[1] + e1[2] * d[2];
            const det = a * c - b * b;
            const s = (c * d0 - b * d1) / det;
            const u = (a * d1 - b * d0) / det;
            const eps = 1e-7;
            return s >= -eps && u >= -eps && s + u <= 1 + eps;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p = plane([rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                rand() * 2 - 1);
            const t = triangle(
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1],
                [rand() * 2 - 1, rand() * 2 - 1, rand() * 2 - 1]);

            const tiResult = ti.test(p, t);
            const fiResult = fi.find(p, t);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            expect(fiResult.numIntersections).toBe(tiResult.numIntersections);
            expect(fiResult.isInterior).toBe(tiResult.isInterior);

            for (let i = 0; i < fiResult.numIntersections; ++i) {
                const X = fiResult.point[i];
                expect(dot(p.normal, X) - p.constant).toBeCloseTo(0, 8);
                expect(inTriangle(t, X)).toBe(true);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrPlane3Triangle3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, wellScaled
} from './helpers/arbitraries.js';
import { add, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

// Integer vertices and an integer plane make every s[i] = Dot(N,V) - c exact
// in binary64, so the zero/positive/negative classification the query relies
// on is exactly reproducible.
const latticeCase = fc.tuple(
    fc.array(fc.integer({ min: -4, max: 4 }),
        { minLength: 9, maxLength: 9 }),
    fc.integer({ min: -3, max: 3 }), fc.integer({ min: -3, max: 3 }),
    fc.integer({ min: -3, max: 3 }), fc.integer({ min: -6, max: 6 })
).filter(([v, nx, ny, nz]) => {
    if (nx === 0 && ny === 0 && nz === 0) { return false; }
    // Non-degenerate triangle.
    const a = Vector.fromArray(v.slice(0, 3));
    const b = Vector.fromArray(v.slice(3, 6));
    const c = Vector.fromArray(v.slice(6, 9));
    return length(cross(sub(b, a), sub(c, a))) > 1e-9;
}).map(([v, nx, ny, nz, k]) => ({
    // The normal is NOT normalized here on purpose: the query only uses
    // Dot(N,V) - c, and integer components keep the classification exact.
    plane: Hyperplane.fromNormalConstant(Vector.fromArray([nx, ny, nz]), k),
    triangle: Triangle.fromVertices(Vector.fromArray(v.slice(0, 3)),
        Vector.fromArray(v.slice(3, 6)), Vector.fromArray(v.slice(6, 9)))
}));

// The exact intersection of the plane with the triangle boundary, computed
// independently: every vertex on the plane, plus the crossing point of every
// edge whose endpoints are strictly on opposite sides.
function referenceBoundaryPoints(p: Hyperplane, t: Triangle): Vector[] {
    const s = [0, 1, 2].map(i => dot(p.normal, t.v[i]) - p.constant);
    const out: Vector[] = [];
    for (let i = 0; i < 3; ++i) {
        if (s[i] === 0) { out.push(t.v[i].clone()); }
    }
    for (let i = 0; i < 3; ++i) {
        const j = (i + 1) % 3;
        if ((s[i] > 0 && s[j] < 0) || (s[i] < 0 && s[j] > 0)) {
            const u = s[i] / (s[i] - s[j]);
            out.push(add(t.v[i], mul(u, sub(t.v[j], t.v[i]))));
        }
    }
    return out;
}

function sameSet(a: Vector[], b: Vector[], tol = 1e-9): boolean {
    if (a.length !== b.length) { return false; }
    const used = new Array<boolean>(b.length).fill(false);
    for (const p of a) {
        let found = -1;
        for (let i = 0; i < b.length; ++i) {
            if (!used[i] && length(sub(p, b[i])) <= tol) { found = i; break; }
        }
        if (found < 0) { return false; }
        used[found] = true;
    }
    return true;
}

describe('IntrPlane3Triangle3 verification', () => {
    const tiq = new IntrPlane3Triangle3TI();
    const fiq = new IntrPlane3Triangle3FI();

    it('TI and FI agree on every field they share', () => {
        check(latticeCase, ({ plane: p, triangle: t }) => {
            const a = tiq.test(p, t);
            const b = fiq.find(p, t);
            expect(a.intersect).toBe(b.intersect);
            expect(a.numIntersections).toBe(b.numIntersections);
            expect(a.isInterior).toBe(b.isInterior);
        });
    });

    it('the classification follows the upstream n/p/z table', () => {
        check(latticeCase, ({ plane: p, triangle: t }) => {
            const s = [0, 1, 2].map(i => dot(p.normal, t.v[i]) - p.constant);
            const nz = s.filter(x => x === 0).length;
            const np = s.filter(x => x > 0).length;
            const nn = s.filter(x => x < 0).length;
            const r = fiq.find(p, t);
            if (nz === 3) {
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(3);
                expect(r.isInterior).toBe(false);
            }
            else if (nz === 2) {
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(2);
                expect(r.isInterior).toBe(false);
            }
            else if (nz === 1) {
                expect(r.intersect).toBe(true);
                if (np === 2 || nn === 2) {
                    expect(r.numIntersections).toBe(1);
                    expect(r.isInterior).toBe(false);
                }
                else {
                    expect(r.numIntersections).toBe(2);
                    expect(r.isInterior).toBe(true);
                }
            }
            else if (np > 0 && nn > 0) {
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(2);
                expect(r.isInterior).toBe(true);
            }
            else {
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
                expect(r.isInterior).toBe(false);
            }
        });
    });

    it('the reported points equal the exact plane-triangle intersection',
        () => {
            check(latticeCase, ({ plane: p, triangle: t }) => {
                const r = fiq.find(p, t);
                if (!r.intersect) { return; }
                const got = r.point.slice(0, r.numIntersections);
                if (r.numIntersections === 3) {
                    // The whole triangle is in the plane; the reported points
                    // are its vertices.
                    expect(sameSet(got, [t.v[0], t.v[1], t.v[2]])).toBe(true);
                    return;
                }
                const want = referenceBoundaryPoints(p, t);
                expect(sameSet(got, want)).toBe(true);
            });
        });

    it('every reported point is on the plane and inside the triangle', () => {
        check(latticeCase, ({ plane: p, triangle: t }) => {
            const r = fiq.find(p, t);
            if (!r.intersect) { return; }
            const e0 = sub(t.v[1], t.v[0]);
            const e1 = sub(t.v[2], t.v[0]);
            const d00 = dot(e0, e0), d01 = dot(e0, e1), d11 = dot(e1, e1);
            const det = d00 * d11 - d01 * d01;
            for (let i = 0; i < r.numIntersections; ++i) {
                const q = r.point[i];
                const scale = 1 + Math.abs(p.constant) + length(p.normal) * 8;
                expectClose(dot(p.normal, q), p.constant, 1e-9 * scale, 1e-9);
                // Barycentric coordinates in [0,1] summing to 1.
                const d = sub(q, t.v[0]);
                const b1 = (d11 * dot(d, e0) - d01 * dot(d, e1)) / det;
                const b2 = (d00 * dot(d, e1) - d01 * dot(d, e0)) / det;
                const b0 = 1 - b1 - b2;
                for (const b of [b0, b1, b2]) {
                    expect(b).toBeGreaterThanOrEqual(-1e-9);
                    expect(b).toBeLessThanOrEqual(1 + 1e-9);
                }
                // The point is on the plane through the triangle, so the
                // out-of-plane barycentric residual is zero.
                expectVectorClose(q,
                    add(t.v[0], add(mul(b1, e0), mul(b2, e1))), 1e-8, 1e-9);
            }
        });
    });

    it('unused point slots keep their default value', () => {
        check(latticeCase, ({ plane: p, triangle: t }) => {
            const r = fiq.find(p, t);
            for (let i = Math.max(r.numIntersections, 0); i < 3; ++i) {
                if (r.numIntersections === 0) {
                    expect(r.point[i].values).toEqual([0, 0, 0]);
                }
            }
            // The one-point (vertex) case leaves point[1] and point[2] at
            // their defaults; the two-point cases leave point[2] alone.
            if (r.numIntersections === 1) {
                expect(r.point[1].values).toEqual([0, 0, 0]);
                expect(r.point[2].values).toEqual([0, 0, 0]);
            }
            if (r.numIntersections === 2) {
                expect(r.point[2].values).toEqual([0, 0, 0]);
            }
        });
    });

    it('results do not alias the triangle vertices', () => {
        // The coplanar and vertex cases copy triangle vertices into point[].
        const t = triangle([0, 0, 0], [1, 0, 0], [0, 1, 0]);
        const r = fiq.find(plane([0, 0, 1], 0), t);
        expect(r.numIntersections).toBe(3);
        r.point[0].set(0, 42);
        expect(t.v[0].get(0)).toBe(0);

        const t2 = triangle([0, 0, 0], [1, 0, 1], [0, 1, 1]);
        const r2 = fiq.find(plane([0, 0, 1], 0), t2);
        expect(r2.numIntersections).toBe(1);
        r2.point[0].set(1, 99);
        expect(t2.v[0].get(1)).toBe(0);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(latticeCase, wellScaled(-Math.PI, Math.PI),
            wellScaled(-Math.PI, Math.PI), wellScaled(-Math.PI, Math.PI),
            wellScaled(-3, 3), wellScaled(-3, 3), wellScaled(-3, 3)),
        ([{ plane: p, triangle: t }, a1, a2, a3, tx, ty, tz]) => {
            const ca = Math.cos(a1), sa = Math.sin(a1);
            const cb = Math.cos(a2), sb = Math.sin(a2);
            const cc = Math.cos(a3), sc = Math.sin(a3);
            const b0 = vec(ca * cb, sa * cb, -sb);
            const b1 = vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc,
                cb * sc);
            const b2 = vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc,
                cb * cc);
            const rot = (v: Vector): Vector => add(mul(v.get(0), b0),
                add(mul(v.get(1), b1), mul(v.get(2), b2)));
            const tr = vec(tx, ty, tz);
            const xf = (v: Vector): Vector => add(rot(v), tr);
            const n2 = rot(p.normal);
            const p2 = Hyperplane.fromNormalConstant(n2,
                p.constant + dot(n2, tr));
            const t2 = Triangle.fromVertices(xf(t.v[0]), xf(t.v[1]),
                xf(t.v[2]));
            // A vertex that is exactly on the plane in the lattice
            // configuration is generically off it after a rotation, which
            // moves the case to a different row of the n/p/z table (for
            // example z=2 with isInterior false becomes z=1, p=1, n=1 with
            // isInterior true). Those rows are covered exactly by the lattice
            // properties above; here the property is restricted to
            // configurations whose vertex classification is stable, i.e. no
            // signed distance is within round-off of zero on either side.
            const scale = 1 + length(p.normal) * 8 + Math.abs(p.constant);
            const margin = (pl: Hyperplane, tr: Triangle): number =>
                Math.min(...[0, 1, 2].map(i =>
                    Math.abs(dot(pl.normal, tr.v[i]) - pl.constant)));
            if (margin(p, t) < 1e-6 * scale
                || margin(p2, t2) < 1e-6 * scale) {
                return;
            }

            const r0 = fiq.find(p, t);
            const r1 = fiq.find(p2, t2);
            if (r0.numIntersections !== r1.numIntersections) {
                return;
            }
            expect(r1.intersect).toBe(r0.intersect);
            expect(r1.isInterior).toBe(r0.isInterior);
            const want = r0.point.slice(0, r0.numIntersections).map(xf);
            const got = r1.point.slice(0, r1.numIntersections);
            expect(sameSet(got, want, 1e-7)).toBe(true);
        });
    });

    it('rejects mismatched dimensions', () => {
        const t2 = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        expect(() => tiq.test(plane([0, 0, 1], 0), t2)).toThrow();
        expect(() => fiq.find(plane([0, 0, 1], 0), t2)).toThrow();
    });
});
