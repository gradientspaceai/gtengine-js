import { describe, expect, it } from 'vitest';
import { DistLine3Triangle3 } from '../src/DistLine3Triangle3.js';
import { Line } from '../src/Line.js';
import { Triangle } from '../src/Triangle.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3Triangle3', () => {
    const query = new DistLine3Triangle3();
    // The triangle in the z = 0 plane with vertices (0,0,0), (1,0,0),
    // (0,1,0).
    const tri = Triangle.fromVertices(v(0, 0, 0), v(1, 0, 0), v(0, 1, 0));

    it('reports zero distance for a line piercing the triangle', () => {
        const result = query.compute(line([0.25, 0.25, -3], [0, 0, 1]), tri);
        expect(result.distance).toBe(0);
        expect(result.parameter).toBeCloseTo(3, 12);
        expect(result.barycentric[0]).toBeCloseTo(0.5, 12);
        expect(result.barycentric[1]).toBeCloseTo(0.25, 12);
        expect(result.barycentric[2]).toBeCloseTo(0.25, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to the plane of the triangle', () => {
        const result = query.compute(line([0, 0, 4], [1, 0, 0]), tri);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a line whose plane intersection is outside the triangle',
        () => {
            // The line pierces the plane at (2,2,0), outside the triangle.
            const result = query.compute(line([2, 2, -1], [0, 0, 1]), tri);
            // The closest triangle point is the midpoint region of the
            // hypotenuse, at (0.5,0.5,0); the distance is sqrt(2*1.5^2).
            expect(result.distance).toBeCloseTo(Math.sqrt(2) * 1.5, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(0.5, 10);
            expect(result.closest[1].values[1]).toBeCloseTo(0.5, 10);
        });

    it('reports the barycentric coordinates of an edge closest point', () => {
        const result = query.compute(line([2, 2, -1], [0, 0, 1]), tri);
        const sum = result.barycentric[0] + result.barycentric[1]
            + result.barycentric[2];
        expect(sum).toBeCloseTo(1, 10);
        for (const b of result.barycentric) {
            expect(b).toBeGreaterThanOrEqual(-1e-12);
        }
        // The closest point reconstructed from the barycentric coordinates.
        const q = add(mul(result.barycentric[0], tri.v[0]),
            add(mul(result.barycentric[1], tri.v[1]),
                mul(result.barycentric[2], tri.v[2])));
        for (let i = 0; i < 3; ++i) {
            expect(q.values[i]).toBeCloseTo(result.closest[1].values[i], 9);
        }
    });

    it('reports zero distance for a line lying in an edge of the triangle',
        () => {
            const result = query.compute(line([0, 0, 0], [1, 0, 0]), tri);
            expect(result.distance).toBeCloseTo(0, 12);
        });

    it('handles a degenerate triangle collapsed to a segment', () => {
        const degenerate = Triangle.fromVertices(v(0, 0, 0), v(1, 0, 0),
            v(2, 0, 0));
        const result = query.compute(line([0, 3, 0], [1, 0, 0]), degenerate);
        expect(result.distance).toBeCloseTo(3, 10);
    });

    it('agrees with a dense sampling of the triangle', () => {
        const rnd = makeRandom(24689);
        const t = Triangle.fromVertices(v(0.5, -1, 0.25), v(2, 0.5, -0.5),
            v(-1, 1.5, 1));

        for (let trial = 0; trial < 30; ++trial) {
            const origin = v(6 * rnd() - 3, 6 * rnd() - 3, 6 * rnd() - 3);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, t);

            // The barycentric coordinates are valid and reproduce closest[1].
            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            for (const bi of b) {
                expect(bi).toBeGreaterThanOrEqual(-1e-8);
            }
            const q = add(mul(b[0], t.v[0]),
                add(mul(b[1], t.v[1]), mul(b[2], t.v[2])));
            for (let i = 0; i < 3; ++i) {
                expect(q.values[i]).toBeCloseTo(result.closest[1].values[i],
                    7);
            }

            // The reported line point matches the reported parameter.
            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onLine.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
            }

            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled triangle point is closer to the line.
            const n = 60;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; i + j <= n; ++j) {
                    const b1 = i / n, b2 = j / n, b0 = 1 - b1 - b2;
                    const p = add(mul(b0, t.v[0]),
                        add(mul(b1, t.v[1]), mul(b2, t.v[2])));
                    const w = sub(p, ln.origin);
                    const s = dot(w, ln.direction) / dd;
                    const f = sub(w, mul(s, ln.direction));
                    best = Math.min(best, dot(f, f));
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine3Triangle3.h.
// ---------------------------------------------------------------------------

describe('DistLine3Triangle3 verification', () => {
    const query = new DistLine3Triangle3();

    // A triangle whose area is bounded away from zero, so the barycentric
    // solve in the port is well conditioned.
    const triArb = fc.tuple(wellScaledVector(3, -5, 5),
        wellScaledVector(3, -5, 5), wellScaledVector(3, -5, 5))
        .filter(([a, b, c]) => {
            const e0 = sub(b, a);
            const e1 = sub(c, a);
            return length(cross(e0, e1)) > 1;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    const lineArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    // Independent distance from a point to a solid triangle in 3D: project
    // into the plane, use the barycentric test, otherwise minimize over the
    // three edges.
    function pointSegmentDistance(p: Vector, p0: Vector, p1: Vector): number {
        const d = sub(p1, p0);
        const dd = dot(d, d);
        const t = dd > 0
            ? Math.min(Math.max(dot(sub(p, p0), d) / dd, 0), 1) : 0;
        return length(sub(p, add(p0, mul(t, d))));
    }

    function pointTriangleDistance(p: Vector, tri: Triangle): number {
        const e1 = sub(tri.v[1], tri.v[0]);
        const e2 = sub(tri.v[2], tri.v[0]);
        const n = cross(e1, e2);
        const nn = dot(n, n);
        const d = sub(p, tri.v[0]);
        const h = dot(n, d) / nn;                 // scaled signed height
        const q = sub(d, mul(h, n));              // in-plane part
        const e1e1 = dot(e1, e1), e1e2 = dot(e1, e2), e2e2 = dot(e2, e2);
        const det = e1e1 * e2e2 - e1e2 * e1e2;
        const b1 = (e2e2 * dot(e1, q) - e1e2 * dot(e2, q)) / det;
        const b2 = (e1e1 * dot(e2, q) - e1e2 * dot(e1, q)) / det;
        if (b1 >= 0 && b2 >= 0 && b1 + b2 <= 1) {
            return Math.abs(h) * Math.sqrt(nn);
        }
        return Math.min(
            pointSegmentDistance(p, tri.v[0], tri.v[1]),
            pointSegmentDistance(p, tri.v[1], tri.v[2]),
            pointSegmentDistance(p, tri.v[2], tri.v[0]));
    }

    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, triArb), ([ln, tri]) => {
            const r = query.compute(ln, tri);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-9, 1e-9);
            expectVectorClose(r.closest[0],
                add(ln.origin, mul(r.parameter, ln.direction)), 1e-9, 1e-9);
            // The barycentric coordinates are a convex combination.
            const b = r.barycentric;
            expectClose(b[0] + b[1] + b[2], 1, 1e-9, 1e-9);
            for (let i = 0; i < 3; ++i) {
                expect(b[i]).toBeGreaterThanOrEqual(-1e-12);
                expect(b[i]).toBeLessThanOrEqual(1 + 1e-12);
            }
            // closest[1] is the point with those barycentric coordinates.
            expectVectorClose(r.closest[1],
                add(add(mul(b[0], tri.v[0]), mul(b[1], tri.v[1])),
                    mul(b[2], tri.v[2])), 1e-8, 1e-8);
        });
    });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(lineArb, triArb), ([ln, tri]) => {
            const r = query.compute(ln, tri);
            const best = ternaryMin(
                t => pointTriangleDistance(
                    add(ln.origin, mul(t, ln.direction)), tri), -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });

    it('reports zero distance when the line crosses the triangle', () => {
        check(fc.tuple(triArb, finite(0.05, 0.9), finite(0.05, 0.9),
            unitVector(3)), ([tri, s, t, dir]) => {
            if (s + t > 0.95) { return; }
            const n = cross(sub(tri.v[1], tri.v[0]), sub(tri.v[2], tri.v[0]));
            normalize(n);
            // A line inside the plane of the triangle can cross it along a
            // whole chord, and the header says only one of the infinitely
            // many closest pairs is returned. Require a transverse direction
            // so that the intersection point is unique.
            if (Math.abs(dot(n, dir)) < 0.1) { return; }
            const target = add(tri.v[0],
                add(mul(s, sub(tri.v[1], tri.v[0])),
                    mul(t, sub(tri.v[2], tri.v[0]))));
            const ln = Line.fromOriginDirection(add(target, mul(4, dir)), dir);
            const r = query.compute(ln, tri);
            expectClose(r.distance, 0, 1e-9, 1e-9);
            expectVectorClose(r.closest[1], target, 1e-7, 1e-7);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, triArb, rotationFrame(3),
            wellScaledVector(3, -6, 6)), ([ln, tri, frame, shift]) => {
            const rot = (p: Vector): Vector =>
                add(add(mul(p.values[0], frame[0]), mul(p.values[1], frame[1])),
                    mul(p.values[2], frame[2]));
            const movedLine = Line.fromOriginDirection(
                add(shift, rot(ln.origin)), rot(ln.direction));
            const movedTri = Triangle.fromVertices(
                add(shift, rot(tri.v[0])), add(shift, rot(tri.v[1])),
                add(shift, rot(tri.v[2])));
            const r0 = query.compute(ln, tri);
            const r1 = query.compute(movedLine, movedTri);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
        });
    });

    it('is invariant to a cyclic relabelling of the vertices', () => {
        check(fc.tuple(lineArb, triArb), ([ln, tri]) => {
            const rotated = Triangle.fromVertices(tri.v[1], tri.v[2],
                tri.v[0]);
            const r0 = query.compute(ln, tri);
            const r1 = query.compute(ln, rotated);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            // Only the distance is relabelling invariant: when the line lies
            // in the plane of the triangle there are infinitely many closest
            // pairs and the edge loop reports whichever edge it visits first,
            // so the barycentric coordinates may legitimately differ. Both
            // results must still be points of the triangle at the reported
            // distance, which the first property checks.
            expectClose(length(sub(r1.closest[0], r1.closest[1])),
                r1.distance, 1e-9, 1e-9);
        });
    });

    it('handles a line parallel to the plane of the triangle', () => {
        check(fc.tuple(triArb, finite(0.1, 5), finite(-Math.PI, Math.PI),
            finite(-3, 3)), ([tri, height, angle, along]) => {
            const e1 = sub(tri.v[1], tri.v[0]);
            const n = cross(e1, sub(tri.v[2], tri.v[0]));
            normalize(n);
            const u = e1.clone();
            normalize(u);
            const w = cross(n, u);
            const dir = add(mul(Math.cos(angle), u), mul(Math.sin(angle), w));
            const perp = add(mul(-Math.sin(angle), u), mul(Math.cos(angle), w));
            const origin = add(tri.v[0],
                add(mul(height, n), mul(along, perp)));
            const ln = Line.fromOriginDirection(origin, dir);
            const r = query.compute(ln, tri);
            const best = ternaryMin(
                t => pointTriangleDistance(
                    add(ln.origin, mul(t, ln.direction)), tri), -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });
});
