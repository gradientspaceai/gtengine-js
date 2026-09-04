import { describe, expect, it } from 'vitest';
import { DistPoint3Frustum3 } from '../src/DistPoint3Frustum3.js';
import { Frustum3 } from '../src/Frustum3.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The canonical frustum: origin at (0,0,0), D = +z, U = +y, R = +x.
function canonical(dMin: number, dMax: number, uBound: number,
    rBound: number): Frustum3 {
    return Frustum3.fromParameters(v(0, 0, 0), v(0, 0, 1), v(0, 1, 0),
        v(1, 0, 0), dMin, dMax, uBound, rBound);
}

describe('DistPoint3Frustum3', () => {
    const query = new DistPoint3Frustum3();

    it('reports zero distance for a point inside the frustum', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, 2), f);
        expect(result.distance).toBe(0);
        expect(result.closest[1].values[2]).toBeCloseTo(2, 12);
    });

    it('measures a point beyond the far plane', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, 7), f);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 12);
    });

    it('measures a point behind the near plane', () => {
        const f = canonical(1, 3, 1, 1);
        const result = query.compute(v(0, 0, -2), f);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('clamps to the far-plane rectangle corner', () => {
        const f = canonical(1, 3, 1, 1);
        // rmax = umax = 3 at the far plane z = 3.
        const result = query.compute(v(6, 7, 3), f);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(3, 12);
    });

    it('is symmetric under reflection in the R and U directions', () => {
        const f = canonical(1, 4, 0.75, 1.5);
        const rnd = makeRandom(555);
        for (let trial = 0; trial < 40; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = query.compute(p, f);
            const r1 = query.compute(
                v(-p.values[0], p.values[1], p.values[2]), f);
            const r2 = query.compute(
                v(p.values[0], -p.values[1], p.values[2]), f);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            expect(r2.distance).toBeCloseTo(r0.distance, 10);
            expect(r1.closest[1].values[0]).toBeCloseTo(
                -r0.closest[1].values[0], 10);
            expect(r2.closest[1].values[1]).toBeCloseTo(
                -r0.closest[1].values[1], 10);
        }
    });

    it('is invariant to a rigid motion of point and frustum', () => {
        const rnd = makeRandom(8080);
        const a = 0.9;
        const ca = Math.cos(a), sa = Math.sin(a);
        const rot = (p: Vector): Vector => v(
            ca * p.values[0] - sa * p.values[2],
            p.values[1],
            sa * p.values[0] + ca * p.values[2]);
        const shift = v(2, -3, 1);

        const f0 = canonical(1, 3, 1, 2);
        const f1 = Frustum3.fromParameters(add(shift, v(0, 0, 0)),
            rot(v(0, 0, 1)), rot(v(0, 1, 0)), rot(v(1, 0, 0)), 1, 3, 1, 2);

        for (let trial = 0; trial < 40; ++trial) {
            const p = v(10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5);
            const r0 = query.compute(p, f0);
            const r1 = query.compute(add(rot(p), shift), f1);
            expect(r1.distance).toBeCloseTo(r0.distance, 10);
            const expected = add(rot(r0.closest[1]), shift);
            for (let i = 0; i < 3; ++i) {
                expect(r1.closest[1].values[i]).toBeCloseTo(
                    expected.values[i], 9);
            }
        }
    });

    it('agrees with a dense sampling of the solid frustum', () => {
        const rnd = makeRandom(13579);
        const origin = v(0.5, -1, 0.25);
        const dVec = v(1, 1, 2);
        normalize(dVec);
        // A right-handed orthonormal frame containing dVec.
        const uVec = v(-1, 1, 0);
        normalize(uVec);
        const rVec = v(
            uVec.values[1] * dVec.values[2] - uVec.values[2] * dVec.values[1],
            uVec.values[2] * dVec.values[0] - uVec.values[0] * dVec.values[2],
            uVec.values[0] * dVec.values[1] - uVec.values[1] * dVec.values[0]);
        normalize(rVec);

        const dMin = 1, dMax = 2.5, uBound = 0.8, rBound = 1.2;
        const f = Frustum3.fromParameters(origin, dVec, uVec, rVec,
            dMin, dMax, uBound, rBound);

        for (let trial = 0; trial < 25; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, f);

            // The reported closest point is in the solid frustum.
            const delta = sub(result.closest[1], origin);
            const z = dot(delta, dVec);
            const x = dot(delta, rVec);
            const y = dot(delta, uVec);
            expect(z).toBeGreaterThanOrEqual(dMin - 1e-9);
            expect(z).toBeLessThanOrEqual(dMax + 1e-9);
            expect(Math.abs(x)).toBeLessThanOrEqual(
                rBound * z / dMin + 1e-9);
            expect(Math.abs(y)).toBeLessThanOrEqual(
                uBound * z / dMin + 1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 9);

            // No sampled solid-frustum point is closer.
            const n = 24;
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                const zz = dMin + (dMax - dMin) * i / n;
                const scale = zz / dMin;
                for (let j = 0; j <= n; ++j) {
                    const xx = rBound * scale * (2 * j / n - 1);
                    for (let k = 0; k <= n; ++k) {
                        const yy = uBound * scale * (2 * k / n - 1);
                        const q = add(origin, add(mul(zz, dVec),
                            add(mul(xx, rVec), mul(yy, uVec))));
                        const g = sub(p, q);
                        best = Math.min(best, dot(g, g));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistPoint3Frustum3.h. The reference computation
// treats the frustum as the convex polyhedron of its eight vertices, which
// is independent of the Voronoi-region case analysis under test.
// ---------------------------------------------------------------------------

describe('DistPoint3Frustum3 verification', () => {
    const query = new DistPoint3Frustum3();

    const frustumArb = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
        finite(0.2, 2), finite(0.2, 4), finite(0.2, 2), finite(0.2, 2))
        .map(([origin, frame, dMin, extra, uBound, rBound]) =>
            Frustum3.fromParameters(origin, frame[0], frame[1], frame[2],
                dMin, dMin + extra, uBound, rBound));

    // The six faces as cyclic vertex index loops. Vertices 0..3 are the near
    // quad and 4..7 the corresponding far quad (see Frustum3.computeVertices).
    const faces: readonly (readonly number[])[] = [
        [0, 1, 2, 3], [4, 5, 6, 7],
        [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]
    ];

    function pointSegmentDistance(p: Vector, a: Vector, b: Vector): number {
        const d = sub(b, a);
        const dd = dot(d, d);
        const t = dd > 0 ? Math.min(Math.max(dot(sub(p, a), d) / dd, 0), 1) : 0;
        return length(sub(p, add(a, mul(t, d))));
    }

    // Outward unit normal of a face, oriented away from the centroid.
    function faceNormal(vertex: Vector[], face: readonly number[],
        centroid: Vector): Vector {
        const a = vertex[face[0]];
        const n = cross(sub(vertex[face[1]], a), sub(vertex[face[2]], a));
        normalize(n);
        return dot(n, sub(centroid, a)) > 0 ? mul(-1, n) : n;
    }

    function pointFaceDistance(p: Vector, vertex: Vector[],
        face: readonly number[], n: Vector): number {
        const a = vertex[face[0]];
        const q = sub(p, mul(dot(n, sub(p, a)), n));
        // The vertex loops are cyclic but their winding relative to the
        // outward normal differs from face to face, so the inside test asks
        // that every edge cross product agree in sign rather than assuming a
        // particular orientation.
        let positive = 0;
        let negative = 0;
        for (let i = 0; i < face.length; ++i) {
            const v0 = vertex[face[i]];
            const v1 = vertex[face[(i + 1) % face.length]];
            const side = dot(cross(sub(v1, v0), sub(q, v0)), n);
            if (side > 0) { ++positive; }
            else if (side < 0) { ++negative; }
        }
        if (positive === 0 || negative === 0) {
            return length(sub(p, q));
        }
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i < face.length; ++i) {
            best = Math.min(best, pointSegmentDistance(p,
                vertex[face[i]], vertex[face[(i + 1) % face.length]]));
        }
        return best;
    }

    function polyhedronDistance(p: Vector, frustum: Frustum3): number {
        const vertex = frustum.computeVertices();
        let centroid = new Vector(3);
        for (const w of vertex) { centroid = add(centroid, w); }
        centroid = mul(1 / 8, centroid);
        let inside = true;
        const normals: Vector[] = [];
        for (const face of faces) {
            const n = faceNormal(vertex, face, centroid);
            normals.push(n);
            if (dot(n, sub(p, vertex[face[0]])) > 0) { inside = false; }
        }
        if (inside) { return 0; }
        let best = Number.POSITIVE_INFINITY;
        for (let i = 0; i < faces.length; ++i) {
            best = Math.min(best,
                pointFaceDistance(p, vertex, faces[i], normals[i]));
        }
        return best;
    }

    function insideFrustum(p: Vector, frustum: Frustum3,
        tol: number): boolean {
        const vertex = frustum.computeVertices();
        let centroid = new Vector(3);
        for (const w of vertex) { centroid = add(centroid, w); }
        centroid = mul(1 / 8, centroid);
        for (const face of faces) {
            const n = faceNormal(vertex, face, centroid);
            if (dot(n, sub(p, vertex[face[0]])) > tol) { return false; }
        }
        return true;
    }

    it('matches the convex-polyhedron distance', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), frustumArb),
            ([p, frustum]) => {
                const r = query.compute(p, frustum);
                expectClose(r.distance, polyhedronDistance(p, frustum),
                    1e-8, 1e-8);
            });
    });

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), frustumArb),
            ([p, frustum]) => {
                const r = query.compute(p, frustum);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                expectVectorClose(r.closest[0], p, 0, 0);
                expect(r.closest[0]).not.toBe(p);
                expectClose(length(sub(r.closest[0], r.closest[1])),
                    r.distance, 1e-8, 1e-8);
                expect(insideFrustum(r.closest[1], frustum, 1e-8)).toBe(true);
            });
    });

    it('reports zero distance for points inside the frustum', () => {
        check(fc.tuple(frustumArb, finite(0.02, 0.98), finite(-0.95, 0.95),
            finite(-0.95, 0.95)), ([frustum, dU, uU, rU]) => {
            const d = frustum.dMin + dU * (frustum.dMax - frustum.dMin);
            const scale = d / frustum.dMin;
            const p = add(frustum.origin,
                add(mul(d, frustum.dVector),
                    add(mul(uU * frustum.uBound * scale, frustum.uVector),
                        mul(rU * frustum.rBound * scale, frustum.rVector))));
            const r = query.compute(p, frustum);
            expectClose(r.distance, 0, 1e-9, 1e-9);
            expectVectorClose(r.closest[1], p, 1e-8, 1e-8);
        });
    });

    it('returns each frustum vertex as its own closest point', () => {
        check(frustumArb, frustum => {
            for (const w of frustum.computeVertices()) {
                const r = query.compute(w, frustum);
                expectClose(r.distance, 0, 1e-8, 1e-8);
            }
        }, 100);
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), frustumArb,
            rotationFrame(3), wellScaledVector(3, -6, 6)),
            ([p, frustum, frame, shift]) => {
                const rot = (q: Vector): Vector =>
                    add(add(mul(q.values[0], frame[0]),
                        mul(q.values[1], frame[1])),
                        mul(q.values[2], frame[2]));
                const moved = Frustum3.fromParameters(
                    add(shift, rot(frustum.origin)), rot(frustum.dVector),
                    rot(frustum.uVector), rot(frustum.rVector),
                    frustum.dMin, frustum.dMax, frustum.uBound,
                    frustum.rBound);
                const r0 = query.compute(p, frustum);
                const r1 = query.compute(add(shift, rot(p)), moved);
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            });
    });

    it('is symmetric under reflection of the R and U coordinates', () => {
        // The algorithm folds the test point into the octant with
        // nonnegative R and U coordinates, so the reflected query must give
        // the same distance and the reflected closest point.
        check(fc.tuple(wellScaledVector(3, -8, 8), frustumArb),
            ([p, frustum]) => {
                const diff = sub(p, frustum.origin);
                const cr = dot(diff, frustum.rVector);
                const cu = dot(diff, frustum.uVector);
                const cd = dot(diff, frustum.dVector);
                const reflected = add(frustum.origin,
                    add(mul(-cr, frustum.rVector),
                        add(mul(-cu, frustum.uVector),
                            mul(cd, frustum.dVector))));
                const r0 = query.compute(p, frustum);
                const r1 = query.compute(reflected, frustum);
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            });
    });
});

// ---------------------------------------------------------------------------
// Regression for the unclamped far-edge coordinate described in
// src/DistPoint3Frustum3.ts.
// ---------------------------------------------------------------------------

describe('DistPoint3Frustum3 far-edge clamping', () => {
    const query = new DistPoint3Frustum3();

    // D is world x, U is world y, R is world z, so a point at frustum
    // coordinates (r,u,d) is the world point (d,u,r).
    const frustum = Frustum3.fromParameters(v(0, 0, 0), v(1, 0, 0),
        v(0, 1, 0), v(0, 0, 1), 1, 2, 1, 1);

    function nearestVertexDistance(p: Vector): number {
        let best = Number.POSITIVE_INFINITY;
        for (const w of frustum.computeVertices()) {
            best = Math.min(best, length(sub(p, w)));
        }
        return best;
    }

    it('clamps the LF-edge u coordinate to the far half-extent', () => {
        // (r,u,d) = (100,3,0). Upstream returns closest = (2,3,2) with
        // distance 98.0204; the frustum only reaches |u| <= 2 at d = 2.
        const p = v(0, 3, 100);
        const r = query.compute(p, frustum);
        expect(r.distance).toBeCloseTo(nearestVertexDistance(p), 12);
        expect(r.distance).toBeCloseTo(Math.sqrt(98 * 98 + 1 + 4), 12);
        expectVectorClose(r.closest[1], v(2, 2, 2), 1e-12, 1e-12);
    });

    it('clamps the UF-edge r coordinate to the far half-extent', () => {
        // The mirror configuration (r,u,d) = (3,100,0).
        const p = v(0, 100, 3);
        const r = query.compute(p, frustum);
        expect(r.distance).toBeCloseTo(nearestVertexDistance(p), 12);
        expectVectorClose(r.closest[1], v(2, 2, 2), 1e-12, 1e-12);
    });

    it('keeps the closest point inside for the whole far-region family', () => {
        for (const dCoord of [0, 0.5, 1, 1.5, 2, 3]) {
            for (const rCoord of [2.5, 5, 20, 100]) {
                for (const uCoord of [2.5, 5, 20, 100]) {
                    for (const sr of [1, -1]) {
                        for (const su of [1, -1]) {
                            const p = v(dCoord, su * uCoord, sr * rCoord);
                            const r = query.compute(p, frustum);
                            // Every far-region point of this family is closest
                            // to a vertex, an edge or the far face; the LUF
                            // vertex bounds it from above.
                            expect(r.distance)
                                .toBeLessThanOrEqual(
                                    nearestVertexDistance(p) + 1e-12);
                            // The closest point must be a frustum point:
                            // |u| <= d and |r| <= d for this frustum, with
                            // 1 <= d <= 2.
                            const cd = r.closest[1].values[0];
                            const cu = Math.abs(r.closest[1].values[1]);
                            const cr = Math.abs(r.closest[1].values[2]);
                            expect(cd).toBeGreaterThanOrEqual(1 - 1e-12);
                            expect(cd).toBeLessThanOrEqual(2 + 1e-12);
                            expect(cu).toBeLessThanOrEqual(cd + 1e-12);
                            expect(cr).toBeLessThanOrEqual(cd + 1e-12);
                        }
                    }
                }
            }
        }
    });
});
