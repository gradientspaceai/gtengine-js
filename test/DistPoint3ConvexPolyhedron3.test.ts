import { describe, expect, it } from 'vitest';
import { ConvexPolyhedron3 } from '../src/ConvexPolyhedron3.js';
import { DistPoint3ConvexPolyhedron3 }
    from '../src/DistPoint3ConvexPolyhedron3.js';
import { Vector, add, dot, hlift, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// The unit cube [0,1]^3 with faces wound counterclockwise when viewed from
// outside, so the generated plane normals point outward.
function cube(): ConvexPolyhedron3 {
    const vertices = [
        v(0, 0, 0), v(1, 0, 0), v(1, 1, 0), v(0, 1, 0),
        v(0, 0, 1), v(1, 0, 1), v(1, 1, 1), v(0, 1, 1)
    ];
    const indices = [
        0, 3, 2, 0, 2, 1,
        4, 5, 6, 4, 6, 7,
        0, 1, 5, 0, 5, 4,
        3, 7, 6, 3, 6, 2,
        0, 4, 7, 0, 7, 3,
        1, 2, 6, 1, 6, 5
    ];
    return new ConvexPolyhedron3(vertices, indices, true, true);
}

function tetrahedron(): ConvexPolyhedron3 {
    const vertices = [v(0, 0, 0), v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
    const indices = [0, 2, 1, 0, 1, 3, 0, 3, 2, 1, 2, 3];
    return new ConvexPolyhedron3(vertices, indices, true, true);
}

// The exact distance from a point to the solid unit cube [0,1]^3.
function cubeDistance(p: Vector): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const d = Math.max(0, -p.values[i], p.values[i] - 1);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

describe('DistPoint3ConvexPolyhedron3', () => {
    it('computes the distance to a cube face', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 0.5, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(1, 8);
        expect(result.closest[0].values).toEqual([2, 0.5, 0.5]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 8);
        expect(result.closest[1].values[1]).toBeCloseTo(0.5, 8);
        expect(result.closest[1].values[2]).toBeCloseTo(0.5, 8);
    });

    it('computes the distance to a cube corner', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 2, 2), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 7);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 7);
        }
    });

    it('computes the distance to a cube edge', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(2, 2, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 7);
    });

    it('returns zero distance for a point inside', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(0.5, 0.25, 0.75), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(0, 8);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i])
                .toBeCloseTo(result.closest[0].values[i], 8);
        }
    });

    it('returns zero distance for a point on the boundary', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const result = query.compute(v(1, 0.5, 0.5), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(0, 8);
    });

    it('handles a tetrahedron', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const t = tetrahedron();
        // The point (1,1,1) is outside the plane x+y+z=1; the closest point
        // is the face centroid (1/3,1/3,1/3) and the distance is
        // (3-1)/sqrt(3).
        const result = query.compute(v(1, 1, 1), t);
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(2 / Math.sqrt(3), 7);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1 / 3, 6);
        }

        const inside = query.compute(v(0.1, 0.1, 0.1), t);
        expect(inside.distance).toBeCloseTo(0, 8);
    });

    it('reports failure for a polyhedron with no planes', () => {
        const query = new DistPoint3ConvexPolyhedron3();
        const empty = new ConvexPolyhedron3();
        const result = query.compute(v(1, 2, 3), empty);
        expect(result.queryIsSuccessful).toBe(false);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(result.numLCPIterations).toBe(0);
        expect(result.closest[0].values).toEqual([0, 0, 0]);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('reuses a preallocated LCP solver', () => {
        // The cube has 12 triangle faces.
        const query = new DistPoint3ConvexPolyhedron3(12);
        const c = cube();
        const r0 = query.compute(v(3, 0.5, 0.5), c);
        const r1 = query.compute(v(-2, 0.5, 0.5), c);
        expect(r0.distance).toBeCloseTo(2, 7);
        expect(r1.distance).toBeCloseTo(2, 7);
        expect(r0.numLCPIterations).toBeGreaterThan(0);
    });

    it('honors setMaxLCPIterations', () => {
        const query = new DistPoint3ConvexPolyhedron3(12);
        query.setMaxLCPIterations(500);
        const result = query.compute(v(2, 2, 2), cube());
        expect(result.queryIsSuccessful).toBe(true);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 7);
    });

    it('agrees with the exact cube distance on random points', () => {
        const query = new DistPoint3ConvexPolyhedron3(12);
        const c = cube();
        let seed = 456789123;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const p = v(4 * rand() - 1.5, 4 * rand() - 1.5, 4 * rand() - 1.5);
            const result = query.compute(p, c);
            expect(result.queryIsSuccessful).toBe(true);
            expect(result.distance).toBeCloseTo(cubeDistance(p), 6);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 6);
            // The closest point is in the cube.
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i])
                    .toBeGreaterThanOrEqual(-1e-6);
                expect(result.closest[1].values[i])
                    .toBeLessThanOrEqual(1 + 1e-6);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistPoint3ConvexPolyhedron3.h.
// ---------------------------------------------------------------------------

// A point is inside the polyhedron when Dot(plane, (X,1)) <= 0 for every
// outward-pointing face plane. The value is normalized by the length of the
// plane normal (ConvexPolyhedron3 does not normalize its generated planes) so
// that the returned quantity is a signed distance.
function maxPlaneValue(poly: ConvexPolyhedron3, x: Vector): number {
    const h = hlift(x, 1);
    let worst = -Number.MAX_VALUE;
    for (const plane of poly.planes) {
        const norm = length(Vector.fromArray(
            [plane.values[0], plane.values[1], plane.values[2]]));
        const value = dot(plane, h) / norm;
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

// Brute-force distance from a point to the solid tetrahedron, by a dense
// barycentric sampling followed by a local pattern-search refinement.
function tetraBruteForce(p: Vector, poly: ConvexPolyhedron3): number {
    const V = poly.vertices;
    const at = (a: number, b: number, c: number): number => {
        const d = 1 - a - b - c;
        const x = add(add(mul(a, V[0]), mul(b, V[1])),
            add(mul(c, V[2]), mul(d, V[3])));
        return length(sub(p, x));
    };
    const n = 24;
    let best = Number.MAX_VALUE;
    let ba = 0;
    let bb = 0;
    let bc = 0;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            for (let k = 0; i + j + k <= n; ++k) {
                const d = at(i / n, j / n, k / n);
                if (d < best) {
                    best = d;
                    ba = i / n;
                    bb = j / n;
                    bc = k / n;
                }
            }
        }
    }
    let h = 1 / n;
    for (let pass = 0; pass < 100; ++pass) {
        for (const [da, db, dc] of [[1, 0, 0], [-1, 0, 0], [0, 1, 0],
            [0, -1, 0], [0, 0, 1], [0, 0, -1], [1, -1, 0], [-1, 1, 0],
            [1, 0, -1], [-1, 0, 1], [0, 1, -1], [0, -1, 1]]) {
            const a = ba + da * h;
            const b = bb + db * h;
            const c = bc + dc * h;
            if (a < 0 || b < 0 || c < 0 || a + b + c > 1) {
                continue;
            }
            const d = at(a, b, c);
            if (d < best) {
                best = d;
                ba = a;
                bb = b;
                bc = c;
            }
        }
        h *= 0.7;
    }
    return best;
}

describe('DistPoint3ConvexPolyhedron3 verification', () => {
    const pointArb = wellScaledVector(3, -3, 3);

    it('result is self consistent and the closest point is in the cube',
        () => {
            const poly = cube();
            check(pointArb, p => {
                const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
                expect(r.queryIsSuccessful).toBe(true);
                expectClose(r.sqrDistance, r.distance * r.distance, 1e-12,
                    1e-12);
                expectVectorClose(r.closest[0], p, 0, 0);
                expectClose(r.distance,
                    length(sub(r.closest[0], r.closest[1])), 1e-9, 1e-9);
                expect(maxPlaneValue(poly, r.closest[1]))
                    .toBeLessThanOrEqual(1e-8);
            }, 60);
        }, 30000);

    it('matches the exact cube distance', () => {
        const poly = cube();
        check(pointArb, p => {
            const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            expectClose(r.distance, cubeDistance(p), 1e-7, 1e-7);
        }, 60);
    }, 30000);

    it('matches a brute-force minimization over the tetrahedron', () => {
        const poly = tetrahedron();
        check(pointArb, p => {
            const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            expect(r.queryIsSuccessful).toBe(true);
            expectClose(r.distance, tetraBruteForce(p, poly), 1e-5, 1e-5);
        }, 25);
    }, 30000);

    it('is invariant under a common translation', () => {
        check(fc.tuple(pointArb, wellScaledVector(3, -4, 4)), ([p, tr]) => {
            const poly = cube();
            const moved = new ConvexPolyhedron3(
                poly.vertices.map(x => add(x, tr)), poly.indices, true, true);
            const r0 = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            const r1 = new DistPoint3ConvexPolyhedron3()
                .compute(add(p, tr), moved);
            expectClose(r0.distance, r1.distance, 1e-6, 1e-6);
        }, 40);
    }, 30000);

    it('a preconstructed solver agrees with the per-query solver', () => {
        // Upstream builds the LCP solver once when numTriangles > 0; the
        // results must not depend on that choice, and the solver must be
        // reusable across calls.
        const poly = cube();
        const shared = new DistPoint3ConvexPolyhedron3(poly.planes.length);
        check(pointArb, p => {
            const r0 = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            const r1 = shared.compute(p, poly);
            expect(r1.queryIsSuccessful).toBe(r0.queryIsSuccessful);
            expectClose(r1.distance, r0.distance, 1e-12, 1e-12);
            expectVectorClose(r1.closest[1], r0.closest[1], 1e-12, 1e-12);
        }, 60);
    }, 30000);

    it('reports zero distance for interior points', () => {
        const poly = cube();
        check(fc.tuple(fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true })),
        ([x, y, z]) => {
            const p = Vector.fromArray([x, y, z]);
            const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            expect(r.queryIsSuccessful).toBe(true);
            expect(r.distance).toBeLessThanOrEqual(1e-7);
        }, 40);
    }, 30000);

    it('is not larger than the distance to any sampled polyhedron point',
        () => {
            const poly = cube();
            check(fc.tuple(pointArb, fc.array(wellScaled(0, 1),
                { minLength: 3, maxLength: 3 })), ([p, s]) => {
                const d = new DistPoint3ConvexPolyhedron3()
                    .compute(p, poly).distance;
                const q = Vector.fromArray([Math.min(1, Math.max(0, s[0])),
                    Math.min(1, Math.max(0, s[1])),
                    Math.min(1, Math.max(0, s[2]))]);
                expect(d).toBeLessThanOrEqual(length(sub(p, q)) + 1e-7);
            }, 60);
        }, 30000);

    it('reports failure with zeroed output when there are no face planes',
        () => {
            const poly = new ConvexPolyhedron3();
            check(pointArb, p => {
                const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
                expect(r.queryIsSuccessful).toBe(false);
                expect(r.distance).toBe(0);
                expect(r.sqrDistance).toBe(0);
                expect(r.numLCPIterations).toBe(0);
                expect(r.closest[0].values).toEqual([0, 0, 0]);
                expect(r.closest[1].values).toEqual([0, 0, 0]);
            }, 20);
        });

    it('does not mutate its inputs', () => {
        const poly = cube();
        check(pointArb, p => {
            const p0 = p.clone();
            const verts = poly.vertices.map(x => x.clone());
            const lo = poly.alignedBox.min.clone();
            const hi = poly.alignedBox.max.clone();
            const r = new DistPoint3ConvexPolyhedron3().compute(p, poly);
            expect(p.values).toEqual(p0.values);
            poly.vertices.forEach((x, i) =>
                expect(x.values).toEqual(verts[i].values));
            expect(poly.alignedBox.min.values).toEqual(lo.values);
            expect(poly.alignedBox.max.values).toEqual(hi.values);
            r.closest[0].values[0] = 987;
            expect(p.values).toEqual(p0.values);
        }, 40);
    }, 30000);
});
