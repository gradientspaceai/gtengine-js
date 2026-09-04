import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2.js';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh.js';
import { Vector } from '../src/Vector.js';
import { fc, check, latticeVector, expectClose } from './helpers/arbitraries.js';
import { orient2 } from './helpers/exact.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function makeTriangulation(points: readonly Vector[]): Delaunay2 {
    const delaunay = new Delaunay2();
    expect(delaunay.compute(points)).toBe(true);
    return delaunay;
}

describe('Delaunay2Mesh', () => {
    it('requires a 2-dimensional triangulation', () => {
        const collinear = new Delaunay2();
        expect(collinear.compute([v2(0, 0), v2(1, 1), v2(2, 2)])).toBe(false);
        expect(() => new Delaunay2Mesh(collinear)).toThrow();
    });

    it('mirrors the accessors of the underlying triangulation', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1), v2(0.4, 0.6)];
        const delaunay = makeTriangulation(points);
        const mesh = new Delaunay2Mesh(delaunay);

        expect(mesh.getNumVertices()).toBe(delaunay.getNumVertices());
        expect(mesh.getNumVertices()).toBe(points.length);
        expect(mesh.getNumTriangles()).toBe(delaunay.getNumTriangles());
        expect(mesh.getVertices()).toBe(delaunay.getVertices());
        expect(mesh.getIndices()).toEqual(delaunay.getIndices());
        expect(mesh.getAdjacencies()).toEqual(delaunay.getAdjacencies());
        expect(mesh.getInvalidIndex()).toBe(-1);

        for (let t = 0; t < mesh.getNumTriangles(); ++t) {
            const indices = mesh.getTriangleIndices(t);
            expect(indices).toEqual(delaunay.getTriangleIndices(t));
            expect(indices).not.toBeNull();
            const adjacencies = mesh.getTriangleAdjacencies(t);
            expect(adjacencies).toEqual(delaunay.getTriangleAdjacencies(t));

            const vertices = mesh.getTriangleVertices(t);
            expect(vertices).not.toBeNull();
            for (let i = 0; i < 3; ++i) {
                expect((vertices as Vector[])[i].equals(
                    points[(indices as number[])[i]])).toBe(true);
            }

            // The compact index array agrees with the per-triangle accessor.
            for (let i = 0; i < 3; ++i) {
                expect((indices as number[])[i]).toBe(mesh.getIndices()[3 * t + i]);
                expect((adjacencies as number[])[i])
                    .toBe(mesh.getAdjacencies()[3 * t + i]);
            }
        }
    });

    it('returns null for out-of-range triangle indices', () => {
        const points = [v2(0, 0), v2(1, 0), v2(0, 1)];
        const mesh = new Delaunay2Mesh(makeTriangulation(points));
        expect(mesh.getNumTriangles()).toBe(1);
        for (const t of [-1, 1, 100]) {
            expect(mesh.getTriangleIndices(t)).toBeNull();
            expect(mesh.getTriangleAdjacencies(t)).toBeNull();
            expect(mesh.getTriangleVertices(t)).toBeNull();
            expect(mesh.getBarycentrics(t, v2(0.25, 0.25))).toBeNull();
        }
    });

    it('computes exact barycentric coordinates', () => {
        const points = [v2(0, 0), v2(4, 0), v2(0, 4)];
        const mesh = new Delaunay2Mesh(makeTriangulation(points));
        const indices = mesh.getTriangleIndices(0) as number[];

        // A point whose barycentric coordinates are exact binary fractions.
        const p = v2(1, 2);
        const bary = mesh.getBarycentrics(0, p) as number[];
        expect(bary).not.toBeNull();
        expect(bary[0] + bary[1] + bary[2]).toBe(1);

        // Reconstruct the point from the barycentric coordinates.
        for (let j = 0; j < 2; ++j) {
            let value = 0;
            for (let i = 0; i < 3; ++i) {
                value += bary[i] * points[indices[i]].values[j];
            }
            expect(value).toBeCloseTo(p.values[j], 12);
        }

        // A vertex of the triangle has a unit barycentric coordinate.
        const atVertex = mesh.getBarycentrics(0, points[indices[1]]) as number[];
        expect(atVertex).toEqual([0, 1, 0]);
    });

    it('locates containing triangles and cross-checks barycentrics', () => {
        const random = makeRandom(9091);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            points.push(v2(random() * 4 - 2, random() * 4 - 2));
        }
        const delaunay = makeTriangulation(points);
        const mesh = new Delaunay2Mesh(delaunay);

        let numFound = 0;
        for (let trial = 0; trial < 40; ++trial) {
            const p = v2(random() * 4 - 2, random() * 4 - 2);
            const t = mesh.getContainingTriangle(p);
            if (t === mesh.getInvalidIndex()) {
                continue;
            }
            ++numFound;
            expect(t).toBeGreaterThanOrEqual(0);
            expect(t).toBeLessThan(mesh.getNumTriangles());

            const bary = mesh.getBarycentrics(t, p) as number[];
            expect(bary).not.toBeNull();
            // The point is inside the triangle, so the barycentrics are
            // nonnegative and sum to one.
            expect(bary[0] + bary[1] + bary[2]).toBeCloseTo(1, 12);
            for (const b of bary) {
                expect(b).toBeGreaterThanOrEqual(-1e-12);
            }

            // The barycentric combination reproduces the query point.
            const indices = mesh.getTriangleIndices(t) as number[];
            for (let j = 0; j < 2; ++j) {
                let value = 0;
                for (let i = 0; i < 3; ++i) {
                    value += bary[i] * points[indices[i]].values[j];
                }
                expect(value).toBeCloseTo(p.values[j], 10);
            }
        }
        expect(numFound).toBeGreaterThan(0);

        // A point far outside the convex hull is not in any triangle.
        expect(mesh.getContainingTriangle(v2(100, 100)))
            .toBe(mesh.getInvalidIndex());
    });
});

// ---------------------------------------------------------------------------
// Verification (V11): property-based cross-checks. All point sets are integer
// lattices, so the containment oracle is exact bigint orientation arithmetic.
// ---------------------------------------------------------------------------

const bigXY = (p: Vector): [bigint, bigint] =>
    [BigInt(p.values[0]), BigInt(p.values[1])];

/**
 * The indices of every triangle of the triangulation that contains P,
 * determined exactly. Triangles of Delaunay2 are counterclockwise, so P is in
 * the closed triangle when all three edge orientations are >= 0.
 */
function containingTrianglesExact(mesh: Delaunay2Mesh, p: Vector): number[] {
    const vertices = mesh.getVertices();
    const [px, py] = bigXY(p);
    const result: number[] = [];
    for (let t = 0; t < mesh.getNumTriangles(); ++t) {
        const idx = mesh.getTriangleIndices(t) as [number, number, number];
        let inside = true;
        for (let i = 0; i < 3; ++i) {
            const [ax, ay] = bigXY(vertices[idx[i]]);
            const [bx, by] = bigXY(vertices[idx[(i + 1) % 3]]);
            if (orient2(ax, ay, bx, by, px, py) < 0) {
                inside = false;
                break;
            }
        }
        if (inside) { result.push(t); }
    }
    return result;
}

const latticeCloud2 = (count: number, range: number): fc.Arbitrary<Vector[]> =>
    fc.array(latticeVector(2, -range, range),
        { minLength: count, maxLength: count });

function buildMesh(points: readonly Vector[]): Delaunay2Mesh | null {
    const delaunay = new Delaunay2();
    if (!delaunay.compute(points) || delaunay.getDimension() !== 2) {
        return null;
    }
    return new Delaunay2Mesh(delaunay);
}

describe('Delaunay2Mesh verification', () => {
    it('locates the containing triangle exactly', () => {
        check(fc.tuple(latticeCloud2(10, 5), latticeVector(2, -6, 6)),
            ([points, query]) => {
                const mesh = buildMesh(points);
                if (mesh === null) { return true; }
                const expected = containingTrianglesExact(mesh, query);
                const t = mesh.getContainingTriangle(query);
                if (expected.length === 0) {
                    // Outside the convex hull.
                    expect(t).toBe(mesh.getInvalidIndex());
                }
                else {
                    // On a shared edge or vertex several triangles contain
                    // the point; the search may return any of them.
                    expect(expected).toContain(t);
                }
                return true;
            }, 150);
    }, 30000);

    it('computes barycentrics that reproduce the query point', () => {
        check(fc.tuple(latticeCloud2(10, 5), latticeVector(2, -6, 6)),
            ([points, query]) => {
                const mesh = buildMesh(points);
                if (mesh === null) { return true; }
                const t = mesh.getContainingTriangle(query);
                if (t === mesh.getInvalidIndex()) { return true; }

                const bary = mesh.getBarycentrics(t, query);
                expect(bary).not.toBeNull();
                const b = bary as [number, number, number];
                const vertices = mesh.getTriangleVertices(t) as
                    [Vector, Vector, Vector];

                // The point is inside the closed triangle, so the
                // barycentrics are in [0, 1] and sum to 1.
                expectClose(b[0] + b[1] + b[2], 1, 1e-12, 1e-12);
                for (const bi of b) {
                    expect(bi).toBeGreaterThanOrEqual(-1e-12);
                    expect(bi).toBeLessThanOrEqual(1 + 1e-12);
                }

                // The barycentric combination reproduces the query point.
                // The coordinates are exact rationals rounded once to
                // double, so the reconstruction is accurate to a few ulps of
                // the coordinate range.
                for (let j = 0; j < 2; ++j) {
                    const value = b[0] * vertices[0].values[j]
                        + b[1] * vertices[1].values[j]
                        + b[2] * vertices[2].values[j];
                    expectClose(value, query.values[j], 1e-11, 1e-11);
                }
                return true;
            }, 150);
    }, 30000);

    it('interpolates a linear function exactly', () => {
        check(fc.tuple(latticeCloud2(10, 5), latticeVector(2, -6, 6),
            latticeVector(3, -4, 4)), ([points, query, coeff]) => {
            const mesh = buildMesh(points);
            if (mesh === null) { return true; }
            const t = mesh.getContainingTriangle(query);
            if (t === mesh.getInvalidIndex()) { return true; }
            const bary = mesh.getBarycentrics(t, query);
            if (bary === null) { return true; }

            const f = (v: Vector): number =>
                coeff.values[0] * v.values[0] + coeff.values[1] * v.values[1]
                + coeff.values[2];
            const vertices = mesh.getTriangleVertices(t) as
                [Vector, Vector, Vector];
            const interpolated = bary[0] * f(vertices[0])
                + bary[1] * f(vertices[1]) + bary[2] * f(vertices[2]);
            expectClose(interpolated, f(query), 1e-10, 1e-10);
            return true;
        }, 150);
    }, 30000);

    it('keeps the accessors consistent with the triangulation', () => {
        check(latticeCloud2(11, 6), points => {
            const mesh = buildMesh(points);
            if (mesh === null) { return true; }
            const vertices = mesh.getVertices();
            const indices = mesh.getIndices();
            const adjacencies = mesh.getAdjacencies();
            expect(indices.length).toBe(3 * mesh.getNumTriangles());
            expect(adjacencies.length).toBe(indices.length);

            for (let t = 0; t < mesh.getNumTriangles(); ++t) {
                const idx = mesh.getTriangleIndices(t) as [number, number, number];
                const adj = mesh.getTriangleAdjacencies(t) as
                    [number, number, number];
                const tv = mesh.getTriangleVertices(t) as [Vector, Vector, Vector];
                for (let j = 0; j < 3; ++j) {
                    expect(idx[j]).toBe(indices[3 * t + j]);
                    expect(adj[j]).toBe(adjacencies[3 * t + j]);
                    expect(tv[j].values).toEqual(vertices[idx[j]].values);
                    // The returned vertices are copies, not aliases of the
                    // triangulation storage (upstream returns by value).
                    expect(tv[j]).not.toBe(vertices[idx[j]]);
                    if (adj[j] !== mesh.getInvalidIndex()) {
                        const back = mesh.getTriangleAdjacencies(adj[j]) as
                            [number, number, number];
                        expect(back).toContain(t);
                    }
                }
            }
            expect(mesh.getTriangleIndices(-1)).toBeNull();
            expect(mesh.getTriangleIndices(mesh.getNumTriangles())).toBeNull();
            expect(mesh.getBarycentrics(mesh.getNumTriangles(),
                points[0])).toBeNull();
            return true;
        }, 100);
    }, 30000);
});
