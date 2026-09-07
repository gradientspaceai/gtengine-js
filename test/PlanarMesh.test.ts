import { describe, it, expect } from 'vitest';
import { PlanarMesh } from '../src/PlanarMesh.js';
import { ETManifoldMesh } from '../src/ETManifoldMesh.js';
import { Vector } from '../src/Vector.js';
import { Delaunay2 } from '../src/Delaunay2.js';
import { IntpLinearNonuniform2 } from '../src/IntpLinearNonuniform2.js';
import type { IntpLinearNonuniform2TriangleMesh } from '../src/IntpLinearNonuniform2.js';
import type { IntpQuadraticNonuniform2TriangleMesh } from '../src/IntpQuadraticNonuniform2.js';
import { check, expectClose, fc, seededRandom } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The unit square split by the diagonal <0,2>, counterclockwise triangles.
const squareVertices = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
const squareIndices = [0, 1, 2, 0, 2, 3];

// A triangulated (n x n) grid of the unit square. The mesh is convex, so the
// linear-walk containment query is guaranteed to succeed for interior points.
function makeGrid(n: number): { vertices: Vector[], indices: number[] } {
    const vertices: Vector[] = [];
    for (let r = 0; r <= n; ++r) {
        for (let c = 0; c <= n; ++c) {
            vertices.push(v2(c / n, r / n));
        }
    }
    const index = (r: number, c: number) => r * (n + 1) + c;
    const indices: number[] = [];
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            indices.push(index(r, c), index(r, c + 1), index(r + 1, c + 1));
            indices.push(index(r, c), index(r + 1, c + 1), index(r + 1, c));
        }
    }
    return { vertices, indices };
}

describe('PlanarMesh', () => {
    it('builds the indices and adjacencies from an index array', () => {
        const mesh = PlanarMesh.fromIndices(squareVertices, squareIndices);
        expect(mesh.getNumVertices()).toBe(4);
        expect(mesh.getNumTriangles()).toBe(2);
        expect(mesh.getVertices()).toBe(squareVertices);
        expect(Array.from(mesh.getIndices())).toEqual(squareIndices);

        // Triangle 0 = <0,1,2> has edges <0,1>, <1,2>, <2,0>; only <2,0> is
        // shared, with triangle 1. Triangle 1 = <0,2,3> has edges <0,2>,
        // <2,3>, <3,0>; only <0,2> is shared, with triangle 0.
        expect(Array.from(mesh.getAdjacencies())).toEqual([-1, -1, 1, 0, -1, -1]);
        expect(mesh.getTriangleAdjacencies(0)).toEqual([-1, -1, 1]);
        expect(mesh.getTriangleAdjacencies(1)).toEqual([0, -1, -1]);

        expect(mesh.getTriangleIndices(0)).toEqual([0, 1, 2]);
        expect(mesh.getTriangleIndices(1)).toEqual([0, 2, 3]);
        expect(mesh.getTriangleVertices(1)).toEqual(
            [squareVertices[0], squareVertices[2], squareVertices[3]]);

        // Out-of-range triangle indices.
        expect(mesh.getTriangleIndices(-1)).toBeNull();
        expect(mesh.getTriangleIndices(2)).toBeNull();
        expect(mesh.getTriangleVertices(2)).toBeNull();
        expect(mesh.getTriangleAdjacencies(2)).toBeNull();
        expect(mesh.getBarycentrics(2, v2(0.5, 0.5))).toBeNull();
    });

    it('builds the same mesh from an ETManifoldMesh', () => {
        const etMesh = new ETManifoldMesh();
        etMesh.insert(0, 1, 2);
        etMesh.insert(0, 2, 3);
        const mesh = PlanarMesh.fromMesh(squareVertices, etMesh);

        // The triangles are ordered by increasing triangle key, which here is
        // the input order.
        expect(mesh.getNumTriangles()).toBe(2);
        expect(Array.from(mesh.getIndices())).toEqual(squareIndices);
        expect(Array.from(mesh.getAdjacencies())).toEqual([-1, -1, 1, 0, -1, -1]);
    });

    it('validates the constructor inputs', () => {
        expect(() => PlanarMesh.fromIndices([v2(0, 0), v2(1, 0)], [0, 1, 0])).toThrow();
        expect(() => PlanarMesh.fromIndices(squareVertices, [0, 1])).toThrow();
        expect(() => PlanarMesh.fromIndices(squareVertices, [0, 1, 2, 3])).toThrow();
        expect(() => PlanarMesh.fromMesh(squareVertices, new ETManifoldMesh())).toThrow();

        // Upstream quirk: a repeated triangle stops the construction silently
        // and the resulting object has no triangles and no vertices.
        const degenerate = PlanarMesh.fromIndices(squareVertices, [0, 1, 2, 1, 2, 0]);
        expect(degenerate.getNumTriangles()).toBe(0);
        expect(degenerate.getNumVertices()).toBe(0);
    });

    it('finds the containing triangle of the two-triangle square', () => {
        const mesh = PlanarMesh.fromIndices(squareVertices, squareIndices);

        // Below the diagonal is triangle 0, above it is triangle 1.
        expect(mesh.getContainingTriangle(v2(0.75, 0.25))).toBe(0);
        expect(mesh.getContainingTriangle(v2(0.25, 0.75))).toBe(1);
        // The walk finds the triangle from either start.
        expect(mesh.getContainingTriangle(v2(0.75, 0.25), 1)).toBe(0);
        expect(mesh.getContainingTriangle(v2(0.25, 0.75), 0)).toBe(1);
        // A point outside the mesh exits through a boundary edge.
        expect(mesh.getContainingTriangle(v2(2, 2))).toBe(-1);
        expect(mesh.getContainingTriangle(v2(-1, 0.5))).toBe(-1);

        // The cycle-trapping overload agrees and reports the visited set.
        const found = mesh.getContainingTriangleVisited(v2(0.25, 0.75), 0);
        expect(found.triangle).toBe(1);
        expect(Array.from(found.visited).sort()).toEqual([0, 1]);
        const missing = mesh.getContainingTriangleVisited(v2(2, 2), 0);
        expect(missing.triangle).toBe(-1);
        expect(missing.visited.has(0)).toBe(true);
    });

    it('tests containment of a specific triangle', () => {
        const mesh = PlanarMesh.fromIndices(squareVertices, squareIndices);
        expect(mesh.contains(0, v2(0.75, 0.25))).toBe(true);
        expect(mesh.contains(1, v2(0.75, 0.25))).toBe(false);
        expect(mesh.contains(1, v2(0.25, 0.75))).toBe(true);
        expect(mesh.contains(0, v2(0.25, 0.75))).toBe(false);
        expect(mesh.contains(0, v2(5, 5))).toBe(false);
        expect(() => mesh.contains(2, v2(0.5, 0.5))).toThrow();
    });

    it('computes barycentric coordinates that reproduce the point', () => {
        const mesh = PlanarMesh.fromIndices(squareVertices, squareIndices);
        const p = v2(0.75, 0.25);
        const bary = mesh.getBarycentrics(0, p);
        expect(bary).not.toBeNull();
        const b = bary!;
        expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 12);
        for (const value of b) {
            expect(value).toBeGreaterThanOrEqual(0);
        }
        const tri = mesh.getTriangleVertices(0)!;
        for (let j = 0; j < 2; ++j) {
            const value = b[0] * tri[0].values[j] + b[1] * tri[1].values[j]
                + b[2] * tri[2].values[j];
            expect(value).toBeCloseTo(p.values[j], 12);
        }

        // A vertex of the triangle has a unit barycentric coordinate.
        expect(mesh.getBarycentrics(0, squareVertices[1])).toEqual([0, 1, 0]);

        // A degenerate triangle has no barycentric coordinates.
        const collinear = [v2(0, 0), v2(1, 0), v2(2, 0), v2(0, 1)];
        const degenerate = PlanarMesh.fromIndices(collinear, [0, 1, 2, 0, 2, 3]);
        expect(degenerate.getBarycentrics(0, v2(0.5, 0))).toBeNull();
    });

    it('agrees with an exhaustive search on a randomized grid mesh', () => {
        const n = 5;
        const { vertices, indices } = makeGrid(n);
        const mesh = PlanarMesh.fromIndices(vertices, indices);
        expect(mesh.getNumTriangles()).toBe(2 * n * n);
        expect(mesh.getNumVertices()).toBe((n + 1) * (n + 1));

        // The adjacency graph is symmetric: if triangle t has neighbor a
        // across some edge, a has t as a neighbor across the same edge.
        const adjacencies = mesh.getAdjacencies();
        let numBoundary = 0;
        for (let t = 0; t < mesh.getNumTriangles(); ++t) {
            for (let i = 0; i < 3; ++i) {
                const a = adjacencies[3 * t + i];
                if (a === -1) {
                    ++numBoundary;
                    continue;
                }
                expect(Array.from(mesh.getTriangleAdjacencies(a)!)).toContain(t);
                // The shared edge has the same unordered vertex pair.
                const ti = mesh.getTriangleIndices(t)!;
                const shared = [ti[i], ti[(i + 1) % 3]].sort((x, y) => x - y);
                const ai = mesh.getTriangleIndices(a)!;
                const j = mesh.getTriangleAdjacencies(a)!.indexOf(t);
                const sharedA = [ai[j], ai[(j + 1) % 3]].sort((x, y) => x - y);
                expect(sharedA).toEqual(shared);
            }
        }
        // The boundary of the unit square has 4*n edges.
        expect(numBoundary).toBe(4 * n);

        // Random interior points: the linear walk from a random start must
        // agree with an exhaustive containment search, and the barycentric
        // coordinates of the reported triangle must be nonnegative and
        // reproduce the point.
        let seed = 13579;
        const nextRandom = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 200; ++trial) {
            const p = v2(nextRandom(), nextRandom());
            const start = Math.floor(nextRandom() * mesh.getNumTriangles());
            const t = mesh.getContainingTriangle(p, start);
            expect(t).toBeGreaterThanOrEqual(0);
            expect(mesh.contains(t, p)).toBe(true);

            const found = mesh.getContainingTriangleVisited(p, start);
            expect(found.triangle).toBe(t);
            expect(found.visited.has(start)).toBe(true);

            const b = mesh.getBarycentrics(t, p)!;
            expect(b).not.toBeNull();
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 12);
            for (const value of b) {
                expect(value).toBeGreaterThan(-1e-12);
            }
            const tri = mesh.getTriangleVertices(t)!;
            for (let j = 0; j < 2; ++j) {
                const value = b[0] * tri[0].values[j] + b[1] * tri[1].values[j]
                    + b[2] * tri[2].values[j];
                expect(value).toBeCloseTo(p.values[j], 12);
            }

            // The exhaustive search agrees; a point can be reported by more
            // than one triangle only when it lies on a shared edge.
            let numContaining = 0;
            for (let u = 0; u < mesh.getNumTriangles(); ++u) {
                if (mesh.contains(u, p)) {
                    ++numContaining;
                }
            }
            expect(numContaining).toBeGreaterThanOrEqual(1);
        }
    });

    it('produces the same mesh from indices and from an equivalent ETManifoldMesh', () => {
        const n = 3;
        const { vertices, indices } = makeGrid(n);
        const fromIndices = PlanarMesh.fromIndices(vertices, indices);

        const etMesh = new ETManifoldMesh();
        for (let t = 0; t < indices.length / 3; ++t) {
            etMesh.insert(indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]);
        }
        const fromMesh = PlanarMesh.fromMesh(vertices, etMesh);

        expect(fromMesh.getNumTriangles()).toBe(fromIndices.getNumTriangles());

        // The triangle orderings differ (key order versus input order), so
        // compare the triangles as sets of index triples and check that the
        // adjacency relations agree.
        const tripleKey = (mesh: PlanarMesh, t: number) => {
            const i = mesh.getTriangleIndices(t)!;
            const rotated = i.slice() as number[];
            const min = Math.min(...rotated);
            while (rotated[0] !== min) {
                rotated.push(rotated.shift()!);
            }
            return rotated.join(',');
        };

        const keysA: string[] = [];
        const keysB: string[] = [];
        for (let t = 0; t < fromIndices.getNumTriangles(); ++t) {
            keysA.push(tripleKey(fromIndices, t));
            keysB.push(tripleKey(fromMesh, t));
        }
        expect(keysB.slice().sort()).toEqual(keysA.slice().sort());

        // The neighbor triples of each triangle agree between the two meshes.
        const neighbors = (mesh: PlanarMesh, keys: string[], t: number) => {
            return mesh.getTriangleAdjacencies(t)!.map(
                a => (a === -1 ? '' : keys[a]));
        };
        for (let t = 0; t < fromIndices.getNumTriangles(); ++t) {
            const u = keysB.indexOf(keysA[t]);
            expect(u).toBeGreaterThanOrEqual(0);
            // Both triangles start at the same smallest vertex index after
            // the rotation used by tripleKey, but the stored triples may be
            // rotations of each other, so compare the neighbor sets.
            expect(neighbors(fromMesh, keysB, u).slice().sort())
                .toEqual(neighbors(fromIndices, keysA, t).slice().sort());
        }

        // Containment queries agree.
        for (let i = 1; i < 10; ++i) {
            const p = v2(i / 10, (10 - i) / 11);
            const ta = fromIndices.getContainingTriangle(p);
            const tb = fromMesh.getContainingTriangle(p);
            expect(ta).toBeGreaterThanOrEqual(0);
            expect(tb).toBeGreaterThanOrEqual(0);
            expect(keysB[tb]).toBe(keysA[ta]);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V37): independent re-check against PlanarMesh.h.
// ---------------------------------------------------------------------------

// A Delaunay triangulation of a set of lattice points, as an index array in
// the input-vertex numbering. Returns null when the points are degenerate
// (fewer than three distinct points, or all collinear).
function delaunayIndices(points: readonly Vector[]): number[] | null {
    const delaunay = new Delaunay2();
    if (!delaunay.compute(points) || delaunay.getDimension() !== 2) {
        return null;
    }
    return Array.from(delaunay.getIndices());
}

// The signed area of the triangle, twice. The lattice coordinates are small
// integers, so this is exact.
function twiceSignedArea(p0: Vector, p1: Vector, p2: Vector): number {
    return (p1.get(0) - p0.get(0)) * (p2.get(1) - p0.get(1))
        - (p1.get(1) - p0.get(1)) * (p2.get(0) - p0.get(0));
}

// The brute-force containment test: the indices of every triangle of the mesh
// whose closed region contains p, using the same convention as upstream's
// PointInPolygon2 (a point on an edge belongs to the triangle).
function bruteForceContainers(mesh: PlanarMesh, p: Vector): number[] {
    const found: number[] = [];
    for (let t = 0; t < mesh.getNumTriangles(); ++t) {
        const tri = mesh.getTriangleVertices(t) as [Vector, Vector, Vector];
        const a0 = twiceSignedArea(p, tri[0], tri[1]);
        const a1 = twiceSignedArea(p, tri[1], tri[2]);
        const a2 = twiceSignedArea(p, tri[2], tri[0]);
        if ((a0 >= 0 && a1 >= 0 && a2 >= 0) || (a0 <= 0 && a1 <= 0 && a2 <= 0)) {
            found.push(t);
        }
    }
    return found;
}

// Distinct lattice points in [-6,6]^2, enough of them for a two-dimensional
// Delaunay triangulation. Integer coordinates keep the PrimalQuery2 sign
// tests exact, so the linear walk is not at the mercy of rounding.
const latticePoints = fc.uniqueArray(
    fc.tuple(fc.integer({ min: -6, max: 6 }), fc.integer({ min: -6, max: 6 })),
    { minLength: 6, maxLength: 14, selector: xy => xy[0] + ':' + xy[1] })
    .map(pairs => pairs.map(xy => v2(xy[0], xy[1])));

describe('PlanarMesh verification', () => {
    it('locates points by linear walk exactly as an exhaustive search does', () => {
        check(fc.tuple(latticePoints, fc.integer({ min: 0, max: 1023 })),
            ([points, salt]) => {
                const indices = delaunayIndices(points);
                if (indices === null) {
                    return;
                }
                const mesh = PlanarMesh.fromIndices(points, indices);
                expect(mesh.getNumTriangles()).toBe(indices.length / 3);

                // The convex hull of a Delaunay triangulation is convex, so
                // the linear walk is guaranteed to succeed for any interior
                // point from any starting triangle.
                const random = seededRandom(salt + 1);
                for (let trial = 0; trial < 12; ++trial) {
                    // A random convex combination of a random triangle's
                    // vertices lies in the hull.
                    const source = Math.floor(random() * mesh.getNumTriangles());
                    const tri = mesh.getTriangleVertices(source) as
                        [Vector, Vector, Vector];
                    let b0 = random();
                    let b1 = random();
                    let b2 = random();
                    const sum = b0 + b1 + b2;
                    if (sum < 1e-6) {
                        continue;
                    }
                    b0 /= sum;
                    b1 /= sum;
                    b2 /= sum;
                    const p = v2(
                        b0 * tri[0].get(0) + b1 * tri[1].get(0) + b2 * tri[2].get(0),
                        b0 * tri[0].get(1) + b1 * tri[1].get(1) + b2 * tri[2].get(1));

                    const containers = bruteForceContainers(mesh, p);
                    expect(containers.length).toBeGreaterThan(0);

                    // The walk finds a containing triangle from every start.
                    for (let start = 0; start < mesh.getNumTriangles(); ++start) {
                        const t = mesh.getContainingTriangle(p, start);
                        expect(containers).toContain(t);
                        expect(mesh.contains(t, p)).toBe(true);

                        // The cycle-trapping overload agrees and its visited
                        // set records the start and the reported triangle.
                        const walk = mesh.getContainingTriangleVisited(p, start);
                        expect(walk.triangle).toBe(t);
                        expect(walk.visited.has(start)).toBe(true);
                        expect(walk.visited.has(t)).toBe(true);

                        // The barycentric coordinates are a partition of
                        // unity that reproduces the query point.
                        const bary = mesh.getBarycentrics(t, p) as
                            [number, number, number];
                        expect(bary).not.toBeNull();
                        expectClose(bary[0] + bary[1] + bary[2], 1, 1e-12, 1e-12);
                        const vtx = mesh.getTriangleVertices(t) as
                            [Vector, Vector, Vector];
                        for (let j = 0; j < 2; ++j) {
                            expectClose(
                                bary[0] * vtx[0].get(j) + bary[1] * vtx[1].get(j)
                                + bary[2] * vtx[2].get(j), p.get(j), 1e-11, 1e-11);
                            expect(bary[j]).toBeGreaterThan(-1e-12);
                        }
                    }
                }
            }, 40);
    }, 30000);

    it('reports the -1 sentinel for points outside the convex hull', () => {
        check(fc.tuple(latticePoints, fc.integer({ min: 0, max: 1023 })),
            ([points, salt]) => {
                const indices = delaunayIndices(points);
                if (indices === null) {
                    return;
                }
                const mesh = PlanarMesh.fromIndices(points, indices);

                // Far outside the [-6,6]^2 lattice: no triangle can contain
                // the point, and the walk must exit through a boundary edge.
                const random = seededRandom(salt + 7);
                for (let trial = 0; trial < 8; ++trial) {
                    const angle = 2 * Math.PI * random();
                    const p = v2(100 * Math.cos(angle), 100 * Math.sin(angle));
                    expect(bruteForceContainers(mesh, p)).toEqual([]);
                    for (let start = 0; start < mesh.getNumTriangles(); ++start) {
                        expect(mesh.getContainingTriangle(p, start)).toBe(-1);
                        expect(mesh.getContainingTriangleVisited(p, start)
                            .triangle).toBe(-1);
                    }
                }
            }, 40);
    }, 30000);

    it('keeps the adjacency graph consistent with the shared edges', () => {
        check(latticePoints, points => {
            const indices = delaunayIndices(points);
            if (indices === null) {
                return;
            }
            const mesh = PlanarMesh.fromIndices(points, indices);
            const adjacencies = mesh.getAdjacencies();
            expect(adjacencies.length).toBe(3 * mesh.getNumTriangles());

            for (let t = 0; t < mesh.getNumTriangles(); ++t) {
                const ti = mesh.getTriangleIndices(t) as [number, number, number];
                // The triangles are counterclockwise (Delaunay2 output).
                const tv = mesh.getTriangleVertices(t) as [Vector, Vector, Vector];
                expect(twiceSignedArea(tv[0], tv[1], tv[2])).toBeGreaterThan(0);

                for (let i = 0; i < 3; ++i) {
                    const a = adjacencies[3 * t + i] as number;
                    if (a === -1) {
                        continue;
                    }
                    // Adjacency i of triangle t sits across the edge
                    // <V[i], V[i+1]>, and the relation is symmetric.
                    const ai = mesh.getTriangleAdjacencies(a) as
                        [number, number, number];
                    const j = ai.indexOf(t);
                    expect(j).toBeGreaterThanOrEqual(0);
                    const aIndices = mesh.getTriangleIndices(a) as
                        [number, number, number];
                    const shared = [ti[i], ti[(i + 1) % 3]].sort((x, y) => x - y);
                    const sharedA = [aIndices[j], aIndices[(j + 1) % 3]]
                        .sort((x, y) => x - y);
                    expect(sharedA).toEqual(shared);
                }
            }
        }, 60);
    }, 30000);

    it('drives IntpLinearNonuniform2 through the duck-typed mesh interface', () => {
        // The interpolator's mesh interface is satisfied by name AND arity;
        // the assignment below is the compile-time half of the check and the
        // exactness assertion is the runtime half. Binding getTriangleIndices
        // to a zero-argument accessor would blend the first triangle's
        // samples for every query point (the V29 defect), which a linear
        // function reproduced exactly cannot hide.
        check(fc.tuple(latticePoints, fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 }), fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: 0, max: 1023 })),
        ([points, a, b, c, salt]) => {
            const indices = delaunayIndices(points);
            if (indices === null) {
                return;
            }
            const planar = PlanarMesh.fromIndices(points, indices);
            const mesh: IntpLinearNonuniform2TriangleMesh = planar;
            const quadMesh: IntpQuadraticNonuniform2TriangleMesh = planar;
            expect(quadMesh.getNumTriangles()).toBe(planar.getNumTriangles());

            // f(x,y) = a*x + b*y + c is reproduced exactly by barycentric
            // interpolation over any triangulation of its samples.
            const F = points.map(p => a * p.get(0) + b * p.get(1) + c);
            const interpolator = new IntpLinearNonuniform2(mesh, F);

            const random = seededRandom(salt + 11);
            for (let trial = 0; trial < 10; ++trial) {
                const source = Math.floor(random() * planar.getNumTriangles());
                const tri = planar.getTriangleVertices(source) as
                    [Vector, Vector, Vector];
                let b0 = random();
                let b1 = random();
                let b2 = random();
                const sum = b0 + b1 + b2;
                if (sum < 1e-6) {
                    continue;
                }
                b0 /= sum;
                b1 /= sum;
                b2 /= sum;
                const p = v2(
                    b0 * tri[0].get(0) + b1 * tri[1].get(0) + b2 * tri[2].get(0),
                    b0 * tri[0].get(1) + b1 * tri[1].get(1) + b2 * tri[2].get(1));

                const result = interpolator.evaluate(p);
                expect(result.valid).toBe(true);
                expectClose(result.F, a * p.get(0) + b * p.get(1) + c,
                    1e-9, 1e-9);
            }

            // Outside the hull the interpolation reports failure.
            const outside = interpolator.evaluate(v2(1000, 1000));
            expect(outside.valid).toBe(false);
        }, 40);
    }, 30000);

    it('gives fromMesh the same triangles as fromIndices', () => {
        check(latticePoints, points => {
            const indices = delaunayIndices(points);
            if (indices === null) {
                return;
            }
            const fromIndices = PlanarMesh.fromIndices(points, indices);

            const etMesh = new ETManifoldMesh();
            for (let t = 0; t < indices.length / 3; ++t) {
                expect(etMesh.insert(indices[3 * t] as number,
                    indices[3 * t + 1] as number,
                    indices[3 * t + 2] as number)).not.toBeNull();
            }
            const fromMesh = PlanarMesh.fromMesh(points, etMesh);
            expect(fromMesh.getNumTriangles()).toBe(fromIndices.getNumTriangles());

            // The same triangles up to rotation of each index triple and to
            // the ordering of the triangles themselves.
            const rotatedKey = (mesh: PlanarMesh, t: number): string => {
                const triple = (mesh.getTriangleIndices(t) as number[]).slice();
                const min = Math.min(...triple);
                while (triple[0] !== min) {
                    triple.push(triple.shift() as number);
                }
                return triple.join(',');
            };
            const keysA: string[] = [];
            const keysB: string[] = [];
            for (let t = 0; t < fromIndices.getNumTriangles(); ++t) {
                keysA.push(rotatedKey(fromIndices, t));
                keysB.push(rotatedKey(fromMesh, t));
            }
            expect(keysB.slice().sort()).toEqual(keysA.slice().sort());

            // And the same point-location answers.
            for (let t = 0; t < fromIndices.getNumTriangles(); ++t) {
                const tri = fromIndices.getTriangleVertices(t) as
                    [Vector, Vector, Vector];
                const p = v2(
                    (tri[0].get(0) + tri[1].get(0) + tri[2].get(0)) / 3,
                    (tri[0].get(1) + tri[1].get(1) + tri[2].get(1)) / 3);
                const ta = fromIndices.getContainingTriangle(p);
                const tb = fromMesh.getContainingTriangle(p);
                expect(ta).toBeGreaterThanOrEqual(0);
                expect(tb).toBeGreaterThanOrEqual(0);
                expect(keysB[tb]).toBe(keysA[ta]);
            }
        }, 40);
    }, 30000);

    it('half-constructs on a duplicated triangle, as upstream does', () => {
        // Upstream issue #256 item 2, preserved: the first constructor
        // returns early when Insert reports a duplicate, leaving the object
        // with no vertices, no triangles and an unset query.
        const mesh = PlanarMesh.fromIndices(squareVertices,
            [0, 1, 2, 0, 2, 3, 0, 1, 2]);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(Array.from(mesh.getIndices())).toEqual([]);
        expect(Array.from(mesh.getAdjacencies())).toEqual([]);
        expect(mesh.getTriangleIndices(0)).toBeNull();
        expect(mesh.getContainingTriangle(v2(0.5, 0.25))).toBe(-1);
    });
});
