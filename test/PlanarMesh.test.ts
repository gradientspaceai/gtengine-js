import { describe, it, expect } from 'vitest';
import { PlanarMesh } from '../src/PlanarMesh';
import { ETManifoldMesh } from '../src/ETManifoldMesh';
import { Vector } from '../src/Vector';

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
