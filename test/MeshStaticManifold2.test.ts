import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { MeshStaticManifold2, MeshStaticManifold2Vertex } from '../src/MeshStaticManifold2.js';

const invalid = MeshStaticManifold2.invalid;

// Two counterclockwise triangles that tile the unit square:
//   3---2
//   | \ |
//   0---1
// with the shared diagonal (0,2).
const squareTriangles = [[0, 1, 2], [0, 2, 3]];

// The closed surface of a tetrahedron: four triangles, no boundary edges.
const tetrahedronSurface = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

// Two disjoint copies of the square mesh.
const twoSquares = [[0, 1, 2], [0, 2, 3], [4, 5, 6], [4, 6, 7]];

describe('MeshStaticManifold2Vertex', () => {
    it('starts empty and records inserted adjacency 4-tuples', () => {
        const vertex = new MeshStaticManifold2Vertex();
        expect(vertex.getNumAdjacents()).toBe(0);
        vertex.insert(1, 2, 7, 0);
        vertex.insert(3, 4, 9, 2);
        expect(vertex.getNumAdjacents()).toBe(2);
        expect(vertex.getAdjacents()).toEqual([[1, 2, 7, 0], [3, 4, 9, 2]]);
    });
});

describe('MeshStaticManifold2 construction', () => {
    it('rejects fewer than three vertices or no triangles', () => {
        expect(() => new MeshStaticManifold2(2, squareTriangles)).toThrow('Invalid input.');
        expect(() => new MeshStaticManifold2(4, [])).toThrow('Invalid input.');
    });

    it('copies the triangles rather than aliasing the input', () => {
        const input = [[0, 1, 2], [0, 2, 3]];
        const mesh = new MeshStaticManifold2(4, input);
        input[0][0] = 99;
        expect(mesh.getTriangles()[0]).toEqual([0, 1, 2]);
    });

    it('accepts and ignores the upstream numThreads argument', () => {
        const single = new MeshStaticManifold2(4, squareTriangles, 1);
        const multi = new MeshStaticManifold2(4, squareTriangles, 4);
        expect(multi.getAdjacents()).toEqual(single.getAdjacents());
    });

    it('reports the triangle counts at the vertices', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        // v0 and v2 are on the diagonal (2 triangles), v1 and v3 are corners.
        expect(mesh.getMinNumTrianglesAtVertex()).toBe(1);
        expect(mesh.getMaxNumTrianglesAtVertex()).toBe(2);

        const closed = new MeshStaticManifold2(4, tetrahedronSurface);
        expect(closed.getMinNumTrianglesAtVertex()).toBe(3);
        expect(closed.getMaxNumTrianglesAtVertex()).toBe(3);
    });
});

describe('MeshStaticManifold2 adjacency', () => {
    it('computes the triangle adjacency of the square', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        // mAdjacents[t][i] is the triangle across the edge opposite the i-th
        // vertex of triangle t. The only shared edge is the diagonal (0,2),
        // opposite vertex 1 of triangle 0 and opposite vertex 3 of triangle 1.
        expect(mesh.getAdjacents()).toEqual([
            [invalid, 1, invalid],
            [invalid, invalid, 0]
        ]);
    });

    it('stores the vertex adjacency 4-tuples', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        const vertices = mesh.getVertices();
        // Vertex 0 belongs to both triangles. The opposite edges (1,2) and
        // (2,3) are both on the boundary.
        expect(vertices[0].getAdjacents()).toEqual([
            [1, 2, 0, invalid],
            [2, 3, 1, invalid]
        ]);
        // Vertex 3 belongs to triangle 1 only; the opposite edge (0,2) is the
        // diagonal, shared with triangle 0.
        expect(vertices[3].getAdjacents()).toEqual([[0, 2, 1, 0]]);
    });

    it('makes every edge of a closed surface interior', () => {
        const mesh = new MeshStaticManifold2(4, tetrahedronSurface);
        for (const adj of mesh.getAdjacents()) {
            for (const a of adj) {
                expect(a).not.toBe(invalid);
            }
        }
        // Adjacency is symmetric: if s is a neighbor of t, then t is a
        // neighbor of s.
        const adjacents = mesh.getAdjacents();
        for (let t = 0; t < adjacents.length; ++t) {
            for (const s of adjacents[t]) {
                expect(adjacents[s]).toContain(t);
            }
        }
    });
});

describe('MeshStaticManifold2.edgeExists', () => {
    it('finds the edges of the square and rejects the others', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        for (const [a, b] of [[0, 1], [1, 2], [2, 0], [2, 3], [3, 0]]) {
            expect(mesh.edgeExists(a, b)).toBe(true);
            expect(mesh.edgeExists(b, a)).toBe(true);
        }
        // The other diagonal is not an edge.
        expect(mesh.edgeExists(1, 3)).toBe(false);
        expect(mesh.edgeExists(3, 1)).toBe(false);
    });

    it('rejects degenerate and out-of-range queries', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        expect(mesh.edgeExists(0, 0)).toBe(false);
        expect(mesh.edgeExists(0, 4)).toBe(false);
        expect(mesh.edgeExists(4, 0)).toBe(false);
        expect(mesh.edgeExists(0, invalid)).toBe(false);
    });
});

describe('MeshStaticManifold2.getAdjacentTriangles', () => {
    it('returns both triangles of an interior edge (documented case 1)', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        // Triangle 1 = <0,2,3> contains the directed edge <0,2>, so it is the
        // L-triangle; triangle 0 = <0,1,2> contains <2,0>.
        expect(mesh.getAdjacentTriangles(0, 2)).toEqual({ adj0: 1, adj1: 0, exists: true });
        expect(mesh.getAdjacentTriangles(2, 0)).toEqual({ adj0: 0, adj1: 1, exists: true });
    });

    it('returns (valid, invalid) and (invalid, valid) for a boundary edge', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        // Triangle 0 = <0,1,2> contains the directed edge <0,1>; <1,0> does
        // not occur.
        expect(mesh.getAdjacentTriangles(0, 1)).toEqual(
            { adj0: 0, adj1: invalid, exists: true });
        expect(mesh.getAdjacentTriangles(1, 0)).toEqual(
            { adj0: invalid, adj1: 0, exists: true });
    });

    it('returns (invalid, invalid) and false for a missing edge (case 4)', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        expect(mesh.getAdjacentTriangles(1, 3)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
        expect(mesh.getAdjacentTriangles(0, 0)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
        expect(mesh.getAdjacentTriangles(0, 7)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
    });

    it('agrees with edgeExists and with the triangle adjacency array', () => {
        const meshes = [
            new MeshStaticManifold2(4, squareTriangles),
            new MeshStaticManifold2(4, tetrahedronSurface),
            new MeshStaticManifold2(8, twoSquares)
        ];
        for (const mesh of meshes) {
            const numVertices = mesh.getVertices().length;
            for (let v0 = 0; v0 < numVertices; ++v0) {
                for (let v1 = 0; v1 < numVertices; ++v1) {
                    const result = mesh.getAdjacentTriangles(v0, v1);
                    expect(result.exists).toBe(mesh.edgeExists(v0, v1));
                    // The two orders of an undirected edge swap the results.
                    const swapped = mesh.getAdjacentTriangles(v1, v0);
                    expect(swapped.adj0).toBe(result.adj1);
                    expect(swapped.adj1).toBe(result.adj0);
                    // A reported triangle really does contain the edge.
                    for (const t of [result.adj0, result.adj1]) {
                        if (t !== invalid) {
                            expect(mesh.getTriangles()[t]).toContain(v0);
                            expect(mesh.getTriangles()[t]).toContain(v1);
                        }
                    }
                }
            }
        }
    });

    it('returns the same pair of triangles as the adjacency array', () => {
        const mesh = new MeshStaticManifold2(4, tetrahedronSurface);
        const triangles = mesh.getTriangles();
        const adjacents = mesh.getAdjacents();
        for (let t = 0; t < triangles.length; ++t) {
            for (let i = 0; i < 3; ++i) {
                // The edge opposite vertex i.
                const a = triangles[t][(i + 1) % 3];
                const b = triangles[t][(i + 2) % 3];
                const result = mesh.getAdjacentTriangles(a, b);
                const pair = [result.adj0, result.adj1].sort((x, y) => x - y);
                const expected = [t, adjacents[t][i]].sort((x, y) => x - y);
                expect(pair).toEqual(expected);
            }
        }
    });
});

describe('MeshStaticManifold2.getComponents', () => {
    it('reports one component for a connected mesh', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        // The depth-first search finishes the deepest triangle first.
        expect(components[0]).toEqual([1, 0]);
    });

    it('reports one component for the closed surface', () => {
        const mesh = new MeshStaticManifold2(4, tetrahedronSurface);
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect([...components[0]].sort((a, b) => a - b)).toEqual([0, 1, 2, 3]);
    });

    it('separates disconnected components', () => {
        const mesh = new MeshStaticManifold2(8, twoSquares);
        const components = mesh.getComponents();
        expect(components.length).toBe(2);
        expect(components[0]).toEqual([1, 0]);
        expect(components[1]).toEqual([3, 2]);
    });

    it('covers every triangle exactly once', () => {
        const mesh = new MeshStaticManifold2(8, twoSquares);
        const all = mesh.getComponents().flat().sort((a, b) => a - b);
        expect(all).toEqual([0, 1, 2, 3]);
    });

    it('handles a single triangle', () => {
        const mesh = new MeshStaticManifold2(3, [[0, 1, 2]]);
        expect(mesh.getComponents()).toEqual([[0]]);
    });
});

describe('MeshStaticManifold2.getBoundaryPolygons', () => {
    it('walks the boundary of the square', () => {
        const mesh = new MeshStaticManifold2(4, squareTriangles);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2, 3]]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([[0, 1, 2, 3, 0]]);
    });

    it('walks the boundary of a single triangle', () => {
        const mesh = new MeshStaticManifold2(3, [[0, 1, 2]]);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2]]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([[0, 1, 2, 0]]);
    });

    it('finds no boundary on a closed surface', () => {
        const mesh = new MeshStaticManifold2(4, tetrahedronSurface);
        expect(mesh.getBoundaryPolygons(false)).toEqual([]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([]);
    });

    it('finds one polygon per component', () => {
        const mesh = new MeshStaticManifold2(8, twoSquares);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2, 3], [4, 5, 6, 7]]);
    });

    it('produces polygons whose consecutive pairs are boundary edges', () => {
        // A triangle fan around the center vertex 0 with boundary
        // 1 -> 2 -> 3 -> 4 -> 1.
        const fan = [[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]];
        const mesh = new MeshStaticManifold2(5, fan);
        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons.length).toBe(1);
        expect(polygons[0].length).toBe(4);
        const polygon = polygons[0];
        for (let i = 0; i < polygon.length; ++i) {
            const a = polygon[i];
            const b = polygon[(i + 1) % polygon.length];
            expect(mesh.edgeExists(a, b)).toBe(true);
            // A boundary edge belongs to exactly one triangle.
            const result = mesh.getAdjacentTriangles(a, b);
            expect(result.adj0 === invalid || result.adj1 === invalid).toBe(true);
        }
        expect([...polygon].sort((x, y) => x - y)).toEqual([1, 2, 3, 4]);
    });

    it('finds the outer and inner boundaries of an annulus', () => {
        // A square annulus: outer vertices 0..3, inner vertices 4..7, eight
        // triangles forming the ring.
        //   3-------2
        //   |  7---6 |
        //   |  |   | |
        //   |  4---5 |
        //   0-------1
        const ring = [
            [0, 1, 5], [0, 5, 4],
            [1, 2, 6], [1, 6, 5],
            [2, 3, 7], [2, 7, 6],
            [3, 0, 4], [3, 4, 7]
        ];
        const mesh = new MeshStaticManifold2(8, ring);
        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons.length).toBe(2);
        const sorted = polygons.map(p => [...p].sort((a, b) => a - b));
        expect(sorted).toContainEqual([0, 1, 2, 3]);
        expect(sorted).toContainEqual([4, 5, 6, 7]);
        // Every polygon edge is a boundary edge of the mesh.
        for (const polygon of polygons) {
            for (let i = 0; i < polygon.length; ++i) {
                const a = polygon[i];
                const b = polygon[(i + 1) % polygon.length];
                const result = mesh.getAdjacentTriangles(a, b);
                expect(result.exists).toBe(true);
                expect(result.adj0 === invalid || result.adj1 === invalid).toBe(true);
            }
        }
        // The mesh is connected even though it has two boundary polygons.
        expect(mesh.getComponents().length).toBe(1);
    });
});

describe('MeshStaticManifold2 verification', () => {
    // ---- random manifold meshes ------------------------------------------
    type Mesh = { numVertices: number; triangles: number[][] };

    /** A CCW-triangulated nx-by-ny grid of quads. */
    function gridMesh(nx: number, ny: number): Mesh {
        const idx = (i: number, j: number) => i + (nx + 1) * j;
        const triangles: number[][] = [];
        for (let j = 0; j < ny; ++j) {
            for (let i = 0; i < nx; ++i) {
                const v00 = idx(i, j), v10 = idx(i + 1, j);
                const v11 = idx(i + 1, j + 1), v01 = idx(i, j + 1);
                triangles.push([v00, v10, v11]);
                triangles.push([v00, v11, v01]);
            }
        }
        return { numVertices: (nx + 1) * (ny + 1), triangles };
    }

    /** The surface of a tetrahedron: closed, every edge interior. */
    const closedTetrahedron: Mesh = {
        numVertices: 4,
        triangles: [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]]
    };

    /** The surface of an octahedron: closed, valence-4 vertices. */
    const closedOctahedron: Mesh = {
        numVertices: 6,
        triangles: [
            [0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1],
            [5, 2, 1], [5, 3, 2], [5, 4, 3], [5, 1, 4]
        ]
    };

    /** Two disjoint copies of a mesh, to exercise multiple components. */
    function disjointUnion(a: Mesh, b: Mesh): Mesh {
        return {
            numVertices: a.numVertices + b.numVertices,
            triangles: [...a.triangles,
            ...b.triangles.map(t => t.map(v => v + a.numVertices))]
        };
    }

    /**
     * Relabel the vertices, reorder the triangles and cyclically rotate each
     * triangle's vertices. All three preserve the mesh and its orientation,
     * so the derived mesh is still a valid input; the port must handle any of
     * them.
     */
    function shuffle(mesh: Mesh, perm: number[], order: number[], rots: number[]): Mesh {
        const triangles = order.map((t, k) => {
            const tri = mesh.triangles[t]!;
            const r = rots[k]! % 3;
            return [perm[tri[r]!]!, perm[tri[(r + 1) % 3]!]!, perm[tri[(r + 2) % 3]!]!];
        });
        return { numVertices: mesh.numVertices, triangles };
    }

    const baseMeshArb: fc.Arbitrary<Mesh> = fc.oneof(
        fc.tuple(fc.integer({ min: 1, max: 3 }), fc.integer({ min: 1, max: 3 }))
            .map(([nx, ny]) => gridMesh(nx, ny)),
        fc.constant(closedTetrahedron),
        fc.constant(closedOctahedron),
        fc.constant(disjointUnion(gridMesh(1, 1), gridMesh(2, 1))),
        fc.constant(disjointUnion(closedTetrahedron, gridMesh(2, 2))));

    const meshArb: fc.Arbitrary<Mesh> = baseMeshArb.chain(mesh =>
        fc.tuple(
            fc.shuffledSubarray([...Array(mesh.numVertices).keys()],
                { minLength: mesh.numVertices, maxLength: mesh.numVertices }),
            fc.shuffledSubarray([...Array(mesh.triangles.length).keys()],
                { minLength: mesh.triangles.length, maxLength: mesh.triangles.length }),
            fc.array(fc.integer({ min: 0, max: 2 }),
                { minLength: mesh.triangles.length, maxLength: mesh.triangles.length }))
            .map(([perm, order, rots]) => shuffle(mesh, perm, order, rots)));

    // ---- independent references -------------------------------------------
    const ekey = (a: number, b: number) => `${a},${b}`;

    /** Map from each directed edge to the triangle that contains it. */
    function directedEdgeMap(triangles: number[][]): Map<string, number> {
        const map = new Map<string, number>();
        triangles.forEach((tri, t) => {
            for (let i = 0; i < 3; ++i) {
                map.set(ekey(tri[i]!, tri[(i + 1) % 3]!), t);
            }
        });
        return map;
    }

    /** Expected triangle adjacency: across the edge opposite tri[i]. */
    function expectedAdjacents(triangles: number[][]): number[][] {
        const map = directedEdgeMap(triangles);
        return triangles.map(tri => {
            const row: number[] = [];
            for (let i = 0; i < 3; ++i) {
                const a = tri[(i + 1) % 3]!, b = tri[(i + 2) % 3]!;
                // The neighbour across (a,b) contains the reversed edge.
                const s = map.get(ekey(b, a));
                row.push(s === undefined ? invalid : s);
            }
            return row;
        });
    }

    // ---- properties -------------------------------------------------------
    it('the triangle adjacency matches a brute-force directed-edge map', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            expect(m.getAdjacents()).toEqual(expectedAdjacents(mesh.triangles));
            return true;
        });
    });

    it('the adjacency relation is symmetric and refers to the shared edge', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const adj = m.getAdjacents();
            const tris = m.getTriangles();
            for (let t = 0; t < adj.length; ++t) {
                for (let i = 0; i < 3; ++i) {
                    const s = adj[t]![i]!;
                    if (s === invalid) { continue; }
                    expect(s).not.toBe(t);
                    // The shared edge is opposite tris[t][i] and opposite
                    // tris[s][j] for exactly one j, and adj[s][j] === t.
                    const edge = [tris[t]![(i + 1) % 3]!, tris[t]![(i + 2) % 3]!];
                    let matches = 0;
                    for (let j = 0; j < 3; ++j) {
                        if (adj[s]![j] === t) {
                            const other = [tris[s]![(j + 1) % 3]!, tris[s]![(j + 2) % 3]!];
                            expect(new Set(other)).toEqual(new Set(edge));
                            ++matches;
                        }
                    }
                    expect(matches).toBe(1);
                }
            }
            return true;
        });
    });

    it('the vertex 4-tuples list the incident triangles in index order', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const adj = expectedAdjacents(mesh.triangles);
            for (let v = 0; v < mesh.numVertices; ++v) {
                const expected: number[][] = [];
                mesh.triangles.forEach((tri, t) => {
                    const p = tri.indexOf(v);
                    if (p >= 0) {
                        expected.push([tri[(p + 1) % 3]!, tri[(p + 2) % 3]!, t, adj[t]![p]!]);
                    }
                });
                expect(m.getVertices()[v]!.getAdjacents()).toEqual(expected);
                expect(m.getVertices()[v]!.getNumAdjacents()).toBe(expected.length);
            }
            return true;
        });
    });

    it('the min and max triangle counts are the extremes over all vertices', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const counts = new Array<number>(mesh.numVertices).fill(0);
            for (const tri of mesh.triangles) { for (const v of tri) { ++counts[v]!; } }
            return m.getMinNumTrianglesAtVertex() === Math.min(...counts)
                && m.getMaxNumTrianglesAtVertex() === Math.max(...counts);
        });
    });

    it('edgeExists is exactly membership in the undirected edge set', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const undirected = new Set<string>();
            for (const tri of mesh.triangles) {
                for (let i = 0; i < 3; ++i) {
                    const a = tri[i]!, b = tri[(i + 1) % 3]!;
                    undirected.add(ekey(Math.min(a, b), Math.max(a, b)));
                }
            }
            for (let a = 0; a < mesh.numVertices; ++a) {
                for (let b = 0; b < mesh.numVertices; ++b) {
                    const expected = a !== b
                        && undirected.has(ekey(Math.min(a, b), Math.max(a, b)));
                    expect(m.edgeExists(a, b)).toBe(expected);
                }
            }
            return true;
        });
    });

    it('getAdjacentTriangles returns the L- and R-triangles of the directed edge', () => {
        // This pins the port's fix of the upstream bug: adj0 is the triangle
        // containing <v0,v1> and adj1 the triangle containing <v1,v0>.
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const map = directedEdgeMap(mesh.triangles);
            for (let a = 0; a < mesh.numVertices; ++a) {
                for (let b = 0; b < mesh.numVertices; ++b) {
                    const r = m.getAdjacentTriangles(a, b);
                    if (a === b) {
                        expect(r).toEqual({ adj0: invalid, adj1: invalid, exists: false });
                        continue;
                    }
                    const l = map.get(ekey(a, b));
                    const rr = map.get(ekey(b, a));
                    expect(r.adj0).toBe(l === undefined ? invalid : l);
                    expect(r.adj1).toBe(rr === undefined ? invalid : rr);
                    expect(r.exists).toBe(l !== undefined || rr !== undefined);
                    expect(r.exists).toBe(m.edgeExists(a, b));
                    // Swapping the arguments swaps the two results.
                    const s = m.getAdjacentTriangles(b, a);
                    expect(s.adj0).toBe(r.adj1);
                    expect(s.adj1).toBe(r.adj0);
                }
            }
            return true;
        });
    });

    it('all four documented cases of getAdjacentTriangles occur', () => {
        // Case 3 (invalid, valid) is unreachable in the upstream code; it is
        // reachable here, which is the point of the fix.
        const m = new MeshStaticManifold2(4, [[0, 1, 2], [0, 2, 3]]);
        // Case 1: the interior edge (0,2).
        const c1 = m.getAdjacentTriangles(0, 2);
        expect(c1.exists).toBe(true);
        expect(c1.adj0).toBe(1);
        expect(c1.adj1).toBe(0);
        // Case 2: <0,1> is directed, <1,0> is not.
        const c2 = m.getAdjacentTriangles(0, 1);
        expect(c2).toEqual({ adj0: 0, adj1: invalid, exists: true });
        // Case 3: <1,0> is not directed but <0,1> is, queried the other way.
        const c3 = m.getAdjacentTriangles(1, 0);
        expect(c3).toEqual({ adj0: invalid, adj1: 0, exists: true });
        // Case 4: the missing diagonal.
        expect(m.getAdjacentTriangles(1, 3))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
    });

    it('rejects negative vertex indices instead of throwing (regression)', () => {
        // Upstream takes size_t parameters, so a negative index becomes a
        // huge unsigned value and the 'v < numVertices' guard rejects it.
        // Before the fix the port indexed the vertex array with the negative
        // value and threw a TypeError.
        const m = new MeshStaticManifold2(4, [[0, 1, 2], [0, 2, 3]]);
        expect(m.edgeExists(-1, 0)).toBe(false);
        expect(m.edgeExists(0, -1)).toBe(false);
        expect(m.edgeExists(-2, -1)).toBe(false);
        expect(m.getAdjacentTriangles(-1, 0))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
        expect(m.getAdjacentTriangles(0, -1))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
    });

    it('getComponents partitions the triangles into adjacency components', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const components = m.getComponents();
            // Every triangle appears exactly once.
            const all = components.flat();
            expect(all.length).toBe(mesh.triangles.length);
            expect(new Set(all).size).toBe(mesh.triangles.length);

            // Independent computation with union-find over the adjacency.
            const parent = [...Array(mesh.triangles.length).keys()];
            const find = (x: number): number =>
                parent[x] === x ? x : (parent[x] = find(parent[x]!));
            const adj = m.getAdjacents();
            for (let t = 0; t < adj.length; ++t) {
                for (const s of adj[t]!) {
                    if (s !== invalid) { parent[find(t)] = find(s); }
                }
            }
            const expectedGroups = new Map<number, Set<number>>();
            for (let t = 0; t < mesh.triangles.length; ++t) {
                const r = find(t);
                if (!expectedGroups.has(r)) { expectedGroups.set(r, new Set()); }
                expectedGroups.get(r)!.add(t);
            }
            expect(components.length).toBe(expectedGroups.size);
            const got = components.map(c => new Set(c));
            for (const group of expectedGroups.values()) {
                expect(got.some(g => g.size === group.size
                    && [...group].every(x => g.has(x)))).toBe(true);
            }
            return true;
        });
    });

    it('getBoundaryPolygons covers every boundary edge exactly once', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            const tris = m.getTriangles();
            const adj = m.getAdjacents();
            // The boundary directed edges, computed independently.
            const boundary = new Set<string>();
            for (let t = 0; t < tris.length; ++t) {
                for (let a = 0; a < 3; ++a) {
                    if (adj[t]![a] === invalid) {
                        boundary.add(ekey(tris[t]![(a + 1) % 3]!, tris[t]![(a + 2) % 3]!));
                    }
                }
            }

            const polygons = m.getBoundaryPolygons(false);
            const walked: string[] = [];
            for (const p of polygons) {
                expect(p.length).toBeGreaterThanOrEqual(3);
                for (let k = 0; k < p.length; ++k) {
                    walked.push(ekey(p[k]!, p[(k + 1) % p.length]!));
                }
            }
            expect(walked.length).toBe(boundary.size);
            expect(new Set(walked)).toEqual(boundary);

            // duplicateEndpoints only appends the starting vertex.
            const closed = m.getBoundaryPolygons(true);
            expect(closed.length).toBe(polygons.length);
            closed.forEach((p, i) => {
                expect(p.length).toBe(polygons[i]!.length + 1);
                expect(p.slice(0, -1)).toEqual(polygons[i]);
                expect(p[p.length - 1]).toBe(p[0]);
            });
            return true;
        });
    });

    it('a closed surface has no boundary polygons', () => {
        for (const mesh of [closedTetrahedron, closedOctahedron]) {
            const m = new MeshStaticManifold2(mesh.numVertices, mesh.triangles);
            expect(m.getBoundaryPolygons(false)).toEqual([]);
            expect(m.getAdjacents().flat().every(a => a !== invalid)).toBe(true);
        }
    });

    it('the constructor copies the triangles and rejects bad sizes', () => {
        check(meshArb, mesh => {
            const input = mesh.triangles.map(t => [...t]);
            const m = new MeshStaticManifold2(mesh.numVertices, input);
            input[0]![0] = 12345;
            expect(m.getTriangles()[0]).toEqual(mesh.triangles[0]);
            return true;
        });
        check(fc.integer({ min: -3, max: 2 }), n => {
            expect(() => new MeshStaticManifold2(n, [[0, 1, 2]])).toThrow();
        });
        expect(() => new MeshStaticManifold2(3, [])).toThrow();
    });
});
