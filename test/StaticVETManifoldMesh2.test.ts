import { describe, it, expect } from 'vitest';
import { StaticVETManifoldMesh2 } from '../src/StaticVETManifoldMesh2.js';
import {
    ETManifoldMesh, ETManifoldMeshTriangle
} from '../src/ETManifoldMesh.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { check, fc } from './helpers/arbitraries.js';

const invalid = StaticVETManifoldMesh2.invalid;

// Two triangles that share the edge (0,2) of the unit square with vertices
// 0=(0,0), 1=(1,0), 2=(1,1), 3=(0,1). Both triangles are counterclockwise.
const squareTriangles: [number, number, number][] = [[0, 1, 2], [0, 2, 3]];

// The four outward-facing triangles of a tetrahedron surface, which is a
// closed mesh with no boundary edges.
const tetrahedronTriangles: [number, number, number][] = [
    [0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]
];

// Build an m-by-n grid of quads, each split into two triangles. The vertex
// at grid location (i,j) has index i + (m + 1) * j.
function makeGrid(m: number, n: number): {
    numVertices: number;
    triangles: [number, number, number][];
} {
    const index = (i: number, j: number) => i + (m + 1) * j;
    const triangles: [number, number, number][] = [];
    for (let j = 0; j < n; ++j) {
        for (let i = 0; i < m; ++i) {
            const v00 = index(i, j);
            const v10 = index(i + 1, j);
            const v11 = index(i + 1, j + 1);
            const v01 = index(i, j + 1);
            triangles.push([v00, v10, v11]);
            triangles.push([v00, v11, v01]);
        }
    }
    return { numVertices: (m + 1) * (n + 1), triangles };
}

describe('StaticVETManifoldMesh2 construction', () => {
    it('rejects invalid input', () => {
        expect(() => new StaticVETManifoldMesh2(2, squareTriangles))
            .toThrow('invalid input');
        expect(() => new StaticVETManifoldMesh2(4, []))
            .toThrow('invalid input');
    });

    it('stores the triangles and the per-vertex triangle counts', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        expect(mesh.getTriangles()).toEqual(squareTriangles);
        // Vertex 0 and vertex 2 belong to both triangles; 1 and 3 to one.
        expect(mesh.getMinNumTrianglesAtVertex()).toBe(1);
        expect(mesh.getMaxNumTrianglesAtVertex()).toBe(2);
        expect(mesh.getVertices().map(v => v.getNumTAdjacents()))
            .toEqual([2, 1, 2, 1]);
    });

    it('computes the triangle-triangle adjacency', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        // mAdjacents[t][a] is the triangle across the directed edge
        // <tri[a], tri[(a+1)%3]>. Only the edge (0,2) is shared.
        expect(mesh.getAdjacents()).toEqual([
            [invalid, invalid, 1],
            [0, invalid, invalid]
        ]);
    });

    it('collects the adjacent vertices without duplicates', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        const vertices = mesh.getVertices();
        // Vertex 0 sees 1 and 2 from triangle 0 and 2 and 3 from triangle 1.
        expect(vertices[0].getVAdjacents()).toEqual([1, 2, 3]);
        expect(vertices[0].getNumVAdjacents()).toBe(3);
        expect(vertices[1].getVAdjacents()).toEqual([2, 0]);
        expect(vertices[2].getVAdjacents()).toEqual([0, 1, 3]);
        expect(vertices[3].getVAdjacents()).toEqual([0, 2]);
    });

    it('stores the outgoing-edge triples <AV,LT,RT>', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        const vertices = mesh.getVertices();
        // Vertex 0 has outgoing edges <0,1> (boundary, from triangle 0) and
        // <0,2> (interior, from triangle 1 with triangle 0 on the right).
        expect(vertices[0].getNumEAdjacents()).toBe(2);
        expect(vertices[0].getEAdjacents()).toEqual([[1, 0, invalid], [2, 1, 0]]);
        // Vertex 2 has the outgoing edge <2,0>, the direction opposite to
        // <0,2>, and the boundary edge <2,3>.
        expect(vertices[2].getEAdjacents())
            .toEqual([[0, 0, 1], [3, 1, invalid]]);
        expect(vertices[1].getEAdjacents()).toEqual([[2, 0, invalid]]);
        expect(vertices[3].getEAdjacents()).toEqual([[0, 1, invalid]]);
    });

    it('builds a closed mesh with no boundary', () => {
        const mesh = new StaticVETManifoldMesh2(4, tetrahedronTriangles);
        expect(mesh.getMinNumTrianglesAtVertex()).toBe(3);
        expect(mesh.getMaxNumTrianglesAtVertex()).toBe(3);
        for (const adjacents of mesh.getAdjacents()) {
            for (const a of adjacents) {
                expect(a).not.toBe(invalid);
            }
        }
        expect(mesh.getBoundaryPolygons(false)).toEqual([]);
    });
});

describe('StaticVETManifoldMesh2.edgeExists', () => {
    it('finds every mesh edge in both orders', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        const meshEdges: [number, number][] =
            [[0, 1], [1, 2], [2, 0], [0, 2], [2, 3], [3, 0]];
        for (const [v0, v1] of meshEdges) {
            expect(mesh.edgeExists(v0, v1)).toBe(true);
            expect(mesh.edgeExists(v1, v0)).toBe(true);
        }
    });

    it('rejects nonedges and out-of-range or degenerate queries', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        expect(mesh.edgeExists(1, 3)).toBe(false);
        expect(mesh.edgeExists(0, 0)).toBe(false);
        expect(mesh.edgeExists(0, 4)).toBe(false);
        expect(mesh.edgeExists(-1, 0)).toBe(false);
    });
});

describe('StaticVETManifoldMesh2.getAdjacentTriangles', () => {
    it('returns both triangles of an interior edge', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        // <0,2> belongs to triangle 1 and <2,0> belongs to triangle 0.
        expect(mesh.getAdjacentTriangles(0, 2))
            .toEqual({ exists: true, adj0: 1, adj1: 0 });
        // Reversing the query swaps the pair.
        expect(mesh.getAdjacentTriangles(2, 0))
            .toEqual({ exists: true, adj0: 0, adj1: 1 });
    });

    it('distinguishes the two boundary-edge cases (upstream bug fixed)', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        // Case 2: <0,1> is outgoing, <1,0> is not.
        expect(mesh.getAdjacentTriangles(0, 1))
            .toEqual({ exists: true, adj0: 0, adj1: invalid });
        // Case 3: <1,0> is not outgoing but <0,1> is. Upstream returns
        // (0, invalid) here, which is the signature it documents for case 2;
        // the port returns the documented (invalid, valid) pair.
        expect(mesh.getAdjacentTriangles(1, 0))
            .toEqual({ exists: true, adj0: invalid, adj1: 0 });

        // The same for the remaining boundary edges.
        expect(mesh.getAdjacentTriangles(2, 3))
            .toEqual({ exists: true, adj0: 1, adj1: invalid });
        expect(mesh.getAdjacentTriangles(3, 2))
            .toEqual({ exists: true, adj0: invalid, adj1: 1 });
    });

    it('returns case 4 for a nonedge', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        expect(mesh.getAdjacentTriangles(1, 3))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
        expect(mesh.getAdjacentTriangles(2, 2))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
        expect(mesh.getAdjacentTriangles(0, 99))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
    });

    it('is antisymmetric and agrees with edgeExists on a grid', () => {
        const { numVertices, triangles } = makeGrid(4, 3);
        const mesh = new StaticVETManifoldMesh2(numVertices, triangles);

        for (const tri of triangles) {
            for (let a = 0; a < 3; ++a) {
                const v0 = tri[a];
                const v1 = tri[(a + 1) % 3];
                const forward = mesh.getAdjacentTriangles(v0, v1);
                const backward = mesh.getAdjacentTriangles(v1, v0);
                expect(forward.exists).toBe(true);
                expect(backward.exists).toBe(true);
                expect(mesh.edgeExists(v0, v1)).toBe(true);
                // Reversing the query swaps the pair.
                expect(backward.adj0).toBe(forward.adj1);
                expect(backward.adj1).toBe(forward.adj0);
                // Exactly one of the documented cases 1, 2 or 3 occurs.
                expect(forward.adj0 !== invalid || forward.adj1 !== invalid)
                    .toBe(true);
            }
        }
    });

    it('agrees with the triangle-triangle adjacency array', () => {
        const { numVertices, triangles } = makeGrid(3, 3);
        const mesh = new StaticVETManifoldMesh2(numVertices, triangles);
        const adjacents = mesh.getAdjacents();

        for (let t = 0; t < triangles.length; ++t) {
            for (let a = 0; a < 3; ++a) {
                const v0 = triangles[t][a];
                const v1 = triangles[t][(a + 1) % 3];
                // <v0,v1> is outgoing in triangle t, so adj0 is t and adj1 is
                // the triangle across the edge, which is adjacents[t][a].
                const { exists, adj0, adj1 } = mesh.getAdjacentTriangles(v0, v1);
                expect(exists).toBe(true);
                expect(adj0).toBe(t);
                expect(adj1).toBe(adjacents[t][a]);
            }
        }
    });
});

describe('StaticVETManifoldMesh2.getComponents', () => {
    it('finds a single component for a connected mesh', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect(components[0].slice().sort((a, b) => a - b)).toEqual([0, 1]);
    });

    it('separates disconnected meshes', () => {
        const triangles: [number, number, number][] =
            [[0, 1, 2], [3, 4, 5], [0, 2, 6]];
        const mesh = new StaticVETManifoldMesh2(7, triangles);
        const components = mesh.getComponents()
            .map(c => c.slice().sort((a, b) => a - b));
        expect(components.length).toBe(2);
        expect(components).toContainEqual([0, 2]);
        expect(components).toContainEqual([1]);
    });

    it('partitions every triangle of a grid into one component', () => {
        const { numVertices, triangles } = makeGrid(5, 4);
        const mesh = new StaticVETManifoldMesh2(numVertices, triangles);
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect(components[0].slice().sort((a, b) => a - b))
            .toEqual(triangles.map((_, t) => t));
    });
});

describe('StaticVETManifoldMesh2.getBoundaryPolygons', () => {
    it('traverses the boundary of two triangles', () => {
        const mesh = new StaticVETManifoldMesh2(4, squareTriangles);
        const open = mesh.getBoundaryPolygons(false);
        expect(open.length).toBe(1);
        expect(open[0]).toEqual([0, 1, 2, 3]);

        const closed = mesh.getBoundaryPolygons(true);
        expect(closed.length).toBe(1);
        expect(closed[0]).toEqual([0, 1, 2, 3, 0]);
    });

    it('returns one polygon per connected component', () => {
        const triangles: [number, number, number][] =
            [[0, 1, 2], [3, 4, 5], [0, 2, 6]];
        const mesh = new StaticVETManifoldMesh2(7, triangles);
        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons.length).toBe(2);
        const sizes = polygons.map(p => p.length).sort((a, b) => a - b);
        expect(sizes).toEqual([3, 4]);
    });

    it('walks the perimeter of a grid', () => {
        const m = 4;
        const n = 3;
        const { numVertices, triangles } = makeGrid(m, n);
        const mesh = new StaticVETManifoldMesh2(numVertices, triangles);
        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons.length).toBe(1);
        // The perimeter of the grid has 2*(m+n) vertices.
        expect(polygons[0].length).toBe(2 * (m + n));
        // Every perimeter vertex is on the grid boundary and occurs once.
        expect(new Set(polygons[0]).size).toBe(polygons[0].length);
        for (const v of polygons[0]) {
            const i = v % (m + 1);
            const j = Math.floor(v / (m + 1));
            expect(i === 0 || i === m || j === 0 || j === n).toBe(true);
        }
        // Consecutive perimeter vertices are mesh edges.
        for (let k = 0; k < polygons[0].length; ++k) {
            const v0 = polygons[0][k];
            const v1 = polygons[0][(k + 1) % polygons[0].length];
            expect(mesh.edgeExists(v0, v1)).toBe(true);
        }
    });
});

describe('StaticVETManifoldMesh2 verification', () => {
    // Random subsets of a consistently oriented grid mesh. Every subset is a
    // manifold mesh with no mixed chirality, which is what the class
    // preconditions require.
    const meshArb = fc.integer({ min: 1, max: 3 }).chain((n) => {
        const grid = makeGrid(n, n);
        return fc.array(fc.boolean(), {
            minLength: grid.triangles.length,
            maxLength: grid.triangles.length
        }).map((mask) => ({
            numVertices: grid.numVertices,
            triangles: grid.triangles.filter((_, i) => mask[i])
        }));
    }).filter(({ triangles }) => triangles.length > 0);

    // The dynamic sibling mesh built from the same triangle soup.
    const dynamicOf = (triangles: [number, number, number][]): ETManifoldMesh => {
        const mesh = new ETManifoldMesh();
        for (const t of triangles) {
            expect(mesh.insert(t[0], t[1], t[2])).not.toBeNull();
        }
        return mesh;
    };

    it('agrees with ETManifoldMesh on edges and triangle adjacency', () => {
        check(meshArb, ({ numVertices, triangles }) => {
            const stat = new StaticVETManifoldMesh2(numVertices, triangles);
            const dyn = dynamicOf(triangles);

            // The same undirected edges.
            const dynEdges = new Set(dyn.getEdgeKeys().map((k) => k.mapKey()));
            const statEdges = new Set<string>();
            for (const t of triangles) {
                for (let i = 0; i < 3; ++i) {
                    const a = t[i];
                    const b = t[(i + 1) % 3];
                    statEdges.add(Math.min(a, b) + ',' + Math.max(a, b));
                    expect(stat.edgeExists(a, b)).toBe(true);
                    expect(stat.edgeExists(b, a)).toBe(true);
                }
            }
            expect(statEdges).toEqual(dynEdges);

            // The same triangle-triangle adjacency. Both classes index the
            // adjacency by the edge (V[i], V[(i+1)%3]).
            const adjacents = stat.getAdjacents();
            for (let t = 0; t < triangles.length; ++t) {
                const tri = triangles[t];
                const dynTri = dyn.getTriangle(tri[0], tri[1], tri[2]);
                expect(dynTri).not.toBeNull();
                for (let i = 0; i < 3; ++i) {
                    // The dynamic triangle may be a rotation of the input, so
                    // look the edge up by its vertices.
                    const dynAdj = (dynTri as ETManifoldMeshTriangle)
                        .getAdjacentOfEdge(tri[i], tri[(i + 1) % 3]);
                    if (adjacents[t][i] === invalid) {
                        expect(dynAdj).toBeNull();
                    } else {
                        expect(dynAdj).not.toBeNull();
                        const other = triangles[adjacents[t][i]];
                        expect(new TriangleKey(true, other[0], other[1], other[2])
                            .mapKey())
                            .toBe(new TriangleKey(true,
                                (dynAdj as ETManifoldMeshTriangle).V[0],
                                (dynAdj as ETManifoldMeshTriangle).V[1],
                                (dynAdj as ETManifoldMeshTriangle).V[2]).mapKey());
                    }
                }
            }
        }, 60);
    });

    it('the packed vertex storage decodes to the mesh adjacency', () => {
        check(meshArb, ({ numVertices, triangles }) => {
            const stat = new StaticVETManifoldMesh2(numVertices, triangles);
            const vertices = stat.getVertices();
            expect(vertices.length).toBe(numVertices);

            // Brute-force per-vertex data.
            const counts: number[] = new Array<number>(numVertices).fill(0);
            const neighbors: Set<number>[] = [];
            const outgoing: Map<number, number>[] = [];
            for (let v = 0; v < numVertices; ++v) {
                neighbors.push(new Set<number>());
                outgoing.push(new Map<number, number>());
            }
            for (let t = 0; t < triangles.length; ++t) {
                const tri = triangles[t];
                for (let i = 0; i < 3; ++i) {
                    ++counts[tri[i]];
                    neighbors[tri[i]].add(tri[(i + 1) % 3]);
                    neighbors[tri[i]].add(tri[(i + 2) % 3]);
                    // The directed edge <tri[i], tri[i+1]> belongs to t.
                    outgoing[tri[i]].set(tri[(i + 1) % 3], t);
                }
            }

            expect(stat.getMinNumTrianglesAtVertex()).toBe(Math.min(...counts));
            expect(stat.getMaxNumTrianglesAtVertex()).toBe(Math.max(...counts));

            for (let v = 0; v < numVertices; ++v) {
                const vertex = vertices[v];
                expect(vertex.getNumTAdjacents()).toBe(counts[v]);
                // One outgoing edge per incident triangle.
                expect(vertex.getNumEAdjacents()).toBe(counts[v]);
                expect(new Set(vertex.getVAdjacents())).toEqual(neighbors[v]);
                expect(vertex.getVAdjacents().length)
                    .toBe(vertex.getNumVAdjacents());
                expect(vertex.getNumVAdjacents())
                    .toBeLessThanOrEqual(2 * vertex.getNumTAdjacents());

                const triples = vertex.getEAdjacents();
                expect(triples.length).toBe(counts[v]);
                expect(new Set(triples.map((e) => e[0])).size)
                    .toBe(triples.length);
                for (let j = 0; j < triples.length; ++j) {
                    const [av, lt, rt] = triples[j];
                    // LT is the triangle containing the directed edge <v,av>.
                    expect(lt).toBe(outgoing[v].get(av));
                    // RT is the triangle containing the reversed directed
                    // edge, or 'invalid' on a boundary edge.
                    const reverse = outgoing[av].get(v);
                    expect(rt).toBe(reverse === undefined ? invalid : reverse);
                    expect(vertex.getEAdjacent(j)).toEqual([av, lt, rt]);
                }
            }
        }, 60);
    });

    it('getAdjacentTriangles reports the documented four cases', () => {
        check(meshArb, ({ numVertices, triangles }) => {
            const stat = new StaticVETManifoldMesh2(numVertices, triangles);
            // Brute force: the triangle containing each directed edge.
            const owner = new Map<string, number>();
            for (let t = 0; t < triangles.length; ++t) {
                for (let i = 0; i < 3; ++i) {
                    owner.set(triangles[t][i] + ',' +
                        triangles[t][(i + 1) % 3], t);
                }
            }

            for (let v0 = 0; v0 < numVertices; ++v0) {
                for (let v1 = 0; v1 < numVertices; ++v1) {
                    const forward = owner.get(v0 + ',' + v1);
                    const backward = owner.get(v1 + ',' + v0);
                    const result = stat.getAdjacentTriangles(v0, v1);
                    if (v0 === v1 || (forward === undefined &&
                        backward === undefined)) {
                        // Case 4.
                        expect(result).toEqual({
                            exists: false, adj0: invalid, adj1: invalid
                        });
                        continue;
                    }
                    expect(result.exists).toBe(true);
                    expect(result.adj0)
                        .toBe(forward === undefined ? invalid : forward);
                    expect(result.adj1)
                        .toBe(backward === undefined ? invalid : backward);
                    // The query is antisymmetric in its arguments.
                    const swapped = stat.getAdjacentTriangles(v1, v0);
                    expect(swapped.adj0).toBe(result.adj1);
                    expect(swapped.adj1).toBe(result.adj0);
                    expect(stat.edgeExists(v0, v1)).toBe(true);
                }
            }
        }, 25);
    });

    it('agrees with ETManifoldMesh on the connected components', () => {
        check(meshArb, ({ numVertices, triangles }) => {
            const stat = new StaticVETManifoldMesh2(numVertices, triangles);
            const dyn = dynamicOf(triangles);

            const key = (t: number[]): string =>
                new TriangleKey(true, t[0], t[1], t[2]).mapKey();
            const statComponents = stat.getComponents()
                .map((c) => c.map((t) => key(triangles[t])).sort().join(';'))
                .sort();
            const dynComponents = dyn.getComponents()
                .map((c) => c.map((t) => key(t.V)).sort().join(';'))
                .sort();
            expect(statComponents).toEqual(dynComponents);

            // The components partition the triangle indices.
            const flat = stat.getComponents().flat();
            expect(flat.length).toBe(triangles.length);
            expect(new Set(flat).size).toBe(triangles.length);
        }, 60);
    });

    it('agrees with ETManifoldMesh on the boundary polygons', () => {
        check(meshArb, ({ numVertices, triangles }) => {
            const stat = new StaticVETManifoldMesh2(numVertices, triangles);
            const dyn = dynamicOf(triangles);

            // Normalize a polygon to its lexicographically smallest rotation
            // so that the two traversals can be compared regardless of where
            // they start.
            const canonical = (polygon: number[]): string => {
                let best: string | null = null;
                for (let s = 0; s < polygon.length; ++s) {
                    const rotated = polygon.slice(s).concat(polygon.slice(0, s));
                    const text = rotated.join(',');
                    if (best === null || text < best) {
                        best = text;
                    }
                }
                return best as string;
            };

            let statPolygons: number[][];
            let dynPolygons: number[][];
            try {
                statPolygons = stat.getBoundaryPolygons(false);
            } catch (e) {
                // A pinched (bow-tie) vertex makes both traversals fail.
                expect(() => dyn.getBoundaryPolygons(false)).toThrow();
                return;
            }
            dynPolygons = dyn.getBoundaryPolygons(false);
            expect(statPolygons.map(canonical).sort())
                .toEqual(dynPolygons.map(canonical).sort());

            // Each boundary edge is traversed exactly once.
            const boundary: string[] = [];
            const adjacents = stat.getAdjacents();
            for (let t = 0; t < triangles.length; ++t) {
                for (let i = 0; i < 3; ++i) {
                    if (adjacents[t][i] === invalid) {
                        boundary.push(triangles[t][i] + ',' +
                            triangles[t][(i + 1) % 3]);
                    }
                }
            }
            const covered: string[] = [];
            for (const polygon of statPolygons) {
                for (let i = 0; i < polygon.length; ++i) {
                    covered.push(polygon[i] + ',' +
                        polygon[(i + 1) % polygon.length]);
                }
            }
            expect(covered.sort()).toEqual(boundary.sort());

            // duplicateEndpoints repeats the first vertex at the end.
            const closed = stat.getBoundaryPolygons(true);
            expect(closed.length).toBe(statPolygons.length);
            for (let i = 0; i < closed.length; ++i) {
                expect(closed[i][closed[i].length - 1]).toBe(closed[i][0]);
                expect(closed[i].slice(0, -1)).toEqual(statPolygons[i]);
            }
        }, 60);
    });

    it('the triangle order does not change the mesh relations', () => {
        check(fc.tuple(meshArb, fc.array(fc.nat(), { maxLength: 16 })),
            ([{ numVertices, triangles }, shuffle]) => {
                const permuted = triangles.slice();
                for (let k = 0; k < shuffle.length; ++k) {
                    const i = shuffle[k] % permuted.length;
                    const j = (shuffle[k] * 5 + k) % permuted.length;
                    const t = permuted[i];
                    permuted[i] = permuted[j];
                    permuted[j] = t;
                }
                const a = new StaticVETManifoldMesh2(numVertices, triangles);
                const b = new StaticVETManifoldMesh2(numVertices, permuted);

                expect(b.getMinNumTrianglesAtVertex())
                    .toBe(a.getMinNumTrianglesAtVertex());
                expect(b.getMaxNumTrianglesAtVertex())
                    .toBe(a.getMaxNumTrianglesAtVertex());
                for (let v0 = 0; v0 < numVertices; ++v0) {
                    for (let v1 = v0 + 1; v1 < numVertices; ++v1) {
                        expect(b.edgeExists(v0, v1)).toBe(a.edgeExists(v0, v1));
                        const ra = a.getAdjacentTriangles(v0, v1);
                        const rb = b.getAdjacentTriangles(v0, v1);
                        expect(rb.exists).toBe(ra.exists);
                        // The adjacency values are triangle indices, so
                        // compare the triangles they name.
                        const nameA = (i: number): string => i === invalid
                            ? 'none' : triangles[i].join(',');
                        const nameB = (i: number): string => i === invalid
                            ? 'none' : permuted[i].join(',');
                        expect(nameB(rb.adj0)).toBe(nameA(ra.adj0));
                        expect(nameB(rb.adj1)).toBe(nameA(ra.adj1));
                    }
                }
            }, 25);
    });

    it('the numThreads argument does not change the result', () => {
        check(fc.tuple(meshArb, fc.integer({ min: 0, max: 8 })),
            ([{ numVertices, triangles }, numThreads]) => {
                const a = new StaticVETManifoldMesh2(numVertices, triangles, 0);
                const b = new StaticVETManifoldMesh2(numVertices, triangles,
                    numThreads);
                expect(b.getAdjacents()).toEqual(a.getAdjacents());
                expect(b.getVertices().map((v) => v.getEAdjacents()))
                    .toEqual(a.getVertices().map((v) => v.getEAdjacents()));
            }, 40);
    });

    it('the constructor copies the caller triangle array', () => {
        const triangles: [number, number, number][] = [[0, 1, 2], [0, 2, 3]];
        const mesh = new StaticVETManifoldMesh2(4, triangles);
        triangles[0][0] = 99;
        expect(mesh.getTriangles()[0]).toEqual([0, 1, 2]);
    });
});
