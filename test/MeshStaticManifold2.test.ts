import { describe, it, expect } from 'vitest';
import { MeshStaticManifold2, MeshStaticManifold2Vertex } from '../src/MeshStaticManifold2';

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
