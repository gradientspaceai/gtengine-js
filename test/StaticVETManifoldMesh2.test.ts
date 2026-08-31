import { describe, it, expect } from 'vitest';
import { StaticVETManifoldMesh2 } from '../src/StaticVETManifoldMesh2';

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
