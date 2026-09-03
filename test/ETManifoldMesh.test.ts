import { describe, it, expect } from 'vitest';
import {
    ETManifoldMesh,
    ETManifoldMeshEdge,
    ETManifoldMeshTriangle
} from '../src/ETManifoldMesh.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { TriangleKey } from '../src/TriangleKey.js';

// The ordered triangle keys of the mesh (or of a component), as arrays of the
// three vertex indices.
function keys(triangles: ETManifoldMeshTriangle[]): number[][] {
    return triangles.map(t => new TriangleKey(true, t.V[0], t.V[1], t.V[2]).V.slice());
}

// The unordered vertex triples, sorted, used where the winding order is not
// the subject of the test.
function unorderedTriples(mesh: ETManifoldMesh): number[][] {
    const triples = mesh.getTriangles().map(
        t => new TriangleKey(false, t.V[0], t.V[1], t.V[2]).V.slice());
    triples.sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]) || (a[2] - b[2]));
    return triples;
}

function edgePairs(mesh: ETManifoldMesh): number[][] {
    return mesh.getEdges().map(e => new EdgeKey(false, e.V[0], e.V[1]).V.slice());
}

// Verify the internal consistency of the mesh: every edge references only
// triangles that are in the mesh, every triangle references the edges of its
// three directed edges, and the triangle-triangle adjacency agrees with the
// edge-triangle adjacency.
function checkInvariants(mesh: ETManifoldMesh): void {
    const triangles = new Set(mesh.getTriangles());
    const edges = new Set(mesh.getEdges());

    for (const edge of mesh.getEdges()) {
        // A one-triangle edge always keeps the reference at index zero.
        expect(edge.T[0]).not.toBeNull();
        for (const tri of edge.T) {
            if (tri !== null) {
                expect(triangles.has(tri)).toBe(true);
                expect(tri.E.indexOf(edge)).toBeGreaterThanOrEqual(0);
            }
        }
        if (edge.T[0] !== null && edge.T[1] !== null) {
            expect(edge.T[0]).not.toBe(edge.T[1]);
        }
    }

    for (const tri of mesh.getTriangles()) {
        expect(mesh.getTriangle(tri.V[0], tri.V[1], tri.V[2])).toBe(tri);
        for (let i = 0; i < 3; ++i) {
            const edge = tri.E[i];
            expect(edge).not.toBeNull();
            const e = edge as ETManifoldMeshEdge;
            expect(edges.has(e)).toBe(true);

            // E[i] is the edge (V[i], V[(i+1)%3]).
            expect(new EdgeKey(false, e.V[0], e.V[1]).mapKey()).toBe(
                new EdgeKey(false, tri.V[i], tri.V[(i + 1) % 3]).mapKey());
            expect(mesh.getEdge(e.V[0], e.V[1])).toBe(e);

            // T[i] is the other triangle sharing E[i].
            const other = e.T[0] === tri ? e.T[1] : e.T[0];
            expect(tri.T[i]).toBe(other);
            if (tri.T[i] !== null) {
                expect((tri.T[i] as ETManifoldMeshTriangle).T.indexOf(tri))
                    .toBeGreaterThanOrEqual(0);
            }
        }
    }
}

// The closed tetrahedron surface with outward-facing triangles.
const tetraSurface: number[][] = [[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]];

// The closed octahedron surface with outward-facing triangles; 4 is the
// "north pole", 5 the "south pole" and (0,2,1,3) the equator.
const octaSurface: number[][] = [
    [0, 2, 4], [2, 1, 4], [1, 3, 4], [3, 0, 4],
    [2, 0, 5], [1, 2, 5], [3, 1, 5], [0, 3, 5]
];

// An annulus (a square ring): the outer square is (0,1,2,3) and the inner
// square is (4,5,6,7), both counterclockwise.
const annulus: number[][] = [
    [0, 1, 5], [0, 5, 4],
    [1, 2, 6], [1, 6, 5],
    [2, 3, 7], [2, 7, 6],
    [3, 0, 4], [3, 4, 7]
];

function build(triples: number[][], throwOnNonmanifold = true): ETManifoldMesh {
    const mesh = new ETManifoldMesh();
    mesh.throwOnNonmanifoldInsertion(throwOnNonmanifold);
    for (const t of triples) {
        mesh.insert(t[0], t[1], t[2]);
    }
    return mesh;
}

// A regular grid of size N x N quads, each quad split into two triangles with
// consistent counterclockwise winding.
function gridTriples(n: number): number[][] {
    const triples: number[][] = [];
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            const v00 = r * (n + 1) + c;
            const v10 = v00 + 1;
            const v01 = v00 + n + 1;
            const v11 = v01 + 1;
            triples.push([v00, v10, v11]);
            triples.push([v00, v11, v01]);
        }
    }
    return triples;
}

describe('ETManifoldMeshTriangle', () => {
    const tri = new ETManifoldMeshTriangle(3, 5, 7);

    it('reports which side of a directed edge it lies on', () => {
        expect(tri.whichSideOfEdge(3, 5)).toBe(+1);
        expect(tri.whichSideOfEdge(5, 7)).toBe(+1);
        expect(tri.whichSideOfEdge(7, 3)).toBe(+1);
        expect(tri.whichSideOfEdge(5, 3)).toBe(-1);
        expect(tri.whichSideOfEdge(7, 5)).toBe(-1);
        expect(tri.whichSideOfEdge(3, 7)).toBe(-1);
        expect(tri.whichSideOfEdge(3, 9)).toBe(0);
        expect(tri.whichSideOfEdge(9, 11)).toBe(0);
    });

    it('finds the vertex opposite an edge', () => {
        expect(tri.getOppositeVertexOfEdge(3, 5)).toEqual({ found: true, uOpposite: 7 });
        expect(tri.getOppositeVertexOfEdge(5, 3)).toEqual({ found: true, uOpposite: 7 });
        expect(tri.getOppositeVertexOfEdge(5, 7)).toEqual({ found: true, uOpposite: 3 });
        expect(tri.getOppositeVertexOfEdge(7, 3)).toEqual({ found: true, uOpposite: 5 });
        expect(tri.getOppositeVertexOfEdge(3, 9)).toEqual({ found: false, uOpposite: -1 });
    });

    it('finds the adjacent triangle of an edge', () => {
        // (0,1,2) and (1,0,3) share the edge {0,1}.
        const mesh = build([[0, 1, 2], [1, 0, 3]]);
        const t0 = mesh.getTriangle(0, 1, 2) as ETManifoldMeshTriangle;
        const t1 = mesh.getTriangle(1, 0, 3) as ETManifoldMeshTriangle;
        expect(t0.getAdjacentOfEdge(0, 1)).toBe(t1);
        expect(t0.getAdjacentOfEdge(1, 0)).toBe(t1);
        expect(t1.getAdjacentOfEdge(0, 1)).toBe(t0);
        expect(t0.getAdjacentOfEdge(1, 2)).toBeNull();
        expect(t0.getAdjacentOfEdge(5, 6)).toBeNull();
    });
});

describe('ETManifoldMesh insertion', () => {
    it('builds a triangle fan with the expected edges and adjacency', () => {
        const mesh = build([[0, 1, 2], [0, 2, 3], [0, 3, 4]]);
        expect(mesh.getNumTriangles()).toBe(3);
        expect(mesh.getNumEdges()).toBe(7);
        expect(edgePairs(mesh)).toEqual([
            [0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [2, 3], [3, 4]
        ]);
        expect(keys(mesh.getTriangles())).toEqual([[0, 1, 2], [0, 2, 3], [0, 3, 4]]);

        const t0 = mesh.getTriangle(0, 1, 2) as ETManifoldMeshTriangle;
        const t1 = mesh.getTriangle(0, 2, 3) as ETManifoldMeshTriangle;
        const t2 = mesh.getTriangle(0, 3, 4) as ETManifoldMeshTriangle;
        // t0 = (0,1,2): E[0] = {0,1}, E[1] = {1,2}, E[2] = {2,0}. Only the
        // last is shared, with t1.
        expect(t0.T).toEqual([null, null, t1]);
        expect(t1.T).toEqual([t0, null, t2]);
        expect(t2.T).toEqual([t1, null, null]);

        const shared = mesh.getEdge(0, 2) as ETManifoldMeshEdge;
        expect(shared.T[0]).toBe(t0);
        expect(shared.T[1]).toBe(t1);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isOriented()).toBe(true);
        checkInvariants(mesh);
    });

    it('rejects a duplicate triangle by returning null', () => {
        const mesh = build([[0, 1, 2]]);
        expect(mesh.insert(0, 1, 2)).toBeNull();
        // The key is the ordered triangle key, so a cyclic rotation is the
        // same triangle.
        expect(mesh.insert(1, 2, 0)).toBeNull();
        expect(mesh.getNumTriangles()).toBe(1);
        // The reversed triangle is a different key and is accepted; the mesh
        // is then closed but not oriented.
        expect(mesh.insert(0, 2, 1)).not.toBeNull();
        expect(mesh.getNumTriangles()).toBe(2);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(false);
    });

    it('builds the closed tetrahedron surface', () => {
        const mesh = build(tetraSurface);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getNumEdges()).toBe(6);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        for (const edge of mesh.getEdges()) {
            expect(edge.T[0]).not.toBeNull();
            expect(edge.T[1]).not.toBeNull();
        }
        for (const tri of mesh.getTriangles()) {
            expect(tri.T.every(t => t !== null)).toBe(true);
        }
        expect(mesh.getBoundaryPolygons(true)).toEqual([]);
        checkInvariants(mesh);
    });

    it('builds the closed octahedron surface', () => {
        const mesh = build(octaSurface);
        expect(mesh.getNumTriangles()).toBe(8);
        expect(mesh.getNumEdges()).toBe(12);
        // Euler characteristic of the sphere: V - E + F = 6 - 12 + 8 = 2.
        expect(6 - mesh.getNumEdges() + mesh.getNumTriangles()).toBe(2);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        checkInvariants(mesh);
    });

    it('builds a disk with boundary from a grid of quads', () => {
        const mesh = build(gridTriples(3));
        expect(mesh.getNumTriangles()).toBe(18);
        // 18*3 edge slots, 12 boundary edges counted once and the rest twice.
        expect(mesh.getNumEdges()).toBe(12 + (18 * 3 - 12) / 2);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isOriented()).toBe(true);
        checkInvariants(mesh);
    });

    it('throws when an insertion would make the mesh nonmanifold', () => {
        const mesh = build([[0, 1, 2], [1, 0, 3]]);
        // The edge {0,1} already has two triangles.
        expect(() => mesh.insert(1, 0, 4)).toThrow('Attempt to create nonmanifold mesh.');
    });

    it('throws when an insertion would make the mesh inconsistently oriented', () => {
        const mesh = build([[0, 1, 2]]);
        // (0,1,3) traverses the shared edge {0,1} in the same direction as
        // (0,1,2), so the insertion is rejected by the orientation check even
        // though the edge has only one triangle.
        expect(() => mesh.insert(0, 1, 3)).toThrow('Attempt to create nonmanifold mesh.');
    });

    // Upstream quirk (the ETManifoldMesh form of upstream issue #73 for
    // VEManifoldMesh::Insert): the graceful rejection path returns null but
    // leaves the edges created earlier in the loop in the edge map, each
    // referencing a triangle that was never added to the triangle map.
    it('leaves the rejected insertion partially in the mesh when not throwing', () => {
        const mesh = build([[0, 1, 2], [1, 0, 3]]);
        expect(mesh.getNumEdges()).toBe(5);

        const previous = mesh.throwOnNonmanifoldInsertion(false);
        expect(previous).toBe(true);

        // The loop creates the edge {1,4} first and then finds the full edge
        // {0,1}, at which point it returns null.
        expect(mesh.insert(1, 0, 4)).toBeNull();
        expect(mesh.getNumTriangles()).toBe(2);
        expect(mesh.getNumEdges()).toBe(6);

        const leaked = mesh.getEdge(1, 4) as ETManifoldMeshEdge;
        expect(leaked).not.toBeNull();
        const phantom = leaked.T[0] as ETManifoldMeshTriangle;
        expect(phantom.V).toEqual([1, 0, 4]);
        // The phantom triangle is not in the triangle map.
        expect(mesh.getTriangle(1, 0, 4)).toBeNull();
        expect(mesh.getTriangles()).toHaveLength(2);

        // Restoring the previous state is reported correctly.
        expect(mesh.throwOnNonmanifoldInsertion(true)).toBe(false);
    });

    it('accepts an inconsistently oriented closed mesh when not throwing', () => {
        const mesh = build([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 3, 2]], false);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(false);
        checkInvariants(mesh);
    });

    it('uses the edge and triangle creator callbacks', () => {
        class MyEdge extends ETManifoldMeshEdge {
            label = 'edge';
        }
        class MyTriangle extends ETManifoldMeshTriangle {
            label = 'triangle';
        }
        const mesh = new ETManifoldMesh(
            (v0, v1) => new MyEdge(v0, v1),
            (v0, v1, v2) => new MyTriangle(v0, v1, v2));
        mesh.insert(0, 1, 2);
        expect(mesh.getEdges().every(e => e instanceof MyEdge)).toBe(true);
        expect(mesh.getTriangles().every(t => t instanceof MyTriangle)).toBe(true);
        // The creators survive a copy.
        const copy = mesh.clone();
        copy.insert(0, 2, 3);
        expect(copy.getTriangles().every(t => t instanceof MyTriangle)).toBe(true);
    });
});

describe('ETManifoldMesh removal', () => {
    it('removes a triangle and restores the adjacency of its neighbors', () => {
        const mesh = build(tetraSurface);
        expect(mesh.remove(1, 2, 3)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(3);
        // The three edges of the removed triangle survive as boundary edges.
        expect(mesh.getNumEdges()).toBe(6);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isOriented()).toBe(true);
        for (const tri of mesh.getTriangles()) {
            expect(tri.T.filter(t => t === null)).toHaveLength(1);
        }
        checkInvariants(mesh);

        // The key is the ordered triangle key, so any cyclic rotation names
        // the same triangle.
        expect(mesh.remove(3, 2, 0)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(2);
        // The edge {2,3} lost its last triangle and is gone.
        expect(mesh.getEdge(2, 3)).toBeNull();
        expect(mesh.getNumEdges()).toBe(5);
        checkInvariants(mesh);
    });

    it('returns false for a triangle that is not in the mesh', () => {
        const mesh = build(tetraSurface);
        expect(mesh.remove(0, 1, 5)).toBe(false);
        // The reversed triangle is a different ordered key.
        expect(mesh.remove(1, 2, 3)).toBe(true);
        expect(mesh.remove(1, 2, 3)).toBe(false);
        expect(mesh.getNumTriangles()).toBe(3);
    });

    it('empties the mesh when all triangles are removed', () => {
        const mesh = build(octaSurface);
        for (const t of octaSurface) {
            expect(mesh.remove(t[0], t[1], t[2])).toBe(true);
        }
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        // A mesh with no edges is vacuously closed and oriented.
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
    });

    it('clears the mesh', () => {
        const mesh = build(octaSurface);
        mesh.clear();
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getTriangles()).toEqual([]);
        expect(mesh.getEdges()).toEqual([]);
    });

    it('supports reinsertion after removal', () => {
        const mesh = build(tetraSurface);
        expect(mesh.remove(1, 2, 3)).toBe(true);
        expect(mesh.insert(1, 2, 3)).not.toBeNull();
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getNumEdges()).toBe(6);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        checkInvariants(mesh);
    });
});

describe('ETManifoldMesh components', () => {
    it('finds a single component for a connected mesh', () => {
        const mesh = build(tetraSurface);
        const components = mesh.getComponents();
        expect(components).toHaveLength(1);
        expect(components[0]).toHaveLength(4);
        expect(mesh.getComponentKeys()).toHaveLength(1);
        expect(mesh.getComponentKeys()[0]).toHaveLength(4);
    });

    it('finds the components of two disjoint tetrahedron surfaces', () => {
        const shifted = tetraSurface.map(t => t.map(v => v + 10));
        const mesh = build(tetraSurface.concat(shifted));
        expect(mesh.getNumTriangles()).toBe(8);

        const components = mesh.getComponents();
        expect(components).toHaveLength(2);
        expect(components.map(c => c.length)).toEqual([4, 4]);
        // The components are seeded in increasing triangle-key order, so the
        // component of the low-index tetrahedron comes first.
        expect(components[0].every(t => t.V.every(v => v < 10))).toBe(true);
        expect(components[1].every(t => t.V.every(v => v >= 10))).toBe(true);

        const keyComponents = mesh.getComponentKeys();
        expect(keyComponents.map(c => c.length)).toEqual([4, 4]);
        for (let i = 0; i < 2; ++i) {
            expect(keyComponents[i].map(k => k.V.slice())).toEqual(keys(components[i]));
        }
    });

    it('finds one component per isolated triangle', () => {
        const mesh = build([[0, 1, 2], [3, 4, 5], [6, 7, 8]]);
        const components = mesh.getComponents();
        expect(components).toHaveLength(3);
        expect(components.map(c => c.length)).toEqual([1, 1, 1]);
        expect(keys(components.map(c => c[0]))).toEqual([[0, 1, 2], [3, 4, 5], [6, 7, 8]]);
    });

    it('has no components for an empty mesh', () => {
        const mesh = new ETManifoldMesh();
        expect(mesh.getComponents()).toEqual([]);
        expect(mesh.getComponentKeys()).toEqual([]);
    });
});

describe('ETManifoldMesh compact graph and chirality', () => {
    it('creates a compact graph with symmetric adjacency', () => {
        const mesh = build(tetraSurface);
        const { triangles, adjacents } = mesh.createCompactGraph();
        expect(triangles).toHaveLength(4);
        expect(adjacents).toHaveLength(4);

        // Every triangle of the closed surface has three neighbors, and the
        // adjacency relation is symmetric.
        for (let i = 0; i < 4; ++i) {
            expect(adjacents[i].every(a => a >= 0 && a < 4)).toBe(true);
            for (const a of adjacents[i]) {
                expect(adjacents[a]).toContain(i);
            }
        }

        // The vertex triples are those of the mesh triangles, in the order in
        // which getTriangles() enumerates them.
        expect(triangles).toEqual(mesh.getTriangles().map(t => t.V.slice()));
    });

    it('marks the boundary edges of a compact graph with -1', () => {
        const mesh = build([[0, 1, 2], [0, 2, 3]]);
        const { triangles, adjacents } = mesh.createCompactGraph();
        expect(triangles).toEqual([[0, 1, 2], [0, 2, 3]]);
        // (0,1,2) has edges (0,1), (1,2), (2,0); only (2,0) is shared.
        expect(adjacents[0]).toEqual([-1, -1, 1]);
        // (0,2,3) has edges (0,2), (2,3), (3,0); only (0,2) is shared.
        expect(adjacents[1]).toEqual([0, -1, -1]);
    });

    it('throws when creating a compact graph of an empty mesh', () => {
        const mesh = new ETManifoldMesh();
        expect(() => mesh.createCompactGraph()).toThrow('Invalid input.');
    });

    it('throws for invalid inputs to getComponentsConsistentChirality', () => {
        expect(() => ETManifoldMesh.getComponentsConsistentChirality([], []))
            .toThrow('Invalid inputs.');
        expect(() => ETManifoldMesh.getComponentsConsistentChirality([[0, 1, 2]], []))
            .toThrow('Invalid inputs.');
    });

    it('makes the chirality consistent for two adjacent triangles', () => {
        // (0,1,2) and (0,3,2) traverse the shared edge {0,2} in the same
        // direction (<2,0>), so the second must be flipped. Each triangle
        // shares the edge at its index 2.
        const triangles = [[0, 1, 2], [0, 3, 2]];
        const adjacents = [[-1, -1, 1], [-1, -1, 0]];
        const { components, numComponentTriangles } =
            ETManifoldMesh.getComponentsConsistentChirality(triangles, adjacents);
        expect(components).toEqual([0, 1]);
        expect(numComponentTriangles).toEqual([2]);
        expect(triangles[0]).toEqual([0, 1, 2]);
        // (0,3,2) becomes (2,3,0), a cyclic rotation of (0,2,3), which
        // traverses the shared edge as <0,2>.
        expect(triangles[1]).toEqual([2, 3, 0]);
    });

    it('reports the sizes of the components of a compact graph', () => {
        const shifted = tetraSurface.map(t => t.map(v => v + 10));
        const mesh = build(tetraSurface.concat(shifted));
        const { triangles, adjacents } = mesh.createCompactGraph();
        const { components, numComponentTriangles } =
            ETManifoldMesh.getComponentsConsistentChirality(triangles, adjacents);
        expect(numComponentTriangles).toEqual([4, 4]);
        expect(components).toHaveLength(8);
        expect(new Set(components).size).toBe(8);
        expect(numComponentTriangles.reduce((a, b) => a + b, 0)).toBe(triangles.length);
    });

    it('repairs a mesh with one flipped triangle', () => {
        // The tetrahedron surface with (1,2,3) flipped to (1,3,2). The
        // insertion is only possible with the nonmanifold check disabled.
        const mesh = build([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 3, 2]], false);
        expect(mesh.isOriented()).toBe(false);

        mesh.makeConsistentChirality();

        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getNumEdges()).toBe(6);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isOriented()).toBe(true);
        expect(unorderedTriples(mesh)).toEqual([
            [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]
        ]);
        checkInvariants(mesh);
    });

    it('leaves an already consistent mesh unchanged', () => {
        const mesh = build(octaSurface);
        const before = keys(mesh.getTriangles());
        mesh.makeConsistentChirality();
        expect(keys(mesh.getTriangles())).toEqual(before);
        expect(mesh.isOriented()).toBe(true);
        checkInvariants(mesh);
    });
});

describe('ETManifoldMesh boundary polygons', () => {
    it('extracts the boundary of a triangle fan', () => {
        const mesh = build([[0, 1, 2], [0, 2, 3], [0, 3, 4]]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([[0, 1, 2, 3, 4, 0]]);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2, 3, 4]]);
    });

    it('extracts the boundary of a single triangle', () => {
        const mesh = build([[0, 1, 2]]);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2]]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([[0, 1, 2, 0]]);
    });

    it('extracts the two boundary polygons of an annulus', () => {
        const mesh = build(annulus);
        expect(mesh.getNumTriangles()).toBe(8);
        expect(mesh.getNumEdges()).toBe(16);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isOriented()).toBe(true);

        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons).toHaveLength(2);
        // The outer boundary is traversed counterclockwise, the inner one
        // clockwise (as seen from the same side).
        expect(polygons[0]).toEqual([0, 1, 2, 3]);
        expect(polygons[1]).toEqual([4, 7, 6, 5]);
        expect(mesh.getBoundaryPolygons(true)).toEqual([
            [0, 1, 2, 3, 0], [4, 7, 6, 5, 4]
        ]);
        checkInvariants(mesh);
    });

    it('extracts the boundary of a grid of quads', () => {
        const n = 3;
        const mesh = build(gridTriples(n));
        const polygons = mesh.getBoundaryPolygons(false);
        expect(polygons).toHaveLength(1);
        expect(polygons[0]).toHaveLength(4 * n);
        // The polygon starts at vertex 0 and traverses the bottom row first.
        expect(polygons[0].slice(0, n + 1)).toEqual([0, 1, 2, 3]);
        // Every boundary vertex is on the border of the grid.
        for (const v of polygons[0]) {
            const r = Math.floor(v / (n + 1));
            const c = v % (n + 1);
            expect(r === 0 || r === n || c === 0 || c === n).toBe(true);
        }
        // The consecutive pairs are edges of the mesh with one triangle.
        for (let i = 0; i < polygons[0].length; ++i) {
            const v0 = polygons[0][i];
            const v1 = polygons[0][(i + 1) % polygons[0].length];
            const edge = mesh.getEdge(v0, v1) as ETManifoldMeshEdge;
            expect(edge).not.toBeNull();
            expect(edge.T[1]).toBeNull();
        }
    });

    it('finds one polygon per boundary component of a disconnected mesh', () => {
        const mesh = build([[0, 1, 2], [3, 4, 5]]);
        expect(mesh.getBoundaryPolygons(false)).toEqual([[0, 1, 2], [3, 4, 5]]);
    });
});

describe('ETManifoldMesh copying', () => {
    it('clones a mesh so that the copies are independent', () => {
        const mesh = build(octaSurface);
        const copy = mesh.clone();

        expect(copy.getNumTriangles()).toBe(mesh.getNumTriangles());
        expect(copy.getNumEdges()).toBe(mesh.getNumEdges());
        expect(keys(copy.getTriangles())).toEqual(keys(mesh.getTriangles()));
        expect(copy.isClosed()).toBe(true);
        expect(copy.isOriented()).toBe(true);
        checkInvariants(copy);

        // No object is shared between the meshes.
        const originalTriangles = new Set<ETManifoldMeshTriangle>(mesh.getTriangles());
        expect(copy.getTriangles().some(t => originalTriangles.has(t))).toBe(false);
        const originalEdges = new Set<ETManifoldMeshEdge>(mesh.getEdges());
        expect(copy.getEdges().some(e => originalEdges.has(e))).toBe(false);

        // Modifying the copy does not affect the original.
        expect(copy.remove(0, 2, 4)).toBe(true);
        expect(copy.getNumTriangles()).toBe(7);
        expect(mesh.getNumTriangles()).toBe(8);
        expect(mesh.isClosed()).toBe(true);
    });

    it('reinserts using the triangle keys, so the vertices may be rotated', () => {
        // The upstream operator= reinserts element.first.V[], the ordered
        // triangle key, not the triangle's own vertex order.
        const mesh = build([[2, 0, 1]]);
        expect((mesh.getTriangles()[0]).V).toEqual([2, 0, 1]);
        const copy = mesh.clone();
        expect((copy.getTriangles()[0]).V).toEqual([0, 1, 2]);
    });

    it('copies the nonmanifold-insertion state', () => {
        const mesh = build(tetraSurface);
        mesh.throwOnNonmanifoldInsertion(false);
        const copy = mesh.clone();
        // The copy does not throw either.
        expect(copy.insert(0, 1, 2)).toBeNull();

        const target = build([[10, 11, 12]]);
        target.assign(mesh);
        expect(target.getNumTriangles()).toBe(4);
        expect(target.getTriangle(10, 11, 12)).toBeNull();
        expect(keys(target.getTriangles())).toEqual(keys(mesh.getTriangles()));
    });
});

describe('ETManifoldMesh randomized consistency', () => {
    // A deterministic pseudorandom generator so the test is reproducible.
    function makeRandom(seed: number): () => number {
        let state = seed >>> 0;
        return () => {
            state = (state * 1664525 + 1013904223) >>> 0;
            return state / 4294967296;
        };
    }

    it('keeps the mesh consistent under random removals and reinsertions', () => {
        const random = makeRandom(12345);
        const triples = gridTriples(4);

        for (let trial = 0; trial < 20; ++trial) {
            const mesh = build(triples);
            const present = new Set<number>(triples.map((_, i) => i));

            // Remove a random subset.
            for (let i = 0; i < triples.length; ++i) {
                if (random() < 0.4) {
                    const t = triples[i];
                    expect(mesh.remove(t[0], t[1], t[2])).toBe(true);
                    present.delete(i);
                }
            }
            expect(mesh.getNumTriangles()).toBe(present.size);
            checkInvariants(mesh);

            // The number of edges equals the number of distinct unordered
            // vertex pairs of the remaining triangles.
            const pairs = new Set<string>();
            for (const i of present) {
                const t = triples[i];
                for (let j = 0; j < 3; ++j) {
                    pairs.add(new EdgeKey(false, t[j], t[(j + 1) % 3]).mapKey());
                }
            }
            expect(mesh.getNumEdges()).toBe(pairs.size);

            // The mesh remains oriented; it is closed only if it is empty.
            expect(mesh.isOriented()).toBe(true);
            expect(mesh.isClosed()).toBe(present.size === 0);

            // The number of components matches an independent union-find over
            // the remaining triangles.
            const indices = Array.from(present).sort((a, b) => a - b);
            const parent = new Map<number, number>(indices.map(i => [i, i]));
            const find = (i: number): number => {
                let r = i;
                while (parent.get(r) !== r) {
                    r = parent.get(r) as number;
                }
                return r;
            };
            const edgeOwner = new Map<string, number>();
            for (const i of indices) {
                const t = triples[i];
                for (let j = 0; j < 3; ++j) {
                    const key = new EdgeKey(false, t[j], t[(j + 1) % 3]).mapKey();
                    const other = edgeOwner.get(key);
                    if (other === undefined) {
                        edgeOwner.set(key, i);
                    } else {
                        parent.set(find(other), find(i));
                    }
                }
            }
            const roots = new Set<number>(indices.map(i => find(i)));
            expect(mesh.getComponents()).toHaveLength(roots.size);
            expect(mesh.getComponents().reduce((sum, c) => sum + c.length, 0))
                .toBe(present.size);

            // Reinsert everything that was removed; the mesh is restored.
            for (let i = 0; i < triples.length; ++i) {
                if (!present.has(i)) {
                    const t = triples[i];
                    expect(mesh.insert(t[0], t[1], t[2])).not.toBeNull();
                }
            }
            expect(mesh.getNumTriangles()).toBe(triples.length);
            expect(keys(mesh.getTriangles())).toEqual(keys(build(triples).getTriangles()));
            checkInvariants(mesh);
        }
    });
});
