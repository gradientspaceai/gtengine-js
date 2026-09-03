import { describe, it, expect } from 'vitest';
import {
    ETNonmanifoldMesh,
    ETNonmanifoldMeshEdge,
    ETNonmanifoldMeshTriangle
} from '../src/ETNonmanifoldMesh.js';
import { TriangleKey } from '../src/TriangleKey.js';

// The keys of the triangles of a mesh (or of a component), as arrays of the
// three vertex indices.
function keys(triangles: ETNonmanifoldMeshTriangle[]): number[][] {
    return triangles.map(t => new TriangleKey(true, t.V[0], t.V[1], t.V[2]).V.slice());
}

function keyArrays(components: TriangleKey[][]): number[][][] {
    return components.map(c => c.map(k => k.V.slice()));
}

// The unordered vertex pairs of the edges, in mesh iteration order.
function edgePairs(mesh: ETNonmanifoldMesh): number[][] {
    return mesh.getEdgeKeys().map(k => k.V.slice());
}

// A triangle fan around vertex 0: (0,1,2), (0,2,3), (0,3,4).
function makeFan(): ETNonmanifoldMesh {
    const mesh = new ETNonmanifoldMesh();
    expect(mesh.insert(0, 1, 2)).not.toBeNull();
    expect(mesh.insert(0, 2, 3)).not.toBeNull();
    expect(mesh.insert(0, 3, 4)).not.toBeNull();
    return mesh;
}

// The closed surface of a tetrahedron with vertices 0,1,2,3.
function makeTetrahedron(): ETNonmanifoldMesh {
    const mesh = new ETNonmanifoldMesh();
    expect(mesh.insert(0, 2, 1)).not.toBeNull();
    expect(mesh.insert(0, 1, 3)).not.toBeNull();
    expect(mesh.insert(0, 3, 2)).not.toBeNull();
    expect(mesh.insert(1, 2, 3)).not.toBeNull();
    return mesh;
}

// Three triangles sharing the edge <0,1>: a genuinely nonmanifold mesh.
function makeBowtieEdge(): ETNonmanifoldMesh {
    const mesh = new ETNonmanifoldMesh();
    expect(mesh.insert(0, 1, 2)).not.toBeNull();
    expect(mesh.insert(0, 1, 3)).not.toBeNull();
    expect(mesh.insert(0, 1, 4)).not.toBeNull();
    return mesh;
}

describe('ETNonmanifoldMesh insertion', () => {
    it('creates the triangle and its three edges', () => {
        const mesh = new ETNonmanifoldMesh();
        const tri = mesh.insert(3, 7, 5);
        expect(tri).not.toBeNull();
        expect(tri!.V).toEqual([3, 7, 5]);
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumEdges()).toBe(3);

        // E[i] is the edge (V[i], V[(i+1)%3]) and stores the vertices in the
        // order in which the edge was first encountered, not sorted.
        expect(tri!.E[0]!.V).toEqual([3, 7]);
        expect(tri!.E[1]!.V).toEqual([7, 5]);
        expect(tri!.E[2]!.V).toEqual([5, 3]);

        // The edge key is unordered, so the lookup order does not matter.
        expect(mesh.getEdge(3, 7)).toBe(tri!.E[0]);
        expect(mesh.getEdge(7, 3)).toBe(tri!.E[0]);
        expect(mesh.getEdge(3, 5)).toBe(tri!.E[2]);
        expect(mesh.getEdge(4, 9)).toBeNull();

        // Every edge knows the triangle.
        for (const edge of mesh.getEdges()) {
            expect(edge.getTriangles()).toEqual([tri]);
        }
    });

    it('finds a triangle by any cyclic rotation but not by the reversal', () => {
        const mesh = new ETNonmanifoldMesh();
        const tri = mesh.insert(3, 7, 5)!;
        expect(mesh.getTriangle(3, 7, 5)).toBe(tri);
        expect(mesh.getTriangle(7, 5, 3)).toBe(tri);
        expect(mesh.getTriangle(5, 3, 7)).toBe(tri);
        // The opposite winding is a different ordered triangle key.
        expect(mesh.getTriangle(3, 5, 7)).toBeNull();
    });

    it('rejects a duplicate triangle, including its cyclic rotations', () => {
        const mesh = new ETNonmanifoldMesh();
        expect(mesh.insert(0, 1, 2)).not.toBeNull();
        expect(mesh.insert(0, 1, 2)).toBeNull();
        expect(mesh.insert(1, 2, 0)).toBeNull();
        expect(mesh.insert(2, 0, 1)).toBeNull();
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumEdges()).toBe(3);

        // The reversed winding is a distinct triangle sharing all edges.
        expect(mesh.insert(0, 2, 1)).not.toBeNull();
        expect(mesh.getNumTriangles()).toBe(2);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isManifold()).toBe(true);
    });

    it('shares an edge between adjacent triangles of a fan', () => {
        const mesh = makeFan();
        expect(mesh.getNumTriangles()).toBe(3);
        // Edges: <0,1>,<0,2>,<0,3>,<0,4>,<1,2>,<2,3>,<3,4>.
        expect(mesh.getNumEdges()).toBe(7);
        expect(edgePairs(mesh)).toEqual([
            [0, 1], [0, 2], [0, 3], [0, 4], [1, 2], [2, 3], [3, 4]
        ]);

        const shared = mesh.getEdge(0, 2)!;
        expect(keys(shared.getTriangles())).toEqual([[0, 1, 2], [0, 2, 3]]);

        const boundary = mesh.getEdge(1, 2)!;
        expect(keys(boundary.getTriangles())).toEqual([[0, 1, 2]]);

        expect(mesh.isManifold()).toBe(true);
        // The fan has boundary edges shared once only.
        expect(mesh.isClosed()).toBe(false);
    });

    it('accepts an arbitrary number of triangles per edge', () => {
        const mesh = makeBowtieEdge();
        expect(mesh.getNumTriangles()).toBe(3);
        // <0,1> plus <0,k> and <1,k> for k = 2,3,4.
        expect(mesh.getNumEdges()).toBe(7);

        const shared = mesh.getEdge(0, 1)!;
        expect(shared.T.size).toBe(3);
        expect(keys(shared.getTriangles())).toEqual([[0, 1, 2], [0, 1, 3], [0, 1, 4]]);

        expect(mesh.isManifold()).toBe(false);
        expect(mesh.isClosed()).toBe(false);
    });
});

describe('ETNonmanifoldMesh closed and manifold queries', () => {
    it('reports the tetrahedron surface as closed and manifold', () => {
        const mesh = makeTetrahedron();
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getNumEdges()).toBe(6);
        expect(edgePairs(mesh)).toEqual([
            [0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]
        ]);
        for (const edge of mesh.getEdges()) {
            expect(edge.T.size).toBe(2);
        }
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isManifold()).toBe(true);

        // Removing one triangle opens the surface but keeps it manifold.
        expect(mesh.remove(1, 2, 3)).toBe(true);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isManifold()).toBe(true);
    });

    it('reports an empty mesh as closed and manifold (vacuously)', () => {
        const mesh = new ETNonmanifoldMesh();
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.isManifold()).toBe(true);
        expect(mesh.getEdges()).toEqual([]);
        expect(mesh.getTriangles()).toEqual([]);
    });

    it('detects the nonmanifold edges of two tetrahedra glued along a face', () => {
        // Two tetrahedra glued along the face (1,2,3) share that face's
        // edges among three triangles each.
        const mesh = makeTetrahedron();
        expect(mesh.insert(1, 3, 4)).not.toBeNull();
        expect(mesh.insert(1, 4, 2)).not.toBeNull();
        expect(mesh.insert(2, 4, 3)).not.toBeNull();
        expect(mesh.getNumTriangles()).toBe(7);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isManifold()).toBe(false);
        expect(mesh.getEdge(1, 2)!.T.size).toBe(3);
        expect(mesh.getEdge(1, 4)!.T.size).toBe(2);
    });
});

describe('ETNonmanifoldMesh removal', () => {
    it('restores the counts when every triangle is removed', () => {
        const mesh = makeFan();
        expect(mesh.remove(0, 2, 3)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(2);
        // The interior edges <0,2> and <0,3> survive because the fan's other
        // triangles still reference them; <2,3> is gone.
        expect(mesh.getNumEdges()).toBe(6);
        expect(mesh.getEdge(2, 3)).toBeNull();
        expect(keys(mesh.getEdge(0, 2)!.getTriangles())).toEqual([[0, 1, 2]]);
        expect(keys(mesh.getEdge(0, 3)!.getTriangles())).toEqual([[0, 3, 4]]);

        expect(mesh.remove(0, 1, 2)).toBe(true);
        expect(mesh.remove(0, 3, 4)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
    });

    it('accepts any cyclic rotation and rejects an absent triangle', () => {
        const mesh = makeTetrahedron();
        // (1,0,2) is the cyclic rotation of the inserted triangle (0,2,1).
        expect(mesh.remove(1, 0, 2)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(3);
        // The same triangle can no longer be removed.
        expect(mesh.remove(0, 2, 1)).toBe(false);
        // The reversed winding was never in the mesh.
        expect(mesh.remove(0, 1, 2)).toBe(false);
        expect(mesh.getNumTriangles()).toBe(3);
    });

    it('removes only the requested triangle from a nonmanifold edge', () => {
        const mesh = makeBowtieEdge();
        expect(mesh.remove(0, 1, 3)).toBe(true);
        const shared = mesh.getEdge(0, 1)!;
        expect(keys(shared.getTriangles())).toEqual([[0, 1, 2], [0, 1, 4]]);
        expect(mesh.getEdge(0, 3)).toBeNull();
        expect(mesh.getEdge(1, 3)).toBeNull();
        expect(mesh.getNumEdges()).toBe(5);
        expect(mesh.isManifold()).toBe(true);
    });

    it('empties the mesh with clear', () => {
        const mesh = makeTetrahedron();
        mesh.clear();
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getComponents()).toEqual([]);
        // Reinsertion after a clear works as for a fresh mesh.
        expect(mesh.insert(0, 2, 1)).not.toBeNull();
        expect(mesh.getNumEdges()).toBe(3);
    });
});

describe('ETNonmanifoldMesh connected components', () => {
    it('finds one component for a connected mesh', () => {
        const mesh = makeTetrahedron();
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect(components[0].length).toBe(4);
        // Every triangle of the mesh appears exactly once.
        const sorted = keys(components[0]).slice().sort();
        expect(sorted).toEqual(keys(mesh.getTriangles()).slice().sort());
    });

    it('finds the separate components of a disconnected mesh', () => {
        const mesh = new ETNonmanifoldMesh();
        // Component A: two triangles sharing edge <0,2>.
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        // Component B: a single triangle on disjoint vertices.
        mesh.insert(10, 11, 12);
        // Component C: two triangles sharing edge <20,22>.
        mesh.insert(20, 21, 22);
        mesh.insert(20, 22, 23);

        const components = mesh.getComponents();
        expect(components.length).toBe(3);
        expect(components.map(c => c.length)).toEqual([2, 1, 2]);
        // The components are discovered in triangle-key order.
        expect(keys(components[0]).slice().sort()).toEqual([[0, 1, 2], [0, 2, 3]]);
        expect(keys(components[1])).toEqual([[10, 11, 12]]);
        expect(keys(components[2]).slice().sort()).toEqual([[20, 21, 22], [20, 22, 23]]);
    });

    it('connects triangles through a nonmanifold edge', () => {
        const mesh = makeBowtieEdge();
        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect(components[0].length).toBe(3);
    });

    it('agrees with the triangle-key overload', () => {
        const mesh = new ETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.insert(10, 11, 12);

        const byRef = mesh.getComponents();
        const byKey = mesh.getComponentKeys();
        expect(keyArrays(byKey)).toEqual(byRef.map(c => keys(c)));

        // The key components survive a clear of the mesh.
        mesh.clear();
        expect(keyArrays(byKey).length).toBe(2);
    });

    it('visits a fan without duplicating triangles', () => {
        // A fan of eight triangles around vertex 0 exercises the preallocated
        // depth-first-search stack.
        const mesh = new ETNonmanifoldMesh();
        for (let i = 1; i <= 8; ++i) {
            mesh.insert(0, i, i === 8 ? 1 : i + 1);
        }
        expect(mesh.getNumTriangles()).toBe(8);
        // Each spoke edge <0,i> is shared by two triangles, but each rim edge
        // <i,i+1> belongs to a single triangle, so the mesh has a boundary.
        expect(mesh.getEdge(0, 3)!.T.size).toBe(2);
        expect(mesh.getEdge(3, 4)!.T.size).toBe(1);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.isManifold()).toBe(true);

        const components = mesh.getComponents();
        expect(components.length).toBe(1);
        expect(components[0].length).toBe(8);
        expect(new Set(components[0]).size).toBe(8);
    });
});

describe('ETNonmanifoldMesh copying', () => {
    it('makes an independent deep copy', () => {
        const mesh = makeTetrahedron();
        const copy = mesh.clone();
        expect(copy.getNumTriangles()).toBe(4);
        expect(copy.getNumEdges()).toBe(6);
        expect(copy.isClosed()).toBe(true);

        // The edge and triangle objects are not shared.
        const original = mesh.getTriangle(0, 2, 1)!;
        const copied = copy.getTriangle(0, 2, 1)!;
        expect(copied).not.toBe(original);
        expect(copy.getEdge(0, 1)).not.toBe(mesh.getEdge(0, 1));

        // Modifying the copy leaves the original alone.
        expect(copy.remove(0, 2, 1)).toBe(true);
        expect(copy.getNumTriangles()).toBe(3);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.isClosed()).toBe(true);
    });

    it('copies the triangles using the triangle keys, as upstream does', () => {
        // Upstream operator= reinserts using the key vertices, so a copied
        // triangle is the cyclic rotation of the original that starts at the
        // smallest vertex index.
        const mesh = new ETNonmanifoldMesh();
        const tri = mesh.insert(7, 3, 5)!;
        expect(tri.V).toEqual([7, 3, 5]);
        const copied = mesh.clone().getTriangle(7, 3, 5)!;
        expect(copied.V).toEqual([3, 5, 7]);
    });

    it('assign replaces the previous contents', () => {
        const target = makeFan();
        const source = makeTetrahedron();
        expect(target.assign(source)).toBe(target);
        expect(target.getNumTriangles()).toBe(4);
        expect(target.getNumEdges()).toBe(6);
        expect(target.getTriangle(0, 1, 2)).toBeNull();
        expect(target.isClosed()).toBe(true);
    });
});

describe('ETNonmanifoldMesh creators', () => {
    it('uses the supplied edge and triangle creators', () => {
        class MyEdge extends ETNonmanifoldMeshEdge {
            label: string;
            constructor(v0: number, v1: number) {
                super(v0, v1);
                this.label = `E${v0}-${v1}`;
            }
        }
        class MyTriangle extends ETNonmanifoldMeshTriangle {
            label: string;
            constructor(v0: number, v1: number, v2: number) {
                super(v0, v1, v2);
                this.label = `T${v0}-${v1}-${v2}`;
            }
        }

        const mesh = new ETNonmanifoldMesh(
            (v0, v1) => new MyEdge(v0, v1),
            (v0, v1, v2) => new MyTriangle(v0, v1, v2));
        const tri = mesh.insert(4, 5, 6)!;
        expect(tri).toBeInstanceOf(MyTriangle);
        expect((tri as MyTriangle).label).toBe('T4-5-6');
        const edge = mesh.getEdge(5, 6)!;
        expect(edge).toBeInstanceOf(MyEdge);
        expect((edge as MyEdge).label).toBe('E5-6');
    });
});

describe('ETNonmanifoldMesh comparisons and ordering', () => {
    it('orders edges and triangles by their feature keys', () => {
        const mesh = new ETNonmanifoldMesh();
        mesh.insert(5, 6, 7);
        mesh.insert(1, 2, 3);
        mesh.insert(2, 4, 9);
        expect(keys(mesh.getTriangles())).toEqual([[1, 2, 3], [2, 4, 9], [5, 6, 7]]);
        expect(edgePairs(mesh)).toEqual([
            [1, 2], [1, 3], [2, 3], [2, 4], [2, 9], [4, 9], [5, 6], [5, 7], [6, 7]
        ]);
    });

    it('compares edges and triangles the way the upstream operator< does', () => {
        const e0 = new ETNonmanifoldMeshEdge(3, 1);
        const e1 = new ETNonmanifoldMeshEdge(1, 5);
        // The comparison uses the unordered key, so e0 is <1,3> and e1 <1,5>.
        expect(e0.lessThan(e1)).toBe(true);
        expect(e1.lessThan(e0)).toBe(false);

        const t0 = new ETNonmanifoldMeshTriangle(2, 0, 1);
        const t1 = new ETNonmanifoldMeshTriangle(0, 1, 3);
        // The comparison uses the ordered key: (0,1,2) < (0,1,3).
        expect(t0.lessThan(t1)).toBe(true);
        expect(t1.lessThan(t0)).toBe(false);
        expect(t0.lessThan(new ETNonmanifoldMeshTriangle(1, 2, 0))).toBe(false);
    });
});

describe('ETNonmanifoldMesh degenerate triangles (upstream behavior)', () => {
    it('accepts a degenerate triangle, which shares one edge object twice', () => {
        // Upstream does not reject a triangle with repeated vertices. The
        // triangle (0,0,1) has edges <1,0>, <0,0> and <0,1>; the first and
        // third have the same unordered key, so E[1] and E[2] are the same
        // edge object.
        const mesh = new ETNonmanifoldMesh();
        const tri = mesh.insert(0, 0, 1);
        expect(tri).not.toBeNull();
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumEdges()).toBe(2);
        expect(tri!.E[1]).toBe(tri!.E[2]);
        expect(mesh.getEdge(0, 0)!.T.size).toBe(1);
        expect(mesh.getEdge(0, 1)!.T.size).toBe(1);
    });

    it('throws when removing a degenerate triangle (upstream bug, preserved)', () => {
        // The removal loop erases the triangle from each of E[0], E[1], E[2].
        // The duplicated edge yields zero erasures on the second visit, which
        // trips the upstream LogAssert(numRemoved > 0, ...).
        const mesh = new ETNonmanifoldMesh();
        expect(mesh.insert(0, 0, 1)).not.toBeNull();
        expect(() => mesh.remove(0, 0, 1)).toThrow('Unexpected condition.');
    });
});
