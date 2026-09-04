import { describe, it, expect } from 'vitest';
import {
    ETNonmanifoldMesh,
    ETNonmanifoldMeshEdge,
    ETNonmanifoldMeshTriangle
} from '../src/ETNonmanifoldMesh.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification block (V16): model-based testing of insert/remove against a
// brute-force adjacency model, over random (mostly nonmanifold) meshes.
// ---------------------------------------------------------------------------

type Tri = [number, number, number];

function triKey(v0: number, v1: number, v2: number): string {
    return new TriangleKey(true, v0, v1, v2).V.join(',');
}

function edgeKey(v0: number, v1: number): string {
    return new EdgeKey(false, v0, v1).V.join(',');
}

// The three unordered edges of a triangle, as keys.
function triEdgeKeys(t: Tri): string[] {
    return [edgeKey(t[2], t[0]), edgeKey(t[0], t[1]), edgeKey(t[1], t[2])];
}

// A brute-force model of the mesh: a set of ordered triangle keys.
class MeshModel {
    readonly triangles = new Map<string, Tri>();

    insert(t: Tri): boolean {
        const key = triKey(t[0], t[1], t[2]);
        if (this.triangles.has(key)) { return false; }
        this.triangles.set(key, t);
        return true;
    }

    remove(t: Tri): boolean {
        return this.triangles.delete(triKey(t[0], t[1], t[2]));
    }

    // edgeKey -> sorted list of the triangle keys sharing that edge.
    edgeMap(): Map<string, string[]> {
        const m = new Map<string, string[]>();
        for (const [key, t] of this.triangles) {
            for (const e of triEdgeKeys(t)) {
                const list = m.get(e);
                if (list === undefined) { m.set(e, [key]); }
                else if (!list.includes(key)) { list.push(key); }
            }
        }
        for (const list of m.values()) { list.sort(compareKeyStrings); }
        return m;
    }

    // Connected components of the triangle graph, triangles adjacent when
    // they share an unordered edge.
    components(): string[][] {
        const byEdge = this.edgeMap();
        const parent = new Map<string, string>();
        const find = (a: string): string => {
            let r = a;
            while (parent.get(r) !== r) { r = parent.get(r) as string; }
            return r;
        };
        for (const key of this.triangles.keys()) { parent.set(key, key); }
        for (const list of byEdge.values()) {
            for (let i = 1; i < list.length; ++i) {
                parent.set(find(list[i]), find(list[0]));
            }
        }
        const groups = new Map<string, string[]>();
        for (const key of this.triangles.keys()) {
            const root = find(key);
            const g = groups.get(root);
            if (g === undefined) { groups.set(root, [key]); } else { g.push(key); }
        }
        return Array.from(groups.values()).map(g => g.slice().sort(compareKeyStrings));
    }
}

// The lexicographic order of the three vertex indices, which is what
// FeatureKey.compare implements and what std::map iteration follows.
function compareKeyStrings(a: string, b: string): number {
    const x = a.split(',').map(Number);
    const y = b.split(',').map(Number);
    for (let i = 0; i < x.length; ++i) {
        if (x[i] !== y[i]) { return x[i] - y[i]; }
    }
    return 0;
}

function meshTriangleKeys(mesh: ETNonmanifoldMesh): string[] {
    return mesh.getTriangleKeys().map(k => k.V.join(','));
}

// Compare the mesh with the model in full.
function expectMeshMatchesModel(mesh: ETNonmanifoldMesh, model: MeshModel): void {
    const modelKeys = Array.from(model.triangles.keys()).sort(compareKeyStrings);
    expect(mesh.getNumTriangles()).toBe(modelKeys.length);
    // getTriangleKeys() must already be in std::map order.
    expect(meshTriangleKeys(mesh)).toEqual(modelKeys);

    const modelEdges = model.edgeMap();
    const edgeKeys = mesh.getEdgeKeys().map(k => k.V.join(','));
    expect(mesh.getNumEdges()).toBe(modelEdges.size);
    expect(edgeKeys).toEqual(
        Array.from(modelEdges.keys()).sort(compareKeyStrings));

    for (const edge of mesh.getEdges()) {
        const key = edgeKey(edge.V[0], edge.V[1]);
        const expected = modelEdges.get(key) as string[];
        expect(expected).toBeDefined();
        expect(edge.T.size).toBe(expected.length);
        // getTriangles() must be in the order of the upstream std::set of
        // weak pointers, which is ordered by TriangleKey<true>.
        expect(edge.getTriangles().map(t => triKey(t.V[0], t.V[1], t.V[2])))
            .toEqual(expected);
    }

    // Every triangle's E[i] is the edge (V[i], V[(i+1)%3]) and points at the
    // object stored in the edge map.
    for (const tri of mesh.getTriangles()) {
        for (let i = 0; i < 3; ++i) {
            const e = tri.E[i];
            expect(e).not.toBeNull();
            const edge = e as ETNonmanifoldMeshEdge;
            expect(edgeKey(edge.V[0], edge.V[1]))
                .toBe(edgeKey(tri.V[i], tri.V[(i + 1) % 3]));
            expect(mesh.getEdge(tri.V[i], tri.V[(i + 1) % 3])).toBe(edge);
        }
        expect(mesh.getTriangle(tri.V[0], tri.V[1], tri.V[2])).toBe(tri);
    }

    const sizes = Array.from(modelEdges.values(), list => list.length);
    expect(mesh.isManifold()).toBe(sizes.every(s => s <= 2));
    expect(mesh.isClosed()).toBe(sizes.every(s => s === 2));
}

// Distinct-vertex triangles over a small pool, so the random meshes are
// densely connected and frequently nonmanifold. Degenerate triangles (two
// equal vertices) are excluded because Remove throws on them, which is the
// upstream defect covered by its own test above.
const triangleArb: fc.Arbitrary<Tri> =
    fc.uniqueArray(fc.integer({ min: 0, max: 5 }),
        { minLength: 3, maxLength: 3 }).map(a => [a[0], a[1], a[2]] as Tri);

const triangleListArb = fc.array(triangleArb, { minLength: 1, maxLength: 14 });

describe('ETNonmanifoldMesh verification', () => {
    it('keeps the edge-triangle adjacency in step with a brute-force model', () => {
        check(fc.tuple(triangleListArb, fc.array(fc.boolean(),
            { minLength: 14, maxLength: 14 })), ([triangles, removeFlags]) => {
            const mesh = new ETNonmanifoldMesh();
            const model = new MeshModel();
            const inserted: Tri[] = [];

            for (let i = 0; i < triangles.length; ++i) {
                const t = triangles[i];
                const created = mesh.insert(t[0], t[1], t[2]);
                const modelCreated = model.insert(t);
                expect(created !== null).toBe(modelCreated);
                if (created !== null) {
                    expect(triKey(created.V[0], created.V[1], created.V[2]))
                        .toBe(triKey(t[0], t[1], t[2]));
                    inserted.push(t);
                }
                expectMeshMatchesModel(mesh, model);

                // Interleave a removal of an earlier triangle.
                if (removeFlags[i] && inserted.length > 0) {
                    const victim = inserted.splice(i % inserted.length, 1)[0];
                    expect(mesh.remove(victim[0], victim[1], victim[2]))
                        .toBe(model.remove(victim));
                    expectMeshMatchesModel(mesh, model);
                    // Removing it again fails on both.
                    expect(mesh.remove(victim[0], victim[1], victim[2]))
                        .toBe(false);
                }
            }

            // Removing everything empties the mesh, edges included.
            for (const t of Array.from(model.triangles.values())) {
                expect(mesh.remove(t[0], t[1], t[2])).toBe(true);
                model.remove(t);
                expectMeshMatchesModel(mesh, model);
            }
            expect(mesh.getNumTriangles()).toBe(0);
            expect(mesh.getNumEdges()).toBe(0);
            expect(mesh.getEdges()).toEqual([]);
            expect(mesh.getTriangles()).toEqual([]);
        });
    });

    it('accepts every cyclic rotation as the same triangle and the reversal as another', () => {
        check(triangleArb, t => {
            const mesh = new ETNonmanifoldMesh();
            expect(mesh.insert(t[0], t[1], t[2])).not.toBeNull();
            expect(mesh.insert(t[1], t[2], t[0])).toBeNull();
            expect(mesh.insert(t[2], t[0], t[1])).toBeNull();
            expect(mesh.getNumTriangles()).toBe(1);
            // The reversal is a different ordered key.
            expect(mesh.insert(t[0], t[2], t[1])).not.toBeNull();
            expect(mesh.getNumTriangles()).toBe(2);
            // ... but shares all three unordered edges.
            expect(mesh.getNumEdges()).toBe(3);
            for (const edge of mesh.getEdges()) { expect(edge.T.size).toBe(2); }
            expect(mesh.isClosed()).toBe(true);
            expect(mesh.remove(t[1], t[2], t[0])).toBe(true);
            expect(mesh.getNumEdges()).toBe(3);
            expect(mesh.remove(t[2], t[1], t[0])).toBe(true);
            expect(mesh.getNumEdges()).toBe(0);
        });
    });

    it('computes the connected components of the edge-triangle graph', () => {
        check(triangleListArb, triangles => {
            const mesh = new ETNonmanifoldMesh();
            const model = new MeshModel();
            for (const t of triangles) {
                mesh.insert(t[0], t[1], t[2]);
                model.insert(t);
            }

            const components = mesh.getComponents().map(
                c => c.map(t => triKey(t.V[0], t.V[1], t.V[2])));
            const keyComponents = mesh.getComponentKeys().map(
                c => c.map(k => k.V.join(',')));
            expect(keyComponents).toEqual(components);

            // Every triangle appears exactly once.
            const all = components.flat();
            expect(all.length).toBe(mesh.getNumTriangles());
            expect(new Set(all).size).toBe(all.length);

            // The partition matches the brute-force union-find.
            const sortedGot = components.map(c => c.slice().sort(compareKeyStrings))
                .sort((a, b) => compareKeyStrings(a[0], b[0]));
            const sortedWant = model.components()
                .sort((a, b) => compareKeyStrings(a[0], b[0]));
            expect(sortedGot).toEqual(sortedWant);

            // Upstream seeds the components by walking mTMap in key order, so
            // the seeds appear in increasing key order.
            const order = meshTriangleKeys(mesh);
            const seeds = components.map(c => Math.min(
                ...c.map(k => order.indexOf(k))));
            for (let i = 1; i < seeds.length; ++i) {
                expect(seeds[i]).toBeGreaterThan(seeds[i - 1]);
            }
        });
    });

    it('deep-copies the triangles through the triangle keys', () => {
        check(triangleListArb, triangles => {
            const mesh = new ETNonmanifoldMesh();
            for (const t of triangles) { mesh.insert(t[0], t[1], t[2]); }

            const copy = mesh.clone();
            expect(meshTriangleKeys(copy)).toEqual(meshTriangleKeys(mesh));
            expect(copy.getNumEdges()).toBe(mesh.getNumEdges());
            // The copy owns its own objects.
            for (const tri of copy.getTriangles()) {
                expect(mesh.getTriangles()).not.toContain(tri);
                // Upstream reinserts using the key vertices, so a copied
                // triangle's V is the canonical rotation of the key.
                expect(tri.V.join(',')).toBe(triKey(tri.V[0], tri.V[1], tri.V[2]));
            }
            for (const edge of copy.getEdges()) {
                expect(mesh.getEdges()).not.toContain(edge);
            }
            // Clearing the source leaves the copy intact.
            const before = meshTriangleKeys(copy);
            mesh.clear();
            expect(mesh.getNumTriangles()).toBe(0);
            expect(meshTriangleKeys(copy)).toEqual(before);
        });
    });

    it('uses the supplied creators for every edge and triangle', () => {
        check(triangleListArb, triangles => {
            let numEdges = 0;
            let numTriangles = 0;
            const mesh = new ETNonmanifoldMesh(
                (v0, v1) => { ++numEdges; return new ETNonmanifoldMeshEdge(v0, v1); },
                (v0, v1, v2) => {
                    ++numTriangles;
                    return new ETNonmanifoldMeshTriangle(v0, v1, v2);
                });
            const model = new MeshModel();
            for (const t of triangles) {
                mesh.insert(t[0], t[1], t[2]);
                model.insert(t);
            }
            // One creator call per distinct triangle key and one per distinct
            // edge that had to be created.
            expect(numTriangles).toBe(model.triangles.size);
            expect(numEdges).toBe(mesh.getNumEdges());
        });
    });
});
