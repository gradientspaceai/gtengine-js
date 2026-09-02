import { describe, it, expect } from 'vitest';
import { VETManifoldMesh, VETManifoldMeshVertex } from '../src/VETManifoldMesh';
import { EdgeKey } from '../src/EdgeKey';
import { TriangleKey } from '../src/TriangleKey';

// A vertex with extra client data, used to exercise the VCreator callback.
class TaggedVertex extends VETManifoldMeshVertex {
    tag: string;

    constructor(vIndex: number) {
        super(vIndex);
        this.tag = 'v' + vIndex;
    }
}

function vIndices(mesh: VETManifoldMesh): number[] {
    return mesh.getVertices().map(v => v.V);
}

function edgeKeys(mesh: VETManifoldMesh, vIndex: number): string[] {
    const vertex = mesh.getVertex(vIndex);
    expect(vertex).not.toBeNull();
    return vertex!.getEAdjacent().map(
        e => new EdgeKey(false, e.V[0], e.V[1]).mapKey()).sort();
}

function triKeys(mesh: VETManifoldMesh, vIndex: number): string[] {
    const vertex = mesh.getVertex(vIndex);
    expect(vertex).not.toBeNull();
    return vertex!.getTAdjacent().map(
        t => new TriangleKey(true, t.V[0], t.V[1], t.V[2]).mapKey()).sort();
}

// The adjacency information a vertex must have, computed directly from the
// list of triangles that are in the mesh.
function expectedAdjacency(triangles: number[][]):
    Map<number, { v: number[], e: string[], t: string[] }> {
    const result = new Map<number, { v: Set<number>, e: Set<string>, t: Set<string> }>();
    const entry = (v: number) => {
        let item = result.get(v);
        if (item === undefined) {
            item = { v: new Set<number>(), e: new Set<string>(), t: new Set<string>() };
            result.set(v, item);
        }
        return item;
    };

    for (const tri of triangles) {
        const tkey = new TriangleKey(true, tri[0], tri[1], tri[2]).mapKey();
        for (let i = 0; i < 3; ++i) {
            const v = tri[i];
            const item = entry(v);
            item.t.add(tkey);
            for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
                const a = tri[j0];
                const b = tri[j1];
                if (a === v || b === v) {
                    item.v.add(a === v ? b : a);
                    item.e.add(new EdgeKey(false, a, b).mapKey());
                }
            }
        }
    }

    const sorted = new Map<number, { v: number[], e: string[], t: string[] }>();
    for (const [key, item] of result) {
        sorted.set(key, {
            v: Array.from(item.v).sort((a, b) => a - b),
            e: Array.from(item.e).sort(),
            t: Array.from(item.t).sort()
        });
    }
    return sorted;
}

function verifyAgainstTriangles(mesh: VETManifoldMesh, triangles: number[][]): void {
    const expected = expectedAdjacency(triangles);
    expect(vIndices(mesh)).toEqual(Array.from(expected.keys()).sort((a, b) => a - b));
    for (const [v, item] of expected) {
        const vertex = mesh.getVertex(v);
        expect(vertex).not.toBeNull();
        expect(vertex!.getVAdjacent()).toEqual(item.v);
        expect(edgeKeys(mesh, v)).toEqual(item.e);
        expect(triKeys(mesh, v)).toEqual(item.t);
    }
}

describe('VETManifoldMesh', () => {
    it('builds vertex adjacency for a single triangle', () => {
        const mesh = new VETManifoldMesh();
        const tri = mesh.insert(0, 1, 2);
        expect(tri).not.toBeNull();

        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.getNumTriangles()).toBe(1);
        expect(vIndices(mesh)).toEqual([0, 1, 2]);

        for (const v of [0, 1, 2]) {
            const vertex = mesh.getVertex(v)!;
            expect(vertex.V).toBe(v);
            expect(vertex.getVAdjacent()).toEqual([0, 1, 2].filter(u => u !== v));
            expect(vertex.EAdjacent.size).toBe(2);
            expect(vertex.getTAdjacent()).toEqual([tri]);
        }
        expect(mesh.getVertex(7)).toBeNull();
    });

    it('shares an interior edge between two triangles', () => {
        // The unit square split by the diagonal <0,2>.
        const mesh = new VETManifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        verifyAgainstTriangles(mesh, [[0, 1, 2], [0, 2, 3]]);

        // Vertex 0 and vertex 2 are on the shared diagonal.
        expect(mesh.getVertex(0)!.getVAdjacent()).toEqual([1, 2, 3]);
        expect(mesh.getVertex(0)!.TAdjacent.size).toBe(2);
        expect(mesh.getVertex(1)!.getVAdjacent()).toEqual([0, 2]);
        expect(mesh.getVertex(1)!.TAdjacent.size).toBe(1);
        expect(mesh.getVertex(3)!.getVAdjacent()).toEqual([0, 2]);
    });

    it('removes a triangle and keeps the shared edge adjacency', () => {
        const mesh = new VETManifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);

        expect(mesh.remove(0, 2, 3)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumEdges()).toBe(3);
        // Vertex 3 is no longer used by any triangle.
        expect(mesh.getVertex(3)).toBeNull();
        verifyAgainstTriangles(mesh, [[0, 1, 2]]);

        // The diagonal <0,2> survives because the remaining triangle uses it.
        expect(mesh.getVertex(0)!.getVAdjacent()).toEqual([1, 2]);

        expect(mesh.remove(0, 1, 2)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });

    it('reports failures for duplicate insertion and missing removal', () => {
        const mesh = new VETManifoldMesh();
        expect(mesh.insert(0, 1, 2)).not.toBeNull();
        // The rotated triple has the same ordered triangle key.
        expect(mesh.insert(1, 2, 0)).toBeNull();
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.remove(5, 6, 7)).toBe(false);
        expect(mesh.getNumTriangles()).toBe(1);
    });

    it('throws when an insertion would make the mesh nonmanifold', () => {
        const mesh = new VETManifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(1, 0, 3);
        // A third triangle on the edge <0,1>.
        expect(() => mesh.insert(1, 0, 4)).toThrow();
    });

    it('clears the vertices along with the edges and triangles', () => {
        const mesh = new VETManifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.clear();
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getVertices()).toEqual([]);
    });

    it('deep copies the mesh with clone and assign', () => {
        const mesh = new VETManifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);

        const copy = mesh.clone();
        expect(copy).not.toBe(mesh);
        verifyAgainstTriangles(copy, [[0, 1, 2], [0, 2, 3]]);
        // No object is shared between the meshes.
        expect(copy.getVertex(0)).not.toBe(mesh.getVertex(0));
        for (const tri of copy.getVertex(0)!.getTAdjacent()) {
            expect(mesh.getVertex(0)!.TAdjacent.has(tri)).toBe(false);
        }

        // Modifying the copy leaves the source alone.
        expect(copy.remove(0, 2, 3)).toBe(true);
        expect(copy.getNumTriangles()).toBe(1);
        expect(mesh.getNumTriangles()).toBe(2);

        const target = new VETManifoldMesh();
        target.insert(10, 11, 12);
        target.assign(mesh);
        verifyAgainstTriangles(target, [[0, 1, 2], [0, 2, 3]]);
        expect(target.getVertex(10)).toBeNull();
    });

    it('uses the vertex creation callback', () => {
        const mesh = new VETManifoldMesh(vIndex => new TaggedVertex(vIndex));
        mesh.insert(0, 1, 2);
        for (const vertex of mesh.getVertices()) {
            expect(vertex).toBeInstanceOf(TaggedVertex);
            expect((vertex as TaggedVertex).tag).toBe('v' + vertex.V);
        }
        // The creator is carried over by clone().
        const copy = mesh.clone();
        expect(copy.getVertex(1)).toBeInstanceOf(TaggedVertex);
    });

    it('maintains the adjacency invariants for a randomized grid mesh', () => {
        // A triangulated (n x n) grid of unit cells is an edge-manifold mesh.
        const n = 6;
        const index = (r: number, c: number) => r * (n + 1) + c;
        const triangles: number[][] = [];
        for (let r = 0; r < n; ++r) {
            for (let c = 0; c < n; ++c) {
                triangles.push([index(r, c), index(r, c + 1), index(r + 1, c + 1)]);
                triangles.push([index(r, c), index(r + 1, c + 1), index(r + 1, c)]);
            }
        }

        const mesh = new VETManifoldMesh();
        for (const tri of triangles) {
            expect(mesh.insert(tri[0], tri[1], tri[2])).not.toBeNull();
        }
        verifyAgainstTriangles(mesh, triangles);

        // Remove the triangles in a deterministic pseudorandom order and
        // verify the adjacency after every removal.
        let seed = 987654321;
        const nextRandom = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const remaining = triangles.slice();
        while (remaining.length > 0) {
            const i = Math.floor(nextRandom() * remaining.length);
            const tri = remaining[i];
            remaining.splice(i, 1);
            expect(mesh.remove(tri[0], tri[1], tri[2])).toBe(true);
            verifyAgainstTriangles(mesh, remaining);
        }

        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });
});
