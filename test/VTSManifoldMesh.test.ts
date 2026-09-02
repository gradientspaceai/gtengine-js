import { describe, it, expect } from 'vitest';
import { VTSManifoldMesh, VTSManifoldMeshVertex } from '../src/VTSManifoldMesh';
import { TetrahedronKey } from '../src/TetrahedronKey';
import { TriangleKey } from '../src/TriangleKey';

// A vertex with extra client data, used to exercise the VCreator callback.
class TaggedVertex extends VTSManifoldMeshVertex {
    tag: string;

    constructor(vIndex: number) {
        super(vIndex);
        this.tag = 'v' + vIndex;
    }
}

function vIndices(mesh: VTSManifoldMesh): number[] {
    return mesh.getVertices().map(v => v.V);
}

// The adjacency information a vertex must have, computed directly from the
// list of tetrahedra that are in the mesh. The faces of <v0,v1,v2,v3> are the
// four triples of distinct vertices.
function expectedAdjacency(tetrahedra: number[][]):
    Map<number, { v: number[], t: string[], s: string[] }> {
    const result = new Map<number, { v: Set<number>, t: Set<string>, s: Set<string> }>();
    const entry = (v: number) => {
        let item = result.get(v);
        if (item === undefined) {
            item = { v: new Set<number>(), t: new Set<string>(), s: new Set<string>() };
            result.set(v, item);
        }
        return item;
    };

    for (const tetra of tetrahedra) {
        const skey = new TetrahedronKey(true, tetra[0], tetra[1], tetra[2],
            tetra[3]).mapKey();
        for (let i = 0; i < 4; ++i) {
            const v = tetra[i];
            const item = entry(v);
            item.s.add(skey);
            for (let j = 0; j < 4; ++j) {
                // The face opposite tetra[j].
                const face = tetra.filter((_, k) => k !== j);
                if (face.indexOf(v) >= 0) {
                    item.t.add(new TriangleKey(false, face[0], face[1], face[2]).mapKey());
                    for (const u of face) {
                        if (u !== v) {
                            item.v.add(u);
                        }
                    }
                }
            }
        }
    }

    const sorted = new Map<number, { v: number[], t: string[], s: string[] }>();
    for (const [key, item] of result) {
        sorted.set(key, {
            v: Array.from(item.v).sort((a, b) => a - b),
            t: Array.from(item.t).sort(),
            s: Array.from(item.s).sort()
        });
    }
    return sorted;
}

function verifyAgainstTetrahedra(mesh: VTSManifoldMesh, tetrahedra: number[][]): void {
    const expected = expectedAdjacency(tetrahedra);
    expect(vIndices(mesh)).toEqual(Array.from(expected.keys()).sort((a, b) => a - b));
    for (const [v, item] of expected) {
        const vertex = mesh.getVertex(v);
        expect(vertex).not.toBeNull();
        expect(vertex!.getVAdjacent()).toEqual(item.v);
        expect(vertex!.getTAdjacent().map(
            t => new TriangleKey(false, t.V[0], t.V[1], t.V[2]).mapKey()).sort())
            .toEqual(item.t);
        expect(vertex!.getSAdjacent().map(
            s => new TetrahedronKey(true, s.V[0], s.V[1], s.V[2], s.V[3]).mapKey()).sort())
            .toEqual(item.s);
    }
}

describe('VTSManifoldMesh', () => {
    it('builds vertex adjacency for a single tetrahedron', () => {
        const mesh = new VTSManifoldMesh();
        const tetra = mesh.insert(0, 1, 2, 3);
        expect(tetra).not.toBeNull();

        expect(mesh.getNumVertices()).toBe(4);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(vIndices(mesh)).toEqual([0, 1, 2, 3]);

        for (const v of [0, 1, 2, 3]) {
            const vertex = mesh.getVertex(v)!;
            expect(vertex.V).toBe(v);
            expect(vertex.getVAdjacent()).toEqual([0, 1, 2, 3].filter(u => u !== v));
            // The vertex is on three of the four faces.
            expect(vertex.TAdjacent.size).toBe(3);
            expect(vertex.getSAdjacent()).toEqual([tetra]);
        }
        expect(mesh.getVertex(9)).toBeNull();
        verifyAgainstTetrahedra(mesh, [[0, 1, 2, 3]]);
    });

    it('shares an interior face between two tetrahedra', () => {
        const mesh = new VTSManifoldMesh();
        mesh.insert(0, 1, 2, 3);
        mesh.insert(0, 1, 2, 4);
        expect(mesh.getNumTetrahedra()).toBe(2);
        // Seven faces: the shared face plus three per tetrahedron.
        expect(mesh.getNumTriangles()).toBe(7);
        verifyAgainstTetrahedra(mesh, [[0, 1, 2, 3], [0, 1, 2, 4]]);

        expect(mesh.getVertex(0)!.getVAdjacent()).toEqual([1, 2, 3, 4]);
        expect(mesh.getVertex(3)!.getVAdjacent()).toEqual([0, 1, 2]);
    });

    it('keeps adjacency supported by surviving faces when a tetrahedron is removed', () => {
        // This is the upstream bug that the port fixes: upstream erases the
        // other two vertices of every destroyed face from VAdjacent, which
        // drops adjacencies that the surviving shared face still supports.
        const mesh = new VTSManifoldMesh();
        mesh.insert(0, 1, 2, 3);
        mesh.insert(0, 1, 2, 4);

        expect(mesh.remove(0, 1, 2, 3)).toBe(true);
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.getVertex(3)).toBeNull();

        // Upstream would report [4] here, losing 1 and 2.
        expect(mesh.getVertex(0)!.getVAdjacent()).toEqual([1, 2, 4]);
        verifyAgainstTetrahedra(mesh, [[0, 1, 2, 4]]);

        expect(mesh.remove(0, 1, 2, 4)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumTetrahedra()).toBe(0);
    });

    it('reports failures for duplicate insertion and missing removal', () => {
        const mesh = new VTSManifoldMesh();
        expect(mesh.insert(0, 1, 2, 3)).not.toBeNull();
        // An even permutation has the same ordered tetrahedron key.
        expect(mesh.insert(1, 0, 3, 2)).toBeNull();
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(mesh.remove(5, 6, 7, 8)).toBe(false);
    });

    it('throws when an insertion would make the mesh nonmanifold', () => {
        const mesh = new VTSManifoldMesh();
        mesh.insert(0, 1, 2, 3);
        mesh.insert(0, 1, 2, 4);
        // A third tetrahedron on the face <0,1,2>.
        expect(() => mesh.insert(0, 1, 2, 5)).toThrow();
    });

    it('clears the vertices along with the faces and tetrahedra', () => {
        const mesh = new VTSManifoldMesh();
        mesh.insert(0, 1, 2, 3);
        mesh.insert(0, 1, 2, 4);
        mesh.clear();
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumTetrahedra()).toBe(0);
        expect(mesh.getVertices()).toEqual([]);
    });

    it('deep copies the mesh with clone and assign', () => {
        const mesh = new VTSManifoldMesh();
        mesh.insert(0, 1, 2, 3);
        mesh.insert(0, 1, 2, 4);

        const copy = mesh.clone();
        expect(copy).not.toBe(mesh);
        verifyAgainstTetrahedra(copy, [[0, 1, 2, 3], [0, 1, 2, 4]]);
        expect(copy.getVertex(0)).not.toBe(mesh.getVertex(0));
        for (const tetra of copy.getVertex(0)!.getSAdjacent()) {
            expect(mesh.getVertex(0)!.SAdjacent.has(tetra)).toBe(false);
        }

        expect(copy.remove(0, 1, 2, 3)).toBe(true);
        expect(copy.getNumTetrahedra()).toBe(1);
        expect(mesh.getNumTetrahedra()).toBe(2);

        const target = new VTSManifoldMesh();
        target.insert(10, 11, 12, 13);
        target.assign(mesh);
        verifyAgainstTetrahedra(target, [[0, 1, 2, 3], [0, 1, 2, 4]]);
        expect(target.getVertex(10)).toBeNull();
    });

    it('uses the vertex creation callback', () => {
        const mesh = new VTSManifoldMesh(vIndex => new TaggedVertex(vIndex));
        mesh.insert(0, 1, 2, 3);
        for (const vertex of mesh.getVertices()) {
            expect(vertex).toBeInstanceOf(TaggedVertex);
            expect((vertex as TaggedVertex).tag).toBe('v' + vertex.V);
        }
        const copy = mesh.clone();
        expect(copy.getVertex(1)).toBeInstanceOf(TaggedVertex);
    });

    it('maintains the adjacency invariants for a randomized chain of tetrahedra', () => {
        // The tetrahedra <i,i+1,i+2,i+3> form a triangle-manifold chain: the
        // face <i+1,i+2,i+3> is shared by consecutive tetrahedra and by no
        // others.
        const numTetrahedra = 12;
        const tetrahedra: number[][] = [];
        for (let i = 0; i < numTetrahedra; ++i) {
            tetrahedra.push([i, i + 1, i + 2, i + 3]);
        }

        const mesh = new VTSManifoldMesh();
        for (const s of tetrahedra) {
            expect(mesh.insert(s[0], s[1], s[2], s[3])).not.toBeNull();
        }
        verifyAgainstTetrahedra(mesh, tetrahedra);

        // Remove the tetrahedra in a deterministic pseudorandom order and
        // verify the adjacency after every removal.
        let seed = 24680135;
        const nextRandom = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const remaining = tetrahedra.slice();
        while (remaining.length > 0) {
            const i = Math.floor(nextRandom() * remaining.length);
            const s = remaining[i];
            remaining.splice(i, 1);
            expect(mesh.remove(s[0], s[1], s[2], s[3])).toBe(true);
            verifyAgainstTetrahedra(mesh, remaining);
        }

        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getNumTetrahedra()).toBe(0);
    });
});
