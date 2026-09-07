import { describe, it, expect } from 'vitest';
import { VTSManifoldMesh, VTSManifoldMeshVertex } from '../src/VTSManifoldMesh.js';
import { TetrahedronKey } from '../src/TetrahedronKey.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V37): independent re-check against VTSManifoldMesh.h.
// ---------------------------------------------------------------------------

// The Kuhn (Freudenthal) subdivision of an (nx x ny x nz) grid of unit cubes
// into 6 tetrahedra per cube. The subdivision is face-conforming, so the
// complex is triangle-manifold and so is any subset of its tetrahedra: a face
// is shared by at most two of them, and an arbitrary insert/remove sequence
// over the subset never trips the nonmanifold rejection. Unlike the chain of
// tetrahedra used above, a vertex here is shared by many tetrahedra and a
// vertex pair by several faces, which is exactly the configuration in which
// upstream's Remove over-erases VAdjacent (upstream issue #256).
function kuhnTetrahedra(nx: number, ny: number, nz: number): number[][] {
    const index = (i: number, j: number, k: number) =>
        (i * (ny + 1) + j) * (nz + 1) + k;
    const permutations = [
        [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
    ];
    const tetrahedra: number[][] = [];
    for (let i = 0; i < nx; ++i) {
        for (let j = 0; j < ny; ++j) {
            for (let k = 0; k < nz; ++k) {
                for (const perm of permutations) {
                    const offset = [0, 0, 0];
                    const tetra = [index(i, j, k)];
                    for (const axis of perm) {
                        offset[axis] = 1;
                        tetra.push(index(i + (offset[0] as number),
                            j + (offset[1] as number),
                            k + (offset[2] as number)));
                    }
                    tetrahedra.push(tetra);
                }
            }
        }
    }
    return tetrahedra;
}

// The face records the base class must hold for a given set of tetrahedra:
// a face key maps to the sorted keys of the tetrahedra that contain it.
function expectedFaces(tetrahedra: number[][]): Map<string, string[]> {
    const faces = new Map<string, Set<string>>();
    for (const tetra of tetrahedra) {
        const skey = new TetrahedronKey(true, tetra[0] as number,
            tetra[1] as number, tetra[2] as number,
            tetra[3] as number).mapKey();
        for (let j = 0; j < 4; ++j) {
            const face = tetra.filter((_, k) => k !== j);
            const tkey = new TriangleKey(false, face[0] as number,
                face[1] as number, face[2] as number).mapKey();
            let s = faces.get(tkey);
            if (s === undefined) {
                s = new Set<string>();
                faces.set(tkey, s);
            }
            s.add(skey);
        }
    }
    const sorted = new Map<string, string[]>();
    for (const [k, s] of faces) {
        sorted.set(k, Array.from(s).sort());
    }
    return sorted;
}

// verifyAgainstTetrahedra plus the base-class face bookkeeping and the
// counts, so that the whole mesh state is pinned.
function verifyFullState(mesh: VTSManifoldMesh, tetrahedra: number[][]): void {
    verifyAgainstTetrahedra(mesh, tetrahedra);

    const faces = expectedFaces(tetrahedra);
    expect(mesh.getNumTetrahedra()).toBe(tetrahedra.length);
    expect(mesh.getNumTriangles()).toBe(faces.size);
    expect(mesh.getTriangles().map(
        t => new TriangleKey(false, t.V[0], t.V[1], t.V[2]).mapKey()).sort())
        .toEqual(Array.from(faces.keys()).sort());

    for (const face of mesh.getTriangles()) {
        const tkey = new TriangleKey(false, face.V[0], face.V[1],
            face.V[2]).mapKey();
        const tetras: string[] = [];
        for (const s of face.S) {
            if (s !== null) {
                tetras.push(new TetrahedronKey(true, s.V[0], s.V[1], s.V[2],
                    s.V[3]).mapKey());
            }
        }
        tetras.sort();
        expect(tetras).toEqual(faces.get(tkey));
        // A one-tetrahedron face always keeps its reference at index 0.
        if (face.S[1] !== null) {
            expect(face.S[0]).not.toBeNull();
        }
    }

    // The vertex adjacency relation is symmetric, and every face a vertex
    // references is still in the mesh.
    const liveFaces = new Set<string>(faces.keys());
    for (const vertex of mesh.getVertices()) {
        for (const w of vertex.getVAdjacent()) {
            const other = mesh.getVertex(w) as VTSManifoldMeshVertex;
            expect(other).not.toBeNull();
            expect(other.getVAdjacent()).toContain(vertex.V);
        }
        for (const face of vertex.getTAdjacent()) {
            expect(liveFaces.has(new TriangleKey(false, face.V[0], face.V[1],
                face.V[2]).mapKey())).toBe(true);
        }
    }
}

describe('VTSManifoldMesh verification', () => {
    const tetrahedra = kuhnTetrahedra(2, 1, 1);

    it('matches a brute-force adjacency model under interleaved insert/remove', () => {
        const op = fc.record({
            doInsert: fc.boolean(),
            s: fc.integer({ min: 0, max: tetrahedra.length - 1 })
        });
        check(fc.array(op, { minLength: 1, maxLength: 24 }), ops => {
            const mesh = new VTSManifoldMesh();
            const present = new Set<number>();
            for (const step of ops) {
                const s = tetrahedra[step.s] as number[];
                if (step.doInsert) {
                    const result = mesh.insert(s[0] as number, s[1] as number,
                        s[2] as number, s[3] as number);
                    if (present.has(step.s)) {
                        // Upstream returns nullptr for a tetrahedron already
                        // in the mesh and leaves the mesh unchanged.
                        expect(result).toBeNull();
                    }
                    else {
                        expect(result).not.toBeNull();
                        present.add(step.s);
                    }
                }
                else {
                    const expected = present.delete(step.s);
                    expect(mesh.remove(s[0] as number, s[1] as number,
                        s[2] as number, s[3] as number)).toBe(expected);
                }
                verifyFullState(mesh,
                    Array.from(present).map(i => tetrahedra[i] as number[]));
            }

            // Removing everything empties every map.
            for (const i of Array.from(present)) {
                const s = tetrahedra[i] as number[];
                expect(mesh.remove(s[0] as number, s[1] as number,
                    s[2] as number, s[3] as number)).toBe(true);
            }
            expect(mesh.getNumVertices()).toBe(0);
            expect(mesh.getNumTriangles()).toBe(0);
            expect(mesh.getNumTetrahedra()).toBe(0);
            expect(mesh.getVertices()).toEqual([]);
        }, 60);
    }, 30000);

    it('rejects a nonmanifold insertion without disturbing the vertex layer', () => {
        check(fc.integer({ min: 0, max: 6 }), k => {
            // Three tetrahedra sharing the face <0,1,2>; the third insertion
            // is nonmanifold.
            const mesh = new VTSManifoldMesh();
            expect(mesh.insert(0, 1, 2, 3)).not.toBeNull();
            expect(mesh.insert(0, 1, 2, 4)).not.toBeNull();
            const good = [[0, 1, 2, 3], [0, 1, 2, 4]];
            verifyFullState(mesh, good);

            // Upstream runs the base-class Insert first and only updates
            // mVMap when it succeeds, so the vertex layer is untouched by a
            // rejection. (The base class leaves behind the faces it created
            // before the failure; that ETManifoldMesh/TSManifoldMesh quirk is
            // preserved by the port and is not part of the vertex layer.)
            const rejected = new TetrahedronKey(true, 0, 1, 2, 5 + k).mapKey();
            const expectVertexLayerUnchanged = (): void => {
                expect(mesh.getNumTetrahedra()).toBe(2);
                expect(mesh.getTetrahedron(0, 1, 2, 5 + k)).toBeNull();
                verifyAgainstTetrahedra(mesh, good);
                for (const vertex of mesh.getVertices()) {
                    const keys = vertex.getSAdjacent().map(
                        s => new TetrahedronKey(true, s.V[0], s.V[1], s.V[2],
                            s.V[3]).mapKey());
                    expect(keys).not.toContain(rejected);
                }
            };

            expect(() => mesh.insert(0, 1, 2, 5 + k)).toThrow();
            expectVertexLayerUnchanged();

            const previous = mesh.throwOnNonmanifoldInsertion(false);
            expect(previous).toBe(true);
            expect(mesh.insert(0, 1, 2, 5 + k)).toBeNull();
            mesh.throwOnNonmanifoldInsertion(previous);
            expectVertexLayerUnchanged();
        }, 7);
    });

    it('clone and assign reproduce the adjacency model', () => {
        const chosenArb = fc.uniqueArray(
            fc.integer({ min: 0, max: tetrahedra.length - 1 }),
            { minLength: 1, maxLength: tetrahedra.length });
        check(chosenArb, chosen => {
            const mesh = new VTSManifoldMesh();
            const chosenTetrahedra = chosen.map(i => tetrahedra[i] as number[]);
            for (const s of chosenTetrahedra) {
                expect(mesh.insert(s[0] as number, s[1] as number,
                    s[2] as number, s[3] as number)).not.toBeNull();
            }

            const copy = mesh.clone();
            verifyFullState(copy, chosenTetrahedra);
            for (const vertex of copy.getVertices()) {
                expect(vertex).not.toBe(mesh.getVertex(vertex.V));
            }

            // Removing from the copy leaves the original alone.
            const first = chosenTetrahedra[0] as number[];
            expect(copy.remove(first[0] as number, first[1] as number,
                first[2] as number, first[3] as number)).toBe(true);
            verifyFullState(mesh, chosenTetrahedra);

            // assign() clears the target first.
            const target = new VTSManifoldMesh();
            expect(target.insert(100, 101, 102, 103)).not.toBeNull();
            target.assign(mesh);
            verifyFullState(target, chosenTetrahedra);
        }, 40);
    }, 30000);
});
