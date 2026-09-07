import { describe, it, expect } from 'vitest';
import {
    TSManifoldMesh,
    TSManifoldMeshTriangle,
    TSManifoldMeshTetrahedron
} from '../src/TSManifoldMesh.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { TetrahedronKey } from '../src/TetrahedronKey.js';
import { check, fc } from './helpers/arbitraries.js';

// The unordered vertex triples of the faces, in mesh iteration order.
function faceTriples(mesh: TSManifoldMesh): number[][] {
    return mesh.getTriangles().map(
        f => new TriangleKey(false, f.V[0], f.V[1], f.V[2]).V.slice());
}

// The ordered tetrahedron keys, in mesh iteration order.
function tetraKeys(mesh: TSManifoldMesh): number[][] {
    return mesh.getTetrahedra().map(
        s => new TetrahedronKey(true, s.V[0], s.V[1], s.V[2], s.V[3]).V.slice());
}

// Verify the internal consistency of the mesh: every face references only
// tetrahedra that are in the mesh, every tetrahedron references the four
// faces opposite its vertices, and the tetrahedron-tetrahedron adjacency
// agrees with the face-tetrahedron adjacency.
function checkInvariants(mesh: TSManifoldMesh): void {
    const tetrahedra = new Set(mesh.getTetrahedra());
    const faces = new Set(mesh.getTriangles());
    const opposite = TetrahedronKey.getOppositeFace();

    for (const face of mesh.getTriangles()) {
        expect(face.S[0]).not.toBeNull();
        for (const tetra of face.S) {
            if (tetra !== null) {
                expect(tetrahedra.has(tetra)).toBe(true);
                expect(tetra.T.indexOf(face)).toBeGreaterThanOrEqual(0);
            }
        }
        if (face.S[0] !== null && face.S[1] !== null) {
            expect(face.S[0]).not.toBe(face.S[1]);
        }
    }

    for (const tetra of mesh.getTetrahedra()) {
        expect(mesh.getTetrahedron(tetra.V[0], tetra.V[1], tetra.V[2], tetra.V[3]))
            .toBe(tetra);
        for (let i = 0; i < 4; ++i) {
            const face = tetra.T[i];
            expect(face).not.toBeNull();
            const f = face as TSManifoldMeshTriangle;
            expect(faces.has(f)).toBe(true);

            // T[i] is the face opposite V[i].
            const o = opposite[i];
            expect(new TriangleKey(false, f.V[0], f.V[1], f.V[2]).mapKey()).toBe(
                new TriangleKey(false, tetra.V[o[0]], tetra.V[o[1]],
                    tetra.V[o[2]]).mapKey());
            expect(mesh.getTriangle(f.V[0], f.V[1], f.V[2])).toBe(f);

            // S[i] is the other tetrahedron sharing T[i].
            const other = f.S[0] === tetra ? f.S[1] : f.S[0];
            expect(tetra.S[i]).toBe(other);
            if (tetra.S[i] !== null) {
                expect((tetra.S[i] as TSManifoldMeshTetrahedron).S.indexOf(tetra))
                    .toBeGreaterThanOrEqual(0);
            }
        }
    }
}

function build(tets: number[][], throwOnNonmanifold = true): TSManifoldMesh {
    const mesh = new TSManifoldMesh();
    mesh.throwOnNonmanifoldInsertion(throwOnNonmanifold);
    for (const s of tets) {
        mesh.insert(s[0], s[1], s[2], s[3]);
    }
    return mesh;
}

// The boundary of the 4-simplex: the five 4-subsets of {0,1,2,3,4}. Every
// 3-subset lies in exactly two of them, so the mesh is closed.
const fiveCell: number[][] = [
    [0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 3, 4], [0, 2, 3, 4], [1, 2, 3, 4]
];

describe('TSManifoldMesh insertion', () => {
    it('builds a mesh with a single tetrahedron', () => {
        const mesh = build([[0, 1, 2, 3]]);
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(mesh.getNumTriangles()).toBe(4);
        expect(mesh.isClosed()).toBe(false);

        const tetra = mesh.getTetrahedron(0, 1, 2, 3) as TSManifoldMeshTetrahedron;
        expect(tetra).not.toBeNull();
        expect(tetra.V).toEqual([0, 1, 2, 3]);
        expect(tetra.S).toEqual([null, null, null, null]);

        // The four faces are those opposite the four vertices, each with a
        // single tetrahedron.
        expect(faceTriples(mesh)).toEqual([
            [0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]
        ]);
        for (const face of mesh.getTriangles()) {
            expect(face.S[0]).toBe(tetra);
            expect(face.S[1]).toBeNull();
        }
        checkInvariants(mesh);
    });

    it('rejects a duplicate tetrahedron by returning null', () => {
        const mesh = build([[0, 1, 2, 3]]);
        expect(mesh.insert(0, 1, 2, 3)).toBeNull();
        // (1,0,3,2) is an even permutation of (0,1,2,3), so it has the same
        // ordered key and names the same tetrahedron.
        expect(mesh.insert(1, 0, 3, 2)).toBeNull();
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(mesh.getNumTriangles()).toBe(4);
    });

    it('finds features by any equivalent vertex ordering', () => {
        const mesh = build([[0, 1, 2, 3]]);
        const tetra = mesh.getTetrahedron(0, 1, 2, 3);
        expect(mesh.getTetrahedron(1, 0, 3, 2)).toBe(tetra);
        expect(mesh.getTetrahedron(2, 3, 0, 1)).toBe(tetra);
        // The face key is unordered, so any permutation names the face.
        const face = mesh.getTriangle(0, 1, 2);
        expect(face).not.toBeNull();
        expect(mesh.getTriangle(2, 1, 0)).toBe(face);
        expect(mesh.getTriangle(1, 2, 0)).toBe(face);
        expect(mesh.getTriangle(0, 1, 5)).toBeNull();
        expect(mesh.getTetrahedron(0, 1, 2, 5)).toBeNull();
    });

    it('links two tetrahedra that share a face', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4]]);
        expect(mesh.getNumTetrahedra()).toBe(2);
        expect(mesh.getNumTriangles()).toBe(7);
        expect(mesh.isClosed()).toBe(false);

        const a = mesh.getTetrahedron(0, 1, 2, 3) as TSManifoldMeshTetrahedron;
        const b = mesh.getTetrahedron(0, 1, 2, 4) as TSManifoldMeshTetrahedron;
        const shared = mesh.getTriangle(0, 1, 2) as TSManifoldMeshTriangle;
        expect(shared.S[0]).toBe(a);
        expect(shared.S[1]).toBe(b);

        // The shared face is opposite V[3] in both tetrahedra.
        expect(a.T[3]).toBe(shared);
        expect(b.T[3]).toBe(shared);
        expect(a.S).toEqual([null, null, null, b]);
        expect(b.S).toEqual([null, null, null, a]);
        checkInvariants(mesh);
    });

    it('links three tetrahedra around a shared edge', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 3, 4], [0, 1, 4, 2]]);
        expect(mesh.getNumTetrahedra()).toBe(3);
        // 12 face slots, 3 of the faces are shared: 12 - 3 = 9 faces.
        expect(mesh.getNumTriangles()).toBe(9);
        expect(mesh.isClosed()).toBe(false);

        const a = mesh.getTetrahedron(0, 1, 2, 3) as TSManifoldMeshTetrahedron;
        const b = mesh.getTetrahedron(0, 1, 3, 4) as TSManifoldMeshTetrahedron;
        const c = mesh.getTetrahedron(0, 1, 4, 2) as TSManifoldMeshTetrahedron;

        // Each tetrahedron is adjacent to the other two; the three faces
        // containing the shared edge {0,1} each have two tetrahedra.
        expect(new Set(a.S.filter(s => s !== null))).toEqual(new Set([b, c]));
        expect(new Set(b.S.filter(s => s !== null))).toEqual(new Set([a, c]));
        expect(new Set(c.S.filter(s => s !== null))).toEqual(new Set([a, b]));

        for (const triple of [[0, 1, 2], [0, 1, 3], [0, 1, 4]]) {
            const face = mesh.getTriangle(triple[0], triple[1], triple[2]);
            expect(face).not.toBeNull();
            expect((face as TSManifoldMeshTriangle).S[1]).not.toBeNull();
        }
        // The faces not containing the shared edge have a single tetrahedron.
        for (const triple of [[0, 2, 3], [1, 2, 3], [0, 3, 4], [1, 3, 4],
            [0, 2, 4], [1, 2, 4]]) {
            const face = mesh.getTriangle(triple[0], triple[1], triple[2]);
            expect(face).not.toBeNull();
            expect((face as TSManifoldMeshTriangle).S[1]).toBeNull();
        }
        checkInvariants(mesh);
    });

    it('builds a closed mesh from the boundary of the 4-simplex', () => {
        const mesh = build(fiveCell);
        expect(mesh.getNumTetrahedra()).toBe(5);
        // Every 3-subset of {0,...,4} is a face and lies in two tetrahedra.
        expect(mesh.getNumTriangles()).toBe(10);
        expect(mesh.isClosed()).toBe(true);
        for (const face of mesh.getTriangles()) {
            expect(face.S[0]).not.toBeNull();
            expect(face.S[1]).not.toBeNull();
        }
        for (const tetra of mesh.getTetrahedra()) {
            expect(tetra.S.every(s => s !== null)).toBe(true);
        }
        expect(tetraKeys(mesh)).toEqual(fiveCell);
        checkInvariants(mesh);
    });

    it('throws when an insertion would make the mesh nonmanifold', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4]]);
        // The face {0,1,2} already has two tetrahedra.
        expect(() => mesh.insert(0, 1, 2, 5)).toThrow('Attempt to create nonmanifold mesh.');
    });

    // Upstream quirk (the TSManifoldMesh form of upstream issue #73 for
    // VEManifoldMesh::Insert): the graceful rejection path returns null but
    // leaves the faces created earlier in the loop in the face map, each
    // referencing a tetrahedron that was never added to the tetrahedron map.
    it('leaves the rejected insertion partially in the mesh when not throwing', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4]]);
        expect(mesh.getNumTriangles()).toBe(7);

        const previous = mesh.throwOnNonmanifoldInsertion(false);
        expect(previous).toBe(true);

        // The loop creates the faces {1,2,5}, {0,2,5} and {0,1,5} before it
        // reaches the full face {0,1,2} and returns null.
        expect(mesh.insert(0, 1, 2, 5)).toBeNull();
        expect(mesh.getNumTetrahedra()).toBe(2);
        expect(mesh.getNumTriangles()).toBe(10);

        const leaked = mesh.getTriangle(1, 2, 5) as TSManifoldMeshTriangle;
        expect(leaked).not.toBeNull();
        const phantom = leaked.S[0] as TSManifoldMeshTetrahedron;
        expect(phantom.V).toEqual([0, 1, 2, 5]);
        // The phantom tetrahedron is not in the tetrahedron map.
        expect(mesh.getTetrahedron(0, 1, 2, 5)).toBeNull();
        expect(mesh.getTetrahedra()).toHaveLength(2);

        expect(mesh.throwOnNonmanifoldInsertion(true)).toBe(false);
    });

    it('uses the triangle and tetrahedron creator callbacks', () => {
        class MyTriangle extends TSManifoldMeshTriangle {
            label = 'face';
        }
        class MyTetrahedron extends TSManifoldMeshTetrahedron {
            label = 'tetrahedron';
        }
        const mesh = new TSManifoldMesh(
            (v0, v1, v2) => new MyTriangle(v0, v1, v2),
            (v0, v1, v2, v3) => new MyTetrahedron(v0, v1, v2, v3));
        mesh.insert(0, 1, 2, 3);
        expect(mesh.getTriangles().every(f => f instanceof MyTriangle)).toBe(true);
        expect(mesh.getTetrahedra().every(s => s instanceof MyTetrahedron)).toBe(true);

        const copy = mesh.clone();
        copy.insert(0, 1, 2, 4);
        expect(copy.getTetrahedra().every(s => s instanceof MyTetrahedron)).toBe(true);
        expect(copy.getTriangles().every(f => f instanceof MyTriangle)).toBe(true);
    });
});

describe('TSManifoldMesh removal', () => {
    it('removes a tetrahedron and restores the adjacency of its neighbor', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4]]);
        expect(mesh.remove(0, 1, 2, 4)).toBe(true);
        expect(mesh.getNumTetrahedra()).toBe(1);
        expect(mesh.getNumTriangles()).toBe(4);

        const a = mesh.getTetrahedron(0, 1, 2, 3) as TSManifoldMeshTetrahedron;
        expect(a.S).toEqual([null, null, null, null]);
        // The shared face keeps its remaining tetrahedron at index zero.
        const shared = mesh.getTriangle(0, 1, 2) as TSManifoldMeshTriangle;
        expect(shared.S[0]).toBe(a);
        expect(shared.S[1]).toBeNull();
        // The faces of the removed tetrahedron that were not shared are gone.
        expect(mesh.getTriangle(0, 1, 4)).toBeNull();
        checkInvariants(mesh);
    });

    it('returns false for a tetrahedron that is not in the mesh', () => {
        const mesh = build([[0, 1, 2, 3]]);
        expect(mesh.remove(0, 1, 2, 4)).toBe(false);
        // An even permutation names the same tetrahedron.
        expect(mesh.remove(1, 0, 3, 2)).toBe(true);
        expect(mesh.remove(0, 1, 2, 3)).toBe(false);
        expect(mesh.getNumTetrahedra()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });

    it('empties the closed mesh when all tetrahedra are removed', () => {
        const mesh = build(fiveCell);
        for (const s of fiveCell) {
            expect(mesh.remove(s[0], s[1], s[2], s[3])).toBe(true);
            checkInvariants(mesh);
        }
        expect(mesh.getNumTetrahedra()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        // A mesh with no faces is vacuously closed.
        expect(mesh.isClosed()).toBe(true);
    });

    it('supports reinsertion after removal', () => {
        const mesh = build(fiveCell);
        expect(mesh.remove(1, 2, 3, 4)).toBe(true);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.insert(1, 2, 3, 4)).not.toBeNull();
        expect(mesh.getNumTetrahedra()).toBe(5);
        expect(mesh.getNumTriangles()).toBe(10);
        expect(mesh.isClosed()).toBe(true);
        checkInvariants(mesh);
    });

    it('clears the mesh', () => {
        const mesh = build(fiveCell);
        mesh.clear();
        expect(mesh.getNumTetrahedra()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getTetrahedra()).toEqual([]);
        expect(mesh.getTriangles()).toEqual([]);
        expect(mesh.isClosed()).toBe(true);
    });

    it('removes the middle tetrahedron of a chain', () => {
        // Three tetrahedra in a row: A-B-C, with B sharing a face with each.
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4], [0, 1, 4, 5]]);
        expect(mesh.getNumTriangles()).toBe(10);
        expect(mesh.remove(0, 1, 2, 4)).toBe(true);
        expect(mesh.getNumTetrahedra()).toBe(2);
        // The two remaining tetrahedra are no longer adjacent.
        for (const tetra of mesh.getTetrahedra()) {
            expect(tetra.S).toEqual([null, null, null, null]);
        }
        expect(mesh.getNumTriangles()).toBe(8);
        checkInvariants(mesh);
    });
});

describe('TSManifoldMesh copying', () => {
    it('clones a mesh so that the copies are independent', () => {
        const mesh = build(fiveCell);
        const copy = mesh.clone();

        expect(copy.getNumTetrahedra()).toBe(5);
        expect(copy.getNumTriangles()).toBe(10);
        expect(tetraKeys(copy)).toEqual(tetraKeys(mesh));
        expect(faceTriples(copy)).toEqual(faceTriples(mesh));
        expect(copy.isClosed()).toBe(true);
        checkInvariants(copy);

        // No object is shared between the meshes.
        const originalTetrahedra = new Set<TSManifoldMeshTetrahedron>(mesh.getTetrahedra());
        expect(copy.getTetrahedra().some(s => originalTetrahedra.has(s))).toBe(false);
        const originalFaces = new Set<TSManifoldMeshTriangle>(mesh.getTriangles());
        expect(copy.getTriangles().some(f => originalFaces.has(f))).toBe(false);

        // Modifying the copy does not affect the original.
        expect(copy.remove(1, 2, 3, 4)).toBe(true);
        expect(copy.getNumTetrahedra()).toBe(4);
        expect(copy.isClosed()).toBe(false);
        expect(mesh.getNumTetrahedra()).toBe(5);
        expect(mesh.isClosed()).toBe(true);
    });

    it('reinserts using the tetrahedron keys, so the vertices may be permuted', () => {
        // The upstream operator= reinserts element.first.V[], the ordered
        // tetrahedron key, not the tetrahedron's own vertex order.
        const mesh = build([[3, 1, 0, 2]]);
        expect(mesh.getTetrahedra()[0].V).toEqual([3, 1, 0, 2]);
        const copy = mesh.clone();
        expect(copy.getTetrahedra()[0].V).toEqual(
            new TetrahedronKey(true, 3, 1, 0, 2).V.slice());
    });

    it('copies the nonmanifold-insertion state and replaces the target', () => {
        const mesh = build([[0, 1, 2, 3], [0, 1, 2, 4]]);
        mesh.throwOnNonmanifoldInsertion(false);
        const copy = mesh.clone();
        // The copy does not throw either.
        expect(copy.insert(0, 1, 2, 5)).toBeNull();

        const target = build([[10, 11, 12, 13]]);
        target.assign(mesh);
        expect(target.getNumTetrahedra()).toBe(2);
        expect(target.getTetrahedron(10, 11, 12, 13)).toBeNull();
        expect(tetraKeys(target)).toEqual(tetraKeys(mesh));
    });
});

describe('TSManifoldMesh randomized consistency', () => {
    // A deterministic pseudorandom generator so the test is reproducible.
    function makeRandom(seed: number): () => number {
        let state = seed >>> 0;
        return () => {
            state = (state * 1664525 + 1013904223) >>> 0;
            return state / 4294967296;
        };
    }

    it('keeps the mesh consistent under random removals and reinsertions', () => {
        const random = makeRandom(2024);

        // A stack of tetrahedra sharing consecutive faces, plus the closed
        // 4-simplex boundary on a disjoint vertex set.
        const chain: number[][] = [];
        for (let i = 0; i < 6; ++i) {
            chain.push([0, 1, 2 + i, 3 + i]);
        }
        const shifted = fiveCell.map(s => s.map(v => v + 20));
        const tets = chain.concat(shifted);

        for (let trial = 0; trial < 20; ++trial) {
            const mesh = build(tets);
            expect(mesh.getNumTetrahedra()).toBe(tets.length);
            checkInvariants(mesh);

            const present = new Set<number>(tets.map((_, i) => i));
            for (let i = 0; i < tets.length; ++i) {
                if (random() < 0.4) {
                    const s = tets[i];
                    expect(mesh.remove(s[0], s[1], s[2], s[3])).toBe(true);
                    present.delete(i);
                }
            }
            expect(mesh.getNumTetrahedra()).toBe(present.size);
            checkInvariants(mesh);

            // The number of faces equals the number of distinct unordered
            // vertex triples of the remaining tetrahedra.
            const opposite = TetrahedronKey.getOppositeFace();
            const counts = new Map<string, number>();
            for (const i of present) {
                const s = tets[i];
                for (let j = 0; j < 4; ++j) {
                    const o = opposite[j];
                    const key = new TriangleKey(false, s[o[0]], s[o[1]], s[o[2]]).mapKey();
                    counts.set(key, (counts.get(key) ?? 0) + 1);
                }
            }
            expect(mesh.getNumTriangles()).toBe(counts.size);

            // The mesh is closed exactly when every face is shared twice.
            const closed = Array.from(counts.values()).every(c => c === 2);
            expect(mesh.isClosed()).toBe(closed);

            // Reinsert everything that was removed; the mesh is restored.
            for (let i = 0; i < tets.length; ++i) {
                if (!present.has(i)) {
                    const s = tets[i];
                    expect(mesh.insert(s[0], s[1], s[2], s[3])).not.toBeNull();
                }
            }
            expect(mesh.getNumTetrahedra()).toBe(tets.length);
            expect(tetraKeys(mesh)).toEqual(tetraKeys(build(tets)));
            checkInvariants(mesh);
        }
    });
});

describe('TSManifoldMesh verification', () => {
    // Kuhn's subdivision of an n x n x n grid of cubes: each cube is split
    // into the six tetrahedra (c, c+e_p0, c+e_p0+e_p1, c+e_p0+e_p1+e_p2) over
    // the permutations p of (0,1,2). The subdivision is face-conforming
    // across cube boundaries, so the whole grid (and every subset of it) is a
    // triangle-tetrahedron manifold mesh.
    const permutations3: number[][] = [
        [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
    ];

    const kuhnTetrahedra = (n: number): number[][] => {
        const index = (x: number, y: number, z: number): number =>
            x + (n + 1) * (y + (n + 1) * z);
        const tets: number[][] = [];
        for (let k = 0; k < n; ++k) {
            for (let j = 0; j < n; ++j) {
                for (let i = 0; i < n; ++i) {
                    for (const p of permutations3) {
                        const c = [i, j, k];
                        const vs = [index(c[0], c[1], c[2])];
                        for (const axis of p) {
                            c[axis] += 1;
                            vs.push(index(c[0], c[1], c[2]));
                        }
                        tets.push(vs);
                    }
                }
            }
        }
        return tets;
    };

    // A brute-force model: the unordered face -> tetrahedra map.
    const faceMapOf = (tets: number[][]): Map<string, number[][]> => {
        const opposite = TetrahedronKey.getOppositeFace();
        const map = new Map<string, number[][]>();
        for (const s of tets) {
            for (let i = 0; i < 4; ++i) {
                const o = opposite[i];
                const key = new TriangleKey(false, s[o[0]], s[o[1]], s[o[2]])
                    .mapKey();
                const list = map.get(key);
                if (list === undefined) {
                    map.set(key, [s]);
                } else {
                    list.push(s);
                }
            }
        }
        return map;
    };

    const skeyOf = (s: number[]): string =>
        new TetrahedronKey(true, s[0], s[1], s[2], s[3]).mapKey();

    const checkAgainstModel = (mesh: TSManifoldMesh, tets: number[][]): void => {
        expect(mesh.getTetrahedronKeys().map((k) => k.mapKey()).sort())
            .toEqual(tets.map(skeyOf).sort());

        const modelFaces = faceMapOf(tets);
        expect(mesh.getTriangleKeys().map((k) => k.mapKey()).sort())
            .toEqual(Array.from(modelFaces.keys()).sort());

        for (const face of mesh.getTriangles()) {
            const key = new TriangleKey(false, face.V[0], face.V[1], face.V[2])
                .mapKey();
            const shared = modelFaces.get(key) as number[][];
            const stored = face.S.filter((s) => s !== null);
            expect(stored.length).toBe(shared.length);
            expect(stored.map((s) => skeyOf(s.V)).sort())
                .toEqual(shared.map(skeyOf).sort());
        }

        expect(mesh.isClosed()).toBe(Array.from(modelFaces.values())
            .every((list) => list.length === 2));
        checkInvariants(mesh);
    };

    const subsetArb = fc.integer({ min: 1, max: 2 }).chain((n) => {
        const all = kuhnTetrahedra(n);
        return fc.array(fc.boolean(),
            { minLength: all.length, maxLength: all.length })
            .map((mask) => all.filter((_, i) => mask[i]));
    }).filter((tets) => tets.length > 0);

    it('agrees with a brute-force face-tetrahedron model', () => {
        check(subsetArb, (tets) => {
            const mesh = new TSManifoldMesh();
            for (const s of tets) {
                expect(mesh.insert(s[0], s[1], s[2], s[3])).not.toBeNull();
            }
            checkAgainstModel(mesh, tets);
        }, 40);
    });

    it('the full Kuhn subdivision of a cube grid has the expected counts', () => {
        check(fc.integer({ min: 1, max: 3 }), (n) => {
            const tets = kuhnTetrahedra(n);
            const mesh = new TSManifoldMesh();
            for (const s of tets) {
                expect(mesh.insert(s[0], s[1], s[2], s[3])).not.toBeNull();
            }
            expect(mesh.getNumTetrahedra()).toBe(6 * n * n * n);
            // Every interior face is shared twice; the boundary faces are the
            // two triangles of each of the 6 * n * n boundary cube faces.
            const boundary = mesh.getTriangles()
                .filter((f) => f.S[1] === null).length;
            expect(boundary).toBe(12 * n * n);
            expect(mesh.isClosed()).toBe(false);
            checkAgainstModel(mesh, tets);
        }, 20);
    });

    it('random insert/remove sequences agree with the model', () => {
        check(fc.array(fc.tuple(fc.boolean(), fc.nat()),
            { minLength: 1, maxLength: 16 }), (ops) => {
            const all = kuhnTetrahedra(2);
            const mesh = new TSManifoldMesh();
            let model: number[][] = [];
            for (const [insertOp, raw] of ops) {
                const s = all[raw % all.length];
                const key = skeyOf(s);
                if (insertOp) {
                    const present = model.some((u) => skeyOf(u) === key);
                    expect(mesh.insert(s[0], s[1], s[2], s[3]) !== null)
                        .toBe(!present);
                    if (!present) {
                        model.push(s);
                    }
                } else {
                    const present = model.some((u) => skeyOf(u) === key);
                    expect(mesh.remove(s[0], s[1], s[2], s[3])).toBe(present);
                    model = model.filter((u) => skeyOf(u) !== key);
                }
                checkAgainstModel(mesh, model);
            }
        }, 30);
    });

    it('the insertion order does not change the mesh contents', () => {
        check(fc.tuple(subsetArb, fc.array(fc.nat(), { maxLength: 16 })),
            ([tets, shuffle]) => {
                const permuted = tets.slice();
                for (let k = 0; k < shuffle.length; ++k) {
                    const i = shuffle[k] % permuted.length;
                    const j = (shuffle[k] * 5 + k) % permuted.length;
                    const t = permuted[i];
                    permuted[i] = permuted[j];
                    permuted[j] = t;
                }
                const a = build(tets);
                const b = build(permuted);
                expect(b.getTetrahedronKeys().map((k) => k.mapKey()))
                    .toEqual(a.getTetrahedronKeys().map((k) => k.mapKey()));
                expect(b.getTriangleKeys().map((k) => k.mapKey()))
                    .toEqual(a.getTriangleKeys().map((k) => k.mapKey()));
                expect(b.isClosed()).toBe(a.isClosed());
            }, 40);
    });

    it('any even permutation of the vertices names the same tetrahedron', () => {
        check(fc.uniqueArray(fc.integer({ min: 0, max: 12 }),
            { minLength: 4, maxLength: 4 }), (v) => {
            const mesh = new TSManifoldMesh();
            const tetra = mesh.insert(v[0], v[1], v[2], v[3]);
            expect(tetra).not.toBeNull();
            // The 12 even permutations of (0,1,2,3).
            const even = [
                [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2],
                [1, 3, 2, 0], [1, 2, 0, 3], [1, 0, 3, 2],
                [2, 3, 0, 1], [2, 0, 1, 3], [2, 1, 3, 0],
                [3, 1, 0, 2], [3, 0, 2, 1], [3, 2, 1, 0]
            ];
            for (const p of even) {
                expect(mesh.getTetrahedron(v[p[0]], v[p[1]], v[p[2]], v[p[3]]))
                    .toBe(tetra);
                expect(mesh.insert(v[p[0]], v[p[1]], v[p[2]], v[p[3]])).toBeNull();
            }
            // An odd permutation is a different key: the reflected
            // tetrahedron. It shares all four faces with the first one, and
            // each of those faces still has a free slot, so the insertion
            // succeeds and the resulting two-tetrahedron mesh is closed.
            expect(mesh.getTetrahedron(v[0], v[1], v[3], v[2])).toBeNull();
            expect(mesh.insert(v[0], v[1], v[3], v[2])).not.toBeNull();
            expect(mesh.getNumTetrahedra()).toBe(2);
            expect(mesh.getNumTriangles()).toBe(4);
            expect(mesh.isClosed()).toBe(true);
            checkInvariants(mesh);
        });
    });

    it('a third tetrahedron on a face is rejected and leaks a phantom face', () => {
        check(fc.integer({ min: 10, max: 20 }), (w) => {
            const mesh = new TSManifoldMesh();
            mesh.insert(0, 1, 2, 3);
            mesh.insert(0, 1, 2, 4);
            expect(mesh.getNumTetrahedra()).toBe(2);

            // A third tetrahedron on the face {0,1,2}.
            expect(() => mesh.insert(0, 1, 2, w))
                .toThrow('Attempt to create nonmanifold mesh.');
            expect(mesh.getNumTetrahedra()).toBe(2);

            // Upstream quirk, preserved: the faces created by the loop
            // iterations that ran before the failure stay in the face map and
            // reference the rejected tetrahedron. (In C++ that Tetrahedron*
            // dangles once the frame's unique_ptr dies.)
            const leaked = mesh.getTriangles()
                .filter((f) => f.S.some((s) => s !== null &&
                    s.V.includes(w)));
            expect(leaked.length).toBeGreaterThan(0);
            for (const face of leaked) {
                for (const s of face.S) {
                    if (s !== null && s.V.includes(w)) {
                        expect(mesh.getTetrahedra()).not.toContain(s);
                    }
                }
            }

            // The graceful path leaks the same way.
            const graceful = new TSManifoldMesh();
            graceful.throwOnNonmanifoldInsertion(false);
            graceful.insert(0, 1, 2, 3);
            graceful.insert(0, 1, 2, 4);
            expect(graceful.insert(0, 1, 2, w)).toBeNull();
            expect(graceful.getNumTetrahedra()).toBe(2);
            expect(graceful.getNumTriangles()).toBeGreaterThan(6);
        });
    });

    it('clone and assign reproduce the mesh', () => {
        check(subsetArb, (tets) => {
            const mesh = build(tets);
            const copy = mesh.clone();
            expect(copy.getTetrahedronKeys().map((k) => k.mapKey()))
                .toEqual(mesh.getTetrahedronKeys().map((k) => k.mapKey()));
            expect(copy.getTriangleKeys().map((k) => k.mapKey()))
                .toEqual(mesh.getTriangleKeys().map((k) => k.mapKey()));
            checkInvariants(copy);

            const first = mesh.getTetrahedra()[0];
            mesh.remove(first.V[0], first.V[1], first.V[2], first.V[3]);
            expect(copy.getNumTetrahedra()).toBe(tets.length);

            const target = build(fiveCell);
            target.assign(copy);
            expect(target.getTetrahedronKeys().map((k) => k.mapKey()))
                .toEqual(copy.getTetrahedronKeys().map((k) => k.mapKey()));
        }, 40);
    });

    it('removing every tetrahedron empties the mesh', () => {
        check(subsetArb, (tets) => {
            const mesh = build(tets);
            for (const s of tets) {
                expect(mesh.remove(s[0], s[1], s[2], s[3])).toBe(true);
            }
            expect(mesh.getNumTetrahedra()).toBe(0);
            expect(mesh.getNumTriangles()).toBe(0);
            expect(mesh.isClosed()).toBe(true);  // vacuously
        }, 40);
    });
});
