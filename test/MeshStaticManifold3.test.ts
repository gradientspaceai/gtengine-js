import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
import { MeshStaticManifold3, MeshStaticManifold3Vertex } from '../src/MeshStaticManifold3.js';

const invalid = MeshStaticManifold3.invalid;

// The canonical tetrahedron <0,1,2,3> with
// V0 = (0,0,0), V1 = (1,0,0), V2 = (0,1,0), V3 = (0,0,1).
const oneTetrahedron = [[0, 1, 2, 3]];

// Two tetrahedra sharing the face {0,1,2}. Tetrahedron 0 = <0,1,2,3> has the
// ordered face <0,2,1> (its face[3]); tetrahedron 1 = <4,0,1,2> has the
// ordered face <0,1,2> (its face[0]). Vertex 4 lies below the plane z = 0, so
// both tetrahedra are positively oriented.
const twoTetrahedra = [[0, 1, 2, 3], [4, 0, 1, 2]];

// The four ordered faces of a tetrahedron <t0,t1,t2,t3>.
function orderedFaces(tet: number[]): number[][] {
    return MeshStaticManifold3.face.map(f => [tet[f[0]], tet[f[1]], tet[f[2]]]);
}

describe('MeshStaticManifold3 face table', () => {
    it('lists the four faces with the opposite vertex omitted', () => {
        expect(MeshStaticManifold3.face.map(f => [...f])).toEqual([
            [1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]
        ]);
        // Face i omits vertex i.
        for (let i = 0; i < 4; ++i) {
            expect(MeshStaticManifold3.face[i]).not.toContain(i);
            expect([...MeshStaticManifold3.face[i]].sort((a, b) => a - b).length).toBe(3);
        }
    });
});

describe('MeshStaticManifold3Vertex', () => {
    it('records inserted adjacency 5-tuples', () => {
        const vertex = new MeshStaticManifold3Vertex();
        expect(vertex.getNumAdjacents()).toBe(0);
        vertex.insert(1, 2, 3, 7, 0);
        expect(vertex.getNumAdjacents()).toBe(1);
        expect(vertex.getAdjacents()).toEqual([[1, 2, 3, 7, 0]]);
    });
});

describe('MeshStaticManifold3 construction', () => {
    it('rejects fewer than four vertices or no tetrahedra', () => {
        expect(() => new MeshStaticManifold3(3, oneTetrahedron)).toThrow('Invalid input.');
        expect(() => new MeshStaticManifold3(4, [])).toThrow('Invalid input.');
    });

    it('copies the tetrahedra rather than aliasing the input', () => {
        const input = [[0, 1, 2, 3]];
        const mesh = new MeshStaticManifold3(4, input);
        input[0][0] = 99;
        expect(mesh.getTetrahedra()[0]).toEqual([0, 1, 2, 3]);
    });

    it('accepts and ignores the upstream numThreads argument', () => {
        const single = new MeshStaticManifold3(5, twoTetrahedra, 1);
        const multi = new MeshStaticManifold3(5, twoTetrahedra, 4);
        expect(multi.getAdjacents()).toEqual(single.getAdjacents());
    });

    it('reports the tetrahedra counts at the vertices', () => {
        const one = new MeshStaticManifold3(4, oneTetrahedron);
        expect(one.getMinNumTetrahedraAtVertex()).toBe(1);
        expect(one.getMaxNumTetrahedraAtVertex()).toBe(1);

        const two = new MeshStaticManifold3(5, twoTetrahedra);
        // Vertices 0, 1, 2 are shared; vertices 3 and 4 are not.
        expect(two.getMinNumTetrahedraAtVertex()).toBe(1);
        expect(two.getMaxNumTetrahedraAtVertex()).toBe(2);
    });

    it('populates one adjacency record per incident tetrahedron', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        const vertices = mesh.getVertices();
        expect(vertices.length).toBe(5);
        expect(vertices[0].getNumAdjacents()).toBe(2);
        expect(vertices[3].getNumAdjacents()).toBe(1);
        expect(vertices[4].getNumAdjacents()).toBe(1);
        // Tetrahedron 0 = <0,1,2,3>: vertex 0 stores {1,2,3,0,a}.
        const record = vertices[0].getAdjacents()[0];
        expect(record.slice(0, 4)).toEqual([1, 2, 3, 0]);
    });
});

describe('MeshStaticManifold3 adjacency', () => {
    it('leaves every face of a lone tetrahedron on the boundary', () => {
        const mesh = new MeshStaticManifold3(4, oneTetrahedron);
        expect(mesh.getAdjacents()).toEqual([[invalid, invalid, invalid, invalid]]);
    });

    it('links the two tetrahedra across their shared face', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        // Tetrahedron 0 = <0,1,2,3>: the shared face {0,1,2} is opposite
        // vertex 3, so it is the location-3 entry.
        // Tetrahedron 1 = <4,0,1,2>: the shared face is opposite vertex 4,
        // the location-0 entry.
        expect(mesh.getAdjacents()).toEqual([
            [invalid, invalid, invalid, 1],
            [0, invalid, invalid, invalid]
        ]);
    });

    it('has symmetric tetrahedron adjacency', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        const adjacents = mesh.getAdjacents();
        for (let t = 0; t < adjacents.length; ++t) {
            for (const s of adjacents[t]) {
                if (s !== invalid) {
                    expect(adjacents[s]).toContain(t);
                }
            }
        }
    });
});

describe('MeshStaticManifold3.faceExists', () => {
    it('finds all four faces of a lone tetrahedron, in either orientation', () => {
        const mesh = new MeshStaticManifold3(4, oneTetrahedron);
        for (const f of orderedFaces([0, 1, 2, 3])) {
            expect(mesh.faceExists(f[0], f[1], f[2])).toBe(true);
            expect(mesh.faceExists(f[0], f[2], f[1])).toBe(true);
            // A cyclic rotation is the same ordered face.
            expect(mesh.faceExists(f[1], f[2], f[0])).toBe(true);
        }
    });

    it('rejects faces that are not in the mesh', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        // Vertices 3 and 4 are on opposite sides and never share a face.
        expect(mesh.faceExists(3, 4, 0)).toBe(false);
        expect(mesh.faceExists(0, 3, 4)).toBe(false);
    });

    it('rejects degenerate and out-of-range queries', () => {
        const mesh = new MeshStaticManifold3(4, oneTetrahedron);
        expect(mesh.faceExists(0, 0, 1)).toBe(false);
        expect(mesh.faceExists(0, 1, 1)).toBe(false);
        expect(mesh.faceExists(0, 1, 4)).toBe(false);
        expect(mesh.faceExists(0, 1, invalid)).toBe(false);
    });

    it('finds the shared face from either side', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        expect(mesh.faceExists(0, 1, 2)).toBe(true);
        expect(mesh.faceExists(0, 2, 1)).toBe(true);
        expect(mesh.faceExists(1, 2, 0)).toBe(true);
    });
});

describe('MeshStaticManifold3.getAdjacentTetrahedra', () => {
    it('returns both tetrahedra of the shared face (documented case 1)', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        // Tetrahedron 1 = <4,0,1,2> owns the ordered face <0,1,2>, so it is
        // the L-tetrahedron for that ordering; tetrahedron 0 owns <0,2,1>.
        expect(mesh.getAdjacentTetrahedra(0, 1, 2)).toEqual(
            { adj0: 1, adj1: 0, exists: true });
        expect(mesh.getAdjacentTetrahedra(0, 2, 1)).toEqual(
            { adj0: 0, adj1: 1, exists: true });
    });

    it('returns (valid, invalid) and (invalid, valid) for a boundary face', () => {
        const mesh = new MeshStaticManifold3(4, oneTetrahedron);
        // <1,2,3> is face[0] of the tetrahedron.
        expect(mesh.getAdjacentTetrahedra(1, 2, 3)).toEqual(
            { adj0: 0, adj1: invalid, exists: true });
        expect(mesh.getAdjacentTetrahedra(1, 3, 2)).toEqual(
            { adj0: invalid, adj1: 0, exists: true });
    });

    it('returns (invalid, invalid) and false for a missing face (case 4)', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        expect(mesh.getAdjacentTetrahedra(0, 3, 4)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
        expect(mesh.getAdjacentTetrahedra(0, 1, 1)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
        expect(mesh.getAdjacentTetrahedra(0, 1, 9)).toEqual(
            { adj0: invalid, adj1: invalid, exists: false });
    });

    it('returns tetrahedron indices, never vertex indices', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        const numTetrahedra = mesh.getTetrahedra().length;
        for (const f of [[0, 1, 2], [0, 2, 1], [1, 2, 3], [0, 1, 3], [0, 1, 4]]) {
            const result = mesh.getAdjacentTetrahedra(f[0], f[1], f[2]);
            for (const t of [result.adj0, result.adj1]) {
                if (t !== invalid) {
                    expect(t).toBeGreaterThanOrEqual(0);
                    expect(t).toBeLessThan(numTetrahedra);
                    // The reported tetrahedron really does contain the face.
                    for (const v of f) {
                        expect(mesh.getTetrahedra()[t]).toContain(v);
                    }
                }
            }
        }
    });

    it('agrees with faceExists and swaps under a face reversal', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        const n = mesh.getVertices().length;
        for (let a = 0; a < n; ++a) {
            for (let b = 0; b < n; ++b) {
                for (let c = 0; c < n; ++c) {
                    const result = mesh.getAdjacentTetrahedra(a, b, c);
                    expect(result.exists).toBe(mesh.faceExists(a, b, c));
                    const reversed = mesh.getAdjacentTetrahedra(a, c, b);
                    expect(reversed.adj0).toBe(result.adj1);
                    expect(reversed.adj1).toBe(result.adj0);
                }
            }
        }
    });

    it('is invariant under a cyclic rotation of the face', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        for (const tet of mesh.getTetrahedra()) {
            for (const f of orderedFaces(tet)) {
                const direct = mesh.getAdjacentTetrahedra(f[0], f[1], f[2]);
                const rotated = mesh.getAdjacentTetrahedra(f[1], f[2], f[0]);
                expect(rotated).toEqual(direct);
            }
        }
    });

    it('reproduces the tetrahedron adjacency array', () => {
        const mesh = new MeshStaticManifold3(5, twoTetrahedra);
        const tetrahedra = mesh.getTetrahedra();
        const adjacents = mesh.getAdjacents();
        for (let t = 0; t < tetrahedra.length; ++t) {
            for (let i = 0; i < 4; ++i) {
                // Face i of tetrahedron t is opposite vertex i, so the
                // neighbor across it is adjacents[t][i].
                const f = MeshStaticManifold3.face[i];
                const tet = tetrahedra[t];
                const result = mesh.getAdjacentTetrahedra(
                    tet[f[0]], tet[f[1]], tet[f[2]]);
                expect(result.exists).toBe(true);
                // The queried ordering belongs to t, so t is the L-tetra.
                expect(result.adj0).toBe(t);
                expect(result.adj1).toBe(adjacents[t][i]);
            }
        }
    });
});

describe('MeshStaticManifold3 on a subdivided cube', () => {
    // The unit cube split into six tetrahedra sharing the diagonal 0-7.
    // Vertex v has coordinates (v & 1, (v >> 1) & 1, (v >> 2) & 1).
    const cube = [
        [0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
        [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]
    ];

    it('builds a mesh in which every interior face has two tetrahedra', () => {
        const mesh = new MeshStaticManifold3(8, cube);
        expect(mesh.getTetrahedra().length).toBe(6);
        // Vertices 0 and 7 belong to all six tetrahedra.
        expect(mesh.getMaxNumTetrahedraAtVertex()).toBe(6);
        expect(mesh.getMinNumTetrahedraAtVertex()).toBe(2);
    });

    it('has symmetric adjacency and consistent queries', () => {
        const mesh = new MeshStaticManifold3(8, cube);
        const adjacents = mesh.getAdjacents();
        let numInterior = 0;
        for (let t = 0; t < adjacents.length; ++t) {
            for (let i = 0; i < 4; ++i) {
                const s = adjacents[t][i];
                if (s !== invalid) {
                    ++numInterior;
                    expect(adjacents[s]).toContain(t);
                    // The shared face is common to both tetrahedra.
                    const f = MeshStaticManifold3.face[i];
                    const tet = mesh.getTetrahedra()[t];
                    const shared = [tet[f[0]], tet[f[1]], tet[f[2]]];
                    for (const v of shared) {
                        expect(mesh.getTetrahedra()[s]).toContain(v);
                    }
                    const result = mesh.getAdjacentTetrahedra(
                        shared[0], shared[1], shared[2]);
                    expect(result.adj0).toBe(t);
                    expect(result.adj1).toBe(s);
                }
            }
        }
        // Six tetrahedra in a fan around the diagonal share six interior
        // faces, counted once from each side.
        expect(numInterior).toBe(12);
    });

    it('reports the cube surface faces as boundary faces', () => {
        const mesh = new MeshStaticManifold3(8, cube);
        // The face <0,1,3> is on the cube surface (the plane z = 0).
        const result = mesh.getAdjacentTetrahedra(0, 1, 3);
        expect(result.exists).toBe(true);
        expect(result.adj0 === invalid || result.adj1 === invalid).toBe(true);
    });
});

describe('MeshStaticManifold3 verification', () => {
    type Mesh = { numVertices: number; tetrahedra: number[][] };

    // The unit cube split into six positively oriented tetrahedra that share
    // the diagonal 0-7.
    const cubeMesh: Mesh = {
        numVertices: 8,
        tetrahedra: [
            [0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
            [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]
        ]
    };
    const oneTet: Mesh = { numVertices: 4, tetrahedra: [[0, 1, 2, 3]] };
    const twoTets: Mesh = { numVertices: 5, tetrahedra: [[0, 1, 2, 3], [4, 0, 1, 2]] };
    // Two cubes side by side, sharing no vertices: several components.
    const twoCubes: Mesh = {
        numVertices: 16,
        tetrahedra: [...cubeMesh.tetrahedra,
        ...cubeMesh.tetrahedra.map(t => t.map(v => v + 8))]
    };

    // The even permutations of (0,1,2,3): relabelling a tetrahedron's own
    // vertices by one of these preserves its orientation.
    const evenPerms: number[][] = [
        [0, 1, 2, 3], [0, 2, 3, 1], [0, 3, 1, 2],
        [1, 0, 3, 2], [1, 2, 0, 3], [1, 3, 2, 0],
        [2, 0, 1, 3], [2, 1, 3, 0], [2, 3, 0, 1],
        [3, 0, 2, 1], [3, 1, 0, 2], [3, 2, 1, 0]
    ];

    function shuffle(mesh: Mesh, perm: number[], order: number[], rots: number[]): Mesh {
        const tetrahedra = order.map((t, k) => {
            const tet = mesh.tetrahedra[t]!;
            const e = evenPerms[rots[k]! % evenPerms.length]!;
            return [perm[tet[e[0]!]!]!, perm[tet[e[1]!]!]!,
            perm[tet[e[2]!]!]!, perm[tet[e[3]!]!]!];
        });
        return { numVertices: mesh.numVertices, tetrahedra };
    }

    const meshArb: fc.Arbitrary<Mesh> =
        fc.constantFrom(oneTet, twoTets, cubeMesh, twoCubes).chain(mesh =>
            fc.tuple(
                fc.shuffledSubarray([...Array(mesh.numVertices).keys()],
                    { minLength: mesh.numVertices, maxLength: mesh.numVertices }),
                fc.shuffledSubarray([...Array(mesh.tetrahedra.length).keys()],
                    { minLength: mesh.tetrahedra.length, maxLength: mesh.tetrahedra.length }),
                fc.array(fc.integer({ min: 0, max: 11 }),
                    { minLength: mesh.tetrahedra.length, maxLength: mesh.tetrahedra.length }))
                .map(([perm, order, rots]) => shuffle(mesh, perm, order, rots)));

    // ---- independent references -------------------------------------------
    /** The canonical key of an oriented triangle, invariant under rotation. */
    function faceKey(a: number, b: number, c: number): string {
        if (a <= b && a <= c) { return `${a},${b},${c}`; }
        if (b <= a && b <= c) { return `${b},${c},${a}`; }
        return `${c},${a},${b}`;
    }

    /** Map from each outward-oriented face to the tetrahedron that has it. */
    function orientedFaceMap(tetrahedra: number[][]): Map<string, number> {
        const map = new Map<string, number>();
        tetrahedra.forEach((tet, t) => {
            for (const f of MeshStaticManifold3.face) {
                map.set(faceKey(tet[f[0]!]!, tet[f[1]!]!, tet[f[2]!]!), t);
            }
        });
        return map;
    }

    /** Expected adjacency: across the face opposite tet[i]. */
    function expectedAdjacents(tetrahedra: number[][]): number[][] {
        const map = orientedFaceMap(tetrahedra);
        return tetrahedra.map(tet => MeshStaticManifold3.face.map(f => {
            const a = tet[f[0]!]!, b = tet[f[1]!]!, c = tet[f[2]!]!;
            // The neighbour has the same face with the opposite orientation.
            const s = map.get(faceKey(a, c, b));
            return s === undefined ? invalid : s;
        }));
    }

    // The 5-tuple prefixes inserted at position p of a tetrahedron, from
    // upstream PopulateVertices.
    const insertOrder: number[][] = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]];

    // ---- properties -------------------------------------------------------
    it('the tetrahedron adjacency matches a brute-force oriented-face map', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            expect(m.getAdjacents()).toEqual(expectedAdjacents(mesh.tetrahedra));
            return true;
        });
    });

    it('the adjacency relation is symmetric and refers to the shared face', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            const adj = m.getAdjacents();
            const tets = m.getTetrahedra();
            for (let t = 0; t < adj.length; ++t) {
                for (let i = 0; i < 4; ++i) {
                    const s = adj[t]![i]!;
                    if (s === invalid) { continue; }
                    expect(s).not.toBe(t);
                    const f = MeshStaticManifold3.face[i]!;
                    const shared = new Set([tets[t]![f[0]!]!, tets[t]![f[1]!]!,
                    tets[t]![f[2]!]!]);
                    let matches = 0;
                    for (let j = 0; j < 4; ++j) {
                        if (adj[s]![j] === t) {
                            const g = MeshStaticManifold3.face[j]!;
                            expect(new Set([tets[s]![g[0]!]!, tets[s]![g[1]!]!,
                            tets[s]![g[2]!]!])).toEqual(shared);
                            ++matches;
                        }
                    }
                    expect(matches).toBe(1);
                }
            }
            return true;
        });
    });

    it('the vertex 5-tuples list the incident tetrahedra in index order', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            const adj = expectedAdjacents(mesh.tetrahedra);
            for (let v = 0; v < mesh.numVertices; ++v) {
                const expected: number[][] = [];
                mesh.tetrahedra.forEach((tet, t) => {
                    const p = tet.indexOf(v);
                    if (p >= 0) {
                        const o = insertOrder[p]!;
                        expected.push([tet[o[0]!]!, tet[o[1]!]!, tet[o[2]!]!, t, adj[t]![p]!]);
                    }
                });
                expect(m.getVertices()[v]!.getAdjacents()).toEqual(expected);
            }
            return true;
        });
    });

    it('the min and max tetrahedra counts are the extremes over all vertices', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            const counts = new Array<number>(mesh.numVertices).fill(0);
            for (const tet of mesh.tetrahedra) { for (const v of tet) { ++counts[v]!; } }
            return m.getMinNumTetrahedraAtVertex() === Math.min(...counts)
                && m.getMaxNumTetrahedraAtVertex() === Math.max(...counts);
        });
    });

    it('faceExists is exactly membership in the unordered face set', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            const unordered = new Set<string>();
            for (const tet of mesh.tetrahedra) {
                for (const f of MeshStaticManifold3.face) {
                    unordered.add([tet[f[0]!]!, tet[f[1]!]!, tet[f[2]!]!]
                        .sort((a, b) => a - b).join(','));
                }
            }
            const n = mesh.numVertices;
            for (let a = 0; a < n; ++a) {
                for (let b = 0; b < n; ++b) {
                    for (let c = 0; c < n; ++c) {
                        const distinct = a !== b && a !== c && b !== c;
                        const expected = distinct
                            && unordered.has([a, b, c].sort((x, y) => x - y).join(','));
                        expect(m.faceExists(a, b, c)).toBe(expected);
                    }
                }
            }
            return true;
        }, 25);
    });

    it('getAdjacentTetrahedra returns the L- and R-tetrahedra of the face', () => {
        // This pins the port's fix of the upstream bug: adj0 is the
        // tetrahedron with the outward face <v0,v1,v2> and adj1 the one with
        // <v0,v2,v1>. Upstream returns a vertex index as adj0.
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            const map = orientedFaceMap(mesh.tetrahedra);
            const n = mesh.numVertices;
            for (let a = 0; a < n; ++a) {
                for (let b = 0; b < n; ++b) {
                    for (let c = 0; c < n; ++c) {
                        const r = m.getAdjacentTetrahedra(a, b, c);
                        if (a === b || a === c || b === c) {
                            expect(r).toEqual({ adj0: invalid, adj1: invalid, exists: false });
                            continue;
                        }
                        const l = map.get(faceKey(a, b, c));
                        const rr = map.get(faceKey(a, c, b));
                        expect(r.adj0).toBe(l === undefined ? invalid : l);
                        expect(r.adj1).toBe(rr === undefined ? invalid : rr);
                        expect(r.exists).toBe(l !== undefined || rr !== undefined);
                        expect(r.exists).toBe(m.faceExists(a, b, c));
                        // Reversing the face swaps the results.
                        const s = m.getAdjacentTetrahedra(a, c, b);
                        expect(s.adj0).toBe(r.adj1);
                        expect(s.adj1).toBe(r.adj0);
                        // A cyclic rotation names the same oriented face.
                        const rot = m.getAdjacentTetrahedra(b, c, a);
                        expect(rot).toEqual(r);
                    }
                }
            }
            return true;
        }, 25);
    });

    it('every returned index is a tetrahedron index, never a vertex index', () => {
        check(meshArb, mesh => {
            const m = new MeshStaticManifold3(mesh.numVertices, mesh.tetrahedra);
            for (const tet of mesh.tetrahedra) {
                for (const f of MeshStaticManifold3.face) {
                    const r = m.getAdjacentTetrahedra(tet[f[0]!]!, tet[f[1]!]!, tet[f[2]!]!);
                    expect(r.exists).toBe(true);
                    for (const a of [r.adj0, r.adj1]) {
                        if (a !== invalid) {
                            expect(Number.isInteger(a)).toBe(true);
                            expect(a).toBeLessThan(mesh.tetrahedra.length);
                        }
                    }
                    // adj0 is the tetrahedron that owns this outward face.
                    expect(r.adj0).not.toBe(invalid);
                    const owner = mesh.tetrahedra[r.adj0]!;
                    for (const v of [tet[f[0]!]!, tet[f[1]!]!, tet[f[2]!]!]) {
                        expect(owner).toContain(v);
                    }
                }
            }
            return true;
        });
    });

    it('all four documented cases of getAdjacentTetrahedra occur', () => {
        const m = new MeshStaticManifold3(5, twoTets.tetrahedra);
        // Case 1: the shared face {0,1,2}. Tetrahedron 1 = <4,0,1,2> has the
        // outward face <0,1,2>; tetrahedron 0 = <0,1,2,3> has <0,2,1>.
        expect(m.getAdjacentTetrahedra(0, 1, 2)).toEqual({ adj0: 1, adj1: 0, exists: true });
        // Case 2: an outward boundary face with no partner.
        const c2 = m.getAdjacentTetrahedra(1, 2, 3);
        expect(c2).toEqual({ adj0: 0, adj1: invalid, exists: true });
        // Case 3: the same face reversed, unreachable in the upstream code.
        const c3 = m.getAdjacentTetrahedra(1, 3, 2);
        expect(c3).toEqual({ adj0: invalid, adj1: 0, exists: true });
        // Case 4: a face that is not in the mesh.
        expect(m.getAdjacentTetrahedra(0, 3, 4))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
    });

    it('rejects negative vertex indices instead of throwing (regression)', () => {
        // Upstream takes size_t parameters, so a negative index becomes a
        // huge unsigned value and the 'v < numVertices' guard rejects it.
        // Before the fix the port indexed the vertex array with the negative
        // value and threw a TypeError.
        const m = new MeshStaticManifold3(5, twoTets.tetrahedra);
        expect(m.faceExists(-1, 0, 1)).toBe(false);
        expect(m.faceExists(0, -1, 1)).toBe(false);
        expect(m.faceExists(0, 1, -1)).toBe(false);
        expect(m.getAdjacentTetrahedra(-1, 0, 1))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
        expect(m.getAdjacentTetrahedra(0, 1, -3))
            .toEqual({ adj0: invalid, adj1: invalid, exists: false });
    });

    it('the constructor copies the tetrahedra and rejects bad sizes', () => {
        check(meshArb, mesh => {
            const input = mesh.tetrahedra.map(t => [...t]);
            const m = new MeshStaticManifold3(mesh.numVertices, input);
            input[0]![0] = 4242;
            expect(m.getTetrahedra()[0]).toEqual(mesh.tetrahedra[0]);
            return true;
        });
        check(fc.integer({ min: -3, max: 3 }), n => {
            expect(() => new MeshStaticManifold3(n, [[0, 1, 2, 3]])).toThrow();
        });
        expect(() => new MeshStaticManifold3(4, [])).toThrow();
    });

    it('the face table gives the four outward faces of the canonical tetrahedron', () => {
        // face[i] omits vertex i and, with V0=(0,0,0), V1=(1,0,0), V2=(0,1,0),
        // V3=(0,0,1), each listed triple is counterclockwise seen from outside;
        // equivalently its normal points away from the opposite vertex.
        const p = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
        for (let i = 0; i < 4; ++i) {
            const f = MeshStaticManifold3.face[i]!;
            expect(f).not.toContain(i);
            const a = p[f[0]!]!, b = p[f[1]!]!, c = p[f[2]!]!;
            const e0 = [b[0]! - a[0]!, b[1]! - a[1]!, b[2]! - a[2]!];
            const e1 = [c[0]! - a[0]!, c[1]! - a[1]!, c[2]! - a[2]!];
            const nrm = [
                e0[1]! * e1[2]! - e0[2]! * e1[1]!,
                e0[2]! * e1[0]! - e0[0]! * e1[2]!,
                e0[0]! * e1[1]! - e0[1]! * e1[0]!
            ];
            const opp = p[i]!;
            const d = nrm[0]! * (opp[0]! - a[0]!) + nrm[1]! * (opp[1]! - a[1]!)
                + nrm[2]! * (opp[2]! - a[2]!);
            expect(d).toBeLessThan(0);
        }
    });
});
