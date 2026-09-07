import { describe, it, expect } from 'vitest';
import { StaticVTSManifoldMesh3 } from '../src/StaticVTSManifoldMesh3.js';
import {
    TSManifoldMesh, TSManifoldMeshTriangle, TSManifoldMeshTetrahedron
} from '../src/TSManifoldMesh.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { TetrahedronKey } from '../src/TetrahedronKey.js';
import { check, fc } from './helpers/arbitraries.js';

const invalid = StaticVTSManifoldMesh3.invalid;

type P3 = [number, number, number];

// Two tetrahedra that share the triangular face with vertices 0, 1 and 2.
// 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(0,0,1), 4=(0,0,-1).
const pairPositions: P3[] = [
    [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, -1]
];
const pairTetrahedra: [number, number, number, number][] = [
    [0, 1, 2, 3],
    [0, 2, 1, 4]
];

function signedVolume(positions: P3[], tetra: readonly number[]): number {
    const p0 = positions[tetra[0]];
    const e = [1, 2, 3].map(i => [
        positions[tetra[i]][0] - p0[0],
        positions[tetra[i]][1] - p0[1],
        positions[tetra[i]][2] - p0[2]
    ]);
    return e[0][0] * (e[1][1] * e[2][2] - e[1][2] * e[2][1])
        - e[0][1] * (e[1][0] * e[2][2] - e[1][2] * e[2][0])
        + e[0][2] * (e[1][0] * e[2][1] - e[1][1] * e[2][0]);
}

// Kuhn's decomposition of the unit cube into six tetrahedra. The cube corner
// (x,y,z) with x, y, z in {0,1} has index x + 2*y + 4*z. Each tetrahedron is
// reordered, if necessary, so that all six have the canonical chirality.
function makeCube(): {
    numVertices: number;
    positions: P3[];
    tetrahedra: [number, number, number, number][];
} {
    const positions: P3[] = [];
    for (let z = 0; z < 2; ++z) {
        for (let y = 0; y < 2; ++y) {
            for (let x = 0; x < 2; ++x) {
                positions.push([x, y, z]);
            }
        }
    }

    const raw: [number, number, number, number][] = [
        [0, 1, 3, 7],
        [0, 1, 5, 7],
        [0, 2, 3, 7],
        [0, 2, 6, 7],
        [0, 4, 5, 7],
        [0, 4, 6, 7]
    ];

    const tetrahedra = raw.map(tetra => {
        if (signedVolume(positions, tetra) < 0) {
            return [tetra[0], tetra[2], tetra[1], tetra[3]] as
                [number, number, number, number];
        }
        return tetra;
    });

    return { numVertices: positions.length, positions, tetrahedra };
}

// The four faces of a tetrahedron, as unordered sorted triples.
function faceKeys(tetra: readonly number[]): string[] {
    return StaticVTSManifoldMesh3.face.map(f => {
        const vs = [tetra[f[0]], tetra[f[1]], tetra[f[2]]];
        vs.sort((a, b) => a - b);
        return vs.join(',');
    });
}

describe('StaticVTSManifoldMesh3.sortFace', () => {
    it('rotates the face so that the first index is the minimum', () => {
        expect(StaticVTSManifoldMesh3.sortFace(0, 1, 2)).toEqual([0, 1, 2]);
        expect(StaticVTSManifoldMesh3.sortFace(1, 2, 0)).toEqual([0, 1, 2]);
        expect(StaticVTSManifoldMesh3.sortFace(2, 0, 1)).toEqual([0, 1, 2]);
        // The opposite orientation stays opposite.
        expect(StaticVTSManifoldMesh3.sortFace(0, 2, 1)).toEqual([0, 2, 1]);
        expect(StaticVTSManifoldMesh3.sortFace(2, 1, 0)).toEqual([0, 2, 1]);
        expect(StaticVTSManifoldMesh3.sortFace(1, 0, 2)).toEqual([0, 2, 1]);
    });

    it('preserves the cyclic order for every triple of distinct indices', () => {
        const values = [3, 11, 7];
        for (let a = 0; a < 3; ++a) {
            for (let b = 0; b < 3; ++b) {
                for (let c = 0; c < 3; ++c) {
                    if (a === b || a === c || b === c) {
                        continue;
                    }
                    const v: [number, number, number] =
                        [values[a], values[b], values[c]];
                    const u = StaticVTSManifoldMesh3.sortFace(v[0], v[1], v[2]);
                    // The minimum comes first.
                    expect(u[0]).toBe(Math.min(v[0], v[1], v[2]));
                    // The result is a rotation of the input.
                    const rotations = [
                        [v[0], v[1], v[2]],
                        [v[1], v[2], v[0]],
                        [v[2], v[0], v[1]]
                    ];
                    expect(rotations).toContainEqual([u[0], u[1], u[2]]);
                }
            }
        }
    });
});

describe('StaticVTSManifoldMesh3 construction', () => {
    it('rejects invalid input', () => {
        expect(() => new StaticVTSManifoldMesh3(3, pairTetrahedra))
            .toThrow('invalid input');
        expect(() => new StaticVTSManifoldMesh3(5, []))
            .toThrow('invalid input');
    });

    it('uses tetrahedra of consistent chirality in the fixtures', () => {
        for (const tetra of pairTetrahedra) {
            expect(signedVolume(pairPositions, tetra)).toBeGreaterThan(0);
        }
        const cube = makeCube();
        for (const tetra of cube.tetrahedra) {
            expect(signedVolume(cube.positions, tetra)).toBeGreaterThan(0);
        }
    });

    it('stores the tetrahedra and the per-vertex counts', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        expect(mesh.getTetrahedra()).toEqual(pairTetrahedra);
        // Vertices 0, 1 and 2 are shared; vertices 3 and 4 are apexes.
        expect(mesh.getMinNumTetrahedraAtVertex()).toBe(1);
        expect(mesh.getMaxNumTetrahedraAtVertex()).toBe(2);
        expect(mesh.getVertices().map(v => v.getNumSAdjacents()))
            .toEqual([2, 2, 2, 1, 1]);
    });

    it('computes the tetrahedron-tetrahedron adjacency', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        // Face 0 of each tetrahedron is the shared face; the other three
        // faces of each are on the boundary.
        expect(mesh.getAdjacents()).toEqual([
            [1, invalid, invalid, invalid],
            [0, invalid, invalid, invalid]
        ]);
    });

    it('collects the adjacent vertices without duplicates', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        const vertices = mesh.getVertices();
        // Vertex 0 is adjacent to every other vertex.
        expect(vertices[0].getVAdjacents().slice().sort((a, b) => a - b))
            .toEqual([1, 2, 3, 4]);
        expect(vertices[3].getVAdjacents().slice().sort((a, b) => a - b))
            .toEqual([0, 1, 2]);
        expect(vertices[4].getVAdjacents().slice().sort((a, b) => a - b))
            .toEqual([0, 1, 2]);
        for (const vertex of vertices) {
            expect(new Set(vertex.getVAdjacents()).size)
                .toBe(vertex.getNumVAdjacents());
        }
    });

    it('stores the outgoing-face quadruples <AV0,AV1,LS,RS>', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        const vertex0 = mesh.getVertices()[0];
        // Each tetrahedron contributes three faces whose minimum vertex is 0.
        expect(vertex0.getNumFAdjacents()).toBe(6);
        const quads = vertex0.getFAdjacents();
        // The shared face <0,2,1> comes from tetrahedron 0; its opposite
        // orientation <0,1,2> comes from tetrahedron 1.
        expect(quads).toContainEqual([2, 1, 0, 1]);
        expect(quads).toContainEqual([1, 2, 1, 0]);
        // The boundary face <0,1,3> of tetrahedron 0 has no right neighbor.
        expect(quads).toContainEqual([1, 3, 0, invalid]);
    });
});

describe('StaticVTSManifoldMesh3.faceExists', () => {
    it('finds every face of the mesh, in either orientation', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        for (const tetra of pairTetrahedra) {
            for (const f of StaticVTSManifoldMesh3.face) {
                const v0 = tetra[f[0]];
                const v1 = tetra[f[1]];
                const v2 = tetra[f[2]];
                expect(mesh.faceExists(v0, v1, v2)).toBe(true);
                expect(mesh.faceExists(v0, v2, v1)).toBe(true);
                expect(mesh.faceExists(v2, v1, v0)).toBe(true);
            }
        }
    });

    it('rejects nonfaces and degenerate or out-of-range queries', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        // Vertices 3 and 4 are never in a common tetrahedron.
        expect(mesh.faceExists(1, 3, 4)).toBe(false);
        expect(mesh.faceExists(0, 3, 4)).toBe(false);
        expect(mesh.faceExists(0, 1, 1)).toBe(false);
        expect(mesh.faceExists(0, 1, 5)).toBe(false);
        expect(mesh.faceExists(-1, 1, 2)).toBe(false);
    });
});

describe('StaticVTSManifoldMesh3.getAdjacentTetrahedra', () => {
    it('returns both tetrahedra of an interior face', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        // <0,1,2> is outgoing from tetrahedron 1 and <0,2,1> from 0.
        expect(mesh.getAdjacentTetrahedra(0, 1, 2))
            .toEqual({ exists: true, adj0: 1, adj1: 0 });
        // Reversing the orientation swaps the pair.
        expect(mesh.getAdjacentTetrahedra(0, 2, 1))
            .toEqual({ exists: true, adj0: 0, adj1: 1 });
        // The query is invariant under rotation of the face.
        expect(mesh.getAdjacentTetrahedra(1, 2, 0))
            .toEqual({ exists: true, adj0: 1, adj1: 0 });
        expect(mesh.getAdjacentTetrahedra(2, 0, 1))
            .toEqual({ exists: true, adj0: 1, adj1: 0 });
    });

    it('distinguishes the two boundary-face cases (upstream bug fixed)', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        // Case 2: <0,1,3> is an outgoing face of tetrahedron 0.
        expect(mesh.getAdjacentTetrahedra(0, 1, 3))
            .toEqual({ exists: true, adj0: 0, adj1: invalid });
        // Case 3: <0,3,1> is not outgoing but its opposite is. Upstream
        // returns (0, invalid) here, the signature it documents for case 2;
        // the port returns the documented (invalid, valid) pair.
        expect(mesh.getAdjacentTetrahedra(0, 3, 1))
            .toEqual({ exists: true, adj0: invalid, adj1: 0 });

        // The same for a boundary face of the second tetrahedron.
        expect(mesh.getAdjacentTetrahedra(0, 2, 4))
            .toEqual({ exists: true, adj0: 1, adj1: invalid });
        expect(mesh.getAdjacentTetrahedra(0, 4, 2))
            .toEqual({ exists: true, adj0: invalid, adj1: 1 });
    });

    it('returns case 4 for a nonface', () => {
        const mesh = new StaticVTSManifoldMesh3(5, pairTetrahedra);
        expect(mesh.getAdjacentTetrahedra(1, 3, 4))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
        expect(mesh.getAdjacentTetrahedra(0, 1, 1))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
        expect(mesh.getAdjacentTetrahedra(0, 1, 99))
            .toEqual({ exists: false, adj0: invalid, adj1: invalid });
    });
});

describe('StaticVTSManifoldMesh3 on a tetrahedralized cube', () => {
    const { numVertices, tetrahedra } = makeCube();

    it('has mutually consistent tetrahedron adjacency', () => {
        const mesh = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
        const adjacents = mesh.getAdjacents();
        expect(adjacents.length).toBe(6);

        let numBoundaryFaces = 0;
        for (let t = 0; t < tetrahedra.length; ++t) {
            for (let i = 0; i < 4; ++i) {
                const s = adjacents[t][i];
                if (s === invalid) {
                    ++numBoundaryFaces;
                    continue;
                }
                expect(s).not.toBe(t);
                // The neighbor names t back across the shared face.
                expect(adjacents[s].filter(a => a === t).length).toBe(1);
                // The shared face is a face of both tetrahedra.
                const f = StaticVTSManifoldMesh3.face[i];
                const key = [
                    tetrahedra[t][f[0]], tetrahedra[t][f[1]], tetrahedra[t][f[2]]
                ].sort((a, b) => a - b).join(',');
                expect(faceKeys(tetrahedra[s])).toContain(key);
            }
        }
        // The cube surface is 6 squares, each split into 2 triangles.
        expect(numBoundaryFaces).toBe(12);
    });

    it('answers face queries consistently with the adjacency array', () => {
        const mesh = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
        const adjacents = mesh.getAdjacents();

        for (let t = 0; t < tetrahedra.length; ++t) {
            for (let i = 0; i < 4; ++i) {
                const f = StaticVTSManifoldMesh3.face[i];
                const v0 = tetrahedra[t][f[0]];
                const v1 = tetrahedra[t][f[1]];
                const v2 = tetrahedra[t][f[2]];

                expect(mesh.faceExists(v0, v1, v2)).toBe(true);
                const forward = mesh.getAdjacentTetrahedra(v0, v1, v2);
                expect(forward.exists).toBe(true);
                // <v0,v1,v2> is outgoing from t, so adj0 is t and adj1 is the
                // neighbor across that face.
                expect(forward.adj0).toBe(t);
                expect(forward.adj1).toBe(adjacents[t][i]);

                // The opposite orientation swaps the returned pair.
                const backward = mesh.getAdjacentTetrahedra(v0, v2, v1);
                expect(backward.adj0).toBe(forward.adj1);
                expect(backward.adj1).toBe(forward.adj0);
            }
        }
    });

    it('bounds the face adjacency counts by three per incident tetrahedron', () => {
        const mesh = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
        let totalFaces = 0;
        for (const vertex of mesh.getVertices()) {
            expect(vertex.getNumFAdjacents())
                .toBeLessThanOrEqual(3 * vertex.getNumSAdjacents());
            expect(vertex.getNumVAdjacents())
                .toBeLessThanOrEqual(3 * vertex.getNumSAdjacents());
            totalFaces += vertex.getNumFAdjacents();
        }
        // Every tetrahedron contributes exactly four outgoing faces.
        expect(totalFaces).toBe(4 * tetrahedra.length);
    });
});

describe('StaticVTSManifoldMesh3 verification', () => {
    // Kuhn's subdivision of an n x n x n grid of cubes, with each tetrahedron
    // reordered if necessary so that all have positive signed volume (the
    // consistent chirality the class requires).
    const kuhnGrid = (n: number): {
        numVertices: number;
        positions: P3[];
        tetrahedra: [number, number, number, number][];
    } => {
        const positions: P3[] = [];
        const index = (x: number, y: number, z: number): number =>
            x + (n + 1) * (y + (n + 1) * z);
        for (let z = 0; z <= n; ++z) {
            for (let y = 0; y <= n; ++y) {
                for (let x = 0; x <= n; ++x) {
                    positions[index(x, y, z)] = [x, y, z];
                }
            }
        }

        const permutations3 = [
            [0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]
        ];
        const tetrahedra: [number, number, number, number][] = [];
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
                        const tetra = vs as [number, number, number, number];
                        tetrahedra.push(signedVolume(positions, tetra) < 0
                            ? [tetra[0], tetra[2], tetra[1], tetra[3]]
                            : tetra);
                    }
                }
            }
        }
        return { numVertices: positions.length, positions, tetrahedra };
    };

    const meshArb = fc.integer({ min: 1, max: 2 }).chain((n) => {
        const grid = kuhnGrid(n);
        return fc.array(fc.boolean(), {
            minLength: grid.tetrahedra.length,
            maxLength: grid.tetrahedra.length
        }).map((mask) => ({
            numVertices: grid.numVertices,
            positions: grid.positions,
            tetrahedra: grid.tetrahedra.filter((_, i) => mask[i])
        }));
    }).filter(({ tetrahedra }) => tetrahedra.length > 0);

    // The sorted (outgoing) orientation of every face of a tetrahedron and
    // the tetrahedron that owns it.
    const outgoingFaces = (tetrahedra: readonly (readonly number[])[]):
        Map<string, number> => {
        const owner = new Map<string, number>();
        for (let t = 0; t < tetrahedra.length; ++t) {
            for (const f of StaticVTSManifoldMesh3.face) {
                const u = StaticVTSManifoldMesh3.sortFace(
                    tetrahedra[t][f[0]], tetrahedra[t][f[1]],
                    tetrahedra[t][f[2]]);
                owner.set(u.join(','), t);
            }
        }
        return owner;
    };

    it('all generated tetrahedra have the canonical chirality', () => {
        check(fc.integer({ min: 1, max: 2 }), (n) => {
            const { positions, tetrahedra } = kuhnGrid(n);
            for (const tetra of tetrahedra) {
                expect(signedVolume(positions, tetra)).toBeGreaterThan(0);
            }
            // Six tetrahedra per cube, and every face orientation is unique
            // (that is the manifold condition the class assumes).
            expect(tetrahedra.length).toBe(6 * n * n * n);
            expect(outgoingFaces(tetrahedra).size).toBe(4 * tetrahedra.length);
        });
    });

    it('agrees with TSManifoldMesh on faces and adjacency', () => {
        check(meshArb, ({ numVertices, tetrahedra }) => {
            const stat = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
            const dyn = new TSManifoldMesh();
            for (const s of tetrahedra) {
                expect(dyn.insert(s[0], s[1], s[2], s[3])).not.toBeNull();
            }

            // The same unordered faces.
            const dynFaces = new Set(dyn.getTriangleKeys().map((k) => k.mapKey()));
            const statFaces = new Set<string>();
            for (const s of tetrahedra) {
                for (const f of StaticVTSManifoldMesh3.face) {
                    statFaces.add(new TriangleKey(false, s[f[0]], s[f[1]],
                        s[f[2]]).mapKey());
                    expect(stat.faceExists(s[f[0]], s[f[1]], s[f[2]])).toBe(true);
                    // Any orientation of the face is found.
                    expect(stat.faceExists(s[f[2]], s[f[1]], s[f[0]])).toBe(true);
                }
            }
            expect(statFaces).toEqual(dynFaces);

            // The same tetrahedron-tetrahedron adjacency. The two classes
            // index the adjacency differently (StaticVTSManifoldMesh3.face[i]
            // versus TetrahedronKey.getOppositeFace()[i]), so compare through
            // the face vertices.
            const adjacents = stat.getAdjacents();
            const skey = (s: readonly number[]): string =>
                new TetrahedronKey(true, s[0], s[1], s[2], s[3]).mapKey();
            for (let t = 0; t < tetrahedra.length; ++t) {
                const s = tetrahedra[t];
                const dynTetra = dyn.getTetrahedron(s[0], s[1], s[2], s[3]);
                expect(dynTetra).not.toBeNull();
                for (let i = 0; i < 4; ++i) {
                    const f = StaticVTSManifoldMesh3.face[i];
                    const faceKey = new TriangleKey(false, s[f[0]], s[f[1]],
                        s[f[2]]).mapKey();
                    // Find the dynamic mesh's adjacency across that face.
                    const dynFace = dyn.getTriangle(s[f[0]], s[f[1]], s[f[2]]);
                    expect(dynFace).not.toBeNull();
                    const face = dynFace as TSManifoldMeshTriangle;
                    expect(new TriangleKey(false, face.V[0], face.V[1],
                        face.V[2]).mapKey()).toBe(faceKey);
                    const other = face.S[0] === dynTetra ? face.S[1] : face.S[0];
                    if (adjacents[t][i] === invalid) {
                        expect(other).toBeNull();
                    } else {
                        expect(other).not.toBeNull();
                        expect(skey(tetrahedra[adjacents[t][i]]))
                            .toBe(skey((other as TSManifoldMeshTetrahedron).V));
                    }
                }
            }
        }, 30);
    });

    it('the packed vertex storage decodes to the mesh adjacency', () => {
        check(meshArb, ({ numVertices, tetrahedra }) => {
            const stat = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
            const vertices = stat.getVertices();
            expect(vertices.length).toBe(numVertices);

            const counts: number[] = new Array<number>(numVertices).fill(0);
            const neighbors: Set<number>[] = [];
            for (let v = 0; v < numVertices; ++v) {
                neighbors.push(new Set<number>());
            }
            for (const s of tetrahedra) {
                for (let i = 0; i < 4; ++i) {
                    ++counts[s[i]];
                    for (let j = 0; j < 4; ++j) {
                        if (j !== i) {
                            neighbors[s[i]].add(s[j]);
                        }
                    }
                }
            }
            expect(stat.getMinNumTetrahedraAtVertex()).toBe(Math.min(...counts));
            expect(stat.getMaxNumTetrahedraAtVertex()).toBe(Math.max(...counts));

            const owner = outgoingFaces(tetrahedra);
            for (let v = 0; v < numVertices; ++v) {
                const vertex = vertices[v];
                expect(vertex.getNumSAdjacents()).toBe(counts[v]);
                expect(new Set(vertex.getVAdjacents())).toEqual(neighbors[v]);
                expect(vertex.getNumVAdjacents())
                    .toBeLessThanOrEqual(3 * vertex.getNumSAdjacents());
                expect(vertex.getNumFAdjacents())
                    .toBeLessThanOrEqual(3 * vertex.getNumSAdjacents());

                const quads = vertex.getFAdjacents();
                // Every outgoing face whose minimum vertex is v is stored
                // exactly once here.
                const expected = Array.from(owner.keys())
                    .filter((k) => Number(k.split(',')[0]) === v);
                expect(quads.length).toBe(expected.length);
                for (let j = 0; j < quads.length; ++j) {
                    const [av0, av1, ls, rs] = quads[j];
                    expect(v).toBeLessThan(av0);
                    expect(v).toBeLessThan(av1);
                    // LS owns the outgoing face <v,av0,av1>.
                    expect(ls).toBe(owner.get([v, av0, av1].join(',')));
                    // RS owns the opposite orientation, or is 'invalid'.
                    const reverse = owner.get([v, av1, av0].join(','));
                    expect(rs).toBe(reverse === undefined ? invalid : reverse);
                    expect(vertex.getFAdjacent(j)).toEqual(quads[j]);
                }
            }
        }, 30);
    });

    it('getAdjacentTetrahedra reports the documented four cases', () => {
        check(meshArb, ({ numVertices, tetrahedra }) => {
            const stat = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
            const owner = outgoingFaces(tetrahedra);

            // Query every face of the mesh in all six vertex orders, plus a
            // few triples that are not faces.
            for (const s of tetrahedra) {
                for (const f of StaticVTSManifoldMesh3.face) {
                    const v = [s[f[0]], s[f[1]], s[f[2]]];
                    const cyclic = [[0, 1, 2], [1, 2, 0], [2, 0, 1]];
                    const anticyclic = [[0, 2, 1], [2, 1, 0], [1, 0, 2]];
                    const u = StaticVTSManifoldMesh3.sortFace(v[0], v[1], v[2]);
                    const forward = owner.get(u.join(','));
                    const backward = owner.get([u[0], u[2], u[1]].join(','));
                    const expected = {
                        exists: true,
                        adj0: forward === undefined ? invalid : forward,
                        adj1: backward === undefined ? invalid : backward
                    };
                    for (const p of cyclic) {
                        expect(stat.getAdjacentTetrahedra(v[p[0]], v[p[1]],
                            v[p[2]])).toEqual(expected);
                    }
                    // Reversing the face orientation swaps the two results.
                    for (const p of anticyclic) {
                        expect(stat.getAdjacentTetrahedra(v[p[0]], v[p[1]],
                            v[p[2]])).toEqual({
                            exists: true,
                            adj0: expected.adj1,
                            adj1: expected.adj0
                        });
                    }
                    // The queried face always exists, and the L-tetrahedron
                    // of the outgoing orientation is a real tetrahedron.
                    expect(stat.faceExists(v[0], v[1], v[2])).toBe(true);
                    expect(expected.adj0 !== invalid ||
                        expected.adj1 !== invalid).toBe(true);
                }
            }

            // Degenerate and out-of-range queries are case 4.
            const none = { exists: false, adj0: invalid, adj1: invalid };
            expect(stat.getAdjacentTetrahedra(0, 0, 1)).toEqual(none);
            expect(stat.getAdjacentTetrahedra(0, 1, 1)).toEqual(none);
            expect(stat.getAdjacentTetrahedra(-1, 0, 1)).toEqual(none);
            expect(stat.getAdjacentTetrahedra(0, 1, numVertices)).toEqual(none);
            expect(stat.faceExists(-1, 0, 1)).toBe(false);
        }, 30);
    });

    it('a boundary face reports exactly one adjacent tetrahedron', () => {
        check(meshArb, ({ numVertices, tetrahedra }) => {
            const stat = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
            const adjacents = stat.getAdjacents();
            for (let t = 0; t < tetrahedra.length; ++t) {
                for (let i = 0; i < 4; ++i) {
                    const f = StaticVTSManifoldMesh3.face[i];
                    const s = tetrahedra[t];
                    const result = stat.getAdjacentTetrahedra(s[f[0]], s[f[1]],
                        s[f[2]]);
                    // The queried orientation is this tetrahedron's outgoing
                    // face, so adj0 is t itself.
                    expect(result.adj0).toBe(t);
                    expect(result.adj1).toBe(adjacents[t][i]);
                    if (adjacents[t][i] === invalid) {
                        // Case 2: a boundary face.
                        expect(result.adj1).toBe(invalid);
                    }
                }
            }
        }, 30);
    });

    it('the tetrahedron order does not change the mesh relations', () => {
        check(fc.tuple(meshArb, fc.array(fc.nat(), { maxLength: 12 })),
            ([{ numVertices, tetrahedra }, shuffle]) => {
                const permuted = tetrahedra.slice();
                for (let k = 0; k < shuffle.length; ++k) {
                    const i = shuffle[k] % permuted.length;
                    const j = (shuffle[k] * 5 + k) % permuted.length;
                    const t = permuted[i];
                    permuted[i] = permuted[j];
                    permuted[j] = t;
                }
                const a = new StaticVTSManifoldMesh3(numVertices, tetrahedra);
                const b = new StaticVTSManifoldMesh3(numVertices, permuted);
                expect(b.getMinNumTetrahedraAtVertex())
                    .toBe(a.getMinNumTetrahedraAtVertex());
                expect(b.getMaxNumTetrahedraAtVertex())
                    .toBe(a.getMaxNumTetrahedraAtVertex());

                const name = (source: readonly (readonly number[])[]) =>
                    (i: number): string => (i === invalid ? 'none'
                        : new TetrahedronKey(true, source[i][0], source[i][1],
                            source[i][2], source[i][3]).mapKey());
                for (const s of tetrahedra) {
                    for (const f of StaticVTSManifoldMesh3.face) {
                        const v = [s[f[0]], s[f[1]], s[f[2]]];
                        const ra = a.getAdjacentTetrahedra(v[0], v[1], v[2]);
                        const rb = b.getAdjacentTetrahedra(v[0], v[1], v[2]);
                        expect(rb.exists).toBe(ra.exists);
                        expect(name(permuted)(rb.adj0))
                            .toBe(name(tetrahedra)(ra.adj0));
                        expect(name(permuted)(rb.adj1))
                            .toBe(name(tetrahedra)(ra.adj1));
                    }
                }
            }, 20);
    });

    it('the numThreads argument does not change the result', () => {
        check(fc.tuple(meshArb, fc.integer({ min: 0, max: 8 })),
            ([{ numVertices, tetrahedra }, numThreads]) => {
                const a = new StaticVTSManifoldMesh3(numVertices, tetrahedra, 0);
                const b = new StaticVTSManifoldMesh3(numVertices, tetrahedra,
                    numThreads);
                expect(b.getAdjacents()).toEqual(a.getAdjacents());
                expect(b.getVertices().map((v) => v.getFAdjacents()))
                    .toEqual(a.getVertices().map((v) => v.getFAdjacents()));
            }, 20);
    });

    it('sortFace keeps the minimum first and the cyclic order', () => {
        check(fc.uniqueArray(fc.integer({ min: 0, max: 20 }),
            { minLength: 3, maxLength: 3 }), (v) => {
            const u = StaticVTSManifoldMesh3.sortFace(v[0], v[1], v[2]);
            expect(u[0]).toBe(Math.min(v[0], v[1], v[2]));
            const rotations = [
                [v[0], v[1], v[2]], [v[1], v[2], v[0]], [v[2], v[0], v[1]]
            ];
            expect(rotations.map((r) => r.join(','))).toContain(u.join(','));
            // Every cyclic rotation of the input sorts to the same triple.
            for (const r of rotations) {
                expect(StaticVTSManifoldMesh3.sortFace(r[0], r[1], r[2]))
                    .toEqual(u);
            }
            // The reversed face sorts to the reversed triple.
            const reversed = StaticVTSManifoldMesh3.sortFace(v[2], v[1], v[0]);
            expect(reversed).toEqual([u[0], u[2], u[1]]);
        });
    });

    it('the constructor copies the caller tetrahedron array', () => {
        const tetrahedra: [number, number, number, number][] = [[0, 1, 2, 3]];
        const mesh = new StaticVTSManifoldMesh3(4, tetrahedra);
        tetrahedra[0][0] = 99;
        expect(mesh.getTetrahedra()[0]).toEqual([0, 1, 2, 3]);
    });
});
