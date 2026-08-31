import { describe, it, expect } from 'vitest';
import { StaticVTSManifoldMesh3 } from '../src/StaticVTSManifoldMesh3';

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
