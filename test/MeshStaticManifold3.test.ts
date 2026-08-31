import { describe, it, expect } from 'vitest';
import { MeshStaticManifold3, MeshStaticManifold3Vertex } from '../src/MeshStaticManifold3';

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
