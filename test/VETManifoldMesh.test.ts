import { describe, it, expect } from 'vitest';
import { VETManifoldMesh, VETManifoldMeshVertex } from '../src/VETManifoldMesh.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V37): independent re-check against VETManifoldMesh.h.
// ---------------------------------------------------------------------------

// The triangles of a triangulated (n x n) grid of unit cells. Any subset of
// these is edge-manifold, because a grid edge is shared by at most two of
// them, so an arbitrary insert/remove sequence over the subset never trips
// the nonmanifold rejection.
function gridTriangles(n: number): [number, number, number][] {
    const index = (r: number, c: number) => r * (n + 1) + c;
    const triangles: [number, number, number][] = [];
    for (let r = 0; r < n; ++r) {
        for (let c = 0; c < n; ++c) {
            triangles.push([index(r, c), index(r, c + 1), index(r + 1, c + 1)]);
            triangles.push([index(r, c), index(r + 1, c + 1), index(r + 1, c)]);
        }
    }
    return triangles;
}

// A brute-force model of the whole mesh, derived only from the set of
// triangles currently in it. The upstream invariants are
//   mVMap keys   = the vertices of the triangles,
//   TAdjacent(v) = the triangles containing v,
//   EAdjacent(v) = the mesh edges containing v,
//   VAdjacent(v) = the other endpoints of those edges,
//   edge.T       = the triangles containing the edge.
interface BruteForceModel {
    vertices: number[];
    vAdjacent: Map<number, number[]>;
    eAdjacent: Map<number, string[]>;
    tAdjacent: Map<number, string[]>;
    edges: string[];
    edgeTriangles: Map<string, string[]>;
}

function bruteForceModel(
    triangles: readonly (readonly number[])[]): BruteForceModel {
    const vAdj = new Map<number, Set<number>>();
    const eAdj = new Map<number, Set<string>>();
    const tAdj = new Map<number, Set<string>>();
    const edgeTris = new Map<string, Set<string>>();
    function touch<K, V>(m: Map<K, Set<V>>, k: K): Set<V> {
        let s = m.get(k);
        if (s === undefined) {
            s = new Set<V>();
            m.set(k, s);
        }
        return s;
    }

    for (const tri of triangles) {
        const t0 = tri[0] as number;
        const t1 = tri[1] as number;
        const t2 = tri[2] as number;
        const tkey = new TriangleKey(true, t0, t1, t2).mapKey();
        const pairs: [number, number][] = [[t0, t1], [t1, t2], [t2, t0]];
        for (const [a, b] of pairs) {
            const ekey = new EdgeKey(false, a, b).mapKey();
            touch(edgeTris, ekey).add(tkey);
            touch(vAdj, a).add(b);
            touch(vAdj, b).add(a);
            touch(eAdj, a).add(ekey);
            touch(eAdj, b).add(ekey);
        }
        touch(tAdj, t0).add(tkey);
        touch(tAdj, t1).add(tkey);
        touch(tAdj, t2).add(tkey);
    }

    const vAdjacent = new Map<number, number[]>();
    for (const [k, s] of vAdj) {
        vAdjacent.set(k, Array.from(s).sort((a, b) => a - b));
    }
    const eAdjacent = new Map<number, string[]>();
    for (const [k, s] of eAdj) {
        eAdjacent.set(k, Array.from(s).sort());
    }
    const tAdjacent = new Map<number, string[]>();
    for (const [k, s] of tAdj) {
        tAdjacent.set(k, Array.from(s).sort());
    }
    const edgeTriangles = new Map<string, string[]>();
    for (const [k, s] of edgeTris) {
        edgeTriangles.set(k, Array.from(s).sort());
    }

    return {
        vertices: Array.from(tAdjacent.keys()).sort((a, b) => a - b),
        vAdjacent,
        eAdjacent,
        tAdjacent,
        edges: Array.from(edgeTriangles.keys()).sort(),
        edgeTriangles
    };
}

function expectMatchesModel(mesh: VETManifoldMesh,
    triangles: readonly (readonly number[])[]): void {
    const model = bruteForceModel(triangles);

    expect(mesh.getNumTriangles()).toBe(triangles.length);
    expect(mesh.getNumVertices()).toBe(model.vertices.length);
    expect(mesh.getNumEdges()).toBe(model.edges.length);
    expect(vIndices(mesh)).toEqual(model.vertices);

    expect(mesh.getEdges().map(e => new EdgeKey(false, e.V[0], e.V[1]).mapKey())
        .sort()).toEqual(model.edges);

    for (const v of model.vertices) {
        const vertex = mesh.getVertex(v);
        expect(vertex).not.toBeNull();
        expect((vertex as VETManifoldMeshVertex).V).toBe(v);
        expect((vertex as VETManifoldMeshVertex).getVAdjacent())
            .toEqual(model.vAdjacent.get(v));
        expect(edgeKeys(mesh, v)).toEqual(model.eAdjacent.get(v));
        expect(triKeys(mesh, v)).toEqual(model.tAdjacent.get(v));
    }

    // The vertex adjacency relation is symmetric.
    for (const v of model.vertices) {
        const vertex = mesh.getVertex(v) as VETManifoldMeshVertex;
        for (const w of vertex.getVAdjacent()) {
            const other = mesh.getVertex(w) as VETManifoldMeshVertex;
            expect(other.getVAdjacent()).toContain(v);
        }
    }

    // The base-class edge records agree with the model.
    for (const edge of mesh.getEdges()) {
        const ekey = new EdgeKey(false, edge.V[0], edge.V[1]).mapKey();
        const tris: string[] = [];
        for (const t of edge.T) {
            if (t !== null) {
                tris.push(new TriangleKey(true, t.V[0], t.V[1], t.V[2]).mapKey());
            }
        }
        tris.sort();
        expect(tris).toEqual(model.edgeTriangles.get(ekey));
    }

    // Every edge referenced by a vertex is still in the mesh.
    const liveEdges = new Set<string>(model.edges);
    for (const v of model.vertices) {
        for (const key of edgeKeys(mesh, v)) {
            expect(liveEdges.has(key)).toBe(true);
        }
    }
}

describe('VETManifoldMesh verification', () => {
    const triangles = gridTriangles(3);

    it('matches a brute-force adjacency model under interleaved insert/remove', () => {
        const op = fc.record({
            doInsert: fc.boolean(),
            t: fc.integer({ min: 0, max: triangles.length - 1 })
        });
        check(fc.array(op, { minLength: 1, maxLength: 30 }), ops => {
            const mesh = new VETManifoldMesh();
            const present = new Set<number>();
            for (const step of ops) {
                const tri = triangles[step.t] as [number, number, number];
                if (step.doInsert) {
                    const result = mesh.insert(tri[0], tri[1], tri[2]);
                    if (present.has(step.t)) {
                        // Upstream returns nullptr for a triangle already in
                        // the mesh and leaves the mesh unchanged.
                        expect(result).toBeNull();
                    }
                    else {
                        expect(result).not.toBeNull();
                        expect(result!.V.slice().sort((a, b) => a - b))
                            .toEqual(tri.slice().sort((a, b) => a - b));
                        present.add(step.t);
                    }
                }
                else {
                    const expected = present.delete(step.t);
                    expect(mesh.remove(tri[0], tri[1], tri[2])).toBe(expected);
                }
                expectMatchesModel(mesh,
                    Array.from(present).map(i => triangles[i] as number[]));
            }

            // Removing everything empties every map.
            for (const t of Array.from(present)) {
                const tri = triangles[t] as [number, number, number];
                expect(mesh.remove(tri[0], tri[1], tri[2])).toBe(true);
            }
            expect(mesh.getNumVertices()).toBe(0);
            expect(mesh.getNumEdges()).toBe(0);
            expect(mesh.getNumTriangles()).toBe(0);
            expect(mesh.getVertices()).toEqual([]);
        }, 100);
    }, 30000);

    it('rejects a nonmanifold insertion without disturbing the adjacency', () => {
        check(fc.integer({ min: 0, max: 6 }), k => {
            // Three triangles sharing the edge <0,1>; the third insertion is
            // nonmanifold.
            const mesh = new VETManifoldMesh();
            expect(mesh.insert(0, 1, 2)).not.toBeNull();
            expect(mesh.insert(1, 0, 3)).not.toBeNull();
            const good: number[][] = [[0, 1, 2], [1, 0, 3]];
            expectMatchesModel(mesh, good);

            // The default behavior throws; the graceful behavior returns
            // null. Either way the vertex layer that VETManifoldMesh adds
            // must be untouched: upstream runs the base-class Insert first
            // and only updates mVMap when it succeeds.
            //
            // Upstream quirk, preserved by the port (documented in
            // ETManifoldMesh.ts): the base-class Insert creates the edges of
            // the rejected triangle one at a time and abandons the ones it
            // already made, so the edge map may grow. Here the loop starts at
            // the edge <4+k,0>, which is new, before it reaches the
            // already-doubly-used edge <0,1> and fails, so exactly one
            // phantom edge is left behind by each rejected insertion.
            const rejected = new TriangleKey(true, 0, 1, 4 + k).mapKey();
            const expectVertexLayerUnchanged = (): void => {
                expect(mesh.getNumTriangles()).toBe(2);
                expect(mesh.getTriangle(0, 1, 4 + k)).toBeNull();
                expect(vIndices(mesh)).toEqual([0, 1, 2, 3]);
                for (const vertex of mesh.getVertices()) {
                    const keys = vertex.getTAdjacent().map(t => new TriangleKey(
                        true, t.V[0], t.V[1], t.V[2]).mapKey());
                    expect(keys).not.toContain(rejected);
                    expect(triKeys(mesh, vertex.V)).toEqual(
                        bruteForceModel(good).tAdjacent.get(vertex.V));
                    expect(vertex.getVAdjacent()).toEqual(
                        bruteForceModel(good).vAdjacent.get(vertex.V));
                }
            };

            expect(() => mesh.insert(0, 1, 4 + k)).toThrow();
            expectVertexLayerUnchanged();
            expect(mesh.getNumEdges()).toBe(6);

            const previous = mesh.throwOnNonmanifoldInsertion(false);
            expect(previous).toBe(true);
            expect(mesh.insert(0, 1, 4 + k)).toBeNull();
            mesh.throwOnNonmanifoldInsertion(previous);
            expectVertexLayerUnchanged();
        }, 7);
    });

    it('clone and assign reproduce the adjacency model', () => {
        const chosenArb = fc.uniqueArray(
            fc.integer({ min: 0, max: triangles.length - 1 }),
            { minLength: 1, maxLength: triangles.length });
        check(chosenArb, chosen => {
            const mesh = new VETManifoldMesh();
            const chosenTriangles = chosen.map(i => triangles[i] as number[]);
            for (const tri of chosenTriangles) {
                expect(mesh.insert(tri[0] as number, tri[1] as number,
                    tri[2] as number)).not.toBeNull();
            }

            const copy = mesh.clone();
            expectMatchesModel(copy, chosenTriangles);

            // The copy shares no vertex object with the original.
            for (const vertex of copy.getVertices()) {
                expect(vertex).not.toBe(mesh.getVertex(vertex.V));
            }

            // Removing from the copy leaves the original alone.
            const first = chosenTriangles[0] as number[];
            expect(copy.remove(first[0] as number, first[1] as number,
                first[2] as number)).toBe(true);
            expectMatchesModel(mesh, chosenTriangles);

            // assign() clears the target first.
            const target = new VETManifoldMesh();
            expect(target.insert(100, 101, 102)).not.toBeNull();
            target.assign(mesh);
            expectMatchesModel(target, chosenTriangles);
        }, 60);
    }, 30000);
});
