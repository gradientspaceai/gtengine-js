import { describe, it, expect } from 'vitest';
import {
    VETNonmanifoldMesh, VETNonmanifoldMeshVertex
} from '../src/VETNonmanifoldMesh.js';
import {
    ETNonmanifoldMeshEdge, ETNonmanifoldMeshTriangle
} from '../src/ETNonmanifoldMesh.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { TriangleKey } from '../src/TriangleKey.js';
import { check, fc } from './helpers/arbitraries.js';

// A vertex with extra client data, used to exercise the VCreator callback.
class TaggedVertex extends VETNonmanifoldMeshVertex {
    tag: string;

    constructor(vIndex: number) {
        super(vIndex);
        this.tag = 'v' + vIndex;
    }
}

function vIndices(mesh: VETNonmanifoldMesh): number[] {
    return mesh.getVertices().map(v => v.V);
}

describe('VETNonmanifoldMesh', () => {
    it('starts empty', () => {
        const mesh = new VETNonmanifoldMesh();
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
        expect(mesh.getVertices()).toEqual([]);
        expect(mesh.getVertex(0)).toBeNull();
    });

    it('creates vertex records when a triangle is inserted', () => {
        const mesh = new VETNonmanifoldMesh();
        const tri = mesh.insert(0, 1, 2);
        expect(tri).not.toBeNull();
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.getNumTriangles()).toBe(1);

        for (const v of [0, 1, 2]) {
            const vertex = mesh.getVertex(v) as VETNonmanifoldMeshVertex;
            expect(vertex).not.toBeNull();
            expect(vertex.V).toBe(v);

            // Each vertex of a lone triangle is adjacent to the other two
            // vertices, to the two incident edges and to the triangle.
            expect(vertex.getVAdjacent()).toEqual([0, 1, 2].filter(w => w !== v));
            expect(vertex.getEAdjacent().length).toBe(2);
            for (const edge of vertex.getEAdjacent()) {
                expect(edge.V[0] === v || edge.V[1] === v).toBe(true);
            }
            expect(vertex.getTAdjacent()).toEqual([tri as ETNonmanifoldMeshTriangle]);
        }
    });

    it('returns null when the triangle is already in the mesh', () => {
        const mesh = new VETNonmanifoldMesh();
        expect(mesh.insert(0, 1, 2)).not.toBeNull();
        expect(mesh.insert(0, 1, 2)).toBeNull();
        // A cyclic rotation is the same ordered triangle key.
        expect(mesh.insert(1, 2, 0)).toBeNull();
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumVertices()).toBe(3);
    });

    it('returns the vertices in increasing index order', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(7, 2, 5);
        mesh.insert(2, 9, 5);
        mesh.insert(0, 2, 7);
        expect(vIndices(mesh)).toEqual([0, 2, 5, 7, 9]);
    });

    it('accumulates adjacency across a fan of triangles', () => {
        // A fan around vertex 0.
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.insert(0, 3, 4);

        const center = mesh.getVertex(0) as VETNonmanifoldMeshVertex;
        expect(center.getVAdjacent()).toEqual([1, 2, 3, 4]);
        expect(center.getEAdjacent().length).toBe(4);
        expect(center.getTAdjacent().length).toBe(3);

        // Vertex 2 is shared by two triangles.
        const v2 = mesh.getVertex(2) as VETNonmanifoldMeshVertex;
        expect(v2.getTAdjacent().length).toBe(2);
        expect(v2.getVAdjacent()).toEqual([0, 1, 3]);

        // Vertex 1 is on the boundary of the fan.
        const v1 = mesh.getVertex(1) as VETNonmanifoldMeshVertex;
        expect(v1.getTAdjacent().length).toBe(1);
        expect(v1.getVAdjacent()).toEqual([0, 2]);
    });

    it('supports a nonmanifold edge shared by three triangles', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 1, 3);
        mesh.insert(0, 1, 4);
        expect(mesh.isManifold()).toBe(false);
        expect(mesh.getNumTriangles()).toBe(3);
        expect(mesh.getNumVertices()).toBe(5);

        const shared = mesh.getEdge(0, 1) as ETNonmanifoldMeshEdge;
        expect(shared.T.size).toBe(3);

        const v0 = mesh.getVertex(0) as VETNonmanifoldMeshVertex;
        expect(v0.getTAdjacent().length).toBe(3);
        expect(v0.getVAdjacent()).toEqual([1, 2, 3, 4]);
        expect(v0.getEAdjacent().length).toBe(4);

        const v2 = mesh.getVertex(2) as VETNonmanifoldMeshVertex;
        expect(v2.getVAdjacent()).toEqual([0, 1]);
    });

    it('removes a triangle and drops the adjacency it contributed', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);

        expect(mesh.remove(0, 2, 3)).toBe(true);
        expect(mesh.getNumTriangles()).toBe(1);
        // Vertex 3 was only in the removed triangle.
        expect(mesh.getVertex(3)).toBeNull();
        expect(vIndices(mesh)).toEqual([0, 1, 2]);

        const v0 = mesh.getVertex(0) as VETNonmanifoldMeshVertex;
        expect(v0.getVAdjacent()).toEqual([1, 2]);
        expect(v0.getEAdjacent().length).toBe(2);
        expect(v0.getTAdjacent().length).toBe(1);

        const v2 = mesh.getVertex(2) as VETNonmanifoldMeshVertex;
        expect(v2.getVAdjacent()).toEqual([0, 1]);
    });

    it('removes the last triangle and empties the vertex map', () => {
        // This is the path guarded by the upstream "Malformed mesh."
        // assertion; the port asserts that the adjacency sets have been
        // emptied rather than the inverted upstream condition.
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        expect(mesh.remove(0, 1, 2)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });

    it('removes triangles from a nonmanifold edge one at a time', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 1, 3);
        mesh.insert(0, 1, 4);

        expect(mesh.remove(0, 1, 4)).toBe(true);
        expect(mesh.getVertex(4)).toBeNull();
        let v0 = mesh.getVertex(0) as VETNonmanifoldMeshVertex;
        expect(v0.getVAdjacent()).toEqual([1, 2, 3]);
        // The shared edge <0,1> survives because two triangles still use it.
        expect((mesh.getEdge(0, 1) as ETNonmanifoldMeshEdge).T.size).toBe(2);
        expect(v0.EAdjacent.has(mesh.getEdge(0, 1) as ETNonmanifoldMeshEdge)).toBe(true);

        expect(mesh.remove(0, 1, 3)).toBe(true);
        expect(mesh.getVertex(3)).toBeNull();
        v0 = mesh.getVertex(0) as VETNonmanifoldMeshVertex;
        expect(v0.getVAdjacent()).toEqual([1, 2]);

        expect(mesh.remove(0, 1, 2)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
    });

    it('returns false when removing a triangle that is not in the mesh', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        expect(mesh.remove(0, 1, 3)).toBe(false);
        expect(mesh.getNumTriangles()).toBe(1);
        expect(mesh.getNumVertices()).toBe(3);
        // A reversed winding is a different ordered triangle key.
        expect(mesh.remove(0, 2, 1)).toBe(false);
    });

    it('clears the vertices along with the edges and triangles', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.clear();
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });

    it('uses the vertex creator callback', () => {
        const mesh = new VETNonmanifoldMesh(v => new TaggedVertex(v));
        mesh.insert(3, 4, 5);
        for (const v of [3, 4, 5]) {
            const vertex = mesh.getVertex(v);
            expect(vertex).toBeInstanceOf(TaggedVertex);
            expect((vertex as TaggedVertex).tag).toBe('v' + v);
        }
    });

    it('makes an independent deep copy with clone()', () => {
        const mesh = new VETNonmanifoldMesh(v => new TaggedVertex(v));
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.insert(0, 1, 4);

        const copy = mesh.clone();
        expect(copy.getNumVertices()).toBe(mesh.getNumVertices());
        expect(copy.getNumEdges()).toBe(mesh.getNumEdges());
        expect(copy.getNumTriangles()).toBe(mesh.getNumTriangles());
        expect(vIndices(copy)).toEqual(vIndices(mesh));

        // The vertex creator is copied, so the copy's vertices are tagged.
        expect(copy.getVertex(0)).toBeInstanceOf(TaggedVertex);

        // No object is shared between the meshes.
        for (const v of vIndices(mesh)) {
            expect(copy.getVertex(v)).not.toBe(mesh.getVertex(v));
        }
        for (const tri of copy.getTriangles()) {
            expect(mesh.getTriangles()).not.toContain(tri);
        }

        // Mutating the original leaves the copy alone.
        mesh.remove(0, 1, 4);
        expect(mesh.getNumTriangles()).toBe(2);
        expect(copy.getNumTriangles()).toBe(3);
        expect(copy.getVertex(4)).not.toBeNull();
        expect(mesh.getVertex(4)).toBeNull();

        // ... and the reverse.
        copy.clear();
        expect(mesh.getNumTriangles()).toBe(2);
        expect(mesh.getNumVertices()).toBe(4);
    });

    it('rebuilds the adjacency of the target when assigning', () => {
        const source = new VETNonmanifoldMesh();
        source.insert(0, 1, 2);
        source.insert(1, 2, 3);

        const target = new VETNonmanifoldMesh();
        target.insert(10, 11, 12);
        target.assign(source);

        expect(vIndices(target)).toEqual([0, 1, 2, 3]);
        expect(target.getNumTriangles()).toBe(2);
        const v1 = target.getVertex(1) as VETNonmanifoldMeshVertex;
        expect(v1.getVAdjacent()).toEqual([0, 2, 3]);
        expect(v1.getTAdjacent().length).toBe(2);
    });

    it('keeps the adjacency consistent through many insert/remove cycles', () => {
        // A strip of triangles plus a nonmanifold flap, built and torn down
        // in a different order, must return to an empty mesh with a
        // consistent adjacency at every step.
        const triangles: Array<[number, number, number]> = [
            [0, 1, 2], [1, 3, 2], [2, 3, 4], [3, 5, 4], [1, 3, 6], [0, 2, 7]
        ];
        const mesh = new VETNonmanifoldMesh();
        for (const [a, b, c] of triangles) {
            expect(mesh.insert(a, b, c)).not.toBeNull();
        }

        const check = (): void => {
            // Every vertex record agrees with the triangles that reference
            // it and with the edges that contain it.
            const expectedT = new Map<number, Set<ETNonmanifoldMeshTriangle>>();
            for (const tri of mesh.getTriangles()) {
                for (const v of tri.V) {
                    if (!expectedT.has(v)) {
                        expectedT.set(v, new Set<ETNonmanifoldMeshTriangle>());
                    }
                    (expectedT.get(v) as Set<ETNonmanifoldMeshTriangle>).add(tri);
                }
            }
            expect(vIndices(mesh)).toEqual(
                Array.from(expectedT.keys()).sort((a, b) => a - b));

            for (const vertex of mesh.getVertices()) {
                const expected = expectedT.get(vertex.V) as Set<ETNonmanifoldMeshTriangle>;
                expect(vertex.getTAdjacent().length).toBe(expected.size);
                for (const tri of vertex.getTAdjacent()) {
                    expect(expected.has(tri)).toBe(true);
                }

                const expectedV = new Set<number>();
                for (const edge of mesh.getEdges()) {
                    if (edge.V[0] === vertex.V) {
                        expectedV.add(edge.V[1]);
                    }
                    else if (edge.V[1] === vertex.V) {
                        expectedV.add(edge.V[0]);
                    }
                }
                expect(vertex.getVAdjacent()).toEqual(
                    Array.from(expectedV).sort((a, b) => a - b));
            }
        };

        check();
        const removalOrder = [2, 0, 5, 1, 4, 3];
        for (const index of removalOrder) {
            const [a, b, c] = triangles[index];
            expect(mesh.remove(a, b, c)).toBe(true);
            check();
        }
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        expect(mesh.getNumTriangles()).toBe(0);
    });

    it('inherits the edge-triangle queries of the base mesh', () => {
        const mesh = new VETNonmanifoldMesh();
        mesh.insert(0, 1, 2);
        mesh.insert(0, 2, 3);
        mesh.insert(10, 11, 12);

        expect(mesh.isManifold()).toBe(true);
        expect(mesh.isClosed()).toBe(false);
        expect(mesh.getComponents().length).toBe(2);
        expect(mesh.getTriangle(0, 1, 2)).not.toBeNull();
        expect(mesh.getTriangle(0, 1, 3)).toBeNull();
    });
});

// ---------------------------------------------------------------------------
// Verification block (V16): model-based testing of the vertex adjacency
// against a brute-force recomputation, over random nonmanifold meshes.
// ---------------------------------------------------------------------------

type Tri = [number, number, number];

function triKeyOf(v0: number, v1: number, v2: number): string {
    return new TriangleKey(true, v0, v1, v2).V.join(',');
}

function edgeKeyOf(v0: number, v1: number): string {
    return new EdgeKey(false, v0, v1).V.join(',');
}

function compareKeyStrings(a: string, b: string): number {
    const x = a.split(',').map(Number);
    const y = b.split(',').map(Number);
    for (let i = 0; i < x.length; ++i) {
        if (x[i] !== y[i]) { return x[i] - y[i]; }
    }
    return 0;
}

// The full adjacency of a set of triangles, computed from scratch.
function bruteForce(triangles: Map<string, Tri>) {
    const vertexTriangles = new Map<number, string[]>();
    const vertexNeighbors = new Map<number, number[]>();
    const vertexEdges = new Map<number, string[]>();
    const edgeTriangles = new Map<string, string[]>();

    const push = <K, V>(m: Map<K, V[]>, k: K, v: V): void => {
        const list = m.get(k);
        if (list === undefined) { m.set(k, [v]); }
        else if (!list.includes(v)) { list.push(v); }
    };

    for (const [key, t] of triangles) {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            const a = t[i0], b = t[i1];
            const e = edgeKeyOf(a, b);
            push(edgeTriangles, e, key);
            push(vertexNeighbors, a, b);
            push(vertexNeighbors, b, a);
            push(vertexEdges, a, e);
            push(vertexEdges, b, e);
        }
        for (const v of t) { push(vertexTriangles, v, key); }
    }
    for (const list of vertexTriangles.values()) { list.sort(compareKeyStrings); }
    for (const list of edgeTriangles.values()) { list.sort(compareKeyStrings); }
    for (const list of vertexNeighbors.values()) { list.sort((a, b) => a - b); }
    for (const list of vertexEdges.values()) { list.sort(compareKeyStrings); }
    return { vertexTriangles, vertexNeighbors, vertexEdges, edgeTriangles };
}

function expectVETMatchesModel(mesh: VETNonmanifoldMesh,
    triangles: Map<string, Tri>): void {
    const model = bruteForce(triangles);
    const expectedVertices = Array.from(model.vertexTriangles.keys())
        .sort((a, b) => a - b);

    expect(mesh.getNumTriangles()).toBe(triangles.size);
    expect(mesh.getNumEdges()).toBe(model.edgeTriangles.size);
    expect(mesh.getNumVertices()).toBe(expectedVertices.length);
    // getVertices() must be in the order of the upstream std::map<int32_t,...>.
    expect(vIndices(mesh)).toEqual(expectedVertices);

    for (const vertex of mesh.getVertices()) {
        const v = vertex.V;
        expect(vertex.getTAdjacent().map(
            t => triKeyOf(t.V[0], t.V[1], t.V[2])))
            .toEqual(model.vertexTriangles.get(v));
        expect(vertex.getVAdjacent()).toEqual(model.vertexNeighbors.get(v));
        expect(vertex.getEAdjacent().map(e => edgeKeyOf(e.V[0], e.V[1])))
            .toEqual(model.vertexEdges.get(v));
        // The adjacent edges and triangles are the mesh's own objects.
        for (const e of vertex.getEAdjacent()) {
            expect(mesh.getEdge(e.V[0], e.V[1])).toBe(e);
        }
        for (const t of vertex.getTAdjacent()) {
            expect(mesh.getTriangle(t.V[0], t.V[1], t.V[2])).toBe(t);
        }
        expect(mesh.getVertex(v)).toBe(vertex);
    }

    // Vertices with no triangle are never kept.
    for (const vertex of mesh.getVertices()) {
        expect(vertex.TAdjacent.size).toBeGreaterThan(0);
    }

    // The inherited edge-triangle adjacency.
    for (const edge of mesh.getEdges()) {
        expect(edge.getTriangles().map(t => triKeyOf(t.V[0], t.V[1], t.V[2])))
            .toEqual(model.edgeTriangles.get(edgeKeyOf(edge.V[0], edge.V[1])));
    }
}

// Distinct-vertex triangles over a small pool. Degenerate triangles are
// excluded: ETNonmanifoldMesh::Remove throws on them (upstream defect).
const triangleArb: fc.Arbitrary<Tri> =
    fc.uniqueArray(fc.integer({ min: 0, max: 5 }),
        { minLength: 3, maxLength: 3 }).map(a => [a[0], a[1], a[2]] as Tri);

const triangleListArb = fc.array(triangleArb, { minLength: 1, maxLength: 12 });

describe('VETNonmanifoldMesh verification', () => {
    it('keeps the vertex adjacency in step with a brute-force recomputation', () => {
        check(fc.tuple(triangleListArb,
            fc.array(fc.boolean(), { minLength: 12, maxLength: 12 })),
        ([triangles, removeFlags]) => {
            const mesh = new VETNonmanifoldMesh();
            const model = new Map<string, Tri>();
            const inserted: Tri[] = [];

            for (let i = 0; i < triangles.length; ++i) {
                const t = triangles[i];
                const key = triKeyOf(t[0], t[1], t[2]);
                const created = mesh.insert(t[0], t[1], t[2]);
                expect(created !== null).toBe(!model.has(key));
                if (created !== null) {
                    model.set(key, t);
                    inserted.push(t);
                }
                expectVETMatchesModel(mesh, model);

                if (removeFlags[i] && inserted.length > 0) {
                    const victim = inserted.splice(i % inserted.length, 1)[0];
                    expect(mesh.remove(victim[0], victim[1], victim[2])).toBe(true);
                    model.delete(triKeyOf(victim[0], victim[1], victim[2]));
                    expectVETMatchesModel(mesh, model);
                    expect(mesh.remove(victim[0], victim[1], victim[2])).toBe(false);
                }
            }

            // Removing every triangle empties the vertex map too. This is also
            // the regression for the upstream 'Malformed mesh' assertion,
            // which is inverted and would fire on the last triangle sharing a
            // vertex (see the port notes in src/VETNonmanifoldMesh.ts).
            for (const t of Array.from(model.values())) {
                expect(mesh.remove(t[0], t[1], t[2])).toBe(true);
                model.delete(triKeyOf(t[0], t[1], t[2]));
                expectVETMatchesModel(mesh, model);
            }
            expect(mesh.getNumVertices()).toBe(0);
            expect(mesh.getNumEdges()).toBe(0);
            expect(mesh.getNumTriangles()).toBe(0);
            expect(mesh.getVertices()).toEqual([]);
        });
    });

    it('removes the only triangle of a mesh without firing the malformed-mesh assertion', () => {
        // Upstream asserts VAdjacent.size() != 0 || EAdjacent.size() != 0 at
        // the moment the vertex is dropped, but by then both sets have been
        // emptied, so the assertion always fires. The port asserts that both
        // are empty instead.
        check(triangleArb, t => {
            const mesh = new VETNonmanifoldMesh();
            expect(mesh.insert(t[0], t[1], t[2])).not.toBeNull();
            expect(() => mesh.remove(t[0], t[1], t[2])).not.toThrow();
            expect(mesh.getNumVertices()).toBe(0);
            expect(mesh.getNumEdges()).toBe(0);
            expect(mesh.getNumTriangles()).toBe(0);
        });
    });

    it('drops a vertex only when its last incident triangle goes away', () => {
        check(triangleListArb, triangles => {
            const mesh = new VETNonmanifoldMesh();
            const model = new Map<string, Tri>();
            for (const t of triangles) {
                if (mesh.insert(t[0], t[1], t[2]) !== null) {
                    model.set(triKeyOf(t[0], t[1], t[2]), t);
                }
            }
            const remaining = Array.from(model.values());
            for (let i = 0; i < remaining.length; ++i) {
                const t = remaining[i];
                mesh.remove(t[0], t[1], t[2]);
                model.delete(triKeyOf(t[0], t[1], t[2]));
                const stillUsed = new Set<number>();
                for (const r of model.values()) {
                    for (const v of r) { stillUsed.add(v); }
                }
                expect(vIndices(mesh)).toEqual(
                    Array.from(stillUsed).sort((a, b) => a - b));
            }
        });
    });

    it('deep-copies the vertices, edges and triangles', () => {
        check(triangleListArb, triangles => {
            const mesh = new VETNonmanifoldMesh();
            const model = new Map<string, Tri>();
            for (const t of triangles) {
                if (mesh.insert(t[0], t[1], t[2]) !== null) {
                    model.set(triKeyOf(t[0], t[1], t[2]), t);
                }
            }

            const copy = mesh.clone();
            expectVETMatchesModel(copy, model);
            for (const vertex of copy.getVertices()) {
                expect(mesh.getVertex(vertex.V)).not.toBe(vertex);
            }
            for (const tri of copy.getTriangles()) {
                expect(mesh.getTriangles()).not.toContain(tri);
            }

            // assign() replaces whatever the target held.
            const target = new VETNonmanifoldMesh();
            expect(target.insert(90, 91, 92)).not.toBeNull();
            target.assign(mesh);
            expectVETMatchesModel(target, model);
            expect(target.getVertex(90)).toBeNull();

            mesh.clear();
            expect(mesh.getNumVertices()).toBe(0);
            expectVETMatchesModel(copy, model);
        });
    });

    it('uses the vertex creator for every distinct vertex index', () => {
        check(triangleListArb, triangles => {
            const created: number[] = [];
            const mesh = new VETNonmanifoldMesh(v => {
                created.push(v);
                return new VETNonmanifoldMeshVertex(v);
            });
            const used = new Set<number>();
            for (const t of triangles) {
                if (mesh.insert(t[0], t[1], t[2]) !== null) {
                    for (const v of t) { used.add(v); }
                }
            }
            expect(created.slice().sort((a, b) => a - b))
                .toEqual(Array.from(used).sort((a, b) => a - b));
        });
    });
});
