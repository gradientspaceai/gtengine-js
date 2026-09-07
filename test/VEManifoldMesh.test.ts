import { describe, it, expect } from 'vitest';
import {
    VEManifoldMesh,
    VEManifoldMeshVertex,
    VEManifoldMeshEdge
} from '../src/VEManifoldMesh.js';
import { check, fc } from './helpers/arbitraries.js';

// Build a closed polyline with the directed edges <0,1>, <1,2>, <2,0>.
function makeTriangleLoop(): VEManifoldMesh {
    const mesh = new VEManifoldMesh();
    expect(mesh.insert(0, 1)).not.toBeNull();
    expect(mesh.insert(1, 2)).not.toBeNull();
    expect(mesh.insert(2, 0)).not.toBeNull();
    return mesh;
}

describe('VEManifoldMesh insertion', () => {
    it('creates the edge and its two vertices', () => {
        const mesh = new VEManifoldMesh();
        const edge = mesh.insert(3, 7);
        expect(edge).not.toBeNull();
        expect(edge!.V).toEqual([3, 7]);
        expect(edge!.E).toEqual([null, null]);

        expect(mesh.getNumVertices()).toBe(2);
        expect(mesh.getNumEdges()).toBe(1);
        expect(mesh.getEdge(3, 7)).toBe(edge);

        // The edge key is the directed pair, so <7,3> is a different edge.
        expect(mesh.getEdge(7, 3)).toBeNull();

        const vertex3 = mesh.getVertex(3);
        const vertex7 = mesh.getVertex(7);
        expect(vertex3).not.toBeNull();
        expect(vertex7).not.toBeNull();
        // A one-edge vertex always has its edge at index zero.
        expect(vertex3!.E).toEqual([edge, null]);
        expect(vertex7!.E).toEqual([edge, null]);
        expect(mesh.isClosed()).toBe(false);
    });

    it('returns null when the directed edge already exists', () => {
        const mesh = new VEManifoldMesh();
        const edge = mesh.insert(0, 1);
        expect(mesh.insert(0, 1)).toBeNull();
        expect(mesh.getNumEdges()).toBe(1);
        expect(mesh.getEdge(0, 1)).toBe(edge);

        // The reversed directed edge is a distinct key and inserts fine, but
        // it makes vertices 0 and 1 shared twice, so the polyline closes.
        expect(mesh.insert(1, 0)).not.toBeNull();
        expect(mesh.getNumEdges()).toBe(2);
        expect(mesh.isClosed()).toBe(true);
    });

    it('links the adjacent edges through the shared vertex', () => {
        const mesh = new VEManifoldMesh();
        const e01 = mesh.insert(0, 1)!;
        const e12 = mesh.insert(1, 2)!;

        // E[i] of an edge points to the edge sharing V[i].
        expect(e01.E).toEqual([null, e12]);
        expect(e12.E).toEqual([e01, null]);

        const vertex1 = mesh.getVertex(1)!;
        expect(vertex1.E).toEqual([e01, e12]);
        expect(mesh.isClosed()).toBe(false);
    });

    it('closes the polyline when every vertex is shared twice', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.isClosed()).toBe(true);

        const e01 = mesh.getEdge(0, 1)!;
        const e12 = mesh.getEdge(1, 2)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(e01.E).toEqual([e20, e12]);
        expect(e12.E).toEqual([e01, e20]);
        expect(e20.E).toEqual([e12, e01]);

        for (const vertex of mesh.getVertices()) {
            expect(vertex.E[0]).not.toBeNull();
            expect(vertex.E[1]).not.toBeNull();
        }
    });

    it('throws when the insertion would make the mesh nonmanifold', () => {
        const mesh = new VEManifoldMesh();
        mesh.insert(0, 1);
        mesh.insert(1, 2);
        // Vertex 1 already has two edges.
        expect(() => mesh.insert(1, 3)).toThrow('The mesh must be manifold.');
    });

    it('returns null instead of throwing when the exception is disabled', () => {
        const mesh = new VEManifoldMesh();
        mesh.insert(0, 1);
        mesh.insert(1, 2);
        mesh.throwOnNonmanifoldInsertion(false);
        expect(mesh.insert(1, 3)).toBeNull();

        // Upstream quirk, preserved by the port: the rejected edge was added
        // to the edge map before the manifold test, so it stays there.
        expect(mesh.getEdge(1, 3)).not.toBeNull();
        // Vertex 3 was never reached, so the loop stopped at the failure.
        expect(mesh.getVertex(3)).toBeNull();
        // Vertex 1 is unchanged: it still refers to the first two edges.
        expect(mesh.getVertex(1)!.E).toEqual([
            mesh.getEdge(0, 1), mesh.getEdge(1, 2)
        ]);
    });
});

describe('VEManifoldMesh removal', () => {
    it('returns false for an edge not in the mesh', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(0, 2)).toBe(false);
        expect(mesh.remove(5, 6)).toBe(false);
        expect(mesh.getNumEdges()).toBe(3);
    });

    it('opens a closed polyline and keeps the endpoints', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(1, 2)).toBe(true);
        expect(mesh.getNumEdges()).toBe(2);
        // Vertices 1 and 2 are still used by the remaining edges.
        expect(mesh.getNumVertices()).toBe(3);
        expect(mesh.isClosed()).toBe(false);

        const e01 = mesh.getEdge(0, 1)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(mesh.getEdge(1, 2)).toBeNull();

        // One-edge vertices always have the edge at index zero.
        expect(mesh.getVertex(1)!.E).toEqual([e01, null]);
        expect(mesh.getVertex(2)!.E).toEqual([e20, null]);
        expect(mesh.getVertex(0)!.E).toEqual([e01, e20]);

        // The adjacency references to the removed edge were cleared.
        expect(e01.E).toEqual([e20, null]);
        expect(e20.E).toEqual([null, e01]);
    });

    it('removes vertices when their last edge is removed', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(0, 1)).toBe(true);
        expect(mesh.remove(1, 2)).toBe(true);
        // Vertex 1 has no remaining edge.
        expect(mesh.getVertex(1)).toBeNull();
        expect(mesh.getNumVertices()).toBe(2);

        expect(mesh.remove(2, 0)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
        // A mesh with no vertices is vacuously closed.
        expect(mesh.isClosed()).toBe(true);
    });

    it('restores the original state after insert then remove', () => {
        const mesh = makeTriangleLoop();
        expect(mesh.remove(1, 2)).toBe(true);
        const e12 = mesh.insert(1, 2)!;
        expect(mesh.isClosed()).toBe(true);
        expect(mesh.getNumEdges()).toBe(3);
        expect(mesh.getNumVertices()).toBe(3);

        const e01 = mesh.getEdge(0, 1)!;
        const e20 = mesh.getEdge(2, 0)!;
        expect(e12.E).toEqual([e01, e20]);
        expect(e01.E).toEqual([e20, e12]);
        expect(e20.E).toEqual([e12, e01]);
    });
});

describe('VEManifoldMesh iteration order', () => {
    it('returns vertices by increasing index and edges lexicographically', () => {
        const mesh = new VEManifoldMesh();
        // Insert out of order; the accessors must sort as std::map iterates.
        mesh.insert(5, 2);
        mesh.insert(2, 9);
        mesh.insert(9, 5);

        expect(mesh.getVertices().map(v => v.V)).toEqual([2, 5, 9]);
        expect(mesh.getEdges().map(e => e.V)).toEqual([[2, 9], [5, 2], [9, 5]]);
    });
});

describe('VEManifoldMesh creators', () => {
    it('uses the supplied vertex and edge creators', () => {
        class TaggedVertex extends VEManifoldMeshVertex {
            tag: string;
            constructor(v: number) {
                super(v);
                this.tag = 'v' + v;
            }
        }
        class TaggedEdge extends VEManifoldMeshEdge {
            tag: string;
            constructor(v0: number, v1: number) {
                super(v0, v1);
                this.tag = 'e' + v0 + '_' + v1;
            }
        }

        const mesh = new VEManifoldMesh(
            v => new TaggedVertex(v),
            (v0, v1) => new TaggedEdge(v0, v1));
        const edge = mesh.insert(4, 6)!;
        expect(edge).toBeInstanceOf(TaggedEdge);
        expect((edge as TaggedEdge).tag).toBe('e4_6');
        expect(mesh.getVertex(4)).toBeInstanceOf(TaggedVertex);
        expect((mesh.getVertex(6) as TaggedVertex).tag).toBe('v6');
    });
});

describe('VEManifoldMesh verification', () => {
    // A brute-force model of the mesh: the set of directed edges that were
    // successfully inserted and are still present.
    class Model {
        edges: [number, number][] = [];

        has(v0: number, v1: number): boolean {
            return this.edges.some((e) => e[0] === v0 && e[1] === v1);
        }

        // The number of model edges incident on a vertex (counting a
        // self-loop <v,v> twice, as the mesh does).
        valence(v: number): number {
            let count = 0;
            for (const e of this.edges) {
                if (e[0] === v) {
                    ++count;
                }
                if (e[1] === v) {
                    ++count;
                }
            }
            return count;
        }

        vertices(): number[] {
            const set = new Set<number>();
            for (const e of this.edges) {
                set.add(e[0]);
                set.add(e[1]);
            }
            return Array.from(set).sort((a, b) => a - b);
        }
    }

    // Check the mesh against the model and against the manifold invariants.
    const checkInvariants = (mesh: VEManifoldMesh, model: Model): void => {
        const edges = mesh.getEdges();
        expect(edges.length).toBe(model.edges.length);
        expect(edges.map((e) => e.V[0] + ',' + e.V[1]).sort())
            .toEqual(model.edges.map((e) => e[0] + ',' + e[1]).sort());
        expect(mesh.getVertices().map((v) => v.V)).toEqual(model.vertices());

        // The number of edge slots incident on a vertex; a self-loop <v,v>
        // occupies two of them.
        const incidence = (v: number): number => edges.reduce((count, e) =>
            count + (e.V[0] === v ? 1 : 0) + (e.V[1] === v ? 1 : 0), 0);

        for (const vertex of mesh.getVertices()) {
            // A one-edge vertex always stores the edge at index 0.
            if (vertex.E[1] !== null) {
                expect(vertex.E[0]).not.toBeNull();
            }
            const stored = vertex.E.filter((e) => e !== null);
            expect(stored.length).toBe(incidence(vertex.V));
            for (const e of stored) {
                expect(edges).toContain(e);
                expect(e.V[0] === vertex.V || e.V[1] === vertex.V).toBe(true);
            }
            // Manifold: a vertex fills at most two edge slots.
            expect(incidence(vertex.V)).toBeLessThanOrEqual(2);
        }

        for (const edge of edges) {
            for (let i = 0; i < 2; ++i) {
                const adjacent = edge.E[i];
                if (adjacent !== null) {
                    // The adjacency is through V[i] and it is symmetric. An
                    // edge is its own neighbour only when it is a self-loop.
                    if (adjacent === edge) {
                        expect(edge.V[0]).toBe(edge.V[1]);
                    }
                    expect(adjacent.V[0] === edge.V[i] ||
                        adjacent.V[1] === edge.V[i]).toBe(true);
                    expect(adjacent.E.includes(edge)).toBe(true);
                } else {
                    // No adjacency means the vertex has only this edge slot.
                    expect(incidence(edge.V[i])).toBe(1);
                }
            }
        }

        // isClosed exactly when every vertex has two incident edges.
        const closed = mesh.getVertices()
            .every((v) => v.E[0] !== null && v.E[1] !== null);
        expect(mesh.isClosed()).toBe(closed);
    };

    it('random insert/remove sequences keep the mesh consistent', () => {
        check(fc.array(fc.tuple(fc.boolean(),
            fc.integer({ min: 0, max: 5 }), fc.integer({ min: 0, max: 5 })),
        { minLength: 1, maxLength: 24 }), (ops) => {
            const mesh = new VEManifoldMesh();
            mesh.throwOnNonmanifoldInsertion(true);
            const model = new Model();
            for (const [insertOp, v0, v1] of ops) {
                if (insertOp) {
                    // Only attempt insertions that keep the mesh manifold,
                    // so the upstream graceful-rejection leak is not
                    // triggered here (it has its own test below).
                    if (model.has(v0, v1)) {
                        expect(mesh.insert(v0, v1)).toBeNull();
                        continue;
                    }
                    const wouldExceed = (v0 === v1)
                        ? model.valence(v0) + 2 > 2
                        : model.valence(v0) + 1 > 2 || model.valence(v1) + 1 > 2;
                    if (wouldExceed) {
                        continue;
                    }
                    const edge = mesh.insert(v0, v1);
                    expect(edge).not.toBeNull();
                    expect((edge as VEManifoldMeshEdge).V).toEqual([v0, v1]);
                    model.edges.push([v0, v1]);
                } else {
                    const present = model.has(v0, v1);
                    expect(mesh.remove(v0, v1)).toBe(present);
                    if (present) {
                        model.edges = model.edges.filter(
                            (e) => !(e[0] === v0 && e[1] === v1));
                    }
                }
                checkInvariants(mesh, model);
            }
        });
    });

    it('a closed polyline is closed and a cut one is not', () => {
        check(fc.integer({ min: 3, max: 12 }), (n) => {
            const mesh = new VEManifoldMesh();
            for (let i = 0; i < n; ++i) {
                expect(mesh.insert(i, (i + 1) % n)).not.toBeNull();
            }
            expect(mesh.getNumEdges()).toBe(n);
            expect(mesh.getNumVertices()).toBe(n);
            expect(mesh.isClosed()).toBe(true);
            // Every edge has both adjacencies set.
            for (const edge of mesh.getEdges()) {
                expect(edge.E[0]).not.toBeNull();
                expect(edge.E[1]).not.toBeNull();
            }

            // Cutting one edge opens the polyline but keeps every vertex.
            expect(mesh.remove(0, 1)).toBe(true);
            expect(mesh.isClosed()).toBe(false);
            expect(mesh.getNumVertices()).toBe(n);
            expect(mesh.getNumEdges()).toBe(n - 1);
            const ends = mesh.getVertices().filter((v) => v.E[1] === null);
            expect(ends.map((v) => v.V).sort((a, b) => a - b)).toEqual([0, 1]);
        });
    });

    it('removing every edge empties the mesh', () => {
        check(fc.integer({ min: 2, max: 10 }), (n) => {
            const mesh = new VEManifoldMesh();
            for (let i = 0; i < n; ++i) {
                mesh.insert(i, (i + 1) % n);
            }
            for (let i = 0; i < n; ++i) {
                expect(mesh.remove(i, (i + 1) % n)).toBe(true);
            }
            expect(mesh.getNumEdges()).toBe(0);
            expect(mesh.getNumVertices()).toBe(0);
            expect(mesh.isClosed()).toBe(true);  // vacuously
        });
    });

    it('the insertion order does not change the final mesh', () => {
        check(fc.tuple(fc.integer({ min: 3, max: 9 }),
            fc.array(fc.nat(), { maxLength: 12 })), ([n, shuffle]) => {
            const edges: [number, number][] = [];
            for (let i = 0; i < n; ++i) {
                edges.push([i, (i + 1) % n]);
            }
            const permuted = edges.slice();
            for (let k = 0; k < shuffle.length; ++k) {
                const i = shuffle[k] % permuted.length;
                const j = (shuffle[k] * 7 + k) % permuted.length;
                const t = permuted[i];
                permuted[i] = permuted[j];
                permuted[j] = t;
            }

            const a = new VEManifoldMesh();
            for (const [v0, v1] of edges) {
                a.insert(v0, v1);
            }
            const b = new VEManifoldMesh();
            for (const [v0, v1] of permuted) {
                b.insert(v0, v1);
            }

            const describeMesh = (mesh: VEManifoldMesh): string[] =>
                mesh.getEdges().map((e) => [
                    e.V.join(','),
                    e.E[0] === null ? '-' : e.E[0].V.join(','),
                    e.E[1] === null ? '-' : e.E[1].V.join(',')
                ].join('|'));
            expect(describeMesh(b)).toEqual(describeMesh(a));
            expect(b.getVertices().map((v) => v.V))
                .toEqual(a.getVertices().map((v) => v.V));
        });
    });

    it('the directed edge key distinguishes <v0,v1> from <v1,v0>', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 4 }),
            fc.integer({ min: 0, max: 4 })).filter(([a, b]) => a !== b),
        ([v0, v1]) => {
            const mesh = new VEManifoldMesh();
            expect(mesh.insert(v0, v1)).not.toBeNull();
            // <v1,v0> is a different edge but would make both vertices
            // have two incident edges, which is still manifold.
            const reversed = mesh.insert(v1, v0);
            expect(reversed).not.toBeNull();
            expect(mesh.getNumEdges()).toBe(2);
            expect(mesh.getNumVertices()).toBe(2);
            expect(mesh.isClosed()).toBe(true);
            expect(mesh.getEdge(v0, v1)).not.toBe(mesh.getEdge(v1, v0));
        });
    });

    it('a third edge at a vertex is rejected (throw and graceful paths)', () => {
        // Throwing path: the mesh must stay usable for the edges inserted
        // before the failure.
        const thrower = new VEManifoldMesh();
        thrower.insert(0, 1);
        thrower.insert(1, 2);
        expect(() => thrower.insert(1, 3)).toThrow('The mesh must be manifold.');

        // Graceful path: upstream returns nullptr but has already inserted
        // the rejected edge into the edge map, so the phantom edge remains.
        // The port preserves this (see the header comment and upstream
        // issue "VEManifoldMesh::Insert leaves rejected edge in map").
        const graceful = new VEManifoldMesh();
        graceful.throwOnNonmanifoldInsertion(false);
        graceful.insert(0, 1);
        graceful.insert(1, 2);
        expect(graceful.insert(1, 3)).toBeNull();
        expect(graceful.getNumEdges()).toBe(3);
        expect(graceful.getEdge(1, 3)).not.toBeNull();
        // Vertex 3 was never created because the loop rejected at i = 0.
        expect(graceful.getVertex(3)).toBeNull();
    });

    it('a self-loop edge consumes both slots of its vertex', () => {
        const mesh = new VEManifoldMesh();
        const edge = mesh.insert(5, 5);
        expect(edge).not.toBeNull();
        const vertex = mesh.getVertex(5) as VEManifoldMeshVertex;
        expect(vertex.E[0]).toBe(edge);
        expect(vertex.E[1]).toBe(edge);
        expect(mesh.isClosed()).toBe(true);
        // The edge is adjacent to itself through V[1].
        expect((edge as VEManifoldMeshEdge).E[1]).toBe(edge);
        expect(mesh.remove(5, 5)).toBe(true);
        expect(mesh.getNumVertices()).toBe(0);
        expect(mesh.getNumEdges()).toBe(0);
    });

    it('the creators are used for every vertex and edge', () => {
        let numVertices = 0;
        let numEdges = 0;
        const mesh = new VEManifoldMesh(
            (v) => {
                ++numVertices;
                return new VEManifoldMeshVertex(v);
            },
            (v0, v1) => {
                ++numEdges;
                return new VEManifoldMeshEdge(v0, v1);
            });
        check(fc.integer({ min: 3, max: 8 }), (n) => {
            numVertices = 0;
            numEdges = 0;
            const fresh = new VEManifoldMesh(
                (v) => {
                    ++numVertices;
                    return new VEManifoldMeshVertex(v);
                },
                (v0, v1) => {
                    ++numEdges;
                    return new VEManifoldMeshEdge(v0, v1);
                });
            for (let i = 0; i < n; ++i) {
                fresh.insert(i, (i + 1) % n);
            }
            expect(numVertices).toBe(n);
            expect(numEdges).toBe(n);
        });
        expect(mesh.getNumEdges()).toBe(0);
    });
});
