import { describe, it, expect } from 'vitest';
import { MinimumSpanningTree, type MSTEdge } from '../src/MinimumSpanningTree.js';
import { check, fc } from './helpers/arbitraries.js';

const nil = MinimumSpanningTree.nil;

// Deterministic pseudorandom generator so failures are reproducible.
function makeRng(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 0x100000000;
    };
}

// An independent minimum-spanning-forest weight, computed by Kruskal's
// algorithm with a union-find structure. Prim's algorithm and Kruskal's
// algorithm can produce different trees when weights tie, but the total
// weight of a minimum spanning forest is unique.
function kruskalWeight(edges: MSTEdge[], weights: number[]): number {
    const parent = new Map<number, number>();
    const find = (v: number): number => {
        let root = v;
        while (parent.get(root) !== root) {
            root = parent.get(root) as number;
        }
        while (parent.get(v) !== root) {
            const next = parent.get(v) as number;
            parent.set(v, root);
            v = next;
        }
        return root;
    };

    for (const edge of edges) {
        for (const v of edge) {
            if (!parent.has(v)) {
                parent.set(v, v);
            }
        }
    }

    const order = edges.map((_, i) => i);
    order.sort((i, j) => weights[i] - weights[j]);

    let total = 0;
    for (const i of order) {
        const r0 = find(edges[i][0]);
        const r1 = find(edges[i][1]);
        if (r0 !== r1) {
            parent.set(r0, r1);
            total += weights[i];
        }
    }
    return total;
}

// Sum the weights of the non-sentinel tree edges.
function treeWeight(tree: MSTEdge[], edges: MSTEdge[], weights: number[]): number {
    const weightOf = new Map<string, number>();
    for (let e = 0; e < edges.length; ++e) {
        const [v0, v1] = edges[e];
        weightOf.set(v0 + ',' + v1, weights[e]);
        weightOf.set(v1 + ',' + v0, weights[e]);
    }

    let total = 0;
    for (const treeEdge of tree) {
        if (treeEdge[0] !== nil) {
            const w = weightOf.get(treeEdge[0] + ',' + treeEdge[1]);
            expect(w).not.toBeUndefined();
            total += w as number;
        }
    }
    return total;
}

// The number of distinct vertices referenced by the edges.
function vertexCount(edges: MSTEdge[]): number {
    const set = new Set<number>();
    for (const edge of edges) {
        set.add(edge[0]);
        set.add(edge[1]);
    }
    return set.size;
}

describe('MinimumSpanningTree on a known graph', () => {
    // A square with one diagonal:
    //   (0,1) w=1, (1,2) w=2, (2,3) w=3, (3,0) w=4, (0,2) w=5
    const edges: MSTEdge[] = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 2]];
    const weights = [1, 2, 3, 4, 5];

    it('produces the unique minimum spanning tree', () => {
        const { minimumSpanningTree, backEdges } =
            MinimumSpanningTree.execute(edges, weights, true);

        // One sentinel plus one edge per remaining vertex.
        expect(minimumSpanningTree.length).toBe(4);
        expect(minimumSpanningTree[0]).toEqual([nil, 0]);

        // The weights 1, 2, 3 are the unique minimum choice; the edges (3,0)
        // and (0,2) close cycles and are more expensive.
        const treeEdges = minimumSpanningTree.slice(1)
            .map(e => [Math.min(e[0], e[1]), Math.max(e[0], e[1])] as MSTEdge);
        treeEdges.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]));
        expect(treeEdges).toEqual([[0, 1], [1, 2], [2, 3]]);
        expect(treeWeight(minimumSpanningTree, edges, weights)).toBe(6);

        // The back edges are the graph edges not in the tree, stored with
        // v0 < v1 and in lexicographic order.
        expect(backEdges).toEqual([[0, 2], [0, 3]]);
    });

    it('gives the same result when validation is skipped', () => {
        const validated = MinimumSpanningTree.execute(edges, weights, true);
        const unvalidated = MinimumSpanningTree.execute(edges, weights, false);
        expect(unvalidated).toEqual(validated);
    });

    it('handles a single edge', () => {
        const { minimumSpanningTree, backEdges } =
            MinimumSpanningTree.execute([[7, 4]], [2.5], true);
        expect(minimumSpanningTree).toEqual([[nil, 7], [7, 4]]);
        expect(backEdges).toEqual([]);
    });

    it('maps non-consecutive vertex indices back to the originals', () => {
        // The vertices are 10, 20, 30; only the edges name them.
        const sparseEdges: MSTEdge[] = [[10, 20], [20, 30], [30, 10]];
        const sparseWeights = [5, 1, 7];
        const { minimumSpanningTree, backEdges } =
            MinimumSpanningTree.execute(sparseEdges, sparseWeights, true);

        expect(minimumSpanningTree.length).toBe(3);
        // The starting vertex is the first vertex encountered, that is, 10.
        expect(minimumSpanningTree[0]).toEqual([nil, 10]);
        expect(treeWeight(minimumSpanningTree, sparseEdges, sparseWeights)).toBe(6);
        expect(backEdges).toEqual([[10, 30]]);
    });
});

describe('MinimumSpanningTree structural properties', () => {
    it('produces a forest with one sentinel per connected component', () => {
        // Two disjoint triangles.
        const edges: MSTEdge[] = [
            [0, 1], [1, 2], [2, 0],
            [3, 4], [4, 5], [5, 3]
        ];
        const weights = [1, 2, 3, 4, 5, 6];
        const { minimumSpanningTree, backEdges } =
            MinimumSpanningTree.execute(edges, weights, true);

        const sentinels = minimumSpanningTree.filter(e => e[0] === nil);
        expect(sentinels.length).toBe(2);
        expect(minimumSpanningTree.length).toBe(6);
        // Each component keeps its two cheapest edges: 1+2 and 4+5.
        expect(treeWeight(minimumSpanningTree, edges, weights)).toBe(12);
        expect(backEdges.length).toBe(2);
    });

    it('has no cycles and spans every vertex', () => {
        const rng = makeRng(0xa5a5f00d);
        for (let trial = 0; trial < 25; ++trial) {
            const numVertices = 4 + Math.floor(8 * rng());
            const edges: MSTEdge[] = [];
            const weights: number[] = [];
            // A random spanning path guarantees connectivity, then random
            // extra edges are added.
            for (let v = 1; v < numVertices; ++v) {
                edges.push([v - 1, v]);
                weights.push(1 + 9 * rng());
            }
            const present = new Set<string>(edges.map(e => e[0] + ',' + e[1]));
            for (let k = 0; k < numVertices; ++k) {
                const a = Math.floor(numVertices * rng());
                const b = Math.floor(numVertices * rng());
                if (a === b) {
                    continue;
                }
                const key = Math.min(a, b) + ',' + Math.max(a, b);
                if (present.has(key)) {
                    continue;
                }
                present.add(key);
                edges.push([a, b]);
                weights.push(1 + 9 * rng());
            }

            const { minimumSpanningTree, backEdges } =
                MinimumSpanningTree.execute(edges, weights, true);

            // Exactly one sentinel for a connected graph, and one output
            // element per vertex.
            expect(minimumSpanningTree.length).toBe(numVertices);
            expect(minimumSpanningTree.filter(e => e[0] === nil).length).toBe(1);

            // The tree reaches every vertex exactly once as a destination.
            const destinations = minimumSpanningTree.map(e => e[1]);
            expect(new Set(destinations).size).toBe(numVertices);

            // Every non-sentinel edge connects to a vertex already reached,
            // so the output has no cycles.
            const reached = new Set<number>();
            for (const treeEdge of minimumSpanningTree) {
                if (treeEdge[0] === nil) {
                    reached.add(treeEdge[1]);
                } else {
                    expect(reached.has(treeEdge[0])).toBe(true);
                    expect(reached.has(treeEdge[1])).toBe(false);
                    reached.add(treeEdge[1]);
                }
            }

            // Tree edges plus back edges account for every input edge.
            expect(minimumSpanningTree.length - 1 + backEdges.length)
                .toBe(edges.length);
            for (const backEdge of backEdges) {
                expect(backEdge[0]).toBeLessThan(backEdge[1]);
                expect(present.has(backEdge[0] + ',' + backEdge[1])).toBe(true);
            }
        }
    });

    it('matches the total weight computed by Kruskal', () => {
        const rng = makeRng(0x1337beef);
        for (let trial = 0; trial < 40; ++trial) {
            const numVertices = 3 + Math.floor(9 * rng());
            const edges: MSTEdge[] = [];
            const weights: number[] = [];
            const present = new Set<string>();
            // A dense random graph, not necessarily connected.
            for (let a = 0; a < numVertices; ++a) {
                for (let b = a + 1; b < numVertices; ++b) {
                    if (rng() < 0.55) {
                        present.add(a + ',' + b);
                        edges.push([a, b]);
                        weights.push(1 + 20 * rng());
                    }
                }
            }
            if (edges.length === 0) {
                continue;
            }

            const { minimumSpanningTree } =
                MinimumSpanningTree.execute(edges, weights, true);
            const prim = treeWeight(minimumSpanningTree, edges, weights);
            expect(prim).toBeCloseTo(kruskalWeight(edges, weights), 12);

            // The output covers exactly the referenced vertices.
            expect(minimumSpanningTree.length).toBe(vertexCount(edges));
        }
    });
});

describe('MinimumSpanningTree input validation', () => {
    it('requires matching edge and weight counts', () => {
        expect(() => MinimumSpanningTree.execute([[0, 1]], [1, 2], true))
            .toThrow('The edges.size() and weights.size() must match.');
    });

    it('rejects negative weights', () => {
        expect(() => MinimumSpanningTree.execute([[0, 1]], [-1], true))
            .toThrow('Encountered a negative weight.');
    });

    it('rejects degenerate edges', () => {
        expect(() => MinimumSpanningTree.execute([[2, 2]], [1], true))
            .toThrow('Encountered a degenerate edge.');
        expect(() => MinimumSpanningTree.execute([[nil, 1]], [1], true))
            .toThrow('Encountered a degenerate edge.');
    });

    it('rejects duplicate edges, in either direction', () => {
        expect(() => MinimumSpanningTree.execute(
            [[0, 1], [1, 2], [1, 0]], [1, 2, 3], true))
            .toThrow('Encountered a duplicate edge.');
        expect(() => MinimumSpanningTree.execute(
            [[0, 1], [0, 1]], [1, 2], true))
            .toThrow('Encountered a duplicate edge.');
    });

    it('detects duplicate edges even without validation', () => {
        // The internal assertions still fire, because the edge map cannot
        // hold two entries for the same remapped directed edge.
        expect(() => MinimumSpanningTree.execute(
            [[0, 1], [1, 0]], [1, 2], false))
            .toThrow('Unexpected result for validated edges.');
    });

    it('accepts zero weights', () => {
        const { minimumSpanningTree } =
            MinimumSpanningTree.execute([[0, 1], [1, 2]], [0, 0], true);
        expect(minimumSpanningTree.length).toBe(3);
    });
});

describe('MinimumSpanningTree verification', () => {
    // A random graph on 'labels' with integer weights. Non-consecutive and
    // out-of-order vertex labels exercise the vertex remapping. Integer
    // weights make ties (and therefore the tie-breaking paths) common.
    const graph = fc.tuple(
        fc.uniqueArray(fc.integer({ min: 0, max: 200 }),
            { minLength: 2, maxLength: 7 }),
        fc.array(fc.tuple(fc.nat(), fc.nat(), fc.integer({ min: 0, max: 4 })),
            { minLength: 1, maxLength: 24 }))
        .map(([labels, raw]) => {
            const edges: MSTEdge[] = [];
            const weights: number[] = [];
            const seen = new Set<string>();
            for (const [a, b, w] of raw) {
                const i = a % labels.length;
                const j = b % labels.length;
                if (i === j) {
                    continue;
                }
                const key = Math.min(i, j) + ',' + Math.max(i, j);
                if (seen.has(key)) {
                    continue;
                }
                seen.add(key);
                edges.push([labels[i], labels[j]]);
                weights.push(w);
            }
            return { edges, weights };
        })
        .filter(({ edges }) => edges.length > 0);

    it('the forest weight equals the Kruskal weight', () => {
        check(graph, ({ edges, weights }) => {
            const { minimumSpanningTree } =
                MinimumSpanningTree.execute(edges, weights, true);
            expect(treeWeight(minimumSpanningTree, edges, weights))
                .toBe(kruskalWeight(edges, weights));
        });
    });

    it('the output is a spanning forest of the referenced vertices', () => {
        check(graph, ({ edges, weights }) => {
            const { minimumSpanningTree } =
                MinimumSpanningTree.execute(edges, weights, true);

            // One output element per referenced vertex, each vertex reached
            // exactly once, and every non-sentinel edge attaches to an
            // already-reached vertex (so the forest is acyclic).
            expect(minimumSpanningTree.length).toBe(vertexCount(edges));
            const reached = new Set<number>();
            for (const treeEdge of minimumSpanningTree) {
                if (treeEdge[0] === nil) {
                    expect(reached.has(treeEdge[1])).toBe(false);
                } else {
                    expect(reached.has(treeEdge[0])).toBe(true);
                    expect(reached.has(treeEdge[1])).toBe(false);
                }
                reached.add(treeEdge[1]);
            }
            const referenced = new Set<number>();
            for (const edge of edges) {
                referenced.add(edge[0]);
                referenced.add(edge[1]);
            }
            expect(reached.size).toBe(referenced.size);
            for (const v of referenced) {
                expect(reached.has(v)).toBe(true);
            }
        });
    });

    it('the sentinel count equals the number of connected components', () => {
        check(graph, ({ edges, weights }) => {
            const parent = new Map<number, number>();
            const find = (v: number): number => {
                while (parent.get(v) !== v) {
                    v = parent.get(v) as number;
                }
                return v;
            };
            for (const edge of edges) {
                for (const v of edge) {
                    if (!parent.has(v)) {
                        parent.set(v, v);
                    }
                }
            }
            for (const [v0, v1] of edges) {
                parent.set(find(v0), find(v1));
            }
            const roots = new Set<number>();
            for (const v of parent.keys()) {
                roots.add(find(v));
            }

            const { minimumSpanningTree } =
                MinimumSpanningTree.execute(edges, weights, true);
            expect(minimumSpanningTree.filter((e) => e[0] === nil).length)
                .toBe(roots.size);
        });
    });

    it('tree edges and back edges partition the input edges', () => {
        check(graph, ({ edges, weights }) => {
            const { minimumSpanningTree, backEdges } =
                MinimumSpanningTree.execute(edges, weights, true);

            const undirected = (e: MSTEdge): string =>
                Math.min(e[0], e[1]) + ',' + Math.max(e[0], e[1]);
            const input = new Set(edges.map(undirected));
            const tree = minimumSpanningTree
                .filter((e) => e[0] !== nil)
                .map(undirected);
            const back = backEdges.map(undirected);

            expect(new Set(tree).size).toBe(tree.length);
            expect(new Set(back).size).toBe(back.length);
            expect(tree.length + back.length).toBe(input.size);
            for (const key of tree.concat(back)) {
                expect(input.has(key)).toBe(true);
            }
        });
    });

    it('back edges are ordered by the internal remapped vertex indices', () => {
        // Upstream keeps the (v0,v1) of a back-edge pair for which v0 < v1
        // in the *remapped* indices (assigned in first-encounter order over
        // edges[]) and emits them in std::map order, then converts back to
        // the original labels. So the reported back edges are neither
        // sorted by the original labels nor necessarily smaller-first.
        check(graph, ({ edges, weights }) => {
            const remap = new Map<number, number>();
            for (const edge of edges) {
                for (const v of edge) {
                    if (!remap.has(v)) {
                        remap.set(v, remap.size);
                    }
                }
            }
            const inverse: number[] = [];
            for (const [label, index] of remap) {
                inverse[index] = label;
            }

            const { minimumSpanningTree, backEdges } =
                MinimumSpanningTree.execute(edges, weights, true);
            const undirected = (e: MSTEdge): string =>
                Math.min(e[0], e[1]) + ',' + Math.max(e[0], e[1]);
            const treeSet = new Set(minimumSpanningTree
                .filter((e) => e[0] !== nil).map(undirected));

            const expected: MSTEdge[] = [];
            for (const edge of edges) {
                if (treeSet.has(undirected(edge))) {
                    continue;
                }
                const r0 = remap.get(edge[0]) as number;
                const r1 = remap.get(edge[1]) as number;
                expected.push([Math.min(r0, r1), Math.max(r0, r1)]);
            }
            expected.sort((a, b) => (a[0] !== b[0] ? a[0] - b[0] : a[1] - b[1]));
            expect(backEdges).toEqual(
                expected.map((e) => [inverse[e[0]], inverse[e[1]]]));
        });
    });

    it('every tree edge is a minimum-weight edge across its own cut', () => {
        check(graph, ({ edges, weights }) => {
            const { minimumSpanningTree } =
                MinimumSpanningTree.execute(edges, weights, true);
            const weightOf = new Map<string, number>();
            for (let e = 0; e < edges.length; ++e) {
                const key = Math.min(edges[e][0], edges[e][1]) + ',' +
                    Math.max(edges[e][0], edges[e][1]);
                weightOf.set(key, weights[e]);
            }

            // Prim grows the tree in output order, so the vertices already
            // added when a tree edge is emitted define one side of the cut.
            const inTree = new Set<number>();
            for (const treeEdge of minimumSpanningTree) {
                if (treeEdge[0] !== nil) {
                    const key = Math.min(treeEdge[0], treeEdge[1]) + ',' +
                        Math.max(treeEdge[0], treeEdge[1]);
                    const chosen = weightOf.get(key) as number;
                    for (const [v0, v1] of edges) {
                        const crosses = (inTree.has(v0) && !inTree.has(v1)) ||
                            (inTree.has(v1) && !inTree.has(v0));
                        if (crosses) {
                            const w = weightOf.get(
                                Math.min(v0, v1) + ',' + Math.max(v0, v1)) as number;
                            expect(chosen).toBeLessThanOrEqual(w);
                        }
                    }
                }
                inTree.add(treeEdge[1]);
            }
        });
    });

    it('validation is optional but does not change the result', () => {
        check(graph, ({ edges, weights }) => {
            const validated = MinimumSpanningTree.execute(edges, weights, true);
            const unvalidated = MinimumSpanningTree.execute(edges, weights, false);
            expect(unvalidated.minimumSpanningTree)
                .toEqual(validated.minimumSpanningTree);
            expect(unvalidated.backEdges).toEqual(validated.backEdges);
        });
    });

    it('execute does not modify its inputs', () => {
        check(graph, ({ edges, weights }) => {
            const edgesCopy = edges.map((e) => e.slice());
            const weightsCopy = weights.slice();
            MinimumSpanningTree.execute(edges, weights, true);
            expect(edges.map((e) => e.slice())).toEqual(edgesCopy);
            expect(weights).toEqual(weightsCopy);
        });
    });

    it('an empty edge list yields empty output instead of upstream UB', () => {
        // Upstream sizes records[] by the number of vertices (zero here) and
        // then writes records[0], which is undefined behavior. The port
        // returns early.
        const result = MinimumSpanningTree.execute([], [], true);
        expect(result.minimumSpanningTree).toEqual([]);
        expect(result.backEdges).toEqual([]);
        const unvalidated = MinimumSpanningTree.execute([], [], false);
        expect(unvalidated.minimumSpanningTree).toEqual([]);
        expect(unvalidated.backEdges).toEqual([]);
    });
});
