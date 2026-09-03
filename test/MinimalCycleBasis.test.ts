import { describe, expect, it } from 'vitest';
import {
    MinimalCycleBasis
} from '../src/MinimalCycleBasis.js';
import type {
    MinimalCycleBasisEdge, MinimalCycleBasisPosition, MinimalCycleBasisTree
} from '../src/MinimalCycleBasis.js';

// A small deterministic pseudorandom generator (mulberry32) so the randomized
// cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function extract(positions: MinimalCycleBasisPosition[],
    edges: MinimalCycleBasisEdge[]): MinimalCycleBasis {
    const mcb = new MinimalCycleBasis();
    mcb.extract(positions, edges, true);
    return mcb;
}

// All non-empty cycles of the forest, in depth-first order.
function collectCycles(trees: MinimalCycleBasisTree[]): number[][] {
    const cycles: number[][] = [];
    const visit = (tree: MinimalCycleBasisTree): void => {
        if (tree.cycle.length > 0) {
            cycles.push(tree.cycle);
        }
        for (const child of tree.children) {
            visit(child);
        }
    };
    for (const tree of trees) {
        visit(tree);
    }
    return cycles;
}

// A rotation- and reflection-invariant key for a cycle whose first and last
// entries are the same vertex.
function cycleKey(cycle: number[]): string {
    expect(cycle.length).toBeGreaterThanOrEqual(4);
    expect(cycle[0]).toBe(cycle[cycle.length - 1]);
    const loop = cycle.slice(0, cycle.length - 1);
    const n = loop.length;
    const candidates: string[] = [];
    for (const seq of [loop, loop.slice().reverse()]) {
        for (let s = 0; s < n; ++s) {
            const rotated: number[] = [];
            for (let i = 0; i < n; ++i) {
                rotated.push(seq[(s + i) % n] as number);
            }
            candidates.push(rotated.join(','));
        }
    }
    candidates.sort();
    return candidates[0] as string;
}

function cycleKeys(trees: MinimalCycleBasisTree[]): string[] {
    const keys = collectCycles(trees).map(cycleKey);
    keys.sort();
    return keys;
}

function edgeName(v0: number, v1: number): string {
    return v0 < v1 ? v0 + ',' + v1 : v1 + ',' + v0;
}

function edgeNameSet(edges: readonly MinimalCycleBasisEdge[]): Set<string> {
    const set = new Set<string>();
    for (const edge of edges) {
        set.add(edgeName(edge[0], edge[1]));
    }
    return set;
}

// The number of independent cycles of a graph is E - V + C, where C counts
// the connected components (an isolated vertex is its own component, and its
// contributions to V and C cancel).
function cycleRank(numVertices: number,
    edges: readonly MinimalCycleBasisEdge[]): number {
    const parent = new Array<number>(numVertices);
    for (let i = 0; i < numVertices; ++i) {
        parent[i] = i;
    }
    const find = (i: number): number => {
        while ((parent[i] as number) !== i) {
            parent[i] = parent[parent[i] as number] as number;
            i = parent[i] as number;
        }
        return i;
    };
    for (const edge of edges) {
        const r0 = find(edge[0]);
        const r1 = find(edge[1]);
        if (r0 !== r1) {
            parent[r0] = r1;
        }
    }
    const roots = new Set<number>();
    for (let i = 0; i < numVertices; ++i) {
        roots.add(find(i));
    }
    return edges.length - numVertices + roots.size;
}

// The m-by-n grid of unit cells; vertices are laid out row by row.
function makeGrid(numCellsX: number, numCellsY: number): {
    positions: MinimalCycleBasisPosition[],
    edges: MinimalCycleBasisEdge[],
    indexOf: (i: number, j: number) => number
} {
    const numX = numCellsX + 1;
    const numY = numCellsY + 1;
    const indexOf = (i: number, j: number): number => j * numX + i;
    const positions: MinimalCycleBasisPosition[] = [];
    for (let j = 0; j < numY; ++j) {
        for (let i = 0; i < numX; ++i) {
            positions.push([i, j]);
        }
    }
    const edges: MinimalCycleBasisEdge[] = [];
    for (let j = 0; j < numY; ++j) {
        for (let i = 0; i + 1 < numX; ++i) {
            edges.push([indexOf(i, j), indexOf(i + 1, j)]);
        }
    }
    for (let j = 0; j + 1 < numY; ++j) {
        for (let i = 0; i < numX; ++i) {
            edges.push([indexOf(i, j), indexOf(i, j + 1)]);
        }
    }
    return { positions, edges, indexOf };
}

describe('MinimalCycleBasis', () => {
    it('reports nothing for empty inputs', () => {
        const mcb = new MinimalCycleBasis();
        mcb.extract([], [], true);
        expect(mcb.getIsolatedVertices()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([]);
        expect(mcb.getForest()).toEqual([]);

        // Upstream returns early when there are no edges, so the vertices are
        // not reported as isolated. The port preserves that behavior.
        mcb.extract([[0, 0], [1, 1]], [], true);
        expect(mcb.getIsolatedVertices()).toEqual([]);
        expect(mcb.getForest()).toEqual([]);
    });

    it('validates the inputs when asked to', () => {
        const mcb = new MinimalCycleBasis();
        expect(() => mcb.extract([[0, 0], [0, 0]], [[0, 1]], true))
            .toThrowError('Input positions must be unique.');
        expect(() => mcb.extract([[0, 0], [1, 0]], [[0, 2]], true))
            .toThrowError('Input index edge[0][1] is out of range.');
        expect(() => mcb.extract([[0, 0], [1, 0]], [[-1, 1]], true))
            .toThrowError('Input index edge[0][0] is out of range.');
        expect(() => mcb.extract([[0, 0], [1, 0]], [[1, 1]], true))
            .toThrowError('Input edge[0] is degenerate.');

        // The same graph is accepted without verification.
        expect(() => mcb.extract([[0, 0], [0, 0]], [[0, 1]], false))
            .not.toThrow();
    });

    it('extracts a single square as one cycle', () => {
        const mcb = extract([[0, 0], [1, 0], [1, 1], [0, 1]],
            [[0, 1], [1, 2], [2, 3], [3, 0]]);
        expect(mcb.getIsolatedVertices()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([]);
        const forest = mcb.getForest();
        expect(forest.length).toBe(1);
        expect((forest[0] as MinimalCycleBasisTree).children).toEqual([]);
        // The closed walk repeats the starting vertex.
        expect((forest[0] as MinimalCycleBasisTree).cycle).toEqual([0, 1, 2, 3, 0]);
    });

    it('extracts a triangle with a spike as one cycle and one filament', () => {
        const mcb = extract([[0, 0], [2, 0], [1, 2], [3, 3]],
            [[0, 1], [1, 2], [2, 0], [2, 3]]);
        expect(mcb.getIsolatedVertices()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([[3, 2]]);
        expect(cycleKeys(mcb.getForest())).toEqual([cycleKey([0, 1, 2, 0])]);
    });

    it('reports isolated vertices', () => {
        const mcb = extract([[0, 0], [1, 0], [1, 1], [0, 1], [5, 5], [7, 7]],
            [[0, 1], [1, 2], [2, 3], [3, 0]]);
        expect(mcb.getIsolatedVertices()).toEqual([4, 5]);
        expect(cycleKeys(mcb.getForest())).toEqual([cycleKey([0, 1, 2, 3, 0])]);
    });

    it('extracts filament chains with branch points', () => {
        // An open polyline has no cycles.
        let mcb = extract([[0, 0], [1, 0], [2, 0], [3, 0]],
            [[0, 1], [1, 2], [2, 3]]);
        expect(mcb.getForest()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([[0, 1, 2, 3]]);

        // A branched polyline: the greedy removal splits it into filaments
        // that together cover every edge exactly once.
        mcb = extract([[0, 0], [1, 0], [2, 0], [2, 1], [2, -1]],
            [[0, 1], [1, 2], [2, 3], [2, 4]]);
        expect(mcb.getForest()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([[0, 1, 2], [3, 2, 4]]);

        // A star: the center is the only vertex with two or more adjacents.
        mcb = extract([[0, 0], [1, 0], [-1, 0], [0, 1]],
            [[0, 1], [0, 2], [0, 3]]);
        expect(mcb.getForest()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([[1, 0], [2, 0, 3]]);

        // Upstream visits only components reachable from a vertex with at
        // least two adjacents, so a lone edge is reported as nothing at all.
        // The port preserves this.
        mcb = extract([[0, 0], [1, 0]], [[0, 1]]);
        expect(mcb.getForest()).toEqual([]);
        expect(mcb.getFilaments()).toEqual([]);
        expect(mcb.getIsolatedVertices()).toEqual([]);
    });

    it('extracts the unit cells of a grid', () => {
        for (const [nx, ny] of [[2, 2], [3, 2], [3, 3]] as [number, number][]) {
            const { positions, edges, indexOf } = makeGrid(nx, ny);
            const mcb = extract(positions, edges);
            expect(mcb.getIsolatedVertices()).toEqual([]);
            expect(mcb.getFilaments()).toEqual([]);

            const expected: string[] = [];
            for (let j = 0; j < ny; ++j) {
                for (let i = 0; i < nx; ++i) {
                    expected.push(cycleKey([
                        indexOf(i, j), indexOf(i + 1, j),
                        indexOf(i + 1, j + 1), indexOf(i, j + 1),
                        indexOf(i, j)
                    ]));
                }
            }
            expected.sort();
            expect(cycleKeys(mcb.getForest())).toEqual(expected);

            // The grid is a single connected component, so the forest holds
            // one tree whose children are the cells (the top-level tree has
            // no cycle of its own when there is more than one cell).
            expect(mcb.getForest().length).toBe(1);
        }
    });

    it('extracts two squares sharing an edge', () => {
        const mcb = extract(
            [[0, 0], [2, 0], [2, 2], [0, 2], [4, 0], [4, 2]],
            [[0, 1], [1, 2], [2, 3], [3, 0], [1, 4], [4, 5], [5, 2]]);
        expect(mcb.getFilaments()).toEqual([]);
        expect(cycleKeys(mcb.getForest())).toEqual([
            cycleKey([0, 1, 2, 3, 0]), cycleKey([1, 4, 5, 2, 1])
        ].sort());
        const forest = mcb.getForest();
        expect(forest.length).toBe(1);
        expect((forest[0] as MinimalCycleBasisTree).cycle).toEqual([]);
        expect((forest[0] as MinimalCycleBasisTree).children.length).toBe(2);
    });

    it('extracts two squares sharing a vertex', () => {
        const mcb = extract(
            [[0, 0], [2, 0], [2, 2], [0, 2], [4, 2], [4, 4], [2, 4]],
            [[0, 1], [1, 2], [2, 3], [3, 0], [2, 4], [4, 5], [5, 6], [6, 2]]);
        expect(mcb.getFilaments()).toEqual([]);
        expect(cycleKeys(mcb.getForest())).toEqual([
            cycleKey([0, 1, 2, 3, 0]), cycleKey([2, 4, 5, 6, 2])
        ].sort());
    });

    it('puts disconnected components in separate trees', () => {
        const mcb = extract(
            [[0, 0], [1, 0], [1, 1], [0, 1], [5, 0], [6, 0], [6, 1], [5, 1]],
            [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]);
        const forest = mcb.getForest();
        expect(forest.length).toBe(2);
        expect((forest[0] as MinimalCycleBasisTree).cycle).toEqual([0, 1, 2, 3, 0]);
        expect((forest[0] as MinimalCycleBasisTree).children).toEqual([]);
        expect((forest[1] as MinimalCycleBasisTree).cycle).toEqual([4, 5, 6, 7, 4]);
        expect((forest[1] as MinimalCycleBasisTree).children).toEqual([]);

        // A square strictly inside another square, but not connected to it,
        // is a separate component and therefore a sibling, not a child.
        const separate = extract(
            [[0, 0], [4, 0], [4, 4], [0, 4], [1, 1], [3, 1], [3, 3], [1, 3]],
            [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4]]);
        expect(separate.getForest().length).toBe(2);
        expect(separate.getForest().every(t => t.children.length === 0))
            .toBe(true);
    });

    it('nests a cycle inside the cycle that surrounds it', () => {
        // An outer square, an inner square, and one edge joining them. The
        // joining edge is a filament and the inner square is a child of the
        // outer one.
        const mcb = extract(
            [[0, 0], [4, 0], [4, 4], [0, 4], [1, 1], [3, 1], [3, 3], [1, 3]],
            [[0, 1], [1, 2], [2, 3], [3, 0], [4, 5], [5, 6], [6, 7], [7, 4],
                [0, 4]]);
        expect(mcb.getFilaments()).toEqual([[0, 4]]);
        const forest = mcb.getForest();
        expect(forest.length).toBe(1);
        const outer = forest[0] as MinimalCycleBasisTree;
        expect(cycleKey(outer.cycle)).toBe(cycleKey([0, 1, 2, 3, 0]));
        expect(outer.children.length).toBe(1);
        const inner = outer.children[0] as MinimalCycleBasisTree;
        expect(cycleKey(inner.cycle)).toBe(cycleKey([4, 5, 6, 7, 4]));
        expect(inner.children).toEqual([]);
    });

    it('nests a doubly nested cycle two levels deep', () => {
        // Regression test for the upstream 'visited' flag bug described in
        // MinimalCycleBasis.ts: with the persistent flags the innermost
        // square is reported as a sibling of the middle square (both under an
        // empty-cycle node) instead of as its child.
        const mcb = extract(
            [[0, 0], [8, 0], [8, 8], [0, 8],
                [1, 1], [7, 1], [7, 7], [1, 7],
                [2, 2], [6, 2], [6, 6], [2, 6]],
            [[0, 1], [1, 2], [2, 3], [3, 0],
                [4, 5], [5, 6], [6, 7], [7, 4],
                [8, 9], [9, 10], [10, 11], [11, 8],
                [0, 4], [4, 8]]);
        expect(mcb.getFilaments()).toEqual([[0, 4], [4, 8]]);
        const forest = mcb.getForest();
        expect(forest.length).toBe(1);
        const outer = forest[0] as MinimalCycleBasisTree;
        expect(cycleKey(outer.cycle)).toBe(cycleKey([0, 1, 2, 3, 0]));
        expect(outer.children.length).toBe(1);
        const middle = outer.children[0] as MinimalCycleBasisTree;
        expect(cycleKey(middle.cycle)).toBe(cycleKey([4, 5, 6, 7, 4]));
        expect(middle.children.length).toBe(1);
        const inner = middle.children[0] as MinimalCycleBasisTree;
        expect(cycleKey(inner.cycle)).toBe(cycleKey([8, 9, 10, 11, 8]));
        expect(inner.children).toEqual([]);
    });

    it('handles collinear vertices along a cycle', () => {
        // A square whose bottom edge is subdivided twice; the extra vertices
        // exercise the exact-arithmetic convexity tests with zero
        // determinants.
        const mcb = extract(
            [[0, 0], [1, 0], [2, 0], [3, 0], [3, 3], [0, 3]],
            [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0]]);
        expect(mcb.getFilaments()).toEqual([]);
        expect(cycleKeys(mcb.getForest()))
            .toEqual([cycleKey([0, 1, 2, 3, 4, 5, 0])]);
    });

    it('agrees with the cycle rank on random grid subgraphs', () => {
        const random = makeRandom(0x5EED1234);
        for (let trial = 0; trial < 40; ++trial) {
            const { positions, edges: allEdges } = makeGrid(3, 3);
            const edges = allEdges.filter(() => random() < 0.75);
            if (edges.length === 0) {
                continue;
            }
            const mcb = extract(positions, edges);

            const names = edgeNameSet(edges);
            const cycles = collectCycles(mcb.getForest());

            // The number of extracted cycles is the cycle rank E - V + C.
            expect(cycles.length).toBe(cycleRank(positions.length, edges));

            // Every cycle is a simple closed walk over graph edges.
            const keys = new Set<string>();
            for (const cycle of cycles) {
                expect(cycle.length).toBeGreaterThanOrEqual(4);
                expect(cycle[0]).toBe(cycle[cycle.length - 1]);
                const seen = new Set<number>();
                for (let i = 0; i + 1 < cycle.length; ++i) {
                    expect(seen.has(cycle[i] as number)).toBe(false);
                    seen.add(cycle[i] as number);
                    expect(names.has(
                        edgeName(cycle[i] as number, cycle[i + 1] as number)))
                        .toBe(true);
                }
                const key = cycleKey(cycle);
                expect(keys.has(key)).toBe(false);
                keys.add(key);
            }

            // Every filament is a polyline over graph edges.
            for (const filament of mcb.getFilaments()) {
                expect(filament.length).toBeGreaterThanOrEqual(2);
                for (let i = 0; i + 1 < filament.length; ++i) {
                    expect(names.has(
                        edgeName(filament[i] as number,
                            filament[i + 1] as number))).toBe(true);
                }
            }

            // The isolated vertices are exactly the vertices with no edge.
            const incident = new Set<number>();
            for (const edge of edges) {
                incident.add(edge[0]);
                incident.add(edge[1]);
            }
            const expectedIsolated: number[] = [];
            for (let i = 0; i < positions.length; ++i) {
                if (!incident.has(i)) {
                    expectedIsolated.push(i);
                }
            }
            expect(mcb.getIsolatedVertices()).toEqual(expectedIsolated);
        }
    });

    it('agrees with the cycle rank on random polygon-with-chords graphs', () => {
        const random = makeRandom(0xC0FFEE);
        for (let trial = 0; trial < 25; ++trial) {
            // Vertices on a circle plus a hub at the center; the spokes are
            // chosen at random. The graph is planar by construction.
            const numRim = 6 + (trial % 5);
            const positions: MinimalCycleBasisPosition[] = [];
            for (let i = 0; i < numRim; ++i) {
                const angle = 2 * Math.PI * i / numRim;
                positions.push([Math.cos(angle), Math.sin(angle)]);
            }
            positions.push([0, 0]);
            const hub = numRim;

            const edges: MinimalCycleBasisEdge[] = [];
            for (let i = 0; i < numRim; ++i) {
                edges.push([i, (i + 1) % numRim]);
            }
            for (let i = 0; i < numRim; ++i) {
                if (random() < 0.6) {
                    edges.push([i, hub]);
                }
            }

            const mcb = extract(positions, edges);
            const cycles = collectCycles(mcb.getForest());
            expect(cycles.length).toBe(cycleRank(positions.length, edges));

            const names = edgeNameSet(edges);
            for (const cycle of cycles) {
                expect(cycle[0]).toBe(cycle[cycle.length - 1]);
                for (let i = 0; i + 1 < cycle.length; ++i) {
                    expect(names.has(
                        edgeName(cycle[i] as number, cycle[i + 1] as number)))
                        .toBe(true);
                }
            }
        }
    });
});
