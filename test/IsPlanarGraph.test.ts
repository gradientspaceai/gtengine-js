import { describe, it, expect } from 'vitest';
import { IsPlanarGraph, OrderedEdge } from '../src/IsPlanarGraph.js';

type P2 = [number, number];
type E2 = [number, number];

describe('OrderedEdge', () => {
    it('stores the pair as (min, max)', () => {
        expect(new OrderedEdge(3, 7).v).toEqual([3, 7]);
        expect(new OrderedEdge(7, 3).v).toEqual([3, 7]);
        expect(new OrderedEdge().v).toEqual([-1, -1]);
    });

    it('orders lexicographically', () => {
        expect(new OrderedEdge(0, 1).lessThan(new OrderedEdge(0, 2))).toBe(true);
        expect(new OrderedEdge(0, 2).lessThan(new OrderedEdge(1, 2))).toBe(true);
        expect(new OrderedEdge(1, 2).lessThan(new OrderedEdge(1, 2))).toBe(false);
        expect(new OrderedEdge(2, 3).lessThan(new OrderedEdge(1, 3))).toBe(false);
    });
});

describe('IsPlanarGraph', () => {
    it('accepts a triangle as planar', () => {
        const positions: P2[] = [[0, 0], [1, 0], [0, 1]];
        const edges: E2[] = [[0, 1], [1, 2], [2, 0]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
        expect(ipg.getInvalidIntersections()).toEqual([]);
    });

    it('accepts a planar straight-line embedding of K4', () => {
        // Outer triangle plus its centroid, all six edges: the standard
        // planar embedding of the complete graph K4.
        const positions: P2[] = [[0, 0], [3, 0], [0, 3], [1, 1]];
        const edges: E2[] = [[0, 1], [1, 2], [2, 0], [0, 3], [1, 3], [2, 3]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
        expect(ipg.getInvalidIntersections()).toEqual([]);
    });

    it('rejects K4 drawn on a convex quadrilateral (crossing diagonals)', () => {
        // The same abstract graph as the planar K4 test, but drawn with all
        // four vertices in convex position: the diagonals (0,2) and (1,3)
        // cross at an edge-interior point.
        const positions: P2[] = [[0, 0], [1, 0], [1, 1], [0, 1]];
        const edges: E2[] = [[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [1, 3]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        const invalid = ipg.getInvalidIntersections();
        expect(invalid.length).toBe(1);
        // Edge indices 4 and 5 are the two diagonals, reported with
        // element.v[0] < element.v[1].
        expect(invalid[0].v).toEqual([4, 5]);
    });

    it('rejects two crossing segments', () => {
        const positions: P2[] = [[0, 0], [1, 1], [0, 1], [1, 0]];
        const edges: E2[] = [[0, 1], [2, 3]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        expect(ipg.getInvalidIntersections().map((e) => e.v)).toEqual([[0, 1]]);
    });

    it('rejects a T-junction (endpoint interior to another edge)', () => {
        // Vertex 2 lies at (1,0), the interior of edge (0,1). The edge
        // (2,3) touches edge (0,1) at that interior point, which is an
        // invalid intersection even though it occurs at an endpoint of one
        // of the segments.
        const positions: P2[] = [[0, 0], [2, 0], [1, 0.0000001], [1, 1]];
        const edges: E2[] = [[0, 1], [2, 3]];
        const ipgSeparated = new IsPlanarGraph();
        expect(ipgSeparated.compute(positions, edges))
            .toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);

        positions[2] = [1, 0];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        expect(ipg.getInvalidIntersections().map((e) => e.v)).toEqual([[0, 1]]);
    });

    it('accepts segments sharing an endpoint', () => {
        const positions: P2[] = [[0, 0], [1, 0], [1, 1], [2, 0]];
        const edges: E2[] = [[0, 1], [1, 2], [1, 3]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
    });

    it('rejects collinear segments overlapping in an interval', () => {
        const positions: P2[] = [[0, 0], [2, 0], [1, 0], [3, 0]];
        const edges: E2[] = [[0, 1], [2, 3]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        expect(ipg.getInvalidIntersections().map((e) => e.v)).toEqual([[0, 1]]);
    });

    it('accepts collinear segments meeting only at a shared vertex', () => {
        const positions: P2[] = [[0, 0], [1, 0], [2, 0]];
        const edges: E2[] = [[0, 1], [1, 2]];
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
    });

    it('flags invalid input sizes', () => {
        const ipg = new IsPlanarGraph();
        expect(ipg.compute([], [])).toBe(IsPlanarGraph.IPG_INVALID_INPUT_SIZES);
        expect(ipg.compute([[0, 0], [1, 0]], []))
            .toBe(IsPlanarGraph.IPG_INVALID_INPUT_SIZES);
        expect(ipg.compute([[0, 0]], [[0, 0]]))
            .toBe(IsPlanarGraph.IPG_INVALID_INPUT_SIZES);
    });

    it('flags duplicated positions and reports the groups in sorted order', () => {
        const positions: P2[] = [[1, 1], [0, 0], [1, 1], [2, 0], [0, 0], [0, 0]];
        const edges: E2[] = [[1, 3]];
        const ipg = new IsPlanarGraph();
        const flags = ipg.compute(positions, edges);
        expect(flags & IsPlanarGraph.IPG_DUPLICATED_POSITIONS)
            .toBe(IsPlanarGraph.IPG_DUPLICATED_POSITIONS);
        // Groups ordered by position: (0,0) before (1,1), indices in input
        // order within a group (the upstream std::map iteration order).
        expect(ipg.getDuplicatedPositions()).toEqual([[1, 4, 5], [0, 2]]);
    });

    it('flags duplicated edges including reversed duplicates', () => {
        const positions: P2[] = [[0, 0], [1, 0], [0, 1]];
        const edges: E2[] = [[0, 1], [1, 2], [1, 0], [2, 1], [2, 0]];
        const ipg = new IsPlanarGraph();
        const flags = ipg.compute(positions, edges);
        expect(flags & IsPlanarGraph.IPG_DUPLICATED_EDGES)
            .toBe(IsPlanarGraph.IPG_DUPLICATED_EDGES);
        // Groups sorted by ordered edge key: (0,1) then (1,2).
        expect(ipg.getDuplicatedEdges()).toEqual([[0, 2], [1, 3]]);
    });

    it('flags degenerate edges', () => {
        const positions: P2[] = [[0, 0], [1, 0]];
        const edges: E2[] = [[0, 1], [1, 1]];
        const ipg = new IsPlanarGraph();
        const flags = ipg.compute(positions, edges);
        expect(flags & IsPlanarGraph.IPG_DEGENERATE_EDGES)
            .toBe(IsPlanarGraph.IPG_DEGENERATE_EDGES);
        expect(ipg.getDegenerateEdges()).toEqual([1]);
    });

    it('flags edges with invalid vertices', () => {
        const positions: P2[] = [[0, 0], [1, 0]];
        const edges: E2[] = [[0, 1], [0, 5], [-1, 1]];
        const ipg = new IsPlanarGraph();
        const flags = ipg.compute(positions, edges);
        expect(flags & IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES)
            .toBe(IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES);
        expect(ipg.getEdgesWithInvalidVertices()).toEqual([1, 2]);
    });

    it('combines multiple failure flags', () => {
        const positions: P2[] = [[0, 0], [1, 1], [0, 1], [1, 0], [1, 1]];
        const edges: E2[] = [[0, 1], [2, 3], [1, 0], [4, 4]];
        const ipg = new IsPlanarGraph();
        const flags = ipg.compute(positions, edges);
        expect(flags).toBe(
            IsPlanarGraph.IPG_DUPLICATED_POSITIONS
            | IsPlanarGraph.IPG_DUPLICATED_EDGES
            | IsPlanarGraph.IPG_DEGENERATE_EDGES
            | IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        // Edge 2 duplicates edge 0 (collinear full overlap, pair (0,2)),
        // and both cross edge 1; the reported pairs are sorted
        // lexicographically.
        expect(ipg.getInvalidIntersections().map((e) => e.v))
            .toEqual([[0, 1], [0, 2], [1, 2]]);
    });

    it('clears state between calls', () => {
        const ipg = new IsPlanarGraph();
        const crossing: P2[] = [[0, 0], [1, 1], [0, 1], [1, 0]];
        expect(ipg.compute(crossing, [[0, 1], [2, 3]]))
            .toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        expect(ipg.getInvalidIntersections().length).toBe(1);

        const triangle: P2[] = [[0, 0], [1, 0], [0, 1]];
        expect(ipg.compute(triangle, [[0, 1], [1, 2], [2, 0]]))
            .toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
        expect(ipg.getInvalidIntersections()).toEqual([]);
    });

    it('handles a larger planar grid graph', () => {
        // 5x5 grid of vertices with horizontal and vertical edges: planar.
        const positions: P2[] = [];
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                positions.push([x, y]);
            }
        }
        const edges: E2[] = [];
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                const i = 5 * y + x;
                if (x < 4) {
                    edges.push([i, i + 1]);
                }
                if (y < 4) {
                    edges.push([i, i + 5]);
                }
            }
        }
        const ipg = new IsPlanarGraph();
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);

        // Adding one long diagonal that crosses interior edges makes it
        // nonplanar as drawn.
        positions.push([0.5, 0.25]);
        positions.push([1.5, 0.75]);
        edges.push([25, 26]);
        expect(ipg.compute(positions, edges)).toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        expect(ipg.getInvalidIntersections().length).toBeGreaterThan(0);
    });
});
