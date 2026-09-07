import { describe, it, expect } from 'vitest';
import { IsPlanarGraph, OrderedEdge } from '../src/IsPlanarGraph.js';
import { check, fc } from './helpers/arbitraries.js';

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

describe('IsPlanarGraph verification', () => {
    type P2 = [number, number];
    type E2 = [number, number];

    // Integer lattice coordinates keep every product in the predicates exact
    // in IEEE double precision, which is what upstream requires of Real.
    const orient = (a: P2, b: P2, c: P2): number =>
        (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);

    const onSegment = (a: P2, b: P2, c: P2): boolean =>
        orient(a, b, c) === 0 &&
        Math.min(a[0], b[0]) <= c[0] && c[0] <= Math.max(a[0], b[0]) &&
        Math.min(a[1], b[1]) <= c[1] && c[1] <= Math.max(a[1], b[1]);

    const strictlyInside = (a: P2, b: P2, c: P2): boolean =>
        onSegment(a, b, c) &&
        !(c[0] === a[0] && c[1] === a[1]) && !(c[0] === b[0] && c[1] === b[1]);

    // An independent predicate for "the two segments share a point that is
    // interior to at least one of them", written with orientation tests
    // rather than the upstream parametric formulation. Both segments must be
    // nondegenerate; the degenerate case is asymmetric upstream and is
    // covered by its own test below.
    const referenceInvalid = (p0: P2, p1: P2, q0: P2, q1: P2): boolean => {
        const o1 = orient(p0, p1, q0);
        const o2 = orient(p0, p1, q1);
        const o3 = orient(q0, q1, p0);
        const o4 = orient(q0, q1, p1);
        if (o1 !== 0 && o2 !== 0 && o3 !== 0 && o4 !== 0) {
            // Interior-interior crossing.
            return (o1 > 0) !== (o2 > 0) && (o3 > 0) !== (o4 > 0);
        }
        if (strictlyInside(p0, p1, q0) || strictlyInside(p0, p1, q1) ||
            strictlyInside(q0, q1, p0) || strictlyInside(q0, q1, p1)) {
            return true;
        }
        // Collinear and equal as point sets (each endpoint of one coincides
        // with an endpoint of the other), so they share interior points.
        return o1 === 0 && o2 === 0 && o3 === 0 && o4 === 0 &&
            onSegment(p0, p1, q0) && onSegment(p0, p1, q1);
    };

    // A graph on a small integer lattice. Positions are unique and edges are
    // unique and nondegenerate with valid indices, so the only flag that can
    // be raised is IPG_INVALID_INTERSECTIONS.
    const latticeGraph = fc.tuple(
        fc.uniqueArray(fc.tuple(fc.integer({ min: 0, max: 4 }),
            fc.integer({ min: 0, max: 4 })),
        { minLength: 2, maxLength: 9, selector: (p) => p[0] + ',' + p[1] }),
        fc.array(fc.tuple(fc.nat(), fc.nat()), { minLength: 1, maxLength: 12 }))
        .map(([points, raw]) => {
            const positions = points.map((p) => [p[0], p[1]] as P2);
            const edges: E2[] = [];
            const seen = new Set<string>();
            for (const [a, b] of raw) {
                const i = a % positions.length;
                const j = b % positions.length;
                if (i === j) {
                    continue;
                }
                const key = Math.min(i, j) + ',' + Math.max(i, j);
                if (seen.has(key)) {
                    continue;
                }
                seen.add(key);
                edges.push([i, j]);
            }
            return { positions, edges };
        })
        .filter(({ edges }) => edges.length > 0);

    it('the sweep finds exactly the brute-force invalid intersections', () => {
        check(latticeGraph, ({ positions, edges }) => {
            const query = new IsPlanarGraph();
            const flags = query.compute(positions, edges);

            const expected: string[] = [];
            for (let i = 0; i < edges.length; ++i) {
                for (let j = i + 1; j < edges.length; ++j) {
                    if (referenceInvalid(
                        positions[edges[i][0]], positions[edges[i][1]],
                        positions[edges[j][0]], positions[edges[j][1]])) {
                        expected.push(i + ',' + j);
                    }
                }
            }
            const actual = query.getInvalidIntersections()
                .map((e) => e.v[0] + ',' + e.v[1]);
            expect(actual).toEqual(expected);
            expect(flags).toBe(expected.length === 0
                ? IsPlanarGraph.IPG_IS_PLANAR_GRAPH
                : IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        });
    });

    it('reported intersection pairs are sorted and have v[0] < v[1]', () => {
        check(latticeGraph, ({ positions, edges }) => {
            const query = new IsPlanarGraph();
            query.compute(positions, edges);
            const pairs = query.getInvalidIntersections();
            for (const pair of pairs) {
                expect(pair.v[0]).toBeLessThan(pair.v[1]);
            }
            for (let i = 0; i + 1 < pairs.length; ++i) {
                expect(pairs[i].lessThan(pairs[i + 1])).toBe(true);
            }
        });
    });

    it('translating and scaling the embedding preserves the verdict', () => {
        // The predicates are exact on the lattice, and an integer similarity
        // transform maps the lattice onto itself, so the flags must agree.
        check(fc.tuple(latticeGraph, fc.integer({ min: -6, max: 6 }),
            fc.integer({ min: -6, max: 6 }), fc.integer({ min: 1, max: 4 })),
        ([{ positions, edges }, dx, dy, s]) => {
            const first = new IsPlanarGraph();
            const flags0 = first.compute(positions, edges);
            const moved = positions.map((p) =>
                [s * p[0] + dx, s * p[1] + dy] as P2);
            const second = new IsPlanarGraph();
            const flags1 = second.compute(moved, edges);
            expect(flags1).toBe(flags0);
            expect(second.getInvalidIntersections().map((e) => e.v))
                .toEqual(first.getInvalidIntersections().map((e) => e.v));
        });
    });

    it('reflecting the embedding preserves the verdict', () => {
        check(latticeGraph, ({ positions, edges }) => {
            const first = new IsPlanarGraph();
            const flags0 = first.compute(positions, edges);
            // Swap x and y: an orientation-reversing lattice isometry, so the
            // signs of every orientation test flip.
            const flipped = positions.map((p) => [p[1], p[0]] as P2);
            const second = new IsPlanarGraph();
            expect(second.compute(flipped, edges)).toBe(flags0);
            expect(second.getInvalidIntersections().map((e) => e.v))
                .toEqual(first.getInvalidIntersections().map((e) => e.v));
        });
    });

    it('reversing an edge direction does not change the verdict', () => {
        check(fc.tuple(latticeGraph, fc.nat()),
            ([{ positions, edges }, which]) => {
                const first = new IsPlanarGraph();
                const flags0 = first.compute(positions, edges);
                const k = which % edges.length;
                const reversed = edges.map((e, i) =>
                    (i === k ? [e[1], e[0]] as E2 : e));
                const second = new IsPlanarGraph();
                expect(second.compute(positions, reversed)).toBe(flags0);
                expect(second.getInvalidIntersections().map((e) => e.v))
                    .toEqual(first.getInvalidIntersections().map((e) => e.v));
            });
    });

    // Points on a parabola are in strictly convex position and have exact
    // integer coordinates, so every predicate in the query is exact. (Points
    // on a circle are NOT usable: with double coordinates a chord and a
    // boundary edge that share a vertex give t-parameters a few ulps away
    // from the endpoint value, which upstream's exact-arithmetic contract
    // rules out but double precision cannot deliver.)
    const convexPolygon = (n: number): P2[] => {
        const positions: P2[] = [];
        for (let i = 0; i < n; ++i) {
            positions.push([i, i * i]);
        }
        return positions;
    };

    // Two chords of a convex polygon whose vertices are cyclically ordered
    // 0..n-1 cross at an interior point exactly when they interleave.
    const interleave = (e0: E2, e1: E2): boolean => {
        const [a, b] = [Math.min(e0[0], e0[1]), Math.max(e0[0], e0[1])];
        const [c, d] = [Math.min(e1[0], e1[1]), Math.max(e1[0], e1[1])];
        return (a < c && c < b && b < d) || (c < a && a < d && d < b);
    };

    it('a triangulated convex polygon fan is always planar', () => {
        check(fc.integer({ min: 3, max: 12 }), (n) => {
            const positions = convexPolygon(n);
            const edges: E2[] = [];
            for (let i = 0; i < n; ++i) {
                edges.push([i, (i + 1) % n]);
            }
            for (let i = 2; i < n - 1; ++i) {
                edges.push([0, i]);
            }
            const query = new IsPlanarGraph();
            expect(query.compute(positions, edges))
                .toBe(IsPlanarGraph.IPG_IS_PLANAR_GRAPH);
        });
    });

    it('the complete graph on a convex polygon reports every crossing', () => {
        // K5 and every larger complete graph is nonplanar, so a
        // straight-line drawing must have crossings; on a convex polygon the
        // crossing pairs are exactly the interleaving chord pairs.
        check(fc.integer({ min: 4, max: 9 }), (n) => {
            const positions = convexPolygon(n);
            const edges: E2[] = [];
            for (let i = 0; i < n; ++i) {
                for (let j = i + 1; j < n; ++j) {
                    edges.push([i, j]);
                }
            }
            let crossings = 0;
            for (let i = 0; i < edges.length; ++i) {
                for (let j = i + 1; j < edges.length; ++j) {
                    if (interleave(edges[i], edges[j])) {
                        ++crossings;
                    }
                }
            }
            // The number of crossings of a convex drawing of K_n is C(n,4).
            expect(crossings).toBe(n * (n - 1) * (n - 2) * (n - 3) / 24);

            const query = new IsPlanarGraph();
            expect(query.compute(positions, edges))
                .toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
            expect(query.getInvalidIntersections().length).toBe(crossings);
        });
    });

    it('a straight-line K3,3 drawn as two rows always has crossings', () => {
        const positions: P2[] = [
            [0, 1], [1, 1], [2, 1],
            [0, 0], [1, 0], [2, 0]
        ];
        const edges: E2[] = [];
        for (let i = 0; i < 3; ++i) {
            for (let j = 3; j < 6; ++j) {
                edges.push([i, j]);
            }
        }
        const query = new IsPlanarGraph();
        expect(query.compute(positions, edges) &
            IsPlanarGraph.IPG_INVALID_INTERSECTIONS).not.toBe(0);
    });

    it('the topology flags are independent of one another', () => {
        // positions[2] duplicates positions[0]; edge 1 duplicates edge 0
        // reversed; edge 2 is degenerate; edge 3 has an invalid vertex.
        const positions: P2[] = [[0, 0], [1, 0], [0, 0]];
        const cases: { edges: E2[], flags: number }[] = [
            { edges: [[0, 1]], flags: IsPlanarGraph.IPG_DUPLICATED_POSITIONS },
            {
                // The duplicated edge is also a self-overlapping segment.
                edges: [[0, 1], [1, 0]],
                flags: IsPlanarGraph.IPG_DUPLICATED_POSITIONS |
                    IsPlanarGraph.IPG_DUPLICATED_EDGES |
                    IsPlanarGraph.IPG_INVALID_INTERSECTIONS
            },
            {
                edges: [[0, 1], [1, 1]],
                flags: IsPlanarGraph.IPG_DUPLICATED_POSITIONS |
                    IsPlanarGraph.IPG_DEGENERATE_EDGES
            },
            {
                edges: [[0, 1], [1, 7]],
                flags: IsPlanarGraph.IPG_DUPLICATED_POSITIONS |
                    IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES
            }
        ];
        for (const { edges, flags } of cases) {
            const query = new IsPlanarGraph();
            expect(query.compute(positions, edges)).toBe(flags);
        }
    });

    it('an out-of-range vertex index returns early instead of reading out of bounds', () => {
        // Upstream indexes positions[] with the invalid index and reads out
        // of bounds; the port returns the accumulated flags first.
        const positions: P2[] = [[0, 0], [1, 0], [0, 1], [1, 1]];
        const edges: E2[] = [[0, 3], [1, 2], [0, 99]];
        const query = new IsPlanarGraph();
        const flags = query.compute(positions, edges);
        expect(flags).toBe(IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES);
        expect(query.getEdgesWithInvalidVertices()).toEqual([2]);
        expect(query.getInvalidIntersections()).toEqual([]);
        // Without the invalid edge, the two diagonals do cross.
        const ok = new IsPlanarGraph();
        expect(ok.compute(positions, [[0, 3], [1, 2]]))
            .toBe(IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
    });

    it('degenerate edges are asymmetric upstream (quirk preserved)', () => {
        // A degenerate edge that lies in the interior of another edge is
        // reported as an invalid intersection only when the degenerate edge
        // has the larger index: upstream's collinear branch parameterizes by
        // the FIRST segment, and a zero-length first segment collapses the
        // overlap test to "no interior overlap".
        const positions: P2[] = [[0, 0], [2, 0], [1, 0]];
        const spanFirst = new IsPlanarGraph();
        expect(spanFirst.compute(positions, [[0, 1], [2, 2]])).toBe(
            IsPlanarGraph.IPG_DEGENERATE_EDGES |
            IsPlanarGraph.IPG_INVALID_INTERSECTIONS);
        const pointFirst = new IsPlanarGraph();
        expect(pointFirst.compute(positions, [[2, 2], [0, 1]]))
            .toBe(IsPlanarGraph.IPG_DEGENERATE_EDGES);
    });

    it('repeated compute calls on one object are independent', () => {
        check(fc.tuple(latticeGraph, latticeGraph), ([g0, g1]) => {
            const shared = new IsPlanarGraph();
            const flags0 = shared.compute(g0.positions, g0.edges);
            const pairs0 = shared.getInvalidIntersections().map((e) => e.v);
            shared.compute(g1.positions, g1.edges);
            expect(shared.compute(g0.positions, g0.edges)).toBe(flags0);
            expect(shared.getInvalidIntersections().map((e) => e.v))
                .toEqual(pairs0);
        });
    });
});
