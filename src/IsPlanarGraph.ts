// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IsPlanarGraph.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Test whether an undirected graph is planar. The input positions must be
// unique and the input edges must be unique, so the number of positions is
// at least 2 and the number of edges is at least one. The elements of the
// edges array must be indices in {0,..,positions.length-1}.
//
// A sort-and-sweep algorithm is used to determine edge-edge intersections.
// If none of the intersections occur at edge-interior points, the graph is
// planar. See Game Physics (2nd edition), Section 6.2.2: Culling with
// Axis-Aligned Bounding Boxes for such an algorithm. The compute function
// returns a combination of IPG_* flags; the combined value 0 indicates the
// graph is planar. If it is not zero, the invalid-intersections list
// contains pairs of edges that intersect at an edge-interior point (that by
// definition is not a graph vertex). Each element is a pair of indices into
// the 'edges' array; the indices are ordered so that element.v[0] <
// element.v[1]. The Real type must be chosen to guarantee exact computation
// of edge-edge intersections; the port computes with IEEE double precision.
//
// Port notes: upstream operator() becomes compute(...). The upstream
// std::map/std::set containers iterate in sorted key order, which
// determines the order of the reported duplicated positions, duplicated
// edges and invalid intersections; the port replicates that order with
// explicit sorting. The upstream mZero/mOne members are replaced by number
// literals.

// A pair of indices (v0,v1) into the positions array is stored as
// (min(v0,v1), max(v0,v1)). This supports sorted containers of edges.
export class OrderedEdge {
    v: [number, number];

    constructor(v0: number = -1, v1: number = -1) {
        if (v0 < v1) {
            // v0 is minimum
            this.v = [v0, v1];
        } else {
            // v1 is minimum
            this.v = [v1, v0];
        }
    }

    // Lexicographical ordering, the port of upstream operator<.
    lessThan(edge: OrderedEdge): boolean {
        if (this.v[0] !== edge.v[0]) {
            return this.v[0] < edge.v[0];
        }
        return this.v[1] < edge.v[1];
    }
}

// The port of the private upstream Endpoint class used by the sweep.
class IPGEndpoint {
    value: number = 0;  // endpoint value
    type: number = 0;   // '0' if interval min, '1' if interval max.
    index: number = 0;  // index of interval containing this endpoint
}

export class IsPlanarGraph {
    static readonly IPG_IS_PLANAR_GRAPH = 0;
    static readonly IPG_INVALID_INPUT_SIZES = 1;
    static readonly IPG_DUPLICATED_POSITIONS = 2;
    static readonly IPG_DUPLICATED_EDGES = 4;
    static readonly IPG_DEGENERATE_EDGES = 8;
    static readonly IPG_EDGES_WITH_INVALID_VERTICES = 16;
    static readonly IPG_INVALID_INTERSECTIONS = 32;

    private mDuplicatedPositions: number[][] = [];
    private mDuplicatedEdges: number[][] = [];
    private mDegenerateEdges: number[] = [];
    private mEdgesWithInvalidVertices: number[] = [];
    private mInvalidIntersections: OrderedEdge[] = [];

    // The function returns a combination of the IPG_* flags listed above. A
    // combined value of 0 indicates the input forms a planar graph. If the
    // combined value is not zero, you may examine the flags for the failure
    // conditions and use the get* member accessors to obtain specific
    // information about the failure. If positions.length < 2 or
    // edges.length === 0, the IPG_INVALID_INPUT_SIZES flag is set.
    compute(positions: [number, number][], edges: [number, number][]): number {
        this.mDuplicatedPositions = [];
        this.mDuplicatedEdges = [];
        this.mDegenerateEdges = [];
        this.mEdgesWithInvalidVertices = [];
        this.mInvalidIntersections = [];

        let flags = this.isValidTopology(positions, edges);
        if (flags === IsPlanarGraph.IPG_INVALID_INPUT_SIZES) {
            return flags;
        }

        // Upstream proceeds to the intersection tests even when an edge
        // references a nonexistent vertex, which reads out of bounds
        // (undefined behavior in C++). The port returns the accumulated
        // flags instead because the sweep cannot be evaluated with invalid
        // vertex indices.
        if ((flags & IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES) !== 0) {
            return flags;
        }

        const overlappingRectangles = this.computeOverlappingRectangles(positions, edges);
        for (const key of overlappingRectangles) {
            // Get the endpoints of the line segments for the edges whose
            // bounding rectangles overlapped. Determine whether the line
            // segments intersect. If they do, determine how they intersect.
            const e0 = edges[key.v[0]];
            const e1 = edges[key.v[1]];
            const p0 = positions[e0[0]];
            const p1 = positions[e0[1]];
            const q0 = positions[e1[0]];
            const q1 = positions[e1[1]];
            if (this.invalidSegmentIntersection(p0, p1, q0, q1)) {
                this.mInvalidIntersections.push(key);
            }
        }

        if (this.mInvalidIntersections.length > 0) {
            flags |= IsPlanarGraph.IPG_INVALID_INTERSECTIONS;
        }

        return flags;
    }

    // Each element lists the indices of positions that share a common value.
    getDuplicatedPositions(): readonly number[][] {
        return this.mDuplicatedPositions;
    }

    // Each element lists the indices of edges that share a common
    // (unordered) vertex pair.
    getDuplicatedEdges(): readonly number[][] {
        return this.mDuplicatedEdges;
    }

    // The indices of edges of the form (v,v).
    getDegenerateEdges(): readonly number[] {
        return this.mDegenerateEdges;
    }

    // The indices of edges that reference nonexistent vertices.
    getEdgesWithInvalidVertices(): readonly number[] {
        return this.mEdgesWithInvalidVertices;
    }

    // Pairs of indices into the 'edges' array for edges that intersect at
    // edge-interior points.
    getInvalidIntersections(): readonly OrderedEdge[] {
        return this.mInvalidIntersections;
    }

    private isValidTopology(positions: [number, number][],
        edges: [number, number][]): number {
        const numPositions = positions.length;
        const numEdges = edges.length;
        if (numPositions < 2 || numEdges === 0) {
            // The graph must have at least one edge.
            return IsPlanarGraph.IPG_INVALID_INPUT_SIZES;
        }

        // The positions must be unique. Upstream uses std::map keyed by the
        // position; the port groups by a canonical string key and sorts the
        // groups lexicographically by position to match the map iteration
        // order. Normalizing -0 to 0 matches the C++ comparison -0 == 0.
        let flags = IsPlanarGraph.IPG_IS_PLANAR_GRAPH;
        const uniquePositions = new Map<string, { p: [number, number], indices: number[] }>();
        for (let i = 0; i < numPositions; ++i) {
            const p: [number, number] = [
                positions[i][0] === 0 ? 0 : positions[i][0],
                positions[i][1] === 0 ? 0 : positions[i][1]
            ];
            const key = p[0] + ' ' + p[1];
            const found = uniquePositions.get(key);
            if (found === undefined) {
                uniquePositions.set(key, { p: p, indices: [i] });
            } else {
                found.indices.push(i);
            }
        }
        if (uniquePositions.size < numPositions) {
            // At least two positions are duplicated.
            const groups = Array.from(uniquePositions.values())
                .filter((g) => g.indices.length > 1)
                .sort((g0, g1) => (g0.p[0] !== g1.p[0] ? g0.p[0] - g1.p[0] : g0.p[1] - g1.p[1]));
            for (const group of groups) {
                this.mDuplicatedPositions.push(group.indices);
            }
            flags |= IsPlanarGraph.IPG_DUPLICATED_POSITIONS;
        }

        // The edges must be unique. Upstream uses std::map keyed by the
        // OrderedEdge; the port groups by key string and sorts the groups.
        const uniqueEdges = new Map<string, { e: OrderedEdge, indices: number[] }>();
        for (let i = 0; i < numEdges; ++i) {
            const e = new OrderedEdge(edges[i][0], edges[i][1]);
            const key = e.v[0] + ' ' + e.v[1];
            const found = uniqueEdges.get(key);
            if (found === undefined) {
                uniqueEdges.set(key, { e: e, indices: [i] });
            } else {
                found.indices.push(i);
            }
        }
        if (uniqueEdges.size < numEdges) {
            // At least two edges are duplicated, possibly even a pair of
            // edges (v0,v1) and (v1,v0) which is not allowed because the
            // graph is undirected.
            const groups = Array.from(uniqueEdges.values())
                .filter((g) => g.indices.length > 1)
                .sort((g0, g1) => (g0.e.lessThan(g1.e) ? -1 : (g1.e.lessThan(g0.e) ? 1 : 0)));
            for (const group of groups) {
                this.mDuplicatedEdges.push(group.indices);
            }
            flags |= IsPlanarGraph.IPG_DUPLICATED_EDGES;
        }

        // The edges are represented as pairs of indices into the 'positions'
        // array. The indices for a single edge must be different (no edges
        // allowed from a vertex to itself) and all indices must be valid.
        for (let i = 0; i < numEdges; ++i) {
            const e = edges[i];
            if (e[0] === e[1]) {
                // The edge is degenerate, originating and terminating at the
                // same vertex.
                this.mDegenerateEdges.push(i);
                flags |= IsPlanarGraph.IPG_DEGENERATE_EDGES;
            }

            if (e[0] < 0 || e[0] >= numPositions || e[1] < 0 || e[1] >= numPositions) {
                // The edge has an index that references a nonexistent
                // vertex.
                this.mEdgesWithInvalidVertices.push(i);
                flags |= IsPlanarGraph.IPG_EDGES_WITH_INVALID_VERTICES;
            }
        }

        return flags;
    }

    // Returns the pairs of edge indices whose bounding rectangles overlap,
    // sorted lexicographically (the upstream std::set iteration order).
    private computeOverlappingRectangles(positions: [number, number][],
        edges: [number, number][]): OrderedEdge[] {
        // Compute axis-aligned bounding rectangles for the edges.
        const numEdges = edges.length;
        const emin: [number, number][] = new Array(numEdges);
        const emax: [number, number][] = new Array(numEdges);
        for (let i = 0; i < numEdges; ++i) {
            const e = edges[i];
            const p0 = positions[e[0]];
            const p1 = positions[e[1]];

            emin[i] = [0, 0];
            emax[i] = [0, 0];
            for (let j = 0; j < 2; ++j) {
                if (p0[j] < p1[j]) {
                    emin[i][j] = p0[j];
                    emax[i][j] = p1[j];
                } else {
                    emin[i][j] = p1[j];
                    emax[i][j] = p0[j];
                }
            }
        }

        // Get the rectangle endpoints. Only the x-endpoints participate in
        // the sweep; the y-overlap test reads emin/emax directly, as
        // upstream does.
        const numEndpoints = 2 * numEdges;
        const xEndpoints: IPGEndpoint[] = new Array(numEndpoints);
        for (let i = 0, j = 0; i < numEdges; ++i) {
            let endpoint = new IPGEndpoint();
            endpoint.type = 0;
            endpoint.value = emin[i][0];
            endpoint.index = i;
            xEndpoints[j++] = endpoint;

            endpoint = new IPGEndpoint();
            endpoint.type = 1;
            endpoint.value = emax[i][0];
            endpoint.index = i;
            xEndpoints[j++] = endpoint;
        }

        // Sort the rectangle endpoints by (value, type).
        xEndpoints.sort((a, b) => (a.value !== b.value ? a.value - b.value : a.type - b.type));

        // Sweep through the endpoints to determine overlapping x-intervals.
        // Use an active set of rectangles to reduce the complexity of the
        // algorithm.
        const overlapping = new Map<string, OrderedEdge>();
        const active = new Set<number>();
        for (let i = 0; i < numEndpoints; ++i) {
            const endpoint = xEndpoints[i];
            const index = endpoint.index;
            if (endpoint.type === 0) {  // an interval 'begin' value
                // In the 1D problem, the current interval overlaps with all
                // the active intervals. In 2D we also need to check for
                // y-overlap.
                for (const activeIndex of active) {
                    // Rectangles activeIndex and index overlap in the
                    // x-dimension. Test for overlap in the y-dimension.
                    const r0min = emin[activeIndex];
                    const r0max = emax[activeIndex];
                    const r1min = emin[index];
                    const r1max = emax[index];
                    if (r0max[1] >= r1min[1] && r0min[1] <= r1max[1]) {
                        const edge = activeIndex < index
                            ? new OrderedEdge(activeIndex, index)
                            : new OrderedEdge(index, activeIndex);
                        overlapping.set(edge.v[0] + ' ' + edge.v[1], edge);
                    }
                }
                active.add(index);
            } else {  // an interval 'end' value
                active.delete(index);
            }
        }

        return Array.from(overlapping.values())
            .sort((e0, e1) => (e0.lessThan(e1) ? -1 : (e1.lessThan(e0) ? 1 : 0)));
    }

    private invalidSegmentIntersection(
        p0: [number, number], p1: [number, number],
        q0: [number, number], q1: [number, number]): boolean {
        // We must solve the two linear equations
        //   p0 + t0 * (p1 - p0) = q0 + t1 * (q1 - q0)
        // for the unknown variables t0 and t1. These may be written as
        //   t0 * (p1 - p0) - t1 * (q1 - q0) = q0 - p0
        // If denom = Dot(p1 - p0, Perp(q1 - q0)) is not zero, then
        //   t0 = Dot(q0 - p0, Perp(q1 - q0)) / denom = numer0 / denom
        //   t1 = Dot(q0 - p0, Perp(p1 - p0)) / denom = numer1 / denom
        // For an invalid intersection, we need (t0,t1) with:
        // ((0 < t0 < 1) and (0 <= t1 <= 1)) or ((0 <= t0 <= 1) and
        // (0 < t1 < 1)).

        const p1mp0: [number, number] = [p1[0] - p0[0], p1[1] - p0[1]];
        const q1mq0: [number, number] = [q1[0] - q0[0], q1[1] - q0[1]];
        const q0mp0: [number, number] = [q0[0] - p0[0], q0[1] - p0[1]];

        const denom = p1mp0[0] * q1mq0[1] - p1mp0[1] * q1mq0[0];
        const numer0 = q0mp0[0] * q1mq0[1] - q0mp0[1] * q1mq0[0];
        const numer1 = q0mp0[0] * p1mp0[1] - q0mp0[1] * p1mp0[0];

        if (denom !== 0) {
            // The lines of the segments are not parallel.
            if (denom > 0) {
                if (0 <= numer0 && numer0 <= denom && 0 <= numer1 && numer1 <= denom) {
                    // The segments intersect.
                    return (numer0 !== 0 && numer0 !== denom) || (numer1 !== 0 && numer1 !== denom);
                } else {
                    return false;
                }
            } else {  // denom < 0
                if (0 >= numer0 && numer0 >= denom && 0 >= numer1 && numer1 >= denom) {
                    // The segments intersect.
                    return (numer0 !== 0 && numer0 !== denom) || (numer1 !== 0 && numer1 !== denom);
                } else {
                    return false;
                }
            }
        } else {
            // The lines of the segments are parallel.
            if (numer0 !== 0 || numer1 !== 0) {
                // The lines of the segments are separated.
                return false;
            } else {
                // The segments lie on the same line. Compute the parameter
                // intervals for the segments in terms of the t0-parameter
                // and determine their overlap (if any).
                const q1mp0: [number, number] = [q1[0] - p0[0], q1[1] - p0[1]];
                const sqrLenP1mP0 = p1mp0[0] * p1mp0[0] + p1mp0[1] * p1mp0[1];
                const value0 = q0mp0[0] * p1mp0[0] + q0mp0[1] * p1mp0[1];
                const value1 = q1mp0[0] * p1mp0[0] + q1mp0[1] * p1mp0[1];
                if ((value0 >= sqrLenP1mP0 && value1 >= sqrLenP1mP0)
                    || (value0 <= 0 && value1 <= 0)) {
                    // If the segments intersect, they must do so at
                    // endpoints of the segments.
                    return false;
                } else {
                    // The segments overlap in a positive-length interval.
                    return true;
                }
            }
        }
    }
}
