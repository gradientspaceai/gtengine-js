// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BSPPolygon2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Boolean operations (intersection, union, difference, exclusive or) on
// planar polygons, implemented with a binary space partitioning (BSP) tree.
// A polygon is a collection of vertices and directed edges; the convention is
// that the interior of the polygon is on the negative side of its directed
// edges, which for a simple polygon means the vertices are ordered
// counterclockwise.
//
// Typical use: insert the vertices and edges of the polygon, call finalize()
// to build the BSP tree, and then apply the Boolean operations. The result of
// a Boolean operation has its tree already built.
//
// Port notes:
// - Upstream 'template <typename Real>' becomes 'number'; Vector2<Real>
//   becomes the runtime-sized Vector with size 2.
// - The nested private class BSPTree2 becomes a module-private class of the
//   same name (nothing is exported for it, so the global export uniqueness
//   rule is respected). The unnamed enum of edge classifications becomes the
//   module-private const object Classification.
// - std::map<Vertex,int32_t> and std::map<Edge,int32_t> become Map objects
//   keyed by a string encoding of the exact vertex coordinates or edge
//   indices; lookups upstream are exact (not epsilon based), so the encoding
//   is faithful. Where upstream relies on the ordered iteration of the edge
//   map (the negation operator), the port sorts the edges lexicographically
//   by (V[0], V[1]) to reproduce the std::map order exactly.
// - The operators map to methods, per PORTING.md: operator~ -> negated(),
//   operator& -> intersection(), operator| -> union(), operator- ->
//   difference(), operator^ -> exclusiveOr(). Copy semantics become clone().
// - The debug-only Print methods (guarded by
//   GTE_BSPPOLYGON2_ENABLE_DEBUG_PRINT) write files and are omitted.
// - splitEdge is private upstream but is called by the nested BSPTree2; since
//   the port's BSPTree2 is a separate module-level class, splitEdge is a
//   public method documented as internal.

import { EdgeKey } from './EdgeKey.js';
import { logAssert } from './Logger.js';
import { Vector, add, dot, mul, sub } from './Vector.js';
import { perp } from './Vector2.js';

// The upstream unnamed enum inside BSPTree2.
const Classification = {
    TRANSVERSE_POSITIVE: 0,
    TRANSVERSE_NEGATIVE: 1,
    ALL_POSITIVE: 2,
    ALL_NEGATIVE: 3,
    COINCIDENT: 4
} as const;

// The key of a vertex in the vertex map. The exact coordinates are the key,
// as they are for std::map<Vector2<Real>, int32_t>; the negative zero is
// normalized to zero because the C++ comparison treats -0 and +0 as equal.
function vertexKey(vertex: Vector): string {
    const x = (vertex.values[0] === 0 ? 0 : vertex.values[0]);
    const y = (vertex.values[1] === 0 ? 0 : vertex.values[1]);
    return `${x},${y}`;
}

// The key of a directed edge in the edge map.
function edgeKey(edge: EdgeKey): string {
    return `${edge.V[0]},${edge.V[1]}`;
}

// A directed edge, the port of 'typedef EdgeKey<true> Edge'.
function makeEdge(v0: number, v1: number): EdgeKey {
    return new EdgeKey(true, v0, v1);
}

// The port of the BSPTree2 class nested privately in BSPPolygon2.
class BSPTree2 {
    private mEpsilon: number;
    private mCoincident: EdgeKey[];
    private mPosChild: BSPTree2 | null;
    private mNegChild: BSPTree2 | null;

    // The empty tree used by getCopy().
    constructor(epsilon: number) {
        this.mEpsilon = (epsilon >= 0 ? epsilon : 0);
        this.mCoincident = [];
        this.mPosChild = null;
        this.mNegChild = null;
    }

    // Build the tree for the specified edges of the polygon. The polygon is
    // modified: edges that cross a splitting line are split, which inserts
    // vertices and edges.
    static build(polygon: BSPPolygon2, edges: readonly EdgeKey[],
        epsilon: number): BSPTree2 {
        const tree = new BSPTree2(epsilon);
        logAssert(edges.length > 0, 'Invalid input.');

        // Construct splitting line from first edge.
        const end0 = polygon.getVertex(edges[0].V[0]);
        const end1 = polygon.getVertex(edges[0].V[1]);

        // Add edge to coincident list.
        tree.mCoincident.push(edges[0]);

        // Split remaining edges.
        const posArray: EdgeKey[] = [];
        const negArray: EdgeKey[] = [];
        for (let i = 1; i < edges.length; ++i) {
            const v0 = edges[i].V[0];
            const v1 = edges[i].V[1];
            const vertex0 = polygon.getVertex(v0);
            const vertex1 = polygon.getVertex(v1);

            const classified = tree.classifyEdge(end0, end1, vertex0, vertex1);
            let vmid: number;

            switch (classified.type) {
                case Classification.TRANSVERSE_POSITIVE:
                    // modify edge <V0,V1> to <V0,I>, add new edge <I,V1>
                    vmid = polygon.insertVertex(classified.intr);
                    if (vmid === v0 || vmid === v1) {
                        // The intersection point is within epsilon of an
                        // endpoint.
                        tree.mCoincident.push(edges[i]);
                    } else {
                        polygon.splitEdge(v0, v1, vmid);
                        posArray.push(makeEdge(vmid, v1));
                        negArray.push(makeEdge(v0, vmid));
                    }
                    break;
                case Classification.TRANSVERSE_NEGATIVE:
                    // modify edge <V0,V1> to <V0,I>, add new edge <I,V1>
                    vmid = polygon.insertVertex(classified.intr);
                    if (vmid === v0 || vmid === v1) {
                        // The intersection point is within epsilon of an
                        // endpoint.
                        tree.mCoincident.push(edges[i]);
                    } else {
                        polygon.splitEdge(v0, v1, vmid);
                        posArray.push(makeEdge(v0, vmid));
                        negArray.push(makeEdge(vmid, v1));
                    }
                    break;
                case Classification.ALL_POSITIVE:
                    posArray.push(edges[i]);
                    break;
                case Classification.ALL_NEGATIVE:
                    negArray.push(edges[i]);
                    break;
                default:  // COINCIDENT
                    tree.mCoincident.push(edges[i]);
                    break;
            }
        }

        if (posArray.length > 0) {
            tree.mPosChild = BSPTree2.build(polygon, posArray, tree.mEpsilon);
        }

        if (negArray.length > 0) {
            tree.mNegChild = BSPTree2.build(polygon, negArray, tree.mEpsilon);
        }

        return tree;
    }

    // A deep copy of the BSP tree (upstream GetCopy).
    getCopy(): BSPTree2 {
        const tree = new BSPTree2(this.mEpsilon);

        tree.mCoincident = this.mCoincident.map(
            edge => makeEdge(edge.V[0], edge.V[1]));

        if (this.mPosChild) {
            tree.mPosChild = this.mPosChild.getCopy();
        }

        if (this.mNegChild) {
            tree.mNegChild = this.mNegChild.getCopy();
        }

        return tree;
    }

    // Polygon Boolean operation support.
    negate(): void {
        // Reverse coincident edge directions.
        for (const edge of this.mCoincident) {
            const save = edge.V[0];
            edge.V[0] = edge.V[1];
            edge.V[1] = save;
        }

        // Swap positive and negative subtrees.
        const saveChild = this.mPosChild;
        this.mPosChild = this.mNegChild;
        this.mNegChild = saveChild;

        if (this.mPosChild) {
            this.mPosChild.negate();
        }

        if (this.mNegChild) {
            this.mNegChild.negate();
        }
    }

    getPartition(polygon: BSPPolygon2, v0: Vector, v1: Vector,
        pos: BSPPolygon2, neg: BSPPolygon2, coSame: BSPPolygon2,
        coDiff: BSPPolygon2): void {
        // Construct splitting line from first coincident edge.
        const end0 = polygon.getVertex(this.mCoincident[0].V[0]);
        const end1 = polygon.getVertex(this.mCoincident[0].V[1]);

        const classified = this.classifyEdge(end0, end1, v0, v1);
        const intr = classified.intr;

        switch (classified.type) {
            case Classification.TRANSVERSE_POSITIVE:
                this.getPosPartition(polygon, intr, v1, pos, neg, coSame, coDiff);
                this.getNegPartition(polygon, v0, intr, pos, neg, coSame, coDiff);
                break;
            case Classification.TRANSVERSE_NEGATIVE:
                this.getPosPartition(polygon, v0, intr, pos, neg, coSame, coDiff);
                this.getNegPartition(polygon, intr, v1, pos, neg, coSame, coDiff);
                break;
            case Classification.ALL_POSITIVE:
                this.getPosPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
                break;
            case Classification.ALL_NEGATIVE:
                this.getNegPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
                break;
            default:  // COINCIDENT
                this.getCoPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
                break;
        }
    }

    // Point-in-polygon support. For a polygon whose interior is on the
    // negative side of its directed edges, the return value is -1 when the
    // point is inside, 0 when it is on the polygon and +1 when it is outside.
    pointLocation(polygon: BSPPolygon2, vertex: Vector): number {
        // Construct splitting line from first coincident edge.
        const end0 = polygon.getVertex(this.mCoincident[0].V[0]);
        const end1 = polygon.getVertex(this.mCoincident[0].V[1]);

        switch (this.classifyVertex(end0, end1, vertex)) {
            case Classification.ALL_POSITIVE:
                if (this.mPosChild) {
                    return this.mPosChild.pointLocation(polygon, vertex);
                }
                return 1;
            case Classification.ALL_NEGATIVE:
                if (this.mNegChild) {
                    return this.mNegChild.pointLocation(polygon, vertex);
                }
                return -1;
            default:  // COINCIDENT
                return this.coPointLocation(polygon, vertex);
        }
    }

    // Classify the edge <v0,v1> against the line through end0 and end1. The
    // 'intr' output of the C++ function becomes a field of the returned
    // object; it is meaningful only for the transverse cases.
    private classifyEdge(end0: Vector, end1: Vector, v0: Vector, v1: Vector):
        { type: number, intr: Vector } {
        const dir = sub(end1, end0);
        const nor = perp(dir);
        const diff0 = sub(v0, end0);
        const diff1 = sub(v1, end0);

        let d0 = dot(nor, diff0);
        let d1 = dot(nor, diff1);

        if (d0 * d1 < 0) {
            // Edge <V0,V1> transversely crosses line. Compute point of
            // intersection I = V0 + t*(V1 - V0).
            const t = d0 / (d0 - d1);
            if (t > this.mEpsilon) {
                if (t < 1 - this.mEpsilon) {
                    const intr = add(v0, mul(sub(v1, v0), t));
                    if (d1 > 0) {
                        return { type: Classification.TRANSVERSE_POSITIVE, intr };
                    }
                    return { type: Classification.TRANSVERSE_NEGATIVE, intr };
                }
                // T is effectively 1 (numerical round-off issue), so set
                // d1 = 0 and go on to other cases.
                d1 = 0;
            } else {
                // T is effectively 0 (numerical round-off issue), so set
                // d0 = 0 and go on to other cases.
                d0 = 0;
            }
        }

        const intr = new Vector(2);
        if (d0 > 0 || d1 > 0) {
            // edge on positive side of line
            return { type: Classification.ALL_POSITIVE, intr };
        }

        if (d0 < 0 || d1 < 0) {
            // edge on negative side of line
            return { type: Classification.ALL_NEGATIVE, intr };
        }

        return { type: Classification.COINCIDENT, intr };
    }

    private getPosPartition(polygon: BSPPolygon2, v0: Vector, v1: Vector,
        pos: BSPPolygon2, neg: BSPPolygon2, coSame: BSPPolygon2,
        coDiff: BSPPolygon2): void {
        if (this.mPosChild) {
            this.mPosChild.getPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        } else {
            const i0 = pos.insertVertex(v0);
            const i1 = pos.insertVertex(v1);
            pos.insertEdge(makeEdge(i0, i1));
        }
    }

    private getNegPartition(polygon: BSPPolygon2, v0: Vector, v1: Vector,
        pos: BSPPolygon2, neg: BSPPolygon2, coSame: BSPPolygon2,
        coDiff: BSPPolygon2): void {
        if (this.mNegChild) {
            this.mNegChild.getPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
        } else {
            const i0 = neg.insertVertex(v0);
            const i1 = neg.insertVertex(v1);
            neg.insertEdge(makeEdge(i0, i1));
        }
    }

    private getCoPartition(polygon: BSPPolygon2, v0: Vector, v1: Vector,
        pos: BSPPolygon2, neg: BSPPolygon2, coSame: BSPPolygon2,
        coDiff: BSPPolygon2): void {
        // Segment the line containing V0 and V1 by the coincident intervals
        // that intersect <V0,V1>.
        const dir = sub(v1, v0);
        const tmax = dot(dir, dir);

        let end0: Vector;
        let end1: Vector;
        let t0: number;
        let t1: number;
        let sameDir: boolean;

        // The C++ std::list<Interval> becomes an array; the list insertions
        // become splice calls at the same positions.
        const intervals: Interval[] = [];

        for (const edge of this.mCoincident) {
            end0 = polygon.getVertex(edge.V[0]);
            end1 = polygon.getVertex(edge.V[1]);

            t0 = dot(dir, sub(end0, v0));
            if (Math.abs(t0) <= this.mEpsilon) {
                t0 = 0;
            } else if (Math.abs(t0 - tmax) <= this.mEpsilon) {
                t0 = tmax;
            }

            t1 = dot(dir, sub(end1, v0));
            if (Math.abs(t1) <= this.mEpsilon) {
                t1 = 0;
            } else if (Math.abs(t1 - tmax) <= this.mEpsilon) {
                t1 = tmax;
            }

            sameDir = (t1 > t0);
            if (!sameDir) {
                const save = t0;
                t0 = t1;
                t1 = save;
            }

            if (t1 > 0 && t0 < tmax) {
                if (intervals.length === 0) {
                    intervals.unshift(new Interval(t0, t1, sameDir, true));
                } else {
                    for (let k = 0; k < intervals.length; ++k) {
                        const iter = intervals[k];
                        if (Math.abs(t1 - iter.t0) <= this.mEpsilon) {
                            t1 = iter.t0;
                        }

                        if (t1 <= iter.t0) {
                            // [t0,t1] is on the left of [I.t0,I.t1]
                            intervals.splice(k, 0, new Interval(t0, t1, sameDir, true));
                            break;
                        }

                        // Theoretically, the intervals are disjoint or
                        // intersect only at an end point. The assert makes
                        // sure that [t0,t1] is to the right of [I.t0,I.t1].
                        if (Math.abs(t0 - iter.t1) <= this.mEpsilon) {
                            t0 = iter.t1;
                        }

                        logAssert(t0 >= iter.t1,
                            'Invalid ordering in BSPTree2::GetCoPartition.');

                        if (k === intervals.length - 1) {
                            intervals.push(new Interval(t0, t1, sameDir, true));
                            break;
                        }
                    }
                }
            }
        }

        if (intervals.length === 0) {
            this.getPosPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
            this.getNegPartition(polygon, v0, v1, pos, neg, coSame, coDiff);
            return;
        }

        // Insert outside intervals between the touching intervals. It is
        // possible that two touching intervals are adjacent, so this is not
        // just a simple alternation of touching and outside intervals.
        const front = intervals[0];
        if (front.t0 > 0) {
            intervals.unshift(new Interval(0, front.t0, front.sameDir, false));
        } else {
            front.t0 = 0;
        }

        const back = intervals[intervals.length - 1];
        if (back.t1 < tmax) {
            intervals.push(new Interval(back.t1, tmax, back.sameDir, false));
        } else {
            back.t1 = tmax;
        }

        // The C++ loop walks two adjacent list iterators and inserts a gap
        // interval before the second one; the insert returns an iterator to
        // the new element, which the loop then advances past.
        let k0 = 0;
        let k1 = 1;
        while (k1 < intervals.length) {
            const a = intervals[k0].t1;
            const b = intervals[k1].t0;
            if (b - a > this.mEpsilon) {
                intervals.splice(k1, 0, new Interval(a, b, true, false));
                k0 = k1 + 1;
                k1 = k1 + 2;
            } else {
                ++k0;
                ++k1;
            }
        }

        // Process the segmentation.
        const invTMax = 1 / tmax;
        end1 = add(v0, mul(dir, intervals[0].t0 * invTMax));
        for (const iter of intervals) {
            end0 = end1;
            end1 = add(v0, mul(dir, iter.t1 * invTMax));

            if (iter.touching) {
                if (iter.sameDir) {
                    const i0 = coSame.insertVertex(end0);
                    const i1 = coSame.insertVertex(end1);
                    if (i0 !== i1) {
                        coSame.insertEdge(makeEdge(i0, i1));
                    }
                } else {
                    const i0 = coDiff.insertVertex(end1);
                    const i1 = coDiff.insertVertex(end0);
                    if (i0 !== i1) {
                        coDiff.insertEdge(makeEdge(i0, i1));
                    }
                }
            } else {
                this.getPosPartition(polygon, end0, end1, pos, neg, coSame, coDiff);
                this.getNegPartition(polygon, end0, end1, pos, neg, coSame, coDiff);
            }
        }
    }

    // Point-in-polygon support: classify a single vertex against the line
    // through end0 and end1.
    private classifyVertex(end0: Vector, end1: Vector, vertex: Vector): number {
        const dir = sub(end1, end0);
        const nor = perp(dir);
        const diff = sub(vertex, end0);
        const c = dot(nor, diff);

        if (c > this.mEpsilon) {
            return Classification.ALL_POSITIVE;
        }

        if (c < -this.mEpsilon) {
            return Classification.ALL_NEGATIVE;
        }

        return Classification.COINCIDENT;
    }

    private coPointLocation(polygon: BSPPolygon2, vertex: Vector): number {
        for (const edge of this.mCoincident) {
            const end0 = polygon.getVertex(edge.V[0]);
            const end1 = polygon.getVertex(edge.V[1]);
            const dir = sub(end1, end0);
            const diff = sub(vertex, end0);
            const tmax = dot(dir, dir);
            const t = dot(dir, diff);

            if (-this.mEpsilon <= t && t <= tmax + this.mEpsilon) {
                return 0;
            }
        }

        // It does not matter which subtree you use.
        if (this.mPosChild) {
            return this.mPosChild.pointLocation(polygon, vertex);
        }

        if (this.mNegChild) {
            return this.mNegChild.pointLocation(polygon, vertex);
        }

        return 0;
    }
}

// A subinterval of the segment processed by BSPTree2.getCoPartition.
class Interval {
    t0: number;
    t1: number;
    sameDir: boolean;
    touching: boolean;

    constructor(inT0: number, inT1: number, inSameDir: boolean, inTouching: boolean) {
        this.t0 = inT0;
        this.t1 = inT1;
        this.sameDir = inSameDir;
        this.touching = inTouching;
    }
}

export class BSPPolygon2 {
    private mEpsilon: number;
    private mVMap: Map<string, number>;
    private mVArray: Vector[];
    private mEMap: Map<string, number>;
    private mEArray: EdgeKey[];
    private mTree: BSPTree2 | null;

    constructor(epsilon: number) {
        this.mEpsilon = (epsilon >= 0 ? epsilon : 0);
        this.mVMap = new Map<string, number>();
        this.mVArray = [];
        this.mEMap = new Map<string, number>();
        this.mEArray = [];
        this.mTree = null;
    }

    // The port of the C++ copy constructor and copy assignment. The BSP tree
    // is deep-copied.
    clone(): BSPPolygon2 {
        const result = new BSPPolygon2(this.mEpsilon);
        result.mVMap = new Map<string, number>(this.mVMap);
        result.mVArray = this.mVArray.map(vertex => vertex.clone());
        result.mEMap = new Map<string, number>(this.mEMap);
        result.mEArray = this.mEArray.map(edge => makeEdge(edge.V[0], edge.V[1]));
        result.mTree = (this.mTree ? this.mTree.getCopy() : null);
        return result;
    }

    // Support for deferred construction.
    insertVertex(vertex: Vector): number {
        logAssert(vertex.values.length === 2, 'BSPPolygon2: vertices must be 2D.');
        const key = vertexKey(vertex);
        const found = this.mVMap.get(key);
        if (found !== undefined) {
            // Vertex already in map, just return its unique index.
            return found;
        }

        // Vertex not in map, insert it and assign it a unique index.
        const i = this.mVArray.length;
        this.mVMap.set(key, i);
        this.mVArray.push(vertex.clone());
        return i;
    }

    insertEdge(edge: EdgeKey): number {
        logAssert(edge.V[0] !== edge.V[1], 'Degenerate edges not allowed.');

        const key = edgeKey(edge);
        const found = this.mEMap.get(key);
        if (found !== undefined) {
            // Edge already in map, just return its unique index.
            return found;
        }

        // Edge not in map, insert it and assign it a unique index.
        const i = this.mEArray.length;
        this.mEMap.set(key, i);
        this.mEArray.push(edge);
        return i;
    }

    // Build the BSP tree. This must be called before the Boolean operations
    // and pointLocation.
    finalize(): void {
        // Upstream note: the BSPTree2 constructor is badly designed. The
        // '*this' is passed as non-const but the previous code in this
        // function passed 'this->mEArray' as const. The mEArray can be
        // modified via '*this' in the BSPTree2 class, which led to an
        // infinite loop in a test data set for a bug report. For now, make a
        // copy of mEArray and pass it.
        const eArray = this.mEArray.map(edge => makeEdge(edge.V[0], edge.V[1]));
        this.mTree = BSPTree2.build(this, eArray, this.mEpsilon);
    }

    // Member access.
    getNumVertices(): number {
        return this.mVMap.size;
    }

    getVertex(i: number): Vector {
        return this.mVArray[i];
    }

    getNumEdges(): number {
        return this.mEMap.size;
    }

    getEdge(i: number): EdgeKey {
        return this.mEArray[i];
    }

    // Negation (upstream operator~). The vertices are shared, the edges are
    // reversed and the BSP tree is copied and negated.
    negated(): BSPPolygon2 {
        logAssert(this.mTree !== null, 'Tree must exist.');

        const neg = new BSPPolygon2(this.mEpsilon);
        neg.mVMap = new Map<string, number>(this.mVMap);
        neg.mVArray = this.mVArray.map(vertex => vertex.clone());

        // Upstream iterates the ordered std::map of edges; the port sorts the
        // edges the same way so that the reversed edges are inserted in the
        // same order.
        for (const edge of this.sortedEdges()) {
            neg.insertEdge(makeEdge(edge.V[1], edge.V[0]));
        }

        neg.mTree = (this.mTree as BSPTree2).getCopy();
        neg.mTree.negate();
        return neg;
    }

    // Intersection (upstream operator&).
    //
    // Upstream quirk, preserved: finalize() is called unconditionally, and
    // the BSP tree builder asserts that the edge list is nonempty. When the
    // two polygons do not overlap in a region of positive area (they are
    // disjoint, or they meet only along oppositely directed coincident
    // edges), the result has no edges and this throws 'Invalid input.'
    // instead of returning the empty polygon.
    intersection(polygon: BSPPolygon2): BSPPolygon2 {
        logAssert(this.mTree !== null, 'Tree must exist.');

        const intersect = new BSPPolygon2(this.mEpsilon);
        this.getInsideEdgesFrom(polygon, intersect);
        polygon.getInsideEdgesFrom(this, intersect);
        intersect.finalize();
        return intersect;
    }

    // Union (upstream operator|).
    union(polygon: BSPPolygon2): BSPPolygon2 {
        return this.negated().intersection(polygon.negated()).negated();
    }

    // Difference (upstream operator-).
    difference(polygon: BSPPolygon2): BSPPolygon2 {
        return this.intersection(polygon.negated());
    }

    // Exclusive or (upstream operator^).
    exclusiveOr(polygon: BSPPolygon2): BSPPolygon2 {
        return this.difference(polygon).union(polygon.difference(this));
    }

    // Point location: -1 when the vertex is inside the polygon, 0 when it is
    // on the polygon and +1 when it is outside. (The convention assumes the
    // polygon interior is on the negative side of its directed edges.)
    pointLocation(vertex: Vector): number {
        logAssert(this.mTree !== null, 'Tree must exist.');
        return (this.mTree as BSPTree2).pointLocation(this, vertex);
    }

    // INTERNAL (private upstream, used by BSPTree2): replace the edge
    // <v0,v1> by the two edges <v0,vmid> and <vmid,v1>.
    splitEdge(v0: number, v1: number, vmid: number): void {
        // Find the edge in the map to get the edge-array index.
        const key = edgeKey(makeEdge(v0, v1));
        const eIndex = this.mEMap.get(key);
        logAssert(eIndex !== undefined, 'Edge does not exist.');

        // Delete edge <V0,V1>.
        this.mEMap.delete(key);

        // Insert edge <V0,VM>. Upstream uses std::map::insert, which does
        // nothing when the key is already present; the port replicates that
        // (Map.set would overwrite the existing index instead).
        const index = eIndex as number;
        this.mEArray[index].V[1] = vmid;
        const splitKey = edgeKey(this.mEArray[index]);
        if (!this.mEMap.has(splitKey)) {
            this.mEMap.set(splitKey, index);
        }

        // Insert edge <VM,V1>.
        this.insertEdge(makeEdge(vmid, v1));
    }

    // The edges in the order the upstream std::map<Edge,int32_t> iterates
    // them: lexicographic by (V[0], V[1]).
    private sortedEdges(): EdgeKey[] {
        // The map values are indices into mEArray, so the map contents are
        // exactly the edges iterated by the upstream std::map.
        const inMap = Array.from(this.mEMap.values(), i => this.mEArray[i]);
        inMap.sort((e0, e1) => {
            if (e0.V[0] !== e1.V[0]) {
                return e0.V[0] - e1.V[0];
            }
            return e0.V[1] - e1.V[1];
        });
        return inMap;
    }

    private getInsideEdgesFrom(polygon: BSPPolygon2, inside: BSPPolygon2): void {
        logAssert(this.mTree !== null, 'Tree must exist.');

        const ignore = new BSPPolygon2(this.mEpsilon);
        const numEdges = polygon.getNumEdges();
        for (let i = 0; i < numEdges; ++i) {
            const v0 = polygon.mEArray[i].V[0];
            const v1 = polygon.mEArray[i].V[1];
            const vertex0 = polygon.mVArray[v0];
            const vertex1 = polygon.mVArray[v1];
            (this.mTree as BSPTree2).getPartition(this, vertex0, vertex1,
                ignore, inside, inside, ignore);
        }
    }
}
