// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TriangulateEC.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Triangulate polygons using ear clipping. The algorithm is described in
// https://www.geometrictools.com/Documentation/TriangulationByEarClipping.pdf
// The algorithm for processing nested polygons involves a division, so the
// upstream ComputeType must be rational-based, say, BSRational. If you process
// only polygons that are simple, you may use BSNumber for the ComputeType.
//
// Port notes:
// - Upstream is 'template <typename InputType, typename ComputeType>'. Per the
//   port's type mapping, both become 'number' (IEEE double). The conversion
//   step from InputType to ComputeType is therefore the identity, but the
//   ConvertPoints machinery is preserved so that a future exact-arithmetic
//   instantiation slots in unchanged. The consequence of a 'number'
//   ComputeType is that the PrimalQuery2 sign computations are exact only when
//   the inputs and the intermediate products and sums are exactly
//   representable in double precision (small integer coordinates, for
//   instance). See the port notes of PrimalQuery2.ts. This port does not offer
//   a BSNumber/BSRational variant because PrimalQuery2.ts is itself
//   number-only.
// - The four upstream operator() overloads become named methods, following
//   the IntrIntervals precedent that the canonical query keeps the short name
//   and the others get descriptive names:
//     operator()()                       -> triangulate()
//     operator()(polygon)                -> triangulatePolygon(polygon)
//     operator()(outer, inner)           -> triangulateWithHole(outer, inner)
//     operator()(outer, inners)          -> triangulateWithHoles(outer, inners)
//     operator()(tree)                   -> triangulateTree(tree)
// - Output reference parameters become returned objects
//   (ComputeNearestOuterPolygonIntersection) or returned arrays
//   (CombineSingle, CombineMultiple, InsertBridge).
// - 'std::numeric_limits<size_t>::max()' becomes Number.MAX_SAFE_INTEGER,
//   matching the PolygonTreeEx.INVALID precedent.
// - 'std::numeric_limits<InputType>::max()' becomes Number.MAX_VALUE.
// - The nested Vertex struct and VertexList class become module-private
//   declarations, because exported names must be unique library-wide.

import { logAssert } from './Logger';
import { PolygonTree } from './PolygonTree';
import { PrimalQuery2 } from './PrimalQuery2';
import { Vector, dot, sub } from './Vector';
import { dotPerp } from './Vector2';

// The fundamental problem is to compute the triangulation of a polygon tree.
// The outer polygons have counterclockwise ordered vertices. The inner
// polygons have clockwise ordered vertices. A polygon is an array of indices
// into the shared point array.
export type TriangulateECPolygon = number[];

// The port of 'std::numeric_limits<size_t>::max()'.
const INVALID = Number.MAX_SAFE_INTEGER;

// A vertex of the doubly linked list used for ear clipping. The vertices are
// specially tagged (convex, reflex, ear).
class ECVertex {
    // Index of the vertex in the points array.
    index: number = -1;

    // Vertex links for the polygon.
    vPrev: number = -1;
    vNext: number = -1;

    // Convex/reflex vertex links (disjoint lists).
    sPrev: number = -1;
    sNext: number = -1;

    // Ear links.
    ePrev: number = -1;
    eNext: number = -1;

    isConvex: boolean = false;
    isEar: boolean = false;
}

// A doubly linked list for storing specially tagged vertices (convex, reflex,
// ear). The vertex list is used for ear clipping.
class ECVertexList {
    private mVertices: ECVertex[] = [];

    // Linear list of convex vertices.
    private mCFirst: number = -1;
    private mCLast: number = -1;

    // Linear list of reflex vertices.
    private mRFirst: number = -1;
    private mRLast: number = -1;

    // Cyclical list of ears.
    private mEFirst: number = -1;
    private mELast: number = -1;

    doEarClipping(polygon: TriangulateECPolygon, computePoints: readonly Vector[],
        query: PrimalQuery2, triangles: [number, number, number][]): void {
        triangles.length = 0;

        // Initialize the vertex list for the incoming polygon. The lists must
        // be cleared in case a single vertex-list object is used two or more
        // times in triangulation queries. This is the case for triangulating a
        // polygon tree. It is also the case if you use a single TriangulateEC
        // object for multiple triangulation queries.
        this.mVertices = new Array<ECVertex>(polygon.length);
        this.mCFirst = -1;
        this.mCLast = -1;
        this.mRFirst = -1;
        this.mRLast = -1;
        this.mEFirst = -1;
        this.mELast = -1;

        // Create a circular list of the polygon vertices for dynamic removal
        // of vertices.
        let numVertices = polygon.length;
        const indices = polygon;
        for (let i = 0, ip1 = 1; ip1 <= numVertices; i = ip1++) {
            const vertex = new ECVertex();
            this.mVertices[i] = vertex;
            vertex.index = indices[i];
            vertex.vPrev = (i > 0 ? i - 1 : numVertices - 1);
            vertex.vNext = (ip1 < numVertices ? ip1 : 0);
            vertex.sPrev = -1;
            vertex.sNext = -1;
            vertex.ePrev = -1;
            vertex.eNext = -1;
            vertex.isConvex = false;
            vertex.isEar = false;
        }

        // Keep track of two linear sublists, one for the convex vertices and
        // one for the reflex vertices. This is an O(N) process where N is the
        // number of polygon vertices.
        for (let i = 0; i < numVertices; ++i) {
            if (this.isConvex(i, query)) {
                this.insertAfterC(i);
            } else {
                this.insertAfterR(i);
            }
        }

        // If the polygon is convex, create a triangle fan.
        if (this.mRFirst === -1) {
            for (let i = 1, ip1 = 2; ip1 < numVertices; i = ip1++) {
                triangles.push([polygon[0], polygon[i], polygon[ip1]]);
            }
            return;
        }

        // Identify the ears and build a circular list of them. Let V0, V1 and
        // V2 be consecutive vertices forming triangle T. The vertex V1 is an
        // ear if no other vertices of the polygon lie inside T. Although it is
        // enough to show that V1 is not an ear by finding at least one other
        // vertex inside T, it is sufficient to search only the reflex
        // vertices. This is an O(C*R) process, where C is the number of convex
        // vertices and R is the number of reflex vertices with N = C+R. The
        // order is O(N^2), for example when C = R = N/2.
        for (let i = this.mCFirst; i !== -1; i = this.V(i).sNext) {
            if (this.isEar(i, computePoints, query)) {
                this.insertEndE(i);
            }
        }
        this.V(this.mEFirst).ePrev = this.mELast;
        this.V(this.mELast).eNext = this.mEFirst;

        // Remove the ears, one at a time.
        let bRemoveAnEar = true;
        while (bRemoveAnEar) {
            // Add the triangle with the ear to the output list of triangles.
            let iVPrev = this.V(this.mEFirst).vPrev;
            let iVNext = this.V(this.mEFirst).vNext;
            triangles.push([this.V(iVPrev).index, this.V(this.mEFirst).index,
                this.V(iVNext).index]);

            // Remove the vertex corresponding to the ear.
            this.removeV(this.mEFirst);
            if (--numVertices === 3) {
                // Only one triangle remains, just remove the ear and copy it.
                this.mEFirst = this.removeE(this.mEFirst);
                iVPrev = this.V(this.mEFirst).vPrev;
                iVNext = this.V(this.mEFirst).vNext;
                triangles.push([this.V(iVPrev).index, this.V(this.mEFirst).index,
                    this.V(iVNext).index]);
                bRemoveAnEar = false;
                continue;
            }

            // Removal of the ear can cause an adjacent vertex to become an ear
            // or to stop being an ear.
            const vPrev = this.V(iVPrev);
            if (vPrev.isEar) {
                if (!this.isEar(iVPrev, computePoints, query)) {
                    this.removeE(iVPrev);
                }
            } else {
                const wasReflex = !vPrev.isConvex;
                if (this.isConvex(iVPrev, query)) {
                    if (wasReflex) {
                        this.removeR(iVPrev);
                    }

                    if (this.isEar(iVPrev, computePoints, query)) {
                        this.insertBeforeE(iVPrev);
                    }
                }
            }

            const vNext = this.V(iVNext);
            if (vNext.isEar) {
                if (!this.isEar(iVNext, computePoints, query)) {
                    this.removeE(iVNext);
                }
            } else {
                const wasReflex = !vNext.isConvex;
                if (this.isConvex(iVNext, query)) {
                    if (wasReflex) {
                        this.removeR(iVNext);
                    }

                    if (this.isEar(iVNext, computePoints, query)) {
                        this.insertAfterE(iVNext);
                    }
                }
            }

            // Remove the ear.
            this.mEFirst = this.removeE(this.mEFirst);
        }
    }

    private V(i: number): ECVertex {
        // If the assertion is triggered, do you have a coincident vertex-edge
        // or edge-edge pair? These violate the assumptions for the algorithm.
        logAssert(0 <= i && i < this.mVertices.length, 'Index out of range..');
        return this.mVertices[i];
    }

    private isConvex(i: number, query: PrimalQuery2): boolean {
        const vertex = this.V(i);
        const curr = vertex.index;
        const prev = this.V(vertex.vPrev).index;
        const next = this.V(vertex.vNext).index;
        vertex.isConvex = (query.toLine(curr, prev, next) > 0);
        return vertex.isConvex;
    }

    private isEar(i: number, computePoints: readonly Vector[],
        query: PrimalQuery2): boolean {
        const vertex = this.V(i);

        if (this.mRFirst === -1) {
            // The remaining polygon is convex.
            vertex.isEar = true;
            return true;
        }

        // Search the reflex vertices and test if any are in the triangle
        // <V[prev],V[curr],V[next]>.
        const prev = this.V(vertex.vPrev).index;
        const curr = vertex.index;
        const next = this.V(vertex.vNext).index;
        vertex.isEar = true;
        for (let j = this.mRFirst; j !== -1; j = this.V(j).sNext) {
            // Check if the test vertex is already one of the triangle
            // vertices.
            if (j === vertex.vPrev || j === i || j === vertex.vNext) {
                continue;
            }

            // V[j] has been ruled out as one of the original vertices of the
            // triangle <V[prev],V[curr],V[next]>. When triangulating polygons
            // with holes, V[j] might be a duplicated vertex, in which case it
            // does not affect the earness of V[curr].
            const testIndex = this.V(j).index;
            const testPoint = computePoints[testIndex];
            if (testPoint.equals(computePoints[prev])
                || testPoint.equals(computePoints[curr])
                || testPoint.equals(computePoints[next])) {
                continue;
            }

            // Test if the vertex is inside or on the triangle. When it is, it
            // causes V[curr] not to be an ear.
            if (query.toTriangle(testIndex, prev, curr, next) <= 0) {
                vertex.isEar = false;
                break;
            }
        }

        return vertex.isEar;
    }

    // Insert a convex vertex.
    private insertAfterC(i: number): void {
        if (this.mCFirst === -1) {
            // Insert the first convex vertex.
            this.mCFirst = i;
        } else {
            this.V(this.mCLast).sNext = i;
            this.V(i).sPrev = this.mCLast;
        }
        this.mCLast = i;
    }

    // Insert a reflex vertex.
    private insertAfterR(i: number): void {
        if (this.mRFirst === -1) {
            // Insert the first reflex vertex.
            this.mRFirst = i;
        } else {
            this.V(this.mRLast).sNext = i;
            this.V(i).sPrev = this.mRLast;
        }
        this.mRLast = i;
    }

    // Insert an ear at the end of the list.
    private insertEndE(i: number): void {
        if (this.mEFirst === -1) {
            // Insert the first ear.
            this.mEFirst = i;
            this.mELast = i;
        }
        this.V(this.mELast).eNext = i;
        this.V(i).ePrev = this.mELast;
        this.mELast = i;
    }

    // Insert an ear after mEFirst.
    private insertAfterE(i: number): void {
        const first = this.V(this.mEFirst);
        const currENext = first.eNext;
        const vertex = this.V(i);
        vertex.ePrev = this.mEFirst;
        vertex.eNext = currENext;
        first.eNext = i;
        this.V(currENext).ePrev = i;
    }

    // Insert an ear before mEFirst.
    private insertBeforeE(i: number): void {
        const first = this.V(this.mEFirst);
        const currEPrev = first.ePrev;
        const vertex = this.V(i);
        vertex.ePrev = currEPrev;
        vertex.eNext = this.mEFirst;
        first.ePrev = i;
        this.V(currEPrev).eNext = i;
    }

    // Remove a vertex.
    private removeV(i: number): void {
        const currVPrev = this.V(i).vPrev;
        const currVNext = this.V(i).vNext;
        this.V(currVPrev).vNext = currVNext;
        this.V(currVNext).vPrev = currVPrev;
    }

    // Remove an ear.
    private removeE(i: number): number {
        const currEPrev = this.V(i).ePrev;
        const currENext = this.V(i).eNext;
        this.V(currEPrev).eNext = currENext;
        this.V(currENext).ePrev = currEPrev;
        return currENext;
    }

    // Remove a reflex vertex.
    private removeR(i: number): void {
        logAssert(this.mRFirst !== -1 && this.mRLast !== -1,
            'Reflex vertices must exist.');

        if (i === this.mRFirst) {
            this.mRFirst = this.V(i).sNext;
            if (this.mRFirst !== -1) {
                this.V(this.mRFirst).sPrev = -1;
            }
            this.V(i).sNext = -1;
        } else if (i === this.mRLast) {
            this.mRLast = this.V(i).sPrev;
            if (this.mRLast !== -1) {
                this.V(this.mRLast).sNext = -1;
            }
            this.V(i).sPrev = -1;
        } else {
            const currSPrev = this.V(i).sPrev;
            const currSNext = this.V(i).sNext;
            this.V(currSPrev).sNext = currSNext;
            this.V(currSNext).sPrev = currSPrev;
            this.V(i).sNext = -1;
            this.V(i).sPrev = -1;
        }
    }
}

// The result of the private nearest-outer-polygon-intersection helper.
interface NearestIntersection {
    intr: Vector;
    v0min: number;
    v1min: number;
    endMin: number;
}

export class TriangulateEC {
    // The input vertex pool.
    private readonly mNumPoints: number;
    private readonly mPoints: readonly Vector[];

    // The output triangulation.
    private mTriangles: [number, number, number][];

    // The array of points used for geometric queries. The InputType points are
    // converted to ComputeType points on demand. The mConverted array keeps
    // track of which input points have been converted.
    private mComputePoints: Vector[];
    private mConverted: boolean[];

    // The object used for toLine and toTriangle queries.
    private mQuery: PrimalQuery2;

    private mVertexList: ECVertexList;

    // The class is a functor to support triangulating multiple polygons that
    // share vertices in a collection of points. The interpretation of
    // 'numPoints' and 'points' is described before each triangulate function.
    // Preconditions are numPoints >= 3 and points has at least numPoints
    // elements. If they are not satisfied, an exception is thrown. Upstream
    // has a second constructor taking only the point array; here 'numPoints'
    // defaults to points.length.
    constructor(points: readonly Vector[], numPoints: number = points.length) {
        this.mNumPoints = numPoints;
        this.mPoints = points;
        this.mTriangles = [];
        this.mVertexList = new ECVertexList();

        logAssert(numPoints >= 3 && points.length >= numPoints, 'Invalid input.');

        this.mComputePoints = new Array<Vector>(this.mNumPoints);
        for (let i = 0; i < this.mNumPoints; ++i) {
            this.mComputePoints[i] = new Vector(2);
        }
        this.mConverted = new Array<boolean>(this.mNumPoints).fill(false);
        this.mQuery = new PrimalQuery2(this.mNumPoints, this.mComputePoints);
    }

    // Access the triangulation after each triangulate call.
    getTriangles(): [number, number, number][] {
        return this.mTriangles;
    }

    // The input 'points' represents an array of vertices for a simple polygon.
    // The vertices are points[0] through points[numPoints-1] and are listed in
    // counterclockwise order.
    triangulate(): void {
        this.mTriangles = [];
        const polygon: TriangulateECPolygon = new Array<number>(this.mNumPoints);
        for (let i = 0; i < this.mNumPoints; ++i) {
            polygon[i] = i;
        }
        this.triangulatePolygon(polygon);
    }

    // The input 'points' represents an array of vertices that contains the
    // vertices of a simple polygon.
    triangulatePolygon(polygon: TriangulateECPolygon): void {
        this.mTriangles = [];

        // Convert InputType polygon vertices to ComputeType.
        this.convertPolygon(polygon);

        // Triangulate the simple polygon using ear clipping.
        this.mVertexList.doEarClipping(polygon, this.mComputePoints, this.mQuery,
            this.mTriangles);
    }

    // The input 'points' is a shared array of vertices that contains the
    // vertices for two simple polygons, an outer polygon and an inner polygon.
    // The inner polygon must be strictly inside the outer polygon.
    triangulateWithHole(outer: TriangulateECPolygon,
        inner: TriangulateECPolygon): void {
        this.mTriangles = [];

        // Convert InputType polygon vertices to ComputeType.
        this.convertPolygon(outer);
        this.convertPolygon(inner);

        // Combine the inner and outer polygon into a pseudosimple polygon.
        const combined = this.combineSingle(outer, inner);

        // Triangulate the pseudosimple polygon using ear clipping.
        this.mVertexList.doEarClipping(combined, this.mComputePoints, this.mQuery,
            this.mTriangles);
    }

    // The input 'points' is a shared array of vertices that contains the
    // vertices for multiple simple polygons, an outer polygon and one or more
    // inner polygons. The inner polygons must be nonoverlapping and strictly
    // inside the outer polygon.
    triangulateWithHoles(outer: TriangulateECPolygon,
        inners: readonly TriangulateECPolygon[]): void {
        this.mTriangles = [];

        // Convert InputType polygon vertices to ComputeType.
        this.convertPolygon(outer);
        for (const inner of inners) {
            this.convertPolygon(inner);
        }

        // Combine the outer polygon and the inner polygons into a pseudosimple
        // polygon using repeated calls to combineSingle.
        const combined = this.combineMultiple(outer, inners);

        // Triangulate the pseudosimple polygon using ear clipping.
        this.mVertexList.doEarClipping(combined, this.mComputePoints, this.mQuery,
            this.mTriangles);
    }

    // The input 'points' is a shared array of vertices that contains the
    // vertices for multiple simple polygons in a tree of polygons.
    triangulateTree(tree: PolygonTree): void {
        this.mTriangles = [];

        // Convert InputType polygon vertices to ComputeType.
        this.convertTree(tree);

        const treeQueue: PolygonTree[] = [tree];
        let head = 0;
        while (head < treeQueue.length) {
            const outer = treeQueue[head++];

            const numChildren = outer.child.length;
            if (numChildren === 0) {
                // The outer polygon is a simple polygon that has no nested
                // inner polygons. Triangulate the pseudosimple polygon using
                // ear clipping.
                const combinedTriangles: [number, number, number][] = [];
                this.mVertexList.doEarClipping(outer.polygon, this.mComputePoints,
                    this.mQuery, combinedTriangles);
                for (const t of combinedTriangles) {
                    this.mTriangles.push(t);
                }
            } else {
                // Place the next level of outer polygon nodes on the queue for
                // triangulation.
                const inners: TriangulateECPolygon[] = new Array<TriangulateECPolygon>(numChildren);
                for (let c = 0; c < numChildren; ++c) {
                    const inner = outer.child[c];
                    inners[c] = inner.polygon;
                    const numGrandChildren = inner.child.length;
                    for (let g = 0; g < numGrandChildren; ++g) {
                        treeQueue.push(inner.child[g]);
                    }
                }

                // Combine the outer polygon and the inner polygons into a
                // pseudosimple polygon using repeated calls to combineSingle.
                const combined = this.combineMultiple(outer.polygon, inners);

                // Triangulate the pseudosimple polygon using ear clipping.
                const combinedTriangles: [number, number, number][] = [];
                this.mVertexList.doEarClipping(combined, this.mComputePoints,
                    this.mQuery, combinedTriangles);
                for (const t of combinedTriangles) {
                    this.mTriangles.push(t);
                }
            }
        }
    }

    // Support for rational arithmetic. The converter transforms points with
    // InputType components to points with ComputeType components. Upstream has
    // four ConvertPoints overloads that differ only in which polygons they
    // enumerate; the port keeps a single per-polygon converter and a tree
    // converter, which enumerate the same indices in the same order.
    private convertPolygon(polygon: TriangulateECPolygon): void {
        for (const index of polygon) {
            if (!this.mConverted[index]) {
                this.mConverted[index] = true;
                for (let j = 0; j < 2; ++j) {
                    this.mComputePoints[index].values[j] = this.mPoints[index].values[j];
                }
            }
        }
    }

    private convertTree(tree: PolygonTree): void {
        const treeQueue: PolygonTree[] = [tree];
        let head = 0;
        while (head < treeQueue.length) {
            // The 'root' is an outer polygon.
            const outer = treeQueue[head++];
            this.convertPolygon(outer.polygon);

            // The grandchildren of the outer polygon are also outer polygons.
            // Insert them into the queue for processing.
            const numChildren = outer.child.length;
            for (let c = 0; c < numChildren; ++c) {
                // The 'child' is an inner polygon.
                const inner = outer.child[c];
                this.convertPolygon(inner.polygon);

                const numGrandChildren = inner.child.length;
                for (let g = 0; g < numGrandChildren; ++g) {
                    treeQueue.push(inner.child[g]);
                }
            }
        }
    }

    // The returned 'xmax' is the maximum x-value of the polygon vertices. The
    // returned 'index' is the index into polygon[] of a vertex that generates
    // a maximum x-value. It is not a problem if the maximum is attained by
    // more than one vertex. It is sufficient to use mPoints directly because
    // the InputType comparisons are exact.
    private getXMaxInfo(polygon: TriangulateECPolygon):
        { xmax: number, index: number } {
        let x = this.mPoints[polygon[0]].values[0];
        const xmaxInfo = { xmax: x, index: 0 };
        for (let i = 1; i < polygon.length; ++i) {
            x = this.mPoints[polygon[i]].values[0];
            if (x > xmaxInfo.xmax) {
                xmaxInfo.xmax = x;
                xmaxInfo.index = i;
            }
        }
        return xmaxInfo;
    }

    // Find the edge whose intersection Intr with the ray M + t * (1,0)
    // minimizes the ray parameter t >= 0.
    private computeNearestOuterPolygonIntersection(M: Vector,
        outer: TriangulateECPolygon): NearestIntersection {
        const cmax = Number.MAX_VALUE;
        const zero = 0;

        const intr = Vector.fromArray([cmax, M.values[1]]);
        let v0min = INVALID;
        let v1min = INVALID;
        let endMin = INVALID;
        let t = cmax;
        const numOuter = outer.length;
        for (let i0 = numOuter - 1, i1 = 0; i1 < numOuter; i0 = i1++) {
            // Consider only edges for which the first vertex is below (or on)
            // the ray and the second vertex is above (or on) the ray.
            let diff0 = sub(this.mComputePoints[outer[i0]], M);
            if (diff0.values[1] > zero) {
                continue;
            }
            let diff1 = sub(this.mComputePoints[outer[i1]], M);
            if (diff1.values[1] < zero) {
                continue;
            }

            // At this time, diff0.y <= 0 and diff1.y >= 0.
            let currentEndMin = INVALID;
            if (diff0.values[1] < zero) {
                if (diff1.values[1] > zero) {
                    // The intersection of the edge and ray occurs at an
                    // interior edge point.
                    const s = diff0.values[1] / (diff0.values[1] - diff1.values[1]);
                    t = diff0.values[0] + s * (diff1.values[0] - diff0.values[0]);
                } else {
                    // diff1.y == 0. The vertex outer[i1] is the intersection
                    // of the edge and the ray.
                    t = diff1.values[0];
                    currentEndMin = i1;
                }
            } else {
                // diff0.y == 0.
                if (diff1.values[1] > zero) {
                    // The vertex outer[i0] is the intersection of the edge and
                    // the ray.
                    t = diff0.values[0];
                    currentEndMin = i0;
                } else {
                    // diff1.y == 0.
                    if (diff0.values[0] < diff1.values[0]) {
                        t = diff0.values[0];
                        currentEndMin = i0;
                    } else {
                        t = diff1.values[0];
                        currentEndMin = i1;
                    }
                }
            }

            if (zero <= t && t < intr.values[0]) {
                intr.values[0] = t;
                v0min = i0;
                v1min = i1;
                if (currentEndMin === INVALID) {
                    // The current closest point is an edge-interior point.
                    endMin = INVALID;
                } else {
                    // The current closest point is a vertex.
                    endMin = currentEndMin;
                }
            } else if (t === intr.values[0]) {
                // The current closest point is a vertex shared by multiple
                // edges; thus, endMin and currentEndMin refer to the same
                // point.
                logAssert(endMin !== INVALID && currentEndMin !== INVALID,
                    'Unexpected condition.');

                // We need to select the edge closest to M. The previous
                // closest edge is <outer[v0min],outer[v1min]>. The current
                // candidate is <outer[i0],outer[i1]>.
                const shared = this.mComputePoints[outer[i1]];

                // For the previous closest edge, endMin refers to a vertex of
                // the edge. Get the index of the other vertex.
                const other = (endMin === v0min ? v1min : v0min);

                // The new edge is closer if the other vertex of the old edge is
                // left-of the new edge.
                diff0 = sub(this.mComputePoints[outer[i0]], shared);
                diff1 = sub(this.mComputePoints[outer[other]], shared);
                const dotperp = dotPerp(diff0, diff1);
                if (dotperp > zero) {
                    // The new edge is closer to M.
                    v0min = i0;
                    v1min = i1;
                    endMin = currentEndMin;
                }
            }
        }

        // The intersection intr[0] stored only the t-value of the ray. The
        // actual point is (mx,my)+t*(1,0), so intr[0] must be adjusted.
        intr.values[0] += M.values[0];
        return { intr, v0min, v1min, endMin };
    }

    private locateOuterVisibleVertex(M: Vector, I: Vector,
        outer: TriangulateECPolygon, v0min: number, v1min: number,
        endMin: number): number {
        // The point mPoints[outer[oVisibleIndex]] maximizes the cosine of the
        // angle between <M,I> and <M,Q> where Q is P or a reflex vertex
        // contained in triangle <M,I,P>.
        let oVisibleIndex = endMin;
        if (endMin === INVALID) {
            // If you reach this assert, there is a good chance that you have
            // two inner polygons that share a vertex or an edge.
            logAssert(v0min !== INVALID && v1min !== INVALID,
                'Is this an invalid nested polygon?');

            // Select mPoints[outer[v0min]] or mPoints[outer[v1min]] that has
            // an x-value larger than M.x, call this vertex P. The triangle
            // <M,I,P> must contain an outer-polygon vertex that is visible to
            // M, which is possibly P itself.
            const triangle: Vector[] = new Array<Vector>(3);
            let pIndex: number;
            if (this.mComputePoints[outer[v0min]].values[0]
                > this.mComputePoints[outer[v1min]].values[0]) {
                const P = this.mComputePoints[outer[v0min]];
                triangle[0] = P;
                triangle[1] = I;
                triangle[2] = M;
                pIndex = v0min;
            } else {
                const P = this.mComputePoints[outer[v1min]];
                triangle[0] = P;
                triangle[1] = M;
                triangle[2] = I;
                pIndex = v1min;
            }

            // If any outer-polygon vertices other than P are inside the
            // triangle <M,I,P>, then at least one of these vertices must be a
            // reflex vertex. It is sufficient to locate the reflex vertex R
            // (if any) in <M,I,P> that minimizes the angle between R-M and
            // (1,0).
            let diff = sub(triangle[0], M);
            let maxSqrLen = dot(diff, diff);
            let maxCos = diff.values[0] * diff.values[0] / maxSqrLen;
            const localQuery = new PrimalQuery2(3, triangle);
            const numOuter = outer.length;
            oVisibleIndex = pIndex;
            for (let i = 0; i < numOuter; ++i) {
                if (i === pIndex) {
                    continue;
                }

                const curr = outer[i];
                const prev = outer[(i + numOuter - 1) % numOuter];
                const next = outer[(i + 1) % numOuter];
                if (this.mQuery.toLine(curr, prev, next) <= 0
                    && localQuery.toTriangle(this.mComputePoints[curr], 0, 1, 2) <= 0) {
                    // The vertex is reflex and inside the <M,I,P> triangle.
                    diff = sub(this.mComputePoints[curr], M);
                    const sqrLen = dot(diff, diff);
                    const cs = diff.values[0] * diff.values[0] / sqrLen;
                    if (cs > maxCos) {
                        // The reflex vertex forms a smaller angle with the
                        // positive x-axis, so it becomes the new visible
                        // candidate.
                        maxSqrLen = sqrLen;
                        maxCos = cs;
                        oVisibleIndex = i;
                    } else if (cs === maxCos && sqrLen < maxSqrLen) {
                        // The reflex vertex has angle equal to the current
                        // minimum but the length is smaller, so it becomes the
                        // new visible candidate.
                        maxSqrLen = sqrLen;
                        oVisibleIndex = i;
                    }
                }
            }
        }

        return oVisibleIndex;
    }

    private combineSingle(outer: TriangulateECPolygon,
        inner: TriangulateECPolygon): TriangulateECPolygon {
        // Get the index into inner[] for the inner-polygon vertex M of maximum
        // x-value.
        const iVisibleIndex = this.getXMaxInfo(inner).index;

        // Get the inner-polygon vertex M of maximum x-value.
        const iVertexIndex = inner[iVisibleIndex];
        const M = this.mComputePoints[iVertexIndex];

        // Compute the closest outer-polygon point I along the ray
        // M + t * (1,0) with t > 0 so that M and I are mutually visible.
        const nearest = this.computeNearestOuterPolygonIntersection(M, outer);

        // Locate Q = mPoints[outer[oVisibleIndex]] so that M and Q are
        // mutually visible.
        const oVisibleIndex = this.locateOuterVisibleVertex(M, nearest.intr,
            outer, nearest.v0min, nearest.v1min, nearest.endMin);

        return TriangulateEC.insertBridge(outer, inner, oVisibleIndex, iVisibleIndex);
    }

    private combineMultiple(outer: TriangulateECPolygon,
        inners: readonly TriangulateECPolygon[]): TriangulateECPolygon {
        // Sort the inner polygons based on maximum x-values.
        const numInners = inners.length;
        const pairs: { xmax: number, index: number }[] =
            new Array<{ xmax: number, index: number }>(numInners);
        for (let p = 0; p < numInners; ++p) {
            const indices = inners[p];
            const numIndices = indices.length;

            let xmax = this.mPoints[indices[0]].values[0];
            for (let i = 1; i < numIndices; ++i) {
                const x = this.mPoints[indices[i]].values[0];
                if (x > xmax) {
                    xmax = x;
                }
            }

            pairs[p] = { xmax, index: p };
        }

        // The port of 'std::sort(..., std::greater<PairType>())', a descending
        // lexicographic order on (xmax, index).
        pairs.sort((a, b) => {
            if (a.xmax !== b.xmax) {
                return b.xmax - a.xmax;
            }
            return b.index - a.index;
        });

        let currentOuter = outer;
        for (const pair of pairs) {
            const inner = inners[pair.index];
            currentOuter = this.combineSingle(currentOuter, inner);
        }
        return currentOuter;
    }

    // The mutually visible vertices are VI = mPoints[inner[iVisibleIndex]] and
    // VO = mPoints[outer[oVisibleIndex]]. Two coincident edges with these
    // endpoints are inserted to connect the outer and inner polygons into a
    // pseudosimple polygon.
    private static insertBridge(outer: TriangulateECPolygon,
        inner: TriangulateECPolygon, oVisibleIndex: number,
        iVisibleIndex: number): TriangulateECPolygon {
        const numOuter = outer.length;
        const numInner = inner.length;
        const combined: TriangulateECPolygon =
            new Array<number>(numOuter + numInner + 2).fill(0);

        // Traverse the outer polygon until the outer polygon bridge point is
        // visited.
        let cIndex = 0;
        for (let i = 0; i <= oVisibleIndex; ++i, ++cIndex) {
            combined[cIndex] = outer[i];
        }

        // Cross the bridge from the outer polygon to the inner polygon.
        // Traverse the inner polygon until the predecessor of the inner
        // polygon bridge point is visited.
        for (let i = 0; i < numInner; ++i, ++cIndex) {
            const j = (iVisibleIndex + i) % numInner;
            combined[cIndex] = inner[j];
        }

        // Inner polygon bridge point.
        combined[cIndex++] = inner[iVisibleIndex];

        // Outer polygon bridge point.
        combined[cIndex++] = outer[oVisibleIndex];

        for (let i = oVisibleIndex + 1; i < numOuter; ++i, ++cIndex) {
            combined[cIndex] = outer[i];
        }

        return combined;
    }
}
