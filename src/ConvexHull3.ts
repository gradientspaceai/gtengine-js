// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ConvexHull3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the convex hull of 3D points using incremental insertion. The only
// way to ensure a correct result for the input vertices is to use an exact
// predicate for computing signs of various expressions. The implementation
// uses interval arithmetic and rational arithmetic for the predicate.
//
// The main cost of the algorithm is testing which side of a plane a point is
// located. This test uses interval arithmetic to determine an exact sign, if
// possible. If that test fails, rational arithmetic is used. For typical
// datasets, the indeterminate sign from interval arithmetic happens rarely.
//
// Port notes:
// * The upstream operator() has an 'lgNumThreads' parameter that selects a
//   multithreaded divide-and-conquer merge of subhulls. JavaScript has no
//   equivalent shared-memory threading, so only the single-threaded path
//   (lgNumThreads = 0) is ported; compute() is the port of that path. The
//   single-threaded and multithreaded paths compute the same hull (upstream
//   notes the triangulation of coplanar faces can differ between them).
// * Upstream selects Rational = BSNumber<UIntegerFP32<27 or 197>>. The port's
//   BSNumber is bigint-backed and grows as needed, so the fixed word count is
//   unnecessary. The exact-arithmetic results are identical.
// * The rational points are memoized exactly as upstream (mRPoints with the
//   mConverted flags becomes an array of nullable rational points).
// * The std::set<Triangle*> 'visited' set and the std::queue 'visible' become
//   a Set and an array used as a FIFO. The seed triangle is selected by
//   iterating the vertex-adjacent triangles in sorted key order so the hull
//   does not depend on pointer-hash ordering.

import { logAssert } from './Logger';
import { Vector } from './Vector';
import { BSNumber } from './BSNumber';
import { SWInterval } from './SWInterval';
import { VETManifoldMesh } from './VETManifoldMesh';
import type { ETManifoldMeshTriangle } from './ETManifoldMesh';
import { ConvexHull2 } from './ConvexHull2';

// A rational 3D point, the port of Vector3<Rational>.
type RationalPoint3 = [BSNumber, BSNumber, BSNumber];

// The port of Cross(v0, v1) for rational vectors.
function rationalCross(v0: RationalPoint3, v1: RationalPoint3): RationalPoint3 {
    return [
        v0[1].mul(v1[2]).sub(v0[2].mul(v1[1])),
        v0[2].mul(v1[0]).sub(v0[0].mul(v1[2])),
        v0[0].mul(v1[1]).sub(v0[1].mul(v1[0]))
    ];
}

// The port of DotCross(v0, v1, v2) = Dot(v0, Cross(v1, v2)) for rational
// vectors.
function rationalDotCross(v0: RationalPoint3, v1: RationalPoint3,
    v2: RationalPoint3): BSNumber {
    const cross = rationalCross(v1, v2);
    let dot = v0[0].mul(cross[0]);
    dot = dot.add(v0[1].mul(cross[1]));
    dot = dot.add(v0[2].mul(cross[2]));
    return dot;
}

// The port of DotCross(v0, v1, v2) for interval vectors.
function intervalDotCross(v0: readonly SWInterval[], v1: readonly SWInterval[],
    v2: readonly SWInterval[]): SWInterval {
    const cross: SWInterval[] = [
        v1[1].mul(v2[2]).sub(v1[2].mul(v2[1])),
        v1[2].mul(v2[0]).sub(v1[0].mul(v2[2])),
        v1[0].mul(v2[1]).sub(v1[1].mul(v2[0]))
    ];
    let dot = v0[0].mul(cross[0]);
    dot = dot.add(v0[1].mul(cross[1]));
    dot = dot.add(v0[2].mul(cross[2]));
    return dot;
}

export class ConvexHull3 {
    // A blend of interval arithmetic and exact arithmetic is used to ensure
    // correctness.
    private mPoints: readonly Vector[];
    private mRPoints: (RationalPoint3 | null)[];

    // The output data.
    private mDimension: number;
    private mVertices: number[];
    private mHull: number[];
    private mHullMesh: VETManifoldMesh;

    // The class is a functor to support computing the convex hull of multiple
    // data sets using the same class object.
    constructor() {
        this.mPoints = [];
        this.mRPoints = [];
        this.mDimension = 0;
        this.mVertices = [];
        this.mHull = [];
        this.mHullMesh = new VETManifoldMesh();
    }

    // Compute the exact convex hull using a blend of interval arithmetic and
    // rational arithmetic.
    compute(points: readonly Vector[]): void {
        logAssert(points.length > 0, 'Invalid argument.');
        for (const point of points) {
            logAssert(point.size === 3, 'ConvexHull3 requires 3D points.');
        }

        // Allocate storage for any rational points that must be computed in
        // the exact sign predicates. The rational points are memoized.
        this.mPoints = points;
        this.mRPoints = new Array<RationalPoint3 | null>(points.length).fill(null);

        // Sort all the points indirectly, then remove duplicates.
        const sorted: number[] = new Array<number>(points.length);
        for (let i = 0; i < points.length; ++i) {
            sorted[i] = i;
        }
        sorted.sort((s0, s1) => {
            const p0 = this.mPoints[s0].values;
            const p1 = this.mPoints[s1].values;
            for (let i = 0; i < 3; ++i) {
                if (p0[i] !== p1[i]) {
                    return p0[i] < p1[i] ? -1 : +1;
                }
            }
            return s0 - s1;
        });
        const unique: number[] = [];
        for (const s of sorted) {
            if (unique.length === 0 ||
                !this.mPoints[unique[unique.length - 1]].equals(this.mPoints[s])) {
                unique.push(s);
            }
        }

        this.computeHull(unique);
    }

    // The dimension is 0 (hull is a single point), 1 (hull is a line
    // segment), 2 (hull is a convex polygon in 3D) or 3 (hull is a convex
    // polyhedron).
    getDimension(): number {
        return this.mDimension;
    }

    // Get the indices into the input 'points[]' that correspond to hull
    // vertices.
    getVertices(): readonly number[] {
        return this.mVertices;
    }

    // Get the indices into the input 'points[]' that correspond to hull
    // vertices. The returned array is organized according to the hull
    // dimension.
    //   0: The hull is a single point. The returned array has size 1 with the
    //      index corresponding to that point.
    //   1: The hull is a line segment. The returned array has size 2 with the
    //      indices corresponding to the segment endpoints.
    //   2: The hull is a convex polygon in 3D. The returned array has size N
    //      with indices corresponding to the ordered polygon vertices.
    //   3: The hull is a convex polyhedron. The returned array has T triples
    //      of indices, each triple corresponding to a triangle face of the
    //      hull. The face vertices are counterclockwise when viewed by an
    //      observer outside the polyhedron. It is possible that some triangle
    //      faces are coplanar.
    // If V is the number of vertices and T is the number of triangles, then
    // the number of edges is E = 3*T/2 and Euler's formula V - E + T = 2 is
    // satisfied. (Upstream documentation bug: the header states E = T/2, which
    // does not satisfy Euler's formula; for a closed triangle mesh every
    // triangle has 3 edges and every edge is shared by 2 triangles, so
    // E = 3*T/2. Example: the hull of the 8 cube corners has V = 8, T = 12 and
    // E = 18, and 8 - 18 + 12 = 2.)
    getHull(): readonly number[] {
        return this.mHull;
    }

    // Get the hull mesh, which is valid only when the dimension is 3. This
    // allows access to the graph of vertices, edges and triangles of the
    // convex (polyhedron) hull.
    getHullMesh(): VETManifoldMesh {
        return this.mHullMesh;
    }

    private computeHull(sorted: readonly number[]): void {
        this.mDimension = 0;
        this.mVertices = [];
        this.mHull = [];
        this.mHullMesh.clear();

        const hull: number[] = [];
        const state = { current: 0 };

        if (this.hull0(hull, sorted, state)) {
            this.mHull = hull;
            this.mVertices = [hull[0]];
            return;
        }

        if (this.hull1(hull, sorted, state)) {
            this.mHull = hull;
            this.mVertices = [hull[0], hull[1]];
            return;
        }

        if (this.hull2(hull, sorted, state)) {
            this.mHull = hull;
            this.mVertices = hull.slice();
            return;
        }

        this.hull3(hull, sorted, state);

        this.mVertices = this.mHullMesh.getVertices().map(vertex => vertex.V);

        const tkeys = this.mHullMesh.getTriangleKeys();
        this.mHull = new Array<number>(3 * tkeys.length);
        let index = 0;
        for (const tkey of tkeys) {
            this.mHull[index++] = tkey.V[0];
            this.mHull[index++] = tkey.V[1];
            this.mHull[index++] = tkey.V[2];
        }
    }

    // Support for computing a 0-dimensional convex hull.
    private hull0(hull: number[], sorted: readonly number[],
        state: { current: number }): boolean {
        hull.push(sorted[state.current]);  // hull[0]
        for (++state.current; state.current < sorted.length; ++state.current) {
            if (!this.colocated(hull[0], sorted[state.current])) {
                this.mDimension = 1;
                break;
            }
        }
        return this.mDimension === 0;
    }

    // Support for computing a 1-dimensional convex hull.
    private hull1(hull: number[], sorted: readonly number[],
        state: { current: number }): boolean {
        hull.push(sorted[state.current]);  // hull[1]
        for (++state.current; state.current < sorted.length; ++state.current) {
            if (!this.colinear(hull[0], hull[1], sorted[state.current])) {
                this.mDimension = 2;
                break;
            }
            hull.push(sorted[state.current]);
        }

        if (hull.length > 2) {
            // Sort the points and choose the extreme points as the endpoints
            // of the line segment that is the convex hull.
            hull.sort((v0, v1) => {
                const p0 = this.mPoints[v0].values;
                const p1 = this.mPoints[v1].values;
                for (let i = 0; i < 3; ++i) {
                    if (p0[i] !== p1[i]) {
                        return p0[i] < p1[i] ? -1 : +1;
                    }
                }
                return v0 - v1;
            });

            const hmin = hull[0];
            const hmax = hull[hull.length - 1];
            hull.length = 0;
            hull.push(hmin);
            hull.push(hmax);
        }

        return this.mDimension === 1;
    }

    // Support for computing a 2-dimensional convex hull.
    private hull2(hull: number[], sorted: readonly number[],
        state: { current: number }): boolean {
        hull.push(sorted[state.current]);  // hull[2]
        for (++state.current; state.current < sorted.length; ++state.current) {
            if (this.toPlane(hull[0], hull[1], hull[2], sorted[state.current]) !== 0) {
                this.mDimension = 3;
                break;
            }
            hull.push(sorted[state.current]);
        }

        if (hull.length > 3) {
            // Compute the planar convex hull of the points. The coplanar
            // points are projected onto a 2D plane determined by the maximum
            // absolute component of the normal of the first triangle. The
            // extreme points of the projected hull generate the extreme
            // points of the planar hull in 3D.
            const rV0 = this.getRationalPoint(hull[0]);
            const rV1 = this.getRationalPoint(hull[1]);
            const rV2 = this.getRationalPoint(hull[2]);
            const rDiff1: RationalPoint3 =
                [rV1[0].sub(rV0[0]), rV1[1].sub(rV0[1]), rV1[2].sub(rV0[2])];
            const rDiff2: RationalPoint3 =
                [rV2[0].sub(rV0[0]), rV2[1].sub(rV0[1]), rV2[2].sub(rV0[2])];
            const rNormal = rationalCross(rDiff1, rDiff2);

            // The signs are used to select 2 of the 3 point components so
            // that when the planar hull is viewed from the side of the plane
            // to which rNormal is directed, the triangles are counterclockwise
            // ordered.
            const sign: [number, number, number] = [0, 0, 0];
            for (let i = 0; i < 3; ++i) {
                sign[i] = rNormal[i].getSign();
                rNormal[i].setSign(Math.abs(sign[i]));
            }

            let c: [number, number];
            if (rNormal[0].greaterThan(rNormal[1])) {
                if (rNormal[0].greaterThan(rNormal[2])) {
                    c = (sign[0] > 0 ? [1, 2] : [2, 1]);
                }
                else {
                    c = (sign[2] > 0 ? [0, 1] : [1, 0]);
                }
            }
            else {
                if (rNormal[1].greaterThan(rNormal[2])) {
                    c = (sign[1] > 0 ? [2, 0] : [0, 2]);
                }
                else {
                    c = (sign[2] > 0 ? [0, 1] : [1, 0]);
                }
            }

            const projections: Vector[] = new Array<Vector>(hull.length);
            for (let i = 0; i < projections.length; ++i) {
                const h = hull[i];
                projections[i] = Vector.fromArray(
                    [this.mPoints[h].values[c[0]], this.mPoints[h].values[c[1]]]);
            }

            const ch2 = new ConvexHull2();
            ch2.compute(projections);
            const hull2 = ch2.getHull();

            const tempHull: number[] = new Array<number>(hull2.length);
            for (let i = 0; i < hull2.length; ++i) {
                tempHull[i] = hull[hull2[i]];
            }
            hull.length = 0;
            for (let i = 0; i < hull2.length; ++i) {
                hull.push(tempHull[i]);
            }
        }

        return this.mDimension === 2;
    }

    // Support for computing a 3-dimensional convex hull.
    private hull3(hull: readonly number[], sorted: readonly number[],
        state: { current: number }): void {
        const hullMesh = this.mHullMesh;

        // The hull points previous to the current one are coplanar and are
        // the vertices of a convex polygon. To initialize the 3D hull, use
        // triangles from a triangle fan of the convex polygon and use
        // triangles connecting the current point to the edges of the convex
        // polygon. The vertex ordering of these triangles depends on whether
        // sorted[current] is on the positive or negative side of the plane
        // determined by hull[0], hull[1] and hull[2].
        let sign = this.toPlane(hull[0], hull[1], hull[2], sorted[state.current]);
        let h0: number, h1: number, h2: number;
        if (sign > 0) {
            h0 = hull[0];
            for (let i1 = 1, i2 = 2; i2 < hull.length; i1 = i2++) {
                h1 = hull[i1];
                h2 = hull[i2];
                const inserted = hullMesh.insert(h0, h2, h1);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }

            h0 = sorted[state.current];
            for (let i1 = hull.length - 1, i2 = 0; i2 < hull.length; i1 = i2++) {
                h1 = hull[i1];
                h2 = hull[i2];
                const inserted = hullMesh.insert(h0, h1, h2);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }
        }
        else {
            h0 = hull[0];
            for (let i1 = 1, i2 = 2; i2 < hull.length; i1 = i2++) {
                h1 = hull[i1];
                h2 = hull[i2];
                const inserted = hullMesh.insert(h0, h1, h2);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }

            h0 = sorted[state.current];
            for (let i1 = hull.length - 1, i2 = 0; i2 < hull.length; i1 = i2++) {
                h1 = hull[i1];
                h2 = hull[i2];
                const inserted = hullMesh.insert(h0, h2, h1);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }
        }

        // The hull is now maintained in mHullMesh, so there is no need to add
        // members to 'hull'. At the time the full hull is known, mHull is
        // assigned the triangle indices.
        for (++state.current; state.current < sorted.length; ++state.current) {
            // The index h0 refers to the previously inserted hull point. The
            // index h1 refers to the current point to be inserted into the
            // hull.
            const vertex = hullMesh.getVertex(h0);
            logAssert(vertex !== null, 'Unexpected condition.');
            h1 = sorted[state.current];

            // The sorting guarantees that the point at h0 is visible to the
            // point at h1. Find the triangles that share h0 and are visible
            // to h1.
            const visible: ETManifoldMeshTriangle[] = [];
            let visibleHead = 0;
            const visited = new Set<ETManifoldMeshTriangle>();
            for (const tri of vertex.getTAdjacent()) {
                sign = this.toPlane(tri.V[0], tri.V[1], tri.V[2], h1);
                if (sign > 0) {
                    visible.push(tri);
                    visited.add(tri);
                    break;
                }
            }
            logAssert(visible.length > 0, 'Unexpected condition.');

            // Remove the connected component of visible triangles. Save the
            // terminator edges for insertion of the new visible set of
            // triangles.
            const terminator: [number, number][] = [];
            while (visibleHead < visible.length) {
                const tri = visible[visibleHead++];
                for (let i = 0; i < 3; ++i) {
                    const adj = tri.T[i];
                    if (adj !== null) {
                        if (this.toPlane(adj.V[0], adj.V[1], adj.V[2], h1) <= 0) {
                            // The shared edge of tri and adj is a terminator.
                            terminator.push([tri.V[i], tri.V[(i + 1) % 3]]);
                        }
                        else {
                            if (!visited.has(adj)) {
                                visible.push(adj);
                                visited.add(adj);
                            }
                        }
                    }
                }
                visited.delete(tri);
                const removed = hullMesh.remove(tri.V[0], tri.V[1], tri.V[2]);
                logAssert(removed, 'Unexpected removal failure.');
            }

            // Insert the new hull triangles.
            for (const edge of terminator) {
                const inserted = hullMesh.insert(edge[0], edge[1], h1);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }

            // The current index h1 becomes the previous index h0 for the next
            // pass of the 'current' loop.
            h0 = h1;
        }
    }

    // Memoized access to the rational representation of the points.
    private getRationalPoint(index: number): RationalPoint3 {
        let rPoint = this.mRPoints[index];
        if (rPoint === null) {
            const values = this.mPoints[index].values;
            rPoint = [BSNumber.fromNumber(values[0]), BSNumber.fromNumber(values[1]),
                BSNumber.fromNumber(values[2])];
            this.mRPoints[index] = rPoint;
        }
        return rPoint;
    }

    private colocated(v0: number, v1: number): boolean {
        const r0 = this.getRationalPoint(v0);
        const r1 = this.getRationalPoint(v1);
        return r0[0].equals(r1[0]) && r0[1].equals(r1[1]) && r0[2].equals(r1[2]);
    }

    private colinear(v0: number, v1: number, v2: number): boolean {
        const rvec0 = this.getRationalPoint(v0);
        const rvec1 = this.getRationalPoint(v1);
        const rvec2 = this.getRationalPoint(v2);
        const rdiff1: RationalPoint3 =
            [rvec1[0].sub(rvec0[0]), rvec1[1].sub(rvec0[1]), rvec1[2].sub(rvec0[2])];
        const rdiff2: RationalPoint3 =
            [rvec2[0].sub(rvec0[0]), rvec2[1].sub(rvec0[1]), rvec2[2].sub(rvec0[2])];
        const rcross = rationalCross(rdiff1, rdiff2);
        return rcross[0].getSign() === 0
            && rcross[1].getSign() === 0
            && rcross[2].getSign() === 0;
    }

    // For a plane with origin V0 and normal N = Cross(V1-V0,V2-V0), toPlane
    // returns
    //   +1, V3 on positive side of plane (side to which N points)
    //   -1, V3 on negative side of plane (side to which -N points)
    //    0, V3 on the plane
    private toPlane(v0: number, v1: number, v2: number, v3: number): number {
        // Attempt to classify the sign using interval arithmetic.
        const p0 = this.mPoints[v0].values;
        const p1 = this.mPoints[v1].values;
        const p2 = this.mPoints[v2].values;
        const p3 = this.mPoints[v3].values;

        const sDiff1: SWInterval[] = [
            SWInterval.sub(p1[0], p0[0]),
            SWInterval.sub(p1[1], p0[1]),
            SWInterval.sub(p1[2], p0[2])
        ];
        const sDiff2: SWInterval[] = [
            SWInterval.sub(p2[0], p0[0]),
            SWInterval.sub(p2[1], p0[1]),
            SWInterval.sub(p2[2], p0[2])
        ];
        const sDiff3: SWInterval[] = [
            SWInterval.sub(p3[0], p0[0]),
            SWInterval.sub(p3[1], p0[1]),
            SWInterval.sub(p3[2], p0[2])
        ];
        const sDet = intervalDotCross(sDiff1, sDiff2, sDiff3);
        if (sDet.get(0) > 0) {
            return +1;
        }
        if (sDet.get(1) < 0) {
            return -1;
        }

        // The sign is indeterminate using interval arithmetic.
        const r0 = this.getRationalPoint(v0);
        const r1 = this.getRationalPoint(v1);
        const r2 = this.getRationalPoint(v2);
        const r3 = this.getRationalPoint(v3);
        const rDiff1: RationalPoint3 =
            [r1[0].sub(r0[0]), r1[1].sub(r0[1]), r1[2].sub(r0[2])];
        const rDiff2: RationalPoint3 =
            [r2[0].sub(r0[0]), r2[1].sub(r0[1]), r2[2].sub(r0[2])];
        const rDiff3: RationalPoint3 =
            [r3[0].sub(r0[0]), r3[1].sub(r0[1]), r3[2].sub(r0[2])];
        const rDet = rationalDotCross(rDiff1, rDiff2, rDiff3);
        return rDet.getSign();
    }
}
