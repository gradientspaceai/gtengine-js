// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OBBTreeOfTriangles.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in OBBTree.ts regarding tree construction.
//
// Port notes:
//   - Upstream 'Create(vertices, triangles, height)' hides the base class
//     'Create(centroids, height)' by C++ name hiding. TypeScript requires an
//     override to be assignable to the base member, so the triangle-specific
//     creation function is named 'createFromTriangles', following the
//     'createFromSegments' precedent in OBBTreeOfSegments.ts. The inherited
//     'create(centroids, height)' remains callable but is not meaningful for
//     a triangle tree.
//   - 'std::array<size_t, 3>' becomes the tuple type
//     [number, number, number].
//   - The upstream overrides ComputeInteriorBox and ComputeLeafBox are
//     'private virtual'; they are 'protected override' here because
//     TypeScript requires an override to be at least as visible as the base
//     member.
//   - The nested struct Intersection becomes the top-level exported class
//     OBBTreeOfTrianglesIntersection (a nested class is not expressible in
//     TypeScript); see the BVTreeOfTrianglesIntersection precedent.
//   - The private static box and triangle query functions are module-private
//     functions here. They are deliberately not exported: BVTreeOfTriangles.ts
//     already owns the global names 'intersectLineTriangle',
//     'intersectRayTriangle' and 'intersectSegmentTriangle'.
//   - 'Execute' becomes 'execute' and returns the sorted intersection array
//     rather than filling an output container.
//   - Upstream initializes pmin and pmax in ComputeInteriorBox to the zero
//     vector rather than to the projections of the first primitive vertex.
//     Combined with the 'if (dot < pmin[j]) ... else if (dot > pmax[j]) ...'
//     update, every projection is still bracketed by [pmin[j], pmax[j]], so
//     the box does contain all triangle vertices; it is merely conservative.
//     The quirk is preserved; see upstream issue #103 (OBBTreeOfSegments has
//     the identical code).
//   - Upstream collects the hits in a 'std::set<Intersection>' whose ordering
//     'Intersection::operator<' compares the parameter alone, so two
//     triangles hit at the same parameter (a linear component through a
//     shared edge or vertex of a closed mesh — the routine case) are
//     set-equivalent and all but one are silently dropped, with the survivor
//     depending on tree partition/traversal order. That corrupts results, so
//     the port FIXES it: hits are ordered lexicographically by
//     (parameter, triangleIndex) and every hit is kept. See upstream issue
//     #167, filed for the identical bug in BVTreeOfTriangles.h, whose port is
//     fixed the same way.

import { IntrLine3OrientedBox3TI } from './IntrLine3OrientedBox3.js';
import { IntrLine3Triangle3FI } from './IntrLine3Triangle3.js';
import { IntrRay3OrientedBox3TI } from './IntrRay3OrientedBox3.js';
import { IntrRay3Triangle3FI } from './IntrRay3Triangle3.js';
import { IntrSegment3OrientedBox3TI } from './IntrSegment3OrientedBox3.js';
import { IntrSegment3Triangle3FI } from './IntrSegment3Triangle3.js';
import { Line } from './Line.js';
import { logAssert } from './Logger.js';
import { OBBNode, OBBTree } from './OBBTree.js';
import { OrientedBox } from './OrientedBox.js';
import { Ray } from './Ray.js';
import { Segment } from './Segment.js';
import { Triangle } from './Triangle.js';
import { Vector, add, div, dot, length, normalize, sub } from './Vector.js';
import { cross, unitCross } from './Vector3.js';

// The port of the upstream function-pointer type BoxQuery.
type BoxQuery = (P: Vector, Q: Vector, box: OrientedBox) => boolean;

// The line is P + t * Q for all real t.
function intersectLineOrientedBox(P: Vector, Q: Vector, box: OrientedBox): boolean {
    const query = new IntrLine3OrientedBox3TI();
    return query.test(Line.fromOriginDirection(P, Q), box).intersect;
}

// The ray is P + t * Q for t >= 0.
function intersectRayOrientedBox(P: Vector, Q: Vector, box: OrientedBox): boolean {
    const query = new IntrRay3OrientedBox3TI();
    return query.test(Ray.fromOriginDirection(P, Q), box).intersect;
}

// The segment has endpoints P and Q.
function intersectSegmentOrientedBox(P: Vector, Q: Vector, box: OrientedBox): boolean {
    const query = new IntrSegment3OrientedBox3TI();
    return query.test(Segment.fromEndpoints(P, Q), box).intersect;
}

// The port of the upstream private struct TriangleResult.
interface TriangleHitResult {
    intersect: boolean;
    point: Vector;
    parameter: number;
}

// The port of the upstream function-pointer type TriangleQuery.
type TriangleQuery = (P: Vector, Q: Vector, triangle: Triangle) => TriangleHitResult;

function intersectLineTri(P: Vector, Q: Vector, triangle: Triangle): TriangleHitResult {
    const query = new IntrLine3Triangle3FI();
    const result = query.find(Line.fromOriginDirection(P, Q), triangle);
    return {
        intersect: result.intersect,
        point: result.point,
        parameter: result.parameter
    };
}

function intersectRayTri(P: Vector, Q: Vector, triangle: Triangle): TriangleHitResult {
    const query = new IntrRay3Triangle3FI();
    const result = query.find(Ray.fromOriginDirection(P, Q), triangle);
    return {
        intersect: result.intersect,
        point: result.point,
        parameter: result.parameter
    };
}

function intersectSegmentTri(P: Vector, Q: Vector, triangle: Triangle): TriangleHitResult {
    const query = new IntrSegment3Triangle3FI();
    const result = query.find(Segment.fromEndpoints(P, Q), triangle);

    // The segment is converted to centered form in the query. That form is
    // C + s * D, where C is the midpoint of the segment, D is a unit-length
    // vector and |s| <= e for segment extent (half length) e. The
    // t-parameter must be converted back to (1-t)*P + t*Q where t in [0,1].
    // Thus, t = (s+e)/(2*e), which is equivalent to s/Length(Q-P) + 1/2.
    return {
        intersect: result.intersect,
        point: result.point,
        parameter: result.parameter / length(sub(Q, P)) + 0.5
    };
}

// The port of the nested struct OBBTreeOfTriangles<T>::Intersection.
export class OBBTreeOfTrianglesIntersection {
    // The port of std::numeric_limits<size_t>::max() used by the default
    // constructor; see the OBBNode.invalid precedent in OBBTree.ts.
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    triangleIndex: number;
    point: Vector;
    parameter: number;

    constructor(triangleIndex: number = OBBTreeOfTrianglesIntersection.invalid,
        point: Vector = Vector.zero(3), parameter: number = 0) {
        this.triangleIndex = triangleIndex;
        this.point = point;
        this.parameter = parameter;
    }

    // The port of 'operator<'. Upstream compares the parameters alone, which
    // makes coincident hits set-equivalent; the port breaks the tie with the
    // triangle index so that no hit is dropped. See the port notes.
    lessThan(other: OBBTreeOfTrianglesIntersection): boolean {
        if (this.parameter < other.parameter) {
            return true;
        }
        if (other.parameter < this.parameter) {
            return false;
        }
        return this.triangleIndex < other.triangleIndex;
    }
}

export class OBBTreeOfTriangles extends OBBTree {
    // These are the queryType inputs to execute(). LINE_QUERY and RAY_QUERY
    // take P as the origin and Q as a unit-length direction; SEGMENT_QUERY
    // takes P and Q as the segment endpoints.
    static readonly LINE_QUERY: number = 0;
    static readonly RAY_QUERY: number = 1;
    static readonly SEGMENT_QUERY: number = 2;

    protected mVertices: Vector[];
    protected mTriangles: [number, number, number][];
    protected mBoxQueries: BoxQuery[];
    protected mTriangleQueries: TriangleQuery[];

    constructor() {
        super();
        this.mVertices = [];
        this.mTriangles = [];
        this.mBoxQueries = [
            intersectLineOrientedBox,
            intersectRayOrientedBox,
            intersectSegmentOrientedBox
        ];
        this.mTriangleQueries = [
            intersectLineTri,
            intersectRayTri,
            intersectSegmentTri
        ];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If OBBTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from triangles.length. If
    // larger than 31, the height is clamped to 31.
    createFromTriangles(vertices: readonly Vector[],
        triangles: readonly (readonly [number, number, number])[],
        height: number = OBBTree.fullHeight): void {
        logAssert(vertices.length >= 3 && triangles.length > 0, 'Invalid input.');

        this.mVertices = new Array<Vector>(vertices.length);
        for (let i = 0; i < vertices.length; ++i) {
            this.mVertices[i] = vertices[i].clone();
        }
        this.mTriangles = new Array<[number, number, number]>(triangles.length);
        for (let i = 0; i < triangles.length; ++i) {
            this.mTriangles[i] = [triangles[i][0], triangles[i][1], triangles[i][2]];
        }

        // Compute the triangle centroids.
        const centroids = new Array<Vector>(this.mTriangles.length);
        const three = 3;
        for (let t = 0; t < this.mTriangles.length; ++t) {
            const tri = this.mTriangles[t];
            centroids[t] = div(add(add(this.mVertices[tri[0]],
                this.mVertices[tri[1]]), this.mVertices[tri[2]]), three);
        }

        // Create the OBB tree for centroids.
        super.create(centroids, height);
    }

    // Member access.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getTriangles(): readonly (readonly [number, number, number])[] {
        return this.mTriangles;
    }

    // Let C be the box center and let U0, U1 and U2 be the box axes. Each
    // input point is of the form X = C + y0*U0 + y1*U1 + y2*U2. The following
    // code computes min(y0), max(y0), min(y1), max(y1), min(y2) and max(y2).
    // The box center is then adjusted to be
    //   C' = C + 0.5*(min(y0)+max(y0))*U0 + 0.5*(min(y1)+max(y1))*U1
    //        + 0.5*(min(y2)+max(y2))*U2
    protected override computeInteriorBox(i0: number, i1: number, box: OrientedBox): void {
        super.computeInteriorBox(i0, i1, box);

        const pmin = Vector.zero(3);
        const pmax = Vector.zero(3);
        for (let i = i0; i <= i1; ++i) {
            const tri = this.mTriangles[this.mPartition[i]];
            for (let k = 0; k < 3; ++k) {
                const diff = sub(this.mVertices[tri[k]], box.center);
                for (let j = 0; j < 3; ++j) {
                    const d = dot(diff, box.axis[j]);
                    if (d < pmin.values[j]) {
                        pmin.values[j] = d;
                    } else if (d > pmax.values[j]) {
                        pmax.values[j] = d;
                    }
                }
            }
        }

        const half = 0.5;
        for (let j = 0; j < 3; ++j) {
            const scalar = half * (pmin.values[j] + pmax.values[j]);
            const axis = box.axis[j];
            for (let k = 0; k < 3; ++k) {
                box.center.values[k] += scalar * axis.values[k];
            }
            box.extent.values[j] = half * (pmax.values[j] - pmin.values[j]);
        }
    }

    protected override computeLeafBox(i: number, box: OrientedBox): void {
        // Create a degenerate box whose center is the midpoint of the
        // triangle primitive, whose axis[0] is the direction of a triangle
        // edge, whose axis[2] is a triangle normal, and whose axis[1] is
        // Cross(axis[2], axis[0]). The extent[0] and extent[1] are chosen so
        // that the box contains the triangle. The extent[2] is zero.
        const tri = this.mTriangles[this.mPartition[i]];
        const edge10 = sub(this.mVertices[tri[1]], this.mVertices[tri[0]]);
        const edge20 = sub(this.mVertices[tri[2]], this.mVertices[tri[0]]);
        normalize(edge10);
        normalize(edge20);
        const normal = unitCross(edge10, edge20);

        box.center = this.mCentroids[this.mPartition[i]].clone();
        box.axis[0] = edge10;
        box.axis[1] = cross(normal, edge10);
        box.axis[2] = normal;

        const V0mC = sub(this.mVertices[tri[0]], box.center);
        const V1mC = sub(this.mVertices[tri[1]], box.center);
        const V2mC = sub(this.mVertices[tri[2]], box.center);
        const ax0 = Math.abs(dot(box.axis[0], V0mC));
        const ax1 = Math.abs(dot(box.axis[0], V1mC));
        const ax2 = Math.abs(dot(box.axis[0], V2mC));
        const ay0 = Math.abs(dot(box.axis[1], V0mC));
        const ay1 = Math.abs(dot(box.axis[1], V1mC));
        const ay2 = Math.abs(dot(box.axis[1], V2mC));
        box.extent.values[0] = Math.max(ax0, Math.max(ax1, ax2));
        box.extent.values[1] = Math.max(ay0, Math.max(ay1, ay2));
        box.extent.values[2] = 0;
    }

    // Generate a list of triangles intersected by a linear component (line,
    // ray or segment). The line is parameterized by P + t * Q, where Q is a
    // unit-length direction and t is any real number. The ray is
    // parameterized by P + t * Q, where Q is a unit-length direction and
    // t >= 0. The segment is parameterized by
    // (1-t) * P + t * Q = P + t * (Q - P), where P and Q are the endpoints of
    // the segment and 0 <= t <= 1.
    //
    // The intersections are sorted by the parameter, ties broken by the
    // triangle index. See the port notes about the upstream std::set.
    execute(queryType: number, P: Vector, Q: Vector): OBBTreeOfTrianglesIntersection[] {
        const invalid = OBBNode.invalid;
        const intersections: OBBTreeOfTrianglesIntersection[] = [];
        const boxQuery = this.mBoxQueries[queryType];
        const triangleQuery = this.mTriangleQueries[queryType];

        const indexStack = new Array<number>(2 * this.mHeight + 1).fill(0);
        let top = 0;
        indexStack[0] = 0;
        while (top >= 0) {
            const nodeIndex = indexStack[top--];
            const node = this.mNodes[nodeIndex];

            // For the balanced tree created by OBBTree, an interior node has
            // two valid children and a leaf node has two invalid children.
            // This is true even if the height passed to OBBTree.create is
            // smaller than the actual height.
            if (node.leftChild !== invalid && node.rightChild !== invalid) {
                // The node is interior.
                if (boxQuery(P, Q, node.box)) {
                    // The linear component intersects the box. Continue the
                    // intersection search to the child nodes.
                    indexStack[++top] = node.rightChild;
                    indexStack[++top] = node.leftChild;
                }
                // Otherwise the linear component does not intersect the box.
                // There are no triangles intersected in the subtree rooted at
                // this node, so do not continue the search to child nodes.
            } else {
                // The node is a leaf.
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const triangleIndex = this.mPartition[i];
                    const tri = this.mTriangles[triangleIndex];
                    const triangle = Triangle.fromVertices(this.mVertices[tri[0]],
                        this.mVertices[tri[1]], this.mVertices[tri[2]]);
                    const triResult = triangleQuery(P, Q, triangle);
                    if (triResult.intersect) {
                        intersections.push(new OBBTreeOfTrianglesIntersection(
                            triangleIndex, triResult.point, triResult.parameter));
                    }
                }
            }
        }

        // The leaf index ranges are disjoint and mPartition is a permutation,
        // so each triangle is tested at most once and no two hits share a
        // triangle index. Sorting by (parameter, triangleIndex) therefore
        // yields the same container as inserting into an ordered set with
        // OBBTreeOfTrianglesIntersection.lessThan.
        intersections.sort((a, b) => (a.lessThan(b) ? -1 : (b.lessThan(a) ? +1 : 0)));
        return intersections;
    }
}
