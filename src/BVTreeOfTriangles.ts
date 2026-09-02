// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BVTreeOfTriangles.h
// Upstream: Geometric Tools Library, https://www.geometrictools.com
// Copyright (c) 2025 Geometric Tools LLC
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Read the comments in BVTree.ts regarding tree construction. Although this
// class appears to be non-abstract upstream, the BoundingVolume type has
// requirements for its interface. In this sense, BVTreeOfTriangles is
// abstract.
//
// Port notes:
//   - See BVTreeOfPoints.ts for the BoundingVolume type-parameter split and
//     for why the class is declared 'abstract'.
//   - Upstream 'Create(vertices, triangles, height)' hides the base class
//     'Create(centroids, height)' by C++ name hiding. TypeScript requires an
//     override to be assignable to the base member, so the triangle-specific
//     creation function is named 'createFromTriangles'.
//   - 'std::array<std::size_t, 3>' becomes the tuple type
//     [number, number, number].
//   - The nested class Intersection becomes the top-level exported class
//     BVTreeOfTrianglesIntersection (a nested class is not expressible in
//     TypeScript); see the BVTreeNode precedent in BVTree.ts.
//   - 'Execute' becomes 'execute' and returns { nodeIndices, intersections }
//     rather than filling output containers. The upstream output container is
//     std::set<Intersection>, whose ordering is by the linear-component
//     parameter alone, so upstream silently drops all but one of the hits
//     that share a parameter (a line through a shared edge or vertex) -- the
//     result-corrupting bug tracked as upstream issue #167. The port instead
//     orders lexicographically by (parameter, triangleIndex) and keeps every
//     hit, matching the OBBTreeOfTriangles port.

import { BVTree } from './BVTree';
import type { BVTreeBoundingVolume, BVTreeVolumeOps } from './BVTree';
import { IntrLine3Triangle3FI } from './IntrLine3Triangle3';
import { IntrRay3Triangle3FI } from './IntrRay3Triangle3';
import { IntrSegment3Triangle3FI } from './IntrSegment3Triangle3';
import { Line } from './Line';
import { logAssert } from './Logger';
import { Ray } from './Ray';
import { Segment } from './Segment';
import { Triangle } from './Triangle';
import { Vector, add, div } from './Vector';

// The result of a linear-component/triangle intersection: whether there is an
// intersection, the intersection point and the parameter of the linear
// component at that point.
export interface LinearTriangleResult {
    intersect: boolean;
    point: Vector;
    parameter: number;
}

// The port of the upstream function-pointer type LinearTriangleQuery.
export type LinearTriangleQuery =
    (P: Vector, Q: Vector, triangle: Triangle) => LinearTriangleResult;

// The line is P + t * Q for all real t.
export function intersectLineTriangle(P: Vector, Q: Vector,
    triangle: Triangle): LinearTriangleResult {
    const query = new IntrLine3Triangle3FI();
    const output = query.find(Line.fromOriginDirection(P, Q), triangle);
    return {
        intersect: output.intersect,
        point: output.point,
        parameter: output.parameter
    };
}

// The ray is P + t * Q for t >= 0.
export function intersectRayTriangle(P: Vector, Q: Vector,
    triangle: Triangle): LinearTriangleResult {
    const query = new IntrRay3Triangle3FI();
    const output = query.find(Ray.fromOriginDirection(P, Q), triangle);
    return {
        intersect: output.intersect,
        point: output.point,
        parameter: output.parameter
    };
}

// The segment has endpoints P and Q.
export function intersectSegmentTriangle(P: Vector, Q: Vector,
    triangle: Triangle): LinearTriangleResult {
    const query = new IntrSegment3Triangle3FI();
    const output = query.find(Segment.fromEndpoints(P, Q), triangle);
    return {
        intersect: output.intersect,
        point: output.point,
        parameter: output.parameter
    };
}

// The port of the nested class
// BVTreeOfTriangles<T, BoundingVolume>::Intersection.
export class BVTreeOfTrianglesIntersection {
    // The port of std::numeric_limits<std::size_t>::max() used by the default
    // constructor; see the same precedent in BVTree.ts.
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    triangleIndex: number;
    point: Vector;
    parameter: number;

    constructor(triangleIndex: number = BVTreeOfTrianglesIntersection.invalid,
        point: Vector = Vector.zero(3), parameter: number = 0) {
        this.triangleIndex = triangleIndex;
        this.point = point;
        this.parameter = parameter;
    }

    // Upstream 'operator<' compares only the parameters, which makes two
    // hits at the same parameter set-equivalent and drops one (upstream
    // issue #167). The port breaks ties by triangleIndex so every hit is
    // kept.
    lessThan(other: BVTreeOfTrianglesIntersection): boolean {
        return this.parameter < other.parameter ||
            (this.parameter === other.parameter &&
                this.triangleIndex < other.triangleIndex);
    }
}

// The result of BVTreeOfTriangles.execute; the port of the two upstream
// output parameters.
export interface BVTreeOfTrianglesExecuteResult {
    nodeIndices: number[];
    intersections: BVTreeOfTrianglesIntersection[];
}

// The port of std::set<Intersection>::insert with the corrected strict weak
// ordering (parameter, then triangleIndex) -- see the header comment and
// upstream issue #167. Only an exact duplicate (same parameter and same
// triangle) is rejected, which cannot occur in execute() since each triangle
// is tested once.
function insertIntersection(intersections: BVTreeOfTrianglesIntersection[],
    item: BVTreeOfTrianglesIntersection): void {
    let lo = 0;
    let hi = intersections.length;
    while (lo < hi) {
        const mid = lo + Math.floor((hi - lo) / 2);
        if (intersections[mid].lessThan(item)) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if (lo < intersections.length && !item.lessThan(intersections[lo])) {
        // An equivalent element is already a member of the set.
        return;
    }
    intersections.splice(lo, 0, item);
}

export abstract class BVTreeOfTriangles<BV extends BVTreeBoundingVolume>
    extends BVTree<BV> {
    protected mVertices: Vector[];
    protected mTriangles: [number, number, number][];
    protected mLinearTriangleQuery: LinearTriangleQuery[];

    protected constructor(ops: BVTreeVolumeOps<BV>) {
        super(ops);
        this.mVertices = [];
        this.mTriangles = [];
        this.mLinearTriangleQuery = [
            intersectLineTriangle,
            intersectRayTriangle,
            intersectSegmentTriangle
        ];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If BVTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from triangles.length. If
    // larger than 31, the height is clamped to 31.
    createFromTriangles(vertices: readonly Vector[],
        triangles: readonly (readonly [number, number, number])[],
        height: number = BVTree.fullHeight): void {
        logAssert(vertices.length >= 3 && triangles.length > 0,
            'Expecting at least 3 vertices and at least 1 triangle.');

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

        // Create the bounding volume tree for centroids.
        super.create(centroids, height);
    }

    // Member access.
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getTriangles(): readonly (readonly [number, number, number])[] {
        return this.mTriangles;
    }

    // Compute intersections of the linear component and triangles. These are
    // sorted by the parameter of the linear component.
    execute(queryType: number, P: Vector, Q: Vector): BVTreeOfTrianglesExecuteResult {
        const nodeIndices = this.getLeafIndices(queryType, P, Q);

        const linearTriangleQuery = this.mLinearTriangleQuery[queryType];
        const intersections: BVTreeOfTrianglesIntersection[] = [];
        for (const leafIndex of nodeIndices) {
            const node = this.mNodes[leafIndex];
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const triangleIndex = this.mPartition[i];
                const tri = this.mTriangles[triangleIndex];
                const triangle = Triangle.fromVertices(this.mVertices[tri[0]],
                    this.mVertices[tri[1]], this.mVertices[tri[2]]);
                const output = linearTriangleQuery(P, Q, triangle);
                if (output.intersect) {
                    insertIntersection(intersections,
                        new BVTreeOfTrianglesIntersection(triangleIndex,
                            output.point, output.parameter));
                }
            }
        }

        return { nodeIndices: nodeIndices, intersections: intersections };
    }
}
