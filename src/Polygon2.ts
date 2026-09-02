// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Polygon2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Polygon2 object represents a simple polygon: no duplicate vertices,
// closed (each vertex is shared by exactly two edges), and no
// self-intersections at interior edge points. The 'vertexPool' array can
// contain more points than needed to define the polygon, which allows the
// vertex pool to have multiple polygons associated with it. Thus, the
// programmer must ensure that the vertex pool persists as long as any
// Polygon2 objects exist that depend on the pool. The number of polygon
// vertices is 'indices.length' and must be 3 or larger. The 'indices' array
// refers to the points in 'vertexPool' that are part of the polygon and must
// have unique elements. The edges of the polygon are pairs of indices into
// 'vertexPool',
//   edge[0] = (indices[0], indices[1])
//   :
//   edge[numIndices-2] = (indices[numIndices-2], indices[numIndices-1])
//   edge[numIndices-1] = (indices[numIndices-1], indices[0])
// The programmer should ensure the polygon is simple. The geometric queries
// are valid regardless of whether the polygon is oriented clockwise or
// counterclockwise.
//
// NOTE: Comparison operators are not provided. The semantics of equal
// polygons is complicated and (at the moment) not useful. The vertex pools
// can be different and indices do not match, but the vertices they reference
// can match. Even with a shared vertex pool, the indices can be "rotated",
// leading to the same polygon abstractly but the data structures do not
// match.
//
// Port notes: the upstream (vertexPool, numIndices, indices) triple becomes
// (vertexPool, indices), with numIndices = indices.length. The vertex pool is
// referenced, not copied, matching the upstream pointer. The upstream
// 'operator bool()' becomes isValid(). The upstream std::set<int32_t> of
// vertex indices becomes a sorted number[], so iteration order (which affects
// the floating-point summation of computeVertexAverage) matches upstream.

import { IntrSegment2Segment2TI } from './IntrSegment2Segment2';
import { Segment } from './Segment';
import { Vector, add, div, length, sub } from './Vector';
import { dotPerp } from './Vector2';

export class Polygon2 {
    private mVertexPool: readonly Vector[] | null;
    private mVertices: number[];
    private mIndices: number[];
    private mCounterClockwise: boolean;

    // Construction. The constructor succeeds when 'indices.length' >= 3, the
    // 'vertexPool' is nonempty, and the indices are unique; we cannot test
    // whether you have a valid number of elements in the input arrays. A copy
    // is made of 'indices', but the 'vertexPool' is not copied. If the
    // constructor fails, the internal vertex pool is set to null, the index
    // array has no elements, and the orientation is set to clockwise.
    constructor(vertexPool: readonly Vector[] | null,
        indices: readonly number[] | null, counterClockwise: boolean) {
        this.mVertexPool = vertexPool;
        this.mVertices = [];
        this.mIndices = [];
        this.mCounterClockwise = counterClockwise;

        const numIndices = (indices !== null ? indices.length : 0);
        if (numIndices >= 3 && vertexPool !== null && vertexPool.length > 0 &&
            indices !== null) {
            const unique = new Set<number>();
            for (let i = 0; i < numIndices; ++i) {
                unique.add(indices[i]);
            }

            if (numIndices === unique.size) {
                this.mVertices = Array.from(unique).sort((a, b) => a - b);
                this.mIndices = indices.slice();
                return;
            }

            // At least one duplicated vertex was encountered, so the polygon
            // is not simple. Fail the constructor call.
        }

        // Invalid input to the Polygon2 constructor.
        this.mVertexPool = null;
        this.mCounterClockwise = false;
    }

    // To validate construction, create an object as shown:
    //     const polygon = new Polygon2(parameters);
    //     if (!polygon.isValid()) { /* constructor failed */ }
    isValid(): boolean {
        return this.mVertexPool !== null;
    }

    // Member access.
    getVertexPool(): readonly Vector[] | null {
        return this.mVertexPool;
    }

    // The sorted, unique vertex indices used by the polygon. The returned
    // array must not be modified.
    getVertices(): readonly number[] {
        return this.mVertices;
    }

    // The polygon indices, in polygon order. The returned array must not be
    // modified.
    getIndices(): readonly number[] {
        return this.mIndices;
    }

    counterClockwise(): boolean {
        return this.mCounterClockwise;
    }

    // Geometric queries.
    computeVertexAverage(): Vector {
        let average = new Vector(2);
        if (this.mVertexPool !== null) {
            for (const index of this.mVertices) {
                average = add(average, this.mVertexPool[index]);
            }
            average = div(average, this.mVertices.length);
        }
        return average;
    }

    computePerimeterLength(): number {
        let perimeter = 0;
        if (this.mVertexPool !== null) {
            const pool = this.mVertexPool;
            let v0 = pool[this.mIndices[this.mIndices.length - 1]];
            for (const index of this.mIndices) {
                const v1 = pool[index];
                perimeter += length(sub(v1, v0));
                v0 = v1;
            }
        }
        return perimeter;
    }

    computeArea(): number {
        let area = 0;
        if (this.mVertexPool !== null) {
            const pool = this.mVertexPool;
            const numIndices = this.mIndices.length;
            let v0 = pool[this.mIndices[numIndices - 2]];
            let v1 = pool[this.mIndices[numIndices - 1]];
            for (const index of this.mIndices) {
                const v2 = pool[index];
                area += v1.values[0] * (v2.values[1] - v0.values[1]);
                v0 = v1;
                v1 = v2;
            }
            area *= 0.5;
        }
        return Math.abs(area);
    }

    // Simple polygons have no self-intersections at interior points of edges.
    // The test is an exhaustive all-pairs intersection test for edges, which
    // is inefficient for polygons with a large number of vertices.
    // TODO (upstream): Provide an efficient algorithm that uses the algorithm
    // of class RectangleManager.h.
    isSimple(): boolean {
        if (this.mVertexPool === null) {
            return false;
        }

        // For mVertexPool to be nonnull, the number of indices is guaranteed
        // to be at least 3.
        if (this.mIndices.length === 3) {
            // The polygon is a triangle.
            return true;
        }

        return this.isSimpleInternal();
    }

    // Convex polygons are simple polygons where the angles between
    // consecutive edges are less than or equal to pi radians.
    isConvex(): boolean {
        if (this.mVertexPool === null) {
            return false;
        }

        // For mVertexPool to be nonnull, the number of indices is guaranteed
        // to be at least 3.
        if (this.mIndices.length === 3) {
            // The polygon is a triangle.
            return true;
        }

        return this.isSimpleInternal() && this.isConvexInternal();
    }

    // These calls have preconditions that mVertexPool is not null and
    // mIndices.length > 3. The heart of the algorithms are implemented here.
    private isSimpleInternal(): boolean {
        const pool = this.mVertexPool as readonly Vector[];
        const query = new IntrSegment2Segment2TI();

        const numIndices = this.mIndices.length;
        for (let i0 = 0; i0 < numIndices; ++i0) {
            const i0p1 = (i0 + 1) % numIndices;
            const seg0 = Segment.fromEndpoints(pool[this.mIndices[i0]],
                pool[this.mIndices[i0p1]]);

            const i1min = (i0 + 2) % numIndices;
            const i1max = (i0 - 2 + numIndices) % numIndices;
            for (let i1 = i1min; i1 <= i1max; ++i1) {
                const i1p1 = (i1 + 1) % numIndices;
                const seg1 = Segment.fromEndpoints(pool[this.mIndices[i1]],
                    pool[this.mIndices[i1p1]]);

                const result = query.test(seg0, seg1);
                if (result.intersect) {
                    return false;
                }
            }
        }
        return true;
    }

    private isConvexInternal(): boolean {
        const pool = this.mVertexPool as readonly Vector[];
        const sign = (this.mCounterClockwise ? 1 : -1);
        const numIndices = this.mIndices.length;
        for (let i = 0; i < numIndices; ++i) {
            const iPrev = (i + numIndices - 1) % numIndices;
            const iNext = (i + 1) % numIndices;
            const vPrev = pool[this.mIndices[iPrev]];
            const vCurr = pool[this.mIndices[i]];
            const vNext = pool[this.mIndices[iNext]];
            const edge0 = sub(vCurr, vPrev);
            const edge1 = sub(vNext, vCurr);
            const test = sign * dotPerp(edge0, edge1);
            if (test < 0) {
                return false;
            }
        }
        return true;
    }
}
