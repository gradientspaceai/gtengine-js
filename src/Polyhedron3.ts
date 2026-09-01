// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Polyhedron3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Polyhedron3 object represents a simple polyhedron. The 'vertexPool'
// array can contain more points than needed to define the polyhedron, which
// allows the vertex pool to have multiple polyhedra associated with it. The
// number of polyhedron indices is 'numIndices' and must be 12 or larger (at
// least four triangles). The 'indices' array refers to the points in
// 'vertexPool' that form the triangle faces, so 'numIndices' must be a
// multiple of 3. The number of vertices is the number of unique elements in
// 'indices' and is determined during construction. The programmer should
// ensure the polyhedron is simple. The geometric queries are valid regardless
// of whether the polyhedron triangles are oriented clockwise or
// counterclockwise.
//
// NOTE: Comparison operators are not provided. The semantics of equal
// polyhedra is complicated and (at the moment) not useful. The vertex pools
// can be different and indices do not match, but the vertices they reference
// can match. Even with a shared vertex pool, the indices can be permuted,
// leading to the same polyhedron abstractly but the data structures do not
// match.
//
// Port notes: the C++ 'std::shared_ptr<std::vector<Vector3<Real>>>' vertex
// pool becomes a plain array reference (null when construction failed); the
// pool is shared, not copied, as upstream. The C++ 'std::set<int32_t>' of
// unique indices becomes a sorted number[] so that iteration order (which
// affects floating-point summation in computeVertexAverage) matches
// upstream. The 'operator bool' validity check becomes 'isValid()'.

import { Vector, add, sub, div, length } from './Vector';
import { dotCross, cross } from './Vector3';

export class Polyhedron3 {
    private mVertexPool: Vector[] | null;
    private mUniqueIndices: number[];
    private mIndices: number[];
    private mCounterClockwise: boolean;

    // Construction. The constructor succeeds when 'numIndices >= 12' (at
    // least 4 triangles), 'numIndices' is a multiple of 3, and 'vertexPool'
    // and 'indices' are not null. A copy is made of 'indices', but the
    // 'vertexPool' is not copied. If the constructor fails, the internal
    // vertex pool is set to null, the number of vertices is set to zero, the
    // index array has no elements, and the triangle face orientation is set
    // to clockwise.
    //
    // Port deviation: upstream cannot test whether 'indices' has at least
    // 'numIndices' elements (it is a raw pointer). Here the array length is
    // known, so a too-short array is also treated as an invalid input rather
    // than reading past the end.
    constructor(vertexPool: Vector[] | null, numIndices: number,
        indices: readonly number[] | null, counterClockwise: boolean) {
        this.mVertexPool = vertexPool;
        this.mCounterClockwise = counterClockwise;
        this.mUniqueIndices = [];
        this.mIndices = [];

        if (vertexPool !== null && indices !== null && numIndices >= 12
            && (numIndices % 3) === 0 && indices.length >= numIndices) {
            const unique = new Set<number>();
            for (let i = 0; i < numIndices; ++i) {
                unique.add(indices[i]);
            }
            // std::set iterates in ascending order; replicate that ordering.
            this.mUniqueIndices = Array.from(unique).sort((a, b) => a - b);
            this.mIndices = indices.slice(0, numIndices);
        } else {
            // Encountered an invalid input.
            this.mVertexPool = null;
            this.mCounterClockwise = false;
        }
    }

    // To validate construction, create an object and check isValid().
    isValid(): boolean {
        return this.mVertexPool !== null;
    }

    // Member access.
    getVertexPool(): Vector[] | null {
        return this.mVertexPool;
    }

    // The port of 'GetVertices', which dereferences the vertex pool. It
    // throws when the polyhedron is invalid (upstream has undefined behavior
    // in that case).
    getVertices(): Vector[] {
        if (this.mVertexPool === null) {
            throw new Error('Polyhedron3: the vertex pool is null.');
        }
        return this.mVertexPool;
    }

    getUniqueIndices(): number[] {
        return this.mUniqueIndices;
    }

    getIndices(): number[] {
        return this.mIndices;
    }

    counterClockwise(): boolean {
        return this.mCounterClockwise;
    }

    // Geometric queries.
    computeVertexAverage(): Vector {
        let average = new Vector(3);
        if (this.mVertexPool !== null) {
            const vertexPool = this.mVertexPool;
            for (const index of this.mUniqueIndices) {
                average = add(average, vertexPool[index]);
            }
            average = div(average, this.mUniqueIndices.length);
        }
        return average;
    }

    computeSurfaceArea(): number {
        let surfaceArea = 0;
        if (this.mVertexPool !== null) {
            const vertexPool = this.mVertexPool;
            const numTriangles = Math.floor(this.mIndices.length / 3);
            for (let t = 0, j = 0; t < numTriangles; ++t) {
                const v0 = this.mIndices[j++];
                const v1 = this.mIndices[j++];
                const v2 = this.mIndices[j++];
                const edge0 = sub(vertexPool[v1], vertexPool[v0]);
                const edge1 = sub(vertexPool[v2], vertexPool[v0]);
                surfaceArea += length(cross(edge0, edge1));
            }
            surfaceArea *= 0.5;
        }
        return surfaceArea;
    }

    computeVolume(): number {
        let volume = 0;
        if (this.mVertexPool !== null) {
            const vertexPool = this.mVertexPool;
            const numTriangles = Math.floor(this.mIndices.length / 3);
            for (let t = 0, j = 0; t < numTriangles; ++t) {
                const v0 = this.mIndices[j++];
                const v1 = this.mIndices[j++];
                const v2 = this.mIndices[j++];
                volume += dotCross(vertexPool[v0], vertexPool[v1],
                    vertexPool[v2]);
            }
            volume /= 6;
        }
        return Math.abs(volume);
    }
}

