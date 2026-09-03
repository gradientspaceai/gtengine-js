// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PlanarMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The planar mesh class is convenient for many applications involving
// searches for triangles containing a specified point.
//
// The input mesh should be consistently oriented, say, the triangles are
// counterclockwise ordered. The vertices should be consistent with this
// ordering. However, floating-point rounding errors in generating the
// vertices can cause apparent fold-over of the mesh; that is, theoretically
// the vertex geometry supports counterclockwise geometry but numerical errors
// cause an inconsistency. This can manifest in the toLine tests whereby
// cycles of triangles occur in the linear walk. When cycles occur,
// getContainingTriangle(P,startTriangle) will iterate numTriangles times
// before reporting that the triangle cannot be found, which is a very slow
// process. The function getContainingTriangleVisited(P,startTriangle) is
// provided to avoid the performance loss, trapping a cycle the first time and
// exiting, but again reporting that the triangle cannot be found. If you know
// that the query should be (theoretically) successful and it fails by
// returning -1, then perform an exhaustive search over the triangles using
// contains(triangle,P), which does not require the triangles to be ordered.
//
// Port notes:
//   - The upstream class template has three numeric parameters, InputType
//     (the type of the input vertices), ComputeType (the type used by the
//     PrimalQuery2 predicates) and RationalType (the type used by the
//     barycentric computation). The port instantiates all three as 'number',
//     which is what the ported PrimalQuery2 supports; mComputeVertices is
//     therefore a copy of the input vertices rather than a conversion of
//     them, and it is kept as a separate array so that the class owns the
//     data handed to PrimalQuery2 and PointInPolygon2.
//   - The two upstream constructors, one taking an index array and one taking
//     an ETManifoldMesh, become the static factories fromIndices() and
//     fromMesh() (the precedent for ambiguous C++ constructor overloads).
//   - The upstream (numVertices, pointer) and (numTriangles, pointer) pairs
//     become arrays. As upstream, the vertex array is referenced, not copied,
//     by getVertices(); the port also stores the separate compute-vertex
//     copy described above.
//   - The overloads GetVertices(t,...), GetIndices(t,...) and
//     GetAdjacencies(t,...), which write to an output array and return a
//     bool, become getTriangleVertices(t), getTriangleIndices(t) and
//     getTriangleAdjacencies(t), which return the triple or null. Likewise
//     GetBarycentrics returns the triple or null.
//   - The GetContainingTriangle overload with the std::set<int32_t>& output
//     parameter becomes getContainingTriangleVisited, which returns
//     { triangle, visited }.
//   - The fromMesh() factory iterates the triangles of the input mesh in
//     increasing triangle-key order (the order in which ETManifoldMesh.
//     getTriangles() reports them), which fixes the triangle ordering of the
//     resulting planar mesh. Upstream iterates an unordered_map, so its
//     ordering is implementation-defined.

import { logAssert } from './Logger.js';
import { PointInPolygon2 } from './ContPointInPolygon2.js';
import { ETManifoldMesh } from './ETManifoldMesh.js';
import { PrimalQuery2 } from './PrimalQuery2.js';
import { TriangleKey } from './TriangleKey.js';
import { Vector } from './Vector.js';
import { computeBarycentrics2 } from './Vector2.js';

export class PlanarMesh {
    private mNumVertices: number;
    private mVertices: readonly Vector[];
    private mNumTriangles: number;
    private mIndices: number[];
    private mMesh: ETManifoldMesh;
    private mTriIndexMap: Map<string, number>;
    private mAdjacencies: number[];
    private mComputeVertices: Vector[];
    private mQuery: PrimalQuery2;

    private constructor() {
        this.mNumVertices = 0;
        this.mVertices = [];
        this.mNumTriangles = 0;
        this.mIndices = [];
        this.mMesh = new ETManifoldMesh();
        this.mTriIndexMap = new Map<string, number>();
        this.mAdjacencies = [];
        this.mComputeVertices = [];
        this.mQuery = new PrimalQuery2();
    }

    // Construction. The inputs must represent a manifold mesh of triangles in
    // the plane. The index array must have 3*numTriangles elements, each
    // triple of indices representing a triangle in the mesh. Each index is
    // into the 'vertices' array.
    //
    // As upstream, if the mesh insertion of a triangle fails (the triangle is
    // a duplicate), the construction stops silently and the returned object
    // has no triangles and no vertices. The insertion throws for nonmanifold
    // input (see ETManifoldMesh.throwOnNonmanifoldInsertion).
    static fromIndices(vertices: readonly Vector[], indices: readonly number[]): PlanarMesh {
        const mesh = new PlanarMesh();

        logAssert(vertices.length >= 3 && indices.length >= 3 && indices.length % 3 === 0,
            'Invalid input.');
        const numTriangles = indices.length / 3;

        // Create a mesh in order to get adjacency information.
        for (let t = 0, current = 0; t < numTriangles; ++t) {
            const v0 = indices[current++];
            const v1 = indices[current++];
            const v2 = indices[current++];
            if (!mesh.mMesh.insert(v0, v1, v2)) {
                // The mesh object throws on nonmanifold inputs; a null return
                // means the triangle is already in the mesh.
                return mesh;
            }
        }

        // We have a valid mesh.
        mesh.createVertices(vertices);

        // Build the adjacency graph using the triangle ordering implied by
        // the indices, not the mesh triangle map, to preserve the triangle
        // ordering of the input indices.
        mesh.mNumTriangles = numTriangles;
        mesh.mIndices = indices.slice();

        for (let t = 0, vIndex = 0; t < numTriangles; ++t) {
            const v0 = indices[vIndex++];
            const v1 = indices[vIndex++];
            const v2 = indices[vIndex++];
            const key = new TriangleKey(true, v0, v1, v2).mapKey();
            // The upstream std::map::insert does not overwrite an existing
            // entry. The keys are distinct here, because a repeated triangle
            // would have caused the early return above.
            if (!mesh.mTriIndexMap.has(key)) {
                mesh.mTriIndexMap.set(key, t);
            }
        }

        mesh.mAdjacencies = new Array<number>(3 * numTriangles).fill(-1);
        for (let t = 0, base = 0; t < numTriangles; ++t, base += 3) {
            const tri = mesh.mMesh.getTriangle(indices[base], indices[base + 1],
                indices[base + 2]);
            logAssert(tri !== null, 'Unexpected condition.');
            for (let i = 0; i < 3; ++i) {
                const adj = tri.T[i];
                if (adj) {
                    const key = new TriangleKey(true, adj.V[0], adj.V[1], adj.V[2]).mapKey();
                    const index = mesh.mTriIndexMap.get(key);
                    logAssert(index !== undefined, 'Unexpected condition.');
                    mesh.mAdjacencies[base + i] = index;
                }
                else {
                    mesh.mAdjacencies[base + i] = -1;
                }
            }
        }

        return mesh;
    }

    // Construction from an existing manifold mesh. The triangle ordering of
    // the planar mesh is the increasing triangle-key ordering of the input
    // mesh.
    static fromMesh(vertices: readonly Vector[], inMesh: ETManifoldMesh): PlanarMesh {
        const mesh = new PlanarMesh();

        logAssert(vertices.length >= 3 && inMesh.getNumTriangles() >= 1,
            'Invalid input in PlanarMesh constructor.');

        // We have a valid mesh.
        mesh.createVertices(vertices);

        // Build the adjacency graph using the triangle ordering implied by
        // the mesh triangle map.
        const triangles = inMesh.getTriangles();
        mesh.mNumTriangles = triangles.length;
        mesh.mIndices = new Array<number>(3 * mesh.mNumTriangles).fill(0);

        let tIndex = 0;
        let vIndex = 0;
        for (const tri of triangles) {
            const key = new TriangleKey(true, tri.V[0], tri.V[1], tri.V[2]).mapKey();
            if (!mesh.mTriIndexMap.has(key)) {
                mesh.mTriIndexMap.set(key, tIndex);
            }
            ++tIndex;
            for (let i = 0; i < 3; ++i, ++vIndex) {
                mesh.mIndices[vIndex] = tri.V[i];
            }
        }

        mesh.mAdjacencies = new Array<number>(3 * mesh.mNumTriangles).fill(-1);
        vIndex = 0;
        for (const tri of triangles) {
            for (let i = 0; i < 3; ++i, ++vIndex) {
                const adj = tri.T[i];
                if (adj) {
                    const key = new TriangleKey(true, adj.V[0], adj.V[1], adj.V[2]).mapKey();
                    const index = mesh.mTriIndexMap.get(key);
                    logAssert(index !== undefined, 'Unexpected condition.');
                    mesh.mAdjacencies[vIndex] = index;
                }
                else {
                    mesh.mAdjacencies[vIndex] = -1;
                }
            }
        }

        return mesh;
    }

    // Mesh information.
    getNumVertices(): number {
        return this.mNumVertices;
    }

    getNumTriangles(): number {
        return this.mNumTriangles;
    }

    // The input vertex array, referenced not copied (as upstream).
    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getIndices(): readonly number[] {
        return this.mIndices;
    }

    getAdjacencies(): readonly number[] {
        return this.mAdjacencies;
    }

    // Containment queries. The function getContainingTriangle works correctly
    // when the planar mesh is a convex set. If the mesh is not convex, it is
    // possible that the linear-walk search algorithm exits the mesh before
    // finding a containing triangle. For example, a C-shaped mesh can contain
    // a point in the top branch of the "C". A starting point in the bottom
    // branch of the "C" will lead to the search exiting the bottom branch and
    // having no path to walk to the top branch. If your mesh is not convex
    // and you want a correct containment query, you will have to append
    // "outside" triangles to your mesh to form a convex set.
    getContainingTriangle(p: Vector, startTriangle: number = 0): number {
        const test = this.toComputeType(p);

        // Use triangle edges as binary separating lines.
        let triangle = startTriangle;
        for (let i = 0; i < this.mNumTriangles; ++i) {
            const ibase = 3 * triangle;
            const v0 = this.mIndices[ibase];
            const v1 = this.mIndices[ibase + 1];
            const v2 = this.mIndices[ibase + 2];

            if (this.mQuery.toLine(test, v0, v1) > 0) {
                triangle = this.mAdjacencies[ibase];
                if (triangle === -1) {
                    return -1;
                }
                continue;
            }

            if (this.mQuery.toLine(test, v1, v2) > 0) {
                triangle = this.mAdjacencies[ibase + 1];
                if (triangle === -1) {
                    return -1;
                }
                continue;
            }

            if (this.mQuery.toLine(test, v2, v0) > 0) {
                triangle = this.mAdjacencies[ibase + 2];
                if (triangle === -1) {
                    return -1;
                }
                continue;
            }

            return triangle;
        }

        return -1;
    }

    // The linear walk that traps cycles. The returned 'visited' set contains
    // the triangles that were visited during the walk.
    getContainingTriangleVisited(p: Vector, startTriangle: number = 0):
        { triangle: number, visited: Set<number> } {
        const test = this.toComputeType(p);
        const visited = new Set<number>();

        // Use triangle edges as binary separating lines.
        let triangle = startTriangle;
        for (let i = 0; i < this.mNumTriangles; ++i) {
            visited.add(triangle);
            const ibase = 3 * triangle;
            const v0 = this.mIndices[ibase];
            const v1 = this.mIndices[ibase + 1];
            const v2 = this.mIndices[ibase + 2];

            if (this.mQuery.toLine(test, v0, v1) > 0) {
                triangle = this.mAdjacencies[ibase];
                if (triangle === -1 || visited.has(triangle)) {
                    return { triangle: -1, visited };
                }
                continue;
            }

            if (this.mQuery.toLine(test, v1, v2) > 0) {
                triangle = this.mAdjacencies[ibase + 1];
                if (triangle === -1 || visited.has(triangle)) {
                    return { triangle: -1, visited };
                }
                continue;
            }

            if (this.mQuery.toLine(test, v2, v0) > 0) {
                triangle = this.mAdjacencies[ibase + 2];
                if (triangle === -1 || visited.has(triangle)) {
                    return { triangle: -1, visited };
                }
                continue;
            }

            return { triangle, visited };
        }

        return { triangle: -1, visited };
    }

    // The vertices of triangle t, or null when t is out of range.
    getTriangleVertices(t: number): [Vector, Vector, Vector] | null {
        const indices = this.getTriangleIndices(t);
        if (indices === null) {
            return null;
        }
        return [this.mVertices[indices[0]], this.mVertices[indices[1]],
            this.mVertices[indices[2]]];
    }

    // The vertex indices of triangle t, or null when t is out of range.
    getTriangleIndices(t: number): [number, number, number] | null {
        if (0 <= t && t < this.mNumTriangles) {
            const vIndex = 3 * t;
            return [this.mIndices[vIndex], this.mIndices[vIndex + 1],
                this.mIndices[vIndex + 2]];
        }
        return null;
    }

    // The adjacent triangles of triangle t, or null when t is out of range.
    // An adjacency of -1 indicates that the corresponding edge is a boundary
    // edge of the mesh.
    getTriangleAdjacencies(t: number): [number, number, number] | null {
        if (0 <= t && t < this.mNumTriangles) {
            const vIndex = 3 * t;
            return [this.mAdjacencies[vIndex], this.mAdjacencies[vIndex + 1],
                this.mAdjacencies[vIndex + 2]];
        }
        return null;
    }

    // The barycentric coordinates of P with respect to triangle t. The result
    // is null when t is out of range or when the triangle is degenerate.
    getBarycentrics(t: number, p: Vector): [number, number, number] | null {
        const indices = this.getTriangleIndices(t);
        if (indices !== null) {
            const rtP = this.toComputeType(p);
            const rtV0 = this.mComputeVertices[indices[0]];
            const rtV1 = this.mComputeVertices[indices[1]];
            const rtV2 = this.mComputeVertices[indices[2]];
            const result = computeBarycentrics2(rtP, rtV0, rtV1, rtV2);
            if (result.valid) {
                return result.bary;
            }
        }
        return null;
    }

    // Test whether P is contained by the specified triangle. The triangles do
    // not have to be ordered for this query.
    contains(triangle: number, p: Vector): boolean {
        const test = this.toComputeType(p);
        // Upstream does no range checking here; the port reports the invalid
        // index rather than reading out of bounds.
        logAssert(0 <= triangle && triangle < this.mNumTriangles,
            'Invalid triangle index.');
        const base = 3 * triangle;
        const v: Vector[] = [
            this.mComputeVertices[this.mIndices[base]],
            this.mComputeVertices[this.mIndices[base + 1]],
            this.mComputeVertices[this.mIndices[base + 2]]
        ];
        const pip = new PointInPolygon2(v);
        return pip.contains(test);
    }

    private createVertices(vertices: readonly Vector[]): void {
        this.mNumVertices = vertices.length;
        this.mVertices = vertices;
        this.mComputeVertices = new Array<Vector>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mComputeVertices[i] = this.toComputeType(vertices[i]);
        }
        this.mQuery.set(this.mNumVertices, this.mComputeVertices);
    }

    // The port of the (ComputeType)P conversion; with ComputeType = number
    // this is a copy of the 2D point.
    private toComputeType(p: Vector): Vector {
        logAssert(p.size === 2, 'PlanarMesh: the points must be 2D.');
        return Vector.fromArray([p.values[0], p.values[1]]);
    }
}
