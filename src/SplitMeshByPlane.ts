// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) SplitMeshByPlane.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The algorithm for splitting a mesh by a plane is described in
// https://www.geometrictools.com/Documentation/ClipMesh.pdf
// Currently, the code here does not include generating a closed mesh (from
// the "positive" and "zero" vertices) by attaching triangulated faces to the
// mesh, where those faces live in the splitting plane. (TODO (upstream): Add
// this code.)
//
// Port notes:
// * The C++ class template is SplitMeshByPlane<Real>; Real becomes 'number'
//   per PORTING.md, so the port is a single non-generic class. 'operator()'
//   becomes compute() and the three output reference parameters become the
//   fields of SplitMeshByPlaneResult.
// * Plane3<Real> is the Hyperplane of dimension 3 (see Hyperplane.ts, which
//   exports the doc-only alias Plane3) and Vector3<Real> is a Vector of size
//   3.
// * 'clipVertices = vertices' copies the input in C++; the port clones each
//   Vector so the caller's vectors are never aliased or mutated.
// * The upstream 'std::map<EdgeKey<false>, std::pair<Vector3<Real>, int32_t>>'
//   becomes a Map keyed by the string form of the unordered EdgeKey. The
//   container is only ever queried by key (never iterated), so no sorted
//   iteration is needed.
// * Upstream reads the map with 'std::map::operator[]', which silently
//   default-constructs a (zero vector, index 0) entry for an absent key. The
//   port asserts instead; see the comment on getSplitIndex.
// * The three copy-pasted sign-change blocks of ClassifyEdges become a loop
//   over the triangle edges (v0,v1), (v1,v2), (v2,v0); the arithmetic for
//   each edge is unchanged.

import { EdgeKey } from './EdgeKey.js';
import { DistPointHyperplane } from './DistPointHyperplane.js';
import type { Hyperplane } from './Hyperplane.js';
import { logAssert } from './Logger.js';
import { Vector, add, mul, sub } from './Vector.js';

export interface SplitMeshByPlaneResult {
    // The input vertices followed by the vertices generated at the
    // edge-plane intersections.
    clipVertices: Vector[];

    // Index triples for the triangles on the negative side of the plane.
    negIndices: number[];

    // Index triples for the triangles on the positive side of the plane.
    posIndices: number[];
}

export class SplitMeshByPlane {
    // Stores the signed distances from the vertices to the plane.
    private mSignedDistances: number[];

    // Stores the edges whose vertices are on opposite sides of the plane.
    // The key is the unordered pair of indices into the vertex array. The
    // value is the point of intersection of the edge with the plane and an
    // index into clipVertices (the index is larger than or equal to the
    // number of vertices of the incoming 'vertices').
    private mEMap: Map<string, { intr: Vector, index: number }>;

    constructor() {
        this.mSignedDistances = [];
        this.mEMap = new Map<string, { intr: Vector, index: number }>();
    }

    // The 'indices' are lookups into the 'vertices' array. The indices
    // represent a triangle mesh. The number of indices must be a multiple of
    // 3, each triple representing a triangle. If t is a triangle index, then
    // the triangle is formed by vertices[indices[3 * t + i]] for 0 <= i <= 2.
    // The outputs 'negIndices' and 'posIndices' are formatted similarly.
    compute(vertices: readonly Vector[], indices: readonly number[],
        plane: Hyperplane): SplitMeshByPlaneResult {
        logAssert(indices.length % 3 === 0,
            'The number of indices must be a multiple of 3.');

        this.mSignedDistances = new Array<number>(vertices.length).fill(0);
        this.mEMap.clear();

        // Make a copy of the incoming vertices. If the mesh intersects the
        // plane, new vertices must be generated. These are appended to the
        // clipVertices array.
        const clipVertices = vertices.map(v => v.clone());
        const negIndices: number[] = [];
        const posIndices: number[] = [];

        this.classifyVertices(clipVertices, plane);
        this.classifyEdges(clipVertices, indices);
        this.classifyTriangles(indices, negIndices, posIndices);

        return { clipVertices, negIndices, posIndices };
    }

    private classifyVertices(clipVertices: readonly Vector[],
        plane: Hyperplane): void {
        const query = new DistPointHyperplane();
        for (let i = 0; i < clipVertices.length; ++i) {
            this.mSignedDistances[i] =
                query.compute(clipVertices[i] as Vector, plane).signedDistance;
        }
    }

    private classifyEdges(clipVertices: Vector[],
        indices: readonly number[]): void {
        const numTriangles = indices.length / 3;
        let nextIndex = clipVertices.length;
        for (let i = 0; i < numTriangles; ++i) {
            const threeI = 3 * i;
            const v: [number, number, number] = [
                indices[threeI + 0] as number,
                indices[threeI + 1] as number,
                indices[threeI + 2] as number
            ];
            const sDist: [number, number, number] = [
                this.mSignedDistances[v[0]] as number,
                this.mSignedDistances[v[1]] as number,
                this.mSignedDistances[v[2]] as number
            ];

            // The change-in-sign tests are structured this way to avoid
            // numerical round-off problems. For example, sDist0 > 0 and
            // sDist1 < 0, but both are very small and sDist0 * sDist1 = 0
            // because of round-off errors. The tests also guarantee
            // consistency between this function and classifyTriangles, the
            // latter function using sign tests only on the individual sDist
            // values.
            for (let j0 = 0; j0 < 3; ++j0) {
                const j1 = (j0 + 1) % 3;
                const sDist0 = sDist[j0] as number;
                const sDist1 = sDist[j1] as number;
                if ((sDist0 > 0 && sDist1 < 0) || (sDist0 < 0 && sDist1 > 0)) {
                    const i0 = v[j0] as number;
                    const i1 = v[j1] as number;
                    const key = SplitMeshByPlane.edgeKeyString(i0, i1);
                    if (!this.mEMap.has(key)) {
                        const t = sDist0 / (sDist0 - sDist1);
                        const diff = sub(clipVertices[i1] as Vector,
                            clipVertices[i0] as Vector);
                        const intr = add(clipVertices[i0] as Vector,
                            mul(t, diff));
                        clipVertices.push(intr);
                        this.mEMap.set(key, { intr, index: nextIndex });
                        ++nextIndex;
                    }
                }
            }
        }
    }

    private static edgeKeyString(v0: number, v1: number): string {
        const key = new EdgeKey(false, v0, v1);
        return key.V[0] + ',' + key.V[1];
    }

    // Upstream uses 'mEMap[EdgeKey<false>(v0, v1)].second', which for an
    // absent key silently inserts a default-constructed value and returns
    // index 0. The classification logic guarantees the key is present, so the
    // port asserts rather than fabricating index 0 (which would corrupt the
    // output mesh).
    private getSplitIndex(v0: number, v1: number): number {
        const value = this.mEMap.get(SplitMeshByPlane.edgeKeyString(v0, v1));
        logAssert(value !== undefined,
            'The edge does not intersect the plane.');
        return (value as { intr: Vector, index: number }).index;
    }

    private classifyTriangles(indices: readonly number[],
        negIndices: number[], posIndices: number[]): void {
        const numTriangles = indices.length / 3;
        for (let i = 0; i < numTriangles; ++i) {
            const threeI = 3 * i;
            const v0 = indices[threeI + 0] as number;
            const v1 = indices[threeI + 1] as number;
            const v2 = indices[threeI + 2] as number;
            const sDist0 = this.mSignedDistances[v0] as number;
            const sDist1 = this.mSignedDistances[v1] as number;
            const sDist2 = this.mSignedDistances[v2] as number;

            if (sDist0 > 0) {
                if (sDist1 > 0) {
                    if (sDist2 > 0) {
                        // +++
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                    else if (sDist2 < 0) {
                        // ++-
                        this.splitTrianglePPM(negIndices, posIndices, v0, v1, v2);
                    }
                    else {
                        // ++0
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                }
                else if (sDist1 < 0) {
                    if (sDist2 > 0) {
                        // +-+
                        this.splitTrianglePPM(negIndices, posIndices, v2, v0, v1);
                    }
                    else if (sDist2 < 0) {
                        // +--
                        this.splitTriangleMMP(negIndices, posIndices, v1, v2, v0);
                    }
                    else {
                        // +-0
                        this.splitTrianglePMZ(negIndices, posIndices, v0, v1, v2);
                    }
                }
                else {
                    if (sDist2 > 0) {
                        // +0+
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                    else if (sDist2 < 0) {
                        // +0-
                        this.splitTriangleMPZ(negIndices, posIndices, v2, v0, v1);
                    }
                    else {
                        // +00
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                }
            }
            else if (sDist0 < 0) {
                if (sDist1 > 0) {
                    if (sDist2 > 0) {
                        // -++
                        this.splitTrianglePPM(negIndices, posIndices, v1, v2, v0);
                    }
                    else if (sDist2 < 0) {
                        // -+-
                        this.splitTriangleMMP(negIndices, posIndices, v2, v0, v1);
                    }
                    else {
                        // -+0
                        this.splitTriangleMPZ(negIndices, posIndices, v0, v1, v2);
                    }
                }
                else if (sDist1 < 0) {
                    if (sDist2 > 0) {
                        // --+
                        this.splitTriangleMMP(negIndices, posIndices, v0, v1, v2);
                    }
                    else if (sDist2 < 0) {
                        // ---
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                    else {
                        // --0
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                }
                else {
                    if (sDist2 > 0) {
                        // -0+
                        this.splitTrianglePMZ(negIndices, posIndices, v2, v0, v1);
                    }
                    else if (sDist2 < 0) {
                        // -0-
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                    else {
                        // -00
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                }
            }
            else {
                if (sDist1 > 0) {
                    if (sDist2 > 0) {
                        // 0++
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                    else if (sDist2 < 0) {
                        // 0+-
                        this.splitTrianglePMZ(negIndices, posIndices, v1, v2, v0);
                    }
                    else {
                        // 0+0
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                }
                else if (sDist1 < 0) {
                    if (sDist2 > 0) {
                        // 0-+
                        this.splitTriangleMPZ(negIndices, posIndices, v1, v2, v0);
                    }
                    else if (sDist2 < 0) {
                        // 0--
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                    else {
                        // 0-0
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                }
                else {
                    if (sDist2 > 0) {
                        // 00+
                        SplitMeshByPlane.appendTriangle(posIndices, v0, v1, v2);
                    }
                    else if (sDist2 < 0) {
                        // 00-
                        SplitMeshByPlane.appendTriangle(negIndices, v0, v1, v2);
                    }
                    else {
                        // 000, reject triangles lying in the plane
                    }
                }
            }
        }
    }

    private static appendTriangle(indices: number[], v0: number, v1: number,
        v2: number): void {
        indices.push(v0);
        indices.push(v1);
        indices.push(v2);
    }

    private splitTrianglePPM(negIndices: number[], posIndices: number[],
        v0: number, v1: number, v2: number): void {
        const v12 = this.getSplitIndex(v1, v2);
        const v20 = this.getSplitIndex(v2, v0);
        posIndices.push(v0);
        posIndices.push(v1);
        posIndices.push(v12);
        posIndices.push(v0);
        posIndices.push(v12);
        posIndices.push(v20);
        negIndices.push(v2);
        negIndices.push(v20);
        negIndices.push(v12);
    }

    private splitTriangleMMP(negIndices: number[], posIndices: number[],
        v0: number, v1: number, v2: number): void {
        const v12 = this.getSplitIndex(v1, v2);
        const v20 = this.getSplitIndex(v2, v0);
        negIndices.push(v0);
        negIndices.push(v1);
        negIndices.push(v12);
        negIndices.push(v0);
        negIndices.push(v12);
        negIndices.push(v20);
        posIndices.push(v2);
        posIndices.push(v20);
        posIndices.push(v12);
    }

    private splitTrianglePMZ(negIndices: number[], posIndices: number[],
        v0: number, v1: number, v2: number): void {
        const v01 = this.getSplitIndex(v0, v1);
        posIndices.push(v2);
        posIndices.push(v0);
        posIndices.push(v01);
        negIndices.push(v2);
        negIndices.push(v01);
        negIndices.push(v1);
    }

    private splitTriangleMPZ(negIndices: number[], posIndices: number[],
        v0: number, v1: number, v2: number): void {
        const v01 = this.getSplitIndex(v0, v1);
        negIndices.push(v2);
        negIndices.push(v0);
        negIndices.push(v01);
        posIndices.push(v2);
        posIndices.push(v01);
        posIndices.push(v1);
    }
}
