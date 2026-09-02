// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Delaunay3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Delaunay tetrahedralization of points (intrinsic dimensionality 3).
//   VQ = number of vertices
//   V  = array of vertices
//   TQ = number of tetrahedra
//   I  = Array of 4-tuples of indices into V that represent the tetrahedra
//        (4*TQ total elements). Access via getIndices().
//   A  = Array of 4-tuples of indices into I that represent the adjacent
//        tetrahedra (4*TQ total elements). Access via getAdjacencies().
// The i-th tetrahedron has vertices
//   vertex[0] = V[I[4*i+0]]
//   vertex[1] = V[I[4*i+1]]
//   vertex[2] = V[I[4*i+2]]
//   vertex[3] = V[I[4*i+3]]
// and face index triples listed below. The face vertex ordering when viewed
// from outside the tetrahedron is counterclockwise.
//   face[0] = <I[4*i+1],I[4*i+2],I[4*i+3]>
//   face[1] = <I[4*i+0],I[4*i+3],I[4*i+2]>
//   face[2] = <I[4*i+0],I[4*i+1],I[4*i+3]>
//   face[3] = <I[4*i+0],I[4*i+2],I[4*i+1]>
// The tetrahedra adjacent to these faces have indices
//   adjacent[j] = A[4*i+j] is the tetrahedron opposite vertex[j], so it is
//                 the tetrahedron sharing face[j].
// If there is no adjacent tetrahedron, the A[*] value is set to -1.
//
// The class uses a blend of interval arithmetic (the fast path) and exact
// rational arithmetic (the fallback when the interval straddles zero) for the
// toPlane and toCircumsphere sign classifications, which is what makes the
// tetrahedralization robust.
//
// Port notes:
// * Upstream Delaunay3.h declares a variadic class template with two
//   specializations: the deprecated Delaunay3<InputType, ComputeType>, which
//   upstream states "will be removed in a future release", and the
//   replacement Delaunay3<T>. Only the replacement is ported. The deprecated
//   class routes its predicates through PrimalQuery3<ComputeType>; the
//   replacement embeds its own interval/rational predicates, which is the
//   behavior ported here (src/PrimalQuery3.ts is number-only and is not used
//   by this file).
// * Upstream selects InputRational = BSNumber<UIntegerFP32<2 or 4>> and
//   ComputeRational = BSNumber<UIntegerFP32<44 or 330>> with a preallocated
//   pool of compute-rational scratch values (mCRPool) and a Copy() helper
//   that widens an input rational into a compute rational. The port's
//   BSNumber is bigint-backed and grows as needed, so the fixed word counts,
//   the two distinct rational types, the pool and Copy() are all unnecessary
//   and are dropped. The exact-arithmetic results are identical.
// * std::unordered_set of directed TriangleKey<true> becomes
//   DirectedTriangleKeySet, a Map keyed by the normalized directed triple and
//   iterated in sorted key order so the tetrahedralization is deterministic.
//   std::unordered_set<Tetrahedron*> becomes a Set, whose JavaScript
//   insertion order is deterministic.
// * size_t 'negOne' (std::numeric_limits<size_t>::max()) becomes -1.
// * The 'bool GetHull(std::vector<size_t>&)' output-parameter form becomes
//   'getHull(): number[]'; the bool was 'true' on every path that returns.
// * The upstream operator() becomes compute(vertices). The overloaded
//   GetIndices(t, array)/GetAdjacencies(t, array) accessors become
//   getTetrahedronIndices(t)/getTetrahedronAdjacencies(t), returning null
//   instead of 'false'.

import { logAssert, logError } from './Logger';
import { Line } from './Line';
import type { Line3 } from './Line';
import { Hyperplane } from './Hyperplane';
import type { Plane3 } from './Hyperplane';
import { Vector } from './Vector';
import { IntrinsicsVector3, unitCross } from './Vector3';
import { BSNumber } from './BSNumber';
import { SWInterval } from './SWInterval';
import { TriangleKey } from './TriangleKey';
import { TetrahedronKey } from './TetrahedronKey';
import { TSManifoldMesh, TSManifoldMeshTetrahedron } from './TSManifoldMesh';

// A rational 3D point, the port of Vector3<InputRational>.
type RationalPoint3 = [BSNumber, BSNumber, BSNumber];

// The sentinel for "the query point is not one of the input vertices"; the
// port of the upstream size_t constant negOne.
const negOne = -1;

// The port of std::unordered_set<TriangleKey<true>, ...>. The keys are
// normalized by TriangleKey<true> (rotated so the smallest index is first,
// which preserves orientation). The port iterates in increasing key order so
// that the constructed tetrahedralization does not depend on hash-table
// ordering.
class DirectedTriangleKeySet {
    private mMap: Map<string, [number, number, number]>;

    constructor() {
        this.mMap = new Map<string, [number, number, number]>();
    }

    insert(v0: number, v1: number, v2: number): void {
        const key = new TriangleKey(true, v0, v1, v2);
        const triple: [number, number, number] = [key.V[0], key.V[1], key.V[2]];
        this.mMap.set(`${triple[0]},${triple[1]},${triple[2]}`, triple);
    }

    get size(): number {
        return this.mMap.size;
    }

    keys(): [number, number, number][] {
        const keys = Array.from(this.mMap.values());
        keys.sort((a, b) => {
            if (a[0] !== b[0]) {
                return a[0] - b[0];
            }
            if (a[1] !== b[1]) {
                return a[1] - b[1];
            }
            return a[2] - b[2];
        });
        return keys;
    }
}

// The port of the SearchInfo struct used by getContainingTetrahedron(). The
// first tetrahedron searched is 'initialTetrahedron'. On return, 'path'
// stores the (ordered) tetrahedron indices visited during the search, with
// 'numPath' valid entries. The last visited tetrahedron has index
// 'finalTetrahedron' and vertex indices finalV[0,1,2,3] stored in volumetric
// counterclockwise order. The last face of the search is
// <finalV[0],finalV[1],finalV[2]>. For spatially coherent query points in
// numerous calls, specify 'finalTetrahedron' of the previous call as
// 'initialTetrahedron' of the next call to reduce search times.
export class Delaunay3SearchInfo {
    initialTetrahedron: number;
    numPath: number;
    finalTetrahedron: number;
    finalV: [number, number, number, number];
    path: number[];

    constructor() {
        this.initialTetrahedron = 0;
        this.numPath = 0;
        this.finalTetrahedron = 0;
        this.finalV = [0, 0, 0, 0];
        this.path = [];
    }
}

export class Delaunay3 {
    // The vertices used for geometric queries. The input vertices are
    // read-only, so they can be represented exactly by rational numbers.
    protected mNumVertices: number;
    protected mVertices: readonly Vector[];
    protected mIRVertices: RationalPoint3[];

    protected mGraph: TSManifoldMesh;

    // If a vertex occurs multiple times in the 'vertices' input, the first
    // processed occurrence of that vertex has an index stored in this array.
    // If there are no duplicates, then mDuplicates[i] = i for all i.
    protected mDuplicates: number[];
    protected mNumUniqueVertices: number;

    // If the intrinsic dimension of the input vertices is 0, 1 or 2,
    // compute() returns early. The caller is responsible for retrieving the
    // dimension and taking an alternate path should the dimension be smaller
    // than 3. If the dimension is 0, all vertices are the same. If the
    // dimension is 1, the vertices lie on a line, in which case the caller
    // can project vertices[] onto the line for further processing. If the
    // dimension is 2, the vertices lie on a plane, in which case the caller
    // can project vertices[] onto the plane and apply Delaunay2.
    protected mDimension: number;
    protected mLine: Line3;
    protected mPlane: Plane3;

    // These are computed by updateIndicesAdjacencies(). They are used for
    // point-containment queries in the tetrahedron mesh.
    protected mNumTetrahedra: number;
    protected mIndices: number[];
    protected mAdjacencies: number[];

    // The query point for update(), getContainingTetrahedron() and
    // getAndRemoveInsertionPolyhedron() when the point is not an input
    // vertex. toPlane() and toCircumsphere() are passed indices into the
    // vertex array. When the index is valid, mVertices[] and mIRVertices[]
    // are used for lookups. When the index is negOne, the query point is
    // used.
    private mQueryPoint: Vector;
    private mIRQueryPoint: RationalPoint3;

    constructor() {
        this.mNumVertices = 0;
        this.mVertices = [];
        this.mIRVertices = [];
        this.mGraph = new TSManifoldMesh();
        this.mDuplicates = [];
        this.mNumUniqueVertices = 0;
        this.mDimension = 0;
        this.mLine = new Line(3);
        this.mLine.direction.makeZero();
        this.mPlane = Hyperplane.fromNormalConstant(new Vector(3), 0);
        this.mNumTetrahedra = 0;
        this.mIndices = [];
        this.mAdjacencies = [];
        this.mQueryPoint = new Vector(3);
        this.mIRQueryPoint = [BSNumber.fromNumber(0), BSNumber.fromNumber(0),
            BSNumber.fromNumber(0)];
    }

    // The input is the array of vertices whose Delaunay tetrahedralization is
    // required. The return value is 'true' if and only if the intrinsic
    // dimension of the points is 3. If the intrinsic dimension is 2, the
    // points lie exactly on a plane which is accessible via getPlane(). If
    // the intrinsic dimension is 1, the points lie exactly on a line which is
    // accessible via getLine(). If the intrinsic dimension is 0, the points
    // are all the same point.
    compute(vertices: readonly Vector[]): boolean {
        // Initialize values in case they were set by a previous call to
        // compute().
        logAssert(vertices.length > 0, 'Invalid argument.');
        for (const vertex of vertices) {
            logAssert(vertex.size === 3, 'Delaunay3 requires 3D vertices.');
        }

        this.mNumVertices = vertices.length;
        this.mVertices = vertices;
        this.mIRVertices = [];
        this.mGraph = new TSManifoldMesh();
        this.mDuplicates = [];
        this.mNumUniqueVertices = 0;
        this.mDimension = 0;
        this.mLine = new Line(3);
        this.mLine.direction.makeZero();
        this.mPlane = Hyperplane.fromNormalConstant(new Vector(3), 0);
        this.mNumTetrahedra = 0;
        this.mIndices = [];
        this.mAdjacencies = [];
        this.mQueryPoint = new Vector(3);
        this.mIRQueryPoint = [BSNumber.fromNumber(0), BSNumber.fromNumber(0),
            BSNumber.fromNumber(0)];

        // Compute the intrinsic dimension and return early if that dimension
        // is 0, 1 or 2.
        const info = new IntrinsicsVector3(vertices, 0);
        if (info.dimension === 0) {
            // The vertices are the same point.
            this.mDimension = 0;
            this.mLine.origin = info.origin.clone();
            return false;
        }

        if (info.dimension === 1) {
            // The vertices are collinear.
            this.mDimension = 1;
            this.mLine.origin = info.origin.clone();
            this.mLine.direction = info.direction[0].clone();
            return false;
        }

        if (info.dimension === 2) {
            // The vertices are coplanar.
            this.mDimension = 2;
            this.mPlane = Hyperplane.fromNormalOrigin(
                unitCross(info.direction[0], info.direction[1]), info.origin);
            return false;
        }

        // The vertices necessarily will have a tetrahedralization.
        this.mDimension = 3;

        // Convert the floating-point inputs to rational type.
        this.mIRVertices = new Array<RationalPoint3>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            const v = vertices[i].values;
            this.mIRVertices[i] = [BSNumber.fromNumber(v[0]),
                BSNumber.fromNumber(v[1]), BSNumber.fromNumber(v[2])];
        }

        // Assume initially the vertices are unique. If duplicates are found
        // during the Delaunay update, mDuplicates[] will be modified
        // accordingly.
        this.mDuplicates = new Array<number>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mDuplicates[i] = i;
        }

        // Insert the nondegenerate tetrahedron constructed by the intrinsics
        // computation. This is necessary for the circumsphere-visibility
        // algorithm to work correctly.
        const extreme: [number, number, number, number] = [info.extreme[0],
            info.extreme[1], info.extreme[2], info.extreme[3]];
        if (!info.extremeCCW) {
            const save = extreme[2];
            extreme[2] = extreme[3];
            extreme[3] = save;
        }

        const inserted = this.mGraph.insert(extreme[0], extreme[1], extreme[2],
            extreme[3]);
        logAssert(inserted !== null, 'The tetrahedron should not be degenerate.');

        // Incrementally update the tetrahedralization. The set of processed
        // points is maintained to eliminate duplicates. The upstream
        // std::unordered_set<ProcessedVertex> is a Map keyed by the exact
        // floating-point coordinate triple; the mapped value is the
        // 'location' member.
        const processed = new Map<string, number>();
        const keyOf = (i: number): string => {
            const v = vertices[i].values;
            return `${v[0]},${v[1]},${v[2]}`;
        };
        for (let i = 0; i < 4; ++i) {
            const j = extreme[i];
            processed.set(keyOf(j), j);
            this.mDuplicates[j] = j;
        }
        for (let i = 0; i < this.mNumVertices; ++i) {
            const key = keyOf(i);
            const location = processed.get(key);
            if (location === undefined) {
                this.update(i);
                processed.set(key, i);
                this.mDuplicates[i] = i;
            }
            else {
                this.mDuplicates[i] = location;
            }
        }
        this.mNumUniqueVertices = processed.size;

        // Assign integer values to the tetrahedra for use by the caller and
        // copy the tetrahedra information to compact arrays mIndices and
        // mAdjacencies.
        this.updateIndicesAdjacencies();

        return true;
    }

    // Dimensional information. If getDimension() returns 1, the points lie on
    // a line P+t*D. You can sort these if you need a polyline output by
    // projecting onto the line each vertex X = P+t*D, where t = Dot(D,X-P).
    // If getDimension() returns 2, the points lie on a plane P+s*U+t*V. You
    // can project each vertex X = P+s*U+t*V, where s = Dot(U,X-P) and
    // t = Dot(V,X-P) and then apply Delaunay2 to the (s,t) tuples.
    getDimension(): number {
        return this.mDimension;
    }

    getLine(): Line3 {
        return this.mLine;
    }

    getPlane(): Plane3 {
        return this.mPlane;
    }

    // Member access.
    //
    // Upstream bug (fixed in the port): GetNumVertices() returns
    // mIRVertices.size(), which is 0 whenever the intrinsic dimension is 0, 1
    // or 2 because operator() clears mIRVertices and returns before resizing
    // it. The port returns mNumVertices, which is the number of input
    // vertices on every path.
    getNumVertices(): number {
        return this.mNumVertices;
    }

    getVertices(): readonly Vector[] {
        return this.mVertices;
    }

    getNumUniqueVertices(): number {
        return this.mNumUniqueVertices;
    }

    // If 'vertices' has no duplicates, getDuplicates()[i] = i for all i. If
    // vertices[i] is the first occurrence of a vertex and if vertices[j] is
    // found later, then getDuplicates()[j] = i.
    getDuplicates(): readonly number[] {
        return this.mDuplicates;
    }

    getNumTetrahedra(): number {
        return this.mNumTetrahedra;
    }

    getGraph(): TSManifoldMesh {
        return this.mGraph;
    }

    getIndices(): readonly number[] {
        return this.mIndices;
    }

    getAdjacencies(): readonly number[] {
        return this.mAdjacencies;
    }

    // Locate those tetrahedra faces that do not share other tetrahedra. The
    // returned array has hull.length = 3*numFaces indices, each triple
    // representing a triangle. The triangles are counterclockwise ordered
    // when viewed from outside the hull. The dimension must be 3.
    getHull(): number[] {
        if (this.mDimension === 3) {
            // Count the number of triangles that are not shared by two
            // tetrahedra.
            let numTriangles = 0;
            for (const adj of this.mAdjacencies) {
                if (adj === -1) {
                    ++numTriangles;
                }
            }

            if (numTriangles > 0) {
                // Enumerate the triangles. The prototypical case is the
                // single tetrahedron V[0] = (0,0,0), V[1] = (1,0,0),
                // V[2] = (0,1,0) and V[3] = (0,0,1) with no adjacent
                // tetrahedra. The mIndices[] array is <0,1,2,3>.
                //   i = 0, face = 0:
                //    skip index 0, <x,1,2,3>, no swap, triangle = <1,2,3>
                //   i = 1, face = 1:
                //    skip index 1, <0,x,2,3>, swap,    triangle = <0,3,2>
                //   i = 2, face = 2:
                //    skip index 2, <0,1,x,3>, no swap, triangle = <0,1,3>
                //   i = 3, face = 3:
                //    skip index 3, <0,1,2,x>, swap,    triangle = <0,2,1>
                // To guarantee counterclockwise order of triangles when
                // viewed outside the tetrahedron, the swap of the last two
                // indices occurs when face is an odd number.
                const hull = new Array<number>(3 * numTriangles);
                let current = 0, i = 0;
                for (const adj of this.mAdjacencies) {
                    if (adj === -1) {
                        const tetra = Math.floor(i / 4), face = i % 4;
                        for (let j = 0; j < 4; ++j) {
                            if (j !== face) {
                                hull[current++] = this.mIndices[4 * tetra + j];
                            }
                        }
                        if ((face % 2) !== 0) {
                            const save = hull[current - 1];
                            hull[current - 1] = hull[current - 2];
                            hull[current - 2] = save;
                        }
                    }
                    ++i;
                }
                return hull;
            }
            else {
                logError('Unexpected condition. There must be at least one tetrahedron.');
            }
        }
        else {
            logError('The dimension must be 3.');
        }
        // Unreachable; logError throws.
        return [];
    }

    // Copy Delaunay tetrahedra to compact arrays mIndices and mAdjacencies.
    updateIndicesAdjacencies(): void {
        // Assign integer values to the tetrahedra for use by the caller.
        const smap = this.mGraph.getTetrahedra();
        const permute = new Map<TSManifoldMeshTetrahedron, number>();
        for (let i = 0; i < smap.length; ++i) {
            permute.set(smap[i], i);
        }

        this.mNumTetrahedra = smap.length;
        const numIndices = 4 * this.mNumTetrahedra;
        if (numIndices > 0) {
            this.mIndices = new Array<number>(numIndices);
            this.mAdjacencies = new Array<number>(numIndices);
            let i = 0;
            for (const tetra of smap) {
                for (let j = 0; j < 4; ++j, ++i) {
                    this.mIndices[i] = tetra.V[j];
                    const adj = tetra.S[j];
                    this.mAdjacencies[i] = (adj !== null ? permute.get(adj) as number : -1);
                }
            }
        }
    }

    // Get the vertex indices for tetrahedron t. The function returns the
    // indices when the dimension is 3 and t is a valid tetrahedron index;
    // otherwise, it returns null.
    getTetrahedronIndices(t: number): [number, number, number, number] | null {
        if (this.mDimension === 3) {
            const numTetrahedra = this.mIndices.length / 4;
            if (0 <= t && t < numTetrahedra) {
                return [this.mIndices[4 * t], this.mIndices[4 * t + 1],
                    this.mIndices[4 * t + 2], this.mIndices[4 * t + 3]];
            }
        }
        return null;
    }

    // Get the indices for tetrahedra adjacent to tetrahedron t. The function
    // returns the adjacencies when the dimension is 3 and t is a valid
    // tetrahedron index; otherwise, it returns null.
    getTetrahedronAdjacencies(t: number): [number, number, number, number] | null {
        if (this.mDimension === 3) {
            const numTetrahedra = this.mIndices.length / 4;
            if (0 <= t && t < numTetrahedra) {
                return [this.mAdjacencies[4 * t], this.mAdjacencies[4 * t + 1],
                    this.mAdjacencies[4 * t + 2], this.mAdjacencies[4 * t + 3]];
            }
        }
        return null;
    }

    // Support for searching the tetrahedralization for a tetrahedron that
    // contains a point. If there is a containing tetrahedron, the returned
    // value is a tetrahedron index t with 0 <= t < getNumTetrahedra(). If
    // there is not a containing tetrahedron, -1 is returned. The computations
    // are performed using exact rational arithmetic.
    getContainingTetrahedron(inP: Vector, info: Delaunay3SearchInfo): number {
        logAssert(this.mDimension === 3, 'Invalid dimension for tetrahedron search.');
        logAssert(inP.size === 3, 'Delaunay3 requires 3D vertices.');

        this.mQueryPoint = inP.clone();
        this.mIRQueryPoint = [BSNumber.fromNumber(inP.values[0]),
            BSNumber.fromNumber(inP.values[1]), BSNumber.fromNumber(inP.values[2])];

        const numTetrahedra = this.mIndices.length / 4;
        info.path = new Array<number>(numTetrahedra).fill(0);
        info.numPath = 0;
        let tetrahedron: number;
        if (0 <= info.initialTetrahedron && info.initialTetrahedron < numTetrahedra) {
            tetrahedron = info.initialTetrahedron;
        }
        else {
            info.initialTetrahedron = 0;
            tetrahedron = 0;
        }

        // Use tetrahedron faces as binary separating planes.
        for (let i = 0; i < numTetrahedra; ++i) {
            const ibase = 4 * tetrahedron;
            const v0 = this.mIndices[ibase];
            const v1 = this.mIndices[ibase + 1];
            const v2 = this.mIndices[ibase + 2];
            const v3 = this.mIndices[ibase + 3];

            info.path[info.numPath++] = tetrahedron;
            info.finalTetrahedron = tetrahedron;
            info.finalV[0] = v0;
            info.finalV[1] = v1;
            info.finalV[2] = v2;
            info.finalV[3] = v3;

            // <V1,V2,V3> counterclockwise when viewed outside tetrahedron.
            if (this.toPlane(negOne, v1, v2, v3) > 0) {
                const adjacent = this.mAdjacencies[ibase];
                if (adjacent === -1) {
                    info.finalV[0] = v1;
                    info.finalV[1] = v2;
                    info.finalV[2] = v3;
                    info.finalV[3] = v0;
                    return negOne;
                }
                tetrahedron = adjacent;
                continue;
            }

            // <V0,V3,V2> counterclockwise when viewed outside tetrahedron.
            if (this.toPlane(negOne, v0, v2, v3) < 0) {
                const adjacent = this.mAdjacencies[ibase + 1];
                if (adjacent === -1) {
                    info.finalV[0] = v0;
                    info.finalV[1] = v2;
                    info.finalV[2] = v3;
                    info.finalV[3] = v1;
                    return negOne;
                }
                tetrahedron = adjacent;
                continue;
            }

            // <V0,V1,V3> counterclockwise when viewed outside tetrahedron.
            if (this.toPlane(negOne, v0, v1, v3) > 0) {
                const adjacent = this.mAdjacencies[ibase + 2];
                if (adjacent === -1) {
                    info.finalV[0] = v0;
                    info.finalV[1] = v1;
                    info.finalV[2] = v3;
                    info.finalV[3] = v2;
                    return negOne;
                }
                tetrahedron = adjacent;
                continue;
            }

            // <V0,V2,V1> counterclockwise when viewed outside tetrahedron.
            if (this.toPlane(negOne, v0, v1, v2) < 0) {
                const adjacent = this.mAdjacencies[ibase + 3];
                if (adjacent === -1) {
                    info.finalV[0] = v0;
                    info.finalV[1] = v1;
                    info.finalV[2] = v2;
                    info.finalV[3] = v3;
                    return negOne;
                }
                tetrahedron = adjacent;
                continue;
            }

            return tetrahedron;
        }

        logError('Unexpected termination of loop while searching for a tetrahedron.');
        return negOne;
    }

    // The floating-point query point for the given index, where negOne
    // selects the stored query point.
    private inPoint(index: number): readonly number[] {
        return (index !== negOne ? this.mVertices[index].values : this.mQueryPoint.values);
    }

    // The rational query point for the given index, where negOne selects the
    // stored query point.
    private irPoint(index: number): RationalPoint3 {
        return (index !== negOne ? this.mIRVertices[index] : this.mIRQueryPoint);
    }

    // Given a plane with origin V0 and normal N = Cross(V1-V0,V2-V0) and
    // given a query point P, toPlane returns
    //   +1, P on positive side of plane (side to which N points)
    //   -1, P on negative side of plane (side to which -N points)
    //    0, P on the plane
    protected toPlane(pIndex: number, v0Index: number, v1Index: number,
        v2Index: number): number {
        // The expression tree has 34 nodes consisting of 12 input leaves and
        // 22 compute nodes.

        // Use interval arithmetic to determine the sign if possible.
        const inP = this.inPoint(pIndex);
        const inV0 = this.mVertices[v0Index].values;
        const inV1 = this.mVertices[v1Index].values;
        const inV2 = this.mVertices[v2Index].values;

        const x0 = SWInterval.sub(inP[0], inV0[0]);
        const y0 = SWInterval.sub(inP[1], inV0[1]);
        const z0 = SWInterval.sub(inP[2], inV0[2]);
        const x1 = SWInterval.sub(inV1[0], inV0[0]);
        const y1 = SWInterval.sub(inV1[1], inV0[1]);
        const z1 = SWInterval.sub(inV1[2], inV0[2]);
        const x2 = SWInterval.sub(inV2[0], inV0[0]);
        const y2 = SWInterval.sub(inV2[1], inV0[1]);
        const z2 = SWInterval.sub(inV2[2], inV0[2]);
        const y0z1 = y0.mul(z1);
        const y0z2 = y0.mul(z2);
        const y1z0 = y1.mul(z0);
        const y1z2 = y1.mul(z2);
        const y2z0 = y2.mul(z0);
        const y2z1 = y2.mul(z1);
        const c0 = y1z2.sub(y2z1);
        const c1 = y2z0.sub(y0z2);
        const c2 = y0z1.sub(y1z0);
        const x0c0 = x0.mul(c0);
        const x1c1 = x1.mul(c1);
        const x2c2 = x2.mul(c2);
        const det = x0c0.add(x1c1).add(x2c2);

        const zero = 0;
        if (det.get(0) > zero) {
            return +1;
        }
        else if (det.get(1) < zero) {
            return -1;
        }

        // The exact sign of the determinant is not known, so compute the
        // determinant using rational arithmetic.
        const irP = this.irPoint(pIndex);
        const irV0 = this.mIRVertices[v0Index];
        const irV1 = this.mIRVertices[v1Index];
        const irV2 = this.mIRVertices[v2Index];

        const crP0 = irP[0], crP1 = irP[1], crP2 = irP[2];
        const crV00 = irV0[0], crV01 = irV0[1], crV02 = irV0[2];
        const crV10 = irV1[0], crV11 = irV1[1], crV12 = irV1[2];
        const crV20 = irV2[0], crV21 = irV2[1], crV22 = irV2[2];

        // Evaluate the expression tree of rational numbers.
        const crX0 = crP0.sub(crV00);
        const crY0 = crP1.sub(crV01);
        const crZ0 = crP2.sub(crV02);
        const crX1 = crV10.sub(crV00);
        const crY1 = crV11.sub(crV01);
        const crZ1 = crV12.sub(crV02);
        const crX2 = crV20.sub(crV00);
        const crY2 = crV21.sub(crV01);
        const crZ2 = crV22.sub(crV02);
        const crY0Z1 = crY0.mul(crZ1);
        const crY0Z2 = crY0.mul(crZ2);
        const crY1Z0 = crY1.mul(crZ0);
        const crY1Z2 = crY1.mul(crZ2);
        const crY2Z0 = crY2.mul(crZ0);
        const crY2Z1 = crY2.mul(crZ1);
        const crC0 = crY1Z2.sub(crY2Z1);
        const crC1 = crY2Z0.sub(crY0Z2);
        const crC2 = crY0Z1.sub(crY1Z0);
        const crX0C0 = crX0.mul(crC0);
        const crX1C1 = crX1.mul(crC1);
        const crX2C2 = crX2.mul(crC2);
        const crDet = crX0C0.add(crX1C1).add(crX2C2);
        return crDet.getSign();
    }

    // For a tetrahedron with vertices ordered as described in TetrahedronKey,
    // toCircumsphere returns
    //   +1, P outside circumsphere of tetrahedron
    //   -1, P inside circumsphere of tetrahedron
    //    0, P on circumsphere of tetrahedron
    protected toCircumsphere(pIndex: number, v0Index: number, v1Index: number,
        v2Index: number, v3Index: number): number {
        // The expression tree has 98 nodes consisting of 15 input leaves and
        // 83 compute nodes.

        // Use interval arithmetic to determine the sign if possible.
        const inP = this.inPoint(pIndex);
        const inV0 = this.mVertices[v0Index].values;
        const inV1 = this.mVertices[v1Index].values;
        const inV2 = this.mVertices[v2Index].values;
        const inV3 = this.mVertices[v3Index].values;

        const x0 = SWInterval.sub(inV0[0], inP[0]);
        const y0 = SWInterval.sub(inV0[1], inP[1]);
        const z0 = SWInterval.sub(inV0[2], inP[2]);
        const s00 = SWInterval.add(inV0[0], inP[0]);
        const s01 = SWInterval.add(inV0[1], inP[1]);
        const s02 = SWInterval.add(inV0[2], inP[2]);
        const x1 = SWInterval.sub(inV1[0], inP[0]);
        const y1 = SWInterval.sub(inV1[1], inP[1]);
        const z1 = SWInterval.sub(inV1[2], inP[2]);
        const s10 = SWInterval.add(inV1[0], inP[0]);
        const s11 = SWInterval.add(inV1[1], inP[1]);
        const s12 = SWInterval.add(inV1[2], inP[2]);
        const x2 = SWInterval.sub(inV2[0], inP[0]);
        const y2 = SWInterval.sub(inV2[1], inP[1]);
        const z2 = SWInterval.sub(inV2[2], inP[2]);
        const s20 = SWInterval.add(inV2[0], inP[0]);
        const s21 = SWInterval.add(inV2[1], inP[1]);
        const s22 = SWInterval.add(inV2[2], inP[2]);
        const x3 = SWInterval.sub(inV3[0], inP[0]);
        const y3 = SWInterval.sub(inV3[1], inP[1]);
        const z3 = SWInterval.sub(inV3[2], inP[2]);
        const s30 = SWInterval.add(inV3[0], inP[0]);
        const s31 = SWInterval.add(inV3[1], inP[1]);
        const s32 = SWInterval.add(inV3[2], inP[2]);
        const t00 = s00.mul(x0);
        const t01 = s01.mul(y0);
        const t02 = s02.mul(z0);
        const t10 = s10.mul(x1);
        const t11 = s11.mul(y1);
        const t12 = s12.mul(z1);
        const t20 = s20.mul(x2);
        const t21 = s21.mul(y2);
        const t22 = s22.mul(z2);
        const t30 = s30.mul(x3);
        const t31 = s31.mul(y3);
        const t32 = s32.mul(z3);
        const w0 = t00.add(t01).add(t02);
        const w1 = t10.add(t11).add(t12);
        const w2 = t20.add(t21).add(t22);
        const w3 = t30.add(t31).add(t32);
        const x0y1 = x0.mul(y1);
        const x0y2 = x0.mul(y2);
        const x0y3 = x0.mul(y3);
        const x1y0 = x1.mul(y0);
        const x1y2 = x1.mul(y2);
        const x1y3 = x1.mul(y3);
        const x2y0 = x2.mul(y0);
        const x2y1 = x2.mul(y1);
        const x2y3 = x2.mul(y3);
        const x3y0 = x3.mul(y0);
        const x3y1 = x3.mul(y1);
        const x3y2 = x3.mul(y2);
        const z0w1 = z0.mul(w1);
        const z0w2 = z0.mul(w2);
        const z0w3 = z0.mul(w3);
        const z1w0 = z1.mul(w0);
        const z1w2 = z1.mul(w2);
        const z1w3 = z1.mul(w3);
        const z2w0 = z2.mul(w0);
        const z2w1 = z2.mul(w1);
        const z2w3 = z2.mul(w3);
        const z3w0 = z3.mul(w0);
        const z3w1 = z3.mul(w1);
        const z3w2 = z3.mul(w2);
        const u0 = x0y1.sub(x1y0);
        const u1 = x0y2.sub(x2y0);
        const u2 = x0y3.sub(x3y0);
        const u3 = x1y2.sub(x2y1);
        const u4 = x1y3.sub(x3y1);
        const u5 = x2y3.sub(x3y2);
        const v0 = z0w1.sub(z1w0);
        const v1 = z0w2.sub(z2w0);
        const v2 = z0w3.sub(z3w0);
        const v3 = z1w2.sub(z2w1);
        const v4 = z1w3.sub(z3w1);
        const v5 = z2w3.sub(z3w2);
        const u0v5 = u0.mul(v5);
        const u1v4 = u1.mul(v4);
        const u2v3 = u2.mul(v3);
        const u3v2 = u3.mul(v2);
        const u4v1 = u4.mul(v1);
        const u5v0 = u5.mul(v0);
        const det = u0v5.sub(u1v4).add(u2v3).add(u3v2).sub(u4v1).add(u5v0);

        const zero = 0;
        if (det.get(0) > zero) {
            return +1;
        }
        else if (det.get(1) < zero) {
            return -1;
        }

        // The exact sign of the determinant is not known, so compute the
        // determinant using rational arithmetic.
        const irP = this.irPoint(pIndex);
        const irV0 = this.mIRVertices[v0Index];
        const irV1 = this.mIRVertices[v1Index];
        const irV2 = this.mIRVertices[v2Index];
        const irV3 = this.mIRVertices[v3Index];

        const crP0 = irP[0], crP1 = irP[1], crP2 = irP[2];
        const crV00 = irV0[0], crV01 = irV0[1], crV02 = irV0[2];
        const crV10 = irV1[0], crV11 = irV1[1], crV12 = irV1[2];
        const crV20 = irV2[0], crV21 = irV2[1], crV22 = irV2[2];
        const crV30 = irV3[0], crV31 = irV3[1], crV32 = irV3[2];

        // Evaluate the expression tree of rational numbers.
        const crX0 = crV00.sub(crP0);
        const crY0 = crV01.sub(crP1);
        const crZ0 = crV02.sub(crP2);
        const crS00 = crV00.add(crP0);
        const crS01 = crV01.add(crP1);
        const crS02 = crV02.add(crP2);
        const crX1 = crV10.sub(crP0);
        const crY1 = crV11.sub(crP1);
        const crZ1 = crV12.sub(crP2);
        const crS10 = crV10.add(crP0);
        const crS11 = crV11.add(crP1);
        const crS12 = crV12.add(crP2);
        const crX2 = crV20.sub(crP0);
        const crY2 = crV21.sub(crP1);
        const crZ2 = crV22.sub(crP2);
        const crS20 = crV20.add(crP0);
        const crS21 = crV21.add(crP1);
        const crS22 = crV22.add(crP2);
        const crX3 = crV30.sub(crP0);
        const crY3 = crV31.sub(crP1);
        const crZ3 = crV32.sub(crP2);
        const crS30 = crV30.add(crP0);
        const crS31 = crV31.add(crP1);
        const crS32 = crV32.add(crP2);
        const crT00 = crS00.mul(crX0);
        const crT01 = crS01.mul(crY0);
        const crT02 = crS02.mul(crZ0);
        const crT10 = crS10.mul(crX1);
        const crT11 = crS11.mul(crY1);
        const crT12 = crS12.mul(crZ1);
        const crT20 = crS20.mul(crX2);
        const crT21 = crS21.mul(crY2);
        const crT22 = crS22.mul(crZ2);
        const crT30 = crS30.mul(crX3);
        const crT31 = crS31.mul(crY3);
        const crT32 = crS32.mul(crZ3);
        const crW0 = crT00.add(crT01).add(crT02);
        const crW1 = crT10.add(crT11).add(crT12);
        const crW2 = crT20.add(crT21).add(crT22);
        const crW3 = crT30.add(crT31).add(crT32);
        const crX0Y1 = crX0.mul(crY1);
        const crX0Y2 = crX0.mul(crY2);
        const crX0Y3 = crX0.mul(crY3);
        const crX1Y0 = crX1.mul(crY0);
        const crX1Y2 = crX1.mul(crY2);
        const crX1Y3 = crX1.mul(crY3);
        const crX2Y0 = crX2.mul(crY0);
        const crX2Y1 = crX2.mul(crY1);
        const crX2Y3 = crX2.mul(crY3);
        const crX3Y0 = crX3.mul(crY0);
        const crX3Y1 = crX3.mul(crY1);
        const crX3Y2 = crX3.mul(crY2);
        const crZ0W1 = crZ0.mul(crW1);
        const crZ0W2 = crZ0.mul(crW2);
        const crZ0W3 = crZ0.mul(crW3);
        const crZ1W0 = crZ1.mul(crW0);
        const crZ1W2 = crZ1.mul(crW2);
        const crZ1W3 = crZ1.mul(crW3);
        const crZ2W0 = crZ2.mul(crW0);
        const crZ2W1 = crZ2.mul(crW1);
        const crZ2W3 = crZ2.mul(crW3);
        const crZ3W0 = crZ3.mul(crW0);
        const crZ3W1 = crZ3.mul(crW1);
        const crZ3W2 = crZ3.mul(crW2);
        const crU0 = crX0Y1.sub(crX1Y0);
        const crU1 = crX0Y2.sub(crX2Y0);
        const crU2 = crX0Y3.sub(crX3Y0);
        const crU3 = crX1Y2.sub(crX2Y1);
        const crU4 = crX1Y3.sub(crX3Y1);
        const crU5 = crX2Y3.sub(crX3Y2);
        const crV0 = crZ0W1.sub(crZ1W0);
        const crV1 = crZ0W2.sub(crZ2W0);
        const crV2 = crZ0W3.sub(crZ3W0);
        const crV3 = crZ1W2.sub(crZ2W1);
        const crV4 = crZ1W3.sub(crZ3W1);
        const crV5 = crZ2W3.sub(crZ3W2);
        const crU0V5 = crU0.mul(crV5);
        const crU1V4 = crU1.mul(crV4);
        const crU2V3 = crU2.mul(crV3);
        const crU3V2 = crU3.mul(crV2);
        const crU4V1 = crU4.mul(crV1);
        const crU5V0 = crU5.mul(crV0);
        const crDet = crU0V5.sub(crU1V4).add(crU2V3).add(crU3V2)
            .sub(crU4V1).add(crU5V0);
        return crDet.getSign();
    }

    // The port of the upstream 'bool GetContainingTetrahedron(pIndex, tetra)'
    // with its Tetrahedron*& in/out parameter; the port returns both.
    private getContainingTetrahedronOfVertex(pIndex: number,
        tetra: TSManifoldMeshTetrahedron):
        { found: boolean; tetra: TSManifoldMeshTetrahedron } {
        const opposite = TetrahedronKey.getOppositeFace();
        const numTetrahedra = this.mGraph.getNumTetrahedra();
        for (let t = 0; t < numTetrahedra; ++t) {
            let j: number;
            for (j = 0; j < 4; ++j) {
                const v0Index = tetra.V[opposite[j][0]];
                const v1Index = tetra.V[opposite[j][1]];
                const v2Index = tetra.V[opposite[j][2]];
                if (this.toPlane(pIndex, v0Index, v1Index, v2Index) > 0) {
                    // The point sees face <v0,v1,v2> from outside the
                    // tetrahedron.
                    const adjTetra = tetra.S[j];
                    if (adjTetra) {
                        // Traverse to the tetrahedron sharing the face.
                        tetra = adjTetra;
                        break;
                    }
                    else {
                        // We reached a hull face, so the point is outside the
                        // hull.
                        return { found: false, tetra };
                    }
                }
            }

            if (j === 4) {
                // The point is inside all four faces, so the point is inside
                // a tetrahedron.
                return { found: true, tetra };
            }
        }

        logError('Unexpected termination of loop while searching for a tetrahedron.');
        return { found: false, tetra };
    }

    private getAndRemoveInsertionPolyhedron(pIndex: number,
        candidates: Set<TSManifoldMeshTetrahedron>,
        boundary: DirectedTriangleKeySet): void {
        const opposite = TetrahedronKey.getOppositeFace();

        // Locate the tetrahedra that make up the insertion polyhedron.
        const polyhedron = new TSManifoldMesh();
        while (candidates.size > 0) {
            const tetra = candidates.values().next().value as TSManifoldMeshTetrahedron;
            candidates.delete(tetra);

            for (let j = 0; j < 4; ++j) {
                const adj = tetra.S[j];
                if (adj && !candidates.has(adj)) {
                    const v0Index = adj.V[0];
                    const v1Index = adj.V[1];
                    const v2Index = adj.V[2];
                    const v3Index = adj.V[3];
                    if (this.toCircumsphere(pIndex, v0Index, v1Index, v2Index,
                        v3Index) <= 0) {
                        // Point P is in the circumsphere.
                        candidates.add(adj);
                    }
                }
            }

            const inserted = polyhedron.insert(tetra.V[0], tetra.V[1], tetra.V[2],
                tetra.V[3]);
            logAssert(inserted !== null, 'Unexpected insertion failure.');
            const removed = this.mGraph.remove(tetra.V[0], tetra.V[1], tetra.V[2],
                tetra.V[3]);
            logAssert(removed, 'Unexpected removal failure.');
        }

        // Get the boundary triangles of the insertion polyhedron.
        for (const tetra of polyhedron.getTetrahedra()) {
            for (let j = 0; j < 4; ++j) {
                if (!tetra.S[j]) {
                    boundary.insert(tetra.V[opposite[j][0]], tetra.V[opposite[j][1]],
                        tetra.V[opposite[j][2]]);
                }
            }
        }
    }

    private update(pIndex: number): void {
        const opposite = TetrahedronKey.getOppositeFace();
        const smap = this.mGraph.getTetrahedra();
        const first = smap[0];
        const containing = this.getContainingTetrahedronOfVertex(pIndex, first);
        if (containing.found) {
            // The point is inside the convex hull. The insertion polyhedron
            // contains only tetrahedra in the current tetrahedralization; the
            // hull does not change.

            // Use a depth-first search for those tetrahedra whose
            // circumspheres contain point P.
            const candidates = new Set<TSManifoldMeshTetrahedron>();
            candidates.add(containing.tetra);

            // Get the boundary of the insertion polyhedron C that contains
            // the tetrahedra whose circumspheres contain point P. Polyhedron
            // C contains this point.
            const boundary = new DirectedTriangleKeySet();
            this.getAndRemoveInsertionPolyhedron(pIndex, candidates, boundary);

            // The insertion polyhedron consists of the tetrahedra formed by
            // point P and the faces of C.
            for (const key of boundary.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                const v2Index = key[2];
                if (this.toPlane(pIndex, v0Index, v1Index, v2Index) < 0) {
                    const inserted = this.mGraph.insert(pIndex, v0Index, v1Index,
                        v2Index);
                    logAssert(inserted !== null, 'Unexpected insertion failure.');
                }
                // else: Point P is on a face of the boundary, so the
                // subdivision would have degenerate tetrahedra. Ignore these.
            }
        }
        else {
            // The point is outside the convex hull. The insertion polyhedron
            // is formed by point P and any tetrahedra in the current
            // tetrahedralization whose circumspheres contain point P.

            // Locate the convex hull of the tetrahedra.
            const hull = new DirectedTriangleKeySet();
            for (const t of smap) {
                for (let j = 0; j < 4; ++j) {
                    if (!t.S[j]) {
                        hull.insert(t.V[opposite[j][0]], t.V[opposite[j][1]],
                            t.V[opposite[j][2]]);
                    }
                }
            }

            // Iterate over all the hull faces and use the ones visible to
            // point P to locate the insertion polyhedron.
            const candidates = new Set<TSManifoldMeshTetrahedron>();
            const visible = new DirectedTriangleKeySet();
            for (const key of hull.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                const v2Index = key[2];
                if (this.toPlane(pIndex, v0Index, v1Index, v2Index) > 0) {
                    const face = this.mGraph.getTriangle(v0Index, v1Index, v2Index);
                    if (face !== null && face.S[1] === null) {
                        const adj = face.S[0];
                        if (adj && !candidates.has(adj)) {
                            const a0Index = adj.V[0];
                            const a1Index = adj.V[1];
                            const a2Index = adj.V[2];
                            const a3Index = adj.V[3];
                            if (this.toCircumsphere(pIndex, a0Index, a1Index, a2Index,
                                a3Index) <= 0) {
                                // Point P is in the circumsphere.
                                candidates.add(adj);
                            }
                            else {
                                // Point P is not in the circumsphere but the
                                // hull face is visible.
                                visible.insert(key[0], key[1], key[2]);
                            }
                        }
                    }
                    else {
                        logError('This condition should not occur for rational arithmetic.');
                    }
                }
            }

            // Get the boundary of the insertion subpolyhedron C that contains
            // the tetrahedra whose circumspheres contain point P.
            const boundary = new DirectedTriangleKeySet();
            this.getAndRemoveInsertionPolyhedron(pIndex, candidates, boundary);

            // The insertion polyhedron P consists of the tetrahedra formed by
            // point P and the back faces of C *and* the visible faces of
            // mGraph-C.
            for (const key of boundary.keys()) {
                const v0Index = key[0];
                const v1Index = key[1];
                const v2Index = key[2];
                if (this.toPlane(pIndex, v0Index, v1Index, v2Index) < 0) {
                    // This is a back face of the boundary.
                    const inserted = this.mGraph.insert(pIndex, v0Index, v1Index,
                        v2Index);
                    logAssert(inserted !== null, 'Unexpected insertion failure.');
                }
            }
            for (const key of visible.keys()) {
                const inserted = this.mGraph.insert(pIndex, key[0], key[2], key[1]);
                logAssert(inserted !== null, 'Unexpected insertion failure.');
            }
        }
    }
}
