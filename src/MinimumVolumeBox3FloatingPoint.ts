// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MinimumVolumeBox3FloatingPoint.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The port of the partial specialization
//   MinimumVolumeBox3<T, IndexType, MVB3FloatingPoint>
// See MinimumVolumeBox3.ts for the documentation of the algorithm and the
// queries. The floating-point pipeline computes the candidate box axes with
// 'number' arithmetic but uses exact (BSNumber) arithmetic for the topology
// decisions that must be exact -- the coplanarity and colinearity tests --
// and exact (BSRational) arithmetic for the final box and volume.

import { Vector, normalize } from './Vector.js';
import { cross, computeOrthogonalComplement3 } from './Vector3.js';
import { OrientedBox, type OrientedBox3 } from './OrientedBox.js';
import { ConvexHull3 } from './ConvexHull3.js';
import { MinimumAreaBox2 } from './MinimumAreaBox2.js';
import { VETManifoldMesh } from './VETManifoldMesh.js';
import { UniqueVerticesSimplices } from './UniqueVerticesSimplices.js';
import { BSNumber } from './BSNumber.js';
import { BSRational } from './BSRational.js';
import { logAssert } from './Logger.js';

// The port of MinimumVolumeBox3::invalidIndex. Upstream uses SIZE_MAX as the
// invalid-index marker for std::size_t members. All valid indices are
// nonnegative, so the port uses -1.
const invalidIndex = -1;

// An edge of the polytope mesh. The v[] are indices into the vertex array and
// the t[] are indices into the triangle array.
export interface MinimumVolumeBox3FloatingPointEdge {
    v: [number, number];
    t: [number, number];
}

// A triangle of the polytope mesh. The v[] are indices into the vertex array,
// the e[] are indices into the edge array and the t[] are indices into the
// triangle array (the adjacent triangles).
export interface MinimumVolumeBox3FloatingPointTriangle {
    v: [number, number, number];
    e: [number, number, number];
    t: [number, number, number];
}

// Information about candidates for the minimum volume box and about that box
// itself.
export interface MinimumVolumeBox3FloatingPointCandidate {
    // Set by processEdgePair.
    edgeIndex: [number, number];
    // Upstream stores copies of the two edges. The port stores references to
    // the immutable mEdges[] elements; the members are never written after
    // extractMeshTopology, so this is behavior preserving.
    edge: [MinimumVolumeBox3FloatingPointEdge, MinimumVolumeBox3FloatingPointEdge];
    N: [Vector, Vector];
    M: [Vector, Vector];
    f00: number;
    f10: number;
    f01: number;
    f11: number;
    levelCurveProcessorIndex: number;

    // Set by pair, minimizerConstantT, minimizerConstantS, minimizerVariableS
    // and minimizerVariableT. The axis[0] and axis[1] are set by the
    // aforementioned functions. The axis[2] is computed by computeVolume.
    axis: [Vector, Vector, Vector];

    // Set by computeVolume.
    minSupportIndex: [number, number, number];
    maxSupportIndex: [number, number, number];
    volume: number;
}

// The output of the queries. The 'dimension' is the dimension of the convex
// hull of the input points; it is 3 for the vertices-and-indices query.
export interface MinimumVolumeBox3FloatingPointResult {
    dimension: number;
    box: OrientedBox3;
    volume: number;
}

// An exact-arithmetic 3-tuple; the port of Vector3<BSNumber<UInteger>>.
type NVector3 = [BSNumber, BSNumber, BSNumber];

function nZero3(): NVector3 {
    return [new BSNumber(), new BSNumber(), new BSNumber()];
}

function nSub3(u: NVector3, v: NVector3): NVector3 {
    return [u[0].sub(v[0]), u[1].sub(v[1]), u[2].sub(v[2])];
}

function nCross3(u: NVector3, v: NVector3): NVector3 {
    return [
        u[1].mul(v[2]).sub(u[2].mul(v[1])),
        u[2].mul(v[0]).sub(u[0].mul(v[2])),
        u[0].mul(v[1]).sub(u[1].mul(v[0]))
    ];
}

function nIsZero3(u: NVector3): boolean {
    return u[0].getSign() === 0 && u[1].getSign() === 0 && u[2].getSign() === 0;
}

function rDot3(u: readonly BSRational[], v: readonly BSRational[]): BSRational {
    return u[0].mul(v[0]).add(u[1].mul(v[1])).add(u[2].mul(v[2]));
}

// Hot-path helpers for 3-tuples stored in the port's Vector.
function dot3(u: Vector, v: Vector): number {
    const a = u.values, b = v.values;
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// The port of 's * u + t * v' for 3-tuples.
function axpby3(s: number, u: Vector, t: number, v: Vector): Vector {
    const a = u.values, b = v.values;
    return Vector.fromArray([
        s * a[0] + t * b[0],
        s * a[1] + t * b[1],
        s * a[2] + t * b[2]
    ]);
}

function makeEdge(): MinimumVolumeBox3FloatingPointEdge {
    return { v: [invalidIndex, invalidIndex], t: [invalidIndex, invalidIndex] };
}

function makeTriangle(): MinimumVolumeBox3FloatingPointTriangle {
    return {
        v: [invalidIndex, invalidIndex, invalidIndex],
        e: [invalidIndex, invalidIndex, invalidIndex],
        t: [invalidIndex, invalidIndex, invalidIndex]
    };
}

function makeCandidate(): MinimumVolumeBox3FloatingPointCandidate {
    return {
        edgeIndex: [invalidIndex, invalidIndex],
        edge: [makeEdge(), makeEdge()],
        N: [Vector.zero(3), Vector.zero(3)],
        M: [Vector.zero(3), Vector.zero(3)],
        f00: 0,
        f10: 0,
        f01: 0,
        f11: 0,
        levelCurveProcessorIndex: invalidIndex,
        axis: [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)],
        minSupportIndex: [invalidIndex, invalidIndex, invalidIndex],
        maxSupportIndex: [invalidIndex, invalidIndex, invalidIndex],
        volume: 0
    };
}

// The port of C++ candidate assignment (value semantics). The vectors are
// cloned; the edges are shared because they are immutable.
function assignCandidate(dst: MinimumVolumeBox3FloatingPointCandidate,
    src: MinimumVolumeBox3FloatingPointCandidate): void {
    dst.edgeIndex = [src.edgeIndex[0], src.edgeIndex[1]];
    dst.edge = [src.edge[0], src.edge[1]];
    dst.N = [src.N[0].clone(), src.N[1].clone()];
    dst.M = [src.M[0].clone(), src.M[1].clone()];
    dst.f00 = src.f00;
    dst.f10 = src.f10;
    dst.f01 = src.f01;
    dst.f11 = src.f11;
    dst.levelCurveProcessorIndex = src.levelCurveProcessorIndex;
    dst.axis = [src.axis[0].clone(), src.axis[1].clone(), src.axis[2].clone()];
    dst.minSupportIndex = [src.minSupportIndex[0], src.minSupportIndex[1], src.minSupportIndex[2]];
    dst.maxSupportIndex = [src.maxSupportIndex[0], src.maxSupportIndex[1], src.maxSupportIndex[2]];
    dst.volume = src.volume;
}

function cloneCandidate(src: MinimumVolumeBox3FloatingPointCandidate):
    MinimumVolumeBox3FloatingPointCandidate {
    const dst = makeCandidate();
    assignCandidate(dst, src);
    return dst;
}

// The type of a level-curve processor. Upstream uses pointers to member
// functions; the port uses functions invoked with 'call' so that a derived
// class sees its own overridden minimizers.
type LevelCurveProcessor = (this: MinimumVolumeBox3FloatingPoint,
    c: MinimumVolumeBox3FloatingPointCandidate,
    mvc: MinimumVolumeBox3FloatingPointCandidate) => void;

export class MinimumVolumeBox3FloatingPoint {
    // The number of threads requested by the caller. The port is
    // single-threaded, so the value is accepted and ignored; it is stored
    // only to report it back. This follows the precedent set by the
    // ConvexHull3 and GenerateMeshUV ports.
    protected mNumThreads: number;

    // A mesh representation of the polytope. These members store topological
    // information.
    protected mEdges: MinimumVolumeBox3FloatingPointEdge[];
    protected mEdgeIndices: Array<[number, number]>;
    protected mTriangles: MinimumVolumeBox3FloatingPointTriangle[];
    protected mAdjacentPool: number[];
    protected mAdjacentPoolLocation: number[];
    protected mVClimbStart: number;

    // These members store geometric information.
    protected mTVertices: Vector[];
    protected mTNormals: Vector[];
    protected mTOrigin: Vector;
    protected mNVertices: NVector3[];
    protected mNNormals: NVector3[];
    protected mNOrigin: NVector3;

    // The axis-aligned bounding box of the vertices is used as the initial
    // candidate for the minimum-volume box.
    protected mAlignedCandidate: MinimumVolumeBox3FloatingPointCandidate;

    // The information for the minimum-volume bounding box of the vertices.
    protected mMinimumVolumeObject: MinimumVolumeBox3FloatingPointCandidate;

    // The maximum sample index used to search each level curve for
    // non-face-supporting boxes (mMaxSample + 1 values). The samples are
    // visited using subdivision of the domain of the level curve. The
    // subdivision information is stored in mDomainIndex(mNumSamples-1).
    protected mMaxSample: number;
    protected mDomainIndex: Array<[number, number, number]>;

    // Each member function a00B10C01D11(*) corresponds to a bilinear function
    // on the domain [0,1]^2. Each corner of the domain has a bilinear
    // function value that is positive, negative or zero, leading to 3^4 = 81
    // possibilities. The 'A', 'B', 'C' and 'D' are in {'P', 'M', 'Z'} [for
    // Plus, Minus, Zero].
    protected mLevelCurveProcessor: Array<LevelCurveProcessor | null>;

    // Construction. To execute in the main thread, set numThreads to 0. The
    // port ignores positive values (it is single-threaded) but accepts them
    // so that calling code is portable.
    constructor(numThreads: number = 0) {
        this.mNumThreads = numThreads;
        this.mEdges = [];
        this.mEdgeIndices = [];
        this.mTriangles = [];
        this.mAdjacentPool = [];
        this.mAdjacentPoolLocation = [];
        this.mVClimbStart = 0;
        this.mTVertices = [];
        this.mTNormals = [];
        this.mTOrigin = Vector.zero(3);
        this.mNVertices = [];
        this.mNNormals = [];
        this.mNOrigin = nZero3();
        this.mAlignedCandidate = makeCandidate();
        this.mMinimumVolumeObject = makeCandidate();
        this.mMaxSample = 0;
        this.mDomainIndex = [];
        this.mLevelCurveProcessor = new Array<LevelCurveProcessor | null>(256).fill(null);
        this.initializeLevelCurveProcessors();
    }

    // The number of threads requested at construction.
    getNumThreads(): number {
        return this.mNumThreads;
    }

    // Compute the minimum volume box for an arbitrary set of points. The
    // returned 'dimension' is the dimension of the convex hull of the points.
    // This is the port of operator()(numPoints, points, lgMaxSample, box,
    // volume) and of its std::vector overload.
    compute(points: readonly Vector[], lgMaxSample: number):
        MinimumVolumeBox3FloatingPointResult {
        // The vertices must be those for a 3-dimensional polytope. The
        // smallest such polytope is a tetrahedron, so there must be at least
        // 4 vertices and 4 triangles. The number of samples must be at least
        // 4.
        logAssert(points.length >= 4 && lgMaxSample >= 2, 'Invalid argument.');

        const hull = this.computeConvexHull(points);
        if (hull.dimension === 3) {
            // Compute the minimum volume box for the 3D convex hull.
            const result = this.computeHull(hull.vertices, hull.indices, lgMaxSample);
            return { dimension: 3, box: result.box, volume: result.volume };
        }

        // else: The minimum volume box has volume zero. The number of
        // zero-valued box extents is 3-dimension.
        return { dimension: hull.dimension, box: hull.box, volume: hull.volume };
    }

    // The minimum volume box algorithm involves processing hyperbolic curves.
    // Each curve has a corresponding parameterization, Volume(s), and the
    // global minimum--or an approximation to it--must be computed. The
    // default method used for minimization is to compute sample points along
    // the curve and choosing that point which provides the minimum among all
    // samples. The number of samples is 2^{lgMaxSample} for lgMaxSample >= 2;
    // this implies there are at least 4 samples. You can override the
    // minimizer functions to use your own minimization algorithm; see the
    // comments before the member function minimizerConstantT.
    //
    // The output is the minimum volume box and its volume, although
    // floating-point rounding errors can lead to a result that is nearly the
    // minimum volume box.
    //
    // This is the port of operator()(numVertices, vertices, numIndices,
    // indices, lgMaxSample, box, volume) and of its std::vector overload.
    computeHull(vertices: readonly Vector[], indices: readonly number[],
        lgMaxSample: number): { box: OrientedBox3, volume: number } {
        // The vertices must be those for a 3-dimensional polytope. The
        // smallest such polytope is a tetrahedron, so there must be at least
        // 4 vertices and 4 triangles. The number of samples must be at least
        // 4.
        logAssert(
            vertices.length >= 4 &&
            indices.length >= 12 && (indices.length % 3) === 0 &&
            lgMaxSample >= 2,
            'Invalid argument.');

        // Generate the 2^{lgMaxSample} sample points for minimizing the
        // volume along hyperbolic curves.
        this.generateSubdivision(lgMaxSample);

        // Create a vertex-edge-triangle graph and extract the topological and
        // geometric information from it.
        const numTriangles = indices.length / 3;
        const mesh = new VETManifoldMesh();
        this.createMeshTopology(numTriangles, indices, mesh);
        this.extractMeshTopology(mesh);
        this.extractVertexAdjacencies(mesh);
        this.extractMeshGeometry(vertices);

        // Given a mesh vertex V, the link polygon is a nonplanar and closed
        // polyline. Its vertices are immediately adjacent to V. Remove those
        // V for which the link polygon and V are coplanar. Effectively, this
        // is an implicit way to merge all triangles that are coplanar and
        // form a non-triangle face of the convex hull. The getExtreme
        // function is a hill-climbing algorithm and is the bottleneck in the
        // MVB3 algorithm. Removing the aforementioned V and disabling the
        // relevant adjacent links makes getExtreme simple to implement. If
        // they are not removed, the logic for getExtreme is more complicated
        // and leads to a massive performance loss for large point sets.
        this.removeCoplanarTriangleAdjacencies();

        // Start the search over pairs of normal vectors for the configuration
        // that leads to the minimum volume box.
        this.computeAlignedCandidate();
        this.getMinimumVolumeCandidate();
        return this.getMinimumVolumeBox();
    }

    // The port of ComputeConvexHull. Upstream writes the hull vertices and
    // indices through reference parameters and writes 'box' and 'volume' for
    // hull dimensions 0, 1 and 2.
    protected computeConvexHull(points: readonly Vector[]):
        { dimension: number, vertices: Vector[], indices: number[], box: OrientedBox3, volume: number } {
        const ch3 = new ConvexHull3();
        ch3.compute(points);
        const dimension = ch3.getDimension();
        const hull = ch3.getHull();
        const box = new OrientedBox(3);
        let hullVertices: Vector[] = [];
        let hullIndices: number[] = [];
        let volume = 0;

        if (dimension === 0) {
            // The points are all the same.
            box.center = points[hull[0]].clone();
            box.axis[0] = Vector.fromArray([1, 0, 0]);
            box.axis[1] = Vector.fromArray([0, 1, 0]);
            box.axis[2] = Vector.fromArray([0, 0, 1]);
            box.extent = Vector.zero(3);
            volume = 0;
        } else if (dimension === 1) {
            // The points lie on a line.
            const p0 = points[hull[0]], p1 = points[hull[1]];
            const direction = Vector.fromArray([
                p1.values[0] - p0.values[0],
                p1.values[1] - p0.values[1],
                p1.values[2] - p0.values[2]
            ]);
            box.center = axpby3(0.5, p0, 0.5, p1);
            box.extent = Vector.fromArray([0.5 * normalize(direction), 0, 0]);
            box.axis[0] = direction;
            computeOrthogonalComplement3(1, box.axis);
            volume = 0;
        } else if (dimension === 2) {
            // The points lie on a plane. Get a coordinate system relative to
            // the plane of the points. Choose the origin to be any of the
            // input points.
            const origin = points[hull[0]];
            let normal = Vector.zero(3);
            const numHull = hull.length;
            // Upstream bug (MinimumVolumeBox3FloatingPoint.h,
            // ComputeConvexHull, dimension == 2): the Newell loop is written
            // as 'for (i0 = numHull - 1, i1 = 1; i1 < numHull; i0 = i1++)',
            // so hull[0] is skipped and the polygon edge list is not closed.
            // For numHull == 3 the sum is Cross(P2,P1) + Cross(P1,P2) = 0,
            // producing a zero normal and, therefore, a garbage basis for
            // every coplanar point set whose hull is a triangle. The port
            // starts the loop at i1 = 0, which closes the polygon.
            for (let i0 = numHull - 1, i1 = 0; i1 < numHull; i0 = i1++) {
                const q0 = points[hull[i0]], q1 = points[hull[i1]];
                const c = cross(q0, q1);
                normal = Vector.fromArray([
                    normal.values[0] + c.values[0],
                    normal.values[1] + c.values[1],
                    normal.values[2] + c.values[2]
                ]);
            }

            const basis: Vector[] = [normal, Vector.zero(3), Vector.zero(3)];
            computeOrthogonalComplement3(1, basis);

            // Project the input points onto the plane.
            const projection: Vector[] = new Array<Vector>(points.length);
            for (let i = 0; i < points.length; ++i) {
                const p = points[i].values, o = origin.values;
                const diff = Vector.fromArray([p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
                projection[i] = Vector.fromArray([dot3(basis[1], diff), dot3(basis[2], diff)]);
            }

            // Compute the minimum area box in 2D.
            const mab2 = new MinimumAreaBox2();
            const rectangle = mab2.compute(projection);

            // Lift the values into 3D.
            const rc = rectangle.center.values;
            box.center = Vector.fromArray([
                origin.values[0] + rc[0] * basis[1].values[0] + rc[1] * basis[2].values[0],
                origin.values[1] + rc[0] * basis[1].values[1] + rc[1] * basis[2].values[1],
                origin.values[2] + rc[0] * basis[1].values[2] + rc[1] * basis[2].values[2]
            ]);
            box.axis[0] = axpby3(rectangle.axis[0].values[0], basis[1],
                rectangle.axis[0].values[1], basis[2]);
            box.axis[1] = axpby3(rectangle.axis[1].values[0], basis[1],
                rectangle.axis[1].values[1], basis[2]);
            box.axis[2] = basis[0].clone();
            box.extent = Vector.fromArray([
                rectangle.extent.values[0], rectangle.extent.values[1], 0]);
            volume = 0;
        } else {
            // dimension == 3. Remove duplicated vertices and reindex them for
            // the polytope.
            const sourceVertices: Vector[] = points.map(p => p.clone());
            const sourceIndices: number[] = Array.from(hull);
            const uvt = new UniqueVerticesSimplices<Vector>(3);
            const output = uvt.removeDuplicateAndUnusedVertices(sourceVertices, sourceIndices);
            hullVertices = output.vertices;
            hullIndices = output.indices;
        }

        return { dimension, vertices: hullVertices, indices: hullIndices, box, volume };
    }

    protected createMeshTopology(numTriangles: number, indices: readonly number[],
        mesh: VETManifoldMesh): void {
        let current = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const v0 = indices[current++];
            const v1 = indices[current++];
            const v2 = indices[current++];
            mesh.insert(v0, v1, v2);
        }

        const numV = mesh.getNumVertices();
        const numE = mesh.getNumEdges();
        const numT = mesh.getNumTriangles();
        this.mEdges = new Array<MinimumVolumeBox3FloatingPointEdge>(numE);
        for (let i = 0; i < numE; ++i) {
            this.mEdges[i] = makeEdge();
        }
        this.mEdgeIndices = [];
        this.mTriangles = new Array<MinimumVolumeBox3FloatingPointTriangle>(numT);
        for (let i = 0; i < numT; ++i) {
            this.mTriangles[i] = makeTriangle();
        }
        this.mTVertices = new Array<Vector>(numV);
        this.mTNormals = new Array<Vector>(numT);
        this.mNVertices = new Array<NVector3>(numV);
        this.mNNormals = new Array<NVector3>(numT);
    }

    protected extractMeshTopology(mesh: VETManifoldMesh): void {
        // The port's getEdges() and getTriangles() return the mesh features
        // in increasing feature-key order, matching the std::map iteration
        // order of the upstream eMap and tMap.
        const eMap = mesh.getEdges();
        const tMap = mesh.getTriangles();

        const edgeIndexMap = new Map<object, number>();
        for (let index = 0; index < eMap.length; ++index) {
            const element = eMap[index];
            edgeIndexMap.set(element, index);
            for (let j = 0; j < 2; ++j) {
                this.mEdges[index].v[j] = element.V[j];
            }
        }

        const triangleIndexMap = new Map<object, number>();
        for (let index = 0; index < tMap.length; ++index) {
            const element = tMap[index];
            triangleIndexMap.set(element, index);
            for (let j = 0; j < 3; ++j) {
                this.mTriangles[index].v[j] = element.V[j];
            }
        }

        for (let index = 0; index < eMap.length; ++index) {
            const element = eMap[index];
            for (let j = 0; j < 2; ++j) {
                const triangle = element.T[j];
                logAssert(triangle !== null, 'Unexpected condition: the mesh must be closed.');
                this.mEdges[index].t[j] = triangleIndexMap.get(triangle) as number;
            }
        }

        for (let index = 0; index < tMap.length; ++index) {
            const element = tMap[index];
            for (let j = 0; j < 3; ++j) {
                const edge = element.E[j];
                logAssert(edge !== null, 'Unexpected condition: missing triangle edge.');
                this.mTriangles[index].e[j] = edgeIndexMap.get(edge) as number;
            }
            for (let j = 0; j < 3; ++j) {
                const triangle = element.T[j];
                logAssert(triangle !== null, 'Unexpected condition: the mesh must be closed.');
                this.mTriangles[index].t[j] = triangleIndexMap.get(triangle) as number;
            }
        }

        for (let e0 = 0; e0 < this.mEdges.length; ++e0) {
            for (let e1 = e0 + 1; e1 < this.mEdges.length; ++e1) {
                this.mEdgeIndices.push([e0, e1]);
            }
        }
    }

    protected extractVertexAdjacencies(mesh: VETManifoldMesh): void {
        // The vertices are stored in a vertex-edge-triangle manifold mesh.
        // Each vertex has a set of adjacent vertices, a set of adjacent edges
        // and a set of adjacent triangles. The adjacent vertices are
        // repackaged into mAdjacentPoolLocation[] and mAdjacentPool[]. For
        // vertex v with n adjacent vertices, mAdjacentPoolLocation[v] is the
        // index into mAdjacentPool[] where the n adjacent vertices are
        // stored. If the adjacent vertices are a[0] through a[n-1], then
        // mAdjacentPool[mAdjacentPoolLocation[v] + i] is a[i] for 0 <= i < n.
        //
        // In the construction of mAdjacentPoolLocation[v]:
        //   (1) the vertex indices v satisfy 0 <= v < N for a mesh of N
        //       vertices and
        //   (2) the vertex map itself is ordered as <0,vertex0>,
        //       <1,vertex1>, ..., <N-1,vertexNm1>.
        // Condition (1) is guaranteed because the input uses the contiguous
        // indices of the position array. The port's getVertices() returns the
        // vertices sorted by vertex index, which satisfies condition (2).
        const sortedVMap = mesh.getVertices();

        let numAdjacentPool = 0;
        for (const element of sortedVMap) {
            numAdjacentPool += element.VAdjacent.size + 1;
        }
        this.mAdjacentPool = new Array<number>(numAdjacentPool).fill(invalidIndex);
        this.mAdjacentPoolLocation = new Array<number>(sortedVMap.length).fill(invalidIndex);
        let apIndex = 0, vaIndex = 0;
        for (const element of sortedVMap) {
            logAssert(element.V === vaIndex,
                'Invalid argument: the vertex indices must be contiguous and start at zero.');
            // The adjacent indices are visited in increasing order, matching
            // the std::set iteration order of the upstream VAdjacent.
            const adjacent = element.getVAdjacent();
            this.mAdjacentPoolLocation[vaIndex++] = apIndex;
            this.mAdjacentPool[apIndex++] = adjacent.length;
            for (const v of adjacent) {
                this.mAdjacentPool[apIndex++] = v;
            }
        }
    }

    protected extractMeshGeometry(vertices: readonly Vector[]): void {
        const numVertices = vertices.length;
        logAssert(numVertices === this.mTVertices.length,
            'Invalid argument: the mesh must use all the vertices exactly once.');

        // Translate the polytope so that vertices[0] becomes the origin. This
        // helps avoid large floating-point rounding errors when the polytope
        // is far away from (0,0,0).
        this.mTOrigin = vertices[0].clone();
        for (let j = 0; j < 3; ++j) {
            this.mNOrigin[j] = BSNumber.fromNumber(vertices[0].values[j]);
        }
        this.mTVertices[0] = Vector.zero(3);
        this.mNVertices[0] = nZero3();
        for (let i = 1; i < numVertices; ++i) {
            const p = vertices[i].values, o = this.mTOrigin.values;
            this.mTVertices[i] = Vector.fromArray([p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
            this.mNVertices[i] = [
                BSNumber.fromNumber(p[0]).sub(this.mNOrigin[0]),
                BSNumber.fromNumber(p[1]).sub(this.mNOrigin[1]),
                BSNumber.fromNumber(p[2]).sub(this.mNOrigin[2])
            ];
        }

        // Create the triangles and normals to the triangles. The normals
        // mNNormals[] are not normalized to avoid floating-point rounding
        // errors. This is necessary for creating the vertex adjacency data
        // structure that supports the hill-climbing algorithm that computes
        // extreme hull points in a specified direction. The mTNormals[] are
        // normalized so that the floating-point-valued dot products in the
        // hill climbing are computed to reduce floating-point rounding
        // errors.
        for (let i = 0; i < this.mTriangles.length; ++i) {
            const tri = this.mTriangles[i];
            const v0 = tri.v[0], v1 = tri.v[1], v2 = tri.v[2];
            const edge10 = nSub3(this.mNVertices[v1], this.mNVertices[v0]);
            const edge20 = nSub3(this.mNVertices[v2], this.mNVertices[v0]);
            this.mNNormals[i] = nCross3(edge20, edge10);
            const normal = Vector.fromArray([
                this.mNNormals[i][0].toNumber(),
                this.mNNormals[i][1].toNumber(),
                this.mNNormals[i][2].toNumber()
            ]);
            normalize(normal);
            this.mTNormals[i] = normal;
        }
    }

    protected insertAdjacent(vertex: number, insertionCandidate: number): void {
        const base = this.mAdjacentPoolLocation[vertex];
        const numAdjacent = this.mAdjacentPool[base] + 1;
        this.mAdjacentPool[base] = numAdjacent;
        this.mAdjacentPool[base + numAdjacent] = insertionCandidate;
    }

    protected removeAdjacent(vertex: number, removalCandidate: number): void {
        const base = this.mAdjacentPoolLocation[vertex];
        const numAdjacent = this.mAdjacentPool[base];
        for (let j = 1; j <= numAdjacent; ++j) {
            if (this.mAdjacentPool[base + j] === removalCandidate) {
                // The adjacent candidate is indeed adjacent to the vertex, so
                // remove it. To maintain a contiguous array of adjacents,
                // move the last element of the array to the location vacated
                // by the adjacent candidate. If the vacated location is
                // already the end of the array, there is nothing to move.
                if (j < numAdjacent) {
                    this.mAdjacentPool[base + j] = this.mAdjacentPool[base + numAdjacent];
                }
                this.mAdjacentPool[base + numAdjacent] = invalidIndex;
                this.mAdjacentPool[base] = numAdjacent - 1;
                return;
            }
        }
    }

    protected removeCoplanarTriangleAdjacencies(): void {
        // Adjacent triangles are coplanar if their unit-length normal vectors
        // are equal. For such triangles, the winding order of the triangles
        // in the manifold mesh guarantees the normals point in the same
        // direction; that is, we cannot have N1 = -N0.
        for (const edge of this.mEdges) {
            const N0 = this.mNNormals[edge.t[0]];
            const N1 = this.mNNormals[edge.t[1]];
            if (nIsZero3(nCross3(N0, N1))) {
                // The triangles sharing the edge are coplanar. Remove the
                // vertex-adjacent information for the edge vertices. This
                // leads to an implied removal of coplanar triangles which
                // then makes the getExtreme hill-climbing algorithm simple to
                // implement by not having to keep track of the bookkeeping
                // while traversing a patch of coplanar vertices.
                this.removeAdjacent(edge.v[0], edge.v[1]);
                this.removeAdjacent(edge.v[1], edge.v[0]);
            }
        }

        // After removing interior edges of a coplanar triangle face, the
        // boundary edges of the face can have colinear vertices. These
        // vertices must be removed so that the face becomes a convex polygon
        // with no colinear vertices.
        for (let v = 0; v < this.mNVertices.length; ++v) {
            const base = this.mAdjacentPoolLocation[v];
            if (this.mAdjacentPool[base] === 2) {
                // Test for colinearity.
                const vPrev = this.mAdjacentPool[base + 1];
                const vNext = this.mAdjacentPool[base + 2];
                const diff0 = nSub3(this.mNVertices[v], this.mNVertices[vPrev]);
                const diff1 = nSub3(this.mNVertices[v], this.mNVertices[vNext]);
                if (nIsZero3(nCross3(diff0, diff1))) {
                    // The points are colinear. Remove the middle point.
                    this.removeAdjacent(v, vPrev);
                    this.removeAdjacent(vPrev, v);
                    this.removeAdjacent(v, vNext);
                    this.removeAdjacent(vNext, v);
                    this.mAdjacentPool[base] = 0;

                    // The endpoints are now adjacent.
                    this.insertAdjacent(vPrev, vNext);
                    this.insertAdjacent(vNext, vPrev);
                }
            }
        }

        // Locate the first nonempty adjacency list and use it to set the
        // initial index into mAdjacentPoolLocation[] for the hill climbing.
        this.mVClimbStart = invalidIndex;
        for (let i = 0; i < this.mNVertices.length; ++i) {
            const numAdjacent = this.mAdjacentPool[this.mAdjacentPoolLocation[i]];
            if (numAdjacent > 0) {
                this.mVClimbStart = i;
                break;
            }
        }

        logAssert(this.mVClimbStart !== invalidIndex,
            'Unexpected condition: At least one adjacency list should be nonempty.');
    }

    protected computeAlignedCandidate(): void {
        const pmin = [0, 0, 0], pmax = [0, 0, 0];
        for (let j = 0; j < 3; ++j) {
            const axis = this.mAlignedCandidate.axis[j];
            const emax = this.getExtreme(axis);
            this.mAlignedCandidate.maxSupportIndex[j] = emax.vMax;
            pmax[j] = emax.dMax;
            const emin = this.getExtreme(Vector.fromArray([
                -axis.values[0], -axis.values[1], -axis.values[2]]));
            this.mAlignedCandidate.minSupportIndex[j] = emin.vMax;
            pmin[j] = -emin.dMax;
        }
        this.mAlignedCandidate.volume =
            (pmax[0] - pmin[0]) * (pmax[1] - pmin[1]) * (pmax[2] - pmin[2]);
    }

    // The hill-climbing algorithm that computes the extreme hull vertex in
    // the specified direction. Upstream returns the vertex index and writes
    // the extreme dot product through a reference parameter.
    protected getExtreme(direction: Vector): { vMax: number, dMax: number } {
        let vMax = this.mVClimbStart;
        let dMax = dot3(direction, this.mTVertices[vMax]);

        for (let i = 0; i < this.mTVertices.length; ++i) {
            let vLocalMax = vMax;
            let dLocalMax = dMax;
            const base = this.mAdjacentPoolLocation[vMax];
            const numAdjacent = this.mAdjacentPool[base];
            for (let j = 1; j <= numAdjacent; ++j) {
                const vCandidate = this.mAdjacentPool[base + j];
                const dCandidate = dot3(direction, this.mTVertices[vCandidate]);
                if (dCandidate > dLocalMax) {
                    vLocalMax = vCandidate;
                    dLocalMax = dCandidate;
                }
            }
            if (vMax !== vLocalMax) {
                vMax = vLocalMax;
                dMax = dLocalMax;
            } else {
                break;
            }
        }

        return { vMax, dMax };
    }

    protected computeVolume(candidate: MinimumVolumeBox3FloatingPointCandidate): void {
        // The last axis is needed only when computing the volume for
        // comparison to the current candidate volume, so compute this axis
        // now.
        candidate.axis[2] = cross(candidate.axis[0], candidate.axis[1]);

        const pmin = [0, 0, 0], pmax = [0, 0, 0];
        candidate.minSupportIndex[0] = this.mEdges[candidate.edgeIndex[0]].v[0];
        pmin[0] = dot3(candidate.axis[0], this.mTVertices[candidate.minSupportIndex[0]]);
        const e0 = this.getExtreme(candidate.axis[0]);
        candidate.maxSupportIndex[0] = e0.vMax;
        pmax[0] = e0.dMax;
        candidate.minSupportIndex[1] = this.mEdges[candidate.edgeIndex[1]].v[0];
        pmin[1] = dot3(candidate.axis[1], this.mTVertices[candidate.minSupportIndex[1]]);
        const e1 = this.getExtreme(candidate.axis[1]);
        candidate.maxSupportIndex[1] = e1.vMax;
        pmax[1] = e1.dMax;
        // Upstream recomputes axis[2] here; the value is identical to the one
        // computed at the start of the function. The port keeps the single
        // computation above.
        const a2 = candidate.axis[2].values;
        const e2min = this.getExtreme(Vector.fromArray([-a2[0], -a2[1], -a2[2]]));
        candidate.minSupportIndex[2] = e2min.vMax;
        pmin[2] = -e2min.dMax;
        const e2max = this.getExtreme(candidate.axis[2]);
        candidate.maxSupportIndex[2] = e2max.vMax;
        pmax[2] = e2max.dMax;
        candidate.volume =
            (pmax[0] - pmin[0]) * (pmax[1] - pmin[1]) * (pmax[2] - pmin[2]) /
            dot3(candidate.axis[2], candidate.axis[2]);
    }

    protected processEdgePair(edgeIndex: readonly [number, number],
        mvCandidate: MinimumVolumeBox3FloatingPointCandidate): void {
        // Examine the zero-valued level curves for
        // F(s,t)
        // = Dot((1-s)*edge0.N0 + s*edge0.N1, (1-t)*edge1.N0 + t*edge1.N1)
        // = (1-s)*(1-t)*Dot(edge0.N0,edge1.N0)
        //   + (1-s)*t*Dot(edge0.N0,edge1.N1)
        //   + s*(1-t)*Dot(edge0.N1,edge1.N0)
        //   + s*t*Dot(edge0.N1,edge1.N1)
        // = (1-s)*(1-t)*f00 + (1-s)*t*f01 + s*(1-t)*f10 + s*t*f11
        // = a00 + a10*s + a01*t + a11*s*t
        // = [(a00*a11 - a01*a10) + (a01 + a11*s)*(a10 + a11*t)]/a11
        // where a00 = f00, a10 = f10-f00, a01 = f01-f00 and
        // a11 = f00-f01-f10+f11. Let d = a00*a11 - a01*a10 =
        // f00*f11 - f01*f10. If d = 0, then the level curves are
        // s = -a01/a11 and t = -a10/a11. If d != 0, then the level curves are
        // hyperbolic curves with asymptotes s = -a01/a11 and t = -a10/a11.

        const candidate = cloneCandidate(this.mAlignedCandidate);
        candidate.edgeIndex = [edgeIndex[0], edgeIndex[1]];
        const edge0 = this.mEdges[candidate.edgeIndex[0]];
        const edge1 = this.mEdges[candidate.edgeIndex[1]];
        candidate.edge[0] = edge0;
        candidate.edge[1] = edge1;
        candidate.N[0] = this.mTNormals[edge0.t[0]];
        candidate.N[1] = this.mTNormals[edge0.t[1]];
        candidate.M[0] = this.mTNormals[edge1.t[0]];
        candidate.M[1] = this.mTNormals[edge1.t[1]];
        candidate.f00 = dot3(candidate.N[0], candidate.M[0]);
        candidate.f10 = dot3(candidate.N[1], candidate.M[0]);
        candidate.f01 = dot3(candidate.N[0], candidate.M[1]);
        candidate.f11 = dot3(candidate.N[1], candidate.M[1]);

        const bits00 = (candidate.f00 > 0 ? 1 : (candidate.f00 < 0 ? 2 : 0));
        const bits10 = (candidate.f10 > 0 ? 1 : (candidate.f10 < 0 ? 2 : 0));
        const bits01 = (candidate.f01 > 0 ? 1 : (candidate.f01 < 0 ? 2 : 0));
        const bits11 = (candidate.f11 > 0 ? 1 : (candidate.f11 < 0 ? 2 : 0));
        const index = bits00 | (bits10 << 2) | (bits01 << 4) | (bits11 << 6);
        if (index !== 0x55 && index !== 0xaa) {
            candidate.levelCurveProcessorIndex = index;
            const processor = this.mLevelCurveProcessor[index];
            logAssert(processor !== null, 'Unexpected condition: missing level-curve processor.');
            processor.call(this, candidate, mvCandidate);
        }
    }

    protected getMinimumVolumeCandidate(): void {
        assignCandidate(this.mMinimumVolumeObject, this.mAlignedCandidate);

        // Upstream splits the edge pairs across std::thread objects when
        // mNumThreads > 0. The port processes the pairs in the calling thread
        // regardless; the result is the same because each thread computes an
        // independent minimum that is then merged.
        for (const edgeIndex of this.mEdgeIndices) {
            this.processEdgePair(edgeIndex, this.mMinimumVolumeObject);
        }
    }

    protected getMinimumVolumeBox(): { box: OrientedBox3, volume: number } {
        const mvc = this.mMinimumVolumeObject;
        const box = new OrientedBox(3);

        // Compute the rational-valued box and volume. Convert this to a
        // floating-point-valued box and volume on return.
        const rCenter: BSRational[] = [new BSRational(), new BSRational(), new BSRational()];
        const rPMin: BSRational[] = [new BSRational(), new BSRational(), new BSRational()];
        const rPMax: BSRational[] = [new BSRational(), new BSRational(), new BSRational()];
        const rAxis: BSRational[][] = [[], [], []];
        const rSqrLengthAxis: BSRational[] = [];
        for (let i = 0; i < 3; ++i) {
            rCenter[i] = BSRational.fromBSNumber(this.mNOrigin[i]);

            const rVertexMin: BSRational[] = [];
            const rVertexMax: BSRational[] = [];
            for (let j = 0; j < 3; ++j) {
                rVertexMin[j] = BSRational.fromBSNumber(
                    this.mNVertices[mvc.minSupportIndex[i]][j]);
                rVertexMax[j] = BSRational.fromBSNumber(
                    this.mNVertices[mvc.maxSupportIndex[i]][j]);
                rAxis[i][j] = BSRational.fromNumber(mvc.axis[i].values[j]);
            }
            rSqrLengthAxis[i] = rDot3(rAxis[i], rAxis[i]);

            rPMin[i] = rDot3(rAxis[i], rVertexMin);
            rPMax[i] = rDot3(rAxis[i], rVertexMax);
        }

        const rHalf = BSRational.fromNumber(0.5);
        const rAverage: BSRational[] = [
            rHalf.mul(rPMax[0].add(rPMin[0])),
            rHalf.mul(rPMax[1].add(rPMin[1])),
            rHalf.mul(rPMax[2].add(rPMin[2]))
        ];
        for (let i = 0; i < 3; ++i) {
            const scale = rAverage[i].div(rSqrLengthAxis[i]);
            for (let j = 0; j < 3; ++j) {
                rCenter[j] = rCenter[j].add(scale.mul(rAxis[i][j]));
            }
        }

        const rDifference: BSRational[] = [
            rPMax[0].sub(rPMin[0]),
            rPMax[1].sub(rPMin[1]),
            rPMax[2].sub(rPMin[2])
        ];
        const rScaledExtent: BSRational[] = [
            rHalf.mul(rDifference[0]),
            rHalf.mul(rDifference[1]),
            rHalf.mul(rDifference[2])
        ];
        const rVolume = rDifference[0].mul(rDifference[1]).mul(rDifference[2])
            .div(rSqrLengthAxis[2]);

        // Compute the floating-point-valued box and volume.
        for (let i = 0; i < 3; ++i) {
            box.center.values[i] = rCenter[i].toNumber();
            const length = Math.sqrt(rSqrLengthAxis[i].toNumber());
            for (let j = 0; j < 3; ++j) {
                box.axis[i].values[j] = rAxis[i][j].toNumber() / length;
            }
            box.extent.values[i] = rScaledExtent[i].toNumber() / length;
        }

        return { box, volume: rVolume.toNumber() };
    }

    protected createDomainIndex(current: { value: number }, end0: number, end1: number): void {
        const mid = Math.floor((end0 + end1) / 2);
        if (mid !== end0 && mid !== end1) {
            this.mDomainIndex[current.value++] = [mid, end0, end1];
            this.createDomainIndex(current, end0, mid);
            this.createDomainIndex(current, mid, end1);
        }
    }

    protected generateSubdivision(lgMaxSample: number): void {
        this.mMaxSample = (1 << lgMaxSample);
        this.mDomainIndex = new Array<[number, number, number]>(this.mMaxSample - 1);
        const current = { value: 0 };
        this.createDomainIndex(current, 0, this.mMaxSample);
    }

    // The minimizers for the queries. The default behavior is to use the
    // built-in minimizers that sample the level curves as a simple search for
    // a minimum volume. However, you can override the minimizers in a derived
    // class and provide a more sophisticated algorithm.
    protected minimizerConstantS(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // Upstream writes 'T const half = static_cast<Number>(0.5)' here,
        // where Number is the BSNumber compute type. The conversion back to T
        // is exact, so the value is 0.5 as in the sibling minimizers.
        const half = 0.5;
        const t = new Array<number>(this.mMaxSample + 1).fill(0);
        t[0] = 0;
        t[this.mMaxSample] = 1;
        for (const item of this.mDomainIndex) {
            t[item[0]] = half * (t[item[1]] + t[item[2]]);
        }

        normalize(c.axis[0]);
        for (let i = 0, j = this.mMaxSample; i <= this.mMaxSample; ++i, --j) {
            c.axis[1] = axpby3(t[j], c.M[0], t[i], c.M[1]);
            normalize(c.axis[1]);
            this.computeVolume(c);
            if (c.volume < mvc.volume) {
                assignCandidate(mvc, c);
            }
        }
    }

    protected minimizerConstantT(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        const half = 0.5;
        const s = new Array<number>(this.mMaxSample + 1).fill(0);
        s[0] = 0;
        s[this.mMaxSample] = 1;
        for (const item of this.mDomainIndex) {
            s[item[0]] = half * (s[item[1]] + s[item[2]]);
        }

        normalize(c.axis[1]);
        for (let i = 0, j = this.mMaxSample; i <= this.mMaxSample; ++i, --j) {
            c.axis[0] = axpby3(s[j], c.N[0], s[i], c.N[1]);
            normalize(c.axis[0]);
            this.computeVolume(c);
            if (c.volume < mvc.volume) {
                assignCandidate(mvc, c);
            }
        }
    }

    protected minimizerVariableS(sminNumer: number, smaxNumer: number, sDenom: number,
        c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        const half = 0.5;
        const s = new Array<number>(this.mMaxSample + 1).fill(0);
        const oms = new Array<number>(this.mMaxSample + 1).fill(0);
        s[0] = sminNumer;
        oms[0] = sDenom - sminNumer;
        s[this.mMaxSample] = smaxNumer;
        oms[this.mMaxSample] = sDenom - smaxNumer;
        for (const item of this.mDomainIndex) {
            s[item[0]] = half * (s[item[1]] + s[item[2]]);
            oms[item[0]] = half * (oms[item[1]] + oms[item[2]]);
        }

        for (let i = 0; i <= this.mMaxSample; ++i) {
            c.axis[0] = axpby3(oms[i], c.N[0], s[i], c.N[1]);
            normalize(c.axis[0]);

            const q0 = oms[i] * c.f00 + s[i] * c.f10;
            const q1 = oms[i] * c.f01 + s[i] * c.f11;
            if (q0 > q1) {
                c.axis[1] = axpby3(q0, c.M[1], -q1, c.M[0]);
            } else {
                c.axis[1] = axpby3(q1, c.M[0], -q0, c.M[1]);
            }
            normalize(c.axis[1]);

            this.computeVolume(c);
            if (c.volume < mvc.volume) {
                assignCandidate(mvc, c);
            }
        }
    }

    protected minimizerVariableT(tminNumer: number, tmaxNumer: number, tDenom: number,
        c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        const half = 0.5;
        const t = new Array<number>(this.mMaxSample + 1).fill(0);
        const omt = new Array<number>(this.mMaxSample + 1).fill(0);
        t[0] = tminNumer;
        omt[0] = tDenom - tminNumer;
        t[this.mMaxSample] = tmaxNumer;
        omt[this.mMaxSample] = tDenom - tmaxNumer;
        for (const item of this.mDomainIndex) {
            t[item[0]] = half * (t[item[1]] + t[item[2]]);
            omt[item[0]] = half * (omt[item[1]] + omt[item[2]]);
        }

        for (let i = 0; i <= this.mMaxSample; ++i) {
            const p0 = omt[i] * c.f00 + t[i] * c.f01;
            const p1 = omt[i] * c.f10 + t[i] * c.f11;
            if (p0 > p1) {
                c.axis[0] = axpby3(p0, c.N[1], -p1, c.N[0]);
            } else {
                c.axis[0] = axpby3(p1, c.N[0], -p0, c.N[1]);
            }
            normalize(c.axis[0]);

            c.axis[1] = axpby3(omt[i], c.M[0], t[i], c.M[1]);
            normalize(c.axis[1]);

            this.computeVolume(c);
            if (c.volume < mvc.volume) {
                assignCandidate(mvc, c);
            }
        }
    }

    protected pair(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        this.computeVolume(c);
        if (c.volume < mvc.volume) {
            assignCandidate(mvc, c);
        }
    }

    // Set axis[0] and axis[1] to the given normals and process the pair. The
    // vectors are cloned because the minimizers normalize the axes in place
    // and must not modify the candidate normals.
    private pairOf(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate, i: number, j: number): void {
        c.axis[0] = c.N[i].clone();
        c.axis[1] = c.M[j].clone();
        this.pair(c, mvc);
    }

    protected initializeLevelCurveProcessors(): void {
        const p = this.mLevelCurveProcessor;
        // The methods are looked up on 'this' so that a derived class that
        // overrides a processor is dispatched to.
        const proto = this;
        p[0x00] = proto.z00Z10Z01Z11;
        p[0x01] = proto.p00Z10Z01Z11;
        p[0x02] = proto.m00Z10Z01Z11;
        p[0x04] = proto.z00P10Z01Z11;
        p[0x05] = proto.p00P10Z01Z11;
        p[0x06] = proto.m00P10Z01Z11;
        p[0x08] = proto.z00M10Z01Z11;
        p[0x09] = proto.p00M10Z01Z11;
        p[0x0a] = proto.m00M10Z01Z11;
        p[0x10] = proto.z00Z10P01Z11;
        p[0x11] = proto.p00Z10P01Z11;
        p[0x12] = proto.m00Z10P01Z11;
        p[0x14] = proto.z00P10P01Z11;
        p[0x15] = proto.p00P10P01Z11;
        p[0x16] = proto.m00P10P01Z11;
        p[0x18] = proto.z00M10P01Z11;
        p[0x19] = proto.p00M10P01Z11;
        p[0x1a] = proto.m00M10P01Z11;
        p[0x20] = proto.z00Z10M01Z11;
        p[0x21] = proto.p00Z10M01Z11;
        p[0x22] = proto.m00Z10M01Z11;
        p[0x24] = proto.z00P10M01Z11;
        p[0x25] = proto.p00P10M01Z11;
        p[0x26] = proto.m00P10M01Z11;
        p[0x28] = proto.z00M10M01Z11;
        p[0x29] = proto.p00M10M01Z11;
        p[0x2a] = proto.m00M10M01Z11;
        p[0x40] = proto.z00Z10Z01P11;
        p[0x41] = proto.p00Z10Z01P11;
        p[0x42] = proto.m00Z10Z01P11;
        p[0x44] = proto.z00P10Z01P11;
        p[0x45] = proto.p00P10Z01P11;
        p[0x46] = proto.m00P10Z01P11;
        p[0x48] = proto.z00M10Z01P11;
        p[0x49] = proto.p00M10Z01P11;
        p[0x4a] = proto.m00M10Z01P11;
        p[0x50] = proto.z00Z10P01P11;
        p[0x51] = proto.p00Z10P01P11;
        p[0x52] = proto.m00Z10P01P11;
        p[0x54] = proto.z00P10P01P11;
        p[0x55] = proto.p00P10P01P11;
        p[0x56] = proto.m00P10P01P11;
        p[0x58] = proto.z00M10P01P11;
        p[0x59] = proto.p00M10P01P11;
        p[0x5a] = proto.m00M10P01P11;
        p[0x60] = proto.z00Z10M01P11;
        p[0x61] = proto.p00Z10M01P11;
        p[0x62] = proto.m00Z10M01P11;
        p[0x64] = proto.z00P10M01P11;
        p[0x65] = proto.p00P10M01P11;
        p[0x66] = proto.m00P10M01P11;
        p[0x68] = proto.z00M10M01P11;
        p[0x69] = proto.p00M10M01P11;
        p[0x6a] = proto.m00M10M01P11;
        p[0x80] = proto.z00Z10Z01M11;
        p[0x81] = proto.p00Z10Z01M11;
        p[0x82] = proto.m00Z10Z01M11;
        p[0x84] = proto.z00P10Z01M11;
        p[0x85] = proto.p00P10Z01M11;
        p[0x86] = proto.m00P10Z01M11;
        p[0x88] = proto.z00M10Z01M11;
        p[0x89] = proto.p00M10Z01M11;
        p[0x8a] = proto.m00M10Z01M11;
        p[0x90] = proto.z00Z10P01M11;
        p[0x91] = proto.p00Z10P01M11;
        p[0x92] = proto.m00Z10P01M11;
        p[0x94] = proto.z00P10P01M11;
        p[0x95] = proto.p00P10P01M11;
        p[0x96] = proto.m00P10P01M11;
        p[0x98] = proto.z00M10P01M11;
        p[0x99] = proto.p00M10P01M11;
        p[0x9a] = proto.m00M10P01M11;
        p[0xa0] = proto.z00Z10M01M11;
        p[0xa1] = proto.p00Z10M01M11;
        p[0xa2] = proto.m00Z10M01M11;
        p[0xa4] = proto.z00P10M01M11;
        p[0xa5] = proto.p00P10M01M11;
        p[0xa6] = proto.m00P10M01M11;
        p[0xa8] = proto.z00M10M01M11;
        p[0xa9] = proto.p00M10M01M11;
        p[0xaa] = proto.m00M10M01M11;
    }

    // ------------------------------------------------------------------
    // The level-curve processors. The comment blocks show the signs of the
    // bilinear function at the corners of the domain [0,1]^2 in the layout
    //   f01 f11
    //   f00 f10
    // ------------------------------------------------------------------

    protected z00Z10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x00
        // 0 0
        // 0 0
        //
        // This case occurs when each edge is shared by two coplanar faces, so
        // we have only two different normals. The normals are perpendicular.
        this.pairOf(c, mvc, 0, 0);
    }

    protected p00Z10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x01
        // 0 0
        // + 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected m00Z10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x02
        // 0 0
        // - 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected z00P10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x04
        // 0 0
        // 0 +

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00P10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x05
        // 0 0
        // + +

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected m00P10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x06
        // 0 0
        // - +

        // tmin = 0, tmax = 1
        // s = -f00 / (f10 - f00), (+)/(+)
        // 1-s = f10 / (f10 - f00), (+)/(+)
        // N = (1-s) * N0 + s * N1, omit denominator
        c.axis[0] = axpby3(c.f10, c.N[0], -c.f00, c.N[1]);
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected z00M10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x08
        // 0 0
        // 0 -

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00M10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x09
        // 0 0
        // + -

        // tmin = 0, tmax = 1
        // s = f00 / (f00 - f10), (+)/(+)
        // 1-s = -f10 / (f00 - f10), (+)/(+)
        // N = s * N1 + (1-s) * N0, omit denominator
        c.axis[0] = axpby3(-c.f10, c.N[0], c.f00, c.N[1]);
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected m00M10Z01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x0a
        // 0 0
        // - -

        // smin = 0, smax = 1, t = 1
        c.axis[1] = c.M[1].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected z00Z10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x10
        // + 0
        // 0 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x11
        // + 0
        // + 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);
    }

    protected m00Z10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x12
        // + 0
        // - 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1
        // t = -f00 / (f01 - f00), (+)/(+)
        // 1-t = f01 / (f01 - f00), (+)/(+)
        // M = (1-t) * M0 + t * M1, omit denominator
        c.axis[1] = axpby3(c.f01, c.M[0], -c.f00, c.M[1]);
        this.minimizerConstantT(c, mvc);
    }

    protected z00P10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x14
        // + 0
        // 0 +
        // It is not possible for a level curve to connect the corners.
        this.pairOf(c, mvc, 0, 0);
        this.pairOf(c, mvc, 1, 1);
    }

    protected p00P10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x15
        // + 0
        // + +
        this.pairOf(c, mvc, 1, 1);
    }

    protected m00P10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x16
        // + 0
        // - +

        // smin = 0
        // smax = -f00 / (f10 - f00), (+)/(+)
        this.minimizerVariableS(0, -c.f00, c.f10 - c.f00, c, mvc);

        this.pairOf(c, mvc, 1, 1);
    }

    protected z00M10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x18
        // + 0
        // 0 -

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected p00M10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x19
        // + 0
        // + -

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected m00M10P01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x1a
        // + 0
        // - -

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected z00Z10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x20
        // - 0
        // 0 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x21
        // - 0
        // + 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1
        // t = f00 / (f00 - f01), (+)/(+)
        // 1-t = -f01 / (f00 - f01), (+)/(+)
        // M = t * M1 + (1-t) * M0, omit denominator
        c.axis[1] = axpby3(-c.f01, c.M[0], c.f00, c.M[1]);
        this.minimizerConstantT(c, mvc);
    }

    protected m00Z10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x22
        // - 0
        // - 0

        // tmin = 0, tmax = 1, s = 1
        c.axis[0] = c.N[1].clone();
        this.minimizerConstantS(c, mvc);
    }

    protected z00P10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x24
        // - 0
        // 0 +

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected p00P10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x25
        // - 0
        // + +

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected m00P10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x26
        // - 0
        // - +

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected z00M10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x28
        // - 0
        // 0 -
        // It is not possible for a level curve to connect the corners.
        this.pairOf(c, mvc, 0, 0);
        this.pairOf(c, mvc, 1, 1);
    }

    protected p00M10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x29
        // - 0
        // + -

        // smin = 0
        // smax = f00 / (f00 - f10), (+)/(+)
        this.minimizerVariableS(0, c.f00, c.f00 - c.f10, c, mvc);

        this.pairOf(c, mvc, 1, 1);
    }

    protected m00M10M01Z11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x2a
        // - 0
        // - -
        this.pairOf(c, mvc, 1, 1);
    }

    protected z00Z10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x40
        // 0 +
        // 0 0

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x41
        // 0 +
        // + 0
        // It is not possible for a level curve to connect the corners.
        this.pairOf(c, mvc, 0, 1);
        this.pairOf(c, mvc, 1, 0);
    }

    protected m00Z10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x42
        // 0 +
        // - 0

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected z00P10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x44
        // 0 +
        // 0 +

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);
    }

    protected p00P10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x45
        // 0 +
        // + +
        this.pairOf(c, mvc, 0, 1);
    }

    protected m00P10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x46
        // 0 +
        // - +

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected z00M10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x48
        // 0 +
        // 0 -

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1
        // t = -f10 / (f11 - f10), (+)/(+)
        // 1-t = f11 / (f11 - f10), (+)/(+)
        // M = (1-t) * M0 + t * M1, omit denominator
        c.axis[1] = axpby3(c.f11, c.M[0], -c.f10, c.M[1]);
        this.minimizerConstantT(c, mvc);
    }

    protected p00M10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x49
        // 0 +
        // + -

        // smin = f00 / (f00 - f10), (+)/(+)
        // smax = 1
        const f00mf10 = c.f00 - c.f10;
        this.minimizerVariableS(c.f00, f00mf10, f00mf10, c, mvc);

        this.pairOf(c, mvc, 0, 1);
    }

    protected m00M10Z01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x4a
        // 0 +
        // - -

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected z00Z10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x50
        // + +
        // 0 0

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x51
        // + +
        // + 0
        this.pairOf(c, mvc, 1, 0);
    }

    protected m00Z10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x52
        // + +
        // - 0

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected z00P10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x54
        // + +
        // 0 +
        this.pairOf(c, mvc, 0, 0);
    }

    protected p00P10P01P11(_c: MinimumVolumeBox3FloatingPointCandidate,
        _mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x55
        // + +
        // + +
        // Nothing to do.
    }

    protected m00P10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x56
        // + +
        // - +

        // smin = 0
        // smax = -f00 / (f10 - f00), (+)/(+)
        this.minimizerVariableS(0, -c.f00, c.f10 - c.f00, c, mvc);
    }

    protected z00M10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x58
        // + +
        // 0 -

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected p00M10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x59
        // + +
        // + -

        // smin = f00 / (f00 - f10), (+)/(+)
        // smax = 1
        const f00mf10 = c.f00 - c.f10;
        this.minimizerVariableS(c.f00, f00mf10, f00mf10, c, mvc);
    }

    protected m00M10P01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x5a
        // + +
        // - -

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected z00Z10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x60
        // - +
        // 0 0

        // tmin = 0, tmax = 1
        // s = -f01 / (f11 - f01), (+)/(+)
        // 1-s = f11 / (f11 - f01), (+)/(+)
        // N = (1-s) * N0 + s * N1, omit denominator
        c.axis[0] = axpby3(c.f11, c.N[0], -c.f01, c.N[1]);
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x61
        // - +
        // + 0

        // smin = 0
        // smax = -f01 / (f11 - f01), (+)/(+)
        this.minimizerVariableS(0, -c.f01, c.f11 - c.f01, c, mvc);

        this.pairOf(c, mvc, 1, 0);
    }

    protected m00Z10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x62
        // - +
        // - 0

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected z00P10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x64
        // - +
        // 0 +

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected p00P10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x65
        // - +
        // + +

        // smin = 0
        // smax = -f01 / (f11 - f01), (+)/(+)
        this.minimizerVariableS(0, -c.f01, c.f11 - c.f01, c, mvc);
    }

    protected m00P10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x66
        // - +
        // - +

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected z00M10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x68
        // - +
        // 0 -

        // smin = -f01 / (f11 - f01), (+)/(+)
        // smax = 1
        const f11mf01 = c.f11 - c.f01;
        this.minimizerVariableS(-c.f01, f11mf01, f11mf01, c, mvc);

        this.pairOf(c, mvc, 0, 0);
    }

    protected p00M10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x69
        // - +
        // + -
        //
        // The level set F = 0 has two hyperbolic curves, each formed by a pair
        // of endpoints in {(0,t0), (s0,0), (s1,1), (1,t1)}, where
        // s0 = -f00 / (f10 - f00), s1 = -f01 / (f11 - f01),
        // t0 = -f00 / (f01 - f00), t1 = -f10 / (f11 - f10), all quantities in
        // (0,1). The two curves are on opposite sides of the asymptotes
        //   sa = (f01 - f00) / ((f01 - f00) + (f10 - f11))
        //   ta = (f10 - f00) / ((f10 - f00) + (f01 - f11))
        // If s0 < sa, one curve has endpoints {(0,t0),(s0,0)} and the other
        // curve has endpoints {(s1,1),(1,t1)}. If s0 > sa, one curve has
        // endpoints {(0,t0),(s1,1)} and the other curve has endpoints
        // {(s0,0),(1,t1)}. If s0 = sa, then segments of the asymptotes are the
        // two curves for the level set. Define d = f00 * f11 - f10 * f01. It
        // can be shown that
        //   s0 - sa = d / ((f10 - f00)((f10 - f00) + (f01 - f11))
        // The denominator is positive, so sign(s0 - sa) = sign(d). A similar
        // argument applies for the comparison between t0 and ta.
        const d = c.f00 * c.f11 - c.f10 * c.f01;
        if (d > 0) {
            // endpoints (s0,0) and (1,t1)
            // smin = f00 / (f00 - f10), (+)/(+)
            // smax = 1
            const f00mf10 = c.f00 - c.f10;
            this.minimizerVariableS(c.f00, f00mf10, f00mf10, c, mvc);

            // endpoints (0,t0) and (s1,1)
            // smin = 0
            // smax = -f01 / (f11 - f01), (+)/(+)
            const f11mf01 = c.f11 - c.f01;
            this.minimizerVariableS(0, -c.f01, f11mf01, c, mvc);
        } else if (d < 0) {
            // endpoints (0,t0) and (s0,0)
            // smin = 0
            // smax = f00 / (f00 - f10), (+)/(+)
            const f00mf10 = c.f00 - c.f10;
            this.minimizerVariableS(0, c.f00, f00mf10, c, mvc);

            // endpoints (s1,1) and (1,t1)
            // smin = -f01 / (f11 - f01), (+)/(+)
            // smax = 1
            const f11mf01 = c.f11 - c.f01;
            this.minimizerVariableS(-c.f01, f11mf01, f11mf01, c, mvc);
        } else {
            // endpoints (sa,0) and (sa,1)
            // sa = (f00 - f01) / ((f00 - f01) + (f11 - f10)), (+)/(+)
            // 1-sa = (f11 - f10) / ((f00 - f01) + (f11 - f10)), (+)/(+)
            // N = (1-sa) * N0 + sa * N1, omit the denominator
            c.axis[0] = axpby3(c.f11 - c.f10, c.N[0], c.f00 - c.f01, c.N[1]);
            this.minimizerConstantS(c, mvc);

            // endpoints (0,ta) and (1,ta)
            // ta = (f00 - f10) / ((f00 - f10) + (f11 - f01)), (+)/(+)
            // 1-ta = (f11 - f01) / ((f00 - f10) + (f11 - f01)), (+)/(+)
            // M = (1-ta) * M0 + ta * M1, omit the denominator
            c.axis[1] = axpby3(c.f11 - c.f01, c.M[0], c.f00 - c.f10, c.M[1]);
            this.minimizerConstantT(c, mvc);
        }
    }

    protected m00M10M01P11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x6a
        // - +
        // - -

        // smin = -f01 / (f11 - f01), (+)/(+)
        // smax = 1
        const f11mf01 = c.f11 - c.f01;
        this.minimizerVariableS(-c.f01, f11mf01, f11mf01, c, mvc);
    }

    protected z00Z10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x80
        // 0 -
        // 0 0

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x81
        // 0 -
        // + 0

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected m00Z10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x82
        // 0 -
        // - 0
        // It is not possible for a level curve to connect the corners.
        this.pairOf(c, mvc, 0, 1);
        this.pairOf(c, mvc, 1, 0);
    }

    protected z00P10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x84
        // 0 -
        // 0 +

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1
        // t = f10 / (f10 - f11), (+)/(+)
        // 1-t = -f11 / (f10 - f11), (+)/(+)
        // M = t * M1 + (1-t) * M0, omit the denominator
        c.axis[1] = axpby3(-c.f11, c.M[0], c.f10, c.M[1]);
        this.minimizerConstantT(c, mvc);
    }

    protected p00P10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x85
        // 0 -
        // + +

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected m00P10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x86
        // 0 -
        // - +

        // smin = -f00 / (f10 - f00), (+)/(+)
        // smax = 1
        const f10mf00 = c.f10 - c.f00;
        this.minimizerVariableS(-c.f00, f10mf00, f10mf00, c, mvc);
    }

    protected z00M10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x88
        // 0 -
        // 0 -

        // tmin = 0, tmax = 1, s = 0
        c.axis[0] = c.N[0].clone();
        this.minimizerConstantS(c, mvc);
    }

    protected p00M10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x89
        // 0 -
        // + -

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected m00M10Z01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x8a
        // 0 -
        // - -
        this.pairOf(c, mvc, 0, 1);
    }

    protected z00Z10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x90
        // + -
        // 0 0

        // tmin = 0, tmax = 1
        // s = f01 / (f01 - f11), (+)/(+)
        // 1-s = -f11 / (f01 - f11), (+)/(+)
        // N = s * N1 + (1-s) * N0, omit the denominator
        c.axis[0] = axpby3(-c.f11, c.N[0], c.f01, c.N[1]);
        this.minimizerConstantS(c, mvc);

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x91
        // + -
        // + 0

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected m00Z10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x92
        // + -
        // - 0

        // smin = 0
        // smax = f01 / (f01 - f11), (+)/(+)
        this.minimizerVariableS(0, c.f01, c.f01 - c.f11, c, mvc);

        this.pairOf(c, mvc, 1, 0);
    }

    protected z00P10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x94
        // + -
        // 0 +

        // smin = f01 / (f01 - f11), (+)/(+)
        // smax = 1
        const f01mf11 = c.f01 - c.f11;
        this.minimizerVariableS(c.f01, f01mf11, f01mf11, c, mvc);

        this.pairOf(c, mvc, 0, 0);
    }

    protected p00P10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x95
        // + -
        // + +

        // smin = f01 / (f01 - f11), (+)/(+)
        // smax = 1
        const f01mf11 = c.f01 - c.f11;
        this.minimizerVariableS(c.f01, f01mf11, f01mf11, c, mvc);
    }

    protected m00P10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x96
        // + -
        // - +
        //
        // The level set F = 0 has two hyperbolic curves, each formed by a pair
        // of endpoints in {(0,t0), (s0,0), (s1,1), (1,t1)}, where
        // s0 = -f00 / (f10 - f00), s1 = -f01 / (f11 - f01),
        // t0 = -f00 / (f01 - f00), t1 = -f10 / (f11 - f10), all quantities in
        // (0,1). The two curves are on opposite sides of the asymptotes
        //   sa = (f01 - f00) / ((f01 - f00) + (f10 - f11))
        //   ta = (f10 - f00) / ((f10 - f00) + (f01 - f11))
        // If s0 < sa, one curve has endpoints {(0,t0),(s0,0)} and the other
        // curve has endpoints {(s1,1),(1,t1)}. If s0 > sa, one curve has
        // endpoints {(0,t0),(s1,1)} and the other curve has endpoints
        // {(s0,0),(1,t1)}. If s0 = sa, then segments of the asymptotes are the
        // two curves for the level set. Define d = f00 * f11 - f10 * f01. It
        // can be shown that
        //   s0 - sa = d / ((f10 - f00)((f10 - f00) + (f01 - f11))
        // The denominator is positive, so sign(s0 - sa) = sign(d). A similar
        // argument applies for the comparison between t0 and ta.
        const d = c.f00 * c.f11 - c.f10 * c.f01;
        if (d > 0) {
            // endpoints (s0,0) and (1,t1)
            // smin = -f00 / (f10 - f00), (+)/(+)
            // smax = 1
            const f10mf00 = c.f10 - c.f00;
            this.minimizerVariableS(-c.f00, f10mf00, f10mf00, c, mvc);

            // endpoints (0,t0) and (s1,1)
            // smin = 0
            // smax = f01 / (f01 - f11)
            const f01mf11 = c.f01 - c.f11;
            this.minimizerVariableS(0, c.f01, f01mf11, c, mvc);
        } else if (d < 0) {
            // endpoints (0,t0) and (s0,0)
            // smin = 0
            // smax = -f00 / (f10 - f00), (+)/(+)
            const f10mf00 = c.f10 - c.f00;
            this.minimizerVariableS(0, -c.f00, f10mf00, c, mvc);

            // endpoints (s1,1) and (1,t1)
            // smin = f01 / (f01 - f11), (+)/(+)
            // smax = 1
            const f01mf11 = c.f01 - c.f11;
            this.minimizerVariableS(c.f01, f01mf11, f01mf11, c, mvc);
        } else {
            // endpoints (sa,0) and (sa,1)
            // sa = (f01 - f00) / ((f01 - f00) + (f10 - f11)), (+)/(+)
            // 1-sa = (f10 - f11) / ((f01 - f00) + (f10 - f11)), (+)/(+)
            // N = (1-sa) * N0 + sa * N1, omit the denominator
            c.axis[0] = axpby3(c.f10 - c.f11, c.N[0], c.f01 - c.f00, c.N[1]);
            this.minimizerConstantS(c, mvc);

            // endpoints (0,ta) and (1,ta)
            // ta = (f10 - f00) / ((f10 - f00) + (f01 - f11)), (+)/(+)
            // 1-ta = (f01 - f11) / ((f10 - f00) + (f01 - f11)), (+)/(+)
            // M = (1-ta) * M0 + ta * M1, omit the denominator
            c.axis[1] = axpby3(c.f01 - c.f11, c.M[0], c.f10 - c.f00, c.M[1]);
            this.minimizerConstantT(c, mvc);
        }
    }

    protected z00M10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x98
        // + -
        // 0 -

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected p00M10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x99
        // + -
        // + -

        // tmin = 0, tmax = 1
        this.minimizerVariableT(0, 1, 1, c, mvc);
    }

    protected m00M10P01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0x9a
        // + -
        // - -

        // smin = 0
        // smax = f01 / (f01 - f11), (+)/(+)
        this.minimizerVariableS(0, c.f01, c.f01 - c.f11, c, mvc);
    }

    protected z00Z10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa0
        // - -
        // 0 0

        // smin = 0, smax = 1, t = 0
        c.axis[1] = c.M[0].clone();
        this.minimizerConstantT(c, mvc);
    }

    protected p00Z10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa1
        // - -
        // + 0

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected m00Z10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa2
        // - -
        // - 0
        this.pairOf(c, mvc, 1, 0);
    }

    protected z00P10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa4
        // - -
        // 0 +

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected p00P10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa5
        // - -
        // + +

        // smin = 0, smax = 1
        this.minimizerVariableS(0, 1, 1, c, mvc);
    }

    protected m00P10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa6
        // - -
        // - +

        // smin = -f00 / (f10 - f00), (+)/(+)
        // smax = 1
        const f10mf00 = c.f10 - c.f00;
        this.minimizerVariableS(-c.f00, f10mf00, f10mf00, c, mvc);
    }

    protected z00M10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa8
        // - -
        // 0 -
        this.pairOf(c, mvc, 0, 0);
    }

    protected p00M10M01M11(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xa9
        // - -
        // + -

        // smin = 0
        // smax = f00 / (f00 - f10), (+)/(+)
        this.minimizerVariableS(0, c.f00, c.f00 - c.f10, c, mvc);
    }

    protected m00M10M01M11(_c: MinimumVolumeBox3FloatingPointCandidate,
        _mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        // index = 0xaa
        // - -
        // - -
        // Nothing to do.
    }
}
