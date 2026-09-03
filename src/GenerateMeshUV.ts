// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) GenerateMeshUV.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// This class is an implementation of the barycentric mapping algorithm
// described in Section 5.3 of the book
//     Polygon Mesh Processing
//     Mario Botsch, Leif Kobbelt, Mark Pauly, Pierre Alliez, Bruno Levy
//     AK Peters, Ltd., Natick MA, 2010
// It uses the mean value weights described in Section 5.3.1 to allow the mesh
// geometry to influence the texture coordinate generation, and it uses
// Gauss-Seidel iteration to solve the sparse linear system. The authors'
// advice is that the Gauss-Seidel approach works well for at most about 5000
// vertices, presumably the convergence rate degrading as the number of
// vertices increases.
//
// The algorithm implemented here has an additional preprocessing step that
// computes a topological distance transform of the vertices. The boundary
// texture coordinates are propagated inward by updating the vertices in
// topological distance order, leading to fast convergence for large numbers
// of vertices.
//
// Port notes:
//   * Upstream 'Vector3<Real>' / 'Vector2<Real>' are runtime-sized Vectors of
//     size 3 / size 2.
//   * The GPU-derived variants (GPUGenerateMeshUV.h, which override
//     SolveSystemInternal with HLSL/GLSL programs) are not ported. The
//     'padding' member of the vertex-graph record, which exists only to match
//     the HLSL/GLSL struct layout, is therefore omitted as well.
//   * The multithreaded CPU path (SolveSystemCPUMultiple) is not ported. It
//     is a pure Jacobi sweep over ping-pong buffers, so it is numerically
//     identical to the single-threaded path; the port keeps the 'numThreads'
//     constructor argument for API compatibility but always runs the
//     single-threaded sweep.
//   * Upstream stores the interior edges in a std::set<Edge*>, whose
//     iteration order is the (nondeterministic) pointer order. The weight
//     computation is independent per edge, so the port instead uses an array
//     in the mesh's edge-key order, which is deterministic.
//   * operator() becomes generate(); the caller supplies the tcoords array,
//     which is filled in place as upstream does.

import { GTE_C_TWO_PI } from './Constants.js';
import { ETManifoldMesh, ETManifoldMeshEdge } from './ETManifoldMesh.js';
import { Vector, length, sub } from './Vector.js';

// The vertex graph record required to set up the sparse linear system of
// equations that determines the texture coordinates. This is the port of the
// nested GenerateMeshUV<Real>::Vertex struct.
export interface GenerateMeshUVVertex {
    // The topological distance from the boundary of the mesh.
    distance: number;

    // The value range0 is the index into the vertex-graph data for the first
    // adjacent vertex. The value range1 is the number of adjacent vertices.
    range0: number;
    range1: number;
}

// One adjacency record: the adjacent vertex index and the mean value weight
// of that adjacency. This is the port of std::pair<int32_t, Real>.
export interface GenerateMeshUVAdjacency {
    first: number;
    second: number;
}

// The port of std::lower_bound on a sorted array of numbers.
function lowerBound(values: readonly number[], value: number): number {
    let lo = 0;
    let hi = values.length;
    while (lo < hi) {
        const mid = lo + ((hi - lo) >> 1);
        if (values[mid] < value) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

export class GenerateMeshUV {
    // Constructor inputs. Set numThreads to 0 when you want the code to run
    // in the main thread of the application. (The port always runs the
    // single-threaded sweep; see the port notes.) Provide a progress callback
    // when you want to monitor each iteration of the uv-solver. The input to
    // the callback is the current iteration; it starts at 1 and increases to
    // the numIterations input of generate().
    protected mNumThreads: number;
    protected mProgress: ((iteration: number) => void) | null;

    // Convenience members that store the input parameters to generate().
    protected mNumVertices: number;
    protected mVertices: readonly Vector[];
    protected mTCoords: Vector[];

    // The edge-triangle manifold graph, where each edge is shared by at most
    // two triangles.
    protected mGraph: ETManifoldMesh;

    // The mVertexInfo array stores -1 for the interior vertices. For a
    // boundary edge <v0,v1> that is counterclockwise, mVertexInfo[v0] = v1,
    // which gives us an ordered boundary polyline.
    protected mVertexInfo: number[];
    protected mNumBoundaryEdges: number;
    protected mBoundaryStart: number;
    protected mInteriorEdges: ETManifoldMeshEdge[];

    // The vertex graph required to set up a sparse linear system of equations
    // to determine the texture coordinates.
    protected mVertexGraph: GenerateMeshUVVertex[];
    protected mVertexGraphData: GenerateMeshUVAdjacency[];

    // The vertices are listed in the order determined by a topological
    // distance transform. Boundary vertices have 'distance' 0. Any vertices
    // that are not boundary vertices but are edge-adjacent to boundary
    // vertices have 'distance' 1. Neighbors of those have distance 2, and so
    // on. The mOrderedVertices array stores distance-0 vertices first,
    // distance-1 vertices second, and so on.
    protected mOrderedVertices: number[];

    constructor(numThreads = 0, progress: ((iteration: number) => void) | null = null) {
        this.mNumThreads = numThreads;
        this.mProgress = progress;
        this.mNumVertices = 0;
        this.mVertices = [];
        this.mTCoords = [];
        this.mGraph = new ETManifoldMesh();
        this.mVertexInfo = [];
        this.mNumBoundaryEdges = 0;
        this.mBoundaryStart = 0;
        this.mInteriorEdges = [];
        this.mVertexGraph = [];
        this.mVertexGraphData = [];
        this.mOrderedVertices = [];
    }

    // The incoming mesh must be edge-triangle manifold and have rectangle
    // topology (simply connected, closed polyline boundary). The arrays
    // 'vertices' and 'tcoords' must both have 'numVertices' elements. Set
    // 'useSquareTopology' to true for the generated coordinates to live in
    // the uv-square [0,1]^2. Set it to false for the generated coordinates to
    // live in a convex polygon that inscribes the uv-disk of center
    // (1/2,1/2) and radius 1/2.
    generate(numIterations: number, useSquareTopology: boolean, numVertices: number,
        vertices: readonly Vector[], numIndices: number, indices: readonly number[],
        tcoords: Vector[]): void {
        // Ensure that numIterations is even, which avoids having a memory
        // copy from the temporary ping-pong buffer to 'tcoords'.
        if (numIterations & 1) {
            ++numIterations;
        }

        this.mNumVertices = numVertices;
        this.mVertices = vertices;
        this.mTCoords = tcoords;

        // The linear system solver has a first pass to initialize the texture
        // coordinates to ensure the Gauss-Seidel iteration converges rapidly.
        // This requires the texture coordinates all start as (-1,-1).
        for (let i = 0; i < numVertices; ++i) {
            this.mTCoords[i].values[0] = -1;
            this.mTCoords[i].values[1] = -1;
        }

        // Create the manifold mesh data structure.
        this.mGraph.clear();
        const numTriangles = Math.trunc(numIndices / 3);
        let index = 0;
        for (let t = 0; t < numTriangles; ++t) {
            const v0 = indices[index++];
            const v1 = indices[index++];
            const v2 = indices[index++];
            this.mGraph.insert(v0, v1, v2);
        }

        this.topologicalVertexDistanceTransform();

        if (useSquareTopology) {
            this.assignBoundaryTextureCoordinatesSquare();
        } else {
            this.assignBoundaryTextureCoordinatesDisk();
        }

        this.computeMeanValueWeights();
        this.solveSystem(numIterations);
    }

    // A CPU-based implementation is provided by this class. Upstream's
    // GPU-derived classes override this function.
    protected solveSystemInternal(numIterations: number): void {
        this.solveSystemCPUSingle(numIterations);
    }

    private topologicalVertexDistanceTransform(): void {
        // Initialize the graph information.
        this.mVertexInfo = new Array<number>(this.mNumVertices).fill(-1);
        const edges = this.mGraph.getEdges();
        this.mVertexGraph = new Array<GenerateMeshUVVertex>(this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mVertexGraph[i] = { distance: 0, range0: 0, range1: 0 };
        }
        this.mVertexGraphData = new Array<GenerateMeshUVAdjacency>(2 * edges.length);
        for (let i = 0; i < this.mVertexGraphData.length; ++i) {
            this.mVertexGraphData[i] = { first: -1, second: -1 };
        }
        this.mOrderedVertices = new Array<number>(this.mNumVertices).fill(0);
        this.mInteriorEdges = [];
        this.mNumBoundaryEdges = 0;
        this.mBoundaryStart = Number.MAX_SAFE_INTEGER;

        // Count the number of adjacent vertices for each vertex. For data
        // sets with a large number of vertices, this is a preprocessing step
        // to avoid a dynamic data structure. Instead, a single array stores
        // all the adjacency information.
        const numAdjacencies = new Array<number>(this.mNumVertices).fill(0);

        for (const edge of edges) {
            ++numAdjacencies[edge.V[0]];
            ++numAdjacencies[edge.V[1]];

            if (edge.T[1]) {
                // This is an interior edge.
                this.mInteriorEdges.push(edge);
            } else {
                // This is a boundary edge. Determine the ordering of the
                // vertex indices to make the edge counterclockwise.
                ++this.mNumBoundaryEdges;
                const v0 = edge.V[0];
                const v1 = edge.V[1];
                const tri = edge.T[0];
                if (tri === null) {
                    continue;
                }
                for (let i = 0; i < 3; ++i) {
                    const v2 = tri.V[i];
                    if (v2 !== v0 && v2 !== v1) {
                        // The vertex is opposite the boundary edge.
                        const b0 = tri.V[(i + 1) % 3];
                        const b1 = tri.V[(i + 2) % 3];
                        this.mVertexInfo[b0] = b1;
                        this.mBoundaryStart = Math.min(this.mBoundaryStart, b0);
                        break;
                    }
                }
            }
        }

        // Set the range data for each vertex.
        for (let vIndex = 0, aIndex = 0; vIndex < this.mNumVertices; ++vIndex) {
            const numAdjacent = numAdjacencies[vIndex];
            this.mVertexGraph[vIndex].range0 = aIndex;
            this.mVertexGraph[vIndex].range1 = numAdjacent;
            aIndex += numAdjacent;
        }

        // Compute a topological distance transform of the vertices.
        const currFrontSet = new Set<number>();
        for (const edge of edges) {
            let v0 = edge.V[0];
            let v1 = edge.V[1];
            for (let i = 0; i < 2; ++i) {
                if (this.mVertexInfo[v0] === -1) {
                    this.mVertexGraph[v0].distance = -1;
                } else {
                    this.mVertexGraph[v0].distance = 0;
                    currFrontSet.add(v0);
                }

                // Insert v1 into the first available slot of the adjacency
                // array.
                const range0 = this.mVertexGraph[v0].range0;
                const range1 = this.mVertexGraph[v0].range1;
                for (let j = 0; j < range1; ++j) {
                    const data = this.mVertexGraphData[range0 + j];
                    if (data.second === -1) {
                        data.first = v1;
                        data.second = 0;
                        break;
                    }
                }

                const swap = v0;
                v0 = v1;
                v1 = swap;
            }
        }

        // Use a breadth-first search to propagate the distance information.
        // Upstream uses std::set<int32_t>, so the fronts are visited in
        // increasing vertex index; the port sorts to match.
        let nextDistance = 1;
        let currFront = Array.from(currFrontSet).sort((a, b) => a - b);
        let numFrontVertices = currFront.length;
        for (let i = 0; i < currFront.length; ++i) {
            this.mOrderedVertices[i] = currFront[i];
        }
        while (currFront.length > 0) {
            const nextFrontSet = new Set<number>();
            for (const v of currFront) {
                const range0 = this.mVertexGraph[v].range0;
                const range1 = this.mVertexGraph[v].range1;
                for (let j = 0; j < range1; ++j) {
                    const a = this.mVertexGraphData[range0 + j].first;
                    if (this.mVertexGraph[a].distance === -1) {
                        this.mVertexGraph[a].distance = nextDistance;
                        nextFrontSet.add(a);
                    }
                }
            }
            const nextFront = Array.from(nextFrontSet).sort((a, b) => a - b);
            for (let i = 0; i < nextFront.length; ++i) {
                this.mOrderedVertices[numFrontVertices + i] = nextFront[i];
            }
            numFrontVertices += nextFront.length;
            currFront = nextFront;
            ++nextDistance;
        }
    }

    private assignBoundaryTextureCoordinatesSquare(): void {
        // Map the boundary of the mesh to the unit square [0,1]^2. The
        // selection of square vertices is such that the relative distances
        // between boundary vertices and the relative distances between
        // polygon vertices is preserved, except that the four corners of the
        // square are required to have boundary points mapped to them. The
        // first boundary point has an implied distance of zero. The value
        // distance[i] is the length of the boundary polyline from vertex 0 to
        // vertex i+1.
        const distance = new Array<number>(this.mNumBoundaryEdges).fill(0);
        let total = 0;
        let v0 = this.mBoundaryStart;
        let v1: number;
        let i: number;
        for (i = 0; i < this.mNumBoundaryEdges; ++i) {
            v1 = this.mVertexInfo[v0];
            total += length(sub(this.mVertices[v1], this.mVertices[v0]));
            distance[i] = total;
            v0 = v1;
        }

        const invTotal = 1 / total;
        for (let k = 0; k < distance.length; ++k) {
            distance[k] *= invTotal;
        }

        // Upstream uses std::lower_bound with no tolerance. When the
        // accumulated arc length rounds so that a boundary vertex sitting
        // exactly on a quarter mark compares slightly less than 0.25 (0.50,
        // 0.75), the forced corner is assigned to the next boundary vertex
        // instead, displacing one vertex along the square boundary. The port
        // preserves this behavior; see the "Upstream bug suspects" notes.
        const endYMin = lowerBound(distance, 0.25);
        const endXMax = lowerBound(distance, 0.50);
        const endYMax = lowerBound(distance, 0.75);
        const endXMin = distance.length - 1;

        // The first polygon vertex is (0,0). The remaining vertices are
        // chosen counterclockwise around the square.
        v0 = this.mBoundaryStart;
        this.mTCoords[v0].values[0] = 0;
        this.mTCoords[v0].values[1] = 0;
        for (i = 0; i < endYMin; ++i) {
            v1 = this.mVertexInfo[v0];
            this.mTCoords[v1].values[0] = distance[i] * 4;
            this.mTCoords[v1].values[1] = 0;
            v0 = v1;
        }

        v1 = this.mVertexInfo[v0];
        this.mTCoords[v1].values[0] = 1;
        this.mTCoords[v1].values[1] = 0;
        v0 = v1;
        for (++i; i < endXMax; ++i) {
            v1 = this.mVertexInfo[v0];
            this.mTCoords[v1].values[0] = 1;
            this.mTCoords[v1].values[1] = distance[i] * 4 - 1;
            v0 = v1;
        }

        v1 = this.mVertexInfo[v0];
        this.mTCoords[v1].values[0] = 1;
        this.mTCoords[v1].values[1] = 1;
        v0 = v1;
        for (++i; i < endYMax; ++i) {
            v1 = this.mVertexInfo[v0];
            this.mTCoords[v1].values[0] = 3 - distance[i] * 4;
            this.mTCoords[v1].values[1] = 1;
            v0 = v1;
        }

        v1 = this.mVertexInfo[v0];
        this.mTCoords[v1].values[0] = 0;
        this.mTCoords[v1].values[1] = 1;
        v0 = v1;
        for (++i; i < endXMin; ++i) {
            v1 = this.mVertexInfo[v0];
            this.mTCoords[v1].values[0] = 0;
            this.mTCoords[v1].values[1] = 4 - distance[i] * 4;
            v0 = v1;
        }
    }

    private assignBoundaryTextureCoordinatesDisk(): void {
        // Map the boundary of the mesh to a convex polygon. The selection of
        // convex polygon vertices is such that the relative distances between
        // boundary vertices and the relative distances between polygon
        // vertices is preserved. The first boundary point has an implied
        // distance of zero. The value distance[i] is the length of the
        // boundary polyline from vertex 0 to vertex i+1.
        const distance = new Array<number>(this.mNumBoundaryEdges).fill(0);
        let total = 0;
        let v0 = this.mBoundaryStart;
        for (let i = 0; i < this.mNumBoundaryEdges; ++i) {
            const v1 = this.mVertexInfo[v0];
            total += length(sub(this.mVertices[v1], this.mVertices[v0]));
            distance[i] = total;
            v0 = v1;
        }

        // The convex polygon lives in [0,1]^2 and inscribes a circle with
        // center (1/2,1/2) and radius 1/2. The polygon center is not
        // necessarily the circle center! This is the case when a boundary
        // edge has length larger than half the total length of the boundary
        // polyline; we do not expect such data for our meshes. The first
        // polygon vertex is (1,1/2), the angle-zero point of the circle. The
        // remaining vertices are chosen counterclockwise around the polygon.
        const multiplier = GTE_C_TWO_PI / total;
        v0 = this.mBoundaryStart;
        this.mTCoords[v0].values[0] = 1;
        this.mTCoords[v0].values[1] = 0.5;
        for (let i = 1, im1 = 0; i < this.mNumBoundaryEdges; ++i, ++im1) {
            const v1 = this.mVertexInfo[v0];
            const angle = multiplier * distance[im1];
            this.mTCoords[v1].values[0] = (Math.cos(angle) + 1) * 0.5;
            this.mTCoords[v1].values[1] = (Math.sin(angle) + 1) * 0.5;
            v0 = v1;
        }
    }

    private computeMeanValueWeights(): void {
        for (const edge of this.mInteriorEdges) {
            let v0 = edge.V[0];
            let v1 = edge.V[1];
            for (let i = 0; i < 2; ++i) {
                // Compute the direction from X0 to X1 and compute the length
                // of the edge (X0,X1).
                const X0 = this.mVertices[v0];
                const X1 = this.mVertices[v1];
                const X1mX0 = sub(X1, X0);
                const x1mx0Length = length(X1mX0);
                let weight: number;
                if (x1mx0Length > 0) {
                    const dirX1mX0 = [
                        X1mX0.values[0] / x1mx0Length,
                        X1mX0.values[1] / x1mx0Length,
                        X1mX0.values[2] / x1mx0Length
                    ];

                    // Compute the weight for X0 associated with X1.
                    weight = 0;
                    for (let j = 0; j < 2; ++j) {
                        // Find the vertex of triangle T[j] opposite edge
                        // <X0,X1>.
                        const tri = edge.T[j];
                        if (tri === null) {
                            continue;
                        }
                        for (let k = 0; k < 3; ++k) {
                            const v2 = tri.V[k];
                            if (v2 !== v0 && v2 !== v1) {
                                const X2 = this.mVertices[v2];
                                const X2mX0 = sub(X2, X0);
                                const x2mx0Length = length(X2mX0);
                                if (x2mx0Length > 0) {
                                    const dot =
                                        (X2mX0.values[0] / x2mx0Length) * dirX1mX0[0] +
                                        (X2mX0.values[1] / x2mx0Length) * dirX1mX0[1] +
                                        (X2mX0.values[2] / x2mx0Length) * dirX1mX0[2];
                                    const cs = Math.min(Math.max(dot, -1), 1);
                                    const angle = Math.acos(cs);
                                    weight += Math.tan(angle * 0.5);
                                } else {
                                    weight += 1;
                                }
                                break;
                            }
                        }
                    }
                    weight /= x1mx0Length;
                } else {
                    weight = 1;
                }

                const range0 = this.mVertexGraph[v0].range0;
                const range1 = this.mVertexGraph[v0].range1;
                for (let j = 0; j < range1; ++j) {
                    const data = this.mVertexGraphData[range0 + j];
                    if (data.first === v1) {
                        data.second = weight;
                    }
                }

                const swap = v0;
                v0 = v1;
                v1 = swap;
            }
        }
    }

    private solveSystem(numIterations: number): void {
        // On the first pass, average only neighbors whose texture coordinates
        // have been computed. This is a good initial guess for the linear
        // system and leads to relatively fast convergence of the Gauss-Seidel
        // iterates.
        for (let i = this.mNumBoundaryEdges; i < this.mNumVertices; ++i) {
            const v0 = this.mOrderedVertices[i];
            const range0 = this.mVertexGraph[v0].range0;
            const range1 = this.mVertexGraph[v0].range1;
            let tcoord0 = 0;
            let tcoord1 = 0;
            let weightSum = 0;
            for (let j = 0; j < range1; ++j) {
                const current = this.mVertexGraphData[range0 + j];
                const v1 = current.first;
                if (this.mTCoords[v1].values[0] !== -1) {
                    const weight = current.second;
                    weightSum += weight;
                    tcoord0 += weight * this.mTCoords[v1].values[0];
                    tcoord1 += weight * this.mTCoords[v1].values[1];
                }
            }
            this.mTCoords[v0].values[0] = tcoord0 / weightSum;
            this.mTCoords[v0].values[1] = tcoord1 / weightSum;
        }

        this.solveSystemInternal(numIterations);
    }

    private solveSystemCPUSingle(numIterations: number): void {
        // Use ping-pong buffers for the texture coordinates.
        let inTCoords = new Array<number>(2 * this.mNumVertices);
        let outTCoords = new Array<number>(2 * this.mNumVertices);
        for (let i = 0; i < this.mNumVertices; ++i) {
            inTCoords[2 * i] = this.mTCoords[i].values[0];
            inTCoords[2 * i + 1] = this.mTCoords[i].values[1];
            outTCoords[2 * i] = inTCoords[2 * i];
            outTCoords[2 * i + 1] = inTCoords[2 * i + 1];
        }

        // The value numIterations is even, so we always swap an even number
        // of times. This ensures that on exit from the loop, inTCoords holds
        // the final texture coordinates.
        for (let i = 1; i <= numIterations; ++i) {
            if (this.mProgress) {
                this.mProgress(i);
            }

            for (let j = this.mNumBoundaryEdges; j < this.mNumVertices; ++j) {
                const v0 = this.mOrderedVertices[j];
                const range0 = this.mVertexGraph[v0].range0;
                const range1 = this.mVertexGraph[v0].range1;
                let tcoord0 = 0;
                let tcoord1 = 0;
                let weightSum = 0;
                for (let k = 0; k < range1; ++k) {
                    const current = this.mVertexGraphData[range0 + k];
                    const v1 = current.first;
                    const weight = current.second;
                    weightSum += weight;
                    tcoord0 += weight * inTCoords[2 * v1];
                    tcoord1 += weight * inTCoords[2 * v1 + 1];
                }
                outTCoords[2 * v0] = tcoord0 / weightSum;
                outTCoords[2 * v0 + 1] = tcoord1 / weightSum;
            }

            const swap = inTCoords;
            inTCoords = outTCoords;
            outTCoords = swap;
        }

        for (let i = 0; i < this.mNumVertices; ++i) {
            this.mTCoords[i].values[0] = inTCoords[2 * i];
            this.mTCoords[i].values[1] = inTCoords[2 * i + 1];
        }
    }
}
