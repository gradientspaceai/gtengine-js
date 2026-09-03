// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) MeshStaticManifold2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// MeshStaticManifold2 represents a vertex-edge-triangle manifold mesh for
// which triangles are provided as a single batch and no mesh modification
// operations are going to be performed on the mesh. It is a requirement that
// the input triangles form a manifold mesh with consistently ordered
// triangles. In most applications, this requirement is already satisfied.
//
// Port notes:
//   - Upstream stores the vertex adjacency 4-tuples in one contiguous block
//     of 12 * numTriangles size_t values and hands each vertex a pointer into
//     that block. The port gives each vertex its own array of 4-tuples; the
//     insertion order, and therefore every observable result, is the same.
//     The single-block layout is a C++ memory-management optimization with no
//     TypeScript counterpart.
//   - Upstream optionally distributes the adjacency update over threads. The
//     port always runs the single-threaded loop, which computes the same
//     values (the multithreaded version is a pure partition of the vertex
//     loop). The numThreads constructor parameter is retained and ignored.
//   - 'invalid' is Number.MAX_SAFE_INTEGER rather than SIZE_MAX, so the
//     upstream 'index < numVertices' guards still reject it.
//   - The depth-first search uses an explicit empty-stack marker instead of
//     relying on unsigned wraparound of the stack top.
//   - The boundary-edge std::map iterates its keys in lexicographic order,
//     which determines the order of the polygons and their starting vertices.
//     The port sorts the keys explicitly to reproduce that order.
//   - getAdjacentTriangles returns { adj0, adj1, exists } instead of using
//     output reference parameters. It also fixes an upstream bug; see the
//     comment in that method.

import { logAssert } from './Logger.js';

// The vertices are stored as an array of MeshStaticManifold2Vertex objects,
// one per vertex index. If triangle[t0] = <v0,v1,v2>, then vertex[v0]
// contains a 4-tuple {v1,v2,t0,a0}. The undirected edge (v1,v2) is opposite
// v0. If there is no adjacent triangle sharing (v1,v2), then a0 is invalid.
// If there is an adjacent triangle, then a0 is the index for that triangle.
// Let triangle[a0] = <v2,v1,v3>; then vertex[v3] contains a 4-tuple
// {v2,v1,a0,t0}.
export class MeshStaticManifold2Vertex {
    // If triangle t0 is <v0,v1,v2> in counterclockwise order, then the
    // corresponding adjacents element is {v1,v2,t0,a0}, where a0 is invalid
    // when (v1,v2) is contained by a single triangle or a0 is the index for
    // the adjacent triangle when (v1,v2) is contained by two triangles.
    //
    // Upstream declares MeshStaticManifold2 a friend so that only it may
    // write the members. TypeScript has no equivalent; treat 'adjacents' and
    // insert() as internal to MeshStaticManifold2.
    readonly adjacents: number[][] = [];

    getNumAdjacents(): number {
        return this.adjacents.length;
    }

    getAdjacents(): number[][] {
        return this.adjacents;
    }

    insert(v1: number, v2: number, t: number, location: number): void {
        this.adjacents.push([v1, v2, t, location]);
    }
}

// A boundary edge of the mesh, the triangle t that contains it and the
// location a of the opposite vertex in that triangle.
class BoundaryEdge {
    constructor(
        public t: number,
        public a: number,
        public visited: boolean) {
    }
}

export class MeshStaticManifold2 {
    // Use the maximum safe integer to denote an invalid index.
    static readonly invalid = Number.MAX_SAFE_INTEGER;

    private mVertices: MeshStaticManifold2Vertex[];
    private mTriangles: number[][];
    private mAdjacents: number[][];
    private mMinTrianglesAtVertex: number;
    private mMaxTrianglesAtVertex: number;

    // Preconditions.
    //   1. The triangles input must have size 1 or larger.
    //   2. The number of vertices must be 3 or larger.
    //   3. The triangles must form a manifold mesh.
    //   4. Each triangle must be nondegenerate; no repeated vertices.
    //   5. The triangles must all be ordered counterclockwise.
    // The numThreads parameter of the upstream constructor is accepted for
    // signature compatibility but ignored; see the port notes.
    constructor(numVertices: number, triangles: number[][], numThreads: number = 0) {
        void numThreads;
        logAssert(numVertices >= 3 && triangles.length > 0, 'Invalid input.');

        this.mVertices = new Array<MeshStaticManifold2Vertex>(numVertices);
        for (let v = 0; v < numVertices; ++v) {
            this.mVertices[v] = new MeshStaticManifold2Vertex();
        }
        this.mTriangles = triangles.map(tri => [tri[0], tri[1], tri[2]]);
        this.mAdjacents = triangles.map(() => [
            MeshStaticManifold2.invalid,
            MeshStaticManifold2.invalid,
            MeshStaticManifold2.invalid
        ]);
        this.mMinTrianglesAtVertex = 0;
        this.mMaxTrianglesAtVertex = 0;

        this.getNumTrianglesAtVertex(numVertices);
        this.populate();
        this.updateAdjacencyForSharedEdges();
    }

    // Member access.
    getVertices(): MeshStaticManifold2Vertex[] {
        return this.mVertices;
    }

    // Each 3-tuple contains indices into the vertices.
    getTriangles(): number[][] {
        return this.mTriangles;
    }

    // Each 3-tuple contains indices into the triangles.
    getAdjacents(): number[][] {
        return this.mAdjacents;
    }

    getMinNumTrianglesAtVertex(): number {
        return this.mMinTrianglesAtVertex;
    }

    getMaxNumTrianglesAtVertex(): number {
        return this.mMaxTrianglesAtVertex;
    }

    // Determine whether or not the undirected edge (v0,v1) exists.
    edgeExists(v0: number, v1: number): boolean {
        if (v0 < this.mVertices.length && v1 < this.mVertices.length && v0 !== v1) {
            return this.getDirectedEdge(v0, v1) !== null
                || this.getDirectedEdge(v1, v0) !== null;
        }
        return false;
    }

    // Get the adjacent triangles for the undirected edge (v0,v1). The
    // returned adjacent triangle indices adj0 and adj1 are the following:
    //
    //   1. <v0,v1> and <v1,v0> are both directed edges, so the edge is shared
    //      by two triangles and both adj0 and adj1 are valid (not equal to
    //      invalid). The index adj0 is the L-triangle for <v0,v1> and the
    //      index adj1 is the R-triangle for <v0,v1>. Equivalently, adj0 is
    //      the R-triangle for <v1,v0> and adj1 is the L-triangle for <v1,v0>.
    //
    //   2. <v0,v1> is directed but <v1,v0> does not exist. The index adj0 is
    //      the L-triangle for <v0,v1> and the index adj1 is invalid (no
    //      R-triangle).
    //
    //   3. <v1,v0> is directed but <v0,v1> does not exist. The index adj0 is
    //      invalid (no L-triangle) for <v1,v0> and the index adj1 is the
    //      R-triangle for <v1,v0>.
    //
    //   4. Neither <v0,v1> nor <v1,v0> exist; that is, the edge does not
    //      occur for any triangle.
    //
    // It is possible to distinguish among the 4 cases by examining the
    // returned indices:
    //   (1) returns (valid, valid) and exists 'true'
    //   (2) returns (valid, invalid) and exists 'true'
    //   (3) returns (invalid, valid) and exists 'true'
    //   (4) returns (invalid, invalid) and exists 'false'
    getAdjacentTriangles(v0: number, v1: number): {
        adj0: number;
        adj1: number;
        exists: boolean;
    } {
        if (v0 < this.mVertices.length && v1 < this.mVertices.length && v0 !== v1) {
            // UPSTREAM BUG (MeshStaticManifold2.h, GetAdjacentTriangles).
            // Upstream looks up only one of the two directed edges and
            // returns components [2] and [3] of the single 4-tuple it finds.
            // The 4-tuple stored at vertex[v0] for triangle <v0,v1,v2> is
            // {v1, v2, t, a}, where 'a' is the triangle adjacent across the
            // edge OPPOSITE v0, that is (v1,v2) -- not across the queried
            // edge (v0,v1). So upstream returns the adjacency of a different
            // edge as adj1. The error is visible in the upstream
            // documentation itself: case 3 is documented to return
            // (invalid, valid), but the upstream code returns
            // (*adjacents1)[2] as adj0, which is a valid triangle index, so
            // case 3 can never be produced.
            //
            // The port computes both directed lookups, as the documented
            // contract requires: adj0 is the triangle containing the directed
            // edge <v0,v1> (the L-triangle) and adj1 is the triangle
            // containing <v1,v0> (the R-triangle).
            const adjacents0 = this.getDirectedEdge(v0, v1);
            const adjacents1 = this.getDirectedEdge(v1, v0);
            if (adjacents0 !== null || adjacents1 !== null) {
                return {
                    adj0: adjacents0 !== null ? adjacents0[2] : MeshStaticManifold2.invalid,
                    adj1: adjacents1 !== null ? adjacents1[2] : MeshStaticManifold2.invalid,
                    exists: true
                };
            }
        }

        return {
            adj0: MeshStaticManifold2.invalid,
            adj1: MeshStaticManifold2.invalid,
            exists: false
        };
    }

    // The connected components of the mesh. Each component is an array of
    // triangle indices.
    getComponents(): number[][] {
        const components: number[][] = [];

        // The values are 0 (unvisited), 1 (discovered), 2 (finished).
        const visited = new Array<number>(this.mTriangles.length).fill(0);

        // Share a stack for the depth-first search. This avoids allocating
        // and deallocating a stack for each call to depthFirstSearch.
        const sharedStack = new Array<number>(this.mTriangles.length).fill(0);

        // The code reserves maximum space for the component in order to avoid
        // allocation costs associated with resizing caused by push_back.
        const sharedComponents = new Array<number>(this.mTriangles.length).fill(0);

        for (let t = 0; t < this.mTriangles.length; ++t) {
            if (visited[t] === 0) {
                const numInserted = this.depthFirstSearch(t, visited, sharedStack,
                    sharedComponents);
                components.push(sharedComponents.slice(0, numInserted));
            }
        }

        return components;
    }

    // Compute the boundary-edge components of the mesh. These are polygons
    // that are simple for the strict definition of manifold mesh that
    // disallows bow-tie configurations. The GTE mesh implementations do allow
    // bow-tie configurations, in which case some polygons might not be
    // simple. If you select duplicateEndpoints to be false, a component has
    // consecutive vertices (v[0], v[1], ..., v[n-1]) and the polygon has
    // edges (v[0],v[1]), (v[1],v[2]), ..., (v[n-2],v[n-1]), (v[n-1],v[0]).
    // If duplicateEndpoints is set to true, a component has consecutive
    // vertices (v[0], v[1], ..., v[n-1], v[0]), emphasizing that the
    // component is closed.
    getBoundaryPolygons(duplicateEndpoints: boolean): number[][] {
        const polygons: number[][] = [];

        // Get the boundary edges.
        const boundaryEdges = new Map<string, BoundaryEdge>();
        for (let t = 0; t < this.mTriangles.length; ++t) {
            const tri = this.mTriangles[t];
            for (let a = 0; a < 3; ++a) {
                if (this.mAdjacents[t][a] === MeshStaticManifold2.invalid) {
                    const directed: [number, number] = [tri[(a + 1) % 3], tri[(a + 2) % 3]];
                    const key = MeshStaticManifold2.edgeKey(directed[0], directed[1]);
                    // std::map::insert does not overwrite an existing key.
                    if (!boundaryEdges.has(key)) {
                        boundaryEdges.set(key, new BoundaryEdge(t, a, false));
                    }
                }
            }
        }

        // Extract the polygons. Each polygon is the boundary for a connected
        // component of the mesh. The keys are visited in the lexicographic
        // order of the std::map that upstream uses.
        const sortedKeys = [...boundaryEdges.keys()].sort((key0, key1) => {
            const e0 = MeshStaticManifold2.parseEdgeKey(key0);
            const e1 = MeshStaticManifold2.parseEdgeKey(key1);
            return e0[0] !== e1[0] ? e0[0] - e1[0] : e0[1] - e1[1];
        });

        for (const key of sortedKeys) {
            const initialEdge = boundaryEdges.get(key) as BoundaryEdge;
            if (!initialEdge.visited) {
                polygons.push(this.getBoundaryPolygon(initialEdge.t, initialEdge.a,
                    boundaryEdges));
            }
        }

        if (!duplicateEndpoints) {
            for (const polygon of polygons) {
                polygon.length = polygon.length - 1;
            }
        }

        return polygons;
    }

    private static edgeKey(v0: number, v1: number): string {
        return `${v0},${v1}`;
    }

    private static parseEdgeKey(key: string): [number, number] {
        const comma = key.indexOf(',');
        return [Number(key.substring(0, comma)), Number(key.substring(comma + 1))];
    }

    // Count the number of triangles sharing each vertex. The total number of
    // indices for triangles adjacent to vertices is 3 * numTriangles. The
    // minimum and maximum triangle counts are for statistical information.
    private getNumTrianglesAtVertex(numVertices: number): void {
        const counts = new Array<number>(numVertices).fill(0);
        for (const tri of this.mTriangles) {
            for (let i = 0; i < 3; ++i) {
                ++counts[tri[i]];
            }
        }

        this.mMinTrianglesAtVertex = counts[0];
        this.mMaxTrianglesAtVertex = counts[0];
        for (const count of counts) {
            if (count < this.mMinTrianglesAtVertex) {
                this.mMinTrianglesAtVertex = count;
            }
            if (count > this.mMaxTrianglesAtVertex) {
                this.mMaxTrianglesAtVertex = count;
            }
        }
    }

    // Populate the adjacency information for the vertices.
    private populate(): void {
        for (let t = 0; t < this.mTriangles.length; ++t) {
            const tri = this.mTriangles[t];
            const v0 = tri[0];
            const v1 = tri[1];
            const v2 = tri[2];

            // The last arguments (i = 0, 1 or 2) are used to set the correct
            // mAdjacents[][i] indices. These arguments are replaced later by
            // the actual indices for adjacent triangles sharing the edge.
            this.mVertices[v0].insert(v1, v2, t, 0);
            this.mVertices[v1].insert(v2, v0, t, 1);
            this.mVertices[v2].insert(v0, v1, t, 2);
        }
    }

    // Update triangle adjacency information for edges that are shared by two
    // triangles.
    private updateAdjacencyForSharedEdges(): void {
        for (let v = 0; v < this.mVertices.length; ++v) {
            this.updateAdjacencyForEdge(v);
        }
    }

    private updateAdjacencyForEdge(v0: number): void {
        const vertex0 = this.mVertices[v0];
        for (const adjacents0 of vertex0.adjacents) {
            const v1 = adjacents0[0];
            const v2 = adjacents0[1];

            // The edge opposite vertex v0 is (v1,v2). We know that vertex[v0]
            // contains a 4-tuple {v1,v2,tri0,loc0}. Determine whether
            // vertex[v2] contains a 4-tuple {v1,v3,adj1,loc1}.
            const adjacents1 = this.getDirectedEdge(v2, v1);
            if (adjacents1 !== null) {
                // The edge <v1,v2> has a triangle adjacent to triangle tri0.
                // Update the vertex adjacency information for triangle tri0
                // at that edge. Triangle a1 adjacency is not updated. It will
                // be updated when <v2,v1> is visited at another time. This
                // avoids two writes of the adjacent triangle indices.
                const tri0 = adjacents0[2];
                const loc0 = adjacents0[3];
                const adj1 = adjacents1[2];
                adjacents0[3] = adj1;
                this.mAdjacents[tri0][loc0] = adj1;
            } else {
                // Replace the mAdjacents[] location value (0, 1 or 2) by an
                // invalid index because edge <v1,v0> does not exist, in which
                // case there is no adjacent triangle to edge <v0,v1>.
                adjacents0[3] = MeshStaticManifold2.invalid;
            }
        }
    }

    private getDirectedEdge(v0: number, v1: number): number[] | null {
        const vertex0 = this.mVertices[v0];
        for (const adjacents0 of vertex0.adjacents) {
            if (adjacents0[0] === v1) {
                return adjacents0;
            }
        }
        return null;
    }

    private depthFirstSearch(tInitial: number, visited: number[], tStack: number[],
        component: number[]): number {
        // Upstream relies on unsigned wraparound of 'top' to denote an empty
        // stack; the port uses -1 for that.
        let top = -1;
        let numInserted = 0;

        tStack[++top] = tInitial;
        while (top !== -1) {
            const t = tStack[top];
            visited[t] = 1;
            let i = 0;
            for (i = 0; i < 3; ++i) {
                const tAdjacent = this.mAdjacents[t][i];
                if (tAdjacent !== MeshStaticManifold2.invalid && visited[tAdjacent] === 0) {
                    tStack[++top] = tAdjacent;
                    break;
                }
            }
            if (i === 3) {
                visited[t] = 2;
                component[numInserted] = t;
                ++numInserted;
                --top;
            }
        }

        return numInserted;
    }

    private getBoundaryPolygon(tStart: number, aStart: number,
        boundaryEdges: Map<string, BoundaryEdge>): number[] {
        const polygon: number[] = [];
        let t = tStart;
        let a = aStart;
        let tri = this.mTriangles[t];
        let i0 = (a + 1) % 3;
        let i1 = (i0 + 1) % 3;
        const vEdge: [number, number] = [tri[i0], tri[i1]];
        polygon.push(vEdge[0]);
        for (; ;) {
            // The upstream std::map::operator[] inserts a default-constructed
            // BoundaryEdge when the key is absent, which happens only for
            // inputs that violate the manifold-mesh precondition.
            let edge = boundaryEdges.get(MeshStaticManifold2.edgeKey(vEdge[0], vEdge[1]));
            if (edge === undefined) {
                edge = new BoundaryEdge(MeshStaticManifold2.invalid,
                    MeshStaticManifold2.invalid, false);
                boundaryEdges.set(MeshStaticManifold2.edgeKey(vEdge[0], vEdge[1]), edge);
            }
            if (edge.visited) {
                break;
            }

            polygon.push(vEdge[1]);
            edge.visited = true;

            // Traverse the triangle strip with vertex at vEdge[1] until the
            // last triangle is encountered. The final edge of the last
            // triangle is the next boundary edge and starts at vEdge[1].
            a = this.mAdjacents[t][i0];
            while (a !== MeshStaticManifold2.invalid) {
                // Get the next triangle in the strip.
                t = a;
                tri = this.mTriangles[t];
                for (i1 = 0; i1 < 3; ++i1) {
                    if (vEdge[1] === tri[i1]) {
                        // Get the next interior edge in the triangle strip,
                        // namely, <tri[i0], tri[i1]>.
                        i0 = (i1 + 2) % 3;
                        a = this.mAdjacents[t][i0];
                        break;
                    }
                }
                logAssert(i1 < 3, 'Unexpected condition.');
            }

            const i2 = (i1 + 1) % 3;
            vEdge[0] = vEdge[1];
            vEdge[1] = tri[i2];
            i0 = i1;
            i1 = (i0 + 1) % 3;
        }

        return polygon;
    }
}
