// gtengine-js: TypeScript port of Geometric Tools Engine (GTE)
// ConstrainedDelaunay2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Compute the Delaunay triangulation of the input points and then insert
// edges that are constrained to be in the triangulation. For each such edge,
// a retriangulation of the triangle strip containing the edge is required.
// NOTE: If two constrained edges overlap at a point that is an interior point
// of each edge, the second insertion will interfere with the retriangulation
// of the first edge. Although the code here will do what is requested, a pair
// of such edges usually indicates the upstream process that generated the
// edges is not doing what it should.
//
// Port notes:
// * Upstream ConstrainedDelaunay2.h declares a variadic class template with
//   two specializations: the deprecated
//   ConstrainedDelaunay2<InputType, ComputeType> and the replacement
//   ConstrainedDelaunay2<InputType>. Only the replacement is ported, matching
//   the Delaunay2 port.
// * Upstream selects ComputeRational = BSNumber<UIntegerFP32<70 or 526>> with
//   a preallocated pool of scratch values (mCRPool) and a Copy() helper that
//   widens an input rational into a compute rational. The port's BSNumber is
//   bigint-backed and grows as needed, so the fixed word counts, the two
//   rational types, the pool and Copy() are unnecessary and are dropped. The
//   exact-arithmetic results are identical.
// * The upstream private ToLine(...) is a verbatim duplicate of the protected
//   Delaunay2::ToLine(...) expression tree, so the port calls the inherited
//   toLine(...) rather than duplicating it.
// * The std::unordered_set<EdgeKey<false>> of inserted edges becomes a Map
//   keyed by the unordered edge; getInsertedEdges() returns the keys in
//   sorted order so the output does not depend on hash ordering.
// * The upstream 'tristrip' scratch ETManifoldMesh is only used as a set of
//   triangle keys, so the port uses a Map keyed by the unordered triangle key.
// * The upstream 'void Insert(edge, std::vector<int32_t>& partitionedEdge)'
//   becomes 'insert(edge): number[]' returning the partitioned edge.

import { logAssert } from './Logger.js';
import { BSNumber } from './BSNumber.js';
import { Delaunay2 } from './Delaunay2.js';
import { EdgeKey } from './EdgeKey.js';
import { TriangleKey } from './TriangleKey.js';
import { FeatureKey } from './FeatureKey.js';
import type { ETManifoldMeshTriangle } from './ETManifoldMesh.js';
import type { VETManifoldMeshVertex } from './VETManifoldMesh.js';
import type { Vector } from './Vector.js';

export class ConstrainedDelaunay2 extends Delaunay2 {
    // All edges inserted via the insert(...) call are stored for use by the
    // caller. If any edge passed to insert(...) is partitioned into subedges,
    // the subedges are inserted into this member.
    private mInsertedEdges: Map<string, EdgeKey>;

    constructor() {
        super();
        this.mInsertedEdges = new Map<string, EdgeKey>();
    }

    // This function computes the Delaunay triangulation only. Edges are
    // inserted later.
    //
    // Upstream bug (fixed in the port): the upstream operator() does not
    // clear mInsertedEdges, so recomputing the triangulation of a second data
    // set with the same object leaves the inserted edges of the first data
    // set in the set returned by GetInsertedEdges(). The class documents
    // itself as "a functor to support computing the constrained Delaunay
    // triangulation of multiple data sets using the same class object", so
    // the stale edges corrupt the reported result. The port clears the set.
    override compute(vertices: readonly Vector[]): boolean {
        this.mInsertedEdges.clear();
        return super.compute(vertices);
    }

    // The 'edge' is the constrained edge to be inserted into the
    // triangulation. If that edge is already in the triangulation, the
    // function returns without any retriangulation and the returned array
    // contains the input 'edge'. If 'edge' is coincident with 1 or more edges
    // already in the triangulation, 'edge' is partitioned into subedges which
    // are then inserted. It is also possible that 'edge' does not overlap
    // already existing edges in the triangulation but has interior points
    // that are vertices in the triangulation; in this case, 'edge' is
    // partitioned and the subedges are inserted. In either case, the returned
    // array is an ordered list of indices into the triangulation vertices
    // that are on the edge. It is guaranteed that the first element is
    // edge[0] and the last element is edge[1] (after replacement of duplicate
    // vertices by their representatives).
    insert(inEdge: readonly [number, number]): number[] {
        logAssert(this.getDimension() === 2,
            'The Delaunay triangulation must have dimension 2.');

        // The upstream code indexes the duplicates array before verifying
        // that the incoming indices are in range, which is undefined behavior
        // for out-of-range input. The port validates the input indices first.
        const numVertices = this.getNumVertices();
        logAssert(
            0 <= inEdge[0] && inEdge[0] < numVertices &&
            0 <= inEdge[1] && inEdge[1] < numVertices,
            'Invalid edge.');

        // Replace duplicates by their vertex indices that represent the
        // equivalence classes. The representatives are guaranteed to be in
        // the mGraph vertex map.
        const duplicates = this.getDuplicates();
        const edge: [number, number] = [duplicates[inEdge[0]], duplicates[inEdge[1]]];

        logAssert(edge[0] !== edge[1], 'Invalid edge.');

        // The 'partition' array stores the endpoints of the incoming edge if
        // that edge does not contain interior points that are vertices of the
        // Delaunay triangulation. If the edge contains one or more vertices
        // in its interior, the edge is partitioned into subedges, each
        // subedge having vertex endpoints but no interior point that is a
        // vertex. The partition is stored in the 'partition' array.
        const partition: [number, number][] = [];

        // When using exact arithmetic, a while(!edgeConsumed) loop suffices.
        // Just in case the code has a bug, guard against an infinite loop.
        let edgeConsumed = false;
        const numTriangles = this.mGraph.getNumTriangles();
        for (let t = 0; t < numTriangles && !edgeConsumed; ++t) {
            if (this.mGraph.getEdge(edge[0], edge[1]) !== null) {
                // The edge already exists in the triangulation.
                this.addInsertedEdge(edge[0], edge[1]);
                partition.push([edge[0], edge[1]]);
                break;
            }

            // Get the link edges for the vertex edge[0]. These edges are
            // opposite the link vertex.
            const linkEdges = this.getLinkEdges(edge[0]);

            // Determine which link triangle contains the to-be-inserted edge.
            for (const linkEdge of linkEdges) {
                // Compute on which side of the to-be-inserted edge the link
                // vertices live. The triangles are not degenerate, so it is
                // not possible for sign0 = sign1 = 0.
                const v0 = linkEdge[0];
                const v1 = linkEdge[1];
                const sign0 = this.toLine(v0, edge[0], edge[1]);
                const sign1 = this.toLine(v1, edge[0], edge[1]);
                if (sign0 >= 0 && sign1 <= 0) {
                    if (sign0 > 0) {
                        if (sign1 < 0) {
                            // The triangle <edge[0], v0, v1> strictly
                            // contains the to-be-inserted edge. Gather the
                            // triangles in the triangle strip containing the
                            // edge.
                            edgeConsumed = this.processTriangleStrip(edge, v0, v1, partition);
                        }
                        else {  // sign1 == 0 && sign0 > 0
                            // The to-be-inserted edge is coincident with the
                            // triangle edge <edge[0], v1>, and it is
                            // guaranteed that the vertex at v1 is an interior
                            // point of <edge[0],edge[1]> because we previously
                            // tested whether edge[] is in the triangulation.
                            edgeConsumed = this.processCoincidentEdge(edge, v1, partition);
                        }
                    }
                    else {  // sign0 == 0 && sign1 < 0
                        // The to-be-inserted edge is coincident with the
                        // triangle edge <edge[0], v0>, and it is guaranteed
                        // that the vertex at v0 is an interior point of
                        // <edge[0],edge[1]> because we previously tested
                        // whether edge[] is in the triangulation.
                        edgeConsumed = this.processCoincidentEdge(edge, v0, partition);
                    }
                    break;
                }
            }
        }

        // If the following assertion is triggered, the linkEdges-loop exited
        // without ever calling processTriangleStrip or processCoincidentEdge.
        logAssert(partition.length > 0, ConstrainedDelaunay2.cdtMessage);

        const partitionedEdge = new Array<number>(partition.length + 1);
        for (let i = 0; i < partition.length; ++i) {
            partitionedEdge[i] = partition[i][0];
        }
        partitionedEdge[partition.length] = partition[partition.length - 1][1];
        return partitionedEdge;
    }

    // All edges inserted via the insert(...) call are stored for use by the
    // caller. If any edge passed to insert(...) is partitioned into subedges,
    // the subedges are stored but not the original edge. The returned array
    // is sorted by edge key.
    getInsertedEdges(): EdgeKey[] {
        const keys = Array.from(this.mInsertedEdges.values());
        keys.sort((e0, e1) => FeatureKey.compare(e0, e1));
        return keys;
    }

    // The interface functions of the base class Delaunay2 are valid, so
    // access to any Delaunay information is allowed. Perhaps the most
    // important member function is getGraph(), which returns the
    // ETManifoldMesh that represents the constrained Delaunay triangulation.
    // NOTE: If you want access to the compact arrays via getIndices(),
    // getAdjacencies(), getTriangleIndices(t) or getTriangleAdjacencies(t),
    // you must first call updateIndicesAdjacencies() to ensure that the
    // compact arrays are synchronized with the Delaunay graph.

    private addInsertedEdge(v0: number, v1: number): void {
        const ekey = new EdgeKey(false, v0, v1);
        this.mInsertedEdges.set(ekey.mapKey(), ekey);
    }

    // For a vertex at index v, return the edges of the adjacent triangles,
    // each triangle having v as a vertex and the returned edge is opposite v.
    private getLinkEdges(v: number): [number, number][] {
        const found = this.mGraph.getVertex(v);
        logAssert(found !== null, 'Failed to find vertex in graph.');
        const vertex = found as VETManifoldMeshVertex;

        const linkEdges: [number, number][] = [];
        for (const linkTri of vertex.getTAdjacent()) {
            let j: number;
            for (j = 0; j < 3; ++j) {
                if (linkTri.V[j] === vertex.V) {
                    linkEdges.push([linkTri.V[(j + 1) % 3], linkTri.V[(j + 2) % 3]]);
                    break;
                }
            }
            logAssert(j < 3, 'Unexpected condition.');
        }
        return linkEdges;
    }

    // The return value is 'true' if the edge did not have to be subdivided
    // because it has an interior point that is a vertex. The return value is
    // 'false' if it does have such a point, in which case edge[0] is updated
    // to the index of that vertex. The caller must process the new edge.
    private processTriangleStrip(edge: [number, number], inV0: number, inV1: number,
        partition: [number, number][]): boolean {
        let v0 = inV0;
        let v1 = inV1;
        let edgeConsumed = true;
        const localEdge: [number, number] = [edge[0], edge[1]];

        // Locate and store the triangles in the triangle strip containing the
        // edge.
        const tristrip = new Map<string, [number, number, number]>();
        const addToTristrip = (a: number, b: number, c: number): void => {
            tristrip.set(new TriangleKey(true, a, b, c).mapKey(), [a, b, c]);
        };
        addToTristrip(localEdge[0], v0, v1);

        const foundTri = this.mGraph.getTriangle(localEdge[0], v0, v1);
        logAssert(foundTri !== null, ConstrainedDelaunay2.cdtMessage);
        let tri = foundTri as ETManifoldMeshTriangle;

        // Keep track of the right and left polylines that bound the triangle
        // strip. These polylines can have coincident edges. In particular,
        // this happens when the current triangle in the strip shares an edge
        // with a previous triangle in the strip and the previous triangle is
        // not the immediate predecessor to the current triangle.
        const rightPolygon: number[] = [];
        const leftPolygon: number[] = [];
        rightPolygon.push(localEdge[0]);
        rightPolygon.push(v0);
        leftPolygon.push(localEdge[0]);
        leftPolygon.push(v1);

        // When using exact arithmetic, an unbounded loop suffices. Guard
        // against an infinite loop.
        const numTriangles = this.mGraph.getNumTriangles();
        let t: number;
        for (t = 0; t < numTriangles; ++t) {
            // The current triangle is tri and has edge <v0,v1>. Get the
            // triangle adj that is adjacent to tri via this edge.
            const foundAdj = tri.getAdjacentOfEdge(v0, v1);
            logAssert(foundAdj !== null, ConstrainedDelaunay2.cdtMessage);
            const adj = foundAdj as ETManifoldMeshTriangle;
            addToTristrip(adj.V[0], adj.V[1], adj.V[2]);

            // Get the vertex of adj that is opposite edge <v0,v1>.
            const opposite = adj.getOppositeVertexOfEdge(v0, v1);
            logAssert(opposite.found, ConstrainedDelaunay2.cdtMessage);
            const vOpposite = opposite.uOpposite;
            if (vOpposite === localEdge[1]) {
                // The triangle strip containing the edge is complete.
                break;
            }

            // The next triangle in the strip depends on whether the opposite
            // vertex is left-of the edge, right-of the edge or on the edge.
            const querySign = this.toLine(vOpposite, localEdge[0], localEdge[1]);
            if (querySign > 0) {
                tri = adj;
                v0 = vOpposite;
                rightPolygon.push(v0);
            }
            else if (querySign < 0) {
                tri = adj;
                v1 = vOpposite;
                leftPolygon.push(v1);
            }
            else {
                // The to-be-inserted edge contains an interior point that is
                // also a vertex in the triangulation. The edge must be
                // subdivided. The first subedge is in a triangle strip that
                // is processed by code below that is outside the loop. The
                // second subedge must be processed by the caller.
                localEdge[1] = vOpposite;
                edge[0] = vOpposite;
                edgeConsumed = false;
                break;
            }
        }
        logAssert(t < numTriangles, ConstrainedDelaunay2.cdtMessage);

        // Insert the final endpoint of the to-be-inserted edge.
        rightPolygon.push(localEdge[1]);
        leftPolygon.push(localEdge[1]);

        // The retriangulation depends on counterclockwise ordering of the
        // boundary right and left polygons. The right polygon is already
        // counterclockwise ordered. The left polygon is clockwise ordered, so
        // reverse it.
        leftPolygon.reverse();

        // Update the inserted edges.
        this.addInsertedEdge(localEdge[0], localEdge[1]);
        partition.push([localEdge[0], localEdge[1]]);

        // Remove the triangle strip from the full triangulation. This must
        // occur before the retriangulation which inserts new triangles into
        // the full triangulation.
        for (const tv of tristrip.values()) {
            this.mGraph.remove(tv[0], tv[1], tv[2]);
        }

        // Retriangulate the tristrip region.
        this.retriangulate(leftPolygon);
        this.retriangulate(rightPolygon);

        return edgeConsumed;
    }

    // Process a to-be-inserted edge that is coincident with an already
    // existing triangulation edge.
    private processCoincidentEdge(edge: [number, number], v: number,
        partition: [number, number][]): boolean {
        this.addInsertedEdge(edge[0], v);
        partition.push([edge[0], v]);
        edge[0] = v;
        return v === edge[1];
    }

    // Retriangulate the polygon via a bisection-like method that finds
    // vertices closest to the current polygon base edge. The function is
    // naturally recursive, but simulated recursion is used to avoid a large
    // program stack by instead using the heap.
    private retriangulate(polygon: readonly number[]): void {
        const stack: [number, number][] = [];
        stack.push([0, polygon.length - 1]);
        while (stack.length > 0) {
            const i = stack.pop() as [number, number];
            if (i[1] > i[0] + 1) {
                // Get the vertex indices for the specified i-values.
                const v0 = polygon[i[0]];
                const v1 = polygon[i[1]];

                // Select isplit in the index range [i[0]+1,i[1]-1] so that
                // the vertex at index polygon[isplit] attains the minimum
                // distance to the edge with vertices at the indices
                // polygon[i[0]] and polygon[i[1]].
                const isplit = this.selectSplit(polygon, i[0], i[1]);
                const vsplit = polygon[isplit];

                // Insert the triangle into the Delaunay graph.
                this.mGraph.insert(v0, vsplit, v1);

                stack.push([i[0], isplit]);
                stack.push([isplit, i[1]]);
            }
        }
    }

    // Determine the polygon vertex with index strictly between i0 and i1 that
    // minimizes the pseudosquared distance from that vertex to the line
    // segment whose endpoints are at indices i0 and i1.
    private selectSplit(polygon: readonly number[], i0: number, i1: number): number {
        let i2: number;
        if (i1 === i0 + 2) {
            // This is the only candidate.
            i2 = i0 + 1;
        }
        else {  // i1 - i0 > 2
            // Select the index i2 in [i0+1,i1-1] for which the distance from
            // the vertex v2 at i2 to the edge <v0,v1> is minimized. To allow
            // exact arithmetic, use a pseudosquared distance that avoids
            // divisions and square roots.
            i2 = i0 + 1;
            const irV0 = this.mIRVertices[polygon[i0]];
            const irV1 = this.mIRVertices[polygon[i1]];

            // Precompute some common values that are used in all calls to
            // computePSD.
            const crV1mV0x = irV1[0].sub(irV0[0]);
            const crV1mV0y = irV1[1].sub(irV0[1]);
            const crSqrLen10 =
                crV1mV0x.mul(crV1mV0x).add(crV1mV0y.mul(crV1mV0y));

            // Locate the minimum pseudosquared distance.
            let crMinPSD = ConstrainedDelaunay2.computePSD(irV0, irV1,
                this.mIRVertices[polygon[i2]], crV1mV0x, crV1mV0y, crSqrLen10);
            for (let i = i2 + 1; i < i1; ++i) {
                const crPSD = ConstrainedDelaunay2.computePSD(irV0, irV1,
                    this.mIRVertices[polygon[i]], crV1mV0x, crV1mV0y, crSqrLen10);
                if (crPSD.lessThan(crMinPSD)) {
                    crMinPSD = crPSD;
                    i2 = i;
                }
            }
        }
        return i2;
    }

    // Compute a pseudosquared distance from the vertex at v2 to the edge
    // <v0,v1>. The result is exact for rational arithmetic and does not
    // involve division.
    private static computePSD(irV0: readonly BSNumber[], irV1: readonly BSNumber[],
        irV2: readonly BSNumber[], crV1mV0x: BSNumber, crV1mV0y: BSNumber,
        crSqrLen10: BSNumber): BSNumber {
        const crV2mV0x = irV2[0].sub(irV0[0]);
        const crV2mV0y = irV2[1].sub(irV0[1]);
        const crDot1020 =
            crV1mV0x.mul(crV2mV0x).add(crV1mV0y.mul(crV2mV0y));

        if (crDot1020.getSign() <= 0) {
            const crSqrLen20 =
                crV2mV0x.mul(crV2mV0x).add(crV2mV0y.mul(crV2mV0y));
            return crSqrLen10.mul(crSqrLen20);
        }

        const crV2mV1x = irV2[0].sub(irV1[0]);
        const crV2mV1y = irV2[1].sub(irV1[1]);
        const crDot1021 =
            crV1mV0x.mul(crV2mV1x).add(crV1mV0y.mul(crV2mV1y));
        if (crDot1021.getSign() >= 0) {
            const crSqrLen21 =
                crV2mV1x.mul(crV2mV1x).add(crV2mV1y.mul(crV2mV1y));
            return crSqrLen10.mul(crSqrLen21);
        }

        const crSqrLen20 =
            crV2mV0x.mul(crV2mV0x).add(crV2mV0y.mul(crV2mV0y));
        return crSqrLen10.mul(crSqrLen20).sub(crDot1020.mul(crDot1020));
    }

    private static readonly cdtMessage =
        'The failed assertion is unexpected when using arbitrary-precision ' +
        'arithmetic. Please file a bug report and provide the input dataset.';
}
