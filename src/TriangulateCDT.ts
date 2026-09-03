// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TriangulateCDT.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The fundamental problem is to compute the triangulation of a polygon tree.
// The outer polygons have counterclockwise ordered vertices. The inner
// polygons have clockwise ordered vertices. The algorithm uses Constrained
// Delaunay Triangulation and the implementation allows polygons to share
// vertices and edges.
//
// The polygons are not required to be simple in the sense that a vertex can
// be shared by an even number of edges, where the number is larger than 2.
// The input points can have duplicates, which the triangulator handles
// correctly. The algorithm supports coincident vertex-edge and coincident
// edge-edge configurations. See the document
//   https://www.geometrictools.com/Documentation/TriangulationByCDT.pdf
// for examples.
//
// If two edges intersect at edge-interior points, the current implementation
// cannot handle this. A pair of such edges cannot simultaneously be inserted
// into the constrained triangulation without affecting each other's local
// re-triangulation.
//
// The input points are a vertex pool. The input tree is a PolygonTree object,
// defined in PolygonTree.ts. Any outer polygon has vertices
// points[outer[0]] through points[outer[outer.length-1]] listed in
// counterclockwise order. Any inner polygon has vertices points[inner[0]]
// through points[inner[inner.length-1]] listed in clockwise order. The output
// tree contains the triangulation of the polygon tree on a per-node basis. If
// coincident vertex-edge or coincident edge-edge configurations exist in the
// polygon tree, the corresponding output polygons differ from the input
// polygons in that they have more vertices due to edge splits. The triangle
// chirality (winding order) is the same as the containing polygon.
//
// Port notes:
// * Upstream TriangulateCDT.h declares a variadic class template with two
//   specializations, the deprecated TriangulateCDT<InputType, ComputeType>
//   and the replacement TriangulateCDT<InputType>. The two bodies are
//   identical apart from the deprecated one's int32_t point counts and its
//   explicit ComputeType for ConstrainedDelaunay2. Only the replacement is
//   ported, matching the ConstrainedDelaunay2 and Delaunay2 ports. The
//   upstream discussion of the worst-case UIntegerFP32<N> word counts does
//   not apply: the port's BSNumber is bigint-backed and grows as needed.
// * The two upstream operator() overloads differ only in how the point array
//   is passed. The port has one compute(...). The 'outputTree' output
//   reference parameter becomes the return value, so each call produces a
//   fresh PolygonTreeEx; upstream clears every member of the output tree at
//   the start anyway, so no state can be carried over from a previous call.
// * std::queue becomes an array plus a read index; std::set<EdgeKey<false>>
//   becomes a Set of EdgeKey.mapKey() strings; std::set<TriangleKey<true>>
//   becomes a Map from TriangleKey.mapKey() to the key, iterated in sorted
//   key order so the extracted triangulations do not depend on insertion
//   order (upstream relies on the std::set ordering).
// * std::map<Vector2<T>, int32_t> in RemapPolygonTree is a Map keyed by the
//   string form of the point coordinates. C++ compares the coordinates with
//   operator< on the doubles, for which -0 and +0 compare equal; the string
//   key uses Number.prototype.toString, which also maps -0 and +0 to "0", so
//   the equivalence classes agree.
//
// Upstream quirk (preserved): in RemapPolygonTree, when a point is a
// duplicate of an earlier point, upstream executes
// 'remapping[iter->second] = node.polygon[i]', which overwrites the
// remapping of the first occurrence with the input index of the later
// duplicate. RestorePolygonTree therefore reports the input index of the
// last duplicate encountered rather than the first. Because the two input
// points have identical coordinates, the geometry of the result is
// unaffected, so the port preserves the behavior.
//
// Upstream quirk (preserved): ConstrainedTriangulate copies the
// ConstrainedDelaunay2 graph into 'graph' and then re-inserts every triangle
// of cdt.GetIndices() into that copy. The copy already contains exactly
// those triangles, so each insertion finds the existing triangle and does
// nothing. The loop is kept because it also builds tree.allTriangles.

import { logAssert } from './Logger';
import { ConstrainedDelaunay2 } from './ConstrainedDelaunay2';
import { EdgeKey } from './EdgeKey';
import { FeatureKey } from './FeatureKey';
import { ETManifoldMesh } from './ETManifoldMesh';
import { PolygonTree, PolygonTreeEx, PolygonTreeExNode } from './PolygonTree';
import { TriangleKey } from './TriangleKey';
import type { Vector } from './Vector';

// The port of std::set<TriangleKey<true>>: membership by the primitive map
// key, with iteration in the sorted order that std::set uses.
class TriangleKeySet {
    private mMap: Map<string, TriangleKey>;

    constructor() {
        this.mMap = new Map<string, TriangleKey>();
    }

    get size(): number {
        return this.mMap.size;
    }

    has(key: TriangleKey): boolean {
        return this.mMap.has(key.mapKey());
    }

    insert(key: TriangleKey): void {
        this.mMap.set(key.mapKey(), key);
    }

    sorted(): TriangleKey[] {
        const keys = Array.from(this.mMap.values());
        keys.sort(FeatureKey.compare);
        return keys;
    }
}

export class TriangulateCDT {
    // The input points are a vertex pool of 2D points; inputPoints.length
    // must be at least 3. The 'inputTree' polygons are indices into
    // 'inputPoints'. The returned tree is the triangulation of the polygon
    // tree, the port of upstream's 'outputTree' parameter.
    compute(inputPoints: readonly Vector[], inputTree: PolygonTree):
        PolygonTreeEx {
        logAssert(inputPoints.length >= 3, 'Invalid argument.');
        for (const point of inputPoints) {
            logAssert(point.size === 2, 'TriangulateCDT requires 2D points.');
        }

        const outputTree = new PolygonTreeEx();
        copyAndCompactify(inputTree, outputTree);
        triangulate(inputPoints, outputTree);
        return outputTree;
    }
}

function copyAndCompactify(input: PolygonTree, output: PolygonTreeEx): void {
    // Count the number of nodes in the tree.
    let numNodes = 1;  // the root node
    let queue: PolygonTree[] = [input];
    for (let head = 0; head < queue.length; ++head) {
        const node = queue[head];
        numNodes += node.child.length;
        for (const child of node.child) {
            queue.push(child);
        }
    }

    // Create the PolygonTreeEx nodes.
    output.nodes = new Array<PolygonTreeExNode>(numNodes);
    for (let i = 0; i < numNodes; ++i) {
        output.nodes[i] = new PolygonTreeExNode();
        output.nodes[i].self = i;
    }
    output.nodes[0].chirality = +1;
    output.nodes[0].parent = PolygonTreeEx.INVALID;

    let current = 0, last = 0, minChild = 1;
    queue = [input];
    for (let head = 0; head < queue.length; ++head) {
        const node = queue[head];
        const exnode = output.nodes[current++];
        exnode.polygon = node.polygon.slice();
        exnode.minChild = minChild;
        exnode.supChild = minChild + node.child.length;
        minChild = exnode.supChild;
        for (const child of node.child) {
            const exchild = output.nodes[++last];
            exchild.chirality = -exnode.chirality;
            exchild.parent = exnode.self;
            queue.push(child);
        }
    }
}

function triangulate(inputPoints: readonly Vector[],
    tree: PolygonTreeEx): void {
    // The constrained Delaunay triangulator will be given the unique points
    // referenced by the polygons in the tree. The tree 'polygon' indices are
    // relative to inputPoints[], but they are temporarily mapped to indices
    // relative to 'points'. Once the triangulation is complete, the indices
    // are restored to those relative to inputPoints[].
    const { points, remapping } = remapPolygonTree(inputPoints, tree);
    logAssert(points.length >= 3, 'Invalid polygon tree.');

    const edges = new Set<string>();
    const graph = constrainedTriangulate(tree, points, edges);
    classifyTriangles(tree, graph, edges);

    restorePolygonTree(tree, remapping);
}

// On return, 'points' are the unique inputPoints[] values referenced by the
// tree. The tree 'polygon' members are modified to be indices into 'points'
// rather than inputPoints[]. The 'remapping' allows us to restore the tree
// 'polygon' members to be indices into inputPoints[] after the triangulation
// is computed.
function remapPolygonTree(inputPoints: readonly Vector[], tree: PolygonTreeEx):
    { points: Vector[]; remapping: number[] } {
    const pointMap = new Map<string, number>();
    const points: Vector[] = [];
    let currentIndex = 0;

    // The remapping is initially the identity, remapping[j] = j.
    const remapping = new Array<number>(inputPoints.length);
    for (let j = 0; j < remapping.length; ++j) {
        remapping[j] = j;
    }

    const queue: number[] = [0];
    for (let head = 0; head < queue.length; ++head) {
        const node = tree.nodes[queue[head]];
        const numIndices = node.polygon.length;
        for (let i = 0; i < numIndices; ++i) {
            const point = inputPoints[node.polygon[i]];
            const key = `${point.values[0]},${point.values[1]}`;
            const found = pointMap.get(key);
            if (found === undefined) {
                // The point is encountered the first time.
                pointMap.set(key, currentIndex);
                remapping[currentIndex] = node.polygon[i];
                node.polygon[i] = currentIndex;
                points.push(point);
                ++currentIndex;
            }
            else {
                // The point is a duplicate. On the remapping, the polygon[]
                // value is set to the index for the first occurrence of the
                // duplicate.
                remapping[found] = node.polygon[i];
                node.polygon[i] = found;
            }
        }

        for (let c = node.minChild; c < node.supChild; ++c) {
            queue.push(c);
        }
    }

    return { points, remapping };
}

function restorePolygonTree(tree: PolygonTreeEx,
    remapping: readonly number[]): void {
    const restore = (triangles: [number, number, number][]): void => {
        for (const tri of triangles) {
            for (let j = 0; j < 3; ++j) {
                tri[j] = remapping[tri[j]];
            }
        }
    };

    const queue: number[] = [0];
    for (let head = 0; head < queue.length; ++head) {
        const node = tree.nodes[queue[head]];

        for (let i = 0; i < node.polygon.length; ++i) {
            node.polygon[i] = remapping[node.polygon[i]];
        }
        restore(node.triangulation);

        for (let c = node.minChild; c < node.supChild; ++c) {
            queue.push(c);
        }
    }

    restore(tree.interiorTriangles);
    restore(tree.exteriorTriangles);
    restore(tree.insideTriangles);
    restore(tree.outsideTriangles);
    restore(tree.allTriangles);
}

// The 'edges' set is filled with all the polygon edges that must be in the
// triangulation. The returned graph is a duplicate of the constrained
// Delaunay graph; it is modified during triangle extraction.
function constrainedTriangulate(tree: PolygonTreeEx,
    points: readonly Vector[], edges: Set<string>): ETManifoldMesh {
    // Use constrained Delaunay triangulation.
    const cdt = new ConstrainedDelaunay2();
    cdt.compute(points);

    const queue: number[] = [0];
    for (let head = 0; head < queue.length; ++head) {
        const node = tree.nodes[queue[head]];

        const replacement: number[] = [];
        let numIndices = node.polygon.length;
        for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
            // Insert the polygon edge into the constrained Delaunay
            // triangulation.
            const outEdge = cdt.insert([node.polygon[i0], node.polygon[i1]]);
            if (outEdge.length > 2) {
                // The polygon edge intersects additional vertices in the
                // triangulation. The outEdge[] edge values are
                // { edge[0], other_vertices, edge[1] } which are ordered
                // along the line segment.
                for (let k = 1; k < outEdge.length; ++k) {
                    replacement.push(outEdge[k]);
                }
            }
            else {
                replacement.push(node.polygon[i1]);
            }
        }
        if (replacement.length > node.polygon.length) {
            node.polygon = replacement;
        }

        numIndices = node.polygon.length;
        for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
            edges.add(new EdgeKey(false, node.polygon[i0],
                node.polygon[i1]).mapKey());
        }

        for (let c = node.minChild; c < node.supChild; ++c) {
            queue.push(c);
        }
    }

    // Copy the graph to the compact arrays mIndices and mAdjacencies for use
    // by the caller.
    cdt.updateIndicesAdjacencies();

    // Duplicate the graph, which will be modified during triangle
    // extraction. The original is kept intact for use by the caller.
    const graph = cdt.getGraph().clone();

    // Store the triangles in allTriangles for potential use by the caller.
    const numTriangles = cdt.getNumTriangles();
    const indices = cdt.getIndices();
    tree.allTriangles = new Array<[number, number, number]>(numTriangles);
    for (let t = 0; t < numTriangles; ++t) {
        const v0 = indices[3 * t];
        const v1 = indices[3 * t + 1];
        const v2 = indices[3 * t + 2];
        graph.insert(v0, v1, v2);
        tree.allTriangles[t] = [v0, v1, v2];
    }

    return graph;
}

function classifyTriangles(tree: PolygonTreeEx, graph: ETManifoldMesh,
    edges: Set<string>): void {
    classifyDFS(tree, 0, graph, edges);
    logAssert(edges.size === 0,
        'The edges should be empty for a correct implementation.');
    getOutsideTriangles(tree, graph);
    getInsideTriangles(tree);
}

function classifyDFS(tree: PolygonTreeEx, index: number, graph: ETManifoldMesh,
    edges: Set<string>): void {
    const node = tree.nodes[index];
    for (let c = node.minChild; c < node.supChild; ++c) {
        classifyDFS(tree, c, graph, edges);
    }

    const region = new TriangleKeySet();
    const numIndices = node.polygon.length;
    for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
        const v0 = node.polygon[i0];
        const v1 = node.polygon[i1];
        const edge = graph.getEdge(v0, v1);
        logAssert(edge !== null, 'Unexpected condition.');
        const tri0 = edge.T[0];
        logAssert(tri0 !== null, 'Unexpected condition.');
        if (tri0.whichSideOfEdge(v0, v1) === node.chirality) {
            region.insert(new TriangleKey(true, tri0.V[0], tri0.V[1], tri0.V[2]));
        }
        else {
            const tri1 = edge.T[1];
            if (tri1 !== null) {
                region.insert(new TriangleKey(true, tri1.V[0], tri1.V[1], tri1.V[2]));
            }
        }
    }

    fillRegion(graph, edges, region);
    extractTriangles(graph, region, node);
    for (let i0 = numIndices - 1, i1 = 0; i1 < numIndices; i0 = i1++) {
        edges.delete(new EdgeKey(false, node.polygon[i0],
            node.polygon[i1]).mapKey());
    }
}

// On input, the set has the initial seeds for the desired region. A
// breadth-first search is performed to find the connected component of the
// seeds. The component is bounded by an outer polygon and the inner polygons
// of its children.
function fillRegion(graph: ETManifoldMesh, edges: ReadonlySet<string>,
    region: TriangleKeySet): void {
    const regionQueue: TriangleKey[] = region.sorted();

    for (let head = 0; head < regionQueue.length; ++head) {
        const tkey = regionQueue[head];
        const tri = graph.getTriangle(tkey.V[0], tkey.V[1], tkey.V[2]);
        logAssert(tri !== null, 'Unexpected condition.');
        for (let j = 0; j < 3; ++j) {
            const edge = tri.E[j];
            if (edge !== null) {
                const ekey = new EdgeKey(false, edge.V[0], edge.V[1]).mapKey();
                if (!edges.has(ekey)) {
                    // The edge is not constrained, so it allows the search to
                    // continue.
                    const adj = tri.T[j];
                    if (adj !== null) {
                        const akey = new TriangleKey(true, adj.V[0], adj.V[1],
                            adj.V[2]);
                        if (!region.has(akey)) {
                            // The adjacent triangle has not yet been visited,
                            // so place it in the queue to continue the search.
                            region.insert(akey);
                            regionQueue.push(akey);
                        }
                    }
                }
            }
        }
    }
}

// Store the region triangles in a triangulation object and remove those
// triangles from the graph in preparation for processing the next layer of
// triangles.
function extractTriangles(graph: ETManifoldMesh, region: TriangleKeySet,
    node: PolygonTreeExNode): void {
    if (node.chirality > 0) {
        for (const tri of region.sorted()) {
            node.triangulation.push([tri.V[0], tri.V[1], tri.V[2]]);
            graph.remove(tri.V[0], tri.V[1], tri.V[2]);
        }
    }
    else {  // node.chirality < 0
        for (const tri of region.sorted()) {
            node.triangulation.push([tri.V[0], tri.V[2], tri.V[1]]);
            graph.remove(tri.V[0], tri.V[1], tri.V[2]);
        }
    }
}

function getOutsideTriangles(tree: PolygonTreeEx, graph: ETManifoldMesh): void {
    const tkeys = graph.getTriangleKeys();
    tree.outsideTriangles = tkeys.map(tkey =>
        [tkey.V[0], tkey.V[1], tkey.V[2]] as [number, number, number]);
    graph.clear();
}

// Get the triangles in the polygon tree, classifying each as interior (in
// the region bounded by an outer polygon and its contained inner polygons)
// or exterior (in the region bounded by an inner polygon and its contained
// outer polygons). The inside triangles are the union of the interior and
// exterior triangles.
function getInsideTriangles(tree: PolygonTreeEx): void {
    for (let nIndex = 0; nIndex < tree.nodes.length; ++nIndex) {
        const node = tree.nodes[nIndex];
        for (const tri of node.triangulation) {
            // Upstream stores std::array<int32_t,3> by value, so each of the
            // classified lists holds an independent copy. The port copies
            // explicitly; otherwise RestorePolygonTree would remap the shared
            // triangle several times.
            if (node.chirality > 0) {
                tree.interiorTriangles.push([tri[0], tri[1], tri[2]]);
                tree.interiorNodeIndices.push(nIndex);
            }
            else {
                tree.exteriorTriangles.push([tri[0], tri[1], tri[2]]);
                tree.exteriorNodeIndices.push(nIndex);
            }

            tree.insideTriangles.push([tri[0], tri[1], tri[2]]);
            tree.insideNodeIndices.push(nIndex);
        }
    }
}
