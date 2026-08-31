// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) OBBTree.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// OBBTree is an abstract class for computing an oriented bounding box tree of
// a collection of primitives. The derived classes are OBBTreeOfPoints (point
// primitives), OBBTreeOfSegments (line segment primitives) and
// OBBTreeOfTriangles (triangle primitives). The derived classes create a box
// for each tree node. The box center is the mean of centroids of the
// primitives that the node represents. The box axis directions are the
// eigenvectors of the covariance matrix of those centroids. The box extents
// are computed to ensure the box contains the primitives represented by the
// node.
//
// The depth of a node in a nonempty tree is the distance from the node to the
// root of the tree. The height is the maximum depth. A tree with a single
// node has height 0. The set of nodes of a tree with the same depth is
// referred to as a level of a tree corresponding to that depth. A complete
// binary tree of height H has 2^{H+1}-1 nodes. The level corresponding to
// depth D has 2^D nodes, in which case the number of leaf nodes (nodes at
// depth H) is 2^H.
//
// The partitioning of primitives between left and right children of a node
// is based on the projection of centroids of the primitives onto a line
// determined by eigenvectors corresponding to the largest eigenvalue of
// covariance matrices. The median of projections is chosen to partition the
// primitives into two subsets of equal size or absolute size difference of 1.
// This leads to a balanced tree, which is helpful for performance of tree
// traversals.
//
// Port notes:
//   - Upstream OBBTree is independent of BVTree; it is the older, oriented-
//     box-specific design and does not derive from BVTree. The port keeps
//     that structure, including the duplicated tree-building code.
//   - 'OBBNode::minIndex/maxIndex/leftChild/rightChild' use the sentinel
//     std::numeric_limits<size_t>::max() upstream; the port uses
//     OBBNode.invalid = Number.MAX_SAFE_INTEGER, following the precedent in
//     MeshStaticManifold2.ts and SingularValueDecomposition.ts.
//   - The 'height' sentinel std::numeric_limits<size_t>::max(), meaning
//     "build the entire tree", is OBBTree.fullHeight (also the default
//     argument of create).
//   - std::nth_element is a partial sort whose order within each of the two
//     resulting partitions is unspecified. The port uses a full sort by
//     projection, which satisfies the postconditions of nth_element and is
//     deterministic (Array.prototype.sort is stable, so centroids with equal
//     projections keep their relative order).
//   - The output reference parameters j0 and j1 of SplitPoints become fields
//     of a returned object.
//   - OrientedBox3<T> is 'new OrientedBox(3)' here (runtime dimension); see
//     OrientedBox.ts.
//   - The eigensolver used is SymmetricEigensolver3x3, whose port fixed an
//     upstream eigenvector bug, so box axes on degenerate (rank-deficient)
//     centroid sets can differ from the buggy upstream axes.

import { BitHacks } from './BitHacks';
import { logAssert } from './Logger';
import { OrientedBox } from './OrientedBox';
import { SymmetricEigensolver3x3 } from './SymmetricEigensolver3x3';
import { Vector, dot, sub } from './Vector';

// The port of the upstream struct OBBNode<T>.
export class OBBNode {
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    box: OrientedBox;
    minIndex: number;
    maxIndex: number;
    leftChild: number;
    rightChild: number;

    constructor() {
        this.box = new OrientedBox(3);
        this.minIndex = OBBNode.invalid;
        this.maxIndex = OBBNode.invalid;
        this.leftChild = OBBNode.invalid;
        this.rightChild = OBBNode.invalid;
    }
}

// The projection of a centroid onto the splitting axis, used for the median
// partitioning of a node's centroids.
interface OBBProjectionInfo {
    pointIndex: number;
    projection: number;
}

// The port of std::nth_element; see the port notes above.
function sortByProjection(info: OBBProjectionInfo[]): void {
    info.sort((info0, info1) => {
        if (info0.projection < info1.projection) {
            return -1;
        }
        if (info1.projection < info0.projection) {
            return +1;
        }
        return 0;
    });
}

export abstract class OBBTree {
    // The 'height' input to create() that requests the entire tree; the port
    // of std::numeric_limits<size_t>::max().
    static readonly fullHeight: number = Number.MAX_SAFE_INTEGER;

    protected mCentroids: Vector[];
    protected mHeight: number;
    protected mNodes: OBBNode[];
    protected mPartition: number[];

    // Abstract base class.
    protected constructor() {
        this.mCentroids = [];
        this.mHeight = 0;
        this.mNodes = [];
        this.mPartition = [];
    }

    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If OBBTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from centroids.length. If
    // larger than 31, the height is clamped to 31.
    create(centroids: readonly Vector[], height: number = OBBTree.fullHeight): void {
        logAssert(centroids.length > 0, 'Invalid input.');

        this.mCentroids = new Array<Vector>(centroids.length);
        for (let i = 0; i < centroids.length; ++i) {
            logAssert(centroids[i].size === 3, 'OBBTree: centroids must be 3D.');
            this.mCentroids[i] = centroids[i].clone();
        }

        if (height === OBBTree.fullHeight) {
            const minPowerOfTwo = BitHacks.roundUpToPowerOfTwo(this.mCentroids.length);
            this.mHeight = BitHacks.log2OfPowerOfTwo(minPowerOfTwo);
        } else {
            this.mHeight = Math.min(height, 31);
        }

        // The tree is built recursively. Preallocate the nodes; a complete
        // binary tree of height H has 2^{H+1}-1 nodes.
        const numNodes = 2 ** (this.mHeight + 1) - 1;
        this.mNodes = new Array<OBBNode>(numNodes);
        for (let i = 0; i < numNodes; ++i) {
            this.mNodes[i] = new OBBNode();
        }

        // The array mPartition stores indices into mCentroids so that at a
        // node, the centroids represented by the node are the indices
        // [mPartition[node.minIndex], mPartition[node.maxIndex]].
        this.mPartition = new Array<number>(this.mCentroids.length);
        for (let i = 0; i < this.mPartition.length; ++i) {
            this.mPartition[i] = i;
        }

        // Build the tree recursively.
        const depth = 0;
        const nodeIndex = 0;
        const i0 = 0;
        const i1 = this.mCentroids.length - 1;
        this.buildTree(depth, nodeIndex, i0, i1);
    }

    // Member access.
    getCentroids(): readonly Vector[] {
        return this.mCentroids;
    }

    getHeight(): number {
        return this.mHeight;
    }

    getNodes(): readonly OBBNode[] {
        return this.mNodes;
    }

    getPartition(): readonly number[] {
        return this.mPartition;
    }

    // The derived classes must override computeInteriorBox, calling the base
    // class function first. They must then compute the box extents to ensure
    // the box contains the primitives represented by the node.
    protected computeInteriorBox(i0: number, i1: number, box: OrientedBox): void {
        // Compute the mean of the centroids.
        box.center = new Vector(3);
        for (let i = i0; i <= i1; ++i) {
            const centroid = this.mCentroids[this.mPartition[i]];
            box.center.values[0] += centroid.values[0];
            box.center.values[1] += centroid.values[1];
            box.center.values[2] += centroid.values[2];
        }
        const denom = i1 - i0 + 1;
        box.center.values[0] /= denom;
        box.center.values[1] /= denom;
        box.center.values[2] /= denom;

        // Compute the covariance matrix of the centroids.
        let covar00 = 0, covar01 = 0, covar02 = 0;
        let covar11 = 0, covar12 = 0, covar22 = 0;
        for (let i = i0; i <= i1; ++i) {
            const diff = sub(this.mCentroids[this.mPartition[i]], box.center);
            covar00 += diff.values[0] * diff.values[0];
            covar01 += diff.values[0] * diff.values[1];
            covar02 += diff.values[0] * diff.values[2];
            covar11 += diff.values[1] * diff.values[1];
            covar12 += diff.values[1] * diff.values[2];
            covar22 += diff.values[2] * diff.values[2];
        }
        covar00 /= denom;
        covar01 /= denom;
        covar02 /= denom;
        covar11 /= denom;
        covar12 /= denom;
        covar22 /= denom;

        // Use the eigenvectors of the covariance matrix for the box axes. The
        // sort type +1 orders the eigenvalues increasingly, so axis[2]
        // corresponds to the largest eigenvalue.
        const es = new SymmetricEigensolver3x3();
        const result = es.solve(covar00, covar01, covar02, covar11, covar12,
            covar22, false, +1);
        for (let i = 0; i < 3; ++i) {
            box.axis[i] = Vector.fromArray(result.evecs[i]);
        }

        // The box.extent values must be computed by the derived classes. For
        // debugging, store the eigenvalues in the extents.
        box.extent = Vector.fromArray(result.evals);
    }

    // The derived classes must override computeLeafBox. The intrinsic box
    // dimension depends on the geometric primitive.
    protected abstract computeLeafBox(i: number, box: OrientedBox): void;

    private buildTree(depth: number, nodeIndex: number, i0: number, i1: number): void {
        const node = this.mNodes[nodeIndex];
        node.minIndex = i0;
        node.maxIndex = i1;

        if (i0 < i1) {
            // The node is interior. Compute an oriented bounding box of
            // centroids, but then with extents modified to ensure the box
            // contains the primitives represented by the node.
            this.computeInteriorBox(i0, i1, node.box);
            if (depth === this.mHeight) {
                // The user-specified height has been reached. Do not continue
                // the recursion past this node.
                return;
            }

            // Use the box axis corresponding to largest extent for the
            // splitting axis. Partition the centroids into two subsets, one
            // for the left child and one for the right child. The subsets
            // have numbers of elements that differ by at most 1, so the tree
            // is balanced.
            const { j0, j1 } = this.splitPoints(i0, i1, node.box.center,
                node.box.axis[2]);

            // Recurse on the two children.
            node.leftChild = 2 * nodeIndex + 1;
            node.rightChild = node.leftChild + 1;
            this.buildTree(depth + 1, node.leftChild, i0, j0);
            this.buildTree(depth + 1, node.rightChild, j1, i1);
        } else {
            // i0 = i1. The node is a leaf. Compute a primitive-dependent
            // oriented bounding box.
            this.computeLeafBox(i0, node.box);
        }
    }

    private splitPoints(i0: number, i1: number, origin: Vector, direction: Vector):
        { j0: number; j1: number } {
        // Project the centroids onto the splitting axis.
        const numProjections = i1 - i0 + 1;
        const info = new Array<OBBProjectionInfo>(numProjections);
        for (let i = i0, j = 0; i <= i1; ++i, ++j) {
            const pointIndex = this.mPartition[i];
            const diff = sub(this.mCentroids[pointIndex], origin);
            info[j] = { pointIndex: pointIndex, projection: dot(direction, diff) };
        }

        // Partition the projections by the median.
        const medianIndex = Math.floor((numProjections - 1) / 2);
        sortByProjection(info);

        // Partition the centroids by the median. When i0 = 0, the initial j0
        // is -1 (upstream, the wrap-around maximum of size_t). However, it is
        // incremented to 0 before the lookup into mPartition[].
        let k = 0;
        let j0 = i0 - 1;
        for (k = 0; k <= medianIndex; ++k) {
            this.mPartition[++j0] = info[k].pointIndex;
        }
        let j1 = i1 + 1;
        for (; k < numProjections; ++k) {
            this.mPartition[--j1] = info[k].pointIndex;
        }

        return { j0: j0, j1: j1 };
    }
}
