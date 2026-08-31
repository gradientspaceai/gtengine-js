// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) BVTree.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// BVTree is an abstract class for computing a bounding volume tree of a
// collection of primitives. The derived classes are BVTreeOfPoints (point
// primitives), BVTreeOfSegments (line segment primitives) and
// BVTreeOfTriangles (triangle primitives). In turn, derived classes of these
// classes create a bounding volume for each tree node.
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
// determined by the bounding volume type. The median of projections is chosen
// to partition the primitives into two subsets of equal size or absolute size
// difference of 1. This leads to a balanced tree, which is helpful for
// performance of tree traversals.
//
// Port notes:
//   - The upstream C++ template parameter 'BoundingVolume' is used both for
//     instance state (the member 'boundingVolume' of a node, and the instance
//     function GetSplittingAxis) and for static functions (IntersectLine,
//     IntersectRay, IntersectSegment) plus default construction. TypeScript
//     has no static-side constraint on a type parameter, so the port splits
//     the requirements: the instance requirement is the interface
//     'BVTreeBoundingVolume' (a type constraint on the parameter BV) and the
//     static requirements are the interface 'BVTreeVolumeOps<BV>', an object
//     passed to the BVTree constructor by the derived class.
//   - 'GetSplittingAxis(origin, direction)' has output reference parameters
//     upstream; the port returns a 'BVTreeSplittingAxis' object.
//   - The nested class 'BVTree<T, BoundingVolume>::Node' becomes the exported
//     top-level generic class 'BVTreeNode<BV>' (a nested generic class is not
//     expressible in TypeScript).
//   - 'Node::invalid' is Number.MAX_SAFE_INTEGER rather than SIZE_MAX; see
//     the same precedent in MeshStaticManifold2.ts and
//     SingularValueDecomposition.ts.
//   - The 'height' sentinel std::numeric_limits<std::size_t>::max(), meaning
//     "build the entire tree", is BVTree.fullHeight (also the default
//     argument of create).
//   - std::nth_element is a partial sort whose order within each of the two
//     resulting partitions is unspecified. The port uses a full sort by
//     projection, which satisfies the postconditions of nth_element and is
//     deterministic (Array.prototype.sort is stable, so centroids with equal
//     projections keep their relative order).
//   - The traversal stack in getLeafIndices uses the size_t wrap-around
//     'top--' at top = 0 as its loop terminator upstream; the port uses the
//     signed value -1 instead.
//   - create() copies the input centroids (upstream moves them, which is the
//     same observable result).

import { BitHacks } from './BitHacks';
import { logAssert } from './Logger';
import { Vector, dot, sub } from './Vector';

// The line/origin pair returned by BoundingVolume::GetSplittingAxis.
export interface BVTreeSplittingAxis {
    origin: Vector;
    direction: Vector;
}

// The instance-side requirement on the BoundingVolume type: it must supply a
// splitting axis, typically one in a direction of largest distribution of
// primitive vertices.
export interface BVTreeBoundingVolume {
    getSplittingAxis(): BVTreeSplittingAxis;
}

// Function signature for {line,ray,segment}-boundingVolume test-intersection
// queries. The line is parameterized by P+t*Q for all real t. The ray is
// parameterized by P+t*Q for nonnegative t. The segment is parameterized by
// (1-t)*P+t*Q for t in [0,1].
export type LinearBoundingVolumeQuery<BV> =
    (P: Vector, Q: Vector, boundingVolume: BV) => boolean;

// The static-side requirement on the BoundingVolume type: default
// construction and the three linear-component intersection queries.
export interface BVTreeVolumeOps<BV extends BVTreeBoundingVolume> {
    // The port of 'BoundingVolume()', used to initialize the nodes of the
    // preallocated node array.
    create(): BV;
    intersectLine(P: Vector, Q: Vector, boundingVolume: BV): boolean;
    intersectRay(P: Vector, Q: Vector, boundingVolume: BV): boolean;
    intersectSegment(P: Vector, Q: Vector, boundingVolume: BV): boolean;
}

// The port of the nested class BVTree<T, BoundingVolume>::Node.
export class BVTreeNode<BV> {
    static readonly invalid: number = Number.MAX_SAFE_INTEGER;

    boundingVolume: BV;
    minIndex: number;
    maxIndex: number;
    leftChild: number;
    rightChild: number;

    constructor(boundingVolume: BV) {
        this.boundingVolume = boundingVolume;
        this.minIndex = BVTreeNode.invalid;
        this.maxIndex = BVTreeNode.invalid;
        this.leftChild = BVTreeNode.invalid;
        this.rightChild = BVTreeNode.invalid;
    }
}

// The projection of a centroid onto the splitting axis, used for the median
// partitioning of a node's centroids.
interface ProjectionInfo {
    centroidIndex: number;
    projection: number;
}

// The port of std::nth_element; see the port notes above.
function sortByProjection(info: ProjectionInfo[]): void {
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

export abstract class BVTree<BV extends BVTreeBoundingVolume> {
    // These are the queryType inputs to the derived classes' execute
    // functions. Generate a list of leaf nodes intersected by a linear
    // component (line, ray or segment). The line is parameterized by
    // P + t * Q, where Q is a unit-length direction and t is any real
    // number. The ray is parameterized by P + t * Q, where Q is a
    // unit-length direction and t >= 0. The segment is parameterized by
    // (1-t) * P + t * Q = P + t * (Q - P), where P and Q are the endpoints
    // of the segment and 0 <= t <= 1.
    static readonly LINE_QUERY: number = 0;
    static readonly RAY_QUERY: number = 1;
    static readonly SEGMENT_QUERY: number = 2;

    // The 'height' input to create() that requests the entire tree; the port
    // of std::numeric_limits<std::size_t>::max().
    static readonly fullHeight: number = Number.MAX_SAFE_INTEGER;

    protected mCentroids: Vector[];
    protected mHeight: number;
    protected mNodes: BVTreeNode<BV>[];
    protected mPartition: number[];
    protected mOps: BVTreeVolumeOps<BV>;
    protected mLinearBoundingVolumeQuery: LinearBoundingVolumeQuery<BV>[];

    // Abstract base class. The derived classes must compute the centroids of
    // the primitives and pass them to the create(...) function.
    protected constructor(ops: BVTreeVolumeOps<BV>) {
        this.mCentroids = [];
        this.mHeight = 0;
        this.mNodes = [];
        this.mPartition = [];
        this.mOps = ops;
        this.mLinearBoundingVolumeQuery = [
            (P, Q, bv) => ops.intersectLine(P, Q, bv),
            (P, Q, bv) => ops.intersectRay(P, Q, bv),
            (P, Q, bv) => ops.intersectSegment(P, Q, bv)
        ];
    }

    // The derived classes must compute the centroids of the primitives and
    // pass them to the create(...) function.
    //
    // The input height specifies the desired height of the tree and must be
    // no larger than 31. If BVTree.fullHeight (the default), the entire tree
    // is built and the actual height is computed from centroids.length. If
    // larger than 31, the height is clamped to 31.
    create(centroids: readonly Vector[], height: number = BVTree.fullHeight): void {
        logAssert(centroids.length > 0,
            'Expecting centroids to create a bounding volume tree.');

        this.mCentroids = new Array<Vector>(centroids.length);
        for (let i = 0; i < centroids.length; ++i) {
            logAssert(centroids[i].size === 3, 'BVTree: centroids must be 3D.');
            this.mCentroids[i] = centroids[i].clone();
        }

        if (height === BVTree.fullHeight) {
            const minPowerOfTwo = BitHacks.roundUpToPowerOfTwo(this.mCentroids.length);
            this.mHeight = BitHacks.log2OfPowerOfTwo(minPowerOfTwo);
        } else {
            this.mHeight = Math.min(height, 31);
        }

        // The tree is built recursively. Preallocate the nodes; a complete
        // binary tree of height H has 2^{H+1}-1 nodes.
        const numNodes = 2 ** (this.mHeight + 1) - 1;
        this.mNodes = new Array<BVTreeNode<BV>>(numNodes);
        for (let i = 0; i < numNodes; ++i) {
            this.mNodes[i] = new BVTreeNode<BV>(this.mOps.create());
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

    getNodes(): readonly BVTreeNode<BV>[] {
        return this.mNodes;
    }

    getPartition(): readonly number[] {
        return this.mPartition;
    }

    // The bounding volume for the primitives' vertices depends on the type of
    // primitive. A derived class representing a primitive tree must implement
    // this. The bounding volume is modified in place, as it is upstream.
    protected abstract computeInteriorBoundingVolume(i0: number, i1: number,
        boundingVolume: BV): void;

    // The bounding volume for a single primitive's vertices depends on the
    // type of primitive. A derived class representing a primitive tree must
    // implement this.
    protected abstract computeLeafBoundingVolume(i: number, boundingVolume: BV): void;

    // Get the node indices for the leaf nodes whose bounding volumes are
    // intersected by the linear component.
    protected getLeafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        const nodeIndices: number[] = [];

        const linearBoundingVolumeQuery = this.mLinearBoundingVolumeQuery[queryType];
        const indexStack = new Array<number>(2 * this.mHeight + 1).fill(0);
        let top = 0;
        indexStack[0] = 0;
        while (top >= 0) {
            const nodeIndex = indexStack[top--];
            const node = this.mNodes[nodeIndex];

            // For the balanced tree created by BVTree, an interior node has
            // two valid children and a leaf node has two invalid children.
            // This is true even if the height passed to create is smaller
            // than the actual height.
            if (node.leftChild !== BVTreeNode.invalid &&
                node.rightChild !== BVTreeNode.invalid) {
                // The node is interior.
                if (linearBoundingVolumeQuery(P, Q, node.boundingVolume)) {
                    // The linear component intersects the bounding volume.
                    // Continue the intersection search to child nodes if they
                    // exist.
                    indexStack[++top] = node.rightChild;
                    indexStack[++top] = node.leftChild;
                }
                // Otherwise the linear component does not intersect the
                // bounding volume, so do not continue the intersection search
                // to child nodes if they exist.
            } else {
                // node.leftChild == invalid && node.rightChild == invalid.
                // NOTE (upstream quirk, preserved): the leaf's own bounding
                // volume is not tested against the linear component, so a
                // leaf is reported whenever its parent's bounding volume is
                // intersected. The result is conservative (a superset of the
                // truly intersected leaves).
                nodeIndices.push(nodeIndex);
            }
        }

        return nodeIndices;
    }

    // Support for tree creation.
    private buildTree(depth: number, nodeIndex: number, i0: number, i1: number): void {
        const node = this.mNodes[nodeIndex];
        node.minIndex = i0;
        node.maxIndex = i1;

        if (i0 < i1) {
            // The node is interior. Compute a bounding volume for the
            // primitives' vertices.
            this.computeInteriorBoundingVolume(i0, i1, node.boundingVolume);
            if (depth === this.mHeight) {
                // The user-specified height has been reached. Do not continue
                // the recursion past this node.
                return;
            }

            // The BoundingVolume type provides a function to access a
            // splitting axis, typically one in a direction of largest
            // distribution of primitive vertices. Use the splitting axis to
            // partition the centroids of the primitives into two subsets, one
            // for the left child and one for the right child. The subsets
            // have numbers of elements that differ by at most 1, so the tree
            // is balanced.
            const { j0, j1 } = this.splitPoints(i0, i1, node.boundingVolume);

            // Recurse on the two children.
            node.leftChild = 2 * nodeIndex + 1;
            node.rightChild = node.leftChild + 1;
            this.buildTree(depth + 1, node.leftChild, i0, j0);
            this.buildTree(depth + 1, node.rightChild, j1, i1);
        } else {
            // i0 = i1. The node is a leaf. Compute a bounding volume for a
            // single primitive's vertices.
            this.computeLeafBoundingVolume(i0, node.boundingVolume);
        }
    }

    private splitPoints(i0: number, i1: number, boundingVolume: BV):
        { j0: number; j1: number } {
        // The direction of the splitting axis is provided by the
        // BoundingVolume type.
        const { origin, direction } = boundingVolume.getSplittingAxis();

        // Project the centroids onto the splitting axis.
        const numProjections = i1 - i0 + 1;
        const info = new Array<ProjectionInfo>(numProjections);
        for (let i = i0, j = 0; i <= i1; ++i, ++j) {
            const centroidIndex = this.mPartition[i];
            const diff = sub(this.mCentroids[centroidIndex], origin);
            info[j] = { centroidIndex: centroidIndex, projection: dot(direction, diff) };
        }

        // Partition the projections by the median.
        const medianIndex = Math.floor((numProjections - 1) / 2);
        sortByProjection(info);

        // Partition the centroids by the median. When i0 = 0, the initial j0
        // is -1 (upstream, the wrap-around maximum of std::size_t). However,
        // it is incremented to 0 before the lookup into mPartition[].
        let k = 0;
        let j0 = i0 - 1;
        for (k = 0; k <= medianIndex; ++k) {
            this.mPartition[++j0] = info[k].centroidIndex;
        }
        let j1 = i1 + 1;
        for (; k < numProjections; ++k) {
            this.mPartition[--j1] = info[k].centroidIndex;
        }

        return { j0: j0, j1: j1 };
    }
}
