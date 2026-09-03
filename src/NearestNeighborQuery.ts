// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) NearestNeighborQuery.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Upstream TODO (preserved): This is not a KD-tree nearest neighbor query.
// Instead, it is an algorithm to get "approximate" nearest neighbors.
// Replace this by the actual KD-tree query.

// Use a kd-tree for sorting used in a query for finding nearest neighbors of
// a point in a space of the specified dimension N. The split order is always
// 0,1,2,...,N-1. The number of sites at a leaf node is controlled by
// 'maxLeafSize' and the maximum level of the tree is controlled by
// 'maxLevel'. The points are of type Vector. The 'Site' is a structure of
// information that minimally implements the function 'getPosition(): Vector'.
// The Site abstraction allows the query to be applied even when it has more
// local information than just point location.
//
// Port notes: the compile-time template parameters <N, T, Site> become
// runtime state; N is inferred from the sites' position dimension. The
// template parameter MaxNeighbors of FindNeighbors becomes the runtime
// argument maxNeighbors, and the (return count, out std::array) pair becomes
// a returned array of indices whose length is the neighbor count. The
// std::pair<Vector<N,T>, int32_t> SortedPoint becomes the
// NearestNeighborSortedPoint interface, and the internal Node struct is
// exported as NearestNeighborNode. std::nth_element is replicated by a
// quickselect that establishes the identical partition postconditions, and
// std::priority_queue<std::pair<T, int32_t>> is replicated by a binary
// max-heap over the same lexicographic pair ordering, so the retained
// neighbor set and the output order are those of the upstream code.

import { logAssert } from './Logger.js';
import { Vector, dot, sub } from './Vector.js';

// The site abstraction: anything that can report its position.
export interface NearestNeighborSite {
    getPosition(): Vector;
}

// Predefined site structs for convenience.
export class PositionSite implements NearestNeighborSite {
    position: Vector;

    constructor(p: Vector) {
        this.position = p.clone();
    }

    getPosition(): Vector {
        return this.position.clone();
    }
}

// Predefined site structs for convenience.
export class PositionDirectionSite implements NearestNeighborSite {
    position: Vector;
    direction: Vector;

    constructor(p: Vector, d: Vector) {
        this.position = p.clone();
        this.direction = d.clone();
    }

    getPosition(): Vector {
        return this.position.clone();
    }
}

// The port of 'using SortedPoint = std::pair<Vector<N, T>, int32_t>'.
export interface NearestNeighborSortedPoint {
    position: Vector;
    index: number;
}

// The port of the nested struct NearestNeighborQuery::Node.
export interface NearestNeighborNode {
    split: number;
    axis: number;
    numSites: number;
    siteOffset: number;
    left: number;
    right: number;
}

// The port of 'using VIPair = std::pair<T, int32_t>'.
interface VIPair {
    first: number;
    second: number;
}

// std::pair<T, int32_t> operator< (lexicographic).
function viPairLess(p0: VIPair, p1: VIPair): boolean {
    return p0.first < p1.first || (p0.first === p1.first && p0.second < p1.second);
}

// The port of std::priority_queue<VIPair> (a binary max-heap ordered by
// std::less on pairs). The site indices in the heap are distinct, so the
// pair ordering is a strict total order and the top/pop sequence is exactly
// that of the upstream priority queue.
class NNPriorityQueue {
    private heap: VIPair[] = [];

    size(): number {
        return this.heap.length;
    }

    empty(): boolean {
        return this.heap.length === 0;
    }

    top(): VIPair {
        return this.heap[0];
    }

    push(element: VIPair): void {
        const heap = this.heap;
        heap.push(element);
        let child = heap.length - 1;
        while (child > 0) {
            const parent = (child - 1) >> 1;
            if (viPairLess(heap[parent], heap[child])) {
                const temp = heap[parent];
                heap[parent] = heap[child];
                heap[child] = temp;
                child = parent;
            } else {
                break;
            }
        }
    }

    pop(): void {
        const heap = this.heap;
        const last = heap.pop() as VIPair;
        if (heap.length === 0) {
            return;
        }
        heap[0] = last;
        let parent = 0;
        for (;;) {
            const left = 2 * parent + 1;
            const right = left + 1;
            let largest = parent;
            if (left < heap.length && viPairLess(heap[largest], heap[left])) {
                largest = left;
            }
            if (right < heap.length && viPairLess(heap[largest], heap[right])) {
                largest = right;
            }
            if (largest === parent) {
                break;
            }
            const temp = heap[parent];
            heap[parent] = heap[largest];
            heap[largest] = temp;
            parent = largest;
        }
    }
}

// The median of three values under 'less'.
function medianOf3<Element>(x: Element, y: Element, z: Element,
    less: (a: Element, b: Element) => boolean): Element {
    if (less(x, y)) {
        if (less(y, z)) {
            return y;
        }
        return less(x, z) ? z : x;
    }
    if (less(x, z)) {
        return x;
    }
    return less(y, z) ? z : y;
}

// The port of std::nth_element on the subrange [first, last) of 'a': after
// the call, a[nth] is the element that would be in that position were the
// subrange sorted, no element of [first, nth) is greater than a[nth], and no
// element of [nth, last) is less than a[nth]. Implemented as a quickselect
// with a median-of-three pivot and a three-way partition; the equal-pivot
// band is nonempty at every step, which guarantees termination.
function nthElement<Element>(a: Element[], first: number, nth: number, last: number,
    less: (p0: Element, p1: Element) => boolean): void {
    while (last - first > 1) {
        const mid = first + ((last - first) >> 1);
        const pivot = medianOf3(a[first], a[mid], a[last - 1], less);

        // Three-way partition: [first, lt) < pivot, [lt, gt) == pivot,
        // [gt, last) > pivot.
        let lt = first, gt = last, i = first;
        while (i < gt) {
            if (less(a[i], pivot)) {
                const temp = a[lt];
                a[lt] = a[i];
                a[i] = temp;
                ++lt;
                ++i;
            } else if (less(pivot, a[i])) {
                --gt;
                const temp = a[gt];
                a[gt] = a[i];
                a[i] = temp;
            } else {
                ++i;
            }
        }

        if (nth < lt) {
            last = lt;
        } else if (nth < gt) {
            return;
        } else {
            first = gt;
        }
    }
}

export class NearestNeighborQuery {
    private mMaxLeafSize: number;
    private mMaxLevel: number;
    private mDimension: number;
    private mSortedPoints: NearestNeighborSortedPoint[];
    private mNodes: NearestNeighborNode[];
    private mDepth: number;
    private mLargestNodeSize: number;

    // Construction. The tree dimension N is the dimension of the site
    // positions.
    constructor(sites: NearestNeighborSite[], maxLeafSize: number, maxLevel: number) {
        this.mMaxLeafSize = maxLeafSize;
        this.mMaxLevel = maxLevel;
        this.mSortedPoints = new Array<NearestNeighborSortedPoint>(sites.length);
        this.mNodes = [];
        this.mDepth = 0;
        this.mLargestNodeSize = 0;

        logAssert(this.mMaxLevel > 0 && this.mMaxLevel <= 32, 'Invalid max level.');

        const numSites = sites.length;
        this.mDimension = numSites > 0 ? sites[0].getPosition().getSize() : 0;
        for (let i = 0; i < numSites; ++i) {
            this.mSortedPoints[i] = { position: sites[i].getPosition(), index: i };
        }

        this.mNodes.push({ split: 0, axis: 0, numSites: 0, siteOffset: 0, left: 0, right: 0 });
        this.build(numSites, 0, 0, 0);
    }

    // Member access.
    getMaxLeafSize(): number {
        return this.mMaxLeafSize;
    }

    getMaxLevel(): number {
        return this.mMaxLevel;
    }

    getDepth(): number {
        return this.mDepth;
    }

    getLargestNodeSize(): number {
        return this.mLargestNodeSize;
    }

    getNumNodes(): number {
        return this.mNodes.length;
    }

    getNodes(): readonly NearestNeighborNode[] {
        return this.mNodes;
    }

    getSortedPoints(): readonly NearestNeighborSortedPoint[] {
        return this.mSortedPoints;
    }

    // Compute up to maxNeighbors nearest neighbors within the specified
    // radius of the point. The returned array stores indices into the sites
    // array passed to the constructor; its length is the number of neighbors
    // found, possibly zero. The indices are in the upstream output order
    // (heap pop order: nonincreasing squared distance, ties by decreasing
    // site index).
    findNeighbors(point: Vector, radius: number, maxNeighbors: number): number[] {
        logAssert(maxNeighbors >= 1, 'Invalid maximum number of neighbors.');

        const sqrRadius = radius * radius;
        const maxHeap = new NNPriorityQueue();

        // The kd-tree construction is recursive, simulated here by using a
        // stack. The maximum depth is limited to 32, because the number of
        // sites is limited to 2^{32} (the number of 32-bit integer indices).
        const stack = new Array<number>(32).fill(0);
        let top = 0;
        stack[0] = 0;

        while (top >= 0) {
            const node = this.mNodes[stack[top--]];

            if (node.siteOffset !== -1) {
                for (let i = 0, j = node.siteOffset; i < node.numSites; ++i, ++j) {
                    const diff = sub(this.mSortedPoints[j].position, point);
                    const sqrLength = dot(diff, diff);
                    if (sqrLength <= sqrRadius) {
                        // Keep track of the nearest neighbors.
                        if (maxHeap.size() < maxNeighbors) {
                            maxHeap.push({ first: sqrLength, second: this.mSortedPoints[j].index });
                        } else if (sqrLength < maxHeap.top().first) {
                            maxHeap.pop();
                            maxHeap.push({ first: sqrLength, second: this.mSortedPoints[j].index });
                        }
                    }
                }
            }

            if (node.left !== -1 && point.get(node.axis) - radius <= node.split) {
                stack[++top] = node.left;
            }

            if (node.right !== -1 && point.get(node.axis) + radius >= node.split) {
                stack[++top] = node.right;
            }
        }

        const neighbors = new Array<number>(maxHeap.size());
        let nidx = 0;
        while (!maxHeap.empty()) {
            neighbors[nidx++] = maxHeap.top().second;
            maxHeap.pop();
        }

        return neighbors;
    }

    // Populate the node so that it contains the points split along the
    // coordinate axes.
    private build(numSites: number, siteOffset: number, nodeIndex: number, level: number): void {
        logAssert(siteOffset !== -1, 'Invalid site offset.');
        logAssert(nodeIndex !== -1, 'Invalid node index.');
        logAssert(numSites > 0, 'Empty point list.');

        this.mDepth = Math.max(this.mDepth, level);

        const node = this.mNodes[nodeIndex];
        node.numSites = numSites;

        if (numSites > this.mMaxLeafSize && level <= this.mMaxLevel) {
            const halfNumSites = Math.trunc(numSites / 2);

            // The point set is too large for a leaf node, so split it at the
            // median. The O(m log m) sort is not needed; rather, we locate
            // the median using an order statistic construction that is
            // expected time O(m).
            const axis = level % this.mDimension;
            const sorter = (p0: NearestNeighborSortedPoint, p1: NearestNeighborSortedPoint) =>
                p0.position.get(axis) < p1.position.get(axis);

            nthElement(this.mSortedPoints, siteOffset, siteOffset + halfNumSites,
                siteOffset + numSites, sorter);

            // Get the median position.
            const index = siteOffset + halfNumSites;
            node.split = this.mSortedPoints[index].position.get(axis);
            node.axis = axis;
            node.siteOffset = -1;

            // Apply a divide-and-conquer step.
            const left = this.mNodes.length;
            const right = left + 1;
            node.left = left;
            node.right = right;
            this.mNodes.push({ split: 0, axis: 0, numSites: 0, siteOffset: 0, left: 0, right: 0 });
            this.mNodes.push({ split: 0, axis: 0, numSites: 0, siteOffset: 0, left: 0, right: 0 });

            const nextLevel = level + 1;
            this.build(halfNumSites, siteOffset, left, nextLevel);
            this.build(numSites - halfNumSites, siteOffset + halfNumSites, right, nextLevel);
        } else {
            // The number of points is small enough or we have run out of
            // depth, so make this node a leaf.
            node.split = Number.MAX_VALUE;
            node.axis = -1;
            node.siteOffset = siteOffset;
            node.left = -1;
            node.right = -1;

            this.mLargestNodeSize = Math.max(this.mLargestNodeSize, node.numSites);
        }
    }
}
