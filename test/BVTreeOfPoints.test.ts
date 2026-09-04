import { describe, it, expect } from 'vitest';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from '../src/BVTree.js';
import { BVTreeOfPoints } from '../src/BVTreeOfPoints.js';
import { Vector, sub } from '../src/Vector.js';
import {
    check, fc, unitVector, vector
} from './helpers/arbitraries.js';

// ---------------------------------------------------------------------------
// A concrete BoundingVolume: an axis-aligned bounding box, the same minimal
// harness used by BVTree.test.ts.
// ---------------------------------------------------------------------------

class AABB implements BVTreeBoundingVolume {
    min: number[] = [Infinity, Infinity, Infinity];
    max: number[] = [-Infinity, -Infinity, -Infinity];

    reset(): void {
        this.min = [Infinity, Infinity, Infinity];
        this.max = [-Infinity, -Infinity, -Infinity];
    }

    grow(p: Vector): void {
        for (let j = 0; j < 3; ++j) {
            if (p.values[j] < this.min[j]) {
                this.min[j] = p.values[j];
            }
            if (p.values[j] > this.max[j]) {
                this.max[j] = p.values[j];
            }
        }
    }

    contains(p: Vector, eps: number = 1e-12): boolean {
        for (let j = 0; j < 3; ++j) {
            if (p.values[j] < this.min[j] - eps || p.values[j] > this.max[j] + eps) {
                return false;
            }
        }
        return true;
    }

    getSplittingAxis(): BVTreeSplittingAxis {
        const origin = Vector.fromArray([
            0.5 * (this.min[0] + this.max[0]),
            0.5 * (this.min[1] + this.max[1]),
            0.5 * (this.min[2] + this.max[2])
        ]);
        let jmax = 0;
        let emax = this.max[0] - this.min[0];
        for (let j = 1; j < 3; ++j) {
            const e = this.max[j] - this.min[j];
            if (e > emax) {
                emax = e;
                jmax = j;
            }
        }
        return { origin: origin, direction: Vector.unit(3, jmax) };
    }
}

// The slab test for a linear component P + t * D restricted to [tmin, tmax].
function intersectSlabs(P: Vector, D: Vector, tmin0: number, tmax0: number,
    box: AABB): boolean {
    let tmin = tmin0;
    let tmax = tmax0;
    for (let j = 0; j < 3; ++j) {
        if (Math.abs(D.values[j]) < 1e-15) {
            if (P.values[j] < box.min[j] || P.values[j] > box.max[j]) {
                return false;
            }
        } else {
            let t0 = (box.min[j] - P.values[j]) / D.values[j];
            let t1 = (box.max[j] - P.values[j]) / D.values[j];
            if (t0 > t1) {
                const t = t0;
                t0 = t1;
                t1 = t;
            }
            tmin = Math.max(tmin, t0);
            tmax = Math.min(tmax, t1);
            if (tmin > tmax) {
                return false;
            }
        }
    }
    return true;
}

const aabbOps: BVTreeVolumeOps<AABB> = {
    create: () => new AABB(),
    intersectLine: (P, Q, bv) => intersectSlabs(P, Q, -Infinity, Infinity, bv),
    intersectRay: (P, Q, bv) => intersectSlabs(P, Q, 0, Infinity, bv),
    intersectSegment: (P, Q, bv) => intersectSlabs(P, sub(Q, P), 0, 1, bv)
};

// The concrete point tree: the bounding volume of a node contains the points
// the node represents.
class AABBTreeOfPoints extends BVTreeOfPoints<AABB> {
    constructor() {
        super(aabbOps);
    }

    protected computeInteriorBoundingVolume(i0: number, i1: number, bv: AABB): void {
        bv.reset();
        for (let i = i0; i <= i1; ++i) {
            bv.grow(this.mVertices[this.mPartition[i]]);
        }
    }

    protected computeLeafBoundingVolume(i: number, bv: AABB): void {
        bv.reset();
        bv.grow(this.mVertices[this.mPartition[i]]);
    }
}

function cubeCorners(): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < 8; ++i) {
        points.push(Vector.fromArray([
            (i & 1) ? 1 : 0,
            (i & 2) ? 2 : 0,
            (i & 4) ? 4 : 0
        ]));
    }
    return points;
}

function isLeaf(node: BVTreeNode<AABB>): boolean {
    return node.leftChild === BVTreeNode.invalid
        && node.rightChild === BVTreeNode.invalid;
}

describe('BVTreeOfPoints', () => {
    it('builds a complete balanced tree for 8 points', () => {
        const tree = new AABBTreeOfPoints();
        const points = cubeCorners();
        tree.create(points);

        // 8 = 2^3 leaves, so the height is 3 and there are 2^4-1 = 15 nodes.
        expect(tree.getHeight()).toBe(3);
        expect(tree.getNodes().length).toBe(15);

        const root = tree.getNodes()[0];
        expect(root.minIndex).toBe(0);
        expect(root.maxIndex).toBe(7);
        expect(isLeaf(root)).toBe(false);

        // The partition is a permutation of the point indices.
        const partition = tree.getPartition();
        expect(partition.length).toBe(8);
        expect([...partition].sort((a, b) => a - b)).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);

        // The leaves are exactly the 8 singleton ranges.
        const leafRanges: number[] = [];
        for (const node of tree.getNodes()) {
            if (isLeaf(node)) {
                expect(node.minIndex).toBe(node.maxIndex);
                leafRanges.push(node.minIndex);
            }
        }
        expect(leafRanges.sort((a, b) => a - b)).toEqual([0, 1, 2, 3, 4, 5, 6, 7]);
    });

    it('the points are also the centroids', () => {
        const tree = new AABBTreeOfPoints();
        const points = cubeCorners();
        tree.create(points);

        const centroids = tree.getCentroids();
        expect(centroids.length).toBe(points.length);
        for (let i = 0; i < points.length; ++i) {
            expect(centroids[i].values).toEqual(points[i].values);
            expect(tree.getVertices()[i].values).toEqual(points[i].values);
        }
    });

    it('copies the input vertices', () => {
        const tree = new AABBTreeOfPoints();
        const points = cubeCorners();
        tree.create(points);
        points[0].values[0] = 100;
        expect(tree.getVertices()[0].values[0]).toBe(0);
    });

    it('each node bounding volume contains the points it represents', () => {
        const tree = new AABBTreeOfPoints();
        // A pseudorandom cloud (a deterministic linear congruential sequence).
        let seed = 12345;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const points: Vector[] = [];
        for (let i = 0; i < 37; ++i) {
            points.push(Vector.fromArray([rand(), 3 * rand(), -2 * rand()]));
        }
        tree.create(points);

        const nodes = tree.getNodes();
        const partition = tree.getPartition();
        for (const node of nodes) {
            if (node.minIndex === BVTreeNode.invalid) {
                // A node of the complete array that the tree does not use.
                continue;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                expect(node.boundingVolume.contains(points[partition[i]])).toBe(true);
            }
        }

        // Sibling ranges partition the parent range.
        for (const node of nodes) {
            if (isLeaf(node) || node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            const left = nodes[node.leftChild];
            const right = nodes[node.rightChild];
            expect(left.minIndex).toBe(node.minIndex);
            expect(right.maxIndex).toBe(node.maxIndex);
            expect(right.minIndex).toBe(left.maxIndex + 1);
        }
    });

    it('execute reports a superset of the leaves whose points are hit', () => {
        const tree = new AABBTreeOfPoints();
        const points = cubeCorners();
        tree.create(points);

        // The line through (0,0,0) and (1,2,4) contains two of the corners.
        const P = Vector.fromArray([0, 0, 0]);
        const Q = Vector.fromArray([1, 2, 4]);
        const nodeIndices = tree.execute(BVTree.LINE_QUERY, P, Q);

        // Every reported index is a leaf and they are distinct.
        expect(new Set(nodeIndices).size).toBe(nodeIndices.length);
        const nodes = tree.getNodes();
        for (const index of nodeIndices) {
            expect(isLeaf(nodes[index])).toBe(true);
        }

        // Superset property: a leaf whose bounding volume is intersected by
        // the line must be reported.
        const partition = tree.getPartition();
        const hit = new Set<number>();
        for (const index of nodeIndices) {
            hit.add(partition[nodes[index].minIndex]);
        }
        for (let i = 0; i < nodes.length; ++i) {
            const node = nodes[i];
            if (!isLeaf(node) || node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            if (aabbOps.intersectLine(P, Q, node.boundingVolume)) {
                expect(hit.has(partition[node.minIndex])).toBe(true);
            }
        }

        // The two corners on the line, (0,0,0) and (1,2,4), are reported.
        expect(hit.has(0)).toBe(true);
        expect(hit.has(7)).toBe(true);
    });

    it('a ray query rejects leaves behind the ray origin', () => {
        const tree = new AABBTreeOfPoints();
        const points = cubeCorners();
        tree.create(points);

        // A ray that starts beyond the cube and points away from it.
        const P = Vector.fromArray([10, 10, 10]);
        const Q = Vector.fromArray([1, 1, 1]);
        expect(tree.execute(BVTree.RAY_QUERY, P, Q).length).toBe(0);

        // The reversed ray hits the cube.
        const Q2 = Vector.fromArray([-1, -1, -1]);
        expect(tree.execute(BVTree.RAY_QUERY, P, Q2).length).toBeGreaterThan(0);

        // A segment that stops short of the cube reports nothing.
        expect(tree.execute(BVTree.SEGMENT_QUERY, P,
            Vector.fromArray([5, 5, 5])).length).toBe(0);
        // A segment that reaches the cube reports leaves.
        expect(tree.execute(BVTree.SEGMENT_QUERY, P,
            Vector.fromArray([0, 0, 0])).length).toBeGreaterThan(0);
    });

    it('handles a single point (height 0)', () => {
        const tree = new AABBTreeOfPoints();
        const point = Vector.fromArray([3, -1, 2]);
        tree.create([point]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(isLeaf(root)).toBe(true);
        expect(root.minIndex).toBe(0);
        expect(root.maxIndex).toBe(0);
        expect(root.boundingVolume.min).toEqual([3, -1, 2]);
        expect(root.boundingVolume.max).toEqual([3, -1, 2]);

        const P = Vector.fromArray([3, -1, 0]);
        const Q = Vector.fromArray([0, 0, 1]);
        expect(tree.execute(BVTree.LINE_QUERY, P, Q)).toEqual([0]);
    });

    it('honors a user-specified height', () => {
        const tree = new AABBTreeOfPoints();
        tree.create(cubeCorners(), 1);

        expect(tree.getHeight()).toBe(1);
        expect(tree.getNodes().length).toBe(3);
        const nodes = tree.getNodes();
        expect(isLeaf(nodes[0])).toBe(false);
        // The children are leaves that each represent 4 points.
        expect(isLeaf(nodes[1])).toBe(true);
        expect(isLeaf(nodes[2])).toBe(true);
        expect(nodes[1].maxIndex - nodes[1].minIndex + 1).toBe(4);
        expect(nodes[2].maxIndex - nodes[2].minIndex + 1).toBe(4);
    });

    it('throws when there are no vertices', () => {
        const tree = new AABBTreeOfPoints();
        expect(() => tree.create([])).toThrow(
            'Expecting vertices to create a bounding volume tree.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of BVTreeOfPoints.h. The
// derived class adds only the vertex copy and the "vertices are already the
// centroids" pass-through, so the properties target those plus the inherited
// Execute/GetLeafIndices contract.
// ---------------------------------------------------------------------------

const bvPointCloud = fc.array(vector(3, -8, 8), { minLength: 1, maxLength: 17 });

function bvpReachable(tree: AABBTreeOfPoints): number[] {
    const nodes = tree.getNodes();
    const list: number[] = [];
    const stack = [0];
    while (stack.length > 0) {
        const index = stack.pop()!;
        list.push(index);
        const node = nodes[index];
        if (!isLeaf(node)) {
            stack.push(node.leftChild, node.rightChild);
        }
    }
    return list;
}

describe('BVTreeOfPoints verification', () => {
    it('uses the vertices as the centroids and copies both', () => {
        check(bvPointCloud, (points) => {
            const tree = new AABBTreeOfPoints();
            tree.create(points);

            const vertices = tree.getVertices();
            const centroids = tree.getCentroids();
            expect(vertices.length).toBe(points.length);
            expect(centroids.length).toBe(points.length);
            for (let i = 0; i < points.length; ++i) {
                expect([...vertices[i].values]).toEqual([...points[i].values]);
                expect([...centroids[i].values]).toEqual([...points[i].values]);
                // Upstream copies the input into mVertices and then copies
                // again into the centroids, so neither aliases the caller's
                // array nor each other.
                expect(vertices[i]).not.toBe(points[i]);
                expect(centroids[i]).not.toBe(vertices[i]);
            }

            for (const p of points) {
                p.values[2] += 100;
            }
            for (let i = 0; i < points.length; ++i) {
                expect(vertices[i].values[2]).not.toBe(points[i].values[2]);
            }
        }, 100);
    });

    it('bounds every node by its own points and tiles the partition', () => {
        check(bvPointCloud, (points) => {
            const tree = new AABBTreeOfPoints();
            tree.create(points);
            const nodes = tree.getNodes();
            const vertices = tree.getVertices();
            const partition = tree.getPartition();

            expect([...partition].sort((a, b) => a - b))
                .toEqual([...Array(points.length).keys()]);

            const covered = new Array<number>(points.length).fill(0);
            for (const index of bvpReachable(tree)) {
                const node = nodes[index];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    expect(node.boundingVolume.contains(vertices[partition[i]]))
                        .toBe(true);
                }
                if (isLeaf(node)) {
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        ++covered[i];
                    }
                    // With the full height every leaf holds exactly one point.
                    expect(node.minIndex).toBe(node.maxIndex);
                }
            }
            expect(covered.every(c => c === 1)).toBe(true);
        }, 100);
    });

    it('execute never misses a leaf whose own volume is hit', () => {
        check(fc.tuple(bvPointCloud, vector(3, -8, 8), unitVector(3),
            vector(3, -8, 8)), (input) => {
                const tree = new AABBTreeOfPoints();
                tree.create(input[0]);
                const P = input[1];
                const D = input[2];
                const Q = input[3];
                const nodes = tree.getNodes();
                const leaves = bvpReachable(tree).filter(i => isLeaf(nodes[i]));

                const cases: Array<{
                    queryType: number; A: Vector; B: Vector;
                    hit: (bv: AABB) => boolean;
                }> = [
                    {
                        queryType: BVTree.LINE_QUERY, A: P, B: D,
                        hit: (bv) => intersectSlabs(P, D, -Infinity, Infinity, bv)
                    },
                    {
                        queryType: BVTree.RAY_QUERY, A: P, B: D,
                        hit: (bv) => intersectSlabs(P, D, 0, Infinity, bv)
                    },
                    {
                        queryType: BVTree.SEGMENT_QUERY, A: P, B: Q,
                        hit: (bv) => intersectSlabs(P, sub(Q, P), 0, 1, bv)
                    }
                ];

                for (const c of cases) {
                    const reported = tree.execute(c.queryType, c.A, c.B);
                    // Every reported index is a leaf, reported once.
                    expect(new Set(reported).size).toBe(reported.length);
                    for (const index of reported) {
                        expect(isLeaf(nodes[index])).toBe(true);
                    }
                    // No false negatives (the set is a superset: upstream
                    // never tests a leaf's own volume, issue #103).
                    for (const leaf of leaves) {
                        if (c.hit(nodes[leaf].boundingVolume)) {
                            expect(reported).toContain(leaf);
                        }
                    }
                }
            }, 60);
    });

    it('honors an explicit height and keeps multi-point leaves bounded', () => {
        check(fc.tuple(bvPointCloud, fc.integer({ min: 0, max: 4 })), (input) => {
            const points = input[0];
            const height = input[1];
            const tree = new AABBTreeOfPoints();
            tree.create(points, height);
            expect(tree.getHeight()).toBe(height);

            const nodes = tree.getNodes();
            const vertices = tree.getVertices();
            const partition = tree.getPartition();
            const covered = new Array<number>(points.length).fill(0);
            for (const index of bvpReachable(tree)) {
                const node = nodes[index];
                if (!isLeaf(node)) {
                    continue;
                }
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    ++covered[i];
                    expect(node.boundingVolume.contains(vertices[partition[i]]))
                        .toBe(true);
                }
                // A truncated leaf holds at most ceil(n / 2^height) points.
                const size = node.maxIndex - node.minIndex + 1;
                expect(size).toBeLessThanOrEqual(
                    Math.ceil(points.length / 2 ** height));
            }
            expect(covered.every(c => c === 1)).toBe(true);
        }, 100);
    });

    it('handles coincident points', () => {
        check(fc.tuple(vector(3, -8, 8), fc.integer({ min: 1, max: 12 })),
            (input) => {
                const points: Vector[] = [];
                for (let i = 0; i < input[1]; ++i) {
                    points.push(input[0].clone());
                }
                const tree = new AABBTreeOfPoints();
                tree.create(points);
                const hits = tree.execute(BVTree.LINE_QUERY, input[0],
                    Vector.unit(3, 0));
                expect(hits.length).toBeGreaterThan(0);
                for (const index of hits) {
                    expect(tree.getNodes()[index].boundingVolume.contains(input[0]))
                        .toBe(true);
                }
            });
    });
});
