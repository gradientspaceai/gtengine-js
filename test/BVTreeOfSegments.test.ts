import { describe, it, expect } from 'vitest';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from '../src/BVTree.js';
import { BVTreeOfSegments } from '../src/BVTreeOfSegments.js';
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

// The concrete segment tree: the bounding volume of a node contains the
// endpoints of the segments the node represents.
class AABBTreeOfSegments extends BVTreeOfSegments<AABB> {
    constructor() {
        super(aabbOps);
    }

    private growSegment(bv: AABB, segmentIndex: number): void {
        const seg = this.mSegments[segmentIndex];
        bv.grow(this.mVertices[seg[0]]);
        bv.grow(this.mVertices[seg[1]]);
    }

    protected computeInteriorBoundingVolume(i0: number, i1: number, bv: AABB): void {
        bv.reset();
        for (let i = i0; i <= i1; ++i) {
            this.growSegment(bv, this.mPartition[i]);
        }
    }

    protected computeLeafBoundingVolume(i: number, bv: AABB): void {
        bv.reset();
        this.growSegment(bv, this.mPartition[i]);
    }
}

function isLeaf(node: BVTreeNode<AABB>): boolean {
    return node.leftChild === BVTreeNode.invalid
        && node.rightChild === BVTreeNode.invalid;
}

// A polyline of n+1 vertices along the x axis, with n segments.
function makePolyline(n: number): {
    vertices: Vector[];
    segments: [number, number][];
} {
    const vertices: Vector[] = [];
    for (let i = 0; i <= n; ++i) {
        vertices.push(Vector.fromArray([i, (i % 2 === 0) ? 0 : 1, 0]));
    }
    const segments: [number, number][] = [];
    for (let i = 0; i < n; ++i) {
        segments.push([i, i + 1]);
    }
    return { vertices: vertices, segments: segments };
}

describe('BVTreeOfSegments', () => {
    it('computes segment midpoints as the centroids', () => {
        const tree = new AABBTreeOfSegments();
        const { vertices, segments } = makePolyline(4);
        tree.createFromSegments(vertices, segments);

        const centroids = tree.getCentroids();
        expect(centroids.length).toBe(4);
        for (let i = 0; i < segments.length; ++i) {
            const seg = segments[i];
            for (let j = 0; j < 3; ++j) {
                expect(centroids[i].values[j]).toBeCloseTo(
                    0.5 * (vertices[seg[0]].values[j] + vertices[seg[1]].values[j]), 15);
            }
        }
    });

    it('builds a balanced tree and copies its inputs', () => {
        const tree = new AABBTreeOfSegments();
        const { vertices, segments } = makePolyline(4);
        tree.createFromSegments(vertices, segments);

        // 4 = 2^2 leaves, so the height is 2 and there are 2^3-1 = 7 nodes.
        expect(tree.getHeight()).toBe(2);
        expect(tree.getNodes().length).toBe(7);
        expect(tree.getNodes()[0].minIndex).toBe(0);
        expect(tree.getNodes()[0].maxIndex).toBe(3);
        expect(tree.getSegments()).toEqual([[0, 1], [1, 2], [2, 3], [3, 4]]);
        expect(tree.getVertices().length).toBe(5);

        // Mutating the inputs afterward does not change the tree.
        vertices[0].values[0] = 100;
        segments[0][1] = 4;
        expect(tree.getVertices()[0].values[0]).toBe(0);
        expect(tree.getSegments()[0][1]).toBe(1);
    });

    it('each node bounding volume contains its segment endpoints', () => {
        const tree = new AABBTreeOfSegments();
        let seed = 987654321;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const vertices: Vector[] = [];
        for (let i = 0; i < 30; ++i) {
            vertices.push(Vector.fromArray([4 * rand(), 2 * rand(), rand() - 0.5]));
        }
        const segments: [number, number][] = [];
        for (let i = 0; i < 29; ++i) {
            segments.push([i, i + 1]);
        }
        tree.createFromSegments(vertices, segments);

        const nodes = tree.getNodes();
        const partition = tree.getPartition();
        expect([...partition].sort((a, b) => a - b)).toEqual(
            segments.map((_s, i) => i));

        for (const node of nodes) {
            if (node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const seg = segments[partition[i]];
                expect(node.boundingVolume.contains(vertices[seg[0]])).toBe(true);
                expect(node.boundingVolume.contains(vertices[seg[1]])).toBe(true);
            }
        }
    });

    it('execute reports the leaves whose bounding volumes are hit', () => {
        const tree = new AABBTreeOfSegments();
        const { vertices, segments } = makePolyline(8);
        tree.createFromSegments(vertices, segments);

        // A line along the polyline's y = 0.5 midline in the x direction.
        const P = Vector.fromArray([-1, 0.5, 0]);
        const Q = Vector.fromArray([1, 0, 0]);
        const nodeIndices = tree.execute(BVTree.LINE_QUERY, P, Q);
        expect(nodeIndices.length).toBe(8);

        const nodes = tree.getNodes();
        for (const index of nodeIndices) {
            expect(isLeaf(nodes[index])).toBe(true);
        }

        // A line far from the polyline reports no leaves.
        const faraway = tree.execute(BVTree.LINE_QUERY,
            Vector.fromArray([-1, 50, 0]), Q);
        expect(faraway.length).toBe(0);

        // Superset property for a line that hits only part of the polyline.
        const P2 = Vector.fromArray([2.5, -1, 0]);
        const Q2 = Vector.fromArray([0, 1, 0]);
        const hitIndices = tree.execute(BVTree.LINE_QUERY, P2, Q2);
        const reported = new Set(hitIndices);
        for (let i = 0; i < nodes.length; ++i) {
            const node = nodes[i];
            if (!isLeaf(node) || node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            if (aabbOps.intersectLine(P2, Q2, node.boundingVolume)) {
                expect(reported.has(i)).toBe(true);
            }
        }
        expect(hitIndices.length).toBeLessThan(8);
    });

    it('handles a single segment (height 0)', () => {
        const tree = new AABBTreeOfSegments();
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([2, 0, 0])
        ];
        tree.createFromSegments(vertices, [[0, 1]]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(isLeaf(root)).toBe(true);
        expect(root.boundingVolume.min).toEqual([0, 0, 0]);
        expect(root.boundingVolume.max).toEqual([2, 0, 0]);
        expect(tree.getCentroids()[0].values).toEqual([1, 0, 0]);

        expect(tree.execute(BVTree.SEGMENT_QUERY, Vector.fromArray([1, -1, 0]),
            Vector.fromArray([1, 1, 0]))).toEqual([0]);
    });

    it('honors a user-specified height', () => {
        const tree = new AABBTreeOfSegments();
        const { vertices, segments } = makePolyline(8);
        tree.createFromSegments(vertices, segments, 1);

        expect(tree.getHeight()).toBe(1);
        expect(tree.getNodes().length).toBe(3);
        const nodes = tree.getNodes();
        expect(isLeaf(nodes[1])).toBe(true);
        expect(isLeaf(nodes[2])).toBe(true);
        expect(nodes[1].maxIndex - nodes[1].minIndex + 1).toBe(4);
        expect(nodes[2].maxIndex - nodes[2].minIndex + 1).toBe(4);
    });

    it('throws when there are no vertices', () => {
        const tree = new AABBTreeOfSegments();
        expect(() => tree.createFromSegments([], [])).toThrow(
            'Expecting vertices to create a bounding volume tree.');
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of BVTreeOfSegments.h. The
// derived class adds the vertex/segment copies and the midpoint centroids;
// the rest is the inherited BVTree contract.
// ---------------------------------------------------------------------------

// A vertex pool with index-pair segments. Indices may repeat, so degenerate
// (zero-length) segments occur, which is what upstream accepts silently.
const segmentSoup = fc.tuple(
    fc.array(vector(3, -8, 8), { minLength: 2, maxLength: 10 }),
    fc.array(fc.tuple(fc.nat({ max: 9 }), fc.nat({ max: 9 })),
        { minLength: 1, maxLength: 13 }))
    .map((input) => ({
        vertices: input[0],
        segments: input[1].map(p => [p[0] % input[0].length,
            p[1] % input[0].length] as [number, number])
    }));

function bvsReachable(tree: AABBTreeOfSegments): number[] {
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

describe('BVTreeOfSegments verification', () => {
    it('computes the exact segment midpoints as the centroids', () => {
        check(segmentSoup, (mesh) => {
            const tree = new AABBTreeOfSegments();
            tree.createFromSegments(mesh.vertices, mesh.segments);

            const centroids = tree.getCentroids();
            expect(centroids.length).toBe(mesh.segments.length);
            for (let i = 0; i < mesh.segments.length; ++i) {
                const seg = mesh.segments[i];
                for (let j = 0; j < 3; ++j) {
                    // Upstream: half * (V[seg[0]] + V[seg[1]]). The port must
                    // use exactly that expression, so this is bit-exact.
                    expect(centroids[i].values[j]).toBe(
                        0.5 * (mesh.vertices[seg[0]].values[j]
                            + mesh.vertices[seg[1]].values[j]));
                }
            }

            // Both inputs are copied, so later mutation cannot reach the tree.
            const vertices = tree.getVertices();
            const segments = tree.getSegments();
            expect(segments.map(s => [s[0], s[1]])).toEqual(
                mesh.segments.map(s => [s[0], s[1]]));
            for (let i = 0; i < mesh.vertices.length; ++i) {
                expect(vertices[i]).not.toBe(mesh.vertices[i]);
                mesh.vertices[i].values[0] += 1000;
                expect(vertices[i].values[0])
                    .not.toBe(mesh.vertices[i].values[0]);
            }
            mesh.segments[0][0] = 0;
        }, 100);
    });

    it('bounds every node by the endpoints of its segments', () => {
        check(segmentSoup, (mesh) => {
            const tree = new AABBTreeOfSegments();
            tree.createFromSegments(mesh.vertices, mesh.segments);
            const nodes = tree.getNodes();
            const partition = tree.getPartition();

            expect([...partition].sort((a, b) => a - b))
                .toEqual([...Array(mesh.segments.length).keys()]);

            const covered = new Array<number>(mesh.segments.length).fill(0);
            for (const index of bvsReachable(tree)) {
                const node = nodes[index];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const seg = mesh.segments[partition[i]];
                    expect(node.boundingVolume.contains(mesh.vertices[seg[0]]))
                        .toBe(true);
                    expect(node.boundingVolume.contains(mesh.vertices[seg[1]]))
                        .toBe(true);
                }
                if (isLeaf(node)) {
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        ++covered[i];
                    }
                }
            }
            expect(covered.every(c => c === 1)).toBe(true);
        }, 100);
    });

    it('execute never misses a leaf whose own volume is hit', () => {
        check(fc.tuple(segmentSoup, vector(3, -8, 8), unitVector(3),
            vector(3, -8, 8)), (input) => {
                const tree = new AABBTreeOfSegments();
                tree.createFromSegments(input[0].vertices, input[0].segments);
                const P = input[1];
                const D = input[2];
                const Q = input[3];
                const nodes = tree.getNodes();
                const leaves = bvsReachable(tree).filter(i => isLeaf(nodes[i]));

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
                    expect(new Set(reported).size).toBe(reported.length);
                    for (const index of reported) {
                        expect(isLeaf(nodes[index])).toBe(true);
                    }
                    for (const leaf of leaves) {
                        if (c.hit(nodes[leaf].boundingVolume)) {
                            expect(reported).toContain(leaf);
                        }
                    }
                }
            }, 60);
    });

    it('accepts degenerate zero-length segments', () => {
        check(fc.tuple(vector(3, -8, 8), fc.integer({ min: 1, max: 8 })),
            (input) => {
                // Every segment collapses to the same vertex.
                const vertices = [input[0].clone(), input[0].clone()];
                const segments: [number, number][] = [];
                for (let i = 0; i < input[1]; ++i) {
                    segments.push([0, 1]);
                }
                const tree = new AABBTreeOfSegments();
                tree.createFromSegments(vertices, segments);
                for (const c of tree.getCentroids()) {
                    for (let j = 0; j < 3; ++j) {
                        expect(c.values[j]).toBe(input[0].values[j]);
                    }
                }
                const hits = tree.execute(BVTree.LINE_QUERY, input[0],
                    Vector.unit(3, 1));
                expect(hits.length).toBeGreaterThan(0);
            });
    });
});
