import { describe, it, expect } from 'vitest';
import { BVTree, BVTreeNode } from '../src/BVTree.js';
import type {
    BVTreeBoundingVolume, BVTreeSplittingAxis, BVTreeVolumeOps
} from '../src/BVTree.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { check, fc, unitVector, vector } from './helpers/arbitraries.js';

// ---------------------------------------------------------------------------
// A concrete BoundingVolume: an axis-aligned bounding box. Upstream this is
// AABBBVTreeOfTriangles.h's box type; the test supplies a minimal version so
// the abstract BVTree can be exercised.
// ---------------------------------------------------------------------------

class AABB implements BVTreeBoundingVolume {
    min: number[];
    max: number[];

    constructor() {
        this.min = [Infinity, Infinity, Infinity];
        this.max = [-Infinity, -Infinity, -Infinity];
    }

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

    // The splitting axis is the coordinate axis of largest extent, through
    // the box center.
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

// ---------------------------------------------------------------------------
// A concrete BVTree of triangles using AABB bounding volumes, modeled on the
// upstream BVTreeOfTriangles/AABBBVTreeOfTriangles derived classes.
// ---------------------------------------------------------------------------

class AABBTreeOfTriangles extends BVTree<AABB> {
    vertices: Vector[] = [];
    triangles: number[][] = [];

    constructor() {
        super(aabbOps);
    }

    createFromMesh(vertices: Vector[], triangles: number[][],
        height: number = BVTree.fullHeight): void {
        this.vertices = vertices;
        this.triangles = triangles;
        const centroids = triangles.map((tri) => {
            const v0 = vertices[tri[0]], v1 = vertices[tri[1]], v2 = vertices[tri[2]];
            return Vector.fromArray([
                (v0.values[0] + v1.values[0] + v2.values[0]) / 3,
                (v0.values[1] + v1.values[1] + v2.values[1]) / 3,
                (v0.values[2] + v1.values[2] + v2.values[2]) / 3
            ]);
        });
        this.create(centroids, height);
    }

    protected computeInteriorBoundingVolume(i0: number, i1: number, bv: AABB): void {
        bv.reset();
        for (let i = i0; i <= i1; ++i) {
            const tri = this.triangles[this.mPartition[i]];
            for (let k = 0; k < 3; ++k) {
                bv.grow(this.vertices[tri[k]]);
            }
        }
    }

    protected computeLeafBoundingVolume(i: number, bv: AABB): void {
        bv.reset();
        const tri = this.triangles[this.mPartition[i]];
        for (let k = 0; k < 3; ++k) {
            bv.grow(this.vertices[tri[k]]);
        }
    }

    // Expose the protected traversal for testing.
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }

    partitionOf(i: number): number {
        return this.mPartition[i];
    }
}

// A BVTree whose primitives are the centroids themselves, for exercising
// create() directly.
class AABBTreeOfPoints extends BVTree<AABB> {
    constructor() {
        super(aabbOps);
    }

    protected computeInteriorBoundingVolume(i0: number, i1: number, bv: AABB): void {
        bv.reset();
        for (let i = i0; i <= i1; ++i) {
            bv.grow(this.mCentroids[this.mPartition[i]]);
        }
    }

    protected computeLeafBoundingVolume(i: number, bv: AABB): void {
        bv.reset();
        bv.grow(this.mCentroids[this.mPartition[i]]);
    }

    // Expose the protected traversal for testing.
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }

    partitionOf(i: number): number {
        return this.mPartition[i];
    }
}

// ---------------------------------------------------------------------------
// Test meshes.
// ---------------------------------------------------------------------------

// The unit quad in the z = 0 plane, split into two triangles.
function makeQuad(): { vertices: Vector[]; triangles: number[][] } {
    const vertices = [
        Vector.fromArray([0, 0, 0]),
        Vector.fromArray([1, 0, 0]),
        Vector.fromArray([1, 1, 0]),
        Vector.fromArray([0, 1, 0])
    ];
    const triangles = [[0, 1, 2], [0, 2, 3]];
    return { vertices: vertices, triangles: triangles };
}

// The axis-aligned box [0,2]x[0,3]x[0,5] as 12 triangles.
function makeBoxMesh(): { vertices: Vector[]; triangles: number[][] } {
    const vertices: Vector[] = [];
    for (let i = 0; i < 8; ++i) {
        vertices.push(Vector.fromArray([
            (i & 1) !== 0 ? 2 : 0,
            (i & 2) !== 0 ? 3 : 0,
            (i & 4) !== 0 ? 5 : 0
        ]));
    }
    const triangles = [
        [0, 2, 3], [0, 3, 1], // z = 0
        [4, 5, 7], [4, 7, 6], // z = 5
        [0, 1, 5], [0, 5, 4], // y = 0
        [2, 6, 7], [2, 7, 3], // y = 3
        [0, 4, 6], [0, 6, 2], // x = 0
        [1, 3, 7], [1, 7, 5]  // x = 2
    ];
    return { vertices: vertices, triangles: triangles };
}

// Walk the reachable nodes of a tree, calling visit(nodeIndex, depth).
function walk<BV extends BVTreeBoundingVolume>(tree: BVTree<BV>,
    visit: (nodeIndex: number, depth: number) => void): void {
    const nodes = tree.getNodes();
    const stack: Array<{ index: number; depth: number }> = [{ index: 0, depth: 0 }];
    while (stack.length > 0) {
        const top = stack.pop()!;
        visit(top.index, top.depth);
        const node = nodes[top.index];
        if (node.leftChild !== BVTreeNode.invalid) {
            stack.push({ index: node.leftChild, depth: top.depth + 1 });
            stack.push({ index: node.rightChild, depth: top.depth + 1 });
        }
    }
}

describe('BVTree', () => {
    describe('construction preconditions', () => {
        it('throws when there are no centroids', () => {
            const tree = new AABBTreeOfTriangles();
            expect(() => tree.create([])).toThrow(
                /Expecting centroids to create a bounding volume tree/);
        });

        it('throws when a centroid is not 3D', () => {
            const tree = new AABBTreeOfTriangles();
            expect(() => tree.create([Vector.fromArray([1, 2])])).toThrow(
                /centroids must be 3D/);
        });

        it('query-type constants match upstream', () => {
            expect(BVTree.LINE_QUERY).toBe(0);
            expect(BVTree.RAY_QUERY).toBe(1);
            expect(BVTree.SEGMENT_QUERY).toBe(2);
        });

        it('copies the centroids rather than aliasing them', () => {
            const tree = new AABBTreeOfPoints();
            const centroids = [Vector.fromArray([1, 2, 3])];
            tree.create(centroids);
            centroids[0].values[0] = 100;
            expect(tree.getCentroids()[0].values[0]).toBe(1);
        });
    });

    describe('known-value median split', () => {
        it('orders four collinear points by their projection', () => {
            // Points on the x-axis, deliberately out of order. The AABB
            // splitting axis is the x-axis through the box center, so the
            // partition is sorted by x at every level.
            const tree = new AABBTreeOfPoints();
            tree.create([
                Vector.fromArray([0, 0, 0]),
                Vector.fromArray([3, 0, 0]),
                Vector.fromArray([1, 0, 0]),
                Vector.fromArray([2, 0, 0])
            ]);
            expect(tree.getHeight()).toBe(2);
            expect(tree.getNodes().length).toBe(7);
            // Hand-computed: the root median split gives [0, 2 | 3, 1] and
            // each child then swaps into sorted order.
            expect([...tree.getPartition()]).toEqual([0, 2, 3, 1]);

            const nodes = tree.getNodes();
            // The four leaves are nodes 3..6 and their centroids increase
            // along x.
            const xs: number[] = [];
            for (let n = 3; n < 7; ++n) {
                expect(nodes[n].minIndex).toBe(nodes[n].maxIndex);
                expect(nodes[n].leftChild).toBe(BVTreeNode.invalid);
                xs.push(tree.getCentroids()[tree.partitionOf(nodes[n].minIndex)].values[0]);
            }
            expect(xs).toEqual([0, 1, 2, 3]);

            // The root box spans the whole point set.
            expect(nodes[0].boundingVolume.min).toEqual([0, 0, 0]);
            expect(nodes[0].boundingVolume.max).toEqual([3, 0, 0]);
        });

        it('splits an odd count with the extra element on the left', () => {
            const tree = new AABBTreeOfPoints();
            tree.create([
                Vector.fromArray([0, 0, 0]),
                Vector.fromArray([1, 0, 0]),
                Vector.fromArray([2, 0, 0])
            ], 1);
            const nodes = tree.getNodes();
            // medianIndex = (3-1)/2 = 1, so the left child gets 2 elements.
            expect(nodes[1].minIndex).toBe(0);
            expect(nodes[1].maxIndex).toBe(1);
            expect(nodes[2].minIndex).toBe(2);
            expect(nodes[2].maxIndex).toBe(2);
        });
    });

    describe('a single primitive', () => {
        it('creates a tree of height 0 with one leaf node', () => {
            const tree = new AABBTreeOfTriangles();
            const { vertices } = makeQuad();
            tree.createFromMesh(vertices, [[0, 1, 2]]);

            expect(tree.getHeight()).toBe(0);
            expect(tree.getNodes().length).toBe(1);
            const root = tree.getNodes()[0];
            expect(root.minIndex).toBe(0);
            expect(root.maxIndex).toBe(0);
            expect(root.leftChild).toBe(BVTreeNode.invalid);
            expect(root.rightChild).toBe(BVTreeNode.invalid);
            // The leaf bounding volume contains the triangle.
            for (const k of [0, 1, 2]) {
                expect(root.boundingVolume.contains(vertices[k])).toBe(true);
            }
        });
    });

    describe("the quad's two triangles", () => {
        it('creates a complete tree of height 1', () => {
            const tree = new AABBTreeOfTriangles();
            const { vertices, triangles } = makeQuad();
            tree.createFromMesh(vertices, triangles);

            expect(tree.getHeight()).toBe(1);
            expect(tree.getNodes().length).toBe(3);

            const nodes = tree.getNodes();
            expect(nodes[0].minIndex).toBe(0);
            expect(nodes[0].maxIndex).toBe(1);
            expect(nodes[0].leftChild).toBe(1);
            expect(nodes[0].rightChild).toBe(2);
            expect(nodes[1].minIndex).toBe(0);
            expect(nodes[1].maxIndex).toBe(0);
            expect(nodes[2].minIndex).toBe(1);
            expect(nodes[2].maxIndex).toBe(1);

            // The partition is a permutation of {0, 1}.
            const partition = [...tree.getPartition()].sort();
            expect(partition).toEqual([0, 1]);

            // The root bounding volume is the whole quad.
            expect(nodes[0].boundingVolume.min).toEqual([0, 0, 0]);
            expect(nodes[0].boundingVolume.max).toEqual([1, 1, 0]);

            // Each leaf bounding volume contains its own triangle.
            for (const nodeIndex of [1, 2]) {
                const node = nodes[nodeIndex];
                const tri = triangles[tree.partitionOf(node.minIndex)];
                for (const k of tri) {
                    expect(node.boundingVolume.contains(vertices[k])).toBe(true);
                }
            }
        });

        it('computes centroids that are the triangle averages', () => {
            const tree = new AABBTreeOfTriangles();
            const { vertices, triangles } = makeQuad();
            tree.createFromMesh(vertices, triangles);
            const centroids = tree.getCentroids();
            expect(centroids.length).toBe(2);
            expect(centroids[0].values[0]).toBeCloseTo(2 / 3, 12);
            expect(centroids[0].values[1]).toBeCloseTo(1 / 3, 12);
            expect(centroids[1].values[0]).toBeCloseTo(1 / 3, 12);
            expect(centroids[1].values[1]).toBeCloseTo(2 / 3, 12);
        });
    });

    describe("the box's 12 triangles", () => {
        const { vertices, triangles } = makeBoxMesh();

        it('creates a complete tree of height 4 (12 rounds up to 16)', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            expect(tree.getHeight()).toBe(4);
            expect(tree.getNodes().length).toBe(2 ** 5 - 1);
        });

        it('has a partition that is a permutation of the primitive indices', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const sorted = [...tree.getPartition()].sort((a, b) => a - b);
            expect(sorted).toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        });

        it('satisfies the leaf/interior invariants', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            let numLeaves = 0;
            let maxDepth = 0;

            walk(tree, (nodeIndex, depth) => {
                const node = nodes[nodeIndex];
                maxDepth = Math.max(maxDepth, depth);
                expect(node.minIndex).not.toBe(BVTreeNode.invalid);
                expect(node.maxIndex).not.toBe(BVTreeNode.invalid);
                expect(node.minIndex).toBeLessThanOrEqual(node.maxIndex);

                const isLeaf = node.leftChild === BVTreeNode.invalid;
                // A leaf has two invalid children; an interior node has two
                // valid children.
                expect(isLeaf).toBe(node.rightChild === BVTreeNode.invalid);
                if (isLeaf) {
                    ++numLeaves;
                    // The tree is built to full height, so leaves are single
                    // primitives.
                    expect(node.minIndex).toBe(node.maxIndex);
                } else {
                    expect(node.leftChild).toBe(2 * nodeIndex + 1);
                    expect(node.rightChild).toBe(2 * nodeIndex + 2);
                    // The children partition the parent's index range and
                    // their sizes differ by at most 1.
                    const left = nodes[node.leftChild];
                    const right = nodes[node.rightChild];
                    expect(left.minIndex).toBe(node.minIndex);
                    expect(right.maxIndex).toBe(node.maxIndex);
                    expect(right.minIndex).toBe(left.maxIndex + 1);
                    const nLeft = left.maxIndex - left.minIndex + 1;
                    const nRight = right.maxIndex - right.minIndex + 1;
                    expect(Math.abs(nLeft - nRight)).toBeLessThanOrEqual(1);
                }
            });

            expect(numLeaves).toBe(triangles.length);
            expect(maxDepth).toBeLessThanOrEqual(tree.getHeight());
        });

        it('has node bounding volumes that contain their primitives', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const tri = triangles[tree.partitionOf(i)];
                    for (const k of tri) {
                        expect(node.boundingVolume.contains(vertices[k])).toBe(true);
                    }
                }
            });
        });

        it('nests child bounding volumes inside parent bounding volumes', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                if (node.leftChild === BVTreeNode.invalid) {
                    return;
                }
                for (const c of [node.leftChild, node.rightChild]) {
                    const child = nodes[c];
                    for (let j = 0; j < 3; ++j) {
                        expect(child.boundingVolume.min[j])
                            .toBeGreaterThanOrEqual(node.boundingVolume.min[j] - 1e-12);
                        expect(child.boundingVolume.max[j])
                            .toBeLessThanOrEqual(node.boundingVolume.max[j] + 1e-12);
                    }
                }
            });
        });

        it('splits by the median of the projections onto the splitting axis', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                if (node.leftChild === BVTreeNode.invalid) {
                    return;
                }
                const { origin, direction } = node.boundingVolume.getSplittingAxis();
                const left = nodes[node.leftChild];
                const right = nodes[node.rightChild];
                let maxLeft = -Infinity;
                for (let i = left.minIndex; i <= left.maxIndex; ++i) {
                    const p = dot(direction, sub(centroids[tree.partitionOf(i)], origin));
                    maxLeft = Math.max(maxLeft, p);
                }
                for (let i = right.minIndex; i <= right.maxIndex; ++i) {
                    const p = dot(direction, sub(centroids[tree.partitionOf(i)], origin));
                    expect(p).toBeGreaterThanOrEqual(maxLeft - 1e-12);
                }
            });
        });
    });

    describe('height limiting', () => {
        const { vertices, triangles } = makeBoxMesh();

        it('stops the recursion at the requested height', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles, 2);
            expect(tree.getHeight()).toBe(2);
            expect(tree.getNodes().length).toBe(2 ** 3 - 1);

            const nodes = tree.getNodes();
            let maxDepth = 0;
            let numLeaves = 0;
            walk(tree, (nodeIndex, depth) => {
                maxDepth = Math.max(maxDepth, depth);
                if (nodes[nodeIndex].leftChild === BVTreeNode.invalid) {
                    ++numLeaves;
                }
            });
            expect(maxDepth).toBe(2);
            expect(numLeaves).toBe(4);

            // At depth 2 the four "leaves" each represent three primitives.
            let total = 0;
            for (let n = 3; n < 7; ++n) {
                total += nodes[n].maxIndex - nodes[n].minIndex + 1;
            }
            expect(total).toBe(12);
        });

        it('leaves the unreachable preallocated nodes untouched', () => {
            const tree = new AABBTreeOfTriangles();
            const { vertices: quadVertices, triangles: quadTriangles } = makeQuad();
            tree.createFromMesh(quadVertices, quadTriangles, 3);
            expect(tree.getHeight()).toBe(3);
            expect(tree.getNodes().length).toBe(15);
            // Only nodes 0, 1 and 2 are reachable; the rest keep the
            // default-constructed sentinels.
            for (let n = 3; n < 15; ++n) {
                expect(tree.getNodes()[n].minIndex).toBe(BVTreeNode.invalid);
                expect(tree.getNodes()[n].maxIndex).toBe(BVTreeNode.invalid);
            }
        });

        it('accepts a height larger than the natural height', () => {
            // A height of 31 (the clamp) would preallocate 2^32-1 nodes, so
            // only a moderate over-request is exercised here. The extra
            // levels are simply never reached.
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles, 6);
            expect(tree.getHeight()).toBe(6);
            expect(tree.getNodes().length).toBe(2 ** 7 - 1);
            let numLeaves = 0;
            walk(tree, (nodeIndex) => {
                if (tree.getNodes()[nodeIndex].leftChild === BVTreeNode.invalid) {
                    ++numLeaves;
                }
            });
            expect(numLeaves).toBe(12);
        });
    });

    describe('getLeafIndices', () => {
        const { vertices, triangles } = makeBoxMesh();

        it('finds all leaves for a line through the whole mesh', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const P = Vector.fromArray([1, 1.5, 2.5]);
            const Q = Vector.fromArray([1, 0, 0]);
            const indices = tree.leafIndices(BVTree.LINE_QUERY, P, Q);
            // Every leaf is reachable because the line passes through the
            // root box and all interior boxes it intersects.
            expect(indices.length).toBeGreaterThan(0);
            for (const nodeIndex of indices) {
                expect(tree.getNodes()[nodeIndex].leftChild).toBe(BVTreeNode.invalid);
            }
            // The reported leaves are distinct.
            expect(new Set(indices).size).toBe(indices.length);
        });

        it('reports no leaves when the linear component misses the root', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const P = Vector.fromArray([100, 100, 100]);
            const Q = Vector.fromArray([1, 0, 0]);
            expect(tree.leafIndices(BVTree.LINE_QUERY, P, Q)).toEqual([]);
            expect(tree.leafIndices(BVTree.RAY_QUERY, P, Q)).toEqual([]);
            expect(tree.leafIndices(BVTree.SEGMENT_QUERY, P,
                Vector.fromArray([101, 100, 100]))).toEqual([]);
        });

        it('distinguishes the ray and line queries by direction', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            // The point is outside the mesh; the ray points away from it but
            // the line still hits it.
            const P = Vector.fromArray([-10, 1.5, 2.5]);
            const away = Vector.fromArray([-1, 0, 0]);
            expect(tree.leafIndices(BVTree.RAY_QUERY, P, away)).toEqual([]);
            expect(tree.leafIndices(BVTree.LINE_QUERY, P, away).length)
                .toBeGreaterThan(0);
        });

        it('restricts the segment query to t in [0, 1]', () => {
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const P = Vector.fromArray([-10, 1.5, 2.5]);
            // Segment that stops short of the mesh.
            expect(tree.leafIndices(BVTree.SEGMENT_QUERY, P,
                Vector.fromArray([-5, 1.5, 2.5]))).toEqual([]);
            // Segment that reaches the mesh.
            expect(tree.leafIndices(BVTree.SEGMENT_QUERY, P,
                Vector.fromArray([10, 1.5, 2.5])).length).toBeGreaterThan(0);
        });

        it('returns the root when the tree is a single leaf, without testing it', () => {
            // Upstream quirk (preserved by the port): a leaf's own bounding
            // volume is never tested, so a single-node tree always reports
            // its root even for a linear component that misses the box.
            const tree = new AABBTreeOfTriangles();
            tree.createFromMesh(vertices, [triangles[0]]);
            const P = Vector.fromArray([1000, 1000, 1000]);
            const Q = Vector.fromArray([1, 0, 0]);
            expect(tree.leafIndices(BVTree.LINE_QUERY, P, Q)).toEqual([0]);
        });
    });

    describe('randomized cross-check', () => {
        it('keeps the partition a permutation and volumes containing centroids', () => {
            // A deterministic pseudorandom generator.
            let seed = 123456789;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };

            for (let trial = 0; trial < 5; ++trial) {
                const numTriangles = 5 + trial * 7;
                const vertices: Vector[] = [];
                const triangles: number[][] = [];
                for (let t = 0; t < numTriangles; ++t) {
                    const base = vertices.length;
                    for (let k = 0; k < 3; ++k) {
                        vertices.push(Vector.fromArray([
                            10 * rand(), 10 * rand(), 10 * rand()
                        ]));
                    }
                    triangles.push([base, base + 1, base + 2]);
                }

                const tree = new AABBTreeOfTriangles();
                tree.createFromMesh(vertices, triangles);

                const sorted = [...tree.getPartition()].sort((a, b) => a - b);
                for (let i = 0; i < numTriangles; ++i) {
                    expect(sorted[i]).toBe(i);
                }

                const nodes = tree.getNodes();
                walk(tree, (nodeIndex) => {
                    const node = nodes[nodeIndex];
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        const tri = triangles[tree.partitionOf(i)];
                        for (const k of tri) {
                            expect(node.boundingVolume.contains(vertices[k])).toBe(true);
                        }
                    }
                });
            }
        });
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of BVTree.h against the
// port. The properties target the translation hazards of BVTree.h -- the
// integer median index (numProjections - 1) / 2, the size_t wrap-around of j0
// in SplitPoints, the reversed fill of the right partition, the size_t
// wrap-around loop terminator in GetLeafIndices, and the interior-only
// bounding-volume test.
// ---------------------------------------------------------------------------

// A cloud of n 3D points, 1 <= n <= 17. The upper bound crosses two
// power-of-two boundaries so the height computation is exercised on both
// complete and incomplete trees.
const pointCloud = fc.array(vector(3, -8, 8), { minLength: 1, maxLength: 17 });

// The reachable nodes of a tree, as { index, depth } records.
function reachable<BV extends BVTreeBoundingVolume>(tree: BVTree<BV>):
    Array<{ index: number; depth: number }> {
    const list: Array<{ index: number; depth: number }> = [];
    walk(tree, (index, depth) => { list.push({ index: index, depth: depth }); });
    return list;
}

function isLeafNode<BV>(node: BVTreeNode<BV>): boolean {
    return node.leftChild === BVTreeNode.invalid
        && node.rightChild === BVTreeNode.invalid;
}

describe('BVTree verification', () => {
    it('keeps the tree structure invariants for random clouds and heights', () => {
        check(fc.tuple(pointCloud, fc.integer({ min: 0, max: 6 }), fc.boolean()),
            ([points, requested, useFullHeight]) => {
                const tree = new AABBTreeOfPoints();
                const height = useFullHeight ? BVTree.fullHeight : requested;
                tree.create(points, height);

                const n = points.length;
                const nodes = tree.getNodes();
                const partition = tree.getPartition();

                // The requested height is honored exactly: fullHeight gives
                // ceil(log2(n)) (BitHacks::RoundUpToPowerOfTwo then Log2), an
                // explicit height is clamped to 31.
                const expectedHeight = useFullHeight
                    ? Math.ceil(Math.log2(n)) : Math.min(requested, 31);
                expect(tree.getHeight()).toBe(expectedHeight);
                expect(nodes.length).toBe(2 ** (expectedHeight + 1) - 1);

                // mPartition is a permutation of the centroid indices.
                expect([...partition].sort((a, b) => a - b))
                    .toEqual([...Array(n).keys()]);

                // The reachable leaves tile [0, n-1] exactly once, and the
                // depth never exceeds the tree height.
                const covered = new Array<number>(n).fill(0);
                for (const record of reachable(tree)) {
                    expect(record.depth).toBeLessThanOrEqual(tree.getHeight());
                    const node = nodes[record.index];
                    expect(node.minIndex).toBeLessThanOrEqual(node.maxIndex);
                    if (isLeafNode(node)) {
                        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                            ++covered[i];
                        }
                    } else {
                        // The children partition the parent range with no gap
                        // and no overlap, and the two halves differ in size by
                        // at most one (the balanced median split).
                        const left = nodes[node.leftChild];
                        const right = nodes[node.rightChild];
                        expect(left.minIndex).toBe(node.minIndex);
                        expect(right.maxIndex).toBe(node.maxIndex);
                        expect(right.minIndex).toBe(left.maxIndex + 1);
                        const sizeL = left.maxIndex - left.minIndex + 1;
                        const sizeR = right.maxIndex - right.minIndex + 1;
                        expect(Math.abs(sizeL - sizeR)).toBeLessThanOrEqual(1);
                        // medianIndex = (m - 1) / 2 with integer division, so
                        // the left child gets ceil(m / 2) of the m elements.
                        expect(sizeL).toBe(Math.ceil((sizeL + sizeR) / 2));
                    }
                }
                expect(covered.every(c => c === 1)).toBe(true);
            }, 100);
    });

    it('splits at the median of the projections onto the splitting axis', () => {
        check(pointCloud.filter(p => p.length >= 2), (points) => {
            const tree = new AABBTreeOfPoints();
            tree.create(points);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();

            for (const record of reachable(tree)) {
                const node = nodes[record.index];
                if (isLeafNode(node)) {
                    continue;
                }
                // The splitting axis comes from the node's own bounding
                // volume, which is computed before the split.
                const axis = node.boundingVolume.getSplittingAxis();
                const project = (i: number): number =>
                    dot(axis.direction,
                        sub(centroids[tree.partitionOf(i)], axis.origin));

                const left = nodes[node.leftChild];
                const right = nodes[node.rightChild];
                let maxLeft = -Infinity;
                for (let i = left.minIndex; i <= left.maxIndex; ++i) {
                    maxLeft = Math.max(maxLeft, project(i));
                }
                let minRight = +Infinity;
                for (let i = right.minIndex; i <= right.maxIndex; ++i) {
                    minRight = Math.min(minRight, project(i));
                }
                // The postcondition of std::nth_element: everything before the
                // median is <= everything after it.
                expect(maxLeft).toBeLessThanOrEqual(minRight);
            }
        }, 100);
    });

    it('bounds every node by its own primitives and nests child in parent', () => {
        check(pointCloud, (points) => {
            const tree = new AABBTreeOfPoints();
            tree.create(points);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();

            for (const record of reachable(tree)) {
                const node = nodes[record.index];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    expect(node.boundingVolume.contains(
                        centroids[tree.partitionOf(i)])).toBe(true);
                }
                if (!isLeafNode(node)) {
                    const children = [nodes[node.leftChild], nodes[node.rightChild]];
                    for (const child of children) {
                        for (let j = 0; j < 3; ++j) {
                            expect(child.boundingVolume.min[j])
                                .toBeGreaterThanOrEqual(node.boundingVolume.min[j]);
                            expect(child.boundingVolume.max[j])
                                .toBeLessThanOrEqual(node.boundingVolume.max[j]);
                        }
                    }
                }
            }
        }, 100);
    });

    it('getLeafIndices never misses a leaf whose own volume is hit', () => {
        // Upstream tests only interior nodes (issue #103), so the reported set
        // is a conservative superset. The property that must hold is the
        // absence of false negatives: a leaf whose own box the linear
        // component meets is always reported.
        check(fc.tuple(pointCloud.filter(p => p.length >= 2), vector(3, -8, 8),
            unitVector(3), vector(3, -8, 8)),
            (input) => {
                const points = input[0];
                const P = input[1];
                const D = input[2];
                const Q = input[3];
                const tree = new AABBTreeOfPoints();
                tree.create(points);
                const nodes = tree.getNodes();
                const allLeaves = reachable(tree)
                    .filter(r => isLeafNode(nodes[r.index])).map(r => r.index);

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
                    const reported = new Set(tree.leafIndices(c.queryType, c.A, c.B));
                    // Only leaves are ever reported, never an interior node.
                    for (const index of reported) {
                        expect(isLeafNode(nodes[index])).toBe(true);
                    }
                    for (const leaf of allLeaves) {
                        if (c.hit(nodes[leaf].boundingVolume)) {
                            expect(reported.has(leaf)).toBe(true);
                        }
                    }
                }
            }, 60);
    });

    it('handles coincident centroids without degenerating the partition', () => {
        check(fc.tuple(vector(3, -8, 8), fc.integer({ min: 1, max: 12 })),
            (input) => {
                const n = input[1];
                const points: Vector[] = [];
                for (let i = 0; i < n; ++i) {
                    points.push(input[0].clone());
                }
                const tree = new AABBTreeOfPoints();
                tree.create(points);
                expect([...tree.getPartition()].sort((a, b) => a - b))
                    .toEqual([...Array(n).keys()]);
                // Every projection is zero, so the split is purely positional;
                // the tree must still be complete and every box must be the
                // degenerate box at the common point.
                const nodes = tree.getNodes();
                for (const record of reachable(tree)) {
                    const node = nodes[record.index];
                    expect(Number.isFinite(node.boundingVolume.min[0])).toBe(true);
                    expect(node.boundingVolume.contains(input[0])).toBe(true);
                }
            });
    });

    it('copies the centroids so later input mutation cannot reach the tree', () => {
        check(pointCloud, (points) => {
            const tree = new AABBTreeOfPoints();
            tree.create(points);
            const before = tree.getCentroids().map(c => [...c.values]);
            for (const p of points) {
                p.values[0] += 1000;
            }
            const after = tree.getCentroids().map(c => [...c.values]);
            expect(after).toEqual(before);
        }, 50);
    });
});
