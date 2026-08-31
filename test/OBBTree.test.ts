import { describe, it, expect } from 'vitest';
import { OBBNode, OBBTree } from '../src/OBBTree';
import { OrientedBox } from '../src/OrientedBox';
import { Vector, dot, normalize, sub } from '../src/Vector';
import { cross, dotCross } from '../src/Vector3';

// ---------------------------------------------------------------------------
// A concrete OBBTree of triangles, modeled on upstream OBBTreeOfTriangles.h.
// The only deviation is that the interior-box projection extremes start at
// +/-infinity rather than at zero (the upstream initialization to zero does
// not guarantee containment).
// ---------------------------------------------------------------------------

class OBBTreeOfTriangles extends OBBTree {
    vertices: Vector[] = [];
    triangles: number[][] = [];

    constructor() {
        super();
    }

    createFromMesh(vertices: Vector[], triangles: number[][],
        height: number = OBBTree.fullHeight): void {
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

    partitionOf(i: number): number {
        return this.mPartition[i];
    }

    protected override computeInteriorBox(i0: number, i1: number, box: OrientedBox): void {
        super.computeInteriorBox(i0, i1, box);

        const pmin = [Infinity, Infinity, Infinity];
        const pmax = [-Infinity, -Infinity, -Infinity];
        for (let i = i0; i <= i1; ++i) {
            const tri = this.triangles[this.mPartition[i]];
            for (let k = 0; k < 3; ++k) {
                const diff = sub(this.vertices[tri[k]], box.center);
                for (let j = 0; j < 3; ++j) {
                    const d = dot(diff, box.axis[j]);
                    if (d < pmin[j]) {
                        pmin[j] = d;
                    }
                    if (d > pmax[j]) {
                        pmax[j] = d;
                    }
                }
            }
        }

        const newCenter = box.center.clone();
        for (let j = 0; j < 3; ++j) {
            const s = 0.5 * (pmin[j] + pmax[j]);
            for (let c = 0; c < 3; ++c) {
                newCenter.values[c] += s * box.axis[j].values[c];
            }
            box.extent.values[j] = 0.5 * (pmax[j] - pmin[j]);
        }
        box.center = newCenter;
    }

    protected computeLeafBox(i: number, box: OrientedBox): void {
        const tri = this.triangles[this.mPartition[i]];
        const edge10 = sub(this.vertices[tri[1]], this.vertices[tri[0]]);
        const edge20 = sub(this.vertices[tri[2]], this.vertices[tri[0]]);
        normalize(edge10);
        normalize(edge20);
        const normal = cross(edge10, edge20);
        normalize(normal);

        box.center = this.mCentroids[this.mPartition[i]].clone();
        box.axis[0] = edge10;
        box.axis[1] = cross(normal, edge10);
        box.axis[2] = normal;

        const V0mC = sub(this.vertices[tri[0]], box.center);
        const V1mC = sub(this.vertices[tri[1]], box.center);
        const V2mC = sub(this.vertices[tri[2]], box.center);
        box.extent.values[0] = Math.max(
            Math.abs(dot(box.axis[0], V0mC)),
            Math.abs(dot(box.axis[0], V1mC)),
            Math.abs(dot(box.axis[0], V2mC)));
        box.extent.values[1] = Math.max(
            Math.abs(dot(box.axis[1], V0mC)),
            Math.abs(dot(box.axis[1], V1mC)),
            Math.abs(dot(box.axis[1], V2mC)));
        box.extent.values[2] = 0;
    }
}

// A minimal derived class that keeps the base-class interior box (eigenvalues
// stored in the extents) so the base ComputeInteriorBox can be tested.
class OBBTreeOfPointsRaw extends OBBTree {
    constructor() {
        super();
    }

    protected computeLeafBox(i: number, box: OrientedBox): void {
        box.center = this.mCentroids[this.mPartition[i]].clone();
        box.extent = new Vector(3);
    }
}

// Is the point inside or on the box (within tolerance)?
function boxContains(box: OrientedBox, p: Vector, eps: number = 1e-10): boolean {
    const diff = sub(p, box.center);
    for (let j = 0; j < 3; ++j) {
        if (Math.abs(dot(diff, box.axis[j])) > box.extent.values[j] + eps) {
            return false;
        }
    }
    return true;
}

function makeQuad(): { vertices: Vector[]; triangles: number[][] } {
    const vertices = [
        Vector.fromArray([0, 0, 0]),
        Vector.fromArray([1, 0, 0]),
        Vector.fromArray([1, 1, 0]),
        Vector.fromArray([0, 1, 0])
    ];
    return { vertices: vertices, triangles: [[0, 1, 2], [0, 2, 3]] };
}

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
        [0, 2, 3], [0, 3, 1],
        [4, 5, 7], [4, 7, 6],
        [0, 1, 5], [0, 5, 4],
        [2, 6, 7], [2, 7, 3],
        [0, 4, 6], [0, 6, 2],
        [1, 3, 7], [1, 7, 5]
    ];
    return { vertices: vertices, triangles: triangles };
}

function walk(tree: OBBTree, visit: (nodeIndex: number, depth: number) => void): void {
    const nodes = tree.getNodes();
    const stack: Array<{ index: number; depth: number }> = [{ index: 0, depth: 0 }];
    while (stack.length > 0) {
        const top = stack.pop()!;
        visit(top.index, top.depth);
        const node = nodes[top.index];
        if (node.leftChild !== OBBNode.invalid) {
            stack.push({ index: node.leftChild, depth: top.depth + 1 });
            stack.push({ index: node.rightChild, depth: top.depth + 1 });
        }
    }
}

describe('OBBNode', () => {
    it('default-constructs with invalid indices and a unit box', () => {
        const node = new OBBNode();
        expect(node.minIndex).toBe(OBBNode.invalid);
        expect(node.maxIndex).toBe(OBBNode.invalid);
        expect(node.leftChild).toBe(OBBNode.invalid);
        expect(node.rightChild).toBe(OBBNode.invalid);
        expect(node.box.dimension).toBe(3);
    });
});

describe('OBBTree', () => {
    describe('construction preconditions', () => {
        it('throws for an empty centroid list', () => {
            const tree = new OBBTreeOfPointsRaw();
            expect(() => tree.create([])).toThrow(/Invalid input/);
        });

        it('throws for non-3D centroids', () => {
            const tree = new OBBTreeOfPointsRaw();
            expect(() => tree.create([Vector.fromArray([1, 2])]))
                .toThrow(/centroids must be 3D/);
        });

        it('copies the centroids rather than aliasing them', () => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids = [Vector.fromArray([1, 2, 3])];
            tree.create(centroids);
            centroids[0].values[1] = -7;
            expect(tree.getCentroids()[0].values[1]).toBe(2);
        });
    });

    describe('height from the number of centroids', () => {
        it.each([
            [1, 0], [2, 1], [3, 2], [4, 2], [5, 3], [8, 3], [9, 4], [12, 4], [16, 4]
        ])('%i centroids gives height %i', (count, expected) => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids: Vector[] = [];
            for (let i = 0; i < count; ++i) {
                centroids.push(Vector.fromArray([i, i * i % 7, (3 * i) % 5]));
            }
            tree.create(centroids);
            expect(tree.getHeight()).toBe(expected);
            expect(tree.getNodes().length).toBe(2 ** (expected + 1) - 1);
        });
    });

    describe('base computeInteriorBox', () => {
        it('centers the box at the mean of the centroids', () => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids = [
                Vector.fromArray([0, 0, 0]),
                Vector.fromArray([4, 0, 0]),
                Vector.fromArray([0, 2, 0]),
                Vector.fromArray([4, 2, 0])
            ];
            tree.create(centroids, 0);
            const root = tree.getNodes()[0];
            expect(root.box.center.values[0]).toBeCloseTo(2, 12);
            expect(root.box.center.values[1]).toBeCloseTo(1, 12);
            expect(root.box.center.values[2]).toBeCloseTo(0, 12);
        });

        it('stores the covariance eigenvalues in the extents, increasing', () => {
            const tree = new OBBTreeOfPointsRaw();
            // Spread 10 along x, 2 along y, 0 along z.
            const centroids = [
                Vector.fromArray([-5, -1, 0]),
                Vector.fromArray([5, -1, 0]),
                Vector.fromArray([-5, 1, 0]),
                Vector.fromArray([5, 1, 0])
            ];
            tree.create(centroids, 0);
            const box = tree.getNodes()[0].box;
            // Covariance is diag(25, 1, 0); eigenvalues sorted increasingly.
            expect(box.extent.values[0]).toBeCloseTo(0, 12);
            expect(box.extent.values[1]).toBeCloseTo(1, 12);
            expect(box.extent.values[2]).toBeCloseTo(25, 12);
        });

        it('recovers the coordinate axes for axis-aligned data', () => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids = [
                Vector.fromArray([-5, -1, -0.5]),
                Vector.fromArray([5, -1, -0.5]),
                Vector.fromArray([-5, 1, -0.5]),
                Vector.fromArray([5, 1, -0.5]),
                Vector.fromArray([-5, -1, 0.5]),
                Vector.fromArray([5, -1, 0.5]),
                Vector.fromArray([-5, 1, 0.5]),
                Vector.fromArray([5, 1, 0.5])
            ];
            tree.create(centroids, 0);
            const box = tree.getNodes()[0].box;
            // Increasing eigenvalues: z (0.25), y (1), x (25). Each axis is a
            // coordinate axis up to sign.
            const expectedAxes = [2, 1, 0];
            for (let j = 0; j < 3; ++j) {
                const axis = box.axis[j];
                for (let c = 0; c < 3; ++c) {
                    const target = c === expectedAxes[j] ? 1 : 0;
                    expect(Math.abs(axis.values[c])).toBeCloseTo(target, 10);
                }
            }
        });

        it('produces a right-handed orthonormal axis frame', () => {
            const tree = new OBBTreeOfPointsRaw();
            let seed = 987654321;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648 - 0.5;
            };
            const centroids: Vector[] = [];
            for (let i = 0; i < 40; ++i) {
                centroids.push(Vector.fromArray([8 * rand(), 3 * rand(), rand()]));
            }
            tree.create(centroids, 0);
            const box = tree.getNodes()[0].box;
            for (let j = 0; j < 3; ++j) {
                expect(dot(box.axis[j], box.axis[j])).toBeCloseTo(1, 10);
                for (let k = j + 1; k < 3; ++k) {
                    expect(dot(box.axis[j], box.axis[k])).toBeCloseTo(0, 10);
                }
            }
            expect(dotCross(box.axis[0], box.axis[1], box.axis[2]))
                .toBeCloseTo(1, 10);
            // The eigenvalues (stored in the extents) are increasing.
            expect(box.extent.values[0]).toBeLessThanOrEqual(box.extent.values[1]);
            expect(box.extent.values[1]).toBeLessThanOrEqual(box.extent.values[2]);
        });
    });

    describe("the quad's two triangles", () => {
        it('creates a complete tree of height 1 with two leaves', () => {
            const tree = new OBBTreeOfTriangles();
            const { vertices, triangles } = makeQuad();
            tree.createFromMesh(vertices, triangles);

            expect(tree.getHeight()).toBe(1);
            expect(tree.getNodes().length).toBe(3);
            const nodes = tree.getNodes();
            expect(nodes[0].leftChild).toBe(1);
            expect(nodes[0].rightChild).toBe(2);
            expect(nodes[1].leftChild).toBe(OBBNode.invalid);
            expect(nodes[2].leftChild).toBe(OBBNode.invalid);
            expect(nodes[1].minIndex).toBe(0);
            expect(nodes[1].maxIndex).toBe(0);
            expect(nodes[2].minIndex).toBe(1);
            expect(nodes[2].maxIndex).toBe(1);
            expect([...tree.getPartition()].sort()).toEqual([0, 1]);
        });

        it('contains the triangles in the node boxes', () => {
            const tree = new OBBTreeOfTriangles();
            const { vertices, triangles } = makeQuad();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const tri = triangles[tree.partitionOf(i)];
                    for (const k of tri) {
                        expect(boxContains(node.box, vertices[k], 1e-9)).toBe(true);
                    }
                }
            });
        });

        it('gives the root a planar box of half-extent 0 in the normal direction', () => {
            const tree = new OBBTreeOfTriangles();
            const { vertices, triangles } = makeQuad();
            tree.createFromMesh(vertices, triangles);
            const box = tree.getNodes()[0].box;
            // All data lies in z = 0, so one extent is zero and the box is
            // centered on the quad center.
            const zeroExtents = box.extent.values.filter(
                (e) => Math.abs(e) < 1e-10).length;
            expect(zeroExtents).toBe(1);
            expect(box.center.values[2]).toBeCloseTo(0, 10);
            expect(box.center.values[0]).toBeCloseTo(0.5, 10);
            expect(box.center.values[1]).toBeCloseTo(0.5, 10);
        });
    });

    describe("the box's 12 triangles", () => {
        const { vertices, triangles } = makeBoxMesh();

        it('creates a complete tree of height 4', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            expect(tree.getHeight()).toBe(4);
            expect(tree.getNodes().length).toBe(31);
        });

        it('has a partition that is a permutation of the triangle indices', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            expect([...tree.getPartition()].sort((a, b) => a - b))
                .toEqual([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]);
        });

        it('satisfies the leaf/interior partition invariants', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            let numLeaves = 0;
            let maxDepth = 0;
            walk(tree, (nodeIndex, depth) => {
                const node = nodes[nodeIndex];
                maxDepth = Math.max(maxDepth, depth);
                const isLeaf = node.leftChild === OBBNode.invalid;
                expect(isLeaf).toBe(node.rightChild === OBBNode.invalid);
                if (isLeaf) {
                    ++numLeaves;
                    expect(node.minIndex).toBe(node.maxIndex);
                } else {
                    expect(node.leftChild).toBe(2 * nodeIndex + 1);
                    expect(node.rightChild).toBe(2 * nodeIndex + 2);
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
            expect(numLeaves).toBe(12);
            expect(maxDepth).toBeLessThanOrEqual(4);
        });

        it('contains every triangle of a node in the node box', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const tri = triangles[tree.partitionOf(i)];
                    for (const k of tri) {
                        expect(boxContains(node.box, vertices[k], 1e-9)).toBe(true);
                    }
                }
            });
        });

        it('splits about the median projection onto the largest-eigenvalue axis', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();
            walk(tree, (nodeIndex) => {
                const node = nodes[nodeIndex];
                if (node.leftChild === OBBNode.invalid) {
                    return;
                }
                // The split used the box state at the time of the split, so
                // recompute the ordering only qualitatively: the two child
                // index ranges are contiguous and cover the parent range.
                const left = nodes[node.leftChild];
                const right = nodes[node.rightChild];
                expect(left.maxIndex + 1).toBe(right.minIndex);
                const nLeft = left.maxIndex - left.minIndex + 1;
                const total = node.maxIndex - node.minIndex + 1;
                expect(nLeft).toBe(Math.floor((total - 1) / 2) + 1);
                expect(centroids.length).toBe(12);
            });
        });

        it('stops the recursion at a requested height', () => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(vertices, triangles, 2);
            expect(tree.getHeight()).toBe(2);
            expect(tree.getNodes().length).toBe(7);
            let maxDepth = 0;
            walk(tree, (_nodeIndex, depth) => {
                maxDepth = Math.max(maxDepth, depth);
            });
            expect(maxDepth).toBe(2);
            // Every primitive is still represented exactly once.
            let total = 0;
            for (let n = 3; n < 7; ++n) {
                const node = tree.getNodes()[n];
                total += node.maxIndex - node.minIndex + 1;
            }
            expect(total).toBe(12);
        });

        it('leaves unreachable preallocated nodes at the sentinels', () => {
            const tree = new OBBTreeOfTriangles();
            const quad = makeQuad();
            tree.createFromMesh(quad.vertices, quad.triangles, 3);
            expect(tree.getNodes().length).toBe(15);
            for (let n = 3; n < 15; ++n) {
                expect(tree.getNodes()[n].minIndex).toBe(OBBNode.invalid);
                expect(tree.getNodes()[n].leftChild).toBe(OBBNode.invalid);
            }
        });
    });

    describe('degenerate inputs', () => {
        it('handles all centroids being identical', () => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids: Vector[] = [];
            for (let i = 0; i < 4; ++i) {
                centroids.push(Vector.fromArray([1, 2, 3]));
            }
            tree.create(centroids);
            // The covariance is the zero matrix; the eigenvalues are zero and
            // the axes form an orthonormal frame.
            const box = tree.getNodes()[0].box;
            expect(box.center.values).toEqual([1, 2, 3]);
            for (let j = 0; j < 3; ++j) {
                expect(box.extent.values[j]).toBeCloseTo(0, 12);
                expect(dot(box.axis[j], box.axis[j])).toBeCloseTo(1, 12);
            }
            expect([...tree.getPartition()].sort((a, b) => a - b))
                .toEqual([0, 1, 2, 3]);
        });

        it('handles collinear centroids', () => {
            const tree = new OBBTreeOfPointsRaw();
            const centroids: Vector[] = [];
            for (let i = 0; i < 5; ++i) {
                centroids.push(Vector.fromArray([i, 2 * i, 3 * i]));
            }
            tree.create(centroids);
            const box = tree.getNodes()[0].box;
            // Two eigenvalues are zero; the largest-eigenvalue axis is the
            // line direction (up to sign).
            expect(box.extent.values[0]).toBeCloseTo(0, 10);
            expect(box.extent.values[1]).toBeCloseTo(0, 10);
            expect(box.extent.values[2]).toBeGreaterThan(1);
            const d = Vector.fromArray([1, 2, 3]);
            normalize(d);
            expect(Math.abs(dot(box.axis[2], d))).toBeCloseTo(1, 10);
        });

        it('builds a single-node tree for one centroid', () => {
            const tree = new OBBTreeOfPointsRaw();
            tree.create([Vector.fromArray([7, 8, 9])]);
            expect(tree.getHeight()).toBe(0);
            expect(tree.getNodes().length).toBe(1);
            const root = tree.getNodes()[0];
            expect(root.minIndex).toBe(0);
            expect(root.maxIndex).toBe(0);
            expect(root.leftChild).toBe(OBBNode.invalid);
            expect(root.box.center.values).toEqual([7, 8, 9]);
        });
    });

    describe('randomized cross-check', () => {
        it('keeps the partition a permutation and boxes containing the triangles', () => {
            let seed = 24680;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };

            for (let trial = 0; trial < 4; ++trial) {
                const numTriangles = 6 + 5 * trial;
                const vertices: Vector[] = [];
                const triangles: number[][] = [];
                for (let t = 0; t < numTriangles; ++t) {
                    const base = vertices.length;
                    const c = [10 * rand(), 10 * rand(), 10 * rand()];
                    for (let k = 0; k < 3; ++k) {
                        vertices.push(Vector.fromArray([
                            c[0] + rand(), c[1] + rand(), c[2] + rand()
                        ]));
                    }
                    triangles.push([base, base + 1, base + 2]);
                }

                const tree = new OBBTreeOfTriangles();
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
                            expect(boxContains(node.box, vertices[k], 1e-8)).toBe(true);
                        }
                    }
                });
            }
        });
    });
});
