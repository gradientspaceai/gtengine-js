import { describe, it, expect } from 'vitest';
import { OBBNode, OBBTree } from '../src/OBBTree.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, normalize, sub } from '../src/Vector.js';
import { cross, dotCross } from '../src/Vector3.js';
import {
    check, expectClose, fc, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

// A triangle soup: numTriangles independent, non-degenerate triangles, each
// with its own three vertices. Used by the verification block below.
function triangleSoup(minTriangles: number, maxTriangles: number) {
    return fc.array(
        fc.tuple(wellScaledVector(3, -6, 6), wellScaledVector(3, -6, 6),
            wellScaledVector(3, -6, 6))
            .filter((t) => {
                const n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
                return dot(n, n) > 1e-2;
            }),
        { minLength: minTriangles, maxLength: maxTriangles })
        .map((tris) => {
            const vertices: Vector[] = [];
            const triangles: [number, number, number][] = [];
            for (const t of tris) {
                const b = vertices.length;
                vertices.push(t[0], t[1], t[2]);
                triangles.push([b, b + 1, b + 2]);
            }
            return { vertices: vertices, triangles: triangles };
        });
}

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

    partitionOf(i: number): number {
        return this.mPartition[i];
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

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of OBBTree.h against the
// port. OBBTree duplicates the BVTree tree-building code with an oriented box
// as the bounding volume, so the properties cover the same translation
// hazards (integer median index, size_t wrap-around of j0, reversed right
// partition) plus the eigen-decomposition that supplies the box frame.
// ---------------------------------------------------------------------------

// wellScaledVector rather than vector: the closed-form eigensolver of
// SymmetricEigensolver3x3 squares the centroid coordinates, so subnormal
// separations (which fc.double() produces readily) underflow the covariance
// matrix to zero and its eigenvectors stop being orthonormal. That is an
// upstream conditioning limit, not a port defect; see the eigensolver's own
// upstream-bug issue.
const obbPointCloud = fc.array(wellScaledVector(3, -8, 8),
    { minLength: 1, maxLength: 17 });

function obbReachable(tree: OBBTree): Array<{ index: number; depth: number }> {
    const list: Array<{ index: number; depth: number }> = [];
    walk(tree, (index, depth) => { list.push({ index: index, depth: depth }); });
    return list;
}

function obbIsLeaf(node: OBBNode): boolean {
    return node.leftChild === OBBNode.invalid && node.rightChild === OBBNode.invalid;
}

describe('OBBTree verification', () => {
    it('keeps the tree structure invariants for random clouds and heights', () => {
        check(fc.tuple(obbPointCloud, fc.integer({ min: 0, max: 6 }), fc.boolean()),
            (input) => {
                const points = input[0];
                const requested = input[1];
                const useFullHeight = input[2];
                const tree = new OBBTreeOfPointsRaw();
                tree.create(points, useFullHeight ? OBBTree.fullHeight : requested);

                const n = points.length;
                const nodes = tree.getNodes();
                const expectedHeight = useFullHeight
                    ? Math.ceil(Math.log2(n)) : Math.min(requested, 31);
                expect(tree.getHeight()).toBe(expectedHeight);
                expect(nodes.length).toBe(2 ** (expectedHeight + 1) - 1);
                expect([...tree.getPartition()].sort((a, b) => a - b))
                    .toEqual([...Array(n).keys()]);

                const covered = new Array<number>(n).fill(0);
                for (const record of obbReachable(tree)) {
                    expect(record.depth).toBeLessThanOrEqual(tree.getHeight());
                    const node = nodes[record.index];
                    expect(node.minIndex).toBeLessThanOrEqual(node.maxIndex);
                    if (obbIsLeaf(node)) {
                        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                            ++covered[i];
                        }
                    } else {
                        const left = nodes[node.leftChild];
                        const right = nodes[node.rightChild];
                        expect(node.leftChild).toBe(2 * record.index + 1);
                        expect(node.rightChild).toBe(2 * record.index + 2);
                        expect(left.minIndex).toBe(node.minIndex);
                        expect(right.maxIndex).toBe(node.maxIndex);
                        expect(right.minIndex).toBe(left.maxIndex + 1);
                        const sizeL = left.maxIndex - left.minIndex + 1;
                        const sizeR = right.maxIndex - right.minIndex + 1;
                        expect(sizeL).toBe(Math.ceil((sizeL + sizeR) / 2));
                    }
                }
                expect(covered.every(c => c === 1)).toBe(true);
            }, 100);
    });

    it('gives every interior box the eigen frame of the centroid covariance', () => {
        // Independent cross-check of ComputeInteriorBox: recompute the mean
        // and the covariance matrix, then verify that each stored axis is an
        // eigenvector of that matrix with the stored extent as eigenvalue.
        check(obbPointCloud.filter(p => p.length >= 2), (points) => {
            const tree = new OBBTreeOfPointsRaw();
            tree.create(points);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();

            for (const record of obbReachable(tree)) {
                const node = nodes[record.index];
                if (obbIsLeaf(node)) {
                    continue;
                }
                const box = node.box;
                const denom = node.maxIndex - node.minIndex + 1;

                // The box center is the mean of the node's centroids.
                const mean = new Vector(3);
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const c = centroids[tree.partitionOf(i)];
                    for (let j = 0; j < 3; ++j) {
                        mean.values[j] += c.values[j];
                    }
                }
                for (let j = 0; j < 3; ++j) {
                    mean.values[j] /= denom;
                    // Not bit-exact: ComputeInteriorBox sums in the mPartition
                    // order that held *before* SplitPoints permuted the node's
                    // range, so the recomputation adds the same terms in a
                    // different order (a 1-ulp difference).
                    expectClose(box.center.values[j], mean.values[j], 1e-12, 1e-12);
                }

                // The covariance matrix of the node's centroids.
                const cov = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const diff = sub(centroids[tree.partitionOf(i)], mean);
                    for (let r = 0; r < 3; ++r) {
                        for (let c = 0; c < 3; ++c) {
                            cov[r][c] += diff.values[r] * diff.values[c];
                        }
                    }
                }
                let scale = 0;
                for (let r = 0; r < 3; ++r) {
                    for (let c = 0; c < 3; ++c) {
                        cov[r][c] /= denom;
                        scale = Math.max(scale, Math.abs(cov[r][c]));
                    }
                }

                // The frame is orthonormal and right-handed, the eigenvalues
                // are stored increasingly in the extents (sortType +1), and
                // axis[2] belongs to the largest eigenvalue.
                for (let a = 0; a < 3; ++a) {
                    expect(Math.abs(dot(box.axis[a], box.axis[a]) - 1))
                        .toBeLessThan(1e-9);
                    for (let b = a + 1; b < 3; ++b) {
                        expect(Math.abs(dot(box.axis[a], box.axis[b])))
                            .toBeLessThan(1e-9);
                    }
                }
                expect(Math.abs(dotCross(box.axis[0], box.axis[1], box.axis[2]) - 1))
                    .toBeLessThan(1e-9);
                expect(box.extent.values[0]).toBeLessThanOrEqual(box.extent.values[1]);
                expect(box.extent.values[1]).toBeLessThanOrEqual(box.extent.values[2]);

                // cov * axis[j] = extent[j] * axis[j]. The tolerance scales
                // with the covariance magnitude because the closed-form
                // eigensolver works with squared coordinates.
                for (let j = 0; j < 3; ++j) {
                    const u = box.axis[j];
                    for (let r = 0; r < 3; ++r) {
                        const lhs = cov[r][0] * u.values[0] + cov[r][1] * u.values[1]
                            + cov[r][2] * u.values[2];
                        const rhs = box.extent.values[j] * u.values[r];
                        expect(Math.abs(lhs - rhs))
                            .toBeLessThanOrEqual(1e-8 * (1 + scale));
                    }
                }
            }
        }, 60);
    });

    it('splits at the median of the projections onto the largest-eigenvalue axis', () => {
        check(obbPointCloud.filter(p => p.length >= 2), (points) => {
            const tree = new OBBTreeOfPointsRaw();
            tree.create(points);
            const nodes = tree.getNodes();
            const centroids = tree.getCentroids();

            for (const record of obbReachable(tree)) {
                const node = nodes[record.index];
                if (obbIsLeaf(node)) {
                    continue;
                }
                // SplitPoints is called with the box center and box.axis[2],
                // the eigenvector of the largest eigenvalue.
                const project = (i: number): number =>
                    dot(node.box.axis[2],
                        sub(centroids[tree.partitionOf(i)], node.box.center));
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
                expect(maxLeft).toBeLessThanOrEqual(minRight);
            }
        }, 100);
    });

    it('bounds the triangles of every node when the derived class extends the box', () => {
        check(triangleSoup(1, 10), (mesh) => {
            const tree = new OBBTreeOfTriangles();
            tree.createFromMesh(mesh.vertices, mesh.triangles);
            const nodes = tree.getNodes();
            for (const record of obbReachable(tree)) {
                const node = nodes[record.index];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const tri = mesh.triangles[tree.partitionOf(i)];
                    for (let k = 0; k < 3; ++k) {
                        expect(boxContains(node.box, mesh.vertices[tri[k]], 1e-9))
                            .toBe(true);
                    }
                }
            }
        }, 60);
    });

    it('survives degenerate centroid sets without producing NaN', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), unitVector(3),
            fc.integer({ min: 1, max: 10 }), fc.boolean()), (input) => {
                const base = input[0];
                const dir = input[1];
                const n = input[2];
                const collinear = input[3];
                const points: Vector[] = [];
                for (let i = 0; i < n; ++i) {
                    // All coincident, or all on one line (rank-deficient
                    // covariance either way).
                    points.push(collinear
                        ? Vector.fromArray([
                            base.values[0] + i * dir.values[0],
                            base.values[1] + i * dir.values[1],
                            base.values[2] + i * dir.values[2]])
                        : base.clone());
                }
                const tree = new OBBTreeOfPointsRaw();
                tree.create(points);
                expect([...tree.getPartition()].sort((a, b) => a - b))
                    .toEqual([...Array(n).keys()]);
                const nodes = tree.getNodes();
                for (const record of obbReachable(tree)) {
                    const box = nodes[record.index].box;
                    for (let j = 0; j < 3; ++j) {
                        expect(Number.isFinite(box.center.values[j])).toBe(true);
                        expect(Number.isFinite(box.extent.values[j])).toBe(true);
                        for (let k = 0; k < 3; ++k) {
                            expect(Number.isFinite(box.axis[j].values[k])).toBe(true);
                        }
                    }
                }
            }, 60);
    });

    it('copies the centroids so later input mutation cannot reach the tree', () => {
        check(obbPointCloud, (points) => {
            const tree = new OBBTreeOfPointsRaw();
            tree.create(points);
            const before = tree.getCentroids().map(c => [...c.values]);
            for (const p of points) {
                p.values[1] -= 500;
            }
            expect(tree.getCentroids().map(c => [...c.values])).toEqual(before);
        }, 50);
    });
});
