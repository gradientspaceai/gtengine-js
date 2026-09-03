import { describe, it, expect } from 'vitest';
import { OrientedBoxTreeOfSegments } from '../src/OrientedBoxTreeOfSegments';
import { OrientedBoxBV } from '../src/OrientedBoxBV';
import { BVTree, BVTreeNode } from '../src/BVTree';
import { Vector, dot, sub } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

type Seg = [number, number];

interface Mesh {
    vertices: Vector[];
    segments: Seg[];
}

function boxContains(bv: OrientedBoxBV, point: Vector, epsilon: number): boolean {
    const box = bv.box;
    const diff = sub(point, box.center);
    for (let k = 0; k < 3; ++k) {
        const y = dot(diff, box.axis[k]);
        if (Math.abs(y) > box.extent.get(k) + epsilon) {
            return false;
        }
    }
    return true;
}

function expectValidBox(bv: OrientedBoxBV): void {
    const box = bv.box;
    for (let i = 0; i < 3; ++i) {
        expect(box.extent.get(i)).toBeGreaterThanOrEqual(0);
        expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 10);
        for (let j = i + 1; j < 3; ++j) {
            expect(dot(box.axis[i], box.axis[j])).toBeCloseTo(0, 10);
        }
    }
}

class TestableTree extends OrientedBoxTreeOfSegments {
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}

function buildTree(mesh: Mesh, height?: number): TestableTree {
    const tree = new TestableTree();
    tree.createFromSegments(mesh.vertices, mesh.segments,
        height === undefined ? BVTree.fullHeight : height);
    return tree;
}

// A polyline of n segments with randomized vertices.
function makePolyline(n: number, rand: () => number): Mesh {
    const vertices: Vector[] = [];
    const segments: Seg[] = [];
    for (let i = 0; i <= n; ++i) {
        vertices.push(v3(6 * rand() - 3, 4 * rand() - 2, 5 * rand() - 1));
    }
    for (let i = 0; i < n; ++i) {
        segments.push([i, i + 1]);
    }
    return { vertices: vertices, segments: segments };
}

// The 12 edges of an axis-aligned box.
function boxEdges(lo: Vector, hi: Vector): Mesh {
    const vertices: Vector[] = [];
    for (let i = 0; i < 8; ++i) {
        vertices.push(v3(
            (i & 1) ? hi.get(0) : lo.get(0),
            (i & 2) ? hi.get(1) : lo.get(1),
            (i & 4) ? hi.get(2) : lo.get(2)));
    }
    const segments: Seg[] = [];
    for (let i = 0; i < 8; ++i) {
        for (const bit of [1, 2, 4]) {
            if ((i & bit) === 0) {
                segments.push([i, i | bit]);
            }
        }
    }
    return { vertices: vertices, segments: segments };
}

function referenceLeaves(tree: TestableTree, queryType: number, P: Vector,
    Q: Vector): number[] {
    const nodes = tree.getNodes();
    const hits = (bv: OrientedBoxBV): boolean =>
        queryType === BVTree.LINE_QUERY ? OrientedBoxBV.intersectLine(P, Q, bv)
            : (queryType === BVTree.RAY_QUERY
                ? OrientedBoxBV.intersectRay(P, Q, bv)
                : OrientedBoxBV.intersectSegment(P, Q, bv));

    const out: number[] = [];
    const visit = (nodeIndex: number): void => {
        const node = nodes[nodeIndex];
        if (node.leftChild !== BVTreeNode.invalid &&
            node.rightChild !== BVTreeNode.invalid) {
            if (hits(node.boundingVolume)) {
                visit(node.leftChild);
                visit(node.rightChild);
            }
        } else {
            out.push(nodeIndex);
        }
    };
    visit(0);
    return out;
}

describe('OrientedBoxTreeOfSegments leaf bounding volumes', () => {
    it('makes a leaf box tight around its segment with two zero extents', () => {
        // A single segment along (1,2,2)/3 with length 6.
        const mesh: Mesh = {
            vertices: [v3(1, 2, 3), v3(3, 6, 7)],
            segments: [[0, 1]]
        };
        const tree = buildTree(mesh);
        expect(tree.getNodes().length).toBe(1);
        const box = tree.getNodes()[0].boundingVolume.box;

        expect(box.center.get(0)).toBeCloseTo(2, 10);
        expect(box.center.get(1)).toBeCloseTo(4, 10);
        expect(box.center.get(2)).toBeCloseTo(5, 10);

        // Exactly two extents are exactly zero and the third is the half
        // length of the segment.
        const zeroCount = box.extent.values.filter(e => e === 0).length;
        expect(zeroCount).toBe(2);
        const maxExtent = Math.max(...box.extent.values);
        expect(maxExtent).toBeCloseTo(0.5 * Math.hypot(2, 4, 4), 10);

        // The endpoints are on the box.
        expect(boxContains(tree.getNodes()[0].boundingVolume, mesh.vertices[0],
            1e-9)).toBe(true);
        expect(boxContains(tree.getNodes()[0].boundingVolume, mesh.vertices[1],
            1e-9)).toBe(true);
    });

    it('zeroes exactly two extents on every leaf of a larger tree', () => {
        const rand = makeRandom(31337);
        for (const mesh of [makePolyline(12, rand), boxEdges(v3(0, 0, 0),
            v3(2, 3, 5))]) {
            const tree = buildTree(mesh);
            const partition = tree.getPartition();
            let leafCount = 0;
            for (const node of tree.getNodes()) {
                if (node.leftChild !== BVTreeNode.invalid ||
                    node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                ++leafCount;
                expect(node.minIndex).toBe(node.maxIndex);
                expectValidBox(node.boundingVolume);
                expect(node.boundingVolume.box.extent.values
                    .filter(e => e === 0).length).toBe(2);
                const seg = mesh.segments[partition[node.minIndex]];
                for (const j of seg) {
                    expect(boxContains(node.boundingVolume, mesh.vertices[j],
                        1e-8)).toBe(true);
                }
            }
            expect(leafCount).toBe(mesh.segments.length);
        }
    });

    it('handles a degenerate (zero-length) segment', () => {
        const mesh: Mesh = {
            vertices: [v3(1, 1, 1), v3(1, 1, 1), v3(4, 0, 0)],
            segments: [[0, 1], [1, 2]]
        };
        const tree = buildTree(mesh);
        const partition = tree.getPartition();
        for (const node of tree.getNodes()) {
            if (node.leftChild !== BVTreeNode.invalid ||
                node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            const seg = mesh.segments[partition[node.minIndex]];
            for (const j of seg) {
                expect(boxContains(node.boundingVolume, mesh.vertices[j],
                    1e-8)).toBe(true);
            }
        }
    });
});

describe('OrientedBoxTreeOfSegments interior bounding volumes', () => {
    it('contains every endpoint of its range in every node box', () => {
        const rand = makeRandom(555);
        const meshes: Mesh[] = [
            makePolyline(9, rand),
            makePolyline(25, rand),
            boxEdges(v3(-1, -2, -3), v3(4, 1, 2))
        ];
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            const partition = tree.getPartition();
            let interior = 0;
            for (const node of tree.getNodes()) {
                if (node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                expectValidBox(node.boundingVolume);
                if (node.maxIndex > node.minIndex) {
                    ++interior;
                }
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    for (const j of mesh.segments[partition[i]]) {
                        expect(boxContains(node.boundingVolume,
                            mesh.vertices[j], 1e-8)).toBe(true);
                    }
                }
            }
            expect(interior).toBeGreaterThan(0);
        }
    });

    it('nests the child boxes of the root inside no larger a point set', () => {
        // The child ranges partition the parent range, so every endpoint of a
        // child range is inside the parent box.
        const rand = makeRandom(8);
        const mesh = makePolyline(20, rand);
        const tree = buildTree(mesh);
        const nodes = tree.getNodes();
        const partition = tree.getPartition();
        let checked = 0;
        for (const node of nodes) {
            if (node.leftChild === BVTreeNode.invalid) {
                continue;
            }
            ++checked;
            for (const c of [node.leftChild, node.rightChild]) {
                const child = nodes[c];
                for (let i = child.minIndex; i <= child.maxIndex; ++i) {
                    for (const j of mesh.segments[partition[i]]) {
                        expect(boxContains(node.boundingVolume,
                            mesh.vertices[j], 1e-8)).toBe(true);
                    }
                }
            }
        }
        expect(checked).toBeGreaterThan(0);
    });

    it('respects a user-specified height', () => {
        const rand = makeRandom(77);
        const mesh = makePolyline(20, rand);
        const shallow = buildTree(mesh, 2);
        expect(shallow.getHeight()).toBe(2);
        expect(buildTree(mesh).getHeight()).toBeGreaterThan(2);

        const partition = shallow.getPartition();
        let multiLeaf = 0;
        for (const node of shallow.getNodes()) {
            if (node.leftChild !== BVTreeNode.invalid ||
                node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            if (node.maxIndex > node.minIndex) {
                ++multiLeaf;
            }
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                for (const j of mesh.segments[partition[i]]) {
                    expect(boxContains(node.boundingVolume, mesh.vertices[j],
                        1e-8)).toBe(true);
                }
            }
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('rejects degenerate creation input', () => {
        const tree = new OrientedBoxTreeOfSegments();
        expect(() => tree.createFromSegments([v3(0, 0, 0), v3(1, 0, 0)], []))
            .toThrow();
    });
});

describe('OrientedBoxTreeOfSegments queries', () => {
    it('agrees with an independent traversal on randomized queries', () => {
        const rand = makeRandom(1357);
        const meshes: Mesh[] = [
            makePolyline(15, rand),
            boxEdges(v3(0, 0, 0), v3(1, 1, 1))
        ];
        for (const mesh of meshes) {
            const tree = buildTree(mesh);
            const nodes = tree.getNodes();
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                for (let trial = 0; trial < 20; ++trial) {
                    const P = v3(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                    const d = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
                    const length = Math.hypot(d.get(0), d.get(1), d.get(2));
                    const Q = queryType === BVTree.SEGMENT_QUERY
                        ? Vector.fromArray(P.values.map((c, k) => c + 8 * d.get(k)))
                        : Vector.fromArray(d.values.map(c => c / length));

                    const leaves = tree.execute(queryType, P, Q);
                    expect(leaves).toEqual(referenceLeaves(tree, queryType, P, Q));
                    expect(tree.leafIndices(queryType, P, Q)).toEqual(leaves);
                    for (const nodeIndex of leaves) {
                        expect(nodes[nodeIndex].leftChild).toBe(BVTreeNode.invalid);
                    }
                }
            }
        }
    });

    it('culls the subtrees the line misses and keeps the one it hits', () => {
        // Four well-separated segments, so the tree is root -> 2 interior ->
        // 4 leaves. A line through the first segment culls the subtree that
        // holds the far cluster. As upstream, getLeafIndices tests only the
        // interior node volumes and reports the reached leaves without
        // testing their own volumes.
        const mesh: Mesh = {
            vertices: [
                v3(0, 0, 0), v3(1, 0, 0),
                v3(0, 1, 0), v3(1, 1, 0),
                v3(0, 100, 0), v3(1, 100, 0),
                v3(0, 101, 0), v3(1, 101, 0)
            ],
            segments: [[0, 1], [2, 3], [4, 5], [6, 7]]
        };
        const tree = buildTree(mesh);
        const partition = tree.getPartition();
        const leaves = tree.execute(BVTree.LINE_QUERY, v3(0.5, 0, -5),
            v3(0, 0, 1));
        expect(leaves.length).toBeGreaterThan(0);
        expect(leaves.length).toBeLessThan(4);
        const reached = new Set<number>();
        for (const nodeIndex of leaves) {
            const node = tree.getNodes()[nodeIndex];
            expect(node.leftChild).toBe(BVTreeNode.invalid);
            reached.add(partition[node.minIndex]);
        }
        // The near cluster is reached; the far cluster is culled.
        expect(reached.has(0)).toBe(true);
        expect(reached.has(2)).toBe(false);
        expect(reached.has(3)).toBe(false);

        // A line far from every segment reports nothing at all.
        expect(tree.execute(BVTree.LINE_QUERY, v3(500, 500, 500), v3(0, 0, 1))
            .length).toBe(0);
    });
});
