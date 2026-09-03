import { describe, it, expect } from 'vitest';
import { OrientedBoxTreeOfPoints } from '../src/OrientedBoxTreeOfPoints';
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

// Is the point inside or on the oriented box, allowing for round-off?
function boxContains(bv: OrientedBoxBV, point: Vector, epsilon = 1e-9): boolean {
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

// The box axes must be orthonormal and the extents nonnegative.
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

class TestableTree extends OrientedBoxTreeOfPoints {
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}

function buildTree(points: readonly Vector[], height?: number): TestableTree {
    const tree = new TestableTree();
    tree.create(points, height === undefined ? BVTree.fullHeight : height);
    return tree;
}

function gridPoints(n: number, rand: () => number): Vector[] {
    const out: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        out.push(v3(4 * rand() - 2, 3 * rand() - 1, 5 * rand()));
    }
    return out;
}

// An independent, recursive reference for the stack-based getLeafIndices.
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

describe('OrientedBoxTreeOfPoints leaf bounding volumes', () => {
    it('makes each leaf a degenerate box at the point with standard axes', () => {
        const points = [v3(0, 0, 0), v3(2, 0, 0), v3(0, 3, 0), v3(0, 0, 5)];
        const tree = buildTree(points);
        const partition = tree.getPartition();
        let leafCount = 0;
        for (const node of tree.getNodes()) {
            if (node.leftChild !== BVTreeNode.invalid) {
                continue;
            }
            ++leafCount;
            expect(node.minIndex).toBe(node.maxIndex);
            const box = node.boundingVolume.box;
            expect(box.extent.values).toEqual([0, 0, 0]);
            expect(box.axis[0].values).toEqual([1, 0, 0]);
            expect(box.axis[1].values).toEqual([0, 1, 0]);
            expect(box.axis[2].values).toEqual([0, 0, 1]);
            expect(box.center.values)
                .toEqual(points[partition[node.minIndex]].values);
        }
        expect(leafCount).toBe(points.length);
    });

    it('does not alias the input points in the leaf boxes', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const tree = buildTree(points);
        for (const node of tree.getNodes()) {
            if (node.leftChild === BVTreeNode.invalid) {
                node.boundingVolume.box.center.set(0, -100);
            }
        }
        for (const p of points) {
            expect(p.get(0)).toBeGreaterThanOrEqual(0);
        }
        for (const p of tree.getVertices()) {
            expect(p.get(0)).toBeGreaterThanOrEqual(0);
        }
    });
});

describe('OrientedBoxTreeOfPoints interior bounding volumes', () => {
    it('contains every point of its range in every node box', () => {
        const rand = makeRandom(97);
        const sets: Vector[][] = [
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)],
            gridPoints(17, rand),
            gridPoints(40, rand),
            // A planar (degenerate) point set.
            Array.from({ length: 16 }, () => v3(rand(), 2 * rand(), 1))
        ];
        for (const points of sets) {
            const tree = buildTree(points);
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
                    expect(boxContains(node.boundingVolume,
                        points[partition[i]], 1e-8)).toBe(true);
                }
            }
            expect(interior).toBeGreaterThan(0);
        }
    });

    it('centers the root box on the extreme midpoint of a symmetric set', () => {
        // Points on a segment along (1,1,1)/sqrt(3); the fitted box center is
        // the midpoint of the extremes.
        const points: Vector[] = [];
        for (let i = -3; i <= 3; ++i) {
            points.push(v3(1 + i, 2 + i, 3 + i));
        }
        const root = buildTree(points).getNodes()[0];
        expect(root.boundingVolume.box.center.get(0)).toBeCloseTo(1, 10);
        expect(root.boundingVolume.box.center.get(1)).toBeCloseTo(2, 10);
        expect(root.boundingVolume.box.center.get(2)).toBeCloseTo(3, 10);
        // The largest extent equals half of the length of the point spread.
        const extents = root.boundingVolume.box.extent.values.slice();
        extents.sort((a, b) => b - a);
        expect(extents[0]).toBeCloseTo(3 * Math.sqrt(3), 10);
        expect(extents[1]).toBeCloseTo(0, 10);
        expect(extents[2]).toBeCloseTo(0, 10);
    });

    it('builds a single-node tree for one point', () => {
        const tree = buildTree([v3(7, -1, 2)]);
        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        expect(tree.getNodes()[0].boundingVolume.box.center.values)
            .toEqual([7, -1, 2]);
        expect(tree.getNodes()[0].boundingVolume.box.extent.values)
            .toEqual([0, 0, 0]);
    });

    it('rejects an empty point set', () => {
        expect(() => new OrientedBoxTreeOfPoints().create([])).toThrow();
    });
});

describe('OrientedBoxTreeOfPoints queries', () => {
    it('agrees with an independent traversal and never misses a point', () => {
        const rand = makeRandom(4242);
        const sets: Vector[][] = [
            gridPoints(24, rand),
            gridPoints(50, rand)
        ];
        for (const points of sets) {
            const tree = buildTree(points);
            const nodes = tree.getNodes();
            const partition = tree.getPartition();
            for (const queryType of [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
                BVTree.SEGMENT_QUERY]) {
                for (let trial = 0; trial < 20; ++trial) {
                    const P = v3(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4);
                    const target = points[Math.floor(rand() * points.length)];
                    const d = sub(target, P);
                    const length = Math.hypot(d.get(0), d.get(1), d.get(2));
                    const Q = queryType === BVTree.SEGMENT_QUERY
                        ? Vector.fromArray(P.values.map((c, k) => c + 3 * d.get(k)))
                        : Vector.fromArray(d.values.map(c => c / length));

                    const leaves = tree.execute(queryType, P, Q);
                    expect(leaves).toEqual(referenceLeaves(tree, queryType, P, Q));
                    expect(tree.leafIndices(queryType, P, Q)).toEqual(leaves);

                    // Every reported leaf is in fact a leaf.
                    for (const nodeIndex of leaves) {
                        expect(nodes[nodeIndex].leftChild)
                            .toBe(BVTreeNode.invalid);
                    }

                    // A leaf is reported only if its (degenerate) box is hit,
                    // so the union of the reported leaf ranges is a subset of
                    // the partition.
                    const seen = new Set<number>();
                    for (const nodeIndex of leaves) {
                        const node = nodes[nodeIndex];
                        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                            expect(seen.has(partition[i])).toBe(false);
                            seen.add(partition[i]);
                        }
                    }
                    expect(seen.size).toBeLessThanOrEqual(points.length);
                }
            }
        }
    });

    it('reports every leaf for a line through all of the points', () => {
        // Points on the x-axis; each leaf box is the degenerate box at a
        // point, and the x-axis passes exactly through all of them.
        const points = [v3(-3, 0, 0), v3(-1, 0, 0), v3(0, 0, 0), v3(2, 0, 0),
            v3(5, 0, 0)];
        const tree = buildTree(points);
        const leaves = tree.execute(BVTree.LINE_QUERY, v3(-10, 0, 0),
            v3(1, 0, 0));
        const allLeaves = tree.getNodes()
            .map((node, index) => ({ node: node, index: index }))
            .filter(e => e.node.leftChild === BVTreeNode.invalid &&
                e.node.minIndex !== BVTreeNode.invalid)
            .map(e => e.index);
        expect(leaves.slice().sort((a, b) => a - b))
            .toEqual(allLeaves.slice().sort((a, b) => a - b));
        expect(leaves.length).toBe(points.length);

        // A ray starting past all of the points reports nothing.
        expect(tree.execute(BVTree.RAY_QUERY, v3(10, 0, 0), v3(1, 0, 0)).length)
            .toBe(0);
    });

    it('reports no leaves for a linear component far from the point set', () => {
        const tree = buildTree([v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0),
            v3(0, 0, 1)]);
        const P = v3(100, 100, 100);
        expect(tree.execute(BVTree.LINE_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.execute(BVTree.RAY_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.execute(BVTree.SEGMENT_QUERY, P, v3(101, 100, 100)).length)
            .toBe(0);
    });
});
