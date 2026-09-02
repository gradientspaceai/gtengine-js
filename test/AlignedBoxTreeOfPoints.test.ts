import { describe, it, expect } from 'vitest';
import { AlignedBoxTreeOfPoints } from '../src/AlignedBoxTreeOfPoints';
import { AlignedBoxBV } from '../src/AlignedBoxBV';
import { BVTree, BVTreeNode } from '../src/BVTree';
import { Vector } from '../src/Vector';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

function randomPoints(seed: number, count: number, scale: number): Vector[] {
    const rand = makeRandom(seed);
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        points.push(v3(scale * (2 * rand() - 1), scale * (2 * rand() - 1),
            scale * (2 * rand() - 1)));
    }
    return points;
}

// The component-wise min/max over a set of points, computed independently of
// the port.
function tightBounds(points: readonly Vector[], indices: readonly number[]):
    { min: number[], max: number[] } {
    const min = [Infinity, Infinity, Infinity];
    const max = [-Infinity, -Infinity, -Infinity];
    for (const i of indices) {
        for (let k = 0; k < 3; ++k) {
            min[k] = Math.min(min[k], points[i].get(k));
            max[k] = Math.max(max[k], points[i].get(k));
        }
    }
    return { min: min, max: max };
}

// Exposes the protected getLeafIndices, the pattern used by
// test/BVTreeOfTriangles.test.ts.
class TestableTree extends AlignedBoxTreeOfPoints {
    leafIndices(queryType: number, P: Vector, Q: Vector): number[] {
        return this.getLeafIndices(queryType, P, Q);
    }
}

function buildTree(points: readonly Vector[], height?: number): TestableTree {
    const tree = new TestableTree();
    tree.create(points, height === undefined ? BVTree.fullHeight : height);
    return tree;
}

// The point indices a tree reports as candidates: the union over the reported
// leaf nodes of their partition ranges.
function candidates(tree: TestableTree, nodeIndices: readonly number[]): Set<number> {
    const nodes = tree.getNodes();
    const partition = tree.getPartition();
    const result = new Set<number>();
    for (const nodeIndex of nodeIndices) {
        const node = nodes[nodeIndex];
        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
            result.add(partition[i]);
        }
    }
    return result;
}

describe('AlignedBoxTreeOfPoints construction', () => {
    it('builds a single-node tree for one point with a degenerate leaf box', () => {
        const tree = buildTree([v3(2, -3, 5)]);
        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(root.leftChild).toBe(BVTreeNode.invalid);
        expect(root.rightChild).toBe(BVTreeNode.invalid);
        expect(root.boundingVolume.box.min.values).toEqual([2, -3, 5]);
        expect(root.boundingVolume.box.max.values).toEqual([2, -3, 5]);
    });

    it('gives the root the tight box of all the points', () => {
        const points = [v3(0, 0, 0), v3(1, -2, 3), v3(-4, 5, 1), v3(2, 2, -6)];
        const tree = buildTree(points);
        const root = tree.getNodes()[0];
        expect(root.boundingVolume.box.min.values).toEqual([-4, -2, -6]);
        expect(root.boundingVolume.box.max.values).toEqual([2, 5, 3]);
    });

    it('gives every visited node the tight box of its point range', () => {
        const sets: Vector[][] = [
            [v3(0, 0, 0), v3(1, 1, 1)],
            randomPoints(101, 9, 3),
            randomPoints(202, 16, 10),
            randomPoints(303, 33, 1)
        ];
        for (const points of sets) {
            const tree = buildTree(points);
            const partition = tree.getPartition();
            let visited = 0;
            let leafCount = 0;
            for (const node of tree.getNodes()) {
                if (node.minIndex === BVTreeNode.invalid) {
                    continue;
                }
                ++visited;
                const range: number[] = [];
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    range.push(partition[i]);
                }
                const expected = tightBounds(points, range);
                expect(node.boundingVolume.box.min.values).toEqual(expected.min);
                expect(node.boundingVolume.box.max.values).toEqual(expected.max);
                if (node.minIndex === node.maxIndex) {
                    ++leafCount;
                    // The leaf box for a point primitive is degenerate.
                    expect(node.boundingVolume.box.min.values)
                        .toEqual(node.boundingVolume.box.max.values);
                }
            }
            expect(leafCount).toBe(points.length);
            expect(visited).toBeGreaterThanOrEqual(points.length);
        }
    });

    it('does not alias the input vertices in the node boxes', () => {
        const points = [v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2)];
        const tree = buildTree(points);
        const root = tree.getNodes()[0];
        root.boundingVolume.box.min.set(0, -100);
        expect(points[0].get(0)).toBe(0);
        expect(tree.getVertices()[0].get(0)).toBe(0);
    });

    it('nests child boxes inside parent boxes', () => {
        const points = randomPoints(404, 24, 7);
        const tree = buildTree(points);
        const nodes = tree.getNodes();
        let interiorCount = 0;
        for (const node of nodes) {
            if (node.leftChild === BVTreeNode.invalid) {
                continue;
            }
            ++interiorCount;
            for (const c of [node.leftChild, node.rightChild]) {
                const child = nodes[c];
                for (let k = 0; k < 3; ++k) {
                    expect(child.boundingVolume.box.min.get(k))
                        .toBeGreaterThanOrEqual(node.boundingVolume.box.min.get(k));
                    expect(child.boundingVolume.box.max.get(k))
                        .toBeLessThanOrEqual(node.boundingVolume.box.max.get(k));
                }
            }
        }
        expect(interiorCount).toBeGreaterThan(0);
    });

    it('respects a user-specified height, leaving multi-point leaves tight', () => {
        const points = randomPoints(505, 20, 4);
        const full = buildTree(points);
        const shallow = buildTree(points, 2);
        expect(shallow.getHeight()).toBe(2);
        expect(full.getHeight()).toBeGreaterThan(2);

        const partition = shallow.getPartition();
        let multiLeaf = 0;
        for (const node of shallow.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid ||
                node.leftChild !== BVTreeNode.invalid) {
                continue;
            }
            if (node.maxIndex > node.minIndex) {
                ++multiLeaf;
            }
            const range: number[] = [];
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                range.push(partition[i]);
            }
            const expected = tightBounds(points, range);
            expect(node.boundingVolume.box.min.values).toEqual(expected.min);
            expect(node.boundingVolume.box.max.values).toEqual(expected.max);
        }
        expect(multiLeaf).toBeGreaterThan(0);
    });

    it('handles coincident points (all boxes degenerate)', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const tree = buildTree(points);
        for (const node of tree.getNodes()) {
            if (node.minIndex === BVTreeNode.invalid) {
                continue;
            }
            expect(node.boundingVolume.box.min.values).toEqual([1, 2, 3]);
            expect(node.boundingVolume.box.max.values).toEqual([1, 2, 3]);
        }
    });

    it('rejects an empty point set', () => {
        const tree = new AlignedBoxTreeOfPoints();
        expect(() => tree.create([])).toThrow();
    });
});

// An independent, recursive reference for the stack-based getLeafIndices: an
// interior node is descended when its box is met by the linear component and
// a leaf is reported without testing its own box (the upstream traversal, see
// BVTree.ts).
function referenceLeaves(tree: TestableTree, queryType: number, P: Vector,
    Q: Vector): number[] {
    const nodes = tree.getNodes();
    const hits = (bv: AlignedBoxBV): boolean =>
        queryType === BVTree.LINE_QUERY ? AlignedBoxBV.intersectLine(P, Q, bv)
            : (queryType === BVTree.RAY_QUERY
                ? AlignedBoxBV.intersectRay(P, Q, bv)
                : AlignedBoxBV.intersectSegment(P, Q, bv));

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

describe('AlignedBoxTreeOfPoints queries', () => {
    it('matches a recursive reference traversal on random queries', () => {
        const points = randomPoints(606, 40, 5);
        const tree = buildTree(points);
        const rand = makeRandom(909);

        const queryTypes = [BVTree.LINE_QUERY, BVTree.RAY_QUERY,
            BVTree.SEGMENT_QUERY];
        let totalReported = 0;
        for (const queryType of queryTypes) {
            for (let trial = 0; trial < 20; ++trial) {
                // Aim through a data point so the queries are not trivial.
                const target = points[Math.floor(rand() * points.length)];
                const P = v3(20 * (2 * rand() - 1), 20 * (2 * rand() - 1),
                    20 * (2 * rand() - 1));
                const dir = Vector.fromArray(
                    target.values.map((c, k) => c - P.get(k)));
                const length = Math.hypot(dir.get(0), dir.get(1), dir.get(2));
                const unit = Vector.fromArray(dir.values.map(c => c / length));
                // For the segment query, Q is the far endpoint.
                const Q = queryType === BVTree.SEGMENT_QUERY
                    ? Vector.fromArray(P.values.map((c, k) => c + 2 * dir.get(k)))
                    : unit;

                const leaves = tree.leafIndices(queryType, P, Q);
                expect(leaves).toEqual(referenceLeaves(tree, queryType, P, Q));
                totalReported += candidates(tree, leaves).size;
            }
        }

        // The queries aim at data points, so the traversals are not trivial.
        // (A single query can still report nothing: the direction is rounded,
        // so the line need not meet the razor-thin box of its target point.)
        expect(totalReported).toBeGreaterThan(0);
    });

    it('reports the leaf of a point hit exactly by an axis-parallel query', () => {
        // The line x = p.x, y = p.y passes exactly through p (no rounding in
        // the coordinates), so p's leaf must be reported, and every other
        // point whose degenerate box the line meets must be reported too.
        const points = randomPoints(606, 40, 5);
        const tree = buildTree(points);
        for (let t = 0; t < points.length; ++t) {
            const p = points[t];
            const P = v3(p.get(0), p.get(1), -100);
            const reported = candidates(tree,
                tree.leafIndices(BVTree.LINE_QUERY, P, v3(0, 0, 1)));
            expect(reported.has(t)).toBe(true);

            for (let i = 0; i < points.length; ++i) {
                const bv = new AlignedBoxBV();
                bv.box.min = points[i].clone();
                bv.box.max = points[i].clone();
                if (AlignedBoxBV.intersectLine(P, v3(0, 0, 1), bv)) {
                    expect(reported.has(i)).toBe(true);
                }
            }

            // The same query as a ray and as a segment reaching past p.
            const rayReported = candidates(tree,
                tree.leafIndices(BVTree.RAY_QUERY, P, v3(0, 0, 1)));
            expect(rayReported.has(t)).toBe(true);
            const segReported = candidates(tree,
                tree.leafIndices(BVTree.SEGMENT_QUERY, P,
                    v3(p.get(0), p.get(1), 100)));
            expect(segReported.has(t)).toBe(true);
        }
    });

    it('reports no leaves for a linear component that misses the root box', () => {
        const points = randomPoints(707, 16, 2);
        const tree = buildTree(points);
        const P = v3(100, 100, 100);
        expect(tree.leafIndices(BVTree.LINE_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.leafIndices(BVTree.RAY_QUERY, P, v3(1, 0, 0)).length).toBe(0);
        expect(tree.leafIndices(BVTree.SEGMENT_QUERY, P, v3(200, 100, 100)).length)
            .toBe(0);
    });

    it('execute() returns the same leaf indices as getLeafIndices', () => {
        const points = randomPoints(808, 12, 3);
        const tree = buildTree(points);
        const P = v3(-10, 0, 0);
        const Q = v3(1, 0, 0);
        expect(tree.execute(BVTree.LINE_QUERY, P, Q))
            .toEqual(tree.leafIndices(BVTree.LINE_QUERY, P, Q));
    });
});
