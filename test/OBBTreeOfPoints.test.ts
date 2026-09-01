import { describe, it, expect } from 'vitest';
import { OBBNode, OBBTree } from '../src/OBBTree';
import { OBBTreeOfPoints } from '../src/OBBTreeOfPoints';
import type { OrientedBox } from '../src/OrientedBox';
import { Vector, dot, sub } from '../src/Vector';

// Whether the point is inside or on the box, |Dot(X-C, U[j])| <= e[j].
function boxContains(box: OrientedBox, p: Vector, eps: number = 1e-9): boolean {
    const diff = sub(p, box.center);
    for (let j = 0; j < 3; ++j) {
        if (Math.abs(dot(diff, box.axis[j])) > box.extent.values[j] + eps) {
            return false;
        }
    }
    return true;
}

function expectOrthonormal(box: OrientedBox): void {
    for (let i = 0; i < 3; ++i) {
        expect(dot(box.axis[i], box.axis[i])).toBeCloseTo(1, 10);
        for (let j = i + 1; j < 3; ++j) {
            expect(dot(box.axis[i], box.axis[j])).toBeCloseTo(0, 10);
        }
    }
}

function isLeaf(node: OBBNode): boolean {
    return node.leftChild === OBBNode.invalid && node.rightChild === OBBNode.invalid;
}

function makeRand(seed0: number): () => number {
    let seed = seed0;
    return () => {
        seed = (1103515245 * seed + 12345) % 2147483648;
        return seed / 2147483648;
    };
}

describe('OBBTreeOfPoints', () => {
    it('uses the points as the centroids', () => {
        const tree = new OBBTreeOfPoints();
        const points = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([2, 1, 0]),
            Vector.fromArray([3, 1, 1])
        ];
        tree.create(points);

        expect(tree.getPoints().length).toBe(4);
        for (let i = 0; i < points.length; ++i) {
            expect(tree.getPoints()[i].values).toEqual(points[i].values);
        }

        // 4 = 2^2 leaves, so the height is 2 and there are 2^3-1 = 7 nodes.
        expect(tree.getHeight()).toBe(2);
        expect(tree.getNodes().length).toBe(7);
        expect(tree.getNodes()[0].minIndex).toBe(0);
        expect(tree.getNodes()[0].maxIndex).toBe(3);
    });

    it('leaf boxes are degenerate boxes at the points', () => {
        const tree = new OBBTreeOfPoints();
        const points = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([-4, 5, 6]),
            Vector.fromArray([7, -8, 9])
        ];
        tree.create(points);

        const partition = tree.getPartition();
        for (const node of tree.getNodes()) {
            if (node.minIndex === OBBNode.invalid || !isLeaf(node)) {
                continue;
            }
            const point = points[partition[node.minIndex]];
            expect(node.box.center.values).toEqual(point.values);
            expect(node.box.extent.values).toEqual([0, 0, 0]);
            expect(node.box.axis[0].values).toEqual([1, 0, 0]);
            expect(node.box.axis[1].values).toEqual([0, 1, 0]);
            expect(node.box.axis[2].values).toEqual([0, 0, 1]);
        }
    });

    it('interior boxes contain the points they represent', () => {
        const tree = new OBBTreeOfPoints();
        const rand = makeRand(20260901);
        const points: Vector[] = [];
        for (let i = 0; i < 40; ++i) {
            // An elongated, rotated cloud so the covariance axes are not the
            // coordinate axes.
            const s = 10 * rand() - 5;
            const u = 0.4 * rand() - 0.2;
            const v = 0.4 * rand() - 0.2;
            points.push(Vector.fromArray([s + u, 2 * s - v, -s + u + v]));
        }
        tree.create(points);

        const partition = tree.getPartition();
        expect([...partition].sort((a, b) => a - b)).toEqual(
            points.map((_p, i) => i));

        let interiorCount = 0;
        for (const node of tree.getNodes()) {
            if (node.minIndex === OBBNode.invalid) {
                continue;
            }
            expectOrthonormal(node.box);
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                expect(boxContains(node.box, points[partition[i]])).toBe(true);
            }
            if (!isLeaf(node)) {
                ++interiorCount;
                // The extents are nonnegative.
                for (let j = 0; j < 3; ++j) {
                    expect(node.box.extent.values[j]).toBeGreaterThanOrEqual(0);
                }
            }
        }
        expect(interiorCount).toBeGreaterThan(0);

        // Sibling ranges partition the parent range.
        const nodes = tree.getNodes();
        for (const node of nodes) {
            if (node.minIndex === OBBNode.invalid || isLeaf(node)) {
                continue;
            }
            expect(nodes[node.leftChild].minIndex).toBe(node.minIndex);
            expect(nodes[node.rightChild].maxIndex).toBe(node.maxIndex);
            expect(nodes[node.rightChild].minIndex).toBe(
                nodes[node.leftChild].maxIndex + 1);
        }
    });

    it('the root box of a symmetric cloud is centered as upstream computes', () => {
        // Eight points at the corners of the box [-1,1]x[-2,2]x[-3,3]. The
        // mean is the origin and the covariance matrix is diagonal, so the
        // box axes are the coordinate axes (in increasing eigenvalue order)
        // and the extents are the half-widths.
        const tree = new OBBTreeOfPoints();
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            points.push(Vector.fromArray([
                (i & 1) ? 1 : -1,
                (i & 2) ? 2 : -2,
                (i & 4) ? 3 : -3
            ]));
        }
        tree.create(points);

        const box = tree.getNodes()[0].box;
        for (let j = 0; j < 3; ++j) {
            expect(box.center.values[j]).toBeCloseTo(0, 12);
        }
        // axis[j] is +/- the coordinate axis of the j-th smallest variance.
        const halfWidths = [1, 2, 3];
        for (let j = 0; j < 3; ++j) {
            let k = 0;
            for (let m = 1; m < 3; ++m) {
                if (Math.abs(box.axis[j].values[m]) > Math.abs(box.axis[j].values[k])) {
                    k = m;
                }
            }
            expect(Math.abs(box.axis[j].values[k])).toBeCloseTo(1, 10);
            expect(box.extent.values[j]).toBeCloseTo(halfWidths[k], 10);
        }
    });

    it('handles a single point (height 0)', () => {
        const tree = new OBBTreeOfPoints();
        tree.create([Vector.fromArray([5, -6, 7])]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(isLeaf(root)).toBe(true);
        expect(root.minIndex).toBe(0);
        expect(root.maxIndex).toBe(0);
        expect(root.box.center.values).toEqual([5, -6, 7]);
        expect(root.box.extent.values).toEqual([0, 0, 0]);
    });

    it('honors a user-specified height', () => {
        const tree = new OBBTreeOfPoints();
        const rand = makeRand(777);
        const points: Vector[] = [];
        for (let i = 0; i < 16; ++i) {
            points.push(Vector.fromArray([rand(), 2 * rand(), 3 * rand()]));
        }
        tree.create(points, 1);

        expect(tree.getHeight()).toBe(1);
        expect(tree.getNodes().length).toBe(3);
        const nodes = tree.getNodes();
        expect(isLeaf(nodes[1])).toBe(true);
        expect(isLeaf(nodes[2])).toBe(true);
        expect(nodes[1].maxIndex - nodes[1].minIndex + 1).toBe(8);
        expect(nodes[2].maxIndex - nodes[2].minIndex + 1).toBe(8);

        // The depth-1 nodes have i0 < i1, so they use the interior-box code
        // and must contain their points.
        const partition = tree.getPartition();
        for (const node of [nodes[1], nodes[2]]) {
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                expect(boxContains(node.box, points[partition[i]])).toBe(true);
            }
        }

        // The full-height default is available as OBBTree.fullHeight.
        expect(OBBTree.fullHeight).toBe(Number.MAX_SAFE_INTEGER);
    });

    it('throws when there are no points', () => {
        const tree = new OBBTreeOfPoints();
        expect(() => tree.create([])).toThrow('Invalid input.');
    });
});
