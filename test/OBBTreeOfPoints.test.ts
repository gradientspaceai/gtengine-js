import { describe, it, expect } from 'vitest';
import { OBBNode, OBBTree } from '../src/OBBTree.js';
import { OBBTreeOfPoints } from '../src/OBBTreeOfPoints.js';
import type { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, sub } from '../src/Vector.js';
import {
    check, fc, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification wave (V07): property-based re-check of OBBTreeOfPoints.h.
//
// wellScaledVector rather than vector: the covariance matrix that
// SymmetricEigensolver3x3 works on squares the coordinates, so subnormal
// point separations (which fc.double() produces readily) underflow it to zero
// and the eigenvectors stop being orthonormal. That is a conditioning limit
// of the eigensolver, not of this file.
// ---------------------------------------------------------------------------

const obbpCloud = fc.array(wellScaledVector(3, -8, 8),
    { minLength: 1, maxLength: 17 });

function obbpReachable(tree: OBBTreeOfPoints): number[] {
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

describe('OBBTreeOfPoints verification', () => {
    it('bounds every node by the points it represents', () => {
        check(fc.tuple(obbpCloud, fc.integer({ min: 0, max: 4 }), fc.boolean()),
            (input) => {
                const points = input[0];
                const tree = new OBBTreeOfPoints();
                tree.create(points, input[2] ? OBBTree.fullHeight : input[1]);

                const nodes = tree.getNodes();
                const stored = tree.getPoints();
                const partition = tree.getPartition();
                expect([...partition].sort((a, b) => a - b))
                    .toEqual([...Array(points.length).keys()]);

                const covered = new Array<number>(points.length).fill(0);
                for (const index of obbpReachable(tree)) {
                    const node = nodes[index];
                    expectOrthonormal(node.box);
                    for (let j = 0; j < 3; ++j) {
                        expect(node.box.extent.values[j])
                            .toBeGreaterThanOrEqual(0);
                    }
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        // The if/else-if update with pmin = pmax = 0 is still
                        // a correct min/max because the projections of the
                        // node's own centroids straddle zero (their mean is
                        // the projection origin), so containment holds.
                        expect(boxContains(node.box, stored[partition[i]], 1e-9))
                            .toBe(true);
                    }
                    if (isLeaf(node)) {
                        for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                            ++covered[i];
                        }
                    }
                }
                expect(covered.every(c => c === 1)).toBe(true);
            }, 60);
    });

    it('gives every leaf the degenerate box at its point', () => {
        check(obbpCloud, (points) => {
            const tree = new OBBTreeOfPoints();
            tree.create(points);
            const nodes = tree.getNodes();
            const stored = tree.getPoints();
            const partition = tree.getPartition();

            for (const index of obbpReachable(tree)) {
                const node = nodes[index];
                if (!isLeaf(node)) {
                    continue;
                }
                expect(node.minIndex).toBe(node.maxIndex);
                const p = stored[partition[node.minIndex]];
                expect([...node.box.center.values]).toEqual([...p.values]);
                expect([...node.box.extent.values]).toEqual([0, 0, 0]);
                for (let j = 0; j < 3; ++j) {
                    expect([...node.box.axis[j].values])
                        .toEqual([...Vector.unit(3, j).values]);
                }
                // The leaf box center must be a copy, not an alias of the
                // stored point: upstream assigns by value.
                expect(node.box.center).not.toBe(p);
            }
        }, 100);
    });

    it('keeps the mean of a node inside its interior box (upstream #103)', () => {
        // Upstream seeds pmin and pmax with the zero vector, so the box is
        // forced to contain the projection origin -- the mean of the node's
        // centroids. That makes the box conservative, and the property pins
        // the quirk rather than a tighter fit.
        check(obbpCloud.filter(p => p.length >= 2), (points) => {
            const tree = new OBBTreeOfPoints();
            tree.create(points);
            const nodes = tree.getNodes();
            const stored = tree.getPoints();
            const partition = tree.getPartition();

            for (const index of obbpReachable(tree)) {
                const node = nodes[index];
                if (isLeaf(node)) {
                    continue;
                }
                const mean = new Vector(3);
                const denom = node.maxIndex - node.minIndex + 1;
                for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                    const p = stored[partition[i]];
                    for (let j = 0; j < 3; ++j) {
                        mean.values[j] += p.values[j];
                    }
                }
                for (let j = 0; j < 3; ++j) {
                    mean.values[j] /= denom;
                }
                expect(boxContains(node.box, mean, 1e-9)).toBe(true);
            }
        }, 60);
    });

    it('rebuilds cleanly and copies the input on every create', () => {
        check(fc.tuple(obbpCloud, obbpCloud), (input) => {
            const tree = new OBBTreeOfPoints();
            tree.create(input[0]);
            tree.create(input[1]);
            expect(tree.getPoints().length).toBe(input[1].length);
            expect(tree.getHeight()).toBe(Math.ceil(Math.log2(input[1].length)));
            expect([...tree.getPartition()].sort((a, b) => a - b))
                .toEqual([...Array(input[1].length).keys()]);
            for (let i = 0; i < input[1].length; ++i) {
                expect(tree.getPoints()[i]).not.toBe(input[1][i]);
                expect([...tree.getPoints()[i].values])
                    .toEqual([...input[1][i].values]);
            }
        }, 60);
    });

    it('survives coincident and collinear point sets', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), unitVector(3),
            fc.integer({ min: 1, max: 10 }), fc.boolean()), (input) => {
                const points: Vector[] = [];
                for (let i = 0; i < input[2]; ++i) {
                    points.push(input[3]
                        ? Vector.fromArray([
                            input[0].values[0] + i * input[1].values[0],
                            input[0].values[1] + i * input[1].values[1],
                            input[0].values[2] + i * input[1].values[2]])
                        : input[0].clone());
                }
                const tree = new OBBTreeOfPoints();
                tree.create(points);
                const nodes = tree.getNodes();
                const partition = tree.getPartition();
                for (const index of obbpReachable(tree)) {
                    const node = nodes[index];
                    for (let j = 0; j < 3; ++j) {
                        expect(Number.isFinite(node.box.center.values[j]))
                            .toBe(true);
                        expect(Number.isFinite(node.box.extent.values[j]))
                            .toBe(true);
                    }
                    for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                        expect(boxContains(node.box,
                            tree.getPoints()[partition[i]], 1e-8)).toBe(true);
                    }
                }
            }, 60);
    });
});
