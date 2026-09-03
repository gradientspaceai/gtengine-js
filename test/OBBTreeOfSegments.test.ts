import { describe, it, expect } from 'vitest';
import { OBBNode } from '../src/OBBTree.js';
import { OBBTreeOfSegments } from '../src/OBBTreeOfSegments.js';
import type { OrientedBox } from '../src/OrientedBox.js';
import { Vector, dot, length, sub } from '../src/Vector.js';

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

// A polyline of n+1 vertices with n segments.
function makePolyline(n: number, rand: () => number): {
    vertices: Vector[];
    segments: [number, number][];
} {
    const vertices: Vector[] = [];
    for (let i = 0; i <= n; ++i) {
        vertices.push(Vector.fromArray([i + rand(), 2 * rand(), 3 * rand() - 1.5]));
    }
    const segments: [number, number][] = [];
    for (let i = 0; i < n; ++i) {
        segments.push([i, i + 1]);
    }
    return { vertices: vertices, segments: segments };
}

describe('OBBTreeOfSegments', () => {
    it('computes segment midpoints as the centroids and copies its inputs', () => {
        const tree = new OBBTreeOfSegments();
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([2, 0, 0]),
            Vector.fromArray([2, 4, 0]),
            Vector.fromArray([2, 4, 6])
        ];
        const segments: [number, number][] = [[0, 1], [1, 2], [2, 3]];
        tree.createFromSegments(vertices, segments);

        const centroids = tree.getCentroids();
        expect(centroids.length).toBe(3);
        expect(centroids[0].values).toEqual([1, 0, 0]);
        expect(centroids[1].values).toEqual([2, 2, 0]);
        expect(centroids[2].values).toEqual([2, 4, 3]);

        expect(tree.getSegments()).toEqual([[0, 1], [1, 2], [2, 3]]);
        expect(tree.getVertices().length).toBe(4);

        // Mutating the inputs afterward does not change the tree.
        vertices[0].values[0] = 100;
        segments[0][1] = 3;
        expect(tree.getVertices()[0].values[0]).toBe(0);
        expect(tree.getSegments()[0][1]).toBe(1);
    });

    it('leaf boxes are degenerate boxes aligned with their segments', () => {
        const tree = new OBBTreeOfSegments();
        const rand = makeRand(31415);
        const { vertices, segments } = makePolyline(8, rand);
        tree.createFromSegments(vertices, segments);

        const partition = tree.getPartition();
        let leafCount = 0;
        for (const node of tree.getNodes()) {
            if (node.minIndex === OBBNode.invalid || !isLeaf(node)) {
                continue;
            }
            ++leafCount;
            const seg = segments[partition[node.minIndex]];
            const v0 = vertices[seg[0]];
            const v1 = vertices[seg[1]];
            const box = node.box;

            // The center is the segment midpoint.
            for (let j = 0; j < 3; ++j) {
                expect(box.center.values[j]).toBeCloseTo(
                    0.5 * (v0.values[j] + v1.values[j]), 12);
            }

            // axis[0] is the unit-length segment direction and extent[0] is
            // half the segment length; the box is degenerate otherwise.
            const halfLength = 0.5 * length(sub(v1, v0));
            expect(box.extent.values[0]).toBeCloseTo(halfLength, 12);
            expect(box.extent.values[1]).toBe(0);
            expect(box.extent.values[2]).toBe(0);
            expectOrthonormal(box);
            const direction = sub(v1, v0);
            expect(Math.abs(dot(box.axis[0], direction))).toBeCloseTo(
                2 * halfLength, 10);

            // The segment endpoints are on the boundary of the leaf box.
            expect(boxContains(box, v0)).toBe(true);
            expect(boxContains(box, v1)).toBe(true);
        }
        expect(leafCount).toBe(8);
    });

    it('interior boxes contain the endpoints of their segments', () => {
        const tree = new OBBTreeOfSegments();
        const rand = makeRand(2718281);
        const { vertices, segments } = makePolyline(21, rand);
        tree.createFromSegments(vertices, segments);

        const partition = tree.getPartition();
        expect([...partition].sort((a, b) => a - b)).toEqual(
            segments.map((_s, i) => i));

        let interiorCount = 0;
        for (const node of tree.getNodes()) {
            if (node.minIndex === OBBNode.invalid) {
                continue;
            }
            expectOrthonormal(node.box);
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const seg = segments[partition[i]];
                expect(boxContains(node.box, vertices[seg[0]])).toBe(true);
                expect(boxContains(node.box, vertices[seg[1]])).toBe(true);
            }
            if (!isLeaf(node)) {
                ++interiorCount;
                for (let j = 0; j < 3; ++j) {
                    expect(node.box.extent.values[j]).toBeGreaterThanOrEqual(0);
                }
            }
        }
        expect(interiorCount).toBeGreaterThan(0);
    });

    it('the root box of a symmetric segment set is centered as upstream computes', () => {
        // Four segments parallel to the x axis whose midpoints are symmetric
        // about the origin in the y and z directions. The mean of the
        // midpoints is the origin.
        const tree = new OBBTreeOfSegments();
        const vertices: Vector[] = [];
        const segments: [number, number][] = [];
        const offsets = [[-1, -2], [-1, 2], [1, -2], [1, 2]];
        for (let i = 0; i < 4; ++i) {
            vertices.push(Vector.fromArray([-2, offsets[i][0], offsets[i][1]]));
            vertices.push(Vector.fromArray([+2, offsets[i][0], offsets[i][1]]));
            segments.push([2 * i, 2 * i + 1]);
        }
        tree.createFromSegments(vertices, segments);

        const box = tree.getNodes()[0].box;
        for (let j = 0; j < 3; ++j) {
            expect(box.center.values[j]).toBeCloseTo(0, 12);
        }
        // The box must contain every endpoint. The centroid variances along
        // x, y and z are 0, 1 and 4, so the eigenvalues are distinct and the
        // box axes are the coordinate axes (up to sign and order). The
        // extents are then the half-widths 2, 1 and 2 of the endpoint set.
        for (const v of vertices) {
            expect(boxContains(box, v)).toBe(true);
        }
        const halfWidths = [2, 1, 2];
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

    it('handles a single segment (height 0)', () => {
        const tree = new OBBTreeOfSegments();
        const vertices = [
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([0, 0, 4])
        ];
        tree.createFromSegments(vertices, [[0, 1]]);

        expect(tree.getHeight()).toBe(0);
        expect(tree.getNodes().length).toBe(1);
        const root = tree.getNodes()[0];
        expect(isLeaf(root)).toBe(true);
        expect(root.box.center.values).toEqual([0, 0, 2]);
        expect(root.box.axis[0].values).toEqual([0, 0, 1]);
        expect(root.box.extent.values[0]).toBeCloseTo(2, 12);
        expect(root.box.extent.values[1]).toBe(0);
        expect(root.box.extent.values[2]).toBe(0);
        expectOrthonormal(root.box);
    });

    it('honors a user-specified height', () => {
        const tree = new OBBTreeOfSegments();
        const rand = makeRand(999983);
        const { vertices, segments } = makePolyline(8, rand);
        tree.createFromSegments(vertices, segments, 1);

        expect(tree.getHeight()).toBe(1);
        expect(tree.getNodes().length).toBe(3);
        const nodes = tree.getNodes();
        expect(isLeaf(nodes[1])).toBe(true);
        expect(isLeaf(nodes[2])).toBe(true);
        expect(nodes[1].maxIndex - nodes[1].minIndex + 1).toBe(4);
        expect(nodes[2].maxIndex - nodes[2].minIndex + 1).toBe(4);

        const partition = tree.getPartition();
        for (const node of [nodes[1], nodes[2]]) {
            for (let i = node.minIndex; i <= node.maxIndex; ++i) {
                const seg = segments[partition[i]];
                expect(boxContains(node.box, vertices[seg[0]])).toBe(true);
                expect(boxContains(node.box, vertices[seg[1]])).toBe(true);
            }
        }
    });

    it('throws for invalid inputs', () => {
        const tree = new OBBTreeOfSegments();
        expect(() => tree.createFromSegments([], [])).toThrow('Invalid input.');
        expect(() => tree.createFromSegments(
            [Vector.zero(3), Vector.unit(3, 0)], [])).toThrow('Invalid input.');
    });
});
