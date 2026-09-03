import { describe, it, expect } from 'vitest';
import {
    PolygonTree, PolygonTreeEx, PolygonTreeExNode
} from '../src/PolygonTree.js';
import { Vector } from '../src/Vector.js';

function pts(points: [number, number][]): Vector[] {
    return points.map(p => Vector.fromArray(p));
}

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// An independent point-in-triangle test using barycentric-style sign tests
// that does not depend on the winding order.
function referenceInTriangle(test: Vector, tri: [number, number, number],
    points: Vector[]): boolean {
    const [a, b, c] = tri.map(i => points[i].values);
    const t = test.values;
    const d1 = (t[0] - b[0]) * (a[1] - b[1]) - (a[0] - b[0]) * (t[1] - b[1]);
    const d2 = (t[0] - c[0]) * (b[1] - c[1]) - (b[0] - c[0]) * (t[1] - c[1]);
    const d3 = (t[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (t[1] - a[1]);
    const hasNeg = d1 < 0 || d2 < 0 || d3 < 0;
    const hasPos = d1 > 0 || d2 > 0 || d3 > 0;
    return !(hasNeg && hasPos);
}

describe('PolygonTree', () => {
    it('default-constructs an empty node', () => {
        const tree = new PolygonTree();
        expect(tree.polygon).toEqual([]);
        expect(tree.child).toEqual([]);
    });

    it('builds a nested hierarchy of outer and inner polygons', () => {
        // Outer square (counterclockwise) with an inner square hole
        // (clockwise) that itself contains a counterclockwise triangle.
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];

        const hole = new PolygonTree();
        hole.polygon = [4, 7, 6, 5];
        root.child.push(hole);

        const island = new PolygonTree();
        island.polygon = [8, 9, 10];
        hole.child.push(island);

        expect(root.child.length).toBe(1);
        expect(root.child[0].polygon).toEqual([4, 7, 6, 5]);
        expect(root.child[0].child[0].polygon).toEqual([8, 9, 10]);
        expect(root.child[0].child[0].child.length).toBe(0);
    });
});

describe('PolygonTreeExNode', () => {
    it('default-constructs with an empty polygon and zero chirality', () => {
        const node = new PolygonTreeExNode();
        expect(node.polygon).toEqual([]);
        expect(node.chirality).toBe(0);
        expect(node.triangulation).toEqual([]);
        expect(node.self).toBe(0);
        expect(node.parent).toBe(0);
        expect(node.minChild).toBe(0);
        expect(node.supChild).toBe(0);
    });
});

// A tree with a counterclockwise outer square [0,4]x[0,4] containing a
// clockwise inner square hole [1,3]x[1,3].
//
//   points 0..3 : outer square corners (CCW)
//   points 4..7 : inner square corners (CW as listed in the node)
//
// The interior region (between the two squares) is triangulated into 8
// triangles; the hole is triangulated into 2 triangles.
function buildTwoLevelTree(): { tree: PolygonTreeEx, points: Vector[] } {
    const points = pts([
        [0, 0], [4, 0], [4, 4], [0, 4],   // 0..3 outer
        [1, 1], [3, 1], [3, 3], [1, 3]    // 4..7 inner
    ]);

    const root = new PolygonTreeExNode();
    root.polygon = [0, 1, 2, 3];
    root.chirality = 1;
    root.self = 0;
    root.parent = 0;
    root.minChild = 1;
    root.supChild = 2;
    // Counterclockwise triangles covering the ring between the squares.
    root.triangulation = [
        [0, 1, 5], [0, 5, 4],
        [1, 2, 6], [1, 6, 5],
        [2, 3, 7], [2, 7, 6],
        [3, 0, 4], [3, 4, 7]
    ];

    const hole = new PolygonTreeExNode();
    hole.polygon = [4, 7, 6, 5];
    hole.chirality = -1;
    hole.self = 1;
    hole.parent = 0;
    hole.minChild = 0;
    hole.supChild = 0;
    // Clockwise triangles covering the hole.
    hole.triangulation = [
        [4, 7, 6], [4, 6, 5]
    ];

    const tree = new PolygonTreeEx();
    tree.nodes = [root, hole];
    tree.interiorTriangles = root.triangulation.slice();
    tree.interiorNodeIndices = root.triangulation.map(() => 0);
    tree.exteriorTriangles = hole.triangulation.slice();
    tree.exteriorNodeIndices = hole.triangulation.map(() => 1);
    tree.insideTriangles = tree.interiorTriangles.concat(tree.exteriorTriangles);
    tree.insideNodeIndices = tree.interiorNodeIndices.concat(tree.exteriorNodeIndices);
    tree.allTriangles = tree.insideTriangles.slice();
    return { tree, points };
}

describe('PolygonTreeEx.getContainingTriangle (tree search)', () => {
    const { tree, points } = buildTwoLevelTree();

    it('finds a triangle of the root node for a point in the ring', () => {
        const result = tree.getContainingTriangle(v2(0.5, 0.5), points);
        expect(result.nIndex).toBe(0);
        expect(result.tIndex).toBeGreaterThanOrEqual(0);
        const tri = tree.nodes[0].triangulation[result.tIndex];
        expect(referenceInTriangle(v2(0.5, 0.5), tri, points)).toBe(true);
    });

    it('finds a triangle of the hole node for a point inside the hole', () => {
        const result = tree.getContainingTriangle(v2(2, 2), points);
        expect(result.nIndex).toBe(1);
        const tri = tree.nodes[1].triangulation[result.tIndex];
        expect(referenceInTriangle(v2(2, 2), tri, points)).toBe(true);
    });

    it('returns the INVALID pair for a point outside the tree', () => {
        const result = tree.getContainingTriangle(v2(10, 10), points);
        expect(result.nIndex).toBe(PolygonTreeEx.INVALID);
        expect(result.tIndex).toBe(PolygonTreeEx.INVALID);
        expect(PolygonTreeEx.INVALID).toBe(Number.MAX_SAFE_INTEGER);
    });

    it('locates a triangle for every interior sample point', () => {
        for (let x = 0.25; x < 4; x += 0.5) {
            for (let y = 0.25; y < 4; y += 0.5) {
                const test = v2(x, y);
                const result = tree.getContainingTriangle(test, points);
                expect(result.nIndex).not.toBe(PolygonTreeEx.INVALID);
                const tri = tree.nodes[result.nIndex].triangulation[result.tIndex];
                expect(referenceInTriangle(test, tri, points)).toBe(true);
            }
        }
    });

    it('reports containment on triangle boundaries and vertices', () => {
        // A point on a shared edge is inside both adjacent triangles; the
        // query reports the first one found.
        const onEdge = tree.getContainingTriangle(v2(2, 0), points);
        expect(onEdge.nIndex).toBe(0);
        const atVertex = tree.getContainingTriangle(v2(0, 0), points);
        expect(atVertex.nIndex).toBe(0);
    });
});

describe('PolygonTreeEx.getContainingTriangleInList', () => {
    const { tree, points } = buildTwoLevelTree();

    it('searches the interior triangles', () => {
        const result = tree.getContainingTriangleInList(v2(0.5, 2),
            tree.interiorTriangles, tree.interiorNodeIndices, points);
        expect(result.nIndex).toBe(0);
        expect(referenceInTriangle(v2(0.5, 2),
            tree.interiorTriangles[result.tIndex], points)).toBe(true);
    });

    it('searches the exterior triangles', () => {
        const result = tree.getContainingTriangleInList(v2(2, 2),
            tree.exteriorTriangles, tree.exteriorNodeIndices, points);
        expect(result.nIndex).toBe(1);
        expect(referenceInTriangle(v2(2, 2),
            tree.exteriorTriangles[result.tIndex], points)).toBe(true);
    });

    it('searches the inside triangles and agrees with the tree search', () => {
        for (let x = 0.25; x < 4; x += 0.5) {
            for (let y = 0.25; y < 4; y += 0.5) {
                const test = v2(x, y);
                const byList = tree.getContainingTriangleInList(test,
                    tree.insideTriangles, tree.insideNodeIndices, points);
                const byTree = tree.getContainingTriangle(test, points);
                expect(byList.nIndex).toBe(byTree.nIndex);
            }
        }
    });

    it('returns the INVALID pair for an outside point', () => {
        const result = tree.getContainingTriangleInList(v2(-1, -1),
            tree.insideTriangles, tree.insideNodeIndices, points);
        expect(result.nIndex).toBe(PolygonTreeEx.INVALID);
        expect(result.tIndex).toBe(PolygonTreeEx.INVALID);
    });

    it('throws when the triangle and node-index arrays differ in length', () => {
        expect(() => tree.getContainingTriangleInList(v2(0, 0),
            tree.insideTriangles, [0], points)).toThrow('Invalid argument.');
    });
});

describe('PolygonTreeEx.getContainingTriangleWithChirality', () => {
    const { tree, points } = buildTwoLevelTree();

    it('finds a counterclockwise triangle by index', () => {
        const t = PolygonTreeEx.getContainingTriangleWithChirality(v2(0.5, 0.5),
            tree.interiorTriangles, 1, points);
        expect(t).not.toBe(PolygonTreeEx.INVALID);
        expect(referenceInTriangle(v2(0.5, 0.5), tree.interiorTriangles[t], points))
            .toBe(true);
    });

    it('finds a clockwise triangle by index', () => {
        const t = PolygonTreeEx.getContainingTriangleWithChirality(v2(2, 2),
            tree.exteriorTriangles, -1, points);
        expect(t).not.toBe(PolygonTreeEx.INVALID);
        expect(referenceInTriangle(v2(2, 2), tree.exteriorTriangles[t], points))
            .toBe(true);
    });

    it('finds nothing when the chirality is wrong for the triangles', () => {
        // The interior triangles are counterclockwise; testing them with
        // chirality -1 inverts every half-plane test, so an interior point
        // is rejected.
        const t = PolygonTreeEx.getContainingTriangleWithChirality(v2(0.5, 0.5),
            tree.interiorTriangles, -1, points);
        expect(t).toBe(PolygonTreeEx.INVALID);
    });

    it('returns INVALID for a point outside every triangle', () => {
        const t = PolygonTreeEx.getContainingTriangleWithChirality(v2(100, 100),
            tree.allTriangles, 1, points);
        expect(t).toBe(PolygonTreeEx.INVALID);
    });

    it('agrees with the reference test on a dense grid', () => {
        for (let x = -1; x <= 5; x += 0.25) {
            for (let y = -1; y <= 5; y += 0.25) {
                const test = v2(x, y);
                const t = PolygonTreeEx.getContainingTriangleWithChirality(test,
                    tree.interiorTriangles, 1, points);
                const anyMatch = tree.interiorTriangles.some(tri =>
                    referenceInTriangle(test, tri, points));
                expect(t !== PolygonTreeEx.INVALID).toBe(anyMatch);
            }
        }
    });
});

describe('PolygonTreeEx default construction', () => {
    it('starts with empty node and triangle arrays', () => {
        const tree = new PolygonTreeEx();
        expect(tree.nodes).toEqual([]);
        expect(tree.interiorTriangles).toEqual([]);
        expect(tree.interiorNodeIndices).toEqual([]);
        expect(tree.exteriorTriangles).toEqual([]);
        expect(tree.exteriorNodeIndices).toEqual([]);
        expect(tree.insideTriangles).toEqual([]);
        expect(tree.insideNodeIndices).toEqual([]);
        expect(tree.outsideTriangles).toEqual([]);
        expect(tree.allTriangles).toEqual([]);
    });
});
