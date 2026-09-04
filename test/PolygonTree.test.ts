import { describe, it, expect } from 'vitest';
import {
    PolygonTree, PolygonTreeEx, PolygonTreeExNode
} from '../src/PolygonTree.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A counterclockwise lattice triangle with a strictly positive area, so the
// sign tests below are exact in binary64.
const ccwLatticeTriangle = fc.tuple(latticeVector(2, -6, 6),
    latticeVector(2, -6, 6), latticeVector(2, -6, 6))
    .map(([a, b, c]) => [a, b, c] as [Vector, Vector, Vector])
    .filter(([a, b, c]) => {
        const area = (b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
            - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0]);
        return area > 0;
    });

// A grid of unit squares over [0,gridN] x [0,gridN], split into two
// counterclockwise triangles each. Coordinates are integers and the query
// points are multiples of 1/2, so every predicate below is exact.
const gridN = 4;
const gridPoints: Vector[] = (() => {
    const p: Vector[] = [];
    for (let j = 0; j <= gridN; ++j) {
        for (let i = 0; i <= gridN; ++i) { p.push(v2(i, j)); }
    }
    return p;
})();
const gridTriangles: [number, number, number][] = (() => {
    const t: [number, number, number][] = [];
    const at = (i: number, j: number) => i + (gridN + 1) * j;
    for (let j = 0; j < gridN; ++j) {
        for (let i = 0; i < gridN; ++i) {
            t.push([at(i, j), at(i + 1, j), at(i + 1, j + 1)]);
            t.push([at(i, j), at(i + 1, j + 1), at(i, j + 1)]);
        }
    }
    return t;
})();

// A point with half-integer coordinates in [-1, gridN + 1], so roughly a
// quarter of the draws land outside the grid.
const halfIntegerPoint = fc.tuple(
    fc.integer({ min: -2, max: 2 * gridN + 2 }),
    fc.integer({ min: -2, max: 2 * gridN + 2 }))
    .map(([a, b]) => v2(a / 2, b / 2));

// Split the grid triangles into a root node and a chain of child nodes; every
// triangle belongs to exactly one node, so the containment answer is unique
// apart from shared edges.
function buildPartitionTree(split: number[]): PolygonTreeEx {
    const tree = new PolygonTreeEx();
    const groups: [number, number, number][][] = [];
    let start = 0;
    for (const s of split.concat([gridTriangles.length])) {
        const end = Math.min(Math.max(s, start), gridTriangles.length);
        groups.push(gridTriangles.slice(start, end));
        start = end;
    }
    tree.nodes = groups.map((g, index) => {
        const node = new PolygonTreeExNode();
        node.chirality = 1;
        node.self = index;
        node.parent = index === 0 ? 0 : index - 1;
        node.minChild = index + 1 < groups.length ? index + 1 : 0;
        node.supChild = index + 1 < groups.length ? index + 2 : 0;
        node.triangulation = g;
        return node;
    });
    for (let index = 0; index < groups.length; ++index) {
        for (const t of groups[index]) {
            tree.insideTriangles.push(t);
            tree.insideNodeIndices.push(index);
        }
    }
    tree.allTriangles = tree.insideTriangles.slice();
    return tree;
}

describe('PolygonTreeEx verification', () => {
    it('matches an exact sign test for a counterclockwise triangle', () => {
        check(fc.tuple(ccwLatticeTriangle, latticeVector(2, -8, 8)),
            ([[a, b, c], test]) => {
                const points = [a, b, c];
                const tri: [number, number, number] = [0, 1, 2];
                const found = PolygonTreeEx.getContainingTriangleWithChirality(
                    test, [tri], 1, points);
                const inside = referenceInTriangle(test, tri, points);
                expect(found === 0).toBe(inside);
                expect(found === PolygonTreeEx.INVALID).toBe(!inside);
            });
    });

    it('finds nothing for a nondegenerate triangle queried with the wrong chirality', () => {
        // With the chirality reversed the test demands that the point be on
        // the outer side of all three edges. The barycentric coordinates sum
        // to one, so no point satisfies that for a triangle of positive area.
        check(fc.tuple(ccwLatticeTriangle, latticeVector(2, -8, 8)),
            ([[a, b, c], test]) => {
                expect(PolygonTreeEx.getContainingTriangleWithChirality(
                    test, [[0, 1, 2]], -1, [a, b, c]))
                    .toBe(PolygonTreeEx.INVALID);
            });
    });

    it('classifies clockwise triangles with chirality -1', () => {
        check(fc.tuple(ccwLatticeTriangle, latticeVector(2, -8, 8)),
            ([[a, b, c], test]) => {
                // Reversing the vertex order makes the triangle clockwise.
                const points = [a, c, b];
                const tri: [number, number, number] = [0, 1, 2];
                const found = PolygonTreeEx.getContainingTriangleWithChirality(
                    test, [tri], -1, points);
                expect(found === 0).toBe(referenceInTriangle(test, tri, points));
            });
    });

    it('reports the first containing triangle of a list, as the linear scan does', () => {
        check(fc.tuple(fc.shuffledSubarray(gridTriangles,
            { minLength: 4, maxLength: gridTriangles.length }),
        halfIntegerPoint), ([triangles, test]) => {
            const found = PolygonTreeEx.getContainingTriangleWithChirality(
                test, triangles, 1, gridPoints);
            let expected = PolygonTreeEx.INVALID;
            for (let t = 0; t < triangles.length; ++t) {
                if (referenceInTriangle(test, triangles[t], gridPoints)) {
                    expected = t;
                    break;
                }
            }
            expect(found).toBe(expected);
        });
    });

    it('locates a point of a partitioned tree in the node that owns it', () => {
        check(fc.tuple(fc.array(fc.integer({ min: 0, max: gridTriangles.length }),
            { minLength: 1, maxLength: 3 }).map(a => a.slice().sort((x, y) => x - y)),
        halfIntegerPoint), ([split, test]) => {
            const tree = buildPartitionTree(split);
            const { nIndex, tIndex } = tree.getContainingTriangle(test, gridPoints);
            const anyContains = gridTriangles.some(
                t => referenceInTriangle(test, t, gridPoints));
            if (!anyContains) {
                expect(nIndex).toBe(PolygonTreeEx.INVALID);
                expect(tIndex).toBe(PolygonTreeEx.INVALID);
                return;
            }
            expect(nIndex).toBeLessThan(tree.nodes.length);
            const triangle = tree.nodes[nIndex].triangulation[tIndex];
            expect(triangle).toBeDefined();
            expect(referenceInTriangle(test, triangle, gridPoints)).toBe(true);
        });
    });

    it('agrees with the list search over the same triangles', () => {
        check(fc.tuple(fc.array(fc.integer({ min: 0, max: gridTriangles.length }),
            { minLength: 1, maxLength: 3 }).map(a => a.slice().sort((x, y) => x - y)),
        halfIntegerPoint), ([split, test]) => {
            const tree = buildPartitionTree(split);
            const inList = tree.getContainingTriangleInList(test,
                tree.insideTriangles, tree.insideNodeIndices, gridPoints);
            const anyContains = gridTriangles.some(
                t => referenceInTriangle(test, t, gridPoints));
            expect(inList.nIndex === PolygonTreeEx.INVALID).toBe(!anyContains);
            if (anyContains) {
                expect(referenceInTriangle(test,
                    tree.insideTriangles[inList.tIndex], gridPoints)).toBe(true);
                expect(tree.insideNodeIndices[inList.tIndex])
                    .toBe(inList.nIndex);
            }
        });
    });

    it('prefers the node it pops first, reproducing the std::stack order', () => {
        // The root is examined first; among the children the last pushed (the
        // largest index) is popped first. Overlapping triangulations pin that
        // order, which the port must reproduce with Array.pop().
        check(fc.tuple(fc.integer({ min: 0, max: gridTriangles.length - 1 }),
            fc.integer({ min: 2, max: 4 })), ([which, numChildren]) => {
            const triangle = gridTriangles[which];
            const inside = v2(
                (gridPoints[triangle[0]].values[0]
                    + gridPoints[triangle[1]].values[0]
                    + gridPoints[triangle[2]].values[0]) / 3,
                (gridPoints[triangle[0]].values[1]
                    + gridPoints[triangle[1]].values[1]
                    + gridPoints[triangle[2]].values[1]) / 3);

            const withRoot = new PolygonTreeEx();
            const root = new PolygonTreeExNode();
            root.chirality = 1;
            root.triangulation = [triangle];
            root.minChild = 1;
            root.supChild = 1 + numChildren;
            withRoot.nodes = [root];
            for (let c = 0; c < numChildren; ++c) {
                const child = new PolygonTreeExNode();
                child.chirality = 1;
                child.self = 1 + c;
                child.triangulation = [triangle];
                withRoot.nodes.push(child);
            }
            // The root owns a containing triangle, so it wins.
            expect(withRoot.getContainingTriangle(inside, gridPoints))
                .toEqual({ nIndex: 0, tIndex: 0 });

            // With an empty root the largest child index is popped first.
            withRoot.nodes[0].triangulation = [];
            expect(withRoot.getContainingTriangle(inside, gridPoints))
                .toEqual({ nIndex: numChildren, tIndex: 0 });
        });
    });

    it('treats edges and vertices of a triangle as inside', () => {
        // Upstream rejects only when sign*(n.d) is strictly positive, so a
        // point on an edge line inside the segment belongs to the triangle.
        check(ccwLatticeTriangle, ([a, b, c]) => {
            const points = [a, b, c];
            const tri: [number, number, number] = [0, 1, 2];
            for (const p of [a, b, c]) {
                expect(PolygonTreeEx.getContainingTriangleWithChirality(
                    p, [tri], 1, points)).toBe(0);
            }
            for (const [p, q] of [[a, b], [b, c], [c, a]]) {
                const mid = v2(0.5 * (p.values[0] + q.values[0]),
                    0.5 * (p.values[1] + q.values[1]));
                expect(PolygonTreeEx.getContainingTriangleWithChirality(
                    mid, [tri], 1, points)).toBe(0);
            }
            const centroid = v2(
                (a.values[0] + b.values[0] + c.values[0]) / 3,
                (a.values[1] + b.values[1] + c.values[1]) / 3);
            expect(PolygonTreeEx.getContainingTriangleWithChirality(
                centroid, [tri], 1, points)).toBe(0);
        });
    });

    it('throws when the triangle and node-index arrays have different lengths', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }),
            fc.integer({ min: 0, max: 5 })), ([numTriangles, numIndices]) => {
            if (numTriangles === numIndices) { return; }
            const tree = buildPartitionTree([gridTriangles.length]);
            const triangles = gridTriangles.slice(0, numTriangles);
            const nodeIndices = new Array<number>(numIndices).fill(0);
            expect(() => tree.getContainingTriangleInList(v2(0.5, 0.5),
                triangles, nodeIndices, gridPoints)).toThrow('Invalid argument.');
        });
    });
});
