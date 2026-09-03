import { describe, expect, it } from 'vitest';
import { PolygonTree, PolygonTreeEx } from '../src/PolygonTree';
import { TriangulateCDT } from '../src/TriangulateCDT';
import { Vector } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

type Tri = [number, number, number];

// Twice the signed area of the triangle, positive for counterclockwise
// ordering.
function signedArea2(tri: Tri, points: readonly Vector[]): number {
    const a = points[tri[0]].values;
    const b = points[tri[1]].values;
    const c = points[tri[2]].values;
    return (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
}

function totalArea(triangles: readonly Tri[], points: readonly Vector[]): number {
    let sum = 0;
    for (const tri of triangles) {
        sum += Math.abs(signedArea2(tri, points)) / 2;
    }
    return sum;
}

function centroid(tri: Tri, points: readonly Vector[]): [number, number] {
    const a = points[tri[0]].values;
    const b = points[tri[1]].values;
    const c = points[tri[2]].values;
    return [(a[0] + b[0] + c[0]) / 3, (a[1] + b[1] + c[1]) / 3];
}

// Even-odd point-in-polygon test on the polygon given by index loop.
function pointInPolygon(p: [number, number], polygon: readonly number[],
    points: readonly Vector[]): boolean {
    let inside = false;
    const n = polygon.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const a = points[polygon[i0]].values;
        const b = points[polygon[i1]].values;
        if ((a[1] > p[1]) !== (b[1] > p[1])) {
            const x = a[0] + (p[1] - a[1]) * (b[0] - a[0]) / (b[1] - a[1]);
            if (p[0] < x) {
                inside = !inside;
            }
        }
    }
    return inside;
}

function edgeSet(triangles: readonly Tri[]): Set<string> {
    const edges = new Set<string>();
    for (const tri of triangles) {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            const a = tri[i0], b = tri[i1];
            edges.add(a < b ? `${a},${b}` : `${b},${a}`);
        }
    }
    return edges;
}

// Every edge of every output polygon must be an edge of the triangulation.
function expectConstraintEdges(tree: PolygonTreeEx): void {
    const edges = edgeSet(tree.allTriangles);
    for (const node of tree.nodes) {
        const n = node.polygon.length;
        for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
            const a = node.polygon[i0], b = node.polygon[i1];
            expect(edges.has(a < b ? `${a},${b}` : `${b},${a}`)).toBe(true);
        }
    }
}

// The classified lists must partition allTriangles and the windings must
// match the chirality of the owning node.
function expectConsistentClassification(tree: PolygonTreeEx,
    points: readonly Vector[]): void {
    expect(tree.insideTriangles.length).toBe(
        tree.interiorTriangles.length + tree.exteriorTriangles.length);
    expect(tree.allTriangles.length).toBe(
        tree.insideTriangles.length + tree.outsideTriangles.length);
    expect(tree.insideNodeIndices.length).toBe(tree.insideTriangles.length);
    expect(tree.interiorNodeIndices.length).toBe(tree.interiorTriangles.length);
    expect(tree.exteriorNodeIndices.length).toBe(tree.exteriorTriangles.length);

    // The union of the per-node triangulations is exactly insideTriangles.
    let nodeCount = 0;
    for (const node of tree.nodes) {
        nodeCount += node.triangulation.length;
        for (const tri of node.triangulation) {
            const area = signedArea2(tri, points);
            expect(area !== 0).toBe(true);
            expect(Math.sign(area)).toBe(node.chirality);
        }
    }
    expect(nodeCount).toBe(tree.insideTriangles.length);

    // The triangles of the whole triangulation are distinct.
    const keys = new Set<string>(tree.allTriangles.map(
        tri => [...tri].sort((a, b) => a - b).join(',')));
    expect(keys.size).toBe(tree.allTriangles.length);
}

describe('TriangulateCDT', () => {
    it('triangulates a convex polygon without holes', () => {
        const points = [v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4)];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];

        const tree = new TriangulateCDT().compute(points, root);
        expect(tree.nodes.length).toBe(1);
        expect(tree.nodes[0].chirality).toBe(+1);
        expect(tree.nodes[0].parent).toBe(PolygonTreeEx.INVALID);
        expect(tree.allTriangles.length).toBe(2);
        expect(tree.interiorTriangles.length).toBe(2);
        expect(tree.exteriorTriangles.length).toBe(0);
        expect(tree.outsideTriangles.length).toBe(0);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(16, 12);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('triangulates a nonconvex polygon and reports outside triangles', () => {
        // An L-shaped polygon of area 3. Its convex hull is the pentagon
        // (0,0),(2,0),(2,1),(1,2),(0,2) of area 3.5, so the outside
        // triangles cover the missing area 0.5.
        const points = [
            v2(0, 0), v2(2, 0), v2(2, 1), v2(1, 1), v2(1, 2), v2(0, 2)
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3, 4, 5];

        const tree = new TriangulateCDT().compute(points, root);
        expect(tree.interiorTriangles.length).toBe(4);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(3, 12);
        expect(tree.outsideTriangles.length).toBeGreaterThan(0);
        expect(totalArea(tree.outsideTriangles, points)).toBeCloseTo(0.5, 12);
        expect(totalArea(tree.allTriangles, points)).toBeCloseTo(3.5, 12);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('triangulates a square with a square hole', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4),   // outer, ccw
            v2(1, 1), v2(1, 3), v2(3, 3), v2(3, 1)    // inner, cw
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole = new PolygonTree();
        hole.polygon = [4, 5, 6, 7];
        root.child = [hole];

        const tree = new TriangulateCDT().compute(points, root);

        expect(tree.nodes.length).toBe(2);
        expect(tree.nodes[0].chirality).toBe(+1);
        expect(tree.nodes[1].chirality).toBe(-1);
        expect(tree.nodes[0].minChild).toBe(1);
        expect(tree.nodes[0].supChild).toBe(2);
        expect(tree.nodes[1].parent).toBe(0);

        // The triangulation of 8 points whose convex hull has 4 vertices has
        // 2*8 - 2 - 4 = 10 triangles. The annulus between the two squares
        // uses 8 + 2*1 - 2 = 8 of them and the hole uses the other 2.
        expect(tree.allTriangles.length).toBe(10);
        expect(tree.interiorTriangles.length).toBe(8);
        expect(tree.exteriorTriangles.length).toBe(2);
        expect(tree.insideTriangles.length).toBe(10);
        expect(tree.outsideTriangles.length).toBe(0);

        // The interior triangles tile the region between the outer boundary
        // and the hole; the exterior triangles tile the hole.
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(12, 12);
        expect(totalArea(tree.exteriorTriangles, points)).toBeCloseTo(4, 12);

        // No interior triangle lies inside the hole.
        for (const tri of tree.interiorTriangles) {
            expect(pointInPolygon(centroid(tri, points), [4, 5, 6, 7], points))
                .toBe(false);
        }
        // Every exterior triangle lies inside the hole.
        for (const tri of tree.exteriorTriangles) {
            expect(pointInPolygon(centroid(tri, points), [4, 5, 6, 7], points))
                .toBe(true);
        }

        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('triangulates a square with two holes', () => {
        const points = [
            v2(0, 0), v2(10, 0), v2(10, 6), v2(0, 6),      // outer, ccw
            v2(1, 1), v2(1, 3), v2(3, 3), v2(3, 1),        // hole 0, cw
            v2(6, 2), v2(6, 5), v2(9, 5), v2(9, 2)         // hole 1, cw
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole0 = new PolygonTree();
        hole0.polygon = [4, 5, 6, 7];
        const hole1 = new PolygonTree();
        hole1.polygon = [8, 9, 10, 11];
        root.child = [hole0, hole1];

        const tree = new TriangulateCDT().compute(points, root);

        expect(tree.nodes.length).toBe(3);
        expect(tree.nodes[0].minChild).toBe(1);
        expect(tree.nodes[0].supChild).toBe(3);
        expect(tree.nodes[1].chirality).toBe(-1);
        expect(tree.nodes[2].chirality).toBe(-1);

        // 12 points, convex hull with 4 vertices: 2*12 - 2 - 4 = 18 triangles.
        expect(tree.allTriangles.length).toBe(18);
        // The interior region has 12 vertices and 2 holes: 12 + 4 - 2 = 14.
        expect(tree.interiorTriangles.length).toBe(14);
        expect(tree.exteriorTriangles.length).toBe(4);
        expect(tree.outsideTriangles.length).toBe(0);

        expect(totalArea(tree.interiorTriangles, points))
            .toBeCloseTo(60 - 4 - 9, 12);
        expect(totalArea(tree.exteriorTriangles, points)).toBeCloseTo(13, 12);

        for (const tri of tree.interiorTriangles) {
            const c = centroid(tri, points);
            expect(pointInPolygon(c, [4, 5, 6, 7], points)).toBe(false);
            expect(pointInPolygon(c, [8, 9, 10, 11], points)).toBe(false);
        }

        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('triangulates a nested island inside a hole', () => {
        const points = [
            v2(0, 0), v2(10, 0), v2(10, 10), v2(0, 10),    // outer, ccw
            v2(2, 2), v2(2, 8), v2(8, 8), v2(8, 2),        // hole, cw
            v2(4, 4), v2(6, 4), v2(6, 6), v2(4, 6)         // island, ccw
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole = new PolygonTree();
        hole.polygon = [4, 5, 6, 7];
        const island = new PolygonTree();
        island.polygon = [8, 9, 10, 11];
        hole.child = [island];
        root.child = [hole];

        const tree = new TriangulateCDT().compute(points, root);

        expect(tree.nodes.length).toBe(3);
        expect(tree.nodes[0].chirality).toBe(+1);
        expect(tree.nodes[1].chirality).toBe(-1);
        expect(tree.nodes[2].chirality).toBe(+1);
        expect(tree.nodes[1].parent).toBe(0);
        expect(tree.nodes[2].parent).toBe(1);

        // 12 points, hull with 4 vertices: 18 triangles. The three regions
        // use 8 (outer minus hole), 8 (hole minus island) and 2 (island).
        expect(tree.allTriangles.length).toBe(18);
        expect(tree.nodes[0].triangulation.length).toBe(8);
        expect(tree.nodes[1].triangulation.length).toBe(8);
        expect(tree.nodes[2].triangulation.length).toBe(2);
        expect(tree.interiorTriangles.length).toBe(10);
        expect(tree.exteriorTriangles.length).toBe(8);
        expect(tree.outsideTriangles.length).toBe(0);

        expect(totalArea(tree.nodes[0].triangulation, points))
            .toBeCloseTo(100 - 36, 12);
        expect(totalArea(tree.nodes[1].triangulation, points))
            .toBeCloseTo(36 - 4, 12);
        expect(totalArea(tree.nodes[2].triangulation, points))
            .toBeCloseTo(4, 12);
        expect(totalArea(tree.allTriangles, points)).toBeCloseTo(100, 12);

        // No triangle of the outermost region is inside the hole, and no
        // triangle of the hole region is inside the island.
        for (const tri of tree.nodes[0].triangulation) {
            expect(pointInPolygon(centroid(tri, points), [4, 5, 6, 7], points))
                .toBe(false);
        }
        for (const tri of tree.nodes[1].triangulation) {
            expect(pointInPolygon(centroid(tri, points), [8, 9, 10, 11], points))
                .toBe(false);
            expect(pointInPolygon(centroid(tri, points), [4, 5, 6, 7], points))
                .toBe(true);
        }
        for (const tri of tree.nodes[2].triangulation) {
            expect(pointInPolygon(centroid(tri, points), [8, 9, 10, 11], points))
                .toBe(true);
        }

        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('handles collinear boundary vertices on the outer polygon', () => {
        // (2,0) is an interior point of the segment from (0,0) to (4,0).
        const points = [
            v2(0, 0), v2(2, 0), v2(4, 0), v2(4, 4), v2(0, 4)
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3, 4];

        const tree = new TriangulateCDT().compute(points, root);
        // 5 points, all 5 on the convex hull boundary: 2*5 - 2 - 5 = 3.
        expect(tree.allTriangles.length).toBe(3);
        expect(tree.interiorTriangles.length).toBe(3);
        expect(tree.outsideTriangles.length).toBe(0);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(16, 12);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('handles collinear boundary vertices on a hole', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4),   // outer, ccw
            v2(1, 1), v2(1, 3), v2(3, 3), v2(3, 1),   // inner, cw
            v2(2, 1)                                  // collinear on (3,1)-(1,1)
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole = new PolygonTree();
        hole.polygon = [4, 5, 6, 7, 8];
        root.child = [hole];

        const tree = new TriangulateCDT().compute(points, root);
        // 9 points, hull with 4 vertices: 2*9 - 2 - 4 = 12 triangles.
        expect(tree.allTriangles.length).toBe(12);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(12, 12);
        expect(totalArea(tree.exteriorTriangles, points)).toBeCloseTo(4, 12);
        expect(tree.outsideTriangles.length).toBe(0);
        for (const tri of tree.interiorTriangles) {
            expect(pointInPolygon(centroid(tri, points), [4, 5, 6, 7], points))
                .toBe(false);
        }
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('handles duplicate points in the vertex pool', () => {
        // Index 4 has the same location as index 0 and index 5 the same
        // location as index 2. The polygon references the duplicates.
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4), v2(0, 0), v2(4, 4)
        ];
        const root = new PolygonTree();
        root.polygon = [4, 1, 5, 3];

        const tree = new TriangulateCDT().compute(points, root);
        expect(tree.allTriangles.length).toBe(2);
        expect(tree.interiorTriangles.length).toBe(2);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(16, 12);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('splits an edge that passes through another polygon vertex', () => {
        // The hole vertex (2,0) is an interior point of the outer edge from
        // (0,0) to (4,0). The outer polygon of the output tree gains that
        // vertex, so the output polygon differs from the input polygon.
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4),   // outer, ccw
            v2(2, 0), v2(1, 2), v2(3, 2)              // inner, cw
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole = new PolygonTree();
        hole.polygon = [4, 5, 6];
        root.child = [hole];

        const tree = new TriangulateCDT().compute(points, root);
        expect(tree.nodes[0].polygon.length).toBe(5);
        expect(tree.nodes[0].polygon).toContain(4);
        expect(tree.nodes[1].polygon).toEqual([4, 5, 6]);

        expect(totalArea(tree.exteriorTriangles, points)).toBeCloseTo(2, 12);
        expect(totalArea(tree.interiorTriangles, points)).toBeCloseTo(14, 12);
        expect(tree.outsideTriangles.length).toBe(0);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
    });

    it('rejects invalid input', () => {
        const root = new PolygonTree();
        root.polygon = [0, 1, 2];
        expect(() => new TriangulateCDT().compute(
            [v2(0, 0), v2(1, 0)], root)).toThrow();
        expect(() => new TriangulateCDT().compute(
            [Vector.fromArray([0, 0, 0]), Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0])], root)).toThrow();

        // Fewer than 3 distinct points referenced by the tree.
        const degenerate = new PolygonTree();
        degenerate.polygon = [0, 1, 0];
        expect(() => new TriangulateCDT().compute(
            [v2(0, 0), v2(1, 0), v2(0, 1)], degenerate)).toThrow();
    });

    it('produces the same output when called twice on the same input', () => {
        const points = [
            v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4),
            v2(1, 1), v2(1, 3), v2(3, 3), v2(3, 1)
        ];
        const makeTree = (): PolygonTree => {
            const root = new PolygonTree();
            root.polygon = [0, 1, 2, 3];
            const hole = new PolygonTree();
            hole.polygon = [4, 5, 6, 7];
            root.child = [hole];
            return root;
        };

        const triangulator = new TriangulateCDT();
        const tree0 = triangulator.compute(points, makeTree());
        const tree1 = triangulator.compute(points, makeTree());
        expect(tree1.allTriangles).toEqual(tree0.allTriangles);
        expect(tree1.interiorTriangles).toEqual(tree0.interiorTriangles);
        expect(tree1.exteriorTriangles).toEqual(tree0.exteriorTriangles);
        expect(tree1.insideNodeIndices).toEqual(tree0.insideNodeIndices);
    });

    it('triangulates randomized star polygons with a hole', () => {
        // The outer boundary is a randomized star-shaped polygon around the
        // origin and the hole is a small square at the center. The area of
        // the interior triangles must equal the area of the outer polygon
        // minus the area of the hole.
        let state = 20260902 >>> 0;
        const rand = (): number => {
            state = (state * 1664525 + 1013904223) >>> 0;
            return state / 4294967296;
        };

        for (let trial = 0; trial < 10; ++trial) {
            const n = 8 + Math.floor(rand() * 8);
            const points: Vector[] = [];
            const outer: number[] = [];
            for (let i = 0; i < n; ++i) {
                const angle = 2 * Math.PI * i / n;
                const radius = 4 + 3 * rand();
                points.push(v2(radius * Math.cos(angle),
                    radius * Math.sin(angle)));
                outer.push(i);
            }
            // Shoelace area of the outer polygon (counterclockwise).
            let outerArea2 = 0;
            for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
                const a = points[i0].values, b = points[i1].values;
                outerArea2 += a[0] * b[1] - b[0] * a[1];
            }

            const hole: number[] = [];
            for (const [x, y] of [[-1, -1], [-1, 1], [1, 1], [1, -1]]) {
                hole.push(points.length);
                points.push(v2(x, y));
            }

            const root = new PolygonTree();
            root.polygon = outer;
            const child = new PolygonTree();
            child.polygon = hole;
            root.child = [child];

            const tree = new TriangulateCDT().compute(points, root);
            expect(totalArea(tree.interiorTriangles, points))
                .toBeCloseTo(outerArea2 / 2 - 4, 9);
            expect(totalArea(tree.exteriorTriangles, points))
                .toBeCloseTo(4, 9);
            expect(tree.insideTriangles.length + tree.outsideTriangles.length)
                .toBe(tree.allTriangles.length);
            expectConstraintEdges(tree);
            expectConsistentClassification(tree, points);
        }
    });
});
