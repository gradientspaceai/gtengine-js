import { describe, expect, it } from 'vitest';
import { ConvexHull2 } from '../src/ConvexHull2.js';
import { PolygonTree, PolygonTreeEx } from '../src/PolygonTree.js';
import { TriangulateCDT } from '../src/TriangulateCDT.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';
import { exactDyadic, inCircle2 } from './helpers/exact.js';

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

// ---------------------------------------------------------------------------
// Verification (V13): property-based exact-cover and Delaunay oracles.
//
// Every generator below places the vertices on an integer lattice with small
// coordinates, so twice-signed-area sums, the point-in-polygon crossing test
// and the separating-axis overlap test are all exact in binary64 and the
// comparisons can be exact equalities rather than tolerances.
// ---------------------------------------------------------------------------

// Twice the signed area of the index polygon; positive for counterclockwise.
function twiceSignedAreaPoly(polygon: readonly number[],
    points: readonly Vector[]): number {
    let area = 0;
    const n = polygon.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const a = points[polygon[i0]].values;
        const b = points[polygon[i1]].values;
        area += a[0] * b[1] - b[0] * a[1];
    }
    return area;
}

function twiceTotalArea(triangles: readonly Tri[],
    points: readonly Vector[]): number {
    let sum = 0;
    for (const tri of triangles) {
        sum += Math.abs(signedArea2(tri, points));
    }
    return sum;
}

// Crossing-number point-in-polygon test in exact integer arithmetic. The point
// (px, py) and the polygon coordinates are both multiplied by 'scale', which
// callers set to 3 so that a triangle centroid becomes an integer. The caller
// must guarantee that the point is not on the polygon boundary, which holds
// for the centroid of a triangle that lies in the region: the open triangle is
// contained in the open region.
function insidePolygonExact(px: number, py: number,
    polygon: readonly number[], points: readonly Vector[],
    scale: number): boolean {
    let inside = false;
    const n = polygon.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const ax = scale * points[polygon[i0]].values[0];
        const ay = scale * points[polygon[i0]].values[1];
        const bx = scale * points[polygon[i1]].values[0];
        const by = scale * points[polygon[i1]].values[1];
        if ((ay > py) !== (by > py)) {
            // Sign of (px - x(py)) * (by - ay), where x(py) is the abscissa of
            // the edge at height py. Multiplying by (by - ay) clears the
            // division, so the test stays in exact integers.
            const t = (bx - ax) * (py - ay) - (px - ax) * (by - ay);
            if (by - ay > 0 ? t > 0 : t < 0) {
                inside = !inside;
            }
        }
    }
    return inside;
}

// Three times the centroid of a triangle, an integer pair on a lattice.
function centroid3(tri: Tri, points: readonly Vector[]): [number, number] {
    const a = points[tri[0]].values;
    const b = points[tri[1]].values;
    const c = points[tri[2]].values;
    return [a[0] + b[0] + c[0], a[1] + b[1] + c[1]];
}

// Exact separating-axis test for two triangles. The return value is true when
// the interiors intersect; triangles that only touch along an edge or at a
// vertex are not overlapping.
function trianglesOverlap(t0: Tri, t1: Tri, points: readonly Vector[]): boolean {
    const corners = (t: Tri): number[][] =>
        [points[t[0]].values, points[t[1]].values, points[t[2]].values];
    const c0 = corners(t0), c1 = corners(t1);
    for (const c of [c0, c1]) {
        for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
            const nx = -(c[i1][1] - c[i0][1]);
            const ny = c[i1][0] - c[i0][0];
            if (nx === 0 && ny === 0) {
                continue;   // degenerate edge
            }
            const project = (cs: number[][]): [number, number] => {
                let lo = Infinity, hi = -Infinity;
                for (const p of cs) {
                    const d = nx * p[0] + ny * p[1];
                    lo = Math.min(lo, d);
                    hi = Math.max(hi, d);
                }
                return [lo, hi];
            };
            const [lo0, hi0] = project(c0);
            const [lo1, hi1] = project(c1);
            if (hi0 <= lo1 || hi1 <= lo0) {
                return false;
            }
        }
    }
    return true;
}

function expectPairwiseDisjoint(triangles: readonly Tri[],
    points: readonly Vector[]): void {
    for (let i = 0; i < triangles.length; ++i) {
        for (let j = i + 1; j < triangles.length; ++j) {
            expect(trianglesOverlap(triangles[i], triangles[j], points),
                'triangles ' + i + ' and ' + j + ' overlap').toBe(false);
        }
    }
}

// A lattice rectangle with a strictly interior rectangular hole. Extra
// collinear vertices are placed on the bottom edge of the outer rectangle, so
// the generator also exercises the non-simple-polygon handling.
interface RectWithHole {
    points: Vector[];
    outer: number[];
    hole: number[];
    outerArea: number;
    holeArea: number;
}

const rectWithHole: fc.Arbitrary<RectWithHole> = fc.tuple(
    fc.integer({ min: 4, max: 9 }),          // width
    fc.integer({ min: 4, max: 9 }),          // height
    fc.integer({ min: 0, max: 4095 }),       // hole position and size
    fc.uniqueArray(fc.integer({ min: 1, max: 8 }),
        { minLength: 0, maxLength: 3 })      // extra bottom-edge abscissas
).map(([w, h, seed, extras]) => {
    // A hole [x0, x0+dw] x [y0, y0+dh] strictly inside the rectangle.
    const x0 = 1 + (seed % (w - 2));
    const y0 = 1 + ((seed >> 3) % (h - 2));
    const dw = 1 + ((seed >> 6) % (w - 1 - x0));
    const dh = 1 + ((seed >> 9) % (h - 1 - y0));

    const points: Vector[] = [];
    const push = (x: number, y: number): number => {
        points.push(v2(x, y));
        return points.length - 1;
    };
    const outer: number[] = [];
    outer.push(push(0, 0));
    for (const x of extras.filter(x => x < w).sort((a, b) => a - b)) {
        outer.push(push(x, 0));
    }
    outer.push(push(w, 0));
    outer.push(push(w, h));
    outer.push(push(0, h));
    // The hole is listed clockwise.
    const hole = [
        push(x0, y0), push(x0, y0 + dh),
        push(x0 + dw, y0 + dh), push(x0 + dw, y0)
    ];
    return { points, outer, hole, outerArea: w * h, holeArea: dw * dh };
});

// A convex lattice polygon: the convex hull of random lattice points, with the
// point pool restricted to the hull vertices so that every input point is a
// polygon vertex.
const convexLatticePolygon:
    fc.Arbitrary<{ points: Vector[]; polygon: number[] }> =
    fc.array(latticeVector(2, -8, 8), { minLength: 4, maxLength: 9 })
        .map(pts => {
            const ch = new ConvexHull2();
            ch.compute(pts);
            if (ch.getDimension() !== 2) {
                return { points: [] as Vector[], polygon: [] as number[] };
            }
            const points = ch.getHull().map(i => pts[i]);
            return { points, polygon: points.map((_unused, i) => i) };
        })
        .filter(c => c.polygon.length >= 3);

describe('TriangulateCDT verification', () => {
    it('exactly covers a lattice rectangle with a rectangular hole', () => {
        check(rectWithHole, c => {
            const root = new PolygonTree();
            root.polygon = c.outer.slice();
            const hole = new PolygonTree();
            hole.polygon = c.hole.slice();
            root.child = [hole];

            const tree = new TriangulateCDT().compute(c.points, root);

            // The rectangle is the convex hull of the points, so the Delaunay
            // triangulation covers exactly the rectangle and nothing is left
            // outside the polygon tree.
            expect(tree.outsideTriangles.length).toBe(0);
            expectConstraintEdges(tree);
            expectConsistentClassification(tree, c.points);

            // Exact areas: the arithmetic is integer, so no tolerance.
            expect(twiceTotalArea(tree.interiorTriangles, c.points))
                .toBe(2 * (c.outerArea - c.holeArea));
            expect(twiceTotalArea(tree.exteriorTriangles, c.points))
                .toBe(2 * c.holeArea);
            expect(twiceTotalArea(tree.allTriangles, c.points))
                .toBe(2 * c.outerArea);

            // The triangles tile the rectangle: pairwise disjoint interiors
            // together with the exact area sum is an exact cover.
            expectPairwiseDisjoint(tree.allTriangles, c.points);

            // Each triangle is attributed to the right region.
            for (const tri of tree.interiorTriangles) {
                const [px, py] = centroid3(tri, c.points);
                expect(insidePolygonExact(px, py, c.outer, c.points, 3))
                    .toBe(true);
                expect(insidePolygonExact(px, py, c.hole, c.points, 3))
                    .toBe(false);
            }
            for (const tri of tree.exteriorTriangles) {
                const [px, py] = centroid3(tri, c.points);
                expect(insidePolygonExact(px, py, c.hole, c.points, 3))
                    .toBe(true);
            }
        }, 60);
    }, 30000);

    it('reproduces the Delaunay triangulation of a convex polygon', () => {
        // With no reflex boundary and no hole, every polygon edge is a hull
        // edge and is already present in the unconstrained triangulation, so
        // ConstrainedDelaunay2 performs no retriangulation and the output must
        // satisfy the empty-circumcircle property. The in-circle predicate is
        // evaluated exactly with bigint arithmetic.
        check(convexLatticePolygon, c => {
            const root = new PolygonTree();
            root.polygon = c.polygon.slice();
            const tree = new TriangulateCDT().compute(c.points, root);

            expect(tree.outsideTriangles.length).toBe(0);
            expect(tree.exteriorTriangles.length).toBe(0);
            expect(twiceTotalArea(tree.interiorTriangles, c.points))
                .toBe(twiceSignedAreaPoly(c.polygon, c.points));
            expectConstraintEdges(tree);
            expectConsistentClassification(tree, c.points);
            expectPairwiseDisjoint(tree.allTriangles, c.points);

            const flat: number[] = [];
            for (const p of c.points) { flat.push(p.values[0], p.values[1]); }
            const e = exactDyadic(flat);
            const xy = (i: number): [bigint, bigint] => [e[2 * i], e[2 * i + 1]];
            for (const tri of tree.allTriangles) {
                const a = tri[0];
                let b = tri[1], d = tri[2];
                if (signedArea2([a, b, d], c.points) < 0) {
                    const swap = b; b = d; d = swap;
                }
                const [ax, ay] = xy(a), [bx, by] = xy(b), [cx, cy] = xy(d);
                for (let k = 0; k < c.points.length; ++k) {
                    if (k === a || k === b || k === d) { continue; }
                    const [kx, ky] = xy(k);
                    expect(inCircle2(ax, ay, bx, by, cx, cy, kx, ky),
                        'point ' + k + ' is inside a circumcircle')
                        .toBeLessThanOrEqual(0);
                }
            }
        }, 60);
    }, 30000);

    it('is unaffected by duplicated points in the vertex pool', () => {
        // Appending copies of polygon vertices and referencing the copies must
        // give the same geometry; the triangulator deduplicates the points
        // before triangulating.
        check(fc.tuple(rectWithHole, fc.array(fc.integer({ min: 0, max: 7 }),
            { minLength: 1, maxLength: 4 })), ([c, picks]) => {
            const points = c.points.slice();
            const outer = c.outer.slice();
            for (const pick of picks) {
                const i = pick % outer.length;
                points.push(v2(points[outer[i]].values[0],
                    points[outer[i]].values[1]));
                outer[i] = points.length - 1;
            }

            const build = (poly: readonly number[]): PolygonTree => {
                const root = new PolygonTree();
                root.polygon = poly.slice();
                const hole = new PolygonTree();
                hole.polygon = c.hole.slice();
                root.child = [hole];
                return root;
            };

            const plain = new TriangulateCDT().compute(c.points, build(c.outer));
            const dupes = new TriangulateCDT().compute(points, build(outer));
            expect(twiceTotalArea(dupes.interiorTriangles, points))
                .toBe(twiceTotalArea(plain.interiorTriangles, c.points));
            expect(twiceTotalArea(dupes.exteriorTriangles, points))
                .toBe(twiceTotalArea(plain.exteriorTriangles, c.points));
            expect(dupes.allTriangles.length).toBe(plain.allTriangles.length);
            expectConstraintEdges(dupes);
            expectConsistentClassification(dupes, points);
        }, 40);
    }, 30000);

    it('preserves the upstream remapping quirk for coincident vertices', () => {
        // The hole shares the location of the outer vertex 0, listed a second
        // time as the input point 4. RemapPolygonTree overwrites the remapping
        // entry of the first occurrence with the input index of the duplicate,
        // so the restored outer polygon reports index 4 where the caller
        // passed index 0. The two input points are equal, so the geometry is
        // unaffected; the port preserves the upstream behavior. See the port
        // notes in src/TriangulateCDT.ts.
        const points = [
            v2(0, 0), v2(6, 0), v2(6, 6), v2(0, 6),
            v2(0, 0), v2(4, 1), v2(1, 4)
        ];
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const hole = new PolygonTree();
        hole.polygon = [4, 6, 5];   // clockwise
        root.child = [hole];

        const tree = new TriangulateCDT().compute(points, root);
        expect(tree.nodes[0].polygon).toEqual([4, 1, 2, 3]);
        expect(tree.nodes[1].polygon).toEqual([4, 6, 5]);

        // The triangulation itself is correct: the outer square has area 36
        // and the hole triangle (0,0), (1,4), (4,1) has area 7.5.
        expect(twiceTotalArea(tree.interiorTriangles, points)).toBe(2 * 28.5);
        expect(twiceTotalArea(tree.exteriorTriangles, points)).toBe(2 * 7.5);
        expect(tree.outsideTriangles.length).toBe(0);
        expectConstraintEdges(tree);
        expectConsistentClassification(tree, points);
        expectPairwiseDisjoint(tree.allTriangles, points);
    });
});
