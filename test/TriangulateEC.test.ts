import { describe, it, expect } from 'vitest';
import { TriangulateEC } from '../src/TriangulateEC.js';
import { PolygonTree } from '../src/PolygonTree.js';
import { Vector } from '../src/Vector.js';

function pts(points: [number, number][]): Vector[] {
    return points.map(p => Vector.fromArray(p));
}

// Twice the signed area of the triangle <a,b,c>; positive when the triangle
// is counterclockwise.
function twiceSignedArea(points: Vector[], tri: [number, number, number]): number {
    const a = points[tri[0]].values;
    const b = points[tri[1]].values;
    const c = points[tri[2]].values;
    return (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
}

// Twice the signed area of a closed polygon given by indices into 'points'.
function twicePolygonArea(points: Vector[], polygon: number[]): number {
    let sum = 0;
    const n = polygon.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const p0 = points[polygon[i0]].values;
        const p1 = points[polygon[i1]].values;
        sum += p0[0] * p1[1] - p1[0] * p0[1];
    }
    return sum;
}

// Every triangle must be counterclockwise and nondegenerate, the total area
// must match, and no triangle may be repeated (up to cyclic rotation).
function validate(points: Vector[], triangles: [number, number, number][],
    expectedTwiceArea: number, expectedCount: number): void {
    expect(triangles.length).toBe(expectedCount);

    let total = 0;
    const seen = new Set<string>();
    for (const tri of triangles) {
        const area2 = twiceSignedArea(points, tri);
        expect(area2).toBeGreaterThan(0);
        total += area2;

        // Canonical form: rotate so the smallest index is first. The winding
        // is preserved, so two triangles with the same canonical form are the
        // same triangle.
        let k = 0;
        for (let i = 1; i < 3; ++i) {
            if (tri[i] < tri[k]) {
                k = i;
            }
        }
        const key = `${tri[k]},${tri[(k + 1) % 3]},${tri[(k + 2) % 3]}`;
        expect(seen.has(key)).toBe(false);
        seen.add(key);
    }

    expect(total).toBeCloseTo(expectedTwiceArea, 9);
}

describe('TriangulateEC', () => {
    it('rejects fewer than three points', () => {
        expect(() => new TriangulateEC(pts([[0, 0], [1, 0]]))).toThrow('Invalid input.');
    });

    it('rejects a numPoints larger than the array', () => {
        expect(() => new TriangulateEC(pts([[0, 0], [1, 0], [0, 1]]), 5))
            .toThrow('Invalid input.');
    });

    it('triangulates a triangle into itself', () => {
        const points = pts([[0, 0], [4, 0], [0, 3]]);
        const tri = new TriangulateEC(points);
        tri.triangulate();
        validate(points, tri.getTriangles(), 12, 1);
        expect(tri.getTriangles()[0]).toEqual([0, 1, 2]);
    });

    it('triangulates a convex square as a fan', () => {
        const points = pts([[0, 0], [4, 0], [4, 4], [0, 4]]);
        const tri = new TriangulateEC(points);
        tri.triangulate();
        expect(tri.getTriangles()).toEqual([[0, 1, 2], [0, 2, 3]]);
        validate(points, tri.getTriangles(), 32, 2);
    });

    it('triangulates a convex hexagon as a fan', () => {
        const points = pts([
            [2, 0], [4, 0], [6, 3], [4, 6], [2, 6], [0, 3]
        ]);
        const tri = new TriangulateEC(points);
        tri.triangulate();
        const triangles = tri.getTriangles();
        expect(triangles).toEqual([[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 5]]);
        validate(points, triangles, twicePolygonArea(points, [0, 1, 2, 3, 4, 5]), 4);
    });

    it('triangulates a concave (arrow) polygon', () => {
        const points = pts([[0, 0], [4, 0], [4, 4], [2, 1], [0, 4]]);
        const tri = new TriangulateEC(points);
        tri.triangulate();
        validate(points, tri.getTriangles(),
            twicePolygonArea(points, [0, 1, 2, 3, 4]), 3);
    });

    it('triangulates a concave comb polygon', () => {
        const points = pts([
            [0, 0], [6, 0], [6, 6], [5, 6], [5, 2], [4, 2],
            [4, 6], [3, 6], [3, 2], [2, 2], [2, 6], [1, 6], [1, 2], [0, 2]
        ]);
        const n = points.length;
        const polygon = Array.from({ length: n }, (_, i) => i);
        expect(twicePolygonArea(points, polygon)).toBeGreaterThan(0);

        const tri = new TriangulateEC(points);
        tri.triangulate();
        validate(points, tri.getTriangles(),
            twicePolygonArea(points, polygon), n - 2);
    });

    it('triangulates an explicit subpolygon of a shared point pool', () => {
        // The point pool has more points than the polygon uses.
        const points = pts([
            [0, 0], [10, 0], [10, 10], [0, 10],
            [100, 100], [200, 100], [200, 200]
        ]);
        const polygon = [0, 1, 2, 3];
        const tri = new TriangulateEC(points);
        tri.triangulatePolygon(polygon);
        validate(points, tri.getTriangles(), twicePolygonArea(points, polygon), 2);
        for (const t of tri.getTriangles()) {
            for (const i of t) {
                expect(i).toBeLessThan(4);
            }
        }
    });

    it('triangulates a square with a square hole', () => {
        const points = pts([
            [0, 0], [10, 0], [10, 10], [0, 10],
            [4, 4], [4, 6], [6, 6], [6, 4]
        ]);
        const outer = [0, 1, 2, 3];
        const inner = [4, 5, 6, 7];
        expect(twicePolygonArea(points, outer)).toBe(200);
        expect(twicePolygonArea(points, inner)).toBe(-8);

        const tri = new TriangulateEC(points);
        tri.triangulateWithHole(outer, inner);

        // n + 2h - 2 with n = 8 total vertices and h = 1 hole.
        validate(points, tri.getTriangles(), 200 - 8, 8);
    });

    it('triangulates a rectangle with two holes', () => {
        const points = pts([
            [0, 0], [20, 0], [20, 10], [0, 10],
            [2, 2], [2, 8], [8, 8], [8, 2],
            [12, 3], [12, 7], [17, 7], [17, 3]
        ]);
        const outer = [0, 1, 2, 3];
        const inner0 = [4, 5, 6, 7];
        const inner1 = [8, 9, 10, 11];
        expect(twicePolygonArea(points, inner0)).toBeLessThan(0);
        expect(twicePolygonArea(points, inner1)).toBeLessThan(0);

        const expectedTwiceArea = twicePolygonArea(points, outer)
            + twicePolygonArea(points, inner0) + twicePolygonArea(points, inner1);

        const tri = new TriangulateEC(points);
        tri.triangulateWithHoles(outer, [inner0, inner1]);

        // n + 2h - 2 with n = 12 and h = 2.
        validate(points, tri.getTriangles(), expectedTwiceArea, 12 + 4 - 2);
    });

    it('gives the same result for triangulateWithHole and triangulateWithHoles', () => {
        const points = pts([
            [0, 0], [10, 0], [10, 10], [0, 10],
            [4, 4], [4, 6], [6, 6], [6, 4]
        ]);
        const outer = [0, 1, 2, 3];
        const inner = [4, 5, 6, 7];

        const a = new TriangulateEC(points);
        a.triangulateWithHole(outer, inner);
        const b = new TriangulateEC(points);
        b.triangulateWithHoles(outer, [inner]);
        expect(b.getTriangles()).toEqual(a.getTriangles());
    });

    it('triangulates a polygon tree with a nested outer polygon', () => {
        const points = pts([
            [0, 0], [20, 0], [20, 20], [0, 20],
            [5, 5], [5, 15], [15, 15], [15, 5],
            [8, 8], [12, 8], [12, 12], [8, 12]
        ]);
        const outer = [0, 1, 2, 3];
        const hole = [4, 5, 6, 7];
        const island = [8, 9, 10, 11];
        expect(twicePolygonArea(points, hole)).toBeLessThan(0);
        expect(twicePolygonArea(points, island)).toBeGreaterThan(0);

        const root = new PolygonTree();
        root.polygon = outer;
        const holeNode = new PolygonTree();
        holeNode.polygon = hole;
        const islandNode = new PolygonTree();
        islandNode.polygon = island;
        holeNode.child = [islandNode];
        root.child = [holeNode];

        const tri = new TriangulateEC(points);
        tri.triangulateTree(root);

        // Region 1: the annulus between 'outer' and 'hole' has 4 + 4 + 2 = 10
        // combined vertices, hence 8 triangles. Region 2: the island is a
        // simple quadrilateral, hence 2 triangles.
        const expectedTwiceArea = twicePolygonArea(points, outer)
            + twicePolygonArea(points, hole) + twicePolygonArea(points, island);
        validate(points, tri.getTriangles(), expectedTwiceArea, 10);
    });

    it('triangulates a tree whose root has no children', () => {
        const points = pts([[0, 0], [4, 0], [4, 4], [2, 1], [0, 4]]);
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3, 4];

        const tri = new TriangulateEC(points);
        tri.triangulateTree(root);
        validate(points, tri.getTriangles(),
            twicePolygonArea(points, root.polygon), 3);
    });

    it('triangulates a tree with two sibling holes', () => {
        const points = pts([
            [0, 0], [20, 0], [20, 10], [0, 10],
            [2, 2], [2, 8], [8, 8], [8, 2],
            [12, 3], [12, 7], [17, 7], [17, 3]
        ]);
        const root = new PolygonTree();
        root.polygon = [0, 1, 2, 3];
        const h0 = new PolygonTree();
        h0.polygon = [4, 5, 6, 7];
        const h1 = new PolygonTree();
        h1.polygon = [8, 9, 10, 11];
        root.child = [h0, h1];

        const tri = new TriangulateEC(points);
        tri.triangulateTree(root);

        const expectedTwiceArea = twicePolygonArea(points, root.polygon)
            + twicePolygonArea(points, h0.polygon)
            + twicePolygonArea(points, h1.polygon);
        validate(points, tri.getTriangles(), expectedTwiceArea, 14);
    });

    it('reuses a single object for several queries', () => {
        const points = pts([
            [0, 0], [10, 0], [10, 10], [0, 10],
            [4, 4], [4, 6], [6, 6], [6, 4]
        ]);
        const tri = new TriangulateEC(points);

        tri.triangulatePolygon([0, 1, 2, 3]);
        validate(points, tri.getTriangles(), 200, 2);

        tri.triangulateWithHole([0, 1, 2, 3], [4, 5, 6, 7]);
        validate(points, tri.getTriangles(), 192, 8);

        // The previous result must be cleared, not appended to.
        tri.triangulatePolygon([0, 1, 2, 3]);
        validate(points, tri.getTriangles(), 200, 2);
    });

    it('is invariant to a cyclic rotation of the polygon indices', () => {
        const points = pts([[0, 0], [4, 0], [4, 4], [2, 1], [0, 4]]);
        const expectedTwiceArea = twicePolygonArea(points, [0, 1, 2, 3, 4]);
        for (let shift = 0; shift < 5; ++shift) {
            const polygon = [0, 1, 2, 3, 4].map(i => (i + shift) % 5);
            const tri = new TriangulateEC(points);
            tri.triangulatePolygon(polygon);
            validate(points, tri.getTriangles(), expectedTwiceArea, 3);
        }
    });

    it('triangulates a randomly generated star-shaped polygon', () => {
        // A star-shaped polygon around the origin: the vertices are sorted by
        // angle so the polygon is simple and counterclockwise. The radii
        // alternate to force many reflex vertices.
        const n = 24;
        const raw: [number, number][] = [];
        for (let i = 0; i < n; ++i) {
            const angle = 2 * Math.PI * i / n;
            const radius = (i % 2 === 0 ? 10 : 4);
            raw.push([radius * Math.cos(angle), radius * Math.sin(angle)]);
        }
        const points = pts(raw);
        const polygon = Array.from({ length: n }, (_, i) => i);
        expect(twicePolygonArea(points, polygon)).toBeGreaterThan(0);

        const tri = new TriangulateEC(points);
        tri.triangulate();
        validate(points, tri.getTriangles(),
            twicePolygonArea(points, polygon), n - 2);
    });

    it('accepts a numPoints smaller than the array length', () => {
        const points = pts([[0, 0], [4, 0], [4, 4], [0, 4], [99, 99]]);
        const tri = new TriangulateEC(points, 4);
        tri.triangulate();
        validate(points, tri.getTriangles(), 32, 2);
    });
});
