import { describe, it, expect } from 'vitest';
import { TriangulateEC } from '../src/TriangulateEC.js';
import { PolygonTree } from '../src/PolygonTree.js';
import { Vector } from '../src/Vector.js';
import { check, fc, scaled, seededRandom } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// TriangulateEC was ported number-only: upstream is
// template <InputType, ComputeType> and expects an exact ComputeType
// (BSNumber/BSRational), but PrimalQuery2 itself is number-only in this port,
// so the sign predicates are ordinary double arithmetic. The generators below
// therefore build configurations in general position (star-shaped polygons
// with distinct radii at well-separated angles), where the predicates are
// unambiguous, rather than lattice configurations full of collinear triples.
//
// The oracle is measure-theoretic and independent of the algorithm: the
// triangles must be positively oriented, their areas must sum to the area of
// the region, and every sample point must be covered by exactly one triangle
// if it is inside the region and by none if it is outside. Together with the
// area identity that rules out both overlaps and gaps.
// ---------------------------------------------------------------------------

const SLOTS = 64;

// Points of a star-shaped polygon: strictly increasing angles about 'center'
// with all angular gaps below pi, so the polygon is simple and
// counterclockwise.
function starPoints(slots: number[], radii: number[], center: [number, number],
    phase: number): [number, number][] {
    return slots.map((slot, k) => {
        const angle = phase + (2 * Math.PI * slot) / SLOTS;
        const radius = radii[k % radii.length];
        return [center[0] + radius * Math.cos(angle),
            center[1] + radius * Math.sin(angle)] as [number, number];
    });
}

const slotArbitrary = (minLength: number, maxLength: number) =>
    fc.uniqueArray(fc.integer({ min: 0, max: SLOTS - 1 }), { minLength, maxLength })
        .map(slots => [...slots].sort((a, b) => a - b))
        .filter(slots => {
            for (let i = 0; i < slots.length; ++i) {
                const next = (i + 1 < slots.length) ? slots[i + 1] : slots[0] + SLOTS;
                if (next - slots[i] >= SLOTS / 2) { return false; }
            }
            return true;
        });

// A simple counterclockwise polygon with reflex vertices (radii vary).
const simplePolygon = fc.record({
    slots: slotArbitrary(3, 14),
    radii: fc.array(fc.integer({ min: 3, max: 10 }), { minLength: 14, maxLength: 14 }),
    phase: scaled(0, 2 * Math.PI, 16)
}).map(({ slots, radii, phase }) =>
    pts(starPoints(slots, radii.map(r => r / 2), [0, 0], phase)));

// The distance from the origin to the closest edge of a polygon: the radius of
// the largest disk about the origin that fits inside a star-shaped polygon.
// The generators use it to guarantee that a hole really is strictly inside its
// container, which the minimum vertex radius alone does not.
function minEdgeDistanceToOrigin(points: [number, number][]): number {
    let best = Number.POSITIVE_INFINITY;
    const n = points.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const a = points[i0], b = points[i1];
        const dx = b[0] - a[0], dy = b[1] - a[1];
        const lengthSquared = dx * dx + dy * dy;
        let t = lengthSquared > 0 ? -(a[0] * dx + a[1] * dy) / lengthSquared : 0;
        t = Math.max(0, Math.min(1, t));
        const x = a[0] + t * dx, y = a[1] + t * dy;
        best = Math.min(best, Math.hypot(x, y));
    }
    return best;
}

function maxRadius(points: [number, number][]): number {
    let best = 0;
    for (const p of points) { best = Math.max(best, Math.hypot(p[0], p[1])); }
    return best;
}

// Crossing-number point-in-polygon test on an index list.
function pointInPolygon(points: Vector[], polygon: number[],
    x: number, y: number): boolean {
    let inside = false;
    const n = polygon.length;
    for (let i0 = n - 1, i1 = 0; i1 < n; i0 = i1++) {
        const p0 = points[polygon[i0]].values;
        const p1 = points[polygon[i1]].values;
        if ((p0[1] > y) !== (p1[1] > y)) {
            const t = (y - p0[1]) / (p1[1] - p0[1]);
            if (x < p0[0] + t * (p1[0] - p0[0])) { inside = !inside; }
        }
    }
    return inside;
}

// The barycentric-sign containment test, returning the smallest of the three
// edge functions scaled by the triangle size so the caller can skip samples
// that are too close to an edge to classify reliably.
function triangleContainment(points: Vector[], tri: [number, number, number],
    x: number, y: number): number {
    const a = points[tri[0]].values, b = points[tri[1]].values, c = points[tri[2]].values;
    const e0 = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]);
    const e1 = (c[0] - b[0]) * (y - b[1]) - (c[1] - b[1]) * (x - b[0]);
    const e2 = (a[0] - c[0]) * (y - c[1]) - (a[1] - c[1]) * (x - c[0]);
    return Math.min(e0, e1, e2);
}

// Every sample point that is safely inside or outside the region must be
// covered by exactly one triangle or by none.
function expectExactCover(points: Vector[], triangles: [number, number, number][],
    inRegion: (x: number, y: number) => boolean, radius: number,
    samples: number): void {
    const random = seededRandom(0x7e51 + triangles.length);
    for (let s = 0; s < samples; ++s) {
        const x = radius * (2 * random() - 1);
        const y = radius * (2 * random() - 1);
        let covered = 0;
        let ambiguous = false;
        for (const tri of triangles) {
            const value = triangleContainment(points, tri, x, y);
            if (Math.abs(value) < 1e-6) { ambiguous = true; break; }
            if (value > 0) { ++covered; }
        }
        if (ambiguous) { continue; }
        expect(covered).toBe(inRegion(x, y) ? 1 : 0);
    }
}

describe('TriangulateEC verification', () => {
    it('triangulates a simple polygon into an exact cover', () => {
        check(simplePolygon, (points) => {
            const n = points.length;
            const polygon = points.map((_, i) => i);
            const query = new TriangulateEC(points);
            query.triangulate();
            const triangles = query.getTriangles();
            validate(points, triangles, twicePolygonArea(points, polygon), n - 2);
            expectExactCover(points, triangles,
                (x, y) => pointInPolygon(points, polygon, x, y), 6, 40);
        }, 60);
    });

    it('triangulate() agrees with triangulatePolygon([0..n-1])', () => {
        check(simplePolygon, (points) => {
            const all = new TriangulateEC(points);
            all.triangulate();
            const explicit = new TriangulateEC(points);
            explicit.triangulatePolygon(points.map((_, i) => i));
            expect(explicit.getTriangles()).toEqual(all.getTriangles());
        }, 60);
    });

    it('a convex polygon is triangulated as the fan from vertex 0', () => {
        check(fc.record({
            slots: slotArbitrary(3, 10),
            radius: scaled(2, 6, 8),
            phase: scaled(0, 2 * Math.PI, 16)
        }), ({ slots, radius, phase }) => {
            // Equal radii put the vertices in convex position.
            const points = pts(starPoints(slots, [radius], [0, 0], phase));
            const query = new TriangulateEC(points);
            query.triangulate();
            const expected: [number, number, number][] = [];
            for (let i = 1; i + 1 < points.length; ++i) {
                expected.push([0, i, i + 1]);
            }
            expect(query.getTriangles()).toEqual(expected);
        }, 60);
    });

    it('triangulates a polygon with one hole', () => {
        const withHole = fc.record({
            outerSlots: slotArbitrary(3, 10),
            outerRadii: fc.array(fc.integer({ min: 14, max: 20 }),
                { minLength: 10, maxLength: 10 }),
            innerSlots: slotArbitrary(3, 8),
            innerRadii: fc.array(fc.integer({ min: 2, max: 6 }),
                { minLength: 8, maxLength: 8 }),
            phase: scaled(0, 2 * Math.PI, 16),
            innerPhase: scaled(0, 2 * Math.PI, 16)
        }).map(({ outerSlots, outerRadii, innerSlots, innerRadii, phase, innerPhase }) => ({
            outerPoints: starPoints(outerSlots, outerRadii.map(r => r / 2), [0, 0], phase),
            // The inner polygon is clockwise, as the algorithm requires.
            innerPoints: starPoints(innerSlots, innerRadii.map(r => r / 2),
                [0, 0], innerPhase).reverse()
        })).filter(({ outerPoints, innerPoints }) =>
            // The hole must be strictly inside the outer polygon; a thin outer
            // polygon can pass close to the origin even when all of its
            // vertices are far away.
            minEdgeDistanceToOrigin(outerPoints) > 1.2 * maxRadius(innerPoints));

        check(withHole, ({ outerPoints, innerPoints }) => {
            const points = pts(outerPoints.concat(innerPoints));
            const outer = outerPoints.map((_, i) => i);
            const inner = innerPoints.map((_, i) => outerPoints.length + i);

            const query = new TriangulateEC(points);
            query.triangulateWithHole(outer, inner);
            const triangles = query.getTriangles();
            // The bridge duplicates two vertices, so the pseudosimple polygon
            // has |outer| + |inner| + 2 vertices and |outer| + |inner|
            // triangles.
            validate(points, triangles,
                twicePolygonArea(points, outer) + twicePolygonArea(points, inner),
                outer.length + inner.length);
            expectExactCover(points, triangles,
                (x, y) => pointInPolygon(points, outer, x, y)
                    && !pointInPolygon(points, inner, x, y), 11, 40);
        }, 50);
    });

    it('triangulates a polygon with several holes', () => {
        const withHoles = fc.record({
            outerSlots: slotArbitrary(3, 10),
            outerRadii: fc.array(fc.integer({ min: 40, max: 52 }),
                { minLength: 10, maxLength: 10 }),
            numHoles: fc.integer({ min: 2, max: 3 }),
            holeSlots: fc.array(slotArbitrary(3, 6), { minLength: 3, maxLength: 3 }),
            holeRadii: fc.array(fc.integer({ min: 2, max: 4 }),
                { minLength: 6, maxLength: 6 }),
            phase: scaled(0, 2 * Math.PI, 16)
        }).map(({ outerSlots, outerRadii, numHoles, holeSlots, holeRadii, phase }) => {
            const outerPoints = starPoints(outerSlots, outerRadii.map(r => r / 2),
                [0, 0], phase);
            const holes: Array<[number, number][]> = [];
            for (let h = 0; h < numHoles; ++h) {
                // Hole centers on a circle of radius 5; each hole has radius
                // at most 2, so the holes are pairwise disjoint.
                const angle = phase + (2 * Math.PI * h) / numHoles;
                const center: [number, number] =
                    [5 * Math.cos(angle), 5 * Math.sin(angle)];
                holes.push(starPoints(holeSlots[h], holeRadii.map(r => r / 2),
                    center, phase + 0.3 * h).reverse());
            }
            return { outerPoints, holes };
        }).filter(({ outerPoints, holes }) =>
            // Every hole must be strictly inside the outer polygon.
            minEdgeDistanceToOrigin(outerPoints)
                > 1.2 * Math.max(...holes.map(maxRadius)));

        check(withHoles, ({ outerPoints, holes }) => {
            const allPoints: [number, number][] = [...outerPoints];
            const outer = outerPoints.map((_, i) => i);
            const inners: number[][] = [];
            for (const holePoints of holes) {
                inners.push(holePoints.map((_, i) => allPoints.length + i));
                allPoints.push(...holePoints);
            }
            const points = pts(allPoints);

            const query = new TriangulateEC(points);
            query.triangulateWithHoles(outer, inners);
            const triangles = query.getTriangles();
            let expectedTwiceArea = twicePolygonArea(points, outer);
            let expectedCount = outer.length - 2;
            for (const inner of inners) {
                expectedTwiceArea += twicePolygonArea(points, inner);
                expectedCount += inner.length + 2;
            }
            validate(points, triangles, expectedTwiceArea, expectedCount);
            expectExactCover(points, triangles, (x, y) =>
                pointInPolygon(points, outer, x, y)
                && !inners.some(inner => pointInPolygon(points, inner, x, y)),
            14, 40);
        }, 40);
    });

    it('triangulates a polygon tree with a nested island', () => {
        const treeArbitrary = fc.record({
            outerSlots: slotArbitrary(3, 8),
            holeSlots: slotArbitrary(3, 8),
            islandSlots: slotArbitrary(3, 6),
            phase: scaled(0, 2 * Math.PI, 16)
        }).map(({ outerSlots, holeSlots, islandSlots, phase }) => ({
            outerPoints: starPoints(outerSlots, [9], [0, 0], phase),
            holePoints: starPoints(holeSlots, [6], [0, 0], phase + 0.11).reverse(),
            islandPoints: starPoints(islandSlots, [3], [0, 0], phase + 0.23)
        })).filter(({ outerPoints, holePoints, islandPoints }) =>
            // The hole must be strictly inside the outer polygon and the
            // island strictly inside the hole.
            minEdgeDistanceToOrigin(outerPoints) > 1.2 * maxRadius(holePoints)
            && minEdgeDistanceToOrigin(holePoints) > 1.2 * maxRadius(islandPoints));

        check(treeArbitrary, ({ outerPoints, holePoints, islandPoints }) => {
            const points = pts([...outerPoints, ...holePoints, ...islandPoints]);
            const outer = outerPoints.map((_, i) => i);
            const hole = holePoints.map((_, i) => outerPoints.length + i);
            const island = islandPoints.map((_, i) =>
                outerPoints.length + holePoints.length + i);

            const root = new PolygonTree();
            root.polygon = outer;
            const holeNode = new PolygonTree();
            holeNode.polygon = hole;
            const islandNode = new PolygonTree();
            islandNode.polygon = island;
            holeNode.child.push(islandNode);
            root.child.push(holeNode);

            const query = new TriangulateEC(points);
            query.triangulateTree(root);
            const triangles = query.getTriangles();
            validate(points, triangles,
                twicePolygonArea(points, outer) + twicePolygonArea(points, hole)
                + twicePolygonArea(points, island),
                (outer.length + hole.length) + (island.length - 2));
            expectExactCover(points, triangles, (x, y) => {
                const inOuter = pointInPolygon(points, outer, x, y);
                const inHole = pointInPolygon(points, hole, x, y);
                const inIsland = pointInPolygon(points, island, x, y);
                return (inOuter && !inHole) || inIsland;
            }, 10, 40);
        }, 40);
    });

    it('a reused query object gives the same answers as fresh ones', () => {
        // The ear-clipping vertex list is a member of TriangulateEC, so a
        // stale list would corrupt the second query. Upstream clears it at the
        // start of DoEarClipping; this pins that behavior. It also pins the
        // port's array semantics: each query builds a fresh triangle array, so
        // an array obtained from an earlier getTriangles() keeps its contents
        // (upstream returns a reference to the member, whose contents change).
        check(fc.tuple(simplePolygon, simplePolygon), ([pointsA, pointsB]) => {
            const combined = [...pointsA, ...pointsB.map(p =>
                Vector.fromArray([p.values[0] + 40, p.values[1]]))];
            const polygonA = pointsA.map((_, i) => i);
            const polygonB = pointsB.map((_, i) => pointsA.length + i);

            const shared = new TriangulateEC(combined);
            shared.triangulatePolygon(polygonA);
            const firstA = shared.getTriangles();
            shared.triangulatePolygon(polygonB);
            const firstB = shared.getTriangles();
            shared.triangulatePolygon(polygonA);
            const secondA = shared.getTriangles();

            const freshA = new TriangulateEC(combined);
            freshA.triangulatePolygon(polygonA);
            const freshB = new TriangulateEC(combined);
            freshB.triangulatePolygon(polygonB);

            expect(firstA).toEqual(freshA.getTriangles());
            expect(firstB).toEqual(freshB.getTriangles());
            expect(secondA).toEqual(freshA.getTriangles());
        }, 40);
    });

    it('is invariant under a cyclic rotation of the polygon indices', () => {
        check(fc.tuple(simplePolygon, fc.integer({ min: 0, max: 13 })),
            ([points, shift]) => {
                const n = points.length;
                const polygon = points.map((_, i) => (i + shift) % n);
                const query = new TriangulateEC(points);
                query.triangulatePolygon(polygon);
                const triangles = query.getTriangles();
                validate(points, triangles,
                    twicePolygonArea(points, polygon), n - 2);
                expectExactCover(points, triangles,
                    (x, y) => pointInPolygon(points, polygon, x, y), 6, 30);
            }, 60);
    });

    it('rejects invalid construction arguments', () => {
        check(fc.integer({ min: 0, max: 5 }), (numPoints) => {
            const points = pts([[0, 0], [4, 0], [4, 4], [0, 4]]);
            if (numPoints >= 3 && numPoints <= points.length) {
                expect(() => new TriangulateEC(points, numPoints)).not.toThrow();
            } else {
                expect(() => new TriangulateEC(points, numPoints))
                    .toThrow('Invalid input.');
            }
        });
    });
});
