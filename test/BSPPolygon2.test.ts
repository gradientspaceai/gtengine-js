import { describe, it, expect } from 'vitest';
import { BSPPolygon2 } from '../src/BSPPolygon2.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { Vector } from '../src/Vector.js';
import { check, fc } from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function vec(x: number, y: number): Vector {
    const v = new Vector(2);
    v.values[0] = x;
    v.values[1] = y;
    return v;
}

// Build a polygon from an ordered loop of vertices. The interior of a BSP
// polygon is on the negative side of its directed edges; for the Classify
// used here (nor = Perp(dir) = (dir.y, -dir.x)), a counterclockwise loop puts
// the interior on the negative side.
function makeLoop(points: Array<[number, number]>, epsilon = 0): BSPPolygon2 {
    const polygon = new BSPPolygon2(epsilon);
    const indices = points.map(p => polygon.insertVertex(vec(p[0], p[1])));
    for (let i = 0; i < indices.length; ++i) {
        const j = (i + 1) % indices.length;
        polygon.insertEdge(new EdgeKey(true, indices[i], indices[j]));
    }
    polygon.finalize();
    return polygon;
}

function square(x0: number, y0: number, x1: number, y1: number,
    epsilon = 0): BSPPolygon2 {
    return makeLoop([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], epsilon);
}

// The total edge length of a polygon, an easily checked invariant of the
// Boolean results.
function totalEdgeLength(polygon: BSPPolygon2): number {
    let sum = 0;
    for (let i = 0; i < polygon.getNumEdges(); ++i) {
        const e = polygon.getEdge(i);
        const p0 = polygon.getVertex(e.V[0]);
        const p1 = polygon.getVertex(e.V[1]);
        sum += Math.hypot(p1.values[0] - p0.values[0], p1.values[1] - p0.values[1]);
    }
    return sum;
}

// The signed area enclosed by the directed edges (the shoelace sum over the
// edge set, which is orientation aware and does not need a loop ordering).
function signedArea(polygon: BSPPolygon2): number {
    let sum = 0;
    for (let i = 0; i < polygon.getNumEdges(); ++i) {
        const e = polygon.getEdge(i);
        const p0 = polygon.getVertex(e.V[0]);
        const p1 = polygon.getVertex(e.V[1]);
        sum += p0.values[0] * p1.values[1] - p1.values[0] * p0.values[1];
    }
    return 0.5 * sum;
}

describe('BSPPolygon2: construction and member access', () => {
    it('assigns unique indices and deduplicates vertices', () => {
        const polygon = new BSPPolygon2(0);
        const i0 = polygon.insertVertex(vec(0, 0));
        const i1 = polygon.insertVertex(vec(1, 0));
        const i2 = polygon.insertVertex(vec(0, 0));
        expect(i0).toBe(0);
        expect(i1).toBe(1);
        expect(i2).toBe(0);
        expect(polygon.getNumVertices()).toBe(2);
        expect(polygon.getVertex(1).values).toEqual([1, 0]);
    });

    it('treats negative zero as zero, as the C++ comparison does', () => {
        const polygon = new BSPPolygon2(0);
        expect(polygon.insertVertex(vec(0, 0))).toBe(0);
        expect(polygon.insertVertex(vec(-0, -0))).toBe(0);
        expect(polygon.getNumVertices()).toBe(1);
    });

    it('copies the inserted vertex', () => {
        const polygon = new BSPPolygon2(0);
        const p = vec(3, 4);
        polygon.insertVertex(p);
        p.values[0] = 99;
        expect(polygon.getVertex(0).values).toEqual([3, 4]);
    });

    it('deduplicates directed edges but distinguishes their direction', () => {
        const polygon = new BSPPolygon2(0);
        polygon.insertVertex(vec(0, 0));
        polygon.insertVertex(vec(1, 0));
        expect(polygon.insertEdge(new EdgeKey(true, 0, 1))).toBe(0);
        expect(polygon.insertEdge(new EdgeKey(true, 0, 1))).toBe(0);
        expect(polygon.insertEdge(new EdgeKey(true, 1, 0))).toBe(1);
        expect(polygon.getNumEdges()).toBe(2);
        expect(polygon.getEdge(1).V).toEqual([1, 0]);
    });

    it('rejects degenerate edges', () => {
        const polygon = new BSPPolygon2(0);
        polygon.insertVertex(vec(0, 0));
        expect(() => polygon.insertEdge(new EdgeKey(true, 0, 0)))
            .toThrow(/Degenerate edges not allowed/);
    });

    it('rejects vertices that are not two dimensional', () => {
        const polygon = new BSPPolygon2(0);
        expect(() => polygon.insertVertex(new Vector(3)))
            .toThrow(/vertices must be 2D/);
    });

    it('clamps a negative epsilon to zero', () => {
        // A negative epsilon would break the classification tests; the
        // constructor clamps it, so a polygon built with -1 behaves as one
        // built with 0.
        const a = square(0, 0, 1, 1, -1);
        expect(a.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(a.pointLocation(vec(2, 2))).toBe(1);
    });

    it('requires the tree before the Boolean operations', () => {
        const polygon = new BSPPolygon2(0);
        polygon.insertVertex(vec(0, 0));
        polygon.insertVertex(vec(1, 0));
        polygon.insertEdge(new EdgeKey(true, 0, 1));
        expect(() => polygon.pointLocation(vec(0, 0))).toThrow(/Tree must exist/);
        expect(() => polygon.negated()).toThrow(/Tree must exist/);
        expect(() => polygon.intersection(polygon)).toThrow(/Tree must exist/);
    });

    it('rejects finalize with no edges', () => {
        const polygon = new BSPPolygon2(0);
        expect(() => polygon.finalize()).toThrow(/Invalid input/);
    });

    it('clone is a deep copy that keeps the tree usable', () => {
        const a = square(0, 0, 2, 2);
        const b = a.clone();
        expect(b.getNumVertices()).toBe(a.getNumVertices());
        expect(b.getNumEdges()).toBe(a.getNumEdges());
        expect(b.pointLocation(vec(1, 1))).toBe(-1);
        // Mutating the copy's vertices does not touch the original.
        b.getVertex(0).values[0] = -50;
        expect(a.getVertex(0).values[0]).toBe(0);
    });
});

describe('BSPPolygon2: point location', () => {
    it('locates points against a counterclockwise square', () => {
        const s = square(0, 0, 2, 2);
        expect(s.pointLocation(vec(1, 1))).toBe(-1);      // inside
        expect(s.pointLocation(vec(0.001, 1))).toBe(-1);
        expect(s.pointLocation(vec(3, 1))).toBe(1);       // outside
        expect(s.pointLocation(vec(-1, -1))).toBe(1);
        expect(s.pointLocation(vec(1, 0))).toBe(0);       // on an edge
        expect(s.pointLocation(vec(0, 0))).toBe(0);       // a corner
        expect(s.pointLocation(vec(2, 2))).toBe(0);
        expect(s.pointLocation(vec(0, 1))).toBe(0);
    });

    it('locates points against a triangle', () => {
        const t = makeLoop([[0, 0], [4, 0], [0, 4]]);
        expect(t.pointLocation(vec(1, 1))).toBe(-1);
        expect(t.pointLocation(vec(3, 3))).toBe(1);
        expect(t.pointLocation(vec(2, 0))).toBe(0);
        expect(t.pointLocation(vec(2, 2))).toBe(0);       // on the hypotenuse
        expect(t.pointLocation(vec(-0.5, 2))).toBe(1);
    });

    it('locates points against a nonconvex L shape', () => {
        // An L-shaped polygon (counterclockwise).
        const l = makeLoop([[0, 0], [3, 0], [3, 1], [1, 1], [1, 3], [0, 3]]);
        expect(l.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(l.pointLocation(vec(2.5, 0.5))).toBe(-1);
        expect(l.pointLocation(vec(0.5, 2.5))).toBe(-1);
        expect(l.pointLocation(vec(2.5, 2.5))).toBe(1);   // in the notch
        expect(l.pointLocation(vec(4, 4))).toBe(1);
        expect(l.pointLocation(vec(1, 0.5))).toBe(-1);   // below the notch
        expect(l.pointLocation(vec(1, 2))).toBe(0);      // on the notch edge
        expect(l.pointLocation(vec(2, 1))).toBe(0);
    });

    it('cross-checks the square against an analytic classification', () => {
        const s = square(-1, -1, 1, 1);
        const rand = makeRandom(20260901);
        for (let trial = 0; trial < 400; ++trial) {
            const x = rand() * 4 - 2;
            const y = rand() * 4 - 2;
            // Skip points near the boundary, where the classification is 0.
            const dx = Math.abs(x) - 1;
            const dy = Math.abs(y) - 1;
            if (Math.abs(dx) < 1e-6 || Math.abs(dy) < 1e-6) {
                continue;
            }
            const inside = (dx < 0 && dy < 0);
            expect(s.pointLocation(vec(x, y))).toBe(inside ? -1 : 1);
        }
    });
});

describe('BSPPolygon2: negation', () => {
    it('reverses every edge and flips the point classification', () => {
        const s = square(0, 0, 2, 2);
        const n = s.negated();
        expect(n.getNumEdges()).toBe(s.getNumEdges());
        expect(n.getNumVertices()).toBe(s.getNumVertices());
        expect(n.pointLocation(vec(1, 1))).toBe(1);
        expect(n.pointLocation(vec(3, 3))).toBe(-1);
        expect(n.pointLocation(vec(1, 0))).toBe(0);
        // The negation reverses the orientation, so the signed area flips.
        expect(signedArea(n)).toBeCloseTo(-signedArea(s), 10);
    });

    it('is an involution on the point classification', () => {
        const t = makeLoop([[0, 0], [4, 0], [0, 4]]);
        const tt = t.negated().negated();
        for (const [x, y] of [[1, 1], [3, 3], [0.5, 0.5], [-1, 2], [2, 0]]) {
            expect(tt.pointLocation(vec(x, y))).toBe(t.pointLocation(vec(x, y)));
        }
    });

    it('does not disturb the original polygon', () => {
        const s = square(0, 0, 2, 2);
        const edgesBefore = s.getNumEdges();
        s.negated();
        expect(s.getNumEdges()).toBe(edgesBefore);
        expect(s.pointLocation(vec(1, 1))).toBe(-1);
    });
});

describe('BSPPolygon2: Boolean operations on axis-aligned squares', () => {
    // Two overlapping squares: [0,2]x[0,2] and [1,3]x[1,3]. The intersection
    // is [1,2]x[1,2] (area 1), the union has area 7 and the difference
    // A - B has area 3.
    const makeA = (): BSPPolygon2 => square(0, 0, 2, 2);
    const makeB = (): BSPPolygon2 => square(1, 1, 3, 3);

    it('computes the intersection', () => {
        const c = makeA().intersection(makeB());
        expect(Math.abs(signedArea(c))).toBeCloseTo(1, 8);
        expect(totalEdgeLength(c)).toBeCloseTo(4, 8);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(-1);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(1);
        expect(c.pointLocation(vec(2.5, 2.5))).toBe(1);
    });

    it('the intersection is commutative in area', () => {
        const c0 = makeA().intersection(makeB());
        const c1 = makeB().intersection(makeA());
        expect(Math.abs(signedArea(c1))).toBeCloseTo(Math.abs(signedArea(c0)), 8);
        expect(c1.pointLocation(vec(1.5, 1.5)))
            .toBe(c0.pointLocation(vec(1.5, 1.5)));
    });

    it('computes the union', () => {
        const c = makeA().union(makeB());
        expect(Math.abs(signedArea(c))).toBeCloseTo(7, 8);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(c.pointLocation(vec(2.5, 2.5))).toBe(-1);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(-1);
        expect(c.pointLocation(vec(2.5, 0.5))).toBe(1);
        expect(c.pointLocation(vec(4, 4))).toBe(1);
    });

    it('computes the difference', () => {
        const c = makeA().difference(makeB());
        expect(Math.abs(signedArea(c))).toBeCloseTo(3, 8);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(1);
        expect(c.pointLocation(vec(2.5, 2.5))).toBe(1);
    });

    it('computes the exclusive or', () => {
        const c = makeA().exclusiveOr(makeB());
        // Symmetric difference: 3 + 3 = 6.
        expect(Math.abs(signedArea(c))).toBeCloseTo(6, 8);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(c.pointLocation(vec(2.5, 2.5))).toBe(-1);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(1);
    });

    it('satisfies |A| + |B| = |A and B| + |A or B|', () => {
        const areaA = Math.abs(signedArea(makeA()));
        const areaB = Math.abs(signedArea(makeB()));
        const areaAnd = Math.abs(signedArea(makeA().intersection(makeB())));
        const areaOr = Math.abs(signedArea(makeA().union(makeB())));
        expect(areaA + areaB).toBeCloseTo(areaAnd + areaOr, 8);
    });

    it('handles disjoint polygons', () => {
        const a = square(0, 0, 1, 1);
        const b = square(5, 5, 6, 6);
        // Upstream's operator& calls Finalize() unconditionally and the
        // BSPTree2 constructor asserts that the edge list is nonempty, so an
        // empty intersection throws rather than producing an empty polygon.
        // The port preserves that behavior.
        expect(() => a.intersection(b)).toThrow(/Invalid input/);
        const u = a.union(b);
        expect(Math.abs(signedArea(u))).toBeCloseTo(2, 8);
        expect(u.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(u.pointLocation(vec(5.5, 5.5))).toBe(-1);
        expect(u.pointLocation(vec(3, 3))).toBe(1);
    });

    it('handles a strictly contained polygon', () => {
        const outer = square(0, 0, 10, 10);
        const inner = square(3, 3, 6, 6);
        const meet = outer.intersection(inner);
        expect(Math.abs(signedArea(meet))).toBeCloseTo(9, 8);
        expect(meet.pointLocation(vec(4, 4))).toBe(-1);
        expect(meet.pointLocation(vec(1, 1))).toBe(1);

        const join = outer.union(inner);
        expect(Math.abs(signedArea(join))).toBeCloseTo(100, 8);

        const diff = outer.difference(inner);
        expect(Math.abs(signedArea(diff))).toBeCloseTo(91, 8);
        expect(diff.pointLocation(vec(1, 1))).toBe(-1);
        expect(diff.pointLocation(vec(4, 4))).toBe(1);
    });

    it('intersects a polygon with itself back to itself', () => {
        const a = square(0, 0, 2, 2);
        const b = square(0, 0, 2, 2);
        const c = a.intersection(b);
        expect(Math.abs(signedArea(c))).toBeCloseTo(4, 8);
        expect(c.pointLocation(vec(1, 1))).toBe(-1);
        expect(c.pointLocation(vec(3, 1))).toBe(1);
    });

    it('subtracting a polygon from itself yields the empty result', () => {
        const a = square(0, 0, 2, 2);
        const b = square(0, 0, 2, 2);
        // Every edge of A is coincident with, and opposite in direction to,
        // an edge of ~B, so both go to the ignored 'coDiff' polygon and the
        // result has no edges; upstream then fails the nonempty assertion in
        // the BSPTree2 constructor.
        expect(() => a.difference(b)).toThrow(/Invalid input/);
    });
});

describe('BSPPolygon2: Boolean operations on non-rectangular polygons', () => {
    it('intersects two triangles that meet only along a shared edge', () => {
        const t0 = makeLoop([[0, 0], [4, 0], [0, 4]]);
        const t1 = makeLoop([[4, 4], [0, 4], [4, 0]]);
        expect(t0.pointLocation(vec(1, 1))).toBe(-1);
        expect(t1.pointLocation(vec(1, 1))).toBe(1);
        expect(t0.pointLocation(vec(3, 3))).toBe(1);
        expect(t1.pointLocation(vec(3, 3))).toBe(-1);
        // The two triangles share only the segment x + y = 4 traversed in
        // opposite directions, so the intersection has no edges and the
        // upstream nonempty assertion fires.
        expect(() => t0.intersection(t1)).toThrow(/Invalid input/);
    });

    it('intersects two overlapping triangles', () => {
        // t0 is the lower-left half of [0,4]^2; t1 is the square [1,5]^2
        // cut to a triangle overlapping t0 in a quadrilateral region.
        const t0 = makeLoop([[0, 0], [4, 0], [0, 4]]);
        const t1 = makeLoop([[1, 1], [5, 1], [1, 5]]);
        const c = t0.intersection(t1);
        // The overlap is {x >= 1, y >= 1, x + y <= 4}, a triangle with legs
        // of length 2, so the area is 2.
        expect(Math.abs(signedArea(c))).toBeCloseTo(2, 8);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(-1);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(1);
        expect(c.pointLocation(vec(3, 3))).toBe(1);
    });

    it('intersects a triangle with a square', () => {
        const t = makeLoop([[0, 0], [4, 0], [0, 4]]);
        const s = square(0, 0, 2, 2);
        const c = t.intersection(s);
        // The square minus the corner triangle above x+y = 4 is the whole
        // square, since the square lies inside the triangle except for the
        // corner (2,2) which is on the hypotenuse.
        expect(Math.abs(signedArea(c))).toBeCloseTo(4, 8);
        expect(c.pointLocation(vec(1, 1))).toBe(-1);
        expect(c.pointLocation(vec(3, 3))).toBe(1);
    });

    it('intersects an L shape with a square', () => {
        const l = makeLoop([[0, 0], [3, 0], [3, 1], [1, 1], [1, 3], [0, 3]]);
        const s = square(0, 0, 2, 2);
        const c = l.intersection(s);
        // The L has area 3*1 + 1*2 = 5; intersecting with [0,2]^2 keeps
        // [0,2]x[0,1] (area 2) and [0,1]x[1,2] (area 1), so 3.
        expect(Math.abs(signedArea(c))).toBeCloseTo(3, 8);
        expect(c.pointLocation(vec(0.5, 0.5))).toBe(-1);
        expect(c.pointLocation(vec(1.5, 0.5))).toBe(-1);
        expect(c.pointLocation(vec(0.5, 1.5))).toBe(-1);
        expect(c.pointLocation(vec(1.5, 1.5))).toBe(1);
    });

    it('unions two triangles into a square', () => {
        const t0 = makeLoop([[0, 0], [2, 0], [2, 2]]);
        const t1 = makeLoop([[0, 0], [2, 2], [0, 2]]);
        const u = t0.union(t1);
        expect(Math.abs(signedArea(u))).toBeCloseTo(4, 8);
        expect(u.pointLocation(vec(1, 0.5))).toBe(-1);
        expect(u.pointLocation(vec(0.5, 1))).toBe(-1);
        expect(u.pointLocation(vec(3, 1))).toBe(1);
    });
});

describe('BSPPolygon2: splitEdge', () => {
    it('replaces an edge by two edges through the midpoint', () => {
        const polygon = new BSPPolygon2(0);
        const i0 = polygon.insertVertex(vec(0, 0));
        const i1 = polygon.insertVertex(vec(2, 0));
        polygon.insertEdge(new EdgeKey(true, i0, i1));
        const imid = polygon.insertVertex(vec(1, 0));

        polygon.splitEdge(i0, i1, imid);
        expect(polygon.getNumEdges()).toBe(2);
        expect(polygon.getEdge(0).V).toEqual([i0, imid]);
        expect(polygon.getEdge(1).V).toEqual([imid, i1]);
    });

    it('rejects splitting an edge that does not exist', () => {
        const polygon = new BSPPolygon2(0);
        polygon.insertVertex(vec(0, 0));
        polygon.insertVertex(vec(1, 0));
        polygon.insertVertex(vec(2, 0));
        expect(() => polygon.splitEdge(0, 1, 2)).toThrow(/Edge does not exist/);
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). The BSP tree, the edge/vertex
// maps and the Boolean operations are cross-checked against independent
// computations: an analytic point-in-rectangle classification, the shoelace
// area of the edge set, and the set identities the Boolean operators must
// satisfy.
// ---------------------------------------------------------------------------

describe('BSPPolygon2 verification', () => {
    // Axis-aligned rectangles on an integer grid, which keeps every
    // intersection point exactly representable and avoids the epsilon paths.
    interface Rect { x0: number; y0: number; x1: number; y1: number }

    const rect: fc.Arbitrary<Rect> =
        fc.tuple(fc.integer({ min: -4, max: 3 }), fc.integer({ min: -4, max: 3 }),
            fc.integer({ min: 1, max: 4 }), fc.integer({ min: 1, max: 4 }))
            .map(([x0, y0, w, h]) => ({ x0, y0, x1: x0 + w, y1: y0 + h }));

    const rectPolygon = (r: Rect): BSPPolygon2 =>
        square(r.x0, r.y0, r.x1, r.y1);

    const rectArea = (r: Rect): number => (r.x1 - r.x0) * (r.y1 - r.y0);

    // The analytic classification of a point against a rectangle: -1 inside,
    // 0 on the boundary, +1 outside.
    function classify(r: Rect, x: number, y: number): number {
        const onX = (x === r.x0 || x === r.x1) && r.y0 <= y && y <= r.y1;
        const onY = (y === r.y0 || y === r.y1) && r.x0 <= x && x <= r.x1;
        if (onX || onY) { return 0; }
        return (r.x0 < x && x < r.x1 && r.y0 < y && y < r.y1) ? -1 : 1;
    }

    // Sample points on a half-integer grid covering the rectangle and a
    // margin around it.
    function samples(r: Rect): Array<[number, number]> {
        const points: Array<[number, number]> = [];
        for (let x = r.x0 - 1; x <= r.x1 + 1; x += 0.5) {
            for (let y = r.y0 - 1; y <= r.y1 + 1; y += 0.5) {
                points.push([x, y]);
            }
        }
        return points;
    }

    it('locates every sample point exactly as the analytic test does', () => {
        check(rect, r => {
            const polygon = rectPolygon(r);
            for (const [x, y] of samples(r)) {
                expect(polygon.pointLocation(vec(x, y))).toBe(classify(r, x, y));
            }
        }, 40);
    });

    it('negation reverses every edge and flips the interior', () => {
        check(rect, r => {
            const polygon = rectPolygon(r);
            const negated = polygon.negated();

            // The reversed edge set is exactly the original one reversed.
            const forward = new Set<string>();
            for (let i = 0; i < polygon.getNumEdges(); ++i) {
                const e = polygon.getEdge(i);
                forward.add(`${e.V[0]},${e.V[1]}`);
            }
            expect(negated.getNumEdges()).toBe(polygon.getNumEdges());
            for (let i = 0; i < negated.getNumEdges(); ++i) {
                const e = negated.getEdge(i);
                expect(forward.has(`${e.V[1]},${e.V[0]}`)).toBe(true);
            }

            // The signed areas are opposite and the interior/exterior swap.
            expect(signedArea(negated)).toBeCloseTo(-signedArea(polygon), 10);
            for (const [x, y] of samples(r)) {
                const inside = classify(r, x, y);
                expect(negated.pointLocation(vec(x, y)))
                    .toBe(inside === 0 ? 0 : -inside);
            }

            // Negation is an involution on the classification, and it does
            // not disturb the original polygon.
            const twice = negated.negated();
            for (const [x, y] of samples(r)) {
                expect(twice.pointLocation(vec(x, y))).toBe(classify(r, x, y));
                expect(polygon.pointLocation(vec(x, y))).toBe(classify(r, x, y));
            }
        }, 30);
    });

    it('the Boolean operations reproduce the set operations on rectangles',
        () => {
            // Two rectangles whose overlap has positive area, so that no
            // Boolean result is the empty polygon (which upstream cannot
            // represent: Finalize asserts a nonempty edge list).
            const overlapping = fc.tuple(rect, rect).filter(([a, b]) =>
                Math.max(a.x0, b.x0) < Math.min(a.x1, b.x1)
                && Math.max(a.y0, b.y0) < Math.min(a.y1, b.y1)
                // Neither contains the other, so the difference is nonempty
                // in both directions.
                && !(a.x0 <= b.x0 && b.x1 <= a.x1 && a.y0 <= b.y0 && b.y1 <= a.y1)
                && !(b.x0 <= a.x0 && a.x1 <= b.x1 && b.y0 <= a.y0 && a.y1 <= b.y1));

            check(overlapping, ([ra, rb]) => {
                const overlap = (Math.min(ra.x1, rb.x1) - Math.max(ra.x0, rb.x0))
                    * (Math.min(ra.y1, rb.y1) - Math.max(ra.y0, rb.y0));
                const areaA = rectArea(ra);
                const areaB = rectArea(rb);

                const inter = rectPolygon(ra).intersection(rectPolygon(rb));
                const uni = rectPolygon(ra).union(rectPolygon(rb));
                const diff = rectPolygon(ra).difference(rectPolygon(rb));
                const xor = rectPolygon(ra).exclusiveOr(rectPolygon(rb));

                expect(Math.abs(signedArea(inter))).toBeCloseTo(overlap, 8);
                expect(Math.abs(signedArea(uni)))
                    .toBeCloseTo(areaA + areaB - overlap, 8);
                expect(Math.abs(signedArea(diff)))
                    .toBeCloseTo(areaA - overlap, 8);
                expect(Math.abs(signedArea(xor)))
                    .toBeCloseTo(areaA + areaB - 2 * overlap, 8);

                // Inclusion-exclusion, from the computed areas alone.
                expect(Math.abs(signedArea(inter)) + Math.abs(signedArea(uni)))
                    .toBeCloseTo(areaA + areaB, 8);

                // Point classification of the intersection agrees with the
                // conjunction of the two analytic tests, away from the
                // boundaries where both polygons disagree by construction.
                for (const [x, y] of samples(ra)) {
                    const ca = classify(ra, x, y);
                    const cb = classify(rb, x, y);
                    if (ca === 0 || cb === 0) { continue; }
                    const expected = (ca < 0 && cb < 0) ? -1 : 1;
                    expect(inter.pointLocation(vec(x, y))).toBe(expected);
                    expect(uni.pointLocation(vec(x, y)))
                        .toBe((ca < 0 || cb < 0) ? -1 : 1);
                    expect(diff.pointLocation(vec(x, y)))
                        .toBe((ca < 0 && cb > 0) ? -1 : 1);
                }
            }, 20);
        });

    it('the vertex and edge maps deduplicate exactly', () => {
        check(fc.array(fc.tuple(fc.integer({ min: -3, max: 3 }),
            fc.integer({ min: -3, max: 3 })), { minLength: 1, maxLength: 24 }),
            points => {
                const polygon = new BSPPolygon2(0);
                const seen = new Map<string, number>();
                for (const [x, y] of points) {
                    const key = `${x},${y}`;
                    const index = polygon.insertVertex(vec(x, y));
                    if (seen.has(key)) {
                        expect(index).toBe(seen.get(key));
                    } else {
                        expect(index).toBe(seen.size);
                        seen.set(key, index);
                    }
                    // The stored vertex is the one that was inserted.
                    expect(polygon.getVertex(index).values[0]).toBe(x);
                    expect(polygon.getVertex(index).values[1]).toBe(y);
                }
                expect(polygon.getNumVertices()).toBe(seen.size);

                // Directed edges: <i,j> and <j,i> are distinct keys.
                const edges = new Map<string, number>();
                for (let i = 0; i + 1 < seen.size; ++i) {
                    for (const [a, b] of [[i, i + 1], [i + 1, i]]) {
                        const key = `${a},${b}`;
                        const index = polygon.insertEdge(new EdgeKey(true, a, b));
                        if (edges.has(key)) {
                            expect(index).toBe(edges.get(key));
                        } else {
                            expect(index).toBe(edges.size);
                            edges.set(key, index);
                        }
                    }
                }
                expect(polygon.getNumEdges()).toBe(edges.size);
            }, 40);
    });

    // Regression for a port defect: insertEdge stored the caller's EdgeKey
    // rather than a copy, so mutating that object after insertion changed
    // the array entry while the map key stayed behind (upstream stores the
    // edge by value in both the map and the vector). splitEdge then mutates
    // array entries in place, which makes the aliasing observable.
    it('copies the inserted edge, as the C++ value semantics require', () => {
        const polygon = new BSPPolygon2(0);
        const i0 = polygon.insertVertex(vec(0, 0));
        const i1 = polygon.insertVertex(vec(1, 0));
        const i2 = polygon.insertVertex(vec(0, 1));

        const edge = new EdgeKey(true, i0, i1);
        const index = polygon.insertEdge(edge);
        edge.V[1] = i2;

        // The stored edge is unchanged, so the array and the map still agree.
        expect(polygon.getEdge(index).V[0]).toBe(i0);
        expect(polygon.getEdge(index).V[1]).toBe(i1);
        // Re-inserting the original key must find the same index rather than
        // appending a duplicate.
        expect(polygon.insertEdge(new EdgeKey(true, i0, i1))).toBe(index);
        expect(polygon.getNumEdges()).toBe(1);
    });

    it('clone is a deep copy in the C++ sense', () => {
        check(rect, r => {
            const polygon = rectPolygon(r);
            const copy = polygon.clone();

            // Mutating the copy's vertices and edges must not disturb the
            // original.
            copy.getVertex(0).values[0] += 100;
            copy.getEdge(0).V[0] = copy.getEdge(0).V[1];
            for (const [x, y] of samples(r)) {
                expect(polygon.pointLocation(vec(x, y))).toBe(classify(r, x, y));
            }

            // The clone keeps a usable tree of its own.
            const fresh = rectPolygon(r).clone();
            for (const [x, y] of samples(r)) {
                expect(fresh.pointLocation(vec(x, y))).toBe(classify(r, x, y));
            }
        }, 25);
    });

    it('classifies points of a convex loop against the half-plane test', () => {
        // Regular n-gons, counterclockwise, so the interior is the
        // intersection of the negative half planes of the directed edges.
        check(fc.tuple(fc.integer({ min: 3, max: 9 }),
            fc.integer({ min: 1, max: 5 })), ([n, radius]) => {
                const points: Array<[number, number]> = [];
                for (let i = 0; i < n; ++i) {
                    const angle = (2 * Math.PI * i) / n;
                    points.push([radius * Math.cos(angle),
                        radius * Math.sin(angle)]);
                }
                const polygon = makeLoop(points);

                const random = makeRandom(1234 + n * 31 + radius);
                for (let k = 0; k < 60; ++k) {
                    const x = (random() * 2 - 1) * (radius + 1);
                    const y = (random() * 2 - 1) * (radius + 1);
                    // Independent test: the point is inside when it is
                    // strictly on the negative side of every directed edge.
                    let minSide = Number.POSITIVE_INFINITY;
                    let maxSide = Number.NEGATIVE_INFINITY;
                    for (let i = 0; i < n; ++i) {
                        const p = points[i];
                        const q = points[(i + 1) % n];
                        const dx = q[0] - p[0], dy = q[1] - p[1];
                        // nor = Perp(dir) = (dir.y, -dir.x)
                        const side = dy * (x - p[0]) - dx * (y - p[1]);
                        minSide = Math.min(minSide, side);
                        maxSide = Math.max(maxSide, side);
                    }
                    // Skip points close to an edge line, where the exact
                    // classification is decided by round-off.
                    if (Math.min(Math.abs(minSide), Math.abs(maxSide)) < 1e-6) {
                        continue;
                    }
                    expect(polygon.pointLocation(vec(x, y)))
                        .toBe(maxSide < 0 ? -1 : 1);
                }
            }, 20);
    });
});
