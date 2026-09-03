import { describe, it, expect } from 'vitest';
import { BSPPolygon2 } from '../src/BSPPolygon2.js';
import { EdgeKey } from '../src/EdgeKey.js';
import { Vector } from '../src/Vector.js';

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
