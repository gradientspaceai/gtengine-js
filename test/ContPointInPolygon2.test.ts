import { describe, it, expect } from 'vitest';
import { PointInPolygon2 } from '../src/ContPointInPolygon2.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Independent reference implementation of the even-odd (crossing number)
// rule, written in the classic PNPOLY form. It is used as a cross-check for
// points that are not on or near the polygon boundary, where the even-odd
// answer is unambiguous.
function referenceInside(polygon: readonly Vector[], p: Vector): boolean {
    const n = polygon.length;
    const x = p.values[0];
    const y = p.values[1];
    let inside = false;
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const xi = polygon[i].values[0], yi = polygon[i].values[1];
        const xj = polygon[j].values[0], yj = polygon[j].values[1];
        if ((yi > y) !== (yj > y)
            && x < (xj - xi) * (y - yi) / (yj - yi) + xi) {
            inside = !inside;
        }
    }
    return inside;
}

// Distance from p to the nearest polygon edge, used to keep randomized
// cross-checks away from boundary ties.
function distanceToBoundary(polygon: readonly Vector[], p: Vector): number {
    const n = polygon.length;
    let best = Number.MAX_VALUE;
    for (let i = 0, j = n - 1; i < n; j = i++) {
        const ax = polygon[j].values[0], ay = polygon[j].values[1];
        const bx = polygon[i].values[0], by = polygon[i].values[1];
        const ex = bx - ax, ey = by - ay;
        const lenSqr = ex * ex + ey * ey;
        let t = lenSqr > 0 ? ((p.values[0] - ax) * ex + (p.values[1] - ay) * ey) / lenSqr : 0;
        t = Math.min(1, Math.max(0, t));
        const dx = p.values[0] - (ax + t * ex);
        const dy = p.values[1] - (ay + t * ey);
        best = Math.min(best, Math.hypot(dx, dy));
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (state * 1103515245 + 12345) & 0x7fffffff;
        return state / 0x7fffffff;
    };
}

// A counterclockwise unit square.
const square = [v(0, 0), v(1, 0), v(1, 1), v(0, 1)];

// A counterclockwise regular hexagon of circumradius 1.
const hexagon: Vector[] = [];
for (let i = 0; i < 6; ++i) {
    const angle = (2 * Math.PI * i) / 6;
    hexagon.push(v(Math.cos(angle), Math.sin(angle)));
}

// A counterclockwise L-shape (concave).
const lShape = [
    v(0, 0), v(3, 0), v(3, 1), v(1, 1), v(1, 3), v(0, 3)
];

describe('PointInPolygon2 construction', () => {
    it('reports the number of vertices', () => {
        expect(new PointInPolygon2(square).numPoints).toBe(4);
    });

    it('requires at least three vertices', () => {
        expect(() => new PointInPolygon2([v(0, 0), v(1, 1)]))
            .toThrow('PointInPolygon2: at least 3 vertices are required.');
    });

    it('requires 2D vertices', () => {
        expect(() => new PointInPolygon2(
            [v(0, 0), v(1, 0), Vector.fromArray([0, 1, 0])]))
            .toThrow('PointInPolygon2: the vertices must be 2D.');
    });

    it('requires a 2D query point', () => {
        const query = new PointInPolygon2(square);
        expect(() => query.contains(Vector.fromArray([0.5, 0.5, 0])))
            .toThrow('PointInPolygon2: the point must be 2D.');
    });
});

describe('PointInPolygon2.contains (simple polygons)', () => {
    it('classifies interior and exterior points of a square', () => {
        const query = new PointInPolygon2(square);
        expect(query.contains(v(0.5, 0.5))).toBe(true);
        expect(query.contains(v(0.01, 0.99))).toBe(true);
        expect(query.contains(v(-0.5, 0.5))).toBe(false);
        expect(query.contains(v(1.5, 0.5))).toBe(false);
        expect(query.contains(v(0.5, -0.5))).toBe(false);
        expect(query.contains(v(0.5, 1.5))).toBe(false);
    });

    it('gives the same answers for a clockwise square (orientation independent)', () => {
        const cw = [...square].reverse();
        const query = new PointInPolygon2(cw);
        expect(query.contains(v(0.5, 0.5))).toBe(true);
        expect(query.contains(v(2, 2))).toBe(false);
    });

    it('uses a half-open boundary rule (crossing counting)', () => {
        // The ray-counting test is exact on the boundary but asymmetric:
        // the bottom and left edges of the unit square count as inside
        // while the top and right edges do not. These are the values the
        // upstream algorithm produces.
        const query = new PointInPolygon2(square);
        expect(query.contains(v(0.5, 0))).toBe(true);   // bottom edge
        expect(query.contains(v(0, 0.5))).toBe(true);   // left edge
        expect(query.contains(v(0.5, 1))).toBe(false);  // top edge
        expect(query.contains(v(1, 0.5))).toBe(false);  // right edge
        expect(query.contains(v(0, 0))).toBe(true);     // min vertex
        expect(query.contains(v(1, 1))).toBe(false);    // max vertex
    });

    it('handles a concave L-shape, including the reflex notch', () => {
        const query = new PointInPolygon2(lShape);
        expect(query.contains(v(0.5, 0.5))).toBe(true);   // in the corner
        expect(query.contains(v(2.5, 0.5))).toBe(true);   // in the arm
        expect(query.contains(v(0.5, 2.5))).toBe(true);   // in the leg
        expect(query.contains(v(2, 2))).toBe(false);      // in the notch
        expect(query.contains(v(2.5, 2.5))).toBe(false);  // in the notch
        expect(query.contains(v(-1, -1))).toBe(false);
    });

    it('applies the even-odd rule to a self-intersecting bowtie', () => {
        // The two lobes of a bowtie are traversed with opposite
        // orientations; even-odd counting reports both as inside and the
        // crossing point region as the shared vertex only.
        const bowtie = [v(0, 0), v(2, 2), v(2, 0), v(0, 2)];
        const query = new PointInPolygon2(bowtie);
        expect(query.contains(v(0.3, 1))).toBe(true);   // left lobe
        expect(query.contains(v(1.7, 1))).toBe(true);   // right lobe
        expect(query.contains(v(1, 0.3))).toBe(false);  // below the crossing
        expect(query.contains(v(1, 1.7))).toBe(false);  // above the crossing
        for (const p of [v(0.3, 1), v(1.7, 1), v(1, 0.3), v(1, 1.7)]) {
            expect(query.contains(p)).toBe(referenceInside(bowtie, p));
        }
    });

    it('applies the even-odd rule to a polygon that touches itself', () => {
        // A square with a slit cut in from the top edge down to the center;
        // the two slit walls touch along x = 1.
        const slit = [
            v(0, 0), v(2, 0), v(2, 2), v(1, 2), v(1, 1), v(1, 2), v(0, 2)
        ];
        const query = new PointInPolygon2(slit);
        expect(query.contains(v(0.5, 0.5))).toBe(true);
        expect(query.contains(v(1.5, 0.5))).toBe(true);
        expect(query.contains(v(0.5, 1.5))).toBe(true);
        expect(query.contains(v(1.5, 1.5))).toBe(true);
        expect(query.contains(v(2.5, 1))).toBe(false);
    });

    it('handles a degenerate triangle with collinear vertices', () => {
        const degenerate = [v(0, 0), v(1, 0), v(2, 0)];
        const query = new PointInPolygon2(degenerate);
        expect(query.contains(v(1, 0.5))).toBe(false);
        expect(query.contains(v(1, -0.5))).toBe(false);
    });

    it('agrees with an independent even-odd implementation (randomized)', () => {
        const rand = makeRandom(9001);
        const polygons = [square, hexagon, lShape,
            [v(0, 0), v(4, 1), v(2, 2), v(4, 3), v(0, 4), v(1, 2)]];
        for (const polygon of polygons) {
            const query = new PointInPolygon2(polygon);
            let tested = 0;
            for (let i = 0; i < 2000; ++i) {
                const p = v(-2 + 8 * rand(), -2 + 8 * rand());
                if (distanceToBoundary(polygon, p) < 1e-6) {
                    continue;
                }
                expect(query.contains(p)).toBe(referenceInside(polygon, p));
                ++tested;
            }
            expect(tested).toBeGreaterThan(1000);
        }
    });
});

describe('PointInPolygon2 convex queries', () => {
    it('classifies points of a counterclockwise square', () => {
        const query = new PointInPolygon2(square);
        for (const p of [v(0.5, 0.5), v(0.1, 0.9)]) {
            expect(query.containsConvexOrderN(p)).toBe(true);
            expect(query.containsConvexOrderLogN(p)).toBe(true);
            expect(query.containsQuadrilateral(p)).toBe(true);
        }
        for (const p of [v(-0.1, 0.5), v(1.1, 0.5), v(0.5, -0.1), v(0.5, 1.1)]) {
            expect(query.containsConvexOrderN(p)).toBe(false);
            expect(query.containsConvexOrderLogN(p)).toBe(false);
            expect(query.containsQuadrilateral(p)).toBe(false);
        }
    });

    it('treats the boundary as inside for the convex queries', () => {
        const query = new PointInPolygon2(square);
        for (const p of [v(0.5, 0), v(0, 0.5), v(0.5, 1), v(1, 0.5),
            v(0, 0), v(1, 1)]) {
            expect(query.containsConvexOrderN(p)).toBe(true);
            expect(query.containsConvexOrderLogN(p)).toBe(true);
            expect(query.containsQuadrilateral(p)).toBe(true);
        }
    });

    it('works for a triangle (the smallest allowed polygon)', () => {
        const triangle = [v(0, 0), v(1, 0), v(0, 1)];
        const query = new PointInPolygon2(triangle);
        expect(query.containsConvexOrderN(v(0.2, 0.2))).toBe(true);
        expect(query.containsConvexOrderLogN(v(0.2, 0.2))).toBe(true);
        expect(query.containsConvexOrderN(v(0.6, 0.6))).toBe(false);
        expect(query.containsConvexOrderLogN(v(0.6, 0.6))).toBe(false);
    });

    it('containsQuadrilateral returns false unless there are four vertices', () => {
        expect(new PointInPolygon2(hexagon).containsQuadrilateral(v(0, 0)))
            .toBe(false);
        expect(new PointInPolygon2([v(0, 0), v(1, 0), v(0, 1)])
            .containsQuadrilateral(v(0.1, 0.1))).toBe(false);
    });

    it('requires counterclockwise input: a clockwise square rejects its interior', () => {
        const cw = [...square].reverse();
        const query = new PointInPolygon2(cw);
        expect(query.containsConvexOrderN(v(0.5, 0.5))).toBe(false);
    });

    it('the O(N), O(log N) and quadrilateral queries agree (randomized)', () => {
        const rand = makeRandom(31337);
        const quad = [v(0, 0), v(3, 1), v(2, 4), v(-1, 2)];
        const convexPolygons = [square, hexagon, quad,
            [v(0, 0), v(2, 0), v(3, 2), v(1, 3), v(-1, 2)]];

        for (const polygon of convexPolygons) {
            const query = new PointInPolygon2(polygon);
            for (let i = 0; i < 1500; ++i) {
                const p = v(-3 + 9 * rand(), -3 + 9 * rand());
                if (distanceToBoundary(polygon, p) < 1e-6) {
                    continue;
                }
                const orderN = query.containsConvexOrderN(p);
                expect(query.containsConvexOrderLogN(p)).toBe(orderN);
                expect(orderN).toBe(referenceInside(polygon, p));
                if (polygon.length === 4) {
                    expect(query.containsQuadrilateral(p)).toBe(orderN);
                }
            }
        }
    });

    it('agrees with the simple-polygon query away from the boundary (randomized)', () => {
        const rand = makeRandom(777);
        const query = new PointInPolygon2(hexagon);
        for (let i = 0; i < 2000; ++i) {
            const p = v(-2 + 4 * rand(), -2 + 4 * rand());
            if (distanceToBoundary(hexagon, p) < 1e-6) {
                continue;
            }
            expect(query.containsConvexOrderN(p)).toBe(query.contains(p));
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks of the port
// against the upstream ContPointInPolygon2.h semantics.
// ---------------------------------------------------------------------------

// Exact CCW convex hull (monotone chain) of small integer points. Collinear
// points are dropped, so the hull is strictly convex, which is what the
// O(log N) bisection query assumes. Coordinates are integers in [-8, 8], so
// every cross product below is exact in binary64.
function convexHullCCW(points: readonly Vector[]): Vector[] {
    const pts = points.map(p => [p.get(0), p.get(1)] as [number, number])
        .sort((a, b) => (a[0] - b[0]) || (a[1] - b[1]))
        .filter((p, i, a) => i === 0 || p[0] !== a[i - 1][0] || p[1] !== a[i - 1][1]);
    if (pts.length < 3) {
        return [];
    }
    const cross = (o: number[], a: number[], b: number[]): number =>
        (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    const lower: number[][] = [];
    for (const p of pts) {
        while (lower.length >= 2
            && cross(lower[lower.length - 2], lower[lower.length - 1], p) <= 0) {
            lower.pop();
        }
        lower.push(p);
    }
    const upper: number[][] = [];
    for (let i = pts.length - 1; i >= 0; --i) {
        const p = pts[i];
        while (upper.length >= 2
            && cross(upper[upper.length - 2], upper[upper.length - 1], p) <= 0) {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    const hull = lower.concat(upper);
    return hull.length >= 3 ? hull.map(p => Vector.fromArray(p)) : [];
}

// Exact side-of-edge classification of p against a CCW convex polygon:
// +1 strictly inside, 0 on the boundary, -1 strictly outside.
function exactConvexClassify(hull: readonly Vector[], p: Vector): number {
    let onBoundary = false;
    for (let i = 0, j = hull.length - 1; i < hull.length; j = i++) {
        const ax = hull[j].get(0), ay = hull[j].get(1);
        const bx = hull[i].get(0), by = hull[i].get(1);
        const s = (bx - ax) * (p.get(1) - ay) - (by - ay) * (p.get(0) - ax);
        if (s < 0) {
            return -1;
        }
        if (s === 0) {
            onBoundary = true;
        }
    }
    return onBoundary ? 0 : 1;
}

describe('ContPointInPolygon2 verification', () => {
    const hullArb = fc.array(latticeVector(2, -8, 8),
        { minLength: 3, maxLength: 12 })
        .map(convexHullCCW)
        .filter(h => h.length >= 3);

    // All four queries must agree with the exact side-of-edge test on a
    // strictly convex lattice polygon. Boundary points are excluded because
    // the even-odd rule of Contains is deliberately half-open there.
    it('all queries agree with an exact convex classification', () => {
        check(fc.tuple(hullArb, latticeVector(2, -10, 10)),
            ([hull, p]: [Vector[], Vector]) => {
                const cls = exactConvexClassify(hull, p);
                if (cls === 0) {
                    return;  // on the boundary: the rules disagree by design
                }
                const inside = cls > 0;
                const query = new PointInPolygon2(hull);
                expect(query.contains(p)).toBe(inside);
                expect(query.containsConvexOrderN(p)).toBe(inside);
                expect(query.containsConvexOrderLogN(p)).toBe(inside);
                if (hull.length === 4) {
                    expect(query.containsQuadrilateral(p)).toBe(inside);
                }
            });
    });

    // The convex queries include the boundary (all side tests use <= 0), so a
    // vertex and an edge midpoint are reported inside.
    it('convex queries include the polygon boundary', () => {
        check(hullArb, (hull: Vector[]) => {
            const query = new PointInPolygon2(hull);
            for (let i = 0; i < hull.length; ++i) {
                expect(query.containsConvexOrderN(hull[i])).toBe(true);
                expect(query.containsConvexOrderLogN(hull[i])).toBe(true);
                const j = (i + 1) % hull.length;
                const mid = Vector.fromArray([
                    0.5 * (hull[i].get(0) + hull[j].get(0)),
                    0.5 * (hull[i].get(1) + hull[j].get(1))]);
                expect(query.containsConvexOrderN(mid)).toBe(true);
                expect(query.containsConvexOrderLogN(mid)).toBe(true);
            }
        });
    });

    // The even-odd ray-casting rule does not depend on the vertex ordering,
    // nor on the starting vertex of the (cyclic) list.
    it('Contains is invariant under reversal and rotation of the vertices', () => {
        check(fc.tuple(hullArb, latticeVector(2, -10, 10)),
            ([hull, p]: [Vector[], Vector]) => {
                if (exactConvexClassify(hull, p) === 0) {
                    return;
                }
                const base = new PointInPolygon2(hull).contains(p);
                const reversed = new PointInPolygon2([...hull].reverse());
                expect(reversed.contains(p)).toBe(base);
                for (let k = 1; k < hull.length; ++k) {
                    const rotated = new PointInPolygon2(
                        [...hull.slice(k), ...hull.slice(0, k)]);
                    expect(rotated.contains(p)).toBe(base);
                }
            });
    });

    // Non-convex case: rectilinear polygons, where every edge is axis
    // aligned, so a test point at half-integer coordinates is never on the
    // boundary and the even-odd answer is unambiguous. Cross-checked against
    // the independent PNPOLY reference above.
    it('Contains matches PNPOLY on non-convex rectilinear polygons', () => {
        const shapes: Vector[][] = [
            // L shape.
            [v(0, 0), v(4, 0), v(4, 2), v(2, 2), v(2, 5), v(0, 5)],
            // Plus shape.
            [v(1, 0), v(2, 0), v(2, 1), v(3, 1), v(3, 2), v(2, 2),
                v(2, 3), v(1, 3), v(1, 2), v(0, 2), v(0, 1), v(1, 1)],
            // Comb with three teeth.
            [v(0, 0), v(6, 0), v(6, 4), v(5, 4), v(5, 1), v(4, 1),
                v(4, 4), v(3, 4), v(3, 1), v(2, 1), v(2, 4), v(1, 4),
                v(1, 1), v(0, 1)],
            // Spiral (strongly non-convex).
            [v(0, 0), v(5, 0), v(5, 5), v(1, 5), v(1, 2), v(3, 2),
                v(3, 3), v(2, 3), v(2, 4), v(4, 4), v(4, 1), v(0, 1)]
        ];
        for (const shape of shapes) {
            const query = new PointInPolygon2(shape);
            const reversed = new PointInPolygon2([...shape].reverse());
            for (let x = -1.5; x <= 7.5; x += 1) {
                for (let y = -1.5; y <= 6.5; y += 1) {
                    const p = v(x, y);
                    const expected = referenceInside(shape, p);
                    expect(query.contains(p)).toBe(expected);
                    expect(reversed.contains(p)).toBe(expected);
                }
            }
        }
    });

    // Degenerate configurations documented by the port.
    it('rejects fewer than three vertices and non-2D vertices', () => {
        expect(() => new PointInPolygon2([v(0, 0), v(1, 0)]))
            .toThrow('PointInPolygon2: at least 3 vertices are required.');
        expect(() => new PointInPolygon2(
            [v(0, 0), v(1, 0), Vector.fromArray([0, 1, 0])]))
            .toThrow('PointInPolygon2: the vertices must be 2D.');
    });

    // Upstream returns false from ContainsQuadrilateral when the polygon does
    // not have exactly four vertices; the return is indistinguishable from
    // "outside". Pinned.
    it('containsQuadrilateral returns false when the count is not four', () => {
        const tri = new PointInPolygon2([v(0, 0), v(4, 0), v(0, 4)]);
        expect(tri.containsQuadrilateral(v(1, 1))).toBe(false);
        expect(tri.containsConvexOrderN(v(1, 1))).toBe(true);
    });
});
