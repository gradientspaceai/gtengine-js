import { describe, it, expect } from 'vitest';
import { ConvexHullSimplePolygon } from '../src/ConvexHullSimplePolygon.js';
import { Vector } from '../src/Vector.js';
import { ConvexHull2 } from '../src/ConvexHull2.js';
import { check, fc, scaled } from './helpers/arbitraries.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Twice the signed area of the polygon; positive for counterclockwise order.
function twiceSignedArea(points: Vector[]): number {
    let area = 0;
    for (let i = 0, j = points.length - 1; i < points.length; j = i++) {
        area += points[j].values[0] * points[i].values[1]
            - points[i].values[0] * points[j].values[1];
    }
    return area;
}

// An independent convex hull (Andrew's monotone chain), returning the hull
// points in counterclockwise order with collinear points removed.
function monotoneChainHull(points: Vector[]): Vector[] {
    const sorted = points.slice().sort((a, b) =>
        (a.values[0] !== b.values[0]
            ? a.values[0] - b.values[0]
            : a.values[1] - b.values[1]));
    const cross = (o: Vector, a: Vector, b: Vector): number =>
        (a.values[0] - o.values[0]) * (b.values[1] - o.values[1])
        - (a.values[1] - o.values[1]) * (b.values[0] - o.values[0]);

    const build = (input: Vector[]): Vector[] => {
        const chain: Vector[] = [];
        for (const p of input) {
            while (chain.length >= 2
                && cross(chain[chain.length - 2], chain[chain.length - 1], p) <= 0) {
                chain.pop();
            }
            chain.push(p);
        }
        chain.pop();
        return chain;
    };

    const lower = build(sorted);
    const upper = build(sorted.slice().reverse());
    return lower.concat(upper);
}

// Compare two counterclockwise cycles of points up to a rotation.
function sameCycle(a: Vector[], b: Vector[]): boolean {
    if (a.length !== b.length) {
        return false;
    }
    const key = (p: Vector): string => `${p.values[0]},${p.values[1]}`;
    for (let shift = 0; shift < a.length; ++shift) {
        let match = true;
        for (let i = 0; i < a.length; ++i) {
            if (key(a[i]) !== key(b[(i + shift) % b.length])) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }
    return false;
}

// Generate a simple (nonself-intersecting) counterclockwise polygon: points
// sampled at increasing angles around the origin with random radii form a
// star-shaped, hence simple, polygon.
function randomStarShapedPolygon(random: () => number, n: number): Vector[] {
    const angles: number[] = [];
    for (let i = 0; i < n; ++i) {
        angles.push(2 * Math.PI * random());
    }
    angles.sort((a, b) => a - b);
    return angles.map(angle => {
        const radius = 1 + 4 * random();
        return v2(radius * Math.cos(angle), radius * Math.sin(angle));
    });
}

describe('ConvexHullSimplePolygon', () => {
    const query = new ConvexHullSimplePolygon();

    it('returns all vertices of a convex counterclockwise polygon', () => {
        const square = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const hull = query.compute(square);
        expect(hull.length).toBe(4);
        expect(hull).toEqual([3, 0, 1, 2]);
        // The hull is counterclockwise ordered.
        expect(twiceSignedArea(hull.map(i => square[i]))).toBeGreaterThan(0);
    });

    it('returns the three vertices of a triangle', () => {
        const triangle = [v2(0, 0), v2(1, 0), v2(0, 1)];
        const hull = query.compute(triangle);
        expect(hull.length).toBe(3);
        expect(new Set(hull)).toEqual(new Set([0, 1, 2]));
        expect(twiceSignedArea(hull.map(i => triangle[i]))).toBeGreaterThan(0);
    });

    it('discards a reflex vertex (hand-worked example)', () => {
        // A square with a notch cut into the top edge. Vertex 3 = (2,2) is
        // reflex and is not on the hull.
        const polygon = [v2(0, 0), v2(4, 0), v2(4, 4), v2(2, 2), v2(0, 4)];
        const hull = query.compute(polygon);
        expect(hull).toEqual([4, 0, 1, 2]);
        expect(hull).not.toContain(3);
        const hullPoints = hull.map(i => polygon[i]);
        expect(twiceSignedArea(hullPoints)).toBeCloseTo(2 * 16, 12);
    });

    it('discards several reflex vertices of a comb-shaped polygon', () => {
        // Counterclockwise comb with two notches in the top edge.
        const polygon = [
            v2(0, 0), v2(6, 0), v2(6, 4),
            v2(5, 1), v2(4, 4),
            v2(3, 1), v2(2, 4),
            v2(0, 4)
        ];
        const hull = query.compute(polygon);
        const hullPoints = hull.map(i => polygon[i]);
        // The hull is the bounding rectangle [0,6]x[0,4].
        expect(twiceSignedArea(hullPoints)).toBeCloseTo(2 * 24, 12);
        expect(hull).toEqual([7, 0, 1, 2]);
        // Vertices 3 and 5 are reflex; 4 and 6 lie on the hull edge from
        // (6,4) to (0,4) and are dropped as collinear.
        for (const i of [3, 4, 5, 6]) {
            expect(hull).not.toContain(i);
        }
    });

    it('keeps a collinear vertex on a hull edge (degenerate input)', () => {
        // The first three vertices are collinear, which is the degenerate
        // start case for Melkman's algorithm. The collinear vertex 1 lies on
        // a hull edge and the port keeps it, as upstream does.
        const polygon = [v2(0, 0), v2(2, 0), v2(4, 0), v2(4, 4), v2(0, 4)];
        const hull = query.compute(polygon);
        expect(hull).toEqual([4, 0, 1, 2, 3]);
        expect(twiceSignedArea(hull.map(i => polygon[i]))).toBeCloseTo(2 * 16, 12);
    });

    it('throws when the polygon has fewer than 3 vertices', () => {
        expect(() => query.compute([v2(0, 0), v2(1, 0)]))
            .toThrow('The input polygon must have at least 3 vertices.');
        expect(() => query.compute([])).toThrow();
    });

    it('agrees with a monotone-chain hull for random simple polygons', () => {
        const random = makeRandom(31415926);
        for (let trial = 0; trial < 200; ++trial) {
            const n = 3 + Math.floor(20 * random());
            const polygon = randomStarShapedPolygon(random, n);
            if (twiceSignedArea(polygon) < 1e-6) {
                // The angular sampling can produce a clockwise (or nearly
                // degenerate) polygon; the query documents counterclockwise
                // input, and reversing a simple polygon keeps it simple.
                polygon.reverse();
            }
            if (twiceSignedArea(polygon) < 1e-6) {
                continue;
            }

            const hull = query.compute(polygon);
            const hullPoints = hull.map(i => polygon[i]);
            const expected = monotoneChainHull(polygon);
            expect(sameCycle(hullPoints, expected)).toBe(true);

            // The hull is counterclockwise ordered and its area is at least
            // that of the polygon.
            expect(twiceSignedArea(hullPoints)).toBeGreaterThan(0);
            expect(twiceSignedArea(hullPoints))
                .toBeGreaterThanOrEqual(twiceSignedArea(polygon) - 1e-9);
        }
    });

    it('contains every polygon vertex inside the computed hull', () => {
        const random = makeRandom(2718281);
        for (let trial = 0; trial < 100; ++trial) {
            const polygon = randomStarShapedPolygon(random, 3 + Math.floor(15 * random()));
            if (twiceSignedArea(polygon) < 1e-6) {
                polygon.reverse();
            }
            if (twiceSignedArea(polygon) < 1e-6) {
                continue;
            }
            const hullPoints = query.compute(polygon).map(i => polygon[i]);
            for (const p of polygon) {
                for (let i = 0, j = hullPoints.length - 1; i < hullPoints.length; j = i++) {
                    const ex = hullPoints[i].values[0] - hullPoints[j].values[0];
                    const ey = hullPoints[i].values[1] - hullPoints[j].values[1];
                    const px = p.values[0] - hullPoints[j].values[0];
                    const py = p.values[1] - hullPoints[j].values[1];
                    // Counterclockwise hull: every point is left of or on
                    // each directed hull edge.
                    expect(ex * py - ey * px).toBeGreaterThan(-1e-9);
                }
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// The reference is ConvexHull2, whose predicates are exact (BSRational), run
// on the same vertex array. Melkman's algorithm drops collinear vertices only
// when they are strictly inside an edge of the hull (WhichSide <= 0 pops), so
// the generator filters out configurations with a nearly collinear triple:
// there the two algorithms may legitimately differ on which vertices lie "on"
// the hull, and the double-precision DotPerp cannot decide the case anyway.
// ---------------------------------------------------------------------------

const SLOT_COUNT = 64;

// A star-shaped (hence simple) counterclockwise polygon: vertices at strictly
// increasing angles about the origin with random radii. Every consecutive
// angular gap must be less than pi, otherwise the origin is not in the
// kernel, the polygon need not be simple and Melkman's algorithm does not
// apply to it.
const starPolygon = fc.record({
    slots: fc.uniqueArray(fc.integer({ min: 0, max: SLOT_COUNT - 1 }),
        { minLength: 3, maxLength: 12 })
        .filter(slots => {
            const sorted = [...slots].sort((a, b) => a - b);
            for (let i = 0; i < sorted.length; ++i) {
                const next = (i + 1 < sorted.length)
                    ? sorted[i + 1] : sorted[0] + SLOT_COUNT;
                if (next - sorted[i] >= SLOT_COUNT / 2) { return false; }
            }
            return true;
        }),
    radii: fc.array(fc.integer({ min: 2, max: 12 }), { minLength: 12, maxLength: 12 }),
    phase: scaled(0, 2 * Math.PI, 32)
}).map(({ slots, radii, phase }) => {
    const sorted = [...slots].sort((a, b) => a - b);
    return sorted.map((slot, k) => {
        const angle = phase + (2 * Math.PI * slot) / SLOT_COUNT;
        const radius = 0.5 * radii[k];
        return v2(radius * Math.cos(angle), radius * Math.sin(angle));
    });
}).filter(polygon => {
    // General position: no triple of vertices is nearly collinear.
    for (let i = 0; i < polygon.length; ++i) {
        for (let j = i + 1; j < polygon.length; ++j) {
            for (let k = j + 1; k < polygon.length; ++k) {
                const ax = polygon[j].values[0] - polygon[i].values[0];
                const ay = polygon[j].values[1] - polygon[i].values[1];
                const bx = polygon[k].values[0] - polygon[i].values[0];
                const by = polygon[k].values[1] - polygon[i].values[1];
                if (Math.abs(ax * by - ay * bx) < 1e-6) { return false; }
            }
        }
    }
    return true;
});

// The counterclockwise hull index list computed by ConvexHull2 (exact
// predicates), used as the independent reference.
function referenceHull(polygon: Vector[]): number[] {
    const hull = new ConvexHull2();
    expect(hull.compute(polygon)).toBe(true);
    return [...hull.getHull()];
}

// Compare two cyclic index sequences up to rotation.
function sameIndexCycle(a: readonly number[], b: readonly number[]): boolean {
    if (a.length !== b.length) { return false; }
    for (let shift = 0; shift < a.length; ++shift) {
        let match = true;
        for (let i = 0; i < a.length; ++i) {
            if (a[i] !== b[(i + shift) % b.length]) { match = false; break; }
        }
        if (match) { return true; }
    }
    return false;
}

describe('ConvexHullSimplePolygon verification', () => {
    const query = new ConvexHullSimplePolygon();

    it('agrees with ConvexHull2 on random simple polygons', () => {
        check(starPolygon, (polygon) => {
            const hull = query.compute(polygon);
            expect(sameIndexCycle(hull, referenceHull(polygon))).toBe(true);
        });
    });

    it('the hull is a strictly convex counterclockwise polygon', () => {
        check(starPolygon, (polygon) => {
            const hull = query.compute(polygon);
            expect(hull.length).toBeGreaterThanOrEqual(3);
            expect(new Set(hull).size).toBe(hull.length);
            const points = hull.map(i => polygon[i]);
            expect(twiceSignedArea(points)).toBeGreaterThan(0);
            for (let i = 0; i < points.length; ++i) {
                const p0 = points[i];
                const p1 = points[(i + 1) % points.length];
                const p2 = points[(i + 2) % points.length];
                const cross = (p1.values[0] - p0.values[0]) * (p2.values[1] - p0.values[1])
                    - (p1.values[1] - p0.values[1]) * (p2.values[0] - p0.values[0]);
                expect(cross).toBeGreaterThan(0);
            }
        });
    });

    it('every polygon vertex lies inside or on the hull', () => {
        check(starPolygon, (polygon) => {
            const points = query.compute(polygon).map(i => polygon[i]);
            for (const q of polygon) {
                for (let i = 0; i < points.length; ++i) {
                    const p0 = points[i];
                    const p1 = points[(i + 1) % points.length];
                    const cross = (p1.values[0] - p0.values[0]) * (q.values[1] - p0.values[1])
                        - (p1.values[1] - p0.values[1]) * (q.values[0] - p0.values[0]);
                    // Left of every directed hull edge, up to round-off in the
                    // difference-of-products (magnitudes are O(100)).
                    expect(cross).toBeGreaterThan(-1e-9);
                }
            }
        });
    });

    it('is invariant under a cyclic rotation of the input polygon', () => {
        check(fc.tuple(starPolygon, fc.integer({ min: 0, max: 11 })),
            ([polygon, shift]) => {
                const s = shift % polygon.length;
                const rotated = polygon.slice(s).concat(polygon.slice(0, s));
                const hull = query.compute(polygon).map(i => polygon[i]);
                const rotatedHull = query.compute(rotated).map(i => rotated[i]);
                expect(sameCycle(rotatedHull, hull)).toBe(true);
            });
    });

    it('is equivariant under rotation of the plane', () => {
        check(fc.tuple(starPolygon, scaled(-Math.PI, Math.PI, 32)),
            ([polygon, angle]) => {
                const c = Math.cos(angle), s = Math.sin(angle);
                const moved = polygon.map(p => v2(
                    c * p.values[0] - s * p.values[1] + 3,
                    s * p.values[0] + c * p.values[1] - 2));
                // A rigid motion maps the hull vertices to the hull vertices,
                // so the index list is unchanged.
                expect(query.compute(moved)).toEqual(query.compute(polygon));
            });
    });

    it('returns every vertex of a convex polygon, keeping the cycle', () => {
        check(fc.record({
            slots: fc.uniqueArray(fc.integer({ min: 0, max: SLOT_COUNT - 1 }),
                { minLength: 3, maxLength: 10 }),
            phase: scaled(0, 2 * Math.PI, 32),
            radius: scaled(1, 6, 8)
        }), ({ slots, phase, radius }) => {
            // Vertices on a circle are in convex position, so the hull is the
            // whole polygon.
            const sorted = [...slots].sort((a, b) => a - b);
            const polygon = sorted.map(slot => {
                const angle = phase + (2 * Math.PI * slot) / SLOT_COUNT;
                return v2(radius * Math.cos(angle), radius * Math.sin(angle));
            });
            const hull = query.compute(polygon);
            expect(hull.length).toBe(polygon.length);
            expect([...hull].sort((a, b) => a - b))
                .toEqual(polygon.map((_, i) => i));
            expect(sameIndexCycle(hull, polygon.map((_, i) => i))).toBe(true);
        });
    });

    it('rejects polygons with fewer than three vertices', () => {
        check(fc.array(fc.tuple(scaled(-5, 5, 32), scaled(-5, 5, 32)),
            { minLength: 0, maxLength: 2 }), (raw) => {
            const polygon = raw.map(([x, y]) => v2(x, y));
            expect(() => query.compute(polygon))
                .toThrow('The input polygon must have at least 3 vertices.');
        });
    });
});
