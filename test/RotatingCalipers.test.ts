import { describe, expect, it } from 'vitest';
import { RotatingCalipers } from '../src/RotatingCalipers.js';
import type { RotatingCalipersAntipode } from '../src/RotatingCalipers.js';
import { ConvexHull2 } from '../src/ConvexHull2.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The indices of the vertices that remain after duplicate and collinear
// vertices are removed: vertex i is kept when the edge arriving at it and
// the edge leaving it are not parallel.
function retainedIndices(vertices: readonly Vector[]): number[] {
    const n = vertices.length;
    const indices: number[] = [];
    for (let i = 0; i < n; ++i) {
        const prev = vertices[(i + n - 1) % n];
        const curr = vertices[i];
        const next = vertices[(i + 1) % n];
        const cross = (curr.get(0) - prev.get(0))
            * (next.get(1) - curr.get(1))
            - (curr.get(1) - prev.get(1)) * (next.get(0) - curr.get(0));
        if (cross !== 0) {
            indices.push(i);
        }
    }
    return indices;
}

// Perpendicular distance from a point to the line through E0 and E1.
function distanceToLine(point: Vector, e0: Vector, e1: Vector): number {
    const dx = e1.get(0) - e0.get(0);
    const dy = e1.get(1) - e0.get(1);
    const len = Math.hypot(dx, dy);
    return Math.abs((point.get(0) - e0.get(0)) * dy
        - (point.get(1) - e0.get(1)) * dx) / len;
}

// The brute-force width of a convex polygon: the minimum over all edges of
// the maximum distance from a vertex to that edge's line.
function bruteForceWidth(vertices: readonly Vector[]): number {
    const n = vertices.length;
    let width = Number.MAX_VALUE;
    for (let i = 0; i < n; ++i) {
        const e0 = vertices[i];
        const e1 = vertices[(i + 1) % n];
        if (e0.equals(e1)) {
            continue;
        }
        let maxDist = 0;
        for (const vertex of vertices) {
            maxDist = Math.max(maxDist, distanceToLine(vertex, e0, e1));
        }
        if (maxDist > 0 && maxDist < width) {
            width = maxDist;
        }
    }
    return width;
}

// Verify that each antipode's vertex really is an extreme point in the
// direction perpendicular to its edge.
function checkAntipodes(vertices: readonly Vector[],
    antipodes: readonly RotatingCalipersAntipode[]): void {
    for (const antipode of antipodes) {
        const e0 = vertices[antipode.edge[0]];
        const e1 = vertices[antipode.edge[1]];
        const vertexDistance = distanceToLine(vertices[antipode.vertex],
            e0, e1);
        let maxDistance = 0;
        for (const vertex of vertices) {
            maxDistance = Math.max(maxDistance,
                distanceToLine(vertex, e0, e1));
        }
        expect(vertexDistance).toBeCloseTo(maxDistance, 9);
    }
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (1664525 * state + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('RotatingCalipers', () => {
    it('computes the antipodes of a unit square', () => {
        const square = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const antipodes = RotatingCalipers.computeAntipodes(square);
        expect(antipodes.length).toBe(4);
        checkAntipodes(square, antipodes);

        // Every reported edge is a polygon edge and the antipodal vertex is
        // diagonally opposite it. Note that for a polygon whose caliper
        // angles all tie (a square), upstream can report the same edge twice
        // and omit another; see the "Upstream bug suspects" note. The pairs
        // that are reported are all valid antipodal pairs.
        for (const antipode of antipodes) {
            expect(antipode.edge[1]).toBe((antipode.edge[0] + 1) % 4);
            // Both vertices of the opposite edge are at the same (maximum)
            // distance from the edge line, so either is a valid antipode.
            expect([(antipode.edge[0] + 2) % 4, (antipode.edge[0] + 3) % 4])
                .toContain(antipode.vertex);
        }
    });

    it('computes the antipodes of a triangle', () => {
        const triangle = [v2(0, 0), v2(4, 0), v2(0, 3)];
        const antipodes = RotatingCalipers.computeAntipodes(triangle);
        expect(antipodes.length).toBe(3);
        checkAntipodes(triangle, antipodes);
        // For a triangle, the antipodal vertex of an edge is the vertex not
        // on the edge.
        for (const antipode of antipodes) {
            expect(antipode.vertex).not.toBe(antipode.edge[0]);
            expect(antipode.vertex).not.toBe(antipode.edge[1]);
        }
    });

    it('removes collinear vertices', () => {
        // A square with an extra vertex at the midpoint of each edge.
        const polygon = [
            v2(0, 0), v2(1, 0), v2(2, 0),
            v2(2, 1), v2(2, 2),
            v2(1, 2), v2(0, 2),
            v2(0, 1)
        ];
        const antipodes = RotatingCalipers.computeAntipodes(polygon);
        expect(antipodes.length).toBe(4);
        // The retained vertices are the corners.
        const used = new Set<number>();
        for (const antipode of antipodes) {
            used.add(antipode.vertex);
            used.add(antipode.edge[0]);
            used.add(antipode.edge[1]);
        }
        expect([...used].sort((a, b) => a - b)).toEqual([0, 2, 4, 6]);
        checkAntipodes(polygon, antipodes);
    });

    it('removes duplicate vertices without losing corners', () => {
        // Upstream drops a genuine corner whenever it is preceded by a
        // duplicate point (the zero-length edge makes DotPerp vanish); the
        // port compares against the most recent nonzero edge instead, so all
        // four corners of this square survive.
        const polygon = [
            v2(0, 0), v2(0, 0), v2(1, 0), v2(1, 1), v2(1, 1), v2(0, 1)
        ];
        const antipodes = RotatingCalipers.computeAntipodes(polygon);
        expect(antipodes.length).toBe(4);
        checkAntipodes(polygon, antipodes);
        // The retained corners are one representative of each duplicate
        // pair plus the two unique corners.
        const used = new Set<number>();
        for (const antipode of antipodes) {
            used.add(antipode.edge[0]);
            used.add(antipode.edge[1]);
        }
        expect([...used].sort((a, b) => a - b)).toEqual([1, 2, 4, 5]);
    });

    it('removes leading and trailing duplicate runs', () => {
        const polygon = [
            v2(0, 0), v2(4, 0), v2(4, 0), v2(4, 0), v2(0, 3), v2(0, 0)
        ];
        const antipodes = RotatingCalipers.computeAntipodes(polygon);
        expect(antipodes.length).toBe(3);
        checkAntipodes(polygon, antipodes);
    });

    it('reproduces the width of a rectangle', () => {
        const rectangle = [v2(0, 0), v2(10, 0), v2(10, 3), v2(0, 3)];
        const antipodes = RotatingCalipers.computeAntipodes(rectangle);
        let width = Number.MAX_VALUE;
        for (const antipode of antipodes) {
            width = Math.min(width, distanceToLine(
                rectangle[antipode.vertex],
                rectangle[antipode.edge[0]],
                rectangle[antipode.edge[1]]));
        }
        expect(width).toBeCloseTo(3, 12);
    });

    it('throws for fewer than three vertices', () => {
        expect(() => RotatingCalipers.computeAntipodes([v2(0, 0), v2(1, 1)]))
            .toThrow('The convex polygon must have at least 3 noncollinear vertices.');
    });

    it('throws for a degenerate (all collinear) polygon', () => {
        expect(() => RotatingCalipers.computeAntipodes(
            [v2(0, 0), v2(1, 1), v2(2, 2), v2(3, 3)]))
            .toThrow('The convex polygon must have at least 3 noncollinear vertices.');
    });

    it('throws for non-2D vertices', () => {
        expect(() => RotatingCalipers.computeAntipodes([
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0])
        ])).toThrow();
    });

    it('matches brute force on random convex polygons', () => {
        const random = makeRandom(1234567);
        for (let trial = 0; trial < 30; ++trial) {
            const numVertices = 3 + (trial % 10);
            // Random distinct angles sorted increasing gives a strictly
            // convex, counterclockwise polygon on an ellipse.
            const angles: number[] = [];
            for (let i = 0; i < numVertices; ++i) {
                angles.push(2 * Math.PI * random());
            }
            angles.sort((a, b) => a - b);
            let degenerate = false;
            for (let i = 1; i < numVertices; ++i) {
                if (angles[i] - angles[i - 1] < 1e-3) {
                    degenerate = true;
                }
            }
            if (degenerate) {
                continue;
            }

            const a = 1 + 3 * random();
            const b = 1 + 3 * random();
            const vertices = angles.map(angle =>
                v2(a * Math.cos(angle), b * Math.sin(angle)));

            const antipodes = RotatingCalipers.computeAntipodes(vertices);
            const retained = retainedIndices(vertices);
            expect(antipodes.length).toBe(retained.length);
            checkAntipodes(vertices, antipodes);

            // The polygon width from the calipers matches the brute-force
            // width.
            let width = Number.MAX_VALUE;
            for (const antipode of antipodes) {
                width = Math.min(width, distanceToLine(
                    vertices[antipode.vertex],
                    vertices[antipode.edge[0]],
                    vertices[antipode.edge[1]]));
            }
            expect(width).toBeCloseTo(bruteForceWidth(vertices), 8);

            // Each polygon edge occurs exactly once as an antipodal edge.
            const edgeKeys = new Set(
                antipodes.map(ap => `${ap.edge[0]},${ap.edge[1]}`));
            expect(edgeKeys.size).toBe(antipodes.length);
            for (const antipode of antipodes) {
                const i = retained.indexOf(antipode.edge[0]);
                expect(i).toBeGreaterThanOrEqual(0);
                expect(antipode.edge[1])
                    .toBe(retained[(i + 1) % retained.length]);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (group V10): property-based cross-checks against exact
// integer brute force. The polygons are convex hulls of small lattice points,
// so every predicate below is evaluated exactly with bigint arithmetic; there
// is no tolerance anywhere in this block.
// ---------------------------------------------------------------------------

// Twice the signed area of the triangle (a, b, c), exactly.
function crossExact(a: Vector, b: Vector, c: Vector): bigint {
    const ax = BigInt(a.get(0)), ay = BigInt(a.get(1));
    const bx = BigInt(b.get(0)), by = BigInt(b.get(1));
    const cx = BigInt(c.get(0)), cy = BigInt(c.get(1));
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

function sqrDistExact(a: Vector, b: Vector): bigint {
    const dx = BigInt(a.get(0)) - BigInt(b.get(0));
    const dy = BigInt(a.get(1)) - BigInt(b.get(1));
    return dx * dx + dy * dy;
}

// A perpendicular distance represented exactly as num / sqrt(den) with
// num >= 0 and den > 0, so distances can be compared without a square root.
interface ExactDistance { num: bigint; den: bigint; }

function distLess(x: ExactDistance, y: ExactDistance): boolean {
    return x.num * x.num * y.den < y.num * y.num * x.den;
}

function distEqual(x: ExactDistance, y: ExactDistance): boolean {
    return x.num * x.num * y.den === y.num * y.num * x.den;
}

// Distance from `point` to the line through (e0, e1), as num / sqrt(den).
function distToLineExact(point: Vector, e0: Vector, e1: Vector): ExactDistance {
    const c = crossExact(e0, e1, point);
    return { num: c < 0n ? -c : c, den: sqrDistExact(e0, e1) };
}

// The strictly convex, counterclockwise hull of `points`, or null when the
// input is degenerate (fewer than three hull vertices) or the hull retains a
// collinear vertex (which the rotating calipers would remove).
function strictHull(points: readonly Vector[]): Vector[] | null {
    const hullQuery = new ConvexHull2();
    if (!hullQuery.compute(points)) {
        return null;
    }
    const hull = hullQuery.getHull().map(i => points[i]);
    if (hull.length < 3) {
        return null;
    }
    for (let i = 0; i < hull.length; ++i) {
        const prev = hull[(i + hull.length - 1) % hull.length];
        const next = hull[(i + 1) % hull.length];
        if (crossExact(prev, hull[i], next) <= 0n) {
            return null;
        }
    }
    return hull;
}

// The minimum, over the reported antipodal pairs, of the distance from the
// antipodal vertex to its edge line: the polygon width.
function calipersWidth(hull: readonly Vector[],
    antipodes: readonly RotatingCalipersAntipode[]): ExactDistance {
    let width: ExactDistance | null = null;
    for (const antipode of antipodes) {
        const d = distToLineExact(hull[antipode.vertex],
            hull[antipode.edge[0]], hull[antipode.edge[1]]);
        if (width === null || distLess(d, width)) {
            width = d;
        }
    }
    return width as ExactDistance;
}

// The maximum, over the reported antipodal pairs, of the squared distance
// between the antipodal vertex and either endpoint of its edge: the polygon
// diameter.
function calipersSqrDiameter(hull: readonly Vector[],
    antipodes: readonly RotatingCalipersAntipode[]): bigint {
    let diameter = 0n;
    for (const antipode of antipodes) {
        for (const e of antipode.edge) {
            const d = sqrDistExact(hull[antipode.vertex], hull[e]);
            if (d > diameter) {
                diameter = d;
            }
        }
    }
    return diameter;
}

const latticeHull = fc.array(latticeVector(2, -8, 8),
    { minLength: 3, maxLength: 14 })
    .map(points => strictHull(points))
    .filter((hull): hull is Vector[] => hull !== null);

// In a counterclockwise convex polygon the edge directions turn monotonically
// through 2*pi, so two edges are parallel only when they are antiparallel.
// Such a pair makes two caliper rotation angles tie, which is where upstream's
// ComputeNextAntipode duplicates one edge and omits another (issue #286).
function hasParallelEdges(hull: readonly Vector[]): boolean {
    const n = hull.length;
    const edge = (i: number): [bigint, bigint] => [
        BigInt(hull[(i + 1) % n].get(0)) - BigInt(hull[i].get(0)),
        BigInt(hull[(i + 1) % n].get(1)) - BigInt(hull[i].get(1))];
    for (let i = 0; i < n; ++i) {
        const a = edge(i);
        for (let j = i + 1; j < n; ++j) {
            const b = edge(j);
            if (a[0] * b[1] - a[1] * b[0] === 0n) {
                return true;
            }
        }
    }
    return false;
}

const latticeHullNoParallelEdges =
    latticeHull.filter(hull => !hasParallelEdges(hull));

describe('RotatingCalipers verification', () => {
    it('reports one antipodal pair per polygon edge, in edge order', () => {
        check(latticeHullNoParallelEdges, hull => {
            const antipodes = RotatingCalipers.computeAntipodes(hull);
            expect(antipodes.length).toBe(hull.length);
            // Consecutive antipodes advance the edge by one; when no two
            // caliper angles tie, every edge of the polygon occurs exactly
            // once.
            const seen = new Set<number>();
            for (const antipode of antipodes) {
                expect(antipode.edge[1])
                    .toBe((antipode.edge[0] + 1) % hull.length);
                seen.add(antipode.edge[0]);
            }
            expect(seen.size).toBe(hull.length);
        });
    });

    it('always emits exactly one antipode per retained vertex', () => {
        // This holds even for the tie case, where the emitted edges repeat.
        check(latticeHull, hull => {
            const antipodes = RotatingCalipers.computeAntipodes(hull);
            expect(antipodes.length).toBe(hull.length);
            for (const antipode of antipodes) {
                expect(antipode.edge[1])
                    .toBe((antipode.edge[0] + 1) % hull.length);
            }
        });
    });

    it('duplicates an edge when the polygon has parallel edges (upstream)', () => {
        // Preserved upstream quirk (issue #286, item 4): when caliper angles
        // tie, one edge is reported twice and another never appears. The unit
        // square is the documented case; this trapezoid with the antiparallel
        // edge pair (1,4) / (-1,-4) is a second, less symmetric one found by
        // the property above. Every reported pair is still genuinely
        // antipodal, which the "exactly extreme" property checks.
        const polygon = [v2(-1, 4), v2(0, -4), v2(1, 0), v2(0, 8)];
        const antipodes = RotatingCalipers.computeAntipodes(polygon);
        expect(antipodes.length).toBe(4);
        const edges = antipodes.map(ap => ap.edge[0]);
        expect(new Set(edges).size).toBe(3);
        checkAntipodes(polygon, antipodes);
    });

    it('each antipodal vertex is exactly extreme for its edge', () => {
        check(latticeHull, hull => {
            const antipodes = RotatingCalipers.computeAntipodes(hull);
            for (const antipode of antipodes) {
                const e0 = hull[antipode.edge[0]];
                const e1 = hull[antipode.edge[1]];
                const chosen = distToLineExact(hull[antipode.vertex], e0, e1);
                for (const vertex of hull) {
                    const other = distToLineExact(vertex, e0, e1);
                    // The antipodal vertex maximizes the distance to the edge
                    // line; ties are allowed (two parallel edges).
                    expect(distLess(chosen, other)).toBe(false);
                }
                // The polygon is counterclockwise, so the extreme vertex lies
                // strictly to the left of the directed edge.
                expect(crossExact(e0, e1, hull[antipode.vertex]) > 0n)
                    .toBe(true);
            }
        });
    });

    it('the calipers width equals the brute-force width', () => {
        check(latticeHull, hull => {
            const antipodes = RotatingCalipers.computeAntipodes(hull);
            // Brute force: minimum over edges of the maximum vertex distance.
            let bruteWidth: ExactDistance | null = null;
            for (let i = 0; i < hull.length; ++i) {
                const e0 = hull[i];
                const e1 = hull[(i + 1) % hull.length];
                let maxDist: ExactDistance =
                    { num: 0n, den: sqrDistExact(e0, e1) };
                for (const vertex of hull) {
                    const d = distToLineExact(vertex, e0, e1);
                    if (distLess(maxDist, d)) {
                        maxDist = d;
                    }
                }
                if (bruteWidth === null || distLess(maxDist, bruteWidth)) {
                    bruteWidth = maxDist;
                }
            }
            expect(distEqual(calipersWidth(hull, antipodes),
                bruteWidth as ExactDistance)).toBe(true);
        });
    });

    it('the calipers diameter equals the brute-force diameter', () => {
        check(latticeHull, hull => {
            const antipodes = RotatingCalipers.computeAntipodes(hull);
            let bruteDiameter = 0n;
            for (let i = 0; i < hull.length; ++i) {
                for (let j = i + 1; j < hull.length; ++j) {
                    const d = sqrDistExact(hull[i], hull[j]);
                    if (d > bruteDiameter) {
                        bruteDiameter = d;
                    }
                }
            }
            expect(calipersSqrDiameter(hull, antipodes) === bruteDiameter)
                .toBe(true);
        });
    });

    it('width and diameter do not depend on the starting vertex', () => {
        check(fc.tuple(latticeHull, fc.nat(20)), ([hull, shift]) => {
            const k = shift % hull.length;
            const rotated = hull.map((_, i) => hull[(i + k) % hull.length]);
            const a0 = RotatingCalipers.computeAntipodes(hull);
            const a1 = RotatingCalipers.computeAntipodes(rotated);
            expect(distEqual(calipersWidth(hull, a0),
                calipersWidth(rotated, a1))).toBe(true);
            expect(calipersSqrDiameter(hull, a0)
                === calipersSqrDiameter(rotated, a1)).toBe(true);
        });
    });

    it('duplicate and collinear vertices do not change the answer', () => {
        check(fc.tuple(latticeHull, fc.nat(1000)), ([hull, seed]) => {
            // Insert a duplicate of one vertex and the midpoint of one edge
            // (the coordinates are doubled so the midpoint stays on the
            // lattice). Neither is a corner, so the calipers must produce the
            // same width and diameter as the clean polygon. Upstream's
            // CreatePolygon drops the corner that follows the duplicate; see
            // the port note in src/RotatingCalipers.ts.
            const scaled = hull.map(v =>
                Vector.fromArray([2 * v.get(0), 2 * v.get(1)]));
            const i = seed % scaled.length;
            const augmented: Vector[] = [];
            for (let k = 0; k < scaled.length; ++k) {
                augmented.push(scaled[k]);
                if (k === i) {
                    augmented.push(scaled[k].clone());
                    const next = scaled[(k + 1) % scaled.length];
                    augmented.push(Vector.fromArray([
                        (scaled[k].get(0) + next.get(0)) / 2,
                        (scaled[k].get(1) + next.get(1)) / 2]));
                }
            }
            const clean = RotatingCalipers.computeAntipodes(scaled);
            const dirty = RotatingCalipers.computeAntipodes(augmented);
            expect(dirty.length).toBe(clean.length);
            expect(distEqual(calipersWidth(scaled, clean),
                calipersWidth(augmented, dirty))).toBe(true);
            expect(calipersSqrDiameter(scaled, clean)
                === calipersSqrDiameter(augmented, dirty)).toBe(true);
        });
    });
});
