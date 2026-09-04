import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2.js';
import {
    IncrementalDelaunay2, IncrementalDelaunay2SearchInfo
} from '../src/IncrementalDelaunay2.js';
import { Vector } from '../src/Vector.js';
import { check, fc, latticeVector } from './helpers/arbitraries.js';
import { inCircle2, orient2, twiceSignedAreaExact } from './helpers/exact.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A deterministic pseudorandom generator, so failures reproduce.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (Math.imul(s, 1664525) + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// The set of triangles as sorted vertex triples, so the comparison ignores
// both the triangle order and the rotation of each triangle.
function triangleSet(triangles: readonly (readonly number[])[]): Set<string> {
    const set = new Set<string>();
    for (const t of triangles) {
        set.add([t[0], t[1], t[2]].slice().sort((a, b) => a - b).join(','));
    }
    return set;
}

// The triangles of a batch Delaunay2 computation, as index triples.
function batchTriangles(points: Vector[]): number[][] {
    const delaunay = new Delaunay2();
    expect(delaunay.compute(points)).toBe(true);
    const triangles: number[][] = [];
    for (let t = 0; t < delaunay.getNumTriangles(); ++t) {
        triangles.push(delaunay.getTriangleIndices(t) as number[]);
    }
    return triangles;
}

// Insert the points into an incremental triangulation, finalize it, and
// return the triangles in terms of the *input* point indices.
function incrementalTriangles(points: Vector[],
    xMin: number, yMin: number, xMax: number, yMax: number): number[][] {
    const incremental = new IncrementalDelaunay2(xMin, yMin, xMax, yMax);
    const toInput = new Map<number, number>();
    for (let i = 0; i < points.length; ++i) {
        const index = incremental.insert(points[i]);
        if (!toInput.has(index)) {
            toInput.set(index, i);
        }
    }
    expect(incremental.finalizeTriangulation()).toBe(true);
    return incremental.getTriangles().map(t => [
        toInput.get(t[0]) as number,
        toInput.get(t[1]) as number,
        toInput.get(t[2]) as number
    ]);
}

// The signed area of the triangle, twice. Positive for counterclockwise.
function signedArea2(a: Vector, b: Vector, c: Vector): number {
    return (b.get(0) - a.get(0)) * (c.get(1) - a.get(1))
        - (c.get(0) - a.get(0)) * (b.get(1) - a.get(1));
}

// The InCircle determinant for a counterclockwise triangle <a,b,c>: positive
// when p is strictly inside the circumcircle.
function inCircle(a: Vector, b: Vector, c: Vector, p: Vector): number {
    const ax = a.get(0) - p.get(0), ay = a.get(1) - p.get(1);
    const bx = b.get(0) - p.get(0), by = b.get(1) - p.get(1);
    const cx = c.get(0) - p.get(0), cy = c.get(1) - p.get(1);
    const a2 = ax * ax + ay * ay;
    const b2 = bx * bx + by * by;
    const c2 = cx * cx + cy * cy;
    return a2 * (bx * cy - by * cx)
        - b2 * (ax * cy - ay * cx)
        + c2 * (ax * by - ay * bx);
}

describe('IncrementalDelaunay2', () => {
    it('rejects an invalid bounding rectangle and out-of-domain points', () => {
        expect(() => new IncrementalDelaunay2(1, 0, 0, 1)).toThrow();
        expect(() => new IncrementalDelaunay2(0, 1, 1, 0)).toThrow();
        expect(() => new IncrementalDelaunay2(0, 0, 0, 1)).toThrow();

        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        // The domain is open: boundary points are rejected.
        expect(() => triangulation.insert(v2(0, 0.5))).toThrow();
        expect(() => triangulation.insert(v2(1, 0.5))).toThrow();
        expect(() => triangulation.insert(v2(0.5, 1))).toThrow();
        expect(() => triangulation.insert(v2(2, 0.5))).toThrow();
        expect(() => triangulation.insert(Vector.fromArray([0.5, 0.5, 0.5]))).toThrow();
    });

    it('starts with the supertriangle and the input rectangle', () => {
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);

        // Three supervertices plus the four rectangle corners.
        expect(triangulation.getNumVertices()).toBe(7);
        const vertices = triangulation.getVertices();
        expect([vertices[3].get(0), vertices[3].get(1)]).toEqual([0, 0]);
        expect([vertices[4].get(0), vertices[4].get(1)]).toEqual([1, 0]);
        expect([vertices[5].get(0), vertices[5].get(1)]).toEqual([0, 1]);
        expect([vertices[6].get(0), vertices[6].get(1)]).toEqual([1, 1]);

        // Nine triangles in the graph, of which the two covering the input
        // rectangle use only Delaunay vertices.
        const graph = triangulation.getTriangulation();
        expect(graph.vertices.length).toBe(7);
        expect(graph.triangles.length).toBe(9);
        expect(triangulation.getNumTriangles()).toBe(2);
        expect(triangleSet(triangulation.getTriangles()))
            .toEqual(triangleSet([[3, 4, 6], [3, 6, 5]]));

        // The hull of the two rectangle triangles is the rectangle itself,
        // counterclockwise.
        const hull = triangulation.getHull();
        expect(hull.length).toBe(4);
        expect(new Set(hull)).toEqual(new Set([3, 4, 6, 5]));
    });

    it('returns the existing index when a point is inserted twice', () => {
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        const p = v2(0.25, 0.75);
        const first = triangulation.insert(p);
        expect(first).toBe(7);
        expect(triangulation.insert(p)).toBe(first);
        expect(triangulation.insert(v2(0.25, 0.75))).toBe(first);
        expect(triangulation.getNumVertices()).toBe(8);

        // Removing a point that is not in the triangulation reports -1.
        expect(triangulation.remove(v2(0.4, 0.4))).toBe(-1);
    });

    it('matches the batch Delaunay2 triangulation on a small point set', () => {
        const points = [
            v2(0.1, 0.1), v2(0.9, 0.1), v2(0.9, 0.9), v2(0.1, 0.9),
            v2(0.5, 0.5), v2(0.3, 0.7), v2(0.7, 0.3)
        ];
        const batch = triangleSet(batchTriangles(points));
        const incremental = triangleSet(incrementalTriangles(points, 0, 0, 1, 1));
        expect(incremental).toEqual(batch);
    });

    it('matches the batch Delaunay2 triangulation on random point sets', () => {
        for (const seed of [12345, 6789, 4242]) {
            const rnd = makeRandom(seed);
            const points: Vector[] = [];
            for (let i = 0; i < 30; ++i) {
                points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
            }
            const batch = triangleSet(batchTriangles(points));
            const incremental = triangleSet(incrementalTriangles(points, 0, 0, 1, 1));
            expect(incremental).toEqual(batch);
        }
    });

    it('produces counterclockwise triangles with the empty-circumcircle property', () => {
        const rnd = makeRandom(9001);
        const points: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
        }

        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        for (const p of points) {
            triangulation.insert(p);
        }
        triangulation.finalizeTriangulation();

        const vertices = triangulation.getVertices();
        const triangles = triangulation.getTriangles();
        expect(triangles.length).toBeGreaterThan(0);

        const used = new Set<number>();
        for (const t of triangles) {
            used.add(t[0]);
            used.add(t[1]);
            used.add(t[2]);
        }
        expect(used.size).toBe(points.length);

        for (const t of triangles) {
            const a = vertices[t[0]];
            const b = vertices[t[1]];
            const c = vertices[t[2]];

            // Counterclockwise and nondegenerate.
            expect(signedArea2(a, b, c)).toBeGreaterThan(0);

            // No other vertex of the Delaunay triangulation is strictly
            // inside the circumcircle. Only the vertices that the triangles
            // actually use qualify: the supervertices and the removed
            // rectangle corners are still in the vertex pool.
            for (const i of used) {
                if (i === t[0] || i === t[1] || i === t[2]) {
                    continue;
                }
                expect(inCircle(a, b, c, vertices[i])).toBeLessThan(1e-12);
            }
        }

        // Euler's formula for a triangulation of points in the plane:
        // numTriangles = 2 * numVertices - numHull - 2.
        const hull = triangulation.getHull();
        expect(triangles.length).toBe(2 * points.length - hull.length - 2);
    });

    it('produces a counterclockwise convex hull that contains every point', () => {
        const rnd = makeRandom(555);
        const points: Vector[] = [];
        for (let i = 0; i < 20; ++i) {
            points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
        }
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        for (const p of points) {
            triangulation.insert(p);
        }
        triangulation.finalizeTriangulation();

        const vertices = triangulation.getVertices();
        const hull = triangulation.getHull();
        expect(hull.length).toBeGreaterThanOrEqual(3);

        // Every point of the triangulation is inside or on the hull, and the
        // hull turns consistently counterclockwise.
        for (let i = 0; i < hull.length; ++i) {
            const a = vertices[hull[i]];
            const b = vertices[hull[(i + 1) % hull.length]];
            const c = vertices[hull[(i + 2) % hull.length]];
            expect(signedArea2(a, b, c)).toBeGreaterThan(-1e-12);
            for (const p of points) {
                expect(signedArea2(a, b, p)).toBeGreaterThan(-1e-12);
            }
        }
    });

    it('supports the containing-triangle search', () => {
        const points = [
            v2(0.1, 0.1), v2(0.9, 0.1), v2(0.9, 0.9), v2(0.1, 0.9), v2(0.5, 0.45)
        ];
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        for (const p of points) {
            triangulation.insert(p);
        }
        triangulation.finalizeTriangulation();

        const vertices = triangulation.getVertices();
        const info = new IncrementalDelaunay2SearchInfo();
        expect(info.initialTriangle).toBe(-1);

        // The centroid of each triangle is inside that triangle.
        const triangles = triangulation.getTriangles();
        for (let t = 0; t < triangles.length; ++t) {
            const tri = triangles[t];
            const centroid = v2(
                (vertices[tri[0]].get(0) + vertices[tri[1]].get(0)
                    + vertices[tri[2]].get(0)) / 3,
                (vertices[tri[0]].get(1) + vertices[tri[1]].get(1)
                    + vertices[tri[2]].get(1)) / 3);
            expect(triangulation.getContainingTriangle(centroid, info)).toBe(t);
            expect(info.numPath).toBeGreaterThan(0);
            expect(info.path[info.numPath - 1]).toBe(t);
            expect(info.finalTriangle).toBe(t);
        }

        // A point well outside the hull has no containing triangle.
        expect(triangulation.getContainingTriangle(v2(-5, -5), info)).toBe(-1);
    });

    it('exposes triangles and adjacencies consistently', () => {
        const rnd = makeRandom(777);
        const points: Vector[] = [];
        for (let i = 0; i < 16; ++i) {
            points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
        }
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        for (const p of points) {
            triangulation.insert(p);
        }
        triangulation.finalizeTriangulation();

        const numTriangles = triangulation.getNumTriangles();
        expect(triangulation.getTriangles().length).toBe(numTriangles);
        expect(triangulation.getAdjacencies().length).toBe(numTriangles);
        expect(triangulation.getTriangle(numTriangles)).toBeNull();
        expect(triangulation.getAdjacent(numTriangles)).toBeNull();
        expect(triangulation.getTriangle(-1)).toBeNull();

        for (let t = 0; t < numTriangles; ++t) {
            const tri = triangulation.getTriangle(t) as [number, number, number];
            const adj = triangulation.getAdjacent(t) as [number, number, number];
            expect(tri).toEqual([...triangulation.getTriangles()[t]]);

            for (let j = 0; j < 3; ++j) {
                const a = adj[j];
                if (a === -1) {
                    continue;
                }
                expect(a).toBeGreaterThanOrEqual(0);
                expect(a).toBeLessThan(numTriangles);

                // The adjacent triangle shares the edge <V[j], V[(j+1)%3]>
                // and lists this triangle as its own neighbor.
                const other = triangulation.getTriangle(a) as [number, number, number];
                const shared = [tri[j], tri[(j + 1) % 3]];
                expect(other.filter(v => shared.includes(v)).length).toBe(2);
                expect(triangulation.getAdjacent(a)).toContain(t);
            }
        }
    });

    it('removes points and reproduces the triangulation of the remaining points', () => {
        const rnd = makeRandom(20260901);
        const points: Vector[] = [];
        for (let i = 0; i < 20; ++i) {
            points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
        }

        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        const indices: number[] = [];
        for (const p of points) {
            indices.push(triangulation.insert(p));
        }

        // Remove three of the points, an interior one and two others.
        const removed = [3, 11, 17];
        for (const i of removed) {
            expect(triangulation.remove(points[i])).toBe(indices[i]);
        }
        triangulation.finalizeTriangulation();

        const survivors: Vector[] = [];
        const toInput = new Map<number, number>();
        for (let i = 0; i < points.length; ++i) {
            if (!removed.includes(i)) {
                toInput.set(indices[i], survivors.length);
                survivors.push(points[i]);
            }
        }

        const incremental = triangleSet(triangulation.getTriangles().map(t => [
            toInput.get(t[0]) as number,
            toInput.get(t[1]) as number,
            toInput.get(t[2]) as number
        ]));
        expect(incremental).toEqual(triangleSet(batchTriangles(survivors)));
    });

    it('handles interleaved insertion and removal', () => {
        const rnd = makeRandom(31337);
        const triangulation = new IncrementalDelaunay2(-1, -1, 1, 1);
        const live = new Map<number, Vector>();
        for (let k = 0; k < 40; ++k) {
            const p = v2(-0.9 + 1.8 * rnd(), -0.9 + 1.8 * rnd());
            const index = triangulation.insert(p);
            live.set(index, p);
            if (k % 5 === 4 && live.size > 4) {
                // Remove the oldest live point.
                const oldest = Math.min(...live.keys());
                expect(triangulation.remove(live.get(oldest) as Vector)).toBe(oldest);
                live.delete(oldest);
            }
        }
        triangulation.finalizeTriangulation();

        const survivors: Vector[] = [];
        const toInput = new Map<number, number>();
        for (const index of Array.from(live.keys()).sort((a, b) => a - b)) {
            toInput.set(index, survivors.length);
            survivors.push(live.get(index) as Vector);
        }

        const incremental = triangleSet(triangulation.getTriangles().map(t => [
            toInput.get(t[0]) as number,
            toInput.get(t[1]) as number,
            toInput.get(t[2]) as number
        ]));
        expect(incremental).toEqual(triangleSet(batchTriangles(survivors)));
    });

    it('handles collinear and duplicated inputs', () => {
        // Collinear points: the triangulation has no interior triangles, so
        // every triangle of the hull-free configuration is degenerate and is
        // rejected; only the supervertex triangles survive.
        const collinear = new IncrementalDelaunay2(0, 0, 1, 1);
        for (let i = 1; i <= 6; ++i) {
            collinear.insert(v2(i / 8, i / 8));
        }
        collinear.finalizeTriangulation();
        expect(collinear.getNumTriangles()).toBe(0);
        // There is no hull polygon for a degenerate triangulation; upstream
        // loops forever here, the port reports the condition.
        expect(() => collinear.getHull()).toThrow(/degenerate/);

        // Duplicates never create new vertices.
        const duplicates = new IncrementalDelaunay2(0, 0, 1, 1);
        const points = [v2(0.2, 0.2), v2(0.8, 0.2), v2(0.5, 0.8)];
        const first = points.map(p => duplicates.insert(p));
        const again = points.map(p => duplicates.insert(p));
        expect(again).toEqual(first);
        expect(duplicates.getNumVertices()).toBe(10);
        duplicates.finalizeTriangulation();
        expect(duplicates.getNumTriangles()).toBe(1);
        expect(triangleSet(duplicates.getTriangles()))
            .toEqual(triangleSet([[first[0], first[1], first[2]]]));

        // Collinear points plus one off the line give a fan.
        const fan = new IncrementalDelaunay2(0, 0, 1, 1);
        for (let i = 1; i <= 4; ++i) {
            fan.insert(v2(i / 8, 0.5));
        }
        fan.insert(v2(0.5, 0.75));
        fan.finalizeTriangulation();
        expect(fan.getNumTriangles()).toBe(3);
    });

    it('is final after finalizeTriangulation', () => {
        const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
        triangulation.insert(v2(0.2, 0.2));
        triangulation.insert(v2(0.8, 0.2));
        triangulation.insert(v2(0.5, 0.8));
        expect(triangulation.finalizeTriangulation()).toBe(true);

        // The second call is a no-op, and insertion/removal is disabled.
        expect(triangulation.finalizeTriangulation()).toBe(false);
        expect(triangulation.insert(v2(0.4, 0.4))).toBe(-1);
        expect(triangulation.remove(v2(0.2, 0.2))).toBe(-1);

        // Points outside the domain no longer throw in remove(), because the
        // rectangle is gone; they simply report -1.
        expect(triangulation.remove(v2(5, 5))).toBe(-1);
    });

    it('is deterministic across runs', () => {
        const rnd = makeRandom(2468);
        const points: Vector[] = [];
        for (let i = 0; i < 25; ++i) {
            points.push(v2(0.02 + 0.96 * rnd(), 0.02 + 0.96 * rnd()));
        }
        const run = () => {
            const triangulation = new IncrementalDelaunay2(0, 0, 1, 1);
            for (const p of points) {
                triangulation.insert(p);
            }
            triangulation.finalizeTriangulation();
            return {
                triangles: triangulation.getTriangles().map(t => [...t]),
                adjacencies: triangulation.getAdjacencies().map(a => [...a]),
                hull: triangulation.getHull()
            };
        };
        expect(run()).toEqual(run());
    });
});

// ---------------------------------------------------------------------------
// Verification pass (group V10). The generators use integer lattice points, so
// every predicate below can be evaluated exactly with bigint arithmetic: no
// tolerance appears anywhere in this block.
//
// Note that the epsilon defect of issue #391 (Delaunay2/Delaunay3 classify the
// dimension of the input with a hardcoded epsilon = 0 against a normalized
// frame) does not apply here: IncrementalDelaunay2 never uses IntrinsicsVector2
// and never classifies the dimension. Degenerate input simply produces a
// triangulation with no Delaunay triangles, which the collinear tests above
// pin.
// ---------------------------------------------------------------------------

// The domain of the incremental triangulation used by the properties. The
// lattice points are strictly inside it, as insert() requires.
const domainMin = 0;
const domainMax = 16;

type Exact2 = [bigint, bigint];

function exactPoint(v: Vector): Exact2 {
    // The generators emit integer coordinates, so BigInt() is exact.
    return [BigInt(v.get(0)), BigInt(v.get(1))];
}

function orientExact(a: Exact2, b: Exact2, c: Exact2): number {
    return orient2(a[0], a[1], b[0], b[1], c[0], c[1]);
}

function twiceAreaExact(a: Exact2, b: Exact2, c: Exact2): bigint {
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
}

// Remove duplicate points, keeping the first occurrence.
function uniquePoints(points: readonly Vector[]): Vector[] {
    const seen = new Set<string>();
    const unique: Vector[] = [];
    for (const p of points) {
        const key = `${p.get(0)},${p.get(1)}`;
        if (!seen.has(key)) {
            seen.add(key);
            unique.push(p);
        }
    }
    return unique;
}

// The strictly convex counterclockwise hull of the points (Andrew's monotone
// chain with the exact orientation predicate), which drops collinear points.
function strictHullExact(points: readonly Vector[]): Exact2[] {
    const sorted = points.map(exactPoint).sort((a, b) =>
        (a[0] !== b[0] ? (a[0] < b[0] ? -1 : 1)
            : (a[1] === b[1] ? 0 : (a[1] < b[1] ? -1 : 1))));
    const build = (input: Exact2[]): Exact2[] => {
        const chain: Exact2[] = [];
        for (const p of input) {
            while (chain.length >= 2
                && orientExact(chain[chain.length - 2],
                    chain[chain.length - 1], p) <= 0) {
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

// Does p lie on the segment [a, b] (endpoints included)?
function onSegmentExact(a: Exact2, b: Exact2, p: Exact2): boolean {
    if (orientExact(a, b, p) !== 0) {
        return false;
    }
    const withinX = (a[0] <= p[0] && p[0] <= b[0])
        || (b[0] <= p[0] && p[0] <= a[0]);
    const withinY = (a[1] <= p[1] && p[1] <= b[1])
        || (b[1] <= p[1] && p[1] <= a[1]);
    return withinX && withinY;
}

// The number of input points on the boundary of the convex hull, collinear
// boundary points included. This is the h of the triangulation count formula
// numTriangles = 2 * n - 2 - h.
function numBoundaryPoints(points: readonly Vector[]): number {
    const hull = strictHullExact(points);
    if (hull.length < 3) {
        return points.length;
    }
    let count = 0;
    for (const p of points.map(exactPoint)) {
        for (let i = 0; i < hull.length; ++i) {
            if (onSegmentExact(hull[i], hull[(i + 1) % hull.length], p)) {
                ++count;
                break;
            }
        }
    }
    return count;
}

function allCollinear(points: readonly Vector[]): boolean {
    const exact = points.map(exactPoint);
    for (let i = 2; i < exact.length; ++i) {
        if (orientExact(exact[0], exact[1], exact[i]) !== 0) {
            return false;
        }
    }
    return true;
}

// True when no four of the points are cocircular and no three are collinear,
// which is exactly the condition under which the Delaunay triangulation of the
// point set is unique.
function inGeneralPosition(points: readonly Vector[]): boolean {
    const exact = points.map(exactPoint);
    const n = exact.length;
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            for (let k = j + 1; k < n; ++k) {
                const sign = orientExact(exact[i], exact[j], exact[k]);
                if (sign === 0) {
                    return false;
                }
                const [a, b, c] = sign > 0
                    ? [exact[i], exact[j], exact[k]]
                    : [exact[i], exact[k], exact[j]];
                for (let l = 0; l < n; ++l) {
                    if (l === i || l === j || l === k) {
                        continue;
                    }
                    if (inCircle2(a[0], a[1], b[0], b[1], c[0], c[1],
                        exact[l][0], exact[l][1]) === 0) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
}

// A finalized incremental triangulation of the points, with the triangles
// expressed as indices into `points`.
interface Finalized {
    triangles: [number, number, number][];
    adjacencies: [number, number, number][];
    hull: number[];
}

function finalize(points: readonly Vector[],
    removeThenReinsert: number = -1): Finalized {
    const query = new IncrementalDelaunay2(domainMin, domainMin,
        domainMax, domainMax);
    const toInput = new Map<number, number>();
    for (let i = 0; i < points.length; ++i) {
        toInput.set(query.insert(points[i]), i);
    }
    if (removeThenReinsert >= 0) {
        const p = points[removeThenReinsert];
        expect(query.remove(p)).toBeGreaterThanOrEqual(0);
        toInput.set(query.insert(p), removeThenReinsert);
    }
    expect(query.finalizeTriangulation()).toBe(true);
    const map = (v: number): number => {
        const i = toInput.get(v);
        expect(i).not.toBeUndefined();
        return i as number;
    };
    return {
        triangles: query.getTriangles().map(t =>
            [map(t[0]), map(t[1]), map(t[2])] as [number, number, number]),
        adjacencies: query.getAdjacencies().map(a => [...a] as
            [number, number, number]),
        hull: query.getHull().map(map)
    };
}

// Check every combinatorial and geometric invariant of a Delaunay
// triangulation of `points`, exactly.
function expectDelaunay(points: readonly Vector[], result: Finalized): void {
    const exact = points.map(exactPoint);
    const n = points.length;

    // Every triangle is counterclockwise and uses distinct input points.
    let twiceArea = 0n;
    for (const t of result.triangles) {
        expect(new Set(t).size).toBe(3);
        const a = exact[t[0]], b = exact[t[1]], c = exact[t[2]];
        expect(orientExact(a, b, c)).toBe(1);
        twiceArea += twiceAreaExact(a, b, c);

        // The empty-circumcircle property: no other input point is strictly
        // inside the circumcircle of the triangle.
        for (let p = 0; p < n; ++p) {
            if (p === t[0] || p === t[1] || p === t[2]) {
                continue;
            }
            expect(inCircle2(a[0], a[1], b[0], b[1], c[0], c[1],
                exact[p][0], exact[p][1])).toBeLessThanOrEqual(0);
        }
    }

    // The triangles tile the convex hull, so their areas add up to the hull
    // area, and Euler's formula fixes their number.
    const hull = strictHullExact(points);
    expect(twiceArea).toBe(twiceSignedAreaExact(hull));
    const h = numBoundaryPoints(points);
    expect(result.triangles.length).toBe(2 * n - 2 - h);

    // Distinct triangles: no triangle is listed twice.
    const keys = new Set(result.triangles.map(t =>
        [...t].sort((x, y) => x - y).join(',')));
    expect(keys.size).toBe(result.triangles.length);

    // The adjacency array is symmetric and describes shared edges.
    for (let t = 0; t < result.triangles.length; ++t) {
        const v = result.triangles[t];
        const adj = result.adjacencies[t];
        for (let j = 0; j < 3; ++j) {
            const a = adj[j];
            if (a === -1) {
                continue;
            }
            expect(a).toBeGreaterThanOrEqual(0);
            expect(a).toBeLessThan(result.triangles.length);
            const w = result.triangles[a];
            // The shared edge is <v[j], v[j+1]> reversed in the neighbor.
            const edge = [v[j], v[(j + 1) % 3]];
            let found = -1;
            for (let k = 0; k < 3; ++k) {
                if (w[k] === edge[1] && w[(k + 1) % 3] === edge[0]) {
                    found = k;
                }
            }
            expect(found).toBeGreaterThanOrEqual(0);
            expect(result.adjacencies[a][found]).toBe(t);
        }
    }

    // The reported hull is the boundary of the triangulation: counterclockwise,
    // with every input point inside it or on it, and with as many vertices as
    // there are boundary points.
    expect(result.hull.length).toBe(h);
    const hullExact = result.hull.map(i => exact[i]);
    for (let i = 0; i < hullExact.length; ++i) {
        const a = hullExact[i];
        const b = hullExact[(i + 1) % hullExact.length];
        for (const p of exact) {
            // Counterclockwise hull: no point is strictly to the right of a
            // hull edge.
            expect(orientExact(a, b, p)).toBeGreaterThanOrEqual(0);
        }
    }
    expect(twiceSignedAreaExact(hullExact)).toBe(twiceSignedAreaExact(hull));
}

const latticeCloud = fc.array(latticeVector(2, 1, domainMax - 1),
    { minLength: 3, maxLength: 9 })
    .map(points => uniquePoints(points))
    .filter(points => points.length >= 3 && !allCollinear(points));

describe('IncrementalDelaunay2 verification', () => {
    it('produces an exact Delaunay triangulation of the inserted points', () => {
        check(latticeCloud, points => {
            expectDelaunay(points, finalize(points));
        }, 200);
    });

    it('is unchanged by removing and re-inserting a point', () => {
        check(fc.tuple(latticeCloud, fc.nat(100)), ([points, seed]) => {
            const index = seed % points.length;
            const roundTrip = finalize(points, index);
            expectDelaunay(points, roundTrip);
            if (inGeneralPosition(points)) {
                // The Delaunay triangulation is unique for a point set in
                // general position, so the round trip must reproduce it
                // exactly.
                const direct = finalize(points);
                const key = (r: Finalized): string =>
                    r.triangles.map(t => [...t].sort((a, b) => a - b).join(','))
                        .sort().join(' ');
                expect(key(roundTrip)).toBe(key(direct));
            }
        }, 100);
    });

    it('does not depend on the order of insertion', () => {
        check(fc.array(latticeVector(2, 1, domainMax - 1),
            { minLength: 3, maxLength: 8 })
            .map(points => uniquePoints(points))
            .filter(points => points.length >= 3 && !allCollinear(points))
            .chain(points => fc.tuple(fc.constant(points),
                fc.shuffledSubarray(points.map((_, i) => i),
                    { minLength: points.length, maxLength: points.length }))),
        ([points, order]) => {
            const shuffled = order.map(i => points[i]);
            const a = finalize(points);
            const b = finalize(shuffled);
            expectDelaunay(shuffled, b);
            if (inGeneralPosition(points)) {
                // Map the shuffled triangles back to the original indices.
                const keyOf = (triangles: [number, number, number][],
                    remap: (i: number) => number): string =>
                    triangles.map(t => t.map(remap).sort((x, y) => x - y)
                        .join(',')).sort().join(' ');
                expect(keyOf(b.triangles, i => order[i]))
                    .toBe(keyOf(a.triangles, i => i));
            }
        }, 100);
    });

    it('finds the triangle containing a lattice query point', () => {
        check(fc.tuple(latticeCloud, latticeVector(2, 1, domainMax - 1)),
            ([points, query]) => {
                const triangulation = new IncrementalDelaunay2(
                    domainMin, domainMin, domainMax, domainMax);
                const toInput = new Map<number, number>();
                for (let i = 0; i < points.length; ++i) {
                    toInput.set(triangulation.insert(points[i]), i);
                }
                expect(triangulation.finalizeTriangulation()).toBe(true);

                const info = new IncrementalDelaunay2SearchInfo();
                const t = triangulation.getContainingTriangle(query, info);
                const exactQuery = exactPoint(query);
                if (t === -1) {
                    // The query point is outside every triangle: it must be
                    // strictly outside the convex hull.
                    const hull = strictHullExact(points);
                    let outside = hull.length < 3;
                    for (let i = 0; i < hull.length && !outside; ++i) {
                        if (orientExact(hull[i], hull[(i + 1) % hull.length],
                            exactQuery) < 0) {
                            outside = true;
                        }
                    }
                    expect(outside).toBe(true);
                    return;
                }
                const tri = triangulation.getTriangle(t);
                expect(tri).not.toBeNull();
                const v = (tri as [number, number, number]).map(i =>
                    exactPoint(points[toInput.get(i) as number]));
                // The point is inside the triangle or on its boundary.
                for (let j = 0; j < 3; ++j) {
                    expect(orientExact(v[j], v[(j + 1) % 3], exactQuery))
                        .toBeGreaterThanOrEqual(0);
                }
                expect(info.numPath).toBeGreaterThan(0);
                expect(info.finalTriangle).toBe(t);
            }, 100);
    });

    it('removing every inserted point restores the rectangle triangulation',
        () => {
            check(latticeCloud, points => {
                const triangulation = new IncrementalDelaunay2(
                    domainMin, domainMin, domainMax, domainMax);
                for (const p of points) {
                    triangulation.insert(p);
                }
                for (const p of points) {
                    expect(triangulation.remove(p)).toBeGreaterThanOrEqual(0);
                }
                // Only the four rectangle corners remain, triangulated by two
                // triangles; the supervertex triangles are not reported.
                expect(triangulation.getNumTriangles()).toBe(2);
                expect(triangulation.getHull().length).toBe(4);
                // Removing a point that is no longer present is a no-op.
                expect(triangulation.remove(points[0])).toBe(-1);
            }, 60);
        });
});
