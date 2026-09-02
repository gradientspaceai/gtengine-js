import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2';
import {
    IncrementalDelaunay2, IncrementalDelaunay2SearchInfo
} from '../src/IncrementalDelaunay2';
import { Vector } from '../src/Vector';

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
