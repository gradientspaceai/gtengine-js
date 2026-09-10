import { describe, it, expect } from 'vitest';
import { Polygon2 } from '../src/Polygon2.js';
import { Vector, length, sub } from '../src/Vector.js';
import { IntrSegment2Segment2TI } from '../src/IntrSegment2Segment2.js';
import { Segment } from '../src/Segment.js';
import { check, expectClose, fc, finite } from './helpers/arbitraries.js';

function V(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// The shoelace area, computed independently of the class.
function shoelaceArea(pool: readonly Vector[],
    indices: readonly number[]): number {
    let sum = 0;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
        const a = pool[indices[i]];
        const b = pool[indices[(i + 1) % n]];
        sum += a.values[0] * b.values[1] - b.values[0] * a.values[1];
    }
    return Math.abs(0.5 * sum);
}

function perimeter(pool: readonly Vector[],
    indices: readonly number[]): number {
    let sum = 0;
    const n = indices.length;
    for (let i = 0; i < n; ++i) {
        sum += length(sub(pool[indices[(i + 1) % n]], pool[indices[i]]));
    }
    return sum;
}

// The unit square, counterclockwise.
const squarePool = [V(0, 0), V(1, 0), V(1, 1), V(0, 1)];
const squareIndices = [0, 1, 2, 3];

describe('Polygon2', () => {
    it('rejects invalid input in the constructor', () => {
        // Fewer than three indices.
        const p0 = new Polygon2(squarePool, [0, 1], true);
        expect(p0.isValid()).toBe(false);
        expect(p0.counterClockwise()).toBe(false);
        expect(p0.getVertexPool()).toBeNull();
        expect(p0.getIndices()).toEqual([]);
        expect(p0.getVertices()).toEqual([]);
        // The queries on a failed polygon return the neutral values.
        expect(p0.computeVertexAverage().values).toEqual([0, 0]);
        expect(p0.computePerimeterLength()).toBe(0);
        expect(p0.computeArea()).toBe(0);
        expect(p0.isSimple()).toBe(false);
        expect(p0.isConvex()).toBe(false);

        // A duplicated index (the polygon is not simple).
        const p1 = new Polygon2(squarePool, [0, 1, 2, 1], true);
        expect(p1.isValid()).toBe(false);

        // A null vertex pool or null indices.
        expect(new Polygon2(null, squareIndices, true).isValid()).toBe(false);
        expect(new Polygon2(squarePool, null, true).isValid()).toBe(false);
        expect(new Polygon2([], squareIndices, true).isValid()).toBe(false);
    });

    it('exposes the pool, sorted vertices and polygon-ordered indices', () => {
        // A polygon that uses a subset of a larger pool, in an order that is
        // not sorted.
        const pool = [V(9, 9), V(0, 0), V(2, 0), V(2, 2), V(0, 2), V(-9, -9)];
        const indices = [4, 1, 2, 3];
        const polygon = new Polygon2(pool, indices, true);
        expect(polygon.isValid()).toBe(true);
        expect(polygon.getVertexPool()).toBe(pool);
        expect(polygon.getIndices()).toEqual([4, 1, 2, 3]);
        // The vertex set is sorted, as the upstream std::set is.
        expect(polygon.getVertices()).toEqual([1, 2, 3, 4]);
        expect(polygon.counterClockwise()).toBe(true);

        // The indices are copied, so mutating the caller's array does not
        // change the polygon.
        indices[0] = 0;
        expect(polygon.getIndices()).toEqual([4, 1, 2, 3]);

        // The square with corners (0,0) and (2,2).
        expect(polygon.computeVertexAverage().values).toEqual([1, 1]);
        expect(polygon.computeArea()).toBeCloseTo(4, 12);
        expect(polygon.computePerimeterLength()).toBeCloseTo(8, 12);
    });

    it('computes area, perimeter and average for the unit square', () => {
        const polygon = new Polygon2(squarePool, squareIndices, true);
        expect(polygon.computeArea()).toBeCloseTo(1, 15);
        expect(polygon.computePerimeterLength()).toBeCloseTo(4, 15);
        expect(polygon.computeVertexAverage().values).toEqual([0.5, 0.5]);

        // The area is orientation independent (it is an absolute value).
        const reversed = new Polygon2(squarePool, [3, 2, 1, 0], false);
        expect(reversed.computeArea()).toBeCloseTo(1, 15);
        expect(reversed.computePerimeterLength()).toBeCloseTo(4, 15);
    });

    it('classifies triangles as simple and convex without further work', () => {
        const pool = [V(0, 0), V(1, 0), V(0, 1)];
        const ccw = new Polygon2(pool, [0, 1, 2], true);
        expect(ccw.isSimple()).toBe(true);
        expect(ccw.isConvex()).toBe(true);
        // The three-index shortcut ignores the orientation flag entirely.
        const cw = new Polygon2(pool, [0, 1, 2], false);
        expect(cw.isConvex()).toBe(true);
        expect(ccw.computeArea()).toBeCloseTo(0.5, 15);
    });

    it('detects simplicity and convexity', () => {
        // The unit square is simple and convex when the orientation flag
        // matches the winding.
        expect(new Polygon2(squarePool, squareIndices, true).isSimple())
            .toBe(true);
        expect(new Polygon2(squarePool, squareIndices, true).isConvex())
            .toBe(true);
        // The same square declared clockwise: the convexity test uses the
        // sign from the flag, so it reports non-convex.
        expect(new Polygon2(squarePool, squareIndices, false).isConvex())
            .toBe(false);
        // The reversed square declared clockwise is convex.
        expect(new Polygon2(squarePool, [3, 2, 1, 0], false).isConvex())
            .toBe(true);

        // A bowtie: the edges (1,0)-(0,1) and (1,1)-(0,0) cross.
        const bowtiePool = [V(0, 0), V(1, 0), V(0, 1), V(1, 1)];
        const bowtie = new Polygon2(bowtiePool, [0, 1, 2, 3], true);
        expect(bowtie.isSimple()).toBe(false);
        expect(bowtie.isConvex()).toBe(false);

        // An L-shaped hexagon: simple but not convex.
        const lPool = [V(0, 0), V(2, 0), V(2, 1), V(1, 1), V(1, 2), V(0, 2)];
        const lIndices = [0, 1, 2, 3, 4, 5];
        const lShape = new Polygon2(lPool, lIndices, true);
        expect(lShape.isSimple()).toBe(true);
        expect(lShape.isConvex()).toBe(false);
        expect(lShape.computeArea()).toBeCloseTo(3, 12);
        expect(lShape.computePerimeterLength()).toBeCloseTo(8, 12);
        expect(lShape.computeVertexAverage().values[0]).toBeCloseTo(1, 12);
        expect(lShape.computeVertexAverage().values[1]).toBeCloseTo(1, 12);

        // A regular hexagon is simple and convex.
        const hexPool: Vector[] = [];
        for (let i = 0; i < 6; ++i) {
            const a = i * Math.PI / 3;
            hexPool.push(V(Math.cos(a), Math.sin(a)));
        }
        const hexagon = new Polygon2(hexPool, [0, 1, 2, 3, 4, 5], true);
        expect(hexagon.isSimple()).toBe(true);
        expect(hexagon.isConvex()).toBe(true);
        expect(hexagon.computeArea()).toBeCloseTo(1.5 * Math.sqrt(3), 12);
        expect(hexagon.computePerimeterLength()).toBeCloseTo(6, 12);
    });

    it('agrees with independent formulas on random convex polygons', () => {
        let seed = 1234567;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 40; ++trial) {
            // Sample angles in increasing order to build a convex polygon
            // inscribed in a circle, which is simple and convex.
            const n = 4 + Math.floor(6 * rand());
            const angles: number[] = [];
            for (let i = 0; i < n; ++i) {
                angles.push(2 * Math.PI * (i + 0.15 + 0.7 * rand()) / n);
            }
            const radius = 0.5 + 2 * rand();
            const cx = 4 * rand() - 2;
            const cy = 4 * rand() - 2;
            const pool = angles.map(a =>
                V(cx + radius * Math.cos(a), cy + radius * Math.sin(a)));
            const indices = pool.map((_, i) => i);

            const polygon = new Polygon2(pool, indices, true);
            expect(polygon.isValid()).toBe(true);
            expect(polygon.computeArea())
                .toBeCloseTo(shoelaceArea(pool, indices), 10);
            expect(polygon.computePerimeterLength())
                .toBeCloseTo(perimeter(pool, indices), 10);

            let sx = 0;
            let sy = 0;
            for (const p of pool) {
                sx += p.values[0];
                sy += p.values[1];
            }
            const average = polygon.computeVertexAverage();
            expect(average.values[0]).toBeCloseTo(sx / n, 10);
            expect(average.values[1]).toBeCloseTo(sy / n, 10);

            expect(polygon.isSimple()).toBe(true);
            expect(polygon.isConvex()).toBe(true);

            // Swapping two nonadjacent vertices makes the polygon
            // self-intersecting.
            if (n >= 5) {
                const swapped = indices.slice();
                const tmp = swapped[1];
                swapped[1] = swapped[3];
                swapped[3] = tmp;
                const bad = new Polygon2(pool, swapped, true);
                expect(bad.isSimple()).toBe(false);
                expect(bad.isConvex()).toBe(false);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// V43 verification: the geometric queries against independent formulas, the
// simple/convex classification on generated star-shaped and convex polygons,
// and the failure semantics of the constructor.
// ---------------------------------------------------------------------------

describe('Polygon2 verification', () => {
    // n points at strictly increasing angles around the origin. The angular
    // jitter is bounded by 0.4 * (2*pi/n), so the angles stay ordered and
    // separated and the polygon (traversed in index order) is star-shaped
    // about the origin, hence simple.
    const starPolygon = (minN: number, maxN: number) =>
        fc.integer({ min: minN, max: maxN }).chain(n =>
            fc.tuple(fc.constant(n),
                fc.array(finite(-0.4, 0.4), { minLength: n, maxLength: n }),
                fc.array(finite(1, 4), { minLength: n, maxLength: n })))
            .map(([n, jitter, radii]) => {
                const step = (2 * Math.PI) / n;
                const pool: Vector[] = [];
                for (let i = 0; i < n; ++i) {
                    const a = (i + jitter[i]) * step;
                    pool.push(Vector.fromArray(
                        [radii[i] * Math.cos(a), radii[i] * Math.sin(a)]));
                }
                return pool;
            });

    // The same construction with all radii equal: the vertices lie on a
    // circle at increasing angles, so the polygon is convex.
    const convexPolygon = (minN: number, maxN: number) =>
        fc.tuple(fc.integer({ min: minN, max: maxN }), finite(1, 4))
            .chain(([n, radius]) =>
                fc.tuple(fc.constant(n), fc.constant(radius),
                    fc.array(finite(-0.4, 0.4),
                        { minLength: n, maxLength: n })))
            .map(([n, radius, jitter]) => {
                const step = (2 * Math.PI) / n;
                const pool: Vector[] = [];
                for (let i = 0; i < n; ++i) {
                    const a = (i + jitter[i]) * step;
                    pool.push(Vector.fromArray(
                        [radius * Math.cos(a), radius * Math.sin(a)]));
                }
                return pool;
            });

    const identity = (n: number) => Array.from({ length: n }, (_, i) => i);

    it('computes area, perimeter and average by the standard formulas', () => {
        check(starPolygon(3, 12), pool => {
            const indices = identity(pool.length);
            const polygon = new Polygon2(pool, indices, true);
            expect(polygon.isValid()).toBe(true);
            // The area matches the shoelace formula. The two summations
            // group the terms differently, so the tolerance is relative to
            // the area scale (radii are at most 4).
            expectClose(polygon.computeArea(), shoelaceArea(pool, indices),
                1e-12, 1e-12);
            expectClose(polygon.computePerimeterLength(),
                perimeter(pool, indices), 1e-12, 1e-12);
            let sx = 0, sy = 0;
            for (const v of pool) { sx += v.get(0); sy += v.get(1); }
            const average = polygon.computeVertexAverage();
            expectClose(average.get(0), sx / pool.length, 1e-12, 1e-12);
            expectClose(average.get(1), sy / pool.length, 1e-12, 1e-12);
        });
    });

    it('reports an orientation-independent absolute area', () => {
        check(starPolygon(3, 10), pool => {
            const indices = identity(pool.length);
            const forward = new Polygon2(pool, indices, true);
            const reversed = new Polygon2(pool, [...indices].reverse(), false);
            // ComputeArea takes std::fabs, so reversing the traversal (which
            // negates the signed area) does not change the result.
            expectClose(forward.computeArea(), reversed.computeArea(),
                1e-12, 1e-12);
            expectClose(forward.computePerimeterLength(),
                reversed.computePerimeterLength(), 1e-12, 1e-12);
            // A cyclic rotation of the indices is the same polygon.
            const rotated = [...indices.slice(1), indices[0]];
            expectClose(new Polygon2(pool, rotated, true).computeArea(),
                forward.computeArea(), 1e-12, 1e-12);
        });
    });

    it('classifies star-shaped polygons as simple', () => {
        check(starPolygon(4, 10), pool => {
            const polygon = new Polygon2(pool, identity(pool.length), true);
            expect(polygon.isSimple()).toBe(true);
        }, 100);
    });

    it('classifies circle polygons as convex in their own orientation', () => {
        check(convexPolygon(4, 10), pool => {
            const indices = identity(pool.length);
            // The angles increase, so the traversal is counterclockwise.
            expect(new Polygon2(pool, indices, true).isConvex()).toBe(true);
            // With the orientation flag reversed, every turn has the wrong
            // sign, so IsConvexInternal rejects it.
            expect(new Polygon2(pool, indices, false).isConvex()).toBe(false);
            // Reversing the traversal makes it clockwise.
            const rev = [...indices].reverse();
            expect(new Polygon2(pool, rev, false).isConvex()).toBe(true);
            expect(new Polygon2(pool, rev, true).isConvex()).toBe(false);
            // Convex implies simple.
            expect(new Polygon2(pool, indices, true).isSimple()).toBe(true);
        }, 100);
    });

    it('matches an all-pairs simplicity test, including swapped polygons',
        () => {
            // IsSimpleInternal does not iterate over every non-adjacent pair
            // of edges: its inner loop runs over the numeric range
            // [(i0+2) % n, (i0-2+n) % n] and is empty when the range wraps.
            // Cross-check it against an exhaustive all-pairs test so a missed
            // pair would show up here.
            const bruteForceSimple = (pool: readonly Vector[],
                indices: readonly number[]): boolean => {
                const n = indices.length;
                const query = new IntrSegment2Segment2TI();
                for (let i = 0; i < n; ++i) {
                    const s0 = Segment.fromEndpoints(pool[indices[i]],
                        pool[indices[(i + 1) % n]]);
                    for (let j = i + 1; j < n; ++j) {
                        // Skip the two adjacent pairs (they always share an
                        // endpoint).
                        if (j === i + 1 || (i === 0 && j === n - 1)) {
                            continue;
                        }
                        const s1 = Segment.fromEndpoints(pool[indices[j]],
                            pool[indices[(j + 1) % n]]);
                        if (query.test(s0, s1).intersect) {
                            return false;
                        }
                    }
                }
                return true;
            };

            check(fc.tuple(convexPolygon(5, 9), fc.nat(), fc.boolean()),
                ([pool, raw, swap]) => {
                    const n = pool.length;
                    const indices = identity(n);
                    if (swap) {
                        // Swapping two non-adjacent vertices of a polygon
                        // whose points are in convex position creates a
                        // crossing.
                        const j = 2 + (raw % (n - 3));
                        [indices[0], indices[j]] = [indices[j], indices[0]];
                    }
                    const polygon = new Polygon2(pool, indices, true);
                    const expected = bruteForceSimple(pool, indices);
                    expect(polygon.isSimple()).toBe(expected);
                    if (swap) {
                        // The crossing must actually be there.
                        expect(expected).toBe(false);
                        expect(polygon.isConvex()).toBe(false);
                    }
                }, 100);

            // The same cross-check on star-shaped polygons, which are simple.
            check(starPolygon(4, 9), pool => {
                const indices = identity(pool.length);
                const polygon = new Polygon2(pool, indices, true);
                expect(polygon.isSimple())
                    .toBe(bruteForceSimple(pool, indices));
            }, 100);
        });

    it('fails construction for short, duplicated or empty input', () => {
        check(fc.tuple(starPolygon(3, 8), fc.nat()), ([pool, raw]) => {
            const n = pool.length;
            const indices = identity(n);

            // Fewer than three indices.
            const short = new Polygon2(pool, indices.slice(0, 2), true);
            expect(short.isValid()).toBe(false);
            expect(short.counterClockwise()).toBe(false);
            expect(short.getIndices().length).toBe(0);
            expect(short.getVertices().length).toBe(0);
            // The queries of an invalid polygon return the zero values.
            expect(short.computeArea()).toBe(0);
            expect(short.computePerimeterLength()).toBe(0);
            expect(short.computeVertexAverage().values).toEqual([0, 0]);
            expect(short.isSimple()).toBe(false);
            expect(short.isConvex()).toBe(false);

            // A duplicated index.
            const dup = [...indices];
            dup[raw % n] = dup[(raw + 1) % n];
            expect(new Polygon2(pool, dup, true).isValid()).toBe(false);

            // A null or empty vertex pool. Upstream can only test the pointer
            // for null; the port also rejects an empty array.
            expect(new Polygon2(null, indices, true).isValid()).toBe(false);
            expect(new Polygon2(pool, null, true).isValid()).toBe(false);
            expect(new Polygon2([], indices, true).isValid()).toBe(false);
        });
    });

    it('stores sorted unique vertices and a copy of the indices', () => {
        check(starPolygon(3, 10), pool => {
            const indices = identity(pool.length);
            const rotated = [...indices.slice(2), ...indices.slice(0, 2)];
            const polygon = new Polygon2(pool, rotated, true);
            // getVertices is the sorted set (the std::set iteration order).
            expect([...polygon.getVertices()]).toEqual(indices);
            // getIndices keeps the polygon order.
            expect([...polygon.getIndices()]).toEqual(rotated);
            // The indices were copied, not aliased.
            rotated[0] = 12345;
            expect(polygon.getIndices()[0]).not.toBe(12345);
            // The vertex pool is referenced, not copied (upstream keeps the
            // caller's pointer).
            expect(polygon.getVertexPool()).toBe(pool);
        });
    });
});
