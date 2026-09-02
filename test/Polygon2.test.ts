import { describe, it, expect } from 'vitest';
import { Polygon2 } from '../src/Polygon2';
import { Vector, length, sub } from '../src/Vector';

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
