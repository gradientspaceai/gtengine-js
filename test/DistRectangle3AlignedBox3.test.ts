import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistRectangle3AlignedBox3 } from '../src/DistRectangle3AlignedBox3.js';
import { Rectangle } from '../src/Rectangle.js';
import { Vector, add, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function rect(center: Vector, a0: Vector, a1: Vector, e0: number, e1: number):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(center, [a0, a1],
        Vector.fromArray([e0, e1]));
}

// The distance from a point to a solid aligned box, computed by clamping.
function distPointBox(p: Vector, box: AlignedBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const x = p.values[i];
        const lo = box.min.values[i], hi = box.max.values[i];
        const d = x < lo ? lo - x : (x > hi ? x - hi : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

// Brute-force minimum of the distance from rectangle points to the box.
function bruteForce(rectangle: Rectangle, box: AlignedBox, n: number): number {
    const e0 = rectangle.extent.values[0], e1 = rectangle.extent.values[1];
    let best = Infinity;
    for (let i = 0; i <= n; ++i) {
        const s0 = -e0 + 2 * e0 * i / n;
        const p0 = add(rectangle.center, mul(s0, rectangle.axis[0]));
        for (let j = 0; j <= n; ++j) {
            const s1 = -e1 + 2 * e1 * j / n;
            const d = distPointBox(add(p0, mul(s1, rectangle.axis[1])), box);
            if (d < best) {
                best = d;
            }
        }
    }
    return best;
}

const unitBox = () => AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));

describe('DistRectangle3AlignedBox3', () => {
    it('reports the separation of a rectangle parallel to a box face', () => {
        // The rectangle lies in the plane z = 3, well inside the box extents
        // in x and y, so the closest box point is directly below it.
        const rectangle = rect(v(0, 0, 3), v(1, 0, 0), v(0, 1, 0), 0.5, 0.5);
        const result = new DistRectangle3AlignedBox3()
            .compute(rectangle, unitBox());
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.sqrDistance).toBeCloseTo(4, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports zero distance for a rectangle that cuts through the box', () => {
        const rectangle = rect(v(0, 0, 0), v(1, 0, 0), v(0, 1, 0), 3, 3);
        const result = new DistRectangle3AlignedBox3()
            .compute(rectangle, unitBox());
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('reports zero distance for a rectangle inside the box', () => {
        const rectangle = rect(v(0, 0, 0), v(1, 0, 0), v(0, 1, 0), 0.25, 0.25);
        const result = new DistRectangle3AlignedBox3()
            .compute(rectangle, unitBox());
        expect(result.distance).toBe(0);
    });

    it('touches the box exactly when the rectangle meets a corner', () => {
        // The rectangle is in the plane x + y + z = 3, which touches the
        // corner (1,1,1) of the box.
        const normal = v(1, 1, 1);
        normalize(normal);
        const a0 = v(1, -1, 0);
        normalize(a0);
        const a1 = cross(normal, a0);
        const rectangle = rect(mul(Math.sqrt(3), normal), a0, a1, 2, 2);
        const result = new DistRectangle3AlignedBox3()
            .compute(rectangle, unitBox());
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('translates the closest points out of the canonical box frame', () => {
        // A box that is not centered at the origin exercises the translation
        // that this file adds on top of the canonical-box query.
        const box = AlignedBox.fromMinMax(v(10, 20, 30), v(12, 24, 36));
        const rectangle = rect(v(11, 22, 40), v(1, 0, 0), v(0, 1, 0), 0.5, 0.5);
        const result = new DistRectangle3AlignedBox3().compute(rectangle, box);
        expect(result.distance).toBeCloseTo(4, 12);
        // The rectangle projects into the top face of the box, so the pair of
        // closest points is not unique; only the plane each lies in is.
        expect(result.closest[0].values[2]).toBeCloseTo(40, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(36, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(
            result.closest[0].values[0], 12);
        expect(result.closest[1].values[1]).toBeCloseTo(
            result.closest[0].values[1], 12);
        expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 12);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 987654321;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistRectangle3AlignedBox3();
        for (let trial = 0; trial < 40; ++trial) {
            const normal = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(normal) < 1e-3) {
                continue;
            }
            normalize(normal);
            let a0 = Math.abs(normal.values[0]) > 0.5
                ? v(-normal.values[1], normal.values[0], 0)
                : v(0, -normal.values[2], normal.values[1]);
            normalize(a0);
            const a1 = cross(normal, a0);
            const rectangle = rect(
                v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3), a0, a1,
                0.2 + rnd(), 0.2 + rnd());
            const lo = v(rnd() - 2, rnd() - 2, rnd() - 2);
            const box = AlignedBox.fromMinMax(lo,
                add(lo, v(0.5 + rnd() * 2, 0.5 + rnd() * 2, 0.5 + rnd() * 2)));

            const result = query.compute(rectangle, box);
            const brute = bruteForce(rectangle, box, 60);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.05);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);

            // The closest points are consistent with the distance and with
            // the reported W-coordinates of the rectangle point.
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            const fromCartesian = add(rectangle.center,
                add(mul(result.cartesian[0], rectangle.axis[0]),
                    mul(result.cartesian[1], rectangle.axis[1])));
            expect(length(sub(fromCartesian, result.closest[0])))
                .toBeCloseTo(0, 8);
            expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 8);
        }
    });
});
