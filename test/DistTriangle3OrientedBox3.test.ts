import { describe, it, expect } from 'vitest';
import { DistTriangle3OrientedBox3 } from '../src/DistTriangle3OrientedBox3';
import { OrientedBox } from '../src/OrientedBox';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function frame(w: Vector): [Vector, Vector, Vector] {
    const n = w.clone();
    normalize(n);
    let u = Math.abs(n.values[0]) > 0.5
        ? v(-n.values[1], n.values[0], 0)
        : v(0, -n.values[2], n.values[1]);
    normalize(u);
    return [u, cross(n, u), n];
}

function distPointBox(p: Vector, box: OrientedBox): number {
    const delta = sub(p, box.center);
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const y = dot(delta, box.axis[i]);
        const e = box.extent.values[i];
        const d = y < -e ? -e - y : (y > e ? y - e : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

function bruteForce(triangle: Triangle, box: OrientedBox, n: number): number {
    let best = Infinity;
    for (let i = 0; i <= n; ++i) {
        for (let j = 0; i + j <= n; ++j) {
            const b0 = i / n, b1 = j / n, b2 = 1 - b0 - b1;
            const p = add(add(mul(b0, triangle.v[0]), mul(b1, triangle.v[1])),
                mul(b2, triangle.v[2]));
            const d = distPointBox(p, box);
            if (d < best) {
                best = d;
            }
        }
    }
    return best;
}

describe('DistTriangle3OrientedBox3', () => {
    it('matches the aligned case when the box axes are the coordinate axes',
        () => {
            const box = new OrientedBox(3);
            box.center = v(0, 0, 0);
            box.extent = v(1, 1, 1);
            const triangle = Triangle.fromVertices(
                v(0, 0, 4), v(0.5, 0, 4), v(0, 0.5, 4));
            const result =
                new DistTriangle3OrientedBox3().compute(triangle, box);
            expect(result.distance).toBeCloseTo(3, 12);
            expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
        });

    it('accounts for the box rotation', () => {
        // A cube of half-width 1 rotated 45 degrees about the z-axis reaches
        // x = sqrt(2), so a triangle vertex at x = 5 on the x-axis is
        // 5 - sqrt(2) away.
        const s = Math.SQRT1_2;
        const box = new OrientedBox(3);
        box.center = v(0, 0, 0);
        box.axis = [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)];
        box.extent = v(1, 1, 1);
        const triangle = Triangle.fromVertices(
            v(5, 0, 0), v(9, 4, 0), v(9, -4, 0));
        const result = new DistTriangle3OrientedBox3().compute(triangle, box);
        expect(result.distance).toBeCloseTo(5 - Math.SQRT2, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT2, 12);
        expect(result.barycentric[0]).toBeCloseTo(1, 10);
    });

    it('reports zero distance when the triangle pierces the box', () => {
        const s = Math.SQRT1_2;
        const box = new OrientedBox(3);
        box.center = v(1, 2, 3);
        box.axis = [v(s, 0, s), v(0, 1, 0), v(-s, 0, s)];
        box.extent = v(1, 2, 0.5);
        const triangle = Triangle.fromVertices(
            v(-3, -1, 3), v(5, -1, 3), v(1, 6, 3));
        const result = new DistTriangle3OrientedBox3().compute(triangle, box);
        expect(result.distance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('is invariant under a rigid motion of both primitives', () => {
        const [u0, u1, u2] = frame(v(1, 2, 3));
        const box = new OrientedBox(3);
        box.center = v(0.5, -0.25, 0.75);
        box.axis = [u0, u1, u2];
        box.extent = v(1, 0.5, 2);
        const triangle = Triangle.fromVertices(
            v(4, 1, -2), v(5, 3, -1), v(3, 2, 2));
        const result = new DistTriangle3OrientedBox3().compute(triangle, box);

        const t = v(-3, 7, 11);
        const rot = (p: Vector) => add(t, add(add(
            mul(p.values[0], u0), mul(p.values[1], u1)),
            mul(p.values[2], u2)));
        const box2 = new OrientedBox(3);
        box2.center = rot(box.center);
        box2.axis = box.axis.map(a => sub(rot(a), t));
        box2.extent = box.extent.clone();
        const triangle2 = Triangle.fromVertices(rot(triangle.v[0]),
            rot(triangle.v[1]), rot(triangle.v[2]));
        const result2 = new DistTriangle3OrientedBox3()
            .compute(triangle2, box2);
        expect(result2.distance).toBeCloseTo(result.distance, 10);
        expect(result2.barycentric[0]).toBeCloseTo(result.barycentric[0], 8);
        expect(result2.barycentric[1]).toBeCloseTo(result.barycentric[1], 8);
        expect(result2.barycentric[2]).toBeCloseTo(result.barycentric[2], 8);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 8675309;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistTriangle3OrientedBox3();
        for (let trial = 0; trial < 40; ++trial) {
            const pt = () => v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3);
            const triangle = Triangle.fromVertices(pt(), pt(), pt());

            const bw = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(bw) < 1e-3) {
                continue;
            }
            const [b0, b1, b2] = frame(bw);
            const box = new OrientedBox(3);
            box.center = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            box.axis = [b0, b1, b2];
            box.extent = v(0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd());

            const result = query.compute(triangle, box);
            const brute = bruteForce(triangle, box, 40);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.2);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 8);

            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            const fromBary = add(add(mul(b[0], triangle.v[0]),
                mul(b[1], triangle.v[1])), mul(b[2], triangle.v[2]));
            expect(length(sub(fromBary, result.closest[0])))
                .toBeCloseTo(0, 8);
        }
    });
});
