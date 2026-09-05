import { describe, it, expect } from 'vitest';
import { DistTriangle3OrientedBox3 } from '../src/DistTriangle3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Triangle } from '../src/Triangle.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistTriangle3CanonicalBox3 } from '../src/DistTriangle3CanonicalBox3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, latticeVector,
    rotationFrame, seededRandom, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

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

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTriangle3OrientedBox3.h.
// ---------------------------------------------------------------------------

describe('DistTriangle3OrientedBox3 verification', () => {
    const query = new DistTriangle3OrientedBox3();
    const cQuery = new DistTriangle3CanonicalBox3();

    const boxArb = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 }))
        .map(([c, R, e]) => OrientedBox.fromCenterAxisExtent(c, R,
            Vector.fromArray([e[0], e[1], e[2]])));

    const triArb = fc.tuple(latticeVector(3, -6, 6), latticeVector(3, -6, 6),
        latticeVector(3, -6, 6))
        .filter(([a, b, c]) => {
            const n = cross(sub(b, a), sub(c, a));
            return dot(n, n) > 4;
        })
        .map(([a, b, c]) => Triangle.fromVertices(a, b, c));

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const res = query.compute(t, b);
            expectClose(res.sqrDistance, res.distance * res.distance,
                1e-12, 1e-12);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(length(sub(res.closest[0], res.closest[1])),
                res.distance, 1e-6, 1e-8);

            // The barycentric coordinates are frame independent, so closest[0]
            // in world coordinates must be their combination of the world
            // vertices. A closest[0] left in box coordinates (or transformed
            // twice) fails this for every non-identity box frame.
            const bary = res.barycentric;
            expectClose(bary[0] + bary[1] + bary[2], 1, 1e-9, 1e-9);
            let rebuilt = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                expect(bary[i]).toBeGreaterThanOrEqual(-1e-8);
                rebuilt = add(rebuilt, mul(bary[i], t.v[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);

            const delta = sub(res.closest[1], b.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(b.axis[i], delta)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-8);
            }
        });
    });

    it('is the canonical-box query composed with the box frame', () => {
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const toBox = (p: Vector) => {
                const d = sub(p, b.center);
                return v(dot(b.axis[0], d), dot(b.axis[1], d),
                    dot(b.axis[2], d));
            };
            const fromBox = (p: Vector) => add(b.center,
                add(add(mul(p.values[0], b.axis[0]),
                    mul(p.values[1], b.axis[1])),
                    mul(p.values[2], b.axis[2])));
            const rc = cQuery.compute(Triangle.fromVertices(toBox(t.v[0]),
                toBox(t.v[1]), toBox(t.v[2])),
                CanonicalBox.fromExtent(b.extent));
            const res = query.compute(t, b);
            expectClose(res.distance, rc.distance, 1e-12, 1e-12);
            expectVectorClose(res.closest[0], fromBox(rc.closest[0]),
                1e-9, 1e-9);
            expectVectorClose(res.closest[1], fromBox(rc.closest[1]),
                1e-9, 1e-9);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(triArb, boxArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([t, b, R, tr]) => {
                const rot = (q: Vector) => add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2]));
                const xf = (q: Vector) => add(rot(q), tr);
                const r0 = query.compute(t, b);
                const r1 = query.compute(
                    Triangle.fromVertices(xf(t.v[0]), xf(t.v[1]), xf(t.v[2])),
                    OrientedBox.fromCenterAxisExtent(xf(b.center),
                        [rot(b.axis[0]), rot(b.axis[1]), rot(b.axis[2])],
                        b.extent));
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectClose(length(sub(r1.closest[0], r1.closest[1])),
                    r0.distance, 1e-7, 1e-7);
            });
    });

    it('never exceeds a barycentric sampling of the triangle', () => {
        const rng = seededRandom(0x77aa33cc);
        for (let k = 0; k < 25; ++k) {
            const w = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(w);
            const f = frame(w);
            const b = OrientedBox.fromCenterAxisExtent(v(0.5, -1, 0.25),
                [f[0], f[1], f[2]], v(1, 2, 0.5));
            const p = () => v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4);
            const t = Triangle.fromVertices(p(), p(), p());
            expect(query.compute(t, b).distance)
                .toBeLessThanOrEqual(bruteForce(t, b, 36) + 1e-9);
        }
    }, 30000);
});
