import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistTriangle3AlignedBox3 } from '../src/DistTriangle3AlignedBox3.js';
import { Triangle } from '../src/Triangle.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistTriangle3CanonicalBox3 } from '../src/DistTriangle3CanonicalBox3.js';
import { DistTriangle3OrientedBox3 } from '../src/DistTriangle3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, latticeVector,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
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

// Brute-force minimum over a barycentric grid on the triangle.
function bruteForce(triangle: Triangle, box: AlignedBox, n: number): number {
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

const unitBox = () => AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));

describe('DistTriangle3AlignedBox3', () => {
    it('reports the separation of a triangle above a box face', () => {
        const triangle = Triangle.fromVertices(
            v(0, 0, 4), v(0.5, 0, 4), v(0, 0.5, 4));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.sqrDistance).toBeCloseTo(9, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(4, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('finds a vertex as the closest triangle point', () => {
        // The triangle is far from the box except for the vertex (2,0,0),
        // whose closest box point is (1,0,0).
        const triangle = Triangle.fromVertices(
            v(2, 0, 0), v(6, 3, 0), v(6, -3, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.closest[0].values).toEqual([2, 0, 0]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.barycentric[0]).toBeCloseTo(1, 10);
        expect(result.barycentric[1]).toBeCloseTo(0, 10);
        expect(result.barycentric[2]).toBeCloseTo(0, 10);
    });

    it('reports zero distance for a triangle that pierces the box', () => {
        const triangle = Triangle.fromVertices(
            v(-3, -3, 0), v(3, -3, 0), v(0, 3, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBe(0);
        expect(length(sub(result.closest[0], result.closest[1])))
            .toBeCloseTo(0, 12);
    });

    it('reports zero distance for a triangle inside the box', () => {
        const triangle = Triangle.fromVertices(
            v(-0.5, -0.5, 0), v(0.5, -0.5, 0), v(0, 0.5, 0));
        const result = new DistTriangle3AlignedBox3()
            .compute(triangle, unitBox());
        expect(result.distance).toBe(0);
    });

    it('translates the closest points out of the canonical box frame', () => {
        const box = AlignedBox.fromMinMax(v(10, 20, 30), v(12, 24, 36));
        const triangle = Triangle.fromVertices(
            v(11, 22, 40), v(11.5, 22, 40), v(11, 22.5, 40));
        const result = new DistTriangle3AlignedBox3().compute(triangle, box);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(40, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(36, 12);
        expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 12);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 13572468;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistTriangle3AlignedBox3();
        for (let trial = 0; trial < 40; ++trial) {
            const pt = () => v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3);
            const triangle = Triangle.fromVertices(pt(), pt(), pt());
            const lo = v(rnd() - 2, rnd() - 2, rnd() - 2);
            const box = AlignedBox.fromMinMax(lo,
                add(lo, v(0.5 + rnd() * 2, 0.5 + rnd() * 2, 0.5 + rnd() * 2)));

            const result = query.compute(triangle, box);
            const brute = bruteForce(triangle, box, 40);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.15);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 8);

            // The barycentric coordinates are a convex combination that
            // reproduces the triangle point.
            const b = result.barycentric;
            expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
            for (const bi of b) {
                expect(bi).toBeGreaterThanOrEqual(-1e-9);
            }
            const fromBary = add(add(mul(b[0], triangle.v[0]),
                mul(b[1], triangle.v[1])), mul(b[2], triangle.v[2]));
            expect(length(sub(fromBary, result.closest[0])))
                .toBeCloseTo(0, 8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistTriangle3AlignedBox3.h.
// ---------------------------------------------------------------------------

describe('DistTriangle3AlignedBox3 verification', () => {
    const query = new DistTriangle3AlignedBox3();
    const cQuery = new DistTriangle3CanonicalBox3();
    const obQuery = new DistTriangle3OrientedBox3();

    const boxArb = fc.tuple(wellScaledVector(3, -4, 4),
        fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            v(c.values[0] - e[0], c.values[1] - e[1], c.values[2] - e[2]),
            v(c.values[0] + e[0], c.values[1] + e[1], c.values[2] + e[2])));

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
            const bary = res.barycentric;
            expectClose(bary[0] + bary[1] + bary[2], 1, 1e-9, 1e-9);
            let rebuilt = new Vector(3);
            for (let i = 0; i < 3; ++i) {
                expect(bary[i]).toBeGreaterThanOrEqual(-1e-8);
                rebuilt = add(rebuilt, mul(bary[i], t.v[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);
            for (let i = 0; i < 3; ++i) {
                expect(res.closest[1].values[i])
                    .toBeGreaterThanOrEqual(b.min.values[i] - 1e-8);
                expect(res.closest[1].values[i])
                    .toBeLessThanOrEqual(b.max.values[i] + 1e-8);
            }
        });
    });

    it('is the canonical-box query composed with the box translation', () => {
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const cf = b.getCenteredForm();
            const xt = Triangle.fromVertices(sub(t.v[0], cf.center),
                sub(t.v[1], cf.center), sub(t.v[2], cf.center));
            const rc = cQuery.compute(xt, CanonicalBox.fromExtent(cf.extent));
            const res = query.compute(t, b);
            expectClose(res.distance, rc.distance, 1e-12, 1e-12);
            expectVectorClose(res.closest[0], add(rc.closest[0], cf.center),
                1e-9, 1e-9);
            expectVectorClose(res.closest[1], add(rc.closest[1], cf.center),
                1e-9, 1e-9);
        });
    });

    it('agrees with the oriented-box query for an identity box frame', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        check(fc.tuple(triArb, boxArb), ([t, b]) => {
            const cf = b.getCenteredForm();
            const ob = OrientedBox.fromCenterAxisExtent(cf.center, axes,
                cf.extent);
            const r0 = query.compute(t, b);
            const r1 = obQuery.compute(t, ob);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectClose(length(sub(r1.closest[0], r1.closest[1])),
                r1.distance, 1e-6, 1e-8);
        });
    });

    it('is invariant under a common translation', () => {
        check(fc.tuple(triArb, boxArb, wellScaledVector(3, -6, 6)),
            ([t, b, tr]) => {
                const r0 = query.compute(t, b);
                const r1 = query.compute(
                    Triangle.fromVertices(add(t.v[0], tr), add(t.v[1], tr),
                        add(t.v[2], tr)),
                    AlignedBox.fromMinMax(add(b.min, tr), add(b.max, tr)));
                // Only the distance is compared: when the two objects touch or
                // several pairs are equidistant, the runs may name different
                // representatives of the same minimum.
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectClose(length(sub(r1.closest[0], r1.closest[1])),
                    r1.distance, 1e-7, 1e-7);
            });
    });

    it('never exceeds a barycentric sampling of the triangle', () => {
        const rng = seededRandom(0xc0ffee11);
        const b = AlignedBox.fromMinMax(v(-1, -2, -0.5), v(1, 2, 0.5));
        for (let k = 0; k < 30; ++k) {
            const p = () => v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4);
            const t = Triangle.fromVertices(p(), p(), p());
            expect(query.compute(t, b).distance)
                .toBeLessThanOrEqual(bruteForce(t, b, 36) + 1e-9);
        }
    }, 30000);

    it('handles a flat (zero-extent) box', () => {
        check(fc.tuple(triArb, boxArb, fc.integer({ min: 0, max: 2 })),
            ([t, b, k]) => {
                const min = b.min.clone(), max = b.max.clone();
                max.values[k] = min.values[k];
                const res = query.compute(t, AlignedBox.fromMinMax(min, max));
                expect(Number.isFinite(res.distance)).toBe(true);
                expect(res.sqrDistance).toBeGreaterThanOrEqual(0);
                expectClose(length(sub(res.closest[0], res.closest[1])),
                    res.distance, 1e-6, 1e-8);
            });
    });
});
