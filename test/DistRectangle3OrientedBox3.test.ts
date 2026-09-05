import { describe, it, expect } from 'vitest';
import { DistRectangle3OrientedBox3 } from '../src/DistRectangle3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Rectangle } from '../src/Rectangle.js';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistRectangle3CanonicalBox3 } from '../src/DistRectangle3CanonicalBox3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function rect(center: Vector, a0: Vector, a1: Vector, e0: number, e1: number):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(center, [a0, a1],
        Vector.fromArray([e0, e1]));
}

// A right-handed orthonormal frame whose third axis is the given direction.
function frame(w: Vector): [Vector, Vector, Vector] {
    const n = w.clone();
    normalize(n);
    let u = Math.abs(n.values[0]) > 0.5
        ? v(-n.values[1], n.values[0], 0)
        : v(0, -n.values[2], n.values[1]);
    normalize(u);
    return [u, cross(n, u), n];
}

// The distance from a point to a solid oriented box, computed by clamping in
// the box frame.
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

function bruteForce(rectangle: Rectangle, box: OrientedBox, n: number): number {
    const e0 = rectangle.extent.values[0], e1 = rectangle.extent.values[1];
    let best = Infinity;
    for (let i = 0; i <= n; ++i) {
        const p0 = add(rectangle.center,
            mul(-e0 + 2 * e0 * i / n, rectangle.axis[0]));
        for (let j = 0; j <= n; ++j) {
            const d = distPointBox(
                add(p0, mul(-e1 + 2 * e1 * j / n, rectangle.axis[1])), box);
            if (d < best) {
                best = d;
            }
        }
    }
    return best;
}

describe('DistRectangle3OrientedBox3', () => {
    it('matches the aligned case when the box axes are the coordinate axes',
        () => {
            const box = new OrientedBox(3);
            box.center = v(0, 0, 0);
            box.extent = v(1, 1, 1);
            const rectangle =
                rect(v(0, 0, 3), v(1, 0, 0), v(0, 1, 0), 0.5, 0.5);
            const result =
                new DistRectangle3OrientedBox3().compute(rectangle, box);
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
        });

    it('accounts for the box rotation', () => {
        // A cube of half-width 1 rotated 45 degrees about the z-axis has its
        // topmost x-extent at sqrt(2), so a rectangle at x = 5 is sqrt(2)
        // away along x.
        const s = Math.SQRT1_2;
        const box = new OrientedBox(3);
        box.center = v(0, 0, 0);
        box.axis = [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)];
        box.extent = v(1, 1, 1);
        const rectangle = rect(v(5, 0, 0), v(0, 1, 0), v(0, 0, 1), 0.5, 0.5);
        const result = new DistRectangle3OrientedBox3()
            .compute(rectangle, box);
        expect(result.distance).toBeCloseTo(5 - Math.SQRT2, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(Math.SQRT2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
    });

    it('reports zero distance when the rectangle intersects the box', () => {
        const s = Math.SQRT1_2;
        const box = new OrientedBox(3);
        box.center = v(1, 2, 3);
        box.axis = [v(s, 0, s), v(0, 1, 0), v(-s, 0, s)];
        box.extent = v(1, 2, 0.5);
        const rectangle = rect(v(1, 2, 3), v(1, 0, 0), v(0, 1, 0), 4, 4);
        const result = new DistRectangle3OrientedBox3()
            .compute(rectangle, box);
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
        const rectangle = rect(v(4, 1, -2), v(1, 0, 0), v(0, 0, 1), 1.5, 0.75);
        const result = new DistRectangle3OrientedBox3()
            .compute(rectangle, box);

        // Rotate everything by the same frame and translate; the distance
        // must not change.
        const t = v(-3, 7, 11);
        const rot = (p: Vector) => add(t, add(add(
            mul(p.values[0], u0), mul(p.values[1], u1)),
            mul(p.values[2], u2)));
        const box2 = new OrientedBox(3);
        box2.center = rot(box.center);
        box2.axis = box.axis.map(a => sub(rot(a), t));
        box2.extent = box.extent.clone();
        const rectangle2 = rect(rot(rectangle.center),
            sub(rot(rectangle.axis[0]), t), sub(rot(rectangle.axis[1]), t),
            rectangle.extent.values[0], rectangle.extent.values[1]);
        const result2 = new DistRectangle3OrientedBox3()
            .compute(rectangle2, box2);
        expect(result2.distance).toBeCloseTo(result.distance, 10);
    });

    it('agrees with brute-force sampling and reports consistent points', () => {
        let seed = 24681357;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistRectangle3OrientedBox3();
        for (let trial = 0; trial < 40; ++trial) {
            const w = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(w) < 1e-3) {
                continue;
            }
            const [a0, a1] = frame(w);
            const rectangle = rect(
                v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3), a0, a1,
                0.2 + rnd(), 0.2 + rnd());

            const bw = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(bw) < 1e-3) {
                continue;
            }
            const [b0, b1, b2] = frame(bw);
            const box = new OrientedBox(3);
            box.center = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            box.axis = [b0, b1, b2];
            box.extent = v(0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd());

            const result = query.compute(rectangle, box);
            const brute = bruteForce(rectangle, box, 60);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.05);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointBox(result.closest[1], box)).toBeCloseTo(0, 8);

            // The rectangle point must lie on the rectangle.
            const delta = sub(result.closest[0], rectangle.center);
            expect(Math.abs(dot(delta, rectangle.axis[0])))
                .toBeLessThanOrEqual(rectangle.extent.values[0] + 1e-8);
            expect(Math.abs(dot(delta, rectangle.axis[1])))
                .toBeLessThanOrEqual(rectangle.extent.values[1] + 1e-8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistRectangle3OrientedBox3.h.
// ---------------------------------------------------------------------------

describe('DistRectangle3OrientedBox3 verification', () => {
    const query = new DistRectangle3OrientedBox3();
    const cQuery = new DistRectangle3CanonicalBox3();

    const boxArb = fc.tuple(wellScaledVector(3, -4, 4), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 }))
        .map(([c, R, e]) => OrientedBox.fromCenterAxisExtent(c, R,
            Vector.fromArray([e[0], e[1], e[2]])));

    const rectArb = fc.tuple(wellScaledVector(3, -6, 6), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 2, maxLength: 2 }))
        .map(([c, R, e]) => Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
            Vector.fromArray([e[0], e[1]])));

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rectArb, boxArb), ([r, b]) => {
            const res = query.compute(r, b);
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

            // closest[0] must be the world-space rectangle point named by the
            // reported W-coordinates. Leaving closest[0] in box coordinates,
            // or transforming it twice, breaks this identity for every box
            // whose frame is not the identity.
            let rebuilt = r.center.clone();
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(res.cartesian[i]))
                    .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
                rebuilt = add(rebuilt, mul(res.cartesian[i], r.axis[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);

            const delta = sub(res.closest[1], b.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(b.axis[i], delta)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-8);
            }
        });
    });

    // The upstream body is "rotate/translate into the box frame, run the
    // canonical query, map both closest points back with the same frame".
    it('is the canonical-box query composed with the box frame', () => {
        check(fc.tuple(rectArb, boxArb), ([r, b]) => {
            const toBox = (p: Vector) => {
                const d = sub(p, b.center);
                return v(dot(b.axis[0], d), dot(b.axis[1], d),
                    dot(b.axis[2], d));
            };
            const dirToBox = (p: Vector) => v(dot(b.axis[0], p),
                dot(b.axis[1], p), dot(b.axis[2], p));
            const fromBox = (p: Vector) => add(b.center,
                add(add(mul(p.values[0], b.axis[0]),
                    mul(p.values[1], b.axis[1])),
                    mul(p.values[2], b.axis[2])));
            const xr = Rectangle.fromCenterAxisExtent(toBox(r.center),
                [dirToBox(r.axis[0]), dirToBox(r.axis[1])], r.extent);
            const rc = cQuery.compute(xr,
                CanonicalBox.fromExtent(b.extent));
            const res = query.compute(r, b);
            expectClose(res.distance, rc.distance, 1e-12, 1e-12);
            expectVectorClose(res.closest[0], fromBox(rc.closest[0]),
                1e-9, 1e-9);
            expectVectorClose(res.closest[1], fromBox(rc.closest[1]),
                1e-9, 1e-9);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(rectArb, boxArb, rotationFrame(3),
            wellScaledVector(3, -4, 4)), ([r, b, R, t]) => {
                const rot = (q: Vector) => add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2]));
                const xf = (q: Vector) => add(rot(q), t);
                const r0 = query.compute(r, b);
                const r1 = query.compute(
                    Rectangle.fromCenterAxisExtent(xf(r.center),
                        [rot(r.axis[0]), rot(r.axis[1])], r.extent),
                    OrientedBox.fromCenterAxisExtent(xf(b.center),
                        [rot(b.axis[0]), rot(b.axis[1]), rot(b.axis[2])],
                        b.extent));
                // Only the distance is compared: when several pairs are
                // equidistant (a rectangle edge parallel to a box face, for
                // example) the two runs may pick different representatives.
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectClose(length(sub(r1.closest[0], r1.closest[1])),
                    r0.distance, 1e-7, 1e-7);
            });
    });

    it('never exceeds a grid sampling of the rectangle', () => {
        const rng = seededRandom(0x5a17ed09);
        for (let k = 0; k < 25; ++k) {
            const w = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(w);
            const f = frame(w);
            const b = OrientedBox.fromCenterAxisExtent(v(0.5, -1, 0.25),
                [f[0], f[1], f[2]], v(1, 2, 0.5));
            const w2 = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(w2);
            const g = frame(w2);
            const r = rect(v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4),
                g[0], g[1], 0.5 + rng(), 0.5 + rng());
            expect(query.compute(r, b).distance)
                .toBeLessThanOrEqual(bruteForce(r, b, 30) + 1e-9);
        }
    }, 30000);
});
