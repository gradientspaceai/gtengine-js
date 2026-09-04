import { describe, expect, it } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistPoint3Circle3 } from '../src/DistPoint3Circle3.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = v(...normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(v(...center), n, radius);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistPoint3Circle3', () => {
    const query = new DistPoint3Circle3();

    it('measures a point in the plane of the circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 2);
        const result = query.compute(v(5, 0, 0), c);
        expect(result.equidistant).toBe(false);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(2, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 12);
    });

    it('measures a point inside the circle in its plane', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 5);
        const result = query.compute(v(1, 0, 0), c);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(5, 12);
    });

    it('measures a point off the plane (3-4-5 configuration)', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 1);
        // Radial distance 4 from the circle, height 3.
        const result = query.compute(v(5, 0, 3), c);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports zero distance for a point on the circle', () => {
        const c = circle([1, 2, 3], [0, 1, 0], 4);
        const result = query.compute(v(5, 2, 3), c);
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('reports equidistance for a point on the normal line', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 3);
        const result = query.compute(v(0, 0, 4), c);
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBeCloseTo(5, 12);
        // The closest point is on the circle.
        const d = sub(result.closest[1], c.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(3, 10);
        expect(Math.abs(dot(d, c.normal))).toBeLessThan(1e-12);
    });

    it('reports equidistance at the circle center', () => {
        const c = circle([2, -1, 4], [1, 1, 1], 2);
        const result = query.compute(v(2, -1, 4), c);
        expect(result.equidistant).toBe(true);
        expect(result.distance).toBeCloseTo(2, 12);
    });

    it('handles a zero-radius circle', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 0);
        const result = query.compute(v(3, 4, 12), c);
        expect(result.distance).toBeCloseTo(13, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(0, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(0, 10);
    });

    it('is accurate for a point nearly on the normal line', () => {
        const c = circle([0, 0, 0], [0, 0, 1], 1);
        const result = query.compute(v(1e-14, 0, 1), c);
        expect(result.equidistant).toBe(false);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 10);
        const d = sub(result.closest[1], c.center);
        expect(Math.sqrt(dot(d, d))).toBeCloseTo(1, 10);
    });

    it('agrees with a dense sampling of the circle', () => {
        const rnd = makeRandom(31415);
        const center = v(0.5, -1, 2);
        const normal = v(1, 2, -0.5);
        normalize(normal);
        const radius = 1.5;
        const c = Circle3.fromCenterNormalRadius(center, normal, radius);

        // An orthonormal basis of the plane of the circle.
        const U = v(-normal.values[1], normal.values[0], 0);
        normalize(U);
        const W = cross(normal, U);

        for (let trial = 0; trial < 50; ++trial) {
            const p = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const result = query.compute(p, c);

            // The reported closest point is on the circle.
            const d = sub(result.closest[1], center);
            expect(Math.sqrt(dot(d, d))).toBeCloseTo(radius, 8);
            expect(Math.abs(dot(d, normal))).toBeLessThan(1e-9);

            // The reported closest point realizes the reported distance.
            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 8);

            // No sampled circle point is closer.
            const n = 5000;
            let best = Number.MAX_VALUE;
            for (let i = 0; i < n; ++i) {
                const t = 2 * Math.PI * i / n;
                const q = add(center, add(mul(radius * Math.cos(t), U),
                    mul(radius * Math.sin(t), W)));
                const f = sub(p, q);
                best = Math.min(best, dot(f, f));
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistPoint3Circle3.h.
// ---------------------------------------------------------------------------

describe('DistPoint3Circle3 verification', () => {
    const query = new DistPoint3Circle3();

    const circleArb = fc.tuple(wellScaledVector(3, -5, 5), unitVector(3),
        finite(0.25, 4)).map(([c, n, r]) =>
            Circle3.fromCenterNormalRadius(c, n, r));

    function planeBasis(circle: Circle3): [Vector, Vector] {
        const n = circle.normal;
        const seed = Math.abs(n.values[0]) < 0.9 ? v(1, 0, 0) : v(0, 1, 0);
        const u = sub(seed, mul(dot(seed, n), n));
        normalize(u);
        return [u, cross(n, u)];
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), circleArb),
            ([p, circle]) => {
                const r = query.compute(p, circle);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                // closest[0] is a copy of the input point.
                expectVectorClose(r.closest[0], p, 0, 0);
                expect(r.closest[0]).not.toBe(p);
                // closest[1] is on the circle.
                const delta = sub(r.closest[1], circle.center);
                expectClose(length(delta), circle.radius, 1e-9, 1e-9);
                expectClose(dot(circle.normal, delta), 0, 1e-9, 1e-9);
                expectClose(length(sub(r.closest[0], r.closest[1])),
                    r.distance, 1e-9, 1e-9);
            });
    });

    it('matches the closed form height^2 + (radial - r)^2', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), circleArb),
            ([p, circle]) => {
                const r = query.compute(p, circle);
                const delta = sub(p, circle.center);
                const height = dot(circle.normal, delta);
                const radial = length(sub(delta, mul(height, circle.normal)));
                const expected = r.equidistant
                    ? Math.sqrt(dot(delta, delta)
                        + circle.radius * circle.radius)
                    : Math.sqrt(height * height
                        + (radial - circle.radius) ** 2);
                expectClose(r.distance, expected, 1e-9, 1e-9);
            });
    });

    it('is the minimum over sampled circle points', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), circleArb),
            ([p, circle]) => {
                const [u, w] = planeBasis(circle);
                const n = 4096;
                let best = Number.POSITIVE_INFINITY;
                for (let i = 0; i < n; ++i) {
                    const a = (2 * Math.PI * i) / n;
                    const q = add(circle.center,
                        add(mul(circle.radius * Math.cos(a), u),
                            mul(circle.radius * Math.sin(a), w)));
                    best = Math.min(best, length(sub(q, p)));
                }
                const r = query.compute(p, circle);
                expect(r.distance).toBeLessThanOrEqual(best + 1e-9);
                // 1-Lipschitz in arc length, so the sampled minimum is at
                // most half an arc-length step above the true minimum.
                const step = (2 * Math.PI * circle.radius) / n;
                expect(r.distance).toBeGreaterThanOrEqual(best - step);
            }, 60);
    }, 30000);

    it('gives sqrt(h^2 + r^2) on the normal line through the center', () => {
        // The 'equidistant' flag is triggered by an exactly zero in-plane
        // projection. For a general unit normal the products that cancel in
        // Dot(U, P-C) are rounded independently, so the projection is a
        // nonzero 1e-17-scale vector and upstream takes the non-equidistant
        // branch. Both branches produce the same distance and a point of the
        // circle, which is what this property pins down; the flag itself is
        // exercised by the axis-aligned fixtures below.
        check(fc.tuple(circleArb, wellScaled(-6, 6)), ([circle, h]) => {
            const p = add(circle.center, mul(h, circle.normal));
            const r = query.compute(p, circle);
            expectClose(r.distance,
                Math.sqrt(h * h + circle.radius * circle.radius), 1e-9, 1e-9);
            const delta = sub(r.closest[1], circle.center);
            expectClose(length(delta), circle.radius, 1e-9, 1e-9);
            expectClose(dot(circle.normal, delta), 0, 1e-9, 1e-9);
        });
    });

    it('sets equidistant for exactly representable axis normals', () => {
        for (const n of [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1),
            v(-1, 0, 0), v(0, -1, 0), v(0, 0, -1)]) {
            for (const h of [-3, -1, 0, 0.5, 4]) {
                const circle = Circle3.fromCenterNormalRadius(v(1, -2, 3),
                    n, 0.75);
                const p = add(circle.center, mul(h, n));
                const r = query.compute(p, circle);
                expect(r.equidistant).toBe(true);
                expect(r.distance).toBeCloseTo(
                    Math.sqrt(h * h + 0.75 * 0.75), 12);
            }
        }
    });

    it('is not equidistant off the normal line', () => {
        check(fc.tuple(circleArb, wellScaled(-6, 6), finite(0.05, 6),
            finite(-Math.PI, Math.PI)), ([circle, h, radial, angle]) => {
            const [u, w] = planeBasis(circle);
            const p = add(circle.center, add(mul(h, circle.normal),
                add(mul(radial * Math.cos(angle), u),
                    mul(radial * Math.sin(angle), w))));
            const r = query.compute(p, circle);
            expect(r.equidistant).toBe(false);
        });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), circleArb,
            rotationFrame(3), wellScaledVector(3, -6, 6)),
            ([p, circle, frame, shift]) => {
                const rot = (q: Vector): Vector =>
                    add(add(mul(q.values[0], frame[0]),
                        mul(q.values[1], frame[1])),
                        mul(q.values[2], frame[2]));
                const moved = Circle3.fromCenterNormalRadius(
                    add(shift, rot(circle.center)), rot(circle.normal),
                    circle.radius);
                const r0 = query.compute(p, circle);
                const r1 = query.compute(add(shift, rot(p)), moved);
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            });
    });

    it('does not depend on the orientation of the circle normal', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), circleArb),
            ([p, circle]) => {
                const flipped = Circle3.fromCenterNormalRadius(circle.center,
                    mul(-1, circle.normal), circle.radius);
                expectClose(query.compute(p, circle).distance,
                    query.compute(p, flipped).distance, 1e-9, 1e-9);
            });
    });
});
