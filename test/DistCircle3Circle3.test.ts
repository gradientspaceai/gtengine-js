import { describe, it, expect } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistCircle3Circle3 } from '../src/DistCircle3Circle3.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame, scaled,
    seededRandom, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function circle(center: Vector, normal: Vector, radius: number): Circle3 {
    const n = normal.clone();
    normalize(n);
    return Circle3.fromCenterNormalRadius(center, n, radius);
}

// An orthonormal pair spanning the plane of the circle.
function planeBasis(n: Vector): [Vector, Vector] {
    let u = Math.abs(n.values[0]) > Math.abs(n.values[1])
        ? v(-n.values[2], 0, n.values[0])
        : v(0, n.values[2], -n.values[1]);
    normalize(u);
    return [u, cross(n, u)];
}

function circlePoint(c: Circle3, u: Vector, w: Vector, angle: number): Vector {
    return add(c.center, add(mul(c.radius * Math.cos(angle), u),
        mul(c.radius * Math.sin(angle), w)));
}

// The exact distance from a point to a circle: split the offset into the
// components along and perpendicular to the circle normal.
function distPointCircle(p: Vector, c: Circle3): number {
    const delta = sub(p, c.center);
    const along = dot(c.normal, delta);
    const inPlane = length(sub(delta, mul(along, c.normal)));
    const radial = inPlane - c.radius;
    return Math.sqrt(along * along + radial * radial);
}

// Brute force: sample circle1 and use the exact point-circle distance to
// circle0. Sampling only one of the two circles is enough because the inner
// minimization is solved in closed form.
function bruteForce(c0: Circle3, c1: Circle3, n: number): number {
    const [u, w] = planeBasis(c1.normal);
    let best = Infinity;
    for (let i = 0; i < n; ++i) {
        const d = distPointCircle(
            circlePoint(c1, u, w, 2 * Math.PI * i / n), c0);
        if (d < best) {
            best = d;
        }
    }
    return best;
}

// A point lies on the circle when it is in the plane and at the radius.
function expectOnCircle(p: Vector, c: Circle3): void {
    const delta = sub(p, c.center);
    expect(Math.abs(dot(c.normal, delta))).toBeLessThan(1e-8);
    expect(length(delta)).toBeCloseTo(c.radius, 8);
}

describe('DistCircle3Circle3', () => {
    describe('coplanar circles (the parallel-planes branch)', () => {
        const inPlane = (cx: number, radius: number) =>
            circle(v(cx, 0, 0), v(0, 0, 1), radius);

        it('measures separated circles', () => {
            const result = new DistCircle3Circle3()
                .compute(inPlane(0, 1), inPlane(5, 2));
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.sqrDistance).toBeCloseTo(4, 12);
            expect(result.numClosestPairs).toBe(1);
            expect(result.equidistant).toBe(false);
            expect(result.circle0Closest[0].values[0]).toBeCloseTo(1, 12);
            expect(result.circle1Closest[0].values[0]).toBeCloseTo(3, 12);
        });

        it('reports zero for externally tangent circles', () => {
            const result = new DistCircle3Circle3()
                .compute(inPlane(0, 1), inPlane(3, 2));
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.numClosestPairs).toBe(1);
            expect(length(sub(result.circle0Closest[0],
                result.circle1Closest[0]))).toBeCloseTo(0, 10);
        });

        it('reports two intersection points for overlapping circles', () => {
            const result = new DistCircle3Circle3()
                .compute(inPlane(0, 1), inPlane(1, 1));
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.numClosestPairs).toBe(2);
            const h = Math.sqrt(3) / 2;
            for (let i = 0; i < 2; ++i) {
                expect(result.circle0Closest[i].values[0]).toBeCloseTo(0.5, 10);
                expect(Math.abs(result.circle0Closest[i].values[1]))
                    .toBeCloseTo(h, 10);
                expect(length(sub(result.circle0Closest[i],
                    result.circle1Closest[i]))).toBeCloseTo(0, 10);
            }
            expect(result.circle0Closest[0].values[1])
                .toBeCloseTo(-result.circle0Closest[1].values[1], 10);
        });

        it('measures one circle nested inside the other', () => {
            // Radii 1 and 3, centers coincident: every pair of points is
            // 2 apart, so the query flags the result as equidistant.
            const result = new DistCircle3Circle3()
                .compute(inPlane(0, 1), inPlane(0, 3));
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.equidistant).toBe(true);
            expect(result.numClosestPairs).toBe(1);
            expectOnCircle(result.circle0Closest[0], inPlane(0, 1));
            expectOnCircle(result.circle1Closest[0], inPlane(0, 3));
        });

        it('measures an off-center nested circle', () => {
            // Radii 1 and 4 with centers 0.5 apart: the closest points are on
            // the side where the gap is smallest, 4 - 0.5 - 1 = 2.5.
            const result = new DistCircle3Circle3()
                .compute(inPlane(0, 1), inPlane(0.5, 4));
            expect(result.distance).toBeCloseTo(2.5, 10);
            expect(result.equidistant).toBe(false);
            expect(result.circle0Closest[0].values[0]).toBeCloseTo(-1, 10);
            expect(result.circle1Closest[0].values[0]).toBeCloseTo(-3.5, 10);
        });

        it('reports zero for identical circles', () => {
            const c = circle(v(1, 2, 3), v(0, 1, 0), 1.5);
            const result = new DistCircle3Circle3().compute(c, c.clone());
            expect(result.distance).toBeCloseTo(0, 12);
            expect(result.equidistant).toBe(true);
        });
    });

    describe('circles in parallel but distinct planes', () => {
        it('measures coaxial circles of equal radius', () => {
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 1);
            const c1 = circle(v(0, 0, 4), v(0, 0, 1), 1);
            const result = new DistCircle3Circle3().compute(c0, c1);
            expect(result.distance).toBeCloseTo(4, 12);
            expect(result.equidistant).toBe(true);
            expectOnCircle(result.circle0Closest[0], c0);
            expectOnCircle(result.circle1Closest[0], c1);
        });

        it('measures coaxial circles of different radii', () => {
            // Radial gap 3 - 1 = 2 and axial gap 2 give sqrt(8).
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 3);
            const c1 = circle(v(0, 0, 2), v(0, 0, 1), 1);
            const result = new DistCircle3Circle3().compute(c0, c1);
            expect(result.distance).toBeCloseTo(Math.sqrt(8), 10);
            expect(result.equidistant).toBe(true);
        });

        it('ignores the sign of the circle normals', () => {
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 1);
            const c1 = circle(v(0, 0, 4), v(0, 0, -1), 1);
            const result = new DistCircle3Circle3().compute(c0, c1);
            expect(result.distance).toBeCloseTo(4, 12);
        });
    });

    describe('circles in non-parallel planes', () => {
        it('measures two perpendicular circles', () => {
            // Circle0 is the unit circle in z = 0; circle1 is the unit circle
            // in x = 3. The configuration is mirror symmetric about y = 0, so
            // there are two closest pairs. The minimum of
            // sqrt((cos a - 3)^2 + (sin a - 1)^2) over a is attained where
            // the point of circle1 is (3,+/-1,0).
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 1);
            const c1 = circle(v(3, 0, 0), v(1, 0, 0), 1);
            const result = new DistCircle3Circle3().compute(c0, c1);
            const brute = bruteForce(c0, c1, 200000);
            expect(result.distance).toBeCloseTo(brute, 6);
            expect(result.distance).toBeCloseTo(2.1622776601683795, 10);
            expect(result.numClosestPairs).toBe(2);
            expectOnCircle(result.circle0Closest[0], c0);
            expectOnCircle(result.circle1Closest[0], c1);
            expect(length(sub(result.circle0Closest[0],
                result.circle1Closest[0]))).toBeCloseTo(result.distance, 8);
        });

        it('measures two linked circles', () => {
            // Circle0 has radius 2 in z = 0 and circle1 has radius 2 in
            // x = 2 centered at (2,0,0). The closest points are at
            // (2,+/-2,0) on circle1 and the nearest point of circle0 is the
            // radial projection at distance 2*sqrt(2) - 2.
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 2);
            const c1 = circle(v(2, 0, 0), v(1, 0, 0), 2);
            const result = new DistCircle3Circle3().compute(c0, c1);
            expect(result.distance).toBeCloseTo(2 * Math.SQRT2 - 2, 10);
            expect(result.numClosestPairs).toBe(2);
        });

        it('reports zero for circles that intersect', () => {
            // Both circles pass through (0,1,0) and (0,-1,0).
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 1);
            const c1 = circle(v(0, 0, 0), v(1, 0, 0), 1);
            const result = new DistCircle3Circle3().compute(c0, c1);
            expect(result.distance).toBeCloseTo(0, 8);
            expect(length(sub(result.circle0Closest[0],
                result.circle1Closest[0]))).toBeCloseTo(0, 8);
        });

        it('is symmetric in its arguments', () => {
            const c0 = circle(v(0.3, 1.1, -0.4), v(0.3, -0.7, 0.5), 1.3);
            const c1 = circle(v(-1.2, 0.5, 2.1), v(-0.2, 0.4, 0.9), 0.7);
            const query = new DistCircle3Circle3();
            const forward = query.compute(c0, c1);
            const backward = query.compute(c1, c0);
            expect(backward.distance).toBeCloseTo(forward.distance, 10);
            expect(length(sub(backward.circle1Closest[0],
                forward.circle0Closest[0]))).toBeCloseTo(0, 6);
            expect(length(sub(backward.circle0Closest[0],
                forward.circle1Closest[0]))).toBeCloseTo(0, 6);
        });

        it('is invariant under translation', () => {
            const t = v(-7, 3.5, 11);
            const c0 = circle(v(0.3, 1.1, -0.4), v(0.3, -0.7, 0.5), 1.3);
            const c1 = circle(v(-1.2, 0.5, 2.1), v(-0.2, 0.4, 0.9), 0.7);
            const query = new DistCircle3Circle3();
            const base = query.compute(c0, c1);
            const moved = query.compute(
                circle(add(c0.center, t), c0.normal, c0.radius),
                circle(add(c1.center, t), c1.normal, c1.radius));
            expect(moved.distance).toBeCloseTo(base.distance, 9);
            expect(length(sub(sub(moved.circle0Closest[0], t),
                base.circle0Closest[0]))).toBeCloseTo(0, 6);
        });
    });

    it('agrees with brute-force sampling on random configurations', () => {
        let seed = 12345;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistCircle3Circle3();
        let tested = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const mk = () => {
                const n = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
                if (length(n) < 1e-2) {
                    return null;
                }
                return circle(v(rnd() * 4 - 2, rnd() * 4 - 2, rnd() * 4 - 2),
                    n, 0.2 + rnd() * 2);
            };
            const c0 = mk(), c1 = mk();
            if (!c0 || !c1) {
                continue;
            }
            ++tested;

            const result = query.compute(c0, c1);
            const brute = Math.min(bruteForce(c0, c1, 4000),
                bruteForce(c1, c0, 4000));
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-6);
            expect(result.distance).toBeGreaterThan(brute - 1e-3);
            expect(result.sqrDistance)
                .toBeCloseTo(result.distance * result.distance, 9);
            expect(result.numClosestPairs).toBeGreaterThanOrEqual(1);
            expect(result.numClosestPairs).toBeLessThanOrEqual(2);

            for (let i = 0; i < result.numClosestPairs; ++i) {
                expect(length(sub(result.circle0Closest[i],
                    result.circle1Closest[i])))
                    .toBeCloseTo(result.distance, 6);
                expectOnCircle(result.circle0Closest[i], c0);
                expectOnCircle(result.circle1Closest[i], c1);
            }
        }
        expect(tested).toBeGreaterThan(100);
    }, 30000);
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistCircle3Circle3.h. The port deviates from upstream
// in three places, all recorded in the source comments: it also feeds the
// roots of p6 into the candidate list (upstream misses the true minimum for
// mirror-symmetric configurations, where round-off splits a double root of
// phi), it counts candidates rather than raw roots when deciding whether to
// report a second closest pair, and it throws on an empty candidate list
// instead of reading a default-constructed ClosestInfo. The properties below
// pin all three down and characterize the inputs that produce an empty
// candidate list.
// ---------------------------------------------------------------------------

describe('DistCircle3Circle3 verification', () => {
    const query = new DistCircle3Circle3();

    const circlePairArb = fc.tuple(
        fc.tuple(wellScaledVector(3, -4, 4), unitVector(3), finite(0.2, 3)),
        fc.tuple(wellScaledVector(3, -4, 4), unitVector(3), finite(0.2, 3)))
        .map(([a, b]) => [Circle3.fromCenterNormalRadius(a[0], a[1], a[2]),
            Circle3.fromCenterNormalRadius(b[0], b[1], b[2])] as
            [Circle3, Circle3]);

    function expectValidResult(c0: Circle3, c1: Circle3, r: {
        distance: number, sqrDistance: number, numClosestPairs: number,
        circle0Closest: [Vector, Vector], circle1Closest: [Vector, Vector]
    }): void {
        expectClose(r.sqrDistance, r.distance * r.distance, 1e-9, 1e-9);
        expect(r.numClosestPairs === 1 || r.numClosestPairs === 2).toBe(true);
        for (let i = 0; i < r.numClosestPairs; ++i) {
            const p0 = r.circle0Closest[i];
            const p1 = r.circle1Closest[i];
            // Each reported point lies in the plane of its circle and at the
            // circle's radius from the center.
            expectClose(dot(c0.normal, sub(p0, c0.center)), 0, 1e-7, 1e-7);
            expectClose(length(sub(p0, c0.center)), c0.radius, 1e-7, 1e-7);
            expectClose(dot(c1.normal, sub(p1, c1.center)), 0, 1e-7, 1e-7);
            expectClose(length(sub(p1, c1.center)), c1.radius, 1e-7, 1e-7);
            // Every reported pair realizes the reported distance. Upstream
            // fabricates a second pair at the origin whenever the unfiltered
            // root count exceeds one but only a single candidate survives;
            // that pair would fail here.
            expectClose(length(sub(p0, p1)), r.distance, 1e-6, 1e-6);
        }
    }

    it('reports valid closest pairs on the two circles', () => {
        check(circlePairArb, ([c0, c1]) => {
            expectValidResult(c0, c1, query.compute(c0, c1));
        }, 100);
    }, 30000);

    it('is symmetric under argument swap', () => {
        check(circlePairArb, ([c0, c1]) => {
            const a = query.compute(c0, c1);
            const b = query.compute(c1, c0);
            // Only the distance is symmetric: when several pairs are
            // equidistant the two argument orders may report a different
            // number of representatives.
            expectClose(a.distance, b.distance, 1e-7, 1e-7);
            expectValidResult(c1, c0, b);
        }, 100);
    }, 30000);

    // PrepareCircles applies a translation, a rotation and a uniform scaling
    // and the query inverts all three at the end, so the result must be
    // equivariant under them.
    it('is equivariant under a rigid motion and uniform scaling', () => {
        check(fc.tuple(circlePairArb, rotationFrame(3),
            wellScaledVector(3, -4, 4), finite(0.25, 4)),
            ([[c0, c1], R, t, s]) => {
                const rot = (p: Vector) => add(add(mul(p.values[0], R[0]),
                    mul(p.values[1], R[1])), mul(p.values[2], R[2]));
                const xf = (c: Circle3) => Circle3.fromCenterNormalRadius(
                    add(mul(s, rot(c.center)), t), rot(c.normal), s * c.radius);
                const a = query.compute(c0, c1);
                const b = query.compute(xf(c0), xf(c1));
                expectClose(s * a.distance, b.distance, 1e-7, 1e-7);
            }, 100);
    }, 30000);

    it('is unaffected by the sign of either circle normal', () => {
        check(circlePairArb, ([c0, c1]) => {
            const flip = (c: Circle3) => Circle3.fromCenterNormalRadius(
                c.center, mul(-1, c.normal), c.radius);
            const a = query.compute(c0, c1);
            expectClose(query.compute(flip(c0), c1).distance, a.distance,
                1e-7, 1e-7);
            expectClose(query.compute(c0, flip(c1)).distance, a.distance,
                1e-7, 1e-7);
            expectClose(query.compute(flip(c0), flip(c1)).distance, a.distance,
                1e-7, 1e-7);
        }, 100);
    }, 30000);

    // Two-sided cross-check. Sampling one circle and using the exact
    // point-to-circle distance to the other gives an upper bound for the true
    // minimum; because that distance is 1-Lipschitz along the sampled circle,
    // the bound exceeds the true minimum by at most half the arc spacing. The
    // reported distance is realized by a genuine pair, so it is never below
    // the true minimum: the two bounds bracket it.
    it('brackets a dense sampling of both circles', () => {
        const rng = seededRandom(0xc1c1e001);
        const N = 4000;
        for (let k = 0; k < 60; ++k) {
            const mk = () => {
                const n = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
                normalize(n);
                return Circle3.fromCenterNormalRadius(
                    v(6 * rng() - 3, 6 * rng() - 3, 6 * rng() - 3), n,
                    0.2 + 2.8 * rng());
            };
            const c0 = mk(), c1 = mk();
            const r = query.compute(c0, c1);
            const sampled = Math.min(bruteForce(c0, c1, N),
                bruteForce(c1, c0, N));
            const bound = Math.PI * Math.max(c0.radius, c1.radius) / N + 1e-9;
            expect(r.distance).toBeLessThanOrEqual(sampled + 1e-7);
            expect(sampled - r.distance).toBeLessThanOrEqual(bound);
        }
    }, 30000);

    // The upstream defect fixed by the port (see the source comments): when
    // the configuration is mirror symmetric, phi = p6^2 - (1-c^2)*p7^2 has a
    // double root that round-off splits into two roots too close together for
    // the bisection to separate, so upstream finds neither and reports a
    // critical point that is not the closest pair. Every configuration below
    // is mirror symmetric about the yz-plane.
    it('finds the minimum for mirror-symmetric configurations', () => {
        const rng = seededRandom(0xc1c1e002);
        const N = 4000;
        for (let k = 0; k < 60; ++k) {
            const tilt = 0.05 + (Math.PI - 0.1) * rng();
            const c1 = Circle3.fromCenterNormalRadius(v(0, 0, 0), v(0, 0, 1),
                1);
            const c0 = Circle3.fromCenterNormalRadius(
                v(0, 4 * rng() - 2, 4 * rng() - 2),
                v(0, -Math.sin(tilt), Math.cos(tilt)), 0.2 + 1.5 * rng());
            const r = query.compute(c0, c1);
            expectValidResult(c0, c1, r);
            const sampled = Math.min(bruteForce(c0, c1, N),
                bruteForce(c1, c0, N));
            const bound = Math.PI * Math.max(c0.radius, c1.radius) / N + 1e-9;
            expect(sampled - r.distance).toBeLessThanOrEqual(bound);
        }
    }, 30000);

    // Coplanar circles take DoQueryParallelPlanes, whose four branches have
    // closed-form answers. 'scaled' keeps the center separation on a uniform
    // grid; a subnormal separation makes PrepareCircles divide by an
    // underflowed length (see the last test in this block).
    it('matches the closed forms for coplanar circles', () => {
        check(fc.tuple(unitVector(3), wellScaledVector(3, -4, 4),
            finite(0.2, 3), finite(0.2, 3), scaled(0, 8)),
            ([n, center, r0, r1, d]) => {
                // Offset the second center by d along a direction lying in
                // the common plane.
                let u = Math.abs(n.values[0]) > Math.abs(n.values[1])
                    ? v(-n.values[2], 0, n.values[0])
                    : v(0, n.values[2], -n.values[1]);
                normalize(u);
                const c0 = Circle3.fromCenterNormalRadius(center, n, r0);
                const c1 = Circle3.fromCenterNormalRadius(
                    add(center, mul(d, u)), n, r1);
                const r = query.compute(c0, c1);
                let expected: number;
                if (d >= r0 + r1) {
                    expected = d - r0 - r1;          // separated
                } else if (d + r1 <= r0) {
                    expected = r0 - d - r1;          // circle1 inside circle0
                } else if (d - r1 <= -r0) {
                    expected = r1 - d - r0;          // circle0 inside circle1
                } else {
                    expected = 0;                    // the circles cross
                }
                expectClose(r.distance, expected, 1e-7, 1e-7);
                expectValidResult(c0, c1, r);
            }, 100);
    });

    it('reports the two crossing points for overlapping coplanar circles',
        () => {
            const c0 = circle(v(0, 0, 0), v(0, 0, 1), 2);
            const c1 = circle(v(3, 0, 0), v(0, 0, 1), 2);
            const r = query.compute(c0, c1);
            expect(r.distance).toBeCloseTo(0, 9);
            expect(r.numClosestPairs).toBe(2);
            expect(r.equidistant).toBe(false);
            // The crossing points are (1.5, +-sqrt(7)/2, 0).
            for (let i = 0; i < 2; ++i) {
                expect(r.circle0Closest[i].values[0]).toBeCloseTo(1.5, 9);
                expect(Math.abs(r.circle0Closest[i].values[1]))
                    .toBeCloseTo(Math.sqrt(7) / 2, 9);
                expectVectorClose(r.circle0Closest[i], r.circle1Closest[i],
                    1e-9, 1e-9);
            }
        });

    it('flags concentric coplanar circles as equidistant', () => {
        const c0 = circle(v(1, 2, 3), v(0, 0, 1), 1);
        const c1 = circle(v(1, 2, 3), v(0, 0, 1), 4);
        const r = query.compute(c0, c1);
        expect(r.distance).toBeCloseTo(3, 9);
        expect(r.numClosestPairs).toBe(1);
        expect(r.equidistant).toBe(true);
        expectValidResult(c0, c1, r);
    });

    it('adds the plane offset for coaxial circles in parallel planes', () => {
        // Coaxial circles of radii 1 and 4 whose planes are 5 apart: the
        // in-plane gap is 3 and the out-of-plane gap is 5.
        const c0 = circle(v(0, 0, 0), v(0, 0, 1), 1);
        const c1 = circle(v(0, 0, 5), v(0, 0, 1), 4);
        const r = query.compute(c0, c1);
        expect(r.distance).toBeCloseTo(Math.sqrt(9 + 25), 9);
        expect(r.equidistant).toBe(true);
        expectValidResult(c0, c1, r);
    });

    it('degenerates to the point-circle distance for a zero radius', () => {
        check(fc.tuple(wellScaledVector(3, -4, 4), unitVector(3),
            wellScaledVector(3, -4, 4), finite(0.2, 3))
            .filter(([p, , c]) => length(sub(p, c)) > 1e-2),
            ([p, n, c, radius]) => {
                const point = Circle3.fromCenterNormalRadius(p, v(0, 0, 1), 0);
                const big = Circle3.fromCenterNormalRadius(c, n, radius);
                const r = query.compute(point, big);
                expectClose(r.distance, distPointCircle(p, big), 1e-7, 1e-7);
                expectVectorClose(r.circle0Closest[0], p, 1e-7, 1e-7);
            }, 100);
    });

    // Port fix 4 (upstream defect). A zero-radius circle whose center is on
    // circle1's axis makes p6 and p7 identically zero -- every point of
    // circle1 is equidistant from it -- and upstream reports its
    // "Unexpected degree for p6" assertion instead of the answer.
    it('answers a point on the circle axis', () => {
        const inv = 1 / Math.sqrt(2);
        for (const dummyNormal of [v(1, 0, 0), v(0, 0, 1), v(0, inv, inv)]) {
            const r = query.compute(
                Circle3.fromCenterNormalRadius(v(0, 0, 5), dummyNormal, 0),
                circle(v(0, 0, 0), v(0, 0, 1), 2));
            expect(r.distance).toBeCloseTo(Math.sqrt(25 + 4), 9);
        }
    });

    // Port fix 1 (upstream defect). PrepareCircles tries to align the two
    // normals by negating any normal whose z-component is negative, which
    // does nothing when both normals are orthogonal to the z-axis. Two
    // circles in parallel planes with anti-parallel normals then reach the
    // 'circle0.normal[2] < 1' test as (0,0,-1), are misclassified as
    // non-parallel, and the general polynomial path answers 3.7385 instead of
    // 5.1851. The port aligns the transformed normal with +z.
    it('handles parallel planes whose normals are anti-parallel', () => {
        const truth = Math.sqrt(
            (3.998533561266931 - 0.2 - 0.20000000000000032) ** 2
            + 3.733128536483254 ** 2);
        const c0 = circle(v(0, -3.733128536483254, 0), v(0, -1, 0), 0.2);
        for (const n1 of [v(0, 1, 0), v(0, -1, 0)]) {
            const c1 = Circle3.fromCenterNormalRadius(
                v(3.998533561266931, 0, 0), n1, 0.20000000000000032);
            const r = query.compute(c0, c1);
            expect(r.distance).toBeCloseTo(truth, 9);
            expectValidResult(c0, c1, r);
        }
    });

    // Port fix 2 (upstream defect). When circle0's center is on circle1's
    // axis, phi and p6 are perfect squares, the sign-change bisection finds
    // no roots, and upstream reads a default-constructed ClosestInfo:
    // distance 0 with both closest points at the shared axis point. The port
    // falls back on the roots of the derivative.
    it('solves concentric and coaxial circles', () => {
        const rng = seededRandom(0xc1c1e003);
        const N = 4000;
        for (let k = 0; k < 60; ++k) {
            const c1 = Circle3.fromCenterNormalRadius(v(0, 0, 0), v(0, 0, 1),
                0.2 + 2.5 * rng());
            const n0 = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(n0);
            // k = 0 .. 29 concentric, 30 .. 59 coaxial.
            const z = k < 30 ? 0 : 6 * rng() - 3;
            const c0 = Circle3.fromCenterNormalRadius(v(0, 0, z), n0,
                0.2 + 2.5 * rng());
            const r = query.compute(c0, c1);
            expectValidResult(c0, c1, r);
            const sampled = Math.min(bruteForce(c0, c1, N),
                bruteForce(c1, c0, N));
            const bound = Math.PI * Math.max(c0.radius, c1.radius) / N + 1e-7;
            expect(r.distance).toBeLessThanOrEqual(sampled + 1e-7);
            expect(sampled - r.distance).toBeLessThanOrEqual(bound);
        }
    }, 30000);

    it('solves the perpendicular-concentric and near-coaxial families', () => {
        // Concentric circles in perpendicular planes make p4^2-(1-cs^2)*p5^2
        // a perfect square as well (p5 or p4 vanishes identically), so the
        // roots of p4 and p5 themselves are needed. Near-coaxial circles are
        // covered by the same candidates because a1 is only approximately
        // zero after the transformation in prepareCircles.
        const rng = seededRandom(0xc1c1e004);
        const N = 4000;
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        const cases: Array<[Circle3, Circle3]> = [];
        for (let k = 0; k < 12; ++k) {
            const r0 = 0.2 + 2 * rng();
            const r1 = k % 3 === 0 ? r0 : 0.2 + 2 * rng();
            cases.push([
                Circle3.fromCenterNormalRadius(v(0, 0, 0), axes[k % 3], r0),
                Circle3.fromCenterNormalRadius(v(0, 0, 0), axes[(k + 1) % 3],
                    r1)]);
        }
        for (let k = 0; k < 12; ++k) {
            const eps = [0, 1e-12, 1e-6, 1e-3][k % 4];
            const n0 = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(n0);
            cases.push([
                Circle3.fromCenterNormalRadius(v(0, 0, 6 * rng() - 3 + eps),
                    n0, 0.2 + 2.5 * rng()),
                Circle3.fromCenterNormalRadius(v(0, 0, 0), v(0, 0, 1),
                    0.2 + 2.5 * rng())]);
        }
        for (const [c0, c1] of cases) {
            const r = query.compute(c0, c1);
            expectValidResult(c0, c1, r);
            const sampled = Math.min(bruteForce(c0, c1, N),
                bruteForce(c1, c0, N));
            const bound = Math.PI * Math.max(c0.radius, c1.radius) / N + 1e-7;
            expect(r.distance).toBeLessThanOrEqual(sampled + 1e-7);
            expect(sampled - r.distance).toBeLessThanOrEqual(bound);
        }
    }, 30000);

    it('reproduces the worked concentric and coaxial examples', () => {
        // Concentric circles of radii 0.2 and 0.5 in perpendicular planes:
        // every point of the small circle is 0.3 from the large one.
        const a = query.compute(circle(v(0, 0, 0), v(0, 0, 1), 0.2),
            circle(v(0, 0, 0), v(0.1, 0.995, 0), 0.5));
        expect(a.distance).toBeCloseTo(0.3, 9);

        // A unit circle in the plane x = 0 centered at (0,0,1) and a circle
        // of radius 2 in the plane z = 0 centered at the origin.
        const b = query.compute(circle(v(0, 0, 1), v(1, 0, 0), 1),
            circle(v(0, 0, 0), v(0, 0, 1), 2));
        expect(b.distance).toBeCloseTo(1.2360679774997898, 9);
    });

    // Port fix 3 (upstream defect). Upstream sets sn = -p6(cs)/p7(cs) for a
    // root cs of phi, which is only on the unit circle when cs is exact. The
    // roots that matter are usually double roots, which the bisection can
    // place only to about sqrt(eps), so the quotient lands off the unit
    // circle and the reported "closest point" is not on the circle at all --
    // here 0.43 inside it, with a distance 0.22 below the true minimum.
    it('never reports a closest point off its circle', () => {
        const c0 = circle(v(0, 0, -0.3195239086635411),
            v(-0.7544658663907117, 0.36320170601171214,
                -0.5466861779865069), 1.5768367294687777);
        const c1 = circle(v(0, 0, 0), v(0, 0, 1), 0.45176617957185955);
        const r = query.compute(c0, c1);
        expectValidResult(c0, c1, r);
        // Upstream (and the port before the fix) reported 0.79957.
        expect(r.distance).toBeCloseTo(1.0234939702204247, 9);
        expect(length(sub(r.circle0Closest[0], c0.center)))
            .toBeCloseTo(c0.radius, 12);
    });

    // Port fix 5 (upstream defect). DoQueryParallelPlanes splits the vector
    // between the two centers into the components along and perpendicular to
    // the common plane normal. For coaxial circles the perpendicular
    // component is mathematically zero, but the normal that reaches the
    // function is a unit vector only up to the rounding of the rotation in
    // PrepareCircles, so the computed component is a residue of size
    // |D|*O(eps) directed along the normal rather than zero. Upstream takes
    // its length as evidence that the circles are not concentric and
    // normalizes it into the radial direction of both closest points, which
    // then land on the common axis instead of on the circles; the reported
    // distance stays correct. Both routes into the parallel-planes code are
    // covered here: the exact 'circle0.normal[2] < 1' test being false, and
    // the p6 == p7 == 0 fallback of port fix 4.
    it('reports coaxial closest points on the circles, not on the axis', () => {
        // The p6 == p7 == 0 fallback. The transformed normal is one ulp
        // below 1, so the polynomial path runs and the fallback fires. Before
        // the fix the closest points were (0,0,-0.2) and (0,0,-0.201), both
        // on the common axis and 0.2 away from their own circles.
        const c0 = Circle3.fromCenterNormalRadius(v(0, 0, 0), v(0, 0, -1), 0.2);
        const c1 = Circle3.fromCenterNormalRadius(v(0, 0, -0.001),
            v(0, 0, 0.9999999999999999), 0.2);
        const r = query.compute(c0, c1);
        expect(r.distance).toBeCloseTo(0.001, 12);
        expectValidResult(c0, c1, r);
        expect(r.equidistant).toBe(true);

        // The parallel-planes branch. The transformed normal rounds to one
        // ulp above 1, so 'circle0.normal[2] < 1' is false. Before the fix
        // both closest points were on the (1,1,1) axis, 2 away from their
        // own circles.
        const n = v(1, 1, 1);
        normalize(n);
        for (const sign of [1, -1]) {
            const a0 = Circle3.fromCenterNormalRadius(v(0, 0, 0), n, 1);
            const a1 = Circle3.fromCenterNormalRadius(mul(0.5, n),
                mul(sign, n), 2);
            const ra = query.compute(a0, a1);
            expect(ra.distance).toBeCloseTo(Math.sqrt(1.25), 12);
            expectValidResult(a0, a1, ra);
            expect(ra.equidistant).toBe(true);
        }

        // A deterministic family of coaxial pairs: separated, nested and
        // equal-radius, on axes whose rotation to the z-axis leaves the
        // transformed normal on either side of 1.
        const axes = [v(1, 1, 1), v(1, 2, 2), v(0, 3, 4), v(1, 0, 1),
            v(2, 3, 6), v(1, 1, 0), v(1, 2, 3), v(0, 0, 1)];
        for (const axis of axes) {
            const u = axis.clone();
            normalize(u);
            for (const offset of [0, 0.25, 0.5, 1, 1.5, 2]) {
                for (const [r0, r1] of [[1, 2], [0.5, 1.5], [1, 1],
                    [2, 0.5]] as Array<[number, number]>) {
                    const b0 = Circle3.fromCenterNormalRadius(v(1, -2, 3), u,
                        r0);
                    const b1 = Circle3.fromCenterNormalRadius(
                        add(v(1, -2, 3), mul(offset, u)), u, r1);
                    const rb = query.compute(b0, b1);
                    expectValidResult(b0, b1, rb);
                    expect(rb.distance).toBeCloseTo(
                        Math.sqrt((r1 - r0) ** 2 + offset ** 2), 12);
                }
            }
        }
    });

    it('records the subnormal-offset breakdown in PrepareCircles', () => {
        // PrepareCircles rotates about the circle1 normal to zero the
        // x-component of the circle0 center, guarding only on
        // 'center[0] != 0'. When the x- and y-components are subnormal the
        // guard passes but sqrt(x*x + y*y) underflows to zero, so sn = x/0
        // and cs = y/0 poison the rotation. Upstream has the same guard and
        // the same division. The port preserves it; update this test if
        // upstream ever tests the length instead.
        const r = query.compute(
            circle(v(5e-324, 0, 0), v(0, 0, 1), 0.2),
            circle(v(0, 0, 0), v(0, 0, 1), 0.2));
        expect(Number.isNaN(r.distance)).toBe(true);

        // The same configuration with a representable offset is fine.
        const ok = query.compute(
            circle(v(1e-6, 0, 0), v(0, 0, 1), 0.2),
            circle(v(0, 0, 0), v(0, 0, 1), 0.2));
        expect(ok.distance).toBeCloseTo(0, 9);
    });
});
