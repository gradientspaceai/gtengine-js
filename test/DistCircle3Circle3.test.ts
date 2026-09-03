import { describe, it, expect } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { DistCircle3Circle3 } from '../src/DistCircle3Circle3.js';
import {
    Vector, add, dot, length, mul, normalize, sub
} from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

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
