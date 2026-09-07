import { describe, it, expect } from 'vitest';
import { Circle3 } from '../src/Circle3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrPlane3Circle3TI,
    IntrPlane3Circle3FI,
    defaultIntrPlane3Circle3FIResult,
    intrPlane3Circle3InfinitePoints
} from '../src/IntrPlane3Circle3.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, fc, plane as arbPlane, positive, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

function circle(center: number[], normal: number[], radius: number): Circle3 {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Circle3.fromCenterNormalRadius(Vector.fromArray(center), n, radius);
}

const ti = new IntrPlane3Circle3TI();
const fi = new IntrPlane3Circle3FI();

describe('IntrPlane3Circle3', () => {
    it('has an invalid circle in the default result', () => {
        const result = defaultIntrPlane3Circle3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.circle.center.values).toEqual([0, 0, 0]);
        expect(result.circle.normal.values).toEqual([0, 0, 0]);
        expect(result.circle.radius).toBe(0);
    });

    it('finds two crossings of a circle cut by a plane', () => {
        // The unit circle in the xy-plane cut by the plane x = 0.5.
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [0.5, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        const y = Math.sqrt(1 - 0.25);
        for (const point of result.point) {
            expect(point.values[0]).toBeCloseTo(0.5, 12);
            expect(Math.abs(point.values[1])).toBeCloseTo(y, 12);
            expect(point.values[2]).toBeCloseTo(0, 12);
        }
        expect(result.point[0].values[1])
            .not.toBeCloseTo(result.point[1].values[1], 6);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports a tangent plane as a single point', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [1, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(result.point[1].values).toEqual(result.point[0].values);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports the whole circle when the planes coincide', () => {
        const C = circle([1, 2, 0], [0, 0, 1], 3);
        const P = plane([0, 0, 1], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(intrPlane3Circle3InfinitePoints);
        expect(result.circle.center.values).toEqual([1, 2, 0]);
        expect(result.circle.radius).toBe(3);
        expect(result.point[0].values).toEqual([0, 0, 0]);
        expect(ti.test(P, C).intersect).toBe(true);
    });

    it('reports no intersection for parallel disjoint planes', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([0, 0, 1], [0, 0, 1]);
        expect(fi.find(P, C).intersect).toBe(false);
        expect(ti.test(P, C).intersect).toBe(false);
    });

    it('reports no intersection when the cut line misses the circle', () => {
        const C = circle([0, 0, 0], [0, 0, 1], 1);
        const P = plane([1, 0, 0], [2, 0, 0]);
        expect(fi.find(P, C).intersect).toBe(false);
        expect(ti.test(P, C).intersect).toBe(false);
    });

    it('agrees with the TI query and geometry on random inputs', () => {
        let state = 31415926;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numTwo = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const C = circle([rand() * 2, rand() * 2, rand() * 2],
                [rand(), rand(), rand() + 0.001],
                0.5 + Math.abs(rand()) * 2);
            const P = plane([rand(), rand(), rand() + 0.001],
                [rand() * 2, rand() * 2, rand() * 2]);
            const result = fi.find(P, C);
            expect(ti.test(P, C).intersect).toBe(result.intersect);

            if (result.numIntersections === 2) {
                ++numTwo;
            }
            if (result.numIntersections === 1
                || result.numIntersections === 2) {
                for (let i = 0; i < 2; ++i) {
                    const X = result.point[i];
                    // The point is on the plane, on the circle plane and at
                    // the circle radius from the circle center.
                    expect(dot(P.normal, X) - P.constant).toBeCloseTo(0, 8);
                    const d = sub(X, C.center);
                    expect(dot(C.normal, d)).toBeCloseTo(0, 8);
                    expect(length(d)).toBeCloseTo(C.radius, 7);
                }
            }
        }
        expect(numTwo).toBeGreaterThan(50);
    });

    it('accepts a circle whose normal is not axis aligned', () => {
        // A circle in the plane x + y + z = 0 cut by the plane z = 0. The
        // intersection points are on the line x + y = 0, z = 0.
        const C = circle([0, 0, 0], [1, 1, 1], 2);
        const P = plane([0, 0, 1], [0, 0, 0]);
        const result = fi.find(P, C);
        expect(result.numIntersections).toBe(2);
        for (const X of result.point) {
            expect(X.values[2]).toBeCloseTo(0, 12);
            expect(X.values[0] + X.values[1]).toBeCloseTo(0, 10);
            expect(length(sub(X, C.center))).toBeCloseTo(2, 10);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrPlane3Circle3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrPlane3Circle3 verification', () => {
    const tiQ = new IntrPlane3Circle3TI();
    const fiQ = new IntrPlane3Circle3FI();

    const arbCircle = fc.tuple(wellScaledVector(3), unitVector(3), positive(5))
        .map(([c, n, r]) => Circle3.fromCenterNormalRadius(c, n, r));

    // IntrPlane3Plane3 decides 'parallel' with the exact test |Dot(N0,N1)|
    // >= 1, so the coplanar and parallel branches are reached only when the
    // shared normal has Dot(N,N) exactly 1. About half of the normalized
    // random vectors do; these are exact.
    const exactUnitNormals = [
        Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
        Vector.fromArray([0, 0, 1]), Vector.fromArray([0.6, 0.8, 0]),
        Vector.fromArray([0, 0.6, 0.8]), Vector.fromArray([0.8, 0, 0.6]),
        Vector.fromArray([-0.6, 0, 0.8])
    ];
    const arbExactCircle = fc.tuple(wellScaledVector(3),
        fc.constantFrom(...exactUnitNormals), positive(5))
        .map(([c, n, r]) => Circle3.fromCenterNormalRadius(c, n, r));

    // The circle plane and the input plane are transverse by a safe margin,
    // which keeps the plane-plane line well conditioned.
    const transverse = (P: Hyperplane, C: Circle3) =>
        Math.abs(dot(P.normal, C.normal)) < 1 - 1e-3;

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbPlane(3), arbCircle), ([P, C]) => {
            expect(fiQ.find(P, C).intersect).toBe(tiQ.test(P, C).intersect);
        });
    });

    it('the FI points lie on the plane and on the circle', () => {
        check(fc.tuple(arbPlane(3), arbCircle), ([P, C]) => {
            if (!transverse(P, C)) { return; }
            const r = fiQ.find(P, C);
            if (!r.intersect
                || r.numIntersections === intrPlane3Circle3InfinitePoints) {
                return;
            }
            const scale = 1 + C.radius + length(C.center);
            for (let i = 0; i < r.numIntersections; ++i) {
                const X = r.point[i];
                // On the input plane.
                expectClose(dot(P.normal, X) - P.constant, 0, 1e-8 * scale, 0);
                // On the circle: in the circle plane and at the radius.
                const d = sub(X, C.center);
                expectClose(dot(C.normal, d), 0, 1e-8 * scale, 0);
                expectClose(length(d), C.radius, 1e-8 * scale, 1e-8);
            }
        });
    });

    it('reports a coplanar circle as the whole circle', () => {
        check(arbExactCircle, C => {
            const P = Hyperplane.fromNormalOrigin(C.normal, C.center);
            const r = fiQ.find(P, C);
            expect(r.intersect).toBe(true);
            expect(r.numIntersections).toBe(intrPlane3Circle3InfinitePoints);
            expect(r.circle.center.values).toEqual(C.center.values);
            expect(r.circle.normal.values).toEqual(C.normal.values);
            expect(r.circle.radius).toBe(C.radius);
            expect(r.point[0].values).toEqual([0, 0, 0]);
            expect(r.point[1].values).toEqual([0, 0, 0]);
            expect(tiQ.test(P, C).intersect).toBe(true);
            // The returned circle is a copy, not an alias of the input.
            r.circle.radius = -1;
            expect(C.radius).not.toBe(-1);
        });
    });

    it('reports a parallel offset plane as empty', () => {
        check(arbExactCircle, C => {
            const P = Hyperplane.fromNormalOrigin(C.normal,
                add(C.center, mul(1.5, C.normal)));
            const r = fiQ.find(P, C);
            expect(r.intersect).toBe(false);
            expect(r.numIntersections).toBe(0);
            expect(tiQ.test(P, C).intersect).toBe(false);
        });
    });

    it('agrees with sampling of the circle when the plane cuts it', () => {
        // The plane cuts the circle when Dot(N,X) - c changes sign along the
        // circle. Sampling proves an intersection; the query must report at
        // least one point in that case.
        check(fc.tuple(arbPlane(3), arbCircle), ([P, C]) => {
            if (!transverse(P, C)) { return; }
            const helper = Math.abs(C.normal.values[0]) < 0.9
                ? Vector.fromArray([1, 0, 0]) : Vector.fromArray([0, 1, 0]);
            const u = cross(C.normal, helper);
            normalize(u);
            const v = cross(C.normal, u);
            let lo = Number.POSITIVE_INFINITY;
            let hi = Number.NEGATIVE_INFINITY;
            for (let i = 0; i < 360; ++i) {
                const a = (2 * Math.PI * i) / 360;
                const X = add(C.center, add(mul(C.radius * Math.cos(a), u),
                    mul(C.radius * Math.sin(a), v)));
                const sd = dot(P.normal, X) - P.constant;
                lo = Math.min(lo, sd);
                hi = Math.max(hi, sd);
            }
            const r = fiQ.find(P, C);
            const scale = 1 + C.radius;
            if (lo < -1e-9 * scale && hi > 1e-9 * scale) {
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(2);
            }
            if (lo > 1e-6 * scale || hi < -1e-6 * scale) {
                expect(r.intersect).toBe(false);
            }
        }, 60);
    });

    it('reports a tangent configuration as a single (doubled) point', () => {
        // The plane z = 0 is tangent to the unit circle in the xz-plane
        // centered at (0,0,1).
        const C = Circle3.fromCenterNormalRadius(Vector.fromArray([0, 0, 1]),
            Vector.fromArray([0, 1, 0]), 1);
        const P = Hyperplane.fromNormalOrigin(Vector.fromArray([0, 0, 1]),
            Vector.zero(3));
        const r = fiQ.find(P, C);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(1);
        expect(r.point[0].values[0]).toBeCloseTo(0, 12);
        expect(r.point[0].values[2]).toBeCloseTo(0, 12);
        expect(r.point[1].values).toEqual(r.point[0].values);
        // point[1] is a copy, not an alias.
        r.point[1].values[0] = 42;
        expect(r.point[0].values[0]).not.toBe(42);
        expect(r.circle.radius).toBe(0);
        expect(tiQ.test(P, C).intersect).toBe(true);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbPlane(3), arbCircle, rotationFrame(3),
            wellScaledVector(3)),
            ([P, C, R, t]) => {
                const rot = (x: Vector) => Vector.fromArray([
                    dot(R[0], x), dot(R[1], x), dot(R[2], x)]);
                const map = (x: Vector) => add(rot(x), t);
                const P2 = Hyperplane.fromNormalOrigin(rot(P.normal),
                    map(P.origin));
                const C2 = Circle3.fromCenterNormalRadius(map(C.center),
                    rot(C.normal), C.radius);
                if (!transverse(P, C)) { return; }
                const r0 = fiQ.find(P, C);
                const r1 = fiQ.find(P2, C2);
                if (r0.numIntersections === 1) { return; }  // tangency
                expect(r1.intersect).toBe(r0.intersect);
                expect(r1.numIntersections).toBe(r0.numIntersections);
            });
    });
    it('pins the near-coplanar conditioning inherited from IntrPlane3Plane3', () => {
        // IntrPlane3Plane3 classifies the planes as parallel only when
        // |Dot(N0,N1)| >= 1 exactly. A normalized normal whose self-dot is
        // one ulp below 1 makes two IDENTICAL planes transverse, and the
        // intersection line is then built with invDet = 1/(1 - d*d) ~ 1e16.
        // This query inherits that: the circle is reported as two points far
        // from the circle rather than as the whole circle. Upstream has the
        // same behaviour; the test pins it so that a future change is noticed.
        const n = Vector.fromArray([0.9995064129969172, 0, 0.03141544807949855]);
        expect(dot(n, n)).toBeLessThan(1);
        const C = Circle3.fromCenterNormalRadius(Vector.zero(3), n, 1);
        const P = Hyperplane.fromNormalOrigin(n, Vector.zero(3));
        const r = fiQ.find(P, C);
        expect(r.numIntersections).not.toBe(intrPlane3Circle3InfinitePoints);
        // With an exactly unit-length normal the same configuration is
        // reported correctly.
        const m = Vector.fromArray([0.6, 0.8, 0]);
        expect(dot(m, m)).toBe(1);
        const C2 = Circle3.fromCenterNormalRadius(Vector.zero(3), m, 1);
        const P2 = Hyperplane.fromNormalOrigin(m, Vector.zero(3));
        expect(fiQ.find(P2, C2).numIntersections)
            .toBe(intrPlane3Circle3InfinitePoints);
    });
});
