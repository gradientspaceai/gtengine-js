import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrAlignedBox2Circle2FI,
    IntrAlignedBox2Circle2FIResultType
} from '../src/IntrAlignedBox2Circle2.js';
import {
    IntrOrientedBox2Circle2TI,
    IntrOrientedBox2Circle2FI,
    defaultIntrOrientedBox2Circle2TIResult
} from '../src/IntrOrientedBox2Circle2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import {
    alignedBox, check, expectClose, expectVectorClose, fc,
    orientedBox as arbOrientedBox, positive, rotationFrame, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function obox(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(c, s), vec(-s, c)], Vector.fromArray(extent));
}

function circle(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

// Rotate a 2D vector by 'angle'.
function rot(v: Vector, angle: number): Vector {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return vec(c * v.values[0] - s * v.values[1],
        s * v.values[0] + c * v.values[1]);
}

const ti = new IntrOrientedBox2Circle2TI();
const fi = new IntrOrientedBox2Circle2FI();

describe('IntrOrientedBox2Circle2', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrOrientedBox2Circle2TIResult().intersect).toBe(false);
    });

    it('tests solid overlap of a rotated box and a circle', () => {
        // A unit square rotated 45 degrees reaches out to sqrt(2) along x.
        const B = obox([0, 0], Math.PI / 4, [1, 1]);
        expect(ti.test(B, circle([2, 0], 0.5)).intersect).toBe(false);
        expect(ti.test(B, circle([2, 0], 0.6)).intersect).toBe(true);
        // The tangent configuration is an intersection.
        expect(ti.test(B, circle([Math.SQRT2 + 1, 0], 1)).intersect).toBe(true);
        // The circle center inside the box is an intersection.
        expect(ti.test(B, circle([0.1, -0.2], 0.01)).intersect).toBe(true);
    });

    it('matches the aligned-box query when the box axes are the identity', () => {
        const abQuery = new IntrAlignedBox2Circle2FI();
        const alignedBox = AlignedBox.fromMinMax(vec(-1, -2), vec(1, 2));
        const orientedBox = obox([0, 0], 0, [1, 2]);
        const C = circle([6, 0], 1);
        const V = vec(-1, 0);
        const expected = abQuery.find(alignedBox, Vector.zero(2), C, V);
        const actual = fi.find(orientedBox, Vector.zero(2), C, V);
        expect(actual.intersectionType).toBe(expected.intersectionType);
        expect(actual.contactTime).toBeCloseTo(expected.contactTime, 12);
        expect(actual.contactPoint.values[0])
            .toBeCloseTo(expected.contactPoint.values[0], 12);
        expect(actual.contactPoint.values[1])
            .toBeCloseTo(expected.contactPoint.values[1], 12);
    });

    it('finds a first contact time on a box face', () => {
        // A circle of radius 1 moving in -x hits the face x = 1 of the box at
        // time 4 (center travels from x = 6 to x = 2).
        const B = obox([0, 0], 0, [1, 2]);
        const result = fi.find(B, Vector.zero(2), circle([6, 0], 1),
            vec(-1, 0));
        expect(result.intersectionType)
            .toBe(IntrAlignedBox2Circle2FIResultType.contact);
        expect(result.contactTime).toBeCloseTo(4, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0, 12);
    });

    it('reports initial overlap and no contact', () => {
        const B = obox([0, 0], 0.3, [1, 1]);
        expect(fi.find(B, Vector.zero(2), circle([0, 0], 0.5), vec(1, 1))
            .intersectionType)
            .toBe(IntrAlignedBox2Circle2FIResultType.initiallyOverlapping);
        expect(fi.find(B, Vector.zero(2), circle([10, 10], 0.5), vec(1, 1))
            .intersectionType)
            .toBe(IntrAlignedBox2Circle2FIResultType.noContact);
    });

    it('is equivariant under rotation of the whole configuration', () => {
        let state = 246810;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numContacts = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const extent = [0.3 + Math.abs(rand()), 0.3 + Math.abs(rand())];
            const boxCenter = [rand() * 2, rand() * 2];
            const boxVelocity = vec(rand(), rand());
            const circleCenter = [rand() * 5, rand() * 5];
            const circleVelocity = vec(rand() * 2, rand() * 2);
            const radius = 0.2 + Math.abs(rand());

            const base = fi.find(obox(boxCenter, 0, extent), boxVelocity,
                circle(circleCenter, radius), circleVelocity);

            // Rotate the box, the circle and both velocities by 'angle'.
            const angle = rand() * Math.PI;
            const rBoxCenter = rot(Vector.fromArray(boxCenter), angle);
            const rBox = OrientedBox.fromCenterAxisExtent(rBoxCenter,
                [rot(vec(1, 0), angle), rot(vec(0, 1), angle)],
                Vector.fromArray(extent));
            const rResult = fi.find(rBox, rot(boxVelocity, angle),
                circle(rot(Vector.fromArray(circleCenter), angle).values.slice(),
                    radius),
                rot(circleVelocity, angle));
            expect(rResult.intersectionType).toBe(base.intersectionType);
            if (base.intersectionType
                !== IntrAlignedBox2Circle2FIResultType.noContact) {
                expect(rResult.contactTime).toBeCloseTo(base.contactTime, 9);
                const expectedPoint = rot(base.contactPoint, angle);
                expect(rResult.contactPoint.values[0])
                    .toBeCloseTo(expectedPoint.values[0], 8);
                expect(rResult.contactPoint.values[1])
                    .toBeCloseTo(expectedPoint.values[1], 8);
            }
            if (base.intersectionType
                === IntrAlignedBox2Circle2FIResultType.contact) {
                ++numContacts;
                // The query works in the frame of the (static) box, with the
                // circle moving at the relative velocity, so the contact
                // point is expressed relative to the initial box position. At
                // the contact time the circle touches the box boundary.
                const relVelocity = sub(circleVelocity, boxVelocity);
                const movedCenter = add(Vector.fromArray(circleCenter),
                    mul(base.contactTime, relVelocity));
                const d = sub(base.contactPoint, movedCenter);
                expect(Math.hypot(d.values[0], d.values[1]))
                    .toBeCloseTo(radius, 8);
                // The contact point is on the box boundary.
                const local = sub(base.contactPoint,
                    Vector.fromArray(boxCenter));
                const ex = Math.abs(local.values[0]) - extent[0];
                const ey = Math.abs(local.values[1]) - extent[1];
                expect(Math.max(ex, ey)).toBeCloseTo(0, 7);
                expect(ex).toBeLessThan(1e-7);
                expect(ey).toBeLessThan(1e-7);
            }
        }
        expect(numContacts).toBeGreaterThan(10);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrOrientedBox2Circle2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrOrientedBox2Circle2 verification', () => {
    const tiQ = new IntrOrientedBox2Circle2TI();
    const fiQ = new IntrOrientedBox2Circle2FI();
    const pbQuery = new DistPointOrientedBox();

    const arbBox = arbOrientedBox(2);
    const arbCircle = fc.tuple(wellScaledVector(2), positive(3))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    // The distance from the moving circle center to the box at time t. It is
    // a convex function of t (distance to a convex set composed with an
    // affine map), which is what makes the ternary search below valid.
    function gap(box: OrientedBox, C: Vector, V: Vector, t: number): number {
        const X = add(C, mul(t, V));
        return pbQuery.compute(X, box).distance;
    }

    it('TI agrees with the point-box distance', () => {
        check(fc.tuple(arbBox, arbCircle), ([box, circle]) => {
            const d = pbQuery.compute(circle.center, box).distance;
            expect(tiQ.test(box, circle).intersect)
                .toBe(d * d <= circle.radius * circle.radius
                    || d <= circle.radius);
        });
    });

    it('the dynamic query agrees with a bisection on the gap function', () => {
        check(fc.tuple(arbBox, wellScaledVector(2), arbCircle,
            wellScaledVector(2)),
            ([box, boxVel, circle, circleVel]) => {
                const V = sub(circleVel, boxVel);
                if (dot(V, V) < 1e-4) { return; }
                const r = fiQ.find(box, boxVel, circle, circleVel);

                if (r.intersectionType
                    === IntrAlignedBox2Circle2FIResultType.initiallyOverlapping) {
                    expect(r.contactTime).toBe(0);
                    expect(gap(box, circle.center, V, 0))
                        .toBeLessThanOrEqual(circle.radius);
                    return;
                }

                // Locate the minimum of the convex gap function on [0, tMax].
                const tMax = 1e4;
                let lo = 0;
                let hi = tMax;
                for (let i = 0; i < 200; ++i) {
                    const m0 = lo + (hi - lo) / 3;
                    const m1 = hi - (hi - lo) / 3;
                    if (gap(box, circle.center, V, m0)
                        <= gap(box, circle.center, V, m1)) {
                        hi = m1;
                    } else {
                        lo = m0;
                    }
                }
                const tMin = (lo + hi) / 2;
                const gMin = gap(box, circle.center, V, tMin);
                const scale = 1 + circle.radius + length(box.extent);

                if (r.intersectionType
                    === IntrAlignedBox2Circle2FIResultType.contact) {
                    expect(r.contactTime).toBeGreaterThan(0);
                    // The reported contact time is a root of gap(t) = radius.
                    expectClose(gap(box, circle.center, V, r.contactTime),
                        circle.radius, 1e-6 * scale, 1e-6);
                    // No earlier contact: the gap is strictly larger before.
                    for (let i = 1; i < 8; ++i) {
                        const t = (r.contactTime * i) / 8;
                        expect(gap(box, circle.center, V, t))
                            .toBeGreaterThan(circle.radius - 1e-6 * scale);
                    }
                } else {
                    // No contact: the smallest gap on [0, tMax] stays above
                    // the radius (the gap grows without bound beyond tMax
                    // because the function is convex and increasing there).
                    expect(gMin).toBeGreaterThan(circle.radius - 1e-7 * scale);
                }
            }, 100);
    });

    it('the contact point is on the circle and on the box', () => {
        check(fc.tuple(arbBox, wellScaledVector(2), arbCircle,
            wellScaledVector(2)),
            ([box, boxVel, circle, circleVel]) => {
                const r = fiQ.find(box, boxVel, circle, circleVel);
                if (r.intersectionType
                    !== IntrAlignedBox2Circle2FIResultType.contact) {
                    return;
                }
                const V = sub(circleVel, boxVel);
                const center = add(circle.center, mul(r.contactTime, V));
                const scale = 1 + circle.radius + length(box.extent)
                    + Math.abs(r.contactTime) * length(V);
                // On the circle at the contact time.
                expectClose(length(sub(r.contactPoint, center)), circle.radius,
                    1e-6 * scale, 1e-6);
                // On the (static) box.
                expect(pbQuery.compute(r.contactPoint, box).distance)
                    .toBeLessThan(1e-6 * scale);
            }, 100);
    });

    it('is equivariant under a rigid motion of the whole configuration', () => {
        check(fc.tuple(arbBox, wellScaledVector(2), arbCircle,
            wellScaledVector(2), rotationFrame(2), wellScaledVector(2)),
            ([box, boxVel, circle, circleVel, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const box2 = OrientedBox.fromCenterAxisExtent(map(box.center),
                    box.axis.map(rot), box.extent);
                const circle2 = Hypersphere.fromCenterRadius(map(circle.center),
                    circle.radius);
                const r0 = fiQ.find(box, boxVel, circle, circleVel);
                const r1 = fiQ.find(box2, rot(boxVel), circle2, rot(circleVel));
                expect(r1.intersectionType).toBe(r0.intersectionType);
                if (r0.intersectionType
                    === IntrAlignedBox2Circle2FIResultType.contact) {
                    const scale = 1 + Math.abs(r0.contactTime);
                    expectClose(r1.contactTime, r0.contactTime, 1e-7 * scale,
                        1e-7);
                    expectVectorClose(r1.contactPoint, map(r0.contactPoint),
                        1e-6, 1e-6);
                }
            }, 100);
    });

    it('matches the aligned-box query when the box axes are the identity', () => {
        check(fc.tuple(alignedBox(2), wellScaledVector(2), arbCircle,
            wellScaledVector(2)),
            ([abox, boxVel, circle, circleVel]) => {
                const center = mul(0.5, add(abox.min, abox.max));
                const extent = mul(0.5, sub(abox.max, abox.min));
                const obox = OrientedBox.fromCenterAxisExtent(center,
                    [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
                    extent);
                const a = new IntrAlignedBox2Circle2FI()
                    .find(abox, boxVel, circle, circleVel);
                const o = fiQ.find(obox, boxVel, circle, circleVel);
                expect(o.intersectionType).toBe(a.intersectionType);
                if (a.intersectionType
                    !== IntrAlignedBox2Circle2FIResultType.noContact) {
                    expectClose(o.contactTime, a.contactTime, 1e-9, 1e-9);
                    expectVectorClose(o.contactPoint, a.contactPoint, 1e-9,
                        1e-9);
                }
            }, 100);
    });

    it('reports an initial overlap with the circle center as contact point', () => {
        const box = OrientedBox.fromCenterAxisExtent(Vector.zero(2),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([2, 1]));
        const circle = Hypersphere.fromCenterRadius(
            Vector.fromArray([0.5, 0.25]), 0.1);
        const r = fiQ.find(box, Vector.zero(2), circle, Vector.zero(2));
        expect(r.intersectionType)
            .toBe(IntrAlignedBox2Circle2FIResultType.initiallyOverlapping);
        expect(r.contactTime).toBe(0);
        expectVectorClose(r.contactPoint, circle.center, 1e-12, 1e-12);
    });

    it('reports the head-on contact time for a rotated box', () => {
        // A unit square rotated by 45 degrees about the origin; a circle of
        // radius 1 approaches along -x from x = 10. The nearest square point
        // on the +x side is the corner at distance sqrt(2)/2 from the center,
        // so contact happens at x = sqrt(2)/2 + 1.
        const s = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(Vector.zero(2),
            [Vector.fromArray([s, s]), Vector.fromArray([-s, s])],
            Vector.fromArray([0.5, 0.5]));
        const circle = Hypersphere.fromCenterRadius(
            Vector.fromArray([10, 0]), 1);
        const r = fiQ.find(box, Vector.zero(2), circle,
            Vector.fromArray([-1, 0]));
        expect(r.intersectionType)
            .toBe(IntrAlignedBox2Circle2FIResultType.contact);
        expect(r.contactTime).toBeCloseTo(10 - (Math.SQRT2 / 2 + 1), 9);
        expect(r.contactPoint.values[0]).toBeCloseTo(Math.SQRT2 / 2, 9);
        expect(r.contactPoint.values[1]).toBeCloseTo(0, 9);
    });
});
