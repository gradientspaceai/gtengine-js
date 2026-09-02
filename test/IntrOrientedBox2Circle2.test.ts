import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Hypersphere } from '../src/Hypersphere';
import {
    IntrAlignedBox2Circle2FI,
    IntrAlignedBox2Circle2FIResultType
} from '../src/IntrAlignedBox2Circle2';
import {
    IntrOrientedBox2Circle2TI,
    IntrOrientedBox2Circle2FI,
    defaultIntrOrientedBox2Circle2TIResult
} from '../src/IntrOrientedBox2Circle2';
import { OrientedBox } from '../src/OrientedBox';
import { Vector, add, mul, sub } from '../src/Vector';

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
