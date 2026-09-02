import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Hypersphere } from '../src/Hypersphere';
import {
    IntrAlignedBox3Sphere3FI,
    IntrAlignedBox3Sphere3FIResultType
} from '../src/IntrAlignedBox3Sphere3';
import {
    IntrOrientedBox3Sphere3TI,
    IntrOrientedBox3Sphere3FI,
    defaultIntrOrientedBox3Sphere3TIResult
} from '../src/IntrOrientedBox3Sphere3';
import { OrientedBox } from '../src/OrientedBox';
import { Vector, add, length, mul, sub } from '../src/Vector';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// A box rotated by 'angle' about the z-axis.
function obox(center: number[], angle: number, extent: number[]): OrientedBox {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(c, s, 0), vec(-s, c, 0), vec(0, 0, 1)],
        Vector.fromArray(extent));
}

function sphere(center: number[], radius: number): Hypersphere {
    return Hypersphere.fromCenterRadius(Vector.fromArray(center), radius);
}

function rotZ(v: Vector, angle: number): Vector {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return vec(c * v.values[0] - s * v.values[1],
        s * v.values[0] + c * v.values[1], v.values[2]);
}

const ti = new IntrOrientedBox3Sphere3TI();
const fi = new IntrOrientedBox3Sphere3FI();

describe('IntrOrientedBox3Sphere3', () => {
    it('defaults to no intersection', () => {
        expect(defaultIntrOrientedBox3Sphere3TIResult().intersect).toBe(false);
    });

    it('tests solid overlap of a rotated box and a sphere', () => {
        // A unit cube rotated 45 degrees about z reaches sqrt(2) along x.
        const B = obox([0, 0, 0], Math.PI / 4, [1, 1, 1]);
        expect(ti.test(B, sphere([2, 0, 0], 0.5)).intersect).toBe(false);
        expect(ti.test(B, sphere([2, 0, 0], 0.6)).intersect).toBe(true);
        // Tangency counts as an intersection.
        expect(ti.test(B, sphere([Math.SQRT2 + 1, 0, 0], 1)).intersect)
            .toBe(true);
        // A sphere strictly inside the box intersects.
        expect(ti.test(B, sphere([0, 0, 0], 0.1)).intersect).toBe(true);
    });

    it('matches the aligned-box query when the box axes are the identity', () => {
        const abQuery = new IntrAlignedBox3Sphere3FI();
        const alignedBox = AlignedBox.fromMinMax(vec(-1, -2, -3),
            vec(1, 2, 3));
        const orientedBox = obox([0, 0, 0], 0, [1, 2, 3]);
        const S = sphere([6, 0, 0], 1);
        const V = vec(-1, 0, 0);
        const expected = abQuery.find(alignedBox, Vector.zero(3), S, V);
        const actual = fi.find(orientedBox, Vector.zero(3), S, V);
        expect(actual.intersectionType).toBe(expected.intersectionType);
        expect(actual.contactTime).toBeCloseTo(expected.contactTime, 12);
        for (let i = 0; i < 3; ++i) {
            expect(actual.contactPoint.values[i])
                .toBeCloseTo(expected.contactPoint.values[i], 12);
        }
    });

    it('finds a first contact time on a box face', () => {
        const B = obox([0, 0, 0], 0, [1, 2, 3]);
        const result = fi.find(B, Vector.zero(3), sphere([6, 0, 0], 1),
            vec(-1, 0, 0));
        expect(result.intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.contact);
        expect(result.contactTime).toBeCloseTo(4, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(0, 12);
    });

    it('reports initial overlap and no contact', () => {
        const B = obox([0, 0, 0], 0.4, [1, 1, 1]);
        expect(fi.find(B, Vector.zero(3), sphere([0, 0, 0], 0.5),
            vec(1, 1, 1)).intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping);
        expect(fi.find(B, Vector.zero(3), sphere([20, 20, 20], 0.5),
            vec(1, 1, 1)).intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.noContact);
    });

    it('is equivariant under rotation and reports valid contacts', () => {
        let state = 777333;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numContacts = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const extent = [0.3 + Math.abs(rand()), 0.3 + Math.abs(rand()),
                0.3 + Math.abs(rand())];
            const boxCenter = [rand() * 2, rand() * 2, rand() * 2];
            const boxVelocity = vec(rand(), rand(), rand());
            const sphereCenter = [rand() * 4, rand() * 4, rand() * 4];
            const sphereVelocity = vec(rand() * 2, rand() * 2, rand() * 2);
            const radius = 0.2 + Math.abs(rand());

            const base = fi.find(obox(boxCenter, 0, extent), boxVelocity,
                sphere(sphereCenter, radius), sphereVelocity);

            const angle = rand() * Math.PI;
            const rBox = OrientedBox.fromCenterAxisExtent(
                rotZ(Vector.fromArray(boxCenter), angle),
                [rotZ(vec(1, 0, 0), angle), rotZ(vec(0, 1, 0), angle),
                    vec(0, 0, 1)],
                Vector.fromArray(extent));
            const rResult = fi.find(rBox, rotZ(boxVelocity, angle),
                Hypersphere.fromCenterRadius(
                    rotZ(Vector.fromArray(sphereCenter), angle), radius),
                rotZ(sphereVelocity, angle));

            expect(rResult.intersectionType).toBe(base.intersectionType);
            if (base.intersectionType
                !== IntrAlignedBox3Sphere3FIResultType.noContact) {
                expect(rResult.contactTime).toBeCloseTo(base.contactTime, 9);
                const expectedPoint = rotZ(base.contactPoint, angle);
                for (let i = 0; i < 3; ++i) {
                    expect(rResult.contactPoint.values[i])
                        .toBeCloseTo(expectedPoint.values[i], 8);
                }
            }

            if (base.intersectionType
                === IntrAlignedBox3Sphere3FIResultType.contact) {
                ++numContacts;
                // The query works in the frame of the (static) box with the
                // sphere moving at the relative velocity, so the contact
                // point is expressed relative to the initial box position.
                const relVelocity = sub(sphereVelocity, boxVelocity);
                const movedCenter = add(Vector.fromArray(sphereCenter),
                    mul(base.contactTime, relVelocity));
                expect(length(sub(base.contactPoint, movedCenter)))
                    .toBeCloseTo(radius, 8);
                // The contact point is on the box boundary.
                const local = sub(base.contactPoint,
                    Vector.fromArray(boxCenter));
                let maxSlack = -Infinity;
                for (let i = 0; i < 3; ++i) {
                    const slack = Math.abs(local.values[i]) - extent[i];
                    expect(slack).toBeLessThan(1e-7);
                    maxSlack = Math.max(maxSlack, slack);
                }
                expect(maxSlack).toBeCloseTo(0, 7);
            }
        }
        expect(numContacts).toBeGreaterThan(5);
    });
});
