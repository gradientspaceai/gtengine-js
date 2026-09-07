import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Hypersphere } from '../src/Hypersphere.js';
import {
    IntrAlignedBox3Sphere3FI,
    IntrAlignedBox3Sphere3FIResultType
} from '../src/IntrAlignedBox3Sphere3.js';
import {
    IntrOrientedBox3Sphere3TI,
    IntrOrientedBox3Sphere3FI,
    defaultIntrOrientedBox3Sphere3TIResult
} from '../src/IntrOrientedBox3Sphere3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { IntrRay3AlignedBox3TI } from '../src/IntrRay3AlignedBox3.js';
import { Ray } from '../src/Ray.js';
import {
    alignedBox, check, expectClose, expectVectorClose, fc,
    orientedBox as arbOrientedBox, positive, rotationFrame, seededRandom,
    wellScaledVector
} from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrOrientedBox3Sphere3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrOrientedBox3Sphere3 verification', () => {
    const tiQ = new IntrOrientedBox3Sphere3TI();
    const fiQ = new IntrOrientedBox3Sphere3FI();
    const pbQuery = new DistPointOrientedBox();

    const arbBox = arbOrientedBox(3);
    const arbSphere = fc.tuple(wellScaledVector(3), positive(3))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    // The distance from the moving sphere center to the box at time t. It is
    // a convex function of t (distance to a convex set composed with an
    // affine map), which is what makes the ternary search below valid.
    function gap(box: OrientedBox, C: Vector, V: Vector, t: number): number {
        return pbQuery.compute(add(C, mul(t, V)), box).distance;
    }

    // The first time in [0, tMax] at which the gap reaches the radius, or
    // null when the gap stays above it. Uses the convexity of the gap.
    function firstContactTime(box: OrientedBox, C: Vector, V: Vector,
        radius: number, tMax: number): number | null {
        let lo = 0;
        let hi = tMax;
        for (let i = 0; i < 200; ++i) {
            const m0 = lo + (hi - lo) / 3;
            const m1 = hi - (hi - lo) / 3;
            if (gap(box, C, V, m0) <= gap(box, C, V, m1)) { hi = m1; }
            else { lo = m0; }
        }
        const tMin = (lo + hi) / 2;
        if (gap(box, C, V, tMin) > radius) { return null; }
        if (gap(box, C, V, 0) <= radius) { return 0; }
        let a = 0;
        let b = tMin;
        for (let i = 0; i < 200; ++i) {
            const m = (a + b) / 2;
            if (gap(box, C, V, m) > radius) { a = m; } else { b = m; }
        }
        return (a + b) / 2;
    }

    it('TI agrees with the point-box distance', () => {
        check(fc.tuple(arbBox, arbSphere), ([box, sphere]) => {
            const d = pbQuery.compute(sphere.center, box).sqrDistance;
            expect(tiQ.test(box, sphere).intersect)
                .toBe(d <= sphere.radius * sphere.radius);
        });
    });

    it('the dynamic query agrees with a bisection on the gap function', () => {
        check(fc.tuple(arbBox, wellScaledVector(3), arbSphere,
            wellScaledVector(3)),
            ([box, boxVel, sphere, sphereVel]) => {
                const V = sub(sphereVel, boxVel);
                if (dot(V, V) < 1e-4) { return; }
                const r = fiQ.find(box, boxVel, sphere, sphereVel);
                const scale = 1 + sphere.radius + length(box.extent);
                const t0 = firstContactTime(box, sphere.center, V,
                    sphere.radius, 1e4);

                if (r.intersectionType
                    === IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping) {
                    expect(r.contactTime).toBe(0);
                    expect(gap(box, sphere.center, V, 0))
                        .toBeLessThanOrEqual(sphere.radius);
                    return;
                }

                if (r.intersectionType
                    === IntrAlignedBox3Sphere3FIResultType.contact) {
                    expect(r.contactTime).toBeGreaterThan(0);
                    expect(t0).not.toBeNull();
                    // Upstream's DoQueryRayRoundedFace accepts the first
                    // rounded-edge probe and never tries the other edge, so
                    // the reported contact can be late (see gtengine-js #458
                    // item 5); the port preserves that. The reported time is
                    // therefore >= the true first contact and the sphere is
                    // in contact or slightly penetrating at it.
                    expect(gap(box, sphere.center, V, r.contactTime))
                        .toBeLessThanOrEqual(sphere.radius + 1e-6 * scale);
                    expect(r.contactTime)
                        .toBeGreaterThanOrEqual((t0 as number) - 1e-6 * scale);
                    return;
                }

                // No contact reported: the gap must stay above the radius on
                // the whole search window.
                expect(t0).toBeNull();
            }, 100);
    });

    it('the added super-box early exit never causes a missed contact', () => {
        // The port delegates to IntrAlignedBox3Sphere3FI, which runs a
        // ray-versus-expanded-box rejection test that upstream's oriented-box
        // path does not have. The expanded box contains the Minkowski sum of
        // the box and the sphere, so the test is conservative; this property
        // checks that whenever the query reports no contact although the gap
        // function reaches the radius, the early exit had ACCEPTED the
        // configuration, that is, the miss comes from the shared DoQuery and
        // not from the delegation. (The known DoQuery miss is pinned below.)
        const rbQuery = new IntrRay3AlignedBox3TI();
        const rnd = seededRandom(0xb073be);
        let numContacts = 0;
        let numMissed = 0;
        for (let iter = 0; iter < 400; ++iter) {
            const axis = [Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0]), Vector.fromArray([0, 0, 1])];
            const box = OrientedBox.fromCenterAxisExtent(
                Vector.fromArray([rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1]),
                axis,
                Vector.fromArray([0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd()]));
            const sphere = Hypersphere.fromCenterRadius(
                Vector.fromArray([4 + rnd() * 2, rnd() * 4 - 2, rnd() * 4 - 2]),
                0.2 + rnd());
            const V = Vector.fromArray([-1, rnd() * 0.4 - 0.2,
                rnd() * 0.4 - 0.2]);
            const r = fiQ.find(box, Vector.zero(3), sphere, V);
            const t0 = firstContactTime(box, sphere.center, V, sphere.radius,
                1e3);
            if (t0 === null) {
                expect(r.intersectionType)
                    .toBe(IntrAlignedBox3Sphere3FIResultType.noContact);
                continue;
            }
            ++numContacts;
            if (r.intersectionType
                !== IntrAlignedBox3Sphere3FIResultType.noContact) {
                continue;
            }
            ++numMissed;
            // The early exit accepted, so the delegation is faithful.
            const superMax = Vector.fromArray([
                box.extent.values[0] + sphere.radius,
                box.extent.values[1] + sphere.radius,
                box.extent.values[2] + sphere.radius]);
            const superBox = AlignedBox.fromMinMax(mul(-1, superMax), superMax);
            const ray = Ray.fromOriginDirection(sub(sphere.center, box.center),
                V);
            expect(rbQuery.test(ray, superBox).intersect).toBe(true);
        }
        expect(numContacts).toBeGreaterThan(50);
        // The shared DoQuery misses a small number of grazing contacts; see
        // the pinned case below.
        expect(numMissed).toBeLessThanOrEqual(2);
    }, 30000);

    it('pins the upstream DoQuery miss inherited from IntrAlignedBox3Sphere3', () => {
        // Upstream's DoQueryRayRoundedFace accepts the first rounded-edge
        // probe and never tries the other edge (gtengine-js #458 item 5).
        // Besides reporting some contacts late, it can miss one entirely: in
        // the configuration below the sphere penetrates the box by 0.0057 at
        // the closest approach yet 'noContact' is reported. The super-box
        // early exit that the port adds ACCEPTS this configuration, so the
        // miss is upstream's, not the port's. Preserved and pinned.
        const axis = [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])];
        const box = OrientedBox.fromCenterAxisExtent(
            Vector.fromArray([0.090417948551476, -0.5309329181909561,
                0.9111482501029968]),
            axis,
            Vector.fromArray([0.3892387102358043, 0.3108668552711606,
                0.5541517764329911]));
        const sphere = Hypersphere.fromCenterRadius(
            Vector.fromArray([5.136491211596876, -1.6951057985424995,
                1.702867484651506]), 0.2129704826977104);
        const V = Vector.fromArray([-1, 0.16651661535725,
            -0.005579612776637072]);

        const r = fiQ.find(box, Vector.zero(3), sphere, V);
        expect(r.intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.noContact);

        // The gap really does drop below the radius.
        const t0 = firstContactTime(box, sphere.center, V, sphere.radius, 1e3);
        expect(t0).not.toBeNull();
        expect(gap(box, sphere.center, V, 5.4365))
            .toBeLessThan(sphere.radius - 5e-3);

        // The port's added early exit is not the cause.
        const superMax = Vector.fromArray([
            box.extent.values[0] + sphere.radius,
            box.extent.values[1] + sphere.radius,
            box.extent.values[2] + sphere.radius]);
        const superBox = AlignedBox.fromMinMax(mul(-1, superMax), superMax);
        expect(new IntrRay3AlignedBox3TI().test(
            Ray.fromOriginDirection(sub(sphere.center, box.center), V),
            superBox).intersect).toBe(true);
    });

    it('the contact point is on the sphere and on the box', () => {
        check(fc.tuple(arbBox, wellScaledVector(3), arbSphere,
            wellScaledVector(3)),
            ([box, boxVel, sphere, sphereVel]) => {
                const r = fiQ.find(box, boxVel, sphere, sphereVel);
                if (r.intersectionType
                    !== IntrAlignedBox3Sphere3FIResultType.contact) {
                    return;
                }
                const V = sub(sphereVel, boxVel);
                const center = add(sphere.center, mul(r.contactTime, V));
                const scale = 1 + sphere.radius + length(box.extent)
                    + Math.abs(r.contactTime) * length(V);
                // The contact point is on the box.
                expect(pbQuery.compute(r.contactPoint, box).distance)
                    .toBeLessThan(1e-6 * scale);
                // ... and no farther from the moved sphere center than the
                // radius (upstream's late contact can put it slightly inside).
                expect(length(sub(r.contactPoint, center)))
                    .toBeLessThanOrEqual(sphere.radius + 1e-6 * scale);
            }, 100);
    });

    it('is equivariant under a rigid motion of the whole configuration', () => {
        check(fc.tuple(arbBox, wellScaledVector(3), arbSphere,
            wellScaledVector(3), rotationFrame(3), wellScaledVector(3)),
            ([box, boxVel, sphere, sphereVel, R, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(R[0], v), dot(R[1], v), dot(R[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const box2 = OrientedBox.fromCenterAxisExtent(map(box.center),
                    box.axis.map(rot), box.extent);
                const sphere2 = Hypersphere.fromCenterRadius(
                    map(sphere.center), sphere.radius);
                const r0 = fiQ.find(box, boxVel, sphere, sphereVel);
                const r1 = fiQ.find(box2, rot(boxVel), sphere2, rot(sphereVel));
                expect(r1.intersectionType).toBe(r0.intersectionType);
                if (r0.intersectionType
                    === IntrAlignedBox3Sphere3FIResultType.contact) {
                    const scale = 1 + Math.abs(r0.contactTime);
                    expectClose(r1.contactTime, r0.contactTime, 1e-6 * scale,
                        1e-6);
                    expectVectorClose(r1.contactPoint, map(r0.contactPoint),
                        1e-5, 1e-5);
                }
            }, 100);
    });

    it('matches the aligned-box query when the box axes are the identity', () => {
        check(fc.tuple(alignedBox(3), wellScaledVector(3), arbSphere,
            wellScaledVector(3)),
            ([abox, boxVel, sphere, sphereVel]) => {
                const center = mul(0.5, add(abox.min, abox.max));
                const extent = mul(0.5, sub(abox.max, abox.min));
                const obox = OrientedBox.fromCenterAxisExtent(center,
                    [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                        Vector.fromArray([0, 0, 1])], extent);
                const a = new IntrAlignedBox3Sphere3FI()
                    .find(abox, boxVel, sphere, sphereVel);
                const o = fiQ.find(obox, boxVel, sphere, sphereVel);
                expect(o.intersectionType).toBe(a.intersectionType);
                if (a.intersectionType
                    !== IntrAlignedBox3Sphere3FIResultType.noContact) {
                    expectClose(o.contactTime, a.contactTime, 1e-9, 1e-9);
                    expectVectorClose(o.contactPoint, a.contactPoint, 1e-9,
                        1e-9);
                }
            }, 100);
    });

    it('reports the head-on contact time for a rotated box', () => {
        // A unit cube rotated 45 degrees about the z-axis; a unit sphere
        // approaches along -x from x = 10. The nearest cube point on the +x
        // side is the edge at x = sqrt(2)/2, so contact is at
        // x = sqrt(2)/2 + 1.
        const s = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([s, s, 0]), Vector.fromArray([-s, s, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([0.5, 0.5, 0.5]));
        const sphere = Hypersphere.fromCenterRadius(
            Vector.fromArray([10, 0, 0]), 1);
        const r = fiQ.find(box, Vector.zero(3), sphere,
            Vector.fromArray([-1, 0, 0]));
        expect(r.intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.contact);
        expect(r.contactTime).toBeCloseTo(10 - (Math.SQRT2 / 2 + 1), 9);
        expect(r.contactPoint.values[0]).toBeCloseTo(Math.SQRT2 / 2, 9);
        expect(r.contactPoint.values[1]).toBeCloseTo(0, 9);
        expect(r.contactPoint.values[2]).toBeCloseTo(0, 9);
    });

    it('reports an initial overlap with the sphere center as contact point', () => {
        const box = OrientedBox.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([2, 1, 1]));
        const sphere = Hypersphere.fromCenterRadius(
            Vector.fromArray([0.5, 0.25, 0]), 0.1);
        const r = fiQ.find(box, Vector.zero(3), sphere, Vector.zero(3));
        expect(r.intersectionType)
            .toBe(IntrAlignedBox3Sphere3FIResultType.initiallyOverlapping);
        expect(r.contactTime).toBe(0);
        expectVectorClose(r.contactPoint, sphere.center, 1e-12, 1e-12);
    });
});
