import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox.js';
import { Vector, add, mul, sub } from '../src/Vector.js';
import {
    IntrAlignedBox3Sphere3TI,
    IntrAlignedBox3Sphere3FI,
    IntrAlignedBox3Sphere3FIResultType as Type
} from '../src/IntrAlignedBox3Sphere3.js';
import { dot } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, positive, seededRandom,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function box(min: Vector, max: Vector): AlignedBox {
    return AlignedBox.fromMinMax(min, max);
}

function sphere(cx: number, cy: number, cz: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v3(cx, cy, cz), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The squared distance from the sphere center to the box at time t, where
// the box moves with velocity bv and the sphere center with velocity sv.
function sqrDistanceAtTime(b: AlignedBox, bv: Vector, s: Hypersphere,
    sv: Vector, t: number): number {
    const movedBox = AlignedBox.fromMinMax(add(b.min, mul(bv, t)),
        add(b.max, mul(bv, t)));
    const movedCenter = add(s.center, mul(sv, t));
    return new DistPointAlignedBox().compute(movedCenter, movedBox).sqrDistance;
}

// The distance from a point to the box (zero when the point is inside).
function distanceToBox(b: AlignedBox, p: Vector): number {
    return new DistPointAlignedBox().compute(p, b).distance;
}

const unitBox = box(v3(-1, -1, -1), v3(1, 1, 1));

describe('IntrAlignedBox3Sphere3TI', () => {
    const query = new IntrAlignedBox3Sphere3TI();

    it('detects a sphere inside the box', () => {
        expect(query.test(unitBox, sphere(0, 0, 0, 0.25)).intersect).toBe(true);
    });

    it('detects a box inside the sphere', () => {
        expect(query.test(unitBox, sphere(0, 0, 0, 10)).intersect).toBe(true);
    });

    it('detects face overlap and separation', () => {
        expect(query.test(unitBox, sphere(1.5, 0, 0, 1)).intersect).toBe(true);
        expect(query.test(unitBox, sphere(2.5, 0, 0, 1)).intersect).toBe(false);
    });

    it('treats tangency as intersection', () => {
        expect(query.test(unitBox, sphere(2, 0, 0, 1)).intersect).toBe(true);
        // The corner (1,1,1) is at distance 5 from (4,5,1) (3-4-5 triangle).
        expect(query.test(unitBox, sphere(4, 5, 1, 5)).intersect).toBe(true);
        expect(query.test(unitBox, sphere(4, 5, 1, 4.999)).intersect).toBe(false);
    });

    it('separates a sphere near a box corner', () => {
        // The corner (1,1,1) is at distance sqrt(3) from (2,2,2).
        expect(query.test(unitBox, sphere(2, 2, 2, 1)).intersect).toBe(false);
        expect(query.test(unitBox, sphere(2, 2, 2, 1.75)).intersect).toBe(true);
    });

    it('separates a sphere near a box edge', () => {
        // The edge x = y = 1 is at distance sqrt(2) from the axis x = y = 2.
        expect(query.test(unitBox, sphere(2, 2, 0, 1.4)).intersect).toBe(false);
        expect(query.test(unitBox, sphere(2, 2, 0, 1.5)).intersect).toBe(true);
    });

    it('handles boxes away from the origin', () => {
        const b = box(v3(10, 20, 30), v3(12, 22, 32));
        expect(query.test(b, sphere(11, 21, 33.5, 1)).intersect).toBe(false);
        expect(query.test(b, sphere(11, 21, 32.5, 1)).intersect).toBe(true);
    });

    it('agrees with the distance query (randomized)', () => {
        const rand = makeRandom(24681357);
        let numIntersect = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const b = box(v3(2 * rand() - 2, 2 * rand() - 2, 2 * rand() - 2),
                v3(2 * rand(), 2 * rand(), 2 * rand()));
            const radius = 0.1 + rand();
            const s = sphere(6 * rand() - 3, 6 * rand() - 3, 6 * rand() - 3,
                radius);
            const expected = distanceToBox(b, s.center) <= radius;
            expect(query.test(b, s).intersect).toBe(expected);
            if (expected) {
                ++numIntersect;
            }
        }
        expect(numIntersect).toBeGreaterThan(50);
        expect(numIntersect).toBeLessThan(450);
    });

    it('throws for mismatched dimensions', () => {
        const b2 = box(Vector.fromArray([-1, -1]), Vector.fromArray([1, 1]));
        expect(() => query.test(b2, sphere(0, 0, 0, 1))).toThrow();
    });
});

describe('IntrAlignedBox3Sphere3FI', () => {
    const query = new IntrAlignedBox3Sphere3FI();
    const zero = v3(0, 0, 0);

    it('reports initial overlap for a sphere inside the box', () => {
        const result = query.find(unitBox, zero, sphere(0.25, -0.5, 0.125, 0.1),
            zero);
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([0.25, -0.5, 0.125]);
    });

    it('reports initial overlap for a sphere crossing a face', () => {
        const result = query.find(unitBox, zero, sphere(1.5, 0.25, -0.5, 1),
            zero);
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 0.25, -0.5]);
    });

    it('reports initial overlap for a sphere covering an edge', () => {
        const result = query.find(unitBox, zero, sphere(1.5, 1.5, 0.25, 1),
            zero);
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1, 0.25]);
    });

    it('reports initial overlap for a sphere covering a vertex', () => {
        const result = query.find(unitBox, zero, sphere(1.5, 1.5, 1.5, 1.5),
            zero);
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1, 1]);
    });

    it('reports contact for exact vertex tangency', () => {
        // The 3-4-5 triangle makes the tangency exact in binary floating
        // point: the corner (1,1,1) is at distance 5 from (4,5,1). Upstream
        // reports the equality case as contact (+1) at time 0.
        const result = query.find(unitBox, zero, sphere(4, 5, 1, 5),
            v3(-1, -1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 1, 1]);
    });

    it('reports no contact for a separated, non-moving sphere', () => {
        const result = query.find(unitBox, zero, sphere(5, 0, 0, 1), zero);
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('reports no contact for a receding sphere', () => {
        const result = query.find(unitBox, zero, sphere(5, 0, 0, 1),
            v3(1, 0, 0));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('reports no contact for a sphere passing by the box', () => {
        const result = query.find(unitBox, zero, sphere(5, -10, 0, 1),
            v3(0, 1, 0));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('computes the contact for approach along the +x face', () => {
        const result = query.find(unitBox, zero, sphere(5, 0.25, -0.5, 1),
            v3(-1, 0, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.25, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(-0.5, 12);
    });

    it('computes the contact for approach along the -z face', () => {
        const result = query.find(unitBox, zero, sphere(0.5, -0.25, -4, 0.5),
            v3(0, 0, 2));
        expect(result.intersectionType).toBe(Type.contact);
        // The sphere boundary reaches z = -1 after moving 2.5 units, so
        // t = 2.5/2 = 1.25.
        expect(result.contactTime).toBeCloseTo(1.25, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(0.5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(-0.25, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(-1, 12);
    });

    it('computes the contact for approach toward an edge', () => {
        // The center travels along the diagonal of the xy-plane toward the
        // box edge x = y = 1. It starts sqrt(2)*3 units from that edge and
        // moves at unit speed, so contact is at t = 3*sqrt(2) - 1 = 3.2426.
        const s = Math.SQRT1_2;
        const result = query.find(unitBox, zero, sphere(4, 4, 0.25, 1),
            v3(-s, -s, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3 * Math.SQRT2 - 1, 10);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 10);
        expect(result.contactPoint.values[1]).toBeCloseTo(1, 10);
        expect(result.contactPoint.values[2]).toBeCloseTo(0.25, 10);
    });

    it('computes the contact for approach toward a vertex', () => {
        // The center starts 4 units (along the diagonal) from the corner
        // (1,1,1) and moves toward it at unit speed, so contact occurs at
        // t = 4 - radius = 3 at the corner.
        const s = 1 / Math.sqrt(3);
        const result = query.find(unitBox, zero,
            sphere(1 + 4 * s, 1 + 4 * s, 1 + 4 * s, 1), v3(-s, -s, -s));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 10);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 10);
        expect(result.contactPoint.values[1]).toBeCloseTo(1, 10);
        expect(result.contactPoint.values[2]).toBeCloseTo(1, 10);
    });

    it('uses the relative velocity of the two objects', () => {
        // Both objects moving identically: no relative motion, no contact.
        const stationary = query.find(unitBox, v3(1, 1, 1),
            sphere(5, 0, 0, 1), v3(1, 1, 1));
        expect(stationary.intersectionType).toBe(Type.noContact);

        // A moving box versus a stationary sphere gives the same contact
        // time as a stationary box versus a sphere with opposite velocity.
        const boxMoves = query.find(unitBox, v3(1, 0, 0),
            sphere(5, 0.25, 0, 1), zero);
        const sphereMoves = query.find(unitBox, zero, sphere(5, 0.25, 0, 1),
            v3(-1, 0, 0));
        expect(boxMoves.intersectionType).toBe(sphereMoves.intersectionType);
        expect(boxMoves.contactTime).toBeCloseTo(sphereMoves.contactTime, 12);
    });

    it('handles boxes that are not centered at the origin', () => {
        const b = box(v3(3, 5, 7), v3(5, 7, 9));  // center (4,6,8)
        const result = query.find(b, zero, sphere(9, 6.25, 8.5, 1),
            v3(-1, 0, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(6.25, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(8.5, 12);
    });

    it('is symmetric under reflection of the configuration', () => {
        const base = query.find(unitBox, zero, sphere(5, 2, 3, 0.75),
            v3(-1, -0.5, -0.75));
        const reflected = query.find(unitBox, zero, sphere(-5, -2, -3, 0.75),
            v3(1, 0.5, 0.75));
        expect(reflected.intersectionType).toBe(base.intersectionType);
        expect(reflected.contactTime).toBeCloseTo(base.contactTime, 12);
        for (let i = 0; i < 3; ++i) {
            expect(reflected.contactPoint.values[i])
                .toBeCloseTo(-base.contactPoint.values[i], 12);
        }
    });

    it('throws for mismatched dimensions', () => {
        const b2 = box(Vector.fromArray([-1, -1]), Vector.fromArray([1, 1]));
        expect(() => query.find(b2, Vector.zero(2),
            Hypersphere.fromCenterRadius(Vector.fromArray([5, 0]), 1),
            Vector.zero(2))).toThrow();
    });

    it('agrees with a sampled simulation (randomized)', () => {
        const rand = makeRandom(11223344);
        let numContact = 0, numOverlap = 0, numNone = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const ex = 0.25 + rand(), ey = 0.25 + rand(), ez = 0.25 + rand();
            const bc = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const b = box(sub(bc, v3(ex, ey, ez)), add(bc, v3(ex, ey, ez)));
            const radius = 0.1 + rand();
            const s = sphere(8 * rand() - 4, 8 * rand() - 4, 8 * rand() - 4,
                radius);
            const bv = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const sv = v3(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            const result = query.find(b, bv, s, sv);
            const rsqr = radius * radius;

            if (result.intersectionType === Type.initiallyOverlapping) {
                ++numOverlap;
                expect(sqrDistanceAtTime(b, bv, s, sv, 0))
                    .toBeLessThanOrEqual(rsqr + 1e-12);
                expect(result.contactTime).toBe(0);
            }
            else if (result.intersectionType === Type.contact) {
                ++numContact;
                const t = result.contactTime;
                expect(t).toBeGreaterThanOrEqual(-1e-15);
                // At the contact time the distance from the moved center to
                // the moved box equals the radius.
                expect(Math.sqrt(sqrDistanceAtTime(b, bv, s, sv, t)))
                    .toBeCloseTo(radius, 8);
                // The contact point is reported in the frame of the box, so
                // it lies on the box at its initial position and is at
                // distance radius from the sphere center after the sphere
                // has moved by the relative velocity for time t.
                expect(distanceToBox(b, result.contactPoint))
                    .toBeLessThan(1e-8);
                const relCenter = add(s.center, mul(sub(sv, bv), t));
                const diff = sub(result.contactPoint, relCenter);
                expect(Math.hypot(diff.values[0], diff.values[1],
                    diff.values[2])).toBeCloseTo(radius, 8);
                // No earlier contact: sample the interval [0,t).
                for (let k = 0; k < 16; ++k) {
                    const u = t * k / 16;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, s, sv, u)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
            else {
                ++numNone;
                // No contact at any sampled future time, out to a horizon
                // beyond which the objects only separate.
                for (let k = 0; k <= 120; ++k) {
                    const u = 20 * k / 120;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, s, sv, u)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
        }
        expect(numContact).toBeGreaterThan(5);
        expect(numOverlap).toBeGreaterThan(5);
        expect(numNone).toBeGreaterThan(5);
    });

    it('reports the contact point in the frame of the box', () => {
        // The box moves and the sphere is stationary. The query works with
        // the relative velocity, so the contact point is on the box at its
        // initial position, not at the position it has at the contact time.
        const s = sphere(5, 0.25, -0.5, 1);
        const result = query.find(unitBox, v3(1, 0, 0), s, zero);
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.25, 12);
        expect(result.contactPoint.values[2]).toBeCloseTo(-0.5, 12);
        expect(distanceToBox(unitBox, result.contactPoint)).toBeLessThan(1e-12);
    });

    it('does not modify the inputs', () => {
        const s = sphere(5, 0.25, 0, 1);
        const sv = v3(-1, 0, 0);
        query.find(unitBox, zero, s, sv);
        expect(s.center.values).toEqual([5, 0.25, 0]);
        expect(sv.values).toEqual([-1, 0, 0]);
        expect(unitBox.min.values).toEqual([-1, -1, -1]);
    });
});

// ---------------------------------------------------------------------------
// Verification (V32): properties cross-checking the port against upstream
// IntrAlignedBox3Sphere3.h.
// ---------------------------------------------------------------------------

describe('IntrAlignedBox3Sphere3 verification', () => {
    const ti = new IntrAlignedBox3Sphere3TI();
    const fi = new IntrAlignedBox3Sphere3FI();

    const arbBox = fc.tuple(wellScaledVector(3, -4, 4),
        fc.array(positive(3, 0.2), { minLength: 3, maxLength: 3 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            Vector.fromArray([c.get(0) - e[0], c.get(1) - e[1],
                c.get(2) - e[2]]),
            Vector.fromArray([c.get(0) + e[0], c.get(1) + e[1],
                c.get(2) + e[2]])));
    const arbSphere = fc.tuple(wellScaledVector(3, -6, 6), positive(2, 0.2))
        .map(([c, r]) => Hypersphere.fromCenterRadius(c, r));

    // The signed gap between the moving sphere and the moving box, from the
    // closed-form point-box distance. Negative when they overlap.
    function gap(b: AlignedBox, bv: Vector, s: Hypersphere, sv: Vector,
        t: number): number {
        let sum = 0;
        for (let i = 0; i < 3; ++i) {
            const lo = b.min.get(i) + t * bv.get(i);
            const hi = b.max.get(i) + t * bv.get(i);
            const p = s.center.get(i) + t * sv.get(i);
            const d = Math.max(lo - p, 0, p - hi);
            sum += d * d;
        }
        return Math.sqrt(sum) - s.radius;
    }

    it('the TI query matches the closed-form point-box distance', () => {
        check(fc.tuple(arbBox, arbSphere), ([b, s]) => {
            const g = gap(b, Vector.zero(3), s, Vector.zero(3), 0);
            // Skip exact tangency, where the closed form and the query's
            // squared comparison can differ by an ulp.
            if (Math.abs(g) > 1e-9 * (1 + s.radius)) {
                expect(ti.test(b, s).intersect).toBe(g <= 0);
            }
        });
    });

    const arbMoving = fc.tuple(arbBox, arbSphere, unitVector(3),
        positive(2, 0.2))
        .map(([b, s, jitter, speed]) => {
            const center = mul(0.5, add(b.min, b.max));
            const aim = add(sub(center, s.center), mul(2.5, jitter));
            const len = Math.sqrt(dot(aim, aim));
            const sv = len > 1e-6 ? mul(speed / len, aim) : mul(speed, jitter);
            return { b, s, sv };
        });

    it('an initial overlap is reported with contact time zero and a point in'
        + ' the box', () => {
        check(arbMoving, ({ b, s, sv }) => {
            if (gap(b, Vector.zero(3), s, sv, 0) >= 0) {
                return;
            }
            const res = fi.find(b, Vector.zero(3), s, sv);
            expect(res.intersectionType).toBe(Type.initiallyOverlapping);
            expect(res.contactTime).toBe(0);
            for (let i = 0; i < 3; ++i) {
                expect(res.contactPoint.get(i))
                    .toBeGreaterThanOrEqual(b.min.get(i) - 1e-9);
                expect(res.contactPoint.get(i))
                    .toBeLessThanOrEqual(b.max.get(i) + 1e-9);
            }
        });
    });

    it('a reported contact time is the first time the objects touch', () => {
        check(arbMoving, ({ b, s, sv }) => {
            const res = fi.find(b, Vector.zero(3), s, sv);
            if (res.intersectionType !== Type.contact
                || res.contactTime === 0) {
                return;
            }
            const t = res.contactTime;
            expect(Number.isFinite(t)).toBe(true);
            expect(t).toBeGreaterThan(0);
            // The objects are never still apart at the reported time: the
            // query does not invent contacts.
            expect(gap(b, Vector.zero(3), s, sv, t)).toBeLessThanOrEqual(1e-7);
            // The reported time is the FIRST contact; that is checked
            // against the Minkowski-sum oracle in the property below, since
            // sampling the gap cannot certify a first crossing. Here the
            // reported point is checked to be a real touch point of the two
            // objects at the reported time.
            // The contact point is on the box and at distance radius from the
            // sphere center at the contact time.
            const p = res.contactPoint;
            const d = sub(p, add(s.center, mul(t, sv)));
            expectClose(Math.sqrt(dot(d, d)), s.radius, 1e-6, 1e-6);
            for (let i = 0; i < 3; ++i) {
                expect(p.get(i)).toBeGreaterThanOrEqual(b.min.get(i) - 1e-7);
                expect(p.get(i)).toBeLessThanOrEqual(b.max.get(i) + 1e-7);
            }
        });
    });

    it('reports the true first contact when the ray leaves the probed'
        + ' rounded edge (#458 item 5)', () => {
        // doQueryRayRoundedFace picks the rounded face the ray enters, and
        // doQueryRayRoundedEdge falls back to the rounded VERTEX of the edge
        // it probes when the cylinder hit is outside the finite cylinder.
        // Upstream stops at that first successful probe, so the other
        // rounded edge of the face - which is where this configuration
        // actually touches first - was never tried and the contact was
        // reported at t = 34.36489152404247 (the box corner), 0.11 late,
        // with the sphere already 0.0139 inside the box.
        const b = box(v3(2.5796760191927124, -0.20000000000000004,
            -0.20000000000000004),
        v3(2.9796760191927127, 0.20000000000000004, 0.20000000000000004));
        const s = sphere(0, -5.576245457464999, -5.370881422178622,
            1.71285599764354);
        const sv = v3(0.05538019511504751, 0.11109696213211791,
            0.15681358038777782);
        const res = fi.find(b, Vector.zero(3), s, sv);
        expect(res.intersectionType).toBe(Type.contact);
        expectClose(res.contactTime, 34.252524073964835, 1e-9, 1e-12);
        expect(res.contactTime).toBeLessThan(34.3);
        // The contact is on the box edge x = b.min[0], y = b.min[1], not at
        // the corner, and the sphere only touches the box there.
        expectVectorClose(res.contactPoint,
            v3(2.5796760191927124, -0.20000000000000004,
                0.00037951517835832504), 1e-9, 1e-9);
        expect(Math.abs(gap(b, Vector.zero(3), s, sv, res.contactTime)))
            .toBeLessThan(1e-9);
        // Strictly separated just before the reported time.
        expect(gap(b, Vector.zero(3), s, sv, res.contactTime - 1e-3))
            .toBeGreaterThan(0);
    });

    it('finds the contact through the rounded edge parallel to the entry'
        + ' face (#465 item 2)', () => {
        // The ray leaves the entry face on both sides, so the first contact
        // is on the rounded edge PARALLEL to the face normal at the overhang
        // corner. That edge bounds neither of the two rounded edges upstream
        // probes, both of whose rounded-vertex fallbacks miss here, so
        // upstream reported noContact although the sphere penetrates the box
        // by 5.7e-3 at the closest approach.
        const c = v3(0.090417948551476, -0.5309329181909561,
            0.9111482501029968);
        const e = v3(0.3892387102358043, 0.3108668552711606,
            0.5541517764329911);
        const b = box(sub(c, e), add(c, e));
        const s = sphere(5.136491211596876, -1.6951057985424995,
            1.702867484651506, 0.2129704826977104);
        const sv = v3(-1, 0.16651661535725, -0.005579612776637072);
        const res = fi.find(b, Vector.zero(3), s, sv);
        expect(res.intersectionType).toBe(Type.contact);
        expectClose(res.contactTime, 4.9164603303241305, 1e-9, 1e-12);
        // The contact point is on the box edge y = b.min[1], z = b.max[2].
        expectVectorClose(res.contactPoint,
            v3(0.22003088127274584, -0.8417997734621168, 1.465300026535988),
            1e-9, 1e-9);
        expect(Math.abs(gap(b, Vector.zero(3), s, sv, res.contactTime)))
            .toBeLessThan(1e-9);
        expect(gap(b, Vector.zero(3), s, sv, res.contactTime - 1e-3))
            .toBeGreaterThan(0);
    });

    it('probes the rounded edges that meet the nearest rounded vertex', () => {
        // The vertexSeparated case of doQuery: the sphere center is in the
        // corner box [K, K + radius] of the first octant, so the nearest
        // piece of the rounded box is the rounded vertex. Upstream probes
        // only that vertex, but the ray can enter through one of the three
        // rounded edges meeting it, either earlier than the vertex or where
        // the vertex is missed altogether. Both configurations were found by
        // the Minkowski-sum oracle below.
        const missed = {
            b: box(v3(-2.8967879977188056, -2.8644249146894176,
                -2.1982821431361246),
            v3(2.8967879977188056, 2.8644249146894176, 2.1982821431361246)),
            s: sphere(7.184056052312964, -4.743390293991202,
                6.546082936666346, 4.44050498629351),
            sv: v3(-4.776835712724244, 5.672260885355733,
                -0.1792759264727523),
            // upstream: noContact
            t: 0.6191840057453503,
            p: v3(2.8967879977188056, -1.2312170773639721, 2.1982821431361246)
        };
        const late = {
            b: box(v3(-0.7669004777317825, -13.33655954226376,
                -5.233129731083094),
            v3(0.7669004777317825, 13.33655954226376, 5.233129731083094)),
            s: sphere(-7.660843902444189, -13.890512876707406,
                -9.942798799926393, 7.294105788085315),
            sv: v3(3.2224692471639997, 10.097123601932271, 5.802238437727227),
            // upstream: contact at t = 0.21189258023326504, 0.031 late
            t: 0.1811483750315202,
            p: v3(-0.7669004777317825, -12.061435343724966,
                -5.233129731083094)
        };
        for (const c of [missed, late]) {
            const res = fi.find(c.b, Vector.zero(3), c.s, c.sv);
            expect(res.intersectionType).toBe(Type.contact);
            expectClose(res.contactTime, c.t, 1e-9, 1e-12);
            expectVectorClose(res.contactPoint, c.p, 1e-9, 1e-9);
            expect(Math.abs(gap(c.b, Vector.zero(3), c.s, c.sv,
                res.contactTime))).toBeLessThan(1e-9);
            expect(gap(c.b, Vector.zero(3), c.s, c.sv,
                res.contactTime * (1 - 1e-4))).toBeGreaterThan(0);
        }
    });

    it('a reported no-contact is confirmed by sampling the motion', () => {
        check(arbMoving, ({ b, s, sv }) => {
            const res = fi.find(b, Vector.zero(3), s, sv);
            if (res.intersectionType !== Type.noContact) {
                return;
            }
            for (let k = 0; k <= 400; ++k) {
                expect(gap(b, Vector.zero(3), s, sv, (20 * k) / 400))
                    .toBeGreaterThan(-1e-7);
            }
        }, 80);
    });

    it('a no-contact result reports either the zero vector or the box center'
        + ' as the contact point (upstream quirk, preserved)', () => {
        // Upstream translates the contact point back by the box center
        // whenever the ray-versus-superbox prefilter succeeds, even when the
        // detailed query reported no contact. The Result documentation claims
        // contactPoint = (0,0,0) in that case.
        check(arbMoving, ({ b, s, sv }) => {
            const res = fi.find(b, Vector.zero(3), s, sv);
            if (res.intersectionType !== Type.noContact) {
                return;
            }
            const center = mul(0.5, add(b.min, b.max));
            const isZero = res.contactPoint.equals(Vector.zero(3));
            const isCenter = res.contactPoint.equals(center);
            expect(isZero || isCenter).toBe(true);
        });
    });

    it('the result depends only on the relative velocity', () => {
        check(fc.tuple(arbMoving, wellScaledVector(3, -2, 2)),
            ([{ b, s, sv }, drift]) => {
                const a = fi.find(b, Vector.zero(3), s, sv);
                const d = fi.find(b, drift, s, add(sv, drift));
                expect(d.intersectionType).toBe(a.intersectionType);
                expectClose(d.contactTime, a.contactTime, 1e-9, 1e-9);
            });
    });

    it('is equivariant under translation and under reflection in the axes',
        () => {
            check(fc.tuple(arbMoving, wellScaledVector(3, -5, 5)),
                ([{ b, s, sv }, shift]) => {
                    const a = fi.find(b, Vector.zero(3), s, sv);

                    const bT = AlignedBox.fromMinMax(add(b.min, shift),
                        add(b.max, shift));
                    const sT = Hypersphere.fromCenterRadius(
                        add(s.center, shift), s.radius);
                    const t = fi.find(bT, Vector.zero(3), sT, sv);
                    expect(t.intersectionType).toBe(a.intersectionType);
                    expectClose(t.contactTime, a.contactTime, 1e-9, 1e-9);
                    if (a.intersectionType === Type.contact
                        || a.intersectionType === Type.initiallyOverlapping) {
                        expectVectorClose(t.contactPoint,
                            add(a.contactPoint, shift), 1e-8, 1e-8);
                    }

                    // Reflect z -> -z; the query mirrors the sphere center
                    // into the first octant, so this exercises the sign
                    // bookkeeping.
                    const flip = (v: Vector): Vector => Vector.fromArray([
                        v.get(0), v.get(1), -v.get(2)]);
                    const bF = AlignedBox.fromMinMax(
                        Vector.fromArray([b.min.get(0), b.min.get(1),
                            -b.max.get(2)]),
                        Vector.fromArray([b.max.get(0), b.max.get(1),
                            -b.min.get(2)]));
                    const sF = Hypersphere.fromCenterRadius(flip(s.center),
                        s.radius);
                    const f = fi.find(bF, Vector.zero(3), sF, flip(sv));
                    expect(f.intersectionType).toBe(a.intersectionType);
                    expectClose(f.contactTime, a.contactTime, 1e-9, 1e-9);
                    if (a.intersectionType === Type.contact
                        || a.intersectionType === Type.initiallyOverlapping) {
                        expectVectorClose(f.contactPoint, flip(a.contactPoint),
                            1e-8, 1e-8);
                    }
                });
        });

    it('reports no NaN in any field', () => {
        check(arbMoving, ({ b, s, sv }) => {
            const res = fi.find(b, Vector.zero(3), s, sv);
            expect(Number.isNaN(res.contactTime)).toBe(false);
            for (let i = 0; i < 3; ++i) {
                expect(Number.isNaN(res.contactPoint.get(i))).toBe(false);
            }
        });
    });

    // -----------------------------------------------------------------
    // The first contact computed independently of the query, as the first
    // intersection of the ray C + t*V (C the sphere center relative to the
    // box center, V the relative velocity) with the Minkowski sum of the box
    // and the ball of radius r. That boundary is made of the six box faces
    // pushed out by r, the twelve cylinders of radius r around the box edges
    // and the eight spheres of radius r at the box vertices. Every piece is
    // a subset of the Minkowski sum, so a ray coming from outside meets each
    // piece no earlier than it enters the sum, and it enters the sum on one
    // of them: the smallest nonnegative parameter over all pieces is the
    // first contact time.
    // -----------------------------------------------------------------
    function minkowskiFirstContact(b: AlignedBox, bv: Vector, s: Hypersphere,
        sv: Vector): { type: Type, t: number, point: Vector } {
        const center = mul(0.5, add(b.min, b.max));
        const K = mul(0.5, sub(b.max, b.min)).values;
        const C = sub(s.center, center).values;
        const V = sub(sv, bv).values;
        const r = s.radius;

        // Already overlapping?
        let sqrLen = 0;
        for (let i = 0; i < 3; ++i) {
            const d = Math.max(Math.abs(C[i]) - K[i], 0);
            sqrLen += d * d;
        }
        if (Math.sqrt(sqrLen) <= r) {
            return {
                type: Type.initiallyOverlapping, t: 0,
                point: add(Vector.fromArray(C), center)
            };
        }

        let best = Number.POSITIVE_INFINITY;
        let bestPoint = [0, 0, 0];
        const at = (t: number): number[] =>
            [C[0] + t * V[0], C[1] + t * V[1], C[2] + t * V[2]];
        const consider = (t: number, q: number[]): void => {
            if (t >= 0 && t < best) { best = t; bestPoint = q; }
        };

        // The six faces, offset by r along their normals.
        for (let i = 0; i < 3; ++i) {
            const i1 = (i + 1) % 3;
            const i2 = (i + 2) % 3;
            for (const sgn of [-1, 1]) {
                if (V[i] === 0) { continue; }
                const t = (sgn * (K[i] + r) - C[i]) / V[i];
                if (t < 0) { continue; }
                const p = at(t);
                if (Math.abs(p[i1]) <= K[i1] && Math.abs(p[i2]) <= K[i2]) {
                    const q = p.slice();
                    q[i] = sgn * K[i];
                    consider(t, q);
                }
            }
        }

        // The twelve cylinders of radius r around the box edges; the edge
        // parallel to axis k passes through (si*K[i], sj*K[j]).
        for (let k = 0; k < 3; ++k) {
            const i = (k + 1) % 3;
            const j = (k + 2) % 3;
            for (const si of [-1, 1]) {
                for (const sj of [-1, 1]) {
                    const ei = C[i] - si * K[i];
                    const ej = C[j] - sj * K[j];
                    const c2 = V[i] * V[i] + V[j] * V[j];
                    if (c2 === 0) { continue; }
                    const c1 = V[i] * ei + V[j] * ej;
                    const c0 = ei * ei + ej * ej - r * r;
                    const discr = c1 * c1 - c2 * c0;
                    if (discr < 0) { continue; }
                    const sq = Math.sqrt(discr);
                    for (const t of [(-c1 - sq) / c2, (-c1 + sq) / c2]) {
                        if (t < 0) { continue; }
                        const p = at(t);
                        if (Math.abs(p[k]) <= K[k]) {
                            const q = [0, 0, 0];
                            q[i] = si * K[i];
                            q[j] = sj * K[j];
                            q[k] = p[k];
                            consider(t, q);
                            break;
                        }
                    }
                }
            }
        }

        // The eight spheres of radius r at the box vertices.
        const a2 = V[0] * V[0] + V[1] * V[1] + V[2] * V[2];
        if (a2 > 0) {
            for (const s0 of [-1, 1]) {
                for (const s1 of [-1, 1]) {
                    for (const s2 of [-1, 1]) {
                        const q = [s0 * K[0], s1 * K[1], s2 * K[2]];
                        const e = [C[0] - q[0], C[1] - q[1], C[2] - q[2]];
                        const a1 = V[0] * e[0] + V[1] * e[1] + V[2] * e[2];
                        const a0 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2]
                            - r * r;
                        const discr = a1 * a1 - a2 * a0;
                        if (discr < 0) { continue; }
                        consider((-a1 - Math.sqrt(discr)) / a2, q);
                    }
                }
            }
        }

        if (!Number.isFinite(best)) {
            return { type: Type.noContact, t: 0, point: Vector.zero(3) };
        }
        return {
            type: Type.contact, t: best,
            point: add(Vector.fromArray(bestPoint), center)
        };
    }

    // The signed gap is convex in t (the distance from a point to a convex
    // set along a line is convex), so a ternary search finds its minimum.
    // A configuration whose minimum gap is within a small margin of zero is
    // tangential: the query and the oracle round the discriminants of the
    // same quadratics differently there and can disagree on whether the
    // objects touch at all, so such configurations are skipped.
    function minimumGap(b: AlignedBox, s: Hypersphere, sv: Vector,
        tHi: number): number {
        let lo = 0;
        let hi = tHi;
        for (let i = 0; i < 100; ++i) {
            const m0 = lo + (hi - lo) / 3;
            const m1 = hi - (hi - lo) / 3;
            if (gap(b, Vector.zero(3), s, sv, m0)
                <= gap(b, Vector.zero(3), s, sv, m1)) {
                hi = m1;
            }
            else {
                lo = m0;
            }
        }
        return gap(b, Vector.zero(3), s, sv, 0.5 * (lo + hi));
    }

    // Compare one configuration against the oracle. Returns false when the
    // configuration is skipped as tangential.
    function checkAgainstOracle(b: AlignedBox, s: Hypersphere,
        sv: Vector): boolean {
        const speed = Math.sqrt(dot(sv, sv));
        if (speed === 0) { return false; }
        const scale = 1 + s.radius + Math.sqrt(dot(b.max, b.max))
            + Math.sqrt(dot(s.center, s.center));
        const oracle = minkowskiFirstContact(b, Vector.zero(3), s, sv);
        const tHi = oracle.type === Type.contact
            ? 2 * oracle.t + scale / speed
            : 100 * scale / speed;
        if (Math.abs(minimumGap(b, s, sv, tHi)) < 1e-7 * scale) {
            return false;
        }
        const res = fi.find(b, Vector.zero(3), s, sv);
        expect(res.intersectionType).toBe(oracle.type);
        if (oracle.type === Type.contact) {
            // The query and the oracle solve the same quadratics, so they
            // agree to the rounding of their coefficients; a time tolerance
            // is a length tolerance divided by the speed.
            expectClose(res.contactTime, oracle.t, 1e-9 * scale / speed,
                1e-9);
        }
        return true;
    }

    it('agrees with the Minkowski-sum first-contact oracle', () => {
        let tested = 0;
        check(arbMoving, ({ b, s, sv }) => {
            if (gap(b, Vector.zero(3), s, sv, 0) <= 0) {
                return;
            }
            if (checkAgainstOracle(b, s, sv)) { ++tested; }
        });
        expect(tested).toBeGreaterThan(50);
    }, 30000);

    it('agrees with the Minkowski-sum oracle on corner-aimed motions', () => {
        // The face and vertex case analyses are exercised by motions aimed
        // at the box corners, where the ray leaves the rounded face it
        // enters; that is the family both #458 item 5 and #465 item 2 come
        // from. The generator spans three decades of scale.
        const rnd = seededRandom(0x5b03e1);
        let contacts = 0;
        let tested = 0;
        for (let iter = 0; iter < 4000; ++iter) {
            const scale = Math.exp(6 * rnd() - 3);
            const K = [0, 1, 2].map(() => scale * (0.05 + 2 * rnd()));
            const radius = scale * (0.02 + 1.5 * rnd());
            const b = box(v3(-K[0], -K[1], -K[2]), v3(K[0], K[1], K[2]));
            let u = [2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1];
            const un = Math.hypot(u[0], u[1], u[2]);
            if (un < 1e-3) { continue; }
            u = u.map(x => x / un);
            const dist = scale * (2 + 10 * rnd());
            const c = u.map(x => x * dist);
            const s = sphere(c[0], c[1], c[2], radius);
            // Aim at a box corner, with a jitter of the order of the corner
            // region so that the ray often just clips or just misses it.
            const target = [0, 1, 2].map(i => K[i] * (rnd() < 0.5 ? -1 : 1));
            const d = [0, 1, 2].map(i => target[i] - c[i]
                + (radius + K[i]) * 0.25 * (2 * rnd() - 1));
            const dn = Math.hypot(d[0], d[1], d[2]);
            if (dn < 1e-6) { continue; }
            const speed = scale * (0.1 + 3 * rnd());
            const sv = v3(d[0] * speed / dn, d[1] * speed / dn,
                d[2] * speed / dn);
            if (gap(b, Vector.zero(3), s, sv, 0) <= 0) { continue; }
            if (checkAgainstOracle(b, s, sv)) {
                ++tested;
                if (fi.find(b, Vector.zero(3), s, sv).intersectionType
                    === Type.contact) {
                    ++contacts;
                }
            }
        }
        expect(tested).toBeGreaterThan(3000);
        expect(contacts).toBeGreaterThan(500);
    }, 30000);

    it('rejects inputs of the wrong dimension', () => {
        const b2 = AlignedBox.fromMinMax(Vector.zero(2), Vector.zero(2));
        const s3 = sphere(0, 0, 0, 1);
        expect(() => ti.test(b2, s3)).toThrow('mismatched sizes');
        expect(() => fi.find(b2, Vector.zero(2), s3, Vector.zero(3)))
            .toThrow('mismatched sizes');
    });
});
