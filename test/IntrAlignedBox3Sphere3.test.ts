import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Hypersphere } from '../src/Hypersphere';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox';
import { Vector, add, mul, sub } from '../src/Vector';
import {
    IntrAlignedBox3Sphere3TI,
    IntrAlignedBox3Sphere3FI,
    IntrAlignedBox3Sphere3FIResultType as Type
} from '../src/IntrAlignedBox3Sphere3';

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
