import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox';
import { Hypersphere } from '../src/Hypersphere';
import { DistPointAlignedBox } from '../src/DistPointAlignedBox';
import { Vector, add, mul } from '../src/Vector';
import {
    IntrAlignedBox2Circle2TI,
    IntrAlignedBox2Circle2FI,
    IntrAlignedBox2Circle2FIResultType as Type
} from '../src/IntrAlignedBox2Circle2';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function box(minx: number, miny: number, maxx: number, maxy: number): AlignedBox {
    return AlignedBox.fromMinMax(v2(minx, miny), v2(maxx, maxy));
}

function circle(cx: number, cy: number, r: number): Hypersphere {
    return Hypersphere.fromCenterRadius(v2(cx, cy), r);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// The distance from the circle center to the box at time t, where the box
// moves with velocity bv and the circle center with velocity cv.
function sqrDistanceAtTime(b: AlignedBox, bv: Vector, c: Hypersphere,
    cv: Vector, t: number): number {
    const movedBox = AlignedBox.fromMinMax(add(b.min, mul(bv, t)),
        add(b.max, mul(bv, t)));
    const movedCenter = add(c.center, mul(cv, t));
    return new DistPointAlignedBox().compute(movedCenter, movedBox).sqrDistance;
}

const unitBox = box(-1, -1, 1, 1);

describe('IntrAlignedBox2Circle2TI', () => {
    const query = new IntrAlignedBox2Circle2TI();

    it('detects a circle inside the box', () => {
        expect(query.test(unitBox, circle(0, 0, 0.25)).intersect).toBe(true);
    });

    it('detects a box inside the circle', () => {
        expect(query.test(unitBox, circle(0, 0, 10)).intersect).toBe(true);
    });

    it('detects face overlap and separation', () => {
        expect(query.test(unitBox, circle(1.5, 0, 1)).intersect).toBe(true);
        expect(query.test(unitBox, circle(2.5, 0, 1)).intersect).toBe(false);
    });

    it('treats tangency as intersection', () => {
        expect(query.test(unitBox, circle(2, 0, 1)).intersect).toBe(true);
        // Corner tangency: distance from (1,1) to the center is exactly 1.
        const s = Math.SQRT1_2;
        expect(query.test(unitBox, circle(1 + s, 1 + s, 1)).intersect).toBe(true);
    });

    it('separates a circle near a box corner', () => {
        // The corner (1,1) is at distance sqrt(2) from (2,2), which exceeds 1.
        expect(query.test(unitBox, circle(2, 2, 1)).intersect).toBe(false);
        expect(query.test(unitBox, circle(2, 2, 1.5)).intersect).toBe(true);
    });
});

describe('IntrAlignedBox2Circle2FI', () => {
    const query = new IntrAlignedBox2Circle2FI();

    it('reports initial overlap for a circle inside the box', () => {
        const result = query.find(unitBox, v2(0, 0), circle(0.25, -0.5, 0.1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([0.25, -0.5]);
    });

    it('reports initial overlap for a circle crossing a face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(1.5, 0.25, 1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 0.25]);
    });

    it('reports initial overlap for a circle covering a vertex', () => {
        const result = query.find(unitBox, v2(0, 0), circle(1.5, 1.5, 1),
            v2(0, 0));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('reports no contact for a separated, non-moving circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0, 1), v2(0, 0));
        expect(result.intersectionType).toBe(Type.noContact);
        expect(result.contactTime).toBe(0);
    });

    it('reports contact at t = 0 for tangency along a face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(2, 0.5, 1),
            v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 0.5]);
    });

    it('reports contact at t = 0 for tangency at a vertex', () => {
        // The 3-4-5 triangle makes the tangency exact in binary floating
        // point: the corner (1,1) is at squared distance 3^2 + 4^2 = 25 from
        // the center (4,5), which equals radius^2 for radius 5. Upstream
        // distinguishes overlap (sqrDistance < sqrRadius) from contact
        // (equality), so this is the +1 case with contactTime 0.
        const result = query.find(unitBox, v2(0, 0), circle(4, 5, 5),
            v2(-1, -1));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBe(0);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('reports overlap when a vertex is strictly inside the circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(4, 5, 5.5),
            v2(-1, -1));
        expect(result.intersectionType).toBe(Type.initiallyOverlapping);
        expect(result.contactPoint.values).toEqual([1, 1]);
    });

    it('computes the contact for approach along the +x face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0.25, 1),
            v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.25, 12);
    });

    it('computes the contact for approach along the -y face', () => {
        const result = query.find(unitBox, v2(0, 0), circle(-0.5, -4, 0.5),
            v2(0, 2));
        expect(result.intersectionType).toBe(Type.contact);
        // The circle boundary reaches y = -1 after moving 2.5 units, so
        // t = 2.5/2 = 1.25.
        expect(result.contactTime).toBeCloseTo(1.25, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(-0.5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(-1, 12);
    });

    it('computes the contact for approach toward a vertex', () => {
        // The center starts 3 units (along the diagonal) from the corner
        // (1,1) and moves toward it at unit speed, so contact occurs at
        // t = 3 - radius = 2 at the corner.
        const s = Math.SQRT1_2;
        const result = query.find(unitBox, v2(0, 0),
            circle(1 + 3 * s, 1 + 3 * s, 1), v2(-s, -s));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(2, 10);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(1, 12);
    });

    it('computes the contact for a vertex approach that grazes a face', () => {
        // The circle is diagonally away from the corner but moves in -x
        // only, so it eventually strikes the +x face (or the corner region).
        const result = query.find(unitBox, v2(0, 0), circle(6, 0.5, 1),
            v2(-2, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(2, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(1, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(0.5, 12);
    });

    it('reports no contact for a receding circle', () => {
        const result = query.find(unitBox, v2(0, 0), circle(5, 0, 1), v2(1, 0));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('reports no contact for a circle passing by the box', () => {
        // The circle travels along x = 5, never nearing the box.
        const result = query.find(unitBox, v2(0, 0), circle(5, -10, 1),
            v2(0, 1));
        expect(result.intersectionType).toBe(Type.noContact);
    });

    it('uses the relative velocity of the two objects', () => {
        // The box chases the circle: with both moving at the same velocity
        // the relative velocity is zero and there is no contact.
        const stationary = query.find(unitBox, v2(1, 1), circle(5, 0, 1),
            v2(1, 1));
        expect(stationary.intersectionType).toBe(Type.noContact);

        // The box moves toward a stationary circle: same answer as a circle
        // moving toward a stationary box with the opposite velocity.
        const boxMoves = query.find(unitBox, v2(1, 0), circle(5, 0.25, 1),
            v2(0, 0));
        const circleMoves = query.find(unitBox, v2(0, 0), circle(5, 0.25, 1),
            v2(-1, 0));
        expect(boxMoves.intersectionType).toBe(circleMoves.intersectionType);
        expect(boxMoves.contactTime).toBeCloseTo(circleMoves.contactTime, 12);
    });

    it('handles boxes that are not centered at the origin', () => {
        const b = box(3, 5, 5, 7);  // center (4,6), extent (1,1)
        const result = query.find(b, v2(0, 0), circle(9, 6.25, 1), v2(-1, 0));
        expect(result.intersectionType).toBe(Type.contact);
        expect(result.contactTime).toBeCloseTo(3, 12);
        expect(result.contactPoint.values[0]).toBeCloseTo(5, 12);
        expect(result.contactPoint.values[1]).toBeCloseTo(6.25, 12);
    });

    it('is symmetric under reflection of the configuration', () => {
        const base = query.find(unitBox, v2(0, 0), circle(5, 2, 0.75),
            v2(-1, -0.5));
        const reflected = query.find(unitBox, v2(0, 0), circle(-5, -2, 0.75),
            v2(1, 0.5));
        expect(reflected.intersectionType).toBe(base.intersectionType);
        expect(reflected.contactTime).toBeCloseTo(base.contactTime, 12);
        expect(reflected.contactPoint.values[0])
            .toBeCloseTo(-base.contactPoint.values[0], 12);
        expect(reflected.contactPoint.values[1])
            .toBeCloseTo(-base.contactPoint.values[1], 12);
    });

    it('agrees with a sampled simulation (randomized)', () => {
        const rand = makeRandom(13572468);
        let numContact = 0, numOverlap = 0, numNone = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const ex = 0.25 + rand(), ey = 0.25 + rand();
            const bcx = 2 * rand() - 1, bcy = 2 * rand() - 1;
            const b = box(bcx - ex, bcy - ey, bcx + ex, bcy + ey);
            const radius = 0.1 + rand();
            const c = circle(8 * rand() - 4, 8 * rand() - 4, radius);
            const bv = v2(2 * rand() - 1, 2 * rand() - 1);
            const cv = v2(2 * rand() - 1, 2 * rand() - 1);
            const result = query.find(b, bv, c, cv);
            const rsqr = radius * radius;

            if (result.intersectionType === Type.initiallyOverlapping) {
                ++numOverlap;
                expect(sqrDistanceAtTime(b, bv, c, cv, 0))
                    .toBeLessThanOrEqual(rsqr + 1e-12);
                expect(result.contactTime).toBe(0);
            }
            else if (result.intersectionType === Type.contact) {
                ++numContact;
                const t = result.contactTime;
                expect(t).toBeGreaterThanOrEqual(0);
                // At the contact time the distance from the moved center to
                // the moved box equals the radius.
                expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, t)))
                    .toBeCloseTo(radius, 8);
                // No earlier contact: sample the interval [0,t).
                for (let k = 0; k < 20; ++k) {
                    const s = t * k / 20;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, s)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
            else {
                ++numNone;
                // No contact at any sampled future time, out to a horizon
                // beyond which the objects only separate.
                for (let k = 0; k <= 200; ++k) {
                    const s = 20 * k / 200;
                    expect(Math.sqrt(sqrDistanceAtTime(b, bv, c, cv, s)))
                        .toBeGreaterThan(radius - 1e-9);
                }
            }
        }
        expect(numContact).toBeGreaterThan(10);
        expect(numOverlap).toBeGreaterThan(10);
        expect(numNone).toBeGreaterThan(10);
    });
});
