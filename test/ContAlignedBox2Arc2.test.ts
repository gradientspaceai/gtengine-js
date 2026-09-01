import { describe, it, expect } from 'vitest';
import { getContainerAlignedBox2Arc2 } from '../src/ContAlignedBox2Arc2';
import { Arc2 } from '../src/Arc2';
import { Vector } from '../src/Vector';

function v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Build the arc on the unit-ish circle (center C, radius r) with endpoints at
// the given angles; counterclockwise from angle0 to angle1.
function arcFromAngles(cx: number, cy: number, r: number, angle0: number,
    angle1: number): Arc2 {
    const e0 = v(cx + r * Math.cos(angle0), cy + r * Math.sin(angle0));
    const e1 = v(cx + r * Math.cos(angle1), cy + r * Math.sin(angle1));
    return Arc2.fromCenterRadiusEnds(v(cx, cy), r, e0, e1);
}

const eps = 1e-12;

describe('getContainerAlignedBox2Arc2', () => {
    it('bounds a full-circle-like arc by the circle bounding box', () => {
        // Endpoints coincide (E0 == E1), so every axis-extreme point passes
        // the side-of-line test (dotPerp(P-E0, E1-E0) = 0 >= 0) and the box
        // is the bounding box of the whole circle.
        const c = v(2, -3);
        const r = 5;
        const arc = Arc2.fromCenterRadiusEnds(c, r, v(7, -3), v(7, -3));
        const box = getContainerAlignedBox2Arc2(arc);
        expect(box.min.values[0]).toBeCloseTo(-3, 12);
        expect(box.min.values[1]).toBeCloseTo(-8, 12);
        expect(box.max.values[0]).toBeCloseTo(7, 12);
        expect(box.max.values[1]).toBeCloseTo(2, 12);
    });

    it('bounds a quarter arc in the first quadrant by its endpoints', () => {
        // From (1,0) counterclockwise to (0,1). Only C+(r,0) = (1,0) and
        // C+(0,r) = (0,1) are on the arc, and both are endpoints, so the box
        // is [0,1]x[0,1].
        const arc = arcFromAngles(0, 0, 1, 0, Math.PI / 2);
        const box = getContainerAlignedBox2Arc2(arc);
        expect(box.min.values[0]).toBeCloseTo(0, 12);
        expect(box.min.values[1]).toBeCloseTo(0, 12);
        expect(box.max.values[0]).toBeCloseTo(1, 12);
        expect(box.max.values[1]).toBeCloseTo(1, 12);
    });

    it('includes an axis-extreme point interior to the arc', () => {
        // From angle -pi/4 to +pi/4 on the unit circle. The point (1,0) is
        // interior to the arc, so max x is 1 rather than cos(pi/4).
        const arc = arcFromAngles(0, 0, 1, -Math.PI / 4, Math.PI / 4);
        const box = getContainerAlignedBox2Arc2(arc);
        const s = Math.SQRT1_2;
        expect(box.max.values[0]).toBeCloseTo(1, 12);
        expect(box.min.values[0]).toBeCloseTo(s, 12);
        expect(box.min.values[1]).toBeCloseTo(-s, 12);
        expect(box.max.values[1]).toBeCloseTo(s, 12);
    });

    it('bounds a three-quarter arc crossing three axis extremes', () => {
        // From angle 0 counterclockwise to 3*pi/2 (i.e. (1,0) to (0,-1)).
        // The arc contains (0,1), (-1,0) and (0,-1); it also contains its
        // endpoints. All four extremes are attained.
        const arc = arcFromAngles(0, 0, 1, 0, 3 * Math.PI / 2);
        const box = getContainerAlignedBox2Arc2(arc);
        expect(box.min.values[0]).toBeCloseTo(-1, 12);
        expect(box.min.values[1]).toBeCloseTo(-1, 12);
        expect(box.max.values[0]).toBeCloseTo(1, 12);
        expect(box.max.values[1]).toBeCloseTo(1, 12);
    });

    it('is translation covariant', () => {
        const arc0 = arcFromAngles(0, 0, 2, 0.3, 2.1);
        const arc1 = arcFromAngles(10, -7, 2, 0.3, 2.1);
        const box0 = getContainerAlignedBox2Arc2(arc0);
        const box1 = getContainerAlignedBox2Arc2(arc1);
        expect(box1.min.values[0]).toBeCloseTo(box0.min.values[0] + 10, 12);
        expect(box1.min.values[1]).toBeCloseTo(box0.min.values[1] - 7, 12);
        expect(box1.max.values[0]).toBeCloseTo(box0.max.values[0] + 10, 12);
        expect(box1.max.values[1]).toBeCloseTo(box0.max.values[1] - 7, 12);
    });

    it('contains sampled points of the arc (randomized)', () => {
        let seed = 12345;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 40; ++trial) {
            const cx = 10 * (rand() - 0.5);
            const cy = 10 * (rand() - 0.5);
            const r = 0.5 + 3 * rand();
            const a0 = 2 * Math.PI * rand();
            const sweep = 0.05 + (2 * Math.PI - 0.1) * rand();
            const a1 = a0 + sweep;
            const arc = arcFromAngles(cx, cy, r, a0, a1);
            const box = getContainerAlignedBox2Arc2(arc);

            for (let k = 0; k <= 64; ++k) {
                const a = a0 + (sweep * k) / 64;
                const x = cx + r * Math.cos(a);
                const y = cy + r * Math.sin(a);
                expect(x).toBeGreaterThanOrEqual(box.min.values[0] - 1e-9);
                expect(x).toBeLessThanOrEqual(box.max.values[0] + 1e-9);
                expect(y).toBeGreaterThanOrEqual(box.min.values[1] - 1e-9);
                expect(y).toBeLessThanOrEqual(box.max.values[1] + 1e-9);
            }
        }
    });

    it('handles a degenerate zero-radius arc', () => {
        const c = v(4, 5);
        const arc = Arc2.fromCenterRadiusEnds(c, 0, c, c);
        const box = getContainerAlignedBox2Arc2(arc);
        expect(box.min.values).toEqual([4, 5]);
        expect(box.max.values).toEqual([4, 5]);
    });

    it('returns a 2D box', () => {
        const box = getContainerAlignedBox2Arc2(new Arc2());
        expect(box.dimension).toBe(2);
        // The default arc has center (0,0), radius 1, ends (1,0) and (0,1),
        // which is the first-quadrant quarter arc.
        expect(Math.abs(box.min.values[0])).toBeLessThan(eps);
        expect(Math.abs(box.min.values[1])).toBeLessThan(eps);
        expect(box.max.values[0]).toBeCloseTo(1, 12);
        expect(box.max.values[1]).toBeCloseTo(1, 12);
    });
});
