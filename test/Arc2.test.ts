import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2';
import { Vector } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

describe('Arc2 construction', () => {
    it('the default constructor is the quarter arc of the unit circle', () => {
        const arc = new Arc2();
        expect(arc.center.values).toEqual([0, 0]);
        expect(arc.radius).toBe(1);
        expect(arc.end[0].values).toEqual([1, 0]);
        expect(arc.end[1].values).toEqual([0, 1]);
    });

    it('fromCenterRadiusEnds copies the inputs', () => {
        const c = v2(1, 2);
        const e0 = v2(3, 2);
        const e1 = v2(1, 4);
        const arc = Arc2.fromCenterRadiusEnds(c, 2, e0, e1);
        c.set(0, 99);
        e0.set(0, 99);
        e1.set(1, 99);
        expect(arc.center.values).toEqual([1, 2]);
        expect(arc.radius).toBe(2);
        expect(arc.end[0].values).toEqual([3, 2]);
        expect(arc.end[1].values).toEqual([1, 4]);
    });

    it('rejects vectors that are not 2D', () => {
        expect(() => Arc2.fromCenterRadiusEnds(Vector.fromArray([0, 0, 0]), 1,
            v2(1, 0), v2(0, 1))).toThrow();
    });

    it('clone is a deep copy', () => {
        const arc = new Arc2();
        const copy = arc.clone();
        copy.center.set(0, 5);
        copy.end[0].set(1, 7);
        expect(arc.center.values).toEqual([0, 0]);
        expect(arc.end[0].values).toEqual([1, 0]);
        expect(copy.center.values).toEqual([5, 0]);
    });
});

describe('Arc2 containment', () => {
    // The default arc runs counterclockwise from (1,0) to (0,1), so it is the
    // first-quadrant quarter of the unit circle.
    const arc = new Arc2();

    it('accepts points on the arc, exactly', () => {
        const s = Math.SQRT1_2;
        expect(arc.contains(v2(1, 0), 0)).toBe(true);
        expect(arc.contains(v2(0, 1), 0)).toBe(true);
        expect(arc.containsOnCircle(v2(s, s))).toBe(true);
    });

    it('rejects points on the circle but off the arc', () => {
        const s = Math.SQRT1_2;
        expect(arc.containsOnCircle(v2(-1, 0))).toBe(false);
        expect(arc.containsOnCircle(v2(0, -1))).toBe(false);
        expect(arc.containsOnCircle(v2(-s, -s))).toBe(false);
        expect(arc.contains(v2(-1, 0), 1e-12)).toBe(false);
    });

    it('rejects points off the circle unless within the tolerance', () => {
        // (0.6, 0.8) is on the unit circle; scale it slightly outward.
        const p = v2(0.6 * 1.001, 0.8 * 1.001);
        expect(arc.contains(p, 0)).toBe(false);
        expect(arc.contains(p, 1e-6)).toBe(false);
        expect(arc.contains(p, 1e-2)).toBe(true);
    });

    it('rejects every point when the tolerance is negative', () => {
        // The upstream comment claims a negative epsilon behaves as zero, but
        // ||P-C| - r| >= 0 > epsilon always fails. Zero works as documented.
        expect(arc.contains(v2(1, 0), -1)).toBe(false);
        expect(arc.contains(v2(1, 0), 0)).toBe(true);
        expect(arc.contains(v2(1.0000001, 0), -1)).toBe(false);
    });

    it('handles arcs subtending more than pi radians', () => {
        // From (1,0) counterclockwise to (0,-1): three quarters of the circle.
        const big = Arc2.fromCenterRadiusEnds(v2(0, 0), 1, v2(1, 0),
            v2(0, -1));
        const s = Math.SQRT1_2;
        expect(big.containsOnCircle(v2(0, 1))).toBe(true);
        expect(big.containsOnCircle(v2(-1, 0))).toBe(true);
        expect(big.containsOnCircle(v2(-s, -s))).toBe(true);
        // The excluded quarter is between (0,-1) and (1,0).
        expect(big.containsOnCircle(v2(s, -s))).toBe(false);
    });

    it('agrees with the side-of-line criterion for a randomized sweep', () => {
        // Arc from angle 0.3 to angle 0.3 + 2.0 radians on a shifted circle.
        const c = v2(2, -1);
        const r = 3;
        const a0 = 0.3;
        const a1 = 0.3 + 2.0;
        const e0 = v2(c.get(0) + r * Math.cos(a0), c.get(1) + r * Math.sin(a0));
        const e1 = v2(c.get(0) + r * Math.cos(a1), c.get(1) + r * Math.sin(a1));
        const arcR = Arc2.fromCenterRadiusEnds(c, r, e0, e1);
        for (let k = 0; k < 64; ++k) {
            const t = (k / 64) * 2 * Math.PI;
            const p = v2(c.get(0) + r * Math.cos(t), c.get(1) + r * Math.sin(t));
            // Point is on the arc when the ccw sweep from a0 to t is <= 2.0.
            let sweep = t - a0;
            while (sweep < 0) {
                sweep += 2 * Math.PI;
            }
            const expected = sweep <= a1 - a0;
            expect(arcR.contains(p, 1e-9)).toBe(expected);
        }
    });
});

describe('Arc2 comparisons', () => {
    const base = new Arc2();

    it('equals compares all members', () => {
        expect(base.equals(new Arc2())).toBe(true);
        expect(base.notEquals(new Arc2())).toBe(false);

        const other = base.clone();
        other.radius = 2;
        expect(base.equals(other)).toBe(false);
        expect(base.notEquals(other)).toBe(true);

        const other2 = base.clone();
        other2.end[1] = v2(0, 2);
        expect(base.equals(other2)).toBe(false);
    });

    it('lessThan orders by center, then radius, then the ends', () => {
        const smallCenter = base.clone();
        smallCenter.center = v2(-1, 0);
        expect(smallCenter.lessThan(base)).toBe(true);
        expect(base.lessThan(smallCenter)).toBe(false);

        const smallRadius = base.clone();
        smallRadius.radius = 0.5;
        expect(smallRadius.lessThan(base)).toBe(true);

        const smallEnd0 = base.clone();
        smallEnd0.end[0] = v2(0.5, 0);
        expect(smallEnd0.lessThan(base)).toBe(true);

        const smallEnd1 = base.clone();
        smallEnd1.end[1] = v2(0, 0.5);
        expect(smallEnd1.lessThan(base)).toBe(true);
        expect(smallEnd1.greaterThan(base)).toBe(false);
    });

    it('the derived comparisons are consistent', () => {
        const other = base.clone();
        other.radius = 2;
        expect(base.lessThanOrEqual(other)).toBe(true);
        expect(base.lessThanOrEqual(base.clone())).toBe(true);
        expect(other.greaterThan(base)).toBe(true);
        expect(other.greaterThanOrEqual(base)).toBe(true);
        expect(base.greaterThanOrEqual(base.clone())).toBe(true);
        expect(base.greaterThan(base.clone())).toBe(false);
    });
});
