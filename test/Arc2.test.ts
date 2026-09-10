import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { Vector, sub, length } from '../src/Vector.js';
import { check, compareKeys, expectStrictWeakOrder, fc, finite, positive,
    scaled, vector } from './helpers/arbitraries.js';

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

describe('Arc2 verification', () => {
    const onCircle = (c: Vector, r: number, t: number) =>
        Vector.fromArray([c.get(0) + r * Math.cos(t),
            c.get(1) + r * Math.sin(t)]);
    const key = (a: Arc2) => [...a.center.values, a.radius,
        ...a.end[0].values, ...a.end[1].values];

    it('containsOnCircle accepts exactly the counterclockwise arc E0 -> E1',
        () => {
            // The chord test dotPerp(P-E0, E1-E0) >= 0 selects the points of
            // the circle on the counterclockwise arc from E0 to E1, for any
            // subtended angle (including angles >= pi).
            // The angles are drawn from a uniform grid: fc.double samples
            // the bit patterns of its range, so tiny magnitudes dominate and
            // the boundary-band precondition below would reject most draws.
            check(fc.tuple(vector(2, -5, 5), positive(5, 0.5),
                scaled(-Math.PI, Math.PI), scaled(0.05, 2 * Math.PI - 0.05),
                scaled(0, 2 * Math.PI)),
                ([c, r, a0, delta, s]) => {
                    // Stay away from the two endpoints where the sign is a
                    // rounding-error coin flip.
                    fc.pre(s > 1e-3 && Math.abs(s - delta) > 1e-3
                        && s < 2 * Math.PI - 1e-3);
                    const arc = Arc2.fromCenterRadiusEnds(c, r,
                        onCircle(c, r, a0), onCircle(c, r, a0 + delta));
                    const p = onCircle(c, r, a0 + s);
                    expect(arc.containsOnCircle(p)).toBe(s < delta);
                });
        });

    it('contains(P, eps) is the circle test and then the arc test', () => {
        check(fc.tuple(vector(2, -5, 5), positive(5, 0.5),
            finite(0.05, 2 * Math.PI - 0.05), finite(0, 2 * Math.PI),
            finite(0.5, 2), positive(0.5, 0.01)),
            ([c, r, delta, s, scale, eps]) => {
                const arc = Arc2.fromCenterRadiusEnds(c, r, onCircle(c, r, 0),
                    onCircle(c, r, delta));
                // A point at radius r*scale, generally off the circle.
                const p = onCircle(c, r * scale, s);
                const off = Math.abs(length(sub(p, c)) - r);
                expect(arc.contains(p, eps))
                    .toBe(off <= eps && arc.containsOnCircle(p));
            });
    });

    it('a negative tolerance rejects every point (upstream quirk #155)', () => {
        // The upstream comment promises a negative epsilon behaves like zero,
        // but ||P-C| - r| <= epsilon is false for every P when epsilon < 0.
        check(fc.tuple(finite(0, 2 * Math.PI), positive(1, 1e-6)),
            ([t, eps]) => {
                const c = Vector.fromArray([0, 0]);
                const arc = Arc2.fromCenterRadiusEnds(c, 1,
                    Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
                const p = onCircle(c, 1, t);
                expect(arc.contains(p, -eps)).toBe(false);
            });
    });

    it('the comparisons follow the (center, radius, end0, end1) order', () => {
        const arc = fc.tuple(vector(2, -2, 2), positive(3),
            vector(2, -2, 2), vector(2, -2, 2))
            .map(([c, r, e0, e1]) => Arc2.fromCenterRadiusEnds(c, r, e0, e1));
        check(fc.tuple(arc, arc), ([a, b]) => {
            const cmp = compareKeys(key(a), key(b));
            expect(a.lessThan(b)).toBe(cmp < 0);
            expect(a.greaterThan(b)).toBe(cmp > 0);
            expect(a.lessThanOrEqual(b)).toBe(cmp <= 0);
            expect(a.greaterThanOrEqual(b)).toBe(cmp >= 0);
            expect(a.equals(b)).toBe(cmp === 0);
        });
        check(fc.array(arc, { minLength: 4, maxLength: 5 }), arcs => {
            expectStrictWeakOrder(arcs, (x, y) => x.lessThan(y));
        }, 30);
    });

    it('equals is element equality, so NaN members break self-equality', () => {
        const arc = Arc2.fromCenterRadiusEnds(Vector.fromArray([0, 0]), 1,
            Vector.fromArray([NaN, 0]), Vector.fromArray([0, 1]));
        expect(arc.equals(arc)).toBe(false);
        expect(arc.lessThan(arc)).toBe(false);
        expect(arc.lessThanOrEqual(arc)).toBe(true);
    });

    it('the factory and clone are independent of their inputs', () => {
        check(fc.tuple(vector(2), vector(2), vector(2)), ([c, e0, e1]) => {
            const arc = Arc2.fromCenterRadiusEnds(c, 1, e0, e1);
            const cloned = arc.clone();
            c.set(0, 999);
            e0.set(1, 888);
            arc.end[1].set(0, 777);
            expect(arc.center.get(0)).not.toBe(999);
            expect(arc.end[0].get(1)).not.toBe(888);
            expect(cloned.end[1].get(0)).not.toBe(777);
        });
    });
});
