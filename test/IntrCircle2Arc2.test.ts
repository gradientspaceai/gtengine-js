import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { Hypersphere } from '../src/Hypersphere.js';
import { Vector, sub, length } from '../src/Vector.js';
import { IntrCircle2Arc2FI } from '../src/IntrCircle2Arc2.js';

const INT32_MAX = 2147483647;

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// A point on the unit circle at the given angle in degrees. The four axis
// directions are exact so that endpoint equality tests are exact.
function pt(deg: number): Vector {
    switch (((deg % 360) + 360) % 360) {
        case 0: return v2(1, 0);
        case 90: return v2(0, 1);
        case 180: return v2(-1, 0);
        case 270: return v2(0, -1);
        default: {
            const a = (deg * Math.PI) / 180;
            return v2(Math.cos(a), Math.sin(a));
        }
    }
}

// The arc of the unit circle from deg0 counterclockwise to deg1.
function arc(deg0: number, deg1: number): Arc2 {
    return Arc2.fromCenterRadiusEnds(v2(0, 0), 1, pt(deg0), pt(deg1));
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

describe('IntrCircle2Arc2', () => {
    const fi = new IntrCircle2Arc2FI();

    it('reports the arc when the circle contains it (cocircular)', () => {
        const a = arc(0, 90);
        const result = fi.find(circle(0, 0, 1), a);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(INT32_MAX);
        expect(result.arc.equals(a)).toBe(true);
    });

    it('reports no intersection for disjoint circles', () => {
        const result = fi.find(circle(5, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('reports no intersection when a nested circle misses', () => {
        const result = fi.find(circle(0, 0, 0.5), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('finds two points when both circle crossings are on the arc', () => {
        // The circle of radius 1 centered at (1,1) meets the unit circle at
        // (1,0) and (0,1), both endpoints of the first-quadrant arc.
        const result = fi.find(circle(1, 1, 1), arc(0, 90));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(length(result.point[i])).toBeCloseTo(1, 12);
        }
    });

    it('finds one point when only one crossing is on the arc', () => {
        // The unit circle centered at (1,1) meets the unit circle at (1,0)
        // and (0,1). The arc from 90 to 270 degrees contains only (0,1).
        const result = fi.find(circle(1, 1, 1), arc(90, 270));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('finds the tangent contact when it lies on the arc', () => {
        // The circle of radius 1 centered at (2,0) is tangent to the unit
        // circle at (1,0), which is the start of the first-quadrant arc.
        const result = fi.find(circle(2, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a tangent contact that is not on the arc', () => {
        // Tangency at (-1,0), which is not on the first-quadrant arc.
        const result = fi.find(circle(-2, 0, 1), arc(0, 90));
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
    });

    it('returns points on both the circle and the arc', () => {
        const rnd = makeRandom(60607);
        let numFound = 0;
        for (let k = 0; k < 400; ++k) {
            const a = arc(rnd() * 360, rnd() * 360);
            const c = circle(rnd() * 4 - 2, rnd() * 4 - 2, 0.2 + rnd() * 2);
            const result = fi.find(c, a);
            if (result.numIntersections === INT32_MAX) {
                continue;
            }
            expect(result.intersect).toBe(result.numIntersections > 0);
            for (let i = 0; i < result.numIntersections; ++i) {
                const p = result.point[i];
                expect(length(sub(p, c.center))).toBeCloseTo(c.radius, 8);
                expect(length(sub(p, a.center))).toBeCloseTo(a.radius, 8);
                expect(a.containsOnCircle(p)).toBe(true);
                ++numFound;
            }
        }
        expect(numFound).toBeGreaterThan(0);
    });

    it('handles a zero-radius circle', () => {
        // A degenerate circle at a point of the arc.
        const result = fi.find(circle(1, 0, 0), arc(0, 90));
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);

        // A degenerate circle away from the unit circle.
        expect(fi.find(circle(0.5, 0, 0), arc(0, 90)).intersect).toBe(false);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, wellScaledVector, expectClose, expectVectorClose
} from './helpers/arbitraries.js';
import { IntrCircle2Circle2FI } from '../src/IntrCircle2Circle2.js';

const INT32_MAX_V = 2147483647;

function v2v(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function arcAtV31(center: Vector, radius: number, a0: number,
    a1: number): Arc2 {
    const on = (a: number): Vector => v2v(center.values[0]
        + radius * Math.cos(a), center.values[1] + radius * Math.sin(a));
    return Arc2.fromCenterRadiusEnds(center, radius, on(a0), on(a1));
}

const circleArc = fc.tuple(wellScaledVector(2, -3, 3),
    fc.double({ min: 0.5, max: 3, noNaN: true }),
    wellScaledVector(2, -3, 3),
    fc.double({ min: 0.5, max: 3, noNaN: true }),
    fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
    fc.double({ min: 0.2, max: 6, noNaN: true }))
    .map(([cc, cr, ac, ar, a0, span]) => ({
        circle: Hypersphere.fromCenterRadius(cc, cr),
        arc: arcAtV31(ac, ar, a0, a0 + span)
    }));

describe('IntrCircle2Arc2 verification', () => {
    const fiq = new IntrCircle2Arc2FI();
    const ccq = new IntrCircle2Circle2FI();

    it('reports exactly the circle-circle points that are on the arc', () => {
        check(circleArc, ({ circle, arc }) => {
            const cc = ccq.find(circle,
                Hypersphere.fromCenterRadius(arc.center, arc.radius));
            const f = fiq.find(circle, arc);
            if (!cc.intersect) {
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
                return;
            }
            if (cc.numIntersections === INT32_MAX_V) {
                // The arc is on the circle.
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(INT32_MAX_V);
                expect(f.arc.equals(arc)).toBe(true);
                return;
            }
            let expected = 0;
            for (let i = 0; i < cc.numIntersections; ++i) {
                if (arc.containsOnCircle(cc.point[i])) {
                    ++expected;
                }
            }
            expect(f.numIntersections).toBe(expected);
            expect(f.intersect).toBe(expected > 0);
        });
    });

    it('the reported points lie on the circle and on the arc', () => {
        check(circleArc, ({ circle, arc }) => {
            const f = fiq.find(circle, arc);
            if (!f.intersect || f.numIntersections === INT32_MAX_V) {
                return;
            }
            for (let i = 0; i < f.numIntersections; ++i) {
                const p = f.point[i];
                expectClose(length(sub(p, circle.center)), circle.radius,
                    1e-7, 1e-9);
                expectClose(length(sub(p, arc.center)), arc.radius,
                    1e-7, 1e-9);
                expect(arc.containsOnCircle(p)).toBe(true);
                for (const x of p.values) {
                    expect(Number.isNaN(x)).toBe(false);
                }
            }
        });
    });

    it('an arc on the circle returns a copy of the arc', () => {
        check(fc.tuple(wellScaledVector(2, -3, 3),
            fc.double({ min: 0.5, max: 3, noNaN: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            fc.double({ min: 0.2, max: 6, noNaN: true })),
            ([c, r, a0, span]) => {
                const arc = arcAtV31(c, r, a0, a0 + span);
                const circle = Hypersphere.fromCenterRadius(c, r);
                const f = fiq.find(circle, arc);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(INT32_MAX_V);
                expect(f.arc.equals(arc)).toBe(true);
                // A copy, not an alias.
                f.arc.center.set(0, f.arc.center.get(0) + 1);
                expect(arc.center.get(0)).not.toBe(f.arc.center.get(0));
            });
    });

    it('a disjoint circle reports no intersection', () => {
        check(fc.tuple(wellScaledVector(2, -3, 3),
            fc.double({ min: 0.5, max: 2, noNaN: true }),
            fc.double({ min: 0.5, max: 2, noNaN: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            fc.double({ min: 0.2, max: 6, noNaN: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true })),
            ([c0, r0, r1, a0, span, dirAngle]) => {
                const gap = r0 + r1 + 1;
                const c1 = v2v(c0.values[0] + gap * Math.cos(dirAngle),
                    c0.values[1] + gap * Math.sin(dirAngle));
                const f = fiq.find(Hypersphere.fromCenterRadius(c0, r0),
                    arcAtV31(c1, r1, a0, a0 + span));
                expect(f.intersect).toBe(false);
                expect(f.numIntersections).toBe(0);
                expectVectorClose(f.point[0], v2v(0, 0), 0, 0);
            });
    });

    it('a circle through a sampled arc point hits that point', () => {
        check(fc.tuple(wellScaledVector(2, -3, 3),
            fc.double({ min: 0.5, max: 3, noNaN: true }),
            fc.double({ min: -Math.PI, max: Math.PI, noNaN: true }),
            fc.double({ min: 0.3, max: 3, noNaN: true }),
            fc.double({ min: 0.15, max: 0.85, noNaN: true }),
            fc.double({ min: 0.3, max: 2, noNaN: true })),
            ([c, r, a0, span, frac, cr]) => {
                const arc = arcAtV31(c, r, a0, a0 + span);
                const target = v2v(c.values[0] + r * Math.cos(a0 + frac * span),
                    c.values[1] + r * Math.sin(a0 + frac * span));
                // A circle of radius cr centred so that 'target' is on it.
                const dir = a0 + frac * span;
                const cc = v2v(target.values[0] + cr * Math.cos(dir + 1),
                    target.values[1] + cr * Math.sin(dir + 1));
                const f = fiq.find(Hypersphere.fromCenterRadius(cc, cr), arc);
                expect(f.intersect).toBe(true);
                let best = Infinity;
                for (let i = 0; i < f.numIntersections; ++i) {
                    best = Math.min(best, length(sub(f.point[i], target)));
                }
                expect(best).toBeLessThan(1e-7);
            });
    });
});
