import { describe, it, expect } from 'vitest';
import { Arc2 } from '../src/Arc2.js';
import { IntrLine2Arc2TI, IntrLine2Arc2FI } from '../src/IntrLine2Arc2.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function line(origin: number[], direction: number[]): Line {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Line.fromOriginDirection(Vector.fromArray(origin), d);
}

// The quarter arc of the unit circle from (1,0) to (0,1), which contains the
// points (cos t, sin t) for t in [0, pi/2].
function quarterArc(): Arc2 {
    return Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0), vec(0, 1));
}

describe('IntrLine2Arc2', () => {
    const ti = new IntrLine2Arc2TI();
    const fi = new IntrLine2Arc2FI();

    it('finds the single intersection of a chord that crosses the arc once', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 0.5], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);

        const expectedX = Math.sqrt(0.75);
        expect(result.point[0].values[0]).toBeCloseTo(expectedX, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0.5, 12);
        expect(result.parameter[0]).toBeCloseTo(expectedX, 12);
        expect(ti.test(line([0, 0.5], [1, 0]), arc).intersect).toBe(true);
    });

    it('finds two intersections when both circle points are on the arc', () => {
        // The arc from (1,0) to (-1,0) going through (0,1) is the upper half
        // circle; the horizontal line y = 0.5 meets it twice.
        const arc = Arc2.fromCenterRadiusEnds(vec(0, 0), 1, vec(1, 0),
            vec(-1, 0));
        const result = fi.find(line([0, 0.5], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        for (let i = 0; i < 2; ++i) {
            expect(result.point[i].values[1]).toBeCloseTo(0.5, 12);
            expect(Math.abs(result.point[i].values[0]))
                .toBeCloseTo(Math.sqrt(0.75), 12);
        }
        // Parameters are ordered as the underlying line-circle query orders
        // them, smallest first.
        expect(result.parameter[0]).toBeLessThan(result.parameter[1]);
    });

    it('accepts an arc endpoint (the on-boundary test is inclusive)', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 0], [0, 1]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('reports no intersection when the circle points miss the arc', () => {
        const arc = quarterArc();
        // The line y = -0.5 meets the circle only in the lower half plane.
        const result = fi.find(line([0, -0.5], [1, 0]), arc);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(line([0, -0.5], [1, 0]), arc).intersect).toBe(false);
    });

    it('reports no intersection when the line misses the circle entirely', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 3], [1, 0]), arc);
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(ti.test(line([0, 3], [1, 0]), arc).intersect).toBe(false);
    });

    it('reports the tangent point when it lies on the arc', () => {
        const arc = quarterArc();
        const result = fi.find(line([0, 1], [1, 0]), arc);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(0, 9);
        expect(result.point[0].values[1]).toBeCloseTo(1, 12);
    });

    it('agrees with a direct arc-membership check on random lines', () => {
        let seed = 24680;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const arc = Arc2.fromCenterRadiusEnds(vec(0.5, -0.25), 2,
            add(vec(0.5, -0.25), mul(2, vec(Math.cos(0.3), Math.sin(0.3)))),
            add(vec(0.5, -0.25), mul(2, vec(Math.cos(2.1), Math.sin(2.1)))));

        for (let trial = 0; trial < 200; ++trial) {
            const l = line([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 2 - 1, rand() * 2 - 1]);
            if (dot(l.direction, l.direction) < 0.5) {
                continue;
            }
            const result = fi.find(l, arc);
            expect(ti.test(l, arc).intersect).toBe(result.intersect);

            for (let i = 0; i < result.numIntersections; ++i) {
                const p = result.point[i];
                // The point is on the line at the reported parameter.
                const q = add(l.origin, mul(result.parameter[i], l.direction));
                expect(sub(p, q).values[0]).toBeCloseTo(0, 9);
                expect(sub(p, q).values[1]).toBeCloseTo(0, 9);

                // The point is on the circle of the arc.
                const r = sub(p, arc.center);
                expect(Math.sqrt(dot(r, r))).toBeCloseTo(arc.radius, 9);

                // The point is on the arc.
                expect(arc.containsOnCircle(p)).toBe(true);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrLine2Arc2.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';
import { IntrLine2Circle2FI } from '../src/IntrLine2Circle2.js';
import { Hypersphere } from '../src/Hypersphere.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so a direction built from
// (cos a, sin a) never has a subnormal component; a subnormal component makes
// the exact DotPerp tests inside the line-line query underflow.
const angle = () => wellScaled(-Math.PI, Math.PI);

// An arc that is a proper sub-arc of its circle: the center is moderately
// scaled, the radius is bounded away from zero, and the traversal span from
// end[0] to end[1] is in (0.2, 2*pi - 0.2) radians, so neither endpoint
// coincides with the other.
const arcArb = fc.tuple(
    wellScaled(-3, 3), wellScaled(-3, 3),
    fc.double({ min: 0.5, max: 3, noNaN: true, noDefaultInfinity: true }),
    angle(),
    fc.double({ min: 0.2, max: 2 * Math.PI - 0.2, noNaN: true,
        noDefaultInfinity: true })
).map(([cx, cy, r, a0, span]) => {
    const c = vec(cx, cy);
    const e0 = vec(cx + r * Math.cos(a0), cy + r * Math.sin(a0));
    const a1 = a0 + span;
    const e1 = vec(cx + r * Math.cos(a1), cy + r * Math.sin(a1));
    return { arc: Arc2.fromCenterRadiusEnds(c, r, e0, e1), a0, span };
});

const lineArb = fc.tuple(wellScaled(-4, 4), wellScaled(-4, 4), angle())
    .map(([ox, oy, a]) =>
        Line.fromOriginDirection(vec(ox, oy),
            vec(Math.cos(a), Math.sin(a))));

const lineArcArb = fc.tuple(lineArb, arcArb)
    .map(([l, a]) => ({ line: l, ...a }));

// A point on the arc at fraction u of the traversal from end[0] to end[1].
function arcPoint(arc: Arc2, a0: number, span: number, u: number): Vector {
    const a = a0 + span * u;
    return add(arc.center,
        vec(arc.radius * Math.cos(a), arc.radius * Math.sin(a)));
}

describe('IntrLine2Arc2 verification', () => {
    const tiq = new IntrLine2Arc2TI();
    const fiq = new IntrLine2Arc2FI();
    const lcq = new IntrLine2Circle2FI();

    it('TI and FI always agree on intersect', () => {
        check(lineArcArb, ({ line: l, arc }) => {
            expect(tiq.test(l, arc).intersect).toBe(fiq.find(l, arc).intersect);
        });
    });

    it('reported points lie on the line, on the circle and on the arc', () => {
        check(lineArcArb, ({ line: l, arc }) => {
            const f = fiq.find(l, arc);
            expect(f.intersect).toBe(f.numIntersections > 0);
            for (let i = 0; i < f.numIntersections; ++i) {
                expectVectorClose(f.point[i],
                    add(l.origin, mul(f.parameter[i], l.direction)), 0, 0);
                // The parameter carries a square root, so the residual of
                // |P - C| = r scales with the conditioning of that root.
                expectClose(length(sub(f.point[i], arc.center)), arc.radius,
                    1e-8, 1e-9);
                expect(arc.containsOnCircle(f.point[i])).toBe(true);
            }
            // Unused slots keep their default values (upstream fills only
            // point[i] for i < numIntersections here).
            for (let i = f.numIntersections; i < 2; ++i) {
                expect(f.parameter[i]).toBe(0);
                expect(f.point[i].get(0)).toBe(0);
                expect(f.point[i].get(1)).toBe(0);
            }
        });
    });

    it('the hits are exactly the line-circle hits that lie on the arc', () => {
        check(lineArcArb, ({ line: l, arc }) => {
            const circle = Hypersphere.fromCenterRadius(arc.center, arc.radius);
            const lc = lcq.find(l, circle);
            const kept: number[] = [];
            for (let i = 0; i < lc.numIntersections; ++i) {
                if (arc.containsOnCircle(lc.point[i])) {
                    kept.push(lc.parameter[i]);
                }
            }
            const f = fiq.find(l, arc);
            expect(f.numIntersections).toBe(kept.length);
            for (let i = 0; i < kept.length; ++i) {
                expect(f.parameter[i]).toBe(kept[i]);
            }
        });
    });

    it('a line through an interior arc point always reports that point', () => {
        check(fc.tuple(arcArb,
            fc.double({ min: 0.15, max: 0.85, noNaN: true,
                noDefaultInfinity: true }),
            angle()),
        ([{ arc, a0, span }, u, ang]) => {
            const p = arcPoint(arc, a0, span, u);
            const l = Line.fromOriginDirection(p,
                vec(Math.cos(ang), Math.sin(ang)));
            const f = fiq.find(l, arc);
            expect(f.intersect).toBe(true);
            let found = false;
            for (let i = 0; i < f.numIntersections; ++i) {
                if (length(sub(f.point[i], p)) <= 1e-7 * (1 + arc.radius)) {
                    found = true;
                }
            }
            expect(found).toBe(true);
        });
    });

    it('a dense sweep of the arc matches the reported crossings', () => {
        const rnd = seededRandom(0x4a2c17);
        for (let trial = 0; trial < 150; ++trial) {
            const cx = rnd() * 4 - 2, cy = rnd() * 4 - 2;
            const r = 0.5 + rnd() * 2;
            const a0 = rnd() * 2 * Math.PI;
            const span = 0.3 + rnd() * (2 * Math.PI - 0.6);
            const arc = Arc2.fromCenterRadiusEnds(vec(cx, cy), r,
                vec(cx + r * Math.cos(a0), cy + r * Math.sin(a0)),
                vec(cx + r * Math.cos(a0 + span),
                    cy + r * Math.sin(a0 + span)));
            const ang = rnd() * 2 * Math.PI;
            const d = vec(Math.cos(ang), Math.sin(ang));
            const l = Line.fromOriginDirection(
                vec(rnd() * 6 - 3, rnd() * 6 - 3), d);

            // Signed distance of a circle point from the line, using the
            // normal n = Perp(D). The arc crossings are the sign changes of g
            // over the traversal interval.
            const n = vec(d.get(1), -d.get(0));
            const c0 = dot(n, l.origin);
            const g = (t: number): number =>
                dot(n, add(arc.center, vec(r * Math.cos(t), r * Math.sin(t))))
                - c0;

            // Skip near-tangency and configurations whose crossing sits
            // within round-off of an arc endpoint; there the discrete count
            // is genuinely ambiguous.
            const centerDist = Math.abs(dot(n, arc.center) - c0);
            if (Math.abs(centerDist - r) < 1e-3) { continue; }
            if (Math.abs(g(a0)) < 1e-3 || Math.abs(g(a0 + span)) < 1e-3) {
                continue;
            }

            const steps = 4000;
            const roots: number[] = [];
            let prev = g(a0);
            for (let k = 1; k <= steps; ++k) {
                const t = a0 + (span * k) / steps;
                const cur = g(t);
                if ((prev < 0 && cur > 0) || (prev > 0 && cur < 0)) {
                    let lo = a0 + (span * (k - 1)) / steps, hi = t;
                    for (let b = 0; b < 60; ++b) {
                        const mid = 0.5 * (lo + hi);
                        if (g(lo) * g(mid) <= 0) { hi = mid; } else { lo = mid; }
                    }
                    roots.push(0.5 * (lo + hi));
                }
                prev = cur;
            }

            const f = fiq.find(l, arc);
            expect(f.numIntersections).toBe(roots.length);
            for (const t of roots) {
                const p = add(arc.center,
                    vec(r * Math.cos(t), r * Math.sin(t)));
                let found = false;
                for (let i = 0; i < f.numIntersections; ++i) {
                    if (length(sub(f.point[i], p)) < 1e-5 * (1 + r)) {
                        found = true;
                    }
                }
                expect(found).toBe(true);
            }
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(lineArcArb, angle(), wellScaled(-2, 2),
            wellScaled(-2, 2)),
        ([{ line: l, arc }, ang, tx, ty]) => {
            const ca = Math.cos(ang), sa = Math.sin(ang);
            const rot = (v: Vector): Vector => vec(
                ca * v.get(0) - sa * v.get(1),
                sa * v.get(0) + ca * v.get(1));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));

            const l2 = Line.fromOriginDirection(xf(l.origin), rot(l.direction));
            const arc2 = Arc2.fromCenterRadiusEnds(xf(arc.center), arc.radius,
                xf(arc.end[0]), xf(arc.end[1]));

            const f0 = fiq.find(l, arc);
            const f1 = fiq.find(l2, arc2);
            if (f0.numIntersections !== f1.numIntersections) {
                // Only a grazing configuration (tangency, or a hit at an arc
                // endpoint) may flip the count under a rotation.
                return;
            }
            for (let i = 0; i < f0.numIntersections; ++i) {
                expectVectorClose(f1.point[i], xf(f0.point[i]), 1e-7, 1e-8);
            }
        });
    });

    it('a line missing the circle leaves every field at its default', () => {
        const arc = quarterArc();
        const f = fiq.find(line([0, 5], [1, 0]), arc);
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(0);
        expect(f.parameter).toEqual([0, 0]);
        expect(f.point[0].values).toEqual([0, 0]);
        expect(f.point[1].values).toEqual([0, 0]);
    });

    it('a line meeting the circle only off the arc reports no hits', () => {
        // The quarter arc spans (cos t, sin t) for t in [0, pi/2]; the line
        // y = -0.5 meets the unit circle only in the lower half plane.
        const f = fiq.find(line([0, -0.5], [1, 0]), quarterArc());
        expect(f.intersect).toBe(false);
        expect(f.numIntersections).toBe(0);
    });
});
