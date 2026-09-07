import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import {
    IntrSegment2AlignedBox2TI,
    IntrSegment2AlignedBox2FI
} from '../src/IntrSegment2AlignedBox2.js';
import {
    IntrSegment2OrientedBox2TI,
    IntrSegment2OrientedBox2FI,
    defaultIntrSegment2OrientedBox2FIResult
} from '../src/IntrSegment2OrientedBox2.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function rotatedBox(center: Vector, angle: number, extent: Vector):
    OrientedBox {
    const c = Math.cos(angle), s = Math.sin(angle);
    return OrientedBox.fromCenterAxisExtent(center, [vec(c, s), vec(-s, c)],
        extent);
}

function boxSignedDepth(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const value = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (value > worst) {
            worst = value;
        }
    }
    return worst;
}

describe('IntrSegment2OrientedBox2', () => {
    const ti = new IntrSegment2OrientedBox2TI();
    const fi = new IntrSegment2OrientedBox2FI();
    const unitAxes = [vec(1, 0), vec(0, 1)];

    it('matches the aligned-box query when the axes are standard', () => {
        const center = vec(1, -2);
        const extent = vec(2, 0.5);
        const obox = OrientedBox.fromCenterAxisExtent(center, unitAxes,
            extent);
        const abox = AlignedBox.fromMinMax(sub(center, extent),
            add(center, extent));
        const aTI = new IntrSegment2AlignedBox2TI();
        const aFI = new IntrSegment2AlignedBox2FI();

        const segments = [
            segment([-10, -2], [10, -2]),
            segment([1, -10], [1, 10]),
            segment([-10, 10], [10, 10]),
            segment([1, -2], [4, -2]),
            segment([-10, -2], [-5, -2])
        ];
        for (const s of segments) {
            expect(ti.test(s, obox).intersect)
                .toBe(aTI.test(s, abox).intersect);
            const o = fi.find(s, obox);
            const a = aFI.find(s, abox);
            expect(o.intersect).toBe(a.intersect);
            expect(o.numIntersections).toBe(a.numIntersections);
            for (let i = 0; i < o.numIntersections; ++i) {
                expect(o.parameter[i]).toBeCloseTo(a.parameter[i], 12);
                expect(o.cdeParameter[i]).toBeCloseTo(a.cdeParameter[i], 12);
                expect(o.point[i].values[0])
                    .toBeCloseTo(a.point[i].values[0], 12);
                expect(o.point[i].values[1])
                    .toBeCloseTo(a.point[i].values[1], 12);
            }
        }
    });

    it('reports correct world-space points for a rotated box', () => {
        // Regression test for the upstream point-reconstruction bug: with a
        // 45-degree rotation the box is a diamond reaching sqrt(2) along x.
        const box = rotatedBox(vec(0, 0), Math.PI / 4, vec(1, 1));
        const result = fi.find(segment([-5, 0], [5, 0]), box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[0]).toBeCloseTo(-Math.SQRT2, 9);
        expect(result.point[0].values[1]).toBeCloseTo(0, 9);
        expect(result.point[1].values[0]).toBeCloseTo(Math.SQRT2, 9);
        expect(result.point[1].values[1]).toBeCloseTo(0, 9);

        // The endpoint-form parameters recover the same points.
        for (let i = 0; i < 2; ++i) {
            const t = result.parameter[i];
            const Q = add(mul(1 - t, vec(-5, 0)), mul(t, vec(5, 0)));
            expect(Q.values[0]).toBeCloseTo(result.point[i].values[0], 9);
            expect(Q.values[1]).toBeCloseTo(result.point[i].values[1], 9);
        }
    });

    it('reports parameters in [0,1] and the centered form in cdeParameter', () => {
        const box = rotatedBox(vec(1, 1), 0.4, vec(1, 0.5));
        const s = segment([-4, 1], [6, 1]);
        const result = fi.find(s, box);
        expect(result.intersect).toBe(true);
        const { center, direction, extent } = s.getCenteredForm();
        for (let i = 0; i < result.numIntersections; ++i) {
            expect(result.parameter[i]).toBeGreaterThanOrEqual(0);
            expect(result.parameter[i]).toBeLessThanOrEqual(1);
            expect(Math.abs(result.cdeParameter[i]))
                .toBeLessThanOrEqual(extent + 1e-9);
            const P = add(center, mul(result.cdeParameter[i], direction));
            expect(P.values[0]).toBeCloseTo(result.point[i].values[0], 9);
            expect(P.values[1]).toBeCloseTo(result.point[i].values[1], 9);
        }
    });

    it('rejects a segment that stops short of the box', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const s = segment([-10, 0], [-5, 0]);
        expect(fi.find(s, box).intersect).toBe(false);
        expect(ti.test(s, box).intersect).toBe(false);
    });

    it('handles a degenerate segment (upstream reports 2 intersections)', () => {
        const box = rotatedBox(vec(0, 0), 0.3, vec(1, 1));
        const inside = segment([0.1, -0.2], [0.1, -0.2]);
        const result = fi.find(inside, box);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter).toEqual([0, 0]);
        expect(result.cdeParameter).toEqual([0, 0]);
        expect(result.point[0].values).toEqual([0.1, -0.2]);
        expect(result.point[1].values).toEqual([0.1, -0.2]);

        const outside = segment([9, 9], [9, 9]);
        expect(fi.find(outside, box).intersect).toBe(false);
    });

    it('has the documented default result', () => {
        const result = defaultIntrSegment2OrientedBox2FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.cdeParameter).toEqual([0, 0]);
    });

    it('rejects non-2D boxes', () => {
        const box3 = new OrientedBox(3);
        const s = segment([0, 0], [1, 1]);
        expect(() => ti.test(s, box3)).toThrow();
        expect(() => fi.find(s, box3)).toThrow();
    });

    it('agrees with brute-force sampling on random segments', () => {
        let seed = 20240901;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = rotatedBox(vec(0.5, -0.25), 0.7, vec(1.5, 0.75));
        for (let trial = 0; trial < 400; ++trial) {
            const s = segment([rand() * 8 - 4, rand() * 8 - 4],
                [rand() * 8 - 4, rand() * 8 - 4]);
            const result = fi.find(s, box);
            expect(ti.test(s, box).intersect).toBe(result.intersect);

            let sampledHit = false;
            for (let k = 0; k <= 400; ++k) {
                const t = k / 400;
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                if (boxSignedDepth(box, P) < -1e-6) {
                    sampledHit = true;
                    break;
                }
            }
            if (sampledHit) {
                expect(result.intersect).toBe(true);
            }

            for (let i = 0; i < result.numIntersections; ++i) {
                const t = result.parameter[i];
                expect(t).toBeGreaterThanOrEqual(-1e-9);
                expect(t).toBeLessThanOrEqual(1 + 1e-9);
                const P = add(s.p[0], mul(t, sub(s.p[1], s.p[0])));
                expect(P.values[0]).toBeCloseTo(result.point[i].values[0], 8);
                expect(P.values[1]).toBeCloseTo(result.point[i].values[1], 8);
                expect(boxSignedDepth(box, result.point[i]))
                    .toBeLessThan(1e-8);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSegment2OrientedBox2.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';
import { IntrLine2OrientedBox2FI } from '../src/IntrLine2OrientedBox2.js';
import { Line } from '../src/Line.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so a direction built from
// (cos a, sin a) never has a subnormal component; a subnormal component makes
// the exact DotPerp tests inside the line-line query underflow.
const angle2 = () => wellScaled(-Math.PI, Math.PI);

const extent2 = () => fc.double(
    { min: 0.25, max: 3, noNaN: true, noDefaultInfinity: true });

const segBox = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-3, 3), wellScaled(-3, 3), angle2(), extent2(), extent2()
).filter(([x0, y0, x1, y1]) =>
    (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0) > 1e-4)
    .map(([x0, y0, x1, y1, cx, cy, ab, e0, e1]) => ({
        segment: Segment.fromEndpoints(vec(x0, y0), vec(x1, y1)),
        box: rotatedBox(vec(cx, cy), ab, vec(e0, e1))
    }));

// The largest signed slab depth: negative strictly inside the box.
function obDepth2(box: OrientedBox, P: Vector): number {
    const diff = sub(P, box.center);
    let worst = -Number.MAX_VALUE;
    for (let i = 0; i < 2; ++i) {
        const v = Math.abs(dot(diff, box.axis[i])) - box.extent.values[i];
        if (v > worst) { worst = v; }
    }
    return worst;
}

// obDepth2 is a maximum of convex functions, so it is convex along the
// segment; a ternary search finds its minimum reliably. A value near zero
// means the segment grazes the box, where the TI slab test and the FI
// interval clip round differently and the discrete answers are ambiguous.
function minSegDepth(box: OrientedBox, s: Segment): number {
    const e = sub(s.p[1], s.p[0]);
    let lo = 0, hi = 1;
    for (let k = 0; k < 200; ++k) {
        const a = lo + (hi - lo) / 3, b = hi - (hi - lo) / 3;
        const fa = obDepth2(box, add(s.p[0], mul(a, e)));
        const fb = obDepth2(box, add(s.p[0], mul(b, e)));
        if (fa < fb) { hi = b; } else { lo = a; }
    }
    return obDepth2(box, add(s.p[0], mul(0.5 * (lo + hi), e)));
}

describe('IntrSegment2OrientedBox2 verification', () => {
    const tiq = new IntrSegment2OrientedBox2TI();
    const fiq = new IntrSegment2OrientedBox2FI();
    const lbq = new IntrLine2OrientedBox2FI();

    it('TI and FI always agree on intersect', () => {
        check(segBox, ({ segment: s, box: b }) => {
            if (Math.abs(minSegDepth(b, s)) < 1e-6) {
                return;   // grazing; the two formulations round differently
            }
            expect(tiq.test(s, b).intersect).toBe(fiq.find(s, b).intersect);
        });
    });

    it('reported points lie on the segment and in the closed box', () => {
        check(segBox, ({ segment: s, box: b }) => {
            const cf = s.getCenteredForm();
            const f = fiq.find(s, b);
            if (!f.intersect) {
                expect(f.numIntersections).toBe(0);
                return;
            }
            // The aligned-box DoQuery promotes a single-point overlap to two
            // reported intersections, so the count is always 2 for a
            // non-degenerate segment that meets the box.
            expect(f.numIntersections).toBe(2);
            for (let i = 0; i < f.numIntersections; ++i) {
                // 'parameter' is the endpoint form; 'cdeParameter' is the
                // centered form. Both must reconstruct the same point.
                expect(f.parameter[i]).toBeGreaterThanOrEqual(-1e-12);
                expect(f.parameter[i]).toBeLessThanOrEqual(1 + 1e-12);
                expect(Math.abs(f.cdeParameter[i]))
                    .toBeLessThanOrEqual(cf.extent + 1e-12);
                expectClose(f.parameter[i],
                    (f.cdeParameter[i] / cf.extent + 1) * 0.5, 1e-12, 1e-12);
                expectVectorClose(f.point[i],
                    add(s.p[0], mul(f.parameter[i], sub(s.p[1], s.p[0]))),
                    1e-9, 1e-9);
                expectVectorClose(f.point[i],
                    add(cf.center, mul(f.cdeParameter[i], cf.direction)),
                    0, 0);
                expect(obDepth2(b, f.point[i])).toBeLessThanOrEqual(1e-9);
            }
        });
    });

    it('the segment hit is the line hit clipped to |t| <= extent', () => {
        check(segBox, ({ segment: s, box: b }) => {
            const cf = s.getCenteredForm();
            const line = Line.fromOriginDirection(cf.center, cf.direction);
            const lb = lbq.find(line, b);
            const f = fiq.find(s, b);
            if (!lb.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const t0 = Math.max(lb.parameter[0], -cf.extent);
            const t1 = Math.min(lb.parameter[1], cf.extent);
            if (t0 > t1) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.cdeParameter[0] + 0).toBe(t0 + 0);
            expect(f.cdeParameter[1] + 0).toBe((t0 < t1 ? t1 : t0) + 0);
        });
    });

    it('a segment strictly inside the box is clipped to its own endpoints',
        () => {
            check(fc.tuple(wellScaled(-3, 3), wellScaled(-3, 3), angle2(),
                extent2(), extent2(),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true })),
            ([cx, cy, ab, e0, e1, a0, a1, b0, b1]) => {
                const box = rotatedBox(vec(cx, cy), ab, vec(e0, e1));
                const p0 = add(box.center, add(mul(a0 * e0, box.axis[0]),
                    mul(a1 * e1, box.axis[1])));
                const p1 = add(box.center, add(mul(b0 * e0, box.axis[0]),
                    mul(b1 * e1, box.axis[1])));
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const s = Segment.fromEndpoints(p0, p1);
                expect(tiq.test(s, box).intersect).toBe(true);
                const f = fiq.find(s, box);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], 0, 1e-12, 1e-12);
                expectClose(f.parameter[1], 1, 1e-12, 1e-12);
                expectVectorClose(f.point[0], p0, 1e-9, 1e-9);
                expectVectorClose(f.point[1], p1, 1e-9, 1e-9);
            });
        });

    it('a fine sweep of the segment agrees with the reported interval', () => {
        const rnd = seededRandom(0x9a3b21d);
        for (let trial = 0; trial < 200; ++trial) {
            const b = rotatedBox(vec(rnd() * 4 - 2, rnd() * 4 - 2),
                rnd() * 2 * Math.PI,
                vec(0.3 + rnd() * 2, 0.3 + rnd() * 2));
            const s = Segment.fromEndpoints(
                vec(rnd() * 10 - 5, rnd() * 10 - 5),
                vec(rnd() * 10 - 5, rnd() * 10 - 5));
            if (length(sub(s.p[1], s.p[0])) < 1e-2) { continue; }
            const f = fiq.find(s, b);
            for (let k = 0; k <= 1000; ++k) {
                const u = k / 1000;
                const p = add(s.p[0], mul(u, sub(s.p[1], s.p[0])));
                if (obDepth2(b, p) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(u).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(u).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('reports world-space points for every rotation of the box', () => {
        // Upstream computes box.center + (segOrigin + t*segDirection) with
        // segOrigin/segDirection in the box frame; the port evaluates the
        // world-space centered form instead (upstream issue #255). This
        // property fails on the upstream formula for any rotated box.
        check(fc.tuple(angle2(), wellScaled(-2, 2), wellScaled(-2, 2),
            extent2(), extent2(), angle2()),
        ([ab, cx, cy, e0, e1, ad]) => {
            const box = rotatedBox(vec(cx, cy), ab, vec(e0, e1));
            const d = vec(Math.cos(ad), Math.sin(ad));
            const s = Segment.fromEndpoints(sub(box.center, mul(20, d)),
                add(box.center, mul(20, d)));
            const f = fiq.find(s, box);
            expect(f.intersect).toBe(true);
            for (let i = 0; i < f.numIntersections; ++i) {
                // On the box boundary: the largest slab depth is zero.
                expectClose(obDepth2(box, f.point[i]), 0, 1e-9, 1e-9);
                expectVectorClose(f.point[i],
                    add(s.p[0], mul(f.parameter[i], sub(s.p[1], s.p[0]))),
                    1e-8, 1e-9);
            }
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(segBox, angle2(), wellScaled(-2, 2), wellScaled(-2, 2)),
            ([{ segment: s, box: b }, ang, tx, ty]) => {
                const ca = Math.cos(ang), sa = Math.sin(ang);
                const rot = (v: Vector): Vector => vec(
                    ca * v.get(0) - sa * v.get(1),
                    sa * v.get(0) + ca * v.get(1));
                const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty));
                const s2 = Segment.fromEndpoints(xf(s.p[0]), xf(s.p[1]));
                const b2 = OrientedBox.fromCenterAxisExtent(xf(b.center),
                    [rot(b.axis[0]), rot(b.axis[1])], b.extent);
                if (Math.abs(minSegDepth(b, s)) < 1e-6) {
                    return;   // grazing corner or edge
                }
                const f0 = fiq.find(s, b);
                const f1 = fiq.find(s2, b2);
                if (f0.intersect !== f1.intersect) {
                    return;
                }
                for (let i = 0; i < f0.numIntersections; ++i) {
                    expectVectorClose(f1.point[i], xf(f0.point[i]),
                        1e-7, 1e-8);
                }
            });
    });

    it('a degenerate segment follows the upstream oriented-box convention',
        () => {
            // Upstream reports numIntersections = 2 here while the
            // aligned-box sibling reports 1 for the same configuration; the
            // quirk is preserved.
            const box = rotatedBox(vec(0, 0), Math.PI / 6, vec(1, 2));
            const inside = segment([0.2, -0.3], [0.2, -0.3]);
            const f = fiq.find(inside, box);
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(2);
            expect(f.parameter).toEqual([0, 0]);
            expect(f.cdeParameter).toEqual([0, 0]);
            expectVectorClose(f.point[0], inside.p[0], 0, 0);
            expectVectorClose(f.point[1], inside.p[1], 0, 0);

            const outside = segment([9, 9], [9, 9]);
            const g = fiq.find(outside, box);
            expect(g.intersect).toBe(false);
            expect(g.numIntersections).toBe(0);
        });
});
