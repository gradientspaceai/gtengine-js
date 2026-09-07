import { describe, it, expect } from 'vitest';
import { Cylinder3 } from '../src/Cylinder3.js';
import {
    IntrSegment3Cylinder3FI,
    defaultIntrSegment3Cylinder3FIResult,
    intrSegment3Cylinder3FIDoQuery
} from '../src/IntrSegment3Cylinder3.js';
import { IntrLine3Cylinder3FI } from '../src/IntrLine3Cylinder3.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function cylinder(origin: number[], direction: number[], radius: number,
    height: number): Cylinder3 {
    const d = Vector.fromArray(direction);
    normalize(d);
    return Cylinder3.fromAxisRadiusHeight(
        Line.fromOriginDirection(Vector.fromArray(origin), d), radius, height);
}

// The signed "outside" amount of a point relative to the solid cylinder:
// negative strictly inside, zero on the boundary, positive outside.
function cylinderSignedDepth(p: Vector, c: Cylinder3): number {
    const diff = sub(p, c.axis.origin);
    const z = dot(diff, c.axis.direction);
    const radial = sub(diff, mul(z, c.axis.direction));
    const rDepth = Math.sqrt(dot(radial, radial)) - c.radius;
    const zDepth = Math.abs(z) - 0.5 * c.height;
    return Math.max(rDepth, zDepth);
}

// The most deeply contained sampled point of the segment.
function bestSampledDepth(s: Segment, c: Cylinder3): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= 4000; ++i) {
        const p = add(s.p[0], mul(i / 4000, sub(s.p[1], s.p[0])));
        const value = cylinderSignedDepth(p, c);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrSegment3Cylinder3', () => {
    const fi = new IntrSegment3Cylinder3FI();
    // A cylinder about the z axis: radius 1, height 2, so |z| <= 1.
    const unit = cylinder([0, 0, 0], [0, 0, 1], 1, 2);

    it('default results report no intersection', () => {
        const result = defaultIntrSegment3Cylinder3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('clips a segment crossing the cylinder wall', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The parameters are relative to the centered form, extent 3.
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('clips a segment crossing both end caps', () => {
        const s = segment([0, 0, -3], [0, 0, 3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.point[0].values[2]).toBeCloseTo(-1, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports the whole segment when it lies inside the cylinder', () => {
        const s = segment([-0.5, 0, 0], [0.5, 0, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-0.5, 12);
        expect(result.parameter[1]).toBeCloseTo(0.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-0.5, 12);
        expect(result.point[1].values[0]).toBeCloseTo(0.5, 12);
    });

    it('clips a segment whose endpoints straddle an end cap', () => {
        const s = segment([0, 0, 0], [0, 0, 4]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The centered form has center (0,0,2), extent 2, so the overlap of
        // the line interval with [-2,2] is [-2,-1].
        expect(result.parameter[0]).toBeCloseTo(-2, 12);
        expect(result.parameter[1]).toBeCloseTo(-1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(0, 12);
        expect(result.point[1].values[2]).toBeCloseTo(1, 12);
    });

    it('reports a tangent contact on the cylinder wall', () => {
        const s = segment([1, -3, 0], [1, 3, 0]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(1, 12);
    });

    it('reports a tangent contact on the cap rim', () => {
        const s = segment([1, -3, 1], [1, 3, 1]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.point[0].values[0]).toBeCloseTo(1, 12);
        expect(result.point[0].values[2]).toBeCloseTo(1, 12);
    });

    it('reports no intersection for a segment that misses the cylinder',
        () => {
            const s = segment([5, 0, 0], [5, 0, 3]);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('reports no intersection when the line hits but the segment does not',
        () => {
            // The line x = 0, y = 0 passes through the cylinder, but the
            // segment lies entirely above the top cap.
            const s = segment([0, 0, 3], [0, 0, 5]);
            expect(fi.find(s, unit).intersect).toBe(false);
        });

    it('clips a diagonal segment symmetrically', () => {
        const s = segment([-3, -3, -3], [3, 3, 3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        // The direction is (1,1,1)/sqrt(3); the wall x^2+y^2 = 1 is met at
        // x = y = +-1/sqrt(2), which is inside |z| <= 1.
        const h = Math.SQRT1_2;
        for (let i = 0; i < 3; ++i) {
            expect(result.point[0].values[i]).toBeCloseTo(-h, 10);
            expect(result.point[1].values[i]).toBeCloseTo(h, 10);
        }
        expect(result.parameter[0]).toBeCloseTo(-h * Math.sqrt(3), 10);
        expect(result.parameter[1]).toBeCloseTo(h * Math.sqrt(3), 10);
    });

    it('handles a tilted, translated cylinder', () => {
        const c = cylinder([1, 2, 3], [1, 1, 0], 0.5, 4);
        // A segment along the cylinder axis, longer than the cylinder.
        const h = Math.SQRT1_2;
        const s = segment([1 - 5 * h, 2 - 5 * h, 3], [1 + 5 * h, 2 + 5 * h, 3]);
        const result = fi.find(s, c);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-2, 10);
        expect(result.parameter[1]).toBeCloseTo(2, 10);
        expect(cylinderSignedDepth(result.point[0], c)).toBeCloseTo(0, 9);
        expect(cylinderSignedDepth(result.point[1], c)).toBeCloseTo(0, 9);
    });

    it('rejects a zero-length segment outside the cylinder', () => {
        const s = segment([5, 5, 5], [5, 5, 5]);
        expect(fi.find(s, unit).intersect).toBe(false);
    });

    it('reports a degenerate zero-length segment inside the cylinder', () => {
        // A zero-length segment has no well-defined centered-form direction
        // (Normalize of the zero vector leaves it zero), so the quadratic in
        // the line-cylinder DoQuery is identically zero and the reported
        // parameters are NaN. This matches the upstream C++ behavior; see
        // the port notes in IntrSegment3Cylinder3.ts.
        const s = segment([0.2, 0.1, 0.3], [0.2, 0.1, 0.3]);
        const result = fi.find(s, unit);
        expect(result.intersect).toBe(true);
        expect(Number.isNaN(result.parameter[0])).toBe(true);
        expect(Number.isNaN(result.parameter[1])).toBe(true);
    });

    it('rejects an infinite cylinder', () => {
        const infinite = cylinder([0, 0, 0], [0, 0, 1], 1, 2);
        infinite.makeInfiniteCylinder();
        expect(() => fi.find(segment([-3, 0, 0], [3, 0, 0]), infinite))
            .toThrow('Infinite cylinders are not yet supported.');
    });

    it('exposes the DoQuery helper, which matches the public query', () => {
        const s = segment([-3, 0, 0.5], [3, 0, 0.5]);
        const { center, direction, extent } = s.getCenteredForm();
        const result = defaultIntrSegment3Cylinder3FIResult();
        intrSegment3Cylinder3FIDoQuery(center, direction, extent, unit,
            result);
        const expected = fi.find(s, unit);
        expect(result.intersect).toBe(expected.intersect);
        expect(result.numIntersections).toBe(expected.numIntersections);
        expect(result.parameter).toEqual(expected.parameter);
        // DoQuery does not compute the points; the public query does.
        expect(result.point[0].values[0]).toBe(0);
        expect(expected.point[0].values[0]).toBeCloseTo(-1, 12);
    });

    it('restricts the line query to the segment t-interval', () => {
        // Whenever the segment query reports an intersection, the line query
        // must report one too, and the segment parameters must lie inside
        // both the line interval and [-extent, extent].
        const lineQuery = new IntrLine3Cylinder3FI();
        const random = makeRandom(90941);
        let intersections = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const c = cylinder(
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                0.2 + random(), 0.4 + 2 * random());
            const s = segment(
                [4 * random() - 2, 4 * random() - 2, 4 * random() - 2],
                [4 * random() - 2, 4 * random() - 2, 4 * random() - 2]);
            const { center, direction, extent } = s.getCenteredForm();
            const segResult = fi.find(s, c);
            const lineResult = lineQuery.find(
                Line.fromOriginDirection(center, direction), c);
            if (segResult.intersect) {
                ++intersections;
                expect(lineResult.intersect).toBe(true);
                for (let i = 0; i < 2; ++i) {
                    expect(segResult.parameter[i])
                        .toBeGreaterThanOrEqual(
                            Math.max(lineResult.parameter[0], -extent) - 1e-9);
                    expect(segResult.parameter[i])
                        .toBeLessThanOrEqual(
                            Math.min(lineResult.parameter[1], extent) + 1e-9);
                }
            }
        }
        expect(intersections).toBeGreaterThan(20);
    });

    it('agrees with sampling on random configurations', () => {
        const random = makeRandom(90942);
        let intersections = 0;
        for (let trial = 0; trial < 400; ++trial) {
            const c = cylinder(
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                [2 * random() - 1, 2 * random() - 1, 2 * random() - 1],
                0.2 + random(), 0.4 + 2 * random());
            const s = segment(
                [5 * random() - 2.5, 5 * random() - 2.5, 5 * random() - 2.5],
                [5 * random() - 2.5, 5 * random() - 2.5, 5 * random() - 2.5]);

            const result = fi.find(s, c);
            const depth = bestSampledDepth(s, c);
            if (depth < -1e-6) {
                expect(result.intersect).toBe(true);
            }
            else if (depth > 1e-2) {
                expect(result.intersect).toBe(false);
            }

            if (result.intersect) {
                ++intersections;
                expect(result.parameter[0])
                    .toBeLessThanOrEqual(result.parameter[1] + 1e-12);
                const { center, direction, extent } = s.getCenteredForm();
                for (let i = 0; i < 2; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(extent + 1e-9);
                    const onSegment = add(center,
                        mul(result.parameter[i], direction));
                    for (let j = 0; j < 3; ++j) {
                        expect(result.point[i].values[j])
                            .toBeCloseTo(onSegment.values[j], 9);
                    }
                    expect(cylinderSignedDepth(result.point[i], c))
                        .toBeLessThan(1e-9);
                }
            }
        }
        expect(intersections).toBeGreaterThan(20);
    });
});

// ---------------------------------------------------------------------------
// Verification (V33): property-based cross-checks against upstream
// IntrSegment3Cylinder3.h.
// ---------------------------------------------------------------------------

import {
    check, fc, expectClose, expectVectorClose, seededRandom, wellScaled
} from './helpers/arbitraries.js';
import { length } from '../src/Vector.js';

// wellScaled snaps |a| < 1e-3 to exactly zero, so no direction component is a
// subnormal that would underflow when squared.
const angle3 = () => wellScaled(-Math.PI, Math.PI);

function unitDir3(th: number, ph: number): Vector {
    const d = vec(Math.cos(th) * Math.cos(ph), Math.sin(th) * Math.cos(ph),
        Math.sin(ph));
    normalize(d);
    return d;
}

const segCylinder = fc.tuple(
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-5, 5), wellScaled(-5, 5), wellScaled(-5, 5),
    wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2),
    angle3(), angle3(),
    fc.double({ min: 0.25, max: 2, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.5, max: 5, noNaN: true, noDefaultInfinity: true })
).filter(([x0, y0, z0, x1, y1, z1]) =>
    (x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2 > 1e-2)
    .map(([x0, y0, z0, x1, y1, z1, cx, cy, cz, ath, aph, r, h]) => ({
        segment: Segment.fromEndpoints(vec(x0, y0, z0), vec(x1, y1, z1)),
        cylinder: Cylinder3.fromAxisRadiusHeight(
            Line.fromOriginDirection(vec(cx, cy, cz), unitDir3(ath, aph)),
            r, h)
    }));

// cylinderSignedDepth is a maximum of convex functions, hence convex along
// the segment; a ternary search finds its minimum reliably.
function minSegCylDepth(c: Cylinder3, s: Segment): number {
    const e = sub(s.p[1], s.p[0]);
    let lo = 0, hi = 1;
    for (let k = 0; k < 200; ++k) {
        const p = lo + (hi - lo) / 3, q = hi - (hi - lo) / 3;
        const fp = cylinderSignedDepth(add(s.p[0], mul(p, e)), c);
        const fq = cylinderSignedDepth(add(s.p[0], mul(q, e)), c);
        if (fp < fq) { hi = q; } else { lo = p; }
    }
    return cylinderSignedDepth(add(s.p[0], mul(0.5 * (lo + hi), e)), c);
}

describe('IntrSegment3Cylinder3 verification', () => {
    const fiq = new IntrSegment3Cylinder3FI();
    const lcq = new IntrLine3Cylinder3FI();

    it('the segment hit is the line hit clipped to |t| <= extent', () => {
        check(segCylinder, ({ segment: s, cylinder: c }) => {
            const cf = s.getCenteredForm();
            const line = Line.fromOriginDirection(cf.center, cf.direction);
            const lc = lcq.find(line, c);
            const f = fiq.find(s, c);
            if (!lc.intersect) {
                expect(f.intersect).toBe(false);
                return;
            }
            const t0 = Math.max(lc.parameter[0], -cf.extent);
            const t1 = Math.min(lc.parameter[1], cf.extent);
            if (t0 > t1) {
                expect(f.intersect).toBe(false);
                return;
            }
            expect(f.intersect).toBe(true);
            expect(f.numIntersections).toBe(t0 < t1 ? 2 : 1);
            // Normalize the -0/+0 tie (toBe uses Object.is).
            expect(f.parameter[0] + 0).toBe(t0 + 0);
        });
    });

    it('reported points lie on the segment and on the cylinder boundary',
        () => {
            check(segCylinder, ({ segment: s, cylinder: c }) => {
                const cf = s.getCenteredForm();
                const f = fiq.find(s, c);
                if (!f.intersect) {
                    expect(f.numIntersections).toBe(0);
                    expect(f.point[0].values).toEqual([0, 0, 0]);
                    expect(f.point[1].values).toEqual([0, 0, 0]);
                    return;
                }
                for (let i = 0; i < 2; ++i) {
                    // The parameter is in the centered convention, C + t*D
                    // with |t| <= e.
                    expect(Math.abs(f.parameter[i]))
                        .toBeLessThanOrEqual(cf.extent + 1e-12);
                    expectVectorClose(f.point[i],
                        add(cf.center, mul(f.parameter[i], cf.direction)),
                        0, 0);
                    const depth = cylinderSignedDepth(f.point[i], c);
                    if (Math.abs(Math.abs(f.parameter[i]) - cf.extent)
                        > 1e-12) {
                        expectClose(depth, 0, 1e-7, 1e-8);
                    }
                    else {
                        expect(depth).toBeLessThanOrEqual(1e-7);
                    }
                }
            });
        });

    it('a segment strictly inside the cylinder is clipped to its endpoints',
        () => {
            check(fc.tuple(wellScaled(-2, 2), wellScaled(-2, 2),
                wellScaled(-2, 2), angle3(), angle3(),
                fc.double({ min: 0.25, max: 2, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: 0.5, max: 5, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true }),
                fc.double({ min: -0.8, max: 0.8, noNaN: true,
                    noDefaultInfinity: true })),
            ([cx, cy, cz, ath, aph, r, h, u0, u1]) => {
                const axis = unitDir3(ath, aph);
                const c = Cylinder3.fromAxisRadiusHeight(
                    Line.fromOriginDirection(vec(cx, cy, cz), axis), r, h);
                // Two points on the axis, strictly between the end disks.
                const p0 = add(c.axis.origin, mul(u0 * 0.5 * h, axis));
                const p1 = add(c.axis.origin, mul(u1 * 0.5 * h, axis));
                if (length(sub(p1, p0)) < 1e-3) { return; }
                const s = Segment.fromEndpoints(p0, p1);
                const cf = s.getCenteredForm();
                const f = fiq.find(s, c);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(2);
                expectClose(f.parameter[0], -cf.extent, 1e-9, 1e-9);
                expectClose(f.parameter[1], cf.extent, 1e-9, 1e-9);
                expectVectorClose(f.point[0], p0, 1e-8, 1e-9);
                expectVectorClose(f.point[1], p1, 1e-8, 1e-9);
            });
        });

    it('a fine sweep of the segment agrees with the reported interval', () => {
        const rnd = seededRandom(0x3b7e12c);
        for (let trial = 0; trial < 150; ++trial) {
            const axis = vec(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (dot(axis, axis) < 1e-4) { continue; }
            normalize(axis);
            const c = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(
                    vec(rnd() * 3 - 1.5, rnd() * 3 - 1.5, rnd() * 3 - 1.5),
                    axis),
                0.3 + rnd() * 1.5, 0.5 + rnd() * 4);
            const s = Segment.fromEndpoints(
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5),
                vec(rnd() * 10 - 5, rnd() * 10 - 5, rnd() * 10 - 5));
            if (length(sub(s.p[1], s.p[0])) < 1e-2) { continue; }
            const cf = s.getCenteredForm();
            const f = fiq.find(s, c);
            for (let k = 0; k <= 500; ++k) {
                const t = -cf.extent + (2 * cf.extent * k) / 500;
                const p = add(cf.center, mul(t, cf.direction));
                if (cylinderSignedDepth(p, c) < -1e-7) {
                    expect(f.intersect).toBe(true);
                    expect(t).toBeGreaterThanOrEqual(f.parameter[0] - 1e-7);
                    expect(t).toBeLessThanOrEqual(f.parameter[1] + 1e-7);
                }
            }
        }
    }, 30000);

    it('the exported DoQuery reproduces the class result', () => {
        check(segCylinder, ({ segment: s, cylinder: c }) => {
            const cf = s.getCenteredForm();
            const result = defaultIntrSegment3Cylinder3FIResult();
            intrSegment3Cylinder3FIDoQuery(cf.center, cf.direction, cf.extent,
                c, result);
            const f = fiq.find(s, c);
            expect(result.intersect).toBe(f.intersect);
            expect(result.numIntersections).toBe(f.numIntersections);
            expect(result.parameter).toEqual(f.parameter);
            expect(result.point[0].values).toEqual([0, 0, 0]);
            expect(result.point[1].values).toEqual([0, 0, 0]);
        });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(segCylinder, angle3(), angle3(), angle3(),
            wellScaled(-2, 2), wellScaled(-2, 2), wellScaled(-2, 2)),
        ([{ segment: s, cylinder: c }, a1, a2, a3, tx, ty, tz]) => {
            if (Math.abs(minSegCylDepth(c, s)) < 1e-6) {
                return;   // grazing the wall or a rim
            }
            const ca = Math.cos(a1), sa = Math.sin(a1);
            const cb = Math.cos(a2), sb = Math.sin(a2);
            const cc = Math.cos(a3), sc = Math.sin(a3);
            const f0v = vec(ca * cb, sa * cb, -sb);
            const f1v = vec(ca * sb * sc - sa * cc, sa * sb * sc + ca * cc,
                cb * sc);
            const f2v = vec(ca * sb * cc + sa * sc, sa * sb * cc - ca * sc,
                cb * cc);
            const rot = (v: Vector): Vector => add(mul(v.get(0), f0v),
                add(mul(v.get(1), f1v), mul(v.get(2), f2v)));
            const xf = (v: Vector): Vector => add(rot(v), vec(tx, ty, tz));
            const s2 = Segment.fromEndpoints(xf(s.p[0]), xf(s.p[1]));
            const c2 = Cylinder3.fromAxisRadiusHeight(
                Line.fromOriginDirection(xf(c.axis.origin),
                    rot(c.axis.direction)), c.radius, c.height);
            const g0 = fiq.find(s, c);
            const g1 = fiq.find(s2, c2);
            expect(g1.intersect).toBe(g0.intersect);
            if (!g0.intersect) { return; }
            expect(g1.numIntersections).toBe(g0.numIntersections);
            for (let i = 0; i < 2; ++i) {
                expectVectorClose(g1.point[i], xf(g0.point[i]), 1e-6, 1e-7);
            }
        });
    });

    it('rejects an infinite cylinder', () => {
        // Port deviation, documented in the source: upstream reads
        // cylinder.height directly and would compute nonsense for the
        // height = -1 infinite sentinel.
        const c = cylinder([0, 0, 0], [0, 0, 1], 1, 4);
        c.makeInfiniteCylinder();
        expect(c.isFinite()).toBe(false);
        expect(() => fiq.find(segment([-5, 0, 0], [5, 0, 0]), c)).toThrow(
            'Infinite cylinders are not yet supported.');
    });
});
