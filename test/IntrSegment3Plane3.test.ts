import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import {
    IntrSegment3Plane3TI,
    IntrSegment3Plane3FI,
    defaultIntrSegment3Plane3FIResult,
    defaultIntrSegment3Plane3TIResult,
    intrSegment3Plane3FIDoQuery
} from '../src/IntrSegment3Plane3.js';
import { IntrLine3Plane3FI } from '../src/IntrLine3Plane3.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame, unitVector,
    wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function plane(normal: number[], origin: number[]): Hyperplane {
    const n = Vector.fromArray(normal);
    normalize(n);
    return Hyperplane.fromNormalOrigin(n, Vector.fromArray(origin));
}

const ti = new IntrSegment3Plane3TI();
const fi = new IntrSegment3Plane3FI();

describe('IntrSegment3Plane3', () => {
    it('default-constructs results as no intersection', () => {
        expect(defaultIntrSegment3Plane3TIResult()).toEqual({ intersect: false });
        const r = defaultIntrSegment3Plane3FIResult();
        expect(r.intersect).toBe(false);
        expect(r.numIntersections).toBe(0);
    });

    it('finds a transverse intersection at a known point', () => {
        const S = segment([0, 0, -1], [0, 0, 3]);
        const P = plane([0, 0, 1], [0, 0, 1]);
        const result = fi.find(S, P);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        // The centered form has center (0,0,1), direction (0,0,1), extent 2,
        // so the intersection parameter is 0.
        expect(result.parameter).toBeCloseTo(0, 12);
        expect(result.point.values[2]).toBeCloseTo(1, 12);
        expect(ti.test(S, P).intersect).toBe(true);
    });

    it('rejects a segment entirely on one side of the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        for (const S of [segment([0, 0, 1], [0, 0, 5]),
            segment([0, 0, -5], [0, 0, -1])]) {
            expect(fi.find(S, P).intersect).toBe(false);
            expect(fi.find(S, P).numIntersections).toBe(0);
            expect(ti.test(S, P).intersect).toBe(false);
        }
    });

    it('accepts an endpoint touching the plane', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        const S0 = segment([1, 2, 0], [1, 2, 4]);
        expect(ti.test(S0, P).intersect).toBe(true);
        const r0 = fi.find(S0, P);
        expect(r0.intersect).toBe(true);
        expect(r0.point.values[2]).toBeCloseTo(0, 12);

        const S1 = segment([1, 2, -4], [1, 2, 0]);
        expect(ti.test(S1, P).intersect).toBe(true);
        expect(fi.find(S1, P).intersect).toBe(true);
    });

    it('reports a segment lying in the plane as infinitely many hits', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        const S = segment([-1, -1, 0], [2, 3, 0]);
        expect(ti.test(S, P).intersect).toBe(true);
        const r = fi.find(S, P);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(2147483647);
    });

    it('reports a parallel disjoint segment as no intersection', () => {
        const P = plane([0, 0, 1], [0, 0, 2]);
        const S = segment([-1, -1, 0], [2, 3, 0]);
        expect(ti.test(S, P).intersect).toBe(false);
        expect(fi.find(S, P).intersect).toBe(false);
    });

    it('exposes the DoQuery helper used by derived queries', () => {
        const P = plane([0, 0, 1], [0, 0, 0]);
        // The line hits the plane at parameter 5, outside the extent 2.
        const outside = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(vec(0, 0, -5), vec(0, 0, 1), 2, P, outside);
        expect(outside.intersect).toBe(false);
        expect(outside.numIntersections).toBe(0);

        const inside = defaultIntrSegment3Plane3FIResult();
        intrSegment3Plane3FIDoQuery(vec(0, 0, -1), vec(0, 0, 1), 2, P, inside);
        expect(inside.intersect).toBe(true);
        expect(inside.parameter).toBeCloseTo(1, 12);
    });

    it('agrees with TI and puts the point on both primitives (randomized)', () => {
        let seed = 24680135;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        const rnd = (): number => 4 * rand() - 2;

        let hits = 0;
        for (let trial = 0; trial < 4000; ++trial) {
            const p0 = [rnd(), rnd(), rnd()];
            const p1 = [rnd(), rnd(), rnd()];
            const S = segment(p0, p1);
            const P = plane([rnd(), rnd(), rnd()], [rnd(), rnd(), rnd()]);
            const fiResult = fi.find(S, P);
            const tiResult = ti.test(S, P);
            expect(fiResult.intersect).toBe(tiResult.intersect);
            if (fiResult.intersect && fiResult.numIntersections === 1) {
                ++hits;
                // The point is on the plane.
                expect(dot(P.normal, fiResult.point) - P.constant)
                    .toBeCloseTo(0, 8);
                // The point is on the segment: X = p0 + s * (p1 - p0) with
                // s in [0,1].
                const d = sub(S.p[1], S.p[0]);
                const s = dot(sub(fiResult.point, S.p[0]), d) / dot(d, d);
                expect(s).toBeGreaterThanOrEqual(-1e-9);
                expect(s).toBeLessThanOrEqual(1 + 1e-9);
                const onSeg = add(S.p[0], mul(s, d));
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.point.values[i]).toBeCloseTo(
                        onSeg.values[i], 8);
                }
            }
        }
        expect(hits).toBeGreaterThan(500);
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): property-based cross-checks against the line query and
// against the geometry of the reported point.
// ---------------------------------------------------------------------------

describe('IntrSegment3Plane3 verification', () => {
    const tiv = new IntrSegment3Plane3TI();
    const fiv = new IntrSegment3Plane3FI();
    const lineFI = new IntrLine3Plane3FI();
    const INT32_MAX = 2147483647;

    const segPlane = fc.tuple(unitVector(3), wellScaledVector(3, -5, 5),
        wellScaledVector(3, -5, 5), wellScaledVector(3, -5, 5))
        .filter(([, , p0, p1]) => length(sub(p1, p0)) > 1e-2)
        .map(([n, po, p0, p1]) => ({
            plane: Hyperplane.fromNormalOrigin(n, po),
            segment: Segment.fromEndpoints(p0, p1)
        }));

    it('TI and FI agree on intersect away from grazing configurations', () => {
        check(segPlane, ({ segment: S, plane: P }) => {
            const sd0 = dot(P.normal, S.p[0]) - P.constant;
            const sd1 = dot(P.normal, S.p[1]) - P.constant;
            // TI decides from the endpoint signs; FI decides from
            // |t| <= extent in the centered form. The two are the same
            // predicate in exact arithmetic, but the FI comparison is a knife
            // edge when the crossing is at an endpoint, so require the
            // endpoints to be off the plane by a relative margin.
            const scale = Math.max(Math.abs(sd0), Math.abs(sd1), 1);
            if (Math.min(Math.abs(sd0), Math.abs(sd1)) < 1e-8 * scale) {
                return;
            }
            expect(fiv.find(S, P).intersect).toBe(tiv.test(S, P).intersect);
        });
    });

    it('equals the line query restricted to |t| <= extent', () => {
        check(segPlane, ({ segment: S, plane: P }) => {
            const cf = S.getCenteredForm();
            const lr = lineFI.find(
                Line.fromOriginDirection(cf.center, cf.direction), P);
            const sr = fiv.find(S, P);
            const expected = lr.intersect
                && !(Math.abs(lr.parameter) > cf.extent);
            expect(sr.intersect).toBe(expected);
            if (expected) {
                expect(sr.numIntersections).toBe(lr.numIntersections);
                expect(sr.parameter).toBe(lr.parameter);
            }
            else {
                expect(sr.numIntersections).toBe(0);
            }
        });
    });

    it('reports a point on the plane and inside the segment', () => {
        check(segPlane, ({ segment: S, plane: P }) => {
            const r = fiv.find(S, P);
            if (!r.intersect) {
                return;
            }
            const cf = S.getCenteredForm();
            for (let i = 0; i < 3; ++i) {
                expect(Number.isNaN(r.point.get(i))).toBe(false);
            }
            if (r.numIntersections === INT32_MAX) {
                // The segment lies on the plane; the reported point is the
                // segment center.
                expect(r.parameter).toBe(0);
                expectVectorClose(r.point, cf.center);
                return;
            }
            expect(r.numIntersections).toBe(1);
            // The result parameter is in the centered form, as documented.
            expect(Math.abs(r.parameter) <= cf.extent).toBe(true);
            expectVectorClose(r.point,
                add(cf.center, mul(r.parameter, cf.direction)), 1e-12, 1e-12);
            const scale = 1 + cf.extent;
            expectClose(dot(P.normal, r.point) - P.constant, 0,
                1e-12 * scale, 1e-12);
        });
    });

    it('finds the crossing of a segment built to straddle the plane', () => {
        check(fc.tuple(unitVector(3), wellScaledVector(3, -5, 5),
            wellScaledVector(3, -5, 5),
            fc.double({ min: 0.2, max: 5, noNaN: true }),
            fc.double({ min: 0.2, max: 5, noNaN: true })),
            ([n, po, tangent, a, b]) => {
                // p0 is at signed distance -a and p1 at +b from the plane, so
                // the segment crosses the plane transversely.
                const base = add(po, sub(tangent, mul(dot(n, tangent), n)));
                const p0 = add(base, mul(-a, n));
                const p1 = add(base, mul(+b, n));
                const P = Hyperplane.fromNormalOrigin(n, po);
                const S = Segment.fromEndpoints(p0, p1);
                expect(tiv.test(S, P).intersect).toBe(true);
                const r = fiv.find(S, P);
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(1);
                // The crossing is at the parameter a/(a+b) of [p0,p1].
                const expected = add(p0, mul(a / (a + b), sub(p1, p0)));
                expectVectorClose(r.point, expected, 1e-9, 1e-9);
            });
    });

    it('reports the whole segment when it lies exactly in the plane', () => {
        // The segment must lie in the plane exactly, not merely to within
        // rounding: for a segment that is only nearly coplanar the FI query
        // divides a tiny signed distance by a tiny Dot(D,N) and the resulting
        // parameter is arbitrary, so the classification is a knife edge. An
        // axis-aligned plane and endpoints sharing that coordinate give exact
        // zeros.
        check(fc.tuple(fc.integer({ min: 0, max: 2 }),
            wellScaledVector(3, -5, 5), wellScaledVector(3, -5, 5),
            wellScaledVector(3, -5, 5)),
            ([k, po, a, b]) => {
                const p0 = a.clone();
                const p1 = b.clone();
                p0.values[k] = po.values[k];
                p1.values[k] = po.values[k];
                if (length(sub(p1, p0)) < 1e-2) {
                    return;
                }
                const P = Hyperplane.fromNormalOrigin(Vector.unit(3, k), po);
                const S = Segment.fromEndpoints(p0, p1);
                expect(tiv.test(S, P).intersect).toBe(true);
                const r = fiv.find(S, P);
                expect(r.intersect).toBe(true);
                expect(r.numIntersections).toBe(INT32_MAX);
                expect(r.parameter).toBe(0);
                expectVectorClose(r.point, S.getCenteredForm().center);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(segPlane, rotationFrame(3), wellScaledVector(3, -5, 5)),
            ([{ segment: S, plane: P }, frame, t]) => {
                const xfDir = (v: Vector): Vector => {
                    const w = new Vector(3);
                    for (let i = 0; i < 3; ++i) {
                        w.values[i] = frame[0].values[i] * v.values[0]
                            + frame[1].values[i] * v.values[1]
                            + frame[2].values[i] * v.values[2];
                    }
                    return w;
                };
                const xf = (v: Vector): Vector => add(xfDir(v), t);
                const sd0 = dot(P.normal, S.p[0]) - P.constant;
                const sd1 = dot(P.normal, S.p[1]) - P.constant;
                const scale = Math.max(Math.abs(sd0), Math.abs(sd1), 1);
                if (Math.min(Math.abs(sd0), Math.abs(sd1)) < 1e-6 * scale) {
                    return;
                }
                const r0 = fiv.find(S, P);
                const r1 = fiv.find(
                    Segment.fromEndpoints(xf(S.p[0]), xf(S.p[1])),
                    Hyperplane.fromNormalOrigin(xfDir(P.normal),
                        xf(mul(P.constant, P.normal))));
                expect(r1.intersect).toBe(r0.intersect);
                if (r0.intersect && r0.numIntersections === 1) {
                    expectVectorClose(r1.point, xf(r0.point), 1e-8, 1e-8);
                }
            });
    });
});
