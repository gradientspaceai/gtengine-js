import { describe, it, expect } from 'vitest';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import {
    IntrSegment3Ellipsoid3TI,
    IntrSegment3Ellipsoid3FI,
    defaultIntrSegment3Ellipsoid3FIResult,
    intrSegment3Ellipsoid3FIDoQuery
} from '../src/IntrSegment3Ellipsoid3.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { mulMatrix } from '../src/Matrix.js';
import { intrLine3Ellipsoid3FIDoQuery } from '../src/IntrLine3Ellipsoid3.js';
import {
    check, expectClose, expectVectorClose, fc, positive, rotationFrame,
    segment as arbSegment, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

function ellipsoid(center: number[], extent: number[]): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(Vector.fromArray(center),
        [vec(1, 0, 0), vec(0, 1, 0), vec(0, 0, 1)],
        Vector.fromArray(extent));
}

// The value of (X-C)^T M (X-C) - 1: zero on the ellipsoid, negative inside.
function level(E: Hyperellipsoid, X: Vector): number {
    const d = sub(X, E.center);
    let sum = -1;
    for (let i = 0; i < 3; ++i) {
        const t = d.values[i] / E.extent.values[i];
        sum += t * t;
    }
    return sum;
}

// An independent test: sample the segment densely and look for a point in the
// solid ellipsoid.
function bruteForceIntersect(S: Segment, E: Hyperellipsoid): boolean {
    const delta = sub(S.p[1], S.p[0]);
    const n = 4000;
    for (let i = 0; i <= n; ++i) {
        if (level(E, add(S.p[0], mul(i / n, delta))) <= 0) {
            return true;
        }
    }
    return false;
}

const ti = new IntrSegment3Ellipsoid3TI();
const fi = new IntrSegment3Ellipsoid3FI();

describe('IntrSegment3Ellipsoid3', () => {
    it('has an empty default result', () => {
        const result = defaultIntrSegment3Ellipsoid3FIResult();
        expect(result.intersect).toBe(false);
        expect(result.numIntersections).toBe(0);
        expect(result.parameter).toEqual([0, 0]);
    });

    it('finds both crossings of a straddling segment', () => {
        // The segment spans x in [-10,10] at y = z = 0; its centered form has
        // origin (0,0,0), direction (1,0,0) and extent 10, so the crossings
        // are at parameters -3 and 3.
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-10, 0, 0], [10, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-3, 12);
        expect(result.parameter[1]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(-3, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('clips a crossing when an endpoint is inside the ellipsoid', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        // The segment from the center to (10,0,0): centered form origin
        // (5,0,0), direction (1,0,0), extent 5, so the ellipsoid crossing at
        // x = 3 is at parameter -2 and the clipped near end is at -5.
        const S = seg([0, 0, 0], [10, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-5, 12);
        expect(result.parameter[1]).toBeCloseTo(-2, 12);
        expect(result.point[0].values[0]).toBeCloseTo(0, 12);
        expect(result.point[1].values[0]).toBeCloseTo(3, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('reports no intersection for a segment that stops short', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-10, 0, 0], [-4, 0, 0]);
        expect(fi.find(S, E).intersect).toBe(false);
        expect(ti.test(S, E).intersect).toBe(false);
        // A segment whose supporting line misses the ellipsoid.
        const M = seg([-10, 5, 0], [10, 5, 0]);
        expect(fi.find(M, E).intersect).toBe(false);
        expect(ti.test(M, E).intersect).toBe(false);
    });

    it('reports a segment fully inside the ellipsoid', () => {
        const E = ellipsoid([0, 0, 0], [3, 2, 1]);
        const S = seg([-1, 0, 0], [1, 0, 0]);
        const result = fi.find(S, E);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-1, 12);
        expect(result.parameter[1]).toBeCloseTo(1, 12);
        expect(ti.test(S, E).intersect).toBe(true);
    });

    it('exposes the DoQuery helper without computing points', () => {
        const result = defaultIntrSegment3Ellipsoid3FIResult();
        intrSegment3Ellipsoid3FIDoQuery(vec(0, 0, 0), vec(1, 0, 0), 10,
            ellipsoid([0, 0, 0], [3, 2, 1]), result);
        expect(result.numIntersections).toBe(2);
        expect(result.parameter[0]).toBeCloseTo(-3, 12);
        expect(result.point[0].values).toEqual([0, 0, 0]);
    });

    it('agrees with the TI query and a dense sampling on random inputs', () => {
        let state = 24680;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        let numHits = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const E = ellipsoid([rand(), rand(), rand()],
                [0.5 + Math.abs(rand()) * 2, 0.5 + Math.abs(rand()) * 2,
                    0.5 + Math.abs(rand()) * 2]);
            const S = seg([rand() * 4, rand() * 4, rand() * 4],
                [rand() * 4, rand() * 4, rand() * 4]);
            const result = fi.find(S, E);
            expect(ti.test(S, E).intersect).toBe(result.intersect);
            expect(result.intersect).toBe(bruteForceIntersect(S, E));

            if (result.intersect) {
                ++numHits;
                const { center: segOrigin, direction: segDirection,
                    extent: segExtent } = S.getCenteredForm();
                for (let i = 0; i < result.numIntersections; ++i) {
                    expect(Math.abs(result.parameter[i]))
                        .toBeLessThanOrEqual(segExtent + 1e-12);
                    const onSegment = add(segOrigin,
                        mul(result.parameter[i], segDirection));
                    expect(length(sub(result.point[i], onSegment)))
                        .toBeCloseTo(0, 10);
                    // The point is on or inside the ellipsoid; it is on the
                    // boundary unless it is a clipped segment endpoint.
                    expect(level(E, result.point[i])).toBeLessThan(1e-8);
                    const clipped =
                        Math.abs(Math.abs(result.parameter[i]) - segExtent)
                        < 1e-14;
                    if (!clipped) {
                        expect(level(E, result.point[i])).toBeCloseTo(0, 8);
                    }
                }
            }
        }
        expect(numHits).toBeGreaterThan(10);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrSegment3Ellipsoid3.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrSegment3Ellipsoid3 verification', () => {
    const tiQ = new IntrSegment3Ellipsoid3TI();
    const fiQ = new IntrSegment3Ellipsoid3FI();

    const arbEllipsoid = fc.tuple(wellScaledVector(3), rotationFrame(3),
        fc.tuple(positive(4, 0.05), positive(4, 0.05), positive(4, 0.05)))
        .map(([c, axis, e]) => Hyperellipsoid.fromCenterAxisExtent(c, axis,
            Vector.fromArray([e[0], e[1], e[2]])));

    function quadratic(E: Hyperellipsoid, X: Vector): number {
        const d = sub(X, E.center);
        return dot(d, mulMatrix(E.getM(), d) as Vector) - 1;
    }

    it('TI and FI agree on intersect', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            expect(tiQ.test(S, E).intersect).toBe(fiQ.find(S, E).intersect);
        });
    });

    it('agrees with a dense sampling of the segment against the solid', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            let minQ = Number.POSITIVE_INFINITY;
            for (let i = 0; i <= 2000; ++i) {
                const t = i / 2000;
                const X = add(mul(1 - t, S.p[0]), mul(t, S.p[1]));
                minQ = Math.min(minQ, quadratic(E, X));
            }
            const got = tiQ.test(S, E).intersect;
            // A sampled interior point proves an intersection with the solid.
            if (minQ < -1e-9) { expect(got).toBe(true); }
            // A comfortable positive minimum proves separation.
            if (minQ > 1e-3) { expect(got).toBe(false); }
        }, 60);
    });

    it('reports a segment fully inside the ellipsoid (upstream defect fixed)', () => {
        // Upstream's TIQuery tests 'qm * qp <= 0' and then concludes with
        // 'qm > 0 && |a1| < a2*e', so a segment strictly inside the solid
        // ellipsoid (both endpoint values negative) was reported as no
        // intersection, contradicting the FI query in the same header
        // (gtengine-js #304 item 2). The port tests 'qm <= 0 || qp <= 0'.
        check(fc.tuple(arbEllipsoid, unitVector(3), unitVector(3),
            fc.double({ min: 0, max: 0.6, noNaN: true }),
            fc.double({ min: 0, max: 0.6, noNaN: true })),
            ([E, u, v, f0, f1]) => {
                const interior = (f: number, w: Vector) => add(E.center, add(
                    mul(f * E.extent.values[0] * w.values[0], E.axis[0]),
                    add(mul(f * E.extent.values[1] * w.values[1], E.axis[1]),
                        mul(f * E.extent.values[2] * w.values[2], E.axis[2]))));
                const p0 = interior(f0, u);
                const p1 = interior(f1, v);
                // A degenerate segment has no centered-form direction.
                if (length(sub(p1, p0)) < 1e-6) { return; }
                expect(quadratic(E, p0)).toBeLessThan(0);
                expect(quadratic(E, p1)).toBeLessThan(0);
                const S = Segment.fromEndpoints(p0, p1);
                expect(tiQ.test(S, E).intersect).toBe(true);
                const r = fiQ.find(S, E);
                expect(r.intersect).toBe(true);
                // The whole segment is the intersection set.
                const cf = S.getCenteredForm();
                expectClose(r.parameter[0], -cf.extent, 1e-12, 1e-12);
                expectClose(r.parameter[1], cf.extent, 1e-12, 1e-12);
            });
    });

    it('pins the contained-segment regression', () => {
        // The unit sphere as an ellipsoid; the segment from (-0.2,0,0) to
        // (0.3,0,0) is strictly inside. Upstream's TI reported false while
        // its FI reported the whole segment.
        const E = Hyperellipsoid.fromCenterAxisExtent(Vector.zero(3),
            [Vector.fromArray([1, 0, 0]), Vector.fromArray([0, 1, 0]),
                Vector.fromArray([0, 0, 1])],
            Vector.fromArray([1, 1, 1]));
        const S = Segment.fromEndpoints(Vector.fromArray([-0.2, 0, 0]),
            Vector.fromArray([0.3, 0, 0]));
        expect(tiQ.test(S, E).intersect).toBe(true);
        const r = fiQ.find(S, E);
        expect(r.intersect).toBe(true);
        expect(r.numIntersections).toBe(2);
        expectVectorClose(r.point[0], S.p[0], 1e-12, 1e-12);
        expectVectorClose(r.point[1], S.p[1], 1e-12, 1e-12);
    });

    it('the FI parameters are in the centered-form domain', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            const r = fiQ.find(S, E);
            if (!r.intersect) {
                expect(r.numIntersections).toBe(0);
                return;
            }
            const cf = S.getCenteredForm();
            const scale = 1 + cf.extent;
            for (let i = 0; i < r.numIntersections; ++i) {
                // The documented convention: C + t*D with |t| <= extent.
                expect(Math.abs(r.parameter[i]))
                    .toBeLessThanOrEqual(cf.extent * (1 + 1e-12));
                expectVectorClose(r.point[i],
                    add(cf.center, mul(r.parameter[i], cf.direction)),
                    1e-12 * scale, 1e-12);
            }
            if (r.numIntersections === 2) {
                expect(r.parameter[0]).toBeLessThanOrEqual(r.parameter[1]);
            }
        });
    });

    it('each reported point is on the ellipsoid or is a clipped endpoint', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            const r = fiQ.find(S, E);
            if (!r.intersect) { return; }
            const cf = S.getCenteredForm();
            const scale = 1 + length(E.extent) + cf.extent;
            for (let i = 0; i < r.numIntersections; ++i) {
                const q = quadratic(E, r.point[i]);
                const clipped = Math.abs(Math.abs(r.parameter[i]) - cf.extent)
                    <= 0;
                if (!clipped) {
                    expectClose(q, 0, 1e-6 * scale, 1e-6);
                } else {
                    expect(q).toBeLessThanOrEqual(1e-9 * scale);
                }
            }
        });
    });

    it('is the line query clipped to [-extent, extent]', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            const cf = S.getCenteredForm();
            const line = defaultIntrSegment3Ellipsoid3FIResult();
            intrLine3Ellipsoid3FIDoQuery(cf.center, cf.direction, E, line);
            const r = fiQ.find(S, E);
            if (!line.intersect) {
                expect(r.intersect).toBe(false);
                return;
            }
            const t0 = line.parameter[0];
            const t1 = line.parameter[1];
            if (t1 < -cf.extent || t0 > cf.extent) {
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
                return;
            }
            expect(r.intersect).toBe(true);
            expect(r.parameter[0]).toBe(Math.max(t0, -cf.extent));
            expect(r.parameter[1]).toBe(Math.min(t1, cf.extent));
        });
    });

    it('reports a segment far from the ellipsoid as empty', () => {
        check(fc.tuple(arbEllipsoid, unitVector(3), positive(4, 1)),
            ([E, d, extra]) => {
                const maxExtent = Math.max(E.extent.values[0],
                    E.extent.values[1], E.extent.values[2]);
                const base = add(E.center, mul(maxExtent + extra + 1, d));
                const S = Segment.fromEndpoints(base, add(base, mul(3, d)));
                expect(tiQ.test(S, E).intersect).toBe(false);
                const r = fiQ.find(S, E);
                expect(r.intersect).toBe(false);
                expect(r.numIntersections).toBe(0);
            });
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid, rotationFrame(3),
            wellScaledVector(3)),
            ([S, E, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v), dot(Rot[2], v)]);
                const map = (v: Vector) => add(rot(v), t);
                const S2 = Segment.fromEndpoints(map(S.p[0]), map(S.p[1]));
                const E2 = Hyperellipsoid.fromCenterAxisExtent(map(E.center),
                    E.axis.map(rot), E.extent);
                const r0 = fiQ.find(S, E);
                const r1 = fiQ.find(S2, E2);
                if (r0.numIntersections === 1 || r1.numIntersections === 1) {
                    return;  // tangency or a clipped coincident endpoint
                }
                expect(r1.intersect).toBe(r0.intersect);
                expect(r1.numIntersections).toBe(r0.numIntersections);
            });
    });

    it('the exported DoQuery agrees with find on the parameters', () => {
        check(fc.tuple(arbSegment(3), arbEllipsoid), ([S, E]) => {
            const cf = S.getCenteredForm();
            const d = defaultIntrSegment3Ellipsoid3FIResult();
            intrSegment3Ellipsoid3FIDoQuery(cf.center, cf.direction, cf.extent,
                E, d);
            const f = fiQ.find(S, E);
            expect(d.intersect).toBe(f.intersect);
            expect(d.numIntersections).toBe(f.numIntersections);
            expect(d.parameter[0]).toBe(f.parameter[0]);
            expect(d.parameter[1]).toBe(f.parameter[1]);
        });
    });
});
