import { describe, it, expect } from 'vitest';
import { Ray } from '../src/Ray.js';
import { Segment } from '../src/Segment.js';
import { Vector, add, mul, normalize, sub, length } from '../src/Vector.js';
import {
    IntrRay2Segment2TI,
    IntrRay2Segment2FI
} from '../src/IntrRay2Segment2.js';

function vec(a: number[]): Vector {
    return Vector.fromArray(a);
}

function ray(p: number[], d: number[]): Ray {
    const dir = vec(d);
    normalize(dir);
    return Ray.fromOriginDirection(vec(p), dir);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(vec(p0), vec(p1));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('IntrRay2Segment2', () => {
    const ti = new IntrRay2Segment2TI();
    const fi = new IntrRay2Segment2FI();

    it('finds the single crossing of a transversal ray and segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, -1], [2, 1]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(2, 12);
        // The segment parameter is relative to the segment center.
        expect(result.segmentParameter[0]).toBeCloseTo(0, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[0].values[1]).toBeCloseTo(0, 12);
    });

    it('rejects a crossing beyond the segment endpoints', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, 2], [2, 3]);
        expect(ti.test(r, s).intersect).toBe(false);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('rejects a crossing behind the ray origin', () => {
        const r = ray([0, 0], [-1, 0]);
        const s = seg([2, -1], [2, 1]);
        expect(ti.test(r, s).intersect).toBe(false);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports no intersection for a parallel, non-collinear segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([1, 1], [5, 1]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports the overlap of a collinear segment', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([2, 0], [5, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(2);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.rayParameter[0]).toBeCloseTo(2, 12);
        expect(result.rayParameter[1]).toBeCloseTo(5, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(-1.5, 12);
        expect(result.segmentParameter[1]).toBeCloseTo(1.5, 12);
        expect(result.point[0].values[0]).toBeCloseTo(2, 12);
        expect(result.point[1].values[0]).toBeCloseTo(5, 12);
    });

    it('clips a collinear segment that straddles the ray origin', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-2, 0], [4, 0]);
        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(2);
        expect(result.rayParameter[0]).toBeCloseTo(0, 12);
        expect(result.rayParameter[1]).toBeCloseTo(4, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(-1, 12);
        expect(result.segmentParameter[1]).toBeCloseTo(3, 12);
    });

    it('reports no intersection for a collinear segment behind the ray', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-5, 0], [-2, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(false);
        expect(tiResult.numIntersections).toBe(0);
        expect(fi.find(r, s).intersect).toBe(false);
    });

    it('reports a single touching point for a collinear segment that ends at the ray origin', () => {
        const r = ray([0, 0], [1, 0]);
        const s = seg([-3, 0], [0, 0]);
        const tiResult = ti.test(r, s);
        expect(tiResult.intersect).toBe(true);
        expect(tiResult.numIntersections).toBe(1);

        const result = fi.find(r, s);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(0, 12);
        expect(result.segmentParameter[0]).toBeCloseTo(1.5, 12);
        expect(result.point[0].values).toEqual([0, 0]);
    });

    it('treats a degenerate zero-length segment as its center point', () => {
        // The centered form of a zero-length segment has a zero direction
        // vector, so the line-line query classifies the two lines as "the
        // same" and only the ray-parameter interval test remains. A
        // degenerate segment on the ray is reported at its parameter.
        const r = ray([0, 0], [1, 0]);
        const onRay = seg([3, 0], [3, 0]);
        const result = fi.find(r, onRay);
        expect(result.intersect).toBe(true);
        expect(result.numIntersections).toBe(1);
        expect(result.rayParameter[0]).toBeCloseTo(3, 12);
        expect(result.point[0].values[0]).toBeCloseTo(3, 12);

        // A degenerate segment behind the ray origin is rejected.
        const behind = seg([-3, 0], [-3, 0]);
        expect(fi.find(r, behind).intersect).toBe(false);
        expect(ti.test(r, behind).intersect).toBe(false);
    });

    it('agrees between the test and find queries on random configurations', () => {
        const rand = makeRandom(24680);
        let single = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const r = ray([6 * rand() - 3, 6 * rand() - 3],
                [2 * rand() - 1, 2 * rand() - 1]);
            const s = seg([6 * rand() - 3, 6 * rand() - 3],
                [6 * rand() - 3, 6 * rand() - 3]);
            const tiResult = ti.test(r, s);
            const fiResult = fi.find(r, s);
            expect(tiResult.intersect).toBe(fiResult.intersect);
            expect(tiResult.numIntersections).toBe(fiResult.numIntersections);

            if (fiResult.numIntersections === 1) {
                ++single;
                expect(fiResult.rayParameter[0]).toBeGreaterThanOrEqual(0);
                const centered = s.getCenteredForm();
                expect(Math.abs(fiResult.segmentParameter[0]))
                    .toBeLessThanOrEqual(centered.extent + 1e-12);
                const p0 = add(r.origin,
                    mul(fiResult.rayParameter[0], r.direction));
                const p1 = add(centered.center,
                    mul(fiResult.segmentParameter[0], centered.direction));
                expect(length(sub(p0, p1))).toBeLessThan(1e-9);
                expect(length(sub(p0, fiResult.point[0]))).toBeLessThan(1e-9);
            }
        }
        expect(single).toBeGreaterThan(20);
    });
});

// ---------------------------------------------------------------------------
// Verification (V31): property-based checks against the upstream header.
// ---------------------------------------------------------------------------
import {
    check, fc, latticeVector, unitVector, wellScaledVector, expectClose,
    expectVectorClose
} from './helpers/arbitraries.js';
import { IntrLine2Line2FI } from '../src/IntrLine2Line2.js';
import { Line } from '../src/Line.js';
import { dotPerp } from '../src/Vector2.js';

// A ray and a segment whose directions are transverse with a comfortable
// margin, so the line-line parameters are well conditioned.
const transverseRaySegment = fc.tuple(wellScaledVector(2), unitVector(2),
    wellScaledVector(2), unitVector(2),
    fc.double({ min: 0.5, max: 5, noNaN: true }))
    .filter(([, d0, , d1]) => Math.abs(dotPerp(d0, d1)) > 0.1)
    .map(([o0, d0, o1, d1, len]) => ({
        ray: Ray.fromOriginDirection(o0, d0),
        segment: Segment.fromEndpoints(o1, add(o1, mul(len, d1)))
    }));

// An exactly representable unit direction, so collinear configurations built
// from integer offsets stay exactly collinear.
const axisDirection2 = fc.integer({ min: 0, max: 3 }).map(i =>
    Vector.fromArray([[1, 0], [0, 1], [-1, 0], [0, -1]][i]));

describe('IntrRay2Segment2 verification', () => {
    const tiq = new IntrRay2Segment2TI();
    const fiq = new IntrRay2Segment2FI();
    const llq = new IntrLine2Line2FI();

    it('TI and FI agree on intersect and numIntersections', () => {
        check(transverseRaySegment, ({ ray: r, segment: s }) => {
            const t = tiq.test(r, s);
            const f = fiq.find(r, s);
            expect(t.intersect).toBe(f.intersect);
            expect(t.numIntersections).toBe(f.numIntersections);
        });
    });

    it('a transverse hit is a line-line hit clipped to ray and segment', () => {
        check(transverseRaySegment, ({ ray: r, segment: s }) => {
            const cf = s.getCenteredForm();
            const ll = llq.find(
                Line.fromOriginDirection(r.origin, r.direction),
                Line.fromOriginDirection(cf.center, cf.direction));
            expect(ll.numIntersections).toBe(1);
            const inside = ll.line0Parameter[0] >= 0
                && Math.abs(ll.line1Parameter[0]) <= cf.extent;
            const f = fiq.find(r, s);
            expect(f.intersect).toBe(inside);
            if (!inside) {
                expect(f.numIntersections).toBe(0);
                return;
            }
            expect(f.numIntersections).toBe(1);
            expect(f.rayParameter[0]).toBeGreaterThanOrEqual(0);
            expect(Math.abs(f.segmentParameter[0]))
                .toBeLessThanOrEqual(cf.extent);
            expectVectorClose(f.point[0],
                add(r.origin, mul(f.rayParameter[0], r.direction)), 0, 0);
            // The segment parameter is the centered-form parameter.
            expectVectorClose(f.point[0],
                add(cf.center, mul(f.segmentParameter[0], cf.direction)),
                1e-9, 1e-12);
        });
    });

    it('a ray fired at a sampled segment point always hits', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 0.05, max: 0.95, noNaN: true }), unitVector(2),
            fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([p0, d, len, u, e, back]) => {
                if (Math.abs(dotPerp(d, e)) <= 0.1) {
                    return;
                }
                const s = Segment.fromEndpoints(p0, add(p0, mul(len, d)));
                const target = add(s.p[0], mul(u, sub(s.p[1], s.p[0])));
                const r = Ray.fromOriginDirection(sub(target, mul(back, e)), e);
                const f = fiq.find(r, s);
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(1);
                expectVectorClose(f.point[0], target, 1e-9, 1e-12);
                // Reversing the ray direction misses.
                const rev = Ray.fromOriginDirection(r.origin, mul(-1, e));
                expect(fiq.find(rev, s).intersect).toBe(false);
            });
    });

    it('a collinear overlap reconstructs its endpoints on the ray', () => {
        check(fc.tuple(latticeVector(2), axisDirection2,
            fc.integer({ min: -6, max: 6 }), fc.integer({ min: 1, max: 5 })),
            ([o, d, k, len]) => {
                const r = Ray.fromOriginDirection(o, d);
                const q0 = add(o, mul(k, d));
                const s = Segment.fromEndpoints(q0, add(q0, mul(len, d)));
                const f = fiq.find(r, s);
                const t = tiq.test(r, s);
                expect(f.intersect).toBe(t.intersect);
                expect(f.numIntersections).toBe(t.numIntersections);
                // The segment spans ray parameters [k, k+len]; the ray covers
                // [0, +infinity), so the overlap is [max(k,0), k+len].
                const lo = Math.max(k, 0), hi = k + len;
                if (hi < 0) {
                    expect(f.intersect).toBe(false);
                    return;
                }
                expect(f.intersect).toBe(true);
                expect(f.numIntersections).toBe(lo < hi ? 2 : 1);
                for (let i = 0; i < f.numIntersections; ++i) {
                    expect(f.rayParameter[i]).toBeGreaterThanOrEqual(0);
                    expectVectorClose(f.point[i],
                        add(r.origin, mul(f.rayParameter[i], r.direction)),
                        0, 0);
                }
                expectClose(f.rayParameter[0], lo, 0, 0);
                if (f.numIntersections === 2) {
                    expectClose(f.rayParameter[1], hi, 0, 0);
                    expect(f.rayParameter[0])
                        .toBeLessThanOrEqual(f.rayParameter[1]);
                    expect(f.segmentParameter[0])
                        .toBeLessThanOrEqual(f.segmentParameter[1]);
                }
            });
    });

    it('a segment parallel to but off the ray never intersects', () => {
        check(fc.tuple(wellScaledVector(2), unitVector(2),
            fc.double({ min: 0.5, max: 5, noNaN: true }),
            fc.double({ min: 0.5, max: 5, noNaN: true })),
            ([o, d, offset, len]) => {
                const n = Vector.fromArray([-d.values[1], d.values[0]]);
                const r = Ray.fromOriginDirection(o, d);
                const q0 = add(o, mul(offset, n));
                const s = Segment.fromEndpoints(q0, add(q0, mul(len, d)));
                expect(fiq.find(r, s).intersect).toBe(false);
                expect(tiq.test(r, s).intersect).toBe(false);
            });
    });

    it('no result field is NaN when the query reports an intersection', () => {
        check(transverseRaySegment, ({ ray: r, segment: s }) => {
            const f = fiq.find(r, s);
            if (!f.intersect) {
                return;
            }
            for (let i = 0; i < f.numIntersections; ++i) {
                expect(Number.isNaN(f.rayParameter[i])).toBe(false);
                expect(Number.isNaN(f.segmentParameter[i])).toBe(false);
                for (const x of f.point[i].values) {
                    expect(Number.isNaN(x)).toBe(false);
                }
            }
        });
    });
});
