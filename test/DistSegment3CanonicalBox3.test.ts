import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { DistSegment3CanonicalBox3 } from '../src/DistSegment3CanonicalBox3.js';
import { Segment } from '../src/Segment.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

// The exact squared distance from a point to the solid canonical box.
function pointBoxSqrDistance(p: Vector, b: CanonicalBox): number {
    let sum = 0;
    for (let i = 0; i < 3; ++i) {
        const e = b.extent.values[i];
        const c = Math.min(Math.max(p.values[i], -e), e);
        sum += (p.values[i] - c) * (p.values[i] - c);
    }
    return sum;
}

// The squared distance from segment(t) to the solid box is a convex function
// of t on [0,1], so a ternary search finds its minimum. A dense sampling is
// used as a second opinion.
function bruteForceSqrDistance(s: Segment, b: CanonicalBox): number {
    const f = (u: number) => pointBoxSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), b);

    let lo = 0;
    let hi = 1;
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }

    let best = f(0.5 * (lo + hi));
    for (let i = 0; i <= 2000; ++i) {
        const value = f(i / 2000);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

// Verify the internal consistency of a result: the closest points lie on
// their primitives and realize the reported distance.
function expectConsistent(
    result: ReturnType<DistSegment3CanonicalBox3['compute']>,
    s: Segment, b: CanonicalBox): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 3; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
        expect(Math.abs(result.closest[1].values[i]))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }

    const diff = sub(result.closest[0], result.closest[1]);
    let sqrLen = 0;
    for (let i = 0; i < 3; ++i) {
        sqrLen += diff.values[i] * diff.values[i];
    }
    expect(Math.sqrt(sqrLen)).toBeCloseTo(result.distance, 6);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);
}

// A small deterministic linear congruential generator.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistSegment3CanonicalBox3', () => {
    const query = new DistSegment3CanonicalBox3();
    const unitBox = box(1, 1, 1);

    it('measures a segment parallel to a face', () => {
        const s = segment([-2, 0, 4], [2, 0, 4]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment passing through the box', () => {
        const s = segment([-3, 0, 0], [3, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance for a segment strictly inside the box', () => {
        const s = segment([-0.5, 0.25, 0.1], [0.5, -0.25, -0.1]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('reports zero distance when an endpoint touches a face', () => {
        const s = segment([1, 0, 0], [4, 0, 0]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            // The line x = t, y = 0, z = 0 has its minimum at t <= 0 relative
            // to p[0] = (3,0,0), so the closest segment point is p[0].
            const s = segment([3, 0, 0], [6, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(2, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
            expectConsistent(result, s, unitBox);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([6, 0, 0], [3, 0, 0]);
            const result = query.compute(s, unitBox);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(2, 10);
            expectConsistent(result, s, unitBox);
        });

    it('measures the distance to a box corner', () => {
        // The segment lies far along the (1,1,1) diagonal direction; the
        // closest box point is the corner (1,1,1).
        const s = segment([3, 3, 3], [4, 5, 6]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(Math.sqrt(12), 10);
        expect(result.parameter).toBe(0);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 10);
        }
        expectConsistent(result, s, unitBox);
    });

    it('gives the same distance for both segment orientations', () => {
        const cases: Array<[number[], number[]]> = [
            [[3, 0, 0], [6, 0, 0]],
            [[-2, 0, 4], [2, 0, 4]],
            [[-3, -3, -3], [3, 3, 3]],
            [[2, 2, -5], [2, 2, 5]],
            [[0.5, 0.5, 0.5], [7, 1, -3]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), unitBox);
            const backward = query.compute(segment(p1, p0), unitBox);
            expect(forward.distance).toBeCloseTo(backward.distance, 10);
        }
    });

    it('handles a degenerate zero-length segment outside the box', () => {
        const s = segment([3, 4, 5], [3, 4, 5]);
        const result = query.compute(s, unitBox);
        const expected = new DistPointCanonicalBox().compute(v(3, 4, 5),
            unitBox);
        expect(result.distance).toBeCloseTo(expected.distance, 12);
        expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
        expectConsistent(result, s, unitBox);
    });

    it('handles a degenerate zero-length segment inside the box', () => {
        const s = segment([0.25, -0.5, 0.75], [0.25, -0.5, 0.75]);
        const result = query.compute(s, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitBox);
    });

    it('handles a degenerate (flat) box', () => {
        const flat = box(2, 3, 0);
        const s = segment([0, 0, 1], [0, 0, 5]);
        const result = query.compute(s, flat);
        expect(result.distance).toBeCloseTo(1, 10);
        expect(result.parameter).toBe(0);
        expectConsistent(result, s, flat);
    });

    it('matches a brute-force minimization on analytic configurations', () => {
        const boxes = [box(1, 1, 1), box(2, 1, 0.5), box(0.5, 3, 2)];
        const cases: Array<[number[], number[]]> = [
            [[-4, 2, 2], [4, 2, 2]],
            [[2, 2, 2], [3, 5, 7]],
            [[-1, -1, -1], [-1, -1, 4]],
            [[0.1, 0.1, 0.1], [0.2, 0.6, 0.9]],
            [[-5, 0.5, 0.5], [-1.5, 0.5, 0.5]],
            [[3, 0, 0], [0, 3, 0]],
            [[0, -6, 0], [0, -2, 0]]
        ];
        for (const b of boxes) {
            for (const [p0, p1] of cases) {
                const s = segment(p0, p1);
                const result = query.compute(s, b);
                const expected = Math.sqrt(bruteForceSqrDistance(s, b));
                expect(result.distance).toBeCloseTo(expected, 6);
                expectConsistent(result, s, b);
            }
        }
    });

    it('matches a brute-force minimization on random configurations', () => {
        const random = makeRandom(987654321);
        const coord = () => 6 * random() - 3;
        for (let trial = 0; trial < 400; ++trial) {
            const b = box(0.25 + 2 * random(), 0.25 + 2 * random(),
                0.25 + 2 * random());
            const s = segment([coord(), coord(), coord()],
                [coord(), coord(), coord()]);
            const result = query.compute(s, b);
            const expected = Math.sqrt(bruteForceSqrDistance(s, b));
            expect(result.distance).toBeCloseTo(expected, 6);
            expectConsistent(result, s, b);

            // Reversing the segment must not change the distance.
            const reversed = query.compute(
                Segment.fromEndpoints(s.p[1], s.p[0]), b);
            expect(reversed.distance).toBeCloseTo(result.distance, 6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistSegment3CanonicalBox3.h.
// ---------------------------------------------------------------------------

describe('DistSegment3CanonicalBox3 verification', () => {
    const query = new DistSegment3CanonicalBox3();
    const lineQuery = new DistLine3CanonicalBox3();
    const pointQuery = new DistPointCanonicalBox();

    const boxArb = fc.array(finite(0, 4), { minLength: 3, maxLength: 3 })
        .map(e => CanonicalBox.fromExtent(Vector.fromArray(e)));

    const segmentArb = fc.tuple(wellScaledVector(3, -8, 8),
        wellScaledVector(3, -8, 8))
        .filter(([p0, p1]) => length(sub(p1, p0)) > 0.5)
        .map(([p0, p1]) => Segment.fromEndpoints(p0, p1));

    function pointBoxDistance(p: Vector, b: CanonicalBox): number {
        let sum = 0;
        for (let i = 0; i < 3; ++i) {
            const over = Math.abs(p.values[i]) - b.extent.values[i];
            if (over > 0) { sum += over * over; }
        }
        return Math.sqrt(sum);
    }

    function ternaryMin(f: (t: number) => number, lo: number,
        hi: number): number {
        let a = lo, b = hi;
        for (let i = 0; i < 200; ++i) {
            const m0 = a + (b - a) / 3;
            const m1 = b - (b - a) / 3;
            if (f(m0) <= f(m1)) { b = m1; } else { a = m0; }
        }
        return f(0.5 * (a + b));
    }

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(segmentArb, boxArb), ([seg, b]) => {
            const r = query.compute(seg, b);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expect(r.parameter).toBeGreaterThanOrEqual(0);
            expect(r.parameter).toBeLessThanOrEqual(1);
            const d = sub(seg.p[1], seg.p[0]);
            expectVectorClose(r.closest[0], add(seg.p[0], mul(r.parameter, d)),
                1e-9, 1e-9);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-8, 1e-8);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(r.closest[1].values[i]))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
            }
        });
    });

    it('matches an independent convex minimization along the segment', () => {
        check(fc.tuple(segmentArb, boxArb), ([seg, b]) => {
            const r = query.compute(seg, b);
            const d = sub(seg.p[1], seg.p[0]);
            const best = ternaryMin(
                t => pointBoxDistance(add(seg.p[0], mul(t, d)), b), 0, 1);
            // See the DistLine3AlignedBox3 note on the incremental
            // canonical-box accumulation.
            expectClose(r.distance, best, 2e-6, 1e-9);
        }, 100);
    });

    it('clamps the line parameter to the segment', () => {
        check(fc.tuple(segmentArb, boxArb), ([seg, b]) => {
            const d = sub(seg.p[1], seg.p[0]);
            const line = Line.fromOriginDirection(seg.p[0], d);
            const rl = lineQuery.compute(line, b);
            const rs = query.compute(seg, b);
            if (rl.parameter >= 0 && rl.parameter <= 1) {
                expectClose(rs.distance, rl.distance, 1e-12, 1e-12);
                expectClose(rs.parameter, rl.parameter, 1e-12, 1e-12);
                expectVectorClose(rs.closest[1], rl.closest[1], 1e-12, 1e-12);
            } else {
                const endpoint = rl.parameter < 0 ? seg.p[0] : seg.p[1];
                expect(rs.parameter).toBe(rl.parameter < 0 ? 0 : 1);
                expectVectorClose(rs.closest[0], endpoint, 1e-12, 1e-12);
                const rp = pointQuery.compute(endpoint, b);
                expectClose(rs.distance, rp.distance, 1e-12, 1e-12);
                expectVectorClose(rs.closest[1], rp.closest[1], 1e-12, 1e-12);
            }
            expect(rs.distance).toBeGreaterThanOrEqual(rl.distance - 1e-12);
        });
    });

    it('is symmetric under swapping the segment endpoints', () => {
        check(fc.tuple(segmentArb, boxArb), ([seg, b]) => {
            const reversed = Segment.fromEndpoints(seg.p[1], seg.p[0]);
            const r0 = query.compute(seg, b);
            const r1 = query.compute(reversed, b);
            expectClose(r0.distance, r1.distance, 2e-6, 1e-9);
        });
    });

    it('reports zero distance for a segment inside the box', () => {
        check(fc.tuple(boxArb, fc.array(finite(-1, 1),
            { minLength: 3, maxLength: 3 }), fc.array(finite(-1, 1),
                { minLength: 3, maxLength: 3 })), ([b, u0, u1]) => {
            const p0 = Vector.fromArray([u0[0] * b.extent.values[0],
                u0[1] * b.extent.values[1], u0[2] * b.extent.values[2]]);
            const p1 = Vector.fromArray([u1[0] * b.extent.values[0],
                u1[1] * b.extent.values[1], u1[2] * b.extent.values[2]]);
            if (length(sub(p1, p0)) < 1e-6) { return; }
            const r = query.compute(Segment.fromEndpoints(p0, p1), b);
            expectClose(r.distance, 0, 1e-9, 1e-9);
        });
    });

    it('is invariant under reflection of any coordinate axis', () => {
        check(fc.tuple(segmentArb, boxArb, fc.nat(2)), ([seg, b, k]) => {
            const flip = (p: Vector): Vector => {
                const q = p.clone();
                q.values[k] = -q.values[k];
                return q;
            };
            const r0 = query.compute(seg, b);
            const r1 = query.compute(
                Segment.fromEndpoints(flip(seg.p[0]), flip(seg.p[1])), b);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
        });
    });
});
