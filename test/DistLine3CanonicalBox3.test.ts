import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistLine3CanonicalBox3 } from '../src/DistLine3CanonicalBox3.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { check, expectClose, expectVectorClose, fc, positive, rotationFrame, seededRandom, wellScaled, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
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

// A ternary search for the minimum of the convex function
// t -> pointBoxSqrDistance(origin + t*direction, box).
function sampledMinimum(l: Line, b: CanonicalBox): number {
    let lo = -1e6;
    let hi = 1e6;
    const f = (t: number) =>
        pointBoxSqrDistance(add(l.origin, mul(t, l.direction)), b);
    for (let i = 0; i < 250; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }
    return f(0.5 * (lo + hi));
}

describe('DistLine3CanonicalBox3', () => {
    it('returns 0, not NaN, for a line through a zero-extent box', () => {
        // The upstream incremental sqrDistance form goes slightly negative
        // here and sqrt() gives NaN (found by V20's equivariance property on
        // a flat OrientedBox). The port clamps negative round-off to 0.
        // Pre-fix: sqrDistance = -1.7763568394002505e-15, distance = NaN.
        const box = new CanonicalBox(3);
        box.extent = Vector.fromArray([0, 3.999999999999937, 3.999999999970607]);
        const line = Line.fromOriginDirection(
            Vector.fromArray([0, -7.499019495845319, 0]),
            Vector.fromArray([-7.55637869396348e-9, -0.9763428607601213,
                0.21622816246442592]));
        const result = new DistLine3CanonicalBox3().compute(line, box);
        expect(Number.isNaN(result.distance)).toBe(false);
        expect(result.sqrDistance).toBeGreaterThanOrEqual(0);
        expect(result.distance).toBeLessThan(1e-6);
    });

    const query = new DistLine3CanonicalBox3();
    const unitBox = box(1, 1, 1);

    it('measures an axis-parallel line above the box', () => {
        const result = query.compute(line([0, 0, 4], [1, 0, 0]), unitBox);
        expect(result.distance).toBeCloseTo(3, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
    });

    it('reports zero distance for a line through the box', () => {
        const result = query.compute(line([0, 0, 0], [1, 2, 3]), unitBox);
        expect(result.distance).toBeCloseTo(0, 10);
        expect(result.sqrDistance).toBeCloseTo(0, 10);
    });

    it('handles a degenerate zero direction as a point query', () => {
        const result = query.compute(line([5, -9, 0.5], [0, 0, 0]), unitBox);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([5, -9, 0.5]);
        expect(result.closest[1].values).toEqual([1, -1, 0.5]);
        expect(result.sqrDistance).toBeCloseTo(16 + 64, 10);
    });

    it('handles a line parallel to an axis but off the box in two axes',
        () => {
            const result = query.compute(line([4, 5, 0], [0, 0, 1]), unitBox);
            expect(result.sqrDistance).toBeCloseTo(9 + 16, 10);
            expect(result.distance).toBeCloseTo(5, 10);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
            expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
        });

    it('measures a diagonal line missing a corner', () => {
        // The line (4,0,0)+t*(-1,1,0) lies in z = 0 and its distance to the
        // unit box is the same as in the 2D case: sqrt(2).
        const result = query.compute(line([4, 0, 0], [-1, 1, 0]), unitBox);
        expect(result.distance).toBeCloseTo(Math.SQRT2, 8);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 8);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 8);
    });

    it('gives the same distance for a line and its reversed direction', () => {
        const b = box(1, 2, 0.5);
        const a = query.compute(line([-3, 4, 5], [1, -2, 3]), b);
        const c = query.compute(line([-3, 4, 5], [-1, 2, -3]), b);
        expect(c.distance).toBeCloseTo(a.distance, 9);
        for (let i = 0; i < 3; ++i) {
            expect(c.closest[1].values[i]).toBeCloseTo(a.closest[1].values[i],
                7);
        }
    });

    it('matches a numeric minimization over the line', () => {
        let seed = 13579;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        for (let trial = 0; trial < 200; ++trial) {
            const b = box(Math.abs(rand()) + 0.1, Math.abs(rand()) + 0.1,
                Math.abs(rand()) + 0.1);
            const dir = [rand(), rand(), rand()];
            if (Math.hypot(dir[0], dir[1], dir[2]) < 0.5) {
                continue;
            }
            const l = line([rand(), rand(), rand()], dir);
            const result = query.compute(l, b);

            expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b), 6);
            expect(result.distance).toBeCloseTo(
                Math.sqrt(result.sqrDistance), 10);

            // The line point matches the reported parameter.
            const onLine = add(l.origin, mul(result.parameter, l.direction));
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[0].values[i]).toBeCloseTo(
                    onLine.values[i], 7);
            }

            // The box point lies in the box.
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(result.closest[1].values[i]))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-7);
            }

            // The reported pair realizes the reported squared distance.
            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 6);
        }
    });

    it('exercises the axis-aligned direction cases', () => {
        let seed = 24601;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const dirs = [
            [1, 0, 0], [0, 1, 0], [0, 0, 1],
            [1, 1, 0], [1, 0, 1], [0, 1, 1],
            [-1, 0, 0], [0, -1, 0], [0, 0, -1],
            [-1, 1, 0], [1, 0, -1], [0, -1, -1]
        ];
        const b = box(1.25, 0.5, 2);
        for (const dir of dirs) {
            for (let trial = 0; trial < 10; ++trial) {
                const l = line([rand(), rand(), rand()], dir);
                const result = query.compute(l, b);
                expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b),
                    6);
                const diff = sub(result.closest[0], result.closest[1]);
                expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 6);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistLine3CanonicalBox3.ts against the upstream header
// DistLine3CanonicalBox3.h.
// ---------------------------------------------------------------------------

function rot3(R: readonly Vector[], p: Vector): Vector {
    let q = mul(p.values[0], R[0]);
    for (let i = 1; i < 3; ++i) {
        q = add(q, mul(p.values[i], R[i]));
    }
    return q;
}

/** Squared distance from a point to the solid canonical box, by clamping. */
function pointCanonicalSqrDistance(p: Vector, extent: Vector): number {
    let sum = 0;
    for (let i = 0; i < p.size; ++i) {
        const over = Math.abs(p.values[i]) - extent.values[i];
        if (over > 0) { sum += over * over; }
    }
    return sum;
}

/**
 * The exact line-box distance, computed independently of the query under
 * test. The function t -> distance(P + t*D, box) is convex (a norm distance
 * to a convex set composed with an affine map), so a ternary search on a
 * bracket containing the minimizer converges to the global minimum. The
 * bracket is |t - t_c| <= (h + R)/|D| around the foot t_c of the
 * perpendicular from the box center, with h the perpendicular distance and R
 * the radius of a ball about the center containing the box; outside it the
 * distance already exceeds the value at t_c.
 */
function referenceLineCanonicalDistance(line: Line, extent: Vector): number {
    const radius = length(extent);
    const dd = dot(line.direction, line.direction);
    const tc = -dot(line.direction, line.origin) / dd;
    const h = length(add(line.origin, mul(tc, line.direction)));
    const half = (h + radius) / Math.sqrt(dd);
    let lo = tc - half - 1;
    let hi = tc + half + 1;
    const f = (t: number): number => pointCanonicalSqrDistance(
        add(line.origin, mul(t, line.direction)), extent);
    for (let k = 0; k < 200; ++k) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) <= f(m1)) { hi = m1; }
        else { lo = m0; }
    }
    return Math.sqrt(Math.min(f(lo), f(hi), f(0.5 * (lo + hi))));
}

const nonUnitLine3 = fc.tuple(wellScaledVector(3, -8, 8),
    wellScaledVector(3, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

const canonicalBox3 = fc.array(positive(5), { minLength: 3, maxLength: 3 })
    .map(e => CanonicalBox.fromExtent(Vector.fromArray(e)));

describe('DistLine3CanonicalBox3 verification', () => {
    const query = new DistLine3CanonicalBox3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(nonUnitLine3, canonicalBox3), ([line, cbox]) => {
                const r = query.compute(line, cbox);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                const diff = sub(r.closest[0], r.closest[1]);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-7, 1e-7);
                expectVectorClose(r.closest[0],
                    add(line.origin, mul(r.parameter, line.direction)), 1e-9,
                    1e-9);
                for (let i = 0; i < 3; ++i) {
                    expect(Math.abs(r.closest[1].values[i]))
                        .toBeLessThanOrEqual(cbox.extent.values[i] + 1e-9);
                }
            });
        });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(nonUnitLine3, canonicalBox3), ([line, cbox]) => {
            expectClose(query.compute(line, cbox).distance,
                referenceLineCanonicalDistance(line, cbox.extent), 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(nonUnitLine3, canonicalBox3), ([line, cbox]) => {
            const o = line.origin.clone();
            const d = line.direction.clone();
            const e = cbox.extent.clone();
            query.compute(line, cbox);
            expectVectorClose(line.origin, o, 0, 0);
            expectVectorClose(line.direction, d, 0, 0);
            expectVectorClose(cbox.extent, e, 0, 0);
        });
    });

    it('is equivariant under the reflections the query applies internally',
        () => {
            check(fc.tuple(nonUnitLine3, canonicalBox3,
                fc.array(fc.boolean(), { minLength: 3, maxLength: 3 })),
            ([line, cbox, flip]) => {
                const refl = (p: Vector): Vector => Vector.fromArray(
                    [0, 1, 2].map(i =>
                        flip[i] ? -p.values[i] : p.values[i]));
                const movedLine = Line.fromOriginDirection(refl(line.origin),
                    refl(line.direction));
                expectClose(query.compute(line, cbox).distance,
                    query.compute(movedLine, cbox).distance, 0, 0);
            });
        });

    it('agrees with the 2D line-box query on a flat embedding', () => {
        // A line in the z = 0 plane sees the canonical box exactly as the 2D
        // aligned box [-e0,e0] x [-e1,e1] (z = 0 is interior to [-e2,e2]).
        check(fc.tuple(wellScaledVector(2, -8, 8), wellScaledVector(2, -3, 3)
            .filter(d => length(d) > 0.25), canonicalBox3),
        ([o2, d2, cbox]) => {
            const line3 = Line.fromOriginDirection(
                Vector.fromArray([o2.values[0], o2.values[1], 0]),
                Vector.fromArray([d2.values[0], d2.values[1], 0]));
            const box2 = AlignedBox.fromMinMax(
                Vector.fromArray([-cbox.extent.values[0],
                    -cbox.extent.values[1]]),
                Vector.fromArray([cbox.extent.values[0],
                    cbox.extent.values[1]]));
            const r3 = query.compute(line3, cbox);
            const r2 = new DistLine2AlignedBox2().compute(
                Line.fromOriginDirection(o2, d2), box2);
            expectClose(r3.distance, r2.distance, 1e-8, 1e-8);
        });
    });

    it('reduces to a point-box distance for a zero direction', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), canonicalBox3),
            ([p, cbox]) => {
                const r = query.compute(
                    Line.fromOriginDirection(p, new Vector(3)), cbox);
                expect(r.parameter).toBe(0);
                expectClose(r.sqrDistance,
                    pointCanonicalSqrDistance(p, cbox.extent), 1e-12, 1e-12);
                const pcb = new DistPointCanonicalBox().compute(p, cbox);
                expectClose(r.distance, pcb.distance, 1e-12, 1e-12);
                expectVectorClose(r.closest[1], pcb.closest[1], 1e-12, 1e-12);
            });
    });

    it('reports zero distance when the line meets the box', () => {
        check(fc.tuple(canonicalBox3, wellScaledVector(3, -1, 1),
            rotationFrame(3), wellScaled(-5, 5)), ([cbox, frac, R, t]) => {
            const q = Vector.fromArray([0, 1, 2].map(i =>
                frac.values[i] * cbox.extent.values[i]));
            const line = Line.fromOriginDirection(add(q, mul(t, R[0])), R[0]);
            expect(query.compute(line, cbox).distance)
                .toBeLessThanOrEqual(1e-9);
        });
    });

    it('is minimal over sampled line/box point pairs', () => {
        const rand = seededRandom(0x51e2);
        check(fc.tuple(nonUnitLine3, canonicalBox3), ([line, cbox]) => {
            const r = query.compute(line, cbox);
            const q = new Vector(3);
            for (let k = 0; k < 24; ++k) {
                const t = 40 * (rand() - 0.5);
                for (let i = 0; i < 3; ++i) {
                    q.values[i] = cbox.extent.values[i] * (2 * rand() - 1);
                }
                const gap = length(
                    sub(add(line.origin, mul(t, line.direction)), q));
                expect(r.distance).toBeLessThanOrEqual(gap + 1e-9 * (1 + gap));
            }
        }, 60);
    });
});
