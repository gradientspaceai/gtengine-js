import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine2AlignedBox2 } from '../src/DistLine2AlignedBox2.js';
import { Line } from '../src/Line.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, rotationFrame, wellScaled, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function box(min: number[], max: number[]): AlignedBox {
    return AlignedBox.fromMinMax(v(...min), v(...max));
}

// The exact squared distance from a point to a solid aligned box.
function pointBoxSqrDistance(p: Vector, b: AlignedBox): number {
    let sum = 0;
    for (let i = 0; i < 2; ++i) {
        const c = Math.min(Math.max(p.values[i], b.min.values[i]),
            b.max.values[i]);
        sum += (p.values[i] - c) * (p.values[i] - c);
    }
    return sum;
}

// A ternary search for the minimum of the convex function
// t -> pointBoxSqrDistance(origin + t*direction, box).
function sampledMinimum(l: Line, b: AlignedBox): number {
    let lo = -1e6;
    let hi = 1e6;
    const f = (t: number) =>
        pointBoxSqrDistance(add(l.origin, mul(t, l.direction)), b);
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
    return f(0.5 * (lo + hi));
}

describe('DistLine2AlignedBox2', () => {
    const query = new DistLine2AlignedBox2();

    it('measures a horizontal line above the box', () => {
        const result = query.compute(line([0, 5], [1, 0]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
    });

    it('measures a vertical line to the right of the box', () => {
        const result = query.compute(line([7, 0], [0, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(6, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[0].values[0]).toBeCloseTo(7, 10);
    });

    it('reports zero distance for a line through the box', () => {
        const result = query.compute(line([0, 0], [1, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(0, 10);
    });

    it('measures a diagonal line missing a corner', () => {
        // The line x + y = 4, i.e. (4,0) + t*(-1,1), and the unit box. The
        // closest box point is the corner (1,1), at distance
        // (4-2)/sqrt(2) = sqrt(2).
        const result = query.compute(line([4, 0], [-1, 1]),
            box([-1, -1], [1, 1]));
        expect(result.distance).toBeCloseTo(Math.SQRT2, 10);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 10);
    });

    it('handles a degenerate zero direction as a point query', () => {
        const result = query.compute(line([5, 9], [0, 0]),
            box([-1, -1], [1, 1]));
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values).toEqual([5, 9]);
        expect(result.closest[1].values).toEqual([1, 1]);
        expect(result.sqrDistance).toBeCloseTo(16 + 64, 10);
    });

    it('handles an off-center box', () => {
        const b = box([10, 20], [12, 26]);
        // The infinite line y = 23 crosses the box, so the distance is 0.
        expect(query.compute(line([0, 23], [1, 0]), b).distance)
            .toBeCloseTo(0, 10);
        // The line y = 30 is above the box.
        const result = query.compute(line([0, 30], [1, 0]), b);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[1].values[1]).toBeCloseTo(26, 10);
    });

    it('gives the same answer for a line and its reversed direction', () => {
        // The line x + y = 12 misses the box, whose far corner is (3,4). The
        // distance is (12-7)/sqrt(2).
        const b = box([-1, -2], [3, 4]);
        const base = line([12, 0], [-1, 1]);
        const flipped = line([12, 0], [1, -1]);
        const a = query.compute(base, b);
        const c = query.compute(flipped, b);
        expect(a.distance).toBeCloseTo(5 / Math.SQRT2, 10);
        expect(c.distance).toBeCloseTo(a.distance, 10);
        expect(a.closest[1].values).toEqual([3, 4]);
        expect(c.closest[1].values[0]).toBeCloseTo(3, 8);
        expect(c.closest[1].values[1]).toBeCloseTo(4, 8);
    });

    it('matches a numeric minimization and reports consistent points', () => {
        let seed = 90210;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 10 - 5;
        };
        for (let trial = 0; trial < 150; ++trial) {
            const lo = [rand(), rand()];
            const b = box(lo, [lo[0] + Math.abs(rand()) + 0.1,
                lo[1] + Math.abs(rand()) + 0.1]);
            const dir = [rand(), rand()];
            if (Math.hypot(dir[0], dir[1]) < 0.25) {
                continue;
            }
            const l = line([rand(), rand()], dir);
            const result = query.compute(l, b);

            expect(result.sqrDistance).toBeCloseTo(sampledMinimum(l, b), 7);

            // The line point matches the reported parameter.
            const onLine = add(l.origin, mul(result.parameter, l.direction));
            expect(result.closest[0].values[0]).toBeCloseTo(onLine.values[0],
                7);
            expect(result.closest[0].values[1]).toBeCloseTo(onLine.values[1],
                7);

            // The box point lies in the box.
            for (let i = 0; i < 2; ++i) {
                expect(result.closest[1].values[i]).toBeGreaterThanOrEqual(
                    b.min.values[i] - 1e-8);
                expect(result.closest[1].values[i]).toBeLessThanOrEqual(
                    b.max.values[i] + 1e-8);
            }

            const diff = sub(result.closest[0], result.closest[1]);
            expect(dot(diff, diff)).toBeCloseTo(result.sqrDistance, 7);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistLine2AlignedBox2.ts against the upstream header DistLine2AlignedBox2.h.
// ---------------------------------------------------------------------------

function rot2(R: readonly Vector[], p: Vector): Vector {
    return add(mul(p.values[0], R[0]), mul(p.values[1], R[1]));
}

/**
 * The exact line-box distance, computed independently of the query under
 * test. The function t -> distance(P + t*D, box) is convex (a norm distance
 * to a convex set composed with an affine map), so a ternary search on a
 * bracket containing the minimizer converges to the global minimum.
 *
 * Bracket: let t_c be the parameter of the foot of the perpendicular from the
 * box center C, h the perpendicular distance and R the radius of a ball about
 * C containing the box. For |t - t_c|*|D| > h + R the distance is at least
 * |t - t_c|*|D| - R > h >= f(t_c), so the minimizer lies within
 * |t - t_c| <= (h + R)/|D|.
 */
function referenceLineBoxDistance(line: Line, box: AlignedBox): number {
    const center = mul(0.5, add(box.min, box.max));
    const radius = 0.5 * length(sub(box.max, box.min));
    const dd = dot(line.direction, line.direction);
    const tc = dot(line.direction, sub(center, line.origin)) / dd;
    const foot = add(line.origin, mul(tc, line.direction));
    const h = length(sub(center, foot));
    const half = (h + radius) / Math.sqrt(dd);
    let lo = tc - half - 1;
    let hi = tc + half + 1;
    const f = (t: number): number =>
        pointBoxSqrDistance(add(line.origin, mul(t, line.direction)), box);
    for (let k = 0; k < 200; ++k) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) <= f(m1)) { hi = m1; }
        else { lo = m0; }
    }
    return Math.sqrt(Math.min(f(lo), f(hi), f(0.5 * (lo + hi))));
}

const nonUnitLine2 = fc.tuple(wellScaledVector(2, -8, 8),
    wellScaledVector(2, -3, 3))
    .filter(([, d]) => length(d) > 0.25)
    .map(([o, d]) => Line.fromOriginDirection(o, d));

const box2 = fc.tuple(wellScaledVector(2, -8, 8), wellScaledVector(2, -8, 8))
    .map(([a, b]) => {
        const lo = Vector.fromArray([Math.min(a.values[0], b.values[0]),
            Math.min(a.values[1], b.values[1])]);
        const hi = Vector.fromArray([Math.max(a.values[0], b.values[0]),
            Math.max(a.values[1], b.values[1])]);
        return AlignedBox.fromMinMax(lo, hi);
    });

describe('DistLine2AlignedBox2 verification', () => {
    const query = new DistLine2AlignedBox2();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(nonUnitLine2, box2), ([line, box]) => {
                const r = query.compute(line, box);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                const diff = sub(r.closest[0], r.closest[1]);
                expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
                expectVectorClose(r.closest[0],
                    add(line.origin, mul(r.parameter, line.direction)), 1e-9,
                    1e-9);
                for (let i = 0; i < 2; ++i) {
                    expect(r.closest[1].values[i])
                        .toBeGreaterThanOrEqual(box.min.values[i] - 1e-9);
                    expect(r.closest[1].values[i])
                        .toBeLessThanOrEqual(box.max.values[i] + 1e-9);
                }
            });
        });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(nonUnitLine2, box2), ([line, box]) => {
            expectClose(query.compute(line, box).distance,
                referenceLineBoxDistance(line, box), 1e-7, 1e-7);
        });
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(nonUnitLine2, box2), ([line, box]) => {
            const o = line.origin.clone();
            const d = line.direction.clone();
            const lo = box.min.clone();
            const hi = box.max.clone();
            query.compute(line, box);
            expectVectorClose(line.origin, o, 0, 0);
            expectVectorClose(line.direction, d, 0, 0);
            expectVectorClose(box.min, lo, 0, 0);
            expectVectorClose(box.max, hi, 0, 0);
        });
    });

    it('is invariant under a common translation', () => {
        check(fc.tuple(nonUnitLine2, box2, wellScaledVector(2, -5, 5)),
            ([line, box, tr]) => {
                const movedLine = Line.fromOriginDirection(
                    add(line.origin, tr), line.direction);
                const movedBox = AlignedBox.fromMinMax(add(box.min, tr),
                    add(box.max, tr));
                expectClose(query.compute(line, box).distance,
                    query.compute(movedLine, movedBox).distance, 1e-8, 1e-8);
            });
    });

    it('is equivariant under the reflections the query applies internally',
        () => {
            // The query reflects the direction into the first quadrant; the
            // distance must be unchanged by reflecting the whole
            // configuration about either coordinate axis.
            check(fc.tuple(nonUnitLine2, box2,
                fc.array(fc.boolean(), { minLength: 2, maxLength: 2 })),
            ([line, box, flip]) => {
                const s = Vector.fromArray([flip[0] ? -1 : 1,
                    flip[1] ? -1 : 1]);
                const refl = (p: Vector): Vector => Vector.fromArray(
                    [s.values[0] * p.values[0], s.values[1] * p.values[1]]);
                const a = refl(box.min);
                const b = refl(box.max);
                const movedBox = AlignedBox.fromMinMax(
                    Vector.fromArray([Math.min(a.values[0], b.values[0]),
                        Math.min(a.values[1], b.values[1])]),
                    Vector.fromArray([Math.max(a.values[0], b.values[0]),
                        Math.max(a.values[1], b.values[1])]));
                const movedLine = Line.fromOriginDirection(refl(line.origin),
                    refl(line.direction));
                expectClose(query.compute(line, box).distance,
                    query.compute(movedLine, movedBox).distance, 0, 0);
            });
        });

    it('reduces to a point-box distance for a zero direction', () => {
        // With D = (0,0) upstream takes the DoQuery0D branch and clamps the
        // line origin to the box.
        check(fc.tuple(wellScaledVector(2, -8, 8), box2), ([p, box]) => {
            const r = query.compute(
                Line.fromOriginDirection(p, new Vector(2)), box);
            expect(r.parameter).toBe(0);
            expectClose(r.sqrDistance, pointBoxSqrDistance(p, box), 1e-12,
                1e-12);
            for (let i = 0; i < 2; ++i) {
                expectClose(r.closest[1].values[i],
                    Math.min(Math.max(p.values[i], box.min.values[i]),
                        box.max.values[i]), 1e-12, 1e-12);
            }
        });
    });

    it('reports zero distance when the line meets the box', () => {
        // Pick a point of the box and a direction through it.
        check(fc.tuple(box2, fc.double({ min: 0, max: 1, noNaN: true }),
            fc.double({ min: 0, max: 1, noNaN: true }), rotationFrame(2),
            wellScaled(-5, 5)), ([box, u, w, R, t]) => {
            const q = Vector.fromArray([
                (1 - u) * box.min.values[0] + u * box.max.values[0],
                (1 - w) * box.min.values[1] + w * box.max.values[1]]);
            const line = Line.fromOriginDirection(
                add(q, mul(t, R[0])), R[0]);
            expect(query.compute(line, box).distance)
                .toBeLessThanOrEqual(1e-9);
        });
    });

    it('is unaffected by the length of the direction vector', () => {
        // Only positive scale factors are covariant in the parameter. A
        // negative factor flips the internal reflections, and when the line
        // meets the box there are infinitely many closest pairs, so upstream
        // may then select the intersection with a different box edge. The
        // distance is unaffected either way.
        check(fc.tuple(nonUnitLine2, box2,
            fc.constantFrom(0.25, 0.5, 1, 2, 4), fc.boolean()),
        ([line, box, k, negate]) => {
            const factor = negate ? -k : k;
            const rescaled = Line.fromOriginDirection(line.origin,
                mul(factor, line.direction));
            const r0 = query.compute(line, box);
            const r1 = query.compute(rescaled, box);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            if (!negate) {
                expectClose(r0.parameter, k * r1.parameter, 1e-9, 1e-8);
            }
        });
    });
});
