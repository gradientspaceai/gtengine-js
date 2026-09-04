import { describe, expect, it } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistLine3AlignedBox3 } from '../src/DistLine3AlignedBox3.js';
import { DistLine3OrientedBox3 } from '../src/DistLine3OrientedBox3.js';
import { Line } from '../src/Line.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function line(origin: number[], direction: number[]): Line {
    return Line.fromOriginDirection(v(...origin), v(...direction));
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistLine3OrientedBox3', () => {
    const query = new DistLine3OrientedBox3();
    const axisAligned = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];

    it('reports zero distance for a line through the box', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0), axisAligned,
            v(1, 1, 1));
        const result = query.compute(line([-5, 0, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('measures a line parallel to a box face', () => {
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0), axisAligned,
            v(1, 1, 1));
        const result = query.compute(line([0, 0, 6], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });

    it('measures a line against a box rotated 45 degrees about z', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(1, 1, 1));
        // A line parallel to the z-axis at (3,0).
        const result = query.compute(line([3, 0, 0], [0, 0, 1]), box);
        expect(result.distance).toBeCloseTo(3 - Math.SQRT2, 9);
    });

    it('handles a degenerate box with zero extents', () => {
        const c = Math.SQRT1_2;
        const box = OrientedBox.fromCenterAxisExtent(v(0, 0, 5),
            [v(c, c, 0), v(-c, c, 0), v(0, 0, 1)], v(0, 0, 0));
        const result = query.compute(line([0, 0, 0], [1, 0, 0]), box);
        expect(result.distance).toBeCloseTo(5, 9);
    });

    it('matches the aligned-box query for an axis-aligned oriented box',
        () => {
            const rnd = makeRandom(60606);
            const abQuery = new DistLine3AlignedBox3();
            const abox = AlignedBox.fromMinMax(v(-1, -2, -3), v(2, 1, 0));
            const { center, extent } = abox.getCenteredForm();
            const obox = OrientedBox.fromCenterAxisExtent(center, axisAligned,
                extent);
            for (let trial = 0; trial < 60; ++trial) {
                const ln = line(
                    [10 * rnd() - 5, 10 * rnd() - 5, 10 * rnd() - 5],
                    [2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1]);
                if (dot(ln.direction, ln.direction) < 1e-4) {
                    continue;
                }
                const r0 = abQuery.compute(ln, abox);
                const r1 = query.compute(ln, obox);
                expect(r1.distance).toBeCloseTo(r0.distance, 9);
            }
        });

    it('agrees with a dense sampling of the box', () => {
        const rnd = makeRandom(70707);
        const c = Math.cos(0.5), s = Math.sin(0.5);
        const axis = [v(c, s, 0), v(-s, c, 0), v(0, 0, 1)];
        const extent = v(1.25, 0.75, 2);
        const box = OrientedBox.fromCenterAxisExtent(v(0.5, -1, 0.25), axis,
            extent);

        for (let trial = 0; trial < 25; ++trial) {
            const origin = v(8 * rnd() - 4, 8 * rnd() - 4, 8 * rnd() - 4);
            const dir = v(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            if (dot(dir, dir) < 1e-4) {
                continue;
            }
            const ln = Line.fromOriginDirection(origin, dir);
            const result = query.compute(ln, box);

            // The reported line point matches the reported parameter.
            const onLine = add(ln.origin, mul(result.parameter, ln.direction));
            for (let i = 0; i < 3; ++i) {
                expect(onLine.values[i]).toBeCloseTo(
                    result.closest[0].values[i], 7);
            }

            // The reported box point is in the box.
            const delta = sub(result.closest[1], box.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(delta, axis[i]))).toBeLessThanOrEqual(
                    extent.values[i] + 1e-8);
            }

            const e = sub(result.closest[0], result.closest[1]);
            expect(Math.sqrt(dot(e, e))).toBeCloseTo(result.distance, 7);

            // No sampled box point is closer to the line.
            const n = 16;
            const dd = dot(ln.direction, ln.direction);
            let best = Number.MAX_VALUE;
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    for (let k = 0; k <= n; ++k) {
                        const q = add(box.center, add(
                            mul((2 * i / n - 1) * extent.values[0], axis[0]),
                            add(
                                mul((2 * j / n - 1) * extent.values[1],
                                    axis[1]),
                                mul((2 * k / n - 1) * extent.values[2],
                                    axis[2]))));
                        const w = sub(q, ln.origin);
                        const t = dot(w, ln.direction) / dd;
                        const f = sub(w, mul(t, ln.direction));
                        best = Math.min(best, dot(f, f));
                    }
                }
            }
            expect(result.sqrDistance).toBeLessThanOrEqual(best + 1e-8);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against the upstream DistLine3OrientedBox3.h.
// ---------------------------------------------------------------------------

describe('DistLine3OrientedBox3 verification', () => {
    const query = new DistLine3OrientedBox3();
    const alignedQuery = new DistLine3AlignedBox3();

    const boxArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0, 4), { minLength: 3, maxLength: 3 }))
        .map(([c, axis, e]) => OrientedBox.fromCenterAxisExtent(c, axis,
            v(e[0], e[1], e[2])));

    const lineArb = fc.tuple(wellScaledVector(3, -8, 8), unitVector(3))
        .map(([o, d]) => Line.fromOriginDirection(o, d));

    function rot3(frame: Vector[], p: Vector): Vector {
        return add(add(mul(p.values[0], frame[0]), mul(p.values[1], frame[1])),
            mul(p.values[2], frame[2]));
    }

    function pointBoxDistance(p: Vector, b: OrientedBox): number {
        const delta = sub(p, b.center);
        let sum = 0;
        for (let i = 0; i < 3; ++i) {
            const over = Math.abs(dot(b.axis[i], delta)) - b.extent.values[i];
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

    // Regression for the upstream defect described in the source comment:
    // DistLine3OrientedBox3.h writes the world-space line point into
    // result.closest[0] and then transforms it a second time as if it were
    // box-frame coordinates. This property fails for the upstream code
    // whenever the box center or axes are not canonical.
    it('reports closest[0] as the line point at the reported parameter', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const r = query.compute(ln, b);
            expectVectorClose(r.closest[0],
                add(ln.origin, mul(r.parameter, ln.direction)), 1e-9, 1e-9);
        });
    });

    it('reproduces the worked example from the source comment', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        const b = OrientedBox.fromCenterAxisExtent(v(10, 0, 0), axes,
            v(1, 1, 1));
        const r = query.compute(line([10, 5, 0], [0, 0, 1]), b);
        expect(r.distance).toBeCloseTo(4, 12);
        expect(r.parameter).toBeCloseTo(1, 12);
        // Upstream would return (20,5,1) for closest[0] here.
        expectVectorClose(r.closest[0], v(10, 5, 1), 1e-12, 1e-12);
        expectVectorClose(r.closest[1], v(10, 1, 1), 1e-12, 1e-12);
    });

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const r = query.compute(ln, b);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-7, 1e-7);
            const delta = sub(r.closest[1], b.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(b.axis[i], delta)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
            }
        });
    });

    it('matches an independent convex minimization along the line', () => {
        check(fc.tuple(lineArb, boxArb), ([ln, b]) => {
            const r = query.compute(ln, b);
            const best = ternaryMin(
                t => pointBoxDistance(add(ln.origin, mul(t, ln.direction)), b),
                -100, 100);
            expectClose(r.distance, best, 1e-7, 1e-7);
        }, 100);
    });

    it('agrees with the aligned-box query for an identity frame', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        check(fc.tuple(lineArb, wellScaledVector(3, -5, 5),
            fc.array(finite(0, 4), { minLength: 3, maxLength: 3 })),
            ([ln, c, e]) => {
                const ext = v(e[0], e[1], e[2]);
                const ob = OrientedBox.fromCenterAxisExtent(c, axes, ext);
                const ab = AlignedBox.fromMinMax(sub(c, ext), add(c, ext));
                const r0 = query.compute(ln, ob);
                const r1 = alignedQuery.compute(ln, ab);
                expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
                expectVectorClose(r0.closest[0], r1.closest[0], 1e-7, 1e-7);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-7, 1e-7);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(lineArb, boxArb, rotationFrame(3),
            wellScaledVector(3, -6, 6)), ([ln, b, frame, shift]) => {
            const movedLine = Line.fromOriginDirection(
                add(shift, rot3(frame, ln.origin)), rot3(frame, ln.direction));
            const movedBox = OrientedBox.fromCenterAxisExtent(
                add(shift, rot3(frame, b.center)),
                [rot3(frame, b.axis[0]), rot3(frame, b.axis[1]),
                    rot3(frame, b.axis[2])], b.extent);
            const r0 = query.compute(ln, b);
            const r1 = query.compute(movedLine, movedBox);
            expectClose(r0.distance, r1.distance, 1e-9, 1e-9);
            expectVectorClose(add(shift, rot3(frame, r0.closest[1])),
                r1.closest[1], 1e-6, 1e-6);
        });
    });
});
