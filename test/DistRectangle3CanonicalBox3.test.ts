import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistRectangle3CanonicalBox3 }
    from '../src/DistRectangle3CanonicalBox3.js';
import type { DistRectangle3CanonicalBox3Result }
    from '../src/DistRectangle3CanonicalBox3.js';
import { Rectangle } from '../src/Rectangle.js';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistRectangle3AlignedBox3 } from '../src/DistRectangle3AlignedBox3.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

// An orthonormal frame parameterized by two angles; the first two vectors
// are used as the rectangle axes.
function frame(a: number, b: number): Vector[] {
    const ca = Math.cos(a), sa = Math.sin(a);
    const cb = Math.cos(b), sb = Math.sin(b);
    return [
        v(ca, sa, 0),
        v(-sa * cb, ca * cb, sb),
        v(sa * sb, -ca * sb, cb)
    ];
}

function rect(center: number[], axis: Vector[], e0: number, e1: number):
    Rectangle {
    return Rectangle.fromCenterAxisExtent(v(...center), [axis[0], axis[1]],
        v(e0, e1));
}

// The exact distance from a point to a solid canonical box.
function pointBoxDistance(p: Vector, box: CanonicalBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const d = Math.abs(p.values[i]) - box.extent.values[i];
        if (d > 0) {
            sqr += d * d;
        }
    }
    return Math.sqrt(sqr);
}

// A grid sampling of the solid rectangle. The distance from each sample to
// the box is exact, so the result is an upper bound whose error is at most
// the grid spacing.
function sampledDistance(r: Rectangle, box: CanonicalBox, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i = 0; i <= n; ++i) {
        const s0 = r.extent.values[0] * (-1 + (2 * i) / n);
        const base = add(r.center, mul(s0, r.axis[0]));
        for (let j = 0; j <= n; ++j) {
            const s1 = r.extent.values[1] * (-1 + (2 * j) / n);
            best = Math.min(best,
                pointBoxDistance(add(base, mul(s1, r.axis[1])), box));
        }
    }
    return best;
}

function expectConsistent(r: Rectangle, box: CanonicalBox,
    result: DistRectangle3CanonicalBox3Result): void {
    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);

    // The cartesian coordinates are within the extents and reproduce
    // closest[0].
    for (let i = 0; i < 2; ++i) {
        expect(Math.abs(result.cartesian[i]))
            .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
    }
    const reconstructed = add(r.center,
        add(mul(result.cartesian[0], r.axis[0]),
            mul(result.cartesian[1], r.axis[1])));
    for (let i = 0; i < 3; ++i) {
        expect(reconstructed.values[i])
            .toBeCloseTo(result.closest[0].values[i], 7);
    }

    // closest[1] is in the box.
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(result.closest[1].values[i]))
            .toBeLessThanOrEqual(box.extent.values[i] + 1e-8);
    }
}

describe('DistRectangle3CanonicalBox3', () => {
    const query = new DistRectangle3CanonicalBox3();
    const unit = CanonicalBox.fromExtent(v(1, 1, 1));

    it('computes the distance for a rectangle parallel to a box face', () => {
        const r = rect([0, 0, 5], frame(0, 0), 2, 3);
        const result = query.compute(r, unit);
        expect(result.distance).toBeCloseTo(4, 10);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 10);
        // There are infinitely many closest-point pairs here; the query
        // returns one of them, so only consistency is asserted.
        expect(result.closest[0].values[2]).toBeCloseTo(5, 10);
        expectConsistent(r, unit, result);
    });

    it('reports zero distance when the rectangle cuts the box', () => {
        const r = rect([0, 0, 0], frame(0, 0), 4, 4);
        const result = query.compute(r, unit);
        expect(result.distance).toBeCloseTo(0, 10);
        expectConsistent(r, unit, result);
    });

    it('uses a rectangle edge when the plane point is outside', () => {
        // The rectangle lies in the plane z = 3 but is offset in x, so the
        // closest rectangle point is on its x = -e0 edge.
        const r = rect([6, 0, 3], frame(0, 0), 1, 1);
        const result = query.compute(r, unit);
        expect(result.cartesian[0]).toBeCloseTo(-1, 8);
        expect(result.distance).toBeCloseTo(Math.sqrt(16 + 4), 8);
        expectConsistent(r, unit, result);
    });

    it('computes the distance for a rectangle beyond a box vertex', () => {
        const r = rect([4, 4, 4], frame(0, 0), 0.25, 0.25);
        const result = query.compute(r, unit);
        // The closest rectangle point is its corner (3.75,3.75,4) and the
        // closest box point is (1,1,1).
        const expected = Math.sqrt(2.75 * 2.75 + 2.75 * 2.75 + 3 * 3);
        expect(result.distance).toBeCloseTo(expected, 8);
        expectConsistent(r, unit, result);
    });

    it('matches a dense sampling of the rectangle on random inputs', () => {
        let seed = 246813579;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = CanonicalBox.fromExtent(v(1, 1.5, 0.75));
        for (let trial = 0; trial < 25; ++trial) {
            const r = rect(
                [rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4],
                frame(rand() * Math.PI, rand() * Math.PI),
                rand() * 2 + 0.1, rand() * 2 + 0.1);
            const result = query.compute(r, box);
            const sampled = sampledDistance(r, box, 60);
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-8);
            expect(sampled - result.distance).toBeLessThan(0.15);
            expectConsistent(r, box, result);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistRectangle3CanonicalBox3.h.
// ---------------------------------------------------------------------------

describe('DistRectangle3CanonicalBox3 verification', () => {
    const query = new DistRectangle3CanonicalBox3();
    const abQuery = new DistRectangle3AlignedBox3();

    const boxArb = fc.array(finite(0.05, 4), { minLength: 3, maxLength: 3 })
        .map(e => CanonicalBox.fromExtent(v(e[0], e[1], e[2])));

    const rectArb = fc.tuple(wellScaledVector(3, -6, 6), rotationFrame(3),
        fc.array(finite(0.05, 4), { minLength: 2, maxLength: 2 }))
        .map(([c, R, e]) =>
            Rectangle.fromCenterAxisExtent(c, [R[0], R[1]], v(e[0], e[1])));

    it('reports consistent distances and on-primitive closest points', () => {
        check(fc.tuple(rectArb, boxArb), ([r, box]) => {
            const res = query.compute(r, box);
            expectClose(res.sqrDistance, res.distance * res.distance,
                1e-12, 1e-12);
            const d = sub(res.closest[0], res.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), res.distance, 1e-8, 1e-8);

            // The W-coordinates are in range and reproduce closest[0].
            let rebuilt = r.center.clone();
            for (let i = 0; i < 2; ++i) {
                expect(Math.abs(res.cartesian[i]))
                    .toBeLessThanOrEqual(r.extent.values[i] + 1e-8);
                rebuilt = add(rebuilt, mul(res.cartesian[i], r.axis[i]));
            }
            expectVectorClose(rebuilt, res.closest[0], 1e-7, 1e-7);

            // closest[1] is in the box.
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(res.closest[1].values[i]))
                    .toBeLessThanOrEqual(box.extent.values[i] + 1e-8);
            }
        });
    });

    // The distance from the rectangle to the box equals the distance from
    // closest[0] to the box (an exact, independent computation), which
    // catches any mismatch between the reported distance and the reported
    // closest points.
    it('agrees with the exact point-box distance at closest[0]', () => {
        check(fc.tuple(rectArb, boxArb), ([r, box]) => {
            const res = query.compute(r, box);
            expectClose(res.distance, pointBoxDistance(res.closest[0], box),
                1e-7, 1e-7);
        });
    });

    it('never exceeds a grid sampling of the rectangle', () => {
        const rng = seededRandom(0x2468ace0);
        const box = CanonicalBox.fromExtent(v(1, 2, 0.5));
        for (let k = 0; k < 40; ++k) {
            const a = 2 * Math.PI * rng(), b = 2 * Math.PI * rng();
            const r = Rectangle.fromCenterAxisExtent(
                v(8 * rng() - 4, 8 * rng() - 4, 8 * rng() - 4),
                [frame(a, b)[0], frame(a, b)[1]],
                v(0.5 + 2 * rng(), 0.5 + 2 * rng()));
            const res = query.compute(r, box);
            expect(res.distance)
                .toBeLessThanOrEqual(sampledDistance(r, box, 24) + 1e-9);
        }
    }, 30000);

    // DistRectangle3AlignedBox3 is "translate, run this query, translate
    // back", so the two must agree for a box centered at the origin.
    it('agrees with the aligned-box query for an origin-centered box', () => {
        check(fc.tuple(rectArb, boxArb), ([r, box]) => {
            const e = box.extent;
            const ab = AlignedBox.fromMinMax(
                v(-e.values[0], -e.values[1], -e.values[2]),
                v(e.values[0], e.values[1], e.values[2]));
            const r0 = query.compute(r, box);
            const r1 = abQuery.compute(r, ab);
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectVectorClose(r0.closest[0], r1.closest[0], 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-12, 1e-12);
        });
    });

    it('maps the cartesian coordinates of the four edges correctly', () => {
        // A rectangle far away and parallel to the xy-plane, positioned so
        // the closest point is a specific corner. This pins down the
        // sign/j0/j1 tables used by the edge branch.
        const axes = [v(1, 0, 0), v(0, 1, 0)];
        const box = CanonicalBox.fromExtent(v(1, 1, 1));
        for (const [cx, cy, sx, sy] of [[10, 10, 1, 1], [-10, 10, -1, 1],
            [10, -10, 1, -1], [-10, -10, -1, -1]]) {
            const r = Rectangle.fromCenterAxisExtent(v(cx, cy, 5), axes,
                v(2, 3));
            const res = query.compute(r, box);
            expect(res.cartesian[0]).toBeCloseTo(-sx * 2, 9);
            expect(res.cartesian[1]).toBeCloseTo(-sy * 3, 9);
            expectVectorClose(res.closest[0],
                v(cx - sx * 2, cy - sy * 3, 5), 1e-9, 1e-9);
        }
    });

    it('degenerates to the point-box distance for a zero-extent rectangle',
        () => {
            check(fc.tuple(wellScaledVector(3, -6, 6), rotationFrame(3),
                boxArb), ([c, R, box]) => {
                    const r = Rectangle.fromCenterAxisExtent(c, [R[0], R[1]],
                        v(0, 0));
                    const res = query.compute(r, box);
                    expectClose(res.distance, pointBoxDistance(c, box),
                        1e-8, 1e-8);
                    expectVectorClose(res.closest[0], c, 1e-9, 1e-9);
                });
        });
});
