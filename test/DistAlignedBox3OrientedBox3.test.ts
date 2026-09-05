import { describe, it, expect } from 'vitest';
import { AlignedBox } from '../src/AlignedBox.js';
import { DistAlignedBox3OrientedBox3 } from '../src/DistAlignedBox3OrientedBox3.js';
import { DistOrientedBox3OrientedBox3 } from '../src/DistOrientedBox3OrientedBox3.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, length, mul, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';
import { cross } from '../src/Vector3.js';

function v(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function frame(w: Vector): [Vector, Vector, Vector] {
    const n = w.clone();
    normalize(n);
    let u = Math.abs(n.values[0]) > 0.5
        ? v(-n.values[1], n.values[0], 0)
        : v(0, -n.values[2], n.values[1]);
    normalize(u);
    return [u, cross(n, u), n];
}

function distPointOrientedBox(p: Vector, b: OrientedBox): number {
    const delta = sub(p, b.center);
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const y = dot(delta, b.axis[i]);
        const e = b.extent.values[i];
        const d = y < -e ? -e - y : (y > e ? y - e : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

function distPointAlignedBox(p: Vector, b: AlignedBox): number {
    let sqr = 0;
    for (let i = 0; i < 3; ++i) {
        const x = p.values[i];
        const lo = b.min.values[i], hi = b.max.values[i];
        const d = x < lo ? lo - x : (x > hi ? x - hi : 0);
        sqr += d * d;
    }
    return Math.sqrt(sqr);
}

// Brute-force minimum from the surface of the aligned box to the solid
// oriented box.
function bruteForce(box0: AlignedBox, box1: OrientedBox, n: number): number {
    let best = Infinity;
    for (let face = 0; face < 3; ++face) {
        const i0 = face, i1 = (face + 1) % 3, i2 = (face + 2) % 3;
        for (const which of [box0.min, box0.max]) {
            for (let i = 0; i <= n; ++i) {
                for (let j = 0; j <= n; ++j) {
                    const p = new Vector(3);
                    p.values[i0] = which.values[i0];
                    p.values[i1] = box0.min.values[i1] + i / n *
                        (box0.max.values[i1] - box0.min.values[i1]);
                    p.values[i2] = box0.min.values[i2] + j / n *
                        (box0.max.values[i2] - box0.min.values[i2]);
                    const d = distPointOrientedBox(p, box1);
                    if (d < best) {
                        best = d;
                    }
                }
            }
        }
    }
    return best;
}

describe('DistAlignedBox3OrientedBox3', () => {
    it('reports the gap between an aligned box and an aligned oriented box',
        () => {
            const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
            const box1 = OrientedBox.fromCenterAxisExtent(v(5, 0, 0),
                [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 1, 1));
            const result =
                new DistAlignedBox3OrientedBox3().compute(box0, box1);
            expect(result.distance).toBeCloseTo(3, 12);
            expect(result.closest[0].values[0]).toBeCloseTo(1, 12);
            expect(result.closest[1].values[0]).toBeCloseTo(4, 12);
        });

    it('handles an aligned box that is not centered at the origin', () => {
        const box0 = AlignedBox.fromMinMax(v(10, 20, 30), v(12, 24, 36));
        const box1 = OrientedBox.fromCenterAxisExtent(v(11, 22, 40),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(0.5, 0.5, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[0].values[2]).toBeCloseTo(36, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(39, 12);
        expect(distPointAlignedBox(result.closest[0], box0)).toBeCloseTo(0, 8);
    });

    it('accounts for the rotation of the oriented box', () => {
        const s = Math.SQRT1_2;
        const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
        const box1 = OrientedBox.fromCenterAxisExtent(v(5, 0, 0),
            [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 12);
    });

    it('reports zero when the boxes overlap', () => {
        const s = Math.SQRT1_2;
        const box0 = AlignedBox.fromMinMax(v(-1, -1, -1), v(1, 1, 1));
        const box1 = OrientedBox.fromCenterAxisExtent(v(1, 1, 0),
            [v(s, s, 0), v(-s, s, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
    });

    it('matches the box-box query on the equivalent oriented box', () => {
        let seed = 246813579;
        const rnd = () => {
            seed = (seed * 1103515245 + 12345) & 0x7fffffff;
            return seed / 0x7fffffff;
        };
        const query = new DistAlignedBox3OrientedBox3();
        const reference = new DistOrientedBox3OrientedBox3();
        for (let trial = 0; trial < 25; ++trial) {
            const lo = v(rnd() * 4 - 3, rnd() * 4 - 3, rnd() * 4 - 3);
            const box0 = AlignedBox.fromMinMax(lo,
                add(lo, v(0.5 + rnd() * 2, 0.5 + rnd() * 2, 0.5 + rnd() * 2)));

            const w = v(rnd() * 2 - 1, rnd() * 2 - 1, rnd() * 2 - 1);
            if (length(w) < 1e-3) {
                continue;
            }
            const [a0, a1, a2] = frame(w);
            const box1 = OrientedBox.fromCenterAxisExtent(
                v(rnd() * 6 - 3, rnd() * 6 - 3, rnd() * 6 - 3), [a0, a1, a2],
                v(0.3 + rnd(), 0.3 + rnd(), 0.3 + rnd()));

            const result = query.compute(box0, box1);

            // The same query, expressed with box0 converted by hand.
            const centered = box0.getCenteredForm();
            const obox0 = OrientedBox.fromCenterAxisExtent(centered.center,
                [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], centered.extent);
            const expected = reference.compute(obox0, box1);
            expect(result.distance).toBe(expected.distance);

            // Cross-check against brute-force sampling of the aligned box.
            const brute = bruteForce(box0, box1, 24);
            expect(result.distance).toBeLessThanOrEqual(brute + 1e-9);
            expect(result.distance).toBeGreaterThan(brute - 0.25);

            expect(length(sub(result.closest[0], result.closest[1])))
                .toBeCloseTo(result.distance, 8);
            expect(distPointAlignedBox(result.closest[0], box0))
                .toBeCloseTo(0, 8);
            expect(distPointOrientedBox(result.closest[1], box1))
                .toBeCloseTo(0, 8);
        }
    });

    it('handles a degenerate aligned box that is a single point', () => {
        // A degenerate aligned box (a point) still produces a valid query.
        const box0 = AlignedBox.fromMinMax(v(2, 2, 2), v(2, 2, 2));
        const box1 = OrientedBox.fromCenterAxisExtent(v(0, 0, 0),
            [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)], v(1, 1, 1));
        const result = new DistAlignedBox3OrientedBox3().compute(box0, box1);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 12);
        expect(result.closest[0].values).toEqual([2, 2, 2]);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(1, 12);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistAlignedBox3OrientedBox3.h.
// ---------------------------------------------------------------------------

describe('DistAlignedBox3OrientedBox3 verification', () => {
    const query = new DistAlignedBox3OrientedBox3();
    const bbQuery = new DistOrientedBox3OrientedBox3();

    const abArb = fc.tuple(wellScaledVector(3, -4, 4),
        fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 }))
        .map(([c, e]) => AlignedBox.fromMinMax(
            v(c.values[0] - e[0], c.values[1] - e[1], c.values[2] - e[2]),
            v(c.values[0] + e[0], c.values[1] + e[1], c.values[2] + e[2])));

    const obArb = fc.tuple(wellScaledVector(3, -5, 5), rotationFrame(3),
        fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 }))
        .map(([c, R, e]) => OrientedBox.fromCenterAxisExtent(c, R,
            Vector.fromArray([e[0], e[1], e[2]])));

    it('reports consistent distances and in-box closest points', () => {
        check(fc.tuple(abArb, obArb), ([a, b]) => {
            const r = query.compute(a, b);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            // The absolute tolerance is 1e-6: these queries accumulate the
            // squared distance while clamping to faces and edges, so a
            // near-touching configuration loses about half the mantissa and
            // the distance carries an absolute error of order sqrt(eps)
            // times the coordinate scale. A translation or frame error
            // would show up as an O(1) discrepancy.
            expectClose(length(sub(r.closest[0], r.closest[1])), r.distance,
                1e-6, 1e-8);
            for (let i = 0; i < 3; ++i) {
                expect(r.closest[0].values[i])
                    .toBeGreaterThanOrEqual(a.min.values[i] - 1e-8);
                expect(r.closest[0].values[i])
                    .toBeLessThanOrEqual(a.max.values[i] + 1e-8);
            }
            const d = sub(r.closest[1], b.center);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(dot(b.axis[i], d)))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-8);
            }
        });
    });

    // The upstream body converts the aligned box to an identity-framed
    // oriented box with center (max+min)/2 and extent (max-min)/2 and
    // delegates. A swapped center/extent or a wrong factor breaks this.
    it('is the box-box query on the centered form of the aligned box', () => {
        const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
        check(fc.tuple(abArb, obArb), ([a, b]) => {
            const cf = a.getCenteredForm();
            const ob = OrientedBox.fromCenterAxisExtent(cf.center, axes,
                cf.extent);
            const r0 = query.compute(a, b);
            const r1 = bbQuery.compute(ob, b);
            expectClose(r0.distance, r1.distance, 1e-12, 1e-12);
            expectVectorClose(r0.closest[0], r1.closest[0], 1e-12, 1e-12);
            expectVectorClose(r0.closest[1], r1.closest[1], 1e-12, 1e-12);
        });
    });

    // With an identity frame on the oriented box the answer is the Euclidean
    // norm of the per-axis interval gaps.
    it('matches the interval-gap formula for an identity oriented frame',
        () => {
            const axes = [v(1, 0, 0), v(0, 1, 0), v(0, 0, 1)];
            check(fc.tuple(abArb, wellScaledVector(3, -5, 5),
                fc.array(finite(0.05, 3), { minLength: 3, maxLength: 3 })),
                ([a, c, e]) => {
                    const b = OrientedBox.fromCenterAxisExtent(c, axes,
                        Vector.fromArray(e));
                    let sqr = 0;
                    for (let i = 0; i < 3; ++i) {
                        const lo = Math.max(a.min.values[i],
                            c.values[i] - e[i]);
                        const hi = Math.min(a.max.values[i],
                            c.values[i] + e[i]);
                        if (lo > hi) { sqr += (lo - hi) * (lo - hi); }
                    }
                    const r = query.compute(a, b);
                    expectClose(r.distance, Math.sqrt(sqr), 1e-8, 1e-8);
                });
        });

    it('never exceeds a face sampling of the aligned box', () => {
        const rng = seededRandom(0x4d1c3b7a);
        for (let k = 0; k < 20; ++k) {
            const c = v(6 * rng() - 3, 6 * rng() - 3, 6 * rng() - 3);
            const e = v(0.3 + rng(), 0.3 + rng(), 0.3 + rng());
            const a = AlignedBox.fromMinMax(sub(c, e), add(c, e));
            const w = v(rng() - 0.5, rng() - 0.5, rng() - 0.5);
            normalize(w);
            const f = frame(w);
            const b = OrientedBox.fromCenterAxisExtent(
                v(6 * rng() - 3, 6 * rng() - 3, 6 * rng() - 3),
                [f[0], f[1], f[2]],
                v(0.3 + rng(), 0.3 + rng(), 0.3 + rng()));
            expect(query.compute(a, b).distance)
                .toBeLessThanOrEqual(bruteForce(a, b, 24) + 1e-9);
        }
    }, 30000);

    it('is invariant under a common translation', () => {
        check(fc.tuple(abArb, obArb, wellScaledVector(3, -5, 5)),
            ([a, b, t]) => {
                const r0 = query.compute(a, b);
                const r1 = query.compute(
                    AlignedBox.fromMinMax(add(a.min, t), add(a.max, t)),
                    OrientedBox.fromCenterAxisExtent(add(b.center, t), b.axis,
                        b.extent));
                // Only the distance is compared: when the boxes touch or
                // overlap there are many equidistant pairs and the two runs
                // may name different representatives.
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectClose(length(sub(r1.closest[0], r1.closest[1])),
                    r0.distance, 1e-7, 1e-7);
            });
    });
});
