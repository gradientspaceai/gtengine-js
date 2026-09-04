import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { Vector, dot, sub } from '../src/Vector.js';
import { check, expectClose, expectVectorClose, fc, positive, seededRandom, wellScaledVector } from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

describe('DistPointCanonicalBox', () => {
    const query = new DistPointCanonicalBox();

    it('reports zero distance for a point inside the box', () => {
        const result = query.compute(v(0.5, -0.5, 0.25), box(1, 1, 1));
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expect(result.closest[1].values).toEqual([0.5, -0.5, 0.25]);
    });

    it('reports zero distance for a point on the boundary', () => {
        const result = query.compute(v(1, 0, 0), box(1, 2, 3));
        expect(result.distance).toBe(0);
    });

    it('measures a face-region distance', () => {
        const result = query.compute(v(4, 0.5, -0.5), box(1, 1, 1));
        expect(result.distance).toBeCloseTo(3, 12);
        expect(result.closest[1].values).toEqual([1, 0.5, -0.5]);
    });

    it('measures an edge-region distance', () => {
        const result = query.compute(v(4, 5, 0), box(1, 1, 1));
        expect(result.distance).toBeCloseTo(5, 12);
        expect(result.closest[1].values).toEqual([1, 1, 0]);
    });

    it('measures a corner-region distance', () => {
        const result = query.compute(v(-3, -5, -13), box(1, 1, 1));
        expect(result.closest[1].values).toEqual([-1, -1, -1]);
        expect(result.sqrDistance).toBeCloseTo(4 + 16 + 144, 12);
    });

    it('works in 2D and 5D', () => {
        const r2 = query.compute(v(3, 4), box(0, 0));
        expect(r2.distance).toBeCloseTo(5, 12);
        const r5 = query.compute(v(2, 2, 2, 2, 2), box(1, 1, 1, 1, 1));
        expect(r5.sqrDistance).toBeCloseTo(5, 12);
    });

    it('leaves the input point untouched in closest[0]', () => {
        const point = v(9, 9, 9);
        const result = query.compute(point, box(1, 1, 1));
        expect(result.closest[0].values).toEqual([9, 9, 9]);
        expect(point.values).toEqual([9, 9, 9]);
    });

    it('agrees with a clamped-coordinate reference', () => {
        let seed = 31337;
        const rand = () => {
            seed = (seed * 1103515245 + 12345) % 2147483648;
            return seed / 2147483648 * 8 - 4;
        };
        const b = box(1.5, 0.5, 2);
        for (let trial = 0; trial < 200; ++trial) {
            const point = v(rand(), rand(), rand());
            const result = query.compute(point, b);
            const clamped = v(
                Math.min(Math.max(point.values[0], -1.5), 1.5),
                Math.min(Math.max(point.values[1], -0.5), 0.5),
                Math.min(Math.max(point.values[2], -2), 2));
            const diff = sub(point, clamped);
            expect(result.sqrDistance).toBeCloseTo(dot(diff, diff), 10);
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i]).toBeCloseTo(
                    clamped.values[i], 12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V19): property-based cross-checks of
// DistPointCanonicalBox.ts against the upstream header
// DistPointCanonicalBox.h.
// ---------------------------------------------------------------------------

const boxOfDim = (n: number): fc.Arbitrary<CanonicalBox> =>
    fc.array(positive(5), { minLength: n, maxLength: n })
        .map(e => CanonicalBox.fromExtent(Vector.fromArray(e)));

describe('DistPointCanonicalBox verification', () => {
    const query = new DistPointCanonicalBox();

    it('matches the closed form sum_i max(0, |p_i| - e_i)^2', () => {
        check(fc.integer({ min: 1, max: 5 }).chain(n =>
            fc.tuple(wellScaledVector(n, -8, 8), boxOfDim(n))),
        ([p, box]) => {
            const r = query.compute(p, box);
            let expected = 0;
            for (let i = 0; i < p.size; ++i) {
                const over = Math.abs(p.values[i]) - box.extent.values[i];
                if (over > 0) { expected += over * over; }
            }
            expectClose(r.sqrDistance, expected, 1e-12, 1e-12);
            expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12, 1e-12);
            const diff = sub(r.closest[0], r.closest[1]);
            expectClose(r.sqrDistance, dot(diff, diff), 1e-12, 1e-12);
        });
    });

    it('returns a closest point inside the box and a copy of the input', () => {
        check(fc.integer({ min: 1, max: 4 }).chain(n =>
            fc.tuple(wellScaledVector(n, -8, 8), boxOfDim(n))),
        ([p, box]) => {
            const r = query.compute(p, box);
            expectVectorClose(r.closest[0], p, 0, 0);
            expect(r.closest[0]).not.toBe(p);
            expect(r.closest[1]).not.toBe(p);
            for (let i = 0; i < p.size; ++i) {
                expect(Math.abs(r.closest[1].values[i]))
                    .toBeLessThanOrEqual(box.extent.values[i] + 1e-12);
            }
        });
    });

    it('reports zero distance exactly for points inside the box', () => {
        check(fc.integer({ min: 1, max: 4 }).chain(n =>
            fc.tuple(fc.array(fc.double({ min: -1, max: 1, noNaN: true }),
                { minLength: n, maxLength: n }), boxOfDim(n))),
        ([frac, box]) => {
            const p = new Vector(box.extent.size);
            for (let i = 0; i < p.size; ++i) {
                p.values[i] = frac[i] * box.extent.values[i];
            }
            const r = query.compute(p, box);
            expect(r.sqrDistance).toBe(0);
            expect(r.distance).toBe(0);
            expectVectorClose(r.closest[1], p, 0, 0);
        });
    });

    it('is minimal over sampled box points', () => {
        const rand = seededRandom(0x51d3);
        check(fc.tuple(wellScaledVector(3, -8, 8), boxOfDim(3)), ([p, box]) => {
            const r = query.compute(p, box);
            const q = new Vector(3);
            for (let k = 0; k < 24; ++k) {
                for (let i = 0; i < 3; ++i) {
                    q.values[i] = box.extent.values[i] * (2 * rand() - 1);
                }
                const diff = sub(p, q);
                const sqr = dot(diff, diff);
                expect(r.sqrDistance).toBeLessThanOrEqual(sqr + 1e-9 * (1 + sqr));
            }
        }, 60);
    });

    it('is equivariant under coordinate reflections', () => {
        check(fc.tuple(wellScaledVector(3, -8, 8), boxOfDim(3),
            fc.array(fc.boolean(), { minLength: 3, maxLength: 3 })),
        ([p, box, flip]) => {
            const q = p.clone();
            for (let i = 0; i < 3; ++i) {
                if (flip[i]) { q.values[i] = -q.values[i]; }
            }
            const r0 = query.compute(p, box);
            const r1 = query.compute(q, box);
            expectClose(r0.sqrDistance, r1.sqrDistance, 1e-12, 1e-12);
            for (let i = 0; i < 3; ++i) {
                const expected = flip[i]
                    ? -r0.closest[1].values[i] : r0.closest[1].values[i];
                expectClose(r1.closest[1].values[i] + 0, expected + 0, 1e-12,
                    1e-12);
            }
        });
    });
});
