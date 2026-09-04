import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPlane3CanonicalBox3 } from '../src/DistPlane3CanonicalBox3.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import {
    check, expectClose, fc, positive, unitVector, wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function box(...extent: number[]): CanonicalBox {
    return CanonicalBox.fromExtent(v(...extent));
}

function plane(normal: number[], constant: number): Hyperplane {
    return Hyperplane.fromNormalConstant(v(...normal), constant);
}

function unitPlane(normal: number[], constant: number): Hyperplane {
    const n = v(...normal);
    normalize(n);
    return Hyperplane.fromNormalConstant(n, constant);
}

// The exact plane-box distance for a unit-length normal. The box is centered
// at the origin, so Dot(N,X) ranges over [-R,R] with R = sum_i e[i]*|n[i]|.
function exactDistance(p: Hyperplane, b: CanonicalBox): number {
    let radius = 0;
    for (let i = 0; i < 3; ++i) {
        radius += b.extent.values[i] * Math.abs(p.normal.values[i]);
    }
    return Math.max(0, Math.abs(p.constant) - radius);
}

function verifyClosest(p: Hyperplane, b: CanonicalBox,
    result: { distance: number, closest: [Vector, Vector] }): void {
    // closest[0] is on the plane.
    expect(dot(p.normal, result.closest[0])).toBeCloseTo(p.constant, 9);
    // closest[1] is in the box.
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(result.closest[1].values[i]))
            .toBeLessThanOrEqual(b.extent.values[i] + 1e-9);
    }
    expect(length(sub(result.closest[0], result.closest[1])))
        .toBeCloseTo(result.distance, 9);
}

describe('DistPlane3CanonicalBox3', () => {
    const query = new DistPlane3CanonicalBox3();
    const unitBox = box(1, 1, 1);

    it('computes the distance for an axis-aligned plane', () => {
        const p = plane([0, 0, 1], 2);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(1, 12);
        expect(result.sqrDistance).toBeCloseTo(1, 12);
        verifyClosest(p, unitBox, result);
    });

    it('reports zero distance when the plane cuts the box', () => {
        const p = plane([0, 0, 1], 0.5);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
        verifyClosest(p, unitBox, result);
    });

    it('handles a plane touching a box face', () => {
        const p = plane([1, 0, 0], 1);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(0, 12);
    });

    it('computes the distance for an oblique plane', () => {
        // The unit normal (1,1,1)/sqrt(3) has box radius sqrt(3).
        const p = unitPlane([1, 1, 1], 3);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(3 - Math.sqrt(3), 10);
        verifyClosest(p, unitBox, result);
    });

    it('handles a normal with a negative component (reflection path)', () => {
        const p = unitPlane([-1, 1, 0], 4);
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(4 - Math.SQRT2, 10);
        verifyClosest(p, unitBox, result);
    });

    it('handles the one-nonzero-component normal path', () => {
        const b = box(1, 2, 3);
        const p = plane([1, 0, 0], 5);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(4, 12);
        // The upstream 1D path reports the plane point (p0,e1,e2) and the box
        // point (clamp(p0),e1,e2).
        expect(result.closest[0].values).toEqual([5, 2, 3]);
        expect(result.closest[1].values).toEqual([1, 2, 3]);
    });

    it('handles a degenerate zero normal', () => {
        // A zero normal is degenerate; upstream falls back to a point-box
        // query using the stored plane origin.
        const p = Hyperplane.fromNormalOrigin(v(0, 0, 0), v(3, 0, 0));
        const result = query.compute(p, unitBox);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[0].values).toEqual([3, 0, 0]);
        expect(result.closest[1].values).toEqual([1, 0, 0]);
    });

    it('handles a degenerate box with zero extents', () => {
        const b = box(0, 0, 0);
        const p = plane([0, 0, 1], 2);
        const result = query.compute(p, b);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values).toEqual([0, 0, 0]);
    });

    it('agrees with the exact formula on random inputs', () => {
        let seed = 1029384756;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 400; ++trial) {
            const n = v(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const p = Hyperplane.fromNormalConstant(n, 6 * rand() - 3);
            const b = box(0.2 + 2 * rand(), 0.2 + 2 * rand(),
                0.2 + 2 * rand());
            const result = query.compute(p, b);
            expect(result.distance).toBeCloseTo(exactDistance(p, b), 9);
            verifyClosest(p, b, result);
        }
    });

    it('is symmetric under negation of the plane', () => {
        let seed = 5647382910 % 2147483648;
        const rand = (): number => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };
        for (let trial = 0; trial < 100; ++trial) {
            const n = v(2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1);
            if (dot(n, n) < 1e-6) {
                continue;
            }
            normalize(n);
            const c = 6 * rand() - 3;
            const b = box(0.5 + rand(), 0.5 + rand(), 0.5 + rand());
            const r0 = query.compute(Hyperplane.fromNormalConstant(n, c), b);
            const negated = v(-n.values[0], -n.values[1], -n.values[2]);
            const r1 = query.compute(
                Hyperplane.fromNormalConstant(negated, -c), b);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistPlane3CanonicalBox3.h.
// ---------------------------------------------------------------------------

// A canonical box with strictly positive extents.
const boxArb = fc.array(positive(5, 1e-2), { minLength: 3, maxLength: 3 })
    .map(e => CanonicalBox.fromExtent(Vector.fromArray(e)));

// A plane with a unit-length normal (the query documents that requirement)
// and a moderate constant, so the plane sometimes misses and sometimes cuts
// the box.
const planeArb = fc.tuple(unitVector(3), wellScaled(-8, 8))
    .map(([n, c]) => Hyperplane.fromNormalConstant(n, c));

// The exact plane-box distance for a unit-length normal. Over the box the
// linear functional Dot(N,X) sweeps [-R,R] with R = sum_i e[i]*|n[i]|, so the
// distance to the level set Dot(N,X) = c is max(0, |c| - R).
function exactPlaneBoxDistance(p: Hyperplane, b: CanonicalBox): number {
    let radius = 0;
    for (let i = 0; i < 3; ++i) {
        radius += b.extent.values[i] * Math.abs(p.normal.values[i]);
    }
    return Math.max(0, Math.abs(p.constant) - radius);
}

describe('DistPlane3CanonicalBox3 verification', () => {
    const query = new DistPlane3CanonicalBox3();

    it('result is self consistent and the points lie on their primitives',
        () => {
            check(fc.tuple(planeArb, boxArb), ([p, b]) => {
                const r = query.compute(p, b);
                expectClose(r.distance, Math.sqrt(r.sqrDistance), 1e-12,
                    1e-12);
                expectClose(r.distance,
                    length(sub(r.closest[0], r.closest[1])), 1e-9, 1e-9);
                // closest[0] is on the plane. Every DoQuery branch reached
                // with a unit-length normal produces an exact plane point.
                expectClose(dot(p.normal, r.closest[0]), p.constant, 1e-9,
                    1e-9);
                // closest[1] is in the box.
                for (let i = 0; i < 3; ++i) {
                    expect(Math.abs(r.closest[1].values[i]))
                        .toBeLessThanOrEqual(b.extent.values[i] + 1e-12);
                }
            });
        });

    it('matches the exact plane-box distance formula', () => {
        check(fc.tuple(planeArb, boxArb), ([p, b]) => {
            expectClose(query.compute(p, b).distance,
                exactPlaneBoxDistance(p, b), 1e-9, 1e-9);
        });
    });

    it('is not larger than the distance to any sampled box point', () => {
        check(fc.tuple(planeArb, boxArb,
            fc.array(wellScaledVector(3, -1, 1), { minLength: 6,
                maxLength: 6 })), ([p, b, dirs]) => {
            const d = query.compute(p, b).distance;
            for (const u of dirs) {
                const x = Vector.fromArray([
                    Math.max(-1, Math.min(1, u.values[0])) * b.extent.values[0],
                    Math.max(-1, Math.min(1, u.values[1])) * b.extent.values[1],
                    Math.max(-1, Math.min(1, u.values[2])) * b.extent.values[2]
                ]);
                // Distance from the sampled box point to the plane.
                const sampled = Math.abs(dot(p.normal, x) - p.constant);
                expect(d).toBeLessThanOrEqual(sampled + 1e-9);
            }
        });
    });

    it('is equivariant under the reflections the query applies internally',
        () => {
            // The query reflects the normal into the first octant. Reflecting
            // the whole configuration about any coordinate plane leaves the
            // canonical box fixed, so the distance must be unchanged.
            check(fc.tuple(planeArb, boxArb,
                fc.array(fc.boolean(), { minLength: 3, maxLength: 3 })),
            ([p, b, flip]) => {
                const s = flip.map(f => (f ? -1 : 1));
                const n = Vector.fromArray([s[0] * p.normal.values[0],
                    s[1] * p.normal.values[1], s[2] * p.normal.values[2]]);
                const moved = Hyperplane.fromNormalConstant(n, p.constant);
                expectClose(query.compute(p, b).distance,
                    query.compute(moved, b).distance, 0, 0);
            });
        });

    it('does not mutate its inputs', () => {
        check(fc.tuple(planeArb, boxArb), ([p, b]) => {
            const n = p.normal.clone();
            const c = p.constant;
            const e = b.extent.clone();
            const r = query.compute(p, b);
            expect(p.normal.values).toEqual(n.values);
            expect(p.constant).toBe(c);
            expect(b.extent.values).toEqual(e.values);
            // The result must not alias the box extent either.
            r.closest[1].values[0] = 12345;
            expect(b.extent.values).toEqual(e.values);
        });
    });

    it('handles axis-aligned normals (the DoQuery1D branches)', () => {
        // These are the branches upstream reaches when two normal components
        // are zero; they are essentially unreachable from random unit
        // normals, so exercise them explicitly.
        check(fc.tuple(fc.integer({ min: 0, max: 2 }), fc.boolean(),
            wellScaled(-8, 8), boxArb), ([axis, negate, c, b]) => {
            const n = new Vector(3);
            n.values[axis] = negate ? -1 : 1;
            const p = Hyperplane.fromNormalConstant(n, c);
            const r = query.compute(p, b);
            expectClose(r.distance,
                Math.max(0, Math.abs(c) - b.extent.values[axis]), 1e-12,
                1e-12);
            expectClose(dot(n, r.closest[0]), c, 1e-12, 1e-12);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(r.closest[1].values[i]))
                    .toBeLessThanOrEqual(b.extent.values[i] + 1e-12);
            }
        });
    });

    it('reports a zero-normal plane through the box origin query (DoQuery0D)',
        () => {
            // With N = (0,0,0) upstream falls through to the point-box query
            // for the stored plane origin, which is (0,0,0) for a zero
            // normal. The closest box point is then the origin itself.
            check(boxArb, b => {
                const p = Hyperplane.fromNormalConstant(new Vector(3), 3);
                const r = query.compute(p, b);
                expect(r.distance).toBe(0);
                expect(r.closest[0].values.map(x => x + 0)).toEqual([0, 0, 0]);
                expect(r.closest[1].values.map(x => x + 0)).toEqual([0, 0, 0]);
            });
        });

    it('reduces to |constant| for a degenerate zero-extent box', () => {
        check(planeArb, p => {
            const b = CanonicalBox.fromExtent(new Vector(3));
            const r = query.compute(p, b);
            expectClose(r.distance, Math.abs(p.constant), 1e-9, 1e-9);
            // '+ 0' normalizes the -0 that clamp(x, -0, +0) returns for
            // negative x; toEqual uses Object.is, which separates -0 from 0.
            expect(r.closest[1].values.map(x => x + 0)).toEqual([0, 0, 0]);
        });
    });
});
