import { describe, expect, it } from 'vitest';
import { CanonicalBox } from '../src/CanonicalBox.js';
import { DistPoint3Parallelepiped3 } from '../src/DistPoint3Parallelepiped3.js';
import type { DistPoint3Parallelepiped3Result }
    from '../src/DistPoint3Parallelepiped3.js';
import { DistPointCanonicalBox } from '../src/DistPointCanonicalBox.js';
import { Matrix, mulMatrix } from '../src/Matrix.js';
import { inverse3x3 } from '../src/Matrix3x3.js';
import { Parallelepiped3 } from '../src/Parallelepiped3.js';
import { DistPointOrientedBox } from '../src/DistPointOrientedBox.js';
import { OrientedBox } from '../src/OrientedBox.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function ppd(center: number[], a0: number[], a1: number[], a2: number[]):
    Parallelepiped3 {
    return Parallelepiped3.fromCenterAxis(v(...center),
        [v(...a0), v(...a1), v(...a2)]);
}

// The parallelepiped coordinates of a point: X = C + sum_i s[i] * A[i].
function coordinates(p: Vector, box: Parallelepiped3): Vector {
    const B = new Matrix(3, 3);
    B.setCol(0, box.axis[0]);
    B.setCol(1, box.axis[1]);
    B.setCol(2, box.axis[2]);
    return mulMatrix(inverse3x3(B).inverse, sub(p, box.center));
}

// A grid sampling of the solid parallelepiped gives an upper bound for the
// distance that converges to the true value as n grows.
function sampledDistance(p: Vector, box: Parallelepiped3, n: number): number {
    let best = Number.MAX_VALUE;
    for (let i0 = 0; i0 <= n; ++i0) {
        const s0 = -1 + (2 * i0) / n;
        const t0 = add(box.center, mul(s0, box.axis[0]));
        for (let i1 = 0; i1 <= n; ++i1) {
            const s1 = -1 + (2 * i1) / n;
            const t1 = add(t0, mul(s1, box.axis[1]));
            for (let i2 = 0; i2 <= n; ++i2) {
                const s2 = -1 + (2 * i2) / n;
                const d = sub(p, add(t1, mul(s2, box.axis[2])));
                best = Math.min(best, dot(d, d));
            }
        }
    }
    return Math.sqrt(best);
}

function expectConsistent(p: Vector, box: Parallelepiped3,
    result: DistPoint3Parallelepiped3Result): void {
    expect(result.closest[0].equals(p)).toBe(true);

    const delta = sub(result.closest[0], result.closest[1]);
    expect(Math.sqrt(dot(delta, delta))).toBeCloseTo(result.distance, 9);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 9);

    // The closest point lies in the parallelepiped.
    const s = coordinates(result.closest[1], box);
    for (let i = 0; i < 3; ++i) {
        expect(Math.abs(s.values[i])).toBeLessThanOrEqual(1 + 1e-9);
    }
}

describe('DistPoint3Parallelepiped3', () => {
    const query = new DistPoint3Parallelepiped3();
    const unit = ppd([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);

    it('reports zero distance for interior points', () => {
        for (const p of [v(0, 0, 0), v(0.5, -0.5, 0.25), v(1, 1, 1),
            v(-1, 0.3, 0.9)]) {
            const result = query.compute(p, unit);
            expect(result.distance).toBeCloseTo(0, 12);
            expectConsistent(p, unit, result);
        }
    });

    it('computes known exterior distances for the unit cube', () => {
        // Face region.
        let result = query.compute(v(3, 0, 0), unit);
        expect(result.distance).toBeCloseTo(2, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(1, 12);

        // Edge region.
        result = query.compute(v(3, 4, 0), unit);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 9), 12);

        // Vertex region.
        result = query.compute(v(-3, -5, -7), unit);
        expect(result.distance).toBeCloseTo(Math.sqrt(4 + 16 + 36), 12);
        expect(result.closest[1].values[0]).toBeCloseTo(-1, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(-1, 12);
        expect(result.closest[1].values[2]).toBeCloseTo(-1, 12);
    });

    it('agrees with the point-canonical-box query for orthogonal axes', () => {
        const extent = v(2, 3, 0.5);
        const cbox = CanonicalBox.fromExtent(extent);
        const box = ppd([0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 0.5]);
        const cbQuery = new DistPointCanonicalBox();

        let seed = 55555;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 300; ++trial) {
            const p = v(rand() * 12 - 6, rand() * 12 - 6, rand() * 12 - 6);
            const expected = cbQuery.compute(p, cbox);
            const result = query.compute(p, box);
            expect(result.distance).toBeCloseTo(expected.distance, 9);
            for (let i = 0; i < 3; ++i) {
                expect(result.closest[1].values[i])
                    .toBeCloseTo(expected.closest[1].values[i], 8);
            }
            expectConsistent(p, box, result);
        }
    });

    it('matches a grid sampling for a sheared parallelepiped', () => {
        let seed = 20260901;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        const box = ppd([0.5, -0.25, 0.75], [1.5, 0.3, -0.2],
            [0.4, 1.2, 0.5], [-0.3, 0.2, 1.1]);
        for (let trial = 0; trial < 30; ++trial) {
            const p = v(rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4);
            const result = query.compute(p, box);
            const sampled = sampledDistance(p, box, 16);
            expect(result.distance).toBeLessThanOrEqual(sampled + 1e-9);
            expect(sampled - result.distance).toBeLessThan(0.25);
            expectConsistent(p, box, result);
        }
    });

    it('is invariant under rigid motion of the point and the solid', () => {
        const c = Math.cos(0.6), s = Math.sin(0.6);
        const rot = (p: Vector) => v(
            c * p.values[0] - s * p.values[2],
            p.values[1],
            s * p.values[0] + c * p.values[2]);
        const shift = v(1, -2, 3);

        const box = ppd([0, 0, 0], [1.2, 0.1, 0], [0.2, 1.1, 0.3],
            [0, 0.1, 0.9]);
        const moved = Parallelepiped3.fromCenterAxis(
            add(rot(box.center), shift),
            [rot(box.axis[0]), rot(box.axis[1]), rot(box.axis[2])]);

        let seed = 3141592;
        const rand = () => {
            seed = (1103515245 * seed + 12345) % 2147483648;
            return seed / 2147483648;
        };

        for (let trial = 0; trial < 100; ++trial) {
            const p = v(rand() * 8 - 4, rand() * 8 - 4, rand() * 8 - 4);
            const r0 = query.compute(p, box);
            const r1 = query.compute(add(rot(p), shift), moved);
            expect(r1.distance).toBeCloseTo(r0.distance, 9);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification wave (see VERIFYING.md): property-based cross-checks of the
// port against upstream DistPoint3Parallelepiped3.h.
// ---------------------------------------------------------------------------

describe('DistPoint3Parallelepiped3 verification', () => {
    const query = new DistPoint3Parallelepiped3();
    const obQuery = new DistPointOrientedBox();

    // A well-conditioned parallelepiped: the axes span a volume bounded away
    // from zero, so Inverse(B) in the query carries significant digits.
    const ppdArb = fc.tuple(wellScaledVector(3, -5, 5), wellScaledVector(3, -3, 3),
        wellScaledVector(3, -3, 3), wellScaledVector(3, -3, 3))
        // Parallelepiped3.fromCenterAxis requires a right-handed basis, so
        // the filter is on the signed volume rather than its magnitude.
        .filter(([, a0, a1, a2]) => dot(a0, cross(a1, a2)) > 0.5)
        .map(([c, a0, a1, a2]) =>
            Parallelepiped3.fromCenterAxis(c, [a0, a1, a2]));

    const pointArb = wellScaledVector(3, -8, 8);

    it('reports consistent distances and in-solid closest points', () => {
        check(fc.tuple(pointArb, ppdArb), ([p, box]) => {
            const r = query.compute(p, box);
            expect(r.closest[0].equals(p)).toBe(true);
            expectClose(r.sqrDistance, r.distance * r.distance, 1e-12, 1e-12);
            const d = sub(r.closest[0], r.closest[1]);
            expectClose(Math.sqrt(dot(d, d)), r.distance, 1e-9, 1e-9);
            const s = coordinates(r.closest[1], box);
            for (let i = 0; i < 3; ++i) {
                expect(Math.abs(s.values[i])).toBeLessThanOrEqual(1 + 1e-8);
            }
        });
    });

    it('reports zero distance and the point itself for interior points',
        () => {
            check(fc.tuple(ppdArb, fc.array(finite(-0.95, 0.95),
                { minLength: 3, maxLength: 3 })), ([box, s]) => {
                    let p = box.center.clone();
                    for (let i = 0; i < 3; ++i) {
                        p = add(p, mul(s[i], box.axis[i]));
                    }
                    const r = query.compute(p, box);
                    expectClose(r.distance, 0, 1e-9, 0);
                    expectVectorClose(r.closest[1], p, 1e-8, 1e-8);
                });
        });

    // With mutually orthogonal axes the parallelepiped is an oriented box
    // whose extents are the axis lengths, so the two queries must agree.
    // This exercises every face/edge/vertex region of GetMinimizer against an
    // independent implementation.
    it('agrees with the point-oriented-box query for orthogonal axes', () => {
        check(fc.tuple(pointArb, wellScaledVector(3, -5, 5), rotationFrame(3),
            fc.array(finite(0.1, 4), { minLength: 3, maxLength: 3 })),
            ([p, c, R, e]) => {
                const axes = [mul(e[0], R[0]), mul(e[1], R[1]),
                    mul(e[2], R[2])];
                const box = Parallelepiped3.fromCenterAxis(c, axes);
                const obox = OrientedBox.fromCenterAxisExtent(c, R,
                    v(e[0], e[1], e[2]));
                const r0 = query.compute(p, box);
                const r1 = obQuery.compute(p, obox);
                expectClose(r0.distance, r1.distance, 1e-7, 1e-7);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-6, 1e-6);
            });
    });

    it('never reports more than a grid sampling of the solid', () => {
        const rng = seededRandom(0x13572468);
        const box = ppd([0.5, -1, 2], [2, 0.5, 0], [-0.5, 1.5, 0.5],
            [0.25, 0, 1.75]);
        for (let k = 0; k < 40; ++k) {
            const p = v(14 * rng() - 7, 14 * rng() - 7, 14 * rng() - 7);
            const r = query.compute(p, box);
            expect(r.distance)
                .toBeLessThanOrEqual(sampledDistance(p, box, 16) + 1e-9);
        }
    }, 30000);

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(pointArb, ppdArb, rotationFrame(3),
            wellScaledVector(3, -5, 5)), ([p, box, R, t]) => {
                const xf = (q: Vector) => add(add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2])), t);
                const rot = (q: Vector) => add(add(
                    mul(q.values[0], R[0]), mul(q.values[1], R[1])),
                    mul(q.values[2], R[2]));
                const r0 = query.compute(p, box);
                const r1 = query.compute(xf(p), Parallelepiped3.fromCenterAxis(
                    xf(box.center), [rot(box.axis[0]), rot(box.axis[1]),
                        rot(box.axis[2])]));
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectVectorClose(xf(r0.closest[1]), r1.closest[1],
                    1e-6, 1e-6);
            });
    });

    it('is invariant under negating two axes', () => {
        // The solid is unchanged when two axes are negated (a rotation by pi
        // about the third), and the basis stays right-handed. This catches a
        // sign slip in the fixed coordinate of the face helpers, where the
        // -1 and +1 faces would swap roles.
        check(fc.tuple(pointArb, ppdArb, fc.integer({ min: 0, max: 2 })),
            ([p, box, k]) => {
                const axes = [box.axis[0].clone(), box.axis[1].clone(),
                    box.axis[2].clone()];
                axes[(k + 1) % 3] = mul(-1, axes[(k + 1) % 3]);
                axes[(k + 2) % 3] = mul(-1, axes[(k + 2) % 3]);
                const r0 = query.compute(p, box);
                const r1 = query.compute(p,
                    Parallelepiped3.fromCenterAxis(box.center, axes));
                expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
                expectVectorClose(r0.closest[1], r1.closest[1], 1e-6, 1e-6);
            });
    });

    it('is invariant under a cyclic permutation of the axes', () => {
        check(fc.tuple(pointArb, ppdArb), ([p, box]) => {
            const r0 = query.compute(p, box);
            const r1 = query.compute(p, Parallelepiped3.fromCenterAxis(
                box.center, [box.axis[1], box.axis[2], box.axis[0]]));
            const r2 = query.compute(p, Parallelepiped3.fromCenterAxis(
                box.center, [box.axis[2], box.axis[0], box.axis[1]]));
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expectClose(r0.distance, r2.distance, 1e-8, 1e-8);
        });
    });
});
