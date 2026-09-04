import { describe, expect, it } from 'vitest';
import { DistPointHyperplane } from '../src/DistPointHyperplane.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, add, dot, length, mul, normalize, sub }
    from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame, unitVector,
    wellScaled, wellScaledVector
} from './helpers/arbitraries.js';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function unit(...values: number[]): Vector {
    const u = v(...values);
    normalize(u);
    return u;
}

describe('DistPointHyperplane', () => {
    const query = new DistPointHyperplane();

    it('computes the signed distance to a 2D line', () => {
        // The line x = 3 with normal (1,0).
        const plane = Hyperplane.fromNormalConstant(v(1, 0), 3);
        const result = query.compute(v(7, -2), plane);
        expect(result.signedDistance).toBeCloseTo(4, 12);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[0].values).toEqual([7, -2]);
        expect(result.closest[1].values[0]).toBeCloseTo(3, 12);
        expect(result.closest[1].values[1]).toBeCloseTo(-2, 12);
    });

    it('reports a negative signed distance on the other side', () => {
        const plane = Hyperplane.fromNormalConstant(v(1, 0), 3);
        const result = query.compute(v(-1, 5), plane);
        expect(result.signedDistance).toBeCloseTo(-4, 12);
        expect(result.distance).toBeCloseTo(4, 12);
        expect(result.closest[1].values[0]).toBeCloseTo(3, 12);
    });

    it('computes the distance to a 3D plane with an oblique normal', () => {
        // The plane x + y + z = 3, whose unit normal is (1,1,1)/sqrt(3) and
        // whose constant is 3/sqrt(3) = sqrt(3).
        const n = unit(1, 1, 1);
        const plane = Hyperplane.fromNormalConstant(n, Math.sqrt(3));
        const result = query.compute(v(0, 0, 0), plane);
        expect(result.signedDistance).toBeCloseTo(-Math.sqrt(3), 12);
        expect(result.distance).toBeCloseTo(Math.sqrt(3), 12);
        // The closest point is (1,1,1).
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(1, 12);
        }
    });

    it('returns zero distance for a point on the hyperplane', () => {
        const n = unit(2, -1, 3);
        const p = v(1, 2, 3);
        const plane = Hyperplane.fromNormalOrigin(n, p);
        const result = query.compute(p, plane);
        expect(result.distance).toBeCloseTo(0, 12);
        expect(result.signedDistance).toBeCloseTo(0, 12);
        for (let i = 0; i < 3; ++i) {
            expect(result.closest[1].values[i]).toBeCloseTo(p.values[i], 12);
        }
    });

    it('places the closest point on the hyperplane and along the normal',
        () => {
            let seed = 918273645;
            const rand = (): number => {
                seed = (1103515245 * seed + 12345) % 2147483648;
                return seed / 2147483648;
            };
            for (let trial = 0; trial < 200; ++trial) {
                const n = 2 + (trial % 4);
                const normalValues: number[] = [];
                for (let i = 0; i < n; ++i) {
                    normalValues.push(2 * rand() - 1);
                }
                const normalVec = v(...normalValues);
                if (dot(normalVec, normalVec) < 1e-6) {
                    continue;
                }
                normalize(normalVec);
                const constant = 4 * rand() - 2;
                const plane = Hyperplane.fromNormalConstant(normalVec,
                    constant);

                const pointValues: number[] = [];
                for (let i = 0; i < n; ++i) {
                    pointValues.push(6 * rand() - 3);
                }
                const point = v(...pointValues);
                const result = query.compute(point, plane);

                // The closest point lies on the hyperplane.
                expect(dot(plane.normal, result.closest[1]))
                    .toBeCloseTo(constant, 10);
                // The displacement is along the normal with the reported
                // signed length.
                for (let i = 0; i < n; ++i) {
                    expect(point.values[i] - result.closest[1].values[i])
                        .toBeCloseTo(
                            result.signedDistance * plane.normal.values[i],
                            10);
                }
                expect(result.distance)
                    .toBeCloseTo(Math.abs(result.signedDistance), 12);
                expect(result.closest[0].values).toEqual(point.values);
            }
        });

    it('throws for a dimension mismatch', () => {
        const plane = Hyperplane.fromNormalConstant(v(1, 0), 1);
        expect(() => query.compute(v(1, 2, 3), plane)).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Independent verification (V21): property-based tests against the upstream
// header DistPointHyperplane.h.
// ---------------------------------------------------------------------------

// The query documents a unit-length normal, so generate one.
const planeArb = (n: number): fc.Arbitrary<Hyperplane> =>
    fc.tuple(unitVector(n), wellScaled(-8, 8))
        .map(([nrm, c]) => Hyperplane.fromNormalConstant(nrm, c));

describe('DistPointHyperplane verification', () => {
    const query = new DistPointHyperplane();

    for (const n of [2, 3, 4, 5]) {
        it(`result is self consistent in ${n}D`, () => {
            check(fc.tuple(wellScaledVector(n, -10, 10), planeArb(n)),
                ([p, plane]) => {
                    const r = query.compute(p, plane);
                    expectClose(r.distance, Math.abs(r.signedDistance), 0, 0);
                    expectVectorClose(r.closest[0], p, 0, 0);
                    // closest[1] is on the hyperplane.
                    expectClose(dot(plane.normal, r.closest[1]),
                        plane.constant, 1e-9, 1e-9);
                    // The reported distance is realized by the pair.
                    expectClose(r.distance,
                        length(sub(r.closest[0], r.closest[1])), 1e-9, 1e-9);
                    // The residual is along the normal.
                    const residual = sub(r.closest[0], r.closest[1]);
                    expectVectorClose(residual,
                        mul(r.signedDistance, plane.normal), 1e-9, 1e-9);
                });
        });
    }

    it('the signed distance has the sign of the side the point is on', () => {
        check(fc.tuple(planeArb(3), wellScaled(-6, 6),
            wellScaledVector(3, -5, 5)), ([plane, h, tangential]) => {
            // Build a point at signed height h above the plane: start from
            // the plane origin, move within the plane, then along the normal.
            const inPlane = sub(tangential,
                mul(dot(tangential, plane.normal), plane.normal));
            const p = add(add(plane.origin, inPlane), mul(h, plane.normal));
            const r = query.compute(p, plane);
            expectClose(r.signedDistance, h, 1e-8, 1e-8);
            expectClose(r.distance, Math.abs(h), 1e-8, 1e-8);
            expectVectorClose(r.closest[1], add(plane.origin, inPlane), 1e-8,
                1e-8);
        });
    });

    it('is not larger than the distance to any sampled plane point', () => {
        check(fc.tuple(wellScaledVector(3, -10, 10), planeArb(3),
            wellScaledVector(3, -8, 8)), ([p, plane, u]) => {
            const q = add(plane.origin,
                sub(u, mul(dot(u, plane.normal), plane.normal)));
            expect(query.compute(p, plane).distance)
                .toBeLessThanOrEqual(length(sub(p, q)) + 1e-9);
        });
    });

    it('flipping the plane orientation negates the signed distance', () => {
        check(fc.tuple(wellScaledVector(3, -10, 10), planeArb(3)),
            ([p, plane]) => {
                const flipped = Hyperplane.fromNormalConstant(
                    mul(-1, plane.normal), -plane.constant);
                const r0 = query.compute(p, plane);
                const r1 = query.compute(p, flipped);
                expectClose(r1.signedDistance, -r0.signedDistance, 0, 0);
                expectClose(r1.distance, r0.distance, 0, 0);
                expectVectorClose(r1.closest[1], r0.closest[1], 1e-12, 1e-12);
            });
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(wellScaledVector(3, -10, 10), planeArb(3),
            rotationFrame(3), wellScaledVector(3, -5, 5)),
        ([p, plane, R, tr]) => {
            const rot = (x: Vector): Vector => {
                let y = new Vector(3);
                for (let i = 0; i < 3; ++i) {
                    y = add(y, mul(x.values[i], R[i]));
                }
                return y;
            };
            const moved = Hyperplane.fromNormalOrigin(rot(plane.normal),
                add(rot(plane.origin), tr));
            const r0 = query.compute(p, plane);
            const r1 = query.compute(add(rot(p), tr), moved);
            expectClose(r0.distance, r1.distance, 1e-8, 1e-8);
            expectVectorClose(add(rot(r0.closest[1]), tr), r1.closest[1],
                1e-8, 1e-8);
        });
    });

    it('returns zero for points on the hyperplane', () => {
        check(fc.tuple(planeArb(4), wellScaledVector(4, -6, 6)),
            ([plane, u]) => {
                const q = add(plane.origin,
                    sub(u, mul(dot(u, plane.normal), plane.normal)));
                expect(query.compute(q, plane).distance)
                    .toBeLessThanOrEqual(1e-9);
            });
    });

    it('throws for dimension mismatches and for dimension 1', () => {
        // Upstream has a static_assert(N >= 2); the port makes it a runtime
        // check. Hyperplane itself rejects dimension 1, so build the 1D case
        // by hand.
        const plane3 = Hyperplane.fromNormalConstant(v(0, 0, 1), 1);
        expect(() => query.compute(v(1, 2), plane3)).toThrow();
        const plane1 = Hyperplane.fromNormalConstant(v(0, 1), 1);
        plane1.normal = v(1);
        expect(() => query.compute(v(1), plane1)).toThrow('Invalid dimension.');
    });

    it('does not mutate its inputs', () => {
        check(fc.tuple(wellScaledVector(3, -10, 10), planeArb(3)),
            ([p, plane]) => {
                const p0 = p.clone();
                const nrm = plane.normal.clone();
                const c = plane.constant;
                const r = query.compute(p, plane);
                expect(p.values).toEqual(p0.values);
                expect(plane.normal.values).toEqual(nrm.values);
                expect(plane.constant).toBe(c);
                r.closest[0].values[0] = 321;
                expect(p.values).toEqual(p0.values);
            });
    });
});
