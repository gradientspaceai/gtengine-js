import { describe, expect, it } from 'vitest';
import { DistPointHyperplane } from '../src/DistPointHyperplane';
import { Hyperplane } from '../src/Hyperplane';
import { Vector, dot, normalize } from '../src/Vector';

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
