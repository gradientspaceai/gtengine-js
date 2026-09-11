import { describe, it, expect } from 'vitest';
import { Hyperplane } from '../src/Hyperplane.js';
import { Vector, dot, sub, add, mul, normalize, length } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, compareKeys, expectClose, expectStrictWeakOrder,
    expectVectorClose, fc, invertibleMatrix, unitVector, wellScaled,
    wellScaledVector
} from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let a = seed >>> 0;
    return () => {
        a = (a + 0x6D2B79F5) >>> 0;
        let t = a;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

function randomVector(n: number, rand: () => number): Vector {
    const v = new Vector(n);
    for (let d = 0; d < n; ++d) {
        v.values[d] = 4 * rand() - 2;
    }
    return v;
}

describe('Hyperplane', () => {
    it('default-constructs the plane with normal (0,...,0,1)', () => {
        for (const n of [2, 3, 4]) {
            const plane = new Hyperplane(n);
            expect(plane.dimension).toBe(n);
            const expected = new Array<number>(n).fill(0);
            expected[n - 1] = 1;
            expect(plane.normal.values).toEqual(expected);
            expect(plane.origin.values).toEqual(new Array<number>(n).fill(0));
            expect(plane.constant).toBe(0);
        }
    });

    it('rejects dimensions less than 2', () => {
        expect(() => new Hyperplane(1)).toThrow();
    });

    it('constructs from a normal and a constant', () => {
        const normal = Vector.fromArray([0, 0, 1]);
        const plane = Hyperplane.fromNormalConstant(normal, 3);
        expect(plane.normal.values).toEqual([0, 0, 1]);
        expect(plane.origin.values).toEqual([0, 0, 3]);
        expect(plane.constant).toBe(3);

        // The vector argument is copied.
        normal.values[2] = 100;
        expect(plane.normal.values[2]).toBe(1);

        // The origin is the point on the plane closest to the origin.
        expect(dot(plane.normal, plane.origin)).toBeCloseTo(plane.constant, 12);
    });

    it('constructs from a normal and an origin', () => {
        const normal = Vector.fromArray([0, 1, 0]);
        const origin = Vector.fromArray([5, 2, -7]);
        const plane = Hyperplane.fromNormalOrigin(normal, origin);
        expect(plane.normal.values).toEqual([0, 1, 0]);
        expect(plane.origin.values).toEqual([5, 2, -7]);
        expect(plane.constant).toBe(2);

        normal.values[1] = 100;
        origin.values[0] = 100;
        expect(plane.normal.values[1]).toBe(1);
        expect(plane.origin.values[0]).toBe(5);
    });

    it('rejects mismatched sizes in fromNormalOrigin', () => {
        expect(() => Hyperplane.fromNormalOrigin(
            Vector.fromArray([1, 0]), Vector.fromArray([1, 0, 0]))).toThrow();
    });

    it('constructs a 3D plane from three points', () => {
        const p = [
            Vector.fromArray([1, 0, 0]),
            Vector.fromArray([0, 1, 0]),
            Vector.fromArray([0, 0, 1])
        ];
        const plane = Hyperplane.fromPoints(p);
        const invSqrt3 = 1 / Math.sqrt(3);
        for (let d = 0; d < 3; ++d) {
            expect(plane.normal.values[d]).toBeCloseTo(invSqrt3, 12);
        }
        expect(plane.constant).toBeCloseTo(invSqrt3, 12);

        // Each input point lies on the plane.
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 12);
        }

        // The normal matches UnitCross(p1-p0, p2-p0) including its sign.
        const c = cross(sub(p[1], p[0]), sub(p[2], p[0]));
        normalize(c);
        expect(plane.normal.values).toEqual(c.values);
    });

    it('constructs a 2D line from two points', () => {
        const p = [Vector.fromArray([1, 1]), Vector.fromArray([3, 1])];
        const plane = Hyperplane.fromPoints(p);
        expect(length(plane.normal)).toBeCloseTo(1, 12);
        // The edge is (2,0), so the normal is +/-(0,1).
        expect(Math.abs(plane.normal.values[0])).toBeCloseTo(0, 12);
        expect(Math.abs(plane.normal.values[1])).toBeCloseTo(1, 12);
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 12);
        }
        expect(dot(plane.normal, plane.origin)).toBeCloseTo(plane.constant, 12);
    });

    it('constructs a 4D hyperplane from four points', () => {
        const p = [
            Vector.fromArray([1, 0, 0, 0]),
            Vector.fromArray([0, 1, 0, 0]),
            Vector.fromArray([0, 0, 1, 0]),
            Vector.fromArray([0, 0, 0, 1])
        ];
        const plane = Hyperplane.fromPoints(p);
        expect(length(plane.normal)).toBeCloseTo(1, 10);
        const half = 0.5;
        for (let d = 0; d < 4; ++d) {
            expect(Math.abs(plane.normal.values[d])).toBeCloseTo(half, 10);
        }
        for (const q of p) {
            expect(dot(plane.normal, q) - plane.constant).toBeCloseTo(0, 10);
        }
    });

    it('produces a unit normal orthogonal to all edges (randomized)', () => {
        const rand = makeRandom(0x91A2E);
        let maxError = 0;
        for (let trial = 0; trial < 150; ++trial) {
            const n = 2 + (trial % 3);  // 2, 3, 4
            const p: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                p.push(randomVector(n, rand));
            }
            const plane = Hyperplane.fromPoints(p);

            maxError = Math.max(maxError,
                Math.abs(length(plane.normal) - 1));
            for (let i = 1; i < n; ++i) {
                maxError = Math.max(maxError,
                    Math.abs(dot(plane.normal, sub(p[i], p[0]))));
            }
            // Every input point is on the plane.
            for (let i = 0; i < n; ++i) {
                maxError = Math.max(maxError,
                    Math.abs(dot(plane.normal, p[i]) - plane.constant));
            }
            // origin = constant * normal and Dot(normal, origin) = constant.
            maxError = Math.max(maxError, Math.abs(
                dot(plane.normal, plane.origin) - plane.constant));
        }
        expect(maxError).toBeLessThan(1e-9);
    });

    it('gives signed point-plane distances with the expected sign', () => {
        // The plane z = 2 with normal (0,0,1).
        const plane = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 2);
        const above = Vector.fromArray([7, -3, 5]);
        const below = Vector.fromArray([7, -3, -1]);
        const on = Vector.fromArray([7, -3, 2]);
        expect(dot(plane.normal, above) - plane.constant).toBeCloseTo(3, 12);
        expect(dot(plane.normal, below) - plane.constant).toBeCloseTo(-3, 12);
        expect(dot(plane.normal, on) - plane.constant).toBeCloseTo(0, 12);

        // The signed distance of a point measured from the plane origin is
        // the same quantity.
        expect(dot(plane.normal, sub(above, plane.origin)))
            .toBeCloseTo(3, 12);
    });

    it('reproduces the signed distance for random planes and points', () => {
        const rand = makeRandom(0x7E571);
        let maxError = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const normal = randomVector(3, rand);
            if (normalize(normal) === 0) {
                continue;
            }
            const origin = randomVector(3, rand);
            const plane = Hyperplane.fromNormalOrigin(normal, origin);
            const t = 4 * rand() - 2;
            // A point at signed distance t from the plane.
            const inPlane = add(origin, mul(
                cross(normal, Vector.fromArray([1, 2, 3])), rand()));
            const X = add(inPlane, mul(plane.normal, t));
            maxError = Math.max(maxError,
                Math.abs((dot(plane.normal, X) - plane.constant) - t));
        }
        expect(maxError).toBeLessThan(1e-10);
    });

    it('yields a zero normal for degenerate point sets', () => {
        // Collinear points in 3D.
        const collinear = Hyperplane.fromPoints([
            Vector.fromArray([0, 0, 0]),
            Vector.fromArray([1, 1, 1]),
            Vector.fromArray([2, 2, 2])
        ]);
        expect(collinear.normal.values).toEqual([0, 0, 0]);
        expect(collinear.constant).toBe(0);
        expect(collinear.origin.values).toEqual([0, 0, 0]);

        // Coincident points in 3D.
        const coincident = Hyperplane.fromPoints([
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([1, 2, 3]),
            Vector.fromArray([1, 2, 3])
        ]);
        expect(coincident.normal.values).toEqual([0, 0, 0]);

        // Coincident points in 2D.
        const coincident2 = Hyperplane.fromPoints([
            Vector.fromArray([4, 5]),
            Vector.fromArray([4, 5])
        ]);
        expect(coincident2.normal.values).toEqual([0, 0]);
        expect(coincident2.constant).toBe(0);
    });

    it('rejects malformed point arrays', () => {
        expect(() => Hyperplane.fromPoints(
            [Vector.fromArray([1, 2, 3])])).toThrow();
        expect(() => Hyperplane.fromPoints([
            Vector.fromArray([1, 2]),
            Vector.fromArray([1, 2, 3])])).toThrow();
    });

    it('clones and compares for sorted containers', () => {
        const a = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 1);
        const b = Hyperplane.fromNormalConstant(
            Vector.fromArray([0, 0, 1]), 2);
        const c = Hyperplane.fromNormalConstant(
            Vector.fromArray([-1, 0, 0]), 1);

        const copy = a.clone();
        expect(copy.equals(a)).toBe(true);
        expect(copy.notEquals(a)).toBe(false);
        copy.constant = 99;
        expect(copy.equals(a)).toBe(false);
        expect(a.constant).toBe(1);

        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(c.lessThan(a)).toBe(true);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(a.greaterThan(c)).toBe(true);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);
        expect(b.greaterThanOrEqual(a)).toBe(true);
    });
});

// ---------------------------------------------------------------------------
// V43 verification: property-based checks of the constructors, of the
// point-set construction in every dimension, and of the signed distance.
// ---------------------------------------------------------------------------

describe('Hyperplane verification', () => {
    // N points in general position: p[0] is arbitrary and the N-1 edges are
    // the leading columns of a well-conditioned invertible matrix, so the
    // simplex is never degenerate.
    const pointsArb = (n: number) =>
        fc.tuple(wellScaledVector(n, -5, 5), invertibleMatrix(n))
            .map(([p0, M]) => {
                const p = [p0];
                for (let i = 1; i < n; ++i) {
                    p.push(add(p0, M.getCol(i - 1)));
                }
                return p;
            });

    it('fromPoints yields a unit normal containing all the points', () => {
        // N = 2 uses ComputeOrthogonalComplement, N = 3 uses UnitCross and
        // N >= 4 uses the SVD, so all three code paths are exercised.
        for (const n of [2, 3, 4, 5]) {
            check(pointsArb(n), p => {
                const plane = Hyperplane.fromPoints(p);
                expect(plane.dimension).toBe(n);
                expectClose(length(plane.normal), 1, 1e-12, 1e-12);
                for (let i = 0; i < n; ++i) {
                    // Dot(U, p[i]) = c for every input point. The dot product
                    // has terms of size |p|, so the tolerance scales with it.
                    const scale = 1 + length(p[i]);
                    expectClose(dot(plane.normal, p[i]), plane.constant,
                        1e-11 * scale, 0);
                }
                // origin = c * U is the point of the plane closest to 0.
                expectVectorClose(plane.origin,
                    mul(plane.normal, plane.constant), 1e-12, 1e-12);
                expectClose(dot(plane.normal, plane.origin), plane.constant,
                    1e-12, 1e-12);
            }, 100);
        }
    });

    it('fromPoints normal is orthogonal to every edge', () => {
        for (const n of [2, 3, 4, 5]) {
            check(pointsArb(n), p => {
                const plane = Hyperplane.fromPoints(p);
                for (let i = 1; i < n; ++i) {
                    const edge = sub(p[i], p[0]);
                    expectClose(dot(plane.normal, edge), 0,
                        1e-11 * (1 + length(edge)), 0);
                }
            }, 100);
        }
    });

    it('fromPoints uses UnitCross orientation in 3D', () => {
        check(pointsArb(3), p => {
            const plane = Hyperplane.fromPoints(p);
            const expected = cross(sub(p[1], p[0]), sub(p[2], p[0]));
            normalize(expected);
            expectVectorClose(plane.normal, expected, 1e-12, 1e-12);
        }, 100);
    });

    it('fromPoints uses the left perpendicular orientation in 2D', () => {
        // Upstream cannot build Hyperplane<2> from points at all (the
        // SingularValueDecomposition(2, 1, 32) constructor asserts), so the
        // port picks ComputeOrthogonalComplement's -Perp(edge), which is the
        // edge rotated by +90 degrees. Pin the convention.
        check(pointsArb(2), p => {
            const edge = sub(p[1], p[0]);
            const expected = Vector.fromArray([-edge.get(1), edge.get(0)]);
            normalize(expected);
            const plane = Hyperplane.fromPoints(p);
            expectVectorClose(plane.normal, expected, 1e-12, 1e-12);
            // DotPerp(edge, normal) > 0: the normal is to the left of the
            // directed edge p[0] -> p[1].
            expect(edge.get(0) * plane.normal.get(1)
                - edge.get(1) * plane.normal.get(0)).toBeGreaterThan(0);
        }, 100);
    });

    it('fromNormalConstant and fromNormalOrigin agree on the plane', () => {
        for (const n of [2, 3, 4]) {
            check(fc.tuple(unitVector(n), wellScaled(-5, 5)), ([u, c]) => {
                const byConstant = Hyperplane.fromNormalConstant(u, c);
                expectClose(byConstant.constant, c, 0, 0);
                expectVectorClose(byConstant.origin, mul(u, c), 0, 0);
                // Rebuilding from (normal, origin) reproduces the constant.
                const byOrigin = Hyperplane.fromNormalOrigin(u,
                    byConstant.origin);
                expectClose(byOrigin.constant, c, 1e-12, 1e-12);
                expectVectorClose(byOrigin.normal, byConstant.normal, 0, 0);
            });
        }
    });

    it('reports the signed distance and its foot point', () => {
        for (const n of [2, 3, 4]) {
            check(fc.tuple(unitVector(n), wellScaled(-5, 5),
                wellScaledVector(n, -5, 5)), ([u, c, X]) => {
                const plane = Hyperplane.fromNormalConstant(u, c);
                const signed = dot(plane.normal, X) - plane.constant;
                // The foot of the perpendicular is on the plane.
                const foot = sub(X, mul(plane.normal, signed));
                expectClose(dot(plane.normal, foot), plane.constant,
                    1e-11 * (1 + length(X)), 0);
                // Moving along +normal increases the signed distance by the
                // step length.
                const moved = add(X, plane.normal);
                expectClose(dot(plane.normal, moved) - plane.constant,
                    signed + 1, 1e-11 * (1 + length(X)), 0);
            });
        }
    });

    it('copies its vector arguments (C++ value semantics)', () => {
        check(fc.tuple(unitVector(3), wellScaledVector(3, -5, 5)),
            ([u, P]) => {
                const plane = Hyperplane.fromNormalOrigin(u, P);
                const c = plane.constant;
                u.set(0, u.get(0) + 1);
                P.set(0, P.get(0) + 1);
                expect(plane.normal.get(0)).not.toBe(u.get(0));
                expect(plane.origin.get(0)).not.toBe(P.get(0));
                expect(plane.constant).toBe(c);
                const copy = plane.clone();
                copy.normal.set(1, 17);
                expect(plane.normal.get(1)).not.toBe(17);
            });
    });

    it('has a trichotomous comparison', () => {
        const planeArb = fc.tuple(unitVector(3), wellScaled(-3, 3))
            .map(([u, c]) => Hyperplane.fromNormalConstant(u, c));
        check(fc.tuple(planeArb, planeArb), ([a, b]) => {
            const lt = a.lessThan(b), gt = a.greaterThan(b), eq = a.equals(b);
            expect([lt, gt, eq].filter(x => x).length).toBe(1);
            expect(a.notEquals(b)).toBe(!eq);
            expect(a.lessThanOrEqual(b)).toBe(!gt);
            expect(a.greaterThanOrEqual(b)).toBe(!lt);
        });
    });

    it('orders by the upstream member sequence, a strict weak ordering', () => {
        const key = (h: Hyperplane): number[] =>
            [...h.normal.values, ...h.origin.values, h.constant];
        const small = fc.tuple(fc.integer({ min: -1, max: 1 }),
            fc.integer({ min: -1, max: 1 }))
            .map(([d, c]) => Hyperplane.fromNormalConstant(
                Vector.fromArray([d, 0, 1 - Math.abs(d)]), c));
        check(fc.tuple(small, small), ([a, b]) => {
            expect(a.lessThan(b)).toBe(compareKeys(key(a), key(b)) < 0);
            expect(a.equals(b)).toBe(compareKeys(key(a), key(b)) === 0);
        });
        check(fc.array(small, { minLength: 3, maxLength: 5 }), items => {
            expectStrictWeakOrder(items, (x, y) => x.lessThan(y));
        }, 50);
    });
});
