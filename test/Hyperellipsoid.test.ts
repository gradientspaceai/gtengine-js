import { describe, it, expect } from 'vitest';
import {
    Hyperellipsoid, hyperellipsoidNumCoefficients
} from '../src/Hyperellipsoid';
import { Matrix, multiplyAB, mulMatrix } from '../src/Matrix';
import { Vector, dot, sub, add, mul, normalize } from '../src/Vector';

// A small deterministic pseudorandom generator (mulberry32) so the randomized
// cross-checks are reproducible.
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

// An orthonormal basis of R^2 obtained by rotating by 'angle'.
function basis2(angle: number): Vector[] {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    return [Vector.fromArray([c, s]), Vector.fromArray([-s, c])];
}

// An orthonormal basis of R^3 built by Gram-Schmidt on three random vectors.
function basis3(rand: () => number): Vector[] {
    const raw: Vector[] = [];
    for (let d = 0; d < 3; ++d) {
        raw.push(Vector.fromArray([
            2 * rand() - 1, 2 * rand() - 1, 2 * rand() - 1]));
    }
    const u: Vector[] = [];
    for (let d = 0; d < 3; ++d) {
        let v = raw[d].clone();
        for (let k = 0; k < d; ++k) {
            v = sub(v, mul(u[k], dot(u[k], v)));
        }
        normalize(v);
        u.push(v);
    }
    return u;
}

function randomHyperellipsoid(n: number, rand: () => number): Hyperellipsoid {
    const axis = (n === 2 ? basis2(2 * Math.PI * rand()) : basis3(rand));
    const center = new Vector(n);
    const extent = new Vector(n);
    for (let d = 0; d < n; ++d) {
        center.values[d] = 4 * rand() - 2;
        // Well-separated extents so the eigenvalue ordering is unambiguous.
        extent.values[d] = 0.5 + 1.5 * d + rand() * 0.5;
    }
    return Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
}

function maxAbsDifference(M0: Matrix, M1: Matrix): number {
    let maxAbs = 0;
    for (let i = 0; i < M0.numElements; ++i) {
        maxAbs = Math.max(maxAbs, Math.abs(M0.values[i] - M1.values[i]));
    }
    return maxAbs;
}

describe('Hyperellipsoid', () => {
    it('default-constructs the unit hyperellipsoid', () => {
        for (const n of [2, 3, 4]) {
            const he = new Hyperellipsoid(n);
            expect(he.dimension).toBe(n);
            expect(he.center.values).toEqual(new Array<number>(n).fill(0));
            expect(he.extent.values).toEqual(new Array<number>(n).fill(1));
            for (let d = 0; d < n; ++d) {
                const unit = new Vector(n);
                unit.makeUnit(d);
                expect(he.axis[d].equals(unit)).toBe(true);
            }
        }
    });

    it('rejects dimensions less than 2', () => {
        expect(() => new Hyperellipsoid(1)).toThrow();
    });

    it('copies the vectors in the factory and in clone', () => {
        const center = Vector.fromArray([1, 2]);
        const axis = [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])];
        const extent = Vector.fromArray([3, 4]);
        const he = Hyperellipsoid.fromCenterAxisExtent(center, axis, extent);
        center.values[0] = 100;
        axis[0].values[0] = 100;
        extent.values[0] = 100;
        expect(he.center.values[0]).toBe(1);
        expect(he.axis[0].values[0]).toBe(1);
        expect(he.extent.values[0]).toBe(3);

        const copy = he.clone();
        copy.center.values[0] = 7;
        expect(he.center.values[0]).toBe(1);
        expect(copy.equals(he)).toBe(false);
        expect(copy.notEquals(he)).toBe(true);
        expect(he.clone().equals(he)).toBe(true);
    });

    it('computes M for an axis-aligned ellipse with known values', () => {
        const he = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([2, 3]));
        const M = he.getM();
        expect(M.get(0, 0)).toBeCloseTo(1 / 4, 12);
        expect(M.get(1, 1)).toBeCloseTo(1 / 9, 12);
        expect(M.get(0, 1)).toBeCloseTo(0, 12);
        expect(M.get(1, 0)).toBeCloseTo(0, 12);

        const MInv = he.getMInverse();
        expect(MInv.get(0, 0)).toBeCloseTo(4, 12);
        expect(MInv.get(1, 1)).toBeCloseTo(9, 12);
    });

    it('has M * MInverse = I for random ellipses and ellipsoids', () => {
        const rand = makeRandom(0x51DE);
        let maxError = 0;
        for (let trial = 0; trial < 200; ++trial) {
            const n = (trial % 2 === 0 ? 2 : 3);
            const he = randomHyperellipsoid(n, rand);
            const M = he.getM();
            const MInv = he.getMInverse();
            maxError = Math.max(maxError,
                maxAbsDifference(multiplyAB(M, MInv), Matrix.identity(n, n)));
            maxError = Math.max(maxError,
                maxAbsDifference(multiplyAB(MInv, M), Matrix.identity(n, n)));
            // M is symmetric.
            for (let r = 0; r < n; ++r) {
                for (let c = 0; c < n; ++c) {
                    maxError = Math.max(maxError,
                        Math.abs(M.get(r, c) - M.get(c, r)));
                }
            }
        }
        expect(maxError).toBeLessThan(1e-10);
    });

    it('satisfies (X-K)^T M (X-K) = 1 for points on the hyperellipsoid', () => {
        const rand = makeRandom(0xA11CE);
        let maxError = 0;
        for (let trial = 0; trial < 100; ++trial) {
            const n = (trial % 2 === 0 ? 2 : 3);
            const he = randomHyperellipsoid(n, rand);
            const M = he.getM();
            for (let k = 0; k < 5; ++k) {
                // Build y with sum (y[d]/e[d])^2 = 1 by normalizing a random
                // direction, then X = K + sum y[d]*U[d].
                const t = new Vector(n);
                for (let d = 0; d < n; ++d) {
                    t.values[d] = 2 * rand() - 1;
                }
                if (normalize(t) === 0) {
                    continue;
                }
                let X = he.center.clone();
                for (let d = 0; d < n; ++d) {
                    X = add(X, mul(he.axis[d],
                        t.values[d] * he.extent.values[d]));
                }
                const diff = sub(X, he.center);
                const q = dot(diff, mulMatrix(M, diff));
                maxError = Math.max(maxError, Math.abs(q - 1));
            }
        }
        expect(maxError).toBeLessThan(1e-10);
    });

    it('evaluates the quadratic equation to zero on the hyperellipsoid', () => {
        // The ellipse (x/2)^2 + (y/3)^2 = 1 centered at (1,-1).
        const he = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([1, -1]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([2, 3]));
        const { A, B, C } = he.toCoefficientsABC();
        for (const angle of [0, 0.7, 1.9, 3.4, 5.1]) {
            const X = Vector.fromArray([
                1 + 2 * Math.cos(angle), -1 + 3 * Math.sin(angle)]);
            const value = C + dot(B, X) + dot(X, mulMatrix(A, X));
            expect(value).toBeCloseTo(0, 10);
        }
    });

    it('packs and unpacks the coefficient array consistently', () => {
        const rand = makeRandom(0xBEE5);
        let maxError = 0;
        for (let trial = 0; trial < 100; ++trial) {
            const n = (trial % 2 === 0 ? 2 : 3);
            const he = randomHyperellipsoid(n, rand);
            const coeff = he.toCoefficients();
            expect(coeff.length).toBe(hyperellipsoidNumCoefficients(n));

            // The largest-magnitude quadratic coefficient is normalized to 1.
            let maxQuadratic = 0;
            for (let i = n + 1; i < coeff.length; ++i) {
                maxQuadratic = Math.max(maxQuadratic, Math.abs(coeff[i]));
            }
            maxError = Math.max(maxError, Math.abs(maxQuadratic - 1));

            // The equation is satisfied at points of the hyperellipsoid; the
            // normalization is an overall scale, so 0 = C + B.X + X.A.X still
            // holds.
            const evaluate = (X: Vector): number => {
                let value = coeff[0];
                for (let d = 0; d < n; ++d) {
                    value += coeff[1 + d] * X.values[d];
                }
                let i = n + 1;
                for (let r = 0; r < n; ++r) {
                    value += coeff[i++] * X.values[r] * X.values[r];
                    for (let c = r + 1; c < n; ++c) {
                        value += coeff[i++] * X.values[r] * X.values[c];
                    }
                }
                return value;
            };
            for (let d = 0; d < n; ++d) {
                for (const sign of [1, -1]) {
                    const X = add(he.center,
                        mul(he.axis[d], sign * he.extent.values[d]));
                    maxError = Math.max(maxError, Math.abs(evaluate(X)));
                }
            }
        }
        expect(maxError).toBeLessThan(1e-9);
    });

    it('round-trips ToCoefficients and FromCoefficients', () => {
        const rand = makeRandom(0xC0FFEE);
        let maxError = 0;
        for (let trial = 0; trial < 120; ++trial) {
            const n = (trial % 2 === 0 ? 2 : 3);
            const he = randomHyperellipsoid(n, rand);
            const coeff = he.toCoefficients();

            const result = new Hyperellipsoid(n);
            expect(result.fromCoefficients(coeff)).toBe(true);

            // The center is recovered exactly (up to round-off).
            for (let d = 0; d < n; ++d) {
                maxError = Math.max(maxError,
                    Math.abs(result.center.values[d] - he.center.values[d]));
            }

            // FromCoefficients returns the eigenvalues in increasing order,
            // so the extents come back in decreasing order, and each axis is
            // determined only up to sign. Compare against the input sorted
            // the same way.
            const order = [...Array(n).keys()].sort(
                (a, b) => he.extent.values[b] - he.extent.values[a]);
            for (let d = 0; d < n; ++d) {
                const src = order[d];
                maxError = Math.max(maxError, Math.abs(
                    result.extent.values[d] - he.extent.values[src]));
                const s = dot(result.axis[d], he.axis[src]);
                maxError = Math.max(maxError, Math.abs(Math.abs(s) - 1));
            }

            // The recovered matrix M matches the original.
            maxError = Math.max(maxError,
                maxAbsDifference(result.getM(), he.getM()));
        }
        expect(maxError).toBeLessThan(1e-8);
    });

    it('round-trips the (A,B,C) form', () => {
        const rand = makeRandom(0xD15EA5E);
        let maxError = 0;
        for (let trial = 0; trial < 60; ++trial) {
            const n = (trial % 2 === 0 ? 2 : 3);
            const he = randomHyperellipsoid(n, rand);
            const { A, B, C } = he.toCoefficientsABC();
            const result = new Hyperellipsoid(n);
            expect(result.fromCoefficientsABC(A, B, C)).toBe(true);
            maxError = Math.max(maxError,
                maxAbsDifference(result.getM(), he.getM()));
            for (let d = 0; d < n; ++d) {
                maxError = Math.max(maxError,
                    Math.abs(result.center.values[d] - he.center.values[d]));
            }
        }
        expect(maxError).toBeLessThan(1e-9);
    });

    it('returns false for coefficients that are not a hyperellipsoid', () => {
        // A is singular (the rank-1 form x^2 = 1), so the center cannot be
        // computed.
        const singular = new Hyperellipsoid(2);
        expect(singular.fromCoefficientsABC(
            Matrix.fromArray(2, 2, [1, 0, 0, 0]),
            Vector.fromArray([0, 0]), -1)).toBe(false);

        // A is invertible but K^T*A*K - C = 0, the degenerate point case.
        const zeroRightSide = new Hyperellipsoid(2);
        expect(zeroRightSide.fromCoefficientsABC(
            Matrix.identity(2, 2), Vector.fromArray([0, 0]), 0)).toBe(false);

        // A has a nonpositive eigenvalue: the hyperbola x^2 - y^2 = 1.
        const hyperbola = new Hyperellipsoid(2);
        expect(hyperbola.fromCoefficientsABC(
            Matrix.fromArray(2, 2, [1, 0, 0, -1]),
            Vector.fromArray([0, 0]), -1)).toBe(false);
    });

    it('handles a zero extent the way upstream does', () => {
        // Upstream computes axis[d]/extent[d]; the GTE Vector division by
        // zero yields the zero vector, so the degenerate direction
        // contributes nothing to M and M is singular. M^{-1} likewise loses
        // the degenerate direction.
        const he = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([2, 0]));
        const M = he.getM();
        expect(M.get(0, 0)).toBeCloseTo(1 / 4, 12);
        expect(M.get(1, 1)).toBe(0);
        const MInv = he.getMInverse();
        expect(MInv.get(0, 0)).toBeCloseTo(4, 12);
        expect(MInv.get(1, 1)).toBe(0);

        // The quadratic form is then not invertible, so FromCoefficients
        // reports failure.
        const { A, B, C } = he.toCoefficientsABC();
        expect(new Hyperellipsoid(2).fromCoefficientsABC(A, B, C)).toBe(false);
    });

    it('produces NaN coefficients when all extents are zero', () => {
        // Upstream divides by the largest-magnitude quadratic coefficient
        // without checking for zero; with M = 0 that divisor is zero. The
        // quirk is preserved.
        const he = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
            Vector.fromArray([0, 0]));
        const coeff = he.toCoefficients();
        expect(coeff[3]).toBeNaN();
        expect(coeff[5]).toBeNaN();
    });

    it('orders hyperellipsoids for sorted containers', () => {
        const make = (cx: number, e0: number): Hyperellipsoid =>
            Hyperellipsoid.fromCenterAxisExtent(
                Vector.fromArray([cx, 0]),
                [Vector.fromArray([1, 0]), Vector.fromArray([0, 1])],
                Vector.fromArray([e0, 1]));

        const a = make(0, 1);
        const b = make(1, 1);
        const c = make(0, 2);
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        expect(a.lessThan(c)).toBe(true);
        expect(a.lessThanOrEqual(a.clone())).toBe(true);
        expect(a.greaterThan(b)).toBe(false);
        expect(b.greaterThan(a)).toBe(true);
        expect(a.greaterThanOrEqual(a.clone())).toBe(true);

        // Axis differences break ties between equal centers.
        const d = Hyperellipsoid.fromCenterAxisExtent(
            Vector.fromArray([0, 0]),
            [Vector.fromArray([0, 1]), Vector.fromArray([-1, 0])],
            Vector.fromArray([1, 1]));
        expect(d.lessThan(a)).toBe(true);
        expect(a.lessThan(d)).toBe(false);
    });

    it('computes the number of coefficients', () => {
        expect(hyperellipsoidNumCoefficients(2)).toBe(6);
        expect(hyperellipsoidNumCoefficients(3)).toBe(10);
        expect(hyperellipsoidNumCoefficients(4)).toBe(15);
    });
});
