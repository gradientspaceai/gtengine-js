import { describe, it, expect } from 'vitest';
import {
    Quaternion, negateQuaternion, addQuaternion, subQuaternion,
    mulQuaternion, divQuaternion, conjugate, inverseQuaternion, rotate,
    slerpQuaternion
} from '../src/Quaternion.js';
import { Vector, dot, length, normalize } from '../src/Vector.js';
import { slerp, slerpUsingCosAngle, slerpUsingMidpoint } from '../src/Slerp.js';
import { Matrix, mulMatrix } from '../src/Matrix.js';

function expectQuaternionClose(actual: Vector, expected: readonly number[],
    tolerance: number = 1e-13): void {
    expect(actual.size).toBe(4);
    for (let i = 0; i < 4; ++i) {
        expect(Math.abs(actual.values[i] - expected[i]))
            .toBeLessThanOrEqual(tolerance);
    }
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296 - 0.5;
    };
}

function randomQuaternion(rand: () => number): Quaternion {
    return new Quaternion(4 * rand(), 4 * rand(), 4 * rand(), 4 * rand());
}

function randomUnitQuaternion(rand: () => number): Quaternion {
    const q = randomQuaternion(rand);
    normalize(q);
    return Quaternion.fromArray(q.values);
}

// The unit quaternion for a rotation of 'angle' radians about the
// unit-length 'axis'.
function axisAngle(axis: readonly number[], angle: number): Quaternion {
    const h = 0.5 * angle;
    const s = Math.sin(h);
    return new Quaternion(s * axis[0], s * axis[1], s * axis[2], Math.cos(h));
}

// An independent rotation matrix from the axis-angle (Rodrigues formula),
// used to cross-check rotate().
function rodrigues(axis: readonly number[], angle: number): Matrix {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const R = new Matrix(3, 3);
    for (let r = 0; r < 3; ++r) {
        for (let col = 0; col < 3; ++col) {
            R.set(r, col, (1 - c) * axis[r] * axis[col] + (r === col ? c : 0));
        }
    }
    R.set(0, 1, R.get(0, 1) - s * axis[2]);
    R.set(0, 2, R.get(0, 2) + s * axis[1]);
    R.set(1, 0, R.get(1, 0) + s * axis[2]);
    R.set(1, 2, R.get(1, 2) - s * axis[0]);
    R.set(2, 0, R.get(2, 0) - s * axis[1]);
    R.set(2, 1, R.get(2, 1) + s * axis[0]);
    return R;
}

describe('Quaternion', () => {
    describe('construction and special quaternions', () => {
        it('stores the tuple in the order (x,y,z,w)', () => {
            const q = new Quaternion(1, 2, 3, 4);
            expect(q.size).toBe(4);
            expect(q.values).toEqual([1, 2, 3, 4]);
            expect(q.get(0)).toBe(1);
            expect(q.get(3)).toBe(4);
            q.set(1, 9);
            expect(q.values[1]).toBe(9);
        });

        it('default-constructs to zero', () => {
            expect(new Quaternion().values).toEqual([0, 0, 0, 0]);
        });

        it('provides the special quaternions', () => {
            expect(Quaternion.zero().values).toEqual([0, 0, 0, 0]);
            expect(Quaternion.i().values).toEqual([1, 0, 0, 0]);
            expect(Quaternion.j().values).toEqual([0, 1, 0, 0]);
            expect(Quaternion.k().values).toEqual([0, 0, 1, 0]);
            expect(Quaternion.identity().values).toEqual([0, 0, 0, 1]);
        });

        it('supports fromArray, clone and the inherited comparisons', () => {
            const q = Quaternion.fromArray([1, 2, 3, 4]);
            expect(q).toBeInstanceOf(Quaternion);
            const c = q.clone();
            expect(c).toBeInstanceOf(Quaternion);
            expect(c.equals(q)).toBe(true);
            c.set(0, 0);
            expect(q.values[0]).toBe(1);
            expect(c.lessThan(q)).toBe(true);
            expect(q.greaterThan(c)).toBe(true);
            expect(q.notEquals(c)).toBe(true);
            expect(() => Quaternion.fromArray([1, 2, 3])).toThrow();
        });

        it('is a Vector, so dot/length/normalize apply', () => {
            const q = new Quaternion(1, 2, 3, 4);
            expect(q).toBeInstanceOf(Vector);
            expect(dot(q, q)).toBe(30);
            expect(length(q)).toBeCloseTo(Math.sqrt(30), 14);
            const n = new Quaternion(0, 0, 0, 5);
            expect(normalize(n)).toBe(5);
            expect(n.values).toEqual([0, 0, 0, 1]);
            const z = new Quaternion(0, 0, 0, 0);
            expect(normalize(z)).toBe(0);
            expect(z.values).toEqual([0, 0, 0, 0]);
        });
    });

    describe('linear algebraic operations', () => {
        it('negates, adds and subtracts componentwise', () => {
            const q0 = new Quaternion(1, 2, 3, 4);
            const q1 = new Quaternion(5, 6, 7, 8);
            expect(negateQuaternion(q0).values).toEqual([-1, -2, -3, -4]);
            expect(addQuaternion(q0, q1).values).toEqual([6, 8, 10, 12]);
            expect(subQuaternion(q1, q0).values).toEqual([4, 4, 4, 4]);
            expect(addQuaternion(q0, q1)).toBeInstanceOf(Quaternion);
        });

        it('multiplies and divides by a scalar', () => {
            const q = new Quaternion(1, 2, 3, 4);
            expect(mulQuaternion(q, 2).values).toEqual([2, 4, 6, 8]);
            expect(mulQuaternion(2, q).values).toEqual([2, 4, 6, 8]);
            expect(divQuaternion(q, 2).values).toEqual([0.5, 1, 1.5, 2]);
            // Division by zero produces the zero quaternion.
            expect(divQuaternion(q, 0).values).toEqual([0, 0, 0, 0]);
        });
    });

    describe('the Hamilton product', () => {
        it('satisfies i*i = j*j = k*k = i*j*k = -1', () => {
            const i = Quaternion.i();
            const j = Quaternion.j();
            const k = Quaternion.k();
            const minusOne = [0, 0, 0, -1];
            expectQuaternionClose(mulQuaternion(i, i), minusOne);
            expectQuaternionClose(mulQuaternion(j, j), minusOne);
            expectQuaternionClose(mulQuaternion(k, k), minusOne);
            expectQuaternionClose(
                mulQuaternion(mulQuaternion(i, j), k), minusOne);
        });

        it('satisfies i*j = k, j*k = i, k*i = j and the anti-commutations',
            () => {
                const i = Quaternion.i();
                const j = Quaternion.j();
                const k = Quaternion.k();
                expectQuaternionClose(mulQuaternion(i, j), [0, 0, 1, 0]);
                expectQuaternionClose(mulQuaternion(j, k), [1, 0, 0, 0]);
                expectQuaternionClose(mulQuaternion(k, i), [0, 1, 0, 0]);
                expectQuaternionClose(mulQuaternion(j, i), [0, 0, -1, 0]);
                expectQuaternionClose(mulQuaternion(k, j), [-1, 0, 0, 0]);
                expectQuaternionClose(mulQuaternion(i, k), [0, -1, 0, 0]);
            });

        it('has 1 as the multiplicative identity and is not commutative',
            () => {
                const one = Quaternion.identity();
                const q = new Quaternion(1, 2, 3, 4);
                expectQuaternionClose(mulQuaternion(q, one), q.values);
                expectQuaternionClose(mulQuaternion(one, q), q.values);
                const p = new Quaternion(-2, 5, 1, 3);
                expect(mulQuaternion(q, p).equals(mulQuaternion(p, q)))
                    .toBe(false);
            });

        it('is associative and distributes over addition', () => {
            const rand = makeRandom(3141);
            for (let trial = 0; trial < 30; ++trial) {
                const a = randomQuaternion(rand);
                const b = randomQuaternion(rand);
                const c = randomQuaternion(rand);
                expectQuaternionClose(
                    mulQuaternion(mulQuaternion(a, b), c),
                    mulQuaternion(a, mulQuaternion(b, c)).values, 1e-12);
                expectQuaternionClose(
                    mulQuaternion(a, addQuaternion(b, c)),
                    addQuaternion(mulQuaternion(a, b),
                        mulQuaternion(a, c)).values, 1e-12);
            }
        });

        it('is multiplicative on the norm: |q0*q1| = |q0|*|q1|', () => {
            const rand = makeRandom(2718);
            for (let trial = 0; trial < 30; ++trial) {
                const a = randomQuaternion(rand);
                const b = randomQuaternion(rand);
                expect(length(mulQuaternion(a, b)))
                    .toBeCloseTo(length(a) * length(b), 12);
            }
        });
    });

    describe('conjugate and inverse', () => {
        it('conjugates the vector part', () => {
            expect(conjugate(new Quaternion(1, 2, 3, 4)).values)
                .toEqual([-1, -2, -3, 4]);
        });

        it('satisfies conj(q0*q1) = conj(q1)*conj(q0)', () => {
            const rand = makeRandom(1618);
            for (let trial = 0; trial < 30; ++trial) {
                const a = randomQuaternion(rand);
                const b = randomQuaternion(rand);
                expectQuaternionClose(conjugate(mulQuaternion(a, b)),
                    mulQuaternion(conjugate(b), conjugate(a)).values, 1e-12);
            }
        });

        it('satisfies q*inv(q) = inv(q)*q = 1', () => {
            const rand = makeRandom(161803);
            for (let trial = 0; trial < 30; ++trial) {
                const q = randomQuaternion(rand);
                const inv = inverseQuaternion(q);
                expectQuaternionClose(mulQuaternion(q, inv), [0, 0, 0, 1],
                    1e-11);
                expectQuaternionClose(mulQuaternion(inv, q), [0, 0, 0, 1],
                    1e-11);
            }
        });

        it('the inverse of a unit quaternion is its conjugate', () => {
            const q = axisAngle([0, 0, 1], 0.7);
            expectQuaternionClose(inverseQuaternion(q), conjugate(q).values);
        });

        it('returns zero for the zero quaternion', () => {
            expect(inverseQuaternion(Quaternion.zero()).values)
                .toEqual([0, 0, 0, 0]);
        });

        it('matches inv(q) = conj(q)/|q|^2 for a known value', () => {
            // |q|^2 = 30 for q = (1,2,3,4).
            expectQuaternionClose(inverseQuaternion(new Quaternion(1, 2, 3, 4)),
                [-1 / 30, -2 / 30, -3 / 30, 4 / 30], 1e-16);
        });
    });

    describe('rotate', () => {
        it('rotates by 90 degrees about z', () => {
            const q = axisAngle([0, 0, 1], Math.PI / 2);
            const v = rotate(q, Vector.fromArray([1, 0, 0]));
            expect(v.values[0]).toBeCloseTo(0, 14);
            expect(v.values[1]).toBeCloseTo(1, 14);
            expect(v.values[2]).toBeCloseTo(0, 14);
        });

        it('handles the 180-degree rotations about each axis exactly', () => {
            // These are the rotations whose matrices have trace -1, i.e. the
            // degenerate branches of a matrix-to-quaternion conversion; each
            // has a zero w component.
            const u = Vector.fromArray([1, 2, 3]);
            const cases: Array<[Quaternion, number[]]> = [
                [Quaternion.i(), [1, -2, -3]],
                [Quaternion.j(), [-1, 2, -3]],
                [Quaternion.k(), [-1, -2, 3]]
            ];
            for (const [q, expected] of cases) {
                const v = rotate(q, u);
                for (let i = 0; i < 3; ++i) {
                    expect(v.values[i]).toBeCloseTo(expected[i], 15);
                }
            }
        });

        it('leaves the identity quaternion acting as the identity', () => {
            const u = Vector.fromArray([1, -2, 3]);
            expect(rotate(Quaternion.identity(), u).values)
                .toEqual([1, -2, 3]);
        });

        it('agrees with the Rodrigues rotation matrix', () => {
            const rand = makeRandom(90210);
            for (let trial = 0; trial < 40; ++trial) {
                const axisV = Vector.fromArray(
                    [4 * rand(), 4 * rand(), 4 * rand()]);
                if (length(axisV) < 1e-3) {
                    continue;
                }
                normalize(axisV);
                const angle = 6 * rand();
                const q = axisAngle(axisV.values, angle);
                const R = rodrigues(axisV.values, angle);
                const u = Vector.fromArray([4 * rand(), 4 * rand(), 4 * rand()]);
                const byQuaternion = rotate(q, u);
                const byMatrix = mulMatrix(R, u);
                for (let i = 0; i < 3; ++i) {
                    expect(byQuaternion.values[i])
                        .toBeCloseTo(byMatrix.values[i], 12);
                }
            }
        });

        it('agrees with q*(0,u)*conj(q), the original definition', () => {
            const rand = makeRandom(24680);
            for (let trial = 0; trial < 30; ++trial) {
                const q = randomUnitQuaternion(rand);
                const u = Vector.fromArray([4 * rand(), 4 * rand(), 4 * rand()]);
                const uq = new Quaternion(u.values[0], u.values[1],
                    u.values[2], 0);
                const p = mulQuaternion(mulQuaternion(q, uq), conjugate(q));
                const v = rotate(q, u);
                for (let i = 0; i < 3; ++i) {
                    expect(v.values[i]).toBeCloseTo(p.values[i], 12);
                }
                expect(p.values[3]).toBeCloseTo(0, 12);
            }
        });

        it('composes: rotate(q0*q1, u) = rotate(q0, rotate(q1, u))', () => {
            const rand = makeRandom(11235);
            for (let trial = 0; trial < 30; ++trial) {
                const q0 = randomUnitQuaternion(rand);
                const q1 = randomUnitQuaternion(rand);
                const u = Vector.fromArray([4 * rand(), 4 * rand(), 4 * rand()]);
                const a = rotate(mulQuaternion(q0, q1), u);
                const b = rotate(q0, rotate(q1, u));
                for (let i = 0; i < 3; ++i) {
                    expect(a.values[i]).toBeCloseTo(b.values[i], 12);
                }
            }
        });

        it('preserves length and gives the same result for q and -q', () => {
            const rand = makeRandom(31415);
            for (let trial = 0; trial < 30; ++trial) {
                const q = randomUnitQuaternion(rand);
                const u = Vector.fromArray([4 * rand(), 4 * rand(), 4 * rand()]);
                const v = rotate(q, u);
                expect(length(v)).toBeCloseTo(length(u), 12);
                const w = rotate(negateQuaternion(q), u);
                for (let i = 0; i < 3; ++i) {
                    expect(v.values[i]).toBeCloseTo(w.values[i], 12);
                }
            }
        });

        it('handles a homogeneous 4D vector, preserving the last component',
            () => {
                const q = axisAngle([0, 0, 1], Math.PI / 2);
                const u4 = Vector.fromArray([1, 0, 0, 0]);
                const v4 = rotate(q, u4);
                expect(v4.size).toBe(4);
                expect(v4.values[0]).toBeCloseTo(0, 14);
                expect(v4.values[1]).toBeCloseTo(1, 14);
                expect(v4.values[2]).toBeCloseTo(0, 14);
                expect(v4.values[3]).toBe(0);

                const rand = makeRandom(4711);
                for (let trial = 0; trial < 20; ++trial) {
                    const p = randomUnitQuaternion(rand);
                    const u = Vector.fromArray(
                        [4 * rand(), 4 * rand(), 4 * rand()]);
                    const h = Vector.fromArray(
                        [u.values[0], u.values[1], u.values[2], 0]);
                    const v3 = rotate(p, u);
                    const vh = rotate(p, h);
                    for (let i = 0; i < 3; ++i) {
                        expect(vh.values[i]).toBeCloseTo(v3.values[i], 13);
                    }
                    expect(vh.values[3]).toBe(0);
                }
            });

        it('rejects dimensions other than 3 and 4', () => {
            expect(() => rotate(Quaternion.identity(), new Vector(2)))
                .toThrow();
        });
    });

    describe('slerpQuaternion', () => {
        it('interpolates the endpoints', () => {
            const q0 = axisAngle([0, 0, 1], 0);
            const q1 = axisAngle([0, 0, 1], 1.0);
            expectQuaternionClose(slerpQuaternion(0, q0, q1), q0.values, 1e-14);
            expectQuaternionClose(slerpQuaternion(1, q0, q1), q1.values, 1e-14);
        });

        it('stays unit length and moves at constant angular speed', () => {
            const q0 = axisAngle([1, 0, 0], 0.2);
            const q1 = axisAngle([0, 1, 0], 1.1);
            const theta = Math.acos(dot(q0, q1));
            for (let k = 0; k <= 10; ++k) {
                const t = k / 10;
                const q = slerpQuaternion(t, q0, q1);
                expect(length(q)).toBeCloseTo(1, 12);
                // acos amplifies round-off near t = 0, so 6 digits.
                expect(Math.acos(Math.min(1, dot(q0, q))))
                    .toBeCloseTo(t * theta, 6);
            }
        });

        it('halfway is the normalized sum of the endpoints', () => {
            const q0 = axisAngle([0, 0, 1], 0.4);
            const q1 = axisAngle([0, 1, 0], 1.3);
            const mid = slerpQuaternion(0.5, q0, q1);
            const sum = addQuaternion(q0, q1);
            normalize(sum);
            expectQuaternionClose(mid, sum.values, 1e-12);
        });

        it('takes the short arc: slerp(t,q0,-q1) = slerp(t,q0,q1)', () => {
            const rand = makeRandom(60613);
            for (let trial = 0; trial < 30; ++trial) {
                const q0 = randomUnitQuaternion(rand);
                const q1 = randomUnitQuaternion(rand);
                for (const t of [0.25, 0.5, 0.75]) {
                    expectQuaternionClose(
                        slerpQuaternion(t, q0, negateQuaternion(q1)),
                        slerpQuaternion(t, q0, q1).values, 1e-13);
                }
            }
        });

        it('agrees with Slerp.ts slerp when the angle is acute', () => {
            const rand = makeRandom(1024);
            for (let trial = 0; trial < 30; ++trial) {
                const q0 = randomUnitQuaternion(rand);
                let q1 = randomUnitQuaternion(rand);
                if (dot(q0, q1) < 0) {
                    q1 = negateQuaternion(q1);
                }
                const t = 0.5 + 0.5 * rand();
                expectQuaternionClose(slerpQuaternion(t, q0, q1),
                    slerp(t, q0.values, q1.values), 1e-13);
                // The restricted/preprocessed forms of Quaternion.h are
                // Slerp.ts's slerpUsingCosAngle and slerpUsingMidpoint.
                const cosA = dot(q0, q1);
                expectQuaternionClose(slerpQuaternion(t, q0, q1),
                    slerpUsingCosAngle(t, q0.values, q1.values, cosA), 1e-13);
                const cosAH = Math.sqrt(0.5 * (1 + cosA));
                const qh = divQuaternion(addQuaternion(q0, q1), 2 * cosAH);
                expectQuaternionClose(slerpQuaternion(t, q0, q1),
                    slerpUsingMidpoint(t, q0.values, q1.values, qh.values,
                        cosAH), 1e-11);
            }
        });

        it('handles identical endpoints', () => {
            const q = axisAngle([0, 1, 0], 0.9);
            expectQuaternionClose(slerpQuaternion(0.3, q, q), q.values, 1e-13);
        });
    });
});
