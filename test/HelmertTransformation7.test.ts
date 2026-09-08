import { describe, it, expect } from 'vitest';
import { HelmertTransformation7 } from '../src/HelmertTransformation7.js';
import { Matrix, mulMatrix, multiplyAB, transpose } from '../src/Matrix.js';
import { determinant3x3 } from '../src/Matrix3x3.js';
import { Vector, add, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

// The rotation parametrization used by HelmertTransformation7 is
// R = Rz(a0) * Ry(a1) * Rx(a2), matching the three Euler-angle update steps.
function rotZ(t: number): Matrix {
    const m = new Matrix(3, 3);
    m.set(0, 0, Math.cos(t)); m.set(0, 1, -Math.sin(t)); m.set(0, 2, 0);
    m.set(1, 0, Math.sin(t)); m.set(1, 1, Math.cos(t)); m.set(1, 2, 0);
    m.set(2, 0, 0); m.set(2, 1, 0); m.set(2, 2, 1);
    return m;
}

function rotY(t: number): Matrix {
    const m = new Matrix(3, 3);
    m.set(0, 0, Math.cos(t)); m.set(0, 1, 0); m.set(0, 2, Math.sin(t));
    m.set(1, 0, 0); m.set(1, 1, 1); m.set(1, 2, 0);
    m.set(2, 0, -Math.sin(t)); m.set(2, 1, 0); m.set(2, 2, Math.cos(t));
    return m;
}

function rotX(t: number): Matrix {
    const m = new Matrix(3, 3);
    m.set(0, 0, 1); m.set(0, 1, 0); m.set(0, 2, 0);
    m.set(1, 0, 0); m.set(1, 1, Math.cos(t)); m.set(1, 2, -Math.sin(t));
    m.set(2, 0, 0); m.set(2, 1, Math.sin(t)); m.set(2, 2, Math.cos(t));
    return m;
}

function makeRotation(a0: number, a1: number, a2: number): Matrix {
    return multiplyAB(multiplyAB(rotZ(a0), rotY(a1)), rotX(a2));
}

function makePRNG(seed: number): () => number {
    let s = seed >>> 0;
    return (): number => {
        s = (s + 0x6d2b79f5) >>> 0;
        let t = s;
        t = Math.imul(t ^ (t >>> 15), t | 1);
        t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
        return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
}

// The 'q' points; the 'p' points are s*R*q + t.
function makeQPoints(count: number, rand: () => number): Vector[] {
    const q: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        q.push(v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2));
    }
    return q;
}

function applyTransform(q: readonly Vector[], scale: number, rotate: Matrix,
    translate: Vector): Vector[] {
    return q.map(v => add(mul(mulMatrix(rotate, v), scale), translate));
}

function maxResidual(p: readonly Vector[], q: readonly Vector[], scale: number,
    rotate: Matrix, translate: Vector): number {
    let maxLen = 0;
    for (let i = 0; i < p.length; ++i) {
        const mapped = add(mul(mulMatrix(rotate, q[i]), scale), translate);
        maxLen = Math.max(maxLen, length(sub(mapped, p[i])));
    }
    return maxLen;
}

describe('HelmertTransformation7', () => {
    it('rejects invalid input', () => {
        const helmert = new HelmertTransformation7();
        const few: Vector[] = [];
        for (let i = 0; i < 6; ++i) {
            few.push(v3(i, 0, 0));
        }
        expect(() => helmert.execute(few, few, 4)).toThrow(/Invalid input/);

        const seven: Vector[] = [];
        for (let i = 0; i < 7; ++i) {
            seven.push(v3(i, 0, 0));
        }
        expect(() => helmert.execute(seven, few, 4)).toThrow(/Invalid input/);
    });

    it('returns the identity transformation when p equals q', () => {
        const rand = makePRNG(7);
        const q = makeQPoints(20, rand);
        const helmert = new HelmertTransformation7();
        const result = helmert.execute(q, q, 64);

        expect(result.scale).toBeCloseTo(1, 10);
        expect(result.functionValue).toBeCloseTo(0, 12);
        for (let i = 0; i < 3; ++i) {
            expect(result.translate.values[i]).toBeCloseTo(0, 8);
            for (let j = 0; j < 3; ++j) {
                expect(result.rotate.get(i, j)).toBeCloseTo(i === j ? 1 : 0, 8);
            }
        }
        // The first sweep cannot improve on F = 0, so the loop breaks
        // immediately.
        expect(result.iterations).toBe(0);
    });

    it('recovers a known scale, rotation and translation', () => {
        const rand = makePRNG(1234);
        const q = makeQPoints(40, rand);
        const scale = 1.75;
        const rotate = makeRotation(0.35, -0.22, 0.17);
        const translate = v3(1.5, -2.25, 0.75);
        const p = applyTransform(q, scale, rotate, translate);

        const helmert = new HelmertTransformation7();
        const result = helmert.execute(p, q, 200);

        expect(result.scale).toBeCloseTo(scale, 8);
        expect(result.functionValue).toBeLessThan(1e-16);
        for (let i = 0; i < 3; ++i) {
            expect(result.translate.values[i]).toBeCloseTo(translate.values[i], 7);
            for (let j = 0; j < 3; ++j) {
                expect(result.rotate.get(i, j)).toBeCloseTo(rotate.get(i, j), 7);
            }
        }
        expect(maxResidual(p, q, result.scale, result.rotate, result.translate))
            .toBeLessThan(1e-7);
    });

    it('produces an orthogonal rotation with determinant 1', () => {
        const rand = makePRNG(99);
        const q = makeQPoints(25, rand);
        const p = applyTransform(q, 0.6, makeRotation(-0.4, 0.3, 0.5), v3(-3, 4, 5));
        const helmert = new HelmertTransformation7();
        const R = helmert.execute(p, q, 150).rotate;

        // R^T * R = I.
        for (let i = 0; i < 3; ++i) {
            for (let j = 0; j < 3; ++j) {
                let sum = 0;
                for (let k = 0; k < 3; ++k) {
                    sum += R.get(k, i) * R.get(k, j);
                }
                expect(sum).toBeCloseTo(i === j ? 1 : 0, 10);
            }
        }
        const det =
            R.get(0, 0) * (R.get(1, 1) * R.get(2, 2) - R.get(1, 2) * R.get(2, 1)) -
            R.get(0, 1) * (R.get(1, 0) * R.get(2, 2) - R.get(1, 2) * R.get(2, 0)) +
            R.get(0, 2) * (R.get(1, 0) * R.get(2, 1) - R.get(1, 1) * R.get(2, 0));
        expect(det).toBeCloseTo(1, 10);
    });

    it('decreases the objective function monotonically with more iterations', () => {
        const rand = makePRNG(2026);
        const q = makeQPoints(30, rand);
        const p = applyTransform(q, 2.5, makeRotation(0.5, 0.4, -0.3), v3(0.5, 0.5, 0.5));

        let previous = Number.MAX_VALUE;
        for (const numIterations of [0, 1, 2, 4, 8, 16, 32, 64]) {
            const result = new HelmertTransformation7().execute(p, q, numIterations);
            expect(result.functionValue).toBeLessThanOrEqual(previous + 1e-15);
            expect(result.iterations).toBeLessThanOrEqual(numIterations);
            previous = result.functionValue;
        }
        expect(previous).toBeLessThan(1e-16);
    });

    it('fits noisy data with a small residual and is invariant to a rigid ' +
        'change of the q-frame', () => {
        const rand = makePRNG(555);
        const q = makeQPoints(50, rand);
        const scale = 1.2;
        const rotate = makeRotation(0.2, 0.1, -0.15);
        const translate = v3(2, -1, 3);
        const noiseScale = 1e-3;
        const p = applyTransform(q, scale, rotate, translate).map(v =>
            add(v, v3(noiseScale * (2 * rand() - 1), noiseScale * (2 * rand() - 1),
                noiseScale * (2 * rand() - 1))));

        const result = new HelmertTransformation7().execute(p, q, 200);
        expect(result.scale).toBeCloseTo(scale, 3);
        expect(Math.sqrt(result.functionValue)).toBeLessThan(1e-2);

        // Pre-composing the q-points with a known rigid motion S must give a
        // transformation whose composition maps the moved points to the same
        // places.
        const S = makeRotation(0.9, -0.6, 0.25);
        const sTranslate = v3(-4, 8, 1.5);
        const qMoved = applyTransform(q, 1, S, sTranslate);
        const moved = new HelmertTransformation7().execute(p, qMoved, 400);
        expect(moved.scale).toBeCloseTo(result.scale, 6);
        for (let i = 0; i < p.length; ++i) {
            const a = add(mul(mulMatrix(result.rotate, q[i]), result.scale),
                result.translate);
            const b = add(mul(mulMatrix(moved.rotate, qMoved[i]), moved.scale),
                moved.translate);
            expect(length(sub(a, b))).toBeLessThan(1e-6);
        }
    });

    it('cross-checks randomized exact-fit cases', () => {
        const rand = makePRNG(31337);
        for (let trial = 0; trial < 12; ++trial) {
            const q = makeQPoints(15, rand);
            const scale = 0.5 + 2 * rand();
            const rotate = makeRotation(0.6 * (2 * rand() - 1), 0.6 * (2 * rand() - 1),
                0.6 * (2 * rand() - 1));
            const translate = v3(10 * rand() - 5, 10 * rand() - 5, 10 * rand() - 5);
            const p = applyTransform(q, scale, rotate, translate);

            const result = new HelmertTransformation7().execute(p, q, 300);
            expect(result.scale).toBeCloseTo(scale, 6);
            expect(maxResidual(p, q, result.scale, result.rotate, result.translate))
                .toBeLessThan(1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V41): properties of the 7-parameter Helmert fit.
// ---------------------------------------------------------------------------

const helmertCoord = () => scaled(-5, 5, 512);

// n >= 7 correspondences of arbitrary (unrelated) points.
const pointSets = (minCount = 7, maxCount = 12) =>
    fc.integer({ min: minCount, max: maxCount }).chain(n => fc.tuple(
        fc.array(fc.array(helmertCoord(), { minLength: 3, maxLength: 3 }),
            { minLength: n, maxLength: n }),
        fc.array(fc.array(helmertCoord(), { minLength: 3, maxLength: 3 }),
            { minLength: n, maxLength: n })))
        .map(([p, q]) => ({
            p: p.map(a => Vector.fromArray(a)),
            q: q.map(a => Vector.fromArray(a))
        }))
        // The scale is numer/denom with denom = sum |q_i - qAverage|^2, so
        // reject the degenerate case of coincident q-points.
        .filter(({ q }) => {
            const c = q.reduce((s, v) => add(s, v), Vector.zero(3));
            const avg = mul(c, 1 / q.length);
            return q.reduce((s, v) => s + length(sub(v, avg)) ** 2, 0) > 1e-3;
        });

// The mean squared residual of the returned similarity, computed
// independently of the class: (1/n) * sum |s*R*q_i + t - p_i|^2.
function meanSquaredResidual(p: readonly Vector[], q: readonly Vector[],
    scale: number, rotate: Matrix, translate: Vector): number {
    let sum = 0;
    for (let i = 0; i < p.length; ++i) {
        const term = sub(add(mul(mulMatrix(rotate, q[i]), scale), translate),
            p[i]);
        sum += term.values[0] ** 2 + term.values[1] ** 2 + term.values[2] ** 2;
    }
    return sum / p.length;
}

function maxAbs(A: Matrix, B: Matrix): number {
    let m = 0;
    for (let i = 0; i < A.numElements; ++i) {
        m = Math.max(m, Math.abs(A.values[i] - B.values[i]));
    }
    return m;
}

describe('HelmertTransformation7 verification', () => {
    it('the reported function value is the mean squared residual of the ' +
        'reported scale, rotation and translation', () => {
            check(fc.tuple(pointSets(), fc.integer({ min: 0, max: 40 })),
                ([{ p, q }, numIterations]) => {
                    const result = new HelmertTransformation7()
                        .execute(p, q, numIterations);
                    const reference = meanSquaredResidual(p, q, result.scale,
                        result.rotate, result.translate);
                    // The magnitudes of the terms set the absolute tolerance:
                    // the residual is a sum of squares of O(10) quantities.
                    let scaleOfTerms = 1;
                    for (let i = 0; i < p.length; ++i) {
                        scaleOfTerms = Math.max(scaleOfTerms,
                            length(p[i]) ** 2,
                            (result.scale * length(q[i])) ** 2);
                    }
                    expectClose(result.functionValue, reference,
                        1e-11 * scaleOfTerms, 1e-9);
                }, 120);
        });

    it('returns a proper rotation (orthogonal with determinant 1)', () => {
        check(fc.tuple(pointSets(), fc.integer({ min: 0, max: 30 })),
            ([{ p, q }, numIterations]) => {
                const R = new HelmertTransformation7()
                    .execute(p, q, numIterations).rotate;
                expect(maxAbs(multiplyAB(R, transpose(R)),
                    Matrix.identity(3, 3))).toBeLessThan(1e-12);
                expectClose(determinant3x3(R), 1, 1e-12, 1e-12);
            }, 120);
    });

    it('the objective decreases monotonically in the iteration count and ' +
        'the iteration count is the loop index at termination', () => {
            check(pointSets(), ({ p, q }) => {
                let previous = Number.POSITIVE_INFINITY;
                let converged = -1;
                for (let k = 0; k <= 12; ++k) {
                    const result = new HelmertTransformation7()
                        .execute(p, q, k);
                    expect(result.iterations).toBeLessThanOrEqual(k);
                    expect(result.functionValue)
                        .toBeLessThanOrEqual(previous + 1e-14);
                    previous = result.functionValue;
                    if (result.iterations < k && converged < 0) {
                        converged = result.iterations;
                    }
                }
                if (converged >= 0) {
                    // Once the coordinate descent stops improving, running
                    // more iterations changes nothing.
                    const a = new HelmertTransformation7()
                        .execute(p, q, converged + 1);
                    const b = new HelmertTransformation7()
                        .execute(p, q, 60);
                    expect(b.iterations).toBe(a.iterations);
                    expect(b.functionValue).toBe(a.functionValue);
                    expect(maxAbs(a.rotate, b.rotate)).toBe(0);
                }
            }, 60);
        });

    it('recovers a known similarity from exact correspondences', () => {
        const angle = () => scaled(-0.5, 0.5, 256);
        const arb = fc.tuple(
            fc.array(fc.array(scaled(-2, 2, 256),
                { minLength: 3, maxLength: 3 }),
                { minLength: 9, maxLength: 9 }),
            angle(), angle(), angle(), scaled(0.4, 2.5, 128),
            fc.array(scaled(-4, 4, 128), { minLength: 3, maxLength: 3 }))
            .filter(([qRaw]) => {
                // Reject nearly collinear/coincident q-point clouds.
                const q = qRaw.map(a => Vector.fromArray(a));
                const c = q.reduce((s, v) => add(s, v), Vector.zero(3));
                const avg = mul(c, 1 / q.length);
                return q.reduce((s, v) => s + length(sub(v, avg)) ** 2, 0) > 1;
            });
        check(arb, ([qRaw, a0, a1, a2, scale, tRaw]) => {
            const q = qRaw.map(a => Vector.fromArray(a));
            const R = makeRotation(a0, a1, a2);
            const t = Vector.fromArray(tRaw);
            const p = applyTransform(q, scale, R, t);
            const result = new HelmertTransformation7().execute(p, q, 400);
            // The fit is exact, so the objective is at machine zero and every
            // parameter is recovered.
            expect(result.functionValue).toBeLessThan(1e-16);
            expectClose(result.scale, scale, 1e-7, 1e-7);
            expect(maxAbs(result.rotate, R)).toBeLessThan(1e-6);
            expect(length(sub(result.translate, t))).toBeLessThan(1e-5);
        }, 60);
    });

    it('is invariant under a permutation of the correspondences', () => {
        const arb = fc.tuple(pointSets(8, 10),
            fc.integer({ min: 0, max: 0xffff }));
        check(arb, ([{ p, q }, seed]) => {
            const rand = makePRNG(seed + 1);
            const order = p.map((_, i) => i);
            for (let i = order.length - 1; i > 0; --i) {
                const j = Math.floor(rand() * (i + 1));
                [order[i], order[j]] = [order[j], order[i]];
            }
            const base = new HelmertTransformation7().execute(p, q, 200);
            const shuffled = new HelmertTransformation7().execute(
                order.map(i => p[i]), order.map(i => q[i]), 200);
            // The attained minimum and the scale agree to a few ulps. The
            // rotation and translation are looser: the descent stops as
            // soon as no Euler-angle update lowers F, so a change
            // of summation order can stop it a half step earlier along a
            // flat direction of the objective (see VERIFYING.md on
            // path-dependent minimizers).
            expectClose(shuffled.scale, base.scale, 1e-12, 1e-12);
            expectClose(shuffled.functionValue, base.functionValue,
                1e-12, 1e-11);
            expect(maxAbs(shuffled.rotate, base.rotate)).toBeLessThan(1e-5);
            expect(length(sub(shuffled.translate, base.translate)))
                .toBeLessThan(1e-4);
        }, 40);
    });

    it('is equivariant under translation of either point set', () => {
        const arb = fc.tuple(pointSets(8, 10),
            fc.array(scaled(-3, 3, 128), { minLength: 3, maxLength: 3 }));
        check(arb, ([{ p, q }, cRaw]) => {
            const c = Vector.fromArray(cRaw);
            const base = new HelmertTransformation7().execute(p, q, 200);

            // p -> p + c: the fit is unchanged except translate -> t + c,
            // because the u-values and the working translation both shift by
            // c and the v-values sum to zero.
            const shiftedP = new HelmertTransformation7()
                .execute(p.map(v => add(v, c)), q, 200);
            expectClose(shiftedP.scale, base.scale, 1e-12, 1e-12);
            expectClose(shiftedP.functionValue, base.functionValue,
                1e-12, 1e-11);
            expect(maxAbs(shiftedP.rotate, base.rotate)).toBeLessThan(1e-5);
            expect(length(sub(shiftedP.translate,
                add(base.translate, c)))).toBeLessThan(1e-4);

            // q -> q + c: translate -> t - s*R*c so that the composed map is
            // unchanged on the shifted points.
            const shiftedQ = new HelmertTransformation7()
                .execute(p, q.map(v => add(v, c)), 200);
            expectClose(shiftedQ.scale, base.scale, 1e-12, 1e-12);
            expectClose(shiftedQ.functionValue, base.functionValue,
                1e-12, 1e-11);
            expect(maxAbs(shiftedQ.rotate, base.rotate)).toBeLessThan(1e-5);
            const expectedT = sub(base.translate,
                mul(mulMatrix(base.rotate, c), base.scale));
            expect(length(sub(shiftedQ.translate, expectedT)))
                .toBeLessThan(1e-4);
        }, 40);
    });

    it('preserves the upstream requirement of at least 7 correspondences',
        () => {
            check(fc.integer({ min: 0, max: 9 }), (n) => {
                const p: Vector[] = [];
                const q: Vector[] = [];
                for (let i = 0; i < n; ++i) {
                    p.push(v3(i, i * i, 1));
                    q.push(v3(2 * i, i, -i));
                }
                const run = () =>
                    new HelmertTransformation7().execute(p, q, 5);
                if (n < 7) {
                    expect(run).toThrow('Invalid input.');
                }
                else {
                    expect(run).not.toThrow();
                }
            });
        });
});
