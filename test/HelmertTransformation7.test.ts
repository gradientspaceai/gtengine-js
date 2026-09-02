import { describe, it, expect } from 'vitest';
import { HelmertTransformation7 } from '../src/HelmertTransformation7';
import { Matrix, mulMatrix, multiplyAB } from '../src/Matrix';
import { Vector, add, length, mul, sub } from '../src/Vector';

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
