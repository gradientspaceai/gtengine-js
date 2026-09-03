import { describe, it, expect } from 'vitest';
import { IntpVectorField2 } from '../src/IntpVectorField2';
import { Vector } from '../src/Vector';

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

function gridPoints(n: number, jitter: number): Vector[] {
    const rand = makeRandom(2024);
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        for (let j = 0; j < n; ++j) {
            const x = i / (n - 1);
            const y = j / (n - 1);
            const onBoundary = (i === 0 || i === n - 1 || j === 0 || j === n - 1);
            const dx = onBoundary ? 0 : jitter * (2 * rand() - 1);
            const dy = onBoundary ? 0 : jitter * (2 * rand() - 1);
            points.push(Vector.fromArray([x + dx, y + dy]));
        }
    }
    return points;
}

// An affine map of the plane; the Cendes-Wong interpolator reproduces affine
// data exactly, so the vector field interpolation is exact for this map.
function affine(p: Vector): Vector {
    const x = p.values[0];
    const y = p.values[1];
    return Vector.fromArray([2 + 3 * x - y, -1 + 0.5 * x + 4 * y]);
}

function inverseAffine(q: Vector): Vector {
    // Invert [[3,-1],[0.5,4]] * (x,y) = (u - 2, v + 1).
    const u = q.values[0] - 2;
    const v = q.values[1] + 1;
    const det = 3 * 4 - (-1) * 0.5;
    return Vector.fromArray([(4 * u + 1 * v) / det, (-0.5 * u + 3 * v) / det]);
}

describe('IntpVectorField2', () => {
    const domain = gridPoints(6, 0.03);
    const range = domain.map(affine);

    it('reproduces an affine vector field exactly', () => {
        const field = new IntpVectorField2(domain, range);
        const rand = makeRandom(31337);
        for (let k = 0; k < 50; ++k) {
            const p = Vector.fromArray([
                0.05 + 0.9 * rand(),
                0.05 + 0.9 * rand()
            ]);
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            const expected = affine(p);
            expect(result.output.values[0]).toBeCloseTo(expected.values[0], 9);
            expect(result.output.values[1]).toBeCloseTo(expected.values[1], 9);
        }
    });

    it('interpolates the data at the domain points', () => {
        // A nonlinear map: the interpolant must still pass through the data.
        const nonlinear = (p: Vector) => Vector.fromArray([
            Math.sin(2 * p.values[0]) + p.values[1],
            p.values[0] * p.values[0] - Math.cos(p.values[1])
        ]);
        const field = new IntpVectorField2(domain, domain.map(nonlinear));
        for (const p of domain) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            const expected = nonlinear(p);
            expect(result.output.values[0]).toBeCloseTo(expected.values[0], 10);
            expect(result.output.values[1]).toBeCloseTo(expected.values[1], 10);
        }
    });

    it('round-trips through the inverse affine field', () => {
        const forward = new IntpVectorField2(domain, range);
        const backward = new IntpVectorField2(range, domain.map(p => p.clone()));
        const rand = makeRandom(5150);
        for (let k = 0; k < 30; ++k) {
            const p = Vector.fromArray([
                0.1 + 0.8 * rand(),
                0.1 + 0.8 * rand()
            ]);
            const q = forward.evaluate(p);
            expect(q.valid).toBe(true);
            // The forward image agrees with the analytic map, and the
            // analytic inverse recovers the input.
            const back = inverseAffine(q.output);
            expect(back.values[0]).toBeCloseTo(p.values[0], 8);
            expect(back.values[1]).toBeCloseTo(p.values[1], 8);

            // The interpolated inverse field also recovers the input.
            const r = backward.evaluate(q.output);
            expect(r.valid).toBe(true);
            expect(r.output.values[0]).toBeCloseTo(p.values[0], 8);
            expect(r.output.values[1]).toBeCloseTo(p.values[1], 8);
        }
    });

    it('reproduces the identity map', () => {
        const field = new IntpVectorField2(domain, domain.map(p => p.clone()));
        for (const p of [
            Vector.fromArray([0.3, 0.7]),
            Vector.fromArray([0.55, 0.21]),
            Vector.fromArray([0.5, 0.5])
        ]) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.output.values[0]).toBeCloseTo(p.values[0], 9);
            expect(result.output.values[1]).toBeCloseTo(p.values[1], 9);
        }
    });

    it('reports invalid outside the convex hull of the domain', () => {
        const field = new IntpVectorField2(domain, range);
        for (const p of [
            Vector.fromArray([-1, 0.5]),
            Vector.fromArray([0.5, 2]),
            Vector.fromArray([5, 5])
        ]) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(false);
        }
    });

    it('validates its input', () => {
        expect(() => new IntpVectorField2([], [])).toThrow();
        expect(() => new IntpVectorField2(domain, range.slice(1))).toThrow();
        expect(() => new IntpVectorField2(
            [Vector.fromArray([0, 0, 0]), Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0])],
            [Vector.fromArray([0, 0, 0]), Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0])])).toThrow();
    });
});
