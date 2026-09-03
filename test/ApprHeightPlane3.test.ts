import { describe, it, expect } from 'vitest';
import { ApprHeightPlane3 } from '../src/ApprHeightPlane3.js';
import { Vector } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprHeightPlane3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprHeightPlane3();
        const p = fitter.getParameters();
        expect(p.average.values).toEqual([0, 0, 0]);
        expect(p.coefficients.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(3);
    });

    it('recovers a plane from points that lie exactly on it', () => {
        // z = 3x - 2y + 5 sampled on a 4x4 grid.
        const points: Vector[] = [];
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 4; ++j) {
                const x = i - 1.5, y = j - 1.5;
                points.push(v3(x, y, 3 * x - 2 * y + 5));
            }
        }
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.average.values[0]).toBeCloseTo(0, 12);
        expect(p.average.values[1]).toBeCloseTo(0, 12);
        expect(p.average.values[2]).toBeCloseTo(5, 12);
        expect(p.coefficients.values[0]).toBeCloseTo(3, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(-2, 10);
        expect(p.coefficients.values[2]).toBe(-1);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 18);
        }
    });

    it('recovers a horizontal plane', () => {
        const points = [v3(0, 0, 4), v3(1, 0, 4), v3(0, 1, 4), v3(3, -2, 4)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.coefficients.values[0]).toBeCloseTo(0, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(0, 10);
        expect(p.average.values[2]).toBeCloseTo(4, 12);
    });

    it('is exact for three noncollinear points', () => {
        const points = [v3(0, 0, 1), v3(1, 0, 4), v3(0, 1, -1)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(true);
        const p = fitter.getParameters();
        expect(p.coefficients.values[0]).toBeCloseTo(3, 10);
        expect(p.coefficients.values[1]).toBeCloseTo(-2, 10);
        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 18);
        }
    });

    it('produces least-squares residuals that satisfy the normal equations', () => {
        const random = makeRandom(777);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            const n = 25;
            for (let i = 0; i < n; ++i) {
                const x = 4 * random() - 2;
                const y = 4 * random() - 2;
                points.push(v3(x, y, 1.5 * x + 0.25 * y - 3
                    + 0.5 * (2 * random() - 1)));
            }
            const fitter = new ApprHeightPlane3();
            expect(fitter.fit(points)).toBe(true);
            const p = fitter.getParameters();
            const a = p.coefficients.values[0], b = p.coefficients.values[1];

            // The residual r = a*dx + b*dy - dz must be orthogonal to dx and
            // to dy (the normal equations of the least-squares problem).
            let dotX = 0, dotY = 0, sumR = 0;
            for (const q of points) {
                const dx = q.values[0] - p.average.values[0];
                const dy = q.values[1] - p.average.values[1];
                const dz = q.values[2] - p.average.values[2];
                const r = a * dx + b * dy - dz;
                dotX += r * dx;
                dotY += r * dy;
                sumR += r;
            }
            expect(dotX).toBeCloseTo(0, 8);
            expect(dotY).toBeCloseTo(0, 8);
            expect(sumR).toBeCloseTo(0, 8);
        }
    });

    it('computes the squared vertical error', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 3), v3(0, 1, -2)];
        const fitter = new ApprHeightPlane3();
        fitter.fit(points);
        // The fitted plane is z = 3x - 2y; at (1,1,5) the deviation is 4.
        expect(fitter.error(v3(1, 1, 5))).toBeCloseTo(16, 8);
    });

    it('fails for collinear (x,y) samples', () => {
        const points = [v3(0, 0, 1), v3(1, 1, 2), v3(2, 2, 3), v3(3, 3, 9)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().average.values).toEqual([0, 0, 0]);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0, 0]);
    });

    it('fails for coincident points', () => {
        const points = [v3(1, 2, 3), v3(1, 2, 3), v3(1, 2, 3)];
        const fitter = new ApprHeightPlane3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().coefficients.values).toEqual([0, 0, 0]);
    });

    it('deep-copies the parameters', () => {
        const source = new ApprHeightPlane3();
        source.fit([v3(0, 0, 1), v3(1, 0, 4), v3(0, 1, -1)]);
        const target = new ApprHeightPlane3();
        target.copyParameters(source);
        expect(target.getParameters().coefficients.values[0]).toBeCloseTo(3, 10);
        source.getParameters().coefficients.values[0] = 99;
        expect(target.getParameters().coefficients.values[0]).toBeCloseTo(3, 10);
    });
});
