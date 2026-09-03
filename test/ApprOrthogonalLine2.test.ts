import { describe, it, expect } from 'vitest';
import { ApprOrthogonalLine2 } from '../src/ApprOrthogonalLine2.js';
import { Vector, dot } from '../src/Vector.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprOrthogonalLine2', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprOrthogonalLine2();
        const line = fitter.getParameters();
        expect(line.origin.values).toEqual([0, 0]);
        expect(line.direction.values).toEqual([0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers a line from points that lie exactly on it', () => {
        // The line through (1,2) with unit direction (3,4)/5.
        const d = v2(0.6, 0.8);
        const points: Vector[] = [];
        for (const t of [-2, -1, 0, 1, 2]) {
            points.push(v2(1 + t * d.values[0], 2 + t * d.values[1]));
        }
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(1, 12);
        expect(line.origin.values[1]).toBeCloseTo(2, 12);
        expect(Math.abs(dot(line.direction, d))).toBeCloseTo(1, 12);
        expect(dot(line.direction, line.direction)).toBeCloseTo(1, 12);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 14);
        }
    });

    it('recovers a vertical line, which the height fitter cannot fit', () => {
        const points = [v2(4, -3), v2(4, 0), v2(4, 1), v2(4, 6)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(4, 12);
        expect(Math.abs(line.direction.values[1])).toBeCloseTo(1, 12);
        expect(line.direction.values[0]).toBeCloseTo(0, 12);
    });

    it('computes the squared orthogonal distance as the error', () => {
        const points = [v2(0, 0), v2(1, 0), v2(2, 0), v2(3, 0)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(true);
        // The fitted line is the x-axis through (1.5,0).
        expect(fitter.error(v2(10, 3))).toBeCloseTo(9, 10);
        expect(fitter.error(v2(-10, 0))).toBeCloseTo(0, 10);
    });

    it('is not unique for a symmetric point set', () => {
        // The covariance matrix is a multiple of the identity, so the
        // maximum eigenvalue has multiplicity 2.
        const points = [v2(1, 0), v2(-1, 0), v2(0, 1), v2(0, -1)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(false);
        // The line parameters are still assigned (the fit always succeeds).
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(0, 12);
        expect(dot(fitter.getParameters().direction,
            fitter.getParameters().direction)).toBeCloseTo(1, 12);
    });

    it('is not unique for coincident points', () => {
        const points = [v2(5, 5), v2(5, 5), v2(5, 5)];
        const fitter = new ApprOrthogonalLine2();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(5, 12);
    });

    it('minimizes the sum of squared orthogonal distances', () => {
        const random = makeRandom(9091);
        for (let trial = 0; trial < 15; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 20; ++i) {
                const t = 6 * random() - 3;
                points.push(v2(t + 0.2 * (2 * random() - 1),
                    2 * t - 1 + 0.2 * (2 * random() - 1)));
            }
            const fitter = new ApprOrthogonalLine2();
            expect(fitter.fit(points)).toBe(true);
            const line = fitter.getParameters();

            let best = 0;
            for (const p of points) {
                best += fitter.error(p);
            }

            // Rotating the fitted direction about the fitted origin must not
            // reduce the sum of squared orthogonal distances.
            for (const delta of [-0.2, -0.05, -0.01, 0.01, 0.05, 0.2]) {
                const angle = Math.atan2(line.direction.values[1],
                    line.direction.values[0]) + delta;
                const u = v2(Math.cos(angle), Math.sin(angle));
                let sum = 0;
                for (const p of points) {
                    const dx = p.values[0] - line.origin.values[0];
                    const dy = p.values[1] - line.origin.values[1];
                    const proj = dx * u.values[0] + dy * u.values[1];
                    sum += dx * dx + dy * dy - proj * proj;
                }
                expect(sum).toBeGreaterThanOrEqual(best - 1e-12);
            }
        }
    });

    it('deep-copies the parameters', () => {
        const source = new ApprOrthogonalLine2();
        source.fit([v2(0, 0), v2(1, 1), v2(2, 2)]);
        const target = new ApprOrthogonalLine2();
        target.copyParameters(source);
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
        source.getParameters().origin.values[0] = 55;
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
    });
});
