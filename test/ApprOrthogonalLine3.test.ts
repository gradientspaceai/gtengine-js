import { describe, it, expect } from 'vitest';
import { ApprOrthogonalLine3 } from '../src/ApprOrthogonalLine3';
import { Vector, dot } from '../src/Vector';

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

describe('ApprOrthogonalLine3', () => {
    it('initializes the parameters to zero', () => {
        const fitter = new ApprOrthogonalLine3();
        const line = fitter.getParameters();
        expect(line.origin.values).toEqual([0, 0, 0]);
        expect(line.direction.values).toEqual([0, 0, 0]);
        expect(fitter.getMinimumRequired()).toBe(2);
    });

    it('recovers a line from points that lie exactly on it', () => {
        const unit = 1 / Math.sqrt(3);
        const d = v3(unit, unit, unit);
        const points: Vector[] = [];
        for (const t of [-3, -1, 0, 1, 3]) {
            points.push(v3(1 + t * d.values[0], -2 + t * d.values[1],
                4 + t * d.values[2]));
        }
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(line.origin.values[0]).toBeCloseTo(1, 12);
        expect(line.origin.values[1]).toBeCloseTo(-2, 12);
        expect(line.origin.values[2]).toBeCloseTo(4, 12);
        expect(Math.abs(dot(line.direction, d))).toBeCloseTo(1, 10);
        expect(dot(line.direction, line.direction)).toBeCloseTo(1, 12);

        for (const point of points) {
            expect(fitter.error(point)).toBeCloseTo(0, 14);
        }
    });

    it('recovers an axis-aligned line', () => {
        const points = [v3(0, 3, -1), v3(1, 3, -1), v3(5, 3, -1), v3(-4, 3, -1)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        const line = fitter.getParameters();
        expect(Math.abs(line.direction.values[0])).toBeCloseTo(1, 10);
        expect(line.origin.values[1]).toBeCloseTo(3, 12);
        expect(line.origin.values[2]).toBeCloseTo(-1, 12);
    });

    it('computes the squared orthogonal distance as the error', () => {
        const points = [v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0), v3(3, 0, 0)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(true);
        // The fitted line is the x-axis through (1.5,0,0).
        expect(fitter.error(v3(0, 3, 4))).toBeCloseTo(25, 8);
        expect(fitter.error(v3(-100, 0, 0))).toBeCloseTo(0, 8);
    });

    it('is not unique for an isotropic point set', () => {
        // The vertices of a regular octahedron have covariance a multiple of
        // the identity, so the maximum eigenvalue has multiplicity 3.
        const points = [v3(1, 0, 0), v3(-1, 0, 0), v3(0, 1, 0),
            v3(0, -1, 0), v3(0, 0, 1), v3(0, 0, -1)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
        expect(dot(fitter.getParameters().direction,
            fitter.getParameters().direction)).toBeCloseTo(1, 10);
    });

    it('is not unique for coincident points', () => {
        const points = [v3(2, 2, 2), v3(2, 2, 2), v3(2, 2, 2)];
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
        expect(fitter.getParameters().origin.values[0]).toBeCloseTo(2, 12);
    });

    it('is not unique for a planar isotropic point set', () => {
        // Points on a circle in the z = 0 plane: the two largest eigenvalues
        // are equal.
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v3(Math.cos(angle), Math.sin(angle), 0));
        }
        const fitter = new ApprOrthogonalLine3();
        expect(fitter.fit(points)).toBe(false);
    });

    it('minimizes the sum of squared orthogonal distances', () => {
        const random = makeRandom(5150);
        for (let trial = 0; trial < 10; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 25; ++i) {
                const t = 8 * random() - 4;
                points.push(v3(
                    t + 0.15 * (2 * random() - 1),
                    2 * t + 1 + 0.15 * (2 * random() - 1),
                    -t + 0.15 * (2 * random() - 1)));
            }
            const fitter = new ApprOrthogonalLine3();
            expect(fitter.fit(points)).toBe(true);
            const line = fitter.getParameters();

            let best = 0;
            for (const p of points) {
                best += fitter.error(p);
            }

            // Perturbing the direction must not reduce the objective.
            for (const eps of [-0.05, -0.01, 0.01, 0.05]) {
                for (let axis = 0; axis < 3; ++axis) {
                    const u = line.direction.clone();
                    u.values[axis] += eps;
                    const len = Math.sqrt(dot(u, u));
                    for (let d = 0; d < 3; ++d) {
                        u.values[d] /= len;
                    }
                    let sum = 0;
                    for (const p of points) {
                        let sqrlen = 0, proj = 0;
                        for (let d = 0; d < 3; ++d) {
                            const diff = p.values[d] - line.origin.values[d];
                            sqrlen += diff * diff;
                            proj += diff * u.values[d];
                        }
                        sum += sqrlen - proj * proj;
                    }
                    expect(sum).toBeGreaterThanOrEqual(best - 1e-10);
                }
            }
        }
    });

    it('deep-copies the parameters', () => {
        const source = new ApprOrthogonalLine3();
        source.fit([v3(0, 0, 0), v3(1, 1, 1), v3(2, 2, 2)]);
        const target = new ApprOrthogonalLine3();
        target.copyParameters(source);
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
        source.getParameters().origin.values[0] = 33;
        expect(target.getParameters().origin.values[0]).toBeCloseTo(1, 12);
    });
});
