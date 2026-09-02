import { describe, it, expect } from 'vitest';
import { ApprEllipsoid3 } from '../src/ApprEllipsoid3';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import { Vector, dot } from '../src/Vector';
import { cross } from '../src/Vector3';

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

// A right-handed orthonormal frame obtained by rotating the standard frame by
// 'yaw' about the z-axis and then by 'pitch' about the rotated y-axis.
function frame(yaw: number, pitch: number): Vector[] {
    const cy = Math.cos(yaw), sy = Math.sin(yaw);
    const cp = Math.cos(pitch), sp = Math.sin(pitch);
    const u0 = v3(cy * cp, sy * cp, -sp);
    const u1 = v3(-sy, cy, 0);
    const u2 = cross(u0, u1);
    return [u0, u1, u2];
}

// Points that lie exactly on the ellipsoid with the given center, frame and
// semi-axis lengths. The sampling uses a spiral so that the points cover the
// whole surface for any 'count'.
function ellipsoidPoints(center: Vector, axis: Vector[], extent: number[],
    count: number, phase: number = 0): Vector[] {
    const points: Vector[] = [];
    const golden = Math.PI * (3 - Math.sqrt(5));
    for (let i = 0; i < count; ++i) {
        const w = 1 - 2 * (i + 0.5) / count;
        const r = Math.sqrt(Math.max(0, 1 - w * w));
        const theta = phase + golden * i;
        const local = [extent[0] * r * Math.cos(theta),
            extent[1] * r * Math.sin(theta), extent[2] * w];
        const p = v3(center.values[0], center.values[1], center.values[2]);
        for (let j = 0; j < 3; ++j) {
            for (let k = 0; k < 3; ++k) {
                p.values[k] += local[j] * axis[j].values[k];
            }
        }
        points.push(p);
    }
    return points;
}

// The implicit ellipsoid value sum_i (Dot(X-C,U[i])/e[i])^2, which is 1 on
// the ellipsoid.
function implicit(ellipsoid: Hyperellipsoid, p: Vector): number {
    let q = 0;
    for (let i = 0; i < 3; ++i) {
        let c = 0;
        for (let k = 0; k < 3; ++k) {
            c += (p.values[k] - ellipsoid.center.values[k]) *
                ellipsoid.axis[i].values[k];
        }
        q += (c * c) / (ellipsoid.extent.values[i] *
            ellipsoid.extent.values[i]);
    }
    return q;
}

function parallelUpToSign(u: Vector, v: Vector): number {
    return Math.abs(dot(u, v));
}

// The extents of the fitted ellipsoid come out sorted decreasing, because the
// eigenvalues 1/e^2 of M are sorted increasing.
function sortedDecreasing(extent: number[]): number[] {
    return [...extent].sort((x, y) => y - x);
}

describe('ApprEllipsoid3', () => {
    it('recovers an axis-aligned ellipsoid from exact points', () => {
        const axis = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const points = ellipsoidPoints(v3(0, 0, 0), axis, [4, 2, 1], 200);
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(points, 512, false,
            ellipsoid);

        // The 2-step gradient descent converges linearly, so after 512
        // iterations the residual is small but not at round-off level.
        expect(error).toBeLessThan(1e-8);
        for (let k = 0; k < 3; ++k) {
            expect(ellipsoid.center.values[k]).toBeCloseTo(0, 4);
        }
        expect(Math.abs(ellipsoid.extent.values[0] - 4)).toBeLessThan(1e-3);
        expect(Math.abs(ellipsoid.extent.values[1] - 2)).toBeLessThan(1e-3);
        expect(Math.abs(ellipsoid.extent.values[2] - 1)).toBeLessThan(1e-3);
        expect(parallelUpToSign(ellipsoid.axis[0], v3(1, 0, 0)))
            .toBeCloseTo(1, 6);
        expect(parallelUpToSign(ellipsoid.axis[1], v3(0, 1, 0)))
            .toBeCloseTo(1, 6);
        expect(parallelUpToSign(ellipsoid.axis[2], v3(0, 0, 1)))
            .toBeCloseTo(1, 6);
    });

    it('recovers a translated and rotated ellipsoid from exact points', () => {
        const axis = frame(0.5, -0.3);
        const center = v3(2, -1, 0.5);
        const extent = [3, 2, 1.25];
        const points = ellipsoidPoints(center, axis, extent, 240);
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(points, 512, false,
            ellipsoid);

        expect(error).toBeLessThan(1e-14);
        for (let k = 0; k < 3; ++k) {
            expect(ellipsoid.center.values[k])
                .toBeCloseTo(center.values[k], 5);
        }
        for (let i = 0; i < 3; ++i) {
            expect(ellipsoid.extent.values[i]).toBeCloseTo(extent[i], 5);
            expect(parallelUpToSign(ellipsoid.axis[i], axis[i]))
                .toBeCloseTo(1, 5);
        }

        // The axes are orthonormal.
        for (let i = 0; i < 3; ++i) {
            expect(dot(ellipsoid.axis[i], ellipsoid.axis[i]))
                .toBeCloseTo(1, 10);
            for (let j = i + 1; j < 3; ++j) {
                expect(dot(ellipsoid.axis[i], ellipsoid.axis[j]))
                    .toBeCloseTo(0, 10);
            }
        }

        // The implicit equation is satisfied by the input points.
        for (const p of points) {
            expect(implicit(ellipsoid, p)).toBeCloseTo(1, 5);
        }
    });

    it('handles the sphere special case', () => {
        const axis = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const points = ellipsoidPoints(v3(-2, 5, 1), axis, [3, 3, 3], 150);
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(points, 256, false,
            ellipsoid);

        expect(error).toBeLessThan(1e-14);
        expect(ellipsoid.center.values[0]).toBeCloseTo(-2, 5);
        expect(ellipsoid.center.values[1]).toBeCloseTo(5, 5);
        expect(ellipsoid.center.values[2]).toBeCloseTo(1, 5);
        for (let i = 0; i < 3; ++i) {
            expect(ellipsoid.extent.values[i]).toBeCloseTo(3, 5);
        }
        // The axis directions are not unique for a sphere, but they are still
        // orthonormal.
        expect(dot(ellipsoid.axis[0], ellipsoid.axis[1])).toBeCloseTo(0, 10);
        expect(dot(ellipsoid.axis[0], ellipsoid.axis[2])).toBeCloseTo(0, 10);
    });

    it('continues the minimization when the ellipsoid is the initial guess', () => {
        const axis = frame(-0.7, 0.4);
        const points = ellipsoidPoints(v3(1, 1, -1), axis, [5, 2, 1.5], 180);
        const fitter = new ApprEllipsoid3();
        const ellipsoid = new Hyperellipsoid(3);
        const error0 = fitter.compute(points, 1, false, ellipsoid);
        const error1 = fitter.compute(points, 3, true, ellipsoid);
        const error2 = fitter.compute(points, 600, true, ellipsoid);

        // Restarting from the previous result continues the descent, so the
        // error never increases (up to round-off once converged).
        expect(error1).toBeLessThanOrEqual(error0 + 1e-20);
        expect(error2).toBeLessThanOrEqual(error1 + 1e-20);
        expect(error2).toBeLessThan(1e-12);
        expect(ellipsoid.center.values[0]).toBeCloseTo(1, 4);
        expect(ellipsoid.center.values[1]).toBeCloseTo(1, 4);
        expect(ellipsoid.center.values[2]).toBeCloseTo(-1, 4);
        expect(sortedDecreasing(ellipsoid.extent.values)[0])
            .toBeCloseTo(5, 4);
        expect(sortedDecreasing(ellipsoid.extent.values)[2])
            .toBeCloseTo(1.5, 4);
    });

    it('reports the initial-guess error when no iterations are requested', () => {
        const axis = [v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1)];
        const points = ellipsoidPoints(v3(0, 0, 0), axis, [4, 2, 1], 120);
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(points, 0, false,
            ellipsoid);

        // The initial guess is the oriented bounding box of the points, so
        // the error is nonzero but the returned ellipsoid is well formed.
        expect(error).toBeGreaterThan(0);
        expect(Number.isFinite(error)).toBe(true);
        for (let i = 0; i < 3; ++i) {
            expect(ellipsoid.extent.values[i]).toBeGreaterThan(0);
        }
        expect(dot(ellipsoid.axis[0], ellipsoid.axis[1])).toBeCloseTo(0, 10);

        // Iterating strictly improves on the initial guess.
        const refined = new Hyperellipsoid(3);
        const refinedError = new ApprEllipsoid3().compute(points, 512, false,
            refined);
        expect(refinedError).toBeLessThan(error);
    });

    it('is robust to small noise', () => {
        const rnd = makeRandom(1234567);
        const axis = frame(0.9, 0.2);
        const center = v3(0.5, -0.25, 0.75);
        const exact = ellipsoidPoints(center, axis, [3, 2, 1.5], 400);
        const noisy = exact.map(p => v3(
            p.values[0] + 1e-3 * (2 * rnd() - 1),
            p.values[1] + 1e-3 * (2 * rnd() - 1),
            p.values[2] + 1e-3 * (2 * rnd() - 1)));
        const ellipsoid = new Hyperellipsoid(3);
        const error = new ApprEllipsoid3().compute(noisy, 512, false,
            ellipsoid);

        expect(error).toBeLessThan(1e-4);
        for (let k = 0; k < 3; ++k) {
            expect(Math.abs(ellipsoid.center.values[k] - center.values[k]))
                .toBeLessThan(1e-2);
        }
        const sorted = sortedDecreasing(ellipsoid.extent.values);
        expect(Math.abs(sorted[0] - 3)).toBeLessThan(1e-2);
        expect(Math.abs(sorted[1] - 2)).toBeLessThan(1e-2);
        expect(Math.abs(sorted[2] - 1.5)).toBeLessThan(1e-2);
    });

    it('recovers random ellipsoids (aggregated)', () => {
        const rnd = makeRandom(24680);
        let worstExtentError = 0;
        let worstCenterError = 0;
        for (let trial = 0; trial < 6; ++trial) {
            const center = v3(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const e0 = 2 + 2 * rnd();
            const e1 = 1 + rnd();
            const e2 = 0.5 + 0.5 * rnd();
            const axis = frame(Math.PI * rnd(), Math.PI * rnd() - 0.5 * Math.PI);
            const points = ellipsoidPoints(center, axis, [e0, e1, e2], 150,
                rnd());
            const ellipsoid = new Hyperellipsoid(3);
            new ApprEllipsoid3().compute(points, 800, false, ellipsoid);

            const expected = sortedDecreasing([e0, e1, e2]);
            const actual = sortedDecreasing(ellipsoid.extent.values);
            for (let i = 0; i < 3; ++i) {
                worstExtentError = Math.max(worstExtentError,
                    Math.abs(actual[i] - expected[i]));
            }
            for (let k = 0; k < 3; ++k) {
                worstCenterError = Math.max(worstCenterError,
                    Math.abs(ellipsoid.center.values[k] - center.values[k]));
            }
        }
        expect(worstExtentError).toBeLessThan(5e-3);
        expect(worstCenterError).toBeLessThan(5e-3);
    });

    it('rejects an empty point set', () => {
        const ellipsoid = new Hyperellipsoid(3);
        expect(() => new ApprEllipsoid3().compute([], 4, false, ellipsoid))
            .toThrow(/no points/);
    });

    it('rejects points that are not 3D', () => {
        const ellipsoid = new Hyperellipsoid(3);
        expect(() => new ApprEllipsoid3().compute(
            [Vector.fromArray([1, 2])], 4, false, ellipsoid)).toThrow(/3D/);
    });
});
