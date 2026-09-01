import { describe, it, expect } from 'vitest';
import { ApprEllipse2 } from '../src/ApprEllipse2';
import { Hyperellipsoid } from '../src/Hyperellipsoid';
import { Vector, dot } from '../src/Vector';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

// Points that lie exactly on the ellipse with the given center, orientation
// angle and semi-axis lengths (a along the rotated x-axis, b along the
// rotated y-axis).
function ellipsePoints(cx: number, cy: number, a: number, b: number,
    angle: number, count: number, phase: number = 0): Vector[] {
    const ca = Math.cos(angle), sa = Math.sin(angle);
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        const t = phase + 2 * Math.PI * i / count;
        const x = a * Math.cos(t), y = b * Math.sin(t);
        points.push(v2(cx + ca * x - sa * y, cy + sa * x + ca * y));
    }
    return points;
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

// |Dot(u,v)| is 1 when the unit-length u and v are parallel up to sign.
function parallelUpToSign(u: Vector, v: Vector): number {
    return Math.abs(dot(u, v));
}

describe('ApprEllipse2', () => {
    it('recovers an axis-aligned ellipse from exact points', () => {
        const points = ellipsePoints(0, 0, 4, 2, 0, 64);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 256, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(0, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(0, 8);
        // The eigenvalues of M are sorted increasing, so extent[0] is the
        // major semi-axis length and extent[1] the minor one.
        expect(ellipse.extent.values[0]).toBeCloseTo(4, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 8);
        expect(parallelUpToSign(ellipse.axis[0], v2(1, 0))).toBeCloseTo(1, 10);
        expect(parallelUpToSign(ellipse.axis[1], v2(0, 1))).toBeCloseTo(1, 10);
    });

    it('recovers a translated and rotated ellipse from exact points', () => {
        const angle = 0.4;
        const points = ellipsePoints(2, -1, 3, 1.5, angle, 64);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 256, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(2, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(-1, 8);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(1.5, 8);

        const major = v2(Math.cos(angle), Math.sin(angle));
        const minor = v2(-Math.sin(angle), Math.cos(angle));
        expect(parallelUpToSign(ellipse.axis[0], major)).toBeCloseTo(1, 10);
        expect(parallelUpToSign(ellipse.axis[1], minor)).toBeCloseTo(1, 10);

        // The unit-length axes are orthogonal.
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
        expect(dot(ellipse.axis[0], ellipse.axis[0])).toBeCloseTo(1, 12);

        // The implicit equation is satisfied by the input points.
        for (const p of points) {
            const d = [p.values[0] - ellipse.center.values[0],
                p.values[1] - ellipse.center.values[1]];
            let q = 0;
            for (let i = 0; i < 2; ++i) {
                const c = d[0] * ellipse.axis[i].values[0] +
                    d[1] * ellipse.axis[i].values[1];
                q += (c * c) / (ellipse.extent.values[i] *
                    ellipse.extent.values[i]);
            }
            expect(q).toBeCloseTo(1, 6);
        }
    });

    it('handles the circle special case', () => {
        const points = ellipsePoints(-2, 5, 3, 3, 0, 48);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 128, false, ellipse);

        expect(error).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(-2, 8);
        expect(ellipse.center.values[1]).toBeCloseTo(5, 8);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 8);
        expect(ellipse.extent.values[1]).toBeCloseTo(3, 8);
        // The axes are still orthonormal even though the axis directions are
        // not unique for a circle.
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
    });

    it('continues the minimization when the ellipse is the initial guess', () => {
        const points = ellipsePoints(1, 1, 5, 2, -0.7, 96);
        const fitter = new ApprEllipse2();
        const ellipse = new Hyperellipsoid(2);
        const error0 = fitter.compute(points, 1, false, ellipse);
        const error1 = fitter.compute(points, 3, true, ellipse);
        const error2 = fitter.compute(points, 200, true, ellipse);

        // Restarting from the previous result continues the descent, so the
        // error never increases (up to floating-point round-off once the
        // minimizer has converged).
        expect(error1).toBeLessThanOrEqual(error0 + 1e-20);
        expect(error2).toBeLessThanOrEqual(error1 + 1e-20);
        expect(error2).toBeLessThan(1e-16);
        expect(ellipse.center.values[0]).toBeCloseTo(1, 6);
        expect(ellipse.center.values[1]).toBeCloseTo(1, 6);
        expect(ellipse.extent.values[0]).toBeCloseTo(5, 6);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 6);
    });

    it('reports the initial-guess error when no iterations are requested', () => {
        const points = ellipsePoints(0, 0, 4, 2, 0, 40);
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(points, 0, false, ellipse);

        // The initial guess is the oriented bounding box of the points, so
        // the error is nonzero but the returned ellipse is well formed.
        expect(error).toBeGreaterThan(0);
        expect(Number.isFinite(error)).toBe(true);
        expect(ellipse.extent.values[0]).toBeGreaterThan(0);
        expect(ellipse.extent.values[1]).toBeGreaterThan(0);
        expect(dot(ellipse.axis[0], ellipse.axis[1])).toBeCloseTo(0, 12);
    });

    it('is robust to small noise', () => {
        const rnd = makeRandom(1234567);
        const exact = ellipsePoints(0.5, -0.25, 3, 2, 0.9, 200);
        const noisy = exact.map(p => v2(
            p.values[0] + 1e-3 * (2 * rnd() - 1),
            p.values[1] + 1e-3 * (2 * rnd() - 1)));
        const ellipse = new Hyperellipsoid(2);
        const error = new ApprEllipse2().compute(noisy, 256, false, ellipse);

        expect(error).toBeLessThan(1e-4);
        expect(ellipse.center.values[0]).toBeCloseTo(0.5, 2);
        expect(ellipse.center.values[1]).toBeCloseTo(-0.25, 2);
        expect(ellipse.extent.values[0]).toBeCloseTo(3, 2);
        expect(ellipse.extent.values[1]).toBeCloseTo(2, 2);
    });

    it('recovers random ellipses (aggregated)', () => {
        const rnd = makeRandom(24680);
        let worstExtentError = 0;
        let worstCenterError = 0;
        for (let trial = 0; trial < 12; ++trial) {
            const cx = 4 * rnd() - 2, cy = 4 * rnd() - 2;
            const a = 1 + 3 * rnd();
            const b = 0.5 + 0.5 * a * rnd();
            const angle = Math.PI * rnd();
            const points = ellipsePoints(cx, cy, a, b, angle, 60, rnd());
            const ellipse = new Hyperellipsoid(2);
            new ApprEllipse2().compute(points, 512, false, ellipse);
            const major = Math.max(a, b), minor = Math.min(a, b);
            worstExtentError = Math.max(worstExtentError,
                Math.abs(ellipse.extent.values[0] - major),
                Math.abs(ellipse.extent.values[1] - minor));
            worstCenterError = Math.max(worstCenterError,
                Math.abs(ellipse.center.values[0] - cx),
                Math.abs(ellipse.center.values[1] - cy));
        }
        expect(worstExtentError).toBeLessThan(1e-5);
        expect(worstCenterError).toBeLessThan(1e-5);
    });

    it('rejects an empty point set', () => {
        const ellipse = new Hyperellipsoid(2);
        expect(() => new ApprEllipse2().compute([], 4, false, ellipse))
            .toThrow(/no points/);
    });

    it('rejects points that are not 2D', () => {
        const ellipse = new Hyperellipsoid(2);
        expect(() => new ApprEllipse2().compute(
            [Vector.fromArray([1, 2, 3])], 4, false, ellipse)).toThrow(/2D/);
    });
});
