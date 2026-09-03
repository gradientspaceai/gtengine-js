import { describe, it, expect } from 'vitest';
import { ApprGreatArc3, ApprGreatCircle3 } from '../src/ApprGreatCircle3.js';
import { Vector, dot, length, normalize } from '../src/Vector.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = v3(x, y, z);
    normalize(v);
    return v;
}

// A right-handed orthonormal basis {N, U, V}.
function basisFromNormal(N: Vector): { U: Vector, V: Vector } {
    const n = N.values;
    const U = Math.abs(n[0]) > Math.abs(n[1]) ?
        unit(-n[2], 0, n[0]) : unit(0, n[2], -n[1]);
    const u = U.values;
    const V = v3(n[1] * u[2] - n[2] * u[1], n[2] * u[0] - n[0] * u[2],
        n[0] * u[1] - n[1] * u[0]);
    return { U, V };
}

// Unit-length points on the great circle with normal N, at the given angles.
function circlePoints(N: Vector, angles: readonly number[]): Vector[] {
    const { U, V } = basisFromNormal(N);
    return angles.map(t => {
        const c = Math.cos(t), s = Math.sin(t);
        return v3(c * U.values[0] + s * V.values[0],
            c * U.values[1] + s * V.values[1],
            c * U.values[2] + s * V.values[2]);
    });
}

function uniformAngles(count: number, span: number = 2 * Math.PI,
    start: number = 0): number[] {
    const angles: number[] = [];
    for (let i = 0; i < count; ++i) {
        angles.push(start + span * i / count);
    }
    return angles;
}

// |Dot(a,b)| is 1 when a and b are parallel up to sign.
function parallelUpToSign(a: Vector, b: Vector): number {
    return Math.abs(dot(a, b));
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

describe('ApprGreatCircle3', () => {
    it('recovers the normal of a great circle up to sign', () => {
        const N = unit(1, -2, 3);
        const points = circlePoints(N, uniformAngles(64));
        const fit = new ApprGreatCircle3().compute(points);
        expect(length(fit)).toBeCloseTo(1, 12);
        expect(parallelUpToSign(fit, N)).toBeCloseTo(1, 12);
        for (const p of points) {
            expect(dot(fit, p)).toBeCloseTo(0, 12);
        }
    });

    it('recovers the normal for the canonical xy great circle', () => {
        const points = circlePoints(v3(0, 0, 1), uniformAngles(37));
        const fit = new ApprGreatCircle3().compute(points);
        expect(Math.abs(fit.values[2])).toBeCloseTo(1, 12);
        expect(fit.values[0]).toBeCloseTo(0, 12);
        expect(fit.values[1]).toBeCloseTo(0, 12);
    });

    it('is nearly exact for noisy samples of a great circle', () => {
        const rnd = makeRandom(20250901);
        const N = unit(-0.3, 0.8, 0.5);
        const points = circlePoints(N, uniformAngles(256)).map(p => {
            const q = v3(p.values[0] + 1e-3 * (2 * rnd() - 1),
                p.values[1] + 1e-3 * (2 * rnd() - 1),
                p.values[2] + 1e-3 * (2 * rnd() - 1));
            normalize(q);
            return q;
        });
        const fit = new ApprGreatCircle3().compute(points);
        expect(length(fit)).toBeCloseTo(1, 12);
        expect(parallelUpToSign(fit, N)).toBeGreaterThan(1 - 1e-6);
    });

    it('returns a unit vector orthogonal to the input when all points are equal', () => {
        // The covariance matrix has rank 1, so the eigenvector for the
        // smallest (zero) eigenvalue is any unit vector orthogonal to the
        // repeated point.
        const p = unit(1, 1, 0);
        const points = new Array<Vector>(10).fill(p);
        const fit = new ApprGreatCircle3().compute(points);
        expect(length(fit)).toBeCloseTo(1, 12);
        expect(dot(fit, p)).toBeCloseTo(0, 12);
    });

    it('rejects an empty point set', () => {
        expect(() => new ApprGreatCircle3().compute([])).toThrow(/no points/);
    });

    it('rejects non-3D points', () => {
        expect(() => new ApprGreatCircle3().compute(
            [Vector.fromArray([1, 0])])).toThrow(/3D/);
    });

    it('recovers random great-circle normals (aggregated)', () => {
        const rnd = makeRandom(97531);
        let worst = 1;
        for (let trial = 0; trial < 50; ++trial) {
            const N = unit(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const angles = uniformAngles(23, 2 * Math.PI, rnd());
            const fit = new ApprGreatCircle3().compute(circlePoints(N, angles));
            worst = Math.min(worst, parallelUpToSign(fit, N));
        }
        expect(worst).toBeGreaterThan(1 - 1e-10);
    });
});

describe('ApprGreatArc3', () => {
    it('recovers the arc endpoints of a quarter circle in the xy-plane', () => {
        // The points span angles [0, pi/2] measured in the xy-plane. The
        // returned normal is +/- (0,0,1) and the arc endpoints are the two
        // extreme points of the sample set.
        const angles = uniformAngles(9, Math.PI / 2);
        const points = circlePoints(v3(0, 0, 1), angles);
        const { normal, arcEnd0, arcEnd1 } = new ApprGreatArc3().compute(points);
        expect(Math.abs(normal.values[2])).toBeCloseTo(1, 12);
        expect(length(arcEnd0)).toBeCloseTo(1, 12);
        expect(length(arcEnd1)).toBeCloseTo(1, 12);
        expect(dot(normal, arcEnd0)).toBeCloseTo(0, 12);
        expect(dot(normal, arcEnd1)).toBeCloseTo(0, 12);

        // The endpoints are the extreme sample points (up to which one is
        // the start, which depends on the sign of the fitted normal).
        const first = points[0];
        const last = points[points.length - 1];
        const matchesDirect = parallelUpToSign(arcEnd0, first) > 1 - 1e-10 &&
            parallelUpToSign(arcEnd1, last) > 1 - 1e-10;
        const matchesSwapped = parallelUpToSign(arcEnd0, last) > 1 - 1e-10 &&
            parallelUpToSign(arcEnd1, first) > 1 - 1e-10;
        expect(matchesDirect || matchesSwapped).toBe(true);

        // Every sample lies inside the sector, so the angle from arcEnd0 to
        // any sample (measured counterclockwise about the normal) is no
        // larger than the angle from arcEnd0 to arcEnd1.
        const sector = Math.acos(Math.min(1, Math.max(-1,
            dot(arcEnd0, arcEnd1))));
        for (const p of points) {
            const angle = Math.acos(Math.min(1, Math.max(-1, dot(arcEnd0, p))));
            expect(angle).toBeLessThanOrEqual(sector + 1e-12);
        }
    });

    it('produces endpoints perpendicular to the fitted normal for a tilted arc', () => {
        const N = unit(1, 2, -1);
        const points = circlePoints(N, uniformAngles(16, 1.1, 0.3));
        const { normal, arcEnd0, arcEnd1 } = new ApprGreatArc3().compute(points);
        expect(parallelUpToSign(normal, N)).toBeCloseTo(1, 10);
        expect(dot(normal, arcEnd0)).toBeCloseTo(0, 10);
        expect(dot(normal, arcEnd1)).toBeCloseTo(0, 10);
        expect(length(arcEnd0)).toBeCloseTo(1, 12);
        expect(length(arcEnd1)).toBeCloseTo(1, 12);

        // The sector angle matches the angular span of the samples.
        const sector = Math.acos(Math.min(1, Math.max(-1,
            dot(arcEnd0, arcEnd1))));
        expect(sector).toBeCloseTo(1.1 * 15 / 16, 8);
    });

    it('rejects an empty point set', () => {
        expect(() => new ApprGreatArc3().compute([])).toThrow(/no points/);
    });
});
