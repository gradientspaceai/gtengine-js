import { describe, expect, it } from 'vitest';
import {
    ApprCone3EllipseAndPoints,
    ApprCone3EllipseAndPointsControl
} from '../src/ApprCone3EllipseAndPoints';
import { Ellipse3 } from '../src/Ellipse3';
import { GTE_C_HALF_PI } from '../src/Constants';
import { Vector, add, dot, mul, normalize, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

const v3 = (x: number, y: number, z: number): Vector =>
    Vector.fromArray([x, y, z]);

const unit = (x: number, y: number, z: number): Vector => {
    const v = v3(x, y, z);
    normalize(v);
    return v;
};

// An orthonormal basis {E0,E1,D} with D the third vector.
function basisFor(D: Vector): { E0: Vector; E1: Vector } {
    const helper = Math.abs(D.get(0)) < 0.9 ? v3(1, 0, 0) : v3(0, 1, 0);
    const E0 = cross(D, helper);
    normalize(E0);
    const E1 = cross(D, E0);
    normalize(E1);
    return { E0, E1 };
}

// Points on the single-sided infinite cone with vertex V, unit axis D and
// half-angle theta, at heights h (measured along the axis).
function conePoints(V: Vector, D: Vector, theta: number,
    heights: readonly number[], numAngles: number): Vector[] {
    const { E0, E1 } = basisFor(D);
    const points: Vector[] = [];
    for (const h of heights) {
        const radius = h * Math.tan(theta);
        for (let i = 0; i < numAngles; ++i) {
            const u = (2 * Math.PI * i) / numAngles;
            const radial = add(mul(radius * Math.cos(u), E0),
                mul(radius * Math.sin(u), E1));
            points.push(add(add(V, mul(h, D)), radial));
        }
    }
    return points;
}

// The exact ellipse cut from the cone (V,D,theta) by the plane
// Dot(N,X) = k. See the derivation in the port notes: with
// cosPhi = Dot(D,N), sinPhi = |D - cosPhi*N|, U = (D - cosPhi*N)/sinPhi and
// h0 = k - Dot(N,V), the section has eccentricity e = sinPhi/cos(theta),
// semimajor a = |h0|*sin(theta)/(cos(theta)*(1-e^2)) along U, semiminor
// b = a*sqrt(1-e^2) along Cross(N,U), and center
// C = V + sinPhi*cosPhi*h0/(cos^2(theta)*(1-e^2))*U + h0*N.
function exactSection(V: Vector, D: Vector, theta: number, N: Vector,
    k: number): Ellipse3 {
    const csPhi = dot(D, N);
    let U = sub(D, mul(csPhi, N));
    const snPhi = normalize(U);
    const csTheta = Math.cos(theta);
    const snTheta = Math.sin(theta);
    const e = snPhi / csTheta;
    expect(e).toBeLessThan(1);
    const omesqr = 1 - e * e;
    const h0 = k - dot(N, V);
    const a = Math.abs(h0) * snTheta / (csTheta * omesqr);
    const b = a * Math.sqrt(omesqr);
    const sc = (snPhi * csPhi * h0) / (csTheta * csTheta * omesqr);
    if (snPhi === 0) {
        // A circular section; any in-plane orthonormal pair works.
        const { E0, E1 } = basisFor(N);
        return Ellipse3.fromCenterNormalAxisExtent(
            add(V, mul(h0, N)), N, [E0, E1], Vector.fromArray([a, b]));
    }
    const C = add(add(V, mul(sc, U)), mul(h0, N));
    const W = cross(N, U);
    return Ellipse3.fromCenterNormalAxisExtent(C, N, [U, W],
        Vector.fromArray([a, b]));
}

// A deterministic pseudorandom generator so the randomized tests repeat.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// The unsigned angular error between the fitted axis and the true axis, plus
// the vertex error and the angle error.
function fitErrors(cone: { ray: { origin: Vector; direction: Vector };
    angle: number }, V: Vector, D: Vector, theta: number) {
    return {
        vertex: Math.hypot(...sub(cone.ray.origin, V).values),
        axis: Math.hypot(...sub(cone.ray.direction, D).values),
        angle: Math.abs(cone.angle - theta)
    };
}

describe('ApprCone3EllipseAndPointsControl', () => {
    it('has the upstream default parameters and validates them', () => {
        const control = new ApprCone3EllipseAndPointsControl();
        expect(control.penalty).toBe(1);
        expect(control.maxSubdivisions).toBe(8);
        expect(control.maxBisections).toBe(64);
        expect(control.epsilon).toBe(1e-08);
        expect(control.tolerance).toBe(1e-04);
        expect(control.padding).toBe(1e-03);
        expect(control.validParameters()).toBe(true);
    });

    it('rejects nonpositive parameters', () => {
        const fields: Array<[keyof ApprCone3EllipseAndPointsControl, number]> = [
            ['penalty', 0], ['maxSubdivisions', 0], ['maxBisections', 0],
            ['epsilon', 0], ['tolerance', 0], ['padding', 0],
            ['penalty', -1], ['padding', -1]
        ];
        for (const [name, value] of fields) {
            const control = new ApprCone3EllipseAndPointsControl();
            (control as unknown as Record<string, number>)[name as string] = value;
            expect(control.validParameters()).toBe(false);
        }
    });
});

describe('ApprCone3EllipseAndPoints.fit', () => {
    it('throws on invalid control parameters', () => {
        const control = new ApprCone3EllipseAndPointsControl();
        control.padding = 0;
        const V = v3(0, 0, 0);
        const D = v3(0, 0, 1);
        const theta = 0.5;
        const ellipse = exactSection(V, D, theta, v3(0, 0, 1), 2);
        const points = conePoints(V, D, theta, [1, 3], 8);
        expect(() => ApprCone3EllipseAndPoints.fit(ellipse, points, control))
            .toThrow('Invalid control parameter.');
    });

    it('throws on an empty point set', () => {
        const V = v3(0, 0, 0);
        const D = v3(0, 0, 1);
        const theta = 0.5;
        const ellipse = exactSection(V, D, theta, v3(0, 0, 1), 2);
        expect(() => ApprCone3EllipseAndPoints.fit(ellipse, []))
            .toThrow('ApprCone3EllipseAndPoints: no points.');
    });

    it('throws on a degenerate ellipse extent', () => {
        const ellipse = new Ellipse3();
        ellipse.extent = Vector.fromArray([1, 0]);
        const points = [v3(1, 0, 1)];
        expect(() => ApprCone3EllipseAndPoints.fit(ellipse, points))
            .toThrow('ellipse extents must be positive');
    });

    it('recovers a cone from a circular cross section (axis-aligned)', () => {
        const V = v3(0, 0, 0);
        const D = v3(0, 0, 1);
        const theta = 0.6;
        const ellipse = exactSection(V, D, theta, v3(0, 0, 1), 2);
        // A circular section: both extents equal h0 * tan(theta).
        expect(ellipse.extent.get(0)).toBeCloseTo(2 * Math.tan(theta), 12);
        expect(ellipse.extent.get(1)).toBeCloseTo(2 * Math.tan(theta), 12);

        const points = conePoints(V, D, theta, [1, 3, 5], 16);
        const cone = ApprCone3EllipseAndPoints.fit(ellipse, points);
        const err = fitErrors(cone, V, D, theta);
        expect(err.angle).toBeLessThan(1e-3);
        expect(err.axis).toBeLessThan(1e-3);
        expect(err.vertex).toBeLessThan(1e-3);
        expect(cone.isInfinite()).toBe(true);
    });

    it('recovers a cone from a tilted elliptical cross section', () => {
        const V = v3(1, -2, 3);
        const D = unit(0.3, -0.4, 1);
        const theta = 0.55;
        const N = unit(0.1, 0.25, 1);
        const k = dot(N, add(V, mul(4, D)));
        const ellipse = exactSection(V, D, theta, N, k);
        expect(ellipse.extent.get(0)).toBeGreaterThan(ellipse.extent.get(1));

        const points = conePoints(V, D, theta, [2, 4, 6], 24);
        const cone = ApprCone3EllipseAndPoints.fit(ellipse, points);
        const err = fitErrors(cone, V, D, theta);
        expect(err.angle).toBeLessThan(2e-3);
        expect(err.axis).toBeLessThan(2e-3);
        expect(err.vertex).toBeLessThan(1e-2);
    });

    it('is insensitive to the sign of the ellipse plane normal', () => {
        const V = v3(-1, 0.5, 0);
        const D = unit(-0.2, 0.3, 1);
        const theta = 0.4;
        const N = unit(0.15, -0.1, 1);
        const k = dot(N, add(V, mul(5, D)));
        const points = conePoints(V, D, theta, [3, 5, 7], 20);

        const ellipsePos = exactSection(V, D, theta, N, k);
        const ellipseNeg = exactSection(V, D, theta, mul(-1, N), -k);
        const conePos = ApprCone3EllipseAndPoints.fit(ellipsePos, points);
        const coneNeg = ApprCone3EllipseAndPoints.fit(ellipseNeg, points);

        for (const cone of [conePos, coneNeg]) {
            const err = fitErrors(cone, V, D, theta);
            expect(err.angle).toBeLessThan(2e-3);
            expect(err.axis).toBeLessThan(2e-3);
            expect(err.vertex).toBeLessThan(1e-2);
        }
    });

    it('reproduces the exact section when the ellipse alone is consistent', () => {
        // The fitted cone must contain the input ellipse: sampling the
        // ellipse and evaluating the cone equation gives (near) zero.
        const V = v3(0.5, 0.5, -1);
        const D = unit(0, 0.5, 1);
        const theta = 0.5;
        const N = unit(0.2, 0, 1);
        const k = dot(N, add(V, mul(6, D)));
        const ellipse = exactSection(V, D, theta, N, k);
        const points = conePoints(V, D, theta, [4, 8], 16);
        const cone = ApprCone3EllipseAndPoints.fit(ellipse, points);

        let maxResidual = 0;
        for (let i = 0; i < 32; ++i) {
            const t = (2 * Math.PI * i) / 32;
            const X = add(add(ellipse.center,
                mul(ellipse.extent.get(0) * Math.cos(t), ellipse.axis[0])),
                mul(ellipse.extent.get(1) * Math.sin(t), ellipse.axis[1]));
            const diff = sub(X, cone.ray.origin);
            const h = dot(cone.ray.direction, diff);
            expect(h).toBeGreaterThan(0);
            const residual = Math.abs(h * h - cone.cosAngleSqr * dot(diff, diff))
                / dot(diff, diff);
            maxResidual = Math.max(maxResidual, residual);
        }
        expect(maxResidual).toBeLessThan(1e-3);
    });

    it('tolerates noise on the sample points', () => {
        const rand = makeRandom(20260902);
        const V = v3(2, 1, -3);
        const D = unit(0.1, -0.2, 1);
        const theta = 0.45;
        const N = unit(-0.2, 0.1, 1);
        const k = dot(N, add(V, mul(5, D)));
        const ellipse = exactSection(V, D, theta, N, k);
        const clean = conePoints(V, D, theta, [3, 5, 7, 9], 24);
        const noisy = clean.map((p) => add(p, v3(
            1e-3 * (2 * rand() - 1), 1e-3 * (2 * rand() - 1),
            1e-3 * (2 * rand() - 1))));

        const cone = ApprCone3EllipseAndPoints.fit(ellipse, noisy);
        const err = fitErrors(cone, V, D, theta);
        expect(err.angle).toBeLessThan(1e-2);
        expect(err.axis).toBeLessThan(1e-2);
        expect(err.vertex).toBeLessThan(5e-2);
    });

    it('recovers randomly generated cones and sections', () => {
        const rand = makeRandom(12345);
        for (let trial = 0; trial < 12; ++trial) {
            const V = v3(4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2);
            const D = unit(rand() - 0.5, rand() - 0.5, 1 + rand());
            // Keep the section elliptic: the plane normal must make an angle
            // phi with the axis satisfying sin(phi) < cos(theta).
            const theta = 0.25 + 0.5 * rand();
            const maxTilt = 0.5 * (GTE_C_HALF_PI - theta);
            const { E0 } = basisFor(D);
            const tilt = maxTilt * rand();
            const N = add(mul(Math.cos(tilt), D), mul(Math.sin(tilt), E0));
            normalize(N);
            const k = dot(N, add(V, mul(3 + 3 * rand(), D)));
            const ellipse = exactSection(V, D, theta, N, k);
            const points = conePoints(V, D, theta, [2, 4, 6], 16);

            const cone = ApprCone3EllipseAndPoints.fit(ellipse, points);
            const err = fitErrors(cone, V, D, theta);
            expect(err.angle).toBeLessThan(5e-3);
            expect(err.axis).toBeLessThan(5e-3);
            expect(err.vertex).toBeLessThan(5e-2);
        }
    });

    it('penalizes points behind the cone vertex plane', () => {
        // The good points lie on the upward cone; a few outliers lie on the
        // reflected (downward) cone. A negligible penalty lets the fit
        // discard the outliers and recover the true cone; a large penalty
        // instead forces a cone whose vertex plane has every point in front
        // of it, trading accuracy for the penalty.
        const V = v3(0, 0, 0);
        const D = v3(0, 0, 1);
        const theta = 0.5;
        const ellipse = exactSection(V, D, theta, v3(0, 0, 1), 3);
        const good = conePoints(V, D, theta, [2, 4], 12);
        const mirrored = good.slice(0, 4).map(
            (p) => v3(p.get(0), p.get(1), -p.get(2)));
        const all = [...good, ...mirrored];

        const numBehind = (cone: { ray: { origin: Vector; direction: Vector } }) =>
            all.filter((p) =>
                dot(cone.ray.direction, sub(p, cone.ray.origin)) < 0).length;

        const cheap = new ApprCone3EllipseAndPointsControl();
        cheap.penalty = 1e-12;
        const cheapCone = ApprCone3EllipseAndPoints.fit(ellipse, all, cheap);
        expect(cheapCone.ray.direction.get(2)).toBeGreaterThan(0.999);
        expect(fitErrors(cheapCone, V, D, theta).angle).toBeLessThan(1e-3);
        expect(numBehind(cheapCone)).toBe(mirrored.length);

        const costly = new ApprCone3EllipseAndPointsControl();
        costly.penalty = 1e4;
        const costlyCone = ApprCone3EllipseAndPoints.fit(ellipse, all, costly);
        expect(costlyCone.ray.direction.get(2)).toBeGreaterThan(0.999);
        expect(numBehind(costlyCone)).toBe(0);
        // The penalty pulled the vertex below all of the points.
        expect(costlyCone.ray.origin.get(2)).toBeLessThan(-4 * Math.tan(theta));
    });

    it('accepts a caller-supplied control that tightens the search', () => {
        const V = v3(0, 0, 0);
        const D = unit(0, 0.2, 1);
        const theta = 0.35;
        const N = unit(0.1, 0, 1);
        const k = dot(N, add(V, mul(5, D)));
        const ellipse = exactSection(V, D, theta, N, k);
        const points = conePoints(V, D, theta, [4, 8], 20);

        const control = new ApprCone3EllipseAndPointsControl();
        control.maxSubdivisions = 16;
        control.tolerance = 1e-08;
        control.epsilon = 1e-12;
        const tight = ApprCone3EllipseAndPoints.fit(ellipse, points, control);
        const loose = ApprCone3EllipseAndPoints.fit(ellipse, points);
        expect(fitErrors(tight, V, D, theta).angle)
            .toBeLessThanOrEqual(fitErrors(loose, V, D, theta).angle + 1e-6);
        expect(fitErrors(tight, V, D, theta).angle).toBeLessThan(1e-4);
    });
});
