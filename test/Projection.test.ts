import { describe, it, expect } from 'vitest';
import { AxisAngle } from '../src/AxisAngle.js';
import { GTE_C_PI } from '../src/Constants.js';
import { Hyperellipsoid } from '../src/Hyperellipsoid.js';
import { Hyperplane } from '../src/Hyperplane.js';
import { Line } from '../src/Line.js';
import { mulMatrix } from '../src/Matrix.js';
import {
    perspectiveProject, projectEllipse2, projectEllipsoid3
} from '../src/Projection.js';
import { Rotation } from '../src/Rotation.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function randomUnitVector3(rand: () => number): Vector {
    const z = 2 * rand() - 1;
    const phi = 2 * GTE_C_PI * rand();
    const r = Math.sqrt(Math.max(0, 1 - z * z));
    return Vector.fromArray([r * Math.cos(phi), r * Math.sin(phi), z]);
}

// A random right-handed orthonormal basis, as the axes of an ellipsoid.
function randomAxes3(rand: () => number): Vector[] {
    const axis = randomUnitVector3(rand);
    const angle = 2 * GTE_C_PI * rand();
    const R = Rotation.fromAxisAngle(new AxisAngle(axis, angle)).toMatrix();
    return [R.getCol(0), R.getCol(1), R.getCol(2)];
}

function makeEllipsoid(center: Vector, axes: readonly Vector[],
    extent: Vector): Hyperellipsoid {
    return Hyperellipsoid.fromCenterAxisExtent(center, axes, extent);
}

function makeSphere3(center: Vector, radius: number): Hyperellipsoid {
    return makeEllipsoid(center,
        [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)],
        Vector.fromArray([radius, radius, radius]));
}

// Points on the ellipsoid boundary, X = C + sum_i e_i * y_i * U_i with
// |y| = 1.
function ellipsoidBoundaryPoints(ellipsoid: Hyperellipsoid,
    steps: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i <= steps; ++i) {
        const theta = (GTE_C_PI * i) / steps;
        const st = Math.sin(theta), ct = Math.cos(theta);
        for (let j = 0; j < 2 * steps; ++j) {
            const phi = (GTE_C_PI * j) / steps;
            const y = [st * Math.cos(phi), st * Math.sin(phi), ct];
            let X = ellipsoid.center.clone();
            for (let d = 0; d < 3; ++d) {
                X = add(X, mul(ellipsoid.axis[d],
                    ellipsoid.extent.values[d] * y[d]));
            }
            points.push(X);
        }
    }
    return points;
}

function ellipseBoundaryPoints(ellipse: Hyperellipsoid,
    steps: number): Vector[] {
    const points: Vector[] = [];
    for (let i = 0; i < steps; ++i) {
        const phi = (2 * GTE_C_PI * i) / steps;
        const y = [Math.cos(phi), Math.sin(phi)];
        let X = ellipse.center.clone();
        for (let d = 0; d < 2; ++d) {
            X = add(X, mul(ellipse.axis[d], ellipse.extent.values[d] * y[d]));
        }
        points.push(X);
    }
    return points;
}

describe('projectEllipsoid3', () => {
    it('projects a sphere onto a line as an interval of the sphere radius',
        () => {
            const rand = makeRandom(1234);
            for (let trial = 0; trial < 100; ++trial) {
                const center = Vector.fromArray([
                    10 * rand() - 5, 10 * rand() - 5, 10 * rand() - 5]);
                const radius = 0.1 + 3 * rand();
                const sphere = makeSphere3(center, radius);
                const origin = Vector.fromArray([
                    10 * rand() - 5, 10 * rand() - 5, 10 * rand() - 5]);
                const line = Line.fromOriginDirection(origin,
                    randomUnitVector3(rand));
                const { smin, smax } = projectEllipsoid3(sphere, line);
                const center1D = dot(line.direction, sub(center, origin));
                expect(Math.abs(smin - (center1D - radius)))
                    .toBeLessThanOrEqual(1e-12);
                expect(Math.abs(smax - (center1D + radius)))
                    .toBeLessThanOrEqual(1e-12);
                expect(smax - smin).toBeCloseTo(2 * radius, 12);
            }
        });

    it('projects an axis-aligned ellipsoid onto the coordinate axes as its '
        + 'extents', () => {
            const extent = Vector.fromArray([3, 5, 0.25]);
            const center = Vector.fromArray([-1, 2, 7]);
            const ellipsoid = makeEllipsoid(center,
                [Vector.unit(3, 0), Vector.unit(3, 1), Vector.unit(3, 2)],
                extent);
            for (let d = 0; d < 3; ++d) {
                const line = Line.fromOriginDirection(new Vector(3),
                    Vector.unit(3, d));
                const { smin, smax } = projectEllipsoid3(ellipsoid, line);
                expect(smin).toBeCloseTo(center.values[d] - extent.values[d],
                    12);
                expect(smax).toBeCloseTo(center.values[d] + extent.values[d],
                    12);
            }
        });

    it('bounds every boundary point and is tight', () => {
        const rand = makeRandom(9090);
        let worstSlack = Number.MAX_VALUE;
        for (let trial = 0; trial < 20; ++trial) {
            const center = Vector.fromArray([
                4 * rand() - 2, 4 * rand() - 2, 4 * rand() - 2]);
            const extent = Vector.fromArray([
                0.2 + 2 * rand(), 0.2 + 2 * rand(), 0.2 + 2 * rand()]);
            const ellipsoid = makeEllipsoid(center, randomAxes3(rand),
                extent);
            const line = Line.fromOriginDirection(
                Vector.fromArray([rand(), rand(), rand()]),
                randomUnitVector3(rand));
            const { smin, smax } = projectEllipsoid3(ellipsoid, line);
            let observedMin = Number.MAX_VALUE;
            let observedMax = -Number.MAX_VALUE;
            for (const X of ellipsoidBoundaryPoints(ellipsoid, 48)) {
                const s = dot(line.direction, sub(X, line.origin));
                observedMin = Math.min(observedMin, s);
                observedMax = Math.max(observedMax, s);
            }
            expect(observedMin).toBeGreaterThanOrEqual(smin - 1e-12);
            expect(observedMax).toBeLessThanOrEqual(smax + 1e-12);
            // The sampled extremes come close to the exact interval.
            worstSlack = Math.min(worstSlack,
                (observedMax - observedMin) / (smax - smin));
        }
        expect(worstSlack).toBeGreaterThan(0.999);
    });

    it('is invariant under translating the line origin along the direction',
        () => {
            const ellipsoid = makeEllipsoid(Vector.fromArray([1, 2, 3]),
                randomAxes3(makeRandom(7)),
                Vector.fromArray([1, 2, 0.5]));
            const D = randomUnitVector3(makeRandom(8));
            const P = Vector.fromArray([0, 0, 0]);
            const a = projectEllipsoid3(ellipsoid,
                Line.fromOriginDirection(P, D));
            const b = projectEllipsoid3(ellipsoid,
                Line.fromOriginDirection(add(P, mul(D, 2.5)), D));
            expect(a.smin - b.smin).toBeCloseTo(2.5, 12);
            expect(a.smax - b.smax).toBeCloseTo(2.5, 12);
        });

    it('throws for the wrong dimension', () => {
        expect(() => projectEllipsoid3(new Hyperellipsoid(2), new Line(2)))
            .toThrow('Projection: expecting dimension 3.');
    });
});

describe('projectEllipse2', () => {
    it('projects a circle onto a line as an interval of the circle radius',
        () => {
            const rand = makeRandom(5150);
            for (let trial = 0; trial < 100; ++trial) {
                const center = Vector.fromArray([6 * rand() - 3,
                    6 * rand() - 3]);
                const radius = 0.1 + 2 * rand();
                const circle = Hyperellipsoid.fromCenterAxisExtent(center,
                    [Vector.unit(2, 0), Vector.unit(2, 1)],
                    Vector.fromArray([radius, radius]));
                const phi = 2 * GTE_C_PI * rand();
                const line = Line.fromOriginDirection(
                    Vector.fromArray([rand(), rand()]),
                    Vector.fromArray([Math.cos(phi), Math.sin(phi)]));
                const { smin, smax } = projectEllipse2(circle, line);
                const center1D = dot(line.direction, sub(center, line.origin));
                expect(smin).toBeCloseTo(center1D - radius, 12);
                expect(smax).toBeCloseTo(center1D + radius, 12);
            }
        });

    it('bounds every boundary point and is tight', () => {
        const rand = makeRandom(4321);
        for (let trial = 0; trial < 20; ++trial) {
            const phi = 2 * GTE_C_PI * rand();
            const u0 = Vector.fromArray([Math.cos(phi), Math.sin(phi)]);
            const u1 = Vector.fromArray([-Math.sin(phi), Math.cos(phi)]);
            const ellipse = Hyperellipsoid.fromCenterAxisExtent(
                Vector.fromArray([4 * rand() - 2, 4 * rand() - 2]),
                [u0, u1],
                Vector.fromArray([0.2 + 2 * rand(), 0.2 + 2 * rand()]));
            const psi = 2 * GTE_C_PI * rand();
            const line = Line.fromOriginDirection(
                Vector.fromArray([rand(), rand()]),
                Vector.fromArray([Math.cos(psi), Math.sin(psi)]));
            const { smin, smax } = projectEllipse2(ellipse, line);
            let observedMin = Number.MAX_VALUE;
            let observedMax = -Number.MAX_VALUE;
            for (const X of ellipseBoundaryPoints(ellipse, 2048)) {
                const s = dot(line.direction, sub(X, line.origin));
                observedMin = Math.min(observedMin, s);
                observedMax = Math.max(observedMax, s);
            }
            expect(observedMin).toBeGreaterThanOrEqual(smin - 1e-12);
            expect(observedMax).toBeLessThanOrEqual(smax + 1e-12);
            expect((observedMax - observedMin) / (smax - smin))
                .toBeGreaterThan(0.999);
        }
    });

    it('throws for the wrong dimension', () => {
        expect(() => projectEllipse2(new Hyperellipsoid(3), new Line(3)))
            .toThrow('Projection: expecting dimension 2.');
    });
});

// The view-plane coordinates of the point where the ray from E through X hits
// the plane Dot(N,P) = Dot(N,E) + n, expressed relative to K = E + n*N.
function projectPointToViewPlane(X: Vector, E: Vector, N: Vector, U: Vector,
    V: Vector, n: number): Vector {
    const dir = sub(X, E);
    const s = n / dot(N, dir);
    const P = add(E, mul(dir, s));
    const K = add(E, mul(N, n));
    const W = sub(P, K);
    return Vector.fromArray([dot(W, U), dot(W, V)]);
}

// (Y - center)^T * M * (Y - center), which is 1 on the ellipse boundary.
function ellipseQuadraticForm(ellipse: Hyperellipsoid, Y: Vector): number {
    const D = sub(Y, ellipse.center);
    return dot(D, mulMatrix(ellipse.getM(), D) as Vector);
}

describe('perspectiveProject', () => {
    it('projects a sphere centered on the view axis to a circle of the known '
        + 'radius', () => {
            // With the eyepoint at the origin, the view direction N = (0,0,1)
            // and a sphere of radius r centered at (0,0,zc), the tangent cone
            // has half-angle asin(r/zc) and meets the plane z = n in a circle
            // of radius n*r/sqrt(zc^2 - r^2) centered at K = (0,0,n).
            const E = new Vector(3);
            const N = Vector.unit(3, 2);
            const U = Vector.unit(3, 0);
            const V = Vector.unit(3, 1);
            for (const [zc, r, n] of [[10, 1, 2], [5, 2, 1], [3, 0.5, 2.5],
                [20, 7, 4]]) {
                const sphere = makeSphere3(Vector.fromArray([0, 0, zc]), r);
                const { ellipse, isEllipse } =
                    perspectiveProject(sphere, E, N, U, V, n);
                expect(isEllipse).toBe(true);
                const expected = (n * r) / Math.sqrt(zc * zc - r * r);
                expect(ellipse.center.values[0]).toBeCloseTo(0, 10);
                expect(ellipse.center.values[1]).toBeCloseTo(0, 10);
                expect(ellipse.extent.values[0]).toBeCloseTo(expected, 10);
                expect(ellipse.extent.values[1]).toBeCloseTo(expected, 10);
            }
        });

    it('grows the projected circle as the sphere approaches the eyepoint',
        () => {
            const E = new Vector(3);
            const N = Vector.unit(3, 2);
            const U = Vector.unit(3, 0);
            const V = Vector.unit(3, 1);
            let previous = 0;
            for (const zc of [20, 15, 10, 6, 4]) {
                const { ellipse } = perspectiveProject(
                    makeSphere3(Vector.fromArray([0, 0, zc]), 1), E, N, U, V,
                    2);
                expect(ellipse.extent.values[0]).toBeGreaterThan(previous);
                previous = ellipse.extent.values[0];
            }
        });

    it('projects every ellipsoid boundary point inside the projected ellipse',
        () => {
            const rand = makeRandom(20260901);
            let maxForm = 0;
            let minMaxForm = Number.MAX_VALUE;
            for (let trial = 0; trial < 25; ++trial) {
                // The eyepoint and a view plane whose normal points away.
                const E = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1]);
                const N = randomUnitVector3(rand);
                // Any U perpendicular to N; V completes the right-handed set
                // so that N = Cross(U, V).
                const seed = Math.abs(N.values[0]) > 0.9
                    ? Vector.unit(3, 1) : Vector.unit(3, 0);
                const Uraw = sub(seed, mul(N, dot(seed, N)));
                const U = mul(Uraw, 1 / Math.sqrt(dot(Uraw, Uraw)));
                const V = cross(N, U);
                const n = 1 + 3 * rand();

                // An ellipsoid strictly between the eyepoint and the plane.
                const extent = Vector.fromArray([0.1 + 0.4 * rand(),
                    0.1 + 0.4 * rand(), 0.1 + 0.4 * rand()]);
                const distance = 3 + 5 * rand();
                const offset = mul(randomUnitVector3(rand), 0.6 * rand());
                const center = add(add(E, mul(N, distance)), offset);
                const ellipsoid = makeEllipsoid(center, randomAxes3(rand),
                    extent);

                // The upstream precondition: the ellipsoid projects onto
                // E + s*N with smin > 0.
                const interval = projectEllipsoid3(ellipsoid,
                    Line.fromOriginDirection(E, N));
                expect(interval.smin).toBeGreaterThan(0);

                const { ellipse, isEllipse } =
                    perspectiveProject(ellipsoid, E, N, U, V, n);
                expect(isEllipse).toBe(true);

                let trialMax = 0;
                for (const X of ellipsoidBoundaryPoints(ellipsoid, 40)) {
                    const Y = projectPointToViewPlane(X, E, N, U, V, n);
                    const form = ellipseQuadraticForm(ellipse, Y);
                    expect(form).toBeLessThanOrEqual(1 + 1e-10);
                    trialMax = Math.max(trialMax, form);
                }
                maxForm = Math.max(maxForm, trialMax);
                minMaxForm = Math.min(minMaxForm, trialMax);
            }
            // The silhouette points lie on the ellipse, so the sampled
            // maximum must approach 1 in every trial.
            expect(maxForm).toBeLessThanOrEqual(1 + 1e-10);
            expect(minMaxForm).toBeGreaterThan(0.99);
        });

    it('agrees with the plane overload up to the choice of the {U,V} frame',
        () => {
            const rand = makeRandom(31415);
            for (let trial = 0; trial < 20; ++trial) {
                const E = Vector.fromArray([2 * rand() - 1, 2 * rand() - 1,
                    2 * rand() - 1]);
                const N = randomUnitVector3(rand);
                const n = 1 + 3 * rand();
                // The plane through K = E + n*N with normal N; its constant
                // is Dot(N,K) = Dot(N,E) + n, so the overload recovers n.
                const plane = Hyperplane.fromNormalConstant(N,
                    dot(N, E) + n);

                const center = add(E, mul(N, 4 + 3 * rand()));
                const ellipsoid = makeEllipsoid(center, randomAxes3(rand),
                    Vector.fromArray([0.2 + rand(), 0.2 + rand(),
                        0.2 + rand()]));

                const viaPlane = perspectiveProject(ellipsoid, E, plane);
                expect(viaPlane.isEllipse).toBe(true);

                // Rebuild the same {U,V} the plane overload uses is not
                // necessary: the ellipse's extents are frame independent.
                const seed = Math.abs(N.values[0]) > 0.9
                    ? Vector.unit(3, 1) : Vector.unit(3, 0);
                const Uraw = sub(seed, mul(N, dot(seed, N)));
                const U = mul(Uraw, 1 / Math.sqrt(dot(Uraw, Uraw)));
                const V = cross(N, U);
                const viaBasis =
                    perspectiveProject(ellipsoid, E, N, U, V, n);

                for (let d = 0; d < 2; ++d) {
                    expect(viaPlane.ellipse.extent.values[d])
                        .toBeCloseTo(viaBasis.ellipse.extent.values[d], 9);
                }
                // The 3D ellipse centers agree: K + y0*U + y1*V.
                const K = add(E, mul(N, n));
                const basisCenter = add(add(K,
                    mul(U, viaBasis.ellipse.center.values[0])),
                    mul(V, viaBasis.ellipse.center.values[1]));
                // The plane overload's frame, recomputed the same way.
                const planeU = Vector.fromArray(
                    orthogonalComplementU(plane.normal));
                const planeV = cross(plane.normal, planeU);
                const planeCenter = add(add(K,
                    mul(planeU, viaPlane.ellipse.center.values[0])),
                    mul(planeV, viaPlane.ellipse.center.values[1]));
                for (let d = 0; d < 3; ++d) {
                    expect(Math.abs(planeCenter.values[d]
                        - basisCenter.values[d])).toBeLessThanOrEqual(1e-9);
                }
            }
        });

    it('throws for the wrong dimension', () => {
        const sphere = makeSphere3(Vector.fromArray([0, 0, 5]), 1);
        expect(() => perspectiveProject(sphere, new Vector(3),
            new Hyperplane(2))).toThrow('Projection: expecting dimension 3.');
        expect(() => perspectiveProject(new Hyperellipsoid(2), new Vector(3),
            Vector.unit(3, 2), Vector.unit(3, 0), Vector.unit(3, 1), 1))
            .toThrow('Projection: expecting dimension 3.');
    });
});

// The {U,V} frame that Projection's plane overload builds, reproduced here
// via the same computeOrthogonalComplement3 rule (v1 chosen from the larger
// of |n0| and |n1|, then Gram-Schmidt).
function orthogonalComplementU(N: Vector): number[] {
    const v: Vector[] = [N.clone(), new Vector(3), new Vector(3)];
    if (Math.abs(v[0].values[0]) > Math.abs(v[0].values[1])) {
        v[1] = Vector.fromArray([-v[0].values[2], 0, +v[0].values[0]]);
    } else {
        v[1] = Vector.fromArray([0, +v[0].values[2], -v[0].values[1]]);
    }
    const inv = 1 / Math.sqrt(dot(v[1], v[1]));
    return [v[1].values[0] * inv, v[1].values[1] * inv, v[1].values[2] * inv];
}
