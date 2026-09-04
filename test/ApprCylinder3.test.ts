import { describe, it, expect } from 'vitest';
import { ApprCylinder3 } from '../src/ApprCylinder3.js';
import { Cylinder3 } from '../src/Cylinder3.js';
import { Vector, dot, length, normalize, sub } from '../src/Vector.js';
import { cross } from '../src/Vector3.js';
import {
    check, expectClose, expectVectorClose, fc, finite, unitVector, vector
} from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function unit(x: number, y: number, z: number): Vector {
    const v = v3(x, y, z);
    normalize(v);
    return v;
}

// A right-handed orthonormal basis {W, U, V}.
function basisFromAxis(W: Vector): { U: Vector, V: Vector } {
    const w = W.values;
    const U = Math.abs(w[0]) > Math.abs(w[1]) ?
        unit(-w[2], 0, w[0]) : unit(0, w[2], -w[1]);
    return { U, V: cross(W, U) };
}

// Points that lie exactly on the cylinder wall with axis C + h*W, radius R,
// h in [-halfHeight, halfHeight].
function cylinderPoints(C: Vector, W: Vector, R: number, halfHeight: number,
    numTheta: number, numH: number): Vector[] {
    const { U, V } = basisFromAxis(W);
    const points: Vector[] = [];
    for (let j = 0; j < numH; ++j) {
        const h = -halfHeight + 2 * halfHeight * j / (numH - 1);
        for (let i = 0; i < numTheta; ++i) {
            const t = 2 * Math.PI * i / numTheta;
            const c = R * Math.cos(t), s = R * Math.sin(t);
            points.push(v3(
                C.values[0] + c * U.values[0] + s * V.values[0] + h * W.values[0],
                C.values[1] + c * U.values[1] + s * V.values[1] + h * W.values[1],
                C.values[2] + c * U.values[2] + s * V.values[2] + h * W.values[2]));
        }
    }
    return points;
}

// The direction on the hemisphere grid used by the search, so that an exact
// fit is attainable.
function gridDirection(i: number, j: number, numTheta: number,
    numPhi: number): Vector {
    const theta = 2 * Math.PI * i / numTheta;
    const phi = (Math.PI / 2) * j / numPhi;
    return v3(Math.cos(theta) * Math.sin(phi), Math.sin(theta) * Math.sin(phi),
        Math.cos(phi));
}

function makeRandom(seed: number): () => number {
    let state = seed;
    return () => {
        state = (1103515245 * state + 12345) % 2147483648;
        return state / 2147483648;
    };
}

// Distance from a point to the cylinder axis line.
function distanceToAxis(p: Vector, cylinder: Cylinder3): number {
    const diff = sub(p, cylinder.axis.origin);
    const h = dot(diff, cylinder.axis.direction);
    return length(sub(diff, Vector.fromArray(
        cylinder.axis.direction.values.map(c => c * h))));
}

// The two axis lines describe the same line when the directions are parallel
// up to sign and the origin difference is along the direction.
function expectSameAxisLine(cylinder: Cylinder3, C: Vector, W: Vector,
    digits: number): void {
    expect(Math.abs(dot(cylinder.axis.direction, W))).toBeCloseTo(1, digits);
    const diff = sub(cylinder.axis.origin, C);
    const h = dot(diff, W);
    const perp = length(sub(diff, Vector.fromArray(W.values.map(c => c * h))));
    expect(perp).toBeCloseTo(0, digits);
}

describe('ApprCylinder3 fitting to points', () => {
    const C = v3(1, -2, 3);
    const W = v3(0, 0, 1);
    const R = 2;
    const halfHeight = 3;
    const axisAlignedPoints = cylinderPoints(C, W, R, halfHeight, 12, 8);

    it('recovers the cylinder from a specified axis direction', () => {
        const cylinder = new Cylinder3();
        // The factory normalizes the input direction.
        const error = ApprCylinder3.fromCylinderAxis(v3(0, 0, 5))
            .compute(axisAlignedPoints, cylinder);

        expect(Math.abs(error)).toBeLessThan(1e-6);
        expectSameAxisLine(cylinder, C, W, 6);
        expect(cylinder.radius).toBeCloseTo(R, 6);
        expect(cylinder.height).toBeCloseTo(2 * halfHeight, 6);
        // The axis origin is the midpoint of the projected h-interval.
        expect(length(sub(cylinder.axis.origin, C))).toBeCloseTo(0, 6);
    });

    it('recovers the cylinder from the covariance eigenvector', () => {
        // For this sample set the axis direction is the eigenvector of the
        // largest eigenvalue, because the spread along the axis exceeds the
        // spread across it.
        const cylinder = new Cylinder3();
        const error = ApprCylinder3.fromEigenIndex(2)
            .compute(axisAlignedPoints, cylinder);

        expect(Math.abs(error)).toBeLessThan(1e-6);
        expectSameAxisLine(cylinder, C, W, 6);
        expect(cylinder.radius).toBeCloseTo(R, 6);
        expect(cylinder.height).toBeCloseTo(2 * halfHeight, 6);
    });

    it('produces a worse fit for the other covariance eigenvectors', () => {
        // Eigenvectors 0 and 1 span the plane orthogonal to the true axis,
        // so the fit is poor. This documents the upstream behavior of the
        // eigenvector-index selection.
        for (const eigenIndex of [0, 1]) {
            const cylinder = new Cylinder3();
            const error = ApprCylinder3.fromEigenIndex(eigenIndex)
                .compute(axisAlignedPoints, cylinder);
            expect(error).toBeGreaterThan(0.1);
            expect(Math.abs(dot(cylinder.axis.direction, W)))
                .toBeLessThan(1e-6);
        }
    });

    it('recovers the cylinder by hemisphere search, single- and multi-threaded', () => {
        const numTheta = 16, numPhi = 8;
        const axis = gridDirection(3, 5, numTheta, numPhi);
        const center = v3(1, -2, 3);
        const points = cylinderPoints(center, axis, 2, 3, 12, 8);

        const single = new Cylinder3();
        const errorSingle = ApprCylinder3
            .fromHemisphereSearch(0, numTheta, numPhi)
            .compute(points, single);
        expect(Math.abs(errorSingle)).toBeLessThan(1e-6);
        expectSameAxisLine(single, center, axis, 6);
        expect(single.radius).toBeCloseTo(2, 6);
        expect(single.height).toBeCloseTo(6, 6);

        // The port runs the "multithreaded" variant sequentially, which
        // partitions the phi samples differently but visits the same grid.
        for (const numThreads of [1, 2, 3, 4]) {
            const multi = new Cylinder3();
            const errorMulti = ApprCylinder3
                .fromHemisphereSearch(numThreads, numTheta, numPhi)
                .compute(points, multi);
            expect(errorMulti).toBeCloseTo(errorSingle, 12);
            expect(multi.axis.direction.values)
                .toEqual(single.axis.direction.values);
            expect(multi.radius).toBeCloseTo(single.radius, 12);
            expect(multi.height).toBeCloseTo(single.height, 12);
        }
    });

    it('finds the north pole axis by hemisphere search', () => {
        const cylinder = new Cylinder3();
        const error = ApprCylinder3.fromHemisphereSearch(0, 32, 16)
            .compute(axisAlignedPoints, cylinder);
        expect(Math.abs(error)).toBeLessThan(1e-6);
        expect(cylinder.axis.direction.values[0]).toBeCloseTo(0, 12);
        expect(cylinder.axis.direction.values[1]).toBeCloseTo(0, 12);
        expect(Math.abs(cylinder.axis.direction.values[2])).toBeCloseTo(1, 12);
        expect(cylinder.radius).toBeCloseTo(R, 10);
    });

    it('gives a small residual for noisy points', () => {
        const rnd = makeRandom(31415926);
        const noisy = axisAlignedPoints.map(p => v3(
            p.values[0] + 1e-3 * (2 * rnd() - 1),
            p.values[1] + 1e-3 * (2 * rnd() - 1),
            p.values[2] + 1e-3 * (2 * rnd() - 1)));
        const cylinder = new Cylinder3();
        const error = ApprCylinder3.fromCylinderAxis(v3(0, 0, 1))
            .compute(noisy, cylinder);

        expect(Math.abs(error)).toBeLessThan(1e-4);
        expect(cylinder.radius).toBeCloseTo(R, 2);
        expectSameAxisLine(cylinder, C, W, 2);
        for (const p of noisy) {
            expect(Math.abs(distanceToAxis(p, cylinder) - cylinder.radius))
                .toBeLessThan(5e-3);
        }
    });

    it('recovers random cylinders from a specified axis (aggregated)', () => {
        const rnd = makeRandom(2718281);
        let worstRadius = 0, worstAxis = 0, worstError = 0;
        for (let trial = 0; trial < 20; ++trial) {
            const axis = unit(2 * rnd() - 1, 2 * rnd() - 1, 2 * rnd() - 1);
            const center = v3(4 * rnd() - 2, 4 * rnd() - 2, 4 * rnd() - 2);
            const radius = 0.5 + 2 * rnd();
            const points = cylinderPoints(center, axis, radius, 2, 11, 7);
            const cylinder = new Cylinder3();
            const error = ApprCylinder3.fromCylinderAxis(axis)
                .compute(points, cylinder);
            worstError = Math.max(worstError, Math.abs(error));
            worstRadius = Math.max(worstRadius,
                Math.abs(cylinder.radius - radius));
            worstAxis = Math.max(worstAxis,
                Math.abs(Math.abs(dot(cylinder.axis.direction, axis)) - 1));
        }
        expect(worstError).toBeLessThan(1e-10);
        expect(worstRadius).toBeLessThan(1e-8);
        expect(worstAxis).toBeLessThan(1e-12);
    });

    it('handles collinear points without throwing', () => {
        // The covariance eigenvector of the largest eigenvalue is the line
        // direction, and the fitted radius degenerates to (nearly) zero.
        const collinear: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            collinear.push(v3(i, 2 * i, 3 * i));
        }
        const cylinder = new Cylinder3();
        const error = ApprCylinder3.fromEigenIndex(2)
            .compute(collinear, cylinder);
        expect(Number.isFinite(error)).toBe(true);
        expect(Math.abs(error)).toBeLessThan(1e-10);
        expect(Math.abs(dot(cylinder.axis.direction, unit(1, 2, 3))))
            .toBeCloseTo(1, 10);
        expect(cylinder.radius).toBeLessThan(1e-6);
    });

    it('handles coplanar points as upstream', () => {
        // Points on an ellipse in the z = 0 plane. The eigenvector for the
        // smallest eigenvalue is the plane normal, so the fitted cylinder has
        // zero height and a radius that is the least-squares circle radius.
        const coplanar: Vector[] = [];
        for (let i = 0; i < 24; ++i) {
            const t = 2 * Math.PI * i / 24;
            coplanar.push(v3(3 * Math.cos(t), 2 * Math.sin(t), 0));
        }
        const cylinder = new Cylinder3();
        const error = ApprCylinder3.fromEigenIndex(0)
            .compute(coplanar, cylinder);
        expect(Number.isFinite(error)).toBe(true);
        expect(error).toBeGreaterThan(0);
        expect(Math.abs(cylinder.axis.direction.values[2])).toBeCloseTo(1, 12);
        expect(cylinder.height).toBeCloseTo(0, 12);
        expect(cylinder.radius).toBeGreaterThan(2);
        expect(cylinder.radius).toBeLessThan(3);
    });

    it('rejects invalid inputs', () => {
        const cylinder = new Cylinder3();
        const fewPoints = axisAlignedPoints.slice(0, 5);
        expect(() => ApprCylinder3.fromEigenIndex(0)
            .compute(fewPoints, cylinder)).toThrow(/at least 6 points/);
        expect(() => ApprCylinder3.fromEigenIndex(3)
            .compute(axisAlignedPoints, cylinder))
            .toThrow(/Eigenvector index is out of range/);
        expect(() => ApprCylinder3.fromCylinderAxis(v3(0, 0, 0))
            .compute(axisAlignedPoints, cylinder))
            .toThrow(/cylinder axis must be nonzero/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 0, 4)
            .compute(axisAlignedPoints, cylinder))
            .toThrow(/theta and psi samples must be positive/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 8, 4, false)
            .compute(axisAlignedPoints, cylinder))
            .toThrow(/computeMesh/);
        expect(() => ApprCylinder3.fromEigenIndex(0).compute(
            [v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(0, 0, 1),
                v3(1, 1, 0), Vector.fromArray([1, 1])], cylinder))
            .toThrow(/3D/);
    });
});

describe('ApprCylinder3 fitting to a mesh', () => {
    // A triangulated cylindrical wall along the z-axis.
    const numTheta = 16, numRings = 5;
    const center = v3(0.5, 0.25, -1);
    const radius = 2, halfHeight = 2;
    const points: Vector[] = [];
    for (let j = 0; j < numRings; ++j) {
        const h = -halfHeight + 2 * halfHeight * j / (numRings - 1);
        for (let i = 0; i < numTheta; ++i) {
            const t = 2 * Math.PI * i / numTheta;
            points.push(v3(center.values[0] + radius * Math.cos(t),
                center.values[1] + radius * Math.sin(t),
                center.values[2] + h));
        }
    }
    const indices: number[] = [];
    for (let j = 0; j < numRings - 1; ++j) {
        for (let i = 0; i < numTheta; ++i) {
            const i1 = (i + 1) % numTheta;
            const a = j * numTheta + i, b = j * numTheta + i1;
            const c = (j + 1) * numTheta + i, d = (j + 1) * numTheta + i1;
            indices.push(a, b, c, b, d, c);
        }
    }

    it('recovers the cylinder from the projected-area measure', () => {
        const cylinder = new Cylinder3();
        ApprCylinder3.fromHemisphereSearch(0, 32, 16, false)
            .computeMesh(points, indices, cylinder);

        expect(Math.abs(cylinder.axis.direction.values[2])).toBeCloseTo(1, 12);
        expect(cylinder.radius).toBeCloseTo(radius, 6);
        expect(cylinder.height).toBeCloseTo(2 * halfHeight, 6);
        expect(cylinder.axis.origin.values[0]).toBeCloseTo(center.values[0], 6);
        expect(cylinder.axis.origin.values[1]).toBeCloseTo(center.values[1], 6);
        expect(cylinder.axis.origin.values[2]).toBeCloseTo(center.values[2], 6);
    });

    it('gives the same mesh fit for the multithreaded partition', () => {
        const single = new Cylinder3();
        ApprCylinder3.fromHemisphereSearch(0, 32, 16, false)
            .computeMesh(points, indices, single);
        for (const numThreads of [1, 3, 4]) {
            const multi = new Cylinder3();
            ApprCylinder3.fromHemisphereSearch(numThreads, 32, 16, false)
                .computeMesh(points, indices, multi);
            expect(multi.axis.direction.values)
                .toEqual(single.axis.direction.values);
            expect(multi.axis.origin.values).toEqual(single.axis.origin.values);
            expect(multi.radius).toBe(single.radius);
            expect(multi.height).toBe(single.height);
        }
    });

    it('rejects invalid mesh inputs', () => {
        const cylinder = new Cylinder3();
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 8, 4)
            .computeMesh(points, indices, cylinder))
            .toThrow(/for fitting to points/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 8, 4, false)
            .computeMesh(points.slice(0, 5), [0, 1, 2, 1, 2, 3], cylinder))
            .toThrow(/at least 6 points and 2 triangles/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 8, 4, false)
            .computeMesh(points, [0, 1, 2], cylinder))
            .toThrow(/at least 6 points and 2 triangles/);
    });
});

// ---------------------------------------------------------------------------
// Property-based verification (VERIFYING.md).

describe('ApprCylinder3 verification', () => {
    // Samples of a cylinder wall, as an arbitrary.
    const config = fc.tuple(vector(3, -5, 5), unitVector(3), finite(0.5, 4),
        finite(0.5, 4), fc.integer({ min: 5, max: 9 }),
        fc.integer({ min: 2, max: 4 }));

    function build(t: [Vector, Vector, number, number, number, number]):
        Vector[] {
        const [C, W, R, halfHeight, numTheta, numH] = t;
        return cylinderPoints(C, W, R, halfHeight, numTheta, numH);
    }

    // The least-squares error function of the cylinder fit, evaluated
    // directly from the returned axis line and radius:
    //   E = (1/n) * sum_i (dist(X[i], axis)^2 - r^2)^2
    // The port computes it through the precomputed moment matrices mF0, mF1,
    // mF2 and mMu in G(...), so this is an independent evaluation.
    function directError(points: readonly Vector[],
        cylinder: Cylinder3): number {
        let sum = 0;
        for (const p of points) {
            const d = distanceToAxis(p, cylinder);
            const term = d * d - cylinder.radius * cylinder.radius;
            sum += term * term;
        }
        return sum / points.length;
    }

    it('reports the error function value of the cylinder it returns', () => {
        check(config, (t) => {
            const points = build(t);
            const cylinder = new Cylinder3();
            const error = ApprCylinder3.fromCylinderAxis(t[1])
                .compute(points, cylinder);
            // The samples lie exactly on a cylinder, so the moments are
            // moderately sized and the two evaluations agree to round-off
            // relative to the moment magnitudes.
            expectClose(error, directError(points, cylinder), 1e-12, 1e-7);
        }, 100);
    });

    it('returns the mean squared axis distance as the squared radius', () => {
        // G(...) sets rsqr = Dot(pVec, mMu) + Dot(PC, PC), which algebraically
        // equals the mean of the squared distances of the samples to the
        // fitted axis line.
        check(config, (t) => {
            const points = build(t);
            const cylinder = new Cylinder3();
            ApprCylinder3.fromCylinderAxis(t[1]).compute(points, cylinder);
            let sum = 0;
            for (const p of points) {
                const d = distanceToAxis(p, cylinder);
                sum += d * d;
            }
            expectClose(cylinder.radius * cylinder.radius, sum / points.length,
                1e-12, 1e-8);
        }, 100);
    });

    it('recovers a cylinder from its exact axis direction', () => {
        check(config, (t) => {
            const [C, W, R, halfHeight] = t;
            const points = build(t);
            const cylinder = new Cylinder3();
            const error = ApprCylinder3.fromCylinderAxis(W)
                .compute(points, cylinder);
            expect(Math.abs(error)).toBeLessThan(1e-8);
            expect(Math.abs(dot(cylinder.axis.direction, W))).toBeGreaterThan(
                1 - 1e-12);
            expectClose(cylinder.radius, R, 1e-8, 1e-8);
            expectClose(cylinder.height, 2 * halfHeight, 1e-7, 1e-8);
            // The axis origin is the midpoint of the projected h-interval,
            // which for a symmetric sample set is the cylinder center.
            expectClose(length(sub(cylinder.axis.origin, C)), 0, 1e-7, 1e-8);
        }, 100);
    });

    it('is equivariant under translation of the samples', () => {
        // Preprocess subtracts the sample average, so a translation only
        // moves the axis origin.
        check(fc.tuple(config, vector(3, -8, 8)), ([t, shift]) => {
            const points = build(t);
            const moved = points.map(p => v3(p.values[0] + shift.values[0],
                p.values[1] + shift.values[1], p.values[2] + shift.values[2]));
            const c0 = new Cylinder3();
            const e0 = ApprCylinder3.fromCylinderAxis(t[1])
                .compute(points, c0);
            const c1 = new Cylinder3();
            const e1 = ApprCylinder3.fromCylinderAxis(t[1])
                .compute(moved, c1);
            expectClose(e0, e1, 1e-9, 1e-7);
            expectClose(c0.radius, c1.radius, 1e-9, 1e-9);
            expectClose(c0.height, c1.height, 1e-9, 1e-9);
            for (let k = 0; k < 3; ++k) {
                expectClose(c0.axis.origin.values[k] + shift.values[k],
                    c1.axis.origin.values[k], 1e-7, 1e-8);
                expectClose(c0.axis.direction.values[k],
                    c1.axis.direction.values[k], 1e-12, 1e-12);
            }
        }, 60);
    });

    it('gives the same hemisphere-search result for every thread count', () => {
        // The port keeps the upstream per-thread partition of the phi samples
        // (numPhiSamples / numThreads, with the last thread taking the
        // remainder plus the north pole) and the order of the final
        // min-reduction. Every partition covers exactly the same direction
        // grid, so the winning sample - and therefore the whole cylinder -
        // must be identical to the single-threaded search.
        check(fc.tuple(config, fc.integer({ min: 4, max: 9 }),
            fc.integer({ min: 2, max: 6 })),
            ([t, numTheta, numPhi]) => {
                const points = build(t);
                const single = new Cylinder3();
                const errorSingle = ApprCylinder3
                    .fromHemisphereSearch(0, numTheta, numPhi)
                    .compute(points, single);
                for (const numThreads of [1, 2, 3, 5, 8]) {
                    const multi = new Cylinder3();
                    const errorMulti = ApprCylinder3
                        .fromHemisphereSearch(numThreads, numTheta, numPhi)
                        .compute(points, multi);
                    expectClose(errorSingle, errorMulti, 0, 1e-14);
                    expectClose(single.radius, multi.radius, 0, 1e-14);
                    expectClose(single.height, multi.height, 0, 1e-14);
                    expectVectorClose(single.axis.origin, multi.axis.origin,
                        1e-12, 1e-12);
                    expectVectorClose(single.axis.direction,
                        multi.axis.direction, 0, 1e-14);
                }
            }, 25);
    });

    it('gives the same mesh fit for every thread count', () => {
        check(fc.tuple(config, fc.integer({ min: 4, max: 8 }),
            fc.integer({ min: 2, max: 5 })), ([t, numTheta, numPhi]) => {
            const points = build(t);
            // A triangle strip over consecutive point triples is enough for
            // the projected-area measure.
            const indices: number[] = [];
            for (let i = 0; i + 2 < points.length; ++i) {
                indices.push(i, i + 1, i + 2);
            }
            const single = new Cylinder3();
            ApprCylinder3.fromHemisphereSearch(0, numTheta, numPhi, false)
                .computeMesh(points, indices, single);
            for (const numThreads of [1, 2, 4]) {
                const multi = new Cylinder3();
                ApprCylinder3
                    .fromHemisphereSearch(numThreads, numTheta, numPhi, false)
                    .computeMesh(points, indices, multi);
                expectClose(single.radius, multi.radius, 0, 1e-14);
                expectClose(single.height, multi.height, 0, 1e-14);
                expectVectorClose(single.axis.direction, multi.axis.direction,
                    0, 1e-14);
                expectVectorClose(single.axis.origin, multi.axis.origin,
                    1e-12, 1e-12);
            }
        }, 20);
    });

    it('rejects the documented degenerate inputs', () => {
        const points = cylinderPoints(v3(0, 0, 0), v3(0, 0, 1), 1, 1, 6, 3);
        const cylinder = new Cylinder3();

        expect(() => ApprCylinder3.fromEigenIndex(0)
            .compute(points.slice(0, 5), cylinder))
            .toThrow(/at least 6 points/);
        expect(() => ApprCylinder3.fromEigenIndex(3).compute(points, cylinder))
            .toThrow(/out of range/);
        expect(() => ApprCylinder3.fromCylinderAxis(v3(0, 0, 0))
            .compute(points, cylinder)).toThrow(/must be nonzero/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 0, 4)
            .compute(points, cylinder)).toThrow(/must be positive/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 4, 4, false)
            .compute(points, cylinder)).toThrow(/fitting to a mesh/);
        expect(() => ApprCylinder3.fromHemisphereSearch(0, 4, 4)
            .computeMesh(points, [0, 1, 2, 1, 2, 3], cylinder))
            .toThrow(/fitting to points/);
    });
});
