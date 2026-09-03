import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh';
import { IntpQuadraticNonuniform2 } from '../src/IntpQuadraticNonuniform2';
import { PlanarMesh } from '../src/PlanarMesh';
import { Vector } from '../src/Vector';

// A deterministic pseudorandom generator so the randomized checks are
// reproducible.
function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// A jittered grid of points over [0,1]^2. The jitter keeps the Delaunay
// triangulation away from the fully cocircular configuration of a regular
// grid.
function gridPoints(n: number, jitter: number): Vector[] {
    const rand = makeRandom(12345);
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

function makeMesh(points: Vector[]): Delaunay2Mesh {
    const delaunay = new Delaunay2();
    const success = delaunay.compute(points);
    expect(success).toBe(true);
    expect(delaunay.getDimension()).toBe(2);
    return new Delaunay2Mesh(delaunay);
}

describe('IntpQuadraticNonuniform2', () => {
    const points = gridPoints(6, 0.03);
    const mesh = makeMesh(points);

    // A linear function; the Cendes-Wong scheme reproduces affine functions
    // exactly when the exact gradients are supplied.
    const a = 0.75, b = -1.5, c = 2.25;
    const linear = (p: Vector) => a + b * p.values[0] + c * p.values[1];
    const F = points.map(linear);
    const FX = points.map(() => b);
    const FY = points.map(() => c);

    it('reproduces the samples at the mesh vertices', () => {
        // A nonlinear sample set: the interpolant must pass through the data.
        const G = points.map(p =>
            Math.sin(3 * p.values[0]) * Math.cos(2 * p.values[1]));
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, G, 1);
        const vertices = mesh.getVertices();
        for (let i = 0; i < vertices.length; ++i) {
            const result = intp.evaluate(vertices[i]);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(G[i], 10);
        }
    });

    it('reproduces an affine function and its gradient exactly '
        + '(explicit derivatives)', () => {
        const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh, F, FX, FY);
        const rand = makeRandom(777);
        for (let k = 0; k < 60; ++k) {
            const p = Vector.fromArray([
                0.05 + 0.9 * rand(),
                0.05 + 0.9 * rand()
            ]);
            const result = intp.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(linear(p), 10);
            expect(result.FX).toBeCloseTo(b, 8);
            expect(result.FY).toBeCloseTo(c, 8);
        }
    });

    it('reproduces an affine function when the derivatives are estimated '
        + 'from mesh normals', () => {
        // For an affine function every triangle normal is the same, so the
        // normal-averaging estimate recovers the exact gradient when
        // spatialDelta is 1.
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, 1);
        const rand = makeRandom(4242);
        for (let k = 0; k < 40; ++k) {
            const p = Vector.fromArray([
                0.05 + 0.9 * rand(),
                0.05 + 0.9 * rand()
            ]);
            const result = intp.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(linear(p), 10);
            expect(result.FX).toBeCloseTo(b, 8);
            expect(result.FY).toBeCloseTo(c, 8);
        }
    });

    it('scales the estimated derivative jets by spatialDelta', () => {
        // The normal-based estimate is FX = spatialDelta * df/dx, so for an
        // affine function the interpolator built with spatialDelta = d must
        // agree with the one built from the explicit jets (d*b, d*c). (The
        // resulting surface is no longer affine for d != 1; only the jets
        // scale.)
        const d = 2;
        const estimated = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, d);
        const explicit = IntpQuadraticNonuniform2.fromDerivatives(mesh, F,
            points.map(() => d * b), points.map(() => d * c));
        for (const p of [
            Vector.fromArray([0.5, 0.5]),
            Vector.fromArray([0.23, 0.71]),
            Vector.fromArray([0.84, 0.36])
        ]) {
            const r0 = estimated.evaluate(p);
            const r1 = explicit.evaluate(p);
            expect(r0.valid).toBe(true);
            expect(r1.valid).toBe(true);
            expect(r0.F).toBeCloseTo(r1.F, 10);
            expect(r0.FX).toBeCloseTo(r1.FX, 10);
            expect(r0.FY).toBeCloseTo(r1.FY, 10);
        }
    });

    it('returns analytic derivatives consistent with finite differences '
        + '(C1 continuity across triangle and subtriangle boundaries)', () => {
        const G = points.map(p =>
            Math.exp(0.5 * p.values[0]) * (1 + p.values[1] * p.values[1]));
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, G, 1);
        const h = 1e-5;
        const rand = makeRandom(99);
        for (let k = 0; k < 40; ++k) {
            const x = 0.15 + 0.7 * rand();
            const y = 0.15 + 0.7 * rand();
            const c0 = intp.evaluate(Vector.fromArray([x, y]));
            const xp = intp.evaluate(Vector.fromArray([x + h, y]));
            const xm = intp.evaluate(Vector.fromArray([x - h, y]));
            const yp = intp.evaluate(Vector.fromArray([x, y + h]));
            const ym = intp.evaluate(Vector.fromArray([x, y - h]));
            expect(c0.valid && xp.valid && xm.valid && yp.valid && ym.valid)
                .toBe(true);
            const dFdx = (xp.F - xm.F) / (2 * h);
            const dFdy = (yp.F - ym.F) / (2 * h);
            expect(Math.abs(dFdx - c0.FX)).toBeLessThan(1e-5);
            expect(Math.abs(dFdy - c0.FY)).toBeLessThan(1e-5);
        }
    });

    it('is continuous across an interior triangle edge', () => {
        const G = points.map(p =>
            p.values[0] * p.values[0] - 2 * p.values[0] * p.values[1]);
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, G, 1);
        // March along a line through the mesh and verify that consecutive
        // samples do not jump; the interpolant is C1, so the increments are
        // bounded by a Lipschitz constant times the step.
        const steps = 400;
        let previous = intp.evaluate(Vector.fromArray([0.05, 0.13]));
        expect(previous.valid).toBe(true);
        for (let k = 1; k <= steps; ++k) {
            const s = k / steps;
            const p = Vector.fromArray([0.05 + 0.9 * s, 0.13 + 0.8 * s]);
            const current = intp.evaluate(p);
            expect(current.valid).toBe(true);
            // Value continuity.
            expect(Math.abs(current.F - previous.F)).toBeLessThan(0.05);
            // Derivative continuity.
            expect(Math.abs(current.FX - previous.FX)).toBeLessThan(0.2);
            expect(Math.abs(current.FY - previous.FY)).toBeLessThan(0.2);
            previous = current;
        }
    });

    it('reports invalid for points outside the convex hull', () => {
        const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh, F, FX, FY);
        const outside = [
            Vector.fromArray([-0.5, 0.5]),
            Vector.fromArray([1.5, 0.5]),
            Vector.fromArray([0.5, -0.25]),
            Vector.fromArray([0.5, 1.25]),
            Vector.fromArray([10, 10])
        ];
        for (const p of outside) {
            const result = intp.evaluate(p);
            expect(result.valid).toBe(false);
            expect(result.F).toBe(0);
        }
    });

    it('works with a PlanarMesh (duck-typed mesh interface)', () => {
        const delaunay = new Delaunay2();
        expect(delaunay.compute(points)).toBe(true);
        const planar = PlanarMesh.fromIndices(delaunay.getVertices(),
            delaunay.getIndices());
        const intp = IntpQuadraticNonuniform2.fromDerivatives(planar, F, FX, FY);
        for (const p of [
            Vector.fromArray([0.3, 0.4]),
            Vector.fromArray([0.62, 0.17]),
            Vector.fromArray([0.5, 0.5])
        ]) {
            const result = intp.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(linear(p), 8);
            expect(result.FX).toBeCloseTo(b, 6);
            expect(result.FY).toBeCloseTo(c, 6);
        }
        expect(intp.evaluate(Vector.fromArray([2, 2])).valid).toBe(false);
    });

    it('handles a single triangle (all edges on the boundary)', () => {
        const tri = [
            Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]),
            Vector.fromArray([0, 1])
        ];
        const triMesh = makeMesh(tri);
        const triF = tri.map(linear);
        const intp = IntpQuadraticNonuniform2.fromDerivatives(triMesh, triF,
            [b, b, b], [c, c, c]);
        for (const p of [
            Vector.fromArray([0.25, 0.25]),
            Vector.fromArray([0.1, 0.7]),
            Vector.fromArray([0.6, 0.2])
        ]) {
            const result = intp.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(linear(p), 10);
            expect(result.FX).toBeCloseTo(b, 8);
            expect(result.FY).toBeCloseTo(c, 8);
        }
    });

    it('validates the input sizes and the point dimension', () => {
        expect(() => IntpQuadraticNonuniform2.fromDerivatives(mesh, [1, 2], FX, FY))
            .toThrow();
        expect(() => IntpQuadraticNonuniform2.fromDerivatives(mesh, F, [1], FY))
            .toThrow();
        const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh, F, FX, FY);
        expect(() => intp.evaluate(Vector.fromArray([0.5, 0.5, 0.5]))).toThrow();
    });
});
