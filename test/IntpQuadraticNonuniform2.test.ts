import { describe, it, expect } from 'vitest';
import { Delaunay2 } from '../src/Delaunay2.js';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh.js';
import { IntpQuadraticNonuniform2 } from '../src/IntpQuadraticNonuniform2.js';
import { PlanarMesh } from '../src/PlanarMesh.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, fc, unitVector }
    from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpQuadraticNonuniform2.h.
// ---------------------------------------------------------------------------

// A small pool of jittered-grid triangulations, built once. Delaunay2Mesh
// computes its barycentrics with exact rational arithmetic and
// IntpQuadraticNonuniform2 calls that four times per triangle, so building a
// fresh mesh inside a fast-check property would dominate the run time (and
// shrinking would rebuild it thousands of times). The properties below vary
// the sample data and the query points instead; each mesh has a different
// size, jitter and connectivity.
function seededGrid(n: number, jitter: number, seed: number): Vector[] {
    const rand = makeRandom(seed);
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        for (let j = 0; j < n; ++j) {
            const x = i / (n - 1);
            const y = j / (n - 1);
            const onBoundary = (i === 0 || i === n - 1 || j === 0 || j === n - 1);
            points.push(Vector.fromArray([
                x + (onBoundary ? 0 : jitter * (2 * rand() - 1)),
                y + (onBoundary ? 0 : jitter * (2 * rand() - 1))]));
        }
    }
    return points;
}

const meshPool = [
    [4, 0.00, 101], [4, 0.05, 202], [5, 0.03, 303], [5, 0.06, 404]
].map(([n, jitter, seed]) => {
    const pts = seededGrid(n, jitter, seed);
    const delaunay = new Delaunay2();
    expect(delaunay.compute(pts)).toBe(true);
    expect(delaunay.getDimension()).toBe(2);
    return { points: pts, mesh: new Delaunay2Mesh(delaunay) };
});

const poolIndex = fc.integer({ min: 0, max: meshPool.length - 1 });

describe('IntpQuadraticNonuniform2 verification', () => {
    // The Cendes-Wong split of each triangle at the incenter, together with
    // the cross-edge points, spans a C1 space that contains the quadratic
    // polynomials. With the exact value and gradient supplied at every vertex
    // the interpolant therefore reproduces any quadratic and its gradient.
    it('reproduces a quadratic polynomial and its gradient', () => {
        let numChecked = 0;
        check(fc.tuple(poolIndex,
            fc.array(fc.integer({ min: -4, max: 4 }),
                { minLength: 6, maxLength: 6 }),
            fc.integer({ min: 0, max: 0xffff })), ([pi, k, seed]) => {
                const { points, mesh } = meshPool[pi];
                const q = (p: Vector) => {
                    const [x, y] = p.values;
                    return k[0] + k[1] * x + k[2] * y + k[3] * x * x
                        + k[4] * x * y + k[5] * y * y;
                };
                const qx = (p: Vector) => {
                    const [x, y] = p.values;
                    return k[1] + 2 * k[3] * x + k[4] * y;
                };
                const qy = (p: Vector) => {
                    const [x, y] = p.values;
                    return k[2] + k[4] * x + 2 * k[5] * y;
                };
                const intp = IntpQuadraticNonuniform2.fromDerivatives(mesh,
                    points.map(q), points.map(qx), points.map(qy));
                const rand = makeRandom(seed + 31);
                for (let n = 0; n < 10; ++n) {
                    const P = Vector.fromArray([rand(), rand()]);
                    const r = intp.evaluate(P);
                    if (!r.valid) {
                        continue;
                    }
                    ++numChecked;
                    expectClose(r.F, q(P), 1e-10, 1e-11);
                    expectClose(r.FX, qx(P), 1e-8, 1e-9);
                    expectClose(r.FY, qy(P), 1e-8, 1e-9);
                }
                return true;
            }, 40);
        expect(numChecked).toBeGreaterThan(100);
    }, 30000);

    // The interpolant reproduces the data at the sample points however the
    // derivatives were obtained.
    it('reproduces the samples at the mesh vertices', () => {
        check(fc.tuple(poolIndex, fc.integer({ min: 0, max: 0xffff })),
            ([pi, seed]) => {
                const { points, mesh } = meshPool[pi];
                const rand = makeRandom(seed + 3);
                const F = points.map(() => Math.round(20 * (2 * rand() - 1)));
                const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, 1);
                const vertices = mesh.getVertices();
                for (let i = 0; i < vertices.length; ++i) {
                    const r = intp.evaluate(vertices[i]);
                    if (!r.valid) {
                        continue;
                    }
                    expectClose(r.F, F[i], 1e-9, 1e-10);
                }
                return true;
            }, 30);
    }, 30000);

    // The whole pipeline - the normal-averaging derivative estimate, the
    // Bezier coefficients and the quadratic evaluation - is linear in the
    // sample data for a fixed mesh: the z-components of the triangle normals
    // and the choice of subtriangle do not depend on F.
    it('is linear in the sample data', () => {
        check(fc.tuple(poolIndex, fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 0, max: 0xffff })), ([pi, a, b, seed]) => {
                const { points, mesh } = meshPool[pi];
                const rand = makeRandom(seed + 17);
                const F1 = points.map(() => Math.round(10 * (2 * rand() - 1)));
                const F2 = points.map(() => Math.round(10 * (2 * rand() - 1)));
                const Fc = F1.map((v, i) => a * v + b * F2[i]);
                const i1 = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F1, 1);
                const i2 = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F2, 1);
                const ic = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, Fc, 1);
                for (let n = 0; n < 6; ++n) {
                    const P = Vector.fromArray([rand(), rand()]);
                    const r1 = i1.evaluate(P);
                    const r2 = i2.evaluate(P);
                    const rc = ic.evaluate(P);
                    if (!(r1.valid && r2.valid && rc.valid)) {
                        continue;
                    }
                    expectClose(a * r1.F + b * r2.F, rc.F, 1e-9, 1e-10);
                    expectClose(a * r1.FX + b * r2.FX, rc.FX, 1e-8, 1e-9);
                    expectClose(a * r1.FY + b * r2.FY, rc.FY, 1e-8, 1e-9);
                }
                return true;
            }, 20);
    }, 30000);

    // The scheme is C1: the value and the gradient are continuous across an
    // interior edge shared by two triangles. Comparing two one-sided
    // evaluations a distance eps apart bounds the jump by O(eps) times the
    // local first/second derivative, so the samples come from a smooth
    // function whose derivatives are O(1) rather than from random data at the
    // grid scale. A gradient discontinuity would show up as a jump of the
    // same order as the gradient itself, far above the bounds below.
    it('has a continuous value and gradient across interior edges', () => {
        let numEdges = 0;
        check(fc.tuple(poolIndex,
            fc.array(fc.integer({ min: -3, max: 3 }),
                { minLength: 5, maxLength: 5 })), ([pi, k]) => {
                const { points, mesh } = meshPool[pi];
                const f = (p: Vector) => {
                    const [x, y] = p.values;
                    return k[0] + k[1] * x + k[2] * y
                        + k[3] * Math.sin(2 * x) + k[4] * Math.cos(2 * y);
                };
                const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh,
                    points.map(f), 1);
                const eps = 1e-6;
                for (let t = 0; t < mesh.getNumTriangles(); ++t) {
                    const V = mesh.getTriangleVertices(t);
                    const adj = mesh.getTriangleAdjacencies(t);
                    if (V === null || adj === null) {
                        continue;
                    }
                    for (let j0 = 2, j1 = 0; j1 < 3; j0 = j1++) {
                        if (adj[j0] < 0) {
                            continue;   // boundary edge
                        }
                        const mx = 0.5 * (V[j0].values[0] + V[j1].values[0]);
                        const my = 0.5 * (V[j0].values[1] + V[j1].values[1]);
                        const ex = V[j1].values[0] - V[j0].values[0];
                        const ey = V[j1].values[1] - V[j0].values[1];
                        const len = Math.hypot(ex, ey);
                        const nx = -ey / len, ny = ex / len;
                        const rp = intp.evaluate(
                            Vector.fromArray([mx + eps * nx, my + eps * ny]));
                        const rm = intp.evaluate(
                            Vector.fromArray([mx - eps * nx, my - eps * ny]));
                        if (!(rp.valid && rm.valid)) {
                            continue;
                        }
                        ++numEdges;
                        const scale = 1 + Math.max(Math.abs(rp.FX),
                            Math.abs(rm.FX), Math.abs(rp.FY), Math.abs(rm.FY));
                        expect(Math.abs(rp.F - rm.F))
                            .toBeLessThan(1e-3 * scale);
                        expect(Math.abs(rp.FX - rm.FX))
                            .toBeLessThan(1e-2 * scale);
                        expect(Math.abs(rp.FY - rm.FY))
                            .toBeLessThan(1e-2 * scale);
                    }
                }
                return true;
            }, 8);
        expect(numEdges).toBeGreaterThan(50);
    }, 30000);

    // The upstream contract: valid if and only if the point is in the convex
    // hull of the mesh vertices.
    it('reports points far outside the hull as invalid', () => {
        const { points, mesh } = meshPool[2];
        const F = points.map((_, i) => i);
        const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh, F, 1);
        check(unitVector(2), dir => {
            const far = Vector.fromArray(
                [1e3 * dir.values[0], 1e3 * dir.values[1]]);
            expect(intp.evaluate(far).valid).toBe(false);
            return true;
        });
    });

    // The mesh normals give the exact gradient of affine data, so the
    // estimated-derivative constructor reproduces an affine function.
    it('reproduces affine data from estimated derivatives', () => {
        check(fc.tuple(poolIndex, fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: 0, max: 0xffff })), ([pi, a, b, c, seed]) => {
                const { points, mesh } = meshPool[pi];
                const f = (p: Vector) => a + b * p.values[0] + c * p.values[1];
                const intp = IntpQuadraticNonuniform2.fromSpatialDelta(mesh,
                    points.map(f), 1);
                const rand = makeRandom(seed + 41);
                for (let n = 0; n < 6; ++n) {
                    const P = Vector.fromArray([rand(), rand()]);
                    const r = intp.evaluate(P);
                    if (!r.valid) {
                        continue;
                    }
                    expectClose(r.F, f(P), 1e-9, 1e-10);
                    expectClose(r.FX, b, 1e-8, 1e-9);
                    expectClose(r.FY, c, 1e-8, 1e-9);
                }
                return true;
            }, 30);
    }, 30000);
});
