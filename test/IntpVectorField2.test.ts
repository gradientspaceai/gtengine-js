import { describe, it, expect } from 'vitest';
import { IntpVectorField2 } from '../src/IntpVectorField2.js';
import { Vector } from '../src/Vector.js';
import { Delaunay2 } from '../src/Delaunay2.js';
import { Delaunay2Mesh } from '../src/Delaunay2Mesh.js';
import { IntpQuadraticNonuniform2 } from '../src/IntpQuadraticNonuniform2.js';
import { check, expectClose, fc, seededRandom, unitVector }
    from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

function gridPoints(n: number, jitter: number): Vector[] {
    const rand = makeRandom(2024);
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

// An affine map of the plane; the Cendes-Wong interpolator reproduces affine
// data exactly, so the vector field interpolation is exact for this map.
function affine(p: Vector): Vector {
    const x = p.values[0];
    const y = p.values[1];
    return Vector.fromArray([2 + 3 * x - y, -1 + 0.5 * x + 4 * y]);
}

function inverseAffine(q: Vector): Vector {
    // Invert [[3,-1],[0.5,4]] * (x,y) = (u - 2, v + 1).
    const u = q.values[0] - 2;
    const v = q.values[1] + 1;
    const det = 3 * 4 - (-1) * 0.5;
    return Vector.fromArray([(4 * u + 1 * v) / det, (-0.5 * u + 3 * v) / det]);
}

describe('IntpVectorField2', () => {
    const domain = gridPoints(6, 0.03);
    const range = domain.map(affine);

    it('reproduces an affine vector field exactly', () => {
        const field = new IntpVectorField2(domain, range);
        const rand = makeRandom(31337);
        for (let k = 0; k < 50; ++k) {
            const p = Vector.fromArray([
                0.05 + 0.9 * rand(),
                0.05 + 0.9 * rand()
            ]);
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            const expected = affine(p);
            expect(result.output.values[0]).toBeCloseTo(expected.values[0], 9);
            expect(result.output.values[1]).toBeCloseTo(expected.values[1], 9);
        }
    });

    it('interpolates the data at the domain points', () => {
        // A nonlinear map: the interpolant must still pass through the data.
        const nonlinear = (p: Vector) => Vector.fromArray([
            Math.sin(2 * p.values[0]) + p.values[1],
            p.values[0] * p.values[0] - Math.cos(p.values[1])
        ]);
        const field = new IntpVectorField2(domain, domain.map(nonlinear));
        for (const p of domain) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            const expected = nonlinear(p);
            expect(result.output.values[0]).toBeCloseTo(expected.values[0], 10);
            expect(result.output.values[1]).toBeCloseTo(expected.values[1], 10);
        }
    });

    it('round-trips through the inverse affine field', () => {
        const forward = new IntpVectorField2(domain, range);
        const backward = new IntpVectorField2(range, domain.map(p => p.clone()));
        const rand = makeRandom(5150);
        for (let k = 0; k < 30; ++k) {
            const p = Vector.fromArray([
                0.1 + 0.8 * rand(),
                0.1 + 0.8 * rand()
            ]);
            const q = forward.evaluate(p);
            expect(q.valid).toBe(true);
            // The forward image agrees with the analytic map, and the
            // analytic inverse recovers the input.
            const back = inverseAffine(q.output);
            expect(back.values[0]).toBeCloseTo(p.values[0], 8);
            expect(back.values[1]).toBeCloseTo(p.values[1], 8);

            // The interpolated inverse field also recovers the input.
            const r = backward.evaluate(q.output);
            expect(r.valid).toBe(true);
            expect(r.output.values[0]).toBeCloseTo(p.values[0], 8);
            expect(r.output.values[1]).toBeCloseTo(p.values[1], 8);
        }
    });

    it('reproduces the identity map', () => {
        const field = new IntpVectorField2(domain, domain.map(p => p.clone()));
        for (const p of [
            Vector.fromArray([0.3, 0.7]),
            Vector.fromArray([0.55, 0.21]),
            Vector.fromArray([0.5, 0.5])
        ]) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(true);
            expect(result.output.values[0]).toBeCloseTo(p.values[0], 9);
            expect(result.output.values[1]).toBeCloseTo(p.values[1], 9);
        }
    });

    it('reports invalid outside the convex hull of the domain', () => {
        const field = new IntpVectorField2(domain, range);
        for (const p of [
            Vector.fromArray([-1, 0.5]),
            Vector.fromArray([0.5, 2]),
            Vector.fromArray([5, 5])
        ]) {
            const result = field.evaluate(p);
            expect(result.valid).toBe(false);
        }
    });

    it('validates its input', () => {
        expect(() => new IntpVectorField2([], [])).toThrow();
        expect(() => new IntpVectorField2(domain, range.slice(1))).toThrow();
        expect(() => new IntpVectorField2(
            [Vector.fromArray([0, 0, 0]), Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0])],
            [Vector.fromArray([0, 0, 0]), Vector.fromArray([1, 0, 0]),
                Vector.fromArray([0, 1, 0])])).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpVectorField2.h.
// ---------------------------------------------------------------------------

// A jittered grid over [0,1]^2 built from a generated seed; the boundary
// points stay on the unit square so the convex hull is fixed.
function seededGrid(n: number, jitter: number, seed: number): Vector[] {
    const rand = seededRandom(seed);
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

const vectorFieldDomain = fc.tuple(
    fc.integer({ min: 4, max: 5 }),
    fc.integer({ min: 1, max: 5 }).map(k => k / 100),
    fc.integer({ min: 0, max: 0xffff })
).map(([n, jitter, seed]) => ({ domain: seededGrid(n, jitter, seed + 1), seed }));

describe('IntpVectorField2 verification', () => {
    // The Cendes-Wong interpolator reproduces affine data exactly when the
    // derivatives it estimates from the mesh normals are exact, which they are
    // for an affine field with spatialDelta = 1. Both components must be
    // reproduced, which catches a swapped x/y repackaging.
    it('reproduces an affine vector field exactly', () => {
        let numChecked = 0;
        check(fc.tuple(vectorFieldDomain,
            fc.array(fc.integer({ min: -4, max: 4 }),
                { minLength: 6, maxLength: 6 })), ([{ domain, seed }, m]) => {
                    const affine = (p: Vector) => Vector.fromArray([
                        m[0] + m[1] * p.values[0] + m[2] * p.values[1],
                        m[3] + m[4] * p.values[0] + m[5] * p.values[1]]);
                    const field = new IntpVectorField2(domain, domain.map(affine));
                    const rand = seededRandom(seed + 11);
                    for (let n = 0; n < 8; ++n) {
                        const P = Vector.fromArray([rand(), rand()]);
                        const r = field.evaluate(P);
                        if (!r.valid) {
                            continue;
                        }
                        ++numChecked;
                        const expected = affine(P);
                        expectClose(r.output.values[0], expected.values[0],
                            1e-9, 1e-10);
                        expectClose(r.output.values[1], expected.values[1],
                            1e-9, 1e-10);
                    }
                    return true;
                }, 20);
        expect(numChecked).toBeGreaterThan(40);
    }, 30000);

    // The interpolant passes through the range data at the domain points.
    it('interpolates the range data at the domain points', () => {
        check(fc.tuple(vectorFieldDomain, fc.integer({ min: 0, max: 0xffff })),
            ([{ domain }, seed]) => {
                const rand = seededRandom(seed + 13);
                const range = domain.map(() => Vector.fromArray(
                    [Math.round(20 * (2 * rand() - 1)),
                        Math.round(20 * (2 * rand() - 1))]));
                const field = new IntpVectorField2(domain, range);
                for (let i = 0; i < domain.length; ++i) {
                    const r = field.evaluate(domain[i]);
                    if (!r.valid) {
                        continue;
                    }
                    expectClose(r.output.values[0], range[i].values[0],
                        1e-8, 1e-9);
                    expectClose(r.output.values[1], range[i].values[1],
                        1e-8, 1e-9);
                }
                return true;
            }, 20);
    }, 30000);

    // The triangulation depends only on the domain, and everything downstream
    // of it - the estimated derivatives, the Bezier coefficients and the
    // quadratic evaluation - is linear in the range data, so the interpolator
    // is linear in the field.
    it('is linear in the range field', () => {
        check(fc.tuple(vectorFieldDomain, fc.integer({ min: -4, max: 4 }),
            fc.integer({ min: -4, max: 4 })), ([{ domain, seed }, a, b]) => {
                const rand = seededRandom(seed + 17);
                const r1 = domain.map(() => Vector.fromArray(
                    [Math.round(10 * (2 * rand() - 1)),
                        Math.round(10 * (2 * rand() - 1))]));
                const r2 = domain.map(() => Vector.fromArray(
                    [Math.round(10 * (2 * rand() - 1)),
                        Math.round(10 * (2 * rand() - 1))]));
                const rc = r1.map((v, i) => Vector.fromArray([
                    a * v.values[0] + b * r2[i].values[0],
                    a * v.values[1] + b * r2[i].values[1]]));
                const f1 = new IntpVectorField2(domain, r1);
                const f2 = new IntpVectorField2(domain, r2);
                const fc2 = new IntpVectorField2(domain, rc);
                for (let n = 0; n < 6; ++n) {
                    const P = Vector.fromArray([rand(), rand()]);
                    const A = f1.evaluate(P);
                    const B = f2.evaluate(P);
                    const C = fc2.evaluate(P);
                    if (!(A.valid && B.valid && C.valid)) {
                        continue;
                    }
                    for (let k = 0; k < 2; ++k) {
                        expectClose(a * A.output.values[k] + b * B.output.values[k],
                            C.output.values[k], 1e-9, 1e-10);
                    }
                }
                return true;
            }, 20);
    }, 30000);

    // The two component interpolators share one triangulation, so the field is
    // exactly the pair of scalar Cendes-Wong interpolants of the components.
    it('equals the pair of component interpolators over the shared mesh', () => {
        check(fc.tuple(vectorFieldDomain, fc.integer({ min: 0, max: 0xffff })),
            ([{ domain }, seed]) => {
                const rand = seededRandom(seed + 19);
                const range = domain.map(() => Vector.fromArray(
                    [Math.round(20 * (2 * rand() - 1)),
                        Math.round(20 * (2 * rand() - 1))]));
                const field = new IntpVectorField2(domain, range);

                const delaunay = new Delaunay2();
                if (!delaunay.compute(domain) || delaunay.getDimension() !== 2) {
                    return true;
                }
                const mesh = new Delaunay2Mesh(delaunay);
                const ix = IntpQuadraticNonuniform2.fromSpatialDelta(mesh,
                    range.map(v => v.values[0]), 1);
                const iy = IntpQuadraticNonuniform2.fromSpatialDelta(mesh,
                    range.map(v => v.values[1]), 1);
                for (let n = 0; n < 6; ++n) {
                    const P = Vector.fromArray([rand(), rand()]);
                    const r = field.evaluate(P);
                    const rx = ix.evaluate(P);
                    const ry = iy.evaluate(P);
                    expect(r.valid).toBe(rx.valid && ry.valid);
                    if (!r.valid) {
                        continue;
                    }
                    expect(r.output.values[0]).toBe(rx.F);
                    expect(r.output.values[1]).toBe(ry.F);
                }
                return true;
            }, 20);
    }, 30000);

    // The upstream contract: valid if and only if the input is in the convex
    // hull of the domain points.
    it('reports points far outside the domain hull as invalid', () => {
        check(fc.tuple(vectorFieldDomain, unitVector(2)), ([{ domain }, dir]) => {
            const range = domain.map(p => Vector.fromArray(
                [p.values[1], -p.values[0]]));
            const field = new IntpVectorField2(domain, range);
            const far = Vector.fromArray(
                [1e3 * dir.values[0], 1e3 * dir.values[1]]);
            expect(field.evaluate(far).valid).toBe(false);
            return true;
        }, 20);
    }, 30000);
});
