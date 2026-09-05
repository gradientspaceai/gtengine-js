import { describe, it, expect } from 'vitest';
import { GTE_C_PI, GTE_C_TWO_PI } from '../src/Constants.js';
import { IntpSphere2 } from '../src/IntpSphere2.js';
import { check, expectClose, fc, finite, seededRandom, unitVector, wellScaled }
    from './helpers/arbitraries.js';

function makeRandom(seed: number): () => number {
    let s = seed >>> 0;
    return () => {
        s = (s * 1664525 + 1013904223) >>> 0;
        return s / 4294967296;
    };
}

// A latitude-longitude sample set. The thetas are in [-pi, pi) and the phis
// are interior; the two antipodal pole points are appended as upstream
// requires for complete spherical coverage.
function latLongSamples(numTheta: number, numPhi: number):
    { theta: number[], phi: number[] } {
    const theta: number[] = [];
    const phi: number[] = [];
    for (let i = 0; i < numTheta; ++i) {
        const t = -GTE_C_PI + i * GTE_C_TWO_PI / numTheta;
        for (let j = 0; j < numPhi; ++j) {
            theta.push(t);
            phi.push(GTE_C_PI * (j + 1) / (numPhi + 1));
        }
    }
    // The poles x = 0, y = 0, |z| = 1.
    theta.push(-GTE_C_PI); phi.push(0);
    theta.push(-GTE_C_PI); phi.push(GTE_C_PI);
    return { theta, phi };
}

describe('IntpSphere2', () => {
    const { theta, phi } = latLongSamples(10, 6);

    describe('getSphericalCoordinates', () => {
        it('converts known unit vectors', () => {
            const north = IntpSphere2.getSphericalCoordinates(0, 0, 1);
            expect(north.theta).toBeCloseTo(-GTE_C_PI, 15);
            expect(north.phi).toBeCloseTo(0, 15);

            const south = IntpSphere2.getSphericalCoordinates(0, 0, -1);
            expect(south.theta).toBeCloseTo(-GTE_C_PI, 15);
            expect(south.phi).toBeCloseTo(GTE_C_PI, 15);

            const xAxis = IntpSphere2.getSphericalCoordinates(1, 0, 0);
            expect(xAxis.theta).toBeCloseTo(0, 15);
            expect(xAxis.phi).toBeCloseTo(GTE_C_PI / 2, 15);

            const yAxis = IntpSphere2.getSphericalCoordinates(0, 1, 0);
            expect(yAxis.theta).toBeCloseTo(GTE_C_PI / 2, 15);
            expect(yAxis.phi).toBeCloseTo(GTE_C_PI / 2, 15);
        });

        it('round-trips through the spherical parameterization', () => {
            const rand = makeRandom(4711);
            for (let k = 0; k < 200; ++k) {
                const t = -GTE_C_PI + GTE_C_TWO_PI * rand();
                const p = GTE_C_PI * rand();
                const x = Math.cos(t) * Math.sin(p);
                const y = Math.sin(t) * Math.sin(p);
                const z = Math.cos(p);
                const sc = IntpSphere2.getSphericalCoordinates(x, y, z);
                expect(sc.phi).toBeCloseTo(p, 10);
                expect(Math.cos(sc.theta) * Math.sin(sc.phi)).toBeCloseTo(x, 10);
                expect(Math.sin(sc.theta) * Math.sin(sc.phi)).toBeCloseTo(y, 10);
                expect(Math.cos(sc.phi)).toBeCloseTo(z, 10);
                expect(-GTE_C_PI <= sc.theta && sc.theta <= GTE_C_PI).toBe(true);
                expect(0 <= sc.phi && sc.phi <= GTE_C_PI).toBe(true);
            }
        });
    });

    it('reproduces a function that is affine in (theta,phi)', () => {
        // f = 1 + 2*phi is affine in the parameter plane and is consistent
        // with the theta-periodic wrap-around copies, so the Cendes-Wong
        // interpolator reproduces it exactly.
        const f = (t: number, p: number) => 1 + 2 * p;
        const F = theta.map((t, i) => f(t, phi[i]));
        const intp = new IntpSphere2(theta, phi, F);
        const rand = makeRandom(8080);
        for (let k = 0; k < 40; ++k) {
            const t = -GTE_C_PI + GTE_C_TWO_PI * rand();
            const p = 0.05 + (GTE_C_PI - 0.1) * rand();
            const result = intp.evaluate(t, p);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(f(t, p), 9);
        }
    });

    it('interpolates the sample data at the sample angles', () => {
        const f = (t: number, p: number) =>
            Math.cos(t) * Math.sin(p) + 0.5 * Math.cos(p);
        const F = theta.map((t, i) => f(t, phi[i]));
        const intp = new IntpSphere2(theta, phi, F);
        for (let i = 0; i < theta.length; ++i) {
            const result = intp.evaluate(theta[i], phi[i]);
            expect(result.valid).toBe(true);
            expect(result.F).toBeCloseTo(F[i], 9);
        }
    });

    it('approximates a smooth spherical function, including across the '
        + 'theta = pi seam', () => {
        // f(x,y,z) = z + 0.25*x is smooth on the sphere and periodic in
        // theta, so the wrap-around copies are consistent.
        const f = (t: number, p: number) =>
            Math.cos(p) + 0.25 * Math.cos(t) * Math.sin(p);
        const build = (numTheta: number, numPhi: number) => {
            const s = latLongSamples(numTheta, numPhi);
            const F = s.theta.map((t, i) => f(t, s.phi[i]));
            return new IntpSphere2(s.theta, s.phi, F);
        };
        const coarse = build(10, 6);
        const intp = build(16, 10);

        const rand = makeRandom(1234);
        let coarseError = 0;
        let maxError = 0;
        for (let k = 0; k < 60; ++k) {
            const t = -GTE_C_PI + GTE_C_TWO_PI * rand();
            const p = 0.2 + (GTE_C_PI - 0.4) * rand();
            const result = intp.evaluate(t, p);
            const coarseResult = coarse.evaluate(t, p);
            expect(result.valid).toBe(true);
            expect(coarseResult.valid).toBe(true);
            maxError = Math.max(maxError, Math.abs(result.F - f(t, p)));
            coarseError =
                Math.max(coarseError, Math.abs(coarseResult.F - f(t, p)));
        }
        // The interpolation converges as the sampling is refined.
        expect(maxError).toBeLessThan(0.03);
        expect(maxError).toBeLessThan(coarseError);

        // Samples straddling the seam theta = +-pi agree with each other,
        // which is what the wrap-around triangulation buys.
        for (const p of [0.5, 1.0, 2.0, 2.7]) {
            const left = intp.evaluate(-GTE_C_PI + 1e-8, p);
            const right = intp.evaluate(GTE_C_PI - 1e-8, p);
            expect(left.valid).toBe(true);
            expect(right.valid).toBe(true);
            // The two evaluations are the same point of the sphere. They
            // agree to interpolation accuracy rather than exactly: the two
            // seam neighborhoods are translates of one another, but the
            // lattice samples are cocircular, so the symbolic perturbation
            // in Delaunay2 can break the ties differently at the two
            // junctions.
            expect(Math.abs(left.F - right.F)).toBeLessThan(5e-3);
            expect(Math.abs(left.F - f(-GTE_C_PI, p))).toBeLessThan(0.03);
        }
    });

    it('reports invalid for angles outside the convex hull of the samples',
        () => {
        const F = theta.map(() => 1);
        const intp = new IntpSphere2(theta, phi, F);
        expect(intp.evaluate(0, -0.5).valid).toBe(false);
        expect(intp.evaluate(0, GTE_C_PI + 0.5).valid).toBe(false);
        expect(intp.evaluate(0, 100).valid).toBe(false);
    });

    it('validates its input', () => {
        expect(() => new IntpSphere2([], [], [])).toThrow();
        expect(() => new IntpSphere2([0, 1], [0, 1], [0])).toThrow();
    });
});

// ---------------------------------------------------------------------------
// Verification pass (V29): property-based checks against upstream
// IntpSphere2.h.
// ---------------------------------------------------------------------------

describe('IntpSphere2 verification', () => {
    // GetSphericalCoordinates inverts the spherical parameterization for every
    // unit vector, with the two documented pole special cases.
    it('round-trips unit vectors through the spherical parameterization', () => {
        check(unitVector(3), v => {
            const [x, y, z] = v.values;
            const { theta: t, phi: p } = IntpSphere2.getSphericalCoordinates(x, y, z);
            expect(t).toBeGreaterThanOrEqual(-GTE_C_PI);
            expect(t).toBeLessThanOrEqual(GTE_C_PI);
            expect(p).toBeGreaterThanOrEqual(0);
            expect(p).toBeLessThanOrEqual(GTE_C_PI);
            // acos is ill-conditioned near |z| = 1: an ulp of z moves phi by
            // about sqrt(2*eps), so the reconstruction of x and y (which are
            // proportional to sin(phi)) carries an absolute error of about
            // eps/sin(phi). The z reconstruction is unaffected.
            const tol = 1e-9 + 1e-15 / Math.max(Math.sin(p), 1e-8);
            expectClose(Math.cos(t) * Math.sin(p), x, tol, 1e-9);
            expectClose(Math.sin(t) * Math.sin(p), y, tol, 1e-9);
            expectClose(Math.cos(p), z, 1e-9, 1e-9);
            return true;
        });
    });

    // The z >= 1 and z <= -1 branches return the poles without calling
    // atan2/acos, so they never produce NaN for a slightly out-of-range z.
    it('returns the poles for |z| at or beyond one', () => {
        check(fc.tuple(finite(-5, 5), finite(-5, 5), finite(1, 10)),
            ([x, y, z]) => {
                const north = IntpSphere2.getSphericalCoordinates(x, y, z);
                expect(north.theta).toBe(-GTE_C_PI);
                expect(north.phi).toBe(0);
                const south = IntpSphere2.getSphericalCoordinates(x, y, -z);
                expect(south.theta).toBe(-GTE_C_PI);
                expect(south.phi).toBe(GTE_C_PI);
                return true;
            });
    });

    // The interpolant passes through the data at every sample angle.
    it('interpolates the sample data at the sample angles', () => {
        check(fc.tuple(fc.integer({ min: 6, max: 10 }),
            fc.integer({ min: 4, max: 6 }), fc.integer({ min: 0, max: 0xffff })),
            ([numTheta, numPhi, seed]) => {
                const s = latLongSamples(numTheta, numPhi);
                const rand = seededRandom(seed + 1);
                const F = s.theta.map(() => Math.round(20 * (2 * rand() - 1)));
                const intp = new IntpSphere2(s.theta, s.phi, F);
                for (let i = 0; i < s.theta.length; ++i) {
                    const r = intp.evaluate(s.theta[i], s.phi[i]);
                    if (!r.valid) {
                        continue;
                    }
                    expectClose(r.F, F[i], 1e-8, 1e-9);
                }
                return true;
            }, 15);
    }, 30000);

    // The construction is equivariant under a rotation about the z-axis: it
    // translates the (theta,phi) sample plane, and both the Delaunay
    // triangulation and the Cendes-Wong scheme are translation invariant. The
    // wrap-around copies at theta +/- 2*pi translate with the data, so the
    // shifted interpolator is the shifted function. (There is no such
    // invariance under a general rotation of the sphere, because the method
    // interpolates in the angle plane rather than on the sphere.)
    //
    // The samples are scattered rather than on a latitude-longitude lattice,
    // and the shift is a dyadic rational: a lattice sample set is massively
    // cocircular, and the round-off of theta + delta then flips some of those
    // ties, which changes the triangulation and hence the interpolant. That
    // is a property of the input, not of the port.
    it('is equivariant under a rotation about the z-axis', () => {
        let numChecked = 0;
        check(fc.tuple(fc.integer({ min: 20, max: 40 }),
            fc.integer({ min: -32, max: 32 }).map(k => k / 64),
            fc.integer({ min: 0, max: 0xffff })),
            ([numPoints, delta, seed]) => {
                const rand = seededRandom(seed + 3);
                const theta: number[] = [];
                const phi: number[] = [];
                const F: number[] = [];
                for (let i = 0; i < numPoints; ++i) {
                    theta.push(-GTE_C_PI + GTE_C_TWO_PI * rand());
                    phi.push(0.05 + (GTE_C_PI - 0.1) * rand());
                    F.push(Math.round(20 * (2 * rand() - 1)));
                }
                theta.push(-GTE_C_PI); phi.push(0); F.push(3);
                theta.push(-GTE_C_PI); phi.push(GTE_C_PI); F.push(-4);
                const base = new IntpSphere2(theta, phi, F);
                const shifted = new IntpSphere2(theta.map(t => t + delta),
                    phi, F);
                for (let k = 0; k < 6; ++k) {
                    const t = -GTE_C_PI + GTE_C_TWO_PI * rand();
                    const p = 0.1 + (GTE_C_PI - 0.2) * rand();
                    const a = base.evaluate(t, p);
                    const b = shifted.evaluate(t + delta, p);
                    if (!a.valid || !b.valid) {
                        continue;
                    }
                    ++numChecked;
                    expectClose(a.F, b.F, 1e-6, 1e-7);
                }
                return true;
            }, 10);
        expect(numChecked).toBeGreaterThan(20);
    }, 30000);

    // A function that is affine in the parameter plane and constant in theta
    // is consistent with the theta-periodic wrap-around copies, so the
    // Cendes-Wong interpolator reproduces it exactly.
    it('reproduces a function affine in phi', () => {
        check(fc.tuple(fc.integer({ min: 6, max: 10 }),
            fc.integer({ min: 4, max: 6 }), fc.integer({ min: -5, max: 5 }),
            fc.integer({ min: -5, max: 5 }), fc.integer({ min: 0, max: 0xffff })),
            ([numTheta, numPhi, a, b, seed]) => {
                const s = latLongSamples(numTheta, numPhi);
                const f = (p: number) => a + b * p;
                const F = s.phi.map(f);
                const intp = new IntpSphere2(s.theta, s.phi, F);
                const rand = seededRandom(seed + 5);
                for (let k = 0; k < 8; ++k) {
                    const t = -GTE_C_PI + GTE_C_TWO_PI * rand();
                    const p = 0.05 + (GTE_C_PI - 0.1) * rand();
                    const r = intp.evaluate(t, p);
                    if (!r.valid) {
                        continue;
                    }
                    expectClose(r.F, f(p), 1e-8, 1e-9);
                }
                return true;
            }, 12);
    }, 30000);

    // The wrap-around copies span theta in [-3*pi, 3*pi] and phi keeps its
    // original range, so a phi far outside [0, pi] is outside the hull.
    it('reports angles far outside the sample hull as invalid', () => {
        const s = latLongSamples(8, 5);
        const F = s.theta.map((t, i) => Math.cos(t) + s.phi[i]);
        const intp = new IntpSphere2(s.theta, s.phi, F);
        check(fc.tuple(finite(-GTE_C_PI, GTE_C_PI), finite(1, 20)),
            ([t, d]) => {
                expect(intp.evaluate(t, GTE_C_PI + d).valid).toBe(false);
                expect(intp.evaluate(t, -d).valid).toBe(false);
                return true;
            });
    });
});
