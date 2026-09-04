import { describe, expect, it } from 'vitest';
import { BasisFunctionInput } from '../src/BasisFunction.js';
import { BSplineGeodesic } from '../src/BSplineGeodesic.js';
import { BSplineSurface } from '../src/BSplineSurface.js';
import { GVector } from '../src/GVector.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function gv(u: number, v: number): GVector {
    return GVector.fromArray([u, v]);
}

// A bilinear (degree 1, two controls per direction) B-spline surface. The
// control points are ordered control[i0 + 2*i1].
function bilinear(c00: Vector, c10: Vector, c01: Vector,
    c11: Vector): BSplineSurface {
    const input = [new BasisFunctionInput(2, 1), new BasisFunctionInput(2, 1)];
    return new BSplineSurface(c00.size, input, [c00, c10, c01, c11]);
}

// The parallelogram patch X(u,v) = P + u*A + v*B, which is a plane.
function planePatch(P: Vector, A: Vector, B: Vector): BSplineSurface {
    return bilinear(P, add(P, A), add(P, B), add(add(P, A), B));
}

function surfacePoint(surface: BSplineSurface, point: GVector): Vector {
    const jet = surface.createJet();
    surface.evaluate(point.values[0], point.values[1], 0, jet);
    return jet[0];
}

function euclideanPathLength(surface: BSplineSurface, path: GVector[],
    quantity: number): number {
    let total = 0;
    for (let i = 1; i < quantity; ++i) {
        total += length(sub(surfacePoint(surface, path[i]),
            surfacePoint(surface, path[i - 1])));
    }
    return total;
}

describe('BSplineGeodesic construction', () => {
    it('requires a 3-dimensional spline surface', () => {
        const input = [new BasisFunctionInput(2, 1),
            new BasisFunctionInput(2, 1)];
        const planar = new BSplineSurface(2, input, [
            Vector.fromArray([0, 0]), Vector.fromArray([1, 0]),
            Vector.fromArray([0, 1]), Vector.fromArray([1, 1])]);
        expect(() => new BSplineGeodesic(planar)).toThrow(/3-dimensional/);

        const spatial = planePatch(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));
        const bg = new BSplineGeodesic(spatial);
        expect(bg.getDimension()).toBe(2);
    });
});

describe('BSplineGeodesic on a planar patch', () => {
    it('reproduces the Euclidean metric of an orthonormal patch', () => {
        // X(u,v) = (u,v,0), so the metric is the identity and lengths are
        // Euclidean distances in the parameter domain.
        const surface = planePatch(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));
        const bg = new BSplineGeodesic(surface);

        expect(bg.computeSegmentLength(gv(0.1, 0.2), gv(0.4, 0.6)))
            .toBeCloseTo(0.5, 10);
        expect(bg.computeSegmentLength(gv(0, 0), gv(1, 1)))
            .toBeCloseTo(Math.SQRT2, 10);
        // A plane has vanishing Christoffel symbols, so any segment has zero
        // geodesic curvature.
        expect(bg.computeSegmentCurvature(gv(0.1, 0.2), gv(0.9, 0.4)))
            .toBeLessThan(1e-10);
    });

    it('reproduces the constant Gram metric of a skew patch', () => {
        // X(u,v) = P + u*A + v*B with non-orthogonal A and B. The metric is
        // the constant Gram matrix of {A,B}, so the length of a segment with
        // parameter difference (du,dv) is |du*A + dv*B|.
        const P = v3(-1, 2, 0.5);
        const A = v3(2, 0, 1);
        const B = v3(0.5, 3, -1);
        const surface = planePatch(P, A, B);
        const bg = new BSplineGeodesic(surface);

        for (const [u0, v0, u1, v1] of [
            [0.0, 0.0, 1.0, 1.0], [0.2, 0.7, 0.9, 0.1], [0.5, 0.5, 0.6, 0.9]]) {
            const expected = length(add(mul(u1 - u0, A), mul(v1 - v0, B)));
            expect(bg.computeSegmentLength(gv(u0, v0), gv(u1, v1)))
                .toBeCloseTo(expected, 8);
        }
    });

    it('keeps a geodesic straight in the parameter domain', () => {
        const surface = planePatch(v3(-1, 2, 0.5), v3(2, 0, 1), v3(0.5, 3, -1));
        const bg = new BSplineGeodesic(surface);
        bg.subdivisions = 3;
        bg.refinements = 3;
        bg.searchSamples = 8;
        bg.updateDerivedParameters();

        const end0 = gv(0.1, 0.1);
        const end1 = gv(0.9, 0.7);
        const result = bg.computeGeodesic(end0, end1);
        expect(result.quantity).toBe(9);

        for (let i = 0; i < result.quantity; ++i) {
            const s = i / (result.quantity - 1);
            expect(result.path[i].values[0]).toBeCloseTo(0.1 + 0.8 * s, 5);
            expect(result.path[i].values[1]).toBeCloseTo(0.1 + 0.6 * s, 5);
        }

        // The polyline length equals the straight-line length, and it equals
        // the Euclidean length of the 3D polyline (the surface is a plane).
        const total = bg.computeTotalLength(result.quantity, result.path);
        const straight = bg.computeSegmentLength(end0, end1);
        expect(total).toBeCloseTo(straight, 6);
        expect(euclideanPathLength(surface, result.path, result.quantity))
            .toBeCloseTo(total, 6);
    });
});

describe('BSplineGeodesic on a curved patch', () => {
    // A hyperbolic paraboloid X(u,v) = (u, v, k*u*v) written as a bilinear
    // patch: only the (1,1) corner is lifted.
    function saddle(k: number): BSplineSurface {
        return bilinear(v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0), v3(1, 1, k));
    }

    it('evaluates the metric of the saddle patch', () => {
        // dX/du = (1,0,k*v) and dX/dv = (0,1,k*u), so
        // g = [[1+k^2 v^2, k^2 u v], [k^2 u v, 1+k^2 u^2]].
        const k = 2;
        const surface = saddle(k);
        const bg = new BSplineGeodesic(surface);

        // Along the boundary v = 0 the metric restricted to du is 1, so the
        // segment length is the parameter difference.
        expect(bg.computeSegmentLength(gv(0.2, 0), gv(0.8, 0)))
            .toBeCloseTo(0.6, 8);

        // Along the boundary u = 1 the surface point is (1, v, k*v), a
        // straight line of speed sqrt(1+k^2).
        expect(bg.computeSegmentLength(gv(1, 0.1), gv(1, 0.9)))
            .toBeCloseTo(0.8 * Math.sqrt(1 + k * k), 8);

        // A boundary line of the patch is a straight line in space, hence a
        // geodesic with zero curvature.
        expect(bg.computeSegmentCurvature(gv(1, 0.1), gv(1, 0.9)))
            .toBeLessThan(1e-8);
    });

    it('shortens a diagonal path and keeps the endpoints', () => {
        const surface = saddle(3);
        const bg = new BSplineGeodesic(surface);
        bg.subdivisions = 3;
        bg.refinements = 3;
        bg.searchSamples = 8;
        bg.updateDerivedParameters();

        const end0 = gv(0.05, 0.2);
        const end1 = gv(0.95, 0.85);
        const straight = bg.computeSegmentLength(end0, end1);
        const result = bg.computeGeodesic(end0, end1);
        const total = bg.computeTotalLength(result.quantity, result.path);

        expect(total).toBeLessThan(straight);

        // The chordal 3D length of the polyline is a lower bound for the
        // metric length of the same polyline.
        const euclid = euclideanPathLength(surface, result.path,
            result.quantity);
        expect(euclid).toBeLessThanOrEqual(total + 1e-9);

        // The endpoints are untouched.
        const last = result.quantity - 1;
        expect(result.path[0].values[0]).toBeCloseTo(0.05, 12);
        expect(result.path[0].values[1]).toBeCloseTo(0.2, 12);
        expect(result.path[last].values[0]).toBeCloseTo(0.95, 12);
        expect(result.path[last].values[1]).toBeCloseTo(0.85, 12);

        // The path stays inside the patch domain.
        for (let i = 0; i < result.quantity; ++i) {
            expect(result.path[i].values[0]).toBeGreaterThanOrEqual(0);
            expect(result.path[i].values[0]).toBeLessThanOrEqual(1);
            expect(result.path[i].values[1]).toBeGreaterThanOrEqual(0);
            expect(result.path[i].values[1]).toBeLessThanOrEqual(1);
        }
    });

    it('has positive curvature for a non-geodesic segment', () => {
        const bg = new BSplineGeodesic(saddle(4));
        // The diagonal u = v lies in a plane of symmetry of the saddle, so
        // it is a geodesic; an asymmetric segment is not.
        expect(bg.computeSegmentCurvature(gv(0.1, 0.1), gv(0.9, 0.9)))
            .toBeLessThan(1e-10);
        expect(bg.computeSegmentCurvature(gv(0.1, 0.2), gv(0.9, 0.7)))
            .toBeGreaterThan(1e-3);
    });

    it('matches a bicubic patch that reproduces the saddle', () => {
        // A degree-1 patch and a degree-3 patch that both represent
        // X(u,v) = (u, v, k*u*v) must give the same metric quantities. The
        // bicubic control net for that surface is C[i0][i1] =
        // (i0/3, i1/3, k*(i0/3)*(i1/3)).
        const k = 1.5;
        const input = [new BasisFunctionInput(4, 3),
            new BasisFunctionInput(4, 3)];
        const controls: Vector[] = [];
        for (let i1 = 0; i1 < 4; ++i1) {
            for (let i0 = 0; i0 < 4; ++i0) {
                const x = i0 / 3, y = i1 / 3;
                controls.push(v3(x, y, k * x * y));
            }
        }
        const bicubic = new BSplineSurface(3, input, controls);
        const bg3 = new BSplineGeodesic(bicubic);
        const bg1 = new BSplineGeodesic(saddle(k));

        for (const [u0, v0, u1, v1] of [
            [0.1, 0.2, 0.7, 0.6], [0.3, 0.8, 0.9, 0.15], [0, 0, 1, 1]]) {
            expect(bg3.computeSegmentLength(gv(u0, v0), gv(u1, v1)))
                .toBeCloseTo(bg1.computeSegmentLength(gv(u0, v0), gv(u1, v1)),
                    8);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). The translation hazard here
// is the mapping of the jet slots to the derivatives: jet[1],jet[2] are
// dX/du,dX/dv and jet[3],jet[4],jet[5] are d2X/du2,d2X/dudv,d2X/dv2, and
// computeChristoffel1 must read them as der00,der01,der11. That mapping is
// pinned by the first-kind Christoffel identity
//   d(g_ij)/du_k = <X_ik, X_j> + <X_i, X_jk> = C[j](i,k) + C[i](j,k),
// checked against central differences of the metric, which uses none of the
// same code.
// ---------------------------------------------------------------------------

// Exposes the protected metric and Christoffel state for the checks below.
class ProbeGeodesic extends BSplineGeodesic {
    metricAt(u: number, v: number): number[][] {
        this.computeMetric(gv(u, v));
        return [[this.mMetric.get(0, 0), this.mMetric.get(0, 1)],
            [this.mMetric.get(1, 0), this.mMetric.get(1, 1)]];
    }

    // The first-kind Christoffel symbols at the point of the previous
    // metricAt call (upstream caches the jet in computeMetric).
    christoffelAt(u: number, v: number): number[][][] {
        this.computeMetric(gv(u, v));
        this.computeChristoffel1(gv(u, v));
        return this.mChristoffel1.map(m => [[m.get(0, 0), m.get(0, 1)],
            [m.get(1, 0), m.get(1, 1)]]);
    }

    // computeChristoffel1 ignores its argument and uses the cached jet.
    christoffelWithCacheFrom(cacheU: number, cacheV: number,
        argU: number, argV: number): number[][][] {
        this.computeMetric(gv(cacheU, cacheV));
        this.computeChristoffel1(gv(argU, argV));
        return this.mChristoffel1.map(m => [[m.get(0, 0), m.get(0, 1)],
            [m.get(1, 0), m.get(1, 1)]]);
    }
}

describe('BSplineGeodesic verification', () => {
    // A bicubic Bezier patch (4x4 controls, degree 3, no interior knots), so
    // the surface is a polynomial and central differences are accurate.
    function bicubic(seed: number): BSplineSurface {
        const rand = (() => {
            let state = seed >>> 0;
            return () => {
                state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
                return state / 4294967296;
            };
        })();
        const input = [new BasisFunctionInput(4, 3),
            new BasisFunctionInput(4, 3)];
        const controls: Vector[] = [];
        for (let i1 = 0; i1 < 4; ++i1) {
            for (let i0 = 0; i0 < 4; ++i0) {
                // A graph over a regular (u,v) grid, so the patch is regular
                // and its metric is well conditioned.
                controls.push(v3(i0 / 3, i1 / 3, rand() * 1.5 - 0.75));
            }
        }
        return new BSplineSurface(3, input, controls);
    }

    const sample = fc.tuple(fc.integer({ min: 1, max: 1 << 20 }),
        fc.integer({ min: 15, max: 85 }), fc.integer({ min: 15, max: 85 }));

    it('reproduces the first fundamental form of the patch', () => {
        check(sample, ([seed, ui, vi]) => {
            const surface = bicubic(seed);
            const probe = new ProbeGeodesic(surface);
            const u = ui / 100, v = vi / 100;
            const g = probe.metricAt(u, v);

            // Independent computation from a freshly evaluated jet.
            const jet = surface.createJet();
            surface.evaluate(u, v, 1, jet);
            const xu = jet[1], xv = jet[2];
            expectClose(g[0][0], dot(xu, xu), 1e-12, 1e-12);
            expectClose(g[0][1], dot(xu, xv), 1e-12, 1e-12);
            expectClose(g[1][0], g[0][1], 0, 0);
            expectClose(g[1][1], dot(xv, xv), 1e-12, 1e-12);

            // The metric of a regular patch is symmetric positive definite.
            expect(g[0][0]).toBeGreaterThan(0);
            expect(g[0][0] * g[1][1] - g[0][1] * g[1][0]).toBeGreaterThan(0);
        }, 40);
    });

    it('the Christoffel symbols differentiate the metric', () => {
        const h = 1e-4;
        check(sample, ([seed, ui, vi]) => {
            const surface = bicubic(seed);
            const probe = new ProbeGeodesic(surface);
            // Keep the stencil inside [0,1].
            const u = 0.15 + (ui - 15) * 0.7 / 70;
            const v = 0.15 + (vi - 15) * 0.7 / 70;

            const C = probe.christoffelAt(u, v);
            // Central differences of the metric in each parameter.
            const dg: number[][][] = [];
            for (let k = 0; k < 2; ++k) {
                const plus = k === 0 ? probe.metricAt(u + h, v)
                    : probe.metricAt(u, v + h);
                const minus = k === 0 ? probe.metricAt(u - h, v)
                    : probe.metricAt(u, v - h);
                dg.push([[(plus[0][0] - minus[0][0]) / (2 * h),
                    (plus[0][1] - minus[0][1]) / (2 * h)],
                [(plus[1][0] - minus[1][0]) / (2 * h),
                    (plus[1][1] - minus[1][1]) / (2 * h)]]);
            }

            for (let i = 0; i < 2; ++i) {
                for (let j = 0; j < 2; ++j) {
                    for (let k = 0; k < 2; ++k) {
                        // d(g_ij)/du_k = C[j](i,k) + C[i](j,k)
                        expectClose(dg[k][i][j], C[j][i][k] + C[i][j][k],
                            1e-5, 1e-5);
                    }
                }
            }

            // Each Christoffel matrix is symmetric (mixed partials commute).
            for (let k = 0; k < 2; ++k) {
                expect(C[k][0][1]).toBe(C[k][1][0]);
            }
        }, 25);
    });

    it('caches the jet from computeMetric, as upstream documents', () => {
        check(fc.tuple(sample, fc.integer({ min: 15, max: 85 }),
            fc.integer({ min: 15, max: 85 })),
            ([[seed, ui, vi], au, av]) => {
                const probe = new ProbeGeodesic(bicubic(seed));
                const u = ui / 100, v = vi / 100;
                const expected = probe.christoffelAt(u, v);
                // The argument of computeChristoffel1 is ignored; the cached
                // jet of the previous computeMetric call decides the result.
                const actual = probe.christoffelWithCacheFrom(u, v,
                    au / 100, av / 100);
                expect(actual).toEqual(expected);
            }, 25);
    });

    it('vanishes on a plane and matches the Gram metric there', () => {
        check(fc.tuple(fc.double({ min: -3, max: 3, noNaN: true }),
            fc.double({ min: -3, max: 3, noNaN: true }),
            fc.double({ min: -3, max: 3, noNaN: true }),
            fc.integer({ min: 0, max: 100 }), fc.integer({ min: 0, max: 100 })),
            ([ax, ay, az, ui, vi]) => {
                const A = v3(1 + Math.abs(ax), ay, az);
                const B = v3(ax, 1 + Math.abs(ay), az);
                const probe = new ProbeGeodesic(
                    planePatch(v3(0.5, -1, 2), A, B));
                const u = ui / 100, v = vi / 100;

                const g = probe.metricAt(u, v);
                expectClose(g[0][0], dot(A, A), 1e-9, 1e-9);
                expectClose(g[0][1], dot(A, B), 1e-9, 1e-9);
                expectClose(g[1][1], dot(B, B), 1e-9, 1e-9);

                // A plane has vanishing second derivatives.
                const C = probe.christoffelAt(u, v);
                for (let k = 0; k < 2; ++k) {
                    for (let i = 0; i < 2; ++i) {
                        for (let j = 0; j < 2; ++j) {
                            expect(Math.abs(C[k][i][j])).toBeLessThan(1e-9);
                        }
                    }
                }

                // The segment length is the Euclidean length of the image.
                const p0 = gv(0.1, 0.2), p1 = gv(0.9, 0.7);
                const expected = length(add(mul(0.8, A), mul(0.5, B)));
                expectClose(probe.computeSegmentLength(p0, p1), expected,
                    1e-6, 1e-6);
            }, 30);
    });
});
