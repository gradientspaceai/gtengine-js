import { describe, expect, it } from 'vitest';
import { BasisFunctionInput } from '../src/BasisFunction';
import { BSplineGeodesic } from '../src/BSplineGeodesic';
import { BSplineSurface } from '../src/BSplineSurface';
import { GVector } from '../src/GVector';
import { Vector, add, length, mul, sub } from '../src/Vector';

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
