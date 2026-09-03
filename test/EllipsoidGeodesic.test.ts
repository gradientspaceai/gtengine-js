import { describe, expect, it } from 'vitest';
import { EllipsoidGeodesic } from '../src/EllipsoidGeodesic.js';
import { GVector } from '../src/GVector.js';
import { Vector, dot, length, sub } from '../src/Vector.js';

function gv(u: number, v: number): GVector {
    return GVector.fromArray([u, v]);
}

// The length of the polyline through the 3D surface points of the path. This
// is an independent check of the metric-based length computed by the class.
function euclideanPathLength(eg: EllipsoidGeodesic, path: GVector[],
    quantity: number): number {
    let total = 0;
    for (let i = 1; i < quantity; ++i) {
        total += length(sub(eg.computePosition(path[i]),
            eg.computePosition(path[i - 1])));
    }
    return total;
}

describe('EllipsoidGeodesic surface points', () => {
    it('evaluates P(u,v) = (a*cos(u)*sin(v), b*sin(u)*sin(v), c*cos(v))', () => {
        const eg = new EllipsoidGeodesic(3, 2, 1);
        const halfPi = Math.PI / 2;

        // The poles v = 0 and v = pi.
        let p = eg.computePosition(gv(0.7, 0));
        expect(p.values[0]).toBeCloseTo(0, 12);
        expect(p.values[1]).toBeCloseTo(0, 12);
        expect(p.values[2]).toBeCloseTo(1, 12);
        p = eg.computePosition(gv(0.7, Math.PI));
        expect(p.values[2]).toBeCloseTo(-1, 12);

        // The equator v = pi/2 traces the ellipse (a*cos(u), b*sin(u), 0).
        p = eg.computePosition(gv(0, halfPi));
        expect(p.values[0]).toBeCloseTo(3, 12);
        expect(p.values[1]).toBeCloseTo(0, 12);
        expect(p.values[2]).toBeCloseTo(0, 12);
        p = eg.computePosition(gv(halfPi, halfPi));
        expect(p.values[0]).toBeCloseTo(0, 12);
        expect(p.values[1]).toBeCloseTo(2, 12);

        // Every point satisfies the implicit ellipsoid equation.
        for (const [u, v] of [[0.3, 0.4], [2.0, 1.1], [5.5, 2.7], [1.0, 3.0]]) {
            const q = eg.computePosition(gv(u, v));
            const f = (q.values[0] / 3) ** 2 + (q.values[1] / 2) ** 2 +
                (q.values[2] / 1) ** 2;
            expect(f).toBeCloseTo(1, 12);
        }
    });
});

describe('EllipsoidGeodesic on the unit sphere', () => {
    it('reproduces the sphere metric through the segment lengths', () => {
        // For a = b = c = 1 the metric is g = diag(sin(v)^2, 1).
        const eg = new EllipsoidGeodesic(1, 1, 1);
        const halfPi = Math.PI / 2;

        // The equator (v = pi/2) is a great circle, so the segment length
        // along it is the difference of the longitudes.
        expect(eg.computeSegmentLength(gv(0.3, halfPi), gv(1.2, halfPi)))
            .toBeCloseTo(0.9, 8);

        // A meridian (u constant) is also a great circle.
        expect(eg.computeSegmentLength(gv(0.5, 0.4), gv(0.5, 1.1)))
            .toBeCloseTo(0.7, 8);

        // A meridian has zero geodesic curvature; a circle of latitude that
        // is not the equator does not.
        expect(eg.computeSegmentCurvature(gv(0.5, 0.4), gv(0.5, 1.1)))
            .toBeLessThan(1e-8);
        expect(eg.computeSegmentCurvature(gv(0.3, 0.6), gv(1.2, 0.6)))
            .toBeGreaterThan(1e-2);
    });

    it('bends a path toward the great circle', () => {
        const eg = new EllipsoidGeodesic(1, 1, 1);
        eg.subdivisions = 4;
        eg.refinements = 4;
        eg.searchSamples = 8;
        eg.updateDerivedParameters();

        const v0 = 1.0;
        const end0 = gv(0, v0);
        const end1 = gv(1, v0);
        const straightLength = eg.computeTotalLength(2, [end0, end1]);

        const result = eg.computeGeodesic(end0, end1);
        expect(result.quantity).toBe(17);
        const total = eg.computeTotalLength(result.quantity, result.path);

        // The exact great-circle distance between the two points.
        const p0 = eg.computePosition(end0);
        const p1 = eg.computePosition(end1);
        const exact = Math.acos(dot(p0, p1));

        expect(total).toBeLessThan(straightLength);
        expect(Math.abs(total - exact)).toBeLessThan(2e-2);

        // The metric length agrees with the Euclidean length of the 3D
        // polyline through the path points (they differ only by the polyline
        // chord-versus-arc discretization).
        const euclid = euclideanPathLength(eg, result.path, result.quantity);
        expect(euclid).toBeLessThanOrEqual(total + 1e-9);
        expect(Math.abs(total - euclid)).toBeLessThan(2e-2);

        // The endpoints are untouched and the path bows toward the pole
        // (smaller v) because the great circle does.
        expect(result.path[0].values[0]).toBeCloseTo(0, 12);
        expect(result.path[0].values[1]).toBeCloseTo(v0, 12);
        expect(result.path[16].values[0]).toBeCloseTo(1, 12);
        expect(result.path[16].values[1]).toBeCloseTo(v0, 12);
        expect(result.path[8].values[1]).toBeLessThan(v0 - 1e-3);

        // The polyline is monotone in the longitude.
        for (let i = 1; i < result.quantity; ++i) {
            expect(result.path[i].values[0])
                .toBeGreaterThan(result.path[i - 1].values[0]);
        }
    });

    it('leaves a meridian alone because it is already a geodesic', () => {
        const eg = new EllipsoidGeodesic(1, 1, 1);
        eg.subdivisions = 3;
        eg.refinements = 3;
        eg.searchSamples = 8;
        eg.updateDerivedParameters();

        const u0 = 0.8;
        const result = eg.computeGeodesic(gv(u0, 0.6), gv(u0, 1.6));
        expect(result.quantity).toBe(9);
        for (let i = 0; i < result.quantity; ++i) {
            expect(result.path[i].values[0]).toBeCloseTo(u0, 6);
        }
        const total = eg.computeTotalLength(result.quantity, result.path);
        expect(Math.abs(total - 1.0)).toBeLessThan(1e-5);
    });
});

describe('EllipsoidGeodesic on a general ellipsoid', () => {
    it('shortens the path and keeps it on the surface', () => {
        const a = 1.5, b = 1.2, c = 1;
        const eg = new EllipsoidGeodesic(a, b, c);
        eg.subdivisions = 3;
        eg.refinements = 3;
        eg.searchSamples = 8;
        eg.updateDerivedParameters();

        const end0 = gv(0.2, 1.0);
        const end1 = gv(1.4, 1.3);
        const straightLength = eg.computeTotalLength(2, [end0, end1]);
        const result = eg.computeGeodesic(end0, end1);
        const total = eg.computeTotalLength(result.quantity, result.path);

        expect(total).toBeLessThan(straightLength);

        // Every path point lies on the ellipsoid, and the chordal 3D length
        // is a lower bound for the metric length.
        for (let i = 0; i < result.quantity; ++i) {
            const p: Vector = eg.computePosition(result.path[i]);
            const f = (p.values[0] / a) ** 2 + (p.values[1] / b) ** 2 +
                (p.values[2] / c) ** 2;
            expect(f).toBeCloseTo(1, 10);
        }
        const euclid = euclideanPathLength(eg, result.path, result.quantity);
        expect(euclid).toBeLessThanOrEqual(total + 1e-9);

        // The endpoints are preserved exactly.
        expect(result.path[0].values[0]).toBeCloseTo(0.2, 12);
        expect(result.path[0].values[1]).toBeCloseTo(1.0, 12);
        expect(result.path[result.quantity - 1].values[0]).toBeCloseTo(1.4, 12);
        expect(result.path[result.quantity - 1].values[1]).toBeCloseTo(1.3, 12);
    });

    it('scales the segment length with the ellipsoid extents', () => {
        // Scaling all extents by s scales the metric by s^2 and hence every
        // length by s.
        const small = new EllipsoidGeodesic(1, 1, 1);
        const big = new EllipsoidGeodesic(3, 3, 3);
        const p0 = gv(0.4, 0.9), p1 = gv(1.3, 1.7);
        expect(big.computeSegmentLength(p0, p1))
            .toBeCloseTo(3 * small.computeSegmentLength(p0, p1), 10);

        // A prolate ellipsoid (c > a = b): a meridian arc through the pole is
        // longer than the corresponding unit-sphere arc.
        const prolate = new EllipsoidGeodesic(1, 1, 2);
        expect(prolate.computeSegmentLength(gv(0.5, 0.2), gv(0.5, 0.9)))
            .toBeGreaterThan(small.computeSegmentLength(gv(0.5, 0.2),
                gv(0.5, 0.9)));

        // Curvature of a non-geodesic latitude circle is positive on any
        // ellipsoid.
        expect(prolate.computeSegmentCurvature(gv(0.2, 0.7), gv(1.0, 0.7)))
            .toBeGreaterThan(0);
    });
});
