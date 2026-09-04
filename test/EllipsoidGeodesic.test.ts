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

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). EllipsoidGeodesic.h was read
// line by line against src/EllipsoidGeodesic.ts. The class contributes only
// ComputePosition, ComputeMetric and ComputeChristoffel1, so the properties
// below check each of those against a definition that does not reuse the
// implementation: the implicit ellipsoid equation, the Gram matrix of the
// numerically differentiated surface, and the textbook formula
// Gamma_{k,ij} = (dg_{ki}/dx_j + dg_{kj}/dx_i - dg_{ij}/dx_k)/2.
import {
    check, fc, expectClose
} from './helpers/arbitraries.js';

// Exposes the protected pipeline so the tensors can be inspected.
class EllipsoidProbe extends EllipsoidGeodesic {
    prepare(point: GVector): void {
        this.computeMetric(point);
        this.computeChristoffel1(point);
        this.computeMetricInverse();
        this.computeChristoffel2();
    }

    metricAt(point: GVector): number[][] {
        this.computeMetric(point);
        return [[this.mMetric.get(0, 0), this.mMetric.get(0, 1)],
            [this.mMetric.get(1, 0), this.mMetric.get(1, 1)]];
    }

    christoffel1(k: number, i: number, j: number): number {
        return this.mChristoffel1[k].get(i, j);
    }

    christoffel2(k: number, i: number, j: number): number {
        return this.mChristoffel2[k].get(i, j);
    }

    metricInverse(i: number, j: number): number {
        return this.mMetricInverse.get(i, j);
    }
}

// Extents bounded away from zero: the metric of a degenerate ellipsoid is
// singular and the Christoffel symbols of the second kind do not exist.
const extents = fc.tuple(
    fc.double({ min: 0.5, max: 3, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.5, max: 3, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.5, max: 3, noNaN: true, noDefaultInfinity: true }));

// (u, v) with v away from the poles, where dP/du degenerates and the metric
// becomes singular.
const uv = fc.tuple(
    fc.double({ min: -3, max: 3, noNaN: true, noDefaultInfinity: true }),
    fc.double({ min: 0.4, max: Math.PI - 0.4, noNaN: true,
        noDefaultInfinity: true }));

describe('EllipsoidGeodesic verification', () => {
    it('places every parameter point on the ellipsoid', () => {
        check(fc.tuple(extents, uv), ([[a, b, c], [u, v]]) => {
            const eg = new EllipsoidGeodesic(a, b, c);
            const p = eg.computePosition(GVector.fromArray([u, v]));
            expectClose((p.values[0] / a) ** 2 + (p.values[1] / b) ** 2 +
                (p.values[2] / c) ** 2, 1, 1e-12, 1e-12);
        });
    });

    it('metric equals the Gram matrix of the surface derivatives', () => {
        // dP/du and dP/dv by centered differences of ComputePosition, an
        // independent route to g_{ij} = Dot(P_i, P_j).
        check(fc.tuple(extents, uv), ([[a, b, c], [u, v]]) => {
            const probe = new EllipsoidProbe(a, b, c);
            const g = probe.metricAt(GVector.fromArray([u, v]));
            const h = 1e-5;
            const der = (k: number): Vector => {
                const plus = [u, v], minus = [u, v];
                plus[k] += h;
                minus[k] -= h;
                const p = probe.computePosition(GVector.fromArray(plus));
                const m = probe.computePosition(GVector.fromArray(minus));
                return Vector.fromArray([0, 1, 2].map(i =>
                    (p.values[i] - m.values[i]) / (2 * h)));
            };
            const d = [der(0), der(1)];
            for (let i = 0; i < 2; ++i) {
                for (let j = 0; j < 2; ++j) {
                    expectClose(g[i][j], dot(d[i], d[j]), 1e-8, 1e-8);
                }
            }
        });
    });

    it('Christoffel symbols of the first kind satisfy the metric identity', () => {
        // Gamma_{k,ij} = (dg_{ki}/dx_j + dg_{kj}/dx_i - dg_{ij}/dx_k)/2, with
        // the metric derivatives taken by centered differences of the metric
        // that ComputeMetric produces.
        check(fc.tuple(extents, uv), ([[a, b, c], [u, v]]) => {
            const probe = new EllipsoidProbe(a, b, c);
            probe.prepare(GVector.fromArray([u, v]));
            const h = 1e-4;
            const dg = (k: number): number[][] => {
                const plus = [u, v], minus = [u, v];
                plus[k] += h;
                minus[k] -= h;
                const gp = new EllipsoidProbe(a, b, c)
                    .metricAt(GVector.fromArray(plus));
                const gm = new EllipsoidProbe(a, b, c)
                    .metricAt(GVector.fromArray(minus));
                return [[(gp[0][0] - gm[0][0]) / (2 * h),
                    (gp[0][1] - gm[0][1]) / (2 * h)],
                [(gp[1][0] - gm[1][0]) / (2 * h),
                    (gp[1][1] - gm[1][1]) / (2 * h)]];
            };
            const d = [dg(0), dg(1)];
            for (let k = 0; k < 2; ++k) {
                for (let i = 0; i < 2; ++i) {
                    for (let j = 0; j < 2; ++j) {
                        expectClose(probe.christoffel1(k, i, j),
                            0.5 * (d[j][k][i] + d[i][k][j] - d[k][i][j]),
                            1e-5, 1e-5);
                    }
                }
            }
        });
    });

    it('Christoffel symbols of the second kind raise the first index', () => {
        // Gamma^{i2}_{i0 i1} = g^{i2 j} Gamma_{j,i0 i1}, with the inverse of
        // the 2x2 metric written out in closed form rather than reusing the
        // GaussianElimination path the class calls.
        check(fc.tuple(extents, uv), ([[a, b, c], [u, v]]) => {
            const probe = new EllipsoidProbe(a, b, c);
            const g = probe.metricAt(GVector.fromArray([u, v]));
            probe.prepare(GVector.fromArray([u, v]));
            const det = g[0][0] * g[1][1] - g[0][1] * g[1][0];
            const ginv = [[g[1][1] / det, -g[0][1] / det],
                [-g[1][0] / det, g[0][0] / det]];
            for (let i = 0; i < 2; ++i) {
                for (let j = 0; j < 2; ++j) {
                    expectClose(probe.metricInverse(i, j), ginv[i][j],
                        1e-9, 1e-9);
                }
            }
            for (let k = 0; k < 2; ++k) {
                for (let i = 0; i < 2; ++i) {
                    for (let j = 0; j < 2; ++j) {
                        const expected = ginv[k][0] * probe.christoffel1(0, i, j)
                            + ginv[k][1] * probe.christoffel1(1, i, j);
                        expectClose(probe.christoffel2(k, i, j), expected,
                            1e-9, 1e-9);
                    }
                }
            }
        });
    });

    it('scales every length with a uniform scaling of the extents', () => {
        check(fc.tuple(extents, uv, uv,
            fc.double({ min: 0.25, max: 4, noNaN: true,
                noDefaultInfinity: true })),
        ([[a, b, c], [u0, v0], [u1, v1], s]) => {
            const p0 = GVector.fromArray([u0, v0]);
            const p1 = GVector.fromArray([u1, v1]);
            if (Math.abs(u1 - u0) + Math.abs(v1 - v0) < 1e-2) { return; }
            const base = new EllipsoidGeodesic(a, b, c);
            const scaled = new EllipsoidGeodesic(s * a, s * b, s * c);
            expectClose(scaled.computeSegmentLength(p0, p1),
                s * base.computeSegmentLength(p0, p1), 1e-10, 1e-10);
        });
    });

    it('gives exact arc lengths along the great circles of a unit sphere', () => {
        // On the unit sphere the metric is diag(sin(v)^2, 1). A meridian
        // segment has integrand |dv| and an equator segment has integrand
        // |du|, both constant, so the trapezoid rule is exact.
        check(fc.tuple(uv, fc.double({ min: -1, max: 1, noNaN: true,
            noDefaultInfinity: true })), ([[u, v], d]) => {
            if (Math.abs(d) < 1e-3) { return; }
            const eg = new EllipsoidGeodesic(1, 1, 1);
            const vEnd = Math.min(Math.max(v + d, 0.05), Math.PI - 0.05);
            expectClose(eg.computeSegmentLength(GVector.fromArray([u, v]),
                GVector.fromArray([u, vEnd])), Math.abs(vEnd - v),
            1e-12, 1e-12);
            const halfPi = Math.PI / 2;
            expectClose(eg.computeSegmentLength(
                GVector.fromArray([u, halfPi]),
                GVector.fromArray([u + d, halfPi])), Math.abs(d),
            1e-12, 1e-12);
        });
    });

    it('bounds the metric length below by the great-circle distance', () => {
        // The metric length of the straight parameter segment is the length of
        // a curve on the unit sphere joining the two surface points, so it is
        // at least the geodesic (great-circle) distance between them.
        check(fc.tuple(uv, uv), ([[u0, v0], [u1, v1]]) => {
            if (Math.abs(u1 - u0) + Math.abs(v1 - v0) < 1e-2) { return; }
            const eg = new EllipsoidGeodesic(1, 1, 1);
            const p0 = GVector.fromArray([u0, v0]);
            const p1 = GVector.fromArray([u1, v1]);
            const q0 = eg.computePosition(p0);
            const q1 = eg.computePosition(p1);
            const great = Math.acos(Math.min(1, Math.max(-1, dot(q0, q1))));
            // The trapezoid rule under-resolves a strongly curved integrand,
            // so allow the quadrature error of 16 samples.
            expect(eg.computeSegmentLength(p0, p1))
                .toBeGreaterThan(great - 1e-2);
        });
    });
});
