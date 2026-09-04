import { describe, expect, it } from 'vitest';
import { GVector } from '../src/GVector.js';
import { RiemannianGeodesic } from '../src/RiemannianGeodesic.js';
import { Vector, sub, length as vectorLength } from '../src/Vector.js';

// A flat (Euclidean) metric: g = I and all Christoffel symbols vanish, so the
// geodesics are straight lines and the arc length is the Euclidean length.
class FlatGeodesic extends RiemannianGeodesic {
    constructor(dimension: number) {
        super(dimension);
    }

    protected computeMetric(_point: GVector): void {
        for (let i = 0; i < this.mDimension; ++i) {
            for (let j = 0; j < this.mDimension; ++j) {
                this.mMetric.set(i, j, i === j ? 1 : 0);
            }
        }
    }

    protected computeChristoffel1(_point: GVector): void {
        for (let k = 0; k < this.mDimension; ++k) {
            this.mChristoffel1[k].makeZero();
        }
    }

    // Test-only accessors for the protected machinery.
    publicComputeIntegrand(pos: GVector, der: Vector): number {
        return this.computeIntegrand(pos, der);
    }

    publicComputeMetricInverse(): boolean {
        this.computeMetric(new GVector(this.mDimension));
        return this.computeMetricInverse();
    }
}

// A non-Euclidean metric of the form g = diag(s(x)^2, s(x)^2, ...) is not
// needed here; instead use the unit sphere parameterized by
//   P(u,v) = (cos(u)*sin(v), sin(u)*sin(v), cos(v)).
// The metric is g = diag(sin(v)^2, 1) and the Christoffel symbols of the
// first kind are Gamma_{k,ij} = Dot(P_ij, P_k), which is how upstream's
// EllipsoidGeodesic computes them.
class SphereGeodesic extends RiemannianGeodesic {
    private mDer0: number[];
    private mDer1: number[];
    private mCos0 = 0;
    private mSin0 = 0;
    private mCos1 = 0;
    private mSin1 = 0;

    constructor() {
        super(2);
        this.mDer0 = [0, 0, 0];
        this.mDer1 = [0, 0, 0];
    }

    static position(u: number, v: number): number[] {
        return [Math.cos(u) * Math.sin(v), Math.sin(u) * Math.sin(v),
            Math.cos(v)];
    }

    protected computeMetric(point: GVector): void {
        this.mCos0 = Math.cos(point.values[0]);
        this.mSin0 = Math.sin(point.values[0]);
        this.mCos1 = Math.cos(point.values[1]);
        this.mSin1 = Math.sin(point.values[1]);

        this.mDer0 = [-this.mSin0 * this.mSin1, this.mCos0 * this.mSin1, 0];
        this.mDer1 = [this.mCos0 * this.mCos1, this.mSin0 * this.mCos1,
            -this.mSin1];

        this.mMetric.set(0, 0, dot3(this.mDer0, this.mDer0));
        this.mMetric.set(0, 1, dot3(this.mDer0, this.mDer1));
        this.mMetric.set(1, 0, this.mMetric.get(0, 1));
        this.mMetric.set(1, 1, dot3(this.mDer1, this.mDer1));
    }

    protected computeChristoffel1(_point: GVector): void {
        const der00 = [-this.mCos0 * this.mSin1, -this.mSin0 * this.mSin1, 0];
        const der01 = [-this.mSin0 * this.mCos1, this.mCos0 * this.mCos1, 0];
        const der11 = [-this.mCos0 * this.mSin1, -this.mSin0 * this.mSin1,
            -this.mCos1];

        this.mChristoffel1[0].set(0, 0, dot3(der00, this.mDer0));
        this.mChristoffel1[0].set(0, 1, dot3(der01, this.mDer0));
        this.mChristoffel1[0].set(1, 0, this.mChristoffel1[0].get(0, 1));
        this.mChristoffel1[0].set(1, 1, dot3(der11, this.mDer0));

        this.mChristoffel1[1].set(0, 0, dot3(der00, this.mDer1));
        this.mChristoffel1[1].set(0, 1, dot3(der01, this.mDer1));
        this.mChristoffel1[1].set(1, 0, this.mChristoffel1[1].get(0, 1));
        this.mChristoffel1[1].set(1, 1, dot3(der11, this.mDer1));
    }
}

function dot3(a: number[], b: number[]): number {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

function gv(...values: number[]): GVector {
    return GVector.fromArray(values);
}

describe('RiemannianGeodesic', () => {
    it('requires a dimension of at least two', () => {
        expect(() => new FlatGeodesic(1)).toThrow();
        expect(new FlatGeodesic(2).getDimension()).toBe(2);
        expect(new FlatGeodesic(4).getDimension()).toBe(4);
    });

    it('exposes the documented default tweakable parameters', () => {
        const g = new FlatGeodesic(2);
        expect(g.integralSamples).toBe(16);
        expect(g.searchSamples).toBe(32);
        expect(g.derivativeStep).toBeCloseTo(1e-4, 12);
        expect(g.subdivisions).toBe(7);
        expect(g.refinements).toBe(8);
        expect(g.searchRadius).toBe(1);
        expect(g.getSubdivisionStep()).toBe(0);
        expect(g.getRefinementStep()).toBe(0);
        expect(g.getCurrentQuantity()).toBe(0);
    });

    it('computes Euclidean lengths for the flat metric', () => {
        const g = new FlatGeodesic(3);
        expect(g.computeSegmentLength(gv(0, 0, 0), gv(3, 4, 0)))
            .toBeCloseTo(5, 10);
        expect(g.computeSegmentLength(gv(1, -2, 5), gv(1, -2, 8)))
            .toBeCloseTo(3, 10);

        const path = [gv(0, 0, 0), gv(1, 0, 0), gv(1, 2, 0), gv(1, 2, 2)];
        expect(g.computeTotalLength(4, path)).toBeCloseTo(5, 10);
        expect(() => g.computeTotalLength(1, path)).toThrow();

        // A degenerate segment has a zero quadratic form, which upstream
        // treats as an unexpected condition.
        expect(() => g.computeSegmentLength(gv(1, 1, 1), gv(1, 1, 1)))
            .toThrow();
    });

    it('produces straight-line geodesics for the flat metric', () => {
        const g = new FlatGeodesic(2);
        g.subdivisions = 3;
        g.refinements = 2;
        g.searchSamples = 8;
        g.updateDerivedParameters();

        const end0 = gv(0, 0);
        const end1 = gv(2, 1);
        const result = g.computeGeodesic(end0, end1);
        expect(result.quantity).toBe(9);
        expect(result.path.length).toBe(9);

        for (let i = 0; i < result.quantity; ++i) {
            const s = i / (result.quantity - 1);
            expect(result.path[i].values[0]).toBeCloseTo(2 * s, 8);
            expect(result.path[i].values[1]).toBeCloseTo(s, 8);
        }

        const total = g.computeTotalLength(result.quantity, result.path);
        expect(total).toBeCloseTo(Math.sqrt(5), 8);

        // The progress counters are reset when the computation ends.
        expect(g.getSubdivisionStep()).toBe(0);
        expect(g.getRefinementStep()).toBe(0);
        expect(g.getCurrentQuantity()).toBe(0);
    });

    it('has zero curvature for the flat metric', () => {
        const g = new FlatGeodesic(2);
        expect(g.computeSegmentCurvature(gv(0, 0), gv(1, 2)))
            .toBeCloseTo(0, 12);
        const path = [gv(0, 0), gv(1, 1), gv(3, 0)];
        expect(g.computeTotalCurvature(3, path)).toBeCloseTo(0, 12);
        expect(() => g.computeTotalCurvature(0, path)).toThrow();
        expect(g.publicComputeIntegrand(gv(0, 0), Vector.fromArray([1, 2])))
            .toBeCloseTo(0, 12);
        expect(g.publicComputeMetricInverse()).toBe(true);
    });

    it('invokes the refine callback and reports progress', () => {
        const g = new FlatGeodesic(2);
        g.subdivisions = 2;
        g.refinements = 1;
        g.searchSamples = 4;
        g.updateDerivedParameters();

        const steps: number[] = [];
        g.refineCallback = () => {
            steps.push(g.getSubdivisionStep());
        };
        g.computeGeodesic(gv(0, 0), gv(1, 1));
        // Subdivide temporarily disables the callback, so only the explicit
        // refinement passes are reported: 1 vertex after the first
        // subdivision and 3 after the second.
        expect(steps).toEqual([1, 2, 2, 2]);
    });

    it('reproduces the sphere metric and its great-circle geodesics', () => {
        const g = new SphereGeodesic();

        // The equator (v = pi/2) is a great circle, so the segment length
        // along it is the difference of the longitudes.
        const halfPi = Math.PI / 2;
        expect(g.computeSegmentLength(gv(0.3, halfPi), gv(1.2, halfPi)))
            .toBeCloseTo(0.9, 8);

        // A meridian (u constant) is also a great circle.
        expect(g.computeSegmentLength(gv(0.5, 0.4), gv(0.5, 1.1)))
            .toBeCloseTo(0.7, 8);

        // A meridian has zero geodesic curvature; a circle of latitude that
        // is not the equator does not.
        expect(g.computeSegmentCurvature(gv(0.5, 0.4), gv(0.5, 1.1)))
            .toBeLessThan(1e-8);
        expect(g.computeSegmentCurvature(gv(0.3, 0.6), gv(1.2, 0.6)))
            .toBeGreaterThan(1e-2);
    });

    it('bends a sphere geodesic toward the great circle', () => {
        const g = new SphereGeodesic();
        g.subdivisions = 4;
        g.refinements = 4;
        g.searchSamples = 8;
        g.updateDerivedParameters();

        const v0 = 1.0;
        const end0 = gv(0, v0);
        const end1 = gv(1, v0);
        const straightLength = g.computeTotalLength(2, [end0, end1]);

        const result = g.computeGeodesic(end0, end1);
        expect(result.quantity).toBe(17);
        const total = g.computeTotalLength(result.quantity, result.path);

        // The exact great-circle distance between the two points.
        const p0 = SphereGeodesic.position(0, v0);
        const p1 = SphereGeodesic.position(1, v0);
        const exact = Math.acos(dot3(p0, p1));
        expect(exact).toBeGreaterThan(0);

        expect(total).toBeLessThan(straightLength);
        expect(Math.abs(total - exact)).toBeLessThan(2e-2);

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

    it('refines a midpoint monotonically', () => {
        const g = new SphereGeodesic();
        g.searchSamples = 16;
        g.updateDerivedParameters();

        const end0 = gv(0, 1.0);
        const end1 = gv(1, 1.0);
        const sub = g.subdivide(end0, end1);
        const naiveMid = gv(0.5, 1.0);
        const naiveLength = g.computeSegmentLength(end0, naiveMid) +
            g.computeSegmentLength(naiveMid, end1);
        const refinedLength = g.computeSegmentLength(end0, sub.mid) +
            g.computeSegmentLength(sub.mid, end1);
        expect(sub.changed).toBe(true);
        expect(refinedLength).toBeLessThan(naiveLength);
        expect(sub.mid.values[1]).toBeLessThan(1.0);

        // refine does not modify its input midpoint (upstream overwrites the
        // reference parameter; the port returns a new vector).
        const mid = gv(0.5, 1.0);
        const refined = g.refine(end0, mid, end1);
        expect(mid.values[0]).toBeCloseTo(0.5, 12);
        expect(mid.values[1]).toBeCloseTo(1.0, 12);
        expect(refined.mid).not.toBe(mid);
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md). RiemannianGeodesic.h was read
// line by line against src/RiemannianGeodesic.ts. The properties below cross
// check the numerical machinery against closed forms: the trapezoid quadrature
// against its exactness on affine integrands (which pins mIntegralStep, the
// 0.5 endpoint weights and the imax bound), the Christoffel symbols of the
// second kind against the textbook values for the unit sphere, and the
// multiresolution solver against the straight lines of a flat metric.
import {
    check, fc, expectClose, wellScaledVector
} from './helpers/arbitraries.js';

// g(x) = f(x)^2 * I with f affine. The length integrand along a segment is
// then f(point0 + t*diff)*|diff|, an affine function of t, and the composite
// trapezoid rule is exact for affine integrands. So the segment length must
// equal |diff| * (f(point0) + f(point1)) / 2 no matter how many integral
// samples are used - an independent closed form for ComputeSegmentLength.
class ConformalGeodesic extends RiemannianGeodesic {
    constructor(dimension: number, private readonly c0: number,
        private readonly c: number[]) {
        super(dimension);
    }

    f(point: GVector): number {
        let value = this.c0;
        for (let i = 0; i < this.mDimension; ++i) {
            value += this.c[i] * point.values[i];
        }
        return value;
    }

    protected computeMetric(point: GVector): void {
        const s = this.f(point) ** 2;
        for (let i = 0; i < this.mDimension; ++i) {
            for (let j = 0; j < this.mDimension; ++j) {
                this.mMetric.set(i, j, i === j ? s : 0);
            }
        }
    }

    protected computeChristoffel1(_point: GVector): void {
        for (let k = 0; k < this.mDimension; ++k) {
            this.mChristoffel1[k].makeZero();
        }
    }
}

// The unit sphere P(u,v) = (cos u sin v, sin u sin v, cos v) with the
// protected machinery exposed for inspection.
class SphereProbe extends RiemannianGeodesic {
    private d0: number[] = [0, 0, 0];
    private d1: number[] = [0, 0, 0];

    constructor() {
        super(2);
    }

    protected computeMetric(point: GVector): void {
        const cu = Math.cos(point.values[0]), su = Math.sin(point.values[0]);
        const cv = Math.cos(point.values[1]), sv = Math.sin(point.values[1]);
        this.d0 = [-su * sv, cu * sv, 0];
        this.d1 = [cu * cv, su * cv, -sv];
        this.mMetric.set(0, 0, dot3(this.d0, this.d0));
        this.mMetric.set(0, 1, dot3(this.d0, this.d1));
        this.mMetric.set(1, 0, this.mMetric.get(0, 1));
        this.mMetric.set(1, 1, dot3(this.d1, this.d1));
    }

    protected computeChristoffel1(point: GVector): void {
        const cu = Math.cos(point.values[0]), su = Math.sin(point.values[0]);
        const cv = Math.cos(point.values[1]), sv = Math.sin(point.values[1]);
        const d00 = [-cu * sv, -su * sv, 0];
        const d01 = [-su * cv, cu * cv, 0];
        const d11 = [-cu * sv, -su * sv, -cv];
        const ders = [this.d0, this.d1];
        for (let k = 0; k < 2; ++k) {
            this.mChristoffel1[k].set(0, 0, dot3(d00, ders[k]));
            this.mChristoffel1[k].set(0, 1, dot3(d01, ders[k]));
            this.mChristoffel1[k].set(1, 0, dot3(d01, ders[k]));
            this.mChristoffel1[k].set(1, 1, dot3(d11, ders[k]));
        }
    }

    // Run the full protected pipeline at a point and expose the results.
    prepare(point: GVector): void {
        this.computeMetric(point);
        this.computeChristoffel1(point);
        this.computeMetricInverse();
        this.computeChristoffel2();
        this.computeMetricDerivative();
    }

    metric(i: number, j: number): number { return this.mMetric.get(i, j); }
    christoffel1(k: number, i: number, j: number): number {
        return this.mChristoffel1[k].get(i, j);
    }
    christoffel2(k: number, i: number, j: number): number {
        return this.mChristoffel2[k].get(i, j);
    }
    metricDerivative(k: number, i: number, j: number): number {
        return this.mMetricDerivative[k].get(i, j);
    }
}

const dimension = fc.integer({ min: 2, max: 5 });
const angle = fc.double({ min: 0.35, max: Math.PI - 0.35, noNaN: true,
    noDefaultInfinity: true });

describe('RiemannianGeodesic verification', () => {
    it('reproduces Euclidean lengths exactly for the flat metric', () => {
        check(dimension.chain(n => fc.tuple(wellScaledVector(n, -8, 8),
            wellScaledVector(n, -8, 8))
            .filter(([a, b]) => vectorLength(sub(b, a)) > 1e-2)),
        ([a, b]) => {
            const g = new FlatGeodesic(a.size);
            const p0 = GVector.fromArray(a.values);
            const p1 = GVector.fromArray(b.values);
            expectClose(g.computeSegmentLength(p0, p1),
                vectorLength(sub(b, a)), 1e-12, 1e-12);
            // The quadrature is symmetric in its endpoints.
            expectClose(g.computeSegmentLength(p1, p0),
                g.computeSegmentLength(p0, p1), 1e-12, 1e-12);
        });
    });

    it('integrates affine integrands exactly for any sample count', () => {
        // The trapezoid rule with mIntegralStep = 1/(integralSamples-1), half
        // weights at the endpoints and interior samples i = 1..integralSamples-2
        // is exact for affine integrands; an off-by-one in the loop bound or a
        // wrong step would break this identity.
        check(fc.tuple(fc.integer({ min: 2, max: 4 }),
            fc.integer({ min: 2, max: 24 }),
            fc.double({ min: 2, max: 6, noNaN: true, noDefaultInfinity: true }),
            fc.array(fc.double({ min: -0.2, max: 0.2, noNaN: true,
                noDefaultInfinity: true }), { minLength: 4, maxLength: 4 }),
            fc.array(fc.double({ min: -1, max: 1, noNaN: true,
                noDefaultInfinity: true }), { minLength: 8, maxLength: 8 })),
        ([n, samples, c0, c, coords]) => {
            const g = new ConformalGeodesic(n, c0, c.slice(0, n));
            g.integralSamples = samples;
            g.updateDerivedParameters();
            const p0 = GVector.fromArray(coords.slice(0, n));
            const p1 = GVector.fromArray(coords.slice(4, 4 + n));
            const diff = sub(p1, p0);
            if (vectorLength(diff) < 1e-2) { return; }
            // f >= c0 - 0.2*4 > 0 on the whole segment by construction.
            const expected = vectorLength(diff) *
                0.5 * (g.f(p0) + g.f(p1));
            expectClose(g.computeSegmentLength(p0, p1), expected, 1e-12, 1e-12);
        });
    });

    it('adds segment lengths to the total length of a path', () => {
        check(fc.tuple(dimension, fc.integer({ min: 2, max: 5 })).chain(
            ([n, q]) => fc.tuple(fc.constant(n), fc.array(
                wellScaledVector(n, -5, 5), { minLength: q, maxLength: q }))),
        ([n, pts]) => {
            const g = new FlatGeodesic(n);
            const path = pts.map(p => GVector.fromArray(p.values));
            for (let i = 1; i < path.length; ++i) {
                if (vectorLength(sub(path[i], path[i - 1])) < 1e-2) { return; }
            }
            let expected = 0;
            for (let i = 1; i < path.length; ++i) {
                expected += g.computeSegmentLength(path[i - 1], path[i]);
            }
            expectClose(g.computeTotalLength(path.length, path), expected,
                1e-12, 1e-12);
        });
    });

    it('produces the uniformly sampled straight line for a flat metric', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 3 }), fc.integer({ min: 1, max: 3 })).chain(
            ([n, s]) => fc.tuple(fc.constant(s), wellScaledVector(n, -5, 5),
                wellScaledVector(n, -5, 5))),
        ([subdivisions, a, b]) => {
            if (vectorLength(sub(b, a)) < 1) { return; }
            const g = new FlatGeodesic(a.size);
            g.subdivisions = subdivisions;
            g.refinements = 1;
            g.searchSamples = 4;
            g.updateDerivedParameters();
            const result = g.computeGeodesic(GVector.fromArray(a.values),
                GVector.fromArray(b.values));
            expect(result.quantity).toBe((1 << subdivisions) + 1);
            expect(result.path.length).toBe(result.quantity);
            for (let i = 0; i < result.quantity; ++i) {
                const s = i / (result.quantity - 1);
                for (let k = 0; k < a.size; ++k) {
                    expectClose(result.path[i].values[k],
                        a.values[k] + s * (b.values[k] - a.values[k]),
                        1e-6, 1e-6);
                }
            }
        }, 30);
    });

    it('leaves the midpoint of a flat-metric segment where it is', () => {
        // F(m) = |m-e0| + |m-e1| attains its minimum on the whole segment, so
        // the estimated gradient at the midpoint vanishes and the line search
        // cannot improve on it.
        check(dimension.chain(n => fc.tuple(wellScaledVector(n, -5, 5),
            wellScaledVector(n, -5, 5))),
        ([a, b]) => {
            if (vectorLength(sub(b, a)) < 1) { return; }
            const g = new FlatGeodesic(a.size);
            g.searchSamples = 8;
            g.updateDerivedParameters();
            const e0 = GVector.fromArray(a.values);
            const e1 = GVector.fromArray(b.values);
            const result = g.subdivide(e0, e1);
            for (let k = 0; k < a.size; ++k) {
                expectClose(result.mid.values[k],
                    0.5 * (a.values[k] + b.values[k]), 1e-6, 1e-6);
            }
            expect(g.computeSegmentLength(e0, result.mid) +
                g.computeSegmentLength(result.mid, e1))
                .toBeLessThanOrEqual(vectorLength(sub(b, a)) * (1 + 1e-9));
        });
    });

    it('matches the textbook Christoffel symbols of the unit sphere', () => {
        // With x = (u,v) and g = diag(sin(v)^2, 1) the nonzero symbols are
        // Gamma_{0,01} = Gamma_{0,10} = sin(v)cos(v), Gamma_{1,00} =
        // -sin(v)cos(v), and of the second kind Gamma^0_{01} = cot(v) and
        // Gamma^1_{00} = -sin(v)cos(v).
        check(fc.tuple(fc.double({ min: -3, max: 3, noNaN: true,
            noDefaultInfinity: true }), angle), ([u, v]) => {
            const probe = new SphereProbe();
            probe.prepare(GVector.fromArray([u, v]));
            const sc = Math.sin(v) * Math.cos(v);
            expectClose(probe.metric(0, 0), Math.sin(v) ** 2, 1e-12, 1e-12);
            expectClose(probe.metric(0, 1), 0, 1e-12, 1e-12);
            expectClose(probe.metric(1, 1), 1, 1e-12, 1e-12);

            expectClose(probe.christoffel1(0, 0, 0), 0, 1e-12, 1e-12);
            expectClose(probe.christoffel1(0, 0, 1), sc, 1e-12, 1e-12);
            expectClose(probe.christoffel1(0, 1, 0), sc, 1e-12, 1e-12);
            expectClose(probe.christoffel1(0, 1, 1), 0, 1e-12, 1e-12);
            expectClose(probe.christoffel1(1, 0, 0), -sc, 1e-12, 1e-12);
            expectClose(probe.christoffel1(1, 0, 1), 0, 1e-12, 1e-12);
            expectClose(probe.christoffel1(1, 1, 1), 0, 1e-12, 1e-12);

            expectClose(probe.christoffel2(0, 0, 1), Math.cos(v) / Math.sin(v),
                1e-9, 1e-9);
            expectClose(probe.christoffel2(0, 1, 0), Math.cos(v) / Math.sin(v),
                1e-9, 1e-9);
            expectClose(probe.christoffel2(0, 0, 0), 0, 1e-9, 1e-9);
            expectClose(probe.christoffel2(1, 0, 0), -sc, 1e-9, 1e-9);
            expectClose(probe.christoffel2(1, 0, 1), 0, 1e-9, 1e-9);
            expectClose(probe.christoffel2(1, 1, 1), 0, 1e-9, 1e-9);
        });
    });

    it('computes twice the Christoffel symbol, not the metric derivative', () => {
        // Upstream issue #295 item 2, preserved: ComputeMetricDerivative sets
        // mMetricDerivative[d](i,j) = Gamma_{d,ij} + Gamma_{d,ji}, which under
        // the Gamma_{k,ij} = Dot(P_ij, P_k) convention is 2*Gamma_{d,ij} and
        // not dg_{ij}/dx_d = Gamma_{i,jd} + Gamma_{j,id}. The quantity is never
        // read inside GTE, so the port keeps it as upstream computes it.
        check(fc.tuple(fc.double({ min: -3, max: 3, noNaN: true,
            noDefaultInfinity: true }), angle), ([u, v]) => {
            const probe = new SphereProbe();
            const point = GVector.fromArray([u, v]);
            probe.prepare(point);
            for (let d = 0; d < 2; ++d) {
                for (let i = 0; i < 2; ++i) {
                    for (let j = 0; j < 2; ++j) {
                        expectClose(probe.metricDerivative(d, i, j),
                            2 * probe.christoffel1(d, i, j), 1e-12, 1e-12);
                    }
                }
            }

            // The true derivative of the metric, by central differences, is
            // Gamma_{i,jd} + Gamma_{j,id}; the two disagree for the sphere.
            const h = 1e-5;
            for (let d = 0; d < 2; ++d) {
                for (let i = 0; i < 2; ++i) {
                    for (let j = 0; j < 2; ++j) {
                        const plus = new SphereProbe();
                        const minus = new SphereProbe();
                        const pp = GVector.fromArray(point.values.slice());
                        const pm = GVector.fromArray(point.values.slice());
                        pp.values[d] += h;
                        pm.values[d] -= h;
                        plus.prepare(pp);
                        minus.prepare(pm);
                        const fd = (plus.metric(i, j) - minus.metric(i, j)) /
                            (2 * h);
                        expectClose(fd, probe.christoffel1(i, j, d) +
                            probe.christoffel1(j, i, d), 1e-6, 1e-6);
                    }
                }
            }
            // sin(v)cos(v) is nonzero away from the equator, where the two
            // quantities differ: dg_{00}/dx_1 = 2 sin v cos v but
            // mMetricDerivative[1](0,0) = 2*Gamma_{1,00} = -2 sin v cos v.
            if (Math.abs(Math.cos(v)) > 0.1) {
                expect(Math.abs(probe.metricDerivative(1, 0, 0) -
                    2 * Math.sin(v) * Math.cos(v))).toBeGreaterThan(1e-3);
            }
        });
    });

    it('reports zero curvature for straight lines in a flat metric', () => {
        check(dimension.chain(n => fc.tuple(wellScaledVector(n, -5, 5),
            wellScaledVector(n, -5, 5))), ([a, b]) => {
            if (vectorLength(sub(b, a)) < 1) { return; }
            const g = new FlatGeodesic(a.size);
            expect(g.computeSegmentCurvature(GVector.fromArray(a.values),
                GVector.fromArray(b.values))).toBeLessThan(1e-12);
        });
    });
});
