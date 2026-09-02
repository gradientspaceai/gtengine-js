import { describe, expect, it } from 'vitest';
import { GVector } from '../src/GVector';
import { RiemannianGeodesic } from '../src/RiemannianGeodesic';
import { Vector, dot } from '../src/Vector';

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
