// Replays oracle/cpp/cases/v18-curves.cpp (verify group 18, curves). Keep the
// two files in the same order.
//
// NaturalSplineCurve is arithmetic only and is compared bit for bit. The
// geodesic algorithm of RiemannianGeodesic is compared bit for bit through
// PolySurfaceGeodesic below, a concrete subclass over a polynomial graph that
// mirrors the C++ one, so that the steepest-descent search, the Trapezoid Rule
// integrators, the metric inverse and the Christoffel symbols of the second
// kind are exercised with no call to the C math library. EllipsoidGeodesic's
// own cases call sin and cos and carry a tolerance; see the C++ file for how
// the libm-decided control flow is kept out of them.
import { describe } from 'vitest';
import { EllipsoidGeodesic } from '../../src/EllipsoidGeodesic.js';
import { GMatrix } from '../../src/GMatrix.js';
import { GVector } from '../../src/GVector.js';
import { NaturalSplineCurve } from '../../src/NaturalSplineCurve.js';
import { RiemannianGeodesic } from '../../src/RiemannianGeodesic.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- NaturalSplineCurve ---------------------------------------------------

interface SplineInput {
    points: Vector[];
    times: number[];
}

// Mirrors MakeSpline: the point count, then that many 3D points, then that
// many times.
function splineInput(io: OracleIO): SplineInput {
    const n = io.integer();
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) { points.push(io.vec(3)); }
    const times = io.reals(n);
    return { points, times };
}

function emitSplineJets(io: OracleIO, curve: NaturalSplineCurve,
    times: readonly number[], order: number): void {
    for (const t of times) {
        const jet: Vector[] = [];
        for (let i = 0; i <= order; ++i) { jet.push(new Vector(3)); }
        curve.evaluate(t, order, jet);
        for (let i = 0; i <= order; ++i) { io.outVecExact(jet[i]); }
    }
}

// ---- the polynomial-surface manifold, mirroring PolySurfaceGeodesic -------
// The graph z = u^2*v + u*v^2. The metric is the Gram matrix of the tangent
// vectors and the Christoffel symbols of the first kind are Dot(P_ij, P_k).

class PolySurfaceGeodesic extends RiemannianGeodesic {
    constructor(dimension: number) {
        super(dimension);
    }

    // Public access to the protected helpers, in the order computeIntegrand
    // calls them, plus computeMetricDerivative.
    probeTensors(point: GVector): boolean {
        this.computeMetric(point);
        this.computeChristoffel1(point);
        const invertible = this.computeMetricInverse();
        this.computeChristoffel2();
        this.computeMetricDerivative();
        return invertible;
    }

    probeIntegrand(pos: GVector, der: GVector): number {
        return this.computeIntegrand(pos, der);
    }

    get metric(): GMatrix { return this.mMetric; }
    get metricInverse(): GMatrix { return this.mMetricInverse; }
    christoffel1(i: number): GMatrix { return this.mChristoffel1[i]; }
    christoffel2(i: number): GMatrix { return this.mChristoffel2[i]; }
    metricDerivative(i: number): GMatrix { return this.mMetricDerivative[i]; }

    protected override computeMetric(point: GVector): void {
        const u = point.values[0];
        const v = point.values[1];
        const du = 2.0 * u * v + v * v;
        const dv = u * u + 2.0 * u * v;
        this.mMetric.set(0, 0, 1.0 + du * du);
        this.mMetric.set(0, 1, du * dv);
        this.mMetric.set(1, 0, this.mMetric.get(0, 1));
        this.mMetric.set(1, 1, 1.0 + dv * dv);
    }

    protected override computeChristoffel1(point: GVector): void {
        const u = point.values[0];
        const v = point.values[1];
        const du = 2.0 * u * v + v * v;
        const dv = u * u + 2.0 * u * v;
        const duu = 2.0 * v;
        const duv = 2.0 * u + 2.0 * v;
        const dvv = 2.0 * u;
        this.mChristoffel1[0].set(0, 0, duu * du);
        this.mChristoffel1[0].set(0, 1, duv * du);
        this.mChristoffel1[0].set(1, 0, this.mChristoffel1[0].get(0, 1));
        this.mChristoffel1[0].set(1, 1, dvv * du);
        this.mChristoffel1[1].set(0, 0, duu * dv);
        this.mChristoffel1[1].set(0, 1, duv * dv);
        this.mChristoffel1[1].set(1, 0, this.mChristoffel1[1].get(0, 1));
        this.mChristoffel1[1].set(1, 1, dvv * dv);
    }
}

function gvec2(io: OracleIO): GVector {
    const x = io.real();
    const y = io.real();
    return GVector.fromArray([x, y]);
}

function gpath(io: OracleIO, n: number): GVector[] {
    const path: GVector[] = [];
    for (let i = 0; i < n; ++i) { path.push(gvec2(io)); }
    return path;
}

function emitMatrix2(io: OracleIO, m: GMatrix): void {
    for (let r = 0; r < 2; ++r) {
        for (let c = 0; c < 2; ++c) { io.outRealExact(m.get(r, c)); }
    }
}

describe('oracle: v18-curves', () => {
    const family = new OracleFamily('v18-curves');

    family.case('NaturalSplineCurve.free.evaluate', (io) => {
        const s = splineInput(io);
        const t = io.reals(5);
        const order = io.integer();
        const curve = NaturalSplineCurve.createFree(s.points, s.times);
        io.outInt(curve.getNumPoints());
        emitSplineJets(io, curve, t, order);
    }, { exact: true });

    family.case('NaturalSplineCurve.closed.evaluate', (io) => {
        const s = splineInput(io);
        const t = io.reals(5);
        const order = io.integer();
        const curve = NaturalSplineCurve.createClosed(s.points, s.times);
        io.outInt(curve.getNumPoints());
        emitSplineJets(io, curve, t, order);
    }, { exact: true });

    family.case('NaturalSplineCurve.closed.wrapRow.deviation', (io) => {
        const s = splineInput(io);
        const t = io.reals(5);
        const order = io.integer();
        const curve = NaturalSplineCurve.createClosed(s.points, s.times);
        io.outInt(curve.getNumPoints());
        emitSplineJets(io, curve, t, order);
    }, { deviation: 'issue #295: NaturalSplineCurve::CreateClosed writes the '
        + 'wrap-around row with three plain assignments whose columns coincide '
        + 'for numPoints 2 and 3' });

    family.case('NaturalSplineCurve.clamped.evaluate', (io) => {
        const s = splineInput(io);
        const t = io.reals(5);
        const order = io.integer();
        const derivative0 = io.vec(3);
        const derivative1 = io.vec(3);
        const curve = NaturalSplineCurve.createClamped(s.points, s.times,
            derivative0, derivative1);
        io.outInt(curve.getNumPoints());
        emitSplineJets(io, curve, t, order);
    }, { exact: true });

    family.case('NaturalSplineCurve.twoPoints.evaluate', (io) => {
        const s = splineInput(io);
        const t = io.reals(5);
        const order = io.integer();
        const isFree = io.boolean();
        const derivative0 = io.vec(3);
        const derivative1 = io.vec(3);
        const curve = isFree
            ? NaturalSplineCurve.createFree(s.points, s.times)
            : NaturalSplineCurve.createClamped(s.points, s.times, derivative0,
                derivative1);
        io.outInt(curve.getNumPoints());
        emitSplineJets(io, curve, t, order);
    }, { exact: true });

    family.case('NaturalSplineCurve.construct.tooFewPoints', (io) => {
        const point = io.vec(3);
        const time = io.real();
        const curve = NaturalSplineCurve.createFree([point], [time]);
        io.outInt(curve.getNumPoints());
    }, { exact: true });

    // ---- RiemannianGeodesic on the polynomial surface ---------------------
    // Polynomial metric and Christoffel symbols, so the whole algorithm is
    // exact IEEE arithmetic, including the comparisons that decide its
    // branches.

    family.case('RiemannianGeodesic.poly.tensors', (io) => {
        const point = gvec2(io);
        const manifold = new PolySurfaceGeodesic(2);
        io.outBool(manifold.probeTensors(point));
        emitMatrix2(io, manifold.metric);
        emitMatrix2(io, manifold.metricInverse);
        emitMatrix2(io, manifold.christoffel1(0));
        emitMatrix2(io, manifold.christoffel1(1));
        emitMatrix2(io, manifold.christoffel2(0));
        emitMatrix2(io, manifold.christoffel2(1));
        emitMatrix2(io, manifold.metricDerivative(0));
        emitMatrix2(io, manifold.metricDerivative(1));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeIntegrand', (io) => {
        const path = gpath(io, 2);
        const manifold = new PolySurfaceGeodesic(2);
        const der = GVector.fromArray([
            path[1].values[0] - path[0].values[0],
            path[1].values[1] - path[0].values[1]
        ]);
        io.outRealExact(manifold.probeIntegrand(path[0], der));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeSegmentLength', (io) => {
        const path = gpath(io, 2);
        const integralSamples = io.integer();
        const manifold = new PolySurfaceGeodesic(2);
        manifold.integralSamples = integralSamples;
        io.outRealExact(manifold.computeSegmentLength(path[0], path[1]));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeTotalLength', (io) => {
        const quantity = io.integer();
        const path = gpath(io, quantity);
        const manifold = new PolySurfaceGeodesic(2);
        io.outRealExact(manifold.computeTotalLength(quantity, path));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeSegmentCurvature', (io) => {
        const path = gpath(io, 2);
        const manifold = new PolySurfaceGeodesic(2);
        io.outRealExact(manifold.computeSegmentCurvature(path[0], path[1]));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeTotalCurvature', (io) => {
        const quantity = io.integer();
        const path = gpath(io, quantity);
        const manifold = new PolySurfaceGeodesic(2);
        io.outRealExact(manifold.computeTotalCurvature(quantity, path));
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.refine', (io) => {
        const path = gpath(io, 2);
        const mid = gvec2(io);
        const searchSamples = io.integer();
        const searchRadius = io.real();
        const derivativeStep = io.real();
        const manifold = new PolySurfaceGeodesic(2);
        manifold.searchSamples = searchSamples;
        manifold.searchRadius = searchRadius;
        manifold.derivativeStep = derivativeStep;
        const r = manifold.refine(path[0], mid, path[1]);
        io.outBool(r.changed);
        io.outVecExact(r.mid);
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.subdivide', (io) => {
        const path = gpath(io, 2);
        const searchSamples = io.integer();
        const searchRadius = io.real();
        const manifold = new PolySurfaceGeodesic(2);
        manifold.searchSamples = searchSamples;
        manifold.searchRadius = searchRadius;
        const r = manifold.subdivide(path[0], path[1]);
        io.outBool(r.changed);
        io.outVecExact(r.mid);
    }, { exact: true });

    family.case('RiemannianGeodesic.poly.computeGeodesic', (io) => {
        const path = gpath(io, 2);
        const subdivisions = io.integer();
        const refinements = io.integer();
        const searchSamples = io.integer();
        const searchRadius = io.real();
        const manifold = new PolySurfaceGeodesic(2);
        manifold.subdivisions = subdivisions;
        manifold.refinements = refinements;
        manifold.searchSamples = searchSamples;
        manifold.searchRadius = searchRadius;
        const r = manifold.computeGeodesic(path[0], path[1]);
        io.outInt(r.quantity);
        for (let i = 0; i < r.quantity; ++i) { io.outVecExact(r.path[i]); }
        io.outRealExact(manifold.computeTotalLength(r.quantity, r.path));
        io.outRealExact(manifold.computeSegmentLength(path[0], path[1]));
        io.outRealExact(manifold.computeTotalCurvature(r.quantity, r.path));
    }, { exact: true, timeout: 600000 });

    // refineCallback and the three progress accessors. subdivide installs a
    // no-op callback for the refine it performs itself, so the trace holds one
    // entry per refine of the refinement loop and none from the subdivision
    // pass.
    family.case('RiemannianGeodesic.poly.refineCallback', (io) => {
        const path = gpath(io, 2);
        const subdivisions = io.integer();
        const refinements = io.integer();
        const searchSamples = io.integer();
        const manifold = new PolySurfaceGeodesic(2);
        manifold.subdivisions = subdivisions;
        manifold.refinements = refinements;
        manifold.searchSamples = searchSamples;
        const trace: number[] = [];
        manifold.refineCallback = () => {
            trace.push(manifold.getSubdivisionStep());
            trace.push(manifold.getRefinementStep());
            trace.push(manifold.getCurrentQuantity());
        };
        manifold.computeGeodesic(path[0], path[1]);
        io.outInt(trace.length);
        for (const value of trace) { io.outInt(value); }
        io.outInt(manifold.getSubdivisionStep());
        io.outInt(manifold.getRefinementStep());
        io.outInt(manifold.getCurrentQuantity());
    }, { exact: true, timeout: 600000 });

    family.case('RiemannianGeodesic.construct.invalidDimension', (io) => {
        const dimension = io.integer();
        const manifold = new PolySurfaceGeodesic(dimension);
        io.outInt(manifold.getDimension());
    }, { exact: true });

    // ---- EllipsoidGeodesic ------------------------------------------------
    // sin and cos in the metric and the Christoffel symbols, so these cases
    // carry the default 1e-12 tolerance.

    family.case('EllipsoidGeodesic.computePosition', (io) => {
        const extents = io.reals(3);
        const point = gvec2(io);
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        io.outVec(eg.computePosition(point));
    });

    family.case('EllipsoidGeodesic.computeSegmentLength', (io) => {
        const extents = io.reals(3);
        const path = gpath(io, 2);
        const integralSamples = io.integer();
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        eg.integralSamples = integralSamples;
        io.outReal(eg.computeSegmentLength(path[0], path[1]));
    });

    family.case('EllipsoidGeodesic.computeTotalLength', (io) => {
        const extents = io.reals(3);
        const quantity = io.integer();
        const path = gpath(io, quantity);
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        io.outReal(eg.computeTotalLength(quantity, path));
    });

    family.case('EllipsoidGeodesic.computeSegmentCurvature', (io) => {
        const extents = io.reals(3);
        const path = gpath(io, 2);
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        io.outReal(eg.computeSegmentCurvature(path[0], path[1]));
    });

    family.case('EllipsoidGeodesic.computeTotalCurvature', (io) => {
        const extents = io.reals(3);
        const quantity = io.integer();
        const path = gpath(io, quantity);
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        io.outReal(eg.computeTotalCurvature(quantity, path));
    });

    // searchSamples = 0 makes the steepest-descent search a single sample at
    // tRay = 0 whose candidate is bit-identical to the midpoint, so no
    // libm-derived value reaches the path and the case is exact.
    family.case('EllipsoidGeodesic.computeGeodesic.noSearch', (io) => {
        const extents = io.reals(3);
        const path = gpath(io, 2);
        const subdivisions = io.integer();
        const refinements = io.integer();
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        eg.subdivisions = subdivisions;
        eg.refinements = refinements;
        eg.searchSamples = 0;
        const r = eg.computeGeodesic(path[0], path[1]);
        io.outInt(r.quantity);
        for (let i = 0; i < r.quantity; ++i) { io.outVecExact(r.path[i]); }
    }, { exact: true, timeout: 600000 });

    // The C++ generator rejected the configurations whose winning sampled
    // length is within a relative 1e-9 of any other sample or of the length at
    // the current midpoint, so both builds select the same sample and the same
    // boolean; the surviving difference is the rounding of the gradient, well
    // inside the default 1e-12 tolerance.
    family.case('EllipsoidGeodesic.refine.separated', (io) => {
        const extents = io.reals(3);
        const end0 = gvec2(io);
        const mid = gvec2(io);
        const end1 = gvec2(io);
        const searchSamples = io.integer();
        const searchRadius = io.real();
        const eg = new EllipsoidGeodesic(extents[0], extents[1], extents[2]);
        eg.searchSamples = searchSamples;
        eg.searchRadius = searchRadius;
        const r = eg.refine(end0, mid, end1);
        io.outBool(r.changed);
        io.outVec(r.mid);
        // The gradient is a centered difference of two nearly equal segment
        // lengths multiplied by mDerivativeFactor = 0.5/1e-4 = 5000, so a
        // 1-ulp libm difference in a length becomes about 5e-12 relative in
        // the gradient before the search scales it by tRay. Measured maximum
        // over a 2000-record deep run: 3.57e-13. The tolerance is still three
        // orders of magnitude below the distance between two search samples,
        // so a flipped argument of the minimum would still fail the case.
    }, { tol: 1e-10, timeout: 600000 });

    family.finish();
});
