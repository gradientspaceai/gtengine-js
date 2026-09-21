// Replays oracle/cpp/cases/v17-curves.cpp. Keep the two files in the same
// order.
//
// Every case is exact unless the comment says otherwise: the only C math
// library calls on these paths are std::pow(x, 1.5) in FrenetFrame2/3 and
// ImplicitCurve2 GetCurvature, and std::acos in SampleCircularArc (which
// only chooses the sample count, pinned by the generator).
import { describe } from 'vitest';
import { Arc2 } from '../../src/Arc2.js';
import { BasisFunction, BasisFunctionInput, UniqueKnot }
    from '../../src/BasisFunction.js';
import { BezierCurve } from '../../src/BezierCurve.js';
import { CLODPolyline } from '../../src/CLODPolyline.js';
import { DarbouxFrame3 } from '../../src/DarbouxFrame.js';
import { FrenetFrame2, FrenetFrame3 } from '../../src/FrenetFrame.js';
import { ImplicitCurve2 } from '../../src/ImplicitCurve2.js';
import { MassSpringCurve } from '../../src/MassSpringCurve.js';
import {
    NURBSCircularArcDegree2, NURBSFullCircleDegree3, NURBSHalfCircleDegree3,
    NURBSQuarterCircleDegree2, NURBSQuarterCircleDegree4
} from '../../src/NURBSCircle.js';
import { NURBSCurve } from '../../src/NURBSCurve.js';
import {
    NURBSEighthSphereDegree4, NURBSFullSphereDegree3, NURBSHalfSphereDegree3,
    createEighthSphereValues
} from '../../src/NURBSSphere.js';
import { NaturalCubicSpline } from '../../src/NaturalCubicSpline.js';
import { NaturalQuinticSpline } from '../../src/NaturalQuinticSpline.js';
import { ParametricCurve } from '../../src/ParametricCurve.js';
import type { ParametricSurface } from '../../src/ParametricSurface.js';
import { PolylineOffset } from '../../src/PolylineOffset.js';
import { ReparameterizeByArclength } from '../../src/ReparameterizeByArclength.js';
import { SampleCircularArc } from '../../src/SampleCircularArc.js';
import { TCBSplineCurve } from '../../src/TCBSplineCurve.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The replay of RecordBasis.
function basisInput(io: OracleIO): BasisFunctionInput {
    const input = new BasisFunctionInput();
    input.numControls = io.integer();
    input.degree = io.integer();
    input.uniform = io.boolean();
    input.periodic = io.boolean();
    input.numUniqueKnots = io.integer();
    input.uniqueKnots = [];
    for (let i = 0; i < input.numUniqueKnots; ++i) {
        const t = io.real();
        input.uniqueKnots.push(new UniqueKnot(t, io.integer()));
    }
    return input;
}

// The replay of RecordControls3 / MakeControls2 / MakeSplineDerivatives.
function points(io: OracleIO, count: number, dimension: number): Vector[] {
    const P = new Array<Vector>(count);
    for (let i = 0; i < count; ++i) { P[i] = io.vec(dimension); }
    return P;
}

// The replay of MakeWeights and MakeTimes (a plain run of recorded reals).
function reals(io: OracleIO, count: number): number[] {
    return io.reals(count);
}

// Emit the leading 'order + 1' entries of a jet.
function outJet(io: OracleIO, jet: Vector[], order: number): void {
    for (let i = 0; i <= order; ++i) { io.outVec(jet[i]); }
}

// The conic of the C++ case file's ConicCurve2, with the same expression
// order.
class ConicCurve2 extends ImplicitCurve2 {
    constructor(private readonly c: number[]) { super(); }

    override f(p: Vector): number {
        const v = p.values;
        return this.c[0] + this.c[1] * v[0] + this.c[2] * v[1]
            + this.c[3] * v[0] * v[0] + this.c[4] * v[0] * v[1]
            + this.c[5] * v[1] * v[1];
    }

    override fx(p: Vector): number {
        const v = p.values;
        return this.c[1] + 2 * this.c[3] * v[0] + this.c[4] * v[1];
    }

    override fy(p: Vector): number {
        const v = p.values;
        return this.c[2] + this.c[4] * v[0] + 2 * this.c[5] * v[1];
    }

    override fxx(_p: Vector): number { return 2 * this.c[3]; }
    override fxy(_p: Vector): number { return this.c[4]; }
    override fyy(_p: Vector): number { return 2 * this.c[5]; }
}

describe('oracle: v17-curves', () => {
    const family = new OracleFamily('v17-curves');

    // ------------------------------------------------------- BasisFunction

    family.case('BasisFunction.create', (io) => {
        const bf = new BasisFunction(basisInput(io));
        io.outInt(bf.getNumControls());
        io.outInt(bf.getDegree());
        io.outInt(bf.getNumUniqueKnots());
        io.outInt(bf.getNumKnots());
        io.outReal(bf.getMinDomain());
        io.outReal(bf.getMaxDomain());
        io.outBool(bf.isOpen());
        io.outBool(bf.isUniform());
        io.outBool(bf.isPeriodic());
        const unique = bf.getUniqueKnots();
        for (let i = 0; i < bf.getNumUniqueKnots(); ++i) {
            io.outReal(unique[i].t);
            io.outInt(unique[i].multiplicity);
        }
        const knots = bf.getKnots();
        for (let i = 0; i < bf.getNumKnots(); ++i) { io.outReal(knots[i]); }
    }, { exact: true });

    family.case('BasisFunction.evaluate', (io) => {
        const bf = new BasisFunction(basisInput(io));
        const order = io.integer();
        const t = io.real();
        const r = bf.evaluate(t, order);
        io.outInt(r.minIndex);
        io.outInt(r.maxIndex);
        const numValues = bf.getNumControls() + bf.getDegree();
        for (let o = 0; o <= order; ++o) {
            for (let i = 0; i < numValues; ++i) { io.outReal(bf.getValue(o, i)); }
        }
    }, { exact: true });

    family.case('BasisFunction.evaluate.repeated', (io) => {
        const bf = new BasisFunction(basisInput(io));
        const t0 = io.real();
        const t1 = io.real();
        const tmin = bf.getMinDomain();
        const tmax = bf.getMaxDomain();
        const u0 = tmin + (tmax - tmin) * t0;
        const u1 = tmin + (tmax - tmin) * t1;
        const r0 = bf.evaluate(u0, 3);
        io.outInt(r0.minIndex);
        io.outInt(r0.maxIndex);
        const r1 = bf.evaluate(u1, 3);
        io.outInt(r1.minIndex);
        io.outInt(r1.maxIndex);
        const numValues = bf.getNumControls() + bf.getDegree();
        for (let o = 0; o <= 3; ++o) {
            for (let i = 0; i < numValues; ++i) { io.outReal(bf.getValue(o, i)); }
        }
    }, { exact: true });

    family.case('BasisFunctionInput.openUniform', (io) => {
        const degree = io.integer();
        // The C++ case records the increment, not numControls itself.
        const numControls = degree + io.integer();
        const input = new BasisFunctionInput(numControls, degree);
        io.outInt(input.numControls);
        io.outInt(input.degree);
        io.outBool(input.uniform);
        io.outBool(input.periodic);
        io.outInt(input.numUniqueKnots);
        for (let i = 0; i < input.numUniqueKnots; ++i) {
            io.outReal(input.uniqueKnots[i].t);
            io.outInt(input.uniqueKnots[i].multiplicity);
        }
    }, { exact: true });

    family.case('BasisFunction.getValue.invalid', (io) => {
        const bf = new BasisFunction(basisInput(io));
        const which = io.integer();
        bf.evaluate(bf.getMinDomain(), 1);
        if (which === 0) {
            io.outReal(bf.getValue(0, bf.getNumControls() + bf.getDegree()));
        }
        else {
            io.outReal(bf.getValue(4, 0));
        }
    }, { exact: true });

    // --------------------------------------------------------- BezierCurve

    family.case('BezierCurve.evaluate', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const order = io.integer();
        const t = io.real();
        const curve = new BezierCurve(3, degree, controls);
        const jet = curve.createJet();
        curve.evaluate(t, order, jet);
        io.outInt(curve.getDegree());
        io.outInt(curve.getNumControls());
        io.outBool(curve.isConstructed());
        outJet(io, jet, order);
        const stored = curve.getControls();
        for (let i = 0; i < curve.getNumControls(); ++i) { io.outVec(stored[i]); }
    }, { exact: true });

    family.case('BezierCurve.evaluate.invalidOrder', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const t = io.real();
        const curve = new BezierCurve(3, degree, controls);
        const jet = curve.createJet();
        curve.evaluate(t, 4, jet);
        for (let i = 0; i < 4; ++i) { io.outVec(jet[i]); }
    }, { exact: true });

    // ----------------------------------------------------- ParametricCurve

    family.case('ParametricCurve.differentialGeometry', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const t = io.real();
        const curve = new BezierCurve(3, degree, controls);
        io.outReal(curve.getTMin());
        io.outReal(curve.getTMax());
        io.outInt(curve.getNumSegments());
        io.outVec(curve.getPosition(t));
        io.outVec(curve.getTangent(t));
        io.outReal(curve.getSpeed(t));
    }, { exact: true });

    family.case('ParametricCurve.setTimeInterval', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const tmin = io.real();
        const tmax = io.real();
        const t = io.real();
        const curve = new BezierCurve(3, degree, controls);
        curve.setTimeInterval(tmin, tmax);
        io.outReal(curve.getTMin());
        io.outReal(curve.getTMax());
        const times = curve.getTimes();
        io.outReal(times[0]);
        io.outReal(times[1]);
        io.outVec(curve.getPosition(t));
    }, { exact: true });

    family.case('ParametricCurve.getLength', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const rombergOrder = io.integer();
        const t0 = io.real();
        const t1 = io.real();
        const curve = new BezierCurve(3, degree, controls);
        curve.setRombergOrder(rombergOrder);
        io.outReal(curve.getLength(t0, t1));
        io.outReal(curve.getTotalLength());
        io.outReal(curve.getLength(t1, t0));
        io.outReal(curve.getTotalLength());
    }, { exact: true });

    family.case('ParametricCurve.getTime', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const rombergOrder = io.integer();
        const maxBisections = io.integer();
        const fraction = io.real();
        const curve = new BezierCurve(3, degree, controls);
        curve.setRombergOrder(rombergOrder);
        curve.setMaxBisections(maxBisections);
        const total = curve.getTotalLength();
        io.outReal(total);
        io.outReal(curve.getTime(fraction * total));
    }, { exact: true });

    family.case('ParametricCurve.subdivideByTime', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const numPoints = io.integer();
        const curve = new BezierCurve(3, degree, controls);
        for (const p of curve.subdivideByTime(numPoints)) { io.outVec(p); }
    }, { exact: true });

    family.case('ParametricCurve.subdivideByLength', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const rombergOrder = io.integer();
        const maxBisections = io.integer();
        const numPoints = io.integer();
        const curve = new BezierCurve(3, degree, controls);
        curve.setRombergOrder(rombergOrder);
        curve.setMaxBisections(maxBisections);
        for (const p of curve.subdivideByLength(numPoints)) { io.outVec(p); }
    }, { exact: true });

    // ---------------------------------------------------------- NURBSCurve

    family.case('NURBSCurve.evaluate', (io) => {
        const input = basisInput(io);
        const controls = points(io, input.numControls, 3);
        const weights = reals(io, input.numControls);
        const order = io.integer();
        const curve = new NURBSCurve(3, input, controls, weights);
        const t = io.real();
        const jet = curve.createJet();
        curve.evaluate(t, order, jet);
        io.outReal(curve.getTMin());
        io.outReal(curve.getTMax());
        io.outInt(curve.getNumControls());
        outJet(io, jet, order);
    }, { exact: true });

    family.case('NURBSCurve.evaluate.invalidOrder', (io) => {
        const input = basisInput(io);
        const controls = points(io, input.numControls, 3);
        const weights = reals(io, input.numControls);
        const curve = new NURBSCurve(3, input, controls, weights);
        const t = io.real();
        const jet = curve.createJet();
        curve.evaluate(t, 4, jet);
        for (let i = 0; i < 4; ++i) { io.outVec(jet[i]); }
    }, { exact: true });

    family.case('NURBSCurve.accessors', (io) => {
        const input = basisInput(io);
        const controls = points(io, input.numControls, 3);
        const weights = reals(io, input.numControls);
        const setIndex = io.integer();
        const newControl = io.vec(3);
        const newWeight = io.real();
        const getIndex = io.integer();
        const curve = new NURBSCurve(3, input, controls, weights);
        curve.setControl(setIndex, newControl);
        curve.setWeight(setIndex, newWeight);
        io.outVec(curve.getControl(getIndex));
        io.outReal(curve.getWeight(getIndex));
        const stored = curve.getControls();
        const storedW = curve.getWeights();
        for (let i = 0; i < curve.getNumControls(); ++i) {
            io.outVec(stored[i]);
            io.outReal(storedW[i]);
        }
    }, { exact: true });

    family.case('NURBSCurve.deferred', (io) => {
        const curve = new NURBSCurve(3, basisInput(io));
        io.outInt(curve.getNumControls());
        const stored = curve.getControls();
        const storedW = curve.getWeights();
        for (let i = 0; i < curve.getNumControls(); ++i) {
            io.outVec(stored[i]);
            io.outReal(storedW[i]);
        }
    }, { exact: true });

    // --------------------------------------------------------- NURBSCircle

    function arcOf(io: OracleIO): Arc2 {
        const center = io.vec(2);
        const radius = io.real();
        const e0 = io.vec(2);
        const e1 = io.vec(2);
        return Arc2.fromCenterRadiusEnds(center, radius, e0, e1);
    }

    function curveJet(io: OracleIO, curve: NURBSCurve, t: number,
        order: number): void {
        const jet = curve.createJet();
        curve.evaluate(t, order, jet);
        io.outReal(curve.getTMin());
        io.outReal(curve.getTMax());
        outJet(io, jet, order);
    }

    family.case('NURBSCircle.quarterDegree2', (io) => {
        const order = io.integer();
        const t = io.real();
        curveJet(io, new NURBSQuarterCircleDegree2(), t, order);
    }, { exact: true });

    family.case('NURBSCircle.quarterDegree4', (io) => {
        const order = io.integer();
        const t = io.real();
        curveJet(io, new NURBSQuarterCircleDegree4(), t, order);
    }, { exact: true });

    family.case('NURBSCircle.halfDegree3', (io) => {
        const order = io.integer();
        const t = io.real();
        curveJet(io, new NURBSHalfCircleDegree3(), t, order);
    }, { exact: true });

    family.case('NURBSCircle.fullDegree3', (io) => {
        const order = io.integer();
        const t = io.real();
        curveJet(io, new NURBSFullCircleDegree3(), t, order);
    }, { exact: true });

    family.case('NURBSCircle.arcDegree2', (io) => {
        const a = arcOf(io);
        const order = io.integer();
        const t = io.real();
        const curve = new NURBSCircularArcDegree2(a);
        io.outReal(curve.getWeight(0));
        io.outReal(curve.getWeight(1));
        io.outReal(curve.getWeight(2));
        for (let i = 0; i < 3; ++i) { io.outVec(curve.getControl(i)); }
        curveJet(io, curve, t, order);
    }, { exact: true });

    // --------------------------------------------------------- NURBSSphere

    family.case('NURBSSphere.eighthDegree4', (io) => {
        const maxOrder = io.integer();
        const u = io.real();
        const v = io.real();
        const patch = new NURBSEighthSphereDegree4();
        const values = createEighthSphereValues();
        patch.evaluate(u, v, maxOrder, values);
        const numValues = (maxOrder === 0 ? 1 : (maxOrder === 1 ? 3 : 6));
        for (let i = 0; i < numValues; ++i) { io.outVec(values[i]); }
    }, { exact: true });

    family.case('NURBSSphere.halfDegree3', (io) => {
        const order = io.integer();
        const u = io.real();
        const v = io.real();
        const surface = new NURBSHalfSphereDegree3();
        const jet = surface.createJet();
        surface.evaluate(u, v, order, jet);
        const numValues = (order === 0 ? 1 : (order === 1 ? 3 : 6));
        for (let i = 0; i < numValues; ++i) { io.outVec(jet[i]); }
    }, { exact: true });

    family.case('NURBSSphere.fullDegree3', (io) => {
        const order = io.integer();
        const u = io.real();
        const v = io.real();
        const surface = new NURBSFullSphereDegree3();
        const jet = surface.createJet();
        surface.evaluate(u, v, order, jet);
        const numValues = (order === 0 ? 1 : (order === 1 ? 3 : 6));
        for (let i = 0; i < numValues; ++i) { io.outVec(jet[i]); }
    }, { exact: true });

    // -------------------------------------------------------- DarbouxFrame

    function sphereSurface(which: number): ParametricSurface {
        return which === 0 ? new NURBSHalfSphereDegree3()
            : new NURBSFullSphereDegree3();
    }

    family.case('DarbouxFrame.compute', (io) => {
        const which = io.integer();
        const u = io.real();
        const v = io.real();
        const frame = new DarbouxFrame3(sphereSurface(which));
        const r = frame.compute(u, v);
        io.outVec(r.position);
        io.outVec(r.tangent0);
        io.outVec(r.tangent1);
        io.outVec(r.normal);
    }, { exact: true });

    family.case('DarbouxFrame.getPrincipalInformation', (io) => {
        const which = io.integer();
        const u = io.real();
        const v = io.real();
        const frame = new DarbouxFrame3(sphereSurface(which));
        const r = frame.getPrincipalInformation(u, v);
        io.outReal(r.curvature0);
        io.outReal(r.curvature1);
        io.outVec(r.direction0);
        io.outVec(r.direction1);
    }, { exact: true });

    // ------------------------------------------------------ TCBSplineCurve

    interface TCBReplay {
        point: Vector[];
        time: number[];
        tension: number[];
        continuity: number[];
        bias: number[];
        lambda: number[];
        firstOut?: Vector;
        lastIn?: Vector;
    }

    function tcbInput(io: OracleIO): TCBReplay {
        const n = io.integer();
        const point = points(io, n, 3);
        const time = reals(io, n);
        const tension = reals(io, n);
        const continuity = reals(io, n);
        const bias = reals(io, n);
        const useLambda = io.integer();
        const lambda = useLambda === 1 ? reals(io, n) : [];
        const useFirst = io.integer();
        const firstOut = useFirst === 1 ? io.vec(3) : undefined;
        const useLast = io.integer();
        const lastIn = useLast === 1 ? io.vec(3) : undefined;
        return {
            point, time, tension, continuity, bias, lambda, firstOut, lastIn
        };
    }

    function tcbCurve(r: TCBReplay): TCBSplineCurve {
        return new TCBSplineCurve(3, r.point, r.time, r.tension, r.continuity,
            r.bias, r.lambda, r.firstOut, r.lastIn);
    }

    family.case('TCBSplineCurve.evaluate', (io) => {
        const r = tcbInput(io);
        const order = io.integer();
        const t = io.real();
        const curve = tcbCurve(r);
        io.outInt(curve.getNumKeyFrames());
        for (const v of curve.getInTangents()) { io.outVec(v); }
        for (const v of curve.getOutTangents()) { io.outVec(v); }
        const jet = curve.createJet();
        curve.evaluate(t, order, jet);
        outJet(io, jet, order);
    }, { exact: true });

    family.case('TCBSplineCurve.degenerateLambda', (io) => {
        const n = io.integer();
        io.integer();
        const point = points(io, n, 3);
        const time = new Array<number>(n);
        for (let i = 0; i < n; ++i) { time[i] = i; }
        const zeros = new Array<number>(n).fill(0);
        const lambda = new Array<number>(n).fill(1);
        const t = io.real();
        const curve = new TCBSplineCurve(3, point, time, zeros, zeros, zeros,
            lambda);
        for (const v of curve.getInTangents()) {
            for (let j = 0; j < 3; ++j) {
                io.outBool(Number.isNaN(v.values[j]));
            }
        }
        for (const v of curve.getOutTangents()) {
            for (let j = 0; j < 3; ++j) {
                io.outBool(Number.isNaN(v.values[j]));
            }
        }
        const jet = curve.createJet();
        curve.evaluate(t, 3, jet);
        for (let i = 0; i < 4; ++i) {
            for (let j = 0; j < 3; ++j) {
                io.outBool(Number.isNaN(jet[i].values[j]));
            }
        }
    }, { exact: true });

    // The port sets mConstructed, upstream never does, so every record
    // deviates (docs/UPSTREAM-FINDINGS.md, issue #182).
    family.case('TCBSplineCurve.isConstructed', (io) => {
        io.outBool(tcbCurve(tcbInput(io)).isConstructed());
    }, { deviation: 'issue #182: TCBSplineCurve never sets mConstructed' });

    // ------------------------------ NaturalCubicSpline / NaturalQuinticSpline

    function emitSpline(io: OracleIO, spline: ParametricCurve,
        polynomials: readonly Vector[][], t: number, order: number): void {
        io.outInt(spline.getNumSegments());
        for (const poly of polynomials) {
            for (const c of poly) { io.outVec(c); }
        }
        const jet = spline.createJet();
        spline.evaluate(t, order, jet);
        outJet(io, jet, order);
    }

    family.case('NaturalCubicSpline.free', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const times = reals(io, n);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalCubicSpline.createFree(f0, times);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    family.case('NaturalCubicSpline.closed', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const times = reals(io, n);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalCubicSpline.createClosed(f0, times);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    family.case('NaturalCubicSpline.clamped', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const times = reals(io, n);
        const d0 = io.vec(3);
        const d1 = io.vec(3);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalCubicSpline.createClamped(f0, times, d0, d1);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    family.case('ParametricCurve.getLength.multiSegment', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const times = reals(io, n);
        const rombergOrder = io.integer();
        const maxBisections = io.integer();
        const a = io.real();
        const b = io.real();
        const spline = NaturalCubicSpline.createFree(f0, times);
        spline.setRombergOrder(rombergOrder);
        spline.setMaxBisections(maxBisections);
        const tmin = spline.getTMin();
        const tmax = spline.getTMax();
        const t0 = tmin + (tmax - tmin) * a;
        const t1 = tmin + (tmax - tmin) * b;
        io.outReal(spline.getLength(t0, t1));
        io.outReal(spline.getTotalLength());
        io.outReal(spline.getTime(0.5 * spline.getTotalLength()));
        io.outReal(spline.getLength(times[0], times[1]));
        io.outReal(spline.getLength(times[1], times[n - 1]));
    }, { exact: true });

    family.case('NaturalQuinticSpline.free', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const f1 = points(io, n, 3);
        const times = reals(io, n);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalQuinticSpline.createFree(f0, f1, times);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    family.case('NaturalQuinticSpline.closed', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const f1 = points(io, n, 3);
        const times = reals(io, n);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalQuinticSpline.createClosed(f0, f1, times);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    family.case('NaturalQuinticSpline.clamped', (io) => {
        const n = io.integer();
        const f0 = points(io, n, 3);
        const f1 = points(io, n, 3);
        const times = reals(io, n);
        const d0 = io.vec(3);
        const d1 = io.vec(3);
        const order = io.integer();
        const t = io.real();
        const spline = NaturalQuinticSpline.createClamped(f0, f1, times, d0, d1);
        emitSpline(io, spline, spline.getPolynomials(), t, order);
    }, { exact: true });

    // --------------------------------------------------------- FrenetFrame

    family.case('FrenetFrame2.compute', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 2);
        const t = io.real();
        const frame = new FrenetFrame2(new BezierCurve(2, degree, controls));
        const r = frame.compute(t);
        io.outVec(r.position);
        io.outVec(r.tangent);
        io.outVec(r.normal);
    }, { exact: true });

    // GetCurvature calls std::pow(speedSqr, 1.5) in C++ and Math.pow here.
    family.case('FrenetFrame2.getCurvature', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 2);
        const t = io.real();
        const frame = new FrenetFrame2(new BezierCurve(2, degree, controls));
        io.outReal(frame.getCurvature(t));
    });

    family.case('FrenetFrame3.compute', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const t = io.real();
        const frame = new FrenetFrame3(new BezierCurve(3, degree, controls));
        const r = frame.compute(t);
        io.outVec(r.position);
        io.outVec(r.tangent);
        io.outVec(r.normal);
        io.outVec(r.binormal);
    }, { exact: true });

    // getCurvature uses Math.pow; getTorsion is arithmetic only and is held
    // to bit-identity.
    family.case('FrenetFrame3.getCurvatureAndTorsion', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const t = io.real();
        const frame = new FrenetFrame3(new BezierCurve(3, degree, controls));
        io.outReal(frame.getCurvature(t));
        io.outRealExact(frame.getTorsion(t));
    });

    // ------------------------------------------------------ ImplicitCurve2

    family.case('ImplicitCurve2.evaluate', (io) => {
        const curve = new ConicCurve2(io.reals(6));
        const p = io.vec(2);
        const epsilon = io.real();
        io.outReal(curve.f(p));
        io.outReal(curve.fx(p));
        io.outReal(curve.fy(p));
        io.outReal(curve.fxx(p));
        io.outReal(curve.fxy(p));
        io.outReal(curve.fyy(p));
        io.outBool(curve.isOnCurve(p, epsilon));
        io.outVec(curve.getGradient(p));
        io.outMat(curve.getHessian(p));
        const frame = curve.getFrame(p);
        io.outVec(frame.tangent);
        io.outVec(frame.normal);
    }, { exact: true });

    // getCurvature uses Math.pow(fx^2 + fy^2, 1.5).
    family.case('ImplicitCurve2.getCurvature', (io) => {
        const curve = new ConicCurve2(io.reals(6));
        const p = io.vec(2);
        const r = curve.getCurvature(p);
        io.outBool(r.valid);
        io.outReal(r.curvature);
    });

    // ------------------------------------------ ReparameterizeByArclength

    family.case('ReparameterizeByArclength.getT', (io) => {
        const degree = io.integer();
        const controls = points(io, degree + 1, 3);
        const rombergOrder = io.integer();
        const useBisection = io.integer();
        const fraction = io.real();
        const curve = new BezierCurve(3, degree, controls);
        curve.setRombergOrder(rombergOrder);
        const reparam = new ReparameterizeByArclength(curve);
        io.outReal(reparam.getTMin());
        io.outReal(reparam.getTMax());
        io.outReal(reparam.getTotalArclength());
        const s = fraction * reparam.getTotalArclength();
        const output = reparam.getT(s, useBisection === 1);
        io.outReal(output.t);
        io.outReal(output.f);
        io.outInt(output.numIterations);
    }, { exact: true });

    // ---------------------------------------------------- SampleCircularArc

    family.case('SampleCircularArc.compute', (io) => {
        const points2 = new SampleCircularArc().compute(arcOf(io));
        io.outInt(points2.length);
        for (const p of points2) { io.outVec(p); }
    }, { exact: true });

    family.case('SampleCircularArc.compute.largeArc', (io) => {
        const points2 = new SampleCircularArc().compute(arcOf(io));
        io.outInt(points2.length);
        for (const p of points2) { io.outVec(p); }
    }, { exact: true });

    // ------------------------------------------------------ PolylineOffset

    family.case('PolylineOffset.execute', (io) => {
        const isOpen = io.integer();
        const n = io.integer();
        const vertices = points(io, n, 2);
        const distance = io.real();
        const which = io.integer();
        const offset = new PolylineOffset(vertices, isOpen === 1);
        const r = offset.execute(distance, which !== 1, which !== 0);
        io.outInt(r.rightPolyline.length);
        for (const p of r.rightPolyline) { io.outVec(p); }
        io.outInt(r.leftPolyline.length);
        for (const p of r.leftPolyline) { io.outVec(p); }
    }, { exact: true });

    family.case('PolylineOffset.throwParity', (io) => {
        const which = io.integer();
        const isOpen = io.integer();
        const n = io.integer();
        const vertices = points(io, n, 2);
        const distance = io.real();
        const offset = new PolylineOffset(vertices, isOpen === 1);
        const r = offset.execute(distance, which !== 1, false);
        io.outInt(r.rightPolyline.length);
        for (const p of r.rightPolyline) { io.outVec(p); }
    }, { exact: true });

    // -------------------------------------------------------- CLODPolyline

    family.case('CLODPolyline.construct', (io) => {
        const closed = io.integer();
        // The C++ case draws the vertex count only for a closed polyline.
        const n = closed === 1 ? io.integer() : 2;
        const vertices = points(io, n, 2);
        const polyline = new CLODPolyline(vertices, closed === 1);
        io.outInt(polyline.getNumVertices());
        io.outBool(polyline.getClosed());
        io.outInt(polyline.getMinLevelOfDetail());
        io.outInt(polyline.getMaxLevelOfDetail());
        io.outInt(polyline.getLevelOfDetail());
        io.outInt(polyline.getNumEdges());
        for (const v of polyline.getVertices()) { io.outVec(v); }
        const edges = polyline.getEdges();
        for (let i = 0; i < 2 * polyline.getNumEdges(); ++i) {
            io.outInt(edges[i]);
        }
    }, { exact: true });

    family.case('CLODPolyline.setLevelOfDetail', (io) => {
        const n = io.integer();
        const vertices = points(io, n, 2);
        const lod = io.integer();
        const lod2 = io.integer();
        const polyline = new CLODPolyline(vertices, true);
        polyline.setLevelOfDetail(lod);
        io.outInt(polyline.getLevelOfDetail());
        io.outInt(polyline.getNumEdges());
        const edges = polyline.getEdges();
        for (let i = 0; i < 2 * polyline.getNumEdges(); ++i) {
            io.outInt(edges[i]);
        }
        polyline.setLevelOfDetail(lod2);
        io.outInt(polyline.getLevelOfDetail());
        io.outInt(polyline.getNumEdges());
        for (let i = 0; i < 2 * polyline.getNumEdges(); ++i) {
            io.outInt(edges[i]);
        }
    }, { exact: true });

    // ----------------------------------------------------- MassSpringCurve

    family.case('MassSpringCurve.update', (io) => {
        const numParticles = io.integer();
        const step = io.real();
        const numSteps = io.integer();
        const immovable = io.integer();
        const mass = io.reals(numParticles);
        const position = new Array<Vector>(numParticles);
        for (let i = 0; i < numParticles; ++i) {
            const dx = io.real();
            const py = io.real();
            const pz = io.real();
            position[i] = Vector.fromArray([i + dx, py, pz]);
        }
        const velocity = points(io, numParticles, 3);
        const numSprings = numParticles - 1;
        const constant = io.reals(numSprings);
        const restLength = io.reals(numSprings);
        const time = io.real();

        const system = new MassSpringCurve(3, numParticles, step);
        for (let i = 0; i < numParticles; ++i) {
            system.setMass(i, i === immovable ? Number.MAX_VALUE : mass[i]);
            system.setPosition(i, position[i]);
            system.setVelocity(i, velocity[i]);
        }
        for (let i = 0; i < numSprings; ++i) {
            system.setConstant(i, constant[i]);
            system.setLength(i, restLength[i]);
        }

        io.outInt(system.getNumParticles());
        io.outInt(system.getNumSprings());
        io.outReal(system.getStep());
        let t = time;
        for (let s = 0; s < numSteps; ++s) {
            system.update(t);
            t += step;
        }
        for (let i = 0; i < numParticles; ++i) {
            io.outReal(system.getMass(i));
            io.outVec(system.getPosition(i));
            io.outVec(system.getVelocity(i));
        }
        for (let i = 0; i < numSprings; ++i) {
            io.outReal(system.getConstant(i));
            io.outReal(system.getLength(i));
        }
    }, { exact: true });

    family.finish();
});
