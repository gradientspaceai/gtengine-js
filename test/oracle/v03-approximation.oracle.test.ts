// Replays oracle/cpp/cases/v03-approximation.cpp. Keep the two files in the
// same order.
//
// Every case is declared exact: no fitting path in this group reaches the C
// math library except through Math.sqrt (see the comment at the top of the
// C++ case file).
import { describe } from 'vitest';
import { ApprCircle2 } from '../../src/ApprCircle2.js';
import { ApprGaussian2 } from '../../src/ApprGaussian2.js';
import { ApprGaussian3 } from '../../src/ApprGaussian3.js';
import { ApprHeightLine2 } from '../../src/ApprHeightLine2.js';
import { ApprHeightPlane3 } from '../../src/ApprHeightPlane3.js';
import { ApprOrthogonalLine2 } from '../../src/ApprOrthogonalLine2.js';
import { ApprOrthogonalLine3 } from '../../src/ApprOrthogonalLine3.js';
import { ApprOrthogonalPlane3 } from '../../src/ApprOrthogonalPlane3.js';
import { ApprPolynomial2 } from '../../src/ApprPolynomial2.js';
import { ApprPolynomial3 } from '../../src/ApprPolynomial3.js';
import { ApprPolynomial4 } from '../../src/ApprPolynomial4.js';
import { ApprPolynomialSpecial2 } from '../../src/ApprPolynomialSpecial2.js';
import { ApprPolynomialSpecial3 } from '../../src/ApprPolynomialSpecial3.js';
import { ApprPolynomialSpecial4 } from '../../src/ApprPolynomialSpecial4.js';
import { ApprQuery } from '../../src/ApprQuery.js';
import { ApprSphere3 } from '../../src/ApprSphere3.js';
import { Hypersphere } from '../../src/Hypersphere.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// The replay of MakePoints2/MakePoints3: the count, the (ignored) generator
// mode, then the coordinates.
function points(io: OracleIO, dimension: number): Vector[] {
    const n = io.integer();
    io.integer();
    const P = new Array<Vector>(n);
    for (let i = 0; i < n; ++i) { P[i] = io.vec(dimension); }
    return P;
}

// The replay of MakeObs2/MakeObs3/MakeObs4: the count, the (ignored)
// generator mode, then 'dimension' doubles per observation.
function observations(io: OracleIO, dimension: number): number[][] {
    const n = io.integer();
    io.integer();
    const obs = new Array<number[]>(n);
    for (let i = 0; i < n; ++i) { obs[i] = io.reals(dimension); }
    return obs;
}

// The replay of MakeDegrees.
function degrees(io: OracleIO): number[] {
    const count = io.integer();
    const values = new Array<number>(count);
    for (let i = 0; i < count; ++i) { values[i] = io.integer(); }
    return values;
}

// The identity degree list that the special-polynomial cases pair with the
// drawn x list.
function identityDegrees(count: number, stride: number): number[] {
    const values = new Array<number>(count);
    for (let i = 0; i < count; ++i) { values[i] = stride * i; }
    return values;
}

// The replay of MakeIndices.
function indices(io: OracleIO): number[] {
    const count = io.integer();
    io.integer();
    const values = new Array<number>(count);
    for (let i = 0; i < count; ++i) { values[i] = io.integer(); }
    return values;
}

describe('oracle: v03-approximation', () => {
    const family = new OracleFamily('v03-approximation');

    // ------------------------------------------------------------ ApprQuery

    function heightLine2(io: OracleIO, fitter: ApprHeightLine2,
        success: boolean): void {
        io.outBool(success);
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
    }

    family.case('ApprQuery.fit.all', (io) => {
        const P = points(io, 2);
        const fitter = new ApprHeightLine2();
        heightLine2(io, fitter, fitter.fit(P));
    }, { exact: true });

    family.case('ApprQuery.fit.range', (io) => {
        const P = points(io, 2);
        const imin = io.integer();
        const imax = io.integer();
        const fitter = new ApprHeightLine2();
        heightLine2(io, fitter, fitter.fit(P, imin, imax));
    }, { exact: true });

    family.case('ApprQuery.fit.indexed', (io) => {
        const P = points(io, 2);
        const subset = indices(io);
        const fitter = new ApprHeightLine2();
        heightLine2(io, fitter, fitter.fit(P, subset));
    }, { exact: true });

    family.case('ApprQuery.fit.numIndices', (io) => {
        const P = points(io, 2);
        const subset = indices(io);
        const numIndices = io.integer();
        const fitter = new ApprHeightLine2();
        heightLine2(io, fitter, fitter.fit(P, subset, numIndices));
    }, { exact: true });

    family.case('ApprQuery.fitIndexed.direct', (io) => {
        const P = points(io, 2);
        const subset = indices(io);
        const probe = io.vec(2);
        const fitter = new ApprHeightLine2();
        const success = fitter.fitIndexed(P, subset);
        io.outBool(success);
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprQuery.ransac.tooFew', (io) => {
        const P = points(io, 2);
        const maxError = io.real();
        const numIterations = io.integer();
        const candidate = new ApprGaussian2();
        const best = new ApprGaussian2();
        const r = ApprQuery.ransac(candidate, P, 1, maxError, numIterations, best);
        io.outBool(r.success);
        io.outInt(r.bestConsensus.length);
        io.outVec(best.getParameters().center);
        io.outVec(best.getParameters().extent);
    }, { exact: true });

    family.case('ApprQuery.ransac.minimum', (io) => {
        const P = points(io, 2);
        const maxError = io.real();
        const numIterations = io.integer();
        const probe = io.vec(2);
        const candidate = new ApprGaussian2();
        const best = new ApprGaussian2();
        const r = ApprQuery.ransac(candidate, P, 2, maxError, numIterations, best);
        io.outBool(r.success);
        io.outInt(r.bestConsensus.length);
        for (const index of r.bestConsensus) { io.outInt(index); }
        io.outVec(best.getParameters().center);
        io.outVec(best.getParameters().axis[0]);
        io.outVec(best.getParameters().axis[1]);
        io.outVec(best.getParameters().extent);
        io.outReal(best.error(probe));
    }, { exact: true });

    // ---------------------------------------------------------- ApprCircle2

    family.case('ApprCircle2.fitUsingSquaredLengths', (io) => {
        const P = points(io, 2);
        const circle = new Hypersphere(2);
        const success = new ApprCircle2().fitUsingSquaredLengths(P, circle);
        io.outBool(success);
        io.outVec(circle.center);
        io.outReal(circle.radius);
    }, { exact: true });

    family.case('ApprCircle2.fitUsingSquaredLengths.cocircular', (io) => {
        const n = io.integer();
        const P = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) { P[i] = io.vec(2); }
        const circle = new Hypersphere(2);
        const success = new ApprCircle2().fitUsingSquaredLengths(P, circle);
        io.outBool(success);
        io.outVec(circle.center);
        io.outReal(circle.radius);
    }, { exact: true });

    family.case('ApprCircle2.fitUsingLengths', (io) => {
        const P = points(io, 2);
        const maxIterations = io.integer();
        const initialCenterIsAverage = io.boolean();
        const circle = Hypersphere.fromCenterRadius(io.vec(2), io.real());
        const iterations = new ApprCircle2().fitUsingLengths(P, maxIterations,
            initialCenterIsAverage, circle);
        io.outInt(iterations);
        io.outVec(circle.center);
        io.outReal(circle.radius);
    }, { exact: true });

    family.case('ApprCircle2.fitUsingLengths.epsilon', (io) => {
        const P = points(io, 2);
        const maxIterations = io.integer();
        const initialCenterIsAverage = io.boolean();
        const circle = Hypersphere.fromCenterRadius(io.vec(2), io.real());
        const epsilon = io.real();
        const iterations = new ApprCircle2().fitUsingLengths(P, maxIterations,
            initialCenterIsAverage, circle, epsilon);
        io.outInt(iterations);
        io.outVec(circle.center);
        io.outReal(circle.radius);
    }, { exact: true });

    // ---------------------------------------------------------- ApprSphere3

    family.case('ApprSphere3.fitUsingSquaredLengths', (io) => {
        const P = points(io, 3);
        const sphere = new Hypersphere(3);
        const success = new ApprSphere3().fitUsingSquaredLengths(P, sphere);
        io.outBool(success);
        io.outVec(sphere.center);
        io.outReal(sphere.radius);
    }, { exact: true });

    family.case('ApprSphere3.fitUsingSquaredLengths.cospherical', (io) => {
        const n = io.integer();
        const P = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) { P[i] = io.vec(3); }
        const sphere = new Hypersphere(3);
        const success = new ApprSphere3().fitUsingSquaredLengths(P, sphere);
        io.outBool(success);
        io.outVec(sphere.center);
        io.outReal(sphere.radius);
    }, { exact: true });

    family.case('ApprSphere3.fitUsingLengths', (io) => {
        const P = points(io, 3);
        const maxIterations = io.integer();
        const initialCenterIsAverage = io.boolean();
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const iterations = new ApprSphere3().fitUsingLengths(P, maxIterations,
            initialCenterIsAverage, sphere);
        io.outInt(iterations);
        io.outVec(sphere.center);
        io.outReal(sphere.radius);
    }, { exact: true });

    family.case('ApprSphere3.fitUsingLengths.epsilon', (io) => {
        const P = points(io, 3);
        const maxIterations = io.integer();
        const initialCenterIsAverage = io.boolean();
        const sphere = Hypersphere.fromCenterRadius(io.vec(3), io.real());
        const epsilon = io.real();
        const iterations = new ApprSphere3().fitUsingLengths(P, maxIterations,
            initialCenterIsAverage, sphere, epsilon);
        io.outInt(iterations);
        io.outVec(sphere.center);
        io.outReal(sphere.radius);
    }, { exact: true });

    // -------------------------------------------------------- ApprGaussian2

    function gaussian2(io: OracleIO, fitter: ApprGaussian2, probe: Vector): void {
        io.outVec(fitter.getParameters().center);
        io.outVec(fitter.getParameters().axis[0]);
        io.outVec(fitter.getParameters().axis[1]);
        io.outVec(fitter.getParameters().extent);
        io.outReal(fitter.error(probe));
    }

    family.case('ApprGaussian2.fit', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const fitter = new ApprGaussian2();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        gaussian2(io, fitter, probe);
    }, { exact: true });

    family.case('ApprGaussian2.fitIndexed', (io) => {
        const P = points(io, 2);
        const subset = indices(io);
        const probe = io.vec(2);
        const fitter = new ApprGaussian2();
        io.outBool(fitter.fitIndexed(P, subset));
        gaussian2(io, fitter, probe);
    }, { exact: true });

    family.case('ApprGaussian2.copyParameters', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const source = new ApprGaussian2();
        const target = new ApprGaussian2();
        source.fit(P);
        target.copyParameters(source);
        gaussian2(io, target, probe);
    }, { exact: true });

    // -------------------------------------------------------- ApprGaussian3

    function gaussian3(io: OracleIO, fitter: ApprGaussian3, probe: Vector): void {
        io.outVec(fitter.getParameters().center);
        io.outVec(fitter.getParameters().axis[0]);
        io.outVec(fitter.getParameters().axis[1]);
        io.outVec(fitter.getParameters().axis[2]);
        io.outVec(fitter.getParameters().extent);
        io.outReal(fitter.error(probe));
    }

    family.case('ApprGaussian3.fit', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const fitter = new ApprGaussian3();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        gaussian3(io, fitter, probe);
    }, { exact: true });

    family.case('ApprGaussian3.fitIndexed', (io) => {
        const P = points(io, 3);
        const subset = indices(io);
        const probe = io.vec(3);
        const fitter = new ApprGaussian3();
        io.outBool(fitter.fitIndexed(P, subset));
        gaussian3(io, fitter, probe);
    }, { exact: true });

    family.case('ApprGaussian3.copyParameters', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const source = new ApprGaussian3();
        const target = new ApprGaussian3();
        source.fit(P);
        target.copyParameters(source);
        gaussian3(io, target, probe);
    }, { exact: true });

    // ------------------------------------------------------ ApprHeightLine2

    family.case('ApprHeightLine2.fit', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const fitter = new ApprHeightLine2();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprHeightLine2.fit.verticalData', (io) => {
        const n = io.integer();
        const P = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) { P[i] = io.vec(2); }
        const probe = io.vec(2);
        const fitter = new ApprHeightLine2();
        io.outBool(fitter.fit(P));
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprHeightLine2.copyParameters', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const source = new ApprHeightLine2();
        const target = new ApprHeightLine2();
        source.fit(P);
        target.copyParameters(source);
        io.outVec(target.getParameters().average);
        io.outVec(target.getParameters().coefficients);
        io.outReal(target.error(probe));
    }, { exact: true });

    // ----------------------------------------------------- ApprHeightPlane3

    family.case('ApprHeightPlane3.fit', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const fitter = new ApprHeightPlane3();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprHeightPlane3.fit.collinearXY', (io) => {
        const n = io.integer();
        const P = new Array<Vector>(n);
        for (let i = 0; i < n; ++i) { P[i] = io.vec(3); }
        const probe = io.vec(3);
        const fitter = new ApprHeightPlane3();
        io.outBool(fitter.fit(P));
        io.outVec(fitter.getParameters().average);
        io.outVec(fitter.getParameters().coefficients);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprHeightPlane3.copyParameters', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const source = new ApprHeightPlane3();
        const target = new ApprHeightPlane3();
        source.fit(P);
        target.copyParameters(source);
        io.outVec(target.getParameters().average);
        io.outVec(target.getParameters().coefficients);
        io.outReal(target.error(probe));
    }, { exact: true });

    // -------------------------------------------------- ApprOrthogonalLine2

    family.case('ApprOrthogonalLine2.fit', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const fitter = new ApprOrthogonalLine2();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().direction);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalLine2.fitIndexed', (io) => {
        const P = points(io, 2);
        const subset = indices(io);
        const probe = io.vec(2);
        const fitter = new ApprOrthogonalLine2();
        io.outBool(fitter.fitIndexed(P, subset));
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().direction);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalLine2.copyParameters', (io) => {
        const P = points(io, 2);
        const probe = io.vec(2);
        const source = new ApprOrthogonalLine2();
        const target = new ApprOrthogonalLine2();
        source.fit(P);
        target.copyParameters(source);
        io.outVec(target.getParameters().origin);
        io.outVec(target.getParameters().direction);
        io.outReal(target.error(probe));
    }, { exact: true });

    // -------------------------------------------------- ApprOrthogonalLine3

    family.case('ApprOrthogonalLine3.fit', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const fitter = new ApprOrthogonalLine3();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().direction);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalLine3.fitIndexed', (io) => {
        const P = points(io, 3);
        const subset = indices(io);
        const probe = io.vec(3);
        const fitter = new ApprOrthogonalLine3();
        io.outBool(fitter.fitIndexed(P, subset));
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().direction);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalLine3.copyParameters', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const source = new ApprOrthogonalLine3();
        const target = new ApprOrthogonalLine3();
        source.fit(P);
        target.copyParameters(source);
        io.outVec(target.getParameters().origin);
        io.outVec(target.getParameters().direction);
        io.outReal(target.error(probe));
    }, { exact: true });

    // ------------------------------------------------- ApprOrthogonalPlane3

    family.case('ApprOrthogonalPlane3.fit', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const fitter = new ApprOrthogonalPlane3();
        io.outBool(fitter.fit(P));
        io.outInt(fitter.getMinimumRequired());
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().normal);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalPlane3.fitIndexed', (io) => {
        const P = points(io, 3);
        const subset = indices(io);
        const probe = io.vec(3);
        const fitter = new ApprOrthogonalPlane3();
        io.outBool(fitter.fitIndexed(P, subset));
        io.outVec(fitter.getParameters().origin);
        io.outVec(fitter.getParameters().normal);
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprOrthogonalPlane3.copyParameters', (io) => {
        const P = points(io, 3);
        const probe = io.vec(3);
        const source = new ApprOrthogonalPlane3();
        const target = new ApprOrthogonalPlane3();
        source.fit(P);
        target.copyParameters(source);
        io.outVec(target.getParameters().origin);
        io.outVec(target.getParameters().normal);
        io.outReal(target.error(probe));
    }, { exact: true });

    // ------------------------------------------------------ ApprPolynomial2

    function coefficients(io: OracleIO, values: readonly number[]): void {
        io.outInt(values.length);
        io.outReals(values);
    }

    family.case('ApprPolynomial2.fit', (io) => {
        const obs = observations(io, 2);
        const degree = io.integer();
        const probeX = io.real();
        const probe = io.reals(2);
        const fitter = new ApprPolynomial2(degree);
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReal(fitter.evaluate(probeX));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomial2.fitIndexed', (io) => {
        const obs = observations(io, 2);
        const subset = indices(io);
        const degree = io.integer();
        const probeX = io.real();
        const fitter = new ApprPolynomial2(degree);
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReal(fitter.evaluate(probeX));
    }, { exact: true });

    family.case('ApprPolynomial2.copyParameters', (io) => {
        const obs = observations(io, 2);
        const degree = io.integer();
        const probeX = io.real();
        const source = new ApprPolynomial2(degree);
        const target = new ApprPolynomial2(degree);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReal(target.evaluate(probeX));
    }, { exact: true });

    // ------------------------------------------------------ ApprPolynomial3

    family.case('ApprPolynomial3.fit', (io) => {
        const obs = observations(io, 3);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const probe = io.reals(3);
        const fitter = new ApprPolynomial3(xDegree, yDegree);
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReal(fitter.evaluate(probeX, probeY));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomial3.fitIndexed', (io) => {
        const obs = observations(io, 3);
        const subset = indices(io);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const fitter = new ApprPolynomial3(xDegree, yDegree);
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReal(fitter.evaluate(probeX, probeY));
    }, { exact: true });

    family.case('ApprPolynomial3.copyParameters', (io) => {
        const obs = observations(io, 3);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const source = new ApprPolynomial3(xDegree, yDegree);
        const target = new ApprPolynomial3(xDegree, yDegree);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReals(target.getYDomain());
        io.outReal(target.evaluate(probeX, probeY));
    }, { exact: true });

    // ------------------------------------------------------ ApprPolynomial4

    family.case('ApprPolynomial4.fit', (io) => {
        const obs = observations(io, 4);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const zDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const probe = io.reals(4);
        const fitter = new ApprPolynomial4(xDegree, yDegree, zDegree);
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReals(fitter.getZDomain());
        io.outReal(fitter.evaluate(probeX, probeY, probeZ));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomial4.fitIndexed', (io) => {
        const obs = observations(io, 4);
        const subset = indices(io);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const zDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const fitter = new ApprPolynomial4(xDegree, yDegree, zDegree);
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReals(fitter.getZDomain());
        io.outReal(fitter.evaluate(probeX, probeY, probeZ));
    }, { exact: true });

    family.case('ApprPolynomial4.copyParameters', (io) => {
        const obs = observations(io, 4);
        const xDegree = io.integer();
        const yDegree = io.integer();
        const zDegree = io.integer();
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const source = new ApprPolynomial4(xDegree, yDegree, zDegree);
        const target = new ApprPolynomial4(xDegree, yDegree, zDegree);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReal(target.evaluate(probeX, probeY, probeZ));
    }, { exact: true });

    // ----------------------------------------------- ApprPolynomialSpecial2

    family.case('ApprPolynomialSpecial2.fit', (io) => {
        const obs = observations(io, 2);
        const list = degrees(io);
        const probeX = io.real();
        const probe = io.reals(2);
        const fitter = new ApprPolynomialSpecial2(list);
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReal(fitter.evaluate(probeX));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomialSpecial2.fitIndexed', (io) => {
        const obs = observations(io, 2);
        const subset = indices(io);
        const list = degrees(io);
        const probeX = io.real();
        const fitter = new ApprPolynomialSpecial2(list);
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReal(fitter.evaluate(probeX));
    }, { exact: true });

    family.case('ApprPolynomialSpecial2.copyParameters', (io) => {
        const obs = observations(io, 2);
        const list = degrees(io);
        const probeX = io.real();
        const source = new ApprPolynomialSpecial2(list);
        const target = new ApprPolynomialSpecial2(list);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReal(target.evaluate(probeX));
    }, { exact: true });

    family.case('ApprPolynomialSpecial2.constructor.assert', (io) => {
        const count = io.integer();
        const list = new Array<number>(count);
        for (let i = 0; i < count; ++i) { list[i] = io.integer(); }
        const fitter = new ApprPolynomialSpecial2(list);
        io.outInt(fitter.getMinimumRequired());
    }, { exact: true });

    // ----------------------------------------------- ApprPolynomialSpecial3

    family.case('ApprPolynomialSpecial3.fit', (io) => {
        const obs = observations(io, 3);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const probe = io.reals(3);
        const fitter = new ApprPolynomialSpecial3(xDegrees,
            identityDegrees(xDegrees.length, 1));
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReal(fitter.evaluate(probeX, probeY));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomialSpecial3.fitIndexed', (io) => {
        const obs = observations(io, 3);
        const subset = indices(io);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const fitter = new ApprPolynomialSpecial3(xDegrees,
            identityDegrees(xDegrees.length, 2));
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReal(fitter.evaluate(probeX, probeY));
    }, { exact: true });

    family.case('ApprPolynomialSpecial3.copyParameters', (io) => {
        const obs = observations(io, 3);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const yDegrees = identityDegrees(xDegrees.length, 1);
        const source = new ApprPolynomialSpecial3(xDegrees, yDegrees);
        const target = new ApprPolynomialSpecial3(xDegrees, yDegrees);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReals(target.getYDomain());
        io.outReal(target.evaluate(probeX, probeY));
    }, { exact: true });

    family.case('ApprPolynomialSpecial3.constructor.affineAssert', (io) => {
        const count = io.integer();
        const xDegrees = new Array<number>(count);
        const yDegrees = new Array<number>(count);
        for (let i = 0; i < count; ++i) {
            xDegrees[i] = io.integer();
            yDegrees[i] = io.integer();
        }
        const fitter = new ApprPolynomialSpecial3(xDegrees, yDegrees);
        io.outInt(fitter.getMinimumRequired());
    }, { exact: true });

    // ----------------------------------------------- ApprPolynomialSpecial4

    family.case('ApprPolynomialSpecial4.fit', (io) => {
        const obs = observations(io, 4);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const probe = io.reals(4);
        const fitter = new ApprPolynomialSpecial4(xDegrees,
            identityDegrees(xDegrees.length, 1),
            identityDegrees(xDegrees.length, 2));
        io.outBool(fitter.fit(obs));
        io.outInt(fitter.getMinimumRequired());
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReals(fitter.getZDomain());
        io.outReal(fitter.evaluate(probeX, probeY, probeZ));
        io.outReal(fitter.error(probe));
    }, { exact: true });

    family.case('ApprPolynomialSpecial4.fitIndexed', (io) => {
        const obs = observations(io, 4);
        const subset = indices(io);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const fitter = new ApprPolynomialSpecial4(xDegrees,
            identityDegrees(xDegrees.length, 2),
            identityDegrees(xDegrees.length, 1));
        io.outBool(fitter.fitIndexed(obs, subset));
        coefficients(io, fitter.getParameters());
        io.outReals(fitter.getXDomain());
        io.outReals(fitter.getYDomain());
        io.outReals(fitter.getZDomain());
        io.outReal(fitter.evaluate(probeX, probeY, probeZ));
    }, { exact: true });

    family.case('ApprPolynomialSpecial4.copyParameters', (io) => {
        const obs = observations(io, 4);
        const xDegrees = degrees(io);
        const probeX = io.real();
        const probeY = io.real();
        const probeZ = io.real();
        const other = identityDegrees(xDegrees.length, 1);
        const source = new ApprPolynomialSpecial4(xDegrees, other, other);
        const target = new ApprPolynomialSpecial4(xDegrees, other, other);
        source.fit(obs);
        target.copyParameters(source);
        coefficients(io, target.getParameters());
        io.outReals(target.getXDomain());
        io.outReal(target.evaluate(probeX, probeY, probeZ));
    }, { exact: true });

    family.case('ApprPolynomialSpecial4.constructor.sizeAssert', (io) => {
        const count = io.integer();
        const extra = io.integer();
        const fitter = new ApprPolynomialSpecial4(identityDegrees(count, 1),
            identityDegrees(count, 1), identityDegrees(count + extra, 1));
        io.outInt(fitter.getMinimumRequired());
    }, { exact: true });

    family.finish();
});
