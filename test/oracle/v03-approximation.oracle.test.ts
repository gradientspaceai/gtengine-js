// Replays oracle/cpp/cases/v03-approximation.cpp. Keep the two files in the
// same order.
//
// Every case is declared exact: no fitting path in this group reaches the C
// math library except through Math.sqrt (see the comment at the top of the
// C++ case file).
import { describe } from 'vitest';
import { ApprCircle2 } from '../../src/ApprCircle2.js';
import { ApprGaussian2 } from '../../src/ApprGaussian2.js';
import { ApprHeightLine2 } from '../../src/ApprHeightLine2.js';
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

    family.finish();
});
