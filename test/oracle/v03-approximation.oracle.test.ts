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

    family.finish();
});
