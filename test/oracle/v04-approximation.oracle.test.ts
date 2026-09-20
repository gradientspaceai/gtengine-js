// Replays oracle/cpp/cases/v04-approximation.cpp (verify group 4,
// approximation). Keep the two files in the same order.
import { describe } from 'vitest';
import { ApprParabola2 } from '../../src/ApprParabola2.js';
import { ApprParaboloid3 } from '../../src/ApprParaboloid3.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Mirrors Points2 / Points3 of the C++ case file: the count is recorded
// first, then that many points, in every generator mode.
function points(io: OracleIO, dimension: number): Vector[] {
    const n = io.integer();
    const list: Vector[] = [];
    for (let i = 0; i < n; ++i) { list.push(io.vec(dimension)); }
    return list;
}

describe('oracle: v04-approximation', () => {
    const family = new OracleFamily('v04-approximation');

    // ---- ApprParabola2 -----------------------------------------------------
    // A 3x3 closed-form inverse and a sum of products: arithmetic only.

    family.case('ApprParabola2.fit', (io) => {
        const r = ApprParabola2.fit(points(io, 2));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fitRobust', (io) => {
        const r = ApprParabola2.fitRobust(points(io, 2));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fit.throw', (io) => {
        const r = ApprParabola2.fit(points(io, 2));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParabola2.fitRobust.throw', (io) => {
        const r = ApprParabola2.fitRobust(points(io, 2));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    // ---- ApprParaboloid3 ---------------------------------------------------
    // An LDL^T decomposition of a 6x6 system: arithmetic only.

    family.case('ApprParaboloid3.fit', (io) => {
        const r = ApprParaboloid3.fit(points(io, 3));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fitRobust', (io) => {
        const r = ApprParaboloid3.fitRobust(points(io, 3));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fit.throw', (io) => {
        const r = ApprParaboloid3.fit(points(io, 3));
        io.outBool(r.success);
        io.outReals(r.u);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.case('ApprParaboloid3.fitRobust.throw', (io) => {
        const r = ApprParaboloid3.fitRobust(points(io, 3));
        io.outBool(r.success);
        io.outVec(r.average);
        io.outReals(r.v);
        io.outReal(r.meanSquareError);
    }, { exact: true });

    family.finish();
});
