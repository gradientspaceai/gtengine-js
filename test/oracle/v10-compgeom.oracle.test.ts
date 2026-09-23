// Replays oracle/cpp/cases/v10-compgeom.cpp. Keep the two files in the same
// order. The C++ file documents what is comparable for each header.
import { describe } from 'vitest';
import { MinimumAreaCircle2 } from '../../src/MinimumAreaCircle2.js';
import { MinimumVolumeSphere3 } from '../../src/MinimumVolumeSphere3.js';
import { RotatingCalipers } from '../../src/RotatingCalipers.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// Read a recorded point list: the count first, then the coordinates.
function readPoints(io: OracleIO, dim: number): Vector[] {
    const count = io.integer();
    const points: Vector[] = [];
    for (let i = 0; i < count; ++i) {
        points.push(io.vec(dim));
    }
    return points;
}

function emitAntipodes(io: OracleIO, polygon: readonly Vector[]): void {
    const antipodes = RotatingCalipers.computeAntipodes(polygon);
    io.outInt(antipodes.length);
    for (const a of antipodes) {
        io.outInt(a.vertex);
        io.outInt(a.edge[0]);
        io.outInt(a.edge[1]);
    }
}

describe('oracle: v10-compgeom', () => {
    const family = new OracleFamily('v10-compgeom');

    // ---- RotatingCalipers ------------------------------------------------

    family.case('RotatingCalipers.computeAntipodes', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { exact: true });

    family.case('RotatingCalipers.computeAntipodes.collinear', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { exact: true });

    // The port compares each candidate corner against the most recent
    // *nonzero* edge, so a duplicated vertex no longer discards the next real
    // corner (issue #286). Upstream keeps one corner fewer, and throws when
    // that leaves fewer than three.
    family.case('RotatingCalipers.computeAntipodes.deviation.duplicate', (io) => {
        emitAntipodes(io, readPoints(io, 2));
    }, { deviation: '#286 (RotatingCalipers::CreatePolygon drops the corner after a duplicate)' });

    // ---- MinimumAreaCircle2 / MinimumVolumeSphere3 -----------------------

    // Only point sets whose result does not depend on the implementation
    // defined std::shuffle permutation reach this case; see the C++ file.
    function emitCircle(io: OracleIO, points: readonly Vector[]): void {
        const query = new MinimumAreaCircle2();
        const { minimal, success } = query.compute(points);
        io.outBool(success);
        io.outReal(minimal.center.get(0));
        io.outReal(minimal.center.get(1));
        io.outReal(minimal.radius);
        io.outInt(query.numSupport);
        const support = query.support.slice(0, query.numSupport).slice().sort(
            (a, b) => a - b);
        for (const s of support) {
            io.outInt(s);
        }
    }

    function emitSphere(io: OracleIO, points: readonly Vector[]): void {
        const query = new MinimumVolumeSphere3();
        const { minimal, success } = query.compute(points);
        io.outBool(success);
        io.outReal(minimal.center.get(0));
        io.outReal(minimal.center.get(1));
        io.outReal(minimal.center.get(2));
        io.outReal(minimal.radius);
        io.outInt(query.numSupport);
        const support = query.support.slice(0, query.numSupport).slice().sort(
            (a, b) => a - b);
        for (const s of support) {
            io.outInt(s);
        }
    }

    family.case('MinimumAreaCircle2.compute', (io) => {
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            points.push(io.vec(2));
        }
        emitCircle(io, points);
    }, { exact: true });

    // The port passes the whole input array to getContainerCircle2 where
    // upstream passes the unique-point count with the full array and so
    // bounds only a prefix (issue #286).
    family.case('MinimumAreaCircle2.compute.deviation.trappedFallback', (io) => {
        emitCircle(io, readPoints(io, 2));
    }, { deviation: '#286 (the trapped-failure fallback bounds only a prefix)' });

    family.case('MinimumAreaCircle2.compute.empty', (io) => {
        io.vec(2);
        const r = new MinimumAreaCircle2().compute([]);
        io.outBool(r.success);
    }, { exact: true });

    family.case('MinimumVolumeSphere3.compute', (io) => {
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            points.push(io.vec(3));
        }
        emitSphere(io, points);
    }, { exact: true });

    family.case('MinimumVolumeSphere3.compute.deviation.trappedFallback', (io) => {
        emitSphere(io, readPoints(io, 3));
    }, { deviation: '#286 (the trapped-failure fallback bounds only a prefix)' });

    family.case('MinimumVolumeSphere3.compute.empty', (io) => {
        io.vec(3);
        const r = new MinimumVolumeSphere3().compute([]);
        io.outBool(r.success);
    }, { exact: true });

    family.finish();
});
