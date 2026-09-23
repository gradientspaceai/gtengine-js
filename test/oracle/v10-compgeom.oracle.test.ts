// Replays oracle/cpp/cases/v10-compgeom.cpp. Keep the two files in the same
// order. The C++ file documents what is comparable for each header.
import { describe } from 'vitest';
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

    family.finish();
});
