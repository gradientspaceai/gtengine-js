// Replays oracle/cpp/cases/v12-compgeom.cpp. Keep the two files in the same
// order; each case reads the inputs the C++ case recorded, in the same order,
// and emits the outputs the C++ case recorded, in the same order.
//
// The C++ side pins the three unspecified enumeration orders of upstream's
// MinimumVolumeBox3 (hash-map edge and triangle numbering, hash-set adjacency
// order, ConvexHull3 triangle order) to sorted-key order, which is the order
// the port uses by construction. So the replay calls the port's public
// queries unchanged, and every case is exact. See the C++ file comment.
import { describe } from 'vitest';
import {
    MinimumVolumeBox3FloatingPoint
} from '../../src/MinimumVolumeBox3FloatingPoint.js';
import { MinimumVolumeBox3Rational } from '../../src/MinimumVolumeBox3Rational.js';
import type { OrientedBox3 } from '../../src/OrientedBox.js';
import type { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

function readCloud(io: OracleIO): Vector[] {
    const n = io.integer();
    const points: Vector[] = [];
    for (let i = 0; i < n; ++i) {
        points.push(io.vec(3));
    }
    return points;
}

function readMesh(io: OracleIO): { vertices: Vector[], indices: number[] } {
    const vertices = readCloud(io);
    const numIndices = io.integer();
    const indices: number[] = [];
    for (let i = 0; i < numIndices; ++i) {
        indices.push(io.integer());
    }
    return { vertices, indices };
}

// The same arithmetic as the C++ ContainmentViolation / PointScale.
function containmentViolation(box: OrientedBox3, points: readonly Vector[]): number {
    let worst = 0;
    for (const p of points) {
        const d = [0, 1, 2].map(j => p.values[j] - box.center.values[j]);
        for (let i = 0; i < 3; ++i) {
            const a = box.axis[i].values;
            const t = Math.abs(d[0] * a[0] + d[1] * a[1] + d[2] * a[2]) - box.extent.values[i];
            if (worst < t) {
                worst = t;
            }
        }
    }
    return worst;
}

function pointScale(points: readonly Vector[]): number {
    let scale = 1;
    for (const p of points) {
        for (let i = 0; i < 3; ++i) {
            const a = Math.abs(p.values[i]);
            if (scale < a) {
                scale = a;
            }
        }
    }
    return scale;
}

function contains(box: OrientedBox3, points: readonly Vector[]): boolean {
    return containmentViolation(box, points) <= 1e-9 * pointScale(points);
}

function outBox(io: OracleIO, box: OrientedBox3, volume: number): void {
    io.outVec(box.center);
    for (let i = 0; i < 3; ++i) {
        io.outVec(box.axis[i]);
    }
    io.outVec(box.extent);
    io.outReal(volume);
}

// Read access to the protected state the C++ Canonical<Base> exposes.
class ProbeFP extends MinimumVolumeBox3FloatingPoint {
    get state() {
        return {
            edges: this.mEdges, numEdgePairs: this.mEdgeIndices.length,
            climbStart: this.mVClimbStart, pool: this.mAdjacentPool,
            poolLocation: this.mAdjacentPoolLocation, maxSample: this.mMaxSample,
            domainIndex: this.mDomainIndex, aligned: this.mAlignedCandidate,
            winner: this.mMinimumVolumeObject
        };
    }
}

class ProbeR extends MinimumVolumeBox3Rational {
    get state() {
        return {
            edges: this.mEdges, numEdgePairs: this.mEdgeIndices.length,
            climbStart: this.mVClimbStart, pool: this.mAdjacentPool,
            poolLocation: this.mAdjacentPoolLocation, maxSample: this.mMaxSample,
            domainIndex: this.mDomainIndex, aligned: this.mAlignedCandidate,
            winner: this.mMinimumVolumeObject
        };
    }
}

type ProbeState = ProbeFP['state'] | ProbeR['state'];

function outTopology(io: OracleIO, s: ProbeState): void {
    io.outInt(s.edges.length);
    for (const e of s.edges) {
        io.outInt(e.v[0]);
        io.outInt(e.v[1]);
        io.outInt(e.t[0]);
        io.outInt(e.t[1]);
    }
    io.outInt(s.numEdgePairs);
    io.outInt(s.climbStart);
    for (let v = 0; v < s.poolLocation.length; ++v) {
        const base = s.poolLocation[v];
        const numAdjacent = s.pool[base];
        io.outInt(numAdjacent);
        for (let j = 1; j <= numAdjacent; ++j) {
            io.outInt(s.pool[base + j]);
        }
    }
    io.outInt(s.maxSample);
    for (const item of s.domainIndex) {
        io.outInt(item[0]);
        io.outInt(item[1]);
        io.outInt(item[2]);
    }
}

function outAlignedFP(io: OracleIO, s: ProbeFP['state']): void {
    const c = s.aligned;
    for (let i = 0; i < 3; ++i) {
        io.outInt(c.minSupportIndex[i]);
        io.outInt(c.maxSupportIndex[i]);
    }
    io.outReal(c.volume);
}

function outWinnerFP(io: OracleIO, s: ProbeFP['state']): void {
    const c = s.winner;
    io.outInt(c.edgeIndex[0]);
    io.outInt(c.edgeIndex[1]);
    io.outInt(c.levelCurveProcessorIndex);
    for (let i = 0; i < 2; ++i) {
        io.outVec(c.N[i]);
        io.outVec(c.M[i]);
    }
    io.outReal(c.f00);
    io.outReal(c.f10);
    io.outReal(c.f01);
    io.outReal(c.f11);
    for (let i = 0; i < 3; ++i) {
        io.outVec(c.axis[i]);
        io.outInt(c.minSupportIndex[i]);
        io.outInt(c.maxSupportIndex[i]);
    }
    io.outReal(c.volume);
}

describe('oracle: v12-compgeom', () => {
    const family = new OracleFamily('v12-compgeom');

    family.case('MinimumVolumeBox3FloatingPoint.compute.canonical', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const numThreads = io.integer();
        const points = readCloud(io);
        io.integer();  // diagnostic RawAgreement (hash-order effect on the raw query)

        const query = new MinimumVolumeBox3FloatingPoint(numThreads);
        const r = query.compute(points, lgMaxSample);
        io.outInt(r.dimension);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, points));
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.computeHull.canonical', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const numThreads = io.integer();
        const { vertices, indices } = readMesh(io);
        io.integer();  // diagnostic RawAgreement

        const query = new ProbeFP(numThreads);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, vertices));
        const s = query.state;
        outTopology(io, s);
        outAlignedFP(io, s);
        outWinnerFP(io, s);
    }, { exact: true });

    family.finish();
});
