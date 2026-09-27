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
    MinimumVolumeBox3FloatingPoint, type MinimumVolumeBox3FloatingPointCandidate
} from '../../src/MinimumVolumeBox3FloatingPoint.js';
import {
    MinimumVolumeBox3Rational, type MinimumVolumeBox3RationalCandidate
} from '../../src/MinimumVolumeBox3Rational.js';
import { BSNumber } from '../../src/BSNumber.js';
import type { OrientedBox3 } from '../../src/OrientedBox.js';
import { VETManifoldMesh } from '../../src/VETManifoldMesh.js';
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

// The C++ BoxIsConsistent: orthonormal axes and volume = 8 * e0 * e1 * e2,
// both to 1e-9.
function boxIsConsistent(box: OrientedBox3, volume: number): boolean {
    for (let i = 0; i < 3; ++i) {
        for (let j = 0; j < 3; ++j) {
            const a = box.axis[i].values, b = box.axis[j].values;
            const d = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
            if (Math.abs(d - (i === j ? 1 : 0)) > 1e-9) {
                return false;
            }
        }
    }
    const e = box.extent.values;
    const product = 8 * e[0] * e[1] * e[2];
    return Math.abs(product - volume) <= 1e-9 * Math.max(1, Math.abs(volume));
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

    // computeHull up to and including removeCoplanarTriangleAdjacencies, the
    // state in which getExtreme is called (the C++ Canonical::Prepare).
    prepare(vertices: readonly Vector[], indices: readonly number[], lgMaxSample: number): void {
        this.generateSubdivision(lgMaxSample);
        const mesh = new VETManifoldMesh();
        this.createMeshTopology(indices.length / 3, indices, mesh);
        this.extractMeshTopology(mesh);
        this.extractVertexAdjacencies(mesh);
        this.extractMeshGeometry(vertices);
        this.removeCoplanarTriangleAdjacencies();
    }

    extreme(direction: Vector): { vMax: number, dMax: number } {
        return this.getExtreme(direction);
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

    prepare(vertices: readonly Vector[], indices: readonly number[], lgMaxSample: number): void {
        this.generateSubdivision(lgMaxSample);
        const mesh = new VETManifoldMesh();
        this.createMeshTopology(indices.length / 3, indices, mesh);
        this.extractMeshTopology(mesh);
        this.extractVertexAdjacencies(mesh);
        this.extractMeshGeometry(vertices);
        this.removeCoplanarTriangleAdjacencies();
    }

    extreme(direction: Vector): { vMax: number, dMax: number } {
        const d = direction.values;
        const r = this.getExtreme([BSNumber.fromNumber(d[0]), BSNumber.fromNumber(d[1]),
            BSNumber.fromNumber(d[2])]);
        return { vMax: r.vMax, dMax: r.dMax.toNumber() };
    }
}

type ProbeState = ProbeFP['state'] | ProbeR['state'];

// The C++ MinimizerProbeFP: records the minimizer calls and can replace any
// minimizer by a no-op (bit mask 1 constantS, 2 constantT, 4 variableS,
// 8 variableT).
class MinimizerProbeFP extends ProbeFP {
    calls = [0, 0, 0, 0];
    argSum = [0, 0, 0, 0, 0, 0];
    processors = new Map<number, number>();

    constructor(readonly disabled: number) {
        super(0);
    }

    private hit(c: MinimumVolumeBox3FloatingPointCandidate): void {
        const p = c.levelCurveProcessorIndex;
        this.processors.set(p, (this.processors.get(p) ?? 0) + 1);
    }

    protected override minimizerConstantS(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        ++this.calls[0];
        this.hit(c);
        if ((this.disabled & 1) === 0) {
            super.minimizerConstantS(c, mvc);
        }
    }

    protected override minimizerConstantT(c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        ++this.calls[1];
        this.hit(c);
        if ((this.disabled & 2) === 0) {
            super.minimizerConstantT(c, mvc);
        }
    }

    protected override minimizerVariableS(sminNumer: number, smaxNumer: number, sDenom: number,
        c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        ++this.calls[2];
        this.hit(c);
        this.argSum[0] += sminNumer;
        this.argSum[1] += smaxNumer;
        this.argSum[2] += sDenom;
        if ((this.disabled & 4) === 0) {
            super.minimizerVariableS(sminNumer, smaxNumer, sDenom, c, mvc);
        }
    }

    protected override minimizerVariableT(tminNumer: number, tmaxNumer: number, tDenom: number,
        c: MinimumVolumeBox3FloatingPointCandidate,
        mvc: MinimumVolumeBox3FloatingPointCandidate): void {
        ++this.calls[3];
        this.hit(c);
        this.argSum[3] += tminNumer;
        this.argSum[4] += tmaxNumer;
        this.argSum[5] += tDenom;
        if ((this.disabled & 8) === 0) {
            super.minimizerVariableT(tminNumer, tmaxNumer, tDenom, c, mvc);
        }
    }
}

function outNVector(io: OracleIO, v: readonly BSNumber[]): void {
    for (let j = 0; j < 3; ++j) {
        io.outReal(v[j].toNumber());
    }
}

function outAlignedR(io: OracleIO, s: ProbeR['state']): void {
    const c = s.aligned;
    for (let i = 0; i < 3; ++i) {
        io.outInt(c.minSupportIndex[i]);
        io.outInt(c.maxSupportIndex[i]);
    }
    io.outReal(c.volume.toNumber());
}

// The winner's exact fields, each rounded to double once.
function outWinnerR(io: OracleIO, s: ProbeR['state']): void {
    const c = s.winner;
    io.outInt(c.edgeIndex[0]);
    io.outInt(c.edgeIndex[1]);
    io.outInt(c.levelCurveProcessorIndex);
    for (let i = 0; i < 2; ++i) {
        outNVector(io, c.N[i]);
        outNVector(io, c.M[i]);
    }
    io.outReal(c.f00.toNumber());
    io.outReal(c.f10.toNumber());
    io.outReal(c.f01.toNumber());
    io.outReal(c.f11.toNumber());
    for (let i = 0; i < 3; ++i) {
        outNVector(io, c.axis[i]);
        io.outInt(c.minSupportIndex[i]);
        io.outInt(c.maxSupportIndex[i]);
    }
    io.outReal(c.volume.toNumber());
}

// The C++ Dimension2Record: mode, lgMaxSample, n, the cloud, the diagnostic.
function dimension2Record(io: OracleIO,
    query: MinimumVolumeBox3FloatingPoint | MinimumVolumeBox3Rational): void {
    io.integer();
    const lgMaxSample = io.integer();
    io.integer();
    const points = readCloud(io);
    io.integer();  // diagnostic: the draw passed the probes
    const r = query.compute(points, lgMaxSample);
    io.outInt(r.dimension);
    outBox(io, r.box, r.volume);
}

// The C++ MinimizerProbeR. Upstream's MinimizerVariableT receives its
// arguments rounded to double (#355); the port's receives them exactly, and
// toNumber() is the same rounding, so the sums agree.
class MinimizerProbeR extends ProbeR {
    calls = [0, 0, 0, 0];
    argSum = [0, 0, 0, 0, 0, 0];
    processors = new Map<number, number>();

    constructor(readonly disabled: number) {
        super(0);
    }

    private hit(c: MinimumVolumeBox3RationalCandidate): void {
        const p = c.levelCurveProcessorIndex;
        this.processors.set(p, (this.processors.get(p) ?? 0) + 1);
    }

    protected override minimizerConstantS(c: MinimumVolumeBox3RationalCandidate,
        mvc: MinimumVolumeBox3RationalCandidate): void {
        ++this.calls[0];
        this.hit(c);
        if ((this.disabled & 1) === 0) {
            super.minimizerConstantS(c, mvc);
        }
    }

    protected override minimizerConstantT(c: MinimumVolumeBox3RationalCandidate,
        mvc: MinimumVolumeBox3RationalCandidate): void {
        ++this.calls[1];
        this.hit(c);
        if ((this.disabled & 2) === 0) {
            super.minimizerConstantT(c, mvc);
        }
    }

    protected override minimizerVariableS(sminNumer: BSNumber, smaxNumer: BSNumber,
        sDenom: BSNumber, c: MinimumVolumeBox3RationalCandidate,
        mvc: MinimumVolumeBox3RationalCandidate): void {
        ++this.calls[2];
        this.hit(c);
        this.argSum[0] += sminNumer.toNumber();
        this.argSum[1] += smaxNumer.toNumber();
        this.argSum[2] += sDenom.toNumber();
        if ((this.disabled & 4) === 0) {
            super.minimizerVariableS(sminNumer, smaxNumer, sDenom, c, mvc);
        }
    }

    protected override minimizerVariableT(tminNumer: BSNumber, tmaxNumer: BSNumber,
        tDenom: BSNumber, c: MinimumVolumeBox3RationalCandidate,
        mvc: MinimumVolumeBox3RationalCandidate): void {
        ++this.calls[3];
        this.hit(c);
        this.argSum[3] += tminNumer.toNumber();
        this.argSum[4] += tmaxNumer.toNumber();
        this.argSum[5] += tDenom.toNumber();
        if ((this.disabled & 8) === 0) {
            super.minimizerVariableT(tminNumer, tmaxNumer, tDenom, c, mvc);
        }
    }
}

// The C++ ReuseRecord: one minimizer-probe functor for mesh A, counters
// cleared, then mesh B.
function reuseRecord(io: OracleIO, query: MinimizerProbeFP | MinimizerProbeR): void {
    const lgMaxSample = io.integer();
    const a = readMesh(io);
    const b = readMesh(io);
    io.integer();  // diagnostic: the stale pairs changed upstream's box
    query.computeHull(a.vertices, a.indices, lgMaxSample);
    query.calls = [0, 0, 0, 0];
    query.argSum = [0, 0, 0, 0, 0, 0];
    query.processors.clear();
    const r = query.computeHull(b.vertices, b.indices, lgMaxSample);
    outMinimizerCalls(io, query);
    outBox(io, r.box, r.volume);
}

// The C++ InvalidArgumentRecord; every record throws.
function invalidArgumentRecord(io: OracleIO,
    query: MinimumVolumeBox3FloatingPoint | MinimumVolumeBox3Rational): void {
    const kind = io.integer();
    io.integer();
    const lgMaxSample = io.integer();
    if (kind <= 1) {
        const n = io.integer();
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            points.push(io.vec(3));
        }
        const r = query.compute(points, lgMaxSample);
        io.outInt(r.dimension);
        outBox(io, r.box, r.volume);
    } else {
        const { vertices, indices } = readMesh(io);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outBox(io, r.box, r.volume);
    }
}

function outMinimizerCalls(io: OracleIO,
    probe: { calls: number[], argSum: number[], processors: Map<number, number> }): void {
    for (let i = 0; i < 4; ++i) {
        io.outInt(probe.calls[i]);
    }
    io.outReals(probe.argSum);
    const keys = [...probe.processors.keys()].sort((a, b) => a - b);
    io.outInt(keys.length);
    for (const k of keys) {
        io.outInt(k);
        io.outInt(probe.processors.get(k)!);
    }
}

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
        io.integer();  // diagnostic: the exact-support separator holds

        const query = new MinimumVolumeBox3FloatingPoint(numThreads);
        const r = query.compute(points, lgMaxSample);
        io.outInt(r.dimension);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, points));
        io.outBool(boxIsConsistent(r.box, r.volume));
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.computeHull.canonical', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const numThreads = io.integer();
        const { vertices, indices } = readMesh(io);
        io.integer();  // diagnostic RawAgreement
        io.integer();  // diagnostic: the exact-support separator holds

        const query = new ProbeFP(numThreads);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, vertices));
        io.outBool(boxIsConsistent(r.box, r.volume));
        const s = query.state;
        outTopology(io, s);
        outAlignedFP(io, s);
        outWinnerFP(io, s);
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.getExtreme', (io) => {
        io.integer();
        io.integer();
        const { vertices, indices } = readMesh(io);
        const query = new ProbeFP(0);
        query.prepare(vertices, indices, 2);
        for (let k = 0; k < 4; ++k) {
            io.integer();
            const r = query.extreme(io.vec(3));
            io.outInt(r.vMax);
            io.outReal(r.dMax);
        }
    }, { exact: true });

    // Deviation (#426): upstream's climb stops short of the exact maximum by
    // more than the rounding bound; the port's plateau traversal does not.
    family.case('MinimumVolumeBox3FloatingPoint.getExtreme.plateau', (io) => {
        io.integer();
        const { vertices, indices } = readMesh(io);
        const direction = io.vec(3);
        io.integer();  // diagnostic: the shortfall exceeds the rounding bound
        const query = new ProbeFP(0);
        query.prepare(vertices, indices, 2);
        const r = query.extreme(direction);
        io.outInt(r.vMax);
        io.outReal(r.dMax);
    }, { exact: true, deviation: '#426 (GetExtreme stalls on a floating-point plateau)' });

    family.case('MinimumVolumeBox3FloatingPoint.minimizers', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const disabled = io.integer();
        const { vertices, indices } = readMesh(io);
        const query = new MinimizerProbeFP(disabled);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outMinimizerCalls(io, query);
        outBox(io, r.box, r.volume);
        outWinnerFP(io, query.state);
    }, { exact: true });

    family.case('MinimumVolumeBox3FloatingPoint.compute.dimension2', (io) => {
        dimension2Record(io, new MinimumVolumeBox3FloatingPoint(0));
    }, { exact: true });

    // Deviation (design): upstream's dimension-2 path uses the floating-point
    // MinimumAreaBox2<T, T>; the port's MinimumAreaBox2 is exact-only (port
    // note in src/MinimumAreaBox2.ts; v11 MinimumAreaBox2.deviation.floatComputeType).
    family.case('MinimumVolumeBox3FloatingPoint.compute.dimension2.floatComputeType', (io) => {
        dimension2Record(io, new MinimumVolumeBox3FloatingPoint(0));
    }, { exact: true, deviation: 'design: MinimumAreaBox2 is exact-only (src/MinimumAreaBox2.ts)' });

    // Deviation (#405, #426): a support vertex of upstream's winner falls
    // short of the exact extreme by more than the rounding bound, so the
    // port's confined fix replaces it.
    family.case('MinimumVolumeBox3FloatingPoint.computeHull.provenViolation', (io) => {
        const lgMaxSample = io.integer();
        io.integer();
        const { vertices, indices } = readMesh(io);
        io.integer();  // diagnostic: the search found such a mesh
        const r = new MinimumVolumeBox3FloatingPoint(0).computeHull(vertices, indices, lgMaxSample);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, vertices));
    }, { exact: true, deviation: '#405, #426 (proven support violations of the winner)' });

    family.case('MinimumVolumeBox3Rational.computeHull.canonical', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const numThreads = io.integer();
        const { vertices, indices } = readMesh(io);
        io.integer();  // diagnostic RawAgreement

        const query = new ProbeR(numThreads);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, vertices));
        io.outBool(boxIsConsistent(r.box, r.volume));
        const s = query.state;
        outTopology(io, s);
        outAlignedR(io, s);
        outWinnerR(io, s);
    }, { exact: true, timeout: 600000 });

    family.case('MinimumVolumeBox3Rational.compute.canonical', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const numThreads = io.integer();
        const points = readCloud(io);

        const query = new MinimumVolumeBox3Rational(numThreads);
        const r = query.compute(points, lgMaxSample);
        io.outInt(r.dimension);
        outBox(io, r.box, r.volume);
        io.outBool(contains(r.box, points));
        io.outBool(boxIsConsistent(r.box, r.volume));
    }, { exact: true, timeout: 600000 });

    family.case('MinimumVolumeBox3Rational.getExtreme', (io) => {
        io.integer();
        io.integer();
        const { vertices, indices } = readMesh(io);
        const query = new ProbeR(0);
        query.prepare(vertices, indices, 2);
        for (let k = 0; k < 4; ++k) {
            io.integer();
            const r = query.extreme(io.vec(3));
            io.outInt(r.vMax);
            io.outReal(r.dMax);
        }
    }, { exact: true, timeout: 600000 });

    family.case('MinimumVolumeBox3Rational.minimizers', (io) => {
        io.integer();
        const lgMaxSample = io.integer();
        io.integer();
        const disabled = 8 | io.integer();
        const { vertices, indices } = readMesh(io);
        const query = new MinimizerProbeR(disabled);
        const r = query.computeHull(vertices, indices, lgMaxSample);
        outMinimizerCalls(io, query);
        outBox(io, r.box, r.volume);
        outWinnerR(io, query.state);
    }, { exact: true, timeout: 600000 });

    family.case('MinimumVolumeBox3Rational.compute.dimension2', (io) => {
        dimension2Record(io, new MinimumVolumeBox3Rational(0));
    }, { exact: true });

    family.case('MinimumVolumeBox3Rational.compute.dimension2.floatComputeType', (io) => {
        dimension2Record(io, new MinimumVolumeBox3Rational(0));
    }, { exact: true, deviation: 'design: MinimumAreaBox2 is exact-only (src/MinimumAreaBox2.ts)' });

    // Deviation (new suspect, group report): upstream appends to
    // mEdgeIndices on every call, so a reused functor processes the previous
    // mesh's edge pairs first; the port resets them.
    family.case('MinimumVolumeBox3FloatingPoint.computeHull.reuse', (io) => {
        reuseRecord(io, new MinimizerProbeFP(0));
    }, { exact: true, deviation: 'v12 suspect: stale mEdgeIndices on functor reuse' });

    family.case('MinimumVolumeBox3Rational.computeHull.reuse', (io) => {
        reuseRecord(io, new MinimizerProbeR(8));
    }, {
        exact: true, timeout: 600000,
        deviation: 'v12 suspect: stale mEdgeIndices on functor reuse'
    });

    family.case('MinimumVolumeBox3FloatingPoint.invalidArgument', (io) => {
        invalidArgumentRecord(io, new MinimumVolumeBox3FloatingPoint(0));
    }, { exact: true });

    family.case('MinimumVolumeBox3Rational.invalidArgument', (io) => {
        invalidArgumentRecord(io, new MinimumVolumeBox3Rational(0));
    }, { exact: true });

    family.finish();
});
