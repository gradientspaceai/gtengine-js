// Replays oracle/cpp/cases/v27-imaging.cpp (AdaptiveSkeletonClimbing3).
// Keep the two files in the same order. Every case is arithmetic-only
// (+ - * /, sqrt, floor, truncating conversions) and compared bit for bit.
import { describe, expect, it } from 'vitest';
import {
    AdaptiveSkeletonClimbing3,
    type AdaptiveSkeletonClimbing3Vertex
} from '../../src/AdaptiveSkeletonClimbing3.js';
import { TriangleKey } from '../../src/TriangleKey.js';
import { OracleFamily, type OracleIO } from './harness.js';

type Vertex = AdaptiveSkeletonClimbing3Vertex;

// FNV-1a over 32-bit words, as the C++ Digest.
class Digest {
    h = 2166136261;

    word(w: number): void {
        for (let i = 0; i < 4; ++i) {
            this.h = (this.h ^ ((w >>> (8 * i)) & 0xFF)) >>> 0;
            this.h = Math.imul(this.h, 16777619) >>> 0;
        }
    }

    integer(v: number): void { this.word(v >>> 0); }

    real(x: number): void {
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(0, x, true);
        this.word(view.getUint32(0, true));
        this.word(view.getUint32(4, true));
    }
}

// Lists of at most this many elements are emitted in full, longer ones as a
// digest (kFullLimit in the case file).
const FULL_LIMIT = 24;

function outVertices(io: OracleIO, v: readonly Vertex[]): void {
    io.outInt(v.length);
    if (v.length <= FULL_LIMIT) {
        for (const p of v) { io.outReal(p[0]); io.outReal(p[1]); io.outReal(p[2]); }
    } else {
        const d = new Digest();
        for (const p of v) { d.real(p[0]); d.real(p[1]); d.real(p[2]); }
        io.outInt(d.h);
    }
}

function outTriangles(io: OracleIO, t: readonly TriangleKey[]): void {
    io.outInt(t.length);
    if (t.length <= FULL_LIMIT) {
        for (const k of t) { io.outInt(k.V[0]); io.outInt(k.V[1]); io.outInt(k.V[2]); }
    } else {
        const d = new Digest();
        for (const k of t) { d.integer(k.V[0]); d.integer(k.V[1]); d.integer(k.V[2]); }
        io.outInt(d.h);
    }
}

type Box = [number, number, number, number, number, number];

function outBoxes(io: OracleIO, boxes: readonly Box[]): void {
    io.outInt(boxes.length);
    if (boxes.length <= FULL_LIMIT) {
        for (const b of boxes) { for (const c of b) { io.outInt(c); } }
    } else {
        const d = new Digest();
        for (const b of boxes) { for (const c of b) { d.integer(c); } }
        io.outInt(d.h);
    }
}

// The monoboxes of the last extract() (upstream reads them back through
// PrintBoxes; the port keeps them in the private mBoxes).
interface OctBoxView { x0: number, y0: number, z0: number, dx: number, dy: number, dz: number }
function getBoxes(asc: AdaptiveSkeletonClimbing3): Box[] {
    const boxes = (asc as unknown as { mBoxes: OctBoxView[] }).mBoxes;
    expect(boxes.length).toBe(asc.getNumBoxes());
    return boxes.map((b) => [b.x0, b.y0, b.z0, b.dx, b.dy, b.dz]);
}

function copyVertices(v: readonly Vertex[]): Vertex[] {
    return v.map((p) => [p[0], p[1], p[2]]);
}

function copyTriangles(t: readonly TriangleKey[]): TriangleKey[] {
    return t.map((k) => {
        const c = new TriangleKey(true);
        c.V[0] = k.V[0];
        c.V[1] = k.V[1];
        c.V[2] = k.V[2];
        return c;
    });
}

interface ASC3Result {
    boxes: Box[];
    v: Vertex[];
    t: TriangleKey[];
    vu: Vertex[];
    tu: TriangleKey[];
    to: TriangleKey[];
    normals: Vertex[];
    saddle: [number, number];
}

function outResult(io: OracleIO, r: ASC3Result): void {
    outBoxes(io, r.boxes);
    outVertices(io, r.v);
    outTriangles(io, r.t);
    outVertices(io, r.vu);
    outTriangles(io, r.tu);
    outTriangles(io, r.to);
    outVertices(io, r.normals);
    io.outInt(r.saddle[0]);
    io.outInt(r.saddle[1]);
}

// Extractions whose outputs the independent checks examine.
interface Checked { voxels: number[], N: number, level: number, depth: number, r: ASC3Result }
const checked: Checked[] = [];

// A finite double as m * 2^e with integer m (exact).
function decompose(x: number): { m: bigint, e: number } {
    const view = new DataView(new ArrayBuffer(8));
    view.setFloat64(0, x);
    const hi = view.getUint32(0), lo = view.getUint32(4);
    const biased = (hi >>> 20) & 0x7FF;
    let m = (BigInt(hi & 0xFFFFF) << 32n) | BigInt(lo);
    let e: number;
    if (biased === 0) {
        e = -1074;
    } else {
        m |= 1n << 52n;
        e = biased - 1075;
    }
    return { m: (hi >>> 31) !== 0 ? -m : m, e };
}

// The exact sign of dg = (f00-L)(f11-L) - (f01-L)(f10-L) for integer corners
// and a finite double level, expanded directly (not through the port's
// det - L*S grouping): with L = m 2^e and s = max(0, -e), 4^s dg is the
// product difference of the integers 2^s f - 2^s L.
function saddleSign(f00: number, f10: number, f01: number, f11: number, level: number): number {
    const { m, e } = decompose(level);
    // Scale by 2^s with s = max(0, -e) so that L 2^s = m 2^(e+s) is an integer.
    const s = Math.max(0, -e);
    const L = e + s >= 0 ? m << BigInt(e + s) : m;
    const k = (f: number): bigint => (BigInt(f) << BigInt(s)) - L;
    const dg = k(f00) * k(f11) - k(f01) * k(f10);
    return dg > 0n ? 1 : (dg < 0n ? -1 : 0);
}

// The independent face-saddle check (SaddleCheck in the case file): every
// four-crossing unit face (finite level, no corner on the level) must be
// cut the way the trilinear interpolant restricted to the face (bilinear)
// is cut. Returns [faces checked, faces contradicting the interpolant].
function saddleCheck(voxels: readonly number[], size: number, level: number,
    vu: readonly Vertex[], tu: readonly TriangleKey[]): [number, number] {
    const result: [number, number] = [0, 0];
    if (!Number.isFinite(level)) { return result; }
    const key = (p: readonly number[]): string => `${p[0]},${p[1]},${p[2]}`;
    const index = new Map<string, number>();
    vu.forEach((p, i) => { if (!index.has(key(p))) { index.set(key(p), i); } });
    const edges = new Set<string>();
    const adjacent = vu.map(() => new Set<number>());
    for (const t of tu) {
        for (let k = 0; k < 3; ++k) {
            const a = t.V[k], b = t.V[(k + 1) % 3];
            edges.add(`${Math.min(a, b)},${Math.max(a, b)}`);
            adjacent[a].add(b);
            adjacent[b].add(a);
        }
    }
    const find = (p: number[]): number => index.get(key(p)) ?? -1;
    const hasEdge = (a: number, b: number): boolean =>
        a >= 0 && b >= 0 && edges.has(`${Math.min(a, b)},${Math.max(a, b)}`);
    const P = (q: number[]): number => voxels[q[0] + size * (q[1] + size * q[2])];
    for (let axis = 0; axis < 3; ++axis) {
        const ku = axis === 0 ? 1 : 0, kv = axis === 2 ? 1 : 2;
        for (let a = 0; a < size; ++a) {
            for (let j = 0; j + 1 < size; ++j) {
                for (let i = 0; i + 1 < size; ++i) {
                    const corner = (du: number, dv: number): number => {
                        const q = [0, 0, 0];
                        q[axis] = a; q[ku] = i + du; q[kv] = j + dv;
                        return P(q);
                    };
                    const f00 = corner(0, 0), f10 = corner(1, 0), f01 = corner(0, 1), f11 = corner(1, 1);
                    const g = [f00, f10, f01, f11].map((f) => f > level);
                    if (!(g[0] === g[3] && g[1] === g[2] && g[0] !== g[1])) { continue; }
                    if ([f00, f10, f01, f11].some((f) => f === level)) { continue; }
                    const crossing = (k: number, base: number, h0: number, h1: number,
                        du: number, dv: number): number => {
                        const p = [0, 0, 0];
                        p[axis] = a; p[ku] = i + du; p[kv] = j + dv;
                        p[k] = base + (level - h0) / (h1 - h0);
                        return find(p);
                    };
                    const ev0 = crossing(ku, i, f00, f10, 0, 0);
                    const ev1 = crossing(ku, i, f01, f11, 0, 1);
                    const eu0 = crossing(kv, j, f00, f01, 0, 0);
                    const eu1 = crossing(kv, j, f10, f11, 1, 0);
                    const sdg = saddleSign(f00, f10, f01, f11, level);
                    let right: boolean;
                    if (sdg > 0) {
                        right = hasEdge(ev0, eu1) && hasEdge(ev1, eu0);
                    } else if (sdg < 0) {
                        right = hasEdge(ev0, eu0) && hasEdge(ev1, eu1);
                    } else {
                        right = ev0 >= 0 && ev1 >= 0 && eu0 >= 0 && eu1 >= 0
                            && [...adjacent[ev0]].some((c) => adjacent[ev1].has(c)
                                && adjacent[eu0].has(c) && adjacent[eu1].has(c));
                    }
                    ++result[0];
                    if (!right) { ++result[1]; }
                }
            }
        }
    }
    return result;
}

function runASC(N: number, voxels: number[], fixBoundary: boolean, levels: readonly number[],
    depths: readonly number[], sameDirs: readonly boolean[]): ASC3Result[] {
    const asc = new AdaptiveSkeletonClimbing3(N, voxels, fixBoundary);
    const results: ASC3Result[] = [];
    for (let k = 0; k < levels.length; ++k) {
        const { vertices, triangles } = asc.extract(levels[k], depths[k]);
        const boxes = getBoxes(asc);
        const vu = copyVertices(vertices);
        const tu = copyTriangles(triangles);
        asc.makeUnique(vu, tu);
        const to = copyTriangles(tu);
        asc.orientTriangles(vu, to, sameDirs[k]);
        const normals = asc.computeNormals(vu, to);
        const saddle = saddleCheck(voxels, (1 << N) + 1, levels[k], vu, tu);
        const r = { boxes, v: vertices, t: triangles, vu, tu, to, normals, saddle };
        results.push(r);
        checked.push({ voxels, N, level: levels[k], depth: depths[k], r });
    }
    return results;
}

// The record layout of the extract, types, saddle, saddlePairing and
// rootMonobox cases: N, fixBoundary, the number of extractions, the voxels,
// then (level, depth, sameDir) per extraction.
function ascCase(io: OracleIO): void {
    const N = io.integer();
    const fixBoundary = io.integer() === 0;
    const numExtracts = io.integer();
    const size = (1 << N) + 1;
    const voxels: number[] = [];
    for (let i = 0; i < size * size * size; ++i) { voxels.push(io.integer()); }
    const levels: number[] = [], depths: number[] = [], sameDirs: boolean[] = [];
    for (let k = 0; k < numExtracts; ++k) {
        levels.push(io.real());
        depths.push(io.integer());
        sameDirs.push(io.boolean());
    }
    for (const r of runASC(N, voxels, fixBoundary, levels, depths, sameDirs)) { outResult(io, r); }
}

function typesCase(io: OracleIO): void {
    io.integer();  // the voxel type; the port has one number type
    ascCase(io);
}

// LargeVoxel of the case file (all values stay far below 2^53).
function largeVoxel(kind: number, p: readonly number[], x: number, y: number, z: number): number {
    const u = x - p[9], v = y - p[10], w = z - p[11];
    const au = Math.abs(u), av = Math.abs(v), aw = Math.abs(w);
    switch (kind) {
        case 0:
            return p[0] * u * u + p[1] * v * v + p[2] * w * w + p[3] * u * v + p[4] * v * w
                + p[5] * w * u + p[6] * x + p[7] * y + p[8] * z;
        case 1:
            return p[0] * au + p[1] * av + p[2] * aw + p[3] * Math.min(au, av)
                + p[4] * Math.max(av, aw) + p[6] * x + p[7] * y + p[8] * z;
        case 2: {
            const s = u * u + v * v + w * w + p[0] * p[0] - p[1] * p[1];
            const a = p[2] === 2 ? v : u, b = p[2] === 0 ? v : w;
            return (s * s - 4 * p[0] * p[0] * (a * a + b * b)) * p[5];
        }
        case 3: {
            const s0 = u * u + v * v + w * w - p[0];
            const s1 = (x - p[6]) ** 2 + (y - p[7]) ** 2 + (z - p[8]) ** 2 - p[1];
            return p[5] * Math.min(s0, s1);
        }
        default:
            return p[0] * u * v * w + p[1] * u * u + p[2] * v * v + p[3] * w * w
                + p[6] * x + p[7] * y + p[8] * z;
    }
}

function largeCase(io: OracleIO): void {
    const N = io.integer();
    const size = (1 << N) + 1;
    const kind = io.integer();
    const p: number[] = [];
    for (let k = 0; k < 12; ++k) { p.push(io.integer()); }
    const level = io.real();
    const depth = io.integer();
    const fixBoundary = io.integer() === 0;
    const sameDir = io.boolean();
    const voxels: number[] = new Array<number>(size * size * size);
    for (let z = 0; z < size; ++z) {
        for (let y = 0; y < size; ++y) {
            for (let x = 0; x < size; ++x) {
                // + 0 turns a -0 product into +0 (the C++ values are int64_t).
                voxels[x + size * (y + size * z)] = largeVoxel(kind, p, x, y, z) + 0;
            }
        }
    }
    const [r] = runASC(N, voxels, fixBoundary, [level], [depth], [sameDir]);
    outResult(io, r);
}

// MakeUnique, OrientTriangles and ComputeNormals on arbitrary meshes (the
// meshOps case).
function meshOpsCase(io: OracleIO): void {
    const N = io.integer();
    const size = (1 << N) + 1;
    const voxels: number[] = [];
    for (let i = 0; i < size * size * size; ++i) { voxels.push(io.integer()); }
    const numVertices = io.integer();
    const vertices: Vertex[] = [];
    for (let i = 0; i < numVertices; ++i) {
        const kind = io.integer();
        if (kind === 0 && i > 0) {
            const j = io.integer();
            const c = vertices[j].map((x) => (x === 0 ? -x : x));
            vertices.push([c[0], c[1], c[2]]);
        } else {
            vertices.push([io.real(), io.real(), io.real()]);
        }
    }
    const numTriangles = io.integer();
    const triangles: TriangleKey[] = [];
    for (let t = 0; t < numTriangles; ++t) {
        const kind = io.integer();
        const key = new TriangleKey(true);
        if (kind === 0 && t > 0) {
            const j = io.integer();
            const rotation = io.integer();
            for (let k = 0; k < 3; ++k) { key.V[k] = triangles[j].V[(k + rotation) % 3]; }
        } else {
            for (let k = 0; k < 3; ++k) { key.V[k] = io.integer(); }
        }
        triangles.push(key);
    }
    const sameDir = io.boolean();
    const asc = new AdaptiveSkeletonClimbing3(N, voxels);
    const vu = copyVertices(vertices);
    const tu = copyTriangles(triangles);
    asc.makeUnique(vu, tu);
    const to = copyTriangles(tu);
    asc.orientTriangles(vu, to, sameDir);
    const normals = asc.computeNormals(vu, to);
    const ta = copyTriangles(triangles);
    asc.orientTriangles(vertices, ta, !sameDir);
    const na = asc.computeNormals(vertices, ta);
    outVertices(io, vu);
    outTriangles(io, tu);
    outTriangles(io, to);
    outVertices(io, normals);
    outTriangles(io, ta);
    outVertices(io, na);
}

function invalidCase(io: OracleIO): void {
    const N = io.integer();
    const useNull = io.boolean();
    const size = (1 << N) + 1;
    const voxels: number[] = [];
    for (let i = 0; i < size * size * size; ++i) { voxels.push(io.integer()); }
    const asc = new AdaptiveSkeletonClimbing3(N,
        useNull ? (null as unknown as number[]) : voxels);
    const { vertices, triangles } = asc.extract(0.5, -1);
    outBoxes(io, getBoxes(asc));
    outVertices(io, vertices);
    outTriangles(io, triangles);
}

describe('oracle: v27-imaging', () => {
    const family = new OracleFamily('v27-imaging');
    const exact = { exact: true };

    family.case('AdaptiveSkeletonClimbing3.extract', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing3.extract.types', typesCase, exact);
    family.case('AdaptiveSkeletonClimbing3.extract.saddle', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing3.extract.saddlePairing', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing3.extract.large', largeCase, { exact: true, timeout: 600000 });
    family.case('AdaptiveSkeletonClimbing3.extract.rootMonobox', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing3.meshOps', meshOpsCase, exact);
    family.case('AdaptiveSkeletonClimbing3.invalid', invalidCase, exact);

    family.finish();
});
