// Replays oracle/cpp/cases/v27-imaging.cpp (AdaptiveSkeletonClimbing3).
// Keep the two files in the same order. Every case is arithmetic-only
// (+ - * /, sqrt, floor, truncating conversions) and compared bit for bit.
import { writeFileSync } from 'node:fs';
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
interface Checked {
    tag: string, voxels: number[], N: number, level: number, depth: number, r: ASC3Result
}
const checked: Checked[] = [];
// The case whose records are being replayed (set by the case wrappers).
let currentTag = '';

// GetZeroBase can return -1 for a CFG_MULT node that HasZeroSubedge admits,
// and Get*EdgesM would then interpolate at grid coordinate -1 (#194, latent).
// Count the interpolations at a negative coordinate the replays perform.
let negativeInterpolations = 0;
{
    const proto = AdaptiveSkeletonClimbing3.prototype as unknown as
        Record<string, (x: number, y: number, z: number) => number>;
    for (const name of ['getXInterp', 'getYInterp', 'getZInterp']) {
        const original = proto[name];
        proto[name] = function (this: unknown, x: number, y: number, z: number): number {
            if (x < 0 || y < 0 || z < 0) { ++negativeInterpolations; }
            return original.call(this, x, y, z);
        };
    }
}

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
        checked.push({ tag: currentTag, voxels, N, level: levels[k], depth: depths[k], r });
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

// ---- independent checks of the extractions (exact arithmetic) ----

// x as the integer x * 2^-E for E <= the exponent of every x passed.
function scaled(x: number, E: number): bigint {
    const { m, e } = decompose(x);
    return m << BigInt(e - E);
}

function minExponent(xs: readonly number[]): number {
    return Math.min(...xs.map((x) => decompose(x).e));
}

interface CheckStats {
    extractions: number;
    edgeVertices: number;
    branchPoints: number;
    centroids: number;
    degenerate: number;
    closedMeshes: number;
    openMeshes: number;
    notClosed: string[];
    orientedMeshes: number;
    misoriented: string[];
    bad: string[];
}

// Checks one extraction; the level must be finite and differ from every
// voxel value (the documented precondition). Every vertex of the unique mesh
// is (a) on a unit grid edge whose end values straddle the level, within
// (a + 4) 2^-53 of the exact root a + (L - f0)/(f1 - f0) (the rounding of
// L - f0, the quotient and the sum; checked exactly), or (b) a plus-sign
// branch point: on a four-crossing unit face whose bilinear saddle value is
// exactly the level, at (the u of the v = 0 crossing, the v of the u = 0
// crossing), or (c) a fan centroid strictly inside a merged box. No
// triangle is degenerate (exact cross product). When the level set does not
// reach the image border (all border voxels on one side) every mesh edge is
// shared by exactly two triangles and, after OrientTriangles, traversed in
// opposite directions by them.
function checkExtraction(c: Checked, stats: CheckStats): void {
    const { voxels, N, level, r } = c;
    const size = (1 << N) + 1;
    const F = (q: readonly number[]): number => voxels[q[0] + size * (q[1] + size * q[2])];
    ++stats.extractions;
    const where = `N=${N} level=${level} depth=${c.depth}`;
    const merged = r.boxes.filter((b) => b[3] > 1 || b[4] > 1 || b[5] > 1);
    for (const p of r.vu) {
        if (onCrossedEdge(p, F, size, level)) { ++stats.edgeVertices; continue; }
        if (isBranchPoint(p, F, size, level)) { ++stats.branchPoints; continue; }
        if (merged.some((b) => [0, 1, 2].every((k) => b[k] < p[k] && p[k] < b[k] + b[k + 3]))) {
            ++stats.centroids;
            continue;
        }
        stats.bad.push(`${where}: vertex (${p.join(', ')}) is on no crossed edge, no saddle face, in no merged box`);
    }
    for (const t of r.tu) {
        const a = r.vu[t.V[0]], b = r.vu[t.V[1]], d = r.vu[t.V[2]];
        const E = minExponent([...a, ...b, ...d]);
        const e1 = [0, 1, 2].map((k) => scaled(b[k], E) - scaled(a[k], E));
        const e2 = [0, 1, 2].map((k) => scaled(d[k], E) - scaled(a[k], E));
        const n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2], e1[0] * e2[1] - e1[1] * e2[0]];
        if (n.every((x) => x === 0n)) { ++stats.degenerate; }
    }
    if (r.tu.length === 0) { return; }
    // Undirected and directed edge uses.
    const undirected = new Map<string, number>();
    const directed = new Map<string, number>();
    for (const t of r.tu) {
        for (let k = 0; k < 3; ++k) {
            const a = t.V[k], b = t.V[(k + 1) % 3];
            const key = `${Math.min(a, b)},${Math.max(a, b)}`;
            undirected.set(key, (undirected.get(key) ?? 0) + 1);
        }
    }
    for (const t of r.to) {
        for (let k = 0; k < 3; ++k) {
            const key = `${t.V[k]},${t.V[(k + 1) % 3]}`;
            directed.set(key, (directed.get(key) ?? 0) + 1);
        }
    }
    let borderSides = 0;
    for (let z = 0; z < size; ++z) {
        for (let y = 0; y < size; ++y) {
            for (let x = 0; x < size; ++x) {
                if ([x, y, z].some((q) => q === 0 || q === size - 1)) {
                    borderSides |= F([x, y, z]) > level ? 1 : 2;
                }
            }
        }
    }
    if (borderSides !== 3) {
        const counts = [...undirected.values()];
        if (counts.every((n) => n === 2)) {
            ++stats.closedMeshes;
        } else {
            stats.notClosed.push(`${where}: edge uses ${[...new Set(counts)].join('/')}`);
        }
        const flipped = [...directed.keys()].filter((key) => {
            const [a, b] = key.split(',');
            return directed.get(key) !== 1 || directed.get(`${b},${a}`) !== 1;
        });
        if (flipped.length === 0) {
            ++stats.orientedMeshes;
        } else {
            stats.misoriented.push(`${where}: ${flipped.length} directed edges`);
        }
    } else {
        ++stats.openMeshes;
    }
}

// (a) of checkExtraction: p lies on a unit grid edge whose end values
// straddle the level, within (a + 4) 2^-53 of the exact root.
function onCrossedEdge(p: Vertex, F: (q: readonly number[]) => number, size: number,
    level: number): boolean {
    for (let k = 0; k < 3; ++k) {
        const others = [0, 1, 2].filter((m) => m !== k);
        if (!others.every((m) => Number.isInteger(p[m]))) { continue; }
        const lows = Number.isInteger(p[k]) ? [p[k] - 1, p[k]] : [Math.floor(p[k])];
        for (const a of lows) {
            if (a < 0 || a + 1 > size - 1) { continue; }
            const q0 = [p[0], p[1], p[2]], q1 = [p[0], p[1], p[2]];
            q0[k] = a;
            q1[k] = a + 1;
            const f0 = F(q0), f1 = F(q1);
            if ((f0 > level) === (f1 > level)) { continue; }
            // |(c - a)(f1 - f0) - (L - f0)| <= (a + 4) 2^-53 |f1 - f0|
            const E = Math.min(minExponent([p[k], level]), -53);
            const lhs = (scaled(p[k], E) - (BigInt(a) << BigInt(-E))) * BigInt(f1 - f0)
                - (scaled(level, E) - (BigInt(f0) << BigInt(-E)));
            const bound = (BigInt(a + 4) << BigInt(-E - 53)) * BigInt(Math.abs(f1 - f0));
            if ((lhs < 0n ? -lhs : lhs) <= bound) { return true; }
        }
    }
    return false;
}

// (b) of checkExtraction.
function isBranchPoint(p: Vertex, F: (q: readonly number[]) => number, size: number,
    level: number): boolean {
    for (let axis = 0; axis < 3; ++axis) {
        if (!Number.isInteger(p[axis])) { continue; }
        const ku = axis === 0 ? 1 : 0, kv = axis === 2 ? 1 : 2;
        const i = Math.floor(p[ku]), j = Math.floor(p[kv]);
        if (i < 0 || j < 0 || i + 1 > size - 1 || j + 1 > size - 1) { continue; }
        const corner = (du: number, dv: number): number => {
            const q = [0, 0, 0];
            q[axis] = p[axis]; q[ku] = i + du; q[kv] = j + dv;
            return F(q);
        };
        const f00 = corner(0, 0), f10 = corner(1, 0), f01 = corner(0, 1), f11 = corner(1, 1);
        const g = [f00, f10, f01, f11].map((f) => f > level);
        if (!(g[0] === g[3] && g[1] === g[2] && g[0] !== g[1])) { continue; }
        if (saddleSign(f00, f10, f01, f11, level) !== 0) { continue; }
        if (p[ku] === i + (level - f00) / (f10 - f00) && p[kv] === j + (level - f00) / (f01 - f00)) {
            return true;
        }
    }
    return false;
}

function tagged(tag: string, body: (io: OracleIO) => void): (io: OracleIO) => void {
    return (io) => {
        currentTag = tag;
        body(io);
    };
}

describe('oracle: v27-imaging', () => {
    const family = new OracleFamily('v27-imaging');
    const exact = { exact: true };

    family.case('AdaptiveSkeletonClimbing3.extract', tagged('extract', ascCase), exact);
    family.case('AdaptiveSkeletonClimbing3.extract.types', tagged('types', typesCase), exact);
    family.case('AdaptiveSkeletonClimbing3.extract.saddle', tagged('saddle', ascCase), exact);
    // Port fix of the face pairing (#544 in 3-D; see the group report and
    // the header of src/AdaptiveSkeletonClimbing3.ts). Every record has a
    // four-crossing face on which upstream contradicts the interpolant: the
    // C++ record's own SaddleCheck counts it, the port's is zero.
    family.case('AdaptiveSkeletonClimbing3.extract.saddlePairing', tagged('saddlePairing', ascCase),
        { exact: true, deviation: '#544 (3-D face cases), v27 report' });
    family.case('AdaptiveSkeletonClimbing3.extract.large', tagged('large', largeCase),
        { exact: true, timeout: 600000 });
    family.case('AdaptiveSkeletonClimbing3.extract.rootMonobox', tagged('rootMonobox', ascCase), exact);
    family.case('AdaptiveSkeletonClimbing3.meshOps', meshOpsCase, exact);
    family.case('AdaptiveSkeletonClimbing3.invalid', invalidCase, exact);

    it('AdaptiveSkeletonClimbing3 meshes against the trilinear interpolant (independent checks)', () => {
        const stats: CheckStats = {
            extractions: 0, edgeVertices: 0, branchPoints: 0, centroids: 0, degenerate: 0,
            closedMeshes: 0, openMeshes: 0, notClosed: [], orientedMeshes: 0, misoriented: [], bad: []
        };
        const faces = new Map<string, [number, number]>();
        let rootDropped = 0;
        let firstBad: unknown;
        const byDepth = new Map<string, number[]>();
        for (const c of checked) {
            const f = faces.get(c.tag) ?? [0, 0];
            f[0] += c.r.saddle[0];
            f[1] += c.r.saddle[1];
            faces.set(c.tag, f);
            const size = (1 << c.N) + 1;
            const sample = c.voxels.includes(c.level);
            if (c.tag === 'rootMonobox') {
                // The level set is not empty (the level is strictly inside
                // the voxel range of a monotone image), the mesh is.
                const crossed = c.voxels.some((v) => v > c.level) && c.voxels.some((v) => v <= c.level);
                if (crossed && c.r.tu.length === 0 && c.r.boxes.length === 0) { ++rootDropped; }
            }
            if (!Number.isFinite(c.level) || sample) { continue; }
            const before = stats.notClosed.length;
            const badBefore = stats.bad.length;
            checkExtraction(c, stats);
            if (firstBad === undefined && stats.bad.length > badBefore) {
                firstBad = { tag: c.tag, voxels: c.voxels, level: c.level, depth: c.depth, boxes: c.r.boxes, vu: c.r.vu,
                    tu: c.r.tu.map((t) => [...t.V]) };
            }
            const key = c.depth > c.N ? 'depth>N' : (c.depth <= 0 ? 'depth<=0' : 'depth1..N');
            const d = byDepth.get(key) ?? [0, 0];
            d[0] += 1;
            d[1] += stats.notClosed.length - before;
            byDepth.set(key, d);
            void size;
        }
        const report = JSON.stringify({
            ...stats, notClosed: stats.notClosed.length, misoriented: stats.misoriented.length,
            bad: stats.bad.length, faces: Object.fromEntries(faces), rootDropped,
            byDepth: Object.fromEntries(byDepth), negativeInterpolations,
            notClosedSample: stats.notClosed.slice(0, 5), misorientedSample: stats.misoriented.slice(0, 5),
            badSample: stats.bad.slice(0, 5), firstBad
        }, null, 1);
        if (process.env['V27_STATS'] !== undefined) { writeFileSync(process.env['V27_STATS'], report); }
    });

    family.finish();
});
