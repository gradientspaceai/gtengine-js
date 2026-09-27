// Replays oracle/cpp/cases/v16-core.cpp. Keep the two files in the same order.
import { describe } from 'vitest';
import { CurveExtractorEdge, CurveExtractorVertex } from '../../src/CurveExtractor.js';
import { CurveExtractorSquares } from '../../src/CurveExtractorSquares.js';
import { CurveExtractorTriangles } from '../../src/CurveExtractorTriangles.js';
import { IEEEBinary16 } from '../../src/IEEEBinary16.js';
import { ImplicitSurface3 } from '../../src/ImplicitSurface3.js';
import { MeshCurvature } from '../../src/MeshCurvature.js';
import { MeshSmoother } from '../../src/MeshSmoother.js';
import { PolygonTreeEx, PolygonTreeExNode } from '../../src/PolygonTree.js';
import { PolygonWindingOrder } from '../../src/PolygonWindingOrder.js';
import { Vector, add, mul } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ------------------------------------------------------------ CurveExtractor

interface ImageDraw {
    xBound: number;
    yBound: number;
    mode: number;
    pixels: number[];
    level: number;
}

function drawImage(io: OracleIO): ImageDraw {
    const mode = io.integer();
    const xBound = io.integer();
    const yBound = io.integer();
    const pixels: number[] = [];
    for (let i = 0; i < xBound * yBound; ++i) { pixels.push(io.integer()); }
    const level = io.integer();
    return { xBound, yBound, mode, pixels, level };
}

function emitRational(io: OracleIO, vertices: CurveExtractorVertex[],
    edges: CurveExtractorEdge[]): void {
    io.outInt(vertices.length);
    for (const v of vertices) {
        io.outInt(v.xNumer); io.outInt(v.xDenom);
        io.outInt(v.yNumer); io.outInt(v.yDenom);
    }
    io.outInt(edges.length);
    for (const e of edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
}

// Independent check (agreement is not correctness): every rational vertex
// of the raw extraction lies on the level set of the interpolant the
// extractor contours, the bilinear interpolant of the square (squares) or
// the linear interpolant of the triangle containing it (triangles, split
// along the diagonal that alternates with the square's parity). Evaluated
// exactly with bigint rationals.
function checkOnLevelSet(d: ImageDraw, vertices: CurveExtractorVertex[],
    triangles: boolean): void {
    const f = (x: number, y: number): bigint =>
        BigInt(d.pixels[x + d.xBound * y] - d.level);
    for (const v of vertices) {
        const xn = BigInt(v.xNumer), xd = BigInt(v.xDenom);
        const yn = BigInt(v.yNumer), yd = BigInt(v.yDenom);
        if (xd <= 0n || yd <= 0n) { throw new Error('nonpositive denominator'); }
        let i = Number(xn / xd), j = Number(yn / yd);
        i = Math.min(Math.max(i, 0), d.xBound - 2);
        j = Math.min(Math.max(j, 0), d.yBound - 2);
        // u = (xn - i*xd)/xd, v = (yn - j*yd)/yd; F scaled by xd*yd.
        const un = xn - BigInt(i) * xd, vn = yn - BigInt(j) * yd;
        const f00 = f(i, j), f10 = f(i + 1, j), f01 = f(i, j + 1), f11 = f(i + 1, j + 1);
        let value: bigint;
        if (!triangles) {
            value = (xd - un) * (yd - vn) * f00 + un * (yd - vn) * f10
                + (xd - un) * vn * f01 + un * vn * f11;
        } else if ((i & 1) === (j & 1)) {
            // Diagonal (i+1,j)-(i,j+1): u + v <= 1 is the triangle of (i,j).
            if (un * yd + vn * xd <= xd * yd) {
                value = f00 * xd * yd + un * yd * (f10 - f00) + vn * xd * (f01 - f00);
            } else {
                value = f11 * xd * yd + (xd - un) * yd * (f01 - f11)
                    + (yd - vn) * xd * (f10 - f11);
            }
        } else {
            // Diagonal (i,j)-(i+1,j+1): v >= u is the triangle of (i,j+1).
            if (vn * xd >= un * yd) {
                value = f00 * xd * yd + vn * xd * (f01 - f00) + un * yd * (f11 - f01);
            } else {
                value = f00 * xd * yd + un * yd * (f10 - f00) + vn * xd * (f11 - f10);
            }
        }
        if (value !== 0n) {
            throw new Error(`vertex ${xn}/${xd}, ${yn}/${yd} is not on the level set`);
        }
    }
}

function runExtractor(io: OracleIO, triangles: boolean): void {
    const d = drawImage(io);
    const make = () => triangles
        ? new CurveExtractorTriangles(d.xBound, d.yBound, d.pixels)
        : new CurveExtractorSquares(d.xBound, d.yBound, d.pixels);
    const extractor = make();
    const { vertices, edges } = extractor.extract(d.level);
    checkOnLevelSet(d, vertices, triangles);
    emitRational(io, vertices, edges);
    extractor.makeUnique(vertices, edges);
    emitRational(io, vertices, edges);
    for (let remove = 0; remove < 2; ++remove) {
        const r = extractor.extractReal(d.level, remove !== 0);
        io.outInt(r.vertices.length);
        for (const v of r.vertices) { io.outReal(v[0]); io.outReal(v[1]); }
        io.outInt(r.edges.length);
        for (const e of r.edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
    }
}

// ------------------------------------------------------------ IEEEBinary16

const f32 = new Float32Array(1);
const u32 = new Uint32Array(f32.buffer);

function floatOfBits(bits: number): number {
    u32[0] = bits >>> 0;
    return f32[0];
}

function bitsOfFloat(value: number): number {
    f32[0] = value;
    return u32[0];
}

// Independent model of the binary32 -> binary16 conversion for non-NaN
// inputs: round to nearest, ties to even, computed on the binary64 value
// (exact for every binary32) by scaling to an integer multiple of the half
// quantum. NaNs follow finding #110: the high 10 trailing bits are kept, so
// a payload below 2^13 becomes infinity.
function referenceConvert32To16(bits: number): number {
    const sign = (bits >>> 16) & 0x8000;
    if ((bits & 0x7F800000) === 0x7F800000) {
        return sign | 0x7C00 | ((bits & 0x007FFFFF) >>> 13);
    }
    const x = Math.abs(floatOfBits(bits));
    const roundEven = (q: number): number => {
        const f = Math.floor(q);
        const r = q - f;
        return (r > 0.5 || (r === 0.5 && (f % 2) === 1)) ? f + 1 : f;
    };
    if (x < 2 ** -14) {
        return sign | roundEven(x * 2 ** 24);
    }
    let e = Math.floor(Math.log2(x));
    if (2 ** e > x) { --e; }
    if (2 ** (e + 1) <= x) { ++e; }
    if (e > 15) { return sign | 0x7C00; }
    const m = roundEven(x * 2 ** (10 - e));
    const enc = ((e + 15) << 10) + (m - 1024);
    return sign | Math.min(enc, 0x7C00);
}

function outHalf(io: OracleIO, h: IEEEBinary16, tol?: number): void {
    io.outReal(h.number, tol);
}

// ------------------------------------------------------------ PolygonWindingOrder

// Exact simplicity test for an integer polygon: distinct vertices, no two
// nonadjacent edges touching, no adjacent edges folding back on each other.
function isSimpleLatticePolygon(polygon: Vector[]): boolean {
    const n = polygon.length;
    const x = polygon.map(v => v.get(0)), y = polygon.map(v => v.get(1));
    const orient = (a: number, b: number, c: number): number =>
        Math.sign((x[b] - x[a]) * (y[c] - y[a]) - (y[b] - y[a]) * (x[c] - x[a]));
    const onSegment = (a: number, b: number, c: number): boolean =>
        Math.min(x[a], x[b]) <= x[c] && x[c] <= Math.max(x[a], x[b])
        && Math.min(y[a], y[b]) <= y[c] && y[c] <= Math.max(y[a], y[b]);
    const intersect = (a: number, b: number, c: number, d: number): boolean => {
        const o1 = orient(a, b, c), o2 = orient(a, b, d);
        const o3 = orient(c, d, a), o4 = orient(c, d, b);
        if (o1 * o2 < 0 && o3 * o4 < 0) { return true; }
        return (o1 === 0 && onSegment(a, b, c)) || (o2 === 0 && onSegment(a, b, d))
            || (o3 === 0 && onSegment(c, d, a)) || (o4 === 0 && onSegment(c, d, b));
    };
    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            if (x[i] === x[j] && y[i] === y[j]) { return false; }
        }
    }
    for (let i = 0; i < n; ++i) {
        const i1 = (i + 1) % n;
        for (let j = i + 1; j < n; ++j) {
            const j1 = (j + 1) % n;
            if (j === i1 || i === j1) {
                // Adjacent edges share one vertex; they must not fold back.
                const shared = (j === i1 ? i1 : i);
                const p = (shared === i1 ? i : i1), q = (shared === i1 ? j1 : j);
                if (orient(shared, p, q) === 0 && (x[p] - x[shared]) * (x[q] - x[shared])
                    + (y[p] - y[shared]) * (y[q] - y[shared]) > 0) { return false; }
                continue;
            }
            if (intersect(i, i1, j, j1)) { return false; }
        }
    }
    return true;
}

// ------------------------------------------------------------ PolygonTree

interface TreeDraw {
    points: Vector[];
    tree: PolygonTreeEx;
}

function drawTree(io: OracleIO): TreeDraw {
    const numPoints = io.integer();
    const points: Vector[] = [];
    for (let i = 0; i < numPoints; ++i) { points.push(io.vec(2)); }
    const numNodes = io.integer();
    const parent: number[] = [0];
    for (let k = 1; k < numNodes; ++k) { parent.push(io.integer()); }
    const tree = new PolygonTreeEx();
    for (let k = 0; k < numNodes; ++k) {
        const node = new PolygonTreeExNode();
        node.self = k;
        node.parent = (k === 0 ? PolygonTreeEx.INVALID : parent[k]);
        node.minChild = numNodes;
        node.supChild = numNodes;
        node.chirality = io.integer();
        const numTriangles = io.integer();
        for (let t = 0; t < numTriangles; ++t) {
            node.triangulation.push([io.integer(), io.integer(), io.integer()]);
        }
        tree.nodes.push(node);
    }
    for (let k = numNodes - 1; k >= 1; --k) {
        const p = tree.nodes[parent[k]];
        p.minChild = k;
        if (p.supChild === numNodes) { p.supChild = k + 1; }
    }
    for (const node of tree.nodes) {
        if (node.minChild === numNodes) { node.supChild = node.minChild; }
    }
    tree.nodes.forEach((node, k) => {
        for (const t of node.triangulation) {
            tree.insideTriangles.push(t);
            tree.insideNodeIndices.push(k);
        }
    });
    return { points, tree };
}

function outIndexOrInvalid(io: OracleIO, i: number): void {
    io.outInt(i === PolygonTreeEx.INVALID ? -1 : i);
}

// ------------------------------------------------------------ surface meshes

interface SurfaceMesh {
    kind: number;
    vertices: Vector[];
    indices: number[];
}

function readSurfaceMesh(io: OracleIO): SurfaceMesh {
    const kind = io.integer();
    const nv = io.integer();
    const vertices: Vector[] = [];
    for (let i = 0; i < nv; ++i) { vertices.push(io.vec(3)); }
    const ni = io.integer();
    const indices: number[] = [];
    for (let i = 0; i < ni; ++i) { indices.push(io.integer()); }
    return { kind, vertices, indices };
}

// The replay of the case file's WeightedSmoother.
class WeightedSmoother extends MeshSmoother {
    protected override vertexInfluenced(i: number, t: number): boolean {
        return (i % 3 !== 1) || t > 0.75;
    }
    protected override getTangentWeight(_i: number, t: number): number {
        return 0.25 + 0.125 * t;
    }
    protected override getNormalWeight(i: number, t: number): number {
        return (i % 2 === 0 ? -0.0625 : 0.03125) * t;
    }
}

// ------------------------------------------------------------ ImplicitSurface3

// The replay of the case file's CubicSurface (same grouping).
class CubicSurface extends ImplicitSurface3 {
    constructor(private readonly c: number[]) { super(); }
    f(p: Vector): number {
        const c = this.c, x = p.get(0), y = p.get(1), z = p.get(2);
        return ((((((((((c[0] * (x * x) + c[1] * (y * y)) + c[2] * (z * z))
            + c[3] * (x * y)) + c[4] * (y * z)) + c[5] * (x * z)) + c[6] * x)
            + c[7] * y) + c[8] * z) + c[9]) + c[10] * ((x * y) * z));
    }
    fx(p: Vector): number {
        const c = this.c;
        return (((((2.0 * c[0]) * p.get(0) + c[3] * p.get(1)) + c[5] * p.get(2)) + c[6])
            + c[10] * (p.get(1) * p.get(2)));
    }
    fy(p: Vector): number {
        const c = this.c;
        return (((((2.0 * c[1]) * p.get(1) + c[3] * p.get(0)) + c[4] * p.get(2)) + c[7])
            + c[10] * (p.get(0) * p.get(2)));
    }
    fz(p: Vector): number {
        const c = this.c;
        return (((((2.0 * c[2]) * p.get(2) + c[4] * p.get(1)) + c[5] * p.get(0)) + c[8])
            + c[10] * (p.get(0) * p.get(1)));
    }
    fxx(_p: Vector): number { return 2.0 * this.c[0]; }
    fxy(p: Vector): number { return this.c[3] + this.c[10] * p.get(2); }
    fxz(p: Vector): number { return this.c[5] + this.c[10] * p.get(1); }
    fyy(_p: Vector): number { return 2.0 * this.c[1]; }
    fyz(p: Vector): number { return this.c[4] + this.c[10] * p.get(0); }
    fzz(_p: Vector): number { return 2.0 * this.c[2]; }
}

describe('oracle: v16-core', () => {
    const family = new OracleFamily('v16-core');

    family.case('CurveExtractorSquares.extract', (io) => runExtractor(io, false),
        { exact: true });

    // Every record plants a '+000' or '00+0' square, for which upstream
    // emits the two edges incident to the nonzero corner; the port emits the
    // zero edges (and passes checkOnLevelSet on every record).
    family.case('CurveExtractorSquares.extract.threeZeroCorners',
        (io) => runExtractor(io, false),
        { exact: true, deviation: 'group 16 finding: CurveExtractorSquares +000/00+0 edges swapped' });

    family.case('CurveExtractorTriangles.extract', (io) => runExtractor(io, true),
        { exact: true });

    family.case('IEEEBinary16.convert32To16', (io) => {
        for (let k = 0; k < 4; ++k) {
            const bits = io.integer();
            const h = IEEEBinary16.convert32To16(bits);
            if (h !== referenceConvert32To16(bits)) {
                throw new Error(`convert32To16(${bits.toString(16)}) = ${h.toString(16)} `
                    + 'disagrees with the round-to-nearest-even model');
            }
            io.outInt(h);
            if (((bits & 0x7FFFFFFF) >>> 0) <= 0x7F800000) {
                io.outInt(IEEEBinary16.fromNumber(floatOfBits(bits)).encoding);
            }
        }
    }, { exact: true });

    family.case('IEEEBinary16.convert16To32', (io) => {
        for (let k = 0; k < 8; ++k) {
            const h = io.integer();
            const bits = IEEEBinary16.convert16To32(h);
            io.outInt(bits);
            if (IEEEBinary16.convert32To16(bits) !== h) {
                throw new Error(`half ${h.toString(16)} does not round-trip`);
            }
            // 'number' is the port of operator float()/operator double().
            // A NaN has no reproducible bits once it is a JS number, so a
            // NaN reports the conversion's bits (compared above).
            const d = new IEEEBinary16(h).number;
            io.outInt(Number.isNaN(d) ? bits : bitsOfFloat(d));
            if (!Number.isNaN(d)) { io.outReal(d); }
        }
    }, { exact: true });

    family.case('IEEEBinary16.fields', (io) => {
        const x = new IEEEBinary16(io.integer());
        io.outInt(x.encoding);
        io.outInt(x.getSign());
        io.outInt(x.getBiased());
        io.outInt(x.getTrailing());
        const enc = x.getEncoding();
        io.outInt(enc.sign);
        io.outInt(enc.biased);
        io.outInt(enc.trailing);
        io.outInt(x.getClassification());
        io.outBool(x.isZero());
        io.outBool(x.isSignMinus());
        io.outBool(x.isSubnormal());
        io.outBool(x.isNormal());
        io.outBool(x.isFinite());
        io.outBool(x.isInfinite());
        io.outBool(x.isNaN());
        io.outBool(x.isQuietNaN());
        io.outBool(x.isSignalingNaN());
        io.outInt(x.getNextUp());
        io.outInt(x.getNextDown());
        const s = io.integer();
        const b = io.integer();
        const t = io.integer();
        const y = new IEEEBinary16(0);
        y.setEncoding(s, b, t);
        io.outInt(y.encoding);
    }, { exact: true });

    family.case('IEEEBinary16.fromNumber', (io) => {
        const d = io.real();
        io.outInt(IEEEBinary16.fromNumber(d).encoding);
    }, { exact: true });

    family.case('IEEEBinary16.compare', (io) => {
        const x = new IEEEBinary16(io.integer());
        const y = new IEEEBinary16(io.integer());
        io.outBool(IEEEBinary16.equals(x, y));
        io.outBool(IEEEBinary16.notEquals(x, y));
        io.outBool(IEEEBinary16.lessThan(x, y));
        io.outBool(IEEEBinary16.lessThanOrEqual(x, y));
        io.outBool(IEEEBinary16.greaterThan(x, y));
        io.outBool(IEEEBinary16.greaterThanOrEqual(x, y));
    }, { exact: true });

    family.case('IEEEBinary16.arithmetic', (io) => {
        const a = io.integer();
        const b = io.integer();
        const f = io.real();
        const H = IEEEBinary16;
        const x = new H(a), y = new H(b);
        io.outInt(H.negate(x).encoding);
        io.outReals([H.add(x, y), H.sub(x, y), H.mul(x, y), H.div(x, y),
            H.add(x, f), H.sub(x, f), H.mul(x, f), H.div(x, f),
            H.add(f, x), H.sub(f, x), H.mul(f, x), H.div(f, x)]);
        const u = [0, 1, 2, 3, 4, 5, 6, 7].map(() => new H(a));
        H.addAssign(u[0], y); H.subAssign(u[1], y); H.mulAssign(u[2], y); H.divAssign(u[3], y);
        H.addAssign(u[4], f); H.subAssign(u[5], f); H.mulAssign(u[6], f); H.divAssign(u[7], f);
        for (const v of u) {
            io.outBool(v.isNaN());
            if (!v.isNaN()) { io.outInt(v.encoding); }
        }
    }, { exact: true });

    family.case('IEEEBinary16.mathExact', (io) => {
        const H = IEEEBinary16;
        const x = new H(io.integer());
        const y = new H(io.integer());
        const z = new H(io.integer());
        const exponent = io.integer();
        outHalf(io, H.ceil(x));
        outHalf(io, H.floor(x));
        outHalf(io, H.fabs(x));
        outHalf(io, H.sqrt(x));
        outHalf(io, H.fmod(x, y));
        const fr = H.frexp(x);
        outHalf(io, fr.result);
        if (x.isFinite()) { io.outInt(fr.exponent); }
        outHalf(io, H.ldexp(x, exponent));
        const yz = H.lessThan(y, z);
        outHalf(io, H.clamp(x, yz ? y : z, yz ? z : y));
        io.outInt(H.isign(x));
        outHalf(io, H.saturate(x));
        outHalf(io, H.sign(x));
        outHalf(io, H.sqr(x));
        outHalf(io, H.invsqrt(x));
    }, { exact: true });

    // Binary32 libm (MSVC) against binary64 V8 libm rounded to binary32:
    // one binary32 ulp moves the half result by at most one half ulp.
    family.case('IEEEBinary16.mathLibm', (io) => {
        const H = IEEEBinary16;
        const x = new H(io.integer());
        const y = new H(io.integer());
        const r = [H.acos(x), H.acosh(x), H.asin(x), H.asinh(x), H.atan(x), H.atanh(x),
            H.atan2(x, y), H.cos(x), H.cosh(x), H.exp(x), H.exp2(x), H.log(x), H.log2(x),
            H.log10(x), H.pow(x, y), H.sin(x), H.sinh(x), H.tan(x), H.tanh(x),
            H.atandivpi(x), H.atan2divpi(x, y), H.cospi(x), H.sinpi(x), H.exp10(x)];
        for (const v of r) { outHalf(io, v, 1e-3); }
    });

    // ------------------------------------------------------------ PolygonWindingOrder

    family.case('PolygonWindingOrder.operator', (io) => {
        const mode = io.integer();
        const n = io.integer();
        const polygon: Vector[] = [];
        if (mode === 3) {
            const base = io.vec(2);
            const dir = io.vec(2);
            for (let i = 0; i < n; ++i) {
                const k = io.real();
                polygon.push(add(base, mul(dir, k)));
            }
        } else {
            for (let i = 0; i < n; ++i) { polygon.push(io.vec(2)); }
        }
        const ccw = new PolygonWindingOrder().isCounterClockwise(polygon);
        // Independent check on the simple (star-shaped) lattice polygons:
        // when the turn at the lexicographically smallest vertex is strict,
        // the answer is the sign of the shoelace area, computed exactly.
        if (mode === 0) {
            let area2 = 0n;
            for (let i = 0, j = n - 1; i < n; j = i++) {
                area2 += BigInt(polygon[j].get(0)) * BigInt(polygon[i].get(1))
                    - BigInt(polygon[i].get(0)) * BigInt(polygon[j].get(1));
            }
            let ll = 0;
            for (let i = 1; i < n; ++i) { if (polygon[i].lessThan(polygon[ll])) { ll = i; } }
            const p = polygon[ll], a = polygon[(ll + 1) % n], b = polygon[(ll + n - 1) % n];
            const turn = (a.get(0) - p.get(0)) * (b.get(1) - p.get(1))
                - (a.get(1) - p.get(1)) * (b.get(0) - p.get(0));
            if (turn !== 0 && isSimpleLatticePolygon(polygon) && (area2 > 0n) !== ccw) {
                throw new Error('winding order disagrees with the shoelace area');
            }
        }
        io.outBool(ccw);
    }, { exact: true });

    // ------------------------------------------------------------ PolygonTree

    family.case('PolygonTree.getContainingTriangle', (io) => {
        const d = drawTree(io);
        const test = io.vec(2);
        const chirality = io.integer();
        const r0 = d.tree.getContainingTriangle(test, d.points);
        const r1 = d.tree.getContainingTriangleInList(test, d.tree.insideTriangles,
            d.tree.insideNodeIndices, d.points);
        const r2 = PolygonTreeEx.getContainingTriangleWithChirality(test,
            d.tree.insideTriangles, chirality, d.points);
        outIndexOrInvalid(io, r0.nIndex);
        outIndexOrInvalid(io, r0.tIndex);
        outIndexOrInvalid(io, r1.nIndex);
        outIndexOrInvalid(io, r1.tIndex);
        outIndexOrInvalid(io, r2);
    }, { exact: true });

    family.case('PolygonTree.getContainingTriangle.invalidArgument', (io) => {
        const d = drawTree(io);
        const test = io.vec(2);
        const drop = io.integer();
        const nodeIndices = [...d.tree.insideNodeIndices];
        if (drop === 1) { nodeIndices.push(0); }
        const r = d.tree.getContainingTriangleInList(test, d.tree.insideTriangles,
            nodeIndices, d.points);
        outIndexOrInvalid(io, r.nIndex);
        outIndexOrInvalid(io, r.tIndex);
    }, { exact: true });

    // ------------------------------------------------------------ MeshSmoother

    family.case('MeshSmoother.update', (io) => {
        const m = readSurfaceMesh(io);
        const weighted = io.boolean();
        const smoother = weighted ? new WeightedSmoother() : new MeshSmoother();
        const vertices = m.vertices.map(v => v.clone());
        smoother.initialize(vertices, m.indices);
        io.outInt(smoother.getNumVertices());
        io.outInt(smoother.getNumTriangles());
        for (const c of smoother.getNeighborCounts()) { io.outInt(c); }
        for (let step = 0; step < 3; ++step) {
            smoother.update(0.5 * step);
            for (const v of smoother.getNormals()) { io.outVec(v); }
            for (const v of smoother.getMeans()) { io.outVec(v); }
            for (const v of vertices) { io.outVec(v); }
        }
    }, { exact: true });

    family.case('MeshSmoother.invalidInput', (io) => {
        const numVertices = io.integer();
        const numIndices = io.integer();
        const vertices: Vector[] = [];
        for (let i = 0; i < numVertices; ++i) { vertices.push(new Vector(3)); }
        const smoother = new MeshSmoother();
        smoother.initialize(vertices, new Array<number>(numIndices).fill(0));
        io.outInt(smoother.getNumTriangles());
    }, { exact: true });

    // ------------------------------------------------------------ MeshCurvature

    family.case('MeshCurvature.compute', (io) => {
        const m = readSurfaceMesh(io);
        const threshold = io.real();
        const curvature = new MeshCurvature();
        curvature.compute(m.vertices, m.indices, threshold);
        const kmin = curvature.getMinCurvatures(), kmax = curvature.getMaxCurvatures();
        // Independent check on the sphere samples of radius R (kind 4) with
        // the default threshold: the estimates are near 1/R.
        // Only the unperturbed samples (vertex 0 exactly on the x-axis): the
        // one-level subdivision is too coarse for tangentially jittered
        // samples (their estimates range over [0.3/R, 1.9/R]).
        const unperturbed = m.vertices[0].get(1) === 0 && m.vertices[0].get(2) === 0;
        if (m.kind === 4 && threshold === 0 && unperturbed) {
            const r = Math.hypot(m.vertices[0].get(0), m.vertices[0].get(1), m.vertices[0].get(2));
            for (let i = 0; i < kmin.length; ++i) {
                if (!(kmin[i] > 0.5 / r && kmax[i] < 1.5 / r && kmin[i] <= kmax[i])) {
                    throw new Error(`curvature ${kmin[i]}, ${kmax[i]} far from 1/R = ${1 / r}`);
                }
            }
        }
        for (const v of curvature.getNormals()) { io.outVec(v); }
        io.outReals(kmin);
        io.outReals(kmax);
        for (const v of curvature.getMinDirections()) { io.outVec(v); }
        for (const v of curvature.getMaxDirections()) { io.outVec(v); }
    }, { exact: true });

    // ------------------------------------------------------------ ImplicitSurface3

    family.case('ImplicitSurface3.queries', (io) => {
        const mode = io.integer();
        const c = io.reals(11);
        const p = io.vec(3);
        const epsilon = io.real();
        const surface = new CubicSurface(c);
        io.outBool(surface.isOnSurface(p, epsilon));
        io.outReal(surface.f(p));
        io.outVec(surface.getGradient(p));
        io.outMat(surface.getHessian(p));
        const frame = surface.getFrame(p);
        io.outVec(frame.tangent0);
        io.outVec(frame.tangent1);
        io.outVec(frame.normal);
        const info = surface.getPrincipalInformation(p);
        // Independent check: both principal curvatures of the sphere of
        // radius R are 1/R (A = Hessian / |gradient| = I / R).
        if (mode === 1) {
            const r = Math.sqrt(-c[9]);
            for (const k of [info.curvature0, info.curvature1]) {
                if (Math.abs(k * r - 1) > 1e-12) {
                    throw new Error(`sphere curvature ${k} is not 1/R = ${1 / r}`);
                }
            }
        }
        io.outBool(info.valid);
        io.outReal(info.curvature0);
        io.outReal(info.curvature1);
        io.outVec(info.direction0);
        io.outVec(info.direction1);
    }, { exact: true });

    family.finish();
});
