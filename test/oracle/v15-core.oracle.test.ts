// Replays oracle/cpp/cases/v15-core.cpp. Keep the two files in the same order.
//
// Bit patterns travel as two 32-bit integers, high half first (u64/outU64).
// Where an independent reference is cheap (BigInt bit arithmetic, an
// ordered-integer model of IEEE neighbours, the heap order, mesh partition
// invariants) the replay checks the port against it too and throws on a
// violation; the harness reports that as "the port threw", so agreement with
// C++ alone never passes a wrong result.
import { describe } from 'vitest';
import { Array2 } from '../../src/Array2.js';
import { Array3 } from '../../src/Array3.js';
import { Array4 } from '../../src/Array4.js';
import { BitHacks } from '../../src/BitHacks.js';
import { BSNumber } from '../../src/BSNumber.js';
import { BSRational } from '../../src/BSRational.js';
import * as C from '../../src/Constants.js';
import {
    CurveExtractor, CurveExtractorEdge, CurveExtractorVertex
} from '../../src/CurveExtractor.js';
import * as F from '../../src/Functions.js';
import { hashCombine, hashValue } from '../../src/HashCombine.js';
import { IEEEBinary32, IEEEBinary64 } from '../../src/IEEEBinary.js';
import { LexicoArray2 } from '../../src/LexicoArray2.js';
import { logAssert, logError } from '../../src/Logger.js';
import { MeshStaticManifold2 } from '../../src/MeshStaticManifold2.js';
import { MeshStaticManifold3 } from '../../src/MeshStaticManifold3.js';
import { MinHeap, type MinHeapRecord } from '../../src/MinHeap.js';
import { QFNumber } from '../../src/QFNumber.js';
import { hasDivisionOperator, isArbitraryPrecision } from '../../src/TypeTraits.js';
import { UniqueVerticesTriangles } from '../../src/UniqueVerticesTriangles.js';
import { OracleFamily, type OracleIO } from './harness.js';

function u64(io: OracleIO): bigint {
    const hi = BigInt(io.integer());
    const lo = BigInt(io.integer());
    return (hi << 32n) | lo;
}

function outU64(io: OracleIO, b: bigint): void {
    io.outInt(Number(b >> 32n));
    io.outInt(Number(b & 0xFFFFFFFFn));
}

// An independent-reference violation. Thrown from a case body it is reported
// as "the port threw ... but C++ returned normally".
function check(condition: boolean, what: string): void {
    if (!condition) { throw new Error(`reference check failed: ${what}`); }
}

// ---------------------------------------------------------------- BitHacks

function refLeading(v: bigint): number {
    return v === 0n ? 0 : v.toString(2).length - 1;
}

function refTrailing(v: bigint): number {
    if (v === 0n) { return 0; }
    let k = 0;
    while (((v >> BigInt(k)) & 1n) === 0n) { ++k; }
    return k;
}

function emitBits32(io: OracleIO, u: number): void {
    const i = u | 0;
    const b = BigInt(u);
    const pow2 = BitHacks.isPowerOfTwo(u);
    check(pow2 === (b !== 0n && (b & (b - 1n)) === 0n), `isPowerOfTwo(${u})`);
    io.outBool(pow2);
    io.outBool(BitHacks.isPowerOfTwo(i));
    const log2 = BitHacks.log2OfPowerOfTwo(u);
    check(!pow2 || log2 === refLeading(b), `log2OfPowerOfTwo(${u})`);
    io.outInt(log2);
    io.outInt(BitHacks.log2OfPowerOfTwo(i));
    const leading = BitHacks.getLeadingBit(u);
    check(leading === refLeading(b), `getLeadingBit(${u})`);
    io.outInt(leading);
    io.outInt(BitHacks.getLeadingBit(i));
    const trailing = BitHacks.getTrailingBit(u);
    check(trailing === refTrailing(b), `getTrailingBit(${u})`);
    io.outInt(trailing);
    io.outInt(BitHacks.getTrailingBit(i));
    const up = BitHacks.roundUpToPowerOfTwo(u);
    const refUp = b === 0n ? 1n
        : ((b & (b - 1n)) === 0n ? b : 1n << BigInt(refLeading(b) + 1));
    check(BigInt(up) === refUp, `roundUpToPowerOfTwo(${u})`);
    io.outInt(up);
    const down = BitHacks.roundDownToPowerOfTwo(u);
    check(down === (u === 0 ? 0 : 2 ** refLeading(b)), `roundDownToPowerOfTwo(${u})`);
    io.outInt(down);
}

function emitBits64(io: OracleIO, u: bigint): void {
    const i = BigInt.asIntN(64, u);
    const leading = BitHacks.getLeadingBit64(u);
    check(leading === refLeading(u), `getLeadingBit64(${u})`);
    io.outInt(leading);
    io.outInt(BitHacks.getLeadingBit64(i));
    const trailing = BitHacks.getTrailingBit64(u);
    check(trailing === refTrailing(u), `getTrailingBit64(${u})`);
    io.outInt(trailing);
    io.outInt(BitHacks.getTrailingBit64(i));
}

// ---------------------------------------------------------------- IEEEBinary

// Independent model of the neighbour functions: map an encoding to an
// integer that is monotone in the represented value (-0 and +0 both to 0),
// step by one and map back (-0 when stepping up from -MIN_SUBNORMAL, +0 when
// stepping down from +MIN_SUBNORMAL, infinities and NaNs as upstream
// documents).
function refNext(bits: bigint, width: 32 | 64, up: boolean): bigint {
    const t = width === 64 ? 52n : 23n;
    const signMask = 1n << BigInt(width - 1);
    const magMask = signMask - 1n;
    const inf = ((1n << (BigInt(width) - 1n - t)) - 1n) << t;
    const mag = bits & magMask;
    const neg = (bits & signMask) !== 0n;
    if (mag > inf) { return bits; }
    if (mag === inf) {
        if (up) { return neg ? signMask | (inf - 1n) : inf; }
        return neg ? signMask | inf : inf - 1n;
    }
    const o = (neg ? -mag : mag) + (up ? 1n : -1n);
    if (mag === 0n) { return up ? 1n : signMask | 1n; }
    if (o === 0n) { return up ? signMask : 0n; }
    return o > 0n ? o : signMask | -o;
}

function refClass(bits: bigint, width: 32 | 64): number {
    // IEEEClassification order: NEG_INFINITY, NEG_SUBNORMAL, NEG_NORMAL,
    // NEG_ZERO, POS_ZERO, POS_SUBNORMAL, POS_NORMAL, POS_INFINITY,
    // QUIET_NAN, SIGNALING_NAN.
    const t = width === 64 ? 52n : 23n;
    const signMask = 1n << BigInt(width - 1);
    const mag = bits & (signMask - 1n);
    const inf = ((1n << (BigInt(width) - 1n - t)) - 1n) << t;
    const neg = (bits & signMask) !== 0n;
    if (mag > inf) { return (mag & (1n << (t - 1n))) !== 0n ? 8 : 9; }
    if (mag === inf) { return neg ? 0 : 7; }
    if (mag === 0n) { return neg ? 3 : 4; }
    if (mag < (1n << t)) { return neg ? 1 : 5; }
    return neg ? 2 : 6;
}

function emitFields64(io: OracleIO, x: IEEEBinary64): void {
    const bits = x.encoding;
    outU64(io, x.encoding);
    io.outInt(x.getSign());
    io.outInt(x.getBiased());
    outU64(io, x.getTrailing());
    const { sign, biased, trailing } = x.getEncoding();
    io.outInt(sign);
    io.outInt(biased);
    outU64(io, trailing);
    const c = x.getClassification();
    check(c === refClass(bits, 64), `classification of ${bits.toString(16)}`);
    io.outInt(c);
    io.outBool(x.isZero());
    io.outBool(x.isSignMinus());
    io.outBool(x.isSubnormal());
    io.outBool(x.isNormal());
    io.outBool(x.isFinite());
    io.outBool(x.isInfinite());
    io.outBool(x.isNaN());
    io.outBool(x.isQuietNaN());
    io.outBool(x.isSignalingNaN());
    const nu = x.getNextUp();
    const nd = x.getNextDown();
    check(nu === refNext(bits, 64, true) && nd === refNext(bits, 64, false),
        `neighbours of ${bits.toString(16)}`);
    outU64(io, nu);
    outU64(io, nd);
    outU64(io, IEEEBinary64.fromParts(sign, biased, trailing).encoding);
}

function emitFields32(io: OracleIO, x: IEEEBinary32): void {
    const bits = BigInt(x.encoding);
    io.outInt(x.encoding);
    io.outInt(x.getSign());
    io.outInt(x.getBiased());
    io.outInt(x.getTrailing());
    const { sign, biased, trailing } = x.getEncoding();
    io.outInt(sign);
    io.outInt(biased);
    io.outInt(trailing);
    const c = x.getClassification();
    check(c === refClass(bits, 32), `classification of ${bits.toString(16)}`);
    io.outInt(c);
    io.outBool(x.isZero());
    io.outBool(x.isSignMinus());
    io.outBool(x.isSubnormal());
    io.outBool(x.isNormal());
    io.outBool(x.isFinite());
    io.outBool(x.isInfinite());
    io.outBool(x.isNaN());
    io.outBool(x.isQuietNaN());
    io.outBool(x.isSignalingNaN());
    const nu = x.getNextUp();
    const nd = x.getNextDown();
    check(BigInt(nu) === refNext(bits, 32, true) && BigInt(nd) === refNext(bits, 32, false),
        `neighbours of ${bits.toString(16)}`);
    io.outInt(nu);
    io.outInt(nd);
    io.outInt(IEEEBinary32.fromParts(sign, biased, trailing).encoding);
}

function emitConstants32(io: OracleIO): void {
    const B = IEEEBinary32;
    for (const v of [B.NUM_ENCODING_BITS, B.NUM_EXPONENT_BITS, B.NUM_SIGNIFICAND_BITS,
        B.NUM_TRAILING_BITS, B.EXPONENT_BIAS, B.MAX_BIASED_EXPONENT, B.MIN_SUB_EXPONENT,
        B.MIN_EXPONENT, B.SIGN_SHIFT, B.SIGN_MASK, B.NOT_SIGN_MASK, B.TRAILING_MASK,
        B.EXPONENT_MASK, B.NAN_QUIET_MASK, B.NAN_PAYLOAD_MASK, B.MAX_TRAILING, B.SUP_TRAILING,
        B.POS_ZERO, B.NEG_ZERO, B.MIN_SUBNORMAL, B.MAX_SUBNORMAL, B.MIN_NORMAL, B.MAX_NORMAL,
        B.POS_INFINITY, B.NEG_INFINITY]) {
        io.outInt(v);
    }
}

function emitConstants64(io: OracleIO): void {
    const B = IEEEBinary64;
    for (const v of [B.NUM_ENCODING_BITS, B.NUM_EXPONENT_BITS, B.NUM_SIGNIFICAND_BITS,
        B.NUM_TRAILING_BITS, B.EXPONENT_BIAS, B.MAX_BIASED_EXPONENT, B.MIN_SUB_EXPONENT,
        B.MIN_EXPONENT, B.SIGN_SHIFT]) {
        io.outInt(v);
    }
    for (const v of [B.SIGN_MASK, B.NOT_SIGN_MASK, B.TRAILING_MASK, B.EXPONENT_MASK,
        B.NAN_QUIET_MASK, B.NAN_PAYLOAD_MASK, B.MAX_TRAILING, B.SUP_TRAILING, B.POS_ZERO,
        B.NEG_ZERO, B.MIN_SUBNORMAL, B.MAX_SUBNORMAL, B.MIN_NORMAL, B.MAX_NORMAL,
        B.POS_INFINITY, B.NEG_INFINITY]) {
        outU64(io, v);
    }
}

// ---------------------------------------------------------------- HashCombine

// The port's documented stand-in for std::hash<double>, recomputed with
// BigInt, and the 32-bit combine step.
function refHashCombine(seed: number, x: number): number {
    const view = new DataView(new ArrayBuffer(8));
    view.setFloat64(0, x === 0 ? 0 : x);
    const b = view.getBigUint64(0);
    const h = (b & 0xFFFFFFFFn) ^ (b >> 32n);
    const s = BigInt(seed);
    const sum = h + 0x9e3779b9n + ((s << 6n) & 0xFFFFFFFFn) + (s >> 2n);
    return Number((s ^ sum) & 0xFFFFFFFFn);
}

// ---------------------------------------------------------------- MinHeap

type HeapRecord = MinHeapRecord<number, number>;

interface HeapInternals {
    mRecords: HeapRecord[];
    mPointers: HeapRecord[];
}

function runHeapSequence(io: OracleIO, allowNaN: boolean): void {
    const maxElements = io.integer();
    const heap = new MinHeap<number, number>(maxElements);
    const internals = heap as unknown as HeapInternals;
    const handles: (HeapRecord | null)[] = [];
    // Independent model of the live multiset of values (for NaN-free
    // sequences): a removed or reported minimum is <= every live value.
    const live = new Map<number, number>();
    // Updating a removed record (a handle Insert returned earlier) sifts it
    // back among the live slots and pushes a live record out, in upstream
    // as in the port; the model stops there.
    let modelValid = !allowNaN;
    const numOps = io.integer();
    for (let op = 0; op < numOps; ++op) {
        const opcode = io.integer();
        if (opcode <= 2) {
            const value = io.real();
            const key = handles.length;
            const record = heap.insert(key, value);
            handles.push(record);
            if (record !== null) { live.set(key, value); }
            io.outInt(record !== null ? record.key : -1);
        } else if (opcode === 3 || opcode === 6) {
            const r = opcode === 3 ? heap.remove() : heap.getMinimum();
            io.outBool(r !== null);
            if (r !== null) {
                if (modelValid) {
                    for (const v of live.values()) {
                        check(r.value <= v, `minimum ${r.value} > live ${v}`);
                    }
                }
                if (opcode === 3) { live.delete(r.key); }
                io.outInt(r.key);
                io.outReal(r.value);
            }
        } else if (handles.length > 0) {
            const which = io.integer();
            const value = io.real();
            const record = handles[which];
            heap.update(record, value);
            if (record !== null && !live.has(record.key)) { modelValid = false; }
            if (record !== null && live.has(record.key)) { live.set(record.key, record.value); }
        }
        check(!modelValid || heap.getNumElements() === live.size, 'live count');
        io.outInt(heap.getNumElements());
        io.outBool(heap.isValid());
    }
    for (let i = 0; i < maxElements; ++i) {
        const record = internals.mPointers[i];
        io.outInt(internals.mRecords.indexOf(record));
        io.outInt(record.index);
        if (i < heap.getNumElements()) {
            io.outInt(record.key);
            io.outReal(record.value);
        }
    }
}

// ---------------------------------------------------------------- meshes

const INVALID2 = MeshStaticManifold2.invalid;
const INVALID3 = MeshStaticManifold3.invalid;

function outIndex(io: OracleIO, i: number, invalid: number): void {
    io.outInt(i === invalid ? -1 : i);
}

interface Mesh {
    numVertices: number;
    simplices: number[][];
}

function readMesh(io: OracleIO, n: number): Mesh {
    const numVertices = io.integer();
    const count = io.integer();
    const simplices: number[][] = [];
    for (let s = 0; s < count; ++s) {
        const simplex: number[] = [];
        for (let k = 0; k < n; ++k) { simplex.push(io.integer()); }
        simplices.push(simplex);
    }
    return { numVertices, simplices };
}

function directedEdgeTriangle(mesh: Mesh, v0: number, v1: number): number {
    for (let t = 0; t < mesh.simplices.length; ++t) {
        const tri = mesh.simplices[t];
        for (let k = 0; k < 3; ++k) {
            if (tri[k] === v0 && tri[(k + 1) % 3] === v1) { return t; }
        }
    }
    return -1;
}

// C++ passes -1 as size_t, which is SIZE_MAX (upstream's 'invalid'). The
// port's guards test 0 <= v explicitly, so -1 reaches the same branch.
function asIndex(i: number): number {
    return i;
}

// ---------------------------------------------------------------- UniqueVerticesTriangles

function readVertexList(io: OracleIO): number[][] {
    const count = io.integer();
    const list: number[][] = [];
    for (let i = 0; i < count; ++i) { list.push(io.reals(3)); }
    return list;
}

function readIndices(io: OracleIO, numVertices: number): number[] {
    const numTriangles = io.integer();
    const indices: number[] = [];
    for (let i = 0; i < 3 * numTriangles; ++i) { indices.push(io.integer()); }
    check(indices.every((i) => i >= 0 && i < numVertices), 'index range');
    return indices;
}

function toTriples(indices: number[]): number[][] {
    const triples: number[][] = [];
    for (let t = 0; t < indices.length / 3; ++t) {
        triples.push([indices[3 * t], indices[3 * t + 1], indices[3 * t + 2]]);
    }
    return triples;
}

// Postconditions of every operation: indices address the output pool, and
// each index maps an input vertex to an equal (component-wise ==) one.
function emitUnique(io: OracleIO, vertices: (number | number[])[], indices: number[]): void {
    check(indices.every((i) => i >= 0 && i < vertices.length), 'output index range');
    io.outInt(vertices.length);
    for (const v of vertices) {
        if (typeof v === 'number') { io.outReal(v); } else { io.outReals(v); }
    }
    io.outInt(indices.length);
    for (const i of indices) { io.outInt(i); }
}

function sameVertex(a: number | number[], b: number | number[]): boolean {
    if (typeof a === 'number' || typeof b === 'number') { return a === b; }
    return a.length === b.length && a.every((x, k) => x === b[k]);
}

// ---------------------------------------------------------------- CurveExtractor

// The replay subclass of the case file's ReplayExtractor: extract() returns
// the recorded rational vertices and edges.
class ReplayExtractor extends CurveExtractor {
    readonly vertices: CurveExtractorVertex[] = [];
    readonly edges: CurveExtractorEdge[] = [];

    constructor(xBound = 2, yBound = 2) {
        super(xBound, yBound, new Array<number>(16).fill(0));
    }

    addRecordedEdge(e: number[]): void {
        this.addEdge(this.vertices, this.edges, e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]);
    }

    addRecordedVertex(v: number[]): void {
        this.addVertex(this.vertices, v[0], v[1], v[2], v[3]);
    }

    extract(_level: number): { vertices: CurveExtractorVertex[]; edges: CurveExtractorEdge[] } {
        return {
            vertices: this.vertices.map((v) => {
                const w = new CurveExtractorVertex();
                w.xNumer = v.xNumer; w.xDenom = v.xDenom; w.yNumer = v.yNumer; w.yDenom = v.yDenom;
                return w;
            }),
            edges: this.edges.map((e) => new CurveExtractorEdge(e.v[0], e.v[1]))
        };
    }
}

describe('oracle: v15-core', () => {
    const family = new OracleFamily('v15-core');

    // ------------------------------------------------------------ Constants

    family.case('Constants.values', (io) => {
        io.outReals([C.GTE_C_PI, C.GTE_C_HALF_PI, C.GTE_C_QUARTER_PI, C.GTE_C_TWO_PI,
            C.GTE_C_INV_PI, C.GTE_C_INV_TWO_PI, C.GTE_C_INV_HALF_PI, C.GTE_C_DEG_TO_RAD,
            C.GTE_C_RAD_TO_DEG, C.GTE_C_SQRT_2, C.GTE_C_INV_SQRT_2, C.GTE_C_LN_2,
            C.GTE_C_INV_LN_2, C.GTE_C_LN_10, C.GTE_C_INV_LN_10]);
    }, { exact: true });

    // ------------------------------------------------------------ Functions

    family.case('Functions.arithmetic', (io) => {
        const x = io.real();
        io.outReal(F.sign(x));
        io.outInt(F.isign(x));
        io.outReal(F.saturate(x));
        io.outReal(F.sqr(x));
        io.outReal(F.invsqrt(x));
    }, { exact: true });

    family.case('Functions.clamp', (io) => {
        const x = io.real();
        const xmin = io.real();
        const xmax = io.real();
        io.outReal(F.clamp(x, xmin, xmax));
    }, { exact: true });

    // Signed-zero ties of clamp; the C++ side compiles Functions.h without
    // optimization because MSVC /O2 turns the inner conditional into a
    // minsd that returns x instead of xmax on such ties (see the case file).
    family.case('Functions.clamp.signedZero', (io) => {
        const x = io.real();
        const xmin = io.real();
        const xmax = io.real();
        io.outReal(F.clamp(x, xmin, xmax));
    }, { exact: true });

    family.case('Functions.stdMaxMin', (io) => {
        const a = io.real();
        const b = io.real();
        io.outReal(F.stdMax(a, b));
        io.outReal(F.stdMin(a, b));
    }, { exact: true });

    // std::atan, std::atan2, std::cos, std::sin, std::exp (MSVC CRT) against
    // Math.* (V8): the default 1e-12 tolerance.
    family.case('Functions.libm', (io) => {
        const x = io.real();
        const y = io.real();
        io.outReal(F.atandivpi(x));
        io.outReal(F.atan2divpi(y, x));
        io.outReal(F.cospi(x));
        io.outReal(F.sinpi(x));
        io.outReal(F.exp10(x));
    });

    for (const name of ['Functions.fma', 'Functions.fma.ties', 'Functions.fma.special']) {
        family.case(name, (io) => {
            const u = io.real();
            const v = io.real();
            const w = io.real();
            io.outReal(F.fma(u, v, w));
        }, { exact: true });
    }

    family.case('Functions.robustSOP', (io) => {
        const [u, v, w, z] = io.reals(4);
        io.outReal(F.robustSOP(u, v, w, z));
    }, { exact: true });

    family.case('Functions.robustDOP', (io) => {
        const [u, v, w, z] = io.reals(4);
        io.outReal(F.robustDOP(u, v, w, z));
    }, { exact: true });

    // ------------------------------------------------------------ BitHacks

    family.case('BitHacks.bits32', (io) => {
        emitBits32(io, io.integer());
    }, { exact: true });

    family.case('BitHacks.bits64', (io) => {
        emitBits64(io, u64(io));
    }, { exact: true });

    family.case('BitHacks.powerSweep', (io) => {
        const part = io.integer();
        for (let k = 4 * part; k < 4 * part + 4; ++k) {
            const p = 1n << BigInt(k);
            for (const v of [p - 1n, p, p + 1n]) {
                emitBits64(io, v);
                emitBits32(io, Number(v & 0xFFFFFFFFn));
            }
        }
    }, { exact: true });

    // ------------------------------------------------------------ IEEEBinary

    family.case('IEEEBinary64.fields', (io) => {
        const x = IEEEBinary64.fromEncoding(u64(io));
        emitFields64(io, x);
        if (!x.isNaN()) { io.outReal(x.number); }
    }, { exact: true });

    family.case('IEEEBinary32.fields', (io) => {
        const x = IEEEBinary32.fromEncoding(io.integer());
        emitFields32(io, x);
        if (!x.isNaN()) { io.outReal(x.number); }
    }, { exact: true });

    family.case('IEEEBinary64.fromNumber', (io) => {
        const x = IEEEBinary64.fromNumber(io.real());
        outU64(io, x.encoding);
        io.outInt(x.getClassification());
    }, { exact: true });

    family.case('IEEEBinary32.fromNumber', (io) => {
        const x = IEEEBinary32.fromNumber(io.real());
        io.outInt(x.encoding);
        io.outInt(x.getClassification());
    }, { exact: true });

    family.case('IEEEBinary64.setEncoding', (io) => {
        const sign = io.integer();
        const biased = io.integer();
        const trailing = u64(io);
        outU64(io, IEEEBinary64.fromParts(sign, biased, trailing).encoding);
    }, { exact: true });

    family.case('IEEEBinary32.setEncoding', (io) => {
        const sign = io.integer();
        const biased = io.integer();
        const trailing = io.integer();
        io.outInt(IEEEBinary32.fromParts(sign, biased, trailing).encoding);
    }, { exact: true });

    family.case('IEEEBinary.constants', (io) => {
        emitConstants32(io);
        emitConstants64(io);
    }, { exact: true });

    // ------------------------------------------------------------ HashCombine

    family.case('HashCombine.combine', (io) => {
        const seed = io.integer();
        const x = io.real();
        const h = hashCombine(seed, x);
        check(h === refHashCombine(seed, x), 'hashCombine reference');
        io.outInt(h);
    }, { exact: true });

    // The inputs include signaling NaNs, so they are read one by one: V8
    // quiets a signaling NaN stored into a JS double array (io.reals builds
    // one), which turned fff7ffffffffffff into ffffffffffffffff, whose
    // 32-bit fold collides with the hash of 0 (deep run, records 1188 and
    // 1546 of the first version of this replay).
    family.case('HashCombine.equalities', (io) => {
        const a = io.real();
        const b = io.real();
        const c = io.real();
        const d = io.real();
        io.outBool(hashValue(a, b) === hashValue(c, d));
        io.outBool(hashValue(a) === hashValue(c));
    }, { exact: true });

    family.case('HashCombine.equalities.signalingNaN', (io) => {
        const a = io.real();
        const c = io.real();
        io.outBool(hashValue(a) === hashValue(c));
        io.outBool(hashValue(a, c) === hashValue(c, a));
    }, { exact: true });

    family.case('HashCombine.hashValue.msvc', (io) => {
        const [a, b, c] = io.reals(3);
        outU64(io, BigInt(hashValue(a, b, c)));
    }, { exact: true, deviation: 'porting-status.json HashCombine.h: 32-bit seeds and the '
        + 'port\'s own std::hash<double> stand-in; no C++ library\'s hash values are reproduced' });

    // ------------------------------------------------------------ TypeTraits

    family.case('TypeTraits.traits', (io) => {
        const values: unknown[] = [1.5, BSNumber.fromNumber(1.5), BSRational.fromNumber(1.5),
            new QFNumber(1, 2, 3), new QFNumber(1, 2, 3)];
        for (const v of values) {
            io.outBool(isArbitraryPrecision(v));
            io.outBool(hasDivisionOperator(v));
        }
    }, { exact: true });

    // ------------------------------------------------------------ Logger

    family.case('Logger.logAssert', (io) => {
        const condition = io.boolean();
        logAssert(condition, 'oracle assertion');
        io.outBool(condition);
    }, { exact: true });

    family.case('Logger.logError', (io) => {
        const unused = io.integer();
        if (unused >= 0) { logError('oracle error'); }
        io.outInt(unused);
    }, { exact: true });

    // ------------------------------------------------------------ MinHeap

    family.case('MinHeap.sequence', (io) => {
        runHeapSequence(io, false);
    }, { exact: true });

    family.case('MinHeap.sequence.nan', (io) => {
        runHeapSequence(io, true);
    }, { exact: true });

    // ------------------------------------------------------------ Array2/3/4

    function readObjects(io: OracleIO, n: number): number[] {
        return io.reals(n);
    }

    family.case('Array2.access', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const objects = readObjects(io, b0 * b1);
        const a = new Array2<number>(b0, b1, objects);
        io.outInt(a.getBound0());
        io.outInt(a.getBound1());
        const numOps = io.integer();
        for (let op = 0; op < numOps; ++op) {
            const i0 = io.integer();
            const i1 = io.integer();
            if (io.boolean()) { a.set(i0, i1, io.real()); } else { io.outReal(a.get(i0, i1)); }
        }
        io.outReals(objects);
    }, { exact: true });

    family.case('Array2.owned', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const a = new Array2<number>(b0, b1);
        for (let i1 = 0; i1 < b1; ++i1) {
            for (let i0 = 0; i0 < b0; ++i0) { a.set(i0, i1, io.real()); }
        }
        const numReads = io.integer();
        for (let r = 0; r < numReads; ++r) {
            const i0 = io.integer();
            const i1 = io.integer();
            io.outReal(a.get(i0, i1));
        }
    }, { exact: true });

    // The port's generic storage is not zero-filled (src/Array2.ts, NOTE in
    // the constructor comment): a never-written element reads undefined.
    family.case('Array2.owned.zeroInit', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const a = new Array2<number>(b0, b1);
        const i0 = io.integer();
        const i1 = io.integer();
        io.outReal(a.get(i0, i1));
    }, { exact: true, deviation: 'src/Array2.ts constructor NOTE: the owned storage of a '
        + 'generic T is not value-initialized; upstream std::vector<double> reads 0' });

    family.case('Array3.access', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const b2 = io.integer();
        const objects = readObjects(io, b0 * b1 * b2);
        const a = new Array3<number>(b0, b1, b2, objects);
        io.outInt(a.getBound0());
        io.outInt(a.getBound1());
        io.outInt(a.getBound2());
        const numOps = io.integer();
        for (let op = 0; op < numOps; ++op) {
            const i0 = io.integer();
            const i1 = io.integer();
            const i2 = io.integer();
            if (io.boolean()) { a.set(i0, i1, i2, io.real()); } else { io.outReal(a.get(i0, i1, i2)); }
        }
        io.outReals(objects);
    }, { exact: true });

    family.case('Array4.access', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const b2 = io.integer();
        const b3 = io.integer();
        const objects = readObjects(io, b0 * b1 * b2 * b3);
        const a = new Array4<number>(b0, b1, b2, b3, objects);
        io.outInt(a.getBound0());
        io.outInt(a.getBound1());
        io.outInt(a.getBound2());
        io.outInt(a.getBound3());
        const numOps = io.integer();
        for (let op = 0; op < numOps; ++op) {
            const i0 = io.integer();
            const i1 = io.integer();
            const i2 = io.integer();
            const i3 = io.integer();
            if (io.boolean()) {
                a.set(i0, i1, i2, i3, io.real());
            } else {
                io.outReal(a.get(i0, i1, i2, i3));
            }
        }
        io.outReals(objects);
    }, { exact: true });

    family.case('Array3.owned', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const b2 = io.integer();
        const a = new Array3<number>(b0, b1, b2);
        for (let i2 = 0; i2 < b2; ++i2) {
            for (let i1 = 0; i1 < b1; ++i1) {
                for (let i0 = 0; i0 < b0; ++i0) { a.set(i0, i1, i2, io.real()); }
            }
        }
        const i0 = io.integer();
        const i1 = io.integer();
        const i2 = io.integer();
        io.outReal(a.get(i0, i1, i2));
    }, { exact: true });

    family.case('Array4.owned', (io) => {
        const b0 = io.integer();
        const b1 = io.integer();
        const b2 = io.integer();
        const b3 = io.integer();
        const a = new Array4<number>(b0, b1, b2, b3);
        for (let i3 = 0; i3 < b3; ++i3) {
            for (let i2 = 0; i2 < b2; ++i2) {
                for (let i1 = 0; i1 < b1; ++i1) {
                    for (let i0 = 0; i0 < b0; ++i0) { a.set(i0, i1, i2, i3, io.real()); }
                }
            }
        }
        const i0 = io.integer();
        const i1 = io.integer();
        const i2 = io.integer();
        const i3 = io.integer();
        io.outReal(a.get(i0, i1, i2, i3));
    }, { exact: true });

    // ------------------------------------------------------------ LexicoArray2

    family.case('LexicoArray2.access', (io) => {
        const mode = io.integer();
        let numRows = mode <= 3 ? 3 : 2;
        let numCols = mode <= 3 ? 4 : 5;
        if (mode <= 1) {
            numRows = io.integer();
            numCols = io.integer();
        }
        const storage = readObjects(io, numRows * numCols);
        const rowMajor = mode % 2 === 0;
        const a = new LexicoArray2(rowMajor, numRows, numCols, storage);
        io.outInt(a.getNumRows());
        io.outInt(a.getNumCols());
        const numOps = io.integer();
        for (let op = 0; op < numOps; ++op) {
            const r = io.integer();
            const c = io.integer();
            if (io.boolean()) { a.set(r, c, io.real()); } else { io.outReal(a.get(r, c)); }
        }
        io.outReals(storage);
    }, { exact: true });

    // ------------------------------------------------------------ MeshStaticManifold2

    family.case('MeshStaticManifold2.construct', (io) => {
        const mesh = readMesh(io, 3);
        const numThreads = io.integer();
        const m = new MeshStaticManifold2(mesh.numVertices, mesh.simplices, numThreads);
        io.outInt(m.getMinNumTrianglesAtVertex());
        io.outInt(m.getMaxNumTrianglesAtVertex());
        for (const tri of m.getTriangles()) { for (const v of tri) { io.outInt(v); } }
        for (const adj of m.getAdjacents()) { for (const a of adj) { outIndex(io, a, INVALID2); } }
        for (const vertex of m.getVertices()) {
            io.outInt(vertex.getNumAdjacents());
            for (const tuple of vertex.getAdjacents()) {
                for (const x of tuple) { outIndex(io, x, INVALID2); }
            }
        }
    }, { exact: true });

    family.case('MeshStaticManifold2.queries', (io) => {
        const mesh = readMesh(io, 3);
        const m = new MeshStaticManifold2(mesh.numVertices, mesh.simplices, 0);
        const numQueries = io.integer();
        for (let q = 0; q < numQueries; ++q) {
            const v0 = asIndex(io.integer());
            const v1 = asIndex(io.integer());
            io.outBool(m.edgeExists(v0, v1));
            const r = m.getAdjacentTriangles(v0, v1);
            io.outBool(r.exists);
            const directed = v0 >= 0 && v1 >= 0 && v0 !== v1 && v0 < mesh.numVertices
                && v1 < mesh.numVertices && directedEdgeTriangle(mesh, v0, v1) >= 0;
            io.outBool(directed);
            if (directed) { outIndex(io, r.adj0, INVALID2); }
        }
    }, { exact: true });

    // Finding #66 (UPSTREAM-FINDINGS, MeshStaticManifold2.h item 1): upstream
    // returns the neighbour across the edge opposite v0 as adj1.
    family.case('MeshStaticManifold2.adjacentTriangles', (io) => {
        const mesh = readMesh(io, 3);
        const m = new MeshStaticManifold2(mesh.numVertices, mesh.simplices, 0);
        const v0 = io.integer();
        const v1 = io.integer();
        const r = m.getAdjacentTriangles(v0, v1);
        const l = directedEdgeTriangle(mesh, v0, v1);
        const rr = directedEdgeTriangle(mesh, v1, v0);
        check(r.adj0 === (l < 0 ? INVALID2 : l) && r.adj1 === (rr < 0 ? INVALID2 : rr),
            'documented adjacency contract');
        io.outBool(r.exists);
        outIndex(io, r.adj0, INVALID2);
        outIndex(io, r.adj1, INVALID2);
    }, { exact: true, deviation: '#66' });

    family.case('MeshStaticManifold2.components', (io) => {
        const mesh = readMesh(io, 3);
        const m = new MeshStaticManifold2(mesh.numVertices, mesh.simplices, 0);
        const components = m.getComponents();
        // Invariants: the components partition the triangles, and the
        // depth-first search from t reaches every valid mAdjacents[t][*], so
        // an adjacent triangle belongs to t's component or to an earlier one
        // (earlier only on the defective meshes, whose adjacency is not
        // symmetric: a duplicated or reversed triangle records a neighbour
        // that does not record it back).
        const owner = new Array<number>(mesh.simplices.length).fill(-1);
        components.forEach((c, i) => c.forEach((t) => {
            check(owner[t] === -1, 'triangle in two components');
            owner[t] = i;
        }));
        check(owner.every((o) => o >= 0), 'triangle in no component');
        m.getAdjacents().forEach((adj, t) => adj.forEach((a) => {
            check(a === INVALID2 || owner[a] <= owner[t], 'adjacency into a later component');
        }));
        io.outInt(components.length);
        for (const c of components) {
            io.outInt(c.length);
            for (const t of c) { io.outInt(t); }
        }
    }, { exact: true });

    family.case('MeshStaticManifold2.boundaryPolygons', (io) => {
        const mesh = readMesh(io, 3);
        const duplicateEndpoints = io.boolean();
        const m = new MeshStaticManifold2(mesh.numVertices, mesh.simplices, 0);
        const polygons = m.getBoundaryPolygons(duplicateEndpoints);
        // Invariant: every boundary edge (a directed edge with no reverse)
        // appears in exactly one polygon, as consecutive vertices.
        const boundary = new Set<string>();
        for (const tri of mesh.simplices) {
            for (let k = 0; k < 3; ++k) {
                const a = tri[k], b = tri[(k + 1) % 3];
                if (directedEdgeTriangle(mesh, b, a) < 0) { boundary.add(`${a},${b}`); }
            }
        }
        let numEdges = 0;
        for (const p of polygons) {
            const n = duplicateEndpoints ? p.length - 1 : p.length;
            for (let i = 0; i < n; ++i) {
                const key = `${p[i]},${p[(i + 1) % p.length]}`;
                check(boundary.has(key),
                    `polygon edge ${key} is not a boundary edge`);
                ++numEdges;
            }
        }
        check(numEdges === boundary.size, 'boundary edge count');
        io.outInt(polygons.length);
        for (const p of polygons) {
            io.outInt(p.length);
            for (const v of p) { io.outInt(v); }
        }
    }, { exact: true });

    family.case('MeshStaticManifold2.invalidInput', (io) => {
        const numVertices = io.integer();
        const empty = io.boolean();
        const triangles = empty ? [] : [[0, 1, 2]];
        const m = new MeshStaticManifold2(numVertices, triangles, 0);
        io.outInt(m.getMaxNumTrianglesAtVertex());
    }, { exact: true });

    // ------------------------------------------------------------ MeshStaticManifold3

    family.case('MeshStaticManifold3.construct', (io) => {
        const mesh = readMesh(io, 4);
        const numThreads = io.integer();
        const m = new MeshStaticManifold3(mesh.numVertices, mesh.simplices, numThreads);
        io.outInt(m.getMinNumTetrahedraAtVertex());
        io.outInt(m.getMaxNumTetrahedraAtVertex());
        for (const tet of m.getTetrahedra()) { for (const v of tet) { io.outInt(v); } }
        for (const adj of m.getAdjacents()) { for (const a of adj) { outIndex(io, a, INVALID3); } }
        for (const vertex of m.getVertices()) {
            io.outInt(vertex.getNumAdjacents());
            for (const tuple of vertex.getAdjacents()) {
                for (const x of tuple) { outIndex(io, x, INVALID3); }
            }
        }
    }, { exact: true });

    family.case('MeshStaticManifold3.queries', (io) => {
        const mesh = readMesh(io, 4);
        const m = new MeshStaticManifold3(mesh.numVertices, mesh.simplices, 0);
        const numQueries = io.integer();
        for (let q = 0; q < numQueries; ++q) {
            const v0 = asIndex(io.integer());
            const v1 = asIndex(io.integer());
            const v2 = asIndex(io.integer());
            io.outBool(m.faceExists(v0, v1, v2));
            io.outBool(m.getAdjacentTetrahedra(v0, v1, v2).exists);
        }
    }, { exact: true });

    // Finding #66 (MeshStaticManifold3.h item 2): upstream returns a vertex
    // index as adj0.
    family.case('MeshStaticManifold3.adjacentTetrahedra', (io) => {
        const mesh = readMesh(io, 4);
        const m = new MeshStaticManifold3(mesh.numVertices, mesh.simplices, 0);
        const tet = mesh.simplices[io.integer()];
        const f = MeshStaticManifold3.face[io.integer()];
        const r = m.getAdjacentTetrahedra(tet[f[0]], tet[f[1]], tet[f[2]]);
        io.outBool(r.exists);
        outIndex(io, r.adj0, INVALID3);
        outIndex(io, r.adj1, INVALID3);
    }, { exact: true, deviation: '#66' });

    family.case('MeshStaticManifold3.invalidInput', (io) => {
        const numVertices = io.integer();
        const empty = io.boolean();
        const tetrahedra = empty ? [] : [[0, 1, 2, 3]];
        const m = new MeshStaticManifold3(numVertices, tetrahedra, 0);
        io.outInt(m.getMaxNumTetrahedraAtVertex());
    }, { exact: true });

    // ------------------------------------------------------------ UniqueVerticesTriangles

    // Every output index maps its input vertex to an equal pool vertex.
    function checkRemap(inVertices: (number | number[])[], inIndices: number[] | null,
        outVertices: (number | number[])[], outIndices: number[]): void {
        outIndices.forEach((o, i) => {
            const source = inIndices === null ? inVertices[i] : inVertices[inIndices[i]];
            check(sameVertex(source, outVertices[o]), 'remapped vertex differs from its source');
        });
    }

    family.case('UniqueVerticesTriangles.generateIndexedTriangles', (io) => {
        const inVertices = readVertexList(io);
        const uvt = new UniqueVerticesTriangles<number[]>();
        const r = uvt.generateIndexedTriangles(inVertices);
        checkRemap(inVertices, null, r.vertices, r.indices);
        emitUnique(io, r.vertices, r.indices);
        const r3 = uvt.generateIndexedTrianglesTriples(inVertices);
        emitUnique(io, r3.vertices, r3.triangles.flat());
    }, { exact: true });

    for (const op of ['removeDuplicateVertices', 'removeUnusedVertices',
        'removeDuplicateAndUnusedVertices'] as const) {
        family.case(`UniqueVerticesTriangles.${op}`, (io) => {
            const inVertices = readVertexList(io);
            const inIndices = readIndices(io, inVertices.length);
            const uvt = new UniqueVerticesTriangles<number[]>();
            const r = uvt[op](inVertices, inIndices);
            checkRemap(inVertices, inIndices, r.vertices, r.indices);
            emitUnique(io, r.vertices, r.indices);
            const r3 = uvt[`${op}Triples`](inVertices, toTriples(inIndices));
            emitUnique(io, r3.vertices, r3.triangles.flat());
        }, { exact: true });
    }

    // The C++ side compiles the checks in with
    // GTL_VALIDATE_UNIQUE_VERTICES_TRIANGLES; the port's switch is 'validate'.
    family.case('UniqueVerticesTriangles.validate', (io) => {
        const op = io.integer();
        io.boolean(); // the generator mode (valid or unconstrained)
        const numVertices = io.integer();
        const inVertices: number[][] = [];
        for (let i = 0; i < numVertices; ++i) { inVertices.push([i, 0, 0]); }
        const numIndices = io.integer();
        const inIndices: number[] = [];
        for (let i = 0; i < numIndices; ++i) { inIndices.push(io.integer()); }
        const uvt = new UniqueVerticesTriangles<number[]>();
        uvt.validate = true;
        let numOut = 0;
        let numOutIndices = 0;
        if (op === 0) {
            const r = uvt.generateIndexedTriangles(inVertices);
            numOut = r.vertices.length; numOutIndices = r.indices.length;
        } else if (op === 1 || op === 2) {
            const r = op === 1 ? uvt.removeDuplicateVertices(inVertices, inIndices)
                : uvt.removeUnusedVertices(inVertices, inIndices);
            numOut = r.vertices.length; numOutIndices = r.indices.length;
        } else {
            const r = op === 3 ? uvt.removeDuplicateVerticesTriples(inVertices, toTriples(inIndices))
                : op === 4 ? uvt.removeUnusedVerticesTriples(inVertices, toTriples(inIndices))
                    : uvt.generateIndexedTrianglesTriples(inVertices);
            numOut = r.vertices.length; numOutIndices = 3 * r.triangles.length;
        }
        io.outInt(numOut);
        io.outInt(numOutIndices);
    }, { exact: true });

    family.case('UniqueVerticesTriangles.scalar', (io) => {
        const numVertices = io.integer();
        const inVertices = io.reals(numVertices);
        const inIndices = readIndices(io, numVertices);
        const uvt = new UniqueVerticesTriangles<number>();
        const r = uvt.removeDuplicateAndUnusedVertices(inVertices, inIndices);
        checkRemap(inVertices, inIndices, r.vertices, r.indices);
        emitUnique(io, r.vertices, r.indices);
        const triple = inVertices.slice(0, 3 * Math.floor(numVertices / 3));
        const g = triple.length > 0 ? uvt.generateIndexedTriangles(triple)
            : { vertices: [], indices: [] };
        emitUnique(io, g.vertices, g.indices);
    }, { exact: true });

    // ------------------------------------------------------------ CurveExtractor

    family.case('CurveExtractor.extract', (io) => {
        const extractor = new ReplayExtractor();
        const numEdges = io.integer();
        for (let e = 0; e < numEdges; ++e) { extractor.addRecordedEdge(io.reals(8)); }
        const numSingles = io.integer();
        for (let s = 0; s < numSingles; ++s) { extractor.addRecordedVertex(io.reals(4)); }
        const removeDuplicates = io.boolean();
        const r = extractor.extractReal(0, removeDuplicates);
        io.outInt(r.vertices.length);
        for (const v of r.vertices) { io.outReal(v[0]); io.outReal(v[1]); }
        io.outInt(r.edges.length);
        for (const e of r.edges) { io.outInt(e.v[0]); io.outInt(e.v[1]); }
    }, { exact: true });

    family.case('CurveExtractor.vertexEdge', (io) => {
        const r = io.reals(8);
        const a = new CurveExtractorVertex(r[0], r[1], r[2], r[3]);
        const b = new CurveExtractorVertex(r[4], r[5], r[6], r[7]);
        io.outInt(a.xNumer); io.outInt(a.xDenom); io.outInt(a.yNumer); io.outInt(a.yDenom);
        io.outInt(b.xNumer); io.outInt(b.xDenom); io.outInt(b.yNumer); io.outInt(b.yDenom);
        // Reference: exact rational comparison with BigInt.
        const cmp = (n0: number, d0: number, n1: number, d1: number): number => {
            const l = BigInt(n0) * BigInt(d1), rr = BigInt(n1) * BigInt(d0);
            return l < rr ? -1 : (l > rr ? 1 : 0);
        };
        const cx = cmp(a.xNumer, a.xDenom, b.xNumer, b.xDenom);
        const cy = cmp(a.yNumer, a.yDenom, b.yNumer, b.yDenom);
        check(a.equals(b) === (cx === 0 && cy === 0), 'Vertex ==');
        check(a.lessThan(b) === (cx < 0 || (cx === 0 && cy < 0)), 'Vertex <');
        io.outBool(a.equals(b));
        io.outBool(a.lessThan(b));
        io.outBool(b.lessThan(a));
        const i0 = io.integer();
        const i1 = io.integer();
        const j0 = io.integer();
        const j1 = io.integer();
        const e0 = new CurveExtractorEdge(i0, i1);
        const e1 = new CurveExtractorEdge(j0, j1);
        io.outInt(e0.v[0]); io.outInt(e0.v[1]);
        io.outBool(e0.equals(e1));
        io.outBool(e0.lessThan(e1));
        io.outBool(e1.lessThan(e0));
    }, { exact: true });

    family.case('CurveExtractor.invalidBounds', (io) => {
        const xBound = io.integer();
        const yBound = io.integer();
        new ReplayExtractor(xBound, yBound);
        io.outInt(xBound * yBound);
    }, { exact: true });

    family.finish();
});
