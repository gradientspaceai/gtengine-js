// Replays oracle/cpp/cases/v06-ap.cpp. Keep the two files in the same order.
//
// BSRational is compared by value and canonical form (OutBSR / outBSR: the
// numerator and the denominator as v05's OutBSN emits a BSNumber - sign,
// biased exponent, exponent, bit count, the 32-bit words of the odd integer
// or a digest of them, both floating-point conversions - followed by the
// rational's own conversions to double and float); the port's bigint storage
// has no word layout of its own, so the words are extracted from the bigint.
//
// Every case is exact except BSRational.std.libm, whose std:: overloads call
// the C math library. Five cases are declared deviations (#95, #168, #280).
import { describe } from 'vitest';
import { APConversion } from '../../src/APConversion.js';
import { APInterval } from '../../src/APInterval.js';
import { BSNumber, BSNumberRoundingMode } from '../../src/BSNumber.js';
import { BSPPolygon2 } from '../../src/BSPPolygon2.js';
import {
    BSRational, convertBSRational, convertBSRationalToBSNumber, convertBSRationalToFloat32,
    convertBSRationalToNumber
} from '../../src/BSRational.js';
import { BasisFunctionInput, UniqueKnot } from '../../src/BasisFunction.js';
import { BSplineCurve } from '../../src/BSplineCurve.js';
import { BSplineGeodesic } from '../../src/BSplineGeodesic.js';
import { BSplineReduction } from '../../src/BSplineReduction.js';
import { BSplineSurface } from '../../src/BSplineSurface.js';
import { EdgeKey } from '../../src/EdgeKey.js';
import { GVector } from '../../src/GVector.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ------------------------------------------------------------ BSRational

const B = (x: number): BSNumber => BSNumber.fromNumber(x);
const R = (x: number): BSRational => BSRational.fromNumber(x);
const MASK64 = (1n << 64n) - 1n;

// The replay of OutBSN (copied from v05-ap.oracle.test.ts).
function outBSN(io: OracleIO, x: BSNumber): void {
    io.outInt(x.getSign());
    io.outInt(x.getBiasedExponent());
    io.outInt(x.getExponent());
    const numBits = x.getNumBits();
    io.outInt(numBits);
    const numWords = Math.floor((numBits + 31) / 32);
    const u = x.getUInteger();
    const words: number[] = [];
    for (let i = 0; i < numWords; ++i) {
        words.push(Number((u >> BigInt(32 * i)) & 0xFFFFFFFFn));
    }
    if (numWords <= 8) {
        for (const w of words) { io.outInt(w); }
    } else {
        io.outInt(words[numWords - 1]);
        io.outInt(words[numWords - 2]);
        io.outInt(words[numWords - 3]);
        io.outInt(words[1]);
        io.outInt(words[0]);
        let h = 0xCBF29CE484222325n;
        for (const w of words) {
            for (let k = 0; k < 4; ++k) {
                h ^= BigInt((w >>> (8 * k)) & 0xFF);
                h = (h * 0x100000001B3n) & MASK64;
            }
        }
        io.outInt(Number(h & ((1n << 48n) - 1n)));
    }
    io.outReal(x.toNumber());
    io.outReal(x.toFloat32());
}

// The replay of OutBSR.
function outBSR(io: OracleIO, r: BSRational): void {
    outBSN(io, r.getNumerator());
    outBSN(io, r.getDenominator());
    io.outReal(r.toNumber());
    io.outReal(r.toFloat32());
}

// The replay of DrawBSR.
function drawBSR(io: OracleIO): BSRational {
    const mode = io.integer();
    switch (mode) {
        case 0:
            return R(io.real());
        case 1: {
            const a = io.real();
            const b = io.real();
            return BSRational.fromBSNumber(B(a), B(b));
        }
        case 2: {
            const a = io.real();
            const b = io.real();
            return R(a).div(R(b));
        }
        case 3: {
            const a = io.real();
            const b = io.real();
            return R(a).add(R(b));
        }
        case 4: {
            const a = io.real();
            const b = io.real();
            const c = io.real();
            return R(a).mul(R(b)).sub(R(c));
        }
        case 5: {
            const p = io.real();
            const q = io.real();
            return BSRational.fromBigInt(BigInt(p), BigInt(q));
        }
        default: {
            const a = io.real();
            const e = io.integer();
            const hs = io.real();
            const j = io.integer();
            const cs = io.real();
            let r = R(a).add(BSRational.ldexp(R(hs), e));
            if (j > 0) {
                r = r.add(BSRational.ldexp(R(cs), e - j).div(R(3)));
            }
            return r;
        }
    }
}

// The replay of DrawPairBSR.
function drawPairBSR(io: OracleIO): [BSRational, BSRational] {
    let x = drawBSR(io);
    const mode = io.integer();
    let y: BSRational;
    switch (mode) {
        case 0: y = drawBSR(io); break;
        case 1: y = x.clone(); break;
        case 2: y = x.negated(); break;
        case 3: y = x.add(R(io.real())); break;
        case 4: y = x.mul(R(io.real())); break;
        case 5: {
            const p = io.real();
            const q = io.real();
            const k = io.real();
            x = BSRational.fromBigInt(BigInt(p * k), BigInt(q * k));
            y = BSRational.fromBigInt(BigInt(p), BigInt(q));
            break;
        }
        default: y = x.add(R(1).div(R(io.real()))); break;
    }
    return [x, y];
}

// The replay of RecordBits: a 64-bit pattern from its two 32-bit words.
function readBits(io: OracleIO): bigint {
    const hi = BigInt(io.integer());
    const lo = BigInt(io.integer());
    return (hi << 32n) | lo;
}

// The value of a pattern for the integer kind (0 int32_t, 1 uint32_t,
// 2 int64_t, 3 uint64_t).
function kindValue(bits: bigint, kind: number): bigint {
    if (kind === 0) { return BigInt.asIntN(32, bits); }
    if (kind === 1) { return BigInt.asUintN(32, bits); }
    if (kind === 2) { return BigInt.asIntN(64, bits); }
    return BigInt.asUintN(64, bits);
}

// The port's constructor for a C++ integer constructor: the 32-bit ones are
// ported by fromNumber (exact), the 64-bit ones by fromBigInt.
function fromKind(kind: number, n: bigint, d?: bigint): BSRational {
    if (kind <= 1) {
        return BSRational.fromNumber(Number(n), d === undefined ? undefined : Number(d));
    }
    return BSRational.fromBigInt(n, d);
}

function readString(io: OracleIO): string {
    const n = io.integer();
    let s = '';
    for (let i = 0; i < n; ++i) { s += String.fromCharCode(io.integer()); }
    return s;
}

const roundingModes: BSNumberRoundingMode[] = [
    BSNumberRoundingMode.FE_TONEAREST, BSNumberRoundingMode.FE_DOWNWARD,
    BSNumberRoundingMode.FE_TOWARDZERO, BSNumberRoundingMode.FE_UPWARD,
    12345 as BSNumberRoundingMode
];

// The replay of DrawBSN (v06): the mode, two doubles.
function drawBSN(io: OracleIO, nonzero: boolean): BSNumber {
    const mode = io.integer();
    const a = io.real();
    const b = io.real();
    if (mode === 0) { return B(a); }
    if (mode === 1) {
        const s = B(a).add(B(b));
        return (nonzero && s.getSign() === 0 ? B(a) : s);
    }
    return B(a).mul(B(b));
}

// ---------------------------------------------------------- APConversion

// The replay of DrawSquare.
function drawSquare(io: OracleIO): BSRational {
    const mode = io.integer();
    const p = io.real();
    const q = io.real();
    const k = io.real();
    const r = R(p).div(R(q));
    switch (mode) {
        case 0: return r.mul(r);
        case 1: return r;
        case 2: return BSRational.ldexp(r, k);
        case 3: {
            const s = BSRational.ldexp(R(p), Math.trunc(k / 2));
            return s.mul(s);
        }
        default: return R(0);
    }
}

// The replay of DrawAmB.
function drawAmB(io: OracleIO): [BSRational, BSRational] {
    const mode = io.integer();
    if (mode === 0) {
        let aSqr = drawSquare(io);
        let bSqr = drawSquare(io);
        if (aSqr.lessThan(bSqr)) { [aSqr, bSqr] = [bSqr, aSqr]; }
        return [aSqr, bSqr];
    }
    if (mode === 1) {
        const P = io.real();
        const Q = io.real();
        return [R(P), R(Q)];
    }
    const aSqr = drawSquare(io);
    const bSqr = (io.boolean() ? aSqr.clone() : R(0));
    return [aSqr, bSqr];
}

function outAmB(io: OracleIO, aSqr: BSRational, bSqr: BSRational): void {
    const precision = io.real();
    const maxIterations = io.real();
    const apc = new APConversion(precision, maxIterations);
    const r = apc.estimateAmB(aSqr, bSqr);
    io.outInt(r.numIterates);
    outBSR(io, r.tMin);
    outBSR(io, r.tMax);
}

// ------------------------------------------------------------ APInterval

// The replay of OutAPE: a finite endpoint in full, an infinity sentinel by
// its sign, its numerator's sign and biased exponent and its denominator.
function outAPE(io: OracleIO, e: BSRational): void {
    const s = e.getSign();
    io.outInt(s);
    if (s >= 2 || s <= -2) {
        io.outInt(e.getNumerator().getSign());
        io.outInt(e.getNumerator().getBiasedExponent());
        outBSN(io, e.getDenominator());
        return;
    }
    outBSR(io, e);
}

function outAPI(io: OracleIO, w: APInterval): void {
    outAPE(io, w.get(0));
    outAPE(io, w.get(1));
}

function drawAPScalar(io: OracleIO): BSRational {
    return (io.integer() === 0 ? R(0) : drawBSR(io));
}

// The replay of DrawAPI.
function drawAPI(io: OracleIO): APInterval {
    const shape = io.integer();
    const a = drawBSR(io);
    const b = drawBSR(io);
    const lo = (a.lessThan(b) ? a : b);
    const up = (a.lessThan(b) ? b : a);
    const fa = BSRational.fabs(a);
    const fb = BSRational.fabs(b).add(R(1));
    const one = R(1);
    switch (shape) {
        case 0: return new APInterval(lo, up);
        case 1: return new APInterval(R(0), fb);
        case 2: return new APInterval(fb.negated(), R(0));
        case 3: return new APInterval(a, a);
        case 4: return new APInterval(fa.negated().sub(one), fb);
        case 5: return new APInterval(R(0), R(0));
        case 6: return new APInterval(fa.add(one), fa.add(fb).add(one));
        default: return new APInterval(fa.negated().sub(fb).sub(one), fa.negated().sub(one));
    }
}

function isZeroInterval(v: APInterval): boolean {
    return v.get(0).getSign() === 0 && v.get(1).getSign() === 0;
}

// ----------------------------------------------------------- BSPPolygon2

interface Soup {
    v: Vector[];
    e: [number, number][];
}

function readSoup(io: OracleIO): Soup {
    const nv = io.integer();
    const v: Vector[] = [];
    for (let i = 0; i < nv; ++i) { v.push(io.vec(2)); }
    const ne = io.integer();
    const e: [number, number][] = [];
    for (let i = 0; i < ne; ++i) {
        const a = io.integer();
        e.push([a, io.integer()]);
    }
    return { v, e };
}

const edge = (v0: number, v1: number): EdgeKey => new EdgeKey(true, v0, v1);

function buildPoly(poly: BSPPolygon2, s: Soup, finalize: boolean): void {
    for (const [a, b] of s.e) {
        const v0 = poly.insertVertex(s.v[a]);
        const v1 = poly.insertVertex(s.v[b]);
        poly.insertEdge(edge(v0, v1));
    }
    if (finalize) { poly.finalize(); }
}

// The private state of the port's BSPPolygon2 and BSPTree2, read the way
// the C++ case reads upstream's (with 'private' opened up).
interface TreeState {
    mCoincident: EdgeKey[];
    mPosChild: TreeState | null;
    mNegChild: TreeState | null;
}

interface PolyState {
    mEArray: EdgeKey[];
    mTree: TreeState | null;
}

function outTree(io: OracleIO, node: TreeState): void {
    io.outInt(node.mCoincident.length);
    for (const e of node.mCoincident) { io.outInt(e.V[0]); io.outInt(e.V[1]); }
    io.outBool(node.mPosChild !== null);
    if (node.mPosChild !== null) { outTree(io, node.mPosChild); }
    io.outBool(node.mNegChild !== null);
    if (node.mNegChild !== null) { outTree(io, node.mNegChild); }
}

function outPoly(io: OracleIO, poly: BSPPolygon2): void {
    io.outInt(poly.getNumVertices());
    for (let i = 0; i < poly.getNumVertices(); ++i) {
        io.outReal(poly.getVertex(i).get(0));
        io.outReal(poly.getVertex(i).get(1));
    }
    io.outInt(poly.getNumEdges());
    const state = poly as unknown as PolyState;
    io.outInt(state.mEArray.length);
    for (const e of state.mEArray) { io.outInt(e.V[0]); io.outInt(e.V[1]); }
    io.outBool(state.mTree !== null);
    if (state.mTree !== null) { outTree(io, state.mTree); }
}

function readQueries(io: OracleIO): Vector[] {
    const pts: Vector[] = [];
    for (let i = 0; i < 8; ++i) { pts.push(io.vec(2)); }
    return pts;
}

function polyOp(P: BSPPolygon2, Q: BSPPolygon2, op: number): BSPPolygon2 {
    switch (op) {
        case 0: return P.negated();
        case 1: return P.intersection(Q);
        case 2: return P.union(Q);
        case 3: return P.difference(Q);
        default: return P.exclusiveOr(Q);
    }
}

// The replay of PolyBoolean.
function polyBoolean(io: OracleIO, op: number): void {
    const eps = io.real();
    const sp = readSoup(io);
    const sq = readSoup(io);
    const pts = readQueries(io);
    const P = new BSPPolygon2(eps);
    const Q = new BSPPolygon2(eps);
    buildPoly(P, sp, true);
    buildPoly(Q, sq, true);
    const result = polyOp(P, Q, op);
    outPoly(io, result);
    for (const x of pts) { io.outInt(result.pointLocation(x)); }
}

// ------------------------------------------------------------- B-splines

// The replay of RecordBasis (copied from v05-ap.oracle.test.ts).
function basisInput(io: OracleIO): BasisFunctionInput {
    const input = new BasisFunctionInput();
    input.numControls = io.integer();
    input.degree = io.integer();
    input.uniform = io.boolean();
    input.periodic = io.boolean();
    input.numUniqueKnots = io.integer();
    input.uniqueKnots = [];
    for (let i = 0; i < input.numUniqueKnots; ++i) {
        const t = io.real();
        input.uniqueKnots.push(new UniqueKnot(t, io.integer()));
    }
    return input;
}

function points(io: OracleIO, count: number, dimension: number): Vector[] {
    const P = new Array<Vector>(count);
    for (let i = 0; i < count; ++i) { P[i] = io.vec(dimension); }
    return P;
}

function gvec2(io: OracleIO): GVector {
    const p = new GVector(2);
    p.values[0] = io.real();
    p.values[1] = io.real();
    return p;
}

function outGVec2(io: OracleIO, v: GVector): void {
    io.outInt(v.size);
    for (let i = 0; i < v.size; ++i) { io.outReal(v.values[i]); }
}

// The replay of DrawGeoSurface: the two bases and the controls.
function geoParts(io: OracleIO): [BasisFunctionInput, BasisFunctionInput, Vector[]] {
    const in0 = basisInput(io);
    const in1 = basisInput(io);
    const controls = points(io, in0.numControls * in1.numControls, 3);
    return [in0, in1, controls];
}

function geoSurface(io: OracleIO): BSplineSurface {
    const [in0, in1, controls] = geoParts(io);
    return new BSplineSurface(3, [in0, in1], controls);
}

// The replay of DrawTuning. The derived steps are not recomputed
// (updateDerivedParameters is not called), as upstream computes them only
// in the constructor.
function tuning(io: OracleIO, geo: BSplineGeodesic): void {
    geo.integralSamples = io.integer();
    geo.searchSamples = io.integer();
    geo.derivativeStep = io.real();
    geo.subdivisions = io.integer();
    geo.refinements = io.integer();
    geo.searchRadius = io.real();
}

describe('oracle: v06-ap', () => {
    const family = new OracleFamily('v06-ap');

    // --------------------------------------------------------- BSRational

    family.case('BSRational.construct.double', (io) => {
        const kind = io.integer();
        const n = io.real();
        if (kind === 0) {
            outBSR(io, R(n));
            return;
        }
        const d = io.real();
        outBSR(io, BSRational.fromNumber(n, d));
    }, { exact: true });

    family.case('BSRational.construct.float', (io) => {
        const kind = io.integer();
        const n = io.real();
        if (kind === 0) {
            outBSR(io, BSRational.fromFloat32(n));
            return;
        }
        const d = io.real();
        outBSR(io, BSRational.fromFloat32(n, d));
    }, { exact: true });

    family.case('BSRational.construct.integer', (io) => {
        const kind = io.integer();
        const pair = io.boolean();
        const n = kindValue(readBits(io), kind);
        if (!pair) {
            outBSR(io, fromKind(kind, n));
            return;
        }
        const d = kindValue(readBits(io), kind);
        outBSR(io, fromKind(kind, n, d));
    }, { exact: true });

    family.case('BSRational.construct.bsnumber', (io) => {
        const kind = io.integer();
        const n = drawBSN(io, kind !== 0);
        if (kind === 0) {
            outBSR(io, BSRational.fromBSNumber(n));
            return;
        }
        const d = (io.integer() === 0 ? new BSNumber() : drawBSN(io, true));
        outBSR(io, BSRational.fromBSNumber(n, d));
    }, { exact: true });

    // Finding #168 (port fixed): every record deviates in the numerator's
    // biased exponent (upstream: minus the denominator's exponent; port: 0).
    family.case('BSRational.construct.bsnumber.zeroNumerator', (io) => {
        const z = io.real();
        const d = io.real();
        outBSR(io, BSRational.fromBSNumber(B(z), B(d)));
    }, { exact: true, deviation: '#168 (UPSTREAM-FINDINGS BSRational.h item 1)' });

    family.case('BSRational.construct.string', (io) => {
        outBSR(io, BSRational.fromString(readString(io)));
    }, { exact: true });

    // Finding #168 (port fixed): upstream gives a zero-valued string a
    // numerator with sign +-1; the port keeps the canonical zero.
    family.case('BSRational.construct.string.zero', (io) => {
        const r = BSRational.fromString(readString(io));
        io.outInt(r.getSign());
        io.outInt(r.getNumerator().getSign());
        io.outInt(r.getNumerator().getNumBits());
    }, { exact: true, deviation: '#168 (UPSTREAM-FINDINGS BSRational.h item 2)' });

    // Finding #95 (port fixed), inherited from BSNumber::ConvertToInteger:
    // upstream accepts a one-character non-digit integer part, the port
    // asserts.
    family.case('BSRational.construct.string.singleChar', (io) => {
        outBSR(io, BSRational.fromString(readString(io)));
    }, { exact: true, deviation: '#95 (UPSTREAM-FINDINGS BSNumber.h)' });

    family.case('BSRational.compare', (io) => {
        const [x, y] = drawPairBSR(io);
        io.outBool(x.equals(y));
        io.outBool(x.notEquals(y));
        io.outBool(x.lessThan(y));
        io.outBool(x.lessThanOrEqual(y));
        io.outBool(x.greaterThan(y));
        io.outBool(x.greaterThanOrEqual(y));
        io.outBool(y.equals(x));
        io.outBool(y.lessThan(x));
        io.outBool(y.lessThanOrEqual(x));
        io.outInt(x.getSign());
        io.outInt(y.getSign());
    }, { exact: true });

    family.case('BSRational.arithmetic', (io) => {
        const [x, y] = drawPairBSR(io);
        outBSR(io, x.clone());
        outBSR(io, x.negated());
        outBSR(io, x.add(y));
        outBSR(io, y.add(x));
        outBSR(io, x.sub(y));
        outBSR(io, y.sub(x));
        outBSR(io, x.mul(y));
        outBSR(io, y.mul(x));
        if (y.getSign() !== 0) { outBSR(io, x.div(y)); }
        if (x.getSign() !== 0) { outBSR(io, y.div(x)); }
        let z = x.clone();
        z = z.add(y);
        outBSR(io, z);
        z = z.sub(x);
        outBSR(io, z);
        z = z.mul(x);
        outBSR(io, z);
        if (y.getSign() !== 0) {
            z = z.div(y);
            outBSR(io, z);
        }
        const w = y.clone();
        w.negate();
        outBSR(io, w);
        if (w.getSign() !== 0) {
            w.setSign(io.integer() !== 0 ? 1 : -1);
            outBSR(io, w);
        }
    }, { exact: true });

    family.case('BSRational.divide.zero', (io) => {
        const x = drawBSR(io);
        const how = io.integer();
        let zero = new BSRational();
        if (how === 1) { zero = R(-0.0); }
        else if (how === 2) { zero = x.sub(x); }
        else if (how === 3) { zero = R(0).mul(x); }
        if (io.boolean()) {
            outBSR(io, x.div(zero));
        } else {
            outBSR(io, x.div(zero));
        }
    }, { exact: true });

    family.case('BSRational.convert', (io) => {
        const x = drawBSR(io);
        const precision = io.integer();
        const mode = roundingModes[io.integer()];
        const outN = convertBSRationalToBSNumber(x, precision, mode);
        const outR = convertBSRational(x, precision, mode);
        const outD = convertBSRationalToNumber(x, mode);
        const outF = convertBSRationalToFloat32(x, mode);
        outBSN(io, outN);
        outBSR(io, outR);
        io.outReal(outD);
        io.outReal(outF);
    }, { exact: true });

    // Upstream's two-step rounding in the subnormal range and beyond the
    // largest finite value is preserved (see the C++ case comment).
    family.case('BSRational.convert.subnormalAndOverflow', (io) => {
        io.integer();  // the range (only the C++ generator uses it)
        const p = io.real();
        const q = io.real();
        const e = io.integer();
        const tail = io.integer();
        const s = io.real();
        const mode = roundingModes[io.integer()];
        let x = BSRational.ldexp(R(p).div(R(q)), e);
        if (tail > 0) { x = x.add(BSRational.ldexp(R(1), e - tail)); }
        x = R(s).mul(x);
        io.outReal(convertBSRationalToNumber(x, mode));
        io.outReal(convertBSRationalToFloat32(x, mode));
        io.outReal(x.toNumber());
        io.outReal(x.toFloat32());
    }, { exact: true });

    family.case('BSRational.std.exact', (io) => {
        const x = drawBSR(io);
        const y = drawBSR(io);
        const z = drawBSR(io);
        const w = drawBSR(io);
        const k = io.integer();
        outBSR(io, BSRational.fabs(x));
        const fr = BSRational.frexp(x);
        outBSR(io, fr.result);
        io.outInt(fr.exponent);
        outBSR(io, BSRational.ldexp(x, k));
        outBSR(io, BSRational.floor(x));
        outBSR(io, BSRational.ceil(x));
        outBSR(io, BSRational.sqrt(x));
        outBSR(io, BSRational.fmod(x, y));
        outBSR(io, BSRational.remainder(x, y));
        outBSR(io, BSRational.clamp(x, y, z));
        outBSR(io, BSRational.invsqrt(x));
        io.outInt(BSRational.isign(x));
        outBSR(io, BSRational.saturate(x));
        outBSR(io, BSRational.sign(x));
        outBSR(io, BSRational.sqr(x));
        outBSR(io, BSRational.fma(x, y, z));
        outBSR(io, BSRational.robustSOP(x, y, z, w));
        outBSR(io, BSRational.robustDOP(x, y, z, w));
    }, { exact: true });

    family.case('BSRational.std.remainder', (io) => {
        const a = io.real();
        const b = io.real();
        const c = io.real();
        const x = R(a).div(R(3));
        const y = R(b);
        const u = R(c);
        outBSR(io, BSRational.remainder(x, y));
        outBSR(io, BSRational.remainder(y, x));
        outBSR(io, BSRational.fmod(x, y));
        outBSR(io, BSRational.remainder(y, u));
        outBSR(io, BSRational.fmod(y, u));
        outBSR(io, BSRational.remainder(y, u.mul(R(4))));
    }, { exact: true });

    // The C math library (std::acos ... std::tanh, gte::atandivpi ...
    // sinpi) of MSVC and V8 may differ by an ulp; default tolerance.
    family.case('BSRational.std.libm', (io) => {
        const v: BSRational[] = [];
        for (let i = 0; i < 2; ++i) {
            const p = io.real();
            const q = io.real();
            v.push(R(p).div(R(q)));
        }
        const [x, y] = v;
        const out = (r: BSRational): void => io.outReal(r.toNumber());
        out(BSRational.acos(x));
        out(BSRational.acosh(x));
        out(BSRational.asin(x));
        out(BSRational.asinh(x));
        out(BSRational.atan(x));
        out(BSRational.atanh(x));
        out(BSRational.atan2(y, x));
        out(BSRational.cos(x));
        out(BSRational.cosh(x));
        out(BSRational.exp(x));
        out(BSRational.exp2(x));
        out(BSRational.log(x));
        out(BSRational.log2(x));
        out(BSRational.log10(x));
        out(BSRational.pow(x, y));
        out(BSRational.sin(x));
        out(BSRational.sinh(x));
        out(BSRational.tan(x));
        out(BSRational.tanh(x));
        out(BSRational.atandivpi(x));
        out(BSRational.atan2divpi(y, x));
        out(BSRational.cospi(x));
        out(BSRational.exp10(x));
        out(BSRational.sinpi(x));
    });

    // ------------------------------------------------------- APConversion

    family.case('APConversion.estimateSqrt', (io) => {
        let aSqr = drawSquare(io);
        if (io.integer() === 0) { aSqr = aSqr.negated().sub(R(1)); }
        const precision = io.real();
        const maxIterations = io.real();
        const apc = new APConversion(precision, maxIterations);
        const r0 = apc.estimateSqrt(aSqr);
        const r1 = apc.estimateSqrtValue(aSqr);
        io.outInt(r0.numIterates);
        outBSR(io, r0.aMin);
        outBSR(io, r0.aMax);
        io.outInt(r1.numIterates);
        outBSR(io, r1.a);
    }, { exact: true });

    family.case('APConversion.estimateApB', (io) => {
        const aSqr = drawSquare(io);
        const bSqr = drawSquare(io);
        const precision = io.real();
        const maxIterations = io.real();
        const apc = new APConversion(precision, maxIterations);
        const r = apc.estimateApB(aSqr, bSqr);
        io.outInt(r.numIterates);
        outBSR(io, r.tMin);
        outBSR(io, r.tMax);
    }, { exact: true });

    family.case('APConversion.estimateAmB', (io) => {
        const [aSqr, bSqr] = drawAmB(io);
        outAmB(io, aSqr, bSqr);
    }, { exact: true });

    family.case('APConversion.estimateAmB.aLessThanB', (io) => {
        const aSqr = drawSquare(io);
        const bSqr = drawSquare(io);
        outAmB(io, aSqr, bSqr);
    }, { exact: true });

    // Finding #280 item 1 (port fixed): after an exhausted bisection upstream
    // can return a bracket that misses a - b; the port returns the bisection
    // bracket exactly then. Every record deviates.
    family.case('APConversion.estimateAmB.bisectionExhausted', (io) => {
        const P = io.real();
        const Q = io.real();
        outAmB(io, R(P), R(Q));
    }, { exact: true, deviation: '#280 (UPSTREAM-FINDINGS APConversion.h item 1)' });

    family.case('APConversion.estimate', (io) => {
        const x = drawBSR(io);
        const ymode = io.integer();
        const y = (ymode === 0 ? R(0) : drawBSR(io));
        const dmode = io.integer();
        let d = (dmode === 0 ? R(0) : drawSquare(io));
        if (dmode === 5) { d = d.negated().sub(R(1)); }
        const precision = io.real();
        const maxIterations = io.real();
        const apc = new APConversion(precision, maxIterations);
        const q = { x: [x, y] as const, d };
        const r0 = apc.estimate(q);
        const r1 = apc.estimateValue(q);
        io.outInt(r0.numIterates);
        outBSR(io, r0.qMin);
        outBSR(io, r0.qMax);
        io.outInt(r1.numIterates);
        outBSR(io, r1.qEstimate);
    }, { exact: true });

    family.case('APConversion.accessors', (io) => {
        const p0 = io.integer();
        const m0 = io.integer();
        const p1 = io.integer();
        const m1 = io.integer();
        const apc = new APConversion(p0, m0);
        const g0 = apc.getPrecision();
        const h0 = apc.getMaxIterations();
        const r0 = apc.estimateSqrt(R(2));
        apc.setPrecision(p1);
        apc.setMaxIterations(m1);
        const r1 = apc.estimateSqrt(R(2));
        io.outInt(g0);
        io.outInt(h0);
        io.outInt(r0.numIterates);
        outBSR(io, r0.aMin);
        outBSR(io, r0.aMax);
        io.outInt(apc.getPrecision());
        io.outInt(apc.getMaxIterations());
        io.outInt(r1.numIterates);
        outBSR(io, r1.aMin);
        outBSR(io, r1.aMax);
    }, { exact: true });

    // --------------------------------------------------------- APInterval

    family.case('APInterval.leaf', (io) => {
        const u = drawAPScalar(io);
        const v = drawAPScalar(io);
        outAPI(io, APInterval.add(u, v));
        outAPI(io, APInterval.sub(u, v));
        outAPI(io, APInterval.mul(u, v));
        outAPI(io, APInterval.div(u, v));
        outAPI(io, APInterval.div(v, u));
    }, { exact: true });

    family.case('APInterval.internal', (io) => {
        const u0 = drawAPScalar(io);
        const u1 = drawAPScalar(io);
        let v0 = drawBSR(io);
        let v1 = drawBSR(io);
        if (v0.getSign() === 0) { v0 = R(3); }
        if (v1.getSign() === 0) { v1 = R(-5); }
        outAPI(io, APInterval.add(u0, u1, v0, v1));
        outAPI(io, APInterval.sub(u0, u1, v0, v1));
        outAPI(io, APInterval.mul(u0, u1, v0, v1));
        outAPI(io, APInterval.mul2(u0, u1, v0, v1));
        outAPI(io, APInterval.div(u0, u1, v0, v1));
        outAPI(io, APInterval.reciprocal(v0, v1));
        outAPI(io, APInterval.reciprocalDown(v0));
        outAPI(io, APInterval.reciprocalUp(v1));
        outAPI(io, APInterval.reals());
    }, { exact: true });

    family.case('APInterval.operators', (io) => {
        const u = drawAPI(io);
        const v = drawAPI(io);
        const s = drawAPScalar(io);
        const vZero = isZeroInterval(v);
        const uZero = isZeroInterval(u);

        outAPI(io, new APInterval());
        outAPI(io, new APInterval(s));
        outAPI(io, APInterval.fromEndpoints(u.getEndpoints()));
        outAPI(io, u.clone());
        outAPE(io, u.get(0));
        outAPE(io, u.get(1));

        outAPI(io, u.clone());
        outAPI(io, u.negate());
        outAPI(io, u.add(s));
        outAPI(io, u.add(s));
        outAPI(io, u.add(v));
        outAPI(io, APInterval.scalarSub(s, u));
        outAPI(io, u.sub(s));
        outAPI(io, u.sub(v));
        outAPI(io, u.mul(s));
        outAPI(io, u.mul(s));
        outAPI(io, u.mul(v));
        outAPI(io, v.mul(u));
        outAPI(io, u.div(s));
        if (!vZero) {
            const q = u.div(v);
            outAPI(io, q);
            outAPI(io, q.negate());
            outAPI(io, APInterval.scalarDiv(s, v));
            outAPI(io, APInterval.scalarDiv(s, v).negate());
        }
        if (!uZero) {
            outAPI(io, v.div(u));
        }

        let w = u.clone();
        w = w.add(s);
        outAPI(io, w);
        w = w.add(v);
        outAPI(io, w);
        w = w.sub(s);
        outAPI(io, w);
        w = w.sub(v);
        outAPI(io, w);
        w = w.mul(s);
        outAPI(io, w);
        w = w.mul(v);
        outAPI(io, w);
        w = w.div(s);
        outAPI(io, w);
        const wFinite = Math.abs(w.get(0).getSign()) < 2 && Math.abs(w.get(1).getSign()) < 2;
        if (!vZero && wFinite) {
            w = w.div(v);
            outAPI(io, w);
        }
    }, { exact: true });

    family.case('APInterval.divide.zeroInterval', (io) => {
        const u = drawAPI(io);
        const s = drawAPScalar(io);
        const how = io.integer();
        const zero = new APInterval(R(0), R(0));
        if (how === 0) { outAPI(io, u.div(zero)); }
        else if (how === 1) { outAPI(io, APInterval.scalarDiv(s, zero)); }
        else { outAPI(io, u.div(zero)); }
    }, { exact: true });

    // -------------------------------------------------------- BSPPolygon2

    family.case('BSPPolygon2.negation', (io) => { polyBoolean(io, 0); }, { exact: true });
    family.case('BSPPolygon2.intersection', (io) => { polyBoolean(io, 1); }, { exact: true });
    family.case('BSPPolygon2.union', (io) => { polyBoolean(io, 2); }, { exact: true });
    family.case('BSPPolygon2.difference', (io) => { polyBoolean(io, 3); }, { exact: true });
    family.case('BSPPolygon2.exclusiveOr', (io) => { polyBoolean(io, 4); }, { exact: true });

    family.case('BSPPolygon2.construct', (io) => {
        const eps = io.real();
        const sp = readSoup(io);
        const pts = readQueries(io);
        const P = new BSPPolygon2(eps);
        const indices: number[] = [];
        for (const [a, b] of sp.e) {
            const v0 = P.insertVertex(sp.v[a]);
            const v1 = P.insertVertex(sp.v[b]);
            indices.push(v0, v1, P.insertEdge(edge(v0, v1)));
        }
        const pre = P.clone();
        P.finalize();
        const locations = pts.map((x) => P.pointLocation(x));
        const C = P.clone();
        const D = C.clone();
        for (const i of indices) { io.outInt(i); }
        outPoly(io, pre);
        outPoly(io, P);
        for (const loc of locations) { io.outInt(loc); }
        outPoly(io, D);
        for (let i = 0; i < D.getNumEdges(); ++i) {
            io.outInt(D.getEdge(i).V[0]);
            io.outInt(D.getEdge(i).V[1]);
        }
    }, { exact: true });

    family.case('BSPPolygon2.splitEdge.unmapped', (io) => {
        const eps = io.real();
        const s = readSoup(io);
        const sq = readSoup(io);
        const pts = readQueries(io);
        const op = io.integer();
        const A = new BSPPolygon2(eps);
        const B2 = new BSPPolygon2(eps);
        buildPoly(A, s, true);
        buildPoly(B2, sq, true);
        outPoly(io, A);
        for (const x of pts) { io.outInt(A.pointLocation(x)); }
        let result: BSPPolygon2 | null = null;
        try {
            switch (op) {
                case 0: result = A.intersection(B2); break;
                case 1: result = A.union(B2); break;
                case 2: result = A.difference(B2); break;
                default: result = B2.difference(A); break;
            }
        } catch {
            result = null;
        }
        io.outBool(result === null);
        if (result !== null) {
            outPoly(io, result);
            for (const x of pts) { io.outInt(result.pointLocation(x)); }
        }
    }, { exact: true });

    family.case('BSPPolygon2.invalid', (io) => {
        const kind = io.integer();
        const sp = readSoup(io);
        const P = new BSPPolygon2(0);
        if (kind === 0) {
            buildPoly(P, sp, false);
            const v = P.insertVertex(sp.v[0]);
            P.insertEdge(edge(v, v));
            io.outInt(v);
        } else if (kind === 1) {
            buildPoly(P, sp, false);
            io.outInt(P.pointLocation(sp.v[0]));
        } else if (kind === 2) {
            buildPoly(P, sp, false);
            io.outInt(P.negated().getNumEdges());
        } else if (kind === 3) {
            P.finalize();
            io.outInt(P.getNumEdges());
        } else {
            buildPoly(P, sp, true);
            const far: Soup = { v: sp.v.map((v) => Vector.fromArray([v.get(0) + 100, v.get(1)])), e: sp.e };
            const F = new BSPPolygon2(0);
            buildPoly(F, far, true);
            io.outInt(P.intersection(F).getNumEdges());
        }
    }, { exact: true });

    // ------------------------------------------------------- BSplineCurve

    family.case('BSplineCurve.evaluate', (io) => {
        const N = io.integer();
        const input = basisInput(io);
        const cmode = io.integer();
        const n = input.numControls;
        const controls = (cmode <= 4 ? points(io, n, N) : undefined);
        const curve = new BSplineCurve(N, input, controls);
        for (let k = 0; k < 3; ++k) {
            const i = io.integer();
            curve.setControl(i, io.vec(N));
        }
        const basis = curve.getBasisFunction();
        io.outBool(curve.isConstructed());
        io.outReal(curve.getTMin());
        io.outReal(curve.getTMax());
        io.outInt(curve.getNumSegments());
        io.outInt(curve.getNumControls());
        io.outInt(basis.getNumKnots());
        io.outInt(basis.getDegree());
        io.outBool(basis.isOpen());
        io.outBool(basis.isUniform());
        io.outBool(basis.isPeriodic());
        io.outReal(basis.getMinDomain());
        io.outReal(basis.getMaxDomain());
        for (let k = 0; k < 3; ++k) {
            io.outVec(curve.getControl(io.integer()));
        }
        io.outVec(curve.getControls()[n - 1]);
        for (let k = 0; k < 5; ++k) {
            const t = io.real();
            const order = io.integer();
            const jet = curve.createJet();
            curve.evaluate(t, order, jet);
            for (let i = 0; i < 4; ++i) { io.outVec(jet[i]); }
        }
        const t = io.real();
        io.outVec(curve.getPosition(t));
        io.outVec(curve.getTangent(t));
        io.outReal(curve.getSpeed(t));
        const tl = io.real();
        io.outReal(curve.getLength(curve.getTMin(), tl));
        const total = curve.getTotalLength();
        io.outReal(total);
        const fraction = io.real();
        io.outReal(curve.getTime(fraction * total));
    }, { exact: true });

    // --------------------------------------------------- BSplineReduction

    family.case('BSplineReduction.reduce', (io) => {
        const N = io.integer();
        const numIn = io.integer();
        const degree = io.integer();
        const fraction = io.real();
        io.integer();  // the control mode (only the C++ generator uses it)
        const input = points(io, numIn, N);
        const reduction = new BSplineReduction();
        const out = reduction.compute(input, degree, fraction);
        io.outInt(out.length);
        for (const c of out) { io.outVec(c); }
        const reversed = [...input].reverse();
        const out2 = reduction.compute(reversed, degree, 0.5);
        io.outInt(out2.length);
        for (const c of out2) { io.outVec(c); }
    }, { exact: true });

    family.case('BSplineReduction.invalid', (io) => {
        const kind = io.integer();
        const numIn = io.integer();
        const degree = io.integer();
        void kind;
        const input = points(io, numIn, 2);
        const out = new BSplineReduction().compute(input, degree, 0.5);
        io.outInt(out.length);
    }, { exact: true });

    // ---------------------------------------------------- BSplineGeodesic

    family.case('BSplineGeodesic.computeGeodesic', (io) => {
        const surface = geoSurface(io);
        const geo = new BSplineGeodesic(surface);
        tuning(io, geo);
        const end0 = gvec2(io);
        const end1 = gvec2(io);
        let calls = 0;
        let progress = 0;
        geo.refineCallback = () => {
            ++calls;
            progress = progress * 7 + geo.getSubdivisionStep() * 100
                + geo.getRefinementStep() * 10 + geo.getCurrentQuantity();
            progress %= 1000000007;
        };
        const { quantity, path } = geo.computeGeodesic(end0, end1);
        const length = geo.computeTotalLength(quantity, path);
        const curvature = geo.computeTotalCurvature(quantity, path);
        io.outInt(geo.getDimension());
        io.outInt(quantity);
        io.outInt(path.length);
        for (const p of path) { outGVec2(io, p); }
        io.outInt(calls);
        io.outInt(progress);
        io.outInt(geo.getSubdivisionStep());
        io.outInt(geo.getRefinementStep());
        io.outInt(geo.getCurrentQuantity());
        io.outReal(length);
        io.outReal(curvature);
    }, { exact: true, timeout: 600000 });

    family.case('BSplineGeodesic.segment', (io) => {
        const surface = geoSurface(io);
        const geo = new BSplineGeodesic(surface);
        tuning(io, geo);
        const p0 = gvec2(io);
        const p1 = gvec2(io);
        const m = gvec2(io);
        let calls = 0;
        geo.refineCallback = () => { ++calls; };
        const length = geo.computeSegmentLength(p0, p1);
        const curvature = geo.computeSegmentCurvature(p0, p1);
        const sub = geo.subdivide(p0, p1);
        const ref = geo.refine(p0, m, p1);
        io.outReal(length);
        io.outReal(curvature);
        io.outBool(sub.changed);
        outGVec2(io, sub.mid);
        io.outBool(ref.changed);
        outGVec2(io, ref.mid);
        io.outInt(calls);
    }, { exact: true, timeout: 600000 });

    family.case('BSplineGeodesic.degenerate', (io) => {
        const [in0, in1, controls] = geoParts(io);
        const kind = io.integer();
        const surface = new BSplineSurface(3, [in0, in1], kind === 0 ? controls : undefined);
        const geo = new BSplineGeodesic(surface);
        const p0 = gvec2(io);
        const p1 = (kind === 0 ? p0.clone() : gvec2(io));
        if (kind === 1 && p1.values[0] === p0.values[0] && p1.values[1] === p0.values[1]) {
            p1.values[0] += 0.125;
        }
        io.outReal(geo.computeSegmentLength(p0, p1));
    }, { exact: true });

    family.finish();
});
