// Replays oracle/cpp/cases/v05-ap.cpp. Keep the two files in the same order.
//
// BSNumber is compared by value and canonical form (OutBSN / outBSN: sign,
// biased exponent, exponent, bit count, the 32-bit words of the odd integer
// or a digest of them, and both floating-point conversions); the port's
// bigint storage has no word layout of its own, so the words are extracted
// from the bigint.
//
// Every case is exact except BSNumber.std.libm, whose std:: overloads call
// the C math library (acos ... tanh, and gte::atandivpi, atan2divpi, cospi,
// exp10, sinpi built on them).
import { describe } from 'vitest';
import { BSNumber, BSNumberRoundingMode, convertBSNumber } from '../../src/BSNumber.js';
import { BSPrecision, BSPrecisionParameters, BSPrecisionType } from '../../src/BSPrecision.js';
import { QFNumber } from '../../src/QFNumber.js';
import { SWInterval } from '../../src/SWInterval.js';
import { BasisFunctionInput, UniqueKnot } from '../../src/BasisFunction.js';
import { BSplineCurveFit } from '../../src/BSplineCurveFit.js';
import { BSplineSurface } from '../../src/BSplineSurface.js';
import { BSplineSurfaceFit } from '../../src/BSplineSurfaceFit.js';
import { BSplineVolume } from '../../src/BSplineVolume.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---------------------------------------------------------------- BSNumber

const B = (x: number): BSNumber => BSNumber.fromNumber(x);
const MASK64 = (1n << 64n) - 1n;

// The replay of OutBSN.
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

// The replay of DrawBSN: the mode, then 1..3 recorded doubles.
function drawBSN(io: OracleIO): BSNumber {
    const mode = io.integer();
    switch (mode) {
        case 0:
            return B(io.real());
        case 1:
        case 4: {
            const a = io.real();
            const b = io.real();
            return B(a).add(B(b));
        }
        case 2:
        case 7: {
            const a = io.real();
            const b = io.real();
            return B(a).mul(B(b));
        }
        case 3: {
            const a = io.real();
            const b = io.real();
            const c = io.real();
            return B(a).mul(B(b)).add(B(c));
        }
        default: {
            // 5, 6, 8: a + h + c
            const a = io.real();
            const h = io.real();
            const c = io.real();
            return B(a).add(B(h)).add(B(c));
        }
    }
}

// The replay of DrawPair.
function drawPair(io: OracleIO): [BSNumber, BSNumber] {
    let x = drawBSN(io);
    const mode = io.integer();
    let y: BSNumber;
    switch (mode) {
        case 0: y = drawBSN(io); break;
        case 1: y = x.clone(); break;
        case 2: y = x.negated(); break;
        case 3: y = x.add(B(io.real())); break;
        case 4: y = x.mul(B(io.real())); break;
        default: {
            const a = io.real();
            const b = io.real();
            const c = io.real();
            x = B(a).mul(B(c));
            y = B(b).mul(B(c));
            break;
        }
    }
    return [x, y];
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

// -------------------------------------------------------------- QFNumber

function readQF(io: OracleIO, n: number): QFNumber {
    if (n === 1) {
        const x0 = io.real();
        const x1 = io.real();
        return new QFNumber(x0, x1, io.real());
    }
    const x0 = readQF(io, n - 1);
    const x1 = readQF(io, n - 1);
    return new QFNumber(x0, x1, io.real());
}

function outQF(io: OracleIO, q: QFNumber, n: number): void {
    if (n === 1) {
        io.outReal(q.x[0] as number);
        io.outReal(q.x[1] as number);
    } else {
        outQF(io, q.x[0] as QFNumber, n - 1);
        outQF(io, q.x[1] as QFNumber, n - 1);
    }
    io.outReal(q.d);
}

// ------------------------------------------------------------ SWInterval

function readSWI(io: OracleIO): SWInterval {
    const e0 = io.real();
    return new SWInterval(e0, io.real());
}

function outSWI(io: OracleIO, w: SWInterval): void {
    io.outReal(w.get(0));
    io.outReal(w.get(1));
}

// --------------------------------------------------------------- BSplines

// The replay of RecordBasis.
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

describe('oracle: v05-ap', () => {
    const family = new OracleFamily('v05-ap');

    // ----------------------------------------------------------- BSNumber

    family.case('BSNumber.construct.double', (io) => {
        outBSN(io, B(io.real()));
    }, { exact: true });

    family.case('BSNumber.construct.float', (io) => {
        outBSN(io, BSNumber.fromFloat32(io.real()));
    }, { exact: true });

    // int32/uint32 are exact doubles (fromNumber); int64/uint64 go through
    // fromBigInt, the port of the 64-bit constructors.
    family.case('BSNumber.construct.integer', (io) => {
        const kind = io.integer();
        const hi = BigInt(io.integer());
        const lo = BigInt(io.integer());
        const bits = (hi << 32n) | lo;
        let x: BSNumber;
        if (kind === 0) {
            x = BSNumber.fromNumber(Number(BigInt.asIntN(32, bits)));
        } else if (kind === 1) {
            x = BSNumber.fromNumber(Number(BigInt.asUintN(32, bits)));
        } else if (kind === 2) {
            x = BSNumber.fromBigInt(BigInt.asIntN(64, bits));
        } else {
            x = BSNumber.fromBigInt(BigInt.asUintN(64, bits));
        }
        outBSN(io, x);
    }, { exact: true });

    family.case('BSNumber.construct.string', (io) => {
        outBSN(io, BSNumber.fromString(readString(io)));
    }, { exact: true });

    family.case('BSNumber.construct.string.singleChar', (io) => {
        outBSN(io, BSNumber.fromString(readString(io)));
    }, { deviation: '#95 (ConvertToInteger validates one-character strings)' });

    family.case('BSNumber.compare', (io) => {
        const [x, y] = drawPair(io);
        io.outBool(x.equals(y));
        io.outBool(x.notEquals(y));
        io.outBool(x.lessThan(y));
        io.outBool(x.lessThanOrEqual(y));
        io.outBool(x.greaterThan(y));
        io.outBool(x.greaterThanOrEqual(y));
        io.outBool(y.lessThan(x));
        io.outBool(y.lessThanOrEqual(x));
    }, { exact: true });

    // Unary + is the identity (not ported); the compound assignments are
    // the binary operators.
    family.case('BSNumber.arithmetic', (io) => {
        const [x, y] = drawPair(io);
        outBSN(io, x.clone());
        outBSN(io, x.negated());
        outBSN(io, x.add(y));
        outBSN(io, x.sub(y));
        outBSN(io, y.sub(x));
        outBSN(io, x.mul(y));
        let z = x.add(y);
        outBSN(io, z);
        z = z.sub(x);
        outBSN(io, z);
        z = z.mul(y);
        outBSN(io, z);
        const w = x.clone();
        w.negate();
        outBSN(io, w);
    }, { exact: true });

    family.case('BSNumber.conversions.boundary', (io) => {
        const x = drawBSN(io);
        outBSN(io, x);
        outBSN(io, x.negated());
    }, { exact: true });

    family.case('BSNumber.accessors', (io) => {
        const x = drawBSN(io);
        const e = io.integer();
        const be = io.integer();
        const y = x.clone();
        y.setExponent(e);
        outBSN(io, y);
        const z = x.clone();
        z.setBiasedExponent(be);
        outBSN(io, z);
        const w = x.clone();
        w.setSign(-x.getSign() | 0);
        outBSN(io, w);
    }, { exact: true });

    family.case('BSNumber.convert', (io) => {
        const x = drawBSN(io);
        const precision = io.integer();
        const mode = roundingModes[io.integer()];
        outBSN(io, convertBSNumber(x, precision, mode));
    }, { exact: true });

    family.case('BSNumber.std.exact', (io) => {
        const x = drawBSN(io);
        const y = drawBSN(io);
        const z = drawBSN(io);
        const w = drawBSN(io);
        const k = io.integer();
        outBSN(io, BSNumber.fabs(x));
        const fr = BSNumber.frexp(x);
        outBSN(io, fr.result);
        io.outInt(fr.exponent);
        outBSN(io, BSNumber.ldexp(x, k));
        outBSN(io, BSNumber.floor(x));
        outBSN(io, BSNumber.ceil(x));
        outBSN(io, BSNumber.sqrt(x));
        outBSN(io, BSNumber.fmod(x, y));
        outBSN(io, BSNumber.remainder(x, y));
        outBSN(io, BSNumber.clamp(x, y, z));
        outBSN(io, BSNumber.invsqrt(x));
        io.outInt(BSNumber.isign(x));
        outBSN(io, BSNumber.saturate(x));
        outBSN(io, BSNumber.sign(x));
        outBSN(io, BSNumber.sqr(x));
        outBSN(io, BSNumber.fma(x, y, z));
        outBSN(io, BSNumber.robustSOP(x, y, z, w));
        outBSN(io, BSNumber.robustDOP(x, y, z, w));
    }, { exact: true });

    family.case('BSNumber.std.remainder', (io) => {
        const x = B(io.real());
        const y = B(io.real());
        outBSN(io, BSNumber.remainder(x, y));
        outBSN(io, BSNumber.remainder(y, x));
        outBSN(io, BSNumber.fmod(x, y));
    }, { exact: true });

    // std::acos ... std::tanh and gte::atandivpi, atan2divpi, cospi, exp10,
    // sinpi: C math library on the MSVC side, V8's Math on this side.
    family.case('BSNumber.std.libm', (io) => {
        const x = B(io.real());
        const y = B(io.real());
        const v = (r: BSNumber): void => { io.outReal(r.toNumber()); };
        v(BSNumber.acos(x));
        v(BSNumber.acosh(x));
        v(BSNumber.asin(x));
        v(BSNumber.asinh(x));
        v(BSNumber.atan(x));
        v(BSNumber.atanh(x));
        v(BSNumber.atan2(y, x));
        v(BSNumber.cos(x));
        v(BSNumber.cosh(x));
        v(BSNumber.exp(x));
        v(BSNumber.exp2(x));
        v(BSNumber.log(x));
        v(BSNumber.log2(x));
        v(BSNumber.log10(x));
        v(BSNumber.pow(x, y));
        v(BSNumber.sin(x));
        v(BSNumber.sinh(x));
        v(BSNumber.tan(x));
        v(BSNumber.tanh(x));
        v(BSNumber.atandivpi(x));
        v(BSNumber.atan2divpi(y, x));
        v(BSNumber.cospi(x));
        v(BSNumber.exp10(x));
        v(BSNumber.sinpi(x));
    });

    // -------------------------------------------------------- BSPrecision

    const outParameters = (io: OracleIO, p: BSPrecisionParameters): void => {
        io.outInt(p.minExponent);
        io.outInt(p.maxExponent);
        io.outInt(p.maxBits);
        io.outInt(p.maxWords);
    };

    const outBSP = (io: OracleIO, p: BSPrecision): void => {
        outParameters(io, p.bsn);
        outParameters(io, p.bsr);
    };

    const drawBSP = (io: OracleIO): BSPrecision => {
        const kind = io.integer();
        if (kind === 0) {
            return new BSPrecision(io.integer() as BSPrecisionType);
        }
        if (kind === 1 || kind === 2) {
            const mn = io.integer();
            const mx = io.integer();
            return new BSPrecision(mn, mx, io.integer());
        }
        switch (io.integer()) {
            case 0: return new BSPrecision(BSPrecisionType.IS_DOUBLE);
            case 1: return new BSPrecision(BSPrecisionType.IS_FLOAT);
            case 2: return new BSPrecision(-17, 20, 24);
            default: return new BSPrecision(-5, 20, 12);
        }
    };

    const drawBSPExpression = (io: OracleIO, depth: number): BSPrecision => {
        const op = (depth === 0 ? 0 : io.integer());
        if (op === 0) { return drawBSP(io); }
        const a = drawBSPExpression(io, depth - 1);
        const b = drawBSPExpression(io, depth - 1);
        switch (op) {
            case 1: return a.add(b);
            case 2: return a.sub(b);
            case 3: return a.mul(b);
            case 4: return a.div(b);
            case 5: return a.equal(b);
            default: return a.lessThan(b);
        }
    };

    family.case('BSPrecision.operators', (io) => {
        const a = drawBSP(io);
        const b = drawBSP(io);
        outBSP(io, a);
        outBSP(io, b);
        outBSP(io, a.add(b));
        outBSP(io, b.add(a));
        outBSP(io, a.sub(b));
        outBSP(io, a.mul(b));
        outBSP(io, a.div(b));
        outBSP(io, a.equal(b));
        outBSP(io, a.notEqual(b));
        outBSP(io, a.lessThan(b));
        outBSP(io, a.lessThanEqual(b));
        outBSP(io, a.greaterThan(b));
        outBSP(io, a.greaterThanEqual(b));
        outBSP(io, new BSPrecision());
        io.outInt(new BSPrecisionParameters().getMaxWords());
    }, { exact: true });

    family.case('BSPrecision.expression', (io) => {
        outBSP(io, drawBSPExpression(io, 3));
    }, { exact: true });

    // ----------------------------------------------------------- QFNumber

    // s + q and s * q are q.add(s) and q.mul(s) (IEEE + and * commute); the
    // compound assignments are the binary operators.
    const qfArithmetic = (io: OracleIO, n: number): void => {
        const q0 = readQF(io, n);
        const q1 = readQF(io, n);
        const s = io.real();
        outQF(io, q0.clone(), n);
        outQF(io, q0.negate(), n);
        outQF(io, q0.add(q1), n);
        outQF(io, q0.add(s), n);
        outQF(io, q0.add(s), n);
        outQF(io, q0.sub(q1), n);
        outQF(io, q0.sub(s), n);
        outQF(io, QFNumber.scalarSub(s, q0), n);
        outQF(io, q0.mul(q1), n);
        outQF(io, q0.mul(s), n);
        outQF(io, q0.mul(s), n);
        outQF(io, q0.div(q1), n);
        outQF(io, q0.div(s), n);
        outQF(io, QFNumber.scalarDiv(s, q0), n);
        let c = q0.add(q1);
        outQF(io, c, n);
        c = c.add(s);
        outQF(io, c, n);
        c = c.sub(q1);
        outQF(io, c, n);
        c = c.sub(s);
        outQF(io, c, n);
        c = c.mul(q1);
        outQF(io, c, n);
        c = c.mul(s);
        outQF(io, c, n);
        c = c.div(q1);
        outQF(io, c, n);
        c = c.div(s);
        outQF(io, c, n);
    };

    const qfCompare = (io: OracleIO, n: number): void => {
        const q0 = readQF(io, n);
        const q1 = readQF(io, n);
        io.real();  // the scalar drawn with the pair, unused here
        io.outBool(q0.equals(q1));
        io.outBool(q0.notEquals(q1));
        io.outBool(q0.lessThan(q1));
        io.outBool(q0.lessThanEqual(q1));
        io.outBool(q0.greaterThan(q1));
        io.outBool(q0.greaterThanEqual(q1));
        io.outBool(q1.lessThan(q0));
        io.outBool(q1.equals(q0));
    };

    family.case('QFNumber.arithmetic.n1', (io) => { qfArithmetic(io, 1); }, { exact: true });
    family.case('QFNumber.arithmetic.n2', (io) => { qfArithmetic(io, 2); }, { exact: true });
    family.case('QFNumber.arithmetic.n3', (io) => { qfArithmetic(io, 3); }, { exact: true });
    family.case('QFNumber.compare.n1', (io) => { qfCompare(io, 1); }, { exact: true });
    family.case('QFNumber.compare.n2', (io) => { qfCompare(io, 2); }, { exact: true });
    family.case('QFNumber.compare.n3', (io) => { qfCompare(io, 3); }, { exact: true });

    family.case('QFNumber.construct', (io) => {
        const d0 = io.real();
        const d1 = io.real();
        const a = io.real();
        const b = io.real();
        outQF(io, new QFNumber(), 1);
        outQF(io, QFNumber.fromD(d0), 1);
        outQF(io, new QFNumber(a, b, d0), 1);
        outQF(io, new QFNumber(b, a, d0), 1);
        const inner = new QFNumber(a, b, d0);
        outQF(io, new QFNumber(inner, inner.negate(), d1), 2);
        outQF(io, new QFNumber(inner.negate(), inner, d1), 2);
    }, { exact: true });

    // --------------------------------------------------------- SWInterval

    family.case('SWInterval.leaf', (io) => {
        const u = io.real();
        const v = io.real();
        outSWI(io, SWInterval.add(u, v));
        outSWI(io, SWInterval.sub(u, v));
        outSWI(io, SWInterval.mul(u, v));
        outSWI(io, SWInterval.div(u, v));
        outSWI(io, SWInterval.div(v, u));
    }, { exact: true });

    family.case('SWInterval.internal', (io) => {
        const u0 = io.real();
        const u1 = io.real();
        const v0 = io.real();
        const v1 = io.real();
        outSWI(io, SWInterval.add(u0, u1, v0, v1));
        outSWI(io, SWInterval.sub(u0, u1, v0, v1));
        outSWI(io, SWInterval.mul(u0, u1, v0, v1));
        outSWI(io, SWInterval.mul2(u0, u1, v0, v1));
        outSWI(io, SWInterval.div(u0, u1, v0, v1));
        outSWI(io, SWInterval.reciprocal(v0, v1));
        outSWI(io, SWInterval.reciprocalDown(v0));
        outSWI(io, SWInterval.reciprocalUp(v1));
        outSWI(io, SWInterval.reals());
    }, { exact: true });

    // s + u and s * u are u.add(s) and u.mul(s): upstream's Add(s, s, u0,
    // u1) and Mul(s, s, ...) form the same commuted products and sums.
    family.case('SWInterval.operators', (io) => {
        const u = readSWI(io);
        const v = readSWI(io);
        const s = io.real();
        outSWI(io, u.clone());
        outSWI(io, u.negate());
        outSWI(io, u.add(v));
        outSWI(io, u.add(s));
        outSWI(io, u.add(s));
        outSWI(io, u.sub(v));
        outSWI(io, u.sub(s));
        outSWI(io, SWInterval.scalarSub(s, u));
        outSWI(io, u.mul(v));
        outSWI(io, v.mul(u));
        outSWI(io, u.mul(s));
        outSWI(io, u.mul(s));
        outSWI(io, u.div(v));
        outSWI(io, u.div(s));
        outSWI(io, SWInterval.scalarDiv(s, u));
        let c = u.add(v);
        outSWI(io, c);
        c = c.add(s);
        outSWI(io, c);
        c = c.sub(v);
        outSWI(io, c);
        c = c.sub(s);
        outSWI(io, c);
        c = c.mul(v);
        outSWI(io, c);
        c = c.mul(s);
        outSWI(io, c);
        c = c.div(v);
        outSWI(io, c);
        c = c.div(s);
        outSWI(io, c);
        const ends = u.getEndpoints();
        io.outReal(ends[0]);
        io.outReal(ends[1]);
        outSWI(io, new SWInterval());
        outSWI(io, new SWInterval(s));
        outSWI(io, SWInterval.fromEndpoints([v.get(0), v.get(1)]));
    }, { exact: true });

    family.case('SWInterval.nextafter', (io) => {
        outSWI(io, SWInterval.mul(io.real(), 1));
    }, { exact: true });

    family.case('SWInterval.overflow', (io) => {
        const a = io.real();
        const b = io.real();
        const t = io.real();
        outSWI(io, SWInterval.add(a, b));
        outSWI(io, SWInterval.mul(a, b));
        outSWI(io, SWInterval.div(a, t));
        outSWI(io, SWInterval.div(8.881784197001252e-16, -Number.MIN_VALUE));
        outSWI(io, new SWInterval(1, 2).div(new SWInterval(0, 4)));
        outSWI(io, new SWInterval(a, a).mul(new SWInterval(b, b)));
    }, { exact: true });

    // ---------------------------------------------------- BSplineCurveFit

    family.case('BSplineCurveFit.fit', (io) => {
        const dimension = io.integer();
        const degree = io.integer();
        const numControls = io.integer();
        const numSamples = io.integer();
        const P = io.reals(numSamples * dimension);
        const fit = new BSplineCurveFit(dimension, numSamples, P, degree, numControls);
        io.outInt(fit.getDimension());
        io.outInt(fit.getNumSamples());
        io.outInt(fit.getDegree());
        io.outInt(fit.getNumControls());
        io.outReal(fit.getSampleData()[0]);
        io.outReals(fit.getControlData().slice(0, dimension * numControls));
        io.outReal(fit.getBasis().getMinDomain());
        io.outReal(fit.getBasis().getMaxDomain());
        io.outInt(fit.getBasis().getNumKnots());
        for (let k = 0; k < 6; ++k) {
            const t = io.real();
            io.outReals(fit.evaluate(t, io.integer()));
        }
        io.outReals(fit.getPosition(io.real()));
    }, { exact: true });

    family.case('BSplineCurveFit.invalid', (io) => {
        io.integer();
        const dimension = io.integer();
        const degree = io.integer();
        const numControls = io.integer();
        const numSamples = io.integer();
        const P = io.reals(Math.max(0, numSamples * dimension));
        const fit = new BSplineCurveFit(dimension, numSamples, P, degree, numControls);
        io.outReal(fit.evaluate(0.5, 4)[0]);
    }, { exact: true });

    // -------------------------------------------------- BSplineSurfaceFit

    family.case('BSplineSurfaceFit.fit', (io) => {
        const degree0 = io.integer();
        const numControls0 = io.integer();
        const numSamples0 = io.integer();
        const degree1 = io.integer();
        const numControls1 = io.integer();
        const numSamples1 = io.integer();
        const samples = points(io, numSamples0 * numSamples1, 3);
        const fit = new BSplineSurfaceFit(degree0, numControls0, numSamples0,
            degree1, numControls1, numSamples1, samples);
        for (let d = 0; d < 2; ++d) {
            io.outInt(fit.getNumSamples(d));
            io.outInt(fit.getDegree(d));
            io.outInt(fit.getNumControls(d));
            io.outReal(fit.getBasis(d).getMinDomain());
            io.outReal(fit.getBasis(d).getMaxDomain());
        }
        io.outVec(fit.getSampleData()[0]);
        const controls = fit.getControlData();
        for (let i = 0; i < numControls0 * numControls1; ++i) { io.outVec(controls[i]); }
        for (let k = 0; k < 5; ++k) {
            const u = io.real();
            io.outVec(fit.getPosition(u, io.real()));
        }
    }, { exact: true });

    family.case('BSplineSurfaceFit.invalid', (io) => {
        io.integer();
        const degree0 = io.integer();
        const numControls0 = io.integer();
        const numSamples0 = io.integer();
        const degree1 = io.integer();
        const numControls1 = io.integer();
        const numSamples1 = io.integer();
        const samples = points(io, numSamples0 * numSamples1, 3);
        const fit = new BSplineSurfaceFit(degree0, numControls0, numSamples0,
            degree1, numControls1, numSamples1, samples);
        io.outInt(fit.getNumControls(0));
    }, { exact: true });

    // ----------------------------------------------------- BSplineSurface

    family.case('BSplineSurface.evaluate', (io) => {
        const N = io.integer();
        const in0 = basisInput(io);
        const in1 = basisInput(io);
        const cmode = io.integer();
        const n0 = in0.numControls, n1 = in1.numControls;
        const controls = (cmode <= 4 ? points(io, n0 * n1, N) : undefined);
        const surface = new BSplineSurface(N, [in0, in1], controls);
        for (let k = 0; k < 3; ++k) {
            const i0 = io.integer();
            const i1 = io.integer();
            surface.setControl(i0, i1, io.vec(N));
        }
        io.outBool(surface.isConstructed());
        io.outReal(surface.getUMin());
        io.outReal(surface.getUMax());
        io.outReal(surface.getVMin());
        io.outReal(surface.getVMax());
        io.outBool(surface.isRectangular());
        io.outInt(surface.getNumControls(0));
        io.outInt(surface.getNumControls(1));
        io.outReal(surface.getBasisFunction(0).getMaxDomain());
        io.outReal(surface.getBasisFunction(1).getMinDomain());
        for (let k = 0; k < 3; ++k) {
            const i0 = io.integer();
            io.outVec(surface.getControl(i0, io.integer()));
        }
        io.outVec(surface.getControls()[n0 * n1 - 1]);
        for (let k = 0; k < 5; ++k) {
            const u = io.real();
            const v = io.real();
            const order = io.integer();
            const jet = surface.createJet();
            surface.evaluate(u, v, order, jet);
            for (let i = 0; i < 6; ++i) { io.outVec(jet[i]); }
        }
        const u = io.real();
        const v = io.real();
        io.outVec(surface.getPosition(u, v));
        io.outVec(surface.getUTangent(u, v));
        io.outVec(surface.getVTangent(u, v));
    }, { exact: true });

    family.case('BSplineSurface.evaluate.invalidOrder', (io) => {
        const in0 = basisInput(io);
        const in1 = basisInput(io);
        const controls = points(io, in0.numControls * in1.numControls, 2);
        const surface = new BSplineSurface(2, [in0, in1], controls);
        const order = io.integer();
        const jet = surface.createJet();
        surface.evaluate(0.5 * (surface.getUMin() + surface.getUMax()), surface.getVMin(),
            order, jet);
        io.outVec(jet[0]);
    }, { exact: true });

    // ------------------------------------------------------ BSplineVolume

    family.case('BSplineVolume.evaluate', (io) => {
        const N = io.integer();
        const inputs = [basisInput(io), basisInput(io), basisInput(io)];
        const cmode = io.integer();
        const n0 = inputs[0].numControls, n1 = inputs[1].numControls, n2 = inputs[2].numControls;
        const controls = (cmode <= 4 ? points(io, n0 * n1 * n2, N) : undefined);
        const volume = new BSplineVolume(N, inputs, controls);
        for (let k = 0; k < 3; ++k) {
            const i0 = io.integer();
            const i1 = io.integer();
            const i2 = io.integer();
            volume.setControl(i0, i1, i2, io.vec(N));
        }
        io.outBool(volume.isConstructed());
        for (let d = 0; d < 3; ++d) {
            io.outReal(volume.getMinDomain(d));
            io.outReal(volume.getMaxDomain(d));
            io.outInt(volume.getNumControls(d));
            io.outInt(volume.getBasisFunction(d).getNumKnots());
        }
        for (let k = 0; k < 3; ++k) {
            const i0 = io.integer();
            const i1 = io.integer();
            io.outVec(volume.getControl(i0, i1, io.integer()));
        }
        io.outVec(volume.getControls()[n0 * n1 * n2 - 1]);
        for (let k = 0; k < 4; ++k) {
            const t = io.reals(3);
            const order = io.integer();
            const jet = volume.createJet();
            volume.evaluate(t[0], t[1], t[2], order, jet);
            for (let i = 0; i < 10; ++i) { io.outVec(jet[i]); }
        }
    }, { exact: true });

    family.case('BSplineVolume.evaluate.invalidOrder', (io) => {
        const inputs = [basisInput(io), basisInput(io), basisInput(io)];
        const controls = points(io,
            inputs[0].numControls * inputs[1].numControls * inputs[2].numControls, 1);
        const volume = new BSplineVolume(1, inputs, controls);
        const order = io.integer();
        const jet = volume.createJet();
        volume.evaluate(volume.getMinDomain(0), volume.getMaxDomain(1), volume.getMinDomain(2),
            order, jet);
        io.outVec(jet[0]);
    }, { exact: true });

    family.finish();
});
