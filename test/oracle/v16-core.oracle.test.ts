// Replays oracle/cpp/cases/v16-core.cpp. Keep the two files in the same order.
import { describe } from 'vitest';
import { CurveExtractorEdge, CurveExtractorVertex } from '../../src/CurveExtractor.js';
import { CurveExtractorSquares } from '../../src/CurveExtractorSquares.js';
import { CurveExtractorTriangles } from '../../src/CurveExtractorTriangles.js';
import { IEEEBinary16 } from '../../src/IEEEBinary16.js';
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

    family.finish();
});
