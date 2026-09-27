// Replays oracle/cpp/cases/v16-core.cpp. Keep the two files in the same order.
import { describe } from 'vitest';
import { CurveExtractorEdge, CurveExtractorVertex } from '../../src/CurveExtractor.js';
import { CurveExtractorSquares } from '../../src/CurveExtractorSquares.js';
import { CurveExtractorTriangles } from '../../src/CurveExtractorTriangles.js';
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

describe('oracle: v16-core', () => {
    const family = new OracleFamily('v16-core');

    family.case('CurveExtractorSquares.extract', (io) => runExtractor(io, false),
        { exact: true });

    family.case('CurveExtractorTriangles.extract', (io) => runExtractor(io, true),
        { exact: true });

    family.finish();
});
