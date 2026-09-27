// Replays oracle/cpp/cases/v26-imaging.cpp. Keep the two files in the same
// order. Every case is integer bookkeeping or IEEE-exact arithmetic (sqrt
// and sqrtf included) and is compared bit for bit. The independent checks
// at the end of each section run on the port's outputs of the replayed
// records.
import { mkdirSync, writeFileSync } from 'node:fs';
import { dirname, join } from 'node:path';
import { fileURLToPath } from 'node:url';
import { describe, expect, it } from 'vitest';
import { Image2 } from '../../src/Image2.js';
import { Image3 } from '../../src/Image3.js';
import { MarchingCubes } from '../../src/MarchingCubes.js';
import { ImageUtility2 } from '../../src/ImageUtility2.js';
import { ImageUtility3 } from '../../src/ImageUtility3.js';
import {
    SurfaceExtractor,
    SurfaceExtractorTriangle,
    SurfaceExtractorVertex
} from '../../src/SurfaceExtractor.js';
import { SurfaceExtractorCubes } from '../../src/SurfaceExtractorCubes.js';
import { SurfaceExtractorMC } from '../../src/SurfaceExtractorMC.js';
import { SurfaceExtractorTetrahedra } from '../../src/SurfaceExtractorTetrahedra.js';
import { Vector } from '../../src/Vector.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- digests and packing (mirroring the C++ helpers) ----

class Digest {
    h = 2166136261;

    word(w: number): void {
        for (let i = 0; i < 4; ++i) {
            this.h = (this.h ^ ((w >>> (8 * i)) & 0xFF)) >>> 0;
            this.h = Math.imul(this.h, 16777619) >>> 0;
        }
    }

    // A signed 64-bit integer (a safe integer here) as two words.
    integer(v: number): void {
        const u = BigInt.asUintN(64, BigInt(v));
        this.word(Number(u & 0xFFFFFFFFn));
        this.word(Number(u >> 32n));
    }

    real(x: number): void {
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(0, x, true);
        this.word(view.getUint32(0, true));
        this.word(view.getUint32(4, true));
    }
}

function pack(values: readonly number[], bits: number): number[] {
    const per = Math.floor(48 / bits);
    const scale = 2 ** bits;
    const out: number[] = [];
    for (let i = 0; i < values.length; i += per) {
        let chunk = 0;
        for (let j = Math.min(values.length, i + per); j > i; --j) {
            chunk = chunk * scale + values[j - 1];
        }
        out.push(chunk);
    }
    return out;
}

function unpack(io: OracleIO, n: number, bits: number): number[] {
    const per = Math.floor(48 / bits);
    const scale = 2 ** bits;
    const out: number[] = [];
    while (out.length < n) {
        let chunk = io.integer();
        for (let j = 0; j < per && out.length < n; ++j) {
            const v = chunk % scale;
            out.push(v);
            chunk = (chunk - v) / scale;
        }
    }
    return out;
}

function outPacked(io: OracleIO, values: readonly number[], bits: number): void {
    for (const c of pack(values, bits)) { io.outReal(c); }
}

function outPoint(io: OracleIO, x: number, y: number, z = 0): void {
    io.outInt((x + 32768) + 65536 * (y + 32768) + 4294967296 * (z + 32768));
}

class PointRecorder {
    readonly points: [number, number, number][] = [];
    readonly cb2 = (x: number, y: number): void => { this.points.push([x, y, 0]); };
    readonly cb3 = (x: number, y: number, z: number): void => { this.points.push([x, y, z]); };

    emit(io: OracleIO, full: boolean): void {
        io.outInt(this.points.length);
        if (full) {
            for (const p of this.points) { outPoint(io, p[0], p[1], p[2]); }
        } else {
            const d = new Digest();
            for (const p of this.points) {
                d.integer(p[0]);
                d.integer(p[1]);
                d.integer(p[2]);
            }
            io.outInt(d.h);
        }
    }
}

// ---- images ----

function makeImage2(d0: number, d1: number, pix: readonly number[]): Image2<number> {
    const image = new Image2<number>(d0, d1);
    const p = image.getPixels();
    for (let i = 0; i < pix.length; ++i) { p[i] = pix[i]; }
    return image;
}

function readImage2(io: OracleIO): Image2<number> {
    const d0 = io.integer();
    const d1 = io.integer();
    return makeImage2(d0, d1, unpack(io, d0 * d1, 2));
}

function outImage2(io: OracleIO, image: Image2<number>, bits: number): void {
    outPacked(io, image.getPixels(), bits);
}

function makeImage3(d0: number, d1: number, d2: number, pix: readonly number[]): Image3<number> {
    const image = new Image3<number>(d0, d1, d2);
    const p = image.getPixels();
    for (let i = 0; i < pix.length; ++i) { p[i] = pix[i]; }
    return image;
}

function readImage3(io: OracleIO): Image3<number> {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    return makeImage3(d0, d1, d2, unpack(io, d0 * d1 * d2, 2));
}

function outImage3(io: OracleIO, image: Image3<number>, bits: number): void {
    outPacked(io, image.getPixels(), bits);
}

// ---- ImageUtility2: components and morphology ----

// Records of the replayed cases, checked independently afterwards.
interface ComponentsRecord {
    input: number[];
    dims: number[];
    // Relative neighbour offsets as coordinate tuples, or null for a
    // caller-specified (possibly asymmetric) list.
    offsets: number[][] | null;
    components: number[][];
    labels: number[];
}
const componentsChecked: ComponentsRecord[] = [];

const NBR4 = [[-1, 0], [1, 0], [0, -1], [0, 1]];
const NBR8 = [...NBR4, [-1, -1], [1, -1], [-1, 1], [1, 1]];

function components2Case(io: OracleIO): void {
    const kind = io.integer();
    const image = readImage2(io);
    const d0 = image.getDimension(0);
    const input = image.getPixels().slice();
    let components: number[][];
    if (kind === 0) {
        components = ImageUtility2.getComponents4(image);
    } else if (kind === 1) {
        components = ImageUtility2.getComponents8(image);
    } else {
        const n = io.integer();
        const nbrs: number[] = [];
        for (let k = 0; k < n; ++k) {
            const dx = io.integer();
            const dy = io.integer();
            nbrs.push(dx + d0 * dy);
        }
        components = ImageUtility2.getComponents(image, nbrs);
    }
    io.outInt(components.length);
    for (let k = 1; k < components.length; ++k) {
        io.outInt(components[k].length);
        outPacked(io, components[k], 8);
    }
    outImage2(io, image, 8);
    componentsChecked.push({
        input, dims: [d0, image.getDimension(1)],
        offsets: kind === 0 ? NBR4 : kind === 1 ? NBR8 : null,
        components, labels: image.getPixels().slice()
    });
}

type Offsets = readonly (readonly number[])[];

function readOffsets(io: OracleIO, dim: number): number[][] {
    const n = io.integer();
    const nbrs: number[][] = [];
    for (let k = 0; k < n; ++k) {
        const d: number[] = [];
        for (let j = 0; j < dim; ++j) { d.push(io.integer()); }
        nbrs.push(d);
    }
    return nbrs;
}

interface MorphRecord {
    op: number;
    zeroExterior: boolean;
    offsets: Offsets;
    dims: number[];
    input: number[];
    output: number[];
}
const morphChecked: MorphRecord[] = [];

function morph2Case(io: OracleIO, op: number): void {
    const kind = io.integer();
    const image = readImage2(io);
    const zeroExterior = (op !== 0 ? io.boolean() : false);
    const nbrs = (kind === 2 ? readOffsets(io, 2) : []);
    const U = ImageUtility2;
    let out: Image2<number>;
    switch (op) {
        case 0:
            out = kind === 0 ? U.dilate4(image) : kind === 1 ? U.dilate8(image) : U.dilate(image, nbrs);
            break;
        case 1:
            out = kind === 0 ? U.erode4(image, zeroExterior) : kind === 1
                ? U.erode8(image, zeroExterior) : U.erode(image, zeroExterior, nbrs);
            break;
        case 2:
            out = kind === 0 ? U.open4(image, zeroExterior) : kind === 1
                ? U.open8(image, zeroExterior) : U.open(image, zeroExterior, nbrs);
            break;
        default:
            out = kind === 0 ? U.close4(image, zeroExterior) : kind === 1
                ? U.close8(image, zeroExterior) : U.close(image, zeroExterior, nbrs);
            break;
    }
    io.outInt(out.getDimension(0));
    io.outInt(out.getDimension(1));
    outImage2(io, out, 2);
    morphChecked.push({
        op, zeroExterior, offsets: kind === 0 ? NBR4 : kind === 1 ? NBR8 : nbrs,
        dims: [image.getDimension(0), image.getDimension(1), 1],
        input: image.getPixels().slice(), output: out.getPixels().slice()
    });
}

interface BoundaryRecord {
    before: number[];
    dims: number[];
    start: number;
    success: boolean;
    boundary: number[];
}
const boundaryChecked: BoundaryRecord[] = [];

function extractBoundaryCase(io: OracleIO): void {
    const image = readImage2(io);
    const calls = io.integer();
    for (let c = 0; c < calls; ++c) {
        const x = io.integer();
        const y = io.integer();
        const before = image.getPixels().slice();
        const r = ImageUtility2.extractBoundary(x, y, image);
        io.outBool(r.success);
        io.outInt(r.boundary.length);
        outPacked(io, r.boundary, 8);
        boundaryChecked.push({
            before, dims: [image.getDimension(0), image.getDimension(1)],
            start: image.getIndex(x, y), success: r.success, boundary: r.boundary
        });
    }
    outImage2(io, image, 2);
}

// ---- ImageUtility2: fills and distance transforms ----

interface FillRecord {
    dims: number[];
    before: number[];
    seed: number[];
    fore: number;
    back: number;
    after: number[];
}
const fillChecked: FillRecord[] = [];

function floodFill2IntCase(io: OracleIO): void {
    const image = readImage2(io);
    const x = io.integer();
    const y = io.integer();
    const back = io.integer();
    const fore = (back + io.integer()) % 4;
    const before = image.getPixels().slice();
    ImageUtility2.floodFill4(image, x, y, fore, back);
    outImage2(io, image, 2);
    fillChecked.push({
        dims: [image.getDimension(0), image.getDimension(1), 1], before, seed: [x, y, 0],
        fore, back, after: image.getPixels().slice()
    });
}

function floodFill2DoubleCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const image = new Image2<number>(d0, d1);
    for (let i = 0; i < d0 * d1; ++i) { image.set(i, io.real()); }
    const x = io.integer();
    const y = io.integer();
    const back = io.real();
    const fore = io.real();
    const before = image.getPixels().slice();
    ImageUtility2.floodFill4(image, x, y, fore, back);
    for (const v of image.getPixels()) { io.outReal(v); }
    fillChecked.push({
        dims: [d0, d1, 1], before, seed: [x, y, 0], fore, back, after: image.getPixels().slice()
    });
}

interface DistanceRecord {
    dims: number[];
    input: number[];
    output: number[];
    l2: boolean;
    max: number;
    arg: number[];
}
const distanceChecked: DistanceRecord[] = [];

function l1DistanceCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    io.boolean();  // whether the generator zeroed the border
    const image = makeImage2(d0, d1, unpack(io, d0 * d1, 2));
    const input = image.getPixels().slice();
    const r = ImageUtility2.getL1Distance(image);
    io.outInt(r.maxDistance);
    io.outInt(r.xMax);
    io.outInt(r.yMax);
    outImage2(io, image, 4);
    distanceChecked.push({
        dims: [image.getDimension(0), image.getDimension(1)], input,
        output: image.getPixels().slice(), l2: false, max: r.maxDistance, arg: [r.xMax, r.yMax]
    });
}

function l2DistanceCase(io: OracleIO): void {
    const image = readImage2(io);
    const r = ImageUtility2.getL2Distance(image);
    io.outReal(r.maxDistance);
    io.outInt(r.xMax);
    io.outInt(r.yMax);
    for (const v of r.transform.getPixels()) { io.outReal(v); }
    distanceChecked.push({
        dims: [image.getDimension(0), image.getDimension(1)], input: image.getPixels().slice(),
        output: r.transform.getPixels().slice(), l2: true, max: r.maxDistance, arg: [r.xMax, r.yMax]
    });
}

function l2DistanceLargeCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const image = new Image2<number>(d0, d1);
    image.getPixels().fill(1);
    if (io.boolean()) {
        for (let y = 0; y < d1; ++y) {
            for (let x = 0; x < d0; ++x) {
                if (x === 0 || y === 0 || x === d0 - 1 || y === d1 - 1) { image.set(x, y, 0); }
            }
        }
    }
    const k = io.integer();
    for (let j = 0; j < k; ++j) { image.set(io.integer(), 0); }
    const r = ImageUtility2.getL2Distance(image);
    io.outReal(r.maxDistance);
    io.outInt(r.xMax);
    io.outInt(r.yMax);
    const d = new Digest();
    for (const v of r.transform.getPixels()) { d.real(v); }
    io.outInt(d.h);
    distanceChecked.push({
        dims: [d0, d1], input: image.getPixels().slice(), output: r.transform.getPixels().slice(),
        l2: true, max: r.maxDistance, arg: [r.xMax, r.yMax]
    });
}

interface SkeletonRecord { dims: number[]; input: number[]; output: number[] }
const skeletonChecked: SkeletonRecord[] = [];

function skeletonCase(io: OracleIO): void {
    const image = readImage2(io);
    const input = image.getPixels().slice();
    ImageUtility2.getSkeleton(image);
    outImage2(io, image, 1);
    skeletonChecked.push({
        dims: [image.getDimension(0), image.getDimension(1)], input, output: image.getPixels().slice()
    });
}

// ---- ImageUtility2: Draw* ----

interface DrawRecord { kind: string; args: number[]; points: [number, number, number][] }
const drawChecked: DrawRecord[] = [];

function drawThickPixelCase(io: OracleIO): void {
    const x = io.integer();
    const y = io.integer();
    const thick = io.integer();
    const rec = new PointRecorder();
    ImageUtility2.drawThickPixel(x, y, thick, rec.cb2);
    rec.emit(io, true);
    drawChecked.push({ kind: 'thick', args: [x, y, thick], points: rec.points });
}

function drawLine2Case(io: OracleIO): void {
    const mode = io.integer();
    const x0 = io.integer();
    const y0 = io.integer();
    let x1 = io.integer();
    let y1 = io.integer();
    const len = io.integer();
    const s = io.integer();
    switch (mode) {
        case 1: if (s >= 0) { y1 = y0; } else { x1 = x0; } break;
        case 2: x1 = x0 + len; y1 = y0 + (s >= 0 ? len : -len); break;
        case 3: x1 = x0; y1 = y0; break;
        case 4: x1 = x0 + len; y1 = y0 + (len >= 0 ? len + s : len - s); break;
        default: break;
    }
    const rec = new PointRecorder();
    ImageUtility2.drawLine(x0, y0, x1, y1, rec.cb2);
    rec.emit(io, mode !== 5);
    drawChecked.push({ kind: 'line', args: [x0, y0, 0, x1, y1, 0], points: rec.points });
}

function drawCircleCase(io: OracleIO): void {
    const mode = io.integer();
    const xc = io.integer();
    const yc = io.integer();
    const solid = io.boolean();
    const radius = io.integer();
    const rec = new PointRecorder();
    ImageUtility2.drawCircle(xc, yc, radius, solid, rec.cb2);
    rec.emit(io, mode !== 3);
    drawChecked.push({ kind: solid ? 'disk' : 'circle', args: [xc, yc, radius], points: rec.points });
}

function drawRectangleCase(io: OracleIO): void {
    const xMin = io.integer();
    const yMin = io.integer();
    const xMax = io.integer();
    const yMax = io.integer();
    const solid = io.boolean();
    const rec = new PointRecorder();
    ImageUtility2.drawRectangle(xMin, yMin, xMax, yMax, solid, rec.cb2);
    rec.emit(io, true);
    drawChecked.push({ kind: solid ? 'box' : 'rect', args: [xMin, yMin, xMax, yMax], points: rec.points });
}

function drawEllipseCase(io: OracleIO): void {
    const mode = io.integer();
    const xc = io.integer();
    const yc = io.integer();
    const a = io.integer();
    let b = io.integer();
    if (mode === 1) { b = a + (b & 1); if (b === 0) { b = 1; } }
    const rec = new PointRecorder();
    ImageUtility2.drawEllipse(xc, yc, a, b, rec.cb2);
    rec.emit(io, mode !== 3);
    drawChecked.push({ kind: 'ellipse', args: [xc, yc, a, b], points: rec.points });
}

const zeroEllipses: { center: number[], points: [number, number, number][] }[] = [];

function drawEllipseZeroCase(io: OracleIO): void {
    const xc = io.integer();
    const yc = io.integer();
    const rec = new PointRecorder();
    ImageUtility2.drawEllipse(xc, yc, 0, 0, rec.cb2);
    zeroEllipses.push({ center: [xc, yc], points: rec.points });
    io.outInt(rec.points.length);
}

function drawFloodFillCase(io: OracleIO): void {
    const image = readImage2(io);
    const d0 = image.getDimension(0);
    const d1 = image.getDimension(1);
    const x = io.integer();
    const y = io.integer();
    const back = io.integer();
    const fore = (back + io.integer()) % 4;
    const before = image.getPixels().slice();
    const sets = new PointRecorder();
    let numGets = 0;
    ImageUtility2.drawFloodFill4(x, y, d0, d1, fore, back,
        (px, py, value) => { sets.cb2(px, py); image.set(px, py, value); },
        (px, py) => { ++numGets; return image.get(px, py); });
    sets.emit(io, true);
    io.outInt(numGets);
    outImage2(io, image, 2);
    fillChecked.push({
        dims: [d0, d1, 1], before, seed: [x, y, 0], fore, back, after: image.getPixels().slice()
    });
}

// ---- ImageUtility3 ----

const NBR6 = [[-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]];
const NBR18: number[][] = [];
const NBR26: number[][] = [];
for (let dz = -1; dz <= 1; ++dz) {
    for (let dy = -1; dy <= 1; ++dy) {
        for (let dx = -1; dx <= 1; ++dx) {
            const m = Math.abs(dx) + Math.abs(dy) + Math.abs(dz);
            if (m > 0) { NBR26.push([dx, dy, dz]); }
            if (m > 0 && m < 3) { NBR18.push([dx, dy, dz]); }
        }
    }
}

function components3Case(io: OracleIO): void {
    const kind = io.integer();
    const image = readImage3(io);
    const d0 = image.getDimension(0);
    const d1 = image.getDimension(1);
    const input = image.getPixels().slice();
    let components: number[][];
    if (kind === 0) {
        components = ImageUtility3.getComponents6(image);
    } else if (kind === 1) {
        components = ImageUtility3.getComponents18(image);
    } else if (kind === 2) {
        components = ImageUtility3.getComponents26(image);
    } else {
        const nbrs = readOffsets(io, 3).map((d) => d[0] + d0 * (d[1] + d1 * d[2]));
        components = ImageUtility3.getComponents(image, nbrs);
    }
    io.outInt(components.length);
    for (let k = 1; k < components.length; ++k) {
        io.outInt(components[k].length);
        outPacked(io, components[k], 8);
    }
    outImage3(io, image, 8);
    componentsChecked.push({
        input, dims: [d0, d1, image.getDimension(2)],
        offsets: kind === 0 ? NBR6 : kind === 1 ? NBR18 : kind === 2 ? NBR26 : null,
        components, labels: image.getPixels().slice()
    });
}

function morph3Case(io: OracleIO, op: number): void {
    const kind = io.integer();
    const image = readImage3(io);
    const zeroExterior = (op !== 0 ? io.boolean() : false);
    const nbrs = (kind === 3 ? readOffsets(io, 3) : []);
    const U = ImageUtility3;
    let out: Image3<number>;
    switch (op) {
        case 0:
            out = kind === 0 ? U.dilate6(image) : kind === 1 ? U.dilate18(image)
                : kind === 2 ? U.dilate26(image) : U.dilate(image, nbrs);
            break;
        case 1:
            out = kind === 0 ? U.erode6(image, zeroExterior) : kind === 1
                ? U.erode18(image, zeroExterior) : kind === 2
                    ? U.erode26(image, zeroExterior) : U.erode(image, zeroExterior, nbrs);
            break;
        case 2:
            out = kind === 0 ? U.open6(image, zeroExterior) : kind === 1
                ? U.open18(image, zeroExterior) : kind === 2
                    ? U.open26(image, zeroExterior) : U.open(image, zeroExterior, nbrs);
            break;
        default:
            out = kind === 0 ? U.close6(image, zeroExterior) : kind === 1
                ? U.close18(image, zeroExterior) : kind === 2
                    ? U.close26(image, zeroExterior) : U.close(image, zeroExterior, nbrs);
            break;
    }
    io.outInt(out.getDimension(0));
    io.outInt(out.getDimension(1));
    io.outInt(out.getDimension(2));
    outImage3(io, out, 2);
    morphChecked.push({
        op, zeroExterior, offsets: kind === 0 ? NBR6 : kind === 1 ? NBR18 : kind === 2 ? NBR26 : nbrs,
        dims: [image.getDimension(0), image.getDimension(1), image.getDimension(2)],
        input: image.getPixels().slice(), output: out.getPixels().slice()
    });
}

interface ConvexRecord { dims: number[]; input: number[]; output: number[] }
const convexChecked: ConvexRecord[] = [];

function cdConvexCase(io: OracleIO): void {
    const image = readImage3(io);
    const input = image.getPixels().slice();
    ImageUtility3.computeCDConvex(image);
    outImage3(io, image, 1);
    convexChecked.push({
        dims: [image.getDimension(0), image.getDimension(1), image.getDimension(2)],
        input, output: image.getPixels().slice()
    });
}

function floodFill3IntCase(io: OracleIO): void {
    const image = readImage3(io);
    const x = io.integer();
    const y = io.integer();
    const z = io.integer();
    const back = io.integer();
    const fore = (back + io.integer()) % 4;
    const before = image.getPixels().slice();
    ImageUtility3.floodFill6(image, x, y, z, fore, back);
    outImage3(io, image, 2);
    fillChecked.push({
        dims: [image.getDimension(0), image.getDimension(1), image.getDimension(2)], before,
        seed: [x, y, z], fore, back, after: image.getPixels().slice()
    });
}

function floodFill3DoubleCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    const image = new Image3<number>(d0, d1, d2);
    for (let i = 0; i < d0 * d1 * d2; ++i) { image.set(i, io.real()); }
    const x = io.integer();
    const y = io.integer();
    const z = io.integer();
    const back = io.real();
    const fore = io.real();
    const before = image.getPixels().slice();
    ImageUtility3.floodFill6(image, x, y, z, fore, back);
    for (const v of image.getPixels()) { io.outReal(v); }
    fillChecked.push({
        dims: [d0, d1, d2], before, seed: [x, y, z], fore, back, after: image.getPixels().slice()
    });
}

function drawLine3Case(io: OracleIO): void {
    const mode = io.integer();
    const x0 = io.integer();
    const y0 = io.integer();
    const z0 = io.integer();
    let x1 = io.integer();
    let y1 = io.integer();
    let z1 = io.integer();
    const len = io.integer();
    const s = io.integer();
    switch (mode) {
        case 1:
            x1 = (s % 3 === 0 ? x1 : x0); y1 = (s % 3 === 1 ? y1 : y0); z1 = (s % 3 === 2 ? z1 : z0);
            break;
        case 2: x1 = x0 + len; y1 = y0 + ((s & 1) ? len : -len); z1 = z0 + ((s & 2) ? len : -len); break;
        case 3: x1 = x0; y1 = y0; z1 = z0; break;
        case 4: {
            const small = (len === 0 ? 0 : (len > 0 ? len - 1 : len + 1));
            const d = [len, (s & 1) ? len : -len, small];
            const k = s % 3;
            const r = [d[k % 3], d[(k + 1) % 3], d[(k + 2) % 3]];
            x1 = x0 + r[0]; y1 = y0 + r[1]; z1 = z0 + r[2];
            break;
        }
        default: break;
    }
    const rec = new PointRecorder();
    ImageUtility3.drawLine(x0, y0, z0, x1, y1, z1, rec.cb3);
    rec.emit(io, mode !== 5);
    drawChecked.push({ kind: 'line', args: [x0, y0, z0, x1, y1, z1], points: rec.points });
}

// ---- independent checks for the image utilities ----

function coords(i: number, dims: readonly number[]): number[] {
    const x = i % dims[0];
    const r = (i - x) / dims[0];
    const y = r % dims[1];
    return [x, y, (r - y) / dims[1]];
}

function indexOf(c: readonly number[], dims: readonly number[]): number {
    return c[0] + dims[0] * (c[1] + dims[1] * (c[2] ?? 0));
}

function inRange(c: readonly number[], dims: readonly number[]): boolean {
    return c.every((v, k) => 0 <= v && v < (dims[k] ?? 1));
}

function neighbour(c: readonly number[], d: readonly number[]): number[] {
    return [c[0] + d[0], c[1] + d[1], c[2] + (d[2] ?? 0)];
}

// Components partition the 1-pixels, the image holds the labels, and for the
// symmetric neighbourhoods each component is connected and no two are
// adjacent.
function checkComponents(r: ComponentsRecord): string[] {
    const bad: string[] = [];
    const n = r.input.length;
    const owner = new Array<number>(n).fill(0);
    for (let k = 1; k < r.components.length; ++k) {
        const c = r.components[k];
        if (c.length === 0) { bad.push(`component ${k} empty`); }
        for (let j = 0; j < c.length; ++j) {
            if (j > 0 && c[j] <= c[j - 1]) { bad.push(`component ${k} not increasing`); }
            if (owner[c[j]] !== 0) { bad.push(`pixel ${c[j]} in two components`); }
            owner[c[j]] = k;
        }
    }
    for (let i = 0; i < n; ++i) {
        if ((r.input[i] === 1) !== (owner[i] !== 0)) { bad.push(`pixel ${i}: foreground/component mismatch`); }
        if (r.labels[i] !== owner[i]) { bad.push(`pixel ${i}: label ${r.labels[i]} vs ${owner[i]}`); }
    }
    if (r.offsets !== null && bad.length === 0) {
        for (let k = 1; k < r.components.length; ++k) {
            const c = r.components[k];
            const seen = new Set<number>([c[0]]);
            const stack = [c[0]];
            while (stack.length > 0) {
                const p = coords(stack.pop()!, r.dims);
                for (const d of r.offsets) {
                    const q = neighbour(p, d);
                    if (!inRange(q, r.dims)) { continue; }
                    const qi = indexOf(q, r.dims);
                    if (owner[qi] !== 0 && owner[qi] !== k) { bad.push(`components ${k} and ${owner[qi]} adjacent`); }
                    if (owner[qi] === k && !seen.has(qi)) { seen.add(qi); stack.push(qi); }
                }
            }
            if (seen.size !== c.length) { bad.push(`component ${k} not connected`); }
        }
    }
    return bad;
}

// Direct definitions of dilation and erosion (the output pixel's view,
// where the library loops over input pixels).
function refDilate(dims: readonly number[], a: readonly number[], nbrs: Offsets): number[] {
    return a.map((v, i) => {
        const p = coords(i, dims);
        for (const d of nbrs) {
            const q = [p[0] - d[0], p[1] - d[1], p[2] - (d[2] ?? 0)];
            if (inRange(q, dims) && a[indexOf(q, dims)] === 1) { return 1; }
        }
        return v;
    });
}

function refErode(dims: readonly number[], a: readonly number[], nbrs: Offsets, zeroExterior: boolean): number[] {
    return a.map((v, i) => {
        if (v !== 1) { return v; }
        const p = coords(i, dims);
        for (const d of nbrs) {
            const q = neighbour(p, d);
            if (inRange(q, dims) ? a[indexOf(q, dims)] === 0 : zeroExterior) { return 0; }
        }
        return v;
    });
}

function checkMorph(r: MorphRecord): string[] {
    const { dims, input, offsets, zeroExterior } = r;
    const ref = r.op === 0 ? refDilate(dims, input, offsets)
        : r.op === 1 ? refErode(dims, input, offsets, zeroExterior)
            : r.op === 2 ? refDilate(dims, refErode(dims, input, offsets, zeroExterior), offsets)
                : refErode(dims, refDilate(dims, input, offsets), offsets, zeroExterior);
    const bad: string[] = [];
    if (ref.some((v, i) => v !== r.output[i])) { bad.push(`op ${r.op}: differs from the definition`); }
    const symmetric = offsets.every((d) => offsets.some((e) => e.every((v, k) => v === -d[k])));
    if (input.every((v) => v <= 1) && symmetric) {
        const sub = (a: readonly number[], b: readonly number[]) => a.every((v, i) => v <= b[i]);
        const ok = r.op === 0 ? sub(input, r.output) : r.op === 1 ? sub(r.output, input)
            : r.op === 2 ? sub(r.output, input) : (zeroExterior || sub(input, r.output));
        if (!ok) { bad.push(`op ${r.op}: set identity fails`); }
    }
    return bad;
}

function checkBoundary(r: BoundaryRecord): string[] {
    const bad: string[] = [];
    let first = r.start;
    while (first < r.before.length && r.before[first] === 0) { ++first; }
    if (r.success !== (first < r.before.length)) { bad.push('success flag'); }
    if (!r.success) { return r.boundary.length === 0 ? bad : ['boundary on failure']; }
    if (r.boundary[0] !== first) { bad.push('first boundary pixel'); }
    for (let j = 0; j < r.boundary.length; ++j) {
        if (r.before[r.boundary[j]] === 0) { bad.push(`boundary pixel ${r.boundary[j]} is background`); }
        const a = coords(r.boundary[j], r.dims);
        const b = coords(r.boundary[(j + 1) % r.boundary.length], r.dims);
        const cheb = Math.max(Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
        if (r.boundary.length > 1 && cheb !== 1) { bad.push(`step ${j} is not an 8-neighbour step`); }
    }
    return bad;
}

function sameValue(a: number, b: number): boolean {
    return Object.is(a, b) || (Number.isNaN(a) && Number.isNaN(b));
}

function checkFill(r: FillRecord): string[] {
    const expected = r.before.slice();
    if (inRange(r.seed, r.dims)) {
        const nbrs = r.dims[2] > 1 || r.seed[2] !== 0 ? NBR6 : NBR4;
        const s = indexOf(r.seed, r.dims);
        const seen = new Set<number>([s]);
        const stack = [s];
        while (stack.length > 0) {
            const i = stack.pop()!;
            expected[i] = r.fore;
            for (const d of nbrs) {
                const q = neighbour(coords(i, r.dims), d);
                if (!inRange(q, r.dims)) { continue; }
                const qi = indexOf(q, r.dims);
                if (!seen.has(qi) && r.before[qi] === r.back) { seen.add(qi); stack.push(qi); }
            }
        }
    }
    return expected.every((v, i) => sameValue(v, r.after[i])) ? [] : ['fill differs from the region definition'];
}

const SQRTF_INT32_MAX = Math.fround(Math.sqrt(2147483648));

// L1: the city-block distance to the nearest background pixel or to the
// image exterior. L2: the exact Euclidean distance to the nearest
// background pixel, rounded to float (every distance here is below 100).
function checkDistance(r: DistanceRecord): string[] {
    const [d0, d1] = r.dims;
    const zeros: number[][] = [];
    r.input.forEach((v, i) => { if (v === 0) { zeros.push(coords(i, r.dims)); } });
    const bad: string[] = [];
    let best = 0;
    let arg = [0, 0];
    for (let i = 0; i < r.input.length; ++i) {
        const [x, y] = coords(i, r.dims);
        let expected: number;
        if (r.l2) {
            let m = Infinity;
            for (const z of zeros) { m = Math.min(m, (x - z[0]) ** 2 + (y - z[1]) ** 2); }
            expected = m === Infinity ? SQRTF_INT32_MAX : Math.fround(Math.sqrt(m));
            if (expected > best) { best = expected; arg = [x, y]; }
        } else {
            let m = r.input[i] === 0 ? 0 : Math.min(x + 1, y + 1, d0 - x, d1 - y);
            for (const z of zeros) { m = Math.min(m, Math.abs(x - z[0]) + Math.abs(y - z[1])); }
            expected = m;
            best = Math.max(best, m);
        }
        if (r.output[i] !== expected) { bad.push(`pixel (${x},${y}): ${r.output[i]} vs ${expected}`); }
    }
    if (r.l2) {
        if (r.max !== best || r.arg[0] !== arg[0] || r.arg[1] !== arg[1]) { bad.push('maximum'); }
    } else {
        const m = Math.max(1, best);
        if (r.max !== m) { bad.push(`maximum ${r.max} vs ${m}`); }
        if (best <= 1 ? (r.arg[0] !== 0 || r.arg[1] !== 0) : r.output[indexOf(r.arg, r.dims)] !== m) {
            bad.push('maximum location');
        }
    }
    return bad;
}

// 8-connected components of the nonzero pixels.
function countComponents8(dims: readonly number[], a: readonly number[]): number {
    const seen = new Array<boolean>(a.length).fill(false);
    let count = 0;
    for (let i = 0; i < a.length; ++i) {
        if (a[i] === 0 || seen[i]) { continue; }
        ++count;
        seen[i] = true;
        const stack = [i];
        while (stack.length > 0) {
            const p = coords(stack.pop()!, dims);
            for (const d of NBR8) {
                const q = neighbour(p, d);
                if (!inRange(q, dims)) { continue; }
                const qi = indexOf(q, dims);
                if (a[qi] !== 0 && !seen[qi]) { seen[qi] = true; stack.push(qi); }
            }
        }
    }
    return count;
}

const skeletonStats = { records: 0, componentsChanged: 0, vanished: 0 };

function checkSkeleton(r: SkeletonRecord): string[] {
    const bad: string[] = [];
    if (r.output.some((v, i) => v !== 0 && (v !== 1 || r.input[i] !== 1))) { bad.push('not a binary subset'); }
    ++skeletonStats.records;
    const before = countComponents8(r.dims, r.input);
    const after = countComponents8(r.dims, r.output);
    if (before !== after) { ++skeletonStats.componentsChanged; }
    if (before > 0 && after === 0) { ++skeletonStats.vanished; }
    return bad;
}

// Bresenham lines: endpoints, one step per major-axis unit, and every minor
// coordinate within 1/2 of the ideal line (exact integer test).
function checkLine(args: readonly number[], pts: readonly number[][]): string[] {
    const p0 = args.slice(0, 3);
    const p1 = args.slice(3, 6);
    const d = p1.map((v, k) => v - p0[k]);
    const m = [0, 1, 2].reduce((a, k) => (Math.abs(d[k]) > Math.abs(d[a]) ? k : a), 0);
    const D = Math.abs(d[m]);
    const bad: string[] = [];
    if (pts.length !== D + 1) { bad.push(`count ${pts.length} vs ${D + 1}`); }
    pts.forEach((p, j) => {
        const t = Math.abs(p[m] - p0[m]);
        if (t !== j) { bad.push(`step ${j} on the major axis`); }
        for (let k = 0; k < 3; ++k) {
            if (Math.abs(2 * (p[k] - p0[k]) * D - 2 * d[k] * t) > D) { bad.push(`point ${j} off the line`); }
        }
    });
    if (pts.length > 0 && pts[pts.length - 1].some((v, k) => v !== p1[k])) { bad.push('last point'); }
    return bad;
}

function checkDraw(r: DrawRecord): string[] {
    const pts = r.points;
    const a = r.args;
    const key = (p: readonly number[]) => `${p[0]},${p[1]}`;
    const set = new Set(pts.map(key));
    switch (r.kind) {
        case 'thick': {
            const t = a[2];
            const n = t < 0 ? 0 : (2 * t + 1) ** 2;
            const ok = pts.length === n && pts.every((p) =>
                Math.abs(p[0] - a[0]) <= t && Math.abs(p[1] - a[1]) <= t) && set.size === n;
            return ok ? [] : ['thick pixel'];
        }
        case 'line':
            return checkLine(a, pts);
        case 'circle':
        case 'disk': {
            const r2 = a[2];
            const bad: string[] = [];
            if (r2 < 0) { return pts.length === 0 ? [] : ['negative radius drew']; }
            for (const p of pts) {
                const d2 = (p[0] - a[0]) ** 2 + (p[1] - a[1]) ** 2;
                if (d2 >= (r2 + 1) ** 2 || (r.kind === 'circle' && r2 > 0 && d2 <= (r2 - 1) ** 2)) {
                    bad.push(`point ${key(p)} at squared distance ${d2}`);
                }
                if (!set.has(key([2 * a[0] - p[0], p[1]])) || !set.has(key([p[0], 2 * a[1] - p[1]]))
                    || !set.has(key([a[0] + p[1] - a[1], a[1] + p[0] - a[0]]))) {
                    bad.push(`point ${key(p)} without its mirror images`);
                }
            }
            if (r.kind === 'disk') {
                for (let y = -r2; y <= r2; ++y) {
                    for (let x = -r2; x <= r2; ++x) {
                        if (x * x + y * y <= (r2 - 1) ** 2 && !set.has(key([a[0] + x, a[1] + y]))) {
                            bad.push(`interior point (${x},${y}) missing`);
                        }
                    }
                }
            }
            return bad.slice(0, 3);
        }
        case 'rect':
        case 'box': {
            const [x0, y0, x1, y1] = a;
            if (x0 > x1 || y0 > y1) { return []; }
            const expected = new Set<string>();
            for (let y = y0; y <= y1; ++y) {
                for (let x = x0; x <= x1; ++x) {
                    if (r.kind === 'box' || x === x0 || x === x1 || y === y0 || y === y1) { expected.add(key([x, y])); }
                }
            }
            const ok = expected.size === set.size && [...set].every((k) => expected.has(k));
            return ok ? [] : [`${r.kind} pixel set`];
        }
        default: {
            // Ellipse: every pixel's 3x3 block straddles the curve
            // b^2 X^2 + a^2 Y^2 = a^2 b^2, and the set is symmetric.
            const [xc, yc, ea, eb] = a;
            if (ea < 0 || eb < 0) { return []; }
            const F = (x: number, y: number) => eb * eb * x * x + ea * ea * y * y - ea * ea * eb * eb;
            const bad: string[] = [];
            for (const p of pts) {
                const X = p[0] - xc;
                const Y = p[1] - yc;
                let lo = Infinity;
                let hi = -Infinity;
                for (let dy = -1; dy <= 1; ++dy) {
                    for (let dx = -1; dx <= 1; ++dx) {
                        lo = Math.min(lo, F(X + dx, Y + dy));
                        hi = Math.max(hi, F(X + dx, Y + dy));
                    }
                }
                if (lo > 0 || hi < 0) { bad.push(`pixel ${key(p)} away from the ellipse`); }
                if (!set.has(key([xc - X, yc + Y])) || !set.has(key([xc + X, yc - Y]))) {
                    bad.push(`pixel ${key(p)} without its mirror images`);
                }
            }
            return bad.slice(0, 3);
        }
    }
}

function checkConvex(r: ConvexRecord): string[] {
    const bad: string[] = [];
    for (let i = 0; i < r.input.length; ++i) {
        const p = coords(i, r.dims);
        let open = false;
        for (const d of NBR6) {
            let q = p;
            let clear = true;
            while (inRange(q, r.dims)) {
                if ((r.input[indexOf(q, r.dims)] & 1) !== 0) { clear = false; break; }
                q = neighbour(q, d);
            }
            open = open || clear;
        }
        if (r.output[i] !== (open ? 0 : 1)) { bad.push(`voxel ${i}`); }
    }
    return bad;
}

// ---- surface extractors ----

function bitsFor(range: number): number {
    let b = 1;
    while (2 ** b <= range) { ++b; }
    return b;
}

function readVoxels(io: OracleIO, n: number, lo: number, hi: number): number[] {
    return unpack(io, n, bitsFor(hi - lo)).map((u) => u + lo);
}

function cloneTriangles(ts: readonly SurfaceExtractorTriangle[]): SurfaceExtractorTriangle[] {
    return ts.map((t) => {
        const c = new SurfaceExtractorTriangle();
        c.v = [t.v[0], t.v[1], t.v[2]];
        return c;
    });
}

function digestRational(d: Digest, vs: readonly SurfaceExtractorVertex[]): void {
    for (const v of vs) {
        d.integer(v.xNumer); d.integer(v.xDenom); d.integer(v.yNumer);
        d.integer(v.yDenom); d.integer(v.zNumer); d.integer(v.zDenom);
    }
}

function digestTriangles(d: Digest, ts: readonly SurfaceExtractorTriangle[]): void {
    for (const t of ts) { d.integer(t.v[0]); d.integer(t.v[1]); d.integer(t.v[2]); }
}

interface SurfaceRecord {
    kind: 'cubes' | 'tetra';
    dims: number[];
    values: number[];
    level: number;
    rational: SurfaceExtractorVertex[];
    triangles: SurfaceExtractorTriangle[];
}
const surfaceChecked: SurfaceRecord[] = [];

// Mirrors EmitSurface.
function emitSurface(io: OracleIO, ex: SurfaceExtractor, level: number, sameDir: boolean,
    allowFull: boolean): void {
    const r = ex.extractRational(level);
    io.outInt(r.vertices.length);
    io.outInt(r.triangles.length);
    const d0 = new Digest();
    digestRational(d0, r.vertices);
    digestTriangles(d0, r.triangles);
    io.outInt(d0.h);

    ex.makeUnique(r.vertices, r.triangles);
    io.outInt(r.vertices.length);
    io.outInt(r.triangles.length);
    const d1 = new Digest();
    digestRational(d1, r.vertices);
    digestTriangles(d1, r.triangles);
    io.outInt(d1.h);

    const a = ex.extract(level, false);
    io.outInt(a.vertices.length);
    const d2 = new Digest();
    for (const v of a.vertices) { d2.real(v[0]); d2.real(v[1]); d2.real(v[2]); }
    digestTriangles(d2, a.triangles);
    io.outInt(d2.h);

    const b = ex.extract(level, true);
    io.outInt(b.vertices.length);
    io.outInt(b.triangles.length);
    const full = allowFull && b.triangles.length <= 16;
    const d3 = new Digest();
    for (const v of b.vertices) {
        for (let k = 0; k < 3; ++k) { if (full) { io.outReal(v[k]); } else { d3.real(v[k]); } }
    }
    for (const t of b.triangles) {
        for (let k = 0; k < 3; ++k) { if (full) { io.outInt(t.v[k]); } else { d3.integer(t.v[k]); } }
    }
    if (!full) { io.outInt(d3.h); }

    const ot = cloneTriangles(b.triangles);
    ex.orientTriangles(b.vertices, ot, sameDir);
    outPacked(io, ot.map((t, i) => (t.v[1] !== b.triangles[i].v[1] ? 1 : 0)), 1);

    const normals = ex.computeNormals(b.vertices, ot);
    const d4 = new Digest();
    for (const n of normals) {
        for (let k = 0; k < 3; ++k) { if (full) { io.outReal(n[k]); } else { d4.real(n[k]); } }
    }
    if (!full) { io.outInt(d4.h); }
}

const TYPE_LO = [-128, -32768, -(2 ** 20), 0, 0, 0];
const TYPE_HI = [127, 32767, 2 ** 20, 255, 65535, 2 ** 20];

function surfaceCase(io: OracleIO, kind: 'cubes' | 'tetra', typed: boolean, lo: number, hi: number,
    allowFull: boolean, check: boolean): void {
    let type = 2;
    if (typed) { type = io.integer(); lo = TYPE_LO[type]; hi = TYPE_HI[type]; }
    const dims = [io.integer(), io.integer(), io.integer()];
    const values = readVoxels(io, dims[0] * dims[1] * dims[2], lo, hi);
    const level = io.integer();
    const sameDir = io.boolean();
    const ex = kind === 'cubes'
        ? new SurfaceExtractorCubes(dims[0], dims[1], dims[2], values)
        : new SurfaceExtractorTetrahedra(dims[0], dims[1], dims[2], values);
    if (check) {
        // Recorded before the outputs are compared: a disagreeing integer
        // output ends the replay of the record.
        const r = ex.extractRational(level);
        ex.makeUnique(r.vertices, r.triangles);
        surfaceChecked.push({ kind, dims, values, level, rational: r.vertices, triangles: r.triangles });
    }
    emitSurface(io, ex, level, sameDir, allowFull);
}

// Mirrors OrientPoints: vertices, index triples (Triangle's rotation), the
// direction; outputs the swap bits and the normals.
function orientPointsCase(io: OracleIO, ex: SurfaceExtractor): void {
    const nv = io.integer();
    io.integer();  // the point mode
    const vs: [number, number, number][] = [];
    for (let i = 0; i < nv; ++i) { vs.push([io.real(), io.real(), io.real()]); }
    const nt = io.integer();
    const ts: SurfaceExtractorTriangle[] = [];
    for (let t = 0; t < nt; ++t) {
        const i0 = io.integer();
        const i1 = io.integer();
        const i2 = io.integer();
        ts.push(new SurfaceExtractorTriangle(i0, i1, i2));
    }
    const sameDir = io.boolean();
    const ot = cloneTriangles(ts);
    ex.orientTriangles(vs, ot, sameDir);
    outPacked(io, ot.map((t, i) => (t.v[1] !== ts[i].v[1] ? 1 : 0)), 1);
    for (const n of ex.computeNormals(vs, ot)) { io.outReal(n[0]); io.outReal(n[1]); io.outReal(n[2]); }
}

function surfacePointsCase(io: OracleIO, kind: 'cubes' | 'tetra'): void {
    const dims = [io.integer(), io.integer(), io.integer()];
    const values = readVoxels(io, dims[0] * dims[1] * dims[2], -20, 20);
    const level = io.integer();
    const ex = kind === 'cubes'
        ? new SurfaceExtractorCubes(dims[0], dims[1], dims[2], values)
        : new SurfaceExtractorTetrahedra(dims[0], dims[1], dims[2], values);
    ex.extractRational(level);
    orientPointsCase(io, ex);
}

function invalidBoundsCase(io: OracleIO): void {
    const which = io.integer();
    const d = [io.integer(), io.integer(), io.integer()];
    const voxels = new Array<number>(64).fill(1);
    if (which === 0) {
        new SurfaceExtractorCubes(d[0], d[1], d[2], voxels);
    } else {
        new SurfaceExtractorTetrahedra(d[0], d[1], d[2], voxels);
    }
    io.outInt(d[0] * d[1] * d[2]);
}

// ---- SurfaceExtractorMC ----

interface MCVoxelRecord { F: number[]; level: number; perturb: number; valid: boolean; vertices: number[][]; pairs: number[][] }
const mcVoxelChecked: MCVoxelRecord[] = [];

function mcVoxelCase(io: OracleIO, withVertices: boolean): void {
    io.integer();  // the value mode
    const F: number[] = [];
    for (let i = 0; i < 8; ++i) { F.push(io.real()); }
    const level = io.real();
    const perturb = io.real();
    const mc = new SurfaceExtractorMC(new Image3<number>(2, 2, 2));
    const { valid, mesh } = mc.extractVoxel(level, perturb, F);
    const t = mesh.topology;
    const vertices = mesh.vertices.slice(0, t.numVertices).map((v) => [...v.values]);
    if (valid) { mcVoxelChecked.push({ F, level, perturb, valid, vertices, pairs: t.vpair.slice(0, t.numVertices) }); }
    io.outBool(valid);
    if (!valid) { return; }
    io.outInt(t.numVertices);
    io.outInt(t.numTriangles);
    const packed: number[] = [];
    for (let i = 0; i < t.numVertices; ++i) { packed.push(t.vpair[i][0], t.vpair[i][1]); }
    for (let i = 0; i < t.numTriangles; ++i) { packed.push(...t.itriple[i]); }
    outPacked(io, packed, 4);
    if (withVertices) {
        for (const v of vertices) { io.outReal(v[0]); io.outReal(v[1]); io.outReal(v[2]); }
    }
}

function readMCImage(io: OracleIO): Image3<number> {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    io.integer();  // the value mode
    const image = new Image3<number>(d0, d1, d2);
    for (let i = 0; i < d0 * d1 * d2; ++i) { image.set(i, io.real()); }
    return image;
}

interface MCImageRecord {
    image: Image3<number>;
    level: number;
    perturb: number;
    vertices: number[][];
    indices: number[];
    unique: number[][];
    uniqueIndices: number[];
}
const mcImageChecked: MCImageRecord[] = [];

function mcImageCase(io: OracleIO, topologyOnly: boolean, allowFull: boolean): void {
    const image = readMCImage(io);
    const level = io.real();
    const perturb = io.real();
    const sameDir = io.boolean();
    const mc = new SurfaceExtractorMC(image);
    const r = mc.extract(level, perturb);
    // Recorded before the outputs are compared (see surfaceCase).
    const rec: MCImageRecord = {
        image, level, perturb, vertices: r.vertices.map((v) => [...v.values]), indices: r.indices,
        unique: [], uniqueIndices: []
    };
    if (!topologyOnly && r.vertices.length > 0) {
        const u0 = mc.makeUnique(r.vertices, r.indices);
        rec.unique = u0.vertices.map((v) => [...v.values]);
        rec.uniqueIndices = u0.indices;
    }
    mcImageChecked.push(rec);
    io.outInt(r.vertices.length);
    io.outInt(r.indices.length);
    const full = allowFull && r.indices.length <= 48;
    const emitMesh = (vs: readonly Vector[], is: readonly number[], withVertices: boolean): void => {
        const d = new Digest();
        if (withVertices) {
            for (const v of vs) {
                for (let k = 0; k < 3; ++k) { if (full) { io.outReal(v.get(k)); } else { d.real(v.get(k)); } }
            }
        }
        for (const i of is) { if (full) { io.outInt(i); } else { d.integer(i); } }
        if (!full) { io.outInt(d.h); }
    };
    emitMesh(r.vertices, r.indices, !topologyOnly);
    if (topologyOnly || r.vertices.length === 0) { return; }
    const u = mc.makeUnique(r.vertices, r.indices);
    io.outInt(u.vertices.length);
    emitMesh(u.vertices, u.indices, true);
    const oriented = u.indices.slice();
    mc.orientTriangles(u.vertices, oriented, sameDir);
    const swapped: number[] = [];
    for (let t = 0; 3 * t < oriented.length; ++t) { swapped.push(oriented[3 * t + 1] !== u.indices[3 * t + 1] ? 1 : 0); }
    outPacked(io, swapped, 1);
    if (oriented.length === 3) {
        for (const n of mc.computeNormals(u.vertices, oriented)) { io.outVec(n); }
    }
}

function mcMakeUniqueCase(io: OracleIO): void {
    const mode = io.integer();
    const nv = (mode === 2 ? 0 : io.integer());
    const vertices: Vector[] = [];
    for (let i = 0; i < nv; ++i) { vertices.push(Vector.fromArray([io.real(), io.real(), io.real()])); }
    const ni = (mode === 3 ? io.integer() : 3 * io.integer());
    const indices: number[] = [];
    for (let i = 0; i < ni; ++i) { indices.push(io.integer()); }
    if (mode === 4) {
        const bad = io.integer();
        indices[bad] = (io.boolean() ? nv : -1);
    }
    const mc = new SurfaceExtractorMC(new Image3<number>(2, 2, 2));
    const u = mc.makeUnique(vertices, indices);
    io.outInt(u.vertices.length);
    for (const v of u.vertices) { io.outVec(v); }
    for (const i of u.indices) { io.outInt(i); }
}

function mcPointsCase(io: OracleIO): void {
    const image = readMCImage(io);
    const nv = io.integer();
    io.integer();  // the point mode
    const vertices: Vector[] = [];
    for (let i = 0; i < nv; ++i) { vertices.push(Vector.fromArray([io.real(), io.real(), io.real()])); }
    const nt = io.integer();
    const indices: number[] = [];
    for (let t = 0; t < 3 * nt; ++t) { indices.push(io.integer()); }
    const sameDir = io.boolean();
    const mc = new SurfaceExtractorMC(image);
    const oriented = indices.slice();
    mc.orientTriangles(vertices, oriented, sameDir);
    const swapped: number[] = [];
    for (let t = 0; t < nt; ++t) { swapped.push(oriented[3 * t + 1] !== indices[3 * t + 1] ? 1 : 0); }
    outPacked(io, swapped, 1);
    if (nt === 1) {
        for (const n of mc.computeNormals(vertices, oriented)) { io.outVec(n); }
    }
}

function mcNormalsCase(io: OracleIO): void {
    const nv = io.integer();
    io.integer();  // the point mode
    const vertices: Vector[] = [];
    for (let i = 0; i < nv; ++i) { vertices.push(Vector.fromArray([io.real(), io.real(), io.real()])); }
    const nt = io.integer();
    const indices: number[] = [];
    for (let t = 0; t < 3 * nt; ++t) { indices.push(io.integer()); }
    const mc = new SurfaceExtractorMC(new Image3<number>(2, 2, 2));
    for (const n of mc.computeNormals(vertices, indices)) { io.outVec(n); }
}

// ---- independent checks for the extracted surfaces (exact, BigInt) ----

function vertexKey(v: SurfaceExtractorVertex): string {
    const g = (a: number, b: number): number => { a = Math.abs(a); while (b !== 0) { [a, b] = [b, a % b]; } return a || 1; };
    const r = (n: number, d: number) => `${n / g(n, d)}/${d / g(n, d)}`;
    return `${r(v.xNumer, v.xDenom)},${r(v.yNumer, v.yDenom)},${r(v.zNumer, v.zDenom)}`;
}

function big(v: SurfaceExtractorVertex): { n: bigint[], d: bigint[] } {
    return {
        n: [BigInt(v.xNumer), BigInt(v.yNumer), BigInt(v.zNumer)],
        d: [BigInt(v.xDenom), BigInt(v.yDenom), BigInt(v.zDenom)]
    };
}

// Cubes: the trilinear interpolant of 2v - 2L - 1 vanishes at every vertex.
function cubesOnLevel(r: SurfaceRecord, v: SurfaceExtractorVertex): boolean {
    const { n, d } = big(v);
    const c = n.map((nk, k) => Math.min(r.dims[k] - 2, Math.max(0, Number(nk / d[k]))));
    const t = n.map((nk, k) => nk - BigInt(c[k]) * d[k]);
    let sum = 0n;
    for (let corner = 0; corner < 8; ++corner) {
        const b = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
        let w = 1n;
        for (let k = 0; k < 3; ++k) { w *= (b[k] === 1 ? t[k] : d[k] - t[k]); }
        const p = [c[0] + b[0], c[1] + b[1], c[2] + b[2]];
        sum += BigInt(2 * r.values[indexOf(p, r.dims)] - 2 * r.level - 1) * w;
    }
    return sum === 0n;
}

// The tetrahedra of the cube at the origin, by parity (x ^ y ^ z odd first),
// as corner offsets, in upstream's Extract order.
const TETRA_ODD = [
    [[1, 0, 0], [1, 1, 0], [0, 0, 0], [1, 0, 1]], [[0, 1, 0], [0, 0, 0], [1, 1, 0], [0, 1, 1]],
    [[0, 0, 1], [0, 1, 1], [1, 0, 1], [0, 0, 0]], [[1, 1, 1], [1, 0, 1], [0, 1, 1], [1, 1, 0]],
    [[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 0]]];
const TETRA_EVEN = [
    [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], [[1, 1, 0], [0, 1, 0], [1, 0, 0], [1, 1, 1]],
    [[1, 0, 1], [0, 0, 1], [1, 1, 1], [1, 0, 0]], [[0, 1, 1], [1, 1, 1], [0, 0, 1], [0, 1, 0]],
    [[1, 1, 1], [0, 1, 0], [1, 0, 0], [0, 0, 1]]];

// Whether the vertex lies in the closed tetrahedron with integer corners
// cs and the linear interpolant of v - L vanishes there.
function inTetraOnLevel(r: SurfaceRecord, cs: number[][], v: SurfaceExtractorVertex): boolean {
    const { n, d } = big(v);
    const D = d[0] * d[1] * d[2];
    const P = [n[0] * d[1] * d[2], n[1] * d[0] * d[2], n[2] * d[0] * d[1]];
    const c0 = cs[0].map(BigInt);
    const M = [1, 2, 3].map((j) => cs[j].map((x, k) => BigInt(x) - c0[k]));  // columns
    const a = (i: number, j: number) => M[j][i];
    const det = a(0, 0) * (a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)) - a(0, 1) * (a(1, 0) * a(2, 2) - a(1, 2) * a(2, 0))
        + a(0, 2) * (a(1, 0) * a(2, 1) - a(1, 1) * a(2, 0));
    const q = P.map((x, k) => x - c0[k] * D);
    // adj(M) q: row i of the adjugate is the cross product of columns j, k.
    const cross = (u: bigint[], w: bigint[]) => [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
    const rows = [cross(M[1], M[2]), cross(M[2], M[0]), cross(M[0], M[1])];
    const L = rows.map((row) => row[0] * q[0] + row[1] * q[1] + row[2] * q[2]);
    const L0 = det * D - L[0] - L[1] - L[2];
    if ([L0, ...L].some((x) => x * det < 0n)) { return false; }
    const f = cs.map((c) => BigInt(r.values[indexOf(c, r.dims)] - r.level));
    return f[0] * L0 + f[1] * L[0] + f[2] * L[1] + f[3] * L[2] === 0n;
}

// Tetrahedra: every triangle lies in one tetrahedron of the decomposition,
// on the zero set of that tetrahedron's linear interpolant.
function tetraTriangleOnLevel(r: SurfaceRecord, tri: SurfaceExtractorVertex[]): boolean {
    const lo = [0, 1, 2].map((k) => Math.min(...tri.map((v) => Number(big(v).n[k] / big(v).d[k]))));
    const ranges = [0, 1, 2].map((k) => [Math.max(0, lo[k] - 1), Math.min(r.dims[k] - 2, lo[k])]);
    for (let z = ranges[2][0]; z <= ranges[2][1]; ++z) {
        for (let y = ranges[1][0]; y <= ranges[1][1]; ++y) {
            for (let x = ranges[0][0]; x <= ranges[0][1]; ++x) {
                const tetras = ((x ^ y ^ z) & 1) ? TETRA_ODD : TETRA_EVEN;
                for (const t of tetras) {
                    const cs = t.map((o) => [x + o[0], y + o[1], z + o[2]]);
                    if (tri.every((v) => inTetraOnLevel(r, cs, v))) { return true; }
                }
            }
        }
    }
    return false;
}

// The crossing of the level on the grid edge from corner a to corner b
// (shifted values fa, fb of opposite signs), as a vertex key.
function crossingKey(a: number[], b: number[], fa: number, fb: number): string {
    const n: number[] = [];
    const d: number[] = [];
    for (let k = 0; k < 3; ++k) {
        // a[k] + (b[k] - a[k]) * fa / (fa - fb)
        n.push(a[k] * (fa - fb) + (b[k] - a[k]) * fa);
        d.push(fa - fb);
    }
    return vertexKey(new SurfaceExtractorVertex(n[0], d[0], n[1], d[1], n[2], d[2]));
}

const surfaceStats = {
    cubesSaddleFaces: 0, cubesPlusFaces: 0, cubesFourCrossingRecords: 0, cubesNonManifoldRecords: 0,
    tetraOpenRecords: 0, tetraZeroRecords: 0, mcOnLevelRecords: 0, mcClosedChecked: 0, mcAmbiguousRecords: 0
};

// Undirected edge -> number of incident triangles.
function edgeCounts(tris: readonly number[][]): Map<string, number> {
    const m = new Map<string, number>();
    for (const t of tris) {
        for (let j = 0; j < 3; ++j) {
            const a = t[j];
            const b = t[(j + 1) % 3];
            const key = a < b ? `${a},${b}` : `${b},${a}`;
            m.set(key, (m.get(key) ?? 0) + 1);
        }
    }
    return m;
}

// Interior edges (not in one image boundary plane) must have exactly two
// triangles; boundary-plane edges one or two.
function openEdges(tris: readonly number[][], onBoundaryPlane: (a: number, b: number) => boolean): number {
    let open = 0;
    for (const [key, count] of edgeCounts(tris)) {
        const [a, b] = key.split(',').map(Number);
        if (count === 2 || (count === 1 && onBoundaryPlane(a, b))) { continue; }
        ++open;
    }
    return open;
}

function checkSurface(r: SurfaceRecord): string[] {
    const bad: string[] = [];
    const keys = r.rational.map(vertexKey);
    const index = new Map<string, number>();
    keys.forEach((k, i) => { if (!index.has(k)) { index.set(k, i); } });
    const tris = r.triangles.map((t) => [t.v[0], t.v[1], t.v[2]]);
    for (const t of tris) {
        if (t[0] === t[1] || t[1] === t[2] || t[2] === t[0]) { bad.push('degenerate triangle'); }
    }
    const boundary = (a: number, b: number): boolean => [0, 1, 2].some((k) => {
        const va = big(r.rational[a]);
        const vb = big(r.rational[b]);
        return [0, r.dims[k] - 1].some((e) => va.n[k] === BigInt(e) * va.d[k] && vb.n[k] === BigInt(e) * vb.d[k]);
    });
    if (r.kind === 'cubes') {
        r.rational.forEach((v, i) => { if (!cubesOnLevel(r, v)) { bad.push(`vertex ${i} off the level set`); } });
        // Saddle faces: the two segments cut off the corners whose sign
        // differs from the bilinear saddle value (sign(det) * sign(f00)).
        const s = (p: number[]) => 2 * r.values[indexOf(p, r.dims)] - 2 * r.level - 1;
        const edges = edgeCounts(tris);
        const hasEdge = (a: string, b: string) => {
            const ia = index.get(a);
            const ib = index.get(b);
            return ia !== undefined && ib !== undefined && edges.has(ia < ib ? `${ia},${ib}` : `${ib},${ia}`);
        };
        let plus = false;
        const saddleCells = new Set<string>();
        for (let axis = 0; axis < 3; ++axis) {
            const u = (axis + 1) % 3;
            const w = (axis + 2) % 3;
            for (let i = 0; i < r.values.length; ++i) {
                const p = coords(i, r.dims);
                if (p[u] + 1 >= r.dims[u] || p[w] + 1 >= r.dims[w]) { continue; }
                const q = [p.slice(), p.slice(), p.slice(), p.slice()];
                q[1][u] += 1; q[2][u] += 1; q[2][w] += 1; q[3][w] += 1;
                const f = q.map(s);
                if (!(f[0] * f[1] < 0 && f[1] * f[2] < 0 && f[2] * f[3] < 0)) { continue; }
                const det = f[0] * f[2] - f[3] * f[1];
                if (det === 0) { ++surfaceStats.cubesPlusFaces; plus = true; continue; }
                ++surfaceStats.cubesSaddleFaces;
                saddleCells.add(`${axis}:${p[axis]}:${[0, 1, 2].filter((j) => j !== axis).map((j) => p[j]).join(',')}`);
                const P = [0, 1, 2, 3].map((j) => crossingKey(q[j], q[(j + 1) % 4], f[j], f[(j + 1) % 4]));
                // P[j] lies on the edge from corner j to corner j+1.
                const ok = det > 0 ? hasEdge(P[0], P[1]) && hasEdge(P[2], P[3])
                    : hasEdge(P[3], P[0]) && hasEdge(P[1], P[2]);
                if (!ok) { bad.push(`saddle face at ${p} (axis ${axis}) paired against the bilinear saddle`); }
            }
        }
        // RemoveTriangles ear-clips each voxel's wireframe without geometry
        // (the lowest-numbered vertex of degree 2 is always the next ear),
        // so on a voxel with a four-crossing face a fan diagonal can lie in
        // that face, and the two voxels sharing the face choose their
        // diagonals independently: edges of three or four triangles and
        // single-triangle edges inside the image (upstream suspect 3 of the
        // v26 report, preserved). The mesh is required to be closed only
        // when no face has four crossings; otherwise the defect is counted.
        const open = openEdges(tris, boundary);
        if (saddleCells.size === 0 && !plus) {
            if (open > 0) { bad.push(`${open} open edges without a four-crossing face`); }
        } else {
            ++surfaceStats.cubesFourCrossingRecords;
            if (open > 0) { ++surfaceStats.cubesNonManifoldRecords; }
        }
    } else {
        r.triangles.forEach((t, i) => {
            if (!tetraTriangleOnLevel(r, t.v.map((j) => r.rational[j]))) { bad.push(`triangle ${i} not on a tetrahedron's level set`); }
        });
        if (r.values.some((v) => v === r.level)) {
            ++surfaceStats.tetraZeroRecords;
        } else if (openEdges(tris, boundary) > 0) {
            ++surfaceStats.tetraOpenRecords;
            bad.push('open edges');
        }
    }
    return bad;
}

// MC vertices on the level set of the linear interpolant along their edge:
// the residual (F0 - L) - t (F0 - F1) is at rounding level.
function mcEdgeResidualOk(F0: number, F1: number, L: number, t: number): boolean {
    const scale = Math.abs(F0) + Math.abs(F1) + Math.abs(L) + Number.MIN_VALUE;
    return Math.abs((F0 - L) - t * (F0 - F1)) <= 2 ** -49 * scale && t >= 0 && t <= 1;
}

function checkMCVoxel(r: MCVoxelRecord): string[] {
    const bad: string[] = [];
    r.vertices.forEach((v, i) => {
        const [j0, j1] = r.pairs[i];
        const k0 = [j0 & 1, (j0 >> 1) & 1, (j0 >> 2) & 1];
        const k1 = [j1 & 1, (j1 >> 1) & 1, (j1 >> 2) & 1];
        for (let k = 0; k < 3; ++k) {
            if (k0[k] === k1[k]) {
                if (v[k] !== k0[k]) { bad.push(`vertex ${i} off its edge`); }
            } else if (!mcEdgeResidualOk(r.F[j0], r.F[j1], r.level, Math.abs(v[k] - k0[k]))) {
                bad.push(`vertex ${i} off the level set`);
            }
        }
    });
    return bad;
}

function checkMCImage(r: MCImageRecord): string[] {
    const bad: string[] = [];
    const dims = [0, 1, 2].map((k) => r.image.getDimension(k));
    const F = (p: number[]) => r.image.get(indexOf(p, dims));
    r.vertices.forEach((v, i) => {
        const frac = [0, 1, 2].filter((k) => !Number.isInteger(v[k]));
        const a = v.map(Math.floor);
        if (frac.length > 1) { bad.push(`vertex ${i} not on a grid edge`); return; }
        const k = frac.length === 1 ? frac[0] : 0;
        const b = a.slice();
        if (b[k] + 1 < dims[k]) { b[k] += 1; } else { a[k] -= 1; }
        const t = v[k] - a[k];
        if (!mcEdgeResidualOk(F(a), F(b), r.level, t)) { bad.push(`vertex ${i} off the level set`); }
    });
    if (r.uniqueIndices.length === 0) { return bad; }
    // Watertightness and consistent orientation of the deduplicated mesh
    // when no voxel lies on the level.
    const cls = r.image.getPixels().map((f) => {
        let g = f - r.level;
        if (g === 0) { g += r.perturb; }
        return g < 0 ? -1 : g > 0 ? 1 : 0;
    });
    // A voxel on the level (perturbed into a class) puts the vertices of
    // all its crossing edges on the voxel corner; MakeUnique merges them
    // and the mesh pinches there, so closedness is not required.
    if (cls.includes(0) || r.image.getPixels().some((f) => f - r.level === 0)) {
        ++surfaceStats.mcOnLevelRecords;
        return bad;
    }
    let ambiguous = false;
    for (let axis = 0; axis < 3 && !ambiguous; ++axis) {
        const u = (axis + 1) % 3;
        const w = (axis + 2) % 3;
        for (let i = 0; i < cls.length && !ambiguous; ++i) {
            const p = coords(i, dims);
            if (p[u] + 1 >= dims[u] || p[w] + 1 >= dims[w]) { continue; }
            const q = [p.slice(), p.slice(), p.slice(), p.slice()];
            q[1][u] += 1; q[2][u] += 1; q[2][w] += 1; q[3][w] += 1;
            const c = q.map((x) => cls[indexOf(x, dims)]);
            ambiguous = c[0] === c[2] && c[1] === c[3] && c[0] !== c[1];
        }
    }
    const tris: number[][] = [];
    for (let t = 0; 3 * t < r.uniqueIndices.length; ++t) { tris.push(r.uniqueIndices.slice(3 * t, 3 * t + 3)); }
    const boundary = (a: number, b: number) => [0, 1, 2].some((k) =>
        [0, dims[k] - 1].some((e) => r.unique[a][k] === e && r.unique[b][k] === e));
    const open = openEdges(tris.filter((t) => t[0] !== t[1] && t[1] !== t[2] && t[2] !== t[0]), boundary);
    // The table resolves every ambiguous face by its signs alone (see the
    // table check below), so the mesh is closed with or without them.
    if (ambiguous) { ++surfaceStats.mcAmbiguousRecords; }
    ++surfaceStats.mcClosedChecked;
    if (open > 0) { bad.push(`${open} open edges${ambiguous ? ' (ambiguous faces)' : ''}`); }
    const directed = new Set<string>();
    for (const t of tris) {
        if (t[0] === t[1] || t[1] === t[2] || t[2] === t[0]) { continue; }
        for (let j = 0; j < 3; ++j) {
            const key = `${t[j]},${t[(j + 1) % 3]}`;
            if (directed.has(key)) { bad.push('inconsistent orientation'); }
            directed.add(key);
        }
    }
    return bad;
}

// Exhaustive check of the MarchingCubes table used by SurfaceExtractorMC:
// on every face whose corner classes alternate (192 face instances over the
// 256 entries), the triangles contain exactly two segments in the face and
// both cut off the negative corners (F < level). The resolution therefore
// depends on the face's signs alone, and two voxels sharing a face agree:
// the table is face-consistent.
function checkMarchingCubesFaces(): { ambiguous: number, bad: string[] } {
    const mc = new MarchingCubes();
    const bad: string[] = [];
    let ambiguous = 0;
    for (let entry = 0; entry < 256; ++entry) {
        const t = mc.getTable(entry);
        const neg = (k: number) => ((entry >> k) & 1) === 1;
        for (let a = 0; a < 3; ++a) {
            for (let side = 0; side < 2; ++side) {
                const [u, w] = [0, 1, 2].filter((k) => k !== a);
                const at = (bu: number, bw: number) => (side << a) | (bu << u) | (bw << w);
                const cyc = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)];
                const s = cyc.map(neg);
                if (!(s[0] === s[2] && s[1] === s[3] && s[0] !== s[1])) { continue; }
                ++ambiguous;
                const onFace = (i: number) => t.vpair[i].every((k) => ((k >> a) & 1) === side);
                const cut: number[] = [];
                for (let tri = 0; tri < t.numTriangles; ++tri) {
                    const v = t.itriple[tri];
                    for (let j = 0; j < 3; ++j) {
                        const p = v[j];
                        const q = v[(j + 1) % 3];
                        if (onFace(p) && onFace(q)) {
                            const shared = t.vpair[p].find((k) => t.vpair[q].includes(k));
                            cut.push(shared === undefined ? -1 : shared);
                        }
                    }
                }
                if (cut.length !== 2 || cut.some((k) => k < 0 || !neg(k))) {
                    bad.push(`entry ${entry}, face ${a}/${side}: cut corners ${cut}`);
                }
            }
        }
    }
    return { ambiguous, bad };
}

function runChecks<T>(records: readonly T[], check: (r: T) => string[]): void {
    const bad: string[] = [];
    records.forEach((r, i) => { for (const m of check(r)) { bad.push(`record ${i}: ${m}`); } });
    expect(records.length).toBeGreaterThan(0);
    expect(bad.slice(0, 10)).toEqual([]);
}

// @@CASES@@

describe('oracle: v26-imaging', () => {
    const family = new OracleFamily('v26-imaging');
    const exact = { exact: true };

    family.case('ImageUtility2.getComponents', components2Case, exact);
    family.case('ImageUtility2.dilate', (io) => morph2Case(io, 0), exact);
    family.case('ImageUtility2.erode', (io) => morph2Case(io, 1), exact);
    family.case('ImageUtility2.open', (io) => morph2Case(io, 2), exact);
    family.case('ImageUtility2.close', (io) => morph2Case(io, 3), exact);
    family.case('ImageUtility2.extractBoundary', extractBoundaryCase, exact);
    family.case('ImageUtility2.floodFill4.int', floodFill2IntCase, exact);
    family.case('ImageUtility2.floodFill4.double', floodFill2DoubleCase, exact);
    family.case('ImageUtility2.getL1Distance', l1DistanceCase, exact);
    family.case('ImageUtility2.getL2Distance', l2DistanceCase, exact);
    family.case('ImageUtility2.getL2Distance.large', l2DistanceLargeCase, exact);
    family.case('ImageUtility2.getSkeleton', skeletonCase, exact);
    family.case('ImageUtility2.drawThickPixel', drawThickPixelCase, exact);
    family.case('ImageUtility2.drawLine', drawLine2Case, exact);
    family.case('ImageUtility2.drawCircle', drawCircleCase, exact);
    family.case('ImageUtility2.drawRectangle', drawRectangleCase, exact);
    family.case('ImageUtility2.drawEllipse', drawEllipseCase, exact);
    // #443: upstream never terminates for zero extents (the C++ callback
    // throws after 1000 calls); the port visits the center once.
    family.case('ImageUtility2.drawEllipse.zeroExtents', drawEllipseZeroCase,
        { exact: true, deviation: '#443 (UPSTREAM-FINDINGS, ImageUtility2 item 3)' });
    family.case('ImageUtility2.drawFloodFill4', drawFloodFillCase, exact);

    family.case('ImageUtility3.getComponents', components3Case, exact);
    family.case('ImageUtility3.dilate', (io) => morph3Case(io, 0), exact);
    family.case('ImageUtility3.erode', (io) => morph3Case(io, 1), exact);
    family.case('ImageUtility3.open', (io) => morph3Case(io, 2), exact);
    family.case('ImageUtility3.close', (io) => morph3Case(io, 3), exact);
    // #129: upstream's Dilate skips x = 0 sources; the port does not.
    family.case('ImageUtility3.morphology.xmin', (io) => {
        const op = io.integer();
        morph3Case(io, op === 0 ? 0 : op === 1 ? 2 : 3);
    }, { exact: true, deviation: '#129 (UPSTREAM-FINDINGS, ImageUtility3 item 1)' });
    family.case('ImageUtility3.computeCDConvex', cdConvexCase, exact);
    family.case('ImageUtility3.floodFill6.int', floodFill3IntCase, exact);
    family.case('ImageUtility3.floodFill6.double', floodFill3DoubleCase, exact);
    family.case('ImageUtility3.drawLine', drawLine3Case, exact);

    it('image utilities satisfy their definitions (independent checks)', () => {
        runChecks(componentsChecked, checkComponents);
        runChecks(morphChecked, checkMorph);
        runChecks(boundaryChecked, checkBoundary);
        runChecks(fillChecked, checkFill);
        runChecks(distanceChecked, checkDistance);
        runChecks(skeletonChecked, checkSkeleton);
        runChecks(drawChecked, checkDraw);
        runChecks(convexChecked, checkConvex);
        runChecks(zeroEllipses, (r) => (r.points.length === 1
            && r.points[0][0] === r.center[0] && r.points[0][1] === r.center[1] ? [] : ['zero ellipse']));
    });

    family.case('SurfaceExtractorCubes.extract',
        (io) => surfaceCase(io, 'cubes', false, -50, 50, true, true), exact);
    family.case('SurfaceExtractorCubes.extract.types',
        (io) => surfaceCase(io, 'cubes', true, 0, 0, false, true), exact);
    family.case('SurfaceExtractorCubes.extract.large',
        (io) => surfaceCase(io, 'cubes', false, -200, 200, false, true), exact);
    family.case('SurfaceExtractorCubes.extract.plusSign',
        (io) => surfaceCase(io, 'cubes', false, -30, 30, true, true), exact);
    // Upstream suspect 1 of the v26 report: the saddle-face pairing.
    family.case('SurfaceExtractorCubes.extract.saddle',
        (io) => surfaceCase(io, 'cubes', false, -20, 20, true, true),
        { exact: true, deviation: 'v26 report, SurfaceExtractorCubes saddle-face pairing' });
    family.case('SurfaceExtractorCubes.orientTriangles.points', (io) => surfacePointsCase(io, 'cubes'), exact);
    // Every record's orientation depends on dividing the gradient sum by 3
    // (upstream) rather than multiplying by 1/3.
    family.case('SurfaceExtractorCubes.orientTriangles.thirds', (io) => surfacePointsCase(io, 'cubes'), exact);
    family.case('SurfaceExtractor.invalidBounds', invalidBoundsCase, exact);

    family.case('SurfaceExtractorMC.extractVoxel', (io) => mcVoxelCase(io, true), exact);
    family.case('SurfaceExtractorMC.extractVoxel.levelTopology', (io) => mcVoxelCase(io, false), exact);
    // #443: upstream's edge interpolation omits the level.
    family.case('SurfaceExtractorMC.extractVoxel.level', (io) => mcVoxelCase(io, true),
        { exact: true, deviation: '#443 (UPSTREAM-FINDINGS, SurfaceExtractorMC item 1)' });
    family.case('SurfaceExtractorMC.extract', (io) => mcImageCase(io, false, true), exact);
    family.case('SurfaceExtractorMC.extract.large', (io) => mcImageCase(io, false, false), exact);
    family.case('SurfaceExtractorMC.extract.levelTopology', (io) => mcImageCase(io, true, true), exact);
    family.case('SurfaceExtractorMC.extract.level', (io) => mcImageCase(io, false, false),
        { exact: true, deviation: '#443 (UPSTREAM-FINDINGS, SurfaceExtractorMC item 1)' });
    family.case('SurfaceExtractorMC.makeUnique', mcMakeUniqueCase, exact);
    family.case('SurfaceExtractorMC.orientTriangles.points', mcPointsCase, exact);
    // Every record's orientation depends on multiplying the gradient sum by
    // 1/3 (upstream's Vector3 operator/) rather than dividing by 3.
    family.case('SurfaceExtractorMC.orientTriangles.thirds', mcPointsCase, exact);
    family.case('SurfaceExtractorMC.extractVoxel.signedZero', (io) => mcVoxelCase(io, true), exact);
    // Upstream suspect 2 of the v26 report: ComputeNormals uses only the
    // first triangle.
    family.case('SurfaceExtractorMC.computeNormals.multiple', mcNormalsCase,
        { exact: true, deviation: 'v26 report, SurfaceExtractorMC ComputeNormals triangle pointer' });

    family.case('SurfaceExtractorTetrahedra.extract',
        (io) => surfaceCase(io, 'tetra', false, -50, 50, true, true), exact);
    family.case('SurfaceExtractorTetrahedra.extract.types',
        (io) => surfaceCase(io, 'tetra', true, 0, 0, false, true), exact);
    family.case('SurfaceExtractorTetrahedra.extract.large',
        (io) => surfaceCase(io, 'tetra', false, -200, 200, false, true), exact);
    family.case('SurfaceExtractorTetrahedra.orientTriangles.points', (io) => surfacePointsCase(io, 'tetra'), exact);
    // #132: GetGradient's central tetrahedron of odd-parity cubes.
    family.case('SurfaceExtractorTetrahedra.orientTriangles.centralTetra', (io) => surfacePointsCase(io, 'tetra'),
        { exact: true, deviation: '#132 (UPSTREAM-FINDINGS, SurfaceExtractorTetrahedra item 3)' });

    it('the MarchingCubes table resolves every ambiguous face by its signs alone (independent check)', () => {
        const { ambiguous, bad } = checkMarchingCubesFaces();
        expect(ambiguous).toBe(192);
        expect(bad).toEqual([]);
    });

    it('extracted surfaces lie on the level set, pair saddle faces by the bilinear saddle and are closed (independent checks)', () => {
        runChecks(surfaceChecked, checkSurface);
        runChecks(mcVoxelChecked, checkMCVoxel);
        runChecks(mcImageChecked, checkMCImage);
        const statsDir = join(dirname(fileURLToPath(import.meta.url)), '..', '..', 'oracle', 'out');
        mkdirSync(statsDir, { recursive: true });
        writeFileSync(join(statsDir, `v26-stats-${surfaceChecked.length}.json`),
            JSON.stringify({ surfaceStats, skeletonStats }, null, 1));
    });

    family.finish();
});
