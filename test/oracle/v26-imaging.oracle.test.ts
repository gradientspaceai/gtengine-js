// Replays oracle/cpp/cases/v26-imaging.cpp. Keep the two files in the same
// order. Every case is integer bookkeeping or IEEE-exact arithmetic (sqrt
// and sqrtf included) and is compared bit for bit. The independent checks
// at the end of each section run on the port's outputs of the replayed
// records.
import { describe, expect, it } from 'vitest';
import { Image2 } from '../../src/Image2.js';
import { Image3 } from '../../src/Image3.js';
import { ImageUtility2 } from '../../src/ImageUtility2.js';
import { ImageUtility3 } from '../../src/ImageUtility3.js';
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

    // @@TESTS@@

    family.finish();
});
