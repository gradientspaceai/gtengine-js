// Replays oracle/cpp/cases/v24-imaging.cpp. Keep the two files in the same
// order. Every case is arithmetic-only or integer bookkeeping and is
// compared bit for bit.
import { describe, expect, it } from 'vitest';
import { AdaptiveSkeletonClimbing2 } from '../../src/AdaptiveSkeletonClimbing2.js';
import { FastGaussianBlur1 } from '../../src/FastGaussianBlur1.js';
import { FastGaussianBlur2 } from '../../src/FastGaussianBlur2.js';
import { FastGaussianBlur3 } from '../../src/FastGaussianBlur3.js';
import { FastMarch } from '../../src/FastMarch.js';
import { Histogram } from '../../src/Histogram.js';
import { Image } from '../../src/Image.js';
import { Image2 } from '../../src/Image2.js';
import { MarchingCubes } from '../../src/MarchingCubes.js';
import { PdeFilter, PdeFilterScaleType } from '../../src/PdeFilter.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- digests (mirroring the C++ Digest) ----

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

function outString(io: OracleIO, s: string): void {
    io.outInt(s.length);
    for (let i = 0; i < s.length; i += 6) {
        let packed = 0;
        for (let j = i; j < Math.min(s.length, i + 6); ++j) {
            packed = packed * 256 + s.charCodeAt(j);
        }
        io.outInt(packed);
    }
}

// ---- FastGaussianBlur ----

type PixelArray = Float64Array | Float32Array | Int32Array | Int16Array;
type PixelCtor = new (n: number) => PixelArray;

function readPixels(io: OracleIO, ctor: PixelCtor, n: number): PixelArray {
    const a = new ctor(n);
    for (let i = 0; i < n; ++i) {
        const v = io.real();
        a[i] = v;
        // The recorded values are already of the pixel type.
        expect(a[i]).toBe(v);
    }
    return a;
}

function blur1Case(io: OracleIO, ctor: PixelCtor): void {
    const xBound = io.integer();
    const input = readPixels(io, ctor, xBound);
    const scale = io.real();
    const logBase = io.real();
    const output = new ctor(xBound);
    new FastGaussianBlur1().execute(xBound, input, output, scale, logBase);
    for (const v of output) { io.outReal(v); }
    checkConstantPreserved(input, output, logBase);
}

function blur2Case(io: OracleIO, ctor: PixelCtor): void {
    const xBound = io.integer();
    const yBound = io.integer();
    const input = readPixels(io, ctor, xBound * yBound);
    const scale = io.real();
    const logBase = io.real();
    const output = new ctor(input.length);
    new FastGaussianBlur2().execute(xBound, yBound, input, output, scale, logBase);
    for (const v of output) { io.outReal(v); }
    checkConstantPreserved(input, output, logBase);
}

function blur3Case(io: OracleIO, ctor: PixelCtor): void {
    const xBound = io.integer();
    const yBound = io.integer();
    const zBound = io.integer();
    const input = readPixels(io, ctor, xBound * yBound * zBound);
    const scale = io.real();
    const logBase = io.real();
    const output = new ctor(input.length);
    new FastGaussianBlur3().execute(xBound, yBound, zBound, input, output, scale, logBase);
    for (const v of output) { io.outReal(v); }
    checkConstantPreserved(input, output, logBase);
}

// Independent reference: every second difference of a constant image is
// exactly zero (the interpolations reproduce the constant), so the blur
// returns the image unchanged for any finite logBase.
let constantImagesChecked = 0;
function checkConstantPreserved(input: PixelArray, output: PixelArray, logBase: number): void {
    if (input.every((v) => v === input[0]) && Number.isFinite(logBase)) {
        for (let i = 0; i < input.length; ++i) {
            expect(Object.is(output[i], input[i]) || output[i] === input[i]).toBe(true);
        }
        ++constantImagesChecked;
    }
}

// ---- Image ----

function outImageShape(io: OracleIO, image: Image<number>): void {
    io.outInt(image.getNumDimensions());
    io.outInt(image.getDimensions().length);
    for (let d = 0; d < image.getNumDimensions(); ++d) {
        io.outInt(image.getDimension(d));
        io.outInt(image.getOffset(d));
        io.outInt(image.getOffsets()[d]);
    }
    io.outInt(image.getNumPixels());
    io.outInt(image.getPixels().length);
}

function imageAccessCase(io: OracleIO): void {
    const numDims = io.integer();
    const dims: number[] = [];
    for (let d = 0; d < numDims; ++d) { dims.push(io.integer()); }
    const how = io.integer();
    let image: Image<number>;
    if (how === 0) {
        image = new Image<number>(dims);
    } else {
        image = new Image<number>();
        if (how === 2) {
            const firstDims = io.integer();
            const first: number[] = [];
            for (let d = 0; d < firstDims; ++d) { first.push(io.integer()); }
            image.reconstruct(first);
            image.set(0, 7.0);
        }
        image.reconstruct(dims);
    }
    outImageShape(io, image);

    const numPixels = image.getNumPixels();
    const nd = image.getNumDimensions();
    const numQueries = io.integer();
    for (let q = 0; q < numQueries; ++q) {
        const coord: number[] = [];
        for (let d = 0; d < Math.max(nd, 1); ++d) { coord.push(io.integer()); }
        io.outInt(image.getIndex(coord));
        if (numPixels > 0) {
            const index = io.integer();
            for (const v of image.getCoordinates(index)) { io.outInt(v); }
            const w = io.integer();
            const value = io.real();
            io.outReal(image.get(w));
            image.set(w, value);
            const r = io.integer();
            io.outReal(image.getClamped(r));
            io.outReal(image.get(w));
        }
    }
}

function image2AccessCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const viaReconstruct = io.boolean();
    let image: Image2<number>;
    if (viaReconstruct) {
        image = new Image2<number>();
        image.reconstruct(3, 2);
        image.reconstruct(d0, d1);
    } else {
        image = new Image2<number>(d0, d1);
    }
    outImageShape(io, image);
    const numPixels = image.getNumPixels();
    if (numPixels === 0) {
        return;
    }
    const dim0 = image.getDimension(0), dim1 = image.getDimension(1);
    for (let y = 0; y < dim1; ++y) {
        for (let x = 0; x < dim0; ++x) {
            if ((x + y) % 2 === 0) { image.set(x, y, x + 10.0 * y + 0.5); }
            else { image.set([x, y], -x - 10.0 * y - 0.25); }
        }
    }
    const numQueries = io.integer();
    for (let q = 0; q < numQueries; ++q) {
        const x = io.integer();
        const y = io.integer();
        io.outInt(image.getIndex(x, y));
        io.outInt(image.getIndex([x, y]));
        io.outReal(image.getClamped(x, y));
        io.outReal(image.getClamped([x, y]));
        if (0 <= x && x < dim0 && 0 <= y && y < dim1) {
            io.outReal(image.get(x, y));
            io.outReal(image.get([x, y]));
        }
        const index = io.integer();
        const c = image.getCoordinates(index);
        io.outInt(c[0]);
        io.outInt(c[1]);
        io.outInt(c[0]);
        io.outInt(c[1]);
        const r = io.integer();
        io.outReal(image.getClamped(r));
    }
}

function outInts(io: OracleIO, a: readonly number[]): void {
    for (const v of a) { io.outInt(v); }
}

function outPairs(io: OracleIO, a: readonly (readonly number[])[]): void {
    for (const v of a) { io.outInt(v[0]); io.outInt(v[1]); }
}

// size_t outputs: C++ emits them as signed 64-bit integers (the main case)
// or as the unsigned value converted to double (the #64 wrap case).
function outAbsolute(io: OracleIO, a: readonly number[], asUnsigned: boolean): void {
    for (const v of a) {
        if (asUnsigned) { io.outReal(v); } else { io.outInt(v); }
    }
}

function neighborhoodCase(io: OracleIO, boundaryOnly: boolean, asUnsigned: boolean): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const image = new Image2<number>(d0, d1);
    if (!boundaryOnly) {
        outInts(io, image.getNeighborhood4());
        outInts(io, image.getNeighborhood8());
        outInts(io, image.getCorners());
        outInts(io, image.getFull());
        outPairs(io, image.getNeighborhood4Coords());
        outPairs(io, image.getNeighborhood8Coords());
        outPairs(io, image.getCornersCoords());
        outPairs(io, image.getFullCoords());
    }
    const count = io.integer();
    for (let k = 0; k < count; ++k) {
        const x = io.integer();
        const y = io.integer();
        outAbsolute(io, image.getNeighborhood4(x, y), asUnsigned);
        outAbsolute(io, image.getNeighborhood8(x, y), asUnsigned);
        outAbsolute(io, image.getCorners(x, y), asUnsigned);
        outAbsolute(io, image.getFull(x, y), asUnsigned);
        for (const t of [image.getNeighborhood4Coords(x, y), image.getNeighborhood8Coords(x, y),
            image.getCornersCoords(x, y), image.getFullCoords(x, y)]) {
            outAbsolute(io, t.flat(), asUnsigned);
        }
    }
}

// ---- Histogram ----

function outHistogram(io: OracleIO, h: Histogram): void {
    const buckets = h.getBuckets();
    io.outInt(buckets.length);
    outInts(io, buckets);
    io.outInt(h.getExcessLess());
    io.outInt(h.getExcessGreater());
    for (let k = 0; k < 3; ++k) {
        const tail = io.real();
        io.outInt(h.getLowerTail(tail));
        io.outInt(h.getUpperTail(tail));
        const tails = h.getTails(tail);
        io.outInt(tails.lower);
        io.outInt(tails.upper);
    }
}
function histogramIntCase(io: OracleIO, noRescaling: boolean): void {
    const numBuckets = io.integer();
    const n = io.integer();
    const samples: number[] = [];
    for (let i = 0; i < n; ++i) { samples.push(io.integer()); }
    const h = Histogram.fromIntegerSamples(numBuckets, samples, noRescaling);
    expect(histogramTotalsAfter(h)).toBe(n);
    outHistogram(io, h);
}

// Independent reference: the buckets and the excess counts account for every
// sample exactly once.
function histogramTotalsAfter(h: Histogram): number {
    return h.getBuckets().reduce((a, b) => a + b, 0) + h.getExcessLess() + h.getExcessGreater();
}

function histogramRealCase(io: OracleIO): void {
    const numBuckets = io.integer();
    const n = io.integer();
    const samples: number[] = [];
    // One by one through io.real(): a NaN sample keeps its bits.
    for (let i = 0; i < n; ++i) { samples.push(io.real()); }
    const h = Histogram.fromRealSamples(numBuckets, samples);
    expect(histogramTotalsAfter(h)).toBe(n);
    outHistogram(io, h);
}

function histogramIncrementalCase(io: OracleIO): void {
    const numBuckets = io.integer();
    const h = new Histogram(numBuckets);
    const numOps = io.integer();
    for (let k = 0; k < numOps; ++k) {
        const checked = io.boolean();
        if (checked) {
            h.insertCheck(io.integer());
        } else {
            h.insert(io.integer());
        }
    }
    expect(histogramTotalsAfter(h)).toBe(numOps);
    outHistogram(io, h);
}

function histogramInvalidCase(io: OracleIO): void {
    const which = io.integer();
    const numBuckets = io.integer();
    const n = io.integer();
    const ones = new Array<number>(n).fill(1);
    let h: Histogram;
    switch (which) {
        case 0: h = Histogram.fromIntegerSamples(numBuckets, ones, true); break;
        case 1: h = Histogram.fromIntegerSamples(numBuckets, ones, false); break;
        case 2: h = Histogram.fromRealSamples(numBuckets, ones); break;
        case 3: h = Histogram.fromRealSamples(numBuckets, ones); break;
        default: h = new Histogram(numBuckets); break;
    }
    io.outInt(h.getBuckets().length);
}

// ---- PdeFilter ----

// The subclass of the C++ case file: exposes the base's protected state and
// logs the order of the three update hooks.
class PdeProbe extends PdeFilter {
    readonly calls: number[] = [];

    constructor(quantity: number, data: ArrayLike<number>, borderValue: number,
        scaleType: PdeFilterScaleType) {
        super(quantity, data, borderValue, scaleType);
    }

    min(): number { return this.mMin; }
    offset(): number { return this.mOffset; }
    scale(): number { return this.mScale; }

    protected onPreUpdate(): void { this.calls.push(1); }
    protected onUpdate(): void { this.calls.push(2); }
    protected onPostUpdate(): void { this.calls.push(3); }
}

function pdeFilterCase(io: OracleIO): void {
    const quantity = io.integer();
    const data = new Float64Array(quantity);
    for (let i = 0; i < quantity; ++i) { data[i] = io.real(); }
    const borderValue = io.real();
    const scaleType = io.integer() as PdeFilterScaleType;
    const filter = new PdeProbe(quantity, data, borderValue, scaleType);
    io.outInt(filter.getQuantity());
    io.outReal(filter.getBorderValue());
    io.outInt(filter.getScaleType());
    io.outReal(filter.min());
    io.outReal(filter.offset());
    io.outReal(filter.scale());
    io.outReal(filter.getTimeStep());
    filter.setTimeStep(io.real());
    io.outReal(filter.getTimeStep());
    const numUpdates = io.integer();
    for (let k = 0; k < numUpdates; ++k) { filter.update(); }
    io.outInt(filter.calls.length);
    outInts(io, filter.calls);
}

// ---- FastMarch ----

// The 1-D subclass of the C++ case file, line for line.
class FastMarch1 extends FastMarch {
    removed = false;
    removedKey = 0;
    removedValue = 0;

    constructor(quantity: number, seeds: readonly number[], speeds: readonly number[] | number) {
        super(quantity, seeds, speeds);
        this.initialize();
    }

    getBoundary(): number[] {
        const boundary: number[] = [];
        for (let i = 0; i < this.mQuantity; ++i) {
            if (this.isBoundary(i)) {
                boundary.push(i);
            }
        }
        return boundary;
    }

    isBoundary(i: number): boolean {
        return this.isValid(i) && !this.isTrial(i)
            && ((i > 0 && this.isTrial(i - 1)) || (i + 1 < this.mQuantity && this.isTrial(i + 1)));
    }

    iterate(): void {
        this.removed = false;
        if (this.mHeap.getNumElements() === 0) {
            return;
        }
        const minimum = this.mHeap.remove()!;
        const i = minimum.key;
        this.removed = true;
        this.removedKey = i;
        this.removedValue = minimum.value;
        this.mTrials[i] = null;
        if (i > 0) { this.visit(i - 1); }
        if (i + 1 < this.mQuantity) { this.visit(i + 1); }
    }

    numTrials(): number { return this.mHeap.getNumElements(); }
    invSpeed(i: number): number { return this.mInvSpeeds[i]; }

    private initialize(): void {
        for (let i = 0; i < this.mQuantity; ++i) {
            if (this.isFar(i)
                && ((i > 0 && this.isValid(i - 1) && !this.isTrial(i - 1))
                    || (i + 1 < this.mQuantity && this.isValid(i + 1) && !this.isTrial(i + 1)))) {
                this.computeTime(i);
                this.mTrials[i] = this.mHeap.insert(i, this.mTimes[i]);
            }
        }
    }

    private visit(j: number): void {
        if (this.isTrial(j)) {
            this.computeTime(j);
            this.mHeap.update(this.mTrials[j], this.mTimes[j]);
        } else if (this.isFar(j)) {
            this.computeTime(j);
            this.mTrials[j] = this.mHeap.insert(j, this.mTimes[j]);
        }
    }

    private computeTime(i: number): void {
        let t = Number.MAX_VALUE;
        if (i > 0 && this.isValid(i - 1) && this.mTimes[i - 1] < t) { t = this.mTimes[i - 1]; }
        if (i + 1 < this.mQuantity && this.isValid(i + 1) && this.mTimes[i + 1] < t) {
            t = this.mTimes[i + 1];
        }
        this.mTimes[i] = t + this.mInvSpeeds[i];
    }
}

function outFastMarchState(io: OracleIO, march: FastMarch1): void {
    const q = march.getQuantity();
    io.outInt(march.numTrials());
    for (let i = 0; i < q; ++i) {
        io.outReal(march.getTime(i));
        io.outInt((march.isValid(i) ? 1 : 0) | (march.isTrial(i) ? 2 : 0)
            | (march.isFar(i) ? 4 : 0) | (march.isZeroSpeed(i) ? 8 : 0)
            | (march.isInterior(i) ? 16 : 0) | (march.isBoundary(i) ? 32 : 0));
    }
    const interior = march.getInterior();
    const boundary = march.getBoundary();
    io.outInt(interior.length);
    outInts(io, interior);
    io.outInt(boundary.length);
    outInts(io, boundary);
    const extremes = march.getTimeExtremes();
    io.outReal(extremes.minValue);
    io.outReal(extremes.maxValue);
}

function drawFastMarch(io: OracleIO): FastMarch1 {
    const quantity = io.integer();
    const numSeeds = io.integer();
    const seeds: number[] = [];
    for (let k = 0; k < numSeeds; ++k) { seeds.push(io.integer()); }
    const perPixel = io.boolean();
    if (perPixel) {
        const speeds: number[] = [];
        for (let i = 0; i < quantity; ++i) { speeds.push(io.real()); }
        return new FastMarch1(quantity, seeds, speeds);
    }
    return new FastMarch1(quantity, seeds, io.real());
}

function fastMarchCase(io: OracleIO): void {
    const march = drawFastMarch(io);
    const q = march.getQuantity();
    for (let i = 0; i < q; ++i) { io.outReal(march.invSpeed(i)); }
    outFastMarchState(io, march);
    const numIterations = io.integer();
    for (let k = 0; k < numIterations; ++k) {
        march.iterate();
        io.outBool(march.removed);
        if (march.removed) {
            io.outInt(march.removedKey);
            io.outReal(march.removedValue);
        }
    }
    outFastMarchState(io, march);
}

function fastMarchAccessorsCase(io: OracleIO): void {
    const march = drawFastMarch(io);
    const numSets = io.integer();
    for (let k = 0; k < numSets; ++k) {
        const i = io.integer();
        const t = io.real();
        march.setTime(i, t);
    }
    outFastMarchState(io, march);
}

// ---- MarchingCubes ----

function marchingCubesTableCase(io: OracleIO): void {
    const first = io.integer();
    const mc = new MarchingCubes();
    const flat = mc.getFlatTable();
    const prebuilt = MarchingCubes.getPrebuiltTable();
    for (let entry = first; entry < first + 13; ++entry) {
        if (entry < 256) {
            const topology = mc.getTable(entry);
            io.outInt(topology.numVertices);
            io.outInt(topology.numTriangles);
            for (const p of topology.vpair) { io.outInt(p[0] * 64 + p[1]); }
            for (const t of topology.itriple) { io.outInt((t[0] * 64 + t[1]) * 64 + t[2]); }
            const flatDigest = new Digest(), prebuiltDigest = new Digest();
            for (let k = 0; k < 41; ++k) {
                flatDigest.integer(flat[41 * entry + k]);
                prebuiltDigest.integer(prebuilt[entry][k]);
            }
            io.outInt(flatDigest.h);
            io.outInt(prebuiltDigest.h);
        }
        outString(io, MarchingCubes.getConfigurationType(entry));
    }
    outString(io, MarchingCubes.getConfigurationType(io.integer()));
}

// ---- AdaptiveSkeletonClimbing2 ----

type Vertex2 = [number, number];
type Edge2 = [number, number];

interface ASCResult {
    v: Vertex2[];
    e: Edge2[];
    vu: Vertex2[];
    eu: Edge2[];
}

// Extractions whose outputs the independent checks below examine.
const ascChecked: { pixels: number[], size: number, level: number, r: ASCResult }[] = [];

function runASC(N: number, pixels: number[], levels: number[], depths: number[]): ASCResult[] {
    const asc = new AdaptiveSkeletonClimbing2(N, pixels);
    const results: ASCResult[] = [];
    for (let k = 0; k < levels.length; ++k) {
        const { vertices, edges } = asc.extract(levels[k], depths[k]);
        const vu = vertices.map((p) => [p[0], p[1]] as Vertex2);
        const eu = edges.map((s) => [s[0], s[1]] as Edge2);
        asc.makeUnique(vu, eu);
        const r = { v: vertices, e: edges, vu, eu };
        results.push(r);
        ascChecked.push({ pixels, size: (1 << N) + 1, level: levels[k], r });
    }
    return results;
}

function outCurves(io: OracleIO, v: readonly Vertex2[], e: readonly Edge2[]): void {
    io.outInt(v.length);
    for (const p of v) { io.outReal(p[0]); io.outReal(p[1]); }
    io.outInt(e.length);
    for (const s of e) { io.outInt(s[0]); io.outInt(s[1]); }
}

function outASCResults(io: OracleIO, results: readonly ASCResult[]): void {
    for (const r of results) {
        outCurves(io, r.v, r.e);
        outCurves(io, r.vu, r.eu);
    }
}

// The record layout of the main, types, saddle and saddlePairing cases:
// N, the number of extractions, the pixels, then (level, depth) pairs.
function ascCase(io: OracleIO): void {
    const N = io.integer();
    const numExtracts = io.integer();
    const size = (1 << N) + 1;
    const pixels: number[] = [];
    for (let i = 0; i < size * size; ++i) { pixels.push(io.integer()); }
    const levels: number[] = [];
    const depths: number[] = [];
    for (let k = 0; k < numExtracts; ++k) {
        levels.push(io.real());
        depths.push(io.integer());
    }
    outASCResults(io, runASC(N, pixels, levels, depths));
}

function largePixel(kind: number, p: readonly number[], x: number, y: number): number {
    const u = x - p[5], v = y - p[6];
    if (kind === 0) {
        return p[0] * u * u + p[1] * v * v + p[2] * u * v + p[3] * x + p[4] * y;
    }
    const au = Math.abs(u), av = Math.abs(v);
    return p[0] * au + p[1] * av + p[2] * Math.min(au, av) + p[3] * x + p[4] * y;
}

function digestCurves(io: OracleIO, v: readonly Vertex2[], e: readonly Edge2[]): void {
    const dv = new Digest(), de = new Digest();
    for (const p of v) { dv.real(p[0]); dv.real(p[1]); }
    for (const s of e) { de.integer(s[0]); de.integer(s[1]); }
    io.outInt(v.length);
    io.outInt(e.length);
    io.outInt(dv.h);
    io.outInt(de.h);
}

function ascLargeCase(io: OracleIO): void {
    const N = io.integer();
    const size = (1 << N) + 1;
    const kind = io.integer();
    const p: number[] = [];
    for (let k = 0; k < 7; ++k) { p.push(io.integer()); }
    const level = io.real();
    const depth = io.integer();
    const pixels: number[] = [];
    for (let y = 0; y < size; ++y) {
        for (let x = 0; x < size; ++x) { pixels.push(largePixel(kind, p, x, y)); }
    }
    const [r] = runASC(N, pixels, [level], [depth]);
    digestCurves(io, r.v, r.e);
    digestCurves(io, r.vu, r.eu);
}

function ascInvalidCase(io: OracleIO): void {
    const N = io.integer();
    const size = (1 << N) + 1;
    const pixels: number[] = [];
    for (let i = 0; i < size * size; ++i) { pixels.push(io.integer()); }
    outASCResults(io, runASC(N, pixels, [0.5], [-1]));
}

// Independent checks of an extraction. Every vertex lies on a grid edge
// (x or y an integer) where the edge's linear interpolant equals the level
// to rounding, except the plus-sign branch points; after MakeUnique every
// vertex on the image border has degree 1 and every interior one degree 2
// (a branch point 4). Returns the violations found.
function checkASC(pixels: readonly number[], size: number, level: number, r: ASCResult): string[] {
    const bad: string[] = [];
    const f = (x: number, y: number): number => pixels[x + size * y];
    for (const [x, y] of r.vu) {
        let value: number | undefined;
        if (Number.isInteger(x)) {
            const y0 = Math.min(Math.floor(y), size - 2);
            value = f(x, y0) + (y - y0) * (f(x, y0 + 1) - f(x, y0));
        } else if (Number.isInteger(y)) {
            const x0 = Math.min(Math.floor(x), size - 2);
            value = f(x0, y) + (x - x0) * (f(x0 + 1, y) - f(x0, y));
        }
        if (value !== undefined) {
            const scale = Math.max(1, Math.abs(level), ...pixels.map(Math.abs));
            if (Math.abs(value - level) > 1e-12 * scale) {
                bad.push(`vertex (${x}, ${y}) interpolates ${value}, level ${level}`);
            }
        }
    }
    const degree = new Array<number>(r.vu.length).fill(0);
    for (const [a, b] of r.eu) { ++degree[a]; ++degree[b]; }
    r.vu.forEach(([x, y], i) => {
        const onBorder = x === 0 || y === 0 || x === size - 1 || y === size - 1;
        const onGrid = Number.isInteger(x) || Number.isInteger(y);
        const expected = onBorder ? [1] : (onGrid ? [2] : [4]);
        if (!expected.includes(degree[i])) {
            bad.push(`vertex (${x}, ${y}) has degree ${degree[i]}`);
        }
    });
    return bad;
}

// Independent checks of the MarchingCubes table. Corner c of the voxel is
// (c & 1, (c >> 1) & 1, (c >> 2) & 1); bit c of the entry set means the
// corner value is negative. The vertex pairs must be exactly the voxel
// edges whose corners differ in sign, and every triangle edge shared by two
// triangles must be traversed in opposite directions (a consistently
// oriented surface). Returns the violations found.
function checkMarchingCubesEntry(entry: number, t: ReturnType<MarchingCubes['getTable']>): string[] {
    const bad: string[] = [];
    const cut = new Set<string>();
    for (let a = 0; a < 8; ++a) {
        for (const bit of [1, 2, 4]) {
            const b = a | bit;
            if (b !== a && (((entry >> a) & 1) !== ((entry >> b) & 1))) { cut.add(`${a},${b}`); }
        }
    }
    const pairs = t.vpair.slice(0, t.numVertices).map((p) => `${p[0]},${p[1]}`);
    if (new Set(pairs).size !== pairs.length || pairs.length !== cut.size
        || !pairs.every((p) => cut.has(p))) {
        bad.push(`entry ${entry}: vertex pairs ${pairs.join(' ')} are not the cut edges`);
    }
    const directed = new Map<string, number>();
    for (const tri of t.itriple.slice(0, t.numTriangles)) {
        for (let k = 0; k < 3; ++k) {
            const a = tri[k], b = tri[(k + 1) % 3];
            if (a >= t.numVertices || b >= t.numVertices) { bad.push(`entry ${entry}: index`); }
            const key = `${a},${b}`;
            directed.set(key, (directed.get(key) ?? 0) + 1);
        }
    }
    for (const [key, count] of directed) {
        if (count > 1) { bad.push(`entry ${entry}: directed edge ${key} used ${count} times`); }
    }
    return bad;
}

describe('oracle: v24-imaging', () => {
    const family = new OracleFamily('v24-imaging');
    const exact = { exact: true };

    family.case('FastGaussianBlur1.execute.double', (io) => blur1Case(io, Float64Array), exact);
    family.case('FastGaussianBlur1.execute.float', (io) => blur1Case(io, Float32Array), exact);
    family.case('FastGaussianBlur1.execute.int32', (io) => blur1Case(io, Int32Array), exact);
    family.case('FastGaussianBlur1.execute.int16', (io) => blur1Case(io, Int16Array), exact);
    family.case('FastGaussianBlur2.execute.double', (io) => blur2Case(io, Float64Array), exact);
    family.case('FastGaussianBlur2.execute.float', (io) => blur2Case(io, Float32Array), exact);
    family.case('FastGaussianBlur2.execute.int32', (io) => blur2Case(io, Int32Array), exact);
    family.case('FastGaussianBlur2.execute.int16', (io) => blur2Case(io, Int16Array), exact);
    family.case('FastGaussianBlur3.execute.double', (io) => blur3Case(io, Float64Array), exact);
    family.case('FastGaussianBlur3.execute.float', (io) => blur3Case(io, Float32Array), exact);
    family.case('FastGaussianBlur3.execute.int32', (io) => blur3Case(io, Int32Array), exact);
    family.case('FastGaussianBlur3.execute.int16', (io) => blur3Case(io, Int16Array), exact);

    it('blurs constant images to themselves (independent check)', () => {
        expect(constantImagesChecked).toBeGreaterThan(0);
    });

    family.case('Image.access', imageAccessCase, exact);
    family.case('Image2.access', image2AccessCase, exact);
    family.case('Image2.neighborhoods', (io) => neighborhoodCase(io, false, false), exact);
    family.case('Image2.neighborhoods.wrap', (io) => neighborhoodCase(io, true, true),
        { exact: true, deviation: '#64 (UPSTREAM-FINDINGS, Image2.h/Image3.h size_t wrap)' });

    family.case('Histogram.int.noRescaling', (io) => histogramIntCase(io, true), exact);
    family.case('Histogram.int.rescaled', (io) => histogramIntCase(io, false), exact);
    family.case('Histogram.double', histogramRealCase, exact);
    family.case('Histogram.float', histogramRealCase, exact);
    family.case('Histogram.incremental', histogramIncrementalCase, exact);
    family.case('Histogram.invalid', histogramInvalidCase, exact);

    family.case('PdeFilter.construct', pdeFilterCase, exact);

    family.case('FastMarch.march1', fastMarchCase, exact);
    family.case('FastMarch.march1.nanSpeed', fastMarchCase, exact);
    family.case('FastMarch.accessors', fastMarchAccessorsCase, exact);

    family.case('MarchingCubes.table', marchingCubesTableCase, exact);

    it('MarchingCubes table: cut edges, consistent orientation, prebuilt table (independent check)', () => {
        const mc = new MarchingCubes();
        const flat = mc.getFlatTable();
        const prebuilt = MarchingCubes.getPrebuiltTable();
        const bad: string[] = [];
        for (let entry = 0; entry < 256; ++entry) {
            bad.push(...checkMarchingCubesEntry(entry, mc.getTable(entry)));
            if (prebuilt[entry].some((v, k) => v !== flat[41 * entry + k])) {
                bad.push(`entry ${entry}: prebuilt row differs from the constructed one`);
            }
        }
        expect(bad).toEqual([]);
    });

    family.case('AdaptiveSkeletonClimbing2.extract', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing2.extract.types', (io) => {
        io.integer();  // the pixel type; the port has one number type
        ascCase(io);
    }, exact);
    family.case('AdaptiveSkeletonClimbing2.extract.saddle', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing2.extract.saddlePairing', ascCase, exact);
    family.case('AdaptiveSkeletonClimbing2.extract.large', ascLargeCase, exact);
    family.case('AdaptiveSkeletonClimbing2.invalid', ascInvalidCase, exact);

    it('AdaptiveSkeletonClimbing2 curves lie on the level set and are closed or end on the border (independent check)', () => {
        const bad: string[] = [];
        for (const c of ascChecked) {
            bad.push(...checkASC(c.pixels, c.size, c.level, c.r));
        }
        expect(ascChecked.length).toBeGreaterThan(0);
        expect(bad.slice(0, 10)).toEqual([]);
    });

    family.finish();
});
