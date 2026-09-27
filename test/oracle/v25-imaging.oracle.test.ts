// Replays oracle/cpp/cases/v25-imaging.cpp. Keep the two files in the same
// order. Every case is arithmetic-only or integer bookkeeping and is
// compared bit for bit, except the updated cells of the GradientAnisotropic
// filters (std::exp; default tolerance 1e-12, everything else exact).
import { describe, expect, it } from 'vitest';
import { Array2 } from '../../src/Array2.js';
import { Array3 } from '../../src/Array3.js';
import { CurvatureFlow2 } from '../../src/CurvatureFlow2.js';
import { CurvatureFlow3 } from '../../src/CurvatureFlow3.js';
import { FastMarch } from '../../src/FastMarch.js';
import { FastMarch2 } from '../../src/FastMarch2.js';
import { FastMarch3 } from '../../src/FastMarch3.js';
import { GaussianBlur2 } from '../../src/GaussianBlur2.js';
import { GaussianBlur3 } from '../../src/GaussianBlur3.js';
import { GradientAnisotropic2 } from '../../src/GradientAnisotropic2.js';
import { GradientAnisotropic3 } from '../../src/GradientAnisotropic3.js';
import { Image } from '../../src/Image.js';
import { Image3 } from '../../src/Image3.js';
import { MinHeap, type MinHeapRecord } from '../../src/MinHeap.js';
import { PdeFilter, PdeFilterScaleType } from '../../src/PdeFilter.js';
import { PdeFilter1 } from '../../src/PdeFilter1.js';
import { PdeFilter2 } from '../../src/PdeFilter2.js';
import { PdeFilter3 } from '../../src/PdeFilter3.js';
import {
    SurfaceExtractor, SurfaceExtractorTriangle, SurfaceExtractorVertex
} from '../../src/SurfaceExtractor.js';
import { TetrahedraRasterizer } from '../../src/TetrahedraRasterizer.js';
import { OracleFamily, type OracleIO } from './harness.js';

// ---- digests (mirroring the C++ Digest; NaN canonicalized) ----

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
        if (Number.isNaN(x)) {
            this.word(0);
            this.word(0x7FF80000);
            return;
        }
        const view = new DataView(new ArrayBuffer(8));
        view.setFloat64(0, x, true);
        this.word(view.getUint32(0, true));
        this.word(view.getUint32(4, true));
    }

    out(io: OracleIO): void { io.outInt(this.h); }
}

// ---- Image3 ----

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

function image3AccessCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    const viaReconstruct = io.boolean();
    let image: Image3<number>;
    if (viaReconstruct) {
        image = new Image3<number>();
        image.reconstruct(2, 3, 2);
        image.reconstruct(d0, d1, d2);
    } else {
        image = new Image3<number>(d0, d1, d2);
    }
    outImageShape(io, image);
    const numPixels = image.getNumPixels();
    if (numPixels === 0) {
        return;
    }
    const dim0 = image.getDimension(0), dim1 = image.getDimension(1);
    const dim2 = image.getDimension(2);
    for (let z = 0; z < dim2; ++z) {
        for (let y = 0; y < dim1; ++y) {
            for (let x = 0; x < dim0; ++x) {
                const v = x + 10.0 * y + 100.0 * z;
                if ((x + y + z) % 2 === 0) { image.set(x, y, z, v + 0.5); }
                else { image.set([x, y, z], -v - 0.25); }
            }
        }
    }
    const numQueries = io.integer();
    for (let q = 0; q < numQueries; ++q) {
        const x = io.integer();
        const y = io.integer();
        const z = io.integer();
        io.outInt(image.getIndex(x, y, z));
        io.outInt(image.getIndex([x, y, z]));
        io.outReal(image.getClamped(x, y, z));
        io.outReal(image.getClamped([x, y, z]));
        if (0 <= x && x < dim0 && 0 <= y && y < dim1 && 0 <= z && z < dim2) {
            io.outReal(image.get(x, y, z));
            io.outReal(image.get([x, y, z]));
        }
        const index = io.integer();
        const c = image.getCoordinates(index);
        io.outInt(c[0]);
        io.outInt(c[1]);
        io.outInt(c[2]);
        io.outInt(c[0]);
        io.outInt(c[1]);
        io.outInt(c[2]);
        const r = io.integer();
        io.outReal(image.getClamped(r));
    }
}

function digestTable(io: OracleIO, a: readonly number[]): void {
    const d = new Digest();
    for (const v of a) { d.integer(v); }
    d.out(io);
}

function digestTuples(io: OracleIO, a: readonly (readonly number[])[]): void {
    const d = new Digest();
    for (const v of a) { d.integer(v[0]); d.integer(v[1]); d.integer(v[2]); }
    d.out(io);
}

function image3NeighborhoodCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    const image = new Image3<number>(d0, d1, d2);
    for (const v of image.getNeighborhood6()) { io.outInt(v); }
    digestTable(io, image.getNeighborhood18());
    digestTable(io, image.getNeighborhood26());
    digestTable(io, image.getCorners8());
    digestTable(io, image.getFull27());
    digestTuples(io, image.getNeighborhood6Coords());
    digestTuples(io, image.getNeighborhood18Coords());
    digestTuples(io, image.getNeighborhood26Coords());
    digestTuples(io, image.getCorners8Coords());
    digestTuples(io, image.getFull27Coords());
    const count = io.integer();
    for (let k = 0; k < count; ++k) {
        const x = io.integer();
        const y = io.integer();
        const z = io.integer();
        for (const v of image.getNeighborhood6(x, y, z)) { io.outInt(v); }
        digestTable(io, image.getNeighborhood18(x, y, z));
        digestTable(io, image.getNeighborhood26(x, y, z));
        digestTable(io, image.getCorners8(x, y, z));
        digestTable(io, image.getFull27(x, y, z));
        digestTuples(io, image.getNeighborhood6Coords(x, y, z));
        digestTuples(io, image.getNeighborhood18Coords(x, y, z));
        digestTuples(io, image.getNeighborhood26Coords(x, y, z));
        digestTuples(io, image.getCorners8Coords(x, y, z));
        digestTuples(io, image.getFull27Coords(x, y, z));
    }
}

// C++ emits the unsigned (wrapped, #64) values; the port's are negative.
function image3WrapCase(io: OracleIO): void {
    const d0 = io.integer();
    const d1 = io.integer();
    const d2 = io.integer();
    const image = new Image3<number>(d0, d1, d2);
    const x = io.integer();
    const y = io.integer();
    const z = io.integer();
    for (const v of image.getNeighborhood6(x, y, z)) { io.outReal(v); }
    for (const v of image.getNeighborhood6Coords(x, y, z)) {
        io.outReal(v[0]);
        io.outReal(v[1]);
        io.outReal(v[2]);
    }
}

// ---- PDE filters ----

// Independent-reference bookkeeping (checked in the last tests).
const refs = {
    blurSumChecked: 0,
    blurSumStaleAfterStep2: 0,
    blurSumStep2Checked: 0,
    flowPlanarChecked: 0,
    flowTvChecked: 0,
    flowTvIncreased: 0
};

function readData(io: OracleIO, n: number): Float64Array {
    const data = new Float64Array(n);
    for (let i = 0; i < n; ++i) { data[i] = io.real(); }
    return data;
}

function readMask(io: OracleIO, n: number): Int32Array | null {
    const mode = io.integer();
    if (mode === 0) {
        return null;
    }
    const mask = new Int32Array(n);
    for (let i = 0; i < n; ++i) { mask[i] = io.integer(); }
    return mask;
}

// The subclass of the C++ case file: the explicit heat equation through
// PdeFilter1's lookUp3.
class Heat1 extends PdeFilter1 {
    constructor(xBound: number, xSpacing: number, data: ArrayLike<number>,
        mask: ArrayLike<number> | null, borderValue: number, scaleType: PdeFilterScaleType) {
        super(xBound, xSpacing, data, mask, borderValue, scaleType);
    }

    min(): number { return this.mMin; }
    offset(): number { return this.mOffset; }
    scale(): number { return this.mScale; }
    src(): number { return this.mSrc; }
    buf(k: number, x: number): number { return this.mBuffer[k][x]; }
    maskAt(x: number): number { return this.mMask[x]; }

    protected onUpdateSingle(x: number): void {
        this.lookUp3(x);
        this.mBuffer[this.mDst][x] = this.mUz + this.mTimeStep * this.mInvDxDx
            * (this.mUp - 2.0 * this.mUz + this.mUm);
    }
}

function outState1(io: OracleIO, f: Heat1): void {
    const xB = f.getXBound();
    io.outInt(f.src());
    for (let k = 0; k < 2; ++k) {
        const d = new Digest();
        for (let x = 0; x < xB + 2; ++x) { d.real(f.buf(k, x)); }
        d.out(io);
    }
    const d = new Digest();
    for (let x = 0; x < xB; ++x) {
        io.outReal(f.getU(x));
        d.real(f.getUx(x));
        d.real(f.getUxx(x));
        d.integer(f.getMask(x));
    }
    for (let x = 0; x < xB + 2; ++x) { d.integer(f.maskAt(x)); }
    d.out(io);
}

function pdeFilter1Case(io: OracleIO): void {
    const xB = io.integer();
    const dx = io.real();
    const data = readData(io, xB);
    const mask = readMask(io, xB);
    const borderValue = io.real();
    const scaleType = io.integer() as PdeFilterScaleType;
    const f = new Heat1(xB, dx, data, mask, borderValue, scaleType);
    f.setTimeStep(io.real());
    io.outInt(f.getQuantity());
    io.outInt(f.getXBound());
    io.outReal(f.getXSpacing());
    io.outReal(f.getBorderValue());
    io.outInt(f.getScaleType());
    io.outReal(f.min());
    io.outReal(f.offset());
    io.outReal(f.scale());
    io.outReal(f.getTimeStep());
    outState1(io, f);
    const numSteps = io.integer();
    for (let k = 0; k < numSteps; ++k) {
        f.update();
        outState1(io, f);
    }
}

// The protected state of the port's filters, read the way the C++ probes
// read upstream's.
interface PdeState {
    mMin: number;
    mOffset: number;
    mScale: number;
    mSrc: number;
    mHasMask: boolean;
    mParameter: number;
    mMHalfParameter: number;
}

interface Pde2State extends PdeState {
    mBuffer: [Array2<number>, Array2<number>];
    mMask: Array2<number>;
}

interface Pde3State extends PdeState {
    mBuffer: [Array3<number>, Array3<number>];
    mMask: Array3<number>;
}

function state2(f: PdeFilter2): Pde2State { return f as unknown as Pde2State; }
function state3(f: PdeFilter3): Pde3State { return f as unknown as Pde3State; }

function outConstantsCommon(io: OracleIO, f: PdeFilter, s: PdeState): void {
    io.outReal(f.getBorderValue());
    io.outInt(f.getScaleType());
    io.outReal(s.mMin);
    io.outReal(s.mOffset);
    io.outReal(s.mScale);
    io.outReal(f.getTimeStep());
}

function digestDerived2(io: OracleIO, f: PdeFilter2): void {
    const xB = f.getXBound(), yB = f.getYBound();
    const d = new Digest();
    for (let y = 0; y < yB; ++y) {
        for (let x = 0; x < xB; ++x) {
            d.real(f.getU(x, y));
            d.real(f.getUx(x, y));
            d.real(f.getUy(x, y));
            d.real(f.getUxx(x, y));
            d.real(f.getUxy(x, y));
            d.real(f.getUyy(x, y));
            d.integer(f.getMask(x, y));
        }
    }
    const s = state2(f);
    for (let y = 0; y < yB + 2; ++y) {
        for (let x = 0; x < xB + 2; ++x) { d.integer(s.mMask.get(x, y)); }
    }
    d.out(io);
}

function outState2(io: OracleIO, f: PdeFilter2): void {
    const xB = f.getXBound(), yB = f.getYBound();
    const s = state2(f);
    io.outInt(s.mSrc);
    for (let k = 0; k < 2; ++k) {
        const d = new Digest();
        for (let y = 0; y < yB + 2; ++y) {
            for (let x = 0; x < xB + 2; ++x) { d.real(s.mBuffer[k].get(x, y)); }
        }
        d.out(io);
    }
    digestDerived2(io, f);
    if (xB * yB <= 9) {
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) { io.outReal(f.getU(x, y)); }
        }
    }
}

// 'exact' holds the construction state (no exp evaluated yet) to bit
// identity.
function outStateAniso2(io: OracleIO, f: GradientAnisotropic2, exact: boolean): void {
    const xB = f.getXBound(), yB = f.getYBound();
    const s = state2(f);
    io.outInt(s.mSrc);
    const d = new Digest();
    for (let k = 0; k < 2; ++k) {
        for (let y = 0; y < yB + 2; ++y) {
            for (let x = 0; x < xB + 2; ++x) {
                const written = 1 <= x && x <= xB && 1 <= y && y <= yB
                    && (!s.mHasMask || s.mMask.get(x, y) !== 0);
                const v = s.mBuffer[k].get(x, y);
                if (!written) { d.real(v); }
                else if (k === s.mSrc) {
                    if (exact) { io.outRealExact(v); } else { io.outReal(v); }
                }
            }
        }
    }
    d.out(io);
    if (exact) {
        io.outRealExact(s.mParameter);
        io.outRealExact(s.mMHalfParameter);
    } else {
        io.outReal(s.mParameter);
        io.outReal(s.mMHalfParameter);
    }
}

type Kind = 'blur' | 'flow' | 'aniso';

// Interior values of the source image, x fastest.
function interior2(f: PdeFilter2): number[] {
    const u: number[] = [];
    for (let y = 0; y < f.getYBound(); ++y) {
        for (let x = 0; x < f.getXBound(); ++x) { u.push(f.getU(x, y)); }
    }
    return u;
}

function interior3(f: PdeFilter3): number[] {
    const u: number[] = [];
    for (let z = 0; z < f.getZBound(); ++z) {
        for (let y = 0; y < f.getYBound(); ++y) {
            for (let x = 0; x < f.getXBound(); ++x) { u.push(f.getU(x, y, z)); }
        }
    }
    return u;
}

// True when the data depend on x only (the planar-level-set mode).
function dependsOnXOnly(data: Float64Array, xB: number): boolean {
    for (let i = 0; i < data.length; ++i) {
        if (data[i] !== data[i % xB]) { return false; }
    }
    return true;
}

function sumAbs(u: readonly number[]): number {
    let s = 0;
    for (const v of u) { s += Math.abs(v); }
    return s;
}

// Total variation of the interior (forward differences along each axis).
function totalVariation(u: readonly number[], dims: readonly number[]): number {
    let tv = 0;
    const [xB, yB, zB] = [dims[0], dims[1] ?? 1, dims[2] ?? 1];
    for (let z = 0; z < zB; ++z) {
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                const i = x + xB * (y + yB * z);
                if (x + 1 < xB) { tv += Math.abs(u[i + 1] - u[i]); }
                if (y + 1 < yB) { tv += Math.abs(u[i + xB] - u[i]); }
                if (z + 1 < zB) { tv += Math.abs(u[i + xB * yB] - u[i]); }
            }
        }
    }
    return tv;
}

// The independent references common to 2D and 3D (see the final tests):
// Neumann Gaussian blur without a mask preserves the interior sum on the
// first step (the ghost cells mirror the boundary, the second differences
// telescope to zero) and, because the ghost cells are never refreshed
// (#60 item 2), not necessarily afterwards; a curvature flow leaves an
// image whose level sets are planes (data depending on x only) unchanged;
// with a small step the curvature flow does not increase the total
// variation (counted, not asserted).
class References {
    private previous: number[];
    private readonly neumannPlain: boolean;
    private readonly planar: boolean;
    private readonly smallStep: boolean;

    constructor(private readonly kind: Kind, private readonly dims: readonly number[],
        u0: number[], border: number, mask: Int32Array | null, data: Float64Array,
        dt: number, minSpacing: number) {
        this.previous = u0;
        this.neumannPlain = border === Number.MAX_VALUE && mask === null
            && u0.every((v) => Number.isFinite(v));
        this.planar = kind === 'flow' && this.neumannPlain && dependsOnXOnly(data, dims[0]);
        this.smallStep = dt > 0 && dt <= 0.125 * minSpacing * minSpacing;
    }

    step(k: number, u: number[], dt: number, invSqrSum: number): void {
        if (this.kind === 'blur' && this.neumannPlain && u.every((v) => Number.isFinite(v))) {
            const s0 = this.previous.reduce((a, b) => a + b, 0);
            const s1 = u.reduce((a, b) => a + b, 0);
            const scale = Math.max(1, sumAbs(this.previous), sumAbs(u))
                * (1 + 4 * Math.abs(dt) * invSqrSum);
            const agrees = Math.abs(s1 - s0) <= 1e-12 * scale;
            if (k === 1) {
                expect(agrees, `blur step 1 changed the Neumann sum from ${s0} to ${s1}`).toBe(true);
                ++refs.blurSumChecked;
            } else if (k === 2) {
                ++refs.blurSumStep2Checked;
                if (!agrees) { ++refs.blurSumStaleAfterStep2; }
            }
        }
        if (this.planar) {
            for (let i = 0; i < u.length; ++i) {
                expect(u[i] === this.previous[i] || (Number.isNaN(u[i])
                    && Number.isNaN(this.previous[i])),
                    `curvature flow changed a planar image at ${i}`).toBe(true);
            }
            ++refs.flowPlanarChecked;
        }
        if (this.kind === 'flow' && this.neumannPlain && this.smallStep && k === 1
            && u.every((v) => Number.isFinite(v))) {
            const tv0 = totalVariation(this.previous, this.dims);
            const tv1 = totalVariation(u, this.dims);
            ++refs.flowTvChecked;
            if (tv1 > tv0 * (1 + 1e-12) + 1e-300) { ++refs.flowTvIncreased; }
        }
        this.previous = u;
    }
}

function filter2Case(io: OracleIO, kind: Kind): void {
    const xB = io.integer();
    const yB = io.integer();
    const dx = io.real();
    const dy = io.real();
    const data = readData(io, xB * yB);
    const mask = readMask(io, xB * yB);
    const border = io.real();
    const scaleType = io.integer() as PdeFilterScaleType;
    let f: PdeFilter2;
    if (kind === 'aniso') {
        f = new GradientAnisotropic2(xB, yB, dx, dy, data, mask, border, scaleType, io.real());
    } else if (kind === 'blur') {
        const blur = new GaussianBlur2(xB, yB, dx, dy, data, mask, border, scaleType);
        io.outReal(blur.getMaximumTimeStep());
        f = blur;
    } else {
        f = new CurvatureFlow2(xB, yB, dx, dy, data, mask, border, scaleType);
    }
    const dt = io.real();
    f.setTimeStep(dt);
    io.outInt(f.getQuantity());
    io.outInt(f.getXBound());
    io.outInt(f.getYBound());
    io.outReal(f.getXSpacing());
    io.outReal(f.getYSpacing());
    outConstantsCommon(io, f, state2(f));
    digestDerived2(io, f);
    const references = new References(kind, [xB, yB], interior2(f), border, mask, data, dt,
        Math.min(dx, dy));
    const numSteps = io.integer();
    for (let k = 0; k <= numSteps; ++k) {
        if (k > 0) {
            f.update();
            references.step(k, interior2(f), dt, 1 / (dx * dx) + 1 / (dy * dy));
        }
        if (kind === 'aniso') { outStateAniso2(io, f as GradientAnisotropic2, k === 0); }
        else { outState2(io, f); }
    }
}

function digestDerived3(io: OracleIO, f: PdeFilter3): void {
    const xB = f.getXBound(), yB = f.getYBound(), zB = f.getZBound();
    const d = new Digest();
    for (let z = 0; z < zB; ++z) {
        for (let y = 0; y < yB; ++y) {
            for (let x = 0; x < xB; ++x) {
                d.real(f.getU(x, y, z));
                d.real(f.getUx(x, y, z));
                d.real(f.getUy(x, y, z));
                d.real(f.getUz(x, y, z));
                d.real(f.getUxx(x, y, z));
                d.real(f.getUxy(x, y, z));
                d.real(f.getUxz(x, y, z));
                d.real(f.getUyy(x, y, z));
                d.real(f.getUyz(x, y, z));
                d.real(f.getUzz(x, y, z));
                d.integer(f.getMask(x, y, z));
            }
        }
    }
    const s = state3(f);
    for (let z = 0; z < zB + 2; ++z) {
        for (let y = 0; y < yB + 2; ++y) {
            for (let x = 0; x < xB + 2; ++x) { d.integer(s.mMask.get(x, y, z)); }
        }
    }
    d.out(io);
}

function outState3(io: OracleIO, f: PdeFilter3): void {
    const xB = f.getXBound(), yB = f.getYBound(), zB = f.getZBound();
    const s = state3(f);
    io.outInt(s.mSrc);
    for (let k = 0; k < 2; ++k) {
        const d = new Digest();
        for (let z = 0; z < zB + 2; ++z) {
            for (let y = 0; y < yB + 2; ++y) {
                for (let x = 0; x < xB + 2; ++x) { d.real(s.mBuffer[k].get(x, y, z)); }
            }
        }
        d.out(io);
    }
    digestDerived3(io, f);
    if (xB * yB * zB <= 8) {
        for (const v of interior3(f)) { io.outReal(v); }
    }
}

function outStateAniso3(io: OracleIO, f: GradientAnisotropic3, exact: boolean): void {
    const xB = f.getXBound(), yB = f.getYBound(), zB = f.getZBound();
    const s = state3(f);
    io.outInt(s.mSrc);
    const d = new Digest();
    for (let k = 0; k < 2; ++k) {
        for (let z = 0; z < zB + 2; ++z) {
            for (let y = 0; y < yB + 2; ++y) {
                for (let x = 0; x < xB + 2; ++x) {
                    const written = 1 <= x && x <= xB && 1 <= y && y <= yB && 1 <= z && z <= zB
                        && (!s.mHasMask || s.mMask.get(x, y, z) !== 0);
                    const v = s.mBuffer[k].get(x, y, z);
                    if (!written) { d.real(v); }
                    else if (k === s.mSrc) {
                        if (exact) { io.outRealExact(v); } else { io.outReal(v); }
                    }
                }
            }
        }
    }
    d.out(io);
    if (exact) {
        io.outRealExact(s.mParameter);
        io.outRealExact(s.mMHalfParameter);
    } else {
        io.outReal(s.mParameter);
        io.outReal(s.mMHalfParameter);
    }
}

function filter3Case(io: OracleIO, kind: Kind): void {
    const xB = io.integer();
    const yB = io.integer();
    const zB = io.integer();
    const dx = io.real();
    const dy = io.real();
    const dz = io.real();
    const n = xB * yB * zB;
    const data = readData(io, n);
    const mask = readMask(io, n);
    const border = io.real();
    const scaleType = io.integer() as PdeFilterScaleType;
    let f: PdeFilter3;
    if (kind === 'aniso') {
        f = new GradientAnisotropic3(xB, yB, zB, dx, dy, dz, data, mask, border, scaleType,
            io.real());
    } else if (kind === 'blur') {
        const blur = new GaussianBlur3(xB, yB, zB, dx, dy, dz, data, mask, border, scaleType);
        io.outReal(blur.getMaximumTimeStep());
        f = blur;
    } else {
        f = new CurvatureFlow3(xB, yB, zB, dx, dy, dz, data, mask, border, scaleType);
    }
    const dt = io.real();
    f.setTimeStep(dt);
    io.outInt(f.getQuantity());
    io.outInt(f.getXBound());
    io.outInt(f.getYBound());
    io.outInt(f.getZBound());
    io.outReal(f.getXSpacing());
    io.outReal(f.getYSpacing());
    io.outReal(f.getZSpacing());
    outConstantsCommon(io, f, state3(f));
    digestDerived3(io, f);
    const references = new References(kind, [xB, yB, zB], interior3(f), border, mask, data, dt,
        Math.min(dx, dy, dz));
    const numSteps = io.integer();
    for (let k = 0; k <= numSteps; ++k) {
        if (k > 0) {
            f.update();
            references.step(k, interior3(f), dt, 1 / (dx * dx) + 1 / (dy * dy) + 1 / (dz * dz));
        }
        if (kind === 'aniso') { outStateAniso3(io, f as GradientAnisotropic3, k === 0); }
        else { outState3(io, f); }
    }
}

// ---- FastMarch2, FastMarch3 ----

interface MarchState {
    mHeap: MinHeap<number, number>;
    mTrials: (MinHeapRecord<number, number> | null)[];
    mInvSpeeds: number[];
}

function marchState(m: FastMarch): MarchState { return m as unknown as MarchState; }

// ComputeTime branch histogram (which upwind terms exist, the sign of the
// discriminant), recomputed from the state before each call.
const branches = { one: 0, twoPos: 0, twoNeg: 0, threePos: 0, threeNeg: 0 };
const marchRefs = { removals: 0, decreases: 0 };

class CountingMarch2 extends FastMarch2 {
    protected override computeTime(i: number): void {
        const xb = this.mXBound;
        const t = this.mTimes, inv = this.mInvSpeeds;
        const term = (a: number, b: number): number | null => {
            if (this.isValid(a)) { return this.isValid(b) && t[b] < t[a] ? t[b] : t[a]; }
            return this.isValid(b) ? t[b] : null;
        };
        const cx = term(i - 1, i + 1), cy = term(i - xb, i + xb);
        if (cx !== null && cy !== null) {
            const diff = cx - cy;
            if (2 * inv[i] * inv[i] - diff * diff >= 0) { ++branches.twoPos; } else { ++branches.twoNeg; }
        } else if (cx !== null || cy !== null) {
            ++branches.one;
        }
        super.computeTime(i);
    }
}

class CountingMarch3 extends FastMarch3 {
    protected override computeTime(i: number): void {
        const t = this.mTimes, inv = this.mInvSpeeds;
        const term = (a: number, b: number): number | null => {
            if (this.isValid(a)) { return this.isValid(b) && t[b] < t[a] ? t[b] : t[a]; }
            return this.isValid(b) ? t[b] : null;
        };
        const c = [term(i - 1, i + 1), term(i - this.mXBound, i + this.mXBound),
            term(i - this.mXYBound, i + this.mXYBound)].filter((v): v is number => v !== null);
        if (c.length === 3) {
            let discr = 3 * inv[i] * inv[i];
            discr -= (c[0] - c[1]) * (c[0] - c[1]);
            discr -= (c[0] - c[2]) * (c[0] - c[2]);
            discr -= (c[1] - c[2]) * (c[1] - c[2]);
            if (discr >= 0) { ++branches.threePos; } else { ++branches.threeNeg; }
        } else if (c.length === 2) {
            const diff = c[0] - c[1];
            if (2 * inv[i] * inv[i] - diff * diff >= 0) { ++branches.twoPos; } else { ++branches.twoNeg; }
        } else if (c.length === 1) {
            ++branches.one;
        }
        super.computeTime(i);
    }
}

function flags(m: FastMarch, i: number): number {
    return (m.isValid(i) ? 1 : 0) | (m.isTrial(i) ? 2 : 0) | (m.isFar(i) ? 4 : 0)
        | (m.isZeroSpeed(i) ? 8 : 0) | (m.isInterior(i) ? 16 : 0) | (m.isBoundary(i) ? 32 : 0);
}

function outMarchState(io: OracleIO, m: FastMarch, include: readonly boolean[]): void {
    const s = marchState(m);
    const q = m.getQuantity();
    io.outInt(s.mHeap.getNumElements());
    const times = new Digest(), fl = new Digest(), speeds = new Digest(), trials = new Digest();
    let numTrials = 0;
    for (let i = 0; i < q; ++i) {
        if (!include[i]) { continue; }
        times.real(m.getTime(i));
        fl.integer(flags(m, i));
        speeds.real(s.mInvSpeeds[i]);
        const record = s.mTrials[i];
        if (record !== null) {
            ++numTrials;
            trials.integer(record.key);
            trials.real(record.value);
            trials.integer(record.index);
        }
    }
    times.out(io);
    fl.out(io);
    speeds.out(io);
    io.outInt(numTrials);
    trials.out(io);
    const interior = m.getInterior();
    const boundary = m.getBoundary();
    io.outInt(interior.length);
    digestTable(io, interior);
    io.outInt(boundary.length);
    digestTable(io, boundary);
    const e = m.getTimeExtremes();
    io.outReal(e.minValue);
    io.outReal(e.maxValue);
}

function march(io: OracleIO, m: FastMarch, maxIter: number, stop: (key: number) => boolean,
    count: boolean): void {
    let last = -Infinity;
    for (let k = 0; k < maxIter; ++k) {
        const heap = marchState(m).mHeap;
        if (heap.getNumElements() === 0) { break; }
        const minimum = heap.getMinimum();
        if (minimum === null || stop(minimum.key)) { break; }
        io.outInt(minimum.key);
        io.outReal(minimum.value);
        if (count && Number.isFinite(minimum.value)) {
            ++marchRefs.removals;
            if (minimum.value < last) { ++marchRefs.decreases; }
            last = minimum.value;
        }
        m.iterate();
    }
    io.outInt(-1);
}

function drawMarch2(io: OracleIO): FastMarch2 {
    const xB = io.integer();
    const yB = io.integer();
    const dx = io.real();
    const dy = io.real();
    const numSeeds = io.integer();
    const seeds: number[] = [];
    for (let k = 0; k < numSeeds; ++k) { seeds.push(io.integer()); }
    const speedMode = io.integer();
    if (speedMode === 0) {
        return new CountingMarch2(xB, yB, dx, dy, seeds, io.real());
    }
    const speeds: number[] = [];
    for (let i = 0; i < xB * yB; ++i) { speeds.push(io.real()); }
    return new CountingMarch2(xB, yB, dx, dy, seeds, speeds);
}

function outMarch2Constants(io: OracleIO, m: FastMarch2): void {
    io.outInt(m.getQuantity());
    io.outInt(m.getXBound());
    io.outInt(m.getYBound());
    io.outReal(m.getXSpacing());
    io.outReal(m.getYSpacing());
    io.outInt(m.index(m.getXBound() - 1, m.getYBound() - 1));
}

function fastMarch2Case(io: OracleIO): void {
    const m = drawMarch2(io);
    const all = new Array<boolean>(m.getQuantity()).fill(true);
    outMarch2Constants(io, m);
    outMarchState(io, m, all);
    march(io, m, io.integer(), () => false, true);
    outMarchState(io, m, all);
}

function fastMarch2AccessorsCase(io: OracleIO): void {
    const m = drawMarch2(io);
    const numSets = io.integer();
    for (let k = 0; k < numSets; ++k) {
        const x = io.integer();
        const y = io.integer();
        m.setTime(m.index(x, y), io.real());
    }
    const all = new Array<boolean>(m.getQuantity()).fill(true);
    outMarchState(io, m, all);
    march(io, m, io.integer(), () => false, false);
    outMarchState(io, m, all);
}

function numExtremes(i: number, xB: number, yB: number, zB: number): number {
    const x = i % xB;
    const y = Math.floor(i / xB) % yB;
    const z = Math.floor(i / (xB * yB));
    return (x === 0 || x === xB - 1 ? 1 : 0) + (y === 0 || y === yB - 1 ? 1 : 0)
        + (z === 0 || z === zB - 1 ? 1 : 0);
}

function outMarch3Constants(io: OracleIO, m: FastMarch3): void {
    io.outInt(m.getQuantity());
    io.outInt(m.getXBound());
    io.outInt(m.getYBound());
    io.outInt(m.getZBound());
    io.outReal(m.getXSpacing());
    io.outReal(m.getYSpacing());
    io.outReal(m.getZSpacing());
    io.outInt(m.index(m.getXBound() - 1, m.getYBound() - 1, m.getZBound() - 1));
}

function readSeeds3(io: OracleIO, xB: number, yB: number, interior: boolean): number[] {
    const numSeeds = io.integer();
    const seeds: number[] = [];
    for (let k = 0; k < numSeeds; ++k) {
        if (interior) {
            const x = io.integer();
            const y = io.integer();
            const z = io.integer();
            seeds.push(x + xB * (y + yB * z));
        } else {
            seeds.push(io.integer());
        }
    }
    return seeds;
}

function fastMarch3Case(io: OracleIO): void {
    const xB = io.integer(), yB = io.integer(), zB = io.integer();
    const dx = io.real(), dy = io.real(), dz = io.real();
    const seeds = readSeeds3(io, xB, yB, false);
    io.boolean();  // the speed mix (recorded for the generator only)
    const speeds: number[] = [];
    for (let i = 0; i < xB * yB * zB; ++i) { speeds.push(io.real()); }
    const m = new CountingMarch3(xB, yB, zB, dx, dy, dz, seeds, speeds);
    const all = new Array<boolean>(m.getQuantity()).fill(true);
    outMarch3Constants(io, m);
    outMarchState(io, m, all);
    march(io, m, io.integer(), () => false, true);
    outMarchState(io, m, all);
}

function fastMarch3ConstantCase(io: OracleIO): void {
    const xB = io.integer(), yB = io.integer(), zB = io.integer();
    const dx = io.real(), dy = io.real(), dz = io.real();
    const seeds = readSeeds3(io, xB, yB, true);
    const m = new CountingMarch3(xB, yB, zB, dx, dy, dz, seeds, io.real());
    const q = m.getQuantity();
    const include: boolean[] = [];
    for (let i = 0; i < q; ++i) { include.push(numExtremes(i, xB, yB, zB) !== 1); }
    outMarch3Constants(io, m);
    outMarchState(io, m, include);
    const hasFaceNeighbor = (key: number): boolean =>
        [key - 1, key + 1, key - xB, key + xB, key - xB * yB, key + xB * yB]
            .some((j) => numExtremes(j, xB, yB, zB) === 1);
    march(io, m, io.integer(), hasFaceNeighbor, true);
    outMarchState(io, m, include);
}

function fastMarch3FacesCase(io: OracleIO): void {
    const xB = io.integer(), yB = io.integer(), zB = io.integer();
    const seeds = readSeeds3(io, xB, yB, true);
    const speeds: number[] = [];
    for (let i = 0; i < xB * yB * zB; ++i) { speeds.push(io.real()); }
    const m = new FastMarch3(xB, yB, zB, 1, 1, 1, seeds, speeds);
    const all = new Array<boolean>(m.getQuantity()).fill(true);
    outMarchState(io, m, all);
    march(io, m, 2 * m.getQuantity(), () => false, false);
    outMarchState(io, m, all);
}

// ---- SurfaceExtractor ----

// The subclass of the C++ case file: recorded rational vertices and
// triangles, and the affine gradient (g0 + g1 x, g2 + g3 y, g4 + g5 z).
class Recorded extends SurfaceExtractor {
    rv: SurfaceExtractorVertex[] = [];
    rt: SurfaceExtractorTriangle[] = [];
    g = [0, 0, 0, 0, 0, 0];

    constructor(xB: number, yB: number, zB: number, voxels: ArrayLike<number>) {
        super(xB, yB, zB, voxels);
    }

    extractRational(): { vertices: SurfaceExtractorVertex[], triangles: SurfaceExtractorTriangle[] } {
        return {
            vertices: this.rv.map((v) => Object.assign(new SurfaceExtractorVertex(), v)),
            triangles: this.rt.map((t) => Object.assign(new SurfaceExtractorTriangle(),
                { v: [...t.v] as [number, number, number] }))
        };
    }

    protected getGradient(p: [number, number, number]): [number, number, number] {
        const g = this.g;
        return [g[0] + g[1] * p[0], g[2] + g[3] * p[1], g[4] + g[5] * p[2]];
    }
}

const voxels8 = new Int32Array(8);

function drawSurface(io: OracleIO, r: Recorded): void {
    const numVertices = io.integer();
    for (let v = 0; v < numVertices; ++v) {
        const n: number[] = [];
        for (let c = 0; c < 3; ++c) {
            n.push(io.integer());
            n.push(io.integer());
        }
        r.rv.push(new SurfaceExtractorVertex(n[0], n[1], n[2], n[3], n[4], n[5]));
    }
    const numTriangles = io.integer();
    for (let t = 0; t < numTriangles; ++t) {
        const v0 = io.integer();
        const v1 = io.integer();
        const v2 = io.integer();
        r.rt.push(new SurfaceExtractorTriangle(v0, v1, v2));
    }
}

function outRationals(io: OracleIO, vertices: readonly SurfaceExtractorVertex[]): void {
    io.outInt(vertices.length);
    for (const v of vertices) {
        io.outInt(v.xNumer); io.outInt(v.xDenom);
        io.outInt(v.yNumer); io.outInt(v.yDenom);
        io.outInt(v.zNumer); io.outInt(v.zDenom);
    }
}

function outTriangles(io: OracleIO, triangles: readonly SurfaceExtractorTriangle[]): void {
    io.outInt(triangles.length);
    for (const t of triangles) { io.outInt(t.v[0]); io.outInt(t.v[1]); io.outInt(t.v[2]); }
}

function outPoints(io: OracleIO, points: readonly (readonly number[])[]): void {
    io.outInt(points.length);
    for (const p of points) { io.outReal(p[0]); io.outReal(p[1]); io.outReal(p[2]); }
}

function vertexTriangleCase(io: OracleIO): void {
    const r = new Recorded(2, 2, 2, voxels8);
    drawSurface(io, r);
    outRationals(io, r.rv);
    for (let k = 0; k + 1 < r.rv.length; ++k) {
        io.outBool(r.rv[k].equals(r.rv[k + 1]));
        io.outBool(r.rv[k].lessThan(r.rv[k + 1]));
        io.outBool(r.rv[k + 1].lessThan(r.rv[k]));
    }
    outTriangles(io, r.rt);
    for (let k = 0; k + 1 < r.rt.length; ++k) {
        io.outBool(r.rt[k].equals(r.rt[k + 1]));
        io.outBool(r.rt[k].lessThan(r.rt[k + 1]));
        io.outBool(r.rt[k + 1].lessThan(r.rt[k]));
    }
}

// Independent reference: after MakeUnique no two vertices are equal as
// rationals, every triangle refers to the same rational points as the
// input triangle it came from, and the triples are pairwise distinct.
let makeUniqueChecked = 0;
let makeUniqueRotatedDuplicates = 0;
function makeUniqueCase(io: OracleIO): void {
    const r = new Recorded(2, 2, 2, voxels8);
    drawSurface(io, r);
    const { vertices: v, triangles: t } = r.extractRational();
    r.makeUnique(v, t);
    outRationals(io, v);
    outTriangles(io, t);
    outPoints(io, r.convert(v));
    if (t.length === 0) {
        return;  // MakeUnique returns at once when there are no triangles
    }
    for (let a = 0; a < v.length; ++a) {
        for (let b = a + 1; b < v.length; ++b) { expect(v[a].equals(v[b])).toBe(false); }
    }
    const keys = new Set(t.map((x) => x.v.join(',')));
    expect(keys.size).toBe(t.length);
    for (const x of t) {
        const source = r.rt.find((y) => [0, 1, 2].every((i) => r.rv[y.v[i]].equals(v[x.v[i]])));
        expect(source, 'a MakeUnique triangle matches an input triangle').toBeDefined();
    }
    // #439 (preserved): the remapped triples are not re-rotated, so two
    // rotations of one triangle can both survive.
    const rotation = (x: SurfaceExtractorTriangle): string =>
        new SurfaceExtractorTriangle(x.v[0], x.v[1], x.v[2]).v.join(',');
    if (new Set(t.map(rotation)).size < t.length) { ++makeUniqueRotatedDuplicates; }
    ++makeUniqueChecked;
}

function extractCase(io: OracleIO): void {
    const r = new Recorded(2, 2, 2, voxels8);
    drawSurface(io, r);
    const unique = io.boolean();
    for (let i = 0; i < 6; ++i) { r.g[i] = io.real(); }
    const sameDir = io.boolean();
    const { vertices, triangles } = r.extract(0, unique);
    outPoints(io, vertices);
    outTriangles(io, triangles);
    r.orientTriangles(vertices, triangles, sameDir);
    outTriangles(io, triangles);
    outPoints(io, r.computeNormals(vertices, triangles));
}

function surfaceInvalidCase(io: OracleIO): void {
    const xB = io.integer();
    const yB = io.integer();
    const zB = io.integer();
    new Recorded(xB, yB, zB, new Int32Array(27));
    io.outInt(1);
}

// ---- TetrahedraRasterizer ----

function outGrid(io: OracleIO, grid: Int32Array): void {
    io.outInt(grid.length);
    const d = new Digest();
    let count = 0;
    for (let j = 0; j < grid.length; ++j) {
        if (grid[j] !== -1) { ++count; d.integer(j); d.integer(grid[j]); }
    }
    io.outInt(count);
    d.out(io);
    if (count <= 40) {
        for (let j = 0; j < grid.length; ++j) {
            if (grid[j] !== -1) { io.outInt(j); io.outInt(grid[j]); }
        }
    }
}

interface Mesh { vertices: number[][]; tetra: number[][]; }

function readMesh(io: OracleIO): Mesh {
    const nv = io.integer();
    const vertices: number[][] = [];
    for (let v = 0; v < nv; ++v) { vertices.push([io.real(), io.real(), io.real()]); }
    const nt = io.integer();
    const tetra: number[][] = [];
    for (let t = 0; t < nt; ++t) {
        tetra.push([io.integer(), io.integer(), io.integer(), io.integer()]);
    }
    return { vertices, tetra };
}

interface Region { rmin: number[]; rmax: number[]; bound: number[]; }

function readRegion(io: OracleIO): Region {
    const rmin = [io.real(), io.real(), io.real()];
    const rmax = [io.real(), io.real(), io.real()];
    const bound = [io.integer(), io.integer(), io.integer()];
    return { rmin, rmax, bound };
}

// Independent reference for the lattice case (grid multiplier exactly 1):
// the voxel (i0, i1, i2) gets the last tetrahedron whose clipped box
// contains it and for which the four signed volumes of upstream's
// PointInTetrahedron are all <= 0, evaluated in BigInt.
let tetraChecked = 0;
function rasterizeReference(mesh: Mesh, region: Region): Int32Array {
    const { rmin, rmax, bound } = region;
    const grid = new Int32Array(bound[0] * bound[1] * bound[2]).fill(-1);
    const sub = (u: bigint[], v: bigint[]): bigint[] => [u[0] - v[0], u[1] - v[1], u[2] - v[2]];
    const dotCross = (u: bigint[], v: bigint[], w: bigint[]): bigint =>
        u[0] * (v[1] * w[2] - v[2] * w[1]) + u[1] * (v[2] * w[0] - v[0] * w[2])
        + u[2] * (v[0] * w[1] - v[1] * w[0]);
    mesh.tetra.forEach((tet, t) => {
        const V = tet.map((k) => mesh.vertices[k].map((c, i) => BigInt(c - rmin[i])));
        const lo: number[] = [], hi: number[] = [];
        for (let i = 0; i < 3; ++i) {
            const cs = tet.map((k) => mesh.vertices[k][i]);
            lo.push(Math.max(Math.min(...cs), rmin[i]) - rmin[i]);
            hi.push(Math.min(Math.max(...cs), rmax[i]) - rmin[i]);
        }
        if (lo.some((v, i) => v > hi[i])) { return; }
        for (let i2 = Math.ceil(lo[2]); i2 <= Math.floor(hi[2]); ++i2) {
            for (let i1 = Math.ceil(lo[1]); i1 <= Math.floor(hi[1]); ++i1) {
                for (let i0 = Math.ceil(lo[0]); i0 <= Math.floor(hi[0]); ++i0) {
                    const P = [BigInt(i0), BigInt(i1), BigInt(i2)];
                    const p0 = sub(P, V[0]), e1 = sub(V[1], V[0]), e2 = sub(V[2], V[0]);
                    const e3 = sub(V[3], V[0]);
                    const inside = dotCross(p0, e2, e1) <= 0n && dotCross(p0, e1, e3) <= 0n
                        && dotCross(p0, e3, e2) <= 0n
                        && dotCross(sub(P, V[1]), sub(V[2], V[1]), sub(V[3], V[1])) <= 0n;
                    if (inside) { grid[i0 + bound[0] * (i1 + bound[1] * i2)] = t; }
                }
            }
        }
    });
    return grid;
}

function tetraCase(io: OracleIO, reference: boolean): void {
    const mesh = readMesh(io);
    const region = readRegion(io);
    io.integer();  // numThreads: 0 or 1 upstream; the port has one code path
    const r = new TetrahedraRasterizer(mesh.vertices, mesh.tetra);
    const grid = r.rasterize(region.rmin, region.rmax, region.bound);
    outGrid(io, grid);
    if (reference) {
        expect([...grid]).toEqual([...rasterizeReference(mesh, region)]);
        ++tetraChecked;
    }
}

function tetraTwiceCase(io: OracleIO): void {
    const mesh = readMesh(io);
    const region0 = readRegion(io);
    const region1 = readRegion(io);
    const r = new TetrahedraRasterizer(mesh.vertices, mesh.tetra);
    outGrid(io, r.rasterize(region0.rmin, region0.rmax, region0.bound));
    outGrid(io, r.rasterize(region1.rmin, region1.rmax, region1.bound));
}

function tetraInvalidCase(io: OracleIO): void {
    const which = io.integer();
    const vertices = which === 0 ? [] : [[0, 0, 0], [2, 0, 0], [0, 2, 0], [0, 0, 2]];
    const tetra = which === 1 ? [] : [[0, 1, 2, 3]];
    const bound = [3, 3, 3];
    if (which === 2) {
        const axis = io.integer();
        bound[axis] = io.integer();
    }
    const r = new TetrahedraRasterizer(vertices, tetra);
    outGrid(io, r.rasterize([0, 0, 0], [2, 2, 2], bound));
}

// ---- registration ----

describe('oracle: v25-imaging', () => {
    const family = new OracleFamily('v25-imaging');
    const exact = { exact: true };

    family.case('Image3.access', image3AccessCase, exact);
    family.case('Image3.neighborhoods', image3NeighborhoodCase, exact);
    family.case('Image3.neighborhoods.wrap', image3WrapCase,
        { exact: true, deviation: '#64 (size_t neighbour coordinates wrap at the xmin/ymin/zmin faces)' });

    family.case('PdeFilter1.heat', pdeFilter1Case, exact);
    family.case('GaussianBlur2.update', (io) => filter2Case(io, 'blur'), exact);
    family.case('CurvatureFlow2.update', (io) => filter2Case(io, 'flow'), exact);
    // std::exp in the conductances: default tolerance 1e-12 on the updated
    // cells after the first step; outRealExact holds the rest.
    family.case('GradientAnisotropic2.update', (io) => filter2Case(io, 'aniso'));
    family.case('GaussianBlur3.update', (io) => filter3Case(io, 'blur'), exact);
    family.case('CurvatureFlow3.update', (io) => filter3Case(io, 'flow'), exact);
    family.case('GradientAnisotropic3.update', (io) => filter3Case(io, 'aniso'));

    family.case('FastMarch2.march', fastMarch2Case, exact);
    family.case('FastMarch2.accessors', fastMarch2AccessorsCase, exact);
    family.case('FastMarch3.march', fastMarch3Case, exact);
    family.case('FastMarch3.march.constantSpeed', fastMarch3ConstantCase, exact);
    family.case('FastMarch3.faces', fastMarch3FacesCase,
        { exact: true, deviation: '#121 (upstream leaves the six boundary faces unmarked)' });

    family.case('SurfaceExtractor.vertexTriangle', vertexTriangleCase, exact);
    family.case('SurfaceExtractor.makeUnique', makeUniqueCase, exact);
    family.case('SurfaceExtractor.extract', extractCase, exact);
    family.case('SurfaceExtractor.invalid', surfaceInvalidCase, exact);

    family.case('TetrahedraRasterizer.lattice', (io) => tetraCase(io, true), exact);
    family.case('TetrahedraRasterizer.general', (io) => tetraCase(io, false), exact);
    family.case('TetrahedraRasterizer.twice', tetraTwiceCase, exact);
    family.case('TetrahedraRasterizer.threads', (io) => tetraCase(io, false), exact);
    family.case('TetrahedraRasterizer.invalid', tetraInvalidCase, exact);

    family.finish();

    it('independent references were exercised', () => {
        console.log('v25 references', JSON.stringify({ refs, branches, marchRefs,
            makeUniqueChecked, makeUniqueRotatedDuplicates, tetraChecked }));
        expect(refs.blurSumChecked).toBeGreaterThan(0);
        expect(refs.flowPlanarChecked).toBeGreaterThan(0);
        expect(makeUniqueChecked).toBeGreaterThan(0);
        expect(tetraChecked).toBeGreaterThan(0);
        expect(branches.one).toBeGreaterThan(0);
        expect(branches.twoPos).toBeGreaterThan(0);
        expect(branches.twoNeg).toBeGreaterThan(0);
        expect(branches.threePos).toBeGreaterThan(0);
    });
});
