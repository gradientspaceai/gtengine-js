import { describe, it, expect } from 'vitest';
import { STLBinaryFile, STLTriangle, type STLTuple3 } from '../src/STLBinaryFile.js';
import { check, fc, finite } from './helpers/arbitraries.js';

function makeTriangle(normal: STLTuple3, v0: STLTuple3, v1: STLTuple3, v2: STLTuple3,
    attributeByteCount: number = 0): STLTriangle {
    const triangle = new STLTriangle();
    triangle.normal = normal;
    triangle.vertex = [v0, v1, v2];
    triangle.attributeByteCount = attributeByteCount;
    return triangle;
}

describe('STLTriangle', () => {
    it('zero-fills on construction', () => {
        const triangle = new STLTriangle();
        expect(triangle.normal).toEqual([0, 0, 0]);
        expect(triangle.vertex).toEqual([[0, 0, 0], [0, 0, 0], [0, 0, 0]]);
        expect(triangle.attributeByteCount).toBe(0);
    });
});

describe('STLBinaryFile', () => {
    it('constructs with a zero-filled header and no triangles', () => {
        const file = new STLBinaryFile();
        expect(file.header.length).toBe(80);
        expect(Array.from(file.header)).toEqual(new Array(80).fill(0));
        expect(file.triangles).toEqual([]);
    });

    it('saves with the exact binary STL layout', () => {
        const file = new STLBinaryFile();
        const headerText = 'gtengine-js unit test';
        for (let i = 0; i < headerText.length; ++i) {
            file.header[i] = headerText.charCodeAt(i);
        }
        file.triangles.push(makeTriangle(
            [0, 0, 1], [0, 0, 0], [1, 0, 0], [0, 1, 0], 0xABCD));
        file.triangles.push(makeTriangle(
            [0, 0, -1], [0, 0, 0], [0, 1, 0], [1, 0, 0]));

        const buffer = file.save();
        expect(buffer.byteLength).toBe(80 + 4 + 2 * 50);

        const bytes = new Uint8Array(buffer);
        const view = new DataView(buffer);

        // Header occupies bytes [0,80): the text then zero padding.
        for (let i = 0; i < 80; ++i) {
            const expected = i < headerText.length ? headerText.charCodeAt(i) : 0;
            expect(bytes[i]).toBe(expected);
        }

        // Triangle count is a little-endian uint32 at byte 80.
        expect(view.getUint32(80, true)).toBe(2);
        expect(Array.from(bytes.subarray(80, 84))).toEqual([2, 0, 0, 0]);

        // Triangle 0 occupies bytes [84,134): normal, three vertices, then
        // the attribute byte count, all little-endian.
        let offset = 84;
        const floats0 = [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0];
        for (const value of floats0) {
            expect(view.getFloat32(offset, true)).toBe(value);
            offset += 4;
        }
        expect(view.getUint16(offset, true)).toBe(0xABCD);
        expect(Array.from(bytes.subarray(offset, offset + 2))).toEqual([0xCD, 0xAB]);
        offset += 2;
        expect(offset).toBe(134);

        // Triangle 1 occupies bytes [134,184).
        const floats1 = [0, 0, -1, 0, 0, 0, 0, 1, 0, 1, 0, 0];
        for (const value of floats1) {
            expect(view.getFloat32(offset, true)).toBe(value);
            offset += 4;
        }
        expect(view.getUint16(offset, true)).toBe(0);
        offset += 2;
        expect(offset).toBe(buffer.byteLength);

        // Spot-check a known float bit pattern: 1.0f is 00 00 80 3F in
        // little-endian order (normal[2] of triangle 0 at byte 92).
        expect(Array.from(bytes.subarray(92, 96))).toEqual([0x00, 0x00, 0x80, 0x3F]);
        // And -1.0f is 00 00 80 BF (normal[2] of triangle 1 at byte 142).
        expect(Array.from(bytes.subarray(142, 146))).toEqual([0x00, 0x00, 0x80, 0xBF]);
    });

    it('round-trips a small triangle soup exactly', () => {
        const file = new STLBinaryFile();
        file.header.set([0x53, 0x54, 0x4C, 0x21]);  // "STL!"
        // Values chosen to be exactly representable in float32 so the
        // round trip is bit-exact.
        file.triangles.push(makeTriangle(
            [0, 0, 1], [0.5, -0.25, 1.5], [2, 0.125, -3], [4.75, 5, 6], 7));
        file.triangles.push(makeTriangle(
            [1, 0, 0], [10, 20, 30], [-40, 50, -60], [0.0625, -0.03125, 2048]));
        file.triangles.push(makeTriangle(
            [0, -1, 0], [1, 1, 1], [2, 2, 2], [3, 3, 3], 65535));

        const buffer = file.save();

        const loaded = new STLBinaryFile();
        expect(loaded.load(buffer)).toBe(true);
        expect(Array.from(loaded.header)).toEqual(Array.from(file.header));
        expect(loaded.triangles).toEqual(file.triangles);

        // Saving the loaded file reproduces the buffer byte-for-byte.
        const resaved = loaded.save();
        expect(Array.from(new Uint8Array(resaved)))
            .toEqual(Array.from(new Uint8Array(buffer)));
    });

    it('converts double-precision inputs through float32 on save', () => {
        const file = new STLBinaryFile();
        file.triangles.push(makeTriangle(
            [0, 0, 1], [0.1, 0.2, 0.3], [0, 0, 0], [0, 0, 0]));

        const loaded = new STLBinaryFile();
        expect(loaded.load(file.save())).toBe(true);
        const v = loaded.triangles[0].vertex[0];
        expect(v[0]).toBe(Math.fround(0.1));
        expect(v[1]).toBe(Math.fround(0.2));
        expect(v[2]).toBe(Math.fround(0.3));
        expect(v[0]).not.toBe(0.1);
    });

    it('loads an empty file (zero triangles)', () => {
        const file = new STLBinaryFile();
        const buffer = file.save();
        expect(buffer.byteLength).toBe(84);

        const loaded = new STLBinaryFile();
        expect(loaded.load(buffer)).toBe(true);
        expect(loaded.triangles).toEqual([]);
    });

    it('rejects truncated buffers like the upstream read failures', () => {
        const file = new STLBinaryFile();
        file.triangles.push(makeTriangle([0, 0, 1], [0, 0, 0], [1, 0, 0], [0, 1, 0]));
        const buffer = file.save();

        const loaded = new STLBinaryFile();
        // Truncated header.
        expect(loaded.load(buffer.slice(0, 79))).toBe(false);
        // Header present, triangle count truncated.
        expect(loaded.load(buffer.slice(0, 82))).toBe(false);
        // Count announces one triangle but its data is truncated.
        expect(loaded.load(buffer.slice(0, 133))).toBe(false);
        // The complete buffer succeeds.
        expect(loaded.load(buffer)).toBe(true);
    });
});

describe('STLBinaryFile verification', () => {
    // float32 values so that save/load is exact rather than rounding.
    const f32 = fc.float({ noNaN: true, min: -1e6, max: 1e6 });
    const tuple3 = fc.tuple(f32, f32, f32);
    const stlTriangle = fc.record({
        normal: tuple3,
        vertex: fc.tuple(tuple3, tuple3, tuple3),
        attributeByteCount: fc.integer({ min: 0, max: 0xffff })
    });

    const build = (header: number[], tris: {
        normal: [number, number, number],
        vertex: [[number, number, number], [number, number, number],
            [number, number, number]],
        attributeByteCount: number
    }[]): STLBinaryFile => {
        const file = new STLBinaryFile();
        file.header = new Uint8Array(80);
        file.header.set(header.slice(0, 80));
        file.triangles = tris.map((t) => {
            const triangle = new STLTriangle();
            triangle.normal = t.normal.slice() as [number, number, number];
            triangle.vertex = t.vertex.map((v) => v.slice()) as
                [[number, number, number], [number, number, number],
                    [number, number, number]];
            triangle.attributeByteCount = t.attributeByteCount;
            return triangle;
        });
        return file;
    };

    it('save then load reproduces header and triangles exactly', () => {
        check(fc.tuple(
            fc.array(fc.integer({ min: 0, max: 255 }),
                { minLength: 80, maxLength: 80 }),
            fc.array(stlTriangle, { maxLength: 6 })),
        ([header, tris]) => {
            const file = build(header, tris as never);
            const buffer = file.save();
            const reloaded = new STLBinaryFile();
            expect(reloaded.load(buffer)).toBe(true);
            expect(Array.from(reloaded.header)).toEqual(header);
            expect(reloaded.triangles.length).toBe(tris.length);
            reloaded.triangles.forEach((t, i) => {
                expect(t.normal).toEqual(file.triangles[i].normal);
                expect(t.vertex).toEqual(file.triangles[i].vertex);
                expect(t.attributeByteCount)
                    .toBe(file.triangles[i].attributeByteCount);
            });
        });
    });

    it('the saved buffer is 84 + 50 * numTriangles bytes', () => {
        check(fc.array(stlTriangle, { maxLength: 8 }), (tris) => {
            const file = build([], tris as never);
            const buffer = file.save();
            expect(buffer.byteLength).toBe(84 + 50 * tris.length);
            // The triangle count field is little-endian at byte 80.
            expect(new DataView(buffer).getUint32(80, true)).toBe(tris.length);
        });
    });

    it('decodes a hand-built little-endian buffer byte for byte', () => {
        check(fc.array(stlTriangle, { minLength: 1, maxLength: 4 }),
            (tris) => {
                // Build the buffer by hand from the format specification.
                const bytes = new Uint8Array(84 + 50 * tris.length);
                const view = new DataView(bytes.buffer);
                for (let i = 0; i < 80; ++i) {
                    bytes[i] = (7 * i + 3) & 0xff;
                }
                // Little-endian uint32 count written one byte at a time.
                const n = tris.length;
                bytes[80] = n & 0xff;
                bytes[81] = (n >>> 8) & 0xff;
                bytes[82] = (n >>> 16) & 0xff;
                bytes[83] = (n >>> 24) & 0xff;
                tris.forEach((t, i) => {
                    let o = 84 + 50 * i;
                    for (const c of t.normal) {
                        view.setFloat32(o, c, true);
                        o += 4;
                    }
                    for (const v of t.vertex) {
                        for (const c of v) {
                            view.setFloat32(o, c, true);
                            o += 4;
                        }
                    }
                    bytes[o] = t.attributeByteCount & 0xff;
                    bytes[o + 1] = (t.attributeByteCount >>> 8) & 0xff;
                });

                const file = new STLBinaryFile();
                expect(file.load(bytes.buffer)).toBe(true);
                for (let i = 0; i < 80; ++i) {
                    expect(file.header[i]).toBe((7 * i + 3) & 0xff);
                }
                expect(file.triangles.length).toBe(tris.length);
                file.triangles.forEach((triangle, i) => {
                    expect(triangle.normal).toEqual(tris[i].normal);
                    expect(triangle.vertex).toEqual(tris[i].vertex);
                    expect(triangle.attributeByteCount)
                        .toBe(tris[i].attributeByteCount);
                });
                // Re-saving reproduces the hand-built bytes exactly.
                expect(Array.from(new Uint8Array(file.save())))
                    .toEqual(Array.from(bytes));
            });
    });

    it('truncation at any byte is rejected', () => {
        check(fc.tuple(fc.array(stlTriangle, { minLength: 1, maxLength: 3 }),
            fc.integer({ min: 0, max: 99 })), ([tris, pct]) => {
            const buffer = build([], tris as never).save();
            const cut = Math.floor(buffer.byteLength * pct / 100);
            const file = new STLBinaryFile();
            expect(file.load(buffer.slice(0, cut))).toBe(false);
        });
    });

    it('load copies the header instead of aliasing the input buffer', () => {
        const source = build([], []).save();
        new Uint8Array(source)[5] = 42;
        const file = new STLBinaryFile();
        expect(file.load(source)).toBe(true);
        expect(file.header[5]).toBe(42);
        new Uint8Array(source)[5] = 7;
        expect(file.header[5]).toBe(42);
    });

    it('non-float32 coordinates round to float32 on save', () => {
        check(fc.tuple(finite(-1e3, 1e3), finite(-1e3, 1e3)),
            ([x, y]) => {
                const file = new STLBinaryFile();
                const triangle = new STLTriangle();
                triangle.normal = [x, y, x + y];
                triangle.vertex = [[x, y, 0], [0, x, y], [y, 0, x]];
                file.triangles = [triangle];
                const reloaded = new STLBinaryFile();
                expect(reloaded.load(file.save())).toBe(true);
                expect(reloaded.triangles[0].normal)
                    .toEqual([Math.fround(x), Math.fround(y),
                        Math.fround(x + y)]);
            });
    });
});
