import { describe, it, expect } from 'vitest';
import { IndexAttribute } from '../src/IndexAttribute.js';
import { check, fc } from './helpers/arbitraries.js';

describe('IndexAttribute', () => {
    it('defaults to a null source with no supported view', () => {
        const attribute = new IndexAttribute();
        expect(attribute.source).toBeNull();
        expect(attribute.size).toBe(0);
        attribute.setTriangle(0, 1, 2, 3);  // no-op
        expect(attribute.getTriangle(0)).toEqual({ v0: 0, v1: 0, v2: 0 });
    });

    it('reads and writes uint32 triangles (size 4)', () => {
        const indices = new Uint32Array(9);
        const attribute = new IndexAttribute(indices, 4);

        attribute.setTriangle(0, 0, 1, 2);
        attribute.setTriangle(1, 2, 1, 3);
        attribute.setTriangle(2, 100000, 200000, 4294967295);

        expect(Array.from(indices)).toEqual([0, 1, 2, 2, 1, 3, 100000, 200000, 4294967295]);
        expect(attribute.getTriangle(0)).toEqual({ v0: 0, v1: 1, v2: 2 });
        expect(attribute.getTriangle(1)).toEqual({ v0: 2, v1: 1, v2: 3 });
        expect(attribute.getTriangle(2)).toEqual({ v0: 100000, v1: 200000, v2: 4294967295 });
    });

    it('reads and writes uint16 triangles (size 2) with 16-bit truncation', () => {
        const indices = new Uint16Array(6);
        const attribute = new IndexAttribute(indices, 2);

        attribute.setTriangle(0, 5, 6, 7);
        // Values are truncated as by the upstream static_cast<uint16_t>.
        attribute.setTriangle(1, 65535, 65536, 65536 + 9);

        expect(Array.from(indices)).toEqual([5, 6, 7, 65535, 0, 9]);
        expect(attribute.getTriangle(0)).toEqual({ v0: 5, v1: 6, v2: 7 });
        expect(attribute.getTriangle(1)).toEqual({ v0: 65535, v1: 0, v2: 9 });
    });

    it('accepts a raw ArrayBuffer as the source', () => {
        const buffer = new ArrayBuffer(6 * 4);
        const attribute = new IndexAttribute(buffer, 4);

        attribute.setTriangle(0, 10, 11, 12);
        attribute.setTriangle(1, 13, 14, 15);

        expect(Array.from(new Uint32Array(buffer))).toEqual([10, 11, 12, 13, 14, 15]);
        expect(attribute.getTriangle(1)).toEqual({ v0: 13, v1: 14, v2: 15 });
    });

    it('respects the byte offset of a view into a shared buffer', () => {
        // 16 bytes of unrelated data followed by two uint32 triangles.
        const buffer = new ArrayBuffer(16 + 6 * 4);
        const view = new Uint32Array(buffer, 16, 6);
        const attribute = new IndexAttribute(view, 4);

        attribute.setTriangle(0, 1, 2, 3);
        attribute.setTriangle(1, 4, 5, 6);

        // The prefix is untouched.
        expect(Array.from(new Uint32Array(buffer, 0, 4))).toEqual([0, 0, 0, 0]);
        expect(Array.from(view)).toEqual([1, 2, 3, 4, 5, 6]);
        expect(attribute.getTriangle(0)).toEqual({ v0: 1, v1: 2, v2: 3 });
    });

    it('treats unsupported index sizes as no-ops that read zeros', () => {
        const indices = new Uint8Array([9, 9, 9]);
        const attribute = new IndexAttribute(indices, 1);

        attribute.setTriangle(0, 1, 2, 3);  // no-op
        expect(Array.from(indices)).toEqual([9, 9, 9]);
        expect(attribute.getTriangle(0)).toEqual({ v0: 0, v1: 0, v2: 0 });
    });
});

describe('IndexAttribute verification', () => {
    const triangles = (numTriangles: number) =>
        fc.array(fc.tuple(
            fc.integer({ min: 0, max: 0xffffffff }),
            fc.integer({ min: 0, max: 0xffffffff }),
            fc.integer({ min: 0, max: 0xffffffff })),
        { minLength: numTriangles, maxLength: numTriangles });

    it('uint32 indices round trip through set/get', () => {
        check(triangles(6), (tris) => {
            const buffer = new ArrayBuffer(4 * 3 * tris.length);
            const attr = new IndexAttribute(buffer, 4);
            tris.forEach(([v0, v1, v2], t) => attr.setTriangle(t, v0, v1, v2));
            tris.forEach(([v0, v1, v2], t) => {
                expect(attr.getTriangle(t)).toEqual({ v0, v1, v2 });
            });
            // The raw layout is 3 consecutive uint32 per triangle.
            const raw = new Uint32Array(buffer);
            tris.forEach(([v0, v1, v2], t) => {
                expect(raw[3 * t]).toBe(v0 >>> 0);
                expect(raw[3 * t + 1]).toBe(v1 >>> 0);
                expect(raw[3 * t + 2]).toBe(v2 >>> 0);
            });
        });
    });

    it('uint16 indices truncate exactly as static_cast<uint16_t>', () => {
        check(triangles(5), (tris) => {
            const buffer = new ArrayBuffer(2 * 3 * tris.length);
            const attr = new IndexAttribute(buffer, 2);
            tris.forEach(([v0, v1, v2], t) => attr.setTriangle(t, v0, v1, v2));
            tris.forEach(([v0, v1, v2], t) => {
                expect(attr.getTriangle(t)).toEqual({
                    v0: v0 & 0xffff,
                    v1: v1 & 0xffff,
                    v2: v2 & 0xffff
                });
            });
        });
    });

    it('writes to one triangle never disturb its neighbours', () => {
        check(fc.tuple(triangles(4), fc.integer({ min: 0, max: 3 }),
            fc.integer({ min: 0, max: 1000 })), ([tris, t, v]) => {
            const attr = new IndexAttribute(new ArrayBuffer(4 * 3 * 4), 4);
            tris.forEach(([v0, v1, v2], i) => attr.setTriangle(i, v0, v1, v2));
            attr.setTriangle(t, v, v, v);
            for (let i = 0; i < 4; ++i) {
                if (i === t) {
                    expect(attr.getTriangle(i)).toEqual({ v0: v, v1: v, v2: v });
                } else {
                    const [v0, v1, v2] = tris[i];
                    expect(attr.getTriangle(i))
                        .toEqual({ v0: v0 >>> 0, v1: v1 >>> 0, v2: v2 >>> 0 });
                }
            }
        });
    });

    it('unsupported index sizes are inert, as upstream documents', () => {
        check(fc.integer({ min: -4, max: 12 })
            .filter((size) => size !== 2 && size !== 4), (size) => {
            const buffer = new ArrayBuffer(64);
            const attr = new IndexAttribute(buffer, size);
            attr.setTriangle(0, 1, 2, 3);
            expect(attr.getTriangle(0)).toEqual({ v0: 0, v1: 0, v2: 0 });
            // Nothing was written to the buffer.
            expect(Array.from(new Uint8Array(buffer)).every((b) => b === 0))
                .toBe(true);
        });
    });

    it('a view source is addressed from its own byteOffset', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 5 }),
            fc.integer({ min: 0, max: 0xffff })), ([offset, v]) => {
            const buffer = new ArrayBuffer(4 * (3 * 4 + 8));
            const view = new Uint32Array(buffer, 4 * offset, 3 * 4);
            const attr = new IndexAttribute(view, 4);
            attr.setTriangle(1, v, v + 1, v + 2);
            const raw = new Uint32Array(buffer);
            expect(raw[offset + 3]).toBe(v);
            expect(raw[offset + 4]).toBe(v + 1);
            expect(raw[offset + 5]).toBe(v + 2);
            expect(attr.getTriangle(1)).toEqual({ v0: v, v1: v + 1, v2: v + 2 });
        });
    });

    it('a null source leaves the attribute inert instead of crashing', () => {
        const attr = new IndexAttribute();
        expect(attr.source).toBeNull();
        expect(attr.size).toBe(0);
        attr.setTriangle(0, 1, 2, 3);
        expect(attr.getTriangle(0)).toEqual({ v0: 0, v1: 0, v2: 0 });
        const sized = new IndexAttribute(null, 4);
        sized.setTriangle(0, 1, 2, 3);
        expect(sized.getTriangle(0)).toEqual({ v0: 0, v1: 0, v2: 0 });
    });
});
