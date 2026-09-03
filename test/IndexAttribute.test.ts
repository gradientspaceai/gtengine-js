import { describe, it, expect } from 'vitest';
import { IndexAttribute } from '../src/IndexAttribute.js';

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
