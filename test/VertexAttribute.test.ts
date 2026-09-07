import { describe, it, expect } from 'vitest';
import { VertexAttribute } from '../src/VertexAttribute.js';
import { check, fc } from './helpers/arbitraries.js';

describe('VertexAttribute', () => {
    it('defaults to an empty semantic with a null source', () => {
        const attribute = new VertexAttribute();
        expect(attribute.semantic).toBe('');
        expect(attribute.source).toBeNull();
        expect(attribute.stride).toBe(0);
    });

    it('stores the semantic, source and stride', () => {
        const positions = new Float32Array(12);
        const attribute = new VertexAttribute('position', positions, 12);
        expect(attribute.semantic).toBe('position');
        expect(attribute.source).toBe(positions);
        expect(attribute.stride).toBe(12);
    });

    it('describes attributes interleaved in a shared vertex buffer', () => {
        // Vertex layout: position (3 floats) then normal (3 floats),
        // stride 24 bytes, two vertices.
        const stride = 24;
        const buffer = new ArrayBuffer(2 * stride);
        const f32 = new Float32Array(buffer);
        f32.set([1, 2, 3, 0, 0, 1], 0);   // vertex 0
        f32.set([4, 5, 6, 0, 1, 0], 6);   // vertex 1

        const position = new VertexAttribute('position', new DataView(buffer, 0), stride);
        const normal = new VertexAttribute('normal', new DataView(buffer, 12), stride);

        // Consume the attributes the way a mesh factory would: source view
        // byteOffset locates the first element, stride advances vertices.
        const readTriple = (attribute: VertexAttribute, i: number): number[] => {
            const view = attribute.source as DataView;
            const base = view.byteOffset + i * attribute.stride;
            return Array.from(new Float32Array(view.buffer, base, 3));
        };

        expect(position.semantic).toBe('position');
        expect(readTriple(position, 0)).toEqual([1, 2, 3]);
        expect(readTriple(position, 1)).toEqual([4, 5, 6]);
        expect(readTriple(normal, 0)).toEqual([0, 0, 1]);
        expect(readTriple(normal, 1)).toEqual([0, 1, 0]);
    });
});

describe('VertexAttribute verification', () => {
    it('the constructor stores its three arguments unchanged', () => {
        check(fc.tuple(fc.string(), fc.integer({ min: 0, max: 64 })),
            ([semantic, stride]) => {
                const buffer = new ArrayBuffer(64);
                const attr = new VertexAttribute(semantic, buffer, stride);
                expect(attr.semantic).toBe(semantic);
                expect(attr.source).toBe(buffer);
                expect(attr.stride).toBe(stride);
            });
    });

    it('stride-based addressing reads back interleaved vertices', () => {
        // Vertex layout: position (3 floats) then normal (3 floats).
        check(fc.array(fc.tuple(
            fc.float({ noNaN: true, min: -100, max: 100 }),
            fc.float({ noNaN: true, min: -100, max: 100 }),
            fc.float({ noNaN: true, min: -100, max: 100 })),
        { minLength: 1, maxLength: 8 }), (positions) => {
            const stride = 6 * 4;
            const buffer = new ArrayBuffer(stride * positions.length);
            const view = new DataView(buffer);
            const position = new VertexAttribute('position', buffer, stride);
            const normal = new VertexAttribute('normal',
                new Uint8Array(buffer, 3 * 4), stride);
            positions.forEach(([x, y, z], i) => {
                const base = i * position.stride;
                view.setFloat32(base, x, true);
                view.setFloat32(base + 4, y, true);
                view.setFloat32(base + 8, z, true);
            });
            positions.forEach(([x, y, z], i) => {
                const base = i * position.stride;
                expect(view.getFloat32(base, true)).toBe(x);
                expect(view.getFloat32(base + 4, true)).toBe(y);
                expect(view.getFloat32(base + 8, true)).toBe(z);
            });
            // The normal attribute names the same buffer at a byte offset
            // that is 4-byte aligned, as the upstream comment requires.
            const source = normal.source as Uint8Array;
            expect(source.byteOffset % 4).toBe(0);
            expect(normal.stride % 4).toBe(0);
        });
    });

    it('the default construction matches upstream defaults', () => {
        const attr = new VertexAttribute();
        expect(attr.semantic).toBe('');
        expect(attr.source).toBeNull();
        expect(attr.stride).toBe(0);
    });
});
