import { describe, it, expect } from 'vitest';
import { VertexAttribute } from '../src/VertexAttribute.js';

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
