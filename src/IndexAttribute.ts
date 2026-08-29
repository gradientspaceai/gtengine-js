// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) IndexAttribute.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The IndexAttribute class represents an array of triples of indices into a
// vertex array for an indexed triangle mesh. For now, the source must be
// either uint16 or uint32 indices.
//
// Port notes: upstream stores a raw 'void* source' and reinterpret_casts it
// to uint16_t* or uint32_t* based on 'size' (the number of bytes per index).
// JavaScript has no raw pointers, so the port stores the source as an
// ArrayBuffer or any ArrayBufferView (typed array or DataView) and creates
// the appropriately typed view over the underlying buffer at construction:
// size 4 selects a Uint32Array view, size 2 selects a Uint16Array view, and
// any other size leaves the attribute without a supported view, in which
// case setTriangle is a no-op and getTriangle returns zeros, matching the
// upstream "Unsupported type" behavior. The upstream 4-byte alignment
// requirement translates to the source byteOffset being a multiple of the
// index size when a view is supplied. The C++ out-parameters of GetTriangle
// become a returned object literal { v0, v1, v2 }.

export class IndexAttribute {
    // The source of the index data. The size is the number of bytes per
    // index (2 for uint16 indices, 4 for uint32 indices).
    readonly source: ArrayBuffer | ArrayBufferView | null;
    readonly size: number;

    // The typed view over the source buffer, selected by 'size'. It is null
    // when the source is null or the size is unsupported.
    private readonly mIndices: Uint16Array | Uint32Array | null;

    // Construction.
    constructor(source: ArrayBuffer | ArrayBufferView | null = null, size: number = 0) {
        this.source = source;
        this.size = size;
        this.mIndices = null;
        if (source !== null) {
            let buffer: ArrayBuffer;
            let byteOffset: number;
            let byteLength: number;
            if (source instanceof ArrayBuffer) {
                buffer = source;
                byteOffset = 0;
                byteLength = source.byteLength;
            } else {
                buffer = source.buffer as ArrayBuffer;
                byteOffset = source.byteOffset;
                byteLength = source.byteLength;
            }

            if (size === 4) {
                this.mIndices = new Uint32Array(buffer, byteOffset, Math.floor(byteLength / 4));
            } else if (size === 2) {
                this.mIndices = new Uint16Array(buffer, byteOffset, Math.floor(byteLength / 2));
            }
            // Otherwise the type is unsupported and mIndices remains null.
        }
    }

    // Triangle access.
    setTriangle(t: number, v0: number, v1: number, v2: number): void {
        const indices = this.mIndices;
        if (indices !== null) {
            const base = 3 * t;
            // Uint16Array assignment truncates to 16 bits, matching the
            // upstream static_cast<uint16_t>.
            indices[base] = v0;
            indices[base + 1] = v1;
            indices[base + 2] = v2;
            return;
        }

        // Unsupported type.
    }

    getTriangle(t: number): { v0: number, v1: number, v2: number } {
        const indices = this.mIndices;
        if (indices !== null) {
            const base = 3 * t;
            return {
                v0: indices[base],
                v1: indices[base + 1],
                v2: indices[base + 2]
            };
        }

        // Unsupported type.
        return { v0: 0, v1: 0, v2: 0 };
    }
}
