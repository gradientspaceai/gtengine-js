// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) VertexAttribute.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: upstream stores a raw 'void* source'. JavaScript has no raw
// pointers, so the port stores the source as an ArrayBuffer or any
// ArrayBufferView (typed array or DataView). The consumer uses 'source' and
// 'stride' (bytes between consecutive vertices) to access the attribute
// data, typically through a DataView or typed array over the underlying
// buffer.

export class VertexAttribute {
    // The 'semantic' string allows you to query for a specific vertex
    // attribute and use the 'source' and 'stride' to access the data of the
    // attribute. For example, you might use the semantics "position"
    // (px,py,pz), "normal" (nx,ny,nz), "tcoord" (texture coordinates (u,v)),
    // "dpdu" (derivative of position with respect to u), or "dpdv"
    // (derivative of position with respect to v) for mesh vertices.
    //
    // The source must be 4-byte aligned; that is, when the source is a view,
    // its byteOffset must be a multiple of 4. The stride must be positive
    // and a multiple of 4. The stride constraint is reasonable given that
    // (usually) geometric attributes are arrays of 32-bit or 64-bit
    // floating-point numbers.

    semantic: string;
    source: ArrayBuffer | ArrayBufferView | null;
    stride: number;

    constructor(semantic: string = '', source: ArrayBuffer | ArrayBufferView | null = null,
        stride: number = 0) {
        this.semantic = semantic;
        this.source = source;
        this.stride = stride;
    }
}
