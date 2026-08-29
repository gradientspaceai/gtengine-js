// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) STLBinaryFile.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// A reader/writer for binary STL files. The file format is described at
// https://en.wikipedia.org/wiki/STL_(file_format).
//
// Port notes: upstream is templated on a Tuple3 type that must represent 3
// contiguous IEEE 32-bit floating-point numbers; the port fixes Tuple3 to
// the triple [number, number, number] (values are converted through 32-bit
// floats when written). Upstream reads/writes std::fstream by filename; the
// port has no file system access, so load takes an ArrayBuffer containing
// the file contents and save returns an ArrayBuffer with the file contents,
// leaving the actual I/O to the caller. The binary layout is identical to
// the on-disk format: an 80-byte header, a uint32 triangle count and 50
// bytes per triangle (twelve float32 values followed by a uint16 attribute
// byte count), all little-endian.

export type STLTuple3 = [number, number, number];

// The port of the upstream nested struct STLBinaryFile::Triangle. The name
// is prefixed to keep the upstream Triangle name free for the port of
// Triangle.h.
export class STLTriangle {
    normal: STLTuple3;
    vertex: [STLTuple3, STLTuple3, STLTuple3];
    attributeByteCount: number;

    constructor() {
        this.normal = [0, 0, 0];
        this.vertex = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
        this.attributeByteCount = 0;
    }
}

export class STLBinaryFile {
    static readonly HEADER_SIZE = 80;
    static readonly TRIANGLE_SIZE = 50;

    header: Uint8Array;
    triangles: STLTriangle[];

    constructor() {
        this.header = new Uint8Array(STLBinaryFile.HEADER_SIZE);
        this.triangles = [];
    }

    // Parse the contents of a binary STL file. The return value is true
    // when the buffer contains a complete file (header, triangle count and
    // all announced triangles); it is false when the buffer is truncated,
    // mirroring the upstream stream-read failure paths.
    load(buffer: ArrayBuffer): boolean {
        const view = new DataView(buffer);
        let offset = 0;

        if (buffer.byteLength < offset + STLBinaryFile.HEADER_SIZE) {
            return false;
        }
        this.header = new Uint8Array(buffer.slice(offset, offset + STLBinaryFile.HEADER_SIZE));
        offset += STLBinaryFile.HEADER_SIZE;

        if (buffer.byteLength < offset + 4) {
            return false;
        }
        const numTriangles = view.getUint32(offset, true);
        offset += 4;

        this.triangles = [];
        for (let t = 0; t < numTriangles; ++t) {
            if (buffer.byteLength < offset + STLBinaryFile.TRIANGLE_SIZE) {
                return false;
            }

            const triangle = new STLTriangle();
            for (let j = 0; j < 3; ++j) {
                triangle.normal[j] = view.getFloat32(offset, true);
                offset += 4;
            }
            for (let v = 0; v < 3; ++v) {
                for (let j = 0; j < 3; ++j) {
                    triangle.vertex[v][j] = view.getFloat32(offset, true);
                    offset += 4;
                }
            }
            triangle.attributeByteCount = view.getUint16(offset, true);
            offset += 2;
            this.triangles.push(triangle);
        }

        return true;
    }

    // The caller is responsible for populating the 'header' and 'triangles'
    // members of the STLBinaryFile object. The returned buffer contains the
    // binary STL file contents. Bytes of 'header' beyond the 80-byte header
    // size are ignored; if 'header' is shorter, the remainder is
    // zero-filled.
    save(): ArrayBuffer {
        const numTriangles = this.triangles.length;
        const buffer = new ArrayBuffer(
            STLBinaryFile.HEADER_SIZE + 4 + STLBinaryFile.TRIANGLE_SIZE * numTriangles);
        const view = new DataView(buffer);
        const bytes = new Uint8Array(buffer);
        let offset = 0;

        const headerLength = Math.min(this.header.length, STLBinaryFile.HEADER_SIZE);
        bytes.set(this.header.subarray(0, headerLength), offset);
        offset += STLBinaryFile.HEADER_SIZE;

        view.setUint32(offset, numTriangles, true);
        offset += 4;

        for (const triangle of this.triangles) {
            for (let j = 0; j < 3; ++j) {
                view.setFloat32(offset, triangle.normal[j], true);
                offset += 4;
            }
            for (let v = 0; v < 3; ++v) {
                for (let j = 0; j < 3; ++j) {
                    view.setFloat32(offset, triangle.vertex[v][j], true);
                    offset += 4;
                }
            }
            view.setUint16(offset, triangle.attributeByteCount, true);
            offset += 2;
        }

        return buffer;
    }
}
