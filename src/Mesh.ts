// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Mesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Mesh class is designed to support triangulations of surfaces of a
// small number of topologies. See the documents
//   https://www.geometrictools.com/MeshDifferentialGeometry.pdf
//   https://www.geometrictools.com/MeshFactory.pdf
// for details.
//
// You must set the vertex attribute sources before calling update().
//
// The semantic "position" is required and its source must be an array of
// floating-point numbers with at least 3 channels so that positions are
// computed as 3-dimensional vectors.
//
// The positions are assumed to be parameterized by texture coordinates
// (u,v); the position is thought of as a function P(u,v). If texture
// coordinates are provided, the semantic must be "tcoord". If texture
// coordinates are not provided, default texture coordinates are computed
// internally as described in the mesh factory document.
//
// The frame for the tangent space is optional. All vectors in the frame
// must have sources that are arrays of floating-point numbers with at least
// 3 channels per attribute. If normal vectors are provided, the semantic
// must be "normal".
//
// Two options are supported for tangent vectors. The first option is that
// the tangents are surface derivatives dP/du and dP/dv, which are not
// necessarily unit length or orthogonal. The semantics must be "dpdu" and
// "dpdv". The second option is that the tangents are unit length and
// orthogonal, with the infrequent possibility that a vertex is degenerate
// in that dP/du and dP/dv are linearly dependent. The semantics must be
// "tangent" and "bitangent".
//
// For each provided vertex attribute, a derived class can initialize that
// attribute by overriding one of the initialize*() functions.
//
// Port notes:
//
// * Upstream stores each bound channel as a raw 'Vector3<Real>*' plus a
//   byte stride and dereferences it with pointer arithmetic. JavaScript has
//   no raw pointers, so the port wraps the (source, stride) pair of a
//   VertexAttribute in a MeshChannel object that owns a typed-array view of
//   the source buffer and performs the same (byteOffset + i*stride)
//   addressing. The channel's element type is float64 unless the supplied
//   source is a Float32Array, in which case it is float32; this is the port
//   of instantiating Mesh<double> versus Mesh<float>. The 'm*Stride' members
//   of upstream are folded into the channel objects (mPositions.stride and
//   so on), so a derived class that substitutes its own storage (as
//   RectangleMesh does for default texture coordinates) assigns a single
//   channel, e.g. 'this.mTCoords = MeshChannel.allocate(numVertices, 2)'.
//
// * Upstream's accessors Position(i), Normal(i), ... return references into
//   the channel memory and are both read and written through. TypeScript
//   cannot return a reference to buffer memory, so each accessor is a pair:
//   position(i) returns a copy as a Vector, and setPosition(i, value)
//   writes. Read-modify-write sequences such as 'Normal(v0) += n' become
//   read, add, write.
//
// * Upstream Mesh copies the MeshDescription by value, so mutations of the
//   internal-use members (constructed, hasTangentSpaceVectors,
//   allowUpdateFrame) are not visible in the caller's object; the caller
//   examines mesh.getDescription() instead. The port preserves this by
//   cloning the description in the constructor. The clone shares the
//   VertexAttribute and IndexAttribute objects, matching the C++ copies that
//   share the underlying buffers.
//
// * The two upstream constructors of MeshDescription are distinguished by
//   the argument types, which erase to 'number' here; the ARBITRARY-topology
//   constructor is therefore the static factory MeshDescription.arbitrary().

import { IndexAttribute } from './IndexAttribute.js';
import { logAssert } from './Logger.js';
import { Matrix, multiplyAB, outerProduct } from './Matrix.js';
import { inverse2x2 } from './Matrix2x2.js';
import { Vector, add, dot, normalize, sub } from './Vector.js';
import { computeOrthogonalComplement3, cross } from './Vector3.js';
import { VertexAttribute } from './VertexAttribute.js';

export enum MeshTopology {
    ARBITRARY,
    RECTANGLE,
    CYLINDER,
    TORUS,
    DISK,
    SPHERE
}

// The port of the (source, stride) channel addressing of upstream's raw
// 'Vector3<Real>*' and 'Vector2<Real>*' members. The stride is measured in
// bytes, as in VertexAttribute.
export class MeshChannel {
    // The source of the attribute data and the number of bytes between
    // consecutive vertices.
    readonly source: ArrayBuffer | ArrayBufferView;
    readonly stride: number;

    // The number of components read and written per vertex (3 for the
    // position and frame channels, 2 for texture coordinates).
    readonly numComponents: number;

    private readonly mData: Float32Array | Float64Array;
    private readonly mElementStride: number;

    constructor(source: ArrayBuffer | ArrayBufferView, stride: number,
        numComponents: number) {
        logAssert(stride > 0, 'The stride must be positive.');
        logAssert(numComponents > 0, 'The number of components must be positive.');

        this.source = source;
        this.stride = stride;
        this.numComponents = numComponents;

        let buffer: ArrayBuffer;
        let byteOffset: number;
        let byteLength: number;
        let isFloat32: boolean;
        if (source instanceof ArrayBuffer) {
            buffer = source;
            byteOffset = 0;
            byteLength = source.byteLength;
            isFloat32 = false;
        } else {
            buffer = source.buffer as ArrayBuffer;
            byteOffset = source.byteOffset;
            byteLength = source.byteLength;
            isFloat32 = (source instanceof Float32Array);
        }

        const bytesPerElement = (isFloat32 ? 4 : 8);
        logAssert(byteOffset % bytesPerElement === 0 &&
            stride % bytesPerElement === 0,
            'The source offset and stride must be multiples of the element size.');

        const numElements = Math.floor(byteLength / bytesPerElement);
        this.mData = (isFloat32
            ? new Float32Array(buffer, byteOffset, numElements)
            : new Float64Array(buffer, byteOffset, numElements));
        this.mElementStride = stride / bytesPerElement;
    }

    // Create a channel backed by tightly packed float64 storage. This is the
    // port of a derived class allocating its own std::vector of vectors and
    // pointing a channel at it (for example, the default texture coordinates
    // of RectangleMesh).
    static allocate(numVertices: number, numComponents: number): MeshChannel {
        const data = new Float64Array(numVertices * numComponents);
        return new MeshChannel(data, 8 * numComponents, numComponents);
    }

    // Create a channel for a vertex attribute whose source is not null.
    static fromAttribute(attribute: VertexAttribute,
        numComponents: number): MeshChannel {
        logAssert(attribute.source !== null, 'The attribute source is null.');
        return new MeshChannel(attribute.source as ArrayBuffer | ArrayBufferView,
            attribute.stride, numComponents);
    }

    // The attribute of vertex i, returned as a copy.
    get(i: number): Vector {
        const result = new Vector(this.numComponents);
        const base = i * this.mElementStride;
        for (let k = 0; k < this.numComponents; ++k) {
            result.values[k] = this.mData[base + k];
        }
        return result;
    }

    // The value of component k of the attribute of vertex i.
    getComponent(i: number, k: number): number {
        return this.mData[i * this.mElementStride + k];
    }

    // Assign the attribute of vertex i. The input must have at least
    // numComponents components; any additional components are ignored, as
    // are any additional channels of the source (upstream writes only the
    // leading components of each vertex record).
    set(i: number, value: Vector): void {
        const base = i * this.mElementStride;
        for (let k = 0; k < this.numComponents; ++k) {
            this.mData[base + k] = value.values[k];
        }
    }

    setComponent(i: number, k: number, value: number): void {
        this.mData[i * this.mElementStride + k] = value;
    }
}

export class MeshDescription {
    topology: MeshTopology;
    numVertices: number;
    numTriangles: number;
    vertexAttributes: VertexAttribute[];
    indexAttribute: IndexAttribute;
    wantDynamicTangentSpaceUpdate: boolean;  // default: false
    wantCCW: boolean;  // default: true

    // For internal use only.
    hasTangentSpaceVectors: boolean;
    allowUpdateFrame: boolean;
    numRows: number;
    numCols: number;
    rMax: number;
    cMax: number;
    rIncrement: number;

    // After an attempt to construct a Mesh or Mesh-derived object, examine
    // this value to determine whether the construction was successful.
    constructed: boolean;

    // Constructor for topologies other than MeshTopology.ARBITRARY. Compute
    // the number of vertices and triangles for the mesh based on the
    // requested number of rows and columns. If the number of rows or columns
    // is invalid for the specified topology, they are modified to be valid,
    // in which case inNumRows/numRows and inNumCols/numCols can differ. If
    // the input topology is MeshTopology.ARBITRARY, then inNumRows and
    // inNumCols are assigned to numVertices and numTriangles, respectively,
    // and numRows and numCols are set to zero. The remaining members must be
    // set explicitly by the client.
    constructor(inTopology: MeshTopology, inNumRows: number, inNumCols: number) {
        this.topology = inTopology;
        this.numVertices = 0;
        this.numTriangles = 0;
        this.vertexAttributes = [];
        this.indexAttribute = new IndexAttribute();
        this.wantDynamicTangentSpaceUpdate = false;
        this.wantCCW = true;
        this.hasTangentSpaceVectors = false;
        this.allowUpdateFrame = false;
        this.numRows = 0;
        this.numCols = 0;
        this.rMax = 0;
        this.cMax = 0;
        this.rIncrement = 0;
        this.constructed = false;

        switch (this.topology) {
            case MeshTopology.ARBITRARY:
                this.numVertices = inNumRows;
                this.numTriangles = inNumCols;
                this.numRows = 0;
                this.numCols = 0;
                this.rMax = 0;
                this.cMax = 0;
                this.rIncrement = 0;
                break;

            case MeshTopology.RECTANGLE:
                this.numRows = Math.max(inNumRows, 2);
                this.numCols = Math.max(inNumCols, 2);
                this.rMax = this.numRows - 1;
                this.cMax = this.numCols - 1;
                this.rIncrement = this.numCols;
                this.numVertices = (this.rMax + 1) * (this.cMax + 1);
                this.numTriangles = 2 * this.rMax * this.cMax;
                break;

            case MeshTopology.CYLINDER:
                this.numRows = Math.max(inNumRows, 2);
                this.numCols = Math.max(inNumCols, 3);
                this.rMax = this.numRows - 1;
                this.cMax = this.numCols;
                this.rIncrement = this.numCols + 1;
                this.numVertices = (this.rMax + 1) * (this.cMax + 1);
                this.numTriangles = 2 * this.rMax * this.cMax;
                break;

            case MeshTopology.TORUS:
                this.numRows = Math.max(inNumRows, 2);
                this.numCols = Math.max(inNumCols, 3);
                this.rMax = this.numRows;
                this.cMax = this.numCols;
                this.rIncrement = this.numCols + 1;
                this.numVertices = (this.rMax + 1) * (this.cMax + 1);
                this.numTriangles = 2 * this.rMax * this.cMax;
                break;

            case MeshTopology.DISK:
                this.numRows = Math.max(inNumRows, 1);
                this.numCols = Math.max(inNumCols, 3);
                this.rMax = this.numRows - 1;
                this.cMax = this.numCols;
                this.rIncrement = this.numCols + 1;
                this.numVertices = (this.rMax + 1) * (this.cMax + 1) + 1;
                this.numTriangles = 2 * this.rMax * this.cMax + this.numCols;
                break;

            case MeshTopology.SPHERE:
                this.numRows = Math.max(inNumRows, 1);
                this.numCols = Math.max(inNumCols, 3);
                this.rMax = this.numRows - 1;
                this.cMax = this.numCols;
                this.rIncrement = this.numCols + 1;
                this.numVertices = (this.rMax + 1) * (this.cMax + 1) + 2;
                this.numTriangles = 2 * this.rMax * this.cMax + 2 * this.numCols;
                break;
        }
    }

    // Constructor for MeshTopology.ARBITRARY. The members topology,
    // numVertices, and numTriangles are set in the obvious manner. The
    // members numRows and numCols are set to zero. The remaining members
    // must be set explicitly by the client.
    static arbitrary(inNumVertices: number, inNumTriangles: number): MeshDescription {
        logAssert(inNumVertices >= 3 && inNumTriangles >= 1, 'Invalid input.');
        const description = new MeshDescription(MeshTopology.ARBITRARY,
            inNumVertices, inNumTriangles);
        return description;
    }

    // The port of the C++ copy constructor. The vertex attributes and the
    // index attribute are shared, matching the C++ copies that refer to the
    // same underlying buffers.
    clone(): MeshDescription {
        const result = new MeshDescription(MeshTopology.ARBITRARY, 0, 0);
        result.topology = this.topology;
        result.numVertices = this.numVertices;
        result.numTriangles = this.numTriangles;
        result.vertexAttributes = this.vertexAttributes.slice();
        result.indexAttribute = this.indexAttribute;
        result.wantDynamicTangentSpaceUpdate = this.wantDynamicTangentSpaceUpdate;
        result.wantCCW = this.wantCCW;
        result.hasTangentSpaceVectors = this.hasTangentSpaceVectors;
        result.allowUpdateFrame = this.allowUpdateFrame;
        result.numRows = this.numRows;
        result.numCols = this.numCols;
        result.rMax = this.rMax;
        result.cMax = this.cMax;
        result.rIncrement = this.rIncrement;
        result.constructed = this.constructed;
        return result;
    }
}

export class Mesh {
    // Constructor inputs.
    // The client requests dynamic tangent-space updates via the description;
    // however, if it is requested and the vertex attributes do not contain
    // entries for "tangent", "bitangent", "dpdu", or "dpdv", then
    // allowUpdateFrame is set to false.
    protected mDescription: MeshDescription;

    // Copied from the vertex attributes when available. Each channel owns
    // the source and the stride of the corresponding attribute.
    protected mPositions: MeshChannel | null;
    protected mNormals: MeshChannel | null;
    protected mTangents: MeshChannel | null;
    protected mBitangents: MeshChannel | null;
    protected mDPDUs: MeshChannel | null;
    protected mDPDVs: MeshChannel | null;
    protected mTCoords: MeshChannel | null;

    // When dynamic tangent-space updates are requested, the update algorithm
    // requires texture coordinates (user-specified or non-local). It is
    // possible to create a vertex-adjacent set (with indices into the vertex
    // array) for each mesh vertex; however, instead we rely on a triangle
    // iteration and incrementally store the information needed for the
    // estimation of the tangent space. Each vertex has associated matrices D
    // and U, but we need to store only U^T*U and D^T*U. See the PDF for
    // details.
    protected mUTU: Matrix[];
    protected mDTU: Matrix[];

    // Construction. This constructor is for ARBITRARY topology. The vertices
    // and indices must already be assigned by the client. Derived classes
    // pass their own set of valid topologies, but assignment of vertices and
    // indices occurs in the derived-class constructors.
    constructor(description: MeshDescription, validTopologies: MeshTopology[]) {
        this.mDescription = description.clone();
        this.mPositions = null;
        this.mNormals = null;
        this.mTangents = null;
        this.mBitangents = null;
        this.mDPDUs = null;
        this.mDPDVs = null;
        this.mTCoords = null;
        this.mUTU = [];
        this.mDTU = [];

        this.mDescription.constructed = false;
        for (const topology of validTopologies) {
            if (this.mDescription.topology === topology) {
                this.mDescription.constructed = true;
                break;
            }
        }

        logAssert(this.mDescription.indexAttribute.source !== null,
            'The mesh needs triangles/indices in Mesh constructor.');

        // Set sources for the requested vertex attributes.
        this.mDescription.hasTangentSpaceVectors = false;
        this.mDescription.allowUpdateFrame =
            this.mDescription.wantDynamicTangentSpaceUpdate;
        for (const attribute of this.mDescription.vertexAttributes) {
            if (attribute.source !== null && attribute.stride > 0) {
                if (attribute.semantic === 'position') {
                    this.mPositions = MeshChannel.fromAttribute(attribute, 3);
                    continue;
                }

                if (attribute.semantic === 'normal') {
                    this.mNormals = MeshChannel.fromAttribute(attribute, 3);
                    continue;
                }

                if (attribute.semantic === 'tangent') {
                    this.mTangents = MeshChannel.fromAttribute(attribute, 3);
                    this.mDescription.hasTangentSpaceVectors = true;
                    continue;
                }

                if (attribute.semantic === 'bitangent') {
                    this.mBitangents = MeshChannel.fromAttribute(attribute, 3);
                    this.mDescription.hasTangentSpaceVectors = true;
                    continue;
                }

                if (attribute.semantic === 'dpdu') {
                    this.mDPDUs = MeshChannel.fromAttribute(attribute, 3);
                    this.mDescription.hasTangentSpaceVectors = true;
                    continue;
                }

                if (attribute.semantic === 'dpdv') {
                    this.mDPDVs = MeshChannel.fromAttribute(attribute, 3);
                    this.mDescription.hasTangentSpaceVectors = true;
                    continue;
                }

                if (attribute.semantic === 'tcoord') {
                    this.mTCoords = MeshChannel.fromAttribute(attribute, 2);
                    continue;
                }
            }
        }

        logAssert(this.mPositions !== null,
            'The mesh needs positions in Mesh constructor.');

        // The initial value of allowUpdateFrame is the client request about
        // wanting dynamic tangent-space updates. If the vertex attributes do
        // not include tangent-space vectors, then dynamic updates are not
        // necessary. If tangent-space vectors are present, the update
        // algorithm requires texture coordinates (mTCoords must be nonnull)
        // or must compute local coordinates (mNormals must be nonnull).
        if (this.mDescription.allowUpdateFrame) {
            if (!this.mDescription.hasTangentSpaceVectors) {
                this.mDescription.allowUpdateFrame = false;
            }

            if (!this.mTCoords && !this.mNormals) {
                this.mDescription.allowUpdateFrame = false;
            }
        }

        if (this.mDescription.allowUpdateFrame) {
            for (let i = 0; i < this.mDescription.numVertices; ++i) {
                this.mUTU.push(new Matrix(2, 2));
                this.mDTU.push(new Matrix(3, 2));
            }
        }
    }

    // Member accessors.
    getDescription(): MeshDescription {
        return this.mDescription;
    }

    // If the underlying geometric data varies dynamically, call this
    // function to update whatever vertex attributes are specified by the
    // vertex pool.
    update(): void {
        logAssert(this.mDescription.constructed,
            'The Mesh object failed the construction.');

        this.updatePositions();

        if (this.mDescription.allowUpdateFrame) {
            this.updateFrame();
        }
        else if (this.mNormals) {
            this.updateNormals();
        }
        // else: The mesh has no frame data, so there is nothing to do.
    }

    // Access the vertex attributes. Upstream returns references into the
    // channel memory; the port returns copies from the get-style accessors
    // and provides matching set-style accessors.
    protected position(i: number): Vector {
        return (this.mPositions as MeshChannel).get(i);
    }

    protected setPosition(i: number, value: Vector): void {
        (this.mPositions as MeshChannel).set(i, value);
    }

    protected normal(i: number): Vector {
        return (this.mNormals as MeshChannel).get(i);
    }

    protected setNormal(i: number, value: Vector): void {
        (this.mNormals as MeshChannel).set(i, value);
    }

    protected tangent(i: number): Vector {
        return (this.mTangents as MeshChannel).get(i);
    }

    protected setTangent(i: number, value: Vector): void {
        (this.mTangents as MeshChannel).set(i, value);
    }

    protected bitangent(i: number): Vector {
        return (this.mBitangents as MeshChannel).get(i);
    }

    protected setBitangent(i: number, value: Vector): void {
        (this.mBitangents as MeshChannel).set(i, value);
    }

    protected dpdu(i: number): Vector {
        return (this.mDPDUs as MeshChannel).get(i);
    }

    protected setDPDU(i: number, value: Vector): void {
        (this.mDPDUs as MeshChannel).set(i, value);
    }

    protected dpdv(i: number): Vector {
        return (this.mDPDVs as MeshChannel).get(i);
    }

    protected setDPDV(i: number, value: Vector): void {
        (this.mDPDVs as MeshChannel).set(i, value);
    }

    protected tcoord(i: number): Vector {
        return (this.mTCoords as MeshChannel).get(i);
    }

    protected setTCoord(i: number, value: Vector): void {
        (this.mTCoords as MeshChannel).set(i, value);
    }

    // Compute the indices for non-arbitrary topologies. This function is
    // called by derived classes.
    protected computeIndices(): void {
        let t = 0;
        for (let r = 0, i = 0; r < this.mDescription.rMax; ++r) {
            let v0 = i, v1 = v0 + 1;
            i += this.mDescription.rIncrement;
            let v2 = i, v3 = v2 + 1;
            for (let c = 0; c < this.mDescription.cMax;
                ++c, ++v0, ++v1, ++v2, ++v3) {
                if (this.mDescription.wantCCW) {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v1, v2);
                    this.mDescription.indexAttribute.setTriangle(t++, v1, v3, v2);
                }
                else {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v2, v1);
                    this.mDescription.indexAttribute.setTriangle(t++, v1, v2, v3);
                }
            }
        }

        if (this.mDescription.topology === MeshTopology.DISK) {
            let v0 = 0, v1 = 1;
            const v2 = this.mDescription.numVertices - 1;
            for (let c = 0; c < this.mDescription.numCols; ++c, ++v0, ++v1) {
                if (this.mDescription.wantCCW) {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v2, v1);
                }
                else {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v1, v2);
                }
            }
        }
        else if (this.mDescription.topology === MeshTopology.SPHERE) {
            let v0 = 0, v1 = 1, v2 = this.mDescription.numVertices - 2;
            for (let c = 0; c < this.mDescription.numCols; ++c, ++v0, ++v1) {
                if (this.mDescription.wantCCW) {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v2, v1);
                }
                else {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v1, v2);
                }
            }

            // Upstream Mesh.h computes the first vertex of the last row as
            // (numRows - 1) * numCols. Each row of a SPHERE mesh has
            // numCols + 1 vertices (the last duplicates the first to close
            // the seam) and rIncrement is numCols + 1, so the first vertex
            // of the last row is rMax * rIncrement. The upstream expression
            // is correct only when numRows is 1; otherwise the south-pole
            // fan references vertices of an interior row and the mesh is
            // corrupted, so the port fixes the index.
            v0 = this.mDescription.rMax * this.mDescription.rIncrement;
            v1 = v0 + 1;
            v2 = this.mDescription.numVertices - 1;
            // Upstream winds this second fan as (v0, v2, v1), the same as the
            // first fan. The first pole sits before row 0 (it plays the role
            // of the grid triangle (v1, v3, v2) with r = -1), but the second
            // pole sits after the last row (the role of (v0, v1, v2) with
            // r = rMax), so its consistent winding is (v0, v1, v2). With the
            // upstream order the last-ring edges are traversed twice in the
            // same direction and the second pole's normal points inward. The
            // port fixes the winding (see the upstream-bug issue for B84).
            for (let c = 0; c < this.mDescription.numCols; ++c, ++v0, ++v1) {
                if (this.mDescription.wantCCW) {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v1, v2);
                }
                else {
                    this.mDescription.indexAttribute.setTriangle(t++, v0, v2, v1);
                }
            }
        }
    }

    // The update() function allows derived classes to use algorithms
    // different from least-squares fitting to compute the normals (when no
    // tangent-space information is requested) or to compute the frame
    // (normals and tangent space). The updatePositions() is a stub; the base
    // class has no knowledge about how positions should be modified. A
    // derived class, however, might choose to use dynamic updating and
    // override updatePositions(). The base-class updateNormals() computes
    // vertex normals as averages of area-weighted triangle normals
    // (nonparametric approach). The base-class updateFrame() uses a
    // least-squares algorithm for estimating the tangent space (parametric
    // approach).
    protected updatePositions(): void {
    }

    protected updateNormals(): void {
        // Compute normal vector as normalized weighted averages of triangle
        // normal vectors.

        // Set the normals to zero to allow accumulation of triangle normals.
        const zero = Vector.fromArray([0, 0, 0]);
        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            this.setNormal(i, zero);
        }

        // Accumulate the triangle normals.
        for (let t = 0; t < this.mDescription.numTriangles; ++t) {
            // Get the positions for the triangle.
            const { v0, v1, v2 } = this.mDescription.indexAttribute.getTriangle(t);
            const P0 = this.position(v0);
            const P1 = this.position(v1);
            const P2 = this.position(v2);

            // Get the edge vectors.
            const E1 = sub(P1, P0);
            const E2 = sub(P2, P0);

            // Compute a triangle normal whose length is twice the area of
            // the triangle.
            const triangleNormal = cross(E1, E2);

            // Accumulate the triangle normals.
            this.setNormal(v0, add(this.normal(v0), triangleNormal));
            this.setNormal(v1, add(this.normal(v1), triangleNormal));
            this.setNormal(v2, add(this.normal(v2), triangleNormal));
        }

        // Normalize the normals.
        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            const n = this.normal(i);
            normalize(n, true);
            this.setNormal(i, n);
        }
    }

    protected updateFrame(): void {
        if (!this.mTCoords) {
            // We need to compute vertex normals first in order to compute
            // local texture coordinates. The vertex normals are recomputed
            // later based on estimated tangent vectors.
            this.updateNormals();
        }

        // Use the least-squares algorithm to estimate the tangent-space
        // vectors and, if requested, normal vectors.
        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            this.mUTU[i] = new Matrix(2, 2);
            this.mDTU[i] = new Matrix(3, 2);
        }

        for (let t = 0; t < this.mDescription.numTriangles; ++t) {
            // Get the positions and differences for the triangle.
            const { v0, v1, v2 } = this.mDescription.indexAttribute.getTriangle(t);
            const P0 = this.position(v0);
            const P1 = this.position(v1);
            const P2 = this.position(v2);
            const D10 = sub(P1, P0);
            const D20 = sub(P2, P0);
            const D21 = sub(P2, P1);

            if (this.mTCoords) {
                // Get the texture coordinates and differences for the
                // triangle.
                const C0 = this.tcoord(v0);
                const C1 = this.tcoord(v1);
                const C2 = this.tcoord(v2);
                const U10 = sub(C1, C0);
                const U20 = sub(C2, C0);
                const U21 = sub(C2, C1);

                // Compute the outer products.
                const outerU10 = outerProduct(U10, U10);
                const outerU20 = outerProduct(U20, U20);
                const outerU21 = outerProduct(U21, U21);
                const outerD10 = outerProduct(D10, U10);
                const outerD20 = outerProduct(D20, U20);
                const outerD21 = outerProduct(D21, U21);

                // Keep a running sum of U^T*U and D^T*U.
                accumulate(this.mUTU[v0], outerU10, outerU20);
                accumulate(this.mUTU[v1], outerU10, outerU21);
                accumulate(this.mUTU[v2], outerU20, outerU21);
                accumulate(this.mDTU[v0], outerD10, outerD20);
                accumulate(this.mDTU[v1], outerD10, outerD21);
                accumulate(this.mDTU[v2], outerD20, outerD21);
            }
            else {
                // Compute local coordinates and differences for the
                // triangle.
                let basis: Vector[] = [this.normal(v0), new Vector(3), new Vector(3)];
                computeOrthogonalComplement3(1, basis, true);
                const U10 = Vector.fromArray(
                    [dot(basis[1], D10), dot(basis[2], D10)]);
                const U20 = Vector.fromArray(
                    [dot(basis[1], D20), dot(basis[2], D20)]);
                accumulate(this.mUTU[v0], outerProduct(U10, U10),
                    outerProduct(U20, U20));
                accumulate(this.mDTU[v0], outerProduct(D10, U10),
                    outerProduct(D20, U20));

                basis = [this.normal(v1), new Vector(3), new Vector(3)];
                computeOrthogonalComplement3(1, basis, true);
                const U01 = Vector.fromArray(
                    [dot(basis[1], D10), dot(basis[2], D10)]);
                const U21 = Vector.fromArray(
                    [dot(basis[1], D21), dot(basis[2], D21)]);
                accumulate(this.mUTU[v1], outerProduct(U01, U01),
                    outerProduct(U21, U21));
                accumulate(this.mDTU[v1], outerProduct(D10, U01),
                    outerProduct(D21, U21));

                basis = [this.normal(v2), new Vector(3), new Vector(3)];
                computeOrthogonalComplement3(1, basis, true);
                const U02 = Vector.fromArray(
                    [dot(basis[1], D20), dot(basis[2], D20)]);
                const U12 = Vector.fromArray(
                    [dot(basis[1], D21), dot(basis[2], D21)]);
                accumulate(this.mUTU[v2], outerProduct(U02, U02),
                    outerProduct(U12, U12));
                accumulate(this.mDTU[v2], outerProduct(D20, U02),
                    outerProduct(D21, U12));
            }
        }

        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            const jacobian = multiplyAB(this.mDTU[i],
                inverse2x2(this.mUTU[i]).inverse);

            const basis: Vector[] = [
                Vector.fromArray(
                    [jacobian.get(0, 0), jacobian.get(1, 0), jacobian.get(2, 0)]),
                Vector.fromArray(
                    [jacobian.get(0, 1), jacobian.get(1, 1), jacobian.get(2, 1)]),
                new Vector(3)
            ];

            if (this.mDPDUs) {
                this.setDPDU(i, basis[0]);
            }
            if (this.mDPDVs) {
                this.setDPDV(i, basis[1]);
            }

            computeOrthogonalComplement3(2, basis, true);

            if (this.mNormals) {
                this.setNormal(i, basis[2]);
            }
            if (this.mTangents) {
                this.setTangent(i, basis[0]);
            }
            if (this.mBitangents) {
                this.setBitangent(i, basis[1]);
            }
        }
    }
}

// The port of 'M += A + B' for the running sums of the least-squares
// algorithm. The matrix M is modified in place.
function accumulate(M: Matrix, A: Matrix, B: Matrix): void {
    for (let i = 0; i < M.values.length; ++i) {
        M.values[i] += A.values[i] + B.values[i];
    }
}
