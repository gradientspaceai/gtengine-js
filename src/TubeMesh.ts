// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) TubeMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Create a mesh (x(u,v),y(u,v),z(u,v)) defined by the specified medial curve
// and radial function. The mesh has torus topology when 'closed' is true and
// has cylinder topology when 'closed' is false. The client is responsible for
// setting the topology correctly in the 'description' input. The rows
// correspond to medial samples and the columns correspond to radial samples.
// The medial curve is sampled according to its natural t-parameter when
// 'sampleByArcLength' is false; otherwise, it is sampled uniformly in
// arclength. TODO (upstream): Allow TORUS and remove the 'closed' input.
//
// Port notes: the upstream std::shared_ptr<ParametricCurve<3,Real>> becomes
// the curve object itself (the B57/B64 precedent). The upstream mTSampler and
// mFSampler std::function objects become private methods. The default texture
// coordinates, which upstream stores in a std::vector and points mTCoords at,
// become a MeshChannel allocated by MeshChannel.allocate.

import { GTE_C_TWO_PI } from './Constants.js';
import { FrenetFrame3 } from './FrenetFrame.js';
import { logAssert } from './Logger.js';
import { Mesh, MeshChannel, MeshDescription, MeshTopology } from './Mesh.js';
import { ParametricCurve } from './ParametricCurve.js';
import { Vector, add, mul } from './Vector.js';
import { unitCross } from './Vector3.js';

// The port of the upstream std::array<Vector3<Real>, 4> frame, which is
// (position, tangent, normal, binormal).
export interface TubeMeshFrame {
    position: Vector;
    tangent: Vector;
    normal: Vector;
    binormal: Vector;
}

export class TubeMesh extends Mesh {
    private mMedial: ParametricCurve | null;
    private mRadial: (t: number) => number;
    private mClosed: boolean;
    private mSampleByArcLength: boolean;
    private mUpVector: Vector;
    private mCosAngle: number[];
    private mSinAngle: number[];
    private mTSamplerFactor: number;
    private mFrenet: FrenetFrame3 | null;

    constructor(description: MeshDescription, medial: ParametricCurve,
        radial: (t: number) => number, closed: boolean,
        sampleByArcLength: boolean, upVector: Vector) {
        super(description, [MeshTopology.CYLINDER]);

        this.mMedial = medial;
        this.mRadial = radial;
        this.mClosed = closed;
        this.mSampleByArcLength = sampleByArcLength;
        this.mUpVector = upVector.clone();
        this.mCosAngle = [];
        this.mSinAngle = [];
        this.mTSamplerFactor = 0;
        this.mFrenet = null;

        if (!this.mDescription.constructed) {
            // The logger system will report these errors in the Mesh
            // constructor.
            this.mMedial = null;
            return;
        }

        logAssert(this.mMedial !== null, 'A nonnull medial curve is required.');
        logAssert(medial.getDimension() === 3,
            'The medial curve must be 3-dimensional.');

        const numCols = this.mDescription.numCols;
        this.mCosAngle = new Array<number>(numCols).fill(0);
        this.mSinAngle = new Array<number>(numCols).fill(0);
        const invRadialSamples = 1 / (numCols - 1);
        for (let i = 0; i < numCols - 1; ++i) {
            const angle = i * invRadialSamples * GTE_C_TWO_PI;
            this.mCosAngle[i] = Math.cos(angle);
            this.mSinAngle[i] = Math.sin(angle);
        }
        this.mCosAngle[numCols - 1] = this.mCosAngle[0];
        this.mSinAngle[numCols - 1] = this.mSinAngle[0];

        let invDenom: number;
        if (this.mClosed) {
            invDenom = 1 / this.mDescription.numRows;
        }
        else {
            invDenom = 1 / (this.mDescription.numRows - 1);
        }

        if (this.mSampleByArcLength) {
            this.mTSamplerFactor = medial.getTotalLength() * invDenom;
        }
        else {
            this.mTSamplerFactor = (medial.getTMax() - medial.getTMin()) * invDenom;
        }

        if (!isZeroVector(this.mUpVector)) {
            this.mFrenet = null;
        }
        else {
            this.mFrenet = new FrenetFrame3(medial);
        }

        if (!this.mTCoords) {
            this.mTCoords = MeshChannel.allocate(this.mDescription.numVertices, 2);

            this.mDescription.allowUpdateFrame =
                this.mDescription.wantDynamicTangentSpaceUpdate;
            if (this.mDescription.allowUpdateFrame) {
                if (!this.mDescription.hasTangentSpaceVectors) {
                    this.mDescription.allowUpdateFrame = false;
                }

                if (!this.mNormals) {
                    this.mDescription.allowUpdateFrame = false;
                }
            }
        }

        this.computeIndices();
        this.initializeTCoords();
        this.updatePositions();
        if (this.mDescription.allowUpdateFrame) {
            this.updateFrame();
        }
        else if (this.mNormals) {
            this.updateNormals();
        }
    }

    // Member access.
    getMedial(): ParametricCurve | null {
        return this.mMedial;
    }

    getRadial(): (t: number) => number {
        return this.mRadial;
    }

    isClosed(): boolean {
        return this.mClosed;
    }

    isSampleByArcLength(): boolean {
        return this.mSampleByArcLength;
    }

    getUpVector(): Vector {
        return this.mUpVector;
    }

    // The port of the upstream mTSampler std::function.
    private tSampler(row: number): number {
        const medial = this.mMedial as ParametricCurve;
        if (this.mSampleByArcLength) {
            return medial.getTime(row * this.mTSamplerFactor);
        }
        else {
            return medial.getTMin() + row * this.mTSamplerFactor;
        }
    }

    // The port of the upstream mFSampler std::function.
    private fSampler(t: number): TubeMeshFrame {
        const medial = this.mMedial as ParametricCurve;
        if (this.mFrenet === null) {
            const position = medial.getPosition(t);
            const tangent = medial.getTangent(t);
            const binormal = unitCross(tangent, this.mUpVector);
            const normal = unitCross(binormal, tangent);
            return { position, tangent, normal, binormal };
        }
        else {
            return this.mFrenet.compute(t);
        }
    }

    private initializeTCoords(): void {
        const description = this.mDescription;
        const tcoord = new Vector(2);
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            tcoord.values[1] = r / description.rMax;
            for (let c = 0; c <= description.numCols; ++c, ++i) {
                tcoord.values[0] = c / description.numCols;
                this.setTCoord(i, tcoord);
            }
        }
    }

    protected override updatePositions(): void {
        if (this.mMedial === null) {
            return;
        }

        const description = this.mDescription;
        let row: number, col: number, v: number, save: number;
        for (row = 0, v = 0; row < description.numRows; ++row, ++v) {
            const t = this.tSampler(row);
            const radius = this.mRadial(t);
            const frame = this.fSampler(t);
            for (col = 0, save = v; col < description.numCols; ++col, ++v) {
                this.setPosition(v, add(frame.position, mul(
                    add(mul(frame.normal, this.mCosAngle[col]),
                        mul(frame.binormal, this.mSinAngle[col])), radius)));
            }
            this.setPosition(v, this.position(save));
        }

        if (this.mClosed) {
            // Upstream TubeMesh.h computes the index of the last-row vertex
            // as 'col + numCols * (numRows - 1)'. The row stride of a
            // CYLINDER mesh is rIncrement = numCols + 1, because each row has
            // an extra vertex that duplicates the first one to close the
            // seam. The upstream expression therefore addresses vertices of
            // an interior row (for numRows > 1) and corrupts the mesh, so the
            // port uses the correct stride. The port also copies the trailing
            // seam vertex of the row (col <= numCols); upstream stops at
            // col < numCols, which leaves the last vertex of the final row
            // with the stale position computed from the wrapped t-value.
            const rIncrement = description.rIncrement;
            for (col = 0; col <= description.numCols; ++col) {
                const i0 = col;
                const i1 = col + rIncrement * (description.numRows - 1);
                this.setPosition(i1, this.position(i0));
            }
        }
    }
}

// The port of the upstream 'mUpVector != Vector3<Real>::Zero()' test.
function isZeroVector(v: Vector): boolean {
    for (let i = 0; i < v.size; ++i) {
        if (v.values[i] !== 0) {
            return false;
        }
    }
    return true;
}
