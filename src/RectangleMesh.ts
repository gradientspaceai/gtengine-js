// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RectangleMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Create a mesh that tessellates a rectangle. The rows correspond to the
// second rectangle axis and the columns correspond to the first.
//
// Port notes: the default texture coordinates, which upstream stores in a
// std::vector and points mTCoords at, become a MeshChannel allocated by
// MeshChannel.allocate (the B77/B84 precedent set by TubeMesh and
// RevolutionMesh).

import { Mesh, MeshChannel, MeshDescription, MeshTopology } from './Mesh.js';
import { Rectangle } from './Rectangle.js';
import { logAssert } from './Logger.js';
import { Vector, add, mul } from './Vector.js';
import { unitCross } from './Vector3.js';

export class RectangleMesh extends Mesh {
    protected mRectangle: Rectangle;

    constructor(description: MeshDescription, rectangle: Rectangle) {
        super(description, [MeshTopology.RECTANGLE]);

        logAssert(rectangle.dimension === 3,
            'RectangleMesh: a 3-dimensional rectangle is required.');
        this.mRectangle = rectangle.clone();

        if (!this.mDescription.constructed) {
            // The logger system will report these errors in the Mesh
            // constructor.
            return;
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
        this.initializePositions();
        if (this.mDescription.allowUpdateFrame) {
            this.initializeFrame();
        }
        else if (this.mNormals) {
            this.initializeNormals();
        }
    }

    // Member access.
    getRectangle(): Rectangle {
        return this.mRectangle;
    }

    protected initializeTCoords(): void {
        const description = this.mDescription;
        const tcoord = new Vector(2);
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            tcoord.values[1] = r / (description.numRows - 1);
            for (let c = 0; c < description.numCols; ++c, ++i) {
                tcoord.values[0] = c / (description.numCols - 1);
                this.setTCoord(i, tcoord);
            }
        }
    }

    protected initializePositions(): void {
        const description = this.mDescription;
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            for (let c = 0; c < description.numCols; ++c, ++i) {
                const tcoord = this.tcoord(i);
                const w0 = (2 * tcoord.values[0] - 1) * this.mRectangle.extent.values[0];
                const w1 = (2 * tcoord.values[1] - 1) * this.mRectangle.extent.values[1];
                this.setPosition(i, add(add(this.mRectangle.center,
                    mul(this.mRectangle.axis[0], w0)),
                    mul(this.mRectangle.axis[1], w1)));
            }
        }
    }

    protected initializeNormals(): void {
        const normal = unitCross(this.mRectangle.axis[0], this.mRectangle.axis[1]);
        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            this.setNormal(i, normal);
        }
    }

    protected initializeFrame(): void {
        // Port fix for an upstream bug. Upstream hardcodes
        //   tangent = (1,0,0), bitangent = (0,1,0)
        // for every vertex, independent of the rectangle, under its own
        // "TODO: Are tangent and bitangent correct?" and a stale comment
        // "bitangent = Cross(normal,tangent)" that the code does not
        // implement. For any rectangle whose axes are not the standard frame,
        // those vectors are not in the tangent plane and are not orthogonal to
        // the normal, so the tangent/bitangent/dpdu/dpdv channels are
        // corrupted. The rectangle's own axes are orthonormal and span its
        // plane, so they are the frame the stale comment intends: the tangent
        // is axis[0], the bitangent is axis[1], and the normal is their cross
        // product (the same value upstream already writes to the normal
        // channel).
        const tangent = this.mRectangle.axis[0];
        const bitangent = this.mRectangle.axis[1];
        const normal = unitCross(tangent, bitangent);
        for (let i = 0; i < this.mDescription.numVertices; ++i) {
            if (this.mNormals) {
                this.setNormal(i, normal);
            }

            if (this.mTangents) {
                this.setTangent(i, tangent);
            }

            if (this.mBitangents) {
                this.setBitangent(i, bitangent);
            }

            if (this.mDPDUs) {
                this.setDPDU(i, tangent);
            }

            if (this.mDPDVs) {
                this.setDPDV(i, bitangent);
            }
        }
    }
}
