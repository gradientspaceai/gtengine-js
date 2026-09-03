// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RectanglePatchMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Create a mesh (x(u,v),y(u,v),z(u,v)) defined by the specified surface. It
// is required that surface.isRectangular() return true.
//
// Port notes: the upstream std::shared_ptr<ParametricSurface<3,Real>> becomes
// the surface object itself (the B57/B64 precedent). The default texture
// coordinates, which upstream stores in a std::vector and points mTCoords at,
// become a MeshChannel allocated by MeshChannel.allocate (the B77/B84
// precedent).

import { logAssert } from './Logger.js';
import { Mesh, MeshChannel, MeshDescription, MeshTopology } from './Mesh.js';
import { ParametricSurface } from './ParametricSurface.js';
import { Vector, normalize } from './Vector.js';
import { computeOrthogonalComplement3, unitCross } from './Vector3.js';

export class RectanglePatchMesh extends Mesh {
    protected mSurface: ParametricSurface | null;

    constructor(description: MeshDescription, surface: ParametricSurface) {
        super(description, [MeshTopology.RECTANGLE]);

        this.mSurface = surface;

        if (!this.mDescription.constructed) {
            // The logger system will report these errors in the Mesh
            // constructor.
            this.mSurface = null;
            return;
        }

        logAssert(this.mSurface !== null && this.mSurface.isRectangular(),
            'A nonnull rectangular surface is required.');
        logAssert(surface.getDimension() === 3,
            'The surface must be 3-dimensional.');

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
    getSurface(): ParametricSurface | null {
        return this.mSurface;
    }

    protected initializeTCoords(): void {
        const surface = this.mSurface as ParametricSurface;
        const description = this.mDescription;
        const uMin = surface.getUMin();
        const uDelta = (surface.getUMax() - uMin) / (description.numCols - 1);
        const vMin = surface.getVMin();
        const vDelta = (surface.getVMax() - vMin) / (description.numRows - 1);
        const tcoord = new Vector(2);
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            tcoord.values[1] = vMin + vDelta * r;
            for (let c = 0; c < description.numCols; ++c, ++i) {
                tcoord.values[0] = uMin + uDelta * c;
                this.setTCoord(i, tcoord);
            }
        }
    }

    protected initializePositions(): void {
        const surface = this.mSurface as ParametricSurface;
        const description = this.mDescription;
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            for (let c = 0; c < description.numCols; ++c, ++i) {
                const tcoord = this.tcoord(i);
                this.setPosition(i,
                    surface.getPosition(tcoord.values[0], tcoord.values[1]));
            }
        }
    }

    protected initializeNormals(): void {
        const surface = this.mSurface as ParametricSurface;
        const description = this.mDescription;
        const jet = surface.createJet();
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            for (let c = 0; c < description.numCols; ++c, ++i) {
                const tcoord = this.tcoord(i);
                surface.evaluate(tcoord.values[0], tcoord.values[1], 1, jet);
                normalize(jet[1], true);
                normalize(jet[2], true);
                this.setNormal(i, unitCross(jet[1], jet[2], true));
            }
        }
    }

    protected initializeFrame(): void {
        const surface = this.mSurface as ParametricSurface;
        const description = this.mDescription;
        const jet = surface.createJet();
        for (let r = 0, i = 0; r < description.numRows; ++r) {
            for (let c = 0; c < description.numCols; ++c, ++i) {
                const tcoord = this.tcoord(i);
                surface.evaluate(tcoord.values[0], tcoord.values[1], 1, jet);
                normalize(jet[1], true);
                normalize(jet[2], true);

                if (this.mDPDUs) {
                    this.setDPDU(i, jet[1]);
                }
                if (this.mDPDVs) {
                    this.setDPDV(i, jet[2]);
                }

                // The upstream call is ComputeOrthogonalComplement(2,
                // &values[1], true), which operates on values[1], values[2]
                // and values[3].
                const basis: Vector[] = [jet[1], jet[2], jet[3]];
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

    protected override updatePositions(): void {
        if (this.mSurface !== null) {
            this.initializePositions();
        }
    }

    protected override updateNormals(): void {
        if (this.mSurface !== null) {
            this.initializeNormals();
        }
    }

    protected override updateFrame(): void {
        if (this.mSurface !== null) {
            this.initializeFrame();
        }
    }
}
