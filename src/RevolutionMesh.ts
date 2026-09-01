// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) RevolutionMesh.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Port notes: the upstream std::shared_ptr<ParametricCurve<2,Real>> becomes
// the curve object itself (the B57/B64 precedent). The upstream mTSampler
// std::function becomes a private method. The default texture coordinates,
// which upstream stores in a std::vector and points mTCoords at, become a
// MeshChannel allocated by MeshChannel.allocate.

import { GTE_C_TWO_PI } from './Constants';
import { logAssert } from './Logger';
import { Mesh, MeshChannel, MeshDescription, MeshTopology } from './Mesh';
import { ParametricCurve } from './ParametricCurve';
import { Vector } from './Vector';

export class RevolutionMesh extends Mesh {
    private mCurve: ParametricCurve | null;
    private mSampleByArcLength: boolean;
    private mCosAngle: number[];
    private mSinAngle: number[];
    private mTSamplerFactor: number;
    private mSamples: Vector[];

    // The axis of revolution is the z-axis. The curve of revolution is
    // p(t) = (x(t),z(t)), where t in [tmin,tmax], x(t) > 0 for t in
    // (tmin,tmax), x(tmin) >= 0, and x(tmax) >= 0. The values tmin and tmax
    // are those for the curve object passed to the constructor. The curve
    // must be non-self-intersecting, except possibly at its endpoints. The
    // curve is closed when p(tmin) = p(tmax), in which case the surface of
    // revolution has torus topology. The curve is open when
    // p(tmin) != p(tmax). For an open curve, define x0 = x(tmin) and
    // x1 = x(tmax). The surface has cylinder topology when x0 > 0 and x1 > 0,
    // disk topology when exactly one of x0 or x1 is zero, or sphere topology
    // when x0 and x1 are both zero. However, to simplify the design, the mesh
    // is always built using cylinder topology. The row samples correspond to
    // curve points and the column samples correspond to the points on the
    // circles of revolution.
    constructor(description: MeshDescription, curve: ParametricCurve,
        sampleByArcLength: boolean = false) {
        super(description, [MeshTopology.CYLINDER, MeshTopology.TORUS,
            MeshTopology.DISK, MeshTopology.SPHERE]);

        this.mCurve = curve;
        this.mSampleByArcLength = sampleByArcLength;
        this.mCosAngle = [];
        this.mSinAngle = [];
        this.mTSamplerFactor = 0;
        this.mSamples = [];

        if (!this.mDescription.constructed) {
            // The logger system will report these errors in the Mesh
            // constructor.
            this.mCurve = null;
            return;
        }

        logAssert(this.mCurve !== null, 'A nonnull revolution curve is required.');
        logAssert(curve.getDimension() === 2,
            'The revolution curve must be 2-dimensional.');

        // The four supported topologies all wrap around in the column
        // direction.
        const numCols = this.mDescription.numCols;
        this.mCosAngle = new Array<number>(numCols + 1).fill(0);
        this.mSinAngle = new Array<number>(numCols + 1).fill(0);
        const invRadialSamples = 1 / numCols;
        for (let c = 0; c < numCols; ++c) {
            const angle = c * invRadialSamples * GTE_C_TWO_PI;
            this.mCosAngle[c] = Math.cos(angle);
            this.mSinAngle[c] = Math.sin(angle);
        }
        this.mCosAngle[numCols] = this.mCosAngle[0];
        this.mSinAngle[numCols] = this.mSinAngle[0];

        this.createSampler();

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
    getCurve(): ParametricCurve | null {
        return this.mCurve;
    }

    isSampleByArcLength(): boolean {
        return this.mSampleByArcLength;
    }

    private createSampler(): void {
        const topology = this.mDescription.topology;
        let numSamples = 0;
        if (topology === MeshTopology.CYLINDER || topology === MeshTopology.TORUS) {
            numSamples = this.mDescription.rMax + 1;
        }
        else if (topology === MeshTopology.DISK) {
            numSamples = this.mDescription.rMax + 2;
        }
        else if (topology === MeshTopology.SPHERE) {
            numSamples = this.mDescription.rMax + 3;
        }

        this.mSamples = new Array<Vector>(numSamples);
        for (let i = 0; i < numSamples; ++i) {
            this.mSamples[i] = new Vector(3);
        }

        const curve = this.mCurve as ParametricCurve;
        const invDenom = 1 / (numSamples - 1);
        if (this.mSampleByArcLength) {
            this.mTSamplerFactor = curve.getTotalLength() * invDenom;
        }
        else {
            this.mTSamplerFactor = (curve.getTMax() - curve.getTMin()) * invDenom;
        }
    }

    // The port of the upstream mTSampler std::function.
    private tSampler(i: number): number {
        const curve = this.mCurve as ParametricCurve;
        if (this.mSampleByArcLength) {
            return curve.getTime(i * this.mTSamplerFactor);
        }
        else {
            return curve.getTMin() + i * this.mTSamplerFactor;
        }
    }

    private initializeTCoords(): void {
        const description = this.mDescription;
        const tcoord = new Vector(2);

        switch (description.topology) {
            case MeshTopology.CYLINDER: {
                for (let r = 0, i = 0; r < description.numRows; ++r) {
                    tcoord.values[1] = r / (description.numRows - 1);
                    for (let c = 0; c <= description.numCols; ++c, ++i) {
                        tcoord.values[0] = c / description.numCols;
                        this.setTCoord(i, tcoord);
                    }
                }
                break;
            }
            case MeshTopology.TORUS: {
                for (let r = 0, i = 0; r <= description.numRows; ++r) {
                    tcoord.values[1] = r / description.numRows;
                    for (let c = 0; c <= description.numCols; ++c, ++i) {
                        tcoord.values[0] = c / description.numCols;
                        this.setTCoord(i, tcoord);
                    }
                }
                break;
            }
            case MeshTopology.DISK: {
                const origin = Vector.fromArray([0.5, 0.5]);
                let i = 0;
                for (let r = 0; r < description.numRows; ++r) {
                    let radius = (r + 1) / (2 * description.numRows);
                    radius = Math.min(radius, 0.5);
                    for (let c = 0; c <= description.numCols; ++c, ++i) {
                        const angle = GTE_C_TWO_PI * c / description.numCols;
                        this.setTCoord(i, Vector.fromArray(
                            [radius * Math.cos(angle), radius * Math.sin(angle)]));
                    }
                }
                this.setTCoord(i, origin);
                break;
            }
            case MeshTopology.SPHERE: {
                let i = 0;
                for (let r = 0; r < description.numRows; ++r) {
                    tcoord.values[1] = r / (description.numRows - 1);
                    for (let c = 0; c <= description.numCols; ++c, ++i) {
                        tcoord.values[0] = c / description.numCols;
                        this.setTCoord(i, tcoord);
                    }
                }
                this.setTCoord(i++, Vector.fromArray([0.5, 0]));
                this.setTCoord(i, Vector.fromArray([0.5, 1]));
                break;
            }
            default:
                // Invalid topology is reported by the Mesh constructor, so
                // there is no need to log a message here.
                break;
        }
    }

    protected override updatePositions(): void {
        if (this.mCurve === null) {
            return;
        }

        const numSamples = this.mSamples.length;
        for (let i = 0; i < numSamples; ++i) {
            const t = this.tSampler(i);
            const position = this.mCurve.getPosition(t);
            this.mSamples[i].values[0] = position.values[0];
            this.mSamples[i].values[1] = 0;
            this.mSamples[i].values[2] = position.values[1];
        }

        switch (this.mDescription.topology) {
            case MeshTopology.CYLINDER:
                this.updateCylinderPositions();
                break;
            case MeshTopology.TORUS:
                this.updateTorusPositions();
                break;
            case MeshTopology.DISK:
                this.updateDiskPositions();
                break;
            case MeshTopology.SPHERE:
                this.updateSpherePositions();
                break;
            default:
                break;
        }
    }

    private updateCylinderPositions(): void {
        for (let r = 0, i = 0; r <= this.mDescription.rMax; ++r) {
            const radius = this.mSamples[r].values[0];
            for (let c = 0; c <= this.mDescription.cMax; ++c, ++i) {
                this.setPosition(i, Vector.fromArray([radius * this.mCosAngle[c],
                    radius * this.mSinAngle[c], this.mSamples[r].values[2]]));
            }
        }
    }

    private updateTorusPositions(): void {
        for (let r = 0, i = 0; r <= this.mDescription.rMax; ++r) {
            const radius = this.mSamples[r].values[0];
            for (let c = 0; c <= this.mDescription.cMax; ++c, ++i) {
                this.setPosition(i, Vector.fromArray([radius * this.mCosAngle[c],
                    radius * this.mSinAngle[c], this.mSamples[r].values[2]]));
            }
        }
    }

    private updateDiskPositions(): void {
        for (let r = 0, rp1 = 1, i = 0; r <= this.mDescription.rMax; ++r, ++rp1) {
            const radius = this.mSamples[rp1].values[0];
            for (let c = 0; c <= this.mDescription.cMax; ++c, ++i) {
                this.setPosition(i, Vector.fromArray([radius * this.mCosAngle[c],
                    radius * this.mSinAngle[c], this.mSamples[rp1].values[2]]));
            }
        }

        this.setPosition(this.mDescription.numVertices - 1,
            Vector.fromArray([0, 0, this.mSamples[0].values[2]]));
    }

    private updateSpherePositions(): void {
        for (let r = 0, rp1 = 1, i = 0; r <= this.mDescription.rMax; ++r, ++rp1) {
            const radius = this.mSamples[rp1].values[0];
            for (let c = 0; c <= this.mDescription.cMax; ++c, ++i) {
                this.setPosition(i, Vector.fromArray([radius * this.mCosAngle[c],
                    radius * this.mSinAngle[c], this.mSamples[rp1].values[2]]));
            }
        }

        this.setPosition(this.mDescription.numVertices - 2,
            Vector.fromArray([0, 0, this.mSamples[0].values[2]]));
        this.setPosition(this.mDescription.numVertices - 1,
            Vector.fromArray([0, 0, this.mSamples[this.mSamples.length - 1].values[2]]));
    }
}
