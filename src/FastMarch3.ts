// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FastMarch3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The topic of fast marching methods are discussed in the book
//   Level Set Methods and Fast Marching Methods:
//     Evolving Interfaces in Computational Geometry, Fluid Mechanics,
//     Computer Vision, and Materials Science
//   J.A. Sethian,
//   Cambridge University Press, 1999
//
// Port notes:
//   - The two upstream constructors (per-voxel speeds versus a single
//     constant speed) are merged into one whose 'speeds' parameter is either
//     an array or a number, matching the merged FastMarch base constructor.
//   - GetBoundary's output parameter becomes the return value of
//     getBoundary().
//   - Upstream Iterate() calls mHeap.Remove(i, value) unconditionally. The
//     ported MinHeap.remove() returns null for an empty heap, so iterate()
//     returns immediately in that case.
//   - As in FastMarch2, ComputeTime solves the Eikonal equation on a unit
//     grid; the spacings are stored and exposed by the accessors but are
//     never used by the numerical method. The quirk is preserved (upstream
//     issue #52, FastMarch dead code).
//
// Upstream bug (fixed here): the upstream Initialize marks only the eight
// box vertices and the twelve box edges with zero speed, not the six box
// faces, even though the comment states that "boundary pixels are marked as
// zero speed to allow us to avoid having to process the boundary pixels
// separately during the iteration" and the 2D code does mark its entire
// border. Voxels in the interior of a face therefore stay "far" and are
// pulled into the marching front, after which their off-grid neighbor
// indices (i - mXYBound for a z = 0 voxel, for instance) wrap in the C++
// size_t arithmetic and read out of bounds. In TypeScript such an index
// yields undefined, and isTrial(undefined index) is true, which would make
// iterate() update a nonexistent heap record. The port marks every voxel of
// the six faces with zero speed, which is what the 2D code and the comment
// intend, and which restores the invariant that valid or trial voxels are
// strictly interior.

import { FastMarch } from './FastMarch';

export class FastMarch3 extends FastMarch {
    protected mXBound: number;
    protected mYBound: number;
    protected mZBound: number;
    protected mXYBound: number;
    protected mXBoundM1: number;
    protected mYBoundM1: number;
    protected mZBoundM1: number;
    protected mXSpacing: number;
    protected mYSpacing: number;
    protected mZSpacing: number;
    protected mInvXSpacing: number;
    protected mInvYSpacing: number;
    protected mInvZSpacing: number;

    constructor(xBound: number, yBound: number, zBound: number,
        xSpacing: number, ySpacing: number, zSpacing: number,
        seeds: readonly number[], speeds: readonly number[] | number) {
        super(xBound * yBound * zBound, seeds, speeds);

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mZBound = zBound;
        this.mXYBound = xBound * yBound;
        this.mXBoundM1 = xBound - 1;
        this.mYBoundM1 = yBound - 1;
        this.mZBoundM1 = zBound - 1;
        this.mXSpacing = xSpacing;
        this.mYSpacing = ySpacing;
        this.mZSpacing = zSpacing;
        this.mInvXSpacing = 1 / xSpacing;
        this.mInvYSpacing = 1 / ySpacing;
        this.mInvZSpacing = 1 / zSpacing;

        this.initialize();
    }

    // Member access.
    getXBound(): number {
        return this.mXBound;
    }

    getYBound(): number {
        return this.mYBound;
    }

    getZBound(): number {
        return this.mZBound;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    getYSpacing(): number {
        return this.mYSpacing;
    }

    getZSpacing(): number {
        return this.mZSpacing;
    }

    index(x: number, y: number, z: number): number {
        return x + this.mXBound * (y + this.mYBound * z);
    }

    // Voxel classification.
    override getBoundary(): number[] {
        const boundary: number[] = [];
        for (let i = 0; i < this.mQuantity; ++i) {
            if (this.isBoundary(i)) {
                boundary.push(i);
            }
        }
        return boundary;
    }

    override isBoundary(i: number): boolean {
        if (this.isValid(i) && !this.isTrial(i)) {
            if (this.isTrial(i - 1)
                || this.isTrial(i + 1)
                || this.isTrial(i - this.mXBound)
                || this.isTrial(i + this.mXBound)
                || this.isTrial(i - this.mXYBound)
                || this.isTrial(i + this.mXYBound)) {
                return true;
            }
        }
        return false;
    }

    // Run one step of the fast marching algorithm.
    override iterate(): void {
        // Remove the minimum trial value from the heap.
        const minimum = this.mHeap.remove();
        if (minimum === null) {
            return;
        }
        const i = minimum.key;

        // Promote the trial value to a known value. The value was negative
        // but is now nonnegative (the heap stores only nonnegative numbers).
        this.mTrials[i] = null;

        // All trial voxels must be updated. All far neighbors must become
        // trial voxels.
        this.updateNeighbor(i - 1);
        this.updateNeighbor(i + 1);
        this.updateNeighbor(i - this.mXBound);
        this.updateNeighbor(i + this.mXBound);
        this.updateNeighbor(i - this.mXYBound);
        this.updateNeighbor(i + this.mXYBound);
    }

    // The common body of the six neighbor updates of iterate(); a
    // behavior-preserving factoring of the repeated upstream blocks.
    private updateNeighbor(j: number): void {
        if (this.isTrial(j)) {
            this.computeTime(j);
            this.mHeap.update(this.mTrials[j], this.mTimes[j]);
        } else if (this.isFar(j)) {
            this.computeTime(j);
            this.mTrials[j] = this.mHeap.insert(j, this.mTimes[j]);
        }
    }

    // Called by the constructor.
    protected initialize(): void {
        // Boundary voxels are marked as zero speed to allow us to avoid
        // having to process the boundary voxels separately during the
        // iteration. See the file comment: upstream marks only the vertices
        // and edges of the box, the port marks the six faces (which include
        // the vertices and edges).
        let x: number, y: number, z: number, i: number;

        // faces z = 0 and z = zmax
        for (y = 0; y < this.mYBound; ++y) {
            for (x = 0; x < this.mXBound; ++x) {
                this.markZeroSpeed(this.index(x, y, 0));
                this.markZeroSpeed(this.index(x, y, this.mZBoundM1));
            }
        }

        // faces y = 0 and y = ymax
        for (z = 0; z < this.mZBound; ++z) {
            for (x = 0; x < this.mXBound; ++x) {
                this.markZeroSpeed(this.index(x, 0, z));
                this.markZeroSpeed(this.index(x, this.mYBoundM1, z));
            }
        }

        // faces x = 0 and x = xmax
        for (z = 0; z < this.mZBound; ++z) {
            for (y = 0; y < this.mYBound; ++y) {
                this.markZeroSpeed(this.index(0, y, z));
                this.markZeroSpeed(this.index(this.mXBoundM1, y, z));
            }
        }

        // Compute the first batch of trial voxels. These are voxels a grid
        // distance of one away from the seed voxels.
        for (z = 1; z < this.mZBoundM1; ++z) {
            for (y = 1; y < this.mYBoundM1; ++y) {
                for (x = 1; x < this.mXBoundM1; ++x) {
                    i = this.index(x, y, z);
                    if (this.isFar(i)) {
                        if (this.isKnown(i - 1)
                            || this.isKnown(i + 1)
                            || this.isKnown(i - this.mXBound)
                            || this.isKnown(i + this.mXBound)
                            || this.isKnown(i - this.mXYBound)
                            || this.isKnown(i + this.mXYBound)) {
                            this.computeTime(i);
                            this.mTrials[i] = this.mHeap.insert(i, this.mTimes[i]);
                        }
                    }
                }
            }
        }
    }

    // A voxel whose crossing time is final (valid but not a trial voxel).
    private isKnown(i: number): boolean {
        return this.isValid(i) && !this.isTrial(i);
    }

    private markZeroSpeed(i: number): void {
        this.mInvSpeeds[i] = Number.MAX_VALUE;
        this.mTimes[i] = -Number.MAX_VALUE;
    }

    // Called by iterate().
    protected computeTime(i: number): void {
        let hasXTerm: boolean;
        let xConst: number;
        if (this.isValid(i - 1)) {
            hasXTerm = true;
            xConst = this.mTimes[i - 1];
            if (this.isValid(i + 1)) {
                if (this.mTimes[i + 1] < xConst) {
                    xConst = this.mTimes[i + 1];
                }
            }
        } else if (this.isValid(i + 1)) {
            hasXTerm = true;
            xConst = this.mTimes[i + 1];
        } else {
            hasXTerm = false;
            xConst = 0;
        }

        let hasYTerm: boolean;
        let yConst: number;
        if (this.isValid(i - this.mXBound)) {
            hasYTerm = true;
            yConst = this.mTimes[i - this.mXBound];
            if (this.isValid(i + this.mXBound)) {
                if (this.mTimes[i + this.mXBound] < yConst) {
                    yConst = this.mTimes[i + this.mXBound];
                }
            }
        } else if (this.isValid(i + this.mXBound)) {
            hasYTerm = true;
            yConst = this.mTimes[i + this.mXBound];
        } else {
            hasYTerm = false;
            yConst = 0;
        }

        let hasZTerm: boolean;
        let zConst: number;
        if (this.isValid(i - this.mXYBound)) {
            hasZTerm = true;
            zConst = this.mTimes[i - this.mXYBound];
            if (this.isValid(i + this.mXYBound)) {
                if (this.mTimes[i + this.mXYBound] < zConst) {
                    zConst = this.mTimes[i + this.mXYBound];
                }
            }
        } else if (this.isValid(i + this.mXYBound)) {
            hasZTerm = true;
            zConst = this.mTimes[i + this.mXYBound];
        } else {
            hasZTerm = false;
            zConst = 0;
        }

        let sum: number, diff: number, discr: number;

        if (hasXTerm) {
            if (hasYTerm) {
                if (hasZTerm) {
                    // xyz
                    sum = xConst + yConst + zConst;
                    discr = 3 * this.mInvSpeeds[i] * this.mInvSpeeds[i];
                    diff = xConst - yConst;
                    discr -= diff * diff;
                    diff = xConst - zConst;
                    discr -= diff * diff;
                    diff = yConst - zConst;
                    discr -= diff * diff;
                    if (discr >= 0) {
                        // The quadratic equation has a real-valued solution.
                        // Choose the largest positive root for the crossing
                        // time.
                        this.mTimes[i] = (sum + Math.sqrt(discr)) / 3;
                    } else {
                        // The quadratic equation does not have a real-valued
                        // solution. This can happen when the speed is so
                        // large that the time gradient has very small length,
                        // which means that the time has not changed
                        // significantly from the neighbors to the current
                        // voxel. Just choose the maximum time of the
                        // neighbors.
                        this.mTimes[i] = xConst;
                        if (yConst > this.mTimes[i]) {
                            this.mTimes[i] = yConst;
                        }
                        if (zConst > this.mTimes[i]) {
                            this.mTimes[i] = zConst;
                        }
                    }
                } else {
                    // xy
                    sum = xConst + yConst;
                    diff = xConst - yConst;
                    discr = 2 * this.mInvSpeeds[i] * this.mInvSpeeds[i] - diff * diff;
                    if (discr >= 0) {
                        this.mTimes[i] = 0.5 * (sum + Math.sqrt(discr));
                    } else {
                        this.mTimes[i] = (diff >= 0 ? xConst : yConst);
                    }
                }
            } else {
                if (hasZTerm) {
                    // xz
                    sum = xConst + zConst;
                    diff = xConst - zConst;
                    discr = 2 * this.mInvSpeeds[i] * this.mInvSpeeds[i] - diff * diff;
                    if (discr >= 0) {
                        this.mTimes[i] = 0.5 * (sum + Math.sqrt(discr));
                    } else {
                        this.mTimes[i] = (diff >= 0 ? xConst : zConst);
                    }
                } else {
                    // x
                    this.mTimes[i] = this.mInvSpeeds[i] + xConst;
                }
            }
        } else {
            if (hasYTerm) {
                if (hasZTerm) {
                    // yz
                    sum = yConst + zConst;
                    diff = yConst - zConst;
                    discr = 2 * this.mInvSpeeds[i] * this.mInvSpeeds[i] - diff * diff;
                    if (discr >= 0) {
                        this.mTimes[i] = 0.5 * (sum + Math.sqrt(discr));
                    } else {
                        this.mTimes[i] = (diff >= 0 ? yConst : zConst);
                    }
                } else {
                    // y
                    this.mTimes[i] = this.mInvSpeeds[i] + yConst;
                }
            } else {
                if (hasZTerm) {
                    // z
                    this.mTimes[i] = this.mInvSpeeds[i] + zConst;
                }
                // else: Assert: The voxel must have at least one valid
                // neighbor.
            }
        }
    }
}
