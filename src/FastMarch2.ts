// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) FastMarch2.h
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
//   - The two upstream constructors (per-pixel speeds versus a single
//     constant speed) are merged into one whose 'speeds' parameter is either
//     an array or a number, matching the merged FastMarch base constructor.
//   - GetBoundary's output parameter becomes the return value of
//     getBoundary().
//   - Upstream Iterate() calls mHeap.Remove(i, value) unconditionally. The
//     ported MinHeap.remove() returns null for an empty heap, so iterate()
//     returns immediately in that case rather than dereferencing null.
//   - The upstream ComputeTime solves the Eikonal equation on a unit grid;
//     mXSpacing/mYSpacing (and their reciprocals) are stored and exposed by
//     the accessors but are never used by the numerical method. That quirk is
//     preserved (see upstream issue #52 for the FastMarch dead-code family).
//   - The neighbor indices i-1, i+1, i-xBound and i+xBound are always in
//     range: Initialize marks every border pixel with zero speed, so a valid
//     or trial pixel is always strictly interior. Upstream relies on the same
//     invariant (its size_t arithmetic would otherwise wrap).

import { FastMarch } from './FastMarch';

export class FastMarch2 extends FastMarch {
    protected mXBound: number;
    protected mYBound: number;
    protected mXBoundM1: number;
    protected mYBoundM1: number;
    protected mXSpacing: number;
    protected mYSpacing: number;
    protected mInvXSpacing: number;
    protected mInvYSpacing: number;

    constructor(xBound: number, yBound: number, xSpacing: number, ySpacing: number,
        seeds: readonly number[], speeds: readonly number[] | number) {
        super(xBound * yBound, seeds, speeds);

        this.mXBound = xBound;
        this.mYBound = yBound;
        this.mXBoundM1 = xBound - 1;
        this.mYBoundM1 = yBound - 1;
        this.mXSpacing = xSpacing;
        this.mYSpacing = ySpacing;
        this.mInvXSpacing = 1 / xSpacing;
        this.mInvYSpacing = 1 / ySpacing;

        this.initialize();
    }

    // Member access.
    getXBound(): number {
        return this.mXBound;
    }

    getYBound(): number {
        return this.mYBound;
    }

    getXSpacing(): number {
        return this.mXSpacing;
    }

    getYSpacing(): number {
        return this.mYSpacing;
    }

    index(x: number, y: number): number {
        return x + this.mXBound * y;
    }

    // Pixel classification.
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
                || this.isTrial(i + this.mXBound)) {
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

        // All trial pixels must be updated. All far neighbors must become
        // trial pixels.
        this.updateNeighbor(i - 1);
        this.updateNeighbor(i + 1);
        this.updateNeighbor(i - this.mXBound);
        this.updateNeighbor(i + this.mXBound);
    }

    // The common body of the four neighbor updates of iterate(); a
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
        // Boundary pixels are marked as zero speed to allow us to avoid
        // having to process the boundary pixels separately during the
        // iteration.
        let x: number, y: number, i: number;

        // vertex (0,0)
        this.markZeroSpeed(this.index(0, 0));

        // vertex (xmax,0)
        this.markZeroSpeed(this.index(this.mXBoundM1, 0));

        // vertex (0,ymax)
        this.markZeroSpeed(this.index(0, this.mYBoundM1));

        // vertex (xmax,ymax)
        this.markZeroSpeed(this.index(this.mXBoundM1, this.mYBoundM1));

        // edges (x,0) and (x,ymax)
        for (x = 0; x < this.mXBound; ++x) {
            this.markZeroSpeed(this.index(x, 0));
            this.markZeroSpeed(this.index(x, this.mYBoundM1));
        }

        // edges (0,y) and (xmax,y)
        for (y = 0; y < this.mYBound; ++y) {
            this.markZeroSpeed(this.index(0, y));
            this.markZeroSpeed(this.index(this.mXBoundM1, y));
        }

        // Compute the first batch of trial pixels. These are pixels a grid
        // distance of one away from the seed pixels.
        for (y = 1; y < this.mYBoundM1; ++y) {
            for (x = 1; x < this.mXBoundM1; ++x) {
                i = this.index(x, y);
                if (this.isFar(i)) {
                    if (this.isKnown(i - 1)
                        || this.isKnown(i + 1)
                        || this.isKnown(i - this.mXBound)
                        || this.isKnown(i + this.mXBound)) {
                        this.computeTime(i);
                        this.mTrials[i] = this.mHeap.insert(i, this.mTimes[i]);
                    }
                }
            }
        }
    }

    // A pixel whose crossing time is final (valid but not a trial pixel).
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

        if (hasXTerm) {
            if (hasYTerm) {
                const sum = xConst + yConst;
                const diff = xConst - yConst;
                const discr = 2 * this.mInvSpeeds[i] * this.mInvSpeeds[i] - diff * diff;
                if (discr >= 0) {
                    // The quadratic equation has a real-valued solution.
                    // Choose the largest positive root for the crossing time.
                    this.mTimes[i] = 0.5 * (sum + Math.sqrt(discr));
                } else {
                    // The quadratic equation does not have a real-valued
                    // solution. This can happen when the speed is so large
                    // that the time gradient has very small length, which
                    // means that the time has not changed significantly from
                    // the neighbors to the current pixel. Just choose the
                    // maximum time of the neighbors.
                    this.mTimes[i] = (diff >= 0 ? xConst : yConst);
                }
            } else {
                // The equation is linear.
                this.mTimes[i] = this.mInvSpeeds[i] + xConst;
            }
        } else if (hasYTerm) {
            // The equation is linear.
            this.mTimes[i] = this.mInvSpeeds[i] + yConst;
        }
        // else: Assert: The pixel must have at least one known neighbor.
    }
}
