// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) PdeFilter.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The abstract base class for all PDE-based image filters (PdeFilter1,
// PdeFilter2, PdeFilter3 build on this).
//
// Port notes: the nested enum PdeFilter::ScaleType is exported as
// PdeFilterScaleType because src/index.ts star-exports every file. The
// template parameter Real is number (IEEE double).

export enum PdeFilterScaleType {
    // The data is processed as is.
    NONE,

    // The data range is d in [min,max].  The scaled values are d'.

    // d' = (d-min)/(max-min) in [0,1]
    UNIT,

    // d' = -1 + 2*(d-min)/(max-min) in [-1,1]
    SYMMETRIC,

    // max > -min:  d' = d/max in [min/max,1]
    // max < -min:  d' = -d/min in [-1,-max/min]
    PRESERVE_ZERO
}

export abstract class PdeFilter {
    // The number of image elements.
    protected mQuantity: number;

    // When set to Number.MAX_VALUE (upstream std::numeric_limits<Real>::max()),
    // Neumann conditions are in use (zero-valued derivatives on the image
    // border). Dirichlet conditions are used, otherwise (image is constant on
    // the border).
    protected mBorderValue: number;

    // These members store how the image data was transformed during the
    // constructor call.
    protected mScaleType: PdeFilterScaleType;
    protected mMin: number;
    protected mOffset: number;
    protected mScale: number;

    // The time step for the PDE solver.  The stability of an algorithm
    // depends on the magnitude of the time step, but the magnitude itself
    // depends on the algorithm.
    protected mTimeStep: number;

    protected constructor(quantity: number, data: ArrayLike<number>,
        borderValue: number, scaleType: PdeFilterScaleType) {
        this.mQuantity = quantity;
        this.mBorderValue = borderValue;
        this.mScaleType = scaleType;
        this.mMin = 0;
        this.mOffset = 0;
        this.mScale = 0;
        this.mTimeStep = 0;

        let maxValue = data[0];
        this.mMin = maxValue;
        for (let i = 1; i < this.mQuantity; i++) {
            const value = data[i];
            if (value < this.mMin) {
                this.mMin = value;
            } else if (value > maxValue) {
                maxValue = value;
            }
        }

        if (this.mMin !== maxValue) {
            switch (this.mScaleType) {
                case PdeFilterScaleType.NONE:
                    this.mOffset = 0;
                    this.mScale = 1;
                    break;
                case PdeFilterScaleType.UNIT:
                    this.mOffset = 0;
                    this.mScale = 1 / (maxValue - this.mMin);
                    break;
                case PdeFilterScaleType.SYMMETRIC:
                    this.mOffset = -1;
                    this.mScale = 2 / (maxValue - this.mMin);
                    break;
                case PdeFilterScaleType.PRESERVE_ZERO:
                    this.mOffset = 0;
                    this.mScale = (maxValue >= -this.mMin ? 1 / maxValue : -1 / this.mMin);
                    this.mMin = 0;
                    break;
            }
        } else {
            this.mOffset = 0;
            this.mScale = 1;
        }
    }

    // Member access.
    getQuantity(): number {
        return this.mQuantity;
    }

    getBorderValue(): number {
        return this.mBorderValue;
    }

    getScaleType(): PdeFilterScaleType {
        return this.mScaleType;
    }

    // Access to the time step for the PDE solver.
    setTimeStep(timeStep: number): void {
        this.mTimeStep = timeStep;
    }

    getTimeStep(): number {
        return this.mTimeStep;
    }

    // This function executes one iteration of the filter.  It calls
    // onPreUpdate, onUpdate and onPostUpdate, in that order.
    update(): void {
        this.onPreUpdate();
        this.onUpdate();
        this.onPostUpdate();
    }

    // The derived classes for 2D and 3D implement this to recompute the
    // boundary values when Neumann conditions are used.  If derived
    // classes built on top of the 2D or 3D classes implement this also,
    // they must call the base-class onPreUpdate first.
    protected abstract onPreUpdate(): void;

    // The derived classes for 2D and 3D implement this to iterate over
    // the image elements, updating an element only if it is not masked
    // out.
    protected abstract onUpdate(): void;

    // The derived classes for 2D and 3D implement this to swap the
    // buffers for the next pass.  If derived classes built on top of the
    // 2D or 3D classes implement this also, they must call the base-class
    // onPostUpdate last.
    protected abstract onPostUpdate(): void;
}
