// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Image.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Image class is a container for an n-dimensional image with pixels of a
// generic type, stored as a flat array in lexicographical order (dimension 0
// varies fastest).
//
// Port notes: upstream 'std::vector<PixelType> mPixels' value-initializes the
// pixels on resize (zero for numeric types, default construction otherwise).
// TypeScript cannot default-construct a generic type, so construction and
// reconstruct() accept an optional 'createPixel' factory; when it is omitted
// the pixels are zero-filled, which matches upstream for the numeric pixel
// types that dominate usage. For object pixel types pass a factory so pixels
// do not alias one another. C++ copy/move machinery is not ported. The
// unchecked operator[] becomes get(i)/set(i, value); the clamping Get()
// becomes getClamped(i)/setClamped(i, value). The out-parameter of
// GetCoordinates becomes the return value.

export class Image<PixelType> {
    protected mDimensions: number[];
    protected mOffsets: number[];
    protected mPixels: PixelType[];

    // Construction. With no arguments the image is empty (zero dimensions,
    // no pixels). Otherwise the image is reconstructed from the dimensions.
    constructor(dimensions?: readonly number[], createPixel?: () => PixelType) {
        this.mDimensions = [];
        this.mOffsets = [];
        this.mPixels = [];
        if (dimensions !== undefined) {
            this.reconstruct(dimensions, createPixel);
        }
    }

    // Support for changing the image dimensions. All pixel data is lost by
    // this operation. As upstream, a nonpositive dimension leaves the image
    // empty.
    reconstruct(dimensions: readonly number[], createPixel?: () => PixelType): void {
        this.mDimensions = [];
        this.mOffsets = [];
        this.mPixels = [];

        if (dimensions.length > 0) {
            for (const dim of dimensions) {
                if (dim <= 0) {
                    return;
                }
            }

            this.mDimensions = dimensions.slice();
            this.mOffsets = new Array<number>(dimensions.length);

            let numPixels = 1;
            for (let d = 0; d < dimensions.length; ++d) {
                numPixels *= this.mDimensions[d];
            }

            this.mOffsets[0] = 1;
            for (let d = 1; d < dimensions.length; ++d) {
                this.mOffsets[d] = this.mDimensions[d - 1] * this.mOffsets[d - 1];
            }

            this.mPixels = new Array<PixelType>(numPixels);
            if (createPixel !== undefined) {
                for (let i = 0; i < numPixels; ++i) {
                    this.mPixels[i] = createPixel();
                }
            } else {
                // Matches upstream value-initialization for numeric types.
                this.mPixels.fill(0 as unknown as PixelType);
            }
        }
    }

    // Access to image data.
    getDimensions(): readonly number[] {
        return this.mDimensions;
    }

    getNumDimensions(): number {
        return this.mDimensions.length;
    }

    getDimension(d: number): number {
        return this.mDimensions[d];
    }

    getOffsets(): readonly number[] {
        return this.mOffsets;
    }

    getOffset(d: number): number {
        return this.mOffsets[d];
    }

    getPixels(): PixelType[] {
        return this.mPixels;
    }

    getNumPixels(): number {
        return this.mPixels.length;
    }

    // Conversions between n-dim and 1-dim structures. The 'coord' array must
    // have getNumDimensions() elements.
    getIndex(coord: readonly number[]): number {
        const numDimensions = this.mDimensions.length;
        let index = coord[0];
        for (let d = 1; d < numDimensions; ++d) {
            index += this.mOffsets[d] * coord[d];
        }
        return index;
    }

    // The upstream out-parameter 'coord' is the return value.
    getCoordinates(index: number): number[] {
        const numDimensions = this.mDimensions.length;
        const coord = new Array<number>(numDimensions);
        for (let d = 0; d < numDimensions; ++d) {
            coord[d] = index % this.mDimensions[d];
            index = Math.floor(index / this.mDimensions[d]);
        }
        return coord;
    }

    // Access the data as a 1-dimensional array. get/set do not test for
    // valid i (upstream operator[]). The clamped accessors test for valid i
    // and clamp to element 0 when invalid; these functions cannot fail
    // (upstream Get()).
    get(i: number): PixelType {
        return this.mPixels[i];
    }

    set(i: number, value: PixelType): void {
        this.mPixels[i] = value;
    }

    getClamped(i: number): PixelType {
        return (0 <= i && i < this.mPixels.length ? this.mPixels[i] : this.mPixels[0]);
    }

    setClamped(i: number, value: PixelType): void {
        if (0 <= i && i < this.mPixels.length) {
            this.mPixels[i] = value;
        } else {
            this.mPixels[0] = value;
        }
    }
}
