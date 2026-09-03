// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Image2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Image2 is a 2-dimensional image with pixels stored in row-major order for
// (x,y): the pixel at (x,y) is at 1-dimensional index x + dimension0 * y.
//
// Port notes: upstream guards coordinate validation behind the compile-time
// define GTE_THROW_ON_IMAGE2_ERRORS; the port implements the unguarded
// (non-throwing) paths. As in the Image base port, the unchecked
// operator()(x, y) becomes get(x, y)/set(x, y, value) overloads and the
// clamping Get(x, y) becomes getClamped(x, y)/setClamped(x, y, value)
// overloads; the std::array<int32_t, 2> coordinate variants accept a
// readonly number[] of length 2. The out-parameter of GetCoordinates
// becomes the return value. The many upstream GetNeighborhood/GetCorners/
// GetFull overloads are distinguished in C++ by the out-parameter array
// type; the port distinguishes them by name and arity:
//   GetNeighborhood(array<int32_t, 4>&)            -> getNeighborhood4()
//   GetNeighborhood(array<int32_t, 8>&)            -> getNeighborhood8()
//   GetCorners(array<int32_t, 4>&)                 -> getCorners()
//   GetFull(array<int32_t, 9>&)                    -> getFull()
//   GetNeighborhood(x, y, array<size_t, 4>&)       -> getNeighborhood4(x, y)
//   GetNeighborhood(x, y, array<size_t, 8>&)       -> getNeighborhood8(x, y)
//   GetCorners(x, y, array<size_t, 4>&)            -> getCorners(x, y)
//   GetFull(x, y, array<size_t, 9>&)               -> getFull(x, y)
//   GetNeighborhood(array<array<int32_t, 2>, 4>&)  -> getNeighborhood4Coords()
//   GetNeighborhood(array<array<int32_t, 2>, 8>&)  -> getNeighborhood8Coords()
//   GetCorners(array<array<int32_t, 2>, 4>&)       -> getCornersCoords()
//   GetFull(array<array<int32_t, 2>, 9>&)          -> getFullCoords()
//   GetNeighborhood(x, y, array<array<size_t, 2>, 4>&) -> getNeighborhood4Coords(x, y)
//   GetNeighborhood(x, y, array<array<size_t, 2>, 8>&) -> getNeighborhood8Coords(x, y)
//   GetCorners(x, y, array<array<size_t, 2>, 4>&)  -> getCornersCoords(x, y)
//   GetFull(x, y, array<array<size_t, 2>, 9>&)     -> getFullCoords(x, y)
// Note that the relative 1-dimensional 8-neighborhood ordering (4-connected
// neighbors first, then diagonal neighbors) differs from the relative
// 2-tuple 8-neighborhood ordering (row-major order without the center);
// both match upstream exactly.

import { Image } from './Image.js';

export class Image2<PixelType> extends Image<PixelType> {
    // Construction. With no arguments the image is empty. Otherwise both
    // dimensions must be positive; a nonpositive dimension leaves the image
    // empty (base class behavior). See the Image port notes for the role of
    // the optional 'createPixel' factory.
    constructor(dimension0?: number, dimension1?: number, createPixel?: () => PixelType) {
        if (dimension0 !== undefined && dimension1 !== undefined) {
            super([dimension0, dimension1], createPixel);
        } else {
            super();
        }
    }

    // Support for changing the image dimensions. All pixel data is lost by
    // this operation. The (dimension0, dimension1) overload is the upstream
    // Image2::Reconstruct; the array overload is the inherited base version.
    override reconstruct(dimensions: readonly number[], createPixel?: () => PixelType): void;
    override reconstruct(dimension0: number, dimension1: number, createPixel?: () => PixelType): void;
    override reconstruct(
        arg0: readonly number[] | number,
        arg1?: number | (() => PixelType),
        arg2?: () => PixelType): void {
        if (typeof arg0 === 'number') {
            super.reconstruct([arg0, arg1 as number], arg2);
        } else {
            super.reconstruct(arg0, arg1 as (() => PixelType) | undefined);
        }
    }

    // Conversion between 1-dimensional indices and 2-dimensional
    // coordinates.
    override getIndex(coord: readonly number[]): number;
    override getIndex(x: number, y: number): number;
    override getIndex(arg0: readonly number[] | number, y?: number): number {
        if (typeof arg0 === 'number') {
            return arg0 + this.mDimensions[0] * (y as number);
        }
        return arg0[0] + this.mDimensions[0] * arg0[1];
    }

    override getCoordinates(index: number): [number, number] {
        return [index % this.mDimensions[0], Math.floor(index / this.mDimensions[0])];
    }

    // Access the data as a 2-dimensional array. The (x, y) and coordinate
    // overloads do not test for valid (x, y) (upstream operator()).
    override get(i: number): PixelType;
    override get(coord: readonly number[]): PixelType;
    override get(x: number, y: number): PixelType;
    override get(arg0: number | readonly number[], y?: number): PixelType {
        if (typeof arg0 === 'number') {
            if (y === undefined) {
                return this.mPixels[arg0];
            }
            return this.mPixels[arg0 + this.mDimensions[0] * y];
        }
        return this.mPixels[arg0[0] + this.mDimensions[0] * arg0[1]];
    }

    override set(i: number, value: PixelType): void;
    override set(coord: readonly number[], value: PixelType): void;
    override set(x: number, y: number, value: PixelType): void;
    override set(arg0: number | readonly number[], arg1: number | PixelType, arg2?: PixelType): void {
        if (typeof arg0 === 'number') {
            if (arguments.length === 3) {
                this.mPixels[arg0 + this.mDimensions[0] * (arg1 as number)] = arg2 as PixelType;
            } else {
                this.mPixels[arg0] = arg1 as PixelType;
            }
        } else {
            this.mPixels[arg0[0] + this.mDimensions[0] * arg0[1]] = arg1 as PixelType;
        }
    }

    // The clamped accessors test for valid (x, y) and clamp each coordinate
    // to the image boundary when invalid; these functions cannot fail
    // (upstream Get(x, y)). The 1-argument flat-index overload is the base
    // class behavior (clamps to element 0).
    override getClamped(i: number): PixelType;
    override getClamped(coord: readonly number[]): PixelType;
    override getClamped(x: number, y: number): PixelType;
    override getClamped(arg0: number | readonly number[], y?: number): PixelType {
        if (typeof arg0 === 'number') {
            if (y === undefined) {
                return super.getClamped(arg0);
            }
            return this.mPixels[this.getClampedIndex(arg0, y)];
        }
        return this.mPixels[this.getClampedIndex(arg0[0], arg0[1])];
    }

    override setClamped(i: number, value: PixelType): void;
    override setClamped(coord: readonly number[], value: PixelType): void;
    override setClamped(x: number, y: number, value: PixelType): void;
    override setClamped(arg0: number | readonly number[], arg1: number | PixelType, arg2?: PixelType): void {
        if (typeof arg0 === 'number') {
            if (arguments.length === 3) {
                this.mPixels[this.getClampedIndex(arg0, arg1 as number)] = arg2 as PixelType;
            } else {
                super.setClamped(arg0, arg1 as PixelType);
            }
        } else {
            this.mPixels[this.getClampedIndex(arg0[0], arg0[1])] = arg1 as PixelType;
        }
    }

    private getClampedIndex(x: number, y: number): number {
        // Clamp to valid (x, y).
        if (x < 0) {
            x = 0;
        } else if (x >= this.mDimensions[0]) {
            x = this.mDimensions[0] - 1;
        }

        if (y < 0) {
            y = 0;
        } else if (y >= this.mDimensions[1]) {
            y = this.mDimensions[1] - 1;
        }

        return x + this.mDimensions[0] * y;
    }

    // In the following discussion, u and v are in {-1,1}. Given a pixel
    // (x,y), the 4-connected neighbors have relative offsets (u,0) and
    // (0,v). The 8-connected neighbors include the 4-connected neighbors
    // and have additional relative offsets (u,v). The corner neighbors
    // have relative offsets (0,0), (1,0), (0,1), and (1,1) in that order.
    // The full neighborhood is the set of 3x3 pixels centered at (x,y).

    // The neighborhoods can be accessed as 1-dimensional indices using
    // these functions. The 0-argument overloads provide 1-dimensional
    // indices relative to any pixel location; these depend only on the
    // image dimensions. The (x, y) overloads provide 1-dimensional indices
    // for the actual pixels in the neighborhood; no clamping is used when
    // (x,y) is on the boundary.
    getNeighborhood4(): number[];
    getNeighborhood4(x: number, y: number): number[];
    getNeighborhood4(x?: number, y?: number): number[] {
        const dim0 = this.mDimensions[0];
        const nbr = [
            -1,     // (x-1,y)
            +1,     // (x+1,y)
            -dim0,  // (x,y-1)
            +dim0   // (x,y+1)
        ];
        if (x !== undefined) {
            const index = this.getIndex(x, y as number);
            for (let i = 0; i < 4; ++i) {
                nbr[i] += index;
            }
        }
        return nbr;
    }

    getNeighborhood8(): number[];
    getNeighborhood8(x: number, y: number): number[];
    getNeighborhood8(x?: number, y?: number): number[] {
        const dim0 = this.mDimensions[0];
        const nbr = [
            -1,         // (x-1,y)
            +1,         // (x+1,y)
            -dim0,      // (x,y-1)
            +dim0,      // (x,y+1)
            -1 - dim0,  // (x-1,y-1)
            +1 - dim0,  // (x+1,y-1)
            -1 + dim0,  // (x-1,y+1)
            +1 + dim0   // (x+1,y+1)
        ];
        if (x !== undefined) {
            const index = this.getIndex(x, y as number);
            for (let i = 0; i < 8; ++i) {
                nbr[i] += index;
            }
        }
        return nbr;
    }

    getCorners(): number[];
    getCorners(x: number, y: number): number[];
    getCorners(x?: number, y?: number): number[] {
        const dim0 = this.mDimensions[0];
        const nbr = [
            0,          // (x,y)
            1,          // (x+1,y)
            dim0,       // (x,y+1)
            dim0 + 1    // (x+1,y+1)
        ];
        if (x !== undefined) {
            const index = this.getIndex(x, y as number);
            for (let i = 0; i < 4; ++i) {
                nbr[i] += index;
            }
        }
        return nbr;
    }

    getFull(): number[];
    getFull(x: number, y: number): number[];
    getFull(x?: number, y?: number): number[] {
        const dim0 = this.mDimensions[0];
        const nbr = [
            -1 - dim0,  // (x-1,y-1)
            -dim0,      // (x,y-1)
            +1 - dim0,  // (x+1,y-1)
            -1,         // (x-1,y)
            0,          // (x,y)
            +1,         // (x+1,y)
            -1 + dim0,  // (x-1,y+1)
            +dim0,      // (x,y+1)
            +1 + dim0   // (x+1,y+1)
        ];
        if (x !== undefined) {
            const index = this.getIndex(x, y as number);
            for (let i = 0; i < 9; ++i) {
                nbr[i] += index;
            }
        }
        return nbr;
    }

    // The neighborhoods can be accessed as 2-tuples using these functions.
    // The 0-argument overloads provide 2-tuples relative to any pixel
    // location; these depend only on the image dimensions. The (x, y)
    // overloads provide 2-tuples for the actual pixels in the neighborhood;
    // no clamping is used when (x,y) is on the boundary.
    getNeighborhood4Coords(): [number, number][];
    getNeighborhood4Coords(x: number, y: number): [number, number][];
    getNeighborhood4Coords(x?: number, y?: number): [number, number][] {
        const nbr: [number, number][] = [
            [-1, 0],
            [+1, 0],
            [0, -1],
            [0, +1]
        ];
        if (x !== undefined) {
            for (let i = 0; i < 4; ++i) {
                nbr[i][0] += x;
                nbr[i][1] += y as number;
            }
        }
        return nbr;
    }

    getNeighborhood8Coords(): [number, number][];
    getNeighborhood8Coords(x: number, y: number): [number, number][];
    getNeighborhood8Coords(x?: number, y?: number): [number, number][] {
        const nbr: [number, number][] = [
            [-1, -1],
            [0, -1],
            [+1, -1],
            [-1, 0],
            [+1, 0],
            [-1, +1],
            [0, +1],
            [+1, +1]
        ];
        if (x !== undefined) {
            for (let i = 0; i < 8; ++i) {
                nbr[i][0] += x;
                nbr[i][1] += y as number;
            }
        }
        return nbr;
    }

    getCornersCoords(): [number, number][];
    getCornersCoords(x: number, y: number): [number, number][];
    getCornersCoords(x?: number, y?: number): [number, number][] {
        const nbr: [number, number][] = [
            [0, 0],
            [1, 0],
            [0, 1],
            [1, 1]
        ];
        if (x !== undefined) {
            for (let i = 0; i < 4; ++i) {
                nbr[i][0] += x;
                nbr[i][1] += y as number;
            }
        }
        return nbr;
    }

    getFullCoords(): [number, number][];
    getFullCoords(x: number, y: number): [number, number][];
    getFullCoords(x?: number, y?: number): [number, number][] {
        const nbr: [number, number][] = [
            [-1, -1],
            [0, -1],
            [+1, -1],
            [-1, 0],
            [0, 0],
            [+1, 0],
            [-1, +1],
            [0, +1],
            [+1, +1]
        ];
        if (x !== undefined) {
            for (let i = 0; i < 9; ++i) {
                nbr[i][0] += x;
                nbr[i][1] += y as number;
            }
        }
        return nbr;
    }
}
