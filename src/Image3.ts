// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) Image3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// The Image3 class is a 3-dimensional image built on the n-dimensional Image
// base. Pixels are stored in a flat array in lexicographical order, so the
// voxel (x,y,z) is at index x + dim0 * (y + dim1 * z).
//
// Port notes:
//   - Upstream overloads Reconstruct/GetIndex/operator()/Get on the argument
//     list; these become TypeScript overload sets on reconstruct/getIndex/
//     get/set (the base-class unchecked operator[]) and getClamped/setClamped
//     (the base-class clamping Get()). C++ copy/move machinery is not ported.
//   - Upstream guards the accessors with GTE_THROW_ON_IMAGE3_ERRORS, which is
//     disabled by default; the port matches the default (unchecked) build.
//   - The many GetNeighborhood/GetCorners/GetFull overloads are distinguished
//     in C++ only by the size and element type of the std::array
//     out-parameter. TypeScript has no such distinction, so (as in the Image2
//     port) the neighborhood size and the result kind are part of the method
//     name and the voxel is an optional argument:
//       GetNeighborhood(array<int32_t, 6>&)             -> getNeighborhood6()
//       GetNeighborhood(array<int32_t, 18>&)            -> getNeighborhood18()
//       GetNeighborhood(array<int32_t, 26>&)            -> getNeighborhood26()
//       GetCorners(array<int32_t, 8>&)                  -> getCorners8()
//       GetFull(array<int32_t, 27>&)                    -> getFull27()
//       GetNeighborhood(x, y, z, array<size_t, 6>&)     -> getNeighborhood6(x, y, z)
//       ... and likewise for 18, 26, GetCorners and GetFull
//       GetNeighborhood(array<array<int32_t, 3>, 6>&)   -> getNeighborhood6Coords()
//       ... and likewise for 18, 26, GetCorners and GetFull
//       GetNeighborhood(x, y, z, array<array<size_t, 3>, 6>&)
//                                                  -> getNeighborhood6Coords(x, y, z)
//       ... and likewise for 18, 26, GetCorners and GetFull
//     Out-parameters become return values.
//   - Upstream forms the absolute 3-tuple neighborhoods as
//     'static_cast<size_t>(x) + inbr[i][0]', which wraps around to a huge
//     unsigned value for a boundary voxel even though the documentation says
//     no clamping is applied. The port returns plain (possibly negative)
//     numbers, which is what the documented semantics call for.

import { Image } from './Image';

export class Image3<PixelType> extends Image<PixelType> {
    // Construction. With no arguments the image is empty. Otherwise the
    // dimensions must be positive; a nonpositive dimension leaves the image
    // empty, matching the Image base. See Image for 'createPixel'.
    constructor(dimension0?: number, dimension1?: number, dimension2?: number,
        createPixel?: () => PixelType) {
        super(dimension0 !== undefined && dimension1 !== undefined && dimension2 !== undefined
            ? [dimension0, dimension1, dimension2] : undefined, createPixel);
    }

    // Support for changing the image dimensions. All pixel data is lost by
    // this operation.
    override reconstruct(dimensions: readonly number[], createPixel?: () => PixelType): void;
    override reconstruct(dimension0: number, dimension1: number, dimension2: number,
        createPixel?: () => PixelType): void;
    override reconstruct(a: number | readonly number[], b?: number | (() => PixelType),
        c?: number, d?: () => PixelType): void {
        if (typeof a === 'number') {
            super.reconstruct([a, b as number, c as number], d);
        } else {
            super.reconstruct(a, b as (() => PixelType) | undefined);
        }
    }

    // Conversion between 1-dimensional indices and 3-dimensional coordinates.
    override getIndex(coord: readonly number[]): number;
    override getIndex(x: number, y: number, z: number): number;
    override getIndex(a: number | readonly number[], y?: number, z?: number): number {
        const dim0 = this.mDimensions[0];
        const dim1 = this.mDimensions[1];
        if (typeof a === 'number') {
            return a + dim0 * ((y as number) + dim1 * (z as number));
        }
        return a[0] + dim0 * (a[1] + dim1 * a[2]);
    }

    // The upstream out-parameters (x,y,z) are the returned 3-tuple.
    override getCoordinates(index: number): number[] {
        const x = index % this.mDimensions[0];
        const rest = Math.floor(index / this.mDimensions[0]);
        const y = rest % this.mDimensions[1];
        const z = Math.floor(rest / this.mDimensions[1]);
        return [x, y, z];
    }

    // Access the data as a 3-dimensional array. get/set do not test for valid
    // (x,y,z) (upstream operator()).
    override get(i: number): PixelType;
    override get(coord: readonly number[]): PixelType;
    override get(x: number, y: number, z: number): PixelType;
    override get(a: number | readonly number[], y?: number, z?: number): PixelType {
        if (typeof a === 'number') {
            if (y === undefined) {
                return this.mPixels[a];
            }
            return this.mPixels[a + this.mDimensions[0] * (y + this.mDimensions[1] * (z as number))];
        }
        return this.mPixels[a[0] + this.mDimensions[0] * (a[1] + this.mDimensions[1] * a[2])];
    }

    override set(i: number, value: PixelType): void;
    override set(coord: readonly number[], value: PixelType): void;
    override set(x: number, y: number, z: number, value: PixelType): void;
    override set(a: number | readonly number[], b: unknown, c?: number, d?: PixelType): void {
        if (typeof a === 'number') {
            if (c === undefined) {
                this.mPixels[a] = b as PixelType;
            } else {
                const i = a + this.mDimensions[0] * ((b as number) + this.mDimensions[1] * c);
                this.mPixels[i] = d as PixelType;
            }
        } else {
            this.mPixels[a[0] + this.mDimensions[0] * (a[1] + this.mDimensions[1] * a[2])] =
                b as PixelType;
        }
    }

    // The clamped accessors test for valid (x,y,z) and clamp when invalid;
    // these functions cannot fail (upstream Get()). The single-argument form
    // is the base-class linear-index version, which clamps to element 0.
    override getClamped(i: number): PixelType;
    override getClamped(coord: readonly number[]): PixelType;
    override getClamped(x: number, y: number, z: number): PixelType;
    override getClamped(a: number | readonly number[], y?: number, z?: number): PixelType {
        if (typeof a === 'number' && y === undefined) {
            return super.getClamped(a);
        }
        const coord = (typeof a === 'number' ? [a, y as number, z as number] : a);
        return this.mPixels[this.clampedIndex(coord)];
    }

    override setClamped(i: number, value: PixelType): void;
    override setClamped(coord: readonly number[], value: PixelType): void;
    override setClamped(x: number, y: number, z: number, value: PixelType): void;
    override setClamped(a: number | readonly number[], b: unknown, c?: number, d?: PixelType): void {
        if (typeof a === 'number' && c === undefined) {
            super.setClamped(a, b as PixelType);
            return;
        }
        const coord = (typeof a === 'number' ? [a, b as number, c as number] : a);
        const value = (typeof a === 'number' ? (d as PixelType) : (b as PixelType));
        this.mPixels[this.clampedIndex(coord)] = value;
    }

    // Clamp a 3-tuple to the valid coordinate range and convert to a
    // 1-dimensional index. The input is not modified (upstream clamps a
    // by-value copy).
    private clampedIndex(coord: readonly number[]): number {
        const clamped = [coord[0], coord[1], coord[2]];
        for (let d = 0; d < 3; ++d) {
            if (clamped[d] < 0) {
                clamped[d] = 0;
            } else if (clamped[d] >= this.mDimensions[d]) {
                clamped[d] = this.mDimensions[d] - 1;
            }
        }
        return clamped[0] + this.mDimensions[0] * (clamped[1] + this.mDimensions[1] * clamped[2]);
    }

    // In the following discussion, u, v and w are in {-1,1}. Given a voxel
    // (x,y,z), the 6-connected neighbors have relative offsets (u,0,0),
    // (0,v,0) and (0,0,w). The 18-connected neighbors include the 6-connected
    // neighbors and have additional relative offsets (u,v,0), (u,0,w) and
    // (0,v,w). The 26-connected neighbors include the 18-connected neighbors
    // and have additional relative offsets (u,v,w). The corner neighbors have
    // offsets (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), (1,0,1), (0,1,1)
    // and (1,1,1) in that order. The full neighborhood is the set of 3x3x3
    // voxels centered at (x,y,z).

    // The neighborhoods can be accessed as 1-dimensional indices using these
    // functions. The 0-argument overloads provide 1-dimensional offsets
    // relative to any voxel location; these depend only on the image
    // dimensions. The (x, y, z) overloads provide 1-dimensional indices for
    // the actual voxels in the neighborhood; no clamping is used when
    // (x,y,z) is on the boundary.
    getNeighborhood6(): number[];
    getNeighborhood6(x: number, y: number, z: number): number[];
    getNeighborhood6(x?: number, y?: number, z?: number): number[] {
        const dim0 = this.mDimensions[0];
        const dim01 = this.mDimensions[0] * this.mDimensions[1];
        const nbr = [
            -1,     // (x-1,y,z)
            +1,     // (x+1,y,z)
            -dim0,  // (x,y-1,z)
            +dim0,  // (x,y+1,z)
            -dim01, // (x,y,z-1)
            +dim01  // (x,y,z+1)
        ];
        return this.offsetBy(nbr, x, y, z);
    }

    getNeighborhood18(): number[];
    getNeighborhood18(x: number, y: number, z: number): number[];
    getNeighborhood18(x?: number, y?: number, z?: number): number[] {
        const dim0 = this.mDimensions[0];
        const dim01 = this.mDimensions[0] * this.mDimensions[1];
        const nbr = [
            ...this.getNeighborhood6(),
            -1 - dim0,      // (x-1,y-1,z)
            +1 - dim0,      // (x+1,y-1,z)
            -1 + dim0,      // (x-1,y+1,z)
            +1 + dim0,      // (x+1,y+1,z)
            -1 + dim01,     // (x-1,y,z+1)
            +1 + dim01,     // (x+1,y,z+1)
            -dim0 + dim01,  // (x,y-1,z+1)
            +dim0 + dim01,  // (x,y+1,z+1)
            -1 - dim01,     // (x-1,y,z-1)
            +1 - dim01,     // (x+1,y,z-1)
            -dim0 - dim01,  // (x,y-1,z-1)
            +dim0 - dim01   // (x,y+1,z-1)
        ];
        return this.offsetBy(nbr, x, y, z);
    }

    getNeighborhood26(): number[];
    getNeighborhood26(x: number, y: number, z: number): number[];
    getNeighborhood26(x?: number, y?: number, z?: number): number[] {
        const dim0 = this.mDimensions[0];
        const dim01 = this.mDimensions[0] * this.mDimensions[1];
        const nbr = [
            ...this.getNeighborhood18(),
            -1 - dim0 - dim01,  // (x-1,y-1,z-1)
            +1 - dim0 - dim01,  // (x+1,y-1,z-1)
            -1 + dim0 - dim01,  // (x-1,y+1,z-1)
            +1 + dim0 - dim01,  // (x+1,y+1,z-1)
            -1 - dim0 + dim01,  // (x-1,y-1,z+1)
            +1 - dim0 + dim01,  // (x+1,y-1,z+1)
            -1 + dim0 + dim01,  // (x-1,y+1,z+1)
            +1 + dim0 + dim01   // (x+1,y+1,z+1)
        ];
        return this.offsetBy(nbr, x, y, z);
    }

    getCorners8(): number[];
    getCorners8(x: number, y: number, z: number): number[];
    getCorners8(x?: number, y?: number, z?: number): number[] {
        const dim0 = this.mDimensions[0];
        const dim01 = this.mDimensions[0] * this.mDimensions[1];
        const nbr = [
            0,                  // (x,y,z)
            1,                  // (x+1,y,z)
            dim0,               // (x,y+1,z)
            dim0 + 1,           // (x+1,y+1,z)
            dim01,              // (x,y,z+1)
            dim01 + 1,          // (x+1,y,z+1)
            dim01 + dim0,       // (x,y+1,z+1)
            dim01 + dim0 + 1    // (x+1,y+1,z+1)
        ];
        return this.offsetBy(nbr, x, y, z);
    }

    getFull27(): number[];
    getFull27(x: number, y: number, z: number): number[];
    getFull27(x?: number, y?: number, z?: number): number[] {
        const dim0 = this.mDimensions[0];
        const dim01 = this.mDimensions[0] * this.mDimensions[1];
        const nbr = [
            -1 - dim0 - dim01,  // (x-1,y-1,z-1)
            -dim0 - dim01,      // (x,  y-1,z-1)
            +1 - dim0 - dim01,  // (x+1,y-1,z-1)
            -1 - dim01,         // (x-1,y,  z-1)
            -dim01,             // (x,  y,  z-1)
            +1 - dim01,         // (x+1,y,  z-1)
            -1 + dim0 - dim01,  // (x-1,y+1,z-1)
            +dim0 - dim01,      // (x,  y+1,z-1)
            +1 + dim0 - dim01,  // (x+1,y+1,z-1)
            -1 - dim0,          // (x-1,y-1,z)
            -dim0,              // (x,  y-1,z)
            +1 - dim0,          // (x+1,y-1,z)
            -1,                 // (x-1,y,  z)
            0,                  // (x,  y,  z)
            +1,                 // (x+1,y,  z)
            -1 + dim0,          // (x-1,y+1,z)
            +dim0,              // (x,  y+1,z)
            +1 + dim0,          // (x+1,y+1,z)
            -1 - dim0 + dim01,  // (x-1,y-1,z+1)
            -dim0 + dim01,      // (x,  y-1,z+1)
            +1 - dim0 + dim01,  // (x+1,y-1,z+1)
            -1 + dim01,         // (x-1,y,  z+1)
            +dim01,             // (x,  y,  z+1)
            +1 + dim01,         // (x+1,y,  z+1)
            -1 + dim0 + dim01,  // (x-1,y+1,z+1)
            +dim0 + dim01,      // (x,  y+1,z+1)
            +1 + dim0 + dim01   // (x+1,y+1,z+1)
        ];
        return this.offsetBy(nbr, x, y, z);
    }

    // Translate relative 1-dimensional offsets to absolute indices at
    // (x,y,z); when no voxel is specified the offsets are returned as is.
    private offsetBy(nbr: number[], x?: number, y?: number, z?: number): number[] {
        if (x !== undefined) {
            const index = this.getIndex(x, y as number, z as number);
            for (let i = 0; i < nbr.length; ++i) {
                nbr[i] += index;
            }
        }
        return nbr;
    }

    // The neighborhoods can be accessed as 3-tuples using these functions.
    // The 0-argument overloads provide 3-tuples relative to any voxel
    // location; these do not depend on the image dimensions. The (x, y, z)
    // overloads provide 3-tuples for the actual voxels in the neighborhood;
    // no clamping is used when (x,y,z) is on the boundary, so the returned
    // coordinates may be negative or exceed the image bounds.
    getNeighborhood6Coords(): [number, number, number][];
    getNeighborhood6Coords(x: number, y: number, z: number): [number, number, number][];
    getNeighborhood6Coords(x?: number, y?: number, z?: number): [number, number, number][] {
        const nbr: [number, number, number][] = [
            [-1, 0, 0],
            [+1, 0, 0],
            [0, -1, 0],
            [0, +1, 0],
            [0, 0, -1],
            [0, 0, +1]
        ];
        return Image3.translateCoords(nbr, x, y, z);
    }

    getNeighborhood18Coords(): [number, number, number][];
    getNeighborhood18Coords(x: number, y: number, z: number): [number, number, number][];
    getNeighborhood18Coords(x?: number, y?: number, z?: number): [number, number, number][] {
        const nbr: [number, number, number][] = [
            ...this.getNeighborhood6Coords(),
            [-1, -1, 0],
            [+1, -1, 0],
            [-1, +1, 0],
            [+1, +1, 0],
            [-1, 0, +1],
            [+1, 0, +1],
            [0, -1, +1],
            [0, +1, +1],
            [-1, 0, -1],
            [+1, 0, -1],
            [0, -1, -1],
            [0, +1, -1]
        ];
        return Image3.translateCoords(nbr, x, y, z);
    }

    getNeighborhood26Coords(): [number, number, number][];
    getNeighborhood26Coords(x: number, y: number, z: number): [number, number, number][];
    getNeighborhood26Coords(x?: number, y?: number, z?: number): [number, number, number][] {
        const nbr: [number, number, number][] = [
            ...this.getNeighborhood18Coords(),
            [-1, -1, -1],
            [+1, -1, -1],
            [-1, +1, -1],
            [+1, +1, -1],
            [-1, -1, +1],
            [+1, -1, +1],
            [-1, +1, +1],
            [+1, +1, +1]
        ];
        return Image3.translateCoords(nbr, x, y, z);
    }

    getCorners8Coords(): [number, number, number][];
    getCorners8Coords(x: number, y: number, z: number): [number, number, number][];
    getCorners8Coords(x?: number, y?: number, z?: number): [number, number, number][] {
        const nbr: [number, number, number][] = [
            [0, 0, 0],
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 0],
            [0, 0, 1],
            [1, 0, 1],
            [0, 1, 1],
            [1, 1, 1]
        ];
        return Image3.translateCoords(nbr, x, y, z);
    }

    getFull27Coords(): [number, number, number][];
    getFull27Coords(x: number, y: number, z: number): [number, number, number][];
    getFull27Coords(x?: number, y?: number, z?: number): [number, number, number][] {
        const nbr: [number, number, number][] = [
            [-1, -1, -1],
            [0, -1, -1],
            [+1, -1, -1],
            [-1, 0, -1],
            [0, 0, -1],
            [+1, 0, -1],
            [-1, +1, -1],
            [0, +1, -1],
            [+1, +1, -1],
            [-1, -1, 0],
            [0, -1, 0],
            [+1, -1, 0],
            [-1, 0, 0],
            [0, 0, 0],
            [+1, 0, 0],
            [-1, +1, 0],
            [0, +1, 0],
            [+1, +1, 0],
            [-1, -1, +1],
            [0, -1, +1],
            [+1, -1, +1],
            [-1, 0, +1],
            [0, 0, +1],
            [+1, 0, +1],
            [-1, +1, +1],
            [0, +1, +1],
            [+1, +1, +1]
        ];
        return Image3.translateCoords(nbr, x, y, z);
    }

    // Translate relative 3-tuples to the absolute coordinates at (x,y,z);
    // when no voxel is specified the relative 3-tuples are returned as is.
    private static translateCoords(coords: [number, number, number][],
        x?: number, y?: number, z?: number): [number, number, number][] {
        if (x !== undefined) {
            for (const coord of coords) {
                coord[0] += x;
                coord[1] += y as number;
                coord[2] += z as number;
            }
        }
        return coords;
    }
}
