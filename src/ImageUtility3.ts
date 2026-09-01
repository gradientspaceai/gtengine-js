// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ImageUtility3.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Image utilities for Image3<number> objects whose voxels are integers.
//
// All but the draw* functions are operations on binary images. Let the image
// have d0 columns, d1 rows, and d2 slices. The input image must have zeros on
// its boundaries x = 0, x = d0-1, y = 0, y = d1-1, z = 0, and z = d2-1. The
// 0-valued voxels are considered to be background. The 1-valued voxels are
// considered to be foreground. In some of the operations, to save memory and
// time the input image is modified by the algorithms. If you need to preserve
// the input image, make a copy of it before calling these functions.
//
// Port notes:
//   - As in the ImageUtility2 port, upstream is a C++ class of static member
//     functions and the port keeps that shape, so the shared names
//     (getComponents, dilate, drawLine, ...) live on the classes rather than
//     the module scope and no 2D/3D suffixes are required for global export
//     uniqueness.
//   - The C++ overloads distinguished by the size N of a std::array template
//     argument (N in {6, 18, 26}) become part of the method name, and the
//     (count, pointer) variants take one array whose length is the count:
//       GetComponents<N>(image, components)   -> getComponents6/18/26(image)
//       GetComponents(image, n, nbrs, comps)  -> getComponents(image, nbrs)
//       Dilate<N>(in, out)                    -> dilate6/18/26(inImage)
//       Dilate(in, n, nbrs, out)              -> dilate(inImage, nbrs)
//       ... and likewise for Erode, Open and Close.
//     The 1-dimensional neighbor offsets come from
//     Image3.getNeighborhood6/18/26(); the 3-tuple offsets used by the
//     morphology operations come from Image3.getNeighborhood6/18/26Coords().
//   - Upstream out-parameters become return values: the morphology functions
//     return a newly allocated output image (so the upstream "input and
//     output must be different objects" assertion holds by construction) and
//     getComponents returns the component index lists. getComponents,
//     computeCDConvex and floodFill6 still modify the input image in place.
//
// Upstream bug suspects (both fixed here; see the PR description):
//   1. ImageUtility3::Dilate iterates 'for (int32_t i0 = 1; i0 < dim0; ++i0)',
//      skipping the x = 0 column as a dilation source. The 2D counterpart
//      starts at i0 = 0 and the documented algorithm places no such
//      restriction, so this port starts at i0 = 0.
//   2. ImageUtility3::Close<N> has the 2D static_assert (N == 4 || N == 8)
//      and constructs its temporary with the two-argument Image3 constructor,
//      which does not exist. The template therefore cannot be instantiated;
//      it is dead code. The port implements close6/close18/close26 with the
//      correct 3D neighborhoods.

import { Image3 } from './Image3';
import { logAssert } from './Logger';

// The type of the callback invoked by drawLine for each visited voxel
// (upstream std::function<void(int32_t, int32_t, int32_t)>).
export type ImageUtility3Callback = (x: number, y: number, z: number) => void;

export class ImageUtility3 {
    // Compute the N-connected components of a binary image (N is 6, 18 or
    // 26). The input image is modified to avoid the cost of making a copy.
    // On output, the image values are the labels for the components. The
    // returned array components[k], k >= 1, contains the indices for the
    // k-th component; components[0] is unused. When there are no components
    // the returned array is empty.
    static getComponents6(image: Image3<number>): number[][] {
        return ImageUtility3.getComponents(image, image.getNeighborhood6());
    }

    static getComponents18(image: Image3<number>): number[][] {
        return ImageUtility3.getComponents(image, image.getNeighborhood18());
    }

    static getComponents26(image: Image3<number>): number[][] {
        return ImageUtility3.getComponents(image, image.getNeighborhood26());
    }

    // Connected component labeling using depth-first search. The
    // neighborhood is specified by the caller as 1-dimensional index offsets
    // into the image.
    static getComponents(image: Image3<number>, neighbors: readonly number[]): number[][] {
        logAssert(neighbors.length > 0, 'Invalid neighbors.');

        const numNeighbors = neighbors.length;
        const numVoxels = image.getNumPixels();
        const voxels = image.getPixels();
        const numElements = new Array<number>(numVoxels).fill(0);
        const vstack = new Array<number>(numVoxels).fill(0);
        let numComponents = 0;
        let label = 2;
        for (let i = 0; i < numVoxels; ++i) {
            if (voxels[i] === 1) {
                let top = -1;
                vstack[++top] = i;

                // Upstream binds a reference 'count' to numElements[k]; the
                // port indexes the array directly.
                const k = numComponents + 1;
                numElements[k] = 0;
                while (top >= 0) {
                    const v = vstack[top];
                    voxels[v] = -1;
                    let j: number;
                    for (j = 0; j < numNeighbors; ++j) {
                        const adj = v + neighbors[j];
                        if (voxels[adj] === 1) {
                            vstack[++top] = adj;
                            break;
                        }
                    }
                    if (j === numNeighbors) {
                        voxels[v] = label;
                        ++numElements[k];
                        --top;
                    }
                }

                ++numComponents;
                ++label;
            }
        }

        const components: number[][] = [];
        if (numComponents > 0) {
            for (let i = 0; i <= numComponents; ++i) {
                components.push([]);
            }
            for (let i = 1; i <= numComponents; ++i) {
                components[i] = new Array<number>(numElements[i]).fill(0);
                numElements[i] = 0;
            }

            for (let i = 0; i < numVoxels; ++i) {
                let value = voxels[i];
                if (value !== 0) {
                    // Labels started at 2 to support the depth-first search,
                    // so they need to be decremented for the correct labels.
                    voxels[i] = --value;
                    components[value][numElements[value]] = i;
                    ++numElements[value];
                }
            }
        }
        return components;
    }

    // Compute a dilation with a structuring element consisting of the
    // N-connected neighbors of each voxel (N is 6, 18 or 26). The input image
    // is binary with 0 for background and 1 for foreground. The dilated image
    // is returned; the input image is not modified.
    static dilate6(inImage: Image3<number>): Image3<number> {
        return ImageUtility3.dilate(inImage, inImage.getNeighborhood6Coords());
    }

    static dilate18(inImage: Image3<number>): Image3<number> {
        return ImageUtility3.dilate(inImage, inImage.getNeighborhood18Coords());
    }

    static dilate26(inImage: Image3<number>): Image3<number> {
        return ImageUtility3.dilate(inImage, inImage.getNeighborhood26Coords());
    }

    // Compute a dilation with a structuring element consisting of neighbors
    // specified by 3-tuple offsets relative to the voxel.
    static dilate(inImage: Image3<number>,
        neighbors: readonly (readonly number[])[]): Image3<number> {
        logAssert(neighbors.length > 0, 'Invalid neighbors.');

        const outImage = copyImage3(inImage);

        // If the voxel at (i0,i1,i2) is 1, then the voxels at
        // (k0,k1,k2) = (i0+nbr0,i1+nbr1,i2+nbr2) are set to 1 where
        // (nbr0,nbr1,nbr2) is in the neighbors array. Boundary testing is
        // used to avoid accessing out-of-range voxels.
        const dim0 = inImage.getDimension(0);
        const dim1 = inImage.getDimension(1);
        const dim2 = inImage.getDimension(2);
        for (let i2 = 0; i2 < dim2; ++i2) {
            for (let i1 = 0; i1 < dim1; ++i1) {
                // Upstream starts this loop at i0 = 1, which skips the x = 0
                // column as a dilation source. See the file header.
                for (let i0 = 0; i0 < dim0; ++i0) {
                    if (inImage.get(i0, i1, i2) === 1) {
                        for (let n = 0; n < neighbors.length; ++n) {
                            const k0 = i0 + neighbors[n][0];
                            const k1 = i1 + neighbors[n][1];
                            const k2 = i2 + neighbors[n][2];
                            if (0 <= k0 && k0 < dim0 &&
                                0 <= k1 && k1 < dim1 &&
                                0 <= k2 && k2 < dim2) {
                                outImage.set(k0, k1, k2, 1);
                            }
                        }
                    }
                }
            }
        }
        return outImage;
    }

    // Compute an erosion with a structuring element consisting of the
    // N-connected neighbors of each voxel (N is 6, 18 or 26). The input image
    // is binary with 0 for background and 1 for foreground. If zeroExterior
    // is true, the image exterior is assumed to be 0, so 1-valued boundary
    // voxels are set to 0; otherwise, boundary voxels are set to 0 only when
    // they have neighboring image voxels that are 0. The eroded image is
    // returned; the input image is not modified.
    static erode6(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode(inImage, zeroExterior, inImage.getNeighborhood6Coords());
    }

    static erode18(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode(inImage, zeroExterior, inImage.getNeighborhood18Coords());
    }

    static erode26(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode(inImage, zeroExterior, inImage.getNeighborhood26Coords());
    }

    // Compute an erosion with a structuring element consisting of neighbors
    // specified by 3-tuple offsets relative to the voxel.
    static erode(inImage: Image3<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image3<number> {
        logAssert(neighbors.length > 0, 'Invalid neighbors.');

        const outImage = copyImage3(inImage);

        // If the voxel at (i0,i1,i2) is 1, it is changed to 0 when at least
        // one neighbor (k0,k1,k2) = (i0+nbr0,i1+nbr1,i2+nbr2) is 0, where
        // (nbr0,nbr1,nbr2) is in the neighbors array.
        const dim0 = inImage.getDimension(0);
        const dim1 = inImage.getDimension(1);
        const dim2 = inImage.getDimension(2);
        for (let i2 = 0; i2 < dim2; ++i2) {
            for (let i1 = 0; i1 < dim1; ++i1) {
                for (let i0 = 0; i0 < dim0; ++i0) {
                    if (inImage.get(i0, i1, i2) === 1) {
                        for (let j = 0; j < neighbors.length; ++j) {
                            const k0 = i0 + neighbors[j][0];
                            const k1 = i1 + neighbors[j][1];
                            const k2 = i2 + neighbors[j][2];
                            if (0 <= k0 && k0 < dim0 &&
                                0 <= k1 && k1 < dim1 &&
                                0 <= k2 && k2 < dim2) {
                                if (inImage.get(k0, k1, k2) === 0) {
                                    outImage.set(i0, i1, i2, 0);
                                    break;
                                }
                            } else if (zeroExterior) {
                                outImage.set(i0, i1, i2, 0);
                                break;
                            }
                        }
                    }
                }
            }
        }
        return outImage;
    }

    // Compute an opening (erosion followed by dilation) with a structuring
    // element consisting of the N-connected neighbors of each voxel (N is 6,
    // 18 or 26). See erode for the meaning of zeroExterior.
    static open6(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.dilate6(ImageUtility3.erode6(inImage, zeroExterior));
    }

    static open18(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.dilate18(ImageUtility3.erode18(inImage, zeroExterior));
    }

    static open26(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.dilate26(ImageUtility3.erode26(inImage, zeroExterior));
    }

    // Compute an opening with a caller-specified structuring element.
    static open(inImage: Image3<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image3<number> {
        const temp = ImageUtility3.erode(inImage, zeroExterior, neighbors);
        return ImageUtility3.dilate(temp, neighbors);
    }

    // Compute a closing (dilation followed by erosion) with a structuring
    // element consisting of the N-connected neighbors of each voxel (N is 6,
    // 18 or 26). See erode for the meaning of zeroExterior.
    static close6(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode6(ImageUtility3.dilate6(inImage), zeroExterior);
    }

    static close18(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode18(ImageUtility3.dilate18(inImage), zeroExterior);
    }

    static close26(inImage: Image3<number>, zeroExterior: boolean): Image3<number> {
        return ImageUtility3.erode26(ImageUtility3.dilate26(inImage), zeroExterior);
    }

    // Compute a closing with a caller-specified structuring element.
    static close(inImage: Image3<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image3<number> {
        const temp = ImageUtility3.dilate(inImage, neighbors);
        return ImageUtility3.erode(temp, zeroExterior, neighbors);
    }

    // Compute the coordinate-directional convex set, in place. For a given
    // coordinate direction (x, y or z), identify the first and last 1-valued
    // voxels on a segment of voxels in that direction. All voxels from first
    // to last are set to 1. This is done for all segments in each of the
    // coordinate directions.
    static computeCDConvex(image: Image3<number>): void {
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);
        const dim2 = image.getDimension(2);

        const temp = copyImage3(image);
        let i0 = 0;
        let i1 = 0;
        let i2 = 0;
        for (i1 = 0; i1 < dim1; ++i1) {
            for (i0 = 0; i0 < dim0; ++i0) {
                let i2min: number;
                for (i2min = 0; i2min < dim2; ++i2min) {
                    if ((temp.get(i0, i1, i2min) & 1) === 0) {
                        temp.set(i0, i1, i2min, temp.get(i0, i1, i2min) | 2);
                    } else {
                        break;
                    }
                }
                if (i2min < dim2) {
                    for (let i2max = dim2 - 1; i2max >= i2min; --i2max) {
                        if ((temp.get(i0, i1, i2max) & 1) === 0) {
                            temp.set(i0, i1, i2max, temp.get(i0, i1, i2max) | 2);
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        for (i2 = 0; i2 < dim2; ++i2) {
            for (i0 = 0; i0 < dim0; ++i0) {
                let i1min: number;
                for (i1min = 0; i1min < dim1; ++i1min) {
                    if ((temp.get(i0, i1min, i2) & 1) === 0) {
                        temp.set(i0, i1min, i2, temp.get(i0, i1min, i2) | 2);
                    } else {
                        break;
                    }
                }
                if (i1min < dim1) {
                    for (let i1max = dim1 - 1; i1max >= i1min; --i1max) {
                        if ((temp.get(i0, i1max, i2) & 1) === 0) {
                            temp.set(i0, i1max, i2, temp.get(i0, i1max, i2) | 2);
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        for (i2 = 0; i2 < dim2; ++i2) {
            for (i1 = 0; i1 < dim1; ++i1) {
                let i0min: number;
                for (i0min = 0; i0min < dim0; ++i0min) {
                    if ((temp.get(i0min, i1, i2) & 1) === 0) {
                        temp.set(i0min, i1, i2, temp.get(i0min, i1, i2) | 2);
                    } else {
                        break;
                    }
                }
                if (i0min < dim0) {
                    for (let i0max = dim0 - 1; i0max >= i0min; --i0max) {
                        if ((temp.get(i0max, i1, i2) & 1) === 0) {
                            temp.set(i0max, i1, i2, temp.get(i0max, i1, i2) | 2);
                        } else {
                            break;
                        }
                    }
                }
            }
        }

        const dst = image.getPixels();
        const src = temp.getPixels();
        for (let i = 0; i < dst.length; ++i) {
            dst[i] = ((src[i] & 2) !== 0 ? 0 : 1);
        }
    }

    // Use a depth-first search for filling a 6-connected region. This is
    // nonrecursive, simulated by using a heap-allocated "stack". The input
    // (x,y,z) is the seed point that starts the fill.
    static floodFill6<PixelType>(image: Image3<PixelType>, x: number, y: number, z: number,
        foreColor: PixelType, backColor: PixelType): void {
        // Test for a valid seed.
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);
        const dim2 = image.getDimension(2);
        if (x < 0 || x >= dim0 || y < 0 || y >= dim1 || z < 0 || z >= dim2) {
            // The seed point is outside the image domain, so there is nothing
            // to fill.
            return;
        }

        // Allocate the maximum amount of space needed for the stack. An empty
        // stack has top == -1.
        const numVoxels = image.getNumPixels();
        const xStack = new Array<number>(numVoxels).fill(0);
        const yStack = new Array<number>(numVoxels).fill(0);
        const zStack = new Array<number>(numVoxels).fill(0);

        // Push the seed point onto the stack. All points pushed onto the
        // stack have the background color backColor.
        let top = 0;
        xStack[top] = x;
        yStack[top] = y;
        zStack[top] = z;

        while (top >= 0) {
            // Read the top of the stack. Do not pop, because we need to
            // return to this top value later to restart the fill in a
            // different direction.
            x = xStack[top];
            y = yStack[top];
            z = zStack[top];

            // Fill the voxel.
            image.set(x, y, z, foreColor);

            const xp1 = x + 1;
            if (xp1 < dim0 && image.get(xp1, y, z) === backColor) {
                ++top;
                xStack[top] = xp1;
                yStack[top] = y;
                zStack[top] = z;
                continue;
            }

            const xm1 = x - 1;
            if (0 <= xm1 && image.get(xm1, y, z) === backColor) {
                ++top;
                xStack[top] = xm1;
                yStack[top] = y;
                zStack[top] = z;
                continue;
            }

            const yp1 = y + 1;
            if (yp1 < dim1 && image.get(x, yp1, z) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = yp1;
                zStack[top] = z;
                continue;
            }

            const ym1 = y - 1;
            if (0 <= ym1 && image.get(x, ym1, z) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = ym1;
                zStack[top] = z;
                continue;
            }

            const zp1 = z + 1;
            if (zp1 < dim2 && image.get(x, y, zp1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = y;
                zStack[top] = zp1;
                continue;
            }

            const zm1 = z - 1;
            if (0 <= zm1 && image.get(x, y, zm1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = y;
                zStack[top] = zm1;
                continue;
            }

            // Done in all directions; pop and return to search the other
            // directions of the predecessor.
            --top;
        }
    }

    // Visit voxels using Bresenham's line drawing algorithm. The callback
    // represents the action you want applied to each voxel as it is visited.
    static drawLine(x0: number, y0: number, z0: number, x1: number, y1: number, z1: number,
        callback: ImageUtility3Callback): void {
        // Starting point of the line.
        let x = x0;
        let y = y0;
        let z = z0;

        // Direction of the line.
        let dx = x1 - x0;
        let dy = y1 - y0;
        let dz = z1 - z0;

        // Increment or decrement depending on the direction of the line.
        const sx = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
        const sy = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));
        const sz = (dz > 0 ? 1 : (dz < 0 ? -1 : 0));

        // Decision parameters for voxel selection.
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        if (dz < 0) {
            dz = -dz;
        }
        const ax = 2 * dx;
        const ay = 2 * dy;
        const az = 2 * dz;
        let decX = 0;
        let decY = 0;
        let decZ = 0;

        // Determine the largest direction component and the single-step
        // related variable.
        let maxValue = dx;
        let variable = 0;
        if (dy > maxValue) {
            maxValue = dy;
            variable = 1;
        }
        if (dz > maxValue) {
            variable = 2;
        }

        // Traverse the Bresenham line.
        if (variable === 0) {
            // Single-step in the x-direction.
            decY = ay - dx;
            decZ = az - dx;
            for (/**/; /**/; x += sx, decY += ay, decZ += az) {
                // Process the voxel.
                callback(x, y, z);

                // Take a Bresenham step.
                if (x === x1) {
                    break;
                }
                if (decY >= 0) {
                    decY -= ax;
                    y += sy;
                }
                if (decZ >= 0) {
                    decZ -= ax;
                    z += sz;
                }
            }
        } else if (variable === 1) {
            // Single-step in the y-direction.
            decX = ax - dy;
            decZ = az - dy;
            for (/**/; /**/; y += sy, decX += ax, decZ += az) {
                // Process the voxel.
                callback(x, y, z);

                // Take a Bresenham step.
                if (y === y1) {
                    break;
                }
                if (decX >= 0) {
                    decX -= ay;
                    x += sx;
                }
                if (decZ >= 0) {
                    decZ -= ay;
                    z += sz;
                }
            }
        } else {
            // Single-step in the z-direction.
            decX = ax - dz;
            decY = ay - dz;
            for (/**/; /**/; z += sz, decX += ax, decY += ay) {
                // Process the voxel.
                callback(x, y, z);

                // Take a Bresenham step.
                if (z === z1) {
                    break;
                }
                if (decX >= 0) {
                    decX -= az;
                    x += sx;
                }
                if (decY >= 0) {
                    decY -= az;
                    y += sy;
                }
            }
        }
    }
}

// Duplicate an Image3<number>; the upstream operations start from
// 'outImage = inImage' or 'Image3<int32_t> temp = image'.
function copyImage3(image: Image3<number>): Image3<number> {
    const out = new Image3<number>(image.getDimension(0), image.getDimension(1),
        image.getDimension(2));
    const src = image.getPixels();
    const dst = out.getPixels();
    for (let i = 0; i < src.length; ++i) {
        dst[i] = src[i];
    }
    return out;
}
