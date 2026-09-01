// gtengine-js: TypeScript port of Geometric Tools Engine (GTE) ImageUtility2.h
// Upstream: David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2026 David Eberly
// Distributed under the Boost Software License, Version 1.0.
// https://www.boost.org/LICENSE_1_0.txt

// Image utilities for Image2<number> objects whose pixels are integers.
//
// All but the draw* functions are operations on binary images. Let the image
// have d0 columns and d1 rows. The input image must have zeros on its
// boundaries x = 0, x = d0-1, y = 0, and y = d1-1. The 0-valued pixels are
// considered to be background. The 1-valued pixels are considered to be
// foreground. In some of the operations, to save memory and time the input
// image is modified by the algorithms. If you need to preserve the input
// image, make a copy of it before calling these functions. Dilation and
// erosion functions do not have the requirement that the boundary pixels of
// the binary image inputs be zero.
//
// Port notes:
//   - Upstream is a C++ class of static member functions; the port keeps that
//     shape (a class named ImageUtility2 whose members are static). This also
//     keeps names such as getComponents, dilate and drawLine off the module
//     scope, so no 2D/3D suffixes are needed to keep global exports unique
//     (ImageUtility3 owns the 3D versions of the same names).
//   - The C++ overloads are distinguished by the size N of a std::array
//     template argument (N in {4, 8}) versus a (count, pointer) pair.
//     TypeScript has no such distinction, so, as in the Image2 port, N is
//     part of the method name and the caller-specified variants take a single
//     array whose length supplies numNeighbors:
//       GetComponents<N>(image, components)   -> getComponents4/8(image)
//       GetComponents(image, n, nbrs, comps)  -> getComponents(image, nbrs)
//       Dilate<N>(in, out)                    -> dilate4/8(inImage)
//       Dilate(in, n, nbrs, out)              -> dilate(inImage, nbrs)
//       ... and likewise for Erode, Open and Close.
//     The 1-dimensional neighbor offsets for getComponents come from
//     Image2.getNeighborhood4()/getNeighborhood8(); the 2-tuple neighbor
//     offsets for dilate/erode/open/close come from
//     Image2.getNeighborhood4Coords()/getNeighborhood8Coords().
//   - Upstream writes results through non-const reference out-parameters. The
//     port returns them instead: the morphology operations return a newly
//     allocated output image (so the upstream "input and output must be
//     different objects" assertion is satisfied by construction), and the
//     multi-valued outputs of extractBoundary, getL1Distance and
//     getL2Distance become returned object literals. Functions documented as
//     modifying the input image in place (getComponents, extractBoundary,
//     floodFill4, getL1Distance, getSkeleton) still do so.
//   - getL2Distance allocates the float transform image itself rather than
//     requiring a preallocated one from the caller, and rounds the stored
//     distances to float precision with Math.fround to match the upstream
//     Image2<float> transform.

import { Image2 } from './Image2';
import { logAssert } from './Logger';

// The type of the callback invoked by the draw* functions for each visited
// pixel (upstream std::function<void(int32_t, int32_t)>).
export type ImageUtility2Callback = (x: number, y: number) => void;

export class ImageUtility2 {
    // Compute the N-connected components of a binary image (N is 4 or 8).
    // The input image is modified to avoid the cost of making a copy. On
    // output, the image values are the labels for the components. The
    // returned array components[k], k >= 1, contains the indices for the
    // k-th component; components[0] is unused. When there are no components
    // the returned array is empty.
    static getComponents4(image: Image2<number>): number[][] {
        return ImageUtility2.getComponents(image, image.getNeighborhood4());
    }

    static getComponents8(image: Image2<number>): number[][] {
        return ImageUtility2.getComponents(image, image.getNeighborhood8());
    }

    // Connected component labeling using depth-first search. The
    // neighborhood is specified by the caller as 1-dimensional index offsets
    // into the image.
    static getComponents(image: Image2<number>, neighbors: readonly number[]): number[][] {
        const numNeighbors = neighbors.length;
        const numPixels = image.getNumPixels();
        const pixels = image.getPixels();
        const numElements = new Array<number>(numPixels).fill(0);
        const vstack = new Array<number>(numPixels).fill(0);
        let numComponents = 0;
        let label = 2;
        for (let i = 0; i < numPixels; ++i) {
            if (pixels[i] === 1) {
                let top = -1;
                vstack[++top] = i;

                // Upstream binds a reference 'count' to numElements[k]; the
                // port indexes the array directly.
                const k = numComponents + 1;
                numElements[k] = 0;
                while (top >= 0) {
                    const v = vstack[top];
                    pixels[v] = -1;
                    let j: number;
                    for (j = 0; j < numNeighbors; ++j) {
                        const adj = v + neighbors[j];
                        if (pixels[adj] === 1) {
                            vstack[++top] = adj;
                            break;
                        }
                    }
                    if (j === numNeighbors) {
                        pixels[v] = label;
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

            for (let i = 0; i < numPixels; ++i) {
                let value = pixels[i];
                if (value !== 0) {
                    // Labels started at 2 to support the depth-first search,
                    // so they need to be decremented for the correct labels.
                    pixels[i] = --value;
                    components[value][numElements[value]] = i;
                    ++numElements[value];
                }
            }
        }
        return components;
    }

    // Compute a dilation with a structuring element consisting of the
    // N-connected neighbors of each pixel (N is 4 or 8). The input image is
    // binary with 0 for background and 1 for foreground. The dilated image
    // is returned; the input image is not modified.
    static dilate4(inImage: Image2<number>): Image2<number> {
        return ImageUtility2.dilate(inImage, inImage.getNeighborhood4Coords());
    }

    static dilate8(inImage: Image2<number>): Image2<number> {
        return ImageUtility2.dilate(inImage, inImage.getNeighborhood8Coords());
    }

    // Compute a dilation with a structuring element consisting of neighbors
    // specified by 2-tuple offsets relative to the pixel.
    static dilate(inImage: Image2<number>,
        neighbors: readonly (readonly number[])[]): Image2<number> {
        logAssert(neighbors.length > 0, 'Invalid neighbors.');

        const outImage = copyImage2(inImage);

        // If the pixel at (i0,i1) is 1, then the pixels at
        // (k0,k1) = (i0+nbr0,i1+nbr1) are set to 1 where (nbr0,nbr1) is in
        // the neighbors array. Boundary testing is used to avoid accessing
        // out-of-range pixels.
        const dim0 = inImage.getDimension(0);
        const dim1 = inImage.getDimension(1);
        for (let i1 = 0; i1 < dim1; ++i1) {
            for (let i0 = 0; i0 < dim0; ++i0) {
                if (inImage.get(i0, i1) === 1) {
                    for (let j = 0; j < neighbors.length; ++j) {
                        const k0 = i0 + neighbors[j][0];
                        const k1 = i1 + neighbors[j][1];
                        if (0 <= k0 && k0 < dim0 && 0 <= k1 && k1 < dim1) {
                            outImage.set(k0, k1, 1);
                        }
                    }
                }
            }
        }
        return outImage;
    }

    // Compute an erosion with a structuring element consisting of the
    // N-connected neighbors of each pixel (N is 4 or 8). The input image is
    // binary with 0 for background and 1 for foreground. If zeroExterior is
    // true, the image exterior is assumed to be 0, so 1-valued boundary
    // pixels are set to 0; otherwise, boundary pixels are set to 0 only when
    // they have neighboring image pixels that are 0. The eroded image is
    // returned; the input image is not modified.
    static erode4(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.erode(inImage, zeroExterior, inImage.getNeighborhood4Coords());
    }

    static erode8(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.erode(inImage, zeroExterior, inImage.getNeighborhood8Coords());
    }

    // Compute an erosion with a structuring element consisting of neighbors
    // specified by 2-tuple offsets relative to the pixel.
    static erode(inImage: Image2<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image2<number> {
        logAssert(neighbors.length > 0, 'Invalid neighbors.');

        const outImage = copyImage2(inImage);

        // If the pixel at (i0,i1) is 1, it is changed to 0 when at least one
        // neighbor (k0,k1) = (i0+nbr0,i1+nbr1) is 0, where (nbr0,nbr1) is in
        // the neighbors array.
        const dim0 = inImage.getDimension(0);
        const dim1 = inImage.getDimension(1);
        for (let i1 = 0; i1 < dim1; ++i1) {
            for (let i0 = 0; i0 < dim0; ++i0) {
                if (inImage.get(i0, i1) === 1) {
                    for (let j = 0; j < neighbors.length; ++j) {
                        const k0 = i0 + neighbors[j][0];
                        const k1 = i1 + neighbors[j][1];
                        if (0 <= k0 && k0 < dim0 && 0 <= k1 && k1 < dim1) {
                            if (inImage.get(k0, k1) === 0) {
                                outImage.set(i0, i1, 0);
                                break;
                            }
                        } else if (zeroExterior) {
                            outImage.set(i0, i1, 0);
                            break;
                        }
                    }
                }
            }
        }
        return outImage;
    }

    // Compute an opening (erosion followed by dilation) with a structuring
    // element consisting of the N-connected neighbors of each pixel (N is 4
    // or 8). See erode for the meaning of zeroExterior.
    static open4(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.dilate4(ImageUtility2.erode4(inImage, zeroExterior));
    }

    static open8(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.dilate8(ImageUtility2.erode8(inImage, zeroExterior));
    }

    // Compute an opening with a caller-specified structuring element.
    static open(inImage: Image2<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image2<number> {
        const temp = ImageUtility2.erode(inImage, zeroExterior, neighbors);
        return ImageUtility2.dilate(temp, neighbors);
    }

    // Compute a closing (dilation followed by erosion) with a structuring
    // element consisting of the N-connected neighbors of each pixel (N is 4
    // or 8). See erode for the meaning of zeroExterior.
    static close4(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.erode4(ImageUtility2.dilate4(inImage), zeroExterior);
    }

    static close8(inImage: Image2<number>, zeroExterior: boolean): Image2<number> {
        return ImageUtility2.erode8(ImageUtility2.dilate8(inImage), zeroExterior);
    }

    // Compute a closing with a caller-specified structuring element.
    static close(inImage: Image2<number>, zeroExterior: boolean,
        neighbors: readonly (readonly number[])[]): Image2<number> {
        const temp = ImageUtility2.dilate(inImage, neighbors);
        return ImageUtility2.erode(temp, zeroExterior, neighbors);
    }

    // Locate a pixel and walk around the edge of a component. The input
    // (x,y) is where the search starts for a nonzero pixel. If (x,y) is
    // outside the component, the walk is around the outside of the
    // component. If the component has a hole and (x,y) is inside that hole,
    // the walk is around the boundary surrounding the hole. The returned
    // 'success' is true for a successful walk; it is false when no boundary
    // was found from the starting (x,y), in which case 'boundary' is empty.
    // Visited pixels of the image are marked with the value 2.
    static extractBoundary(x: number, y: number, image: Image2<number>):
        { success: boolean, boundary: number[] } {
        const boundary: number[] = [];

        // Find a first boundary pixel.
        const numPixels = image.getNumPixels();
        let i = image.getIndex(x, y);
        for (/**/; i < numPixels; ++i) {
            if (image.get(i) !== 0) {
                break;
            }
        }
        if (i === numPixels) {
            // No boundary pixel found.
            return { success: false, boundary };
        }

        const dx = [-1, 0, +1, +1, +1, 0, -1, -1];
        const dy = [-1, -1, -1, 0, +1, +1, +1, 0];

        // Create a new point list that contains the first boundary point.
        boundary.push(i);

        // The direction from background 0 to boundary pixel 1 is
        // (dx[7],dy[7]).
        const coord = image.getCoordinates(i);
        const x0 = coord[0];
        const y0 = coord[1];
        let cx = x0;
        let cy = y0;
        let nx = x0 - 1;
        let ny = y0;
        let dir = 7;

        // Traverse the boundary in clockwise order. Mark visited pixels as 2.
        image.set(cx, cy, 2);
        let notDone = true;
        while (notDone) {
            let j = 0;
            let nbr = dir;
            for (j = 0, nbr = dir; j < 8; ++j, nbr = (nbr + 1) % 8) {
                nx = cx + dx[nbr];
                ny = cy + dy[nbr];
                if (image.get(nx, ny) !== 0) {
                    // Next boundary pixel found.
                    break;
                }
            }

            if (j === 8) {
                // (cx,cy) is isolated.
                notDone = false;
                continue;
            }

            if (nx === x0 && ny === y0) {
                // Boundary traversal completed.
                notDone = false;
                continue;
            }

            // (nx,ny) is the next boundary point; add it to the list.
            boundary.push(image.getIndex(nx, ny));

            // Mark visited pixels as 2.
            image.set(nx, ny, 2);

            // Start the search for the next point.
            cx = nx;
            cy = ny;
            dir = (j + 5 + dir) % 8;
        }

        return { success: true, boundary };
    }

    // Use a depth-first search for filling a 4-connected region. This is
    // nonrecursive, simulated by using a heap-allocated "stack". The input
    // (x,y) is the seed point that starts the fill.
    static floodFill4<PixelType>(image: Image2<PixelType>, x: number, y: number,
        foreColor: PixelType, backColor: PixelType): void {
        // Test for a valid seed.
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);
        if (x < 0 || x >= dim0 || y < 0 || y >= dim1) {
            // The seed point is outside the image domain, so there is
            // nothing to fill.
            return;
        }

        // Allocate the maximum amount of space needed for the stack. An
        // empty stack has top == -1.
        const numPixels = image.getNumPixels();
        const xStack = new Array<number>(numPixels).fill(0);
        const yStack = new Array<number>(numPixels).fill(0);

        // Push the seed point onto the stack. All points pushed onto the
        // stack have the background color backColor.
        let top = 0;
        xStack[top] = x;
        yStack[top] = y;

        while (top >= 0) {
            // Read the top of the stack. Do not pop, because we need to
            // return to this top value later to restart the fill in a
            // different direction.
            x = xStack[top];
            y = yStack[top];

            // Fill the pixel.
            image.set(x, y, foreColor);

            const xp1 = x + 1;
            if (xp1 < dim0 && image.get(xp1, y) === backColor) {
                ++top;
                xStack[top] = xp1;
                yStack[top] = y;
                continue;
            }

            const xm1 = x - 1;
            if (0 <= xm1 && image.get(xm1, y) === backColor) {
                ++top;
                xStack[top] = xm1;
                yStack[top] = y;
                continue;
            }

            const yp1 = y + 1;
            if (yp1 < dim1 && image.get(x, yp1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = yp1;
                continue;
            }

            const ym1 = y - 1;
            if (0 <= ym1 && image.get(x, ym1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = ym1;
                continue;
            }

            // Done in all directions; pop and return to search the other
            // directions of the predecessor.
            --top;
        }
    }

    // Compute the L1-distance transform of the binary image, in place. The
    // function returns the maximum distance and a point at which the maximum
    // distance is attained.
    static getL1Distance(image: Image2<number>):
        { maxDistance: number, xMax: number, yMax: number } {
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);
        const dim0m1 = dim0 - 1;
        const dim1m1 = dim1 - 1;

        // Use a grass-fire approach, computing the distance from boundary to
        // interior one pass at a time.
        let changeMade = true;
        let distance = 1;
        let xMax = 0;
        let yMax = 0;
        for (/**/; changeMade; ++distance) {
            changeMade = false;
            const distanceP1 = distance + 1;
            for (let y = 1; y < dim1m1; ++y) {
                for (let x = 1; x < dim0m1; ++x) {
                    if (image.get(x, y) === distance) {
                        if (image.get(x - 1, y) >= distance
                            && image.get(x + 1, y) >= distance
                            && image.get(x, y - 1) >= distance
                            && image.get(x, y + 1) >= distance) {
                            image.set(x, y, distanceP1);
                            xMax = x;
                            yMax = y;
                            changeMade = true;
                        }
                    }
                }
            }
        }

        return { maxDistance: --distance, xMax, yMax };
    }

    // Compute the L2-distance transform of the binary image. The maximum
    // distance should not be larger than 100, so you have to ensure this is
    // the case for the input image. The function returns the maximum
    // distance, a point at which the maximum distance is attained, and the
    // float-valued distance transform image. The input image is not
    // modified.
    //
    // This computes the Euclidean distance transform of a binary input
    // image. The adaptive algorithm is guaranteed to give exact distances
    // for all distances < 100. The algorithm was provided by John Gauch at
    // the University of Kansas. The following is a quote:
    //
    // The basic idea is similar to a EDT described recently in PAMI by
    // Laymarie from McGill. By keeping the dx and dy offset to the nearest
    // edge (feature) point in the image, we can search to see which dx dy is
    // closest to a given point by examining a set of neighbors. The Laymarie
    // method (and Borgfors) look at a fixed 3x3 or 5x5 neighborhood and call
    // it a day. What we did was calculate (painfully) what neighborhoods you
    // need to look at to guarantee that the exact distance is obtained.
    // Thus, you will see in the code, that we check the current distance and
    // depending on what we have so far, we extend the search region. Since
    // our algorithm for checking the exactness of each neighborhood is on the
    // order N^4, we have only gone to N=100. In theory, you could make this
    // large enough to get all distances exact. We have implemented the
    // algorithm to get all distances < 100 to be exact.
    static getL2Distance(image: Image2<number>):
        { maxDistance: number, xMax: number, yMax: number, transform: Image2<number> } {
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);
        const dim0m1 = dim0 - 1;
        const dim1m1 = dim1 - 1;
        let x = 0;
        let y = 0;
        let distance = 0;

        // Create and initialize intermediate images.
        const xNear = new Image2<number>(dim0, dim1);
        const yNear = new Image2<number>(dim0, dim1);
        const dist = new Image2<number>(dim0, dim1);
        for (y = 0; y < dim1; ++y) {
            for (x = 0; x < dim0; ++x) {
                if (image.get(x, y) !== 0) {
                    xNear.set(x, y, 0);
                    yNear.set(x, y, 0);
                    dist.set(x, y, INT32_MAX);
                } else {
                    xNear.set(x, y, x);
                    yNear.set(x, y, y);
                    dist.set(x, y, 0);
                }
            }
        }

        const K1 = 1;
        const K2 = 169;   // 13^2
        const K3 = 961;   // 31^2
        const K4 = 2401;  // 49^2
        const K5 = 5184;  // 72^2

        // Pass in the ++ direction.
        for (y = 0; y < dim1; ++y) {
            for (x = 0; x < dim0; ++x) {
                distance = dist.get(x, y);
                if (distance > K1) {
                    l2Check(x, y, -1, 0, xNear, yNear, dist);
                    l2Check(x, y, -1, -1, xNear, yNear, dist);
                    l2Check(x, y, 0, -1, xNear, yNear, dist);
                }
                if (distance > K2) {
                    l2Check(x, y, -2, -1, xNear, yNear, dist);
                    l2Check(x, y, -1, -2, xNear, yNear, dist);
                }
                if (distance > K3) {
                    l2Check(x, y, -3, -1, xNear, yNear, dist);
                    l2Check(x, y, -3, -2, xNear, yNear, dist);
                    l2Check(x, y, -2, -3, xNear, yNear, dist);
                    l2Check(x, y, -1, -3, xNear, yNear, dist);
                }
                if (distance > K4) {
                    l2Check(x, y, -4, -1, xNear, yNear, dist);
                    l2Check(x, y, -4, -3, xNear, yNear, dist);
                    l2Check(x, y, -3, -4, xNear, yNear, dist);
                    l2Check(x, y, -1, -4, xNear, yNear, dist);
                }
                if (distance > K5) {
                    l2Check(x, y, -5, -1, xNear, yNear, dist);
                    l2Check(x, y, -5, -2, xNear, yNear, dist);
                    l2Check(x, y, -5, -3, xNear, yNear, dist);
                    l2Check(x, y, -5, -4, xNear, yNear, dist);
                    l2Check(x, y, -4, -5, xNear, yNear, dist);
                    l2Check(x, y, -2, -5, xNear, yNear, dist);
                    l2Check(x, y, -3, -5, xNear, yNear, dist);
                    l2Check(x, y, -1, -5, xNear, yNear, dist);
                }
            }
        }

        // Pass in the -- direction.
        for (y = dim1m1; y >= 0; --y) {
            for (x = dim0m1; x >= 0; --x) {
                distance = dist.get(x, y);
                if (distance > K1) {
                    l2Check(x, y, 1, 0, xNear, yNear, dist);
                    l2Check(x, y, 1, 1, xNear, yNear, dist);
                    l2Check(x, y, 0, 1, xNear, yNear, dist);
                }
                if (distance > K2) {
                    l2Check(x, y, 2, 1, xNear, yNear, dist);
                    l2Check(x, y, 1, 2, xNear, yNear, dist);
                }
                if (distance > K3) {
                    l2Check(x, y, 3, 1, xNear, yNear, dist);
                    l2Check(x, y, 3, 2, xNear, yNear, dist);
                    l2Check(x, y, 2, 3, xNear, yNear, dist);
                    l2Check(x, y, 1, 3, xNear, yNear, dist);
                }
                if (distance > K4) {
                    l2Check(x, y, 4, 1, xNear, yNear, dist);
                    l2Check(x, y, 4, 3, xNear, yNear, dist);
                    l2Check(x, y, 3, 4, xNear, yNear, dist);
                    l2Check(x, y, 1, 4, xNear, yNear, dist);
                }
                if (distance > K5) {
                    l2Check(x, y, 5, 1, xNear, yNear, dist);
                    l2Check(x, y, 5, 2, xNear, yNear, dist);
                    l2Check(x, y, 5, 3, xNear, yNear, dist);
                    l2Check(x, y, 5, 4, xNear, yNear, dist);
                    l2Check(x, y, 4, 5, xNear, yNear, dist);
                    l2Check(x, y, 2, 5, xNear, yNear, dist);
                    l2Check(x, y, 3, 5, xNear, yNear, dist);
                    l2Check(x, y, 1, 5, xNear, yNear, dist);
                }
            }
        }

        // Pass in the +- direction.
        for (y = dim1m1; y >= 0; --y) {
            for (x = 0; x < dim0; ++x) {
                distance = dist.get(x, y);
                if (distance > K1) {
                    l2Check(x, y, -1, 0, xNear, yNear, dist);
                    l2Check(x, y, -1, 1, xNear, yNear, dist);
                    l2Check(x, y, 0, 1, xNear, yNear, dist);
                }
                if (distance > K2) {
                    l2Check(x, y, -2, 1, xNear, yNear, dist);
                    l2Check(x, y, -1, 2, xNear, yNear, dist);
                }
                if (distance > K3) {
                    l2Check(x, y, -3, 1, xNear, yNear, dist);
                    l2Check(x, y, -3, 2, xNear, yNear, dist);
                    l2Check(x, y, -2, 3, xNear, yNear, dist);
                    l2Check(x, y, -1, 3, xNear, yNear, dist);
                }
                if (distance > K4) {
                    l2Check(x, y, -4, 1, xNear, yNear, dist);
                    l2Check(x, y, -4, 3, xNear, yNear, dist);
                    l2Check(x, y, -3, 4, xNear, yNear, dist);
                    l2Check(x, y, -1, 4, xNear, yNear, dist);
                }
                if (distance > K5) {
                    l2Check(x, y, -5, 1, xNear, yNear, dist);
                    l2Check(x, y, -5, 2, xNear, yNear, dist);
                    l2Check(x, y, -5, 3, xNear, yNear, dist);
                    l2Check(x, y, -5, 4, xNear, yNear, dist);
                    l2Check(x, y, -4, 5, xNear, yNear, dist);
                    l2Check(x, y, -2, 5, xNear, yNear, dist);
                    l2Check(x, y, -3, 5, xNear, yNear, dist);
                    l2Check(x, y, -1, 5, xNear, yNear, dist);
                }
            }
        }

        // Pass in the -+ direction.
        for (y = 0; y < dim1; ++y) {
            for (x = dim0m1; x >= 0; --x) {
                distance = dist.get(x, y);
                if (distance > K1) {
                    l2Check(x, y, 1, 0, xNear, yNear, dist);
                    l2Check(x, y, 1, -1, xNear, yNear, dist);
                    l2Check(x, y, 0, -1, xNear, yNear, dist);
                }
                if (distance > K2) {
                    l2Check(x, y, 2, -1, xNear, yNear, dist);
                    l2Check(x, y, 1, -2, xNear, yNear, dist);
                }
                if (distance > K3) {
                    l2Check(x, y, 3, -1, xNear, yNear, dist);
                    l2Check(x, y, 3, -2, xNear, yNear, dist);
                    l2Check(x, y, 2, -3, xNear, yNear, dist);
                    l2Check(x, y, 1, -3, xNear, yNear, dist);
                }
                if (distance > K4) {
                    l2Check(x, y, 4, -1, xNear, yNear, dist);
                    l2Check(x, y, 4, -3, xNear, yNear, dist);
                    l2Check(x, y, 3, -4, xNear, yNear, dist);
                    l2Check(x, y, 1, -4, xNear, yNear, dist);
                }
                if (distance > K5) {
                    l2Check(x, y, 5, -1, xNear, yNear, dist);
                    l2Check(x, y, 5, -2, xNear, yNear, dist);
                    l2Check(x, y, 5, -3, xNear, yNear, dist);
                    l2Check(x, y, 5, -4, xNear, yNear, dist);
                    l2Check(x, y, 4, -5, xNear, yNear, dist);
                    l2Check(x, y, 2, -5, xNear, yNear, dist);
                    l2Check(x, y, 3, -5, xNear, yNear, dist);
                    l2Check(x, y, 1, -5, xNear, yNear, dist);
                }
            }
        }

        const transform = new Image2<number>(dim0, dim1);
        let xMax = 0;
        let yMax = 0;
        let maxDistance = 0;
        for (y = 0; y < dim1; ++y) {
            for (x = 0; x < dim0; ++x) {
                const fdistance = Math.fround(Math.sqrt(dist.get(x, y)));
                if (fdistance > maxDistance) {
                    maxDistance = fdistance;
                    xMax = x;
                    yMax = y;
                }
                transform.set(x, y, fdistance);
            }
        }

        return { maxDistance, xMax, yMax, transform };
    }

    // Compute a skeleton of a binary image. Boundary pixels are trimmed from
    // the object one layer at a time based on their adjacency to interior
    // pixels. At each step the connectivity and cycles of the object are
    // preserved. The skeleton overwrites the contents of the input image.
    static getSkeleton(image: Image2<number>): void {
        const dim0 = image.getDimension(0);
        const dim1 = image.getDimension(1);

        // Trim pixels, mark interior as 4.
        let notDone = true;
        while (notDone) {
            if (markInterior(image, 4, interior4)) {
                // No interior pixels; the trimmed set is at most 2 pixels
                // thick.
                notDone = false;
                continue;
            }

            if (clearInteriorAdjacent(image, 4)) {
                // All remaining interior pixels are either articulation
                // points or part of blobs whose boundary pixels are all
                // articulation points. An example of the latter case is shown
                // below. The background pixels are marked with '.' rather
                // than '0' for readability. The interior pixels are marked
                // with '4' and the boundary pixels are marked with '1'.
                //
                //   .........
                //   .....1...
                //   ..1.1.1..
                //   .1.141...
                //   ..14441..
                //   ..1441.1.
                //   .1.11.1..
                //   ..1..1...
                //   .........
                //
                // This is a pathological problem where there are many small
                // holes (a 0-pixel with north, south, west and east neighbors
                // all 1-pixels) that your application can try to avoid by an
                // initial pass over the image to fill in such holes. Of
                // course, you do have problems with checkerboard patterns.
                notDone = false;
                continue;
            }
        }

        // Trim pixels, mark interior as 3.
        notDone = true;
        while (notDone) {
            if (markInterior(image, 3, interior3)) {
                // No interior pixels; the trimmed set is at most 2 pixels
                // thick.
                notDone = false;
                continue;
            }

            if (clearInteriorAdjacent(image, 3)) {
                // All remaining 3-values can be safely removed since they are
                // not articulation points and the removal will not cause new
                // holes.
                for (let y = 0; y < dim1; ++y) {
                    for (let x = 0; x < dim0; ++x) {
                        if (image.get(x, y) === 3 && !isArticulation(image, x, y)) {
                            image.set(x, y, 0);
                        }
                    }
                }
                notDone = false;
                continue;
            }
        }

        // Trim pixels, mark interior as 2.
        notDone = true;
        while (notDone) {
            if (markInterior(image, 2, interior2)) {
                // No interior pixels; the trimmed set is at most 1 pixel
                // thick. Call it a skeleton.
                notDone = false;
                continue;
            }

            if (clearInteriorAdjacent(image, 2)) {
                // Remove 2-values that are not articulation points.
                for (let y = 0; y < dim1; ++y) {
                    for (let x = 0; x < dim0; ++x) {
                        if (image.get(x, y) === 2 && !isArticulation(image, x, y)) {
                            image.set(x, y, 0);
                        }
                    }
                }
                notDone = false;
                continue;
            }
        }

        // Make the skeleton a binary image.
        const pixels = image.getPixels();
        for (let i = 0; i < pixels.length; ++i) {
            if (pixels[i] !== 0) {
                pixels[i] = 1;
            }
        }
    }

    // In the remaining member functions, the callback represents the action
    // you want applied to each pixel as it is visited.

    // Visit pixels in a (2*thick+1)x(2*thick+1) square centered at (x,y).
    static drawThickPixel(x: number, y: number, thick: number,
        callback: ImageUtility2Callback): void {
        for (let dy = -thick; dy <= thick; ++dy) {
            for (let dx = -thick; dx <= thick; ++dx) {
                callback(x + dx, y + dy);
            }
        }
    }

    // Visit pixels using Bresenham's line drawing algorithm.
    static drawLine(x0: number, y0: number, x1: number, y1: number,
        callback: ImageUtility2Callback): void {
        // Starting point of the line.
        let x = x0;
        let y = y0;

        // Direction of the line.
        let dx = x1 - x0;
        let dy = y1 - y0;

        // Increment or decrement depending on the direction of the line.
        const sx = (dx > 0 ? 1 : (dx < 0 ? -1 : 0));
        const sy = (dy > 0 ? 1 : (dy < 0 ? -1 : 0));

        // Decision parameters for pixel selection.
        if (dx < 0) {
            dx = -dx;
        }
        if (dy < 0) {
            dy = -dy;
        }
        const ax = 2 * dx;
        const ay = 2 * dy;

        // Determine the largest direction component and the single-step
        // related variable.
        const maxValue = dx;
        let variable = 0;
        if (dy > maxValue) {
            variable = 1;
        }

        // Traverse the Bresenham line.
        if (variable === 0) {
            // Single-step in the x-direction.
            let decY = ay - dx;
            for (/**/; /**/; x += sx, decY += ay) {
                callback(x, y);

                // Take a Bresenham step.
                if (x === x1) {
                    break;
                }
                if (decY >= 0) {
                    decY -= ax;
                    y += sy;
                }
            }
        } else {
            // Single-step in the y-direction.
            let decX = ax - dy;
            for (/**/; /**/; y += sy, decX += ax) {
                callback(x, y);

                // Take a Bresenham step.
                if (y === y1) {
                    break;
                }
                if (decX >= 0) {
                    decX -= ay;
                    x += sx;
                }
            }
        }
    }

    // Visit pixels using Bresenham's circle drawing algorithm. Set 'solid' to
    // false for drawing only the circle. Set 'solid' to true to draw all
    // pixels on and inside the circle.
    static drawCircle(xCenter: number, yCenter: number, radius: number, solid: boolean,
        callback: ImageUtility2Callback): void {
        if (solid) {
            for (let x = 0, y = radius, dec = 3 - 2 * radius; x <= y; ++x) {
                let xMin = xCenter - x;
                let xMax = xCenter + x;
                let yValue = yCenter - y;
                for (let xValue = xMin; xValue <= xMax; ++xValue) {
                    callback(xValue, yValue);
                }

                yValue = yCenter + y;
                for (let xValue = xMin; xValue <= xMax; ++xValue) {
                    callback(xValue, yValue);
                }

                xMin = xCenter - y;
                xMax = xCenter + y;
                yValue = yCenter - x;
                for (let xValue = xMin; xValue <= xMax; ++xValue) {
                    callback(xValue, yValue);
                }

                yValue = yCenter + x;
                for (let xValue = xMin; xValue <= xMax; ++xValue) {
                    callback(xValue, yValue);
                }

                if (dec >= 0) {
                    dec += -4 * y + 4;
                    --y;
                }
                dec += 4 * x + 6;
            }
        } else {
            for (let x = 0, y = radius, dec = 3 - 2 * radius; x <= y; ++x) {
                callback(xCenter + x, yCenter + y);
                callback(xCenter + x, yCenter - y);
                callback(xCenter - x, yCenter + y);
                callback(xCenter - x, yCenter - y);
                callback(xCenter + y, yCenter + x);
                callback(xCenter + y, yCenter - x);
                callback(xCenter - y, yCenter + x);
                callback(xCenter - y, yCenter - x);

                if (dec >= 0) {
                    dec += -4 * y + 4;
                    --y;
                }
                dec += 4 * x + 6;
            }
        }
    }

    // Visit pixels in a rectangle of the specified dimensions. Set 'solid' to
    // false for drawing only the rectangle. Set 'solid' to true to draw all
    // pixels on and inside the rectangle.
    static drawRectangle(xMin: number, yMin: number, xMax: number, yMax: number,
        solid: boolean, callback: ImageUtility2Callback): void {
        if (solid) {
            for (let y = yMin; y <= yMax; ++y) {
                for (let x = xMin; x <= xMax; ++x) {
                    callback(x, y);
                }
            }
        } else {
            for (let x = xMin; x <= xMax; ++x) {
                callback(x, yMin);
                callback(x, yMax);
            }
            for (let y = yMin + 1; y <= yMax - 1; ++y) {
                callback(xMin, y);
                callback(xMax, y);
            }
        }
    }

    // Visit the pixels using Bresenham's algorithm for the axis-aligned
    // ellipse ((x-xc)/a)^2 + ((y-yc)/b)^2 = 1, where xCenter is xc, yCenter
    // is yc, xExtent is a, and yExtent is b.
    static drawEllipse(xCenter: number, yCenter: number, xExtent: number, yExtent: number,
        callback: ImageUtility2Callback): void {
        const xExtSqr = xExtent * xExtent;
        const yExtSqr = yExtent * yExtent;
        let x = 0;
        let y = yExtent;
        let dec = 2 * yExtSqr + xExtSqr * (1 - 2 * yExtent);
        for (/**/; yExtSqr * x <= xExtSqr * y; ++x) {
            callback(xCenter + x, yCenter + y);
            callback(xCenter - x, yCenter + y);
            callback(xCenter + x, yCenter - y);
            callback(xCenter - x, yCenter - y);

            if (dec >= 0) {
                dec += 4 * xExtSqr * (1 - y);
                --y;
            }
            dec += yExtSqr * (4 * x + 6);
        }
        if (y === 0 && x < xExtent) {
            // The discretization caused us to reach the y-axis before the
            // x-values reached the ellipse vertices. Draw a solid line along
            // the x-axis to those vertices.
            for (/**/; x <= xExtent; ++x) {
                callback(xCenter + x, yCenter);
                callback(xCenter - x, yCenter);
            }
            return;
        }

        x = xExtent;
        y = 0;
        dec = 2 * xExtSqr + yExtSqr * (1 - 2 * xExtent);
        for (/**/; xExtSqr * y <= yExtSqr * x; ++y) {
            callback(xCenter + x, yCenter + y);
            callback(xCenter - x, yCenter + y);
            callback(xCenter + x, yCenter - y);
            callback(xCenter - x, yCenter - y);

            if (dec >= 0) {
                dec += 4 * yExtSqr * (1 - x);
                --x;
            }
            dec += xExtSqr * (4 * y + 6);
        }
        if (x === 0 && y < yExtent) {
            // The discretization caused us to reach the x-axis before the
            // y-values reached the ellipse vertices. Draw a solid line along
            // the y-axis to those vertices.
            for (/**/; y <= yExtent; ++y) {
                callback(xCenter, yCenter + y);
                callback(xCenter, yCenter - y);
            }
        }
    }

    // Use a depth-first search for filling a 4-connected region. This is
    // nonrecursive, simulated by using a heap-allocated "stack". The input
    // (x,y) is the seed point that starts the fill. The x-value is in
    // {0..xSize-1} and the y-value is in {0..ySize-1}. Pixel reads and writes
    // are performed by the caller-supplied callbacks.
    static drawFloodFill4<PixelType>(x: number, y: number, xSize: number, ySize: number,
        foreColor: PixelType, backColor: PixelType,
        setCallback: (x: number, y: number, value: PixelType) => void,
        getCallback: (x: number, y: number) => PixelType): void {
        // Test for a valid seed.
        if (x < 0 || x >= xSize || y < 0 || y >= ySize) {
            // The seed point is outside the image domain, so nothing to fill.
            return;
        }

        // Allocate the maximum amount of space needed for the stack. An empty
        // stack has top == -1.
        const numPixels = xSize * ySize;
        const xStack = new Array<number>(numPixels).fill(0);
        const yStack = new Array<number>(numPixels).fill(0);

        // Push the seed point onto the stack. All points pushed onto the
        // stack have the background color backColor.
        let top = 0;
        xStack[top] = x;
        yStack[top] = y;

        while (top >= 0) {
            // Read the top of the stack. Do not pop, because we need to
            // return to this top value later to restart the fill in a
            // different direction.
            x = xStack[top];
            y = yStack[top];

            // Fill the pixel.
            setCallback(x, y, foreColor);

            const xp1 = x + 1;
            if (xp1 < xSize && getCallback(xp1, y) === backColor) {
                ++top;
                xStack[top] = xp1;
                yStack[top] = y;
                continue;
            }

            const xm1 = x - 1;
            if (0 <= xm1 && getCallback(xm1, y) === backColor) {
                ++top;
                xStack[top] = xm1;
                yStack[top] = y;
                continue;
            }

            const yp1 = y + 1;
            if (yp1 < ySize && getCallback(x, yp1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = yp1;
                continue;
            }

            const ym1 = y - 1;
            if (0 <= ym1 && getCallback(x, ym1) === backColor) {
                ++top;
                xStack[top] = x;
                yStack[top] = ym1;
                continue;
            }

            // Done in all directions; pop and return to search the other
            // directions of the predecessor.
            --top;
        }
    }
}

// The value of std::numeric_limits<int32_t>::max(), used to initialize the
// distances of foreground pixels in getL2Distance.
const INT32_MAX = 2147483647;

// Duplicate an Image2<number>; the upstream morphology operations start from
// 'outImage = inImage'.
function copyImage2(image: Image2<number>): Image2<number> {
    const out = new Image2<number>(image.getDimension(0), image.getDimension(1));
    const src = image.getPixels();
    const dst = out.getPixels();
    for (let i = 0; i < src.length; ++i) {
        dst[i] = src[i];
    }
    return out;
}

// Support for getL2Distance.
function l2Check(x: number, y: number, dx: number, dy: number,
    xNear: Image2<number>, yNear: Image2<number>, dist: Image2<number>): void {
    const dim0 = dist.getDimension(0);
    const dim1 = dist.getDimension(1);
    const xp = x + dx;
    const yp = y + dy;
    if (0 <= xp && xp < dim0 && 0 <= yp && yp < dim1) {
        if (dist.get(xp, yp) < dist.get(x, y)) {
            const dx0 = xNear.get(xp, yp) - x;
            const dy0 = yNear.get(xp, yp) - y;
            const newDist = dx0 * dx0 + dy0 * dy0;
            if (newDist < dist.get(x, y)) {
                xNear.set(x, y, xNear.get(xp, yp));
                yNear.set(x, y, yNear.get(xp, yp));
                dist.set(x, y, newDist);
            }
        }
    }
}

// Support for getSkeleton. The 'interior' predicates and the callers below
// read the 8 neighbors of (x,y) without range testing, which is safe because
// the input image is required to have zero boundary pixels.
type InteriorPredicate = (image: Image2<number>, x: number, y: number) => boolean;

function interior2(image: Image2<number>, x: number, y: number): boolean {
    const b1 = (image.get(x, y - 1) !== 0);
    const b3 = (image.get(x + 1, y) !== 0);
    const b5 = (image.get(x, y + 1) !== 0);
    const b7 = (image.get(x - 1, y) !== 0);
    return (b1 && b3) || (b3 && b5) || (b5 && b7) || (b7 && b1);
}

function interior3(image: Image2<number>, x: number, y: number): boolean {
    let numNeighbors = 0;
    if (image.get(x - 1, y) !== 0) {
        ++numNeighbors;
    }
    if (image.get(x + 1, y) !== 0) {
        ++numNeighbors;
    }
    if (image.get(x, y - 1) !== 0) {
        ++numNeighbors;
    }
    if (image.get(x, y + 1) !== 0) {
        ++numNeighbors;
    }
    return numNeighbors === 3;
}

function interior4(image: Image2<number>, x: number, y: number): boolean {
    return image.get(x - 1, y) !== 0
        && image.get(x + 1, y) !== 0
        && image.get(x, y - 1) !== 0
        && image.get(x, y + 1) !== 0;
}

function markInterior(image: Image2<number>, value: number,
    predicate: InteriorPredicate): boolean {
    const dim0 = image.getDimension(0);
    const dim1 = image.getDimension(1);
    let noInterior = true;
    for (let y = 0; y < dim1; ++y) {
        for (let x = 0; x < dim0; ++x) {
            if (image.get(x, y) > 0) {
                if (predicate(image, x, y)) {
                    image.set(x, y, value);
                    noInterior = false;
                } else {
                    image.set(x, y, 1);
                }
            }
        }
    }
    return noInterior;
}

// A lookup table indexed by the 8-bit neighbor mask of a pixel; the entry is
// 1 when removing the pixel would disconnect its neighborhood.
const ARTICULATION: readonly number[] = [
    0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
];

function isArticulation(image: Image2<number>, x: number, y: number): boolean {
    // Convert the 8 neighbors of pixel (x,y) to an 8-bit value; a bit is 1
    // if and only if the corresponding pixel is set.
    let byteMask = 0;
    if (image.get(x - 1, y - 1) !== 0) {
        byteMask |= 0x01;
    }
    if (image.get(x, y - 1) !== 0) {
        byteMask |= 0x02;
    }
    if (image.get(x + 1, y - 1) !== 0) {
        byteMask |= 0x04;
    }
    if (image.get(x + 1, y) !== 0) {
        byteMask |= 0x08;
    }
    if (image.get(x + 1, y + 1) !== 0) {
        byteMask |= 0x10;
    }
    if (image.get(x, y + 1) !== 0) {
        byteMask |= 0x20;
    }
    if (image.get(x - 1, y + 1) !== 0) {
        byteMask |= 0x40;
    }
    if (image.get(x - 1, y) !== 0) {
        byteMask |= 0x80;
    }

    return ARTICULATION[byteMask] === 1;
}

function clearInteriorAdjacent(image: Image2<number>, value: number): boolean {
    const dim0 = image.getDimension(0);
    const dim1 = image.getDimension(1);
    let noRemoval = true;
    for (let y = 0; y < dim1; ++y) {
        for (let x = 0; x < dim0; ++x) {
            if (image.get(x, y) === 1) {
                const interiorAdjacent =
                    image.get(x - 1, y - 1) === value ||
                    image.get(x, y - 1) === value ||
                    image.get(x + 1, y - 1) === value ||
                    image.get(x + 1, y) === value ||
                    image.get(x + 1, y + 1) === value ||
                    image.get(x, y + 1) === value ||
                    image.get(x - 1, y + 1) === value ||
                    image.get(x - 1, y) === value;

                if (interiorAdjacent && !isArticulation(image, x, y)) {
                    image.set(x, y, 0);
                    noRemoval = false;
                }
            }
        }
    }
    return noRemoval;
}
