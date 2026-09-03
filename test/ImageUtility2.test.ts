import { describe, it, expect } from 'vitest';
import { Image2 } from '../src/Image2.js';
import { ImageUtility2 } from '../src/ImageUtility2.js';

// Build an Image2<number> of the given size whose pixels are 1 exactly at the
// listed (x,y) coordinates.
function makeImage(dim0: number, dim1: number, ones: readonly (readonly [number, number])[]): Image2<number> {
    const image = new Image2<number>(dim0, dim1);
    for (const [x, y] of ones) {
        image.set(x, y, 1);
    }
    return image;
}

// The set of (x,y) with a nonzero pixel, as sorted "x,y" strings.
function nonzeroCoords(image: Image2<number>): string[] {
    const result: string[] = [];
    for (let y = 0; y < image.getDimension(1); ++y) {
        for (let x = 0; x < image.getDimension(0); ++x) {
            if (image.get(x, y) !== 0) {
                result.push(`${x},${y}`);
            }
        }
    }
    return result.sort();
}

function coordKeys(pairs: readonly (readonly [number, number])[]): string[] {
    return pairs.map(([x, y]) => `${x},${y}`).sort();
}

// A simple xorshift generator so the randomized tests are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state ^= state << 13; state >>>= 0;
        state ^= state >>> 17;
        state ^= state << 5; state >>>= 0;
        return state / 4294967296;
    };
}

// Reference connected-component labeling by breadth-first search over the
// 2D neighborhood, used to cross-check the ported depth-first algorithm.
function referenceComponents(pixels: readonly number[], dim0: number, dim1: number,
    offsets: readonly (readonly [number, number])[]): number[][] {
    const label = new Array<number>(dim0 * dim1).fill(0);
    const components: number[][] = [[]];
    for (let i = 0; i < pixels.length; ++i) {
        if (pixels[i] !== 1 || label[i] !== 0) {
            continue;
        }
        const current = components.length;
        const member: number[] = [];
        const queue = [i];
        label[i] = current;
        while (queue.length > 0) {
            const v = queue.pop() as number;
            member.push(v);
            const x = v % dim0;
            const y = Math.floor(v / dim0);
            for (const [dx, dy] of offsets) {
                const nx = x + dx;
                const ny = y + dy;
                if (0 <= nx && nx < dim0 && 0 <= ny && ny < dim1) {
                    const w = nx + dim0 * ny;
                    if (pixels[w] === 1 && label[w] === 0) {
                        label[w] = current;
                        queue.push(w);
                    }
                }
            }
        }
        member.sort((a, b) => a - b);
        components.push(member);
    }
    return components;
}

describe('ImageUtility2.getComponents', () => {
    it('separates diagonal pixels with 4-connectivity', () => {
        // (1,1) and (2,2) touch only diagonally; (5,5) and (5,6) are a
        // separate 4-connected pair.
        const image = makeImage(8, 8, [[1, 1], [2, 2], [5, 5], [5, 6]]);
        const components = ImageUtility2.getComponents4(image);
        expect(components.length).toBe(4);
        expect(components[0]).toEqual([]);
        expect(components[1]).toEqual([1 + 8 * 1]);
        expect(components[2]).toEqual([2 + 8 * 2]);
        expect(components[3]).toEqual([5 + 8 * 5, 5 + 8 * 6]);

        // The image now holds the component labels.
        expect(image.get(1, 1)).toBe(1);
        expect(image.get(2, 2)).toBe(2);
        expect(image.get(5, 5)).toBe(3);
        expect(image.get(5, 6)).toBe(3);
        expect(image.get(4, 4)).toBe(0);
    });

    it('merges diagonal pixels with 8-connectivity', () => {
        const image = makeImage(8, 8, [[1, 1], [2, 2], [5, 5], [5, 6]]);
        const components = ImageUtility2.getComponents8(image);
        expect(components.length).toBe(3);
        expect(components[1]).toEqual([1 + 8 * 1, 2 + 8 * 2]);
        expect(components[2]).toEqual([5 + 8 * 5, 5 + 8 * 6]);
        expect(image.get(1, 1)).toBe(1);
        expect(image.get(2, 2)).toBe(1);
    });

    it('returns an empty array for an all-background image', () => {
        const image = new Image2<number>(6, 6);
        expect(ImageUtility2.getComponents4(image)).toEqual([]);
    });

    it('labels a single component of every foreground pixel', () => {
        // A solid 4x4 block inside a zero boundary is one component.
        const ones: [number, number][] = [];
        for (let y = 1; y <= 4; ++y) {
            for (let x = 1; x <= 4; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(6, 6, ones);
        const components = ImageUtility2.getComponents4(image);
        expect(components.length).toBe(2);
        expect(components[1].length).toBe(16);
        // The indices are produced in increasing order by the final pass.
        const sorted = components[1].slice().sort((a, b) => a - b);
        expect(components[1]).toEqual(sorted);
    });

    it('accepts caller-specified 1-dimensional neighbor offsets', () => {
        // A horizontal-only structuring element connects (1,1)-(2,1) but not
        // (1,1)-(1,2).
        const image = makeImage(6, 6, [[1, 1], [2, 1], [1, 2]]);
        const components = ImageUtility2.getComponents(image, [-1, +1]);
        expect(components.length).toBe(3);
        expect(components[1]).toEqual([1 + 6 * 1, 2 + 6 * 1]);
        expect(components[2]).toEqual([1 + 6 * 2]);
    });

    it('agrees with a breadth-first reference on random images', () => {
        const random = makeRandom(20250830);
        const dim0 = 14;
        const dim1 = 11;
        for (let trial = 0; trial < 20; ++trial) {
            const original = new Array<number>(dim0 * dim1).fill(0);
            for (let y = 1; y < dim1 - 1; ++y) {
                for (let x = 1; x < dim0 - 1; ++x) {
                    original[x + dim0 * y] = (random() < 0.4 ? 1 : 0);
                }
            }

            for (const connectivity of [4, 8] as const) {
                const image = new Image2<number>(dim0, dim1);
                for (let i = 0; i < original.length; ++i) {
                    image.set(i, original[i]);
                }
                const offsets = (connectivity === 4
                    ? image.getNeighborhood4Coords()
                    : image.getNeighborhood8Coords()) as [number, number][];
                const expected = referenceComponents(original, dim0, dim1, offsets);
                const actual = connectivity === 4
                    ? ImageUtility2.getComponents4(image)
                    : ImageUtility2.getComponents8(image);

                // Both algorithms visit pixels in increasing index order, so
                // the components are labeled in the same order.
                if (expected.length === 1) {
                    expect(actual).toEqual([]);
                } else {
                    expect(actual.length).toBe(expected.length);
                    for (let k = 1; k < expected.length; ++k) {
                        expect(actual[k]).toEqual(expected[k]);
                    }
                    // The image labels agree with the component lists.
                    for (let k = 1; k < expected.length; ++k) {
                        for (const index of expected[k]) {
                            expect(image.get(index)).toBe(k);
                        }
                    }
                }
            }
        }
    });
});

describe('ImageUtility2.dilate and erode', () => {
    it('dilates a single pixel to the 4- and 8-neighborhoods', () => {
        const image = makeImage(7, 7, [[3, 3]]);
        expect(nonzeroCoords(ImageUtility2.dilate4(image))).toEqual(
            coordKeys([[3, 3], [2, 3], [4, 3], [3, 2], [3, 4]]));

        const block: [number, number][] = [];
        for (let y = 2; y <= 4; ++y) {
            for (let x = 2; x <= 4; ++x) {
                block.push([x, y]);
            }
        }
        expect(nonzeroCoords(ImageUtility2.dilate8(image))).toEqual(coordKeys(block));

        // The input image is not modified.
        expect(nonzeroCoords(image)).toEqual(['3,3']);
    });

    it('clips dilation at the image boundary', () => {
        const image = makeImage(3, 3, [[0, 0]]);
        expect(nonzeroCoords(ImageUtility2.dilate4(image))).toEqual(
            coordKeys([[0, 0], [1, 0], [0, 1]]));
    });

    it('erodes a solid block to its center', () => {
        const ones: [number, number][] = [];
        for (let y = 2; y <= 4; ++y) {
            for (let x = 2; x <= 4; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(7, 7, ones);
        expect(nonzeroCoords(ImageUtility2.erode4(image, true))).toEqual(['3,3']);
        expect(nonzeroCoords(ImageUtility2.erode8(image, true))).toEqual(['3,3']);
    });

    it('honors zeroExterior when the foreground touches the border', () => {
        const ones: [number, number][] = [];
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(5, 5, ones);

        // With a zero exterior, the border pixels are eroded away.
        const interior: [number, number][] = [];
        for (let y = 1; y <= 3; ++y) {
            for (let x = 1; x <= 3; ++x) {
                interior.push([x, y]);
            }
        }
        expect(nonzeroCoords(ImageUtility2.erode4(image, true))).toEqual(coordKeys(interior));

        // Without it, no in-range neighbor is background so nothing erodes.
        expect(nonzeroCoords(ImageUtility2.erode4(image, false))).toEqual(coordKeys(ones));
    });

    it('rejects an empty structuring element', () => {
        const image = makeImage(4, 4, [[1, 1]]);
        expect(() => ImageUtility2.dilate(image, [])).toThrow('Invalid neighbors.');
        expect(() => ImageUtility2.erode(image, true, [])).toThrow('Invalid neighbors.');
    });

    it('dilation and erosion are dual under complementation', () => {
        // erode(A) == complement(dilate(complement(A))) for a symmetric
        // structuring element, when the exterior is treated consistently.
        const random = makeRandom(777);
        const dim0 = 9;
        const dim1 = 8;
        const image = new Image2<number>(dim0, dim1);
        const complement = new Image2<number>(dim0, dim1);
        for (let i = 0; i < dim0 * dim1; ++i) {
            const value = (random() < 0.5 ? 1 : 0);
            image.set(i, value);
            complement.set(i, 1 - value);
        }

        const eroded = ImageUtility2.erode8(image, false);
        const dilatedComplement = ImageUtility2.dilate8(complement);
        for (let i = 0; i < dim0 * dim1; ++i) {
            expect(eroded.get(i)).toBe(1 - dilatedComplement.get(i));
        }
    });
});

describe('ImageUtility2.open and close', () => {
    it('opening removes an isolated pixel and rounds a block', () => {
        const ones: [number, number][] = [[8, 8]];
        for (let y = 2; y <= 4; ++y) {
            for (let x = 2; x <= 4; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(11, 11, ones);
        // Erosion leaves only the block center, and dilation grows it back
        // into a plus.
        expect(nonzeroCoords(ImageUtility2.open4(image, true))).toEqual(
            coordKeys([[3, 3], [2, 3], [4, 3], [3, 2], [3, 4]]));
        // With the 8-neighborhood, the block is restored exactly and the
        // isolated pixel is still removed.
        const block = ones.slice(1);
        expect(nonzeroCoords(ImageUtility2.open8(image, true))).toEqual(coordKeys(block));
    });

    it('closing fills a one-pixel hole', () => {
        const ones: [number, number][] = [];
        for (let y = 1; y <= 5; ++y) {
            for (let x = 1; x <= 5; ++x) {
                if (x !== 3 || y !== 3) {
                    ones.push([x, y]);
                }
            }
        }
        const image = makeImage(7, 7, ones);
        const closed = ImageUtility2.close4(image, true);
        const filled = ones.concat([[3, 3]]);
        expect(nonzeroCoords(closed)).toEqual(coordKeys(filled));
    });

    it('caller-specified structuring elements match the 4-neighborhood forms', () => {
        const image = makeImage(9, 9, [[3, 3], [4, 3], [3, 4], [6, 6]]);
        const nbrs = image.getNeighborhood4Coords();
        expect(nonzeroCoords(ImageUtility2.open(image, true, nbrs)))
            .toEqual(nonzeroCoords(ImageUtility2.open4(image, true)));
        expect(nonzeroCoords(ImageUtility2.close(image, true, nbrs)))
            .toEqual(nonzeroCoords(ImageUtility2.close4(image, true)));
    });
});

describe('ImageUtility2.extractBoundary', () => {
    it('walks the ring of a solid 3x3 block', () => {
        const ones: [number, number][] = [];
        for (let y = 2; y <= 4; ++y) {
            for (let x = 2; x <= 4; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(7, 7, ones);
        const { success, boundary } = ImageUtility2.extractBoundary(0, 0, image);
        expect(success).toBe(true);
        // Clockwise traversal starting at the first nonzero pixel (2,2).
        expect(boundary).toEqual([16, 17, 18, 25, 32, 31, 30, 23]);

        // Visited boundary pixels are marked 2 and the center is untouched.
        for (const index of boundary) {
            expect(image.get(index)).toBe(2);
        }
        expect(image.get(3, 3)).toBe(1);
    });

    it('reports failure when there is no foreground pixel', () => {
        const image = new Image2<number>(5, 5);
        const { success, boundary } = ImageUtility2.extractBoundary(0, 0, image);
        expect(success).toBe(false);
        expect(boundary).toEqual([]);
    });

    it('stops immediately at an isolated pixel', () => {
        const image = makeImage(5, 5, [[2, 2]]);
        const { success, boundary } = ImageUtility2.extractBoundary(0, 0, image);
        expect(success).toBe(true);
        expect(boundary).toEqual([2 + 5 * 2]);
        expect(image.get(2, 2)).toBe(2);
    });
});

describe('ImageUtility2.floodFill4 and drawFloodFill4', () => {
    it('fills only the 4-connected region on one side of a wall', () => {
        const image = new Image2<number>(5, 5);
        for (let y = 0; y < 5; ++y) {
            image.set(2, y, 5);
        }
        ImageUtility2.floodFill4(image, 0, 0, 1, 0);
        for (let y = 0; y < 5; ++y) {
            expect(image.get(0, y)).toBe(1);
            expect(image.get(1, y)).toBe(1);
            expect(image.get(2, y)).toBe(5);
            expect(image.get(3, y)).toBe(0);
            expect(image.get(4, y)).toBe(0);
        }
    });

    it('does nothing for a seed outside the image', () => {
        const image = new Image2<number>(4, 4);
        ImageUtility2.floodFill4(image, -1, 0, 1, 0);
        ImageUtility2.floodFill4(image, 0, 4, 1, 0);
        expect(image.getPixels().every((v) => v === 0)).toBe(true);
    });

    it('drawFloodFill4 matches floodFill4 through callbacks', () => {
        const image = new Image2<number>(6, 6);
        const mirror = new Image2<number>(6, 6);
        for (let y = 0; y < 6; ++y) {
            image.set(3, y, 7);
            mirror.set(3, y, 7);
        }
        ImageUtility2.floodFill4(image, 5, 5, 1, 0);
        ImageUtility2.drawFloodFill4<number>(5, 5, 6, 6, 1, 0,
            (x, y, value) => mirror.set(x, y, value),
            (x, y) => mirror.get(x, y));
        expect(mirror.getPixels()).toEqual(image.getPixels());
        // Only the region to the right of the wall is filled.
        expect(image.get(4, 0)).toBe(1);
        expect(image.get(2, 0)).toBe(0);
    });
});

describe('ImageUtility2.getL1Distance', () => {
    it('computes the city-block distance to the background', () => {
        const dim = 9;
        const ones: [number, number][] = [];
        for (let y = 1; y <= dim - 2; ++y) {
            for (let x = 1; x <= dim - 2; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(dim, dim, ones);
        const { maxDistance, xMax, yMax } = ImageUtility2.getL1Distance(image);
        expect(maxDistance).toBe(4);
        expect(xMax).toBe(4);
        expect(yMax).toBe(4);
        for (let y = 1; y <= dim - 2; ++y) {
            for (let x = 1; x <= dim - 2; ++x) {
                expect(image.get(x, y)).toBe(Math.min(x, y, dim - 1 - x, dim - 1 - y));
            }
        }
    });

    it('leaves a one-pixel-thick object unchanged', () => {
        const image = makeImage(9, 9, [[2, 4], [3, 4], [4, 4], [5, 4], [6, 4]]);
        const { maxDistance } = ImageUtility2.getL1Distance(image);
        expect(maxDistance).toBe(1);
        for (let x = 2; x <= 6; ++x) {
            expect(image.get(x, 4)).toBe(1);
        }
    });
});

describe('ImageUtility2.getL2Distance', () => {
    it('matches the exact Euclidean distance to the nearest background pixel', () => {
        const random = makeRandom(31337);
        const dim0 = 13;
        const dim1 = 10;
        for (let trial = 0; trial < 8; ++trial) {
            const image = new Image2<number>(dim0, dim1);
            const background: [number, number][] = [];
            for (let y = 0; y < dim1; ++y) {
                for (let x = 0; x < dim0; ++x) {
                    if (random() < 0.75) {
                        image.set(x, y, 1);
                    } else {
                        background.push([x, y]);
                    }
                }
            }
            if (background.length === 0) {
                continue;
            }

            const { maxDistance, transform } = ImageUtility2.getL2Distance(image);

            let expectedMax = 0;
            for (let y = 0; y < dim1; ++y) {
                for (let x = 0; x < dim0; ++x) {
                    let best = Number.MAX_VALUE;
                    for (const [u, v] of background) {
                        const sq = (x - u) * (x - u) + (y - v) * (y - v);
                        if (sq < best) {
                            best = sq;
                        }
                    }
                    const expected = Math.fround(Math.sqrt(best));
                    expect(transform.get(x, y)).toBe(expected);
                    if (expected > expectedMax) {
                        expectedMax = expected;
                    }
                }
            }
            expect(maxDistance).toBe(expectedMax);
        }
    });

    it('is zero everywhere for an all-background image', () => {
        const image = new Image2<number>(5, 5);
        const { maxDistance, xMax, yMax, transform } = ImageUtility2.getL2Distance(image);
        expect(maxDistance).toBe(0);
        expect(xMax).toBe(0);
        expect(yMax).toBe(0);
        expect(transform.getPixels().every((v) => v === 0)).toBe(true);
    });
});

describe('ImageUtility2.getSkeleton', () => {
    it('leaves a one-pixel-thick segment unchanged', () => {
        const ones: [number, number][] = [[2, 4], [3, 4], [4, 4], [5, 4], [6, 4]];
        const image = makeImage(9, 9, ones);
        ImageUtility2.getSkeleton(image);
        expect(nonzeroCoords(image)).toEqual(coordKeys(ones));
    });

    it('thins a solid block to a binary subset of the original', () => {
        const ones: [number, number][] = [];
        for (let y = 2; y <= 6; ++y) {
            for (let x = 2; x <= 8; ++x) {
                ones.push([x, y]);
            }
        }
        const image = makeImage(11, 9, ones);
        const before = nonzeroCoords(image);
        ImageUtility2.getSkeleton(image);
        const after = nonzeroCoords(image);

        // The skeleton is binary, nonempty, and a proper subset.
        expect(image.getPixels().every((v) => v === 0 || v === 1)).toBe(true);
        expect(after.length).toBeGreaterThan(0);
        expect(after.length).toBeLessThan(before.length);
        for (const key of after) {
            expect(before).toContain(key);
        }

        // The skeleton of a solid convex block stays 8-connected.
        const components = ImageUtility2.getComponents8(image);
        expect(components.length).toBe(2);
    });

    it('leaves an empty image empty', () => {
        const image = new Image2<number>(6, 6);
        ImageUtility2.getSkeleton(image);
        expect(image.getPixels().every((v) => v === 0)).toBe(true);
    });
});

describe('ImageUtility2 drawing primitives', () => {
    it('drawThickPixel visits a square of the requested half-width', () => {
        const visited: [number, number][] = [];
        ImageUtility2.drawThickPixel(5, 7, 1, (x, y) => visited.push([x, y]));
        expect(visited.length).toBe(9);
        expect(coordKeys(visited)).toEqual(coordKeys([
            [4, 6], [5, 6], [6, 6],
            [4, 7], [5, 7], [6, 7],
            [4, 8], [5, 8], [6, 8]
        ]));

        const single: [number, number][] = [];
        ImageUtility2.drawThickPixel(0, 0, 0, (x, y) => single.push([x, y]));
        expect(single).toEqual([[0, 0]]);
    });

    it('drawLine reproduces the Bresenham pixels', () => {
        const visited: [number, number][] = [];
        ImageUtility2.drawLine(0, 0, 4, 2, (x, y) => visited.push([x, y]));
        expect(visited).toEqual([[0, 0], [1, 1], [2, 1], [3, 2], [4, 2]]);
    });

    it('drawLine handles steep, axis-aligned and degenerate lines', () => {
        const steep: [number, number][] = [];
        ImageUtility2.drawLine(0, 0, 2, 4, (x, y) => steep.push([x, y]));
        expect(steep).toEqual([[0, 0], [1, 1], [1, 2], [2, 3], [2, 4]]);

        const horizontal: [number, number][] = [];
        ImageUtility2.drawLine(2, 3, 5, 3, (x, y) => horizontal.push([x, y]));
        expect(horizontal).toEqual([[2, 3], [3, 3], [4, 3], [5, 3]]);

        const vertical: [number, number][] = [];
        ImageUtility2.drawLine(2, 3, 2, 0, (x, y) => vertical.push([x, y]));
        expect(vertical).toEqual([[2, 3], [2, 2], [2, 1], [2, 0]]);

        const point: [number, number][] = [];
        ImageUtility2.drawLine(4, 4, 4, 4, (x, y) => point.push([x, y]));
        expect(point).toEqual([[4, 4]]);
    });

    it('drawLine visits max(|dx|,|dy|)+1 pixels and ends at the endpoint', () => {
        const random = makeRandom(9182736);
        for (let trial = 0; trial < 40; ++trial) {
            const x0 = Math.floor(random() * 21) - 10;
            const y0 = Math.floor(random() * 21) - 10;
            const x1 = Math.floor(random() * 21) - 10;
            const y1 = Math.floor(random() * 21) - 10;
            const visited: [number, number][] = [];
            ImageUtility2.drawLine(x0, y0, x1, y1, (x, y) => visited.push([x, y]));
            const expectedCount = Math.max(Math.abs(x1 - x0), Math.abs(y1 - y0)) + 1;
            expect(visited.length).toBe(expectedCount);
            expect(visited[0]).toEqual([x0, y0]);
            expect(visited[visited.length - 1]).toEqual([x1, y1]);
            // Consecutive pixels are 8-connected.
            for (let i = 1; i < visited.length; ++i) {
                expect(Math.abs(visited[i][0] - visited[i - 1][0])).toBeLessThanOrEqual(1);
                expect(Math.abs(visited[i][1] - visited[i - 1][1])).toBeLessThanOrEqual(1);
            }
        }
    });

    it('drawCircle visits the Bresenham circle', () => {
        const visited = new Set<string>();
        ImageUtility2.drawCircle(0, 0, 3, false, (x, y) => visited.add(`${x},${y}`));
        const expected = new Set<string>();
        for (const [a, b] of [[0, 3], [1, 3], [2, 2], [3, 1], [3, 0]] as [number, number][]) {
            for (const [sx, sy] of [[1, 1], [1, -1], [-1, 1], [-1, -1]] as [number, number][]) {
                expected.add(`${sx * a},${sy * b}`);
            }
        }
        expect([...visited].sort()).toEqual([...expected].sort());
    });

    it('drawCircle solid covers the disk bounded by the outline', () => {
        const outline = new Set<string>();
        ImageUtility2.drawCircle(10, 10, 5, false, (x, y) => outline.add(`${x},${y}`));
        const solid = new Set<string>();
        ImageUtility2.drawCircle(10, 10, 5, true, (x, y) => solid.add(`${x},${y}`));

        // Every outline pixel is in the solid disk.
        for (const key of outline) {
            expect(solid.has(key)).toBe(true);
        }
        // The solid disk is a set of full rows spanning the outline.
        for (const key of solid) {
            const [x, y] = key.split(',').map(Number);
            expect((x - 10) * (x - 10) + (y - 10) * (y - 10)).toBeLessThanOrEqual(5 * 5 + 5);
        }
        expect(solid.has('10,10')).toBe(true);
        expect(solid.size).toBeGreaterThan(outline.size);
    });

    it('drawCircle with radius zero visits the center only', () => {
        const visited = new Set<string>();
        ImageUtility2.drawCircle(2, 3, 0, false, (x, y) => visited.add(`${x},${y}`));
        expect([...visited]).toEqual(['2,3']);
    });

    it('drawRectangle visits the outline or the filled interior', () => {
        const outline: [number, number][] = [];
        ImageUtility2.drawRectangle(1, 2, 3, 4, false, (x, y) => outline.push([x, y]));
        expect(coordKeys(outline)).toEqual(coordKeys([
            [1, 2], [2, 2], [3, 2],
            [1, 4], [2, 4], [3, 4],
            [1, 3], [3, 3]
        ]));

        const solid: [number, number][] = [];
        ImageUtility2.drawRectangle(1, 2, 3, 4, true, (x, y) => solid.push([x, y]));
        expect(solid.length).toBe(9);
        expect(coordKeys(solid)).toContain('2,3');
    });

    it('drawEllipse is symmetric about its center and hits the vertices', () => {
        const visited = new Set<string>();
        ImageUtility2.drawEllipse(0, 0, 6, 3, (x, y) => visited.add(`${x},${y}`));
        for (const key of visited) {
            const [x, y] = key.split(',').map(Number);
            expect(visited.has(`${-x},${y}`)).toBe(true);
            expect(visited.has(`${x},${-y}`)).toBe(true);
        }
        expect(visited.has('6,0')).toBe(true);
        expect(visited.has('-6,0')).toBe(true);
        expect(visited.has('0,3')).toBe(true);
        expect(visited.has('0,-3')).toBe(true);
        // Every visited pixel is near the ellipse.
        for (const key of visited) {
            const [x, y] = key.split(',').map(Number);
            const value = (x / 6) * (x / 6) + (y / 3) * (y / 3);
            expect(value).toBeGreaterThan(0.3);
            expect(value).toBeLessThan(2.0);
        }
    });

    it('drawEllipse degenerates to the vertices when an extent is zero', () => {
        const visited = new Set<string>();
        ImageUtility2.drawEllipse(0, 0, 4, 0, (x, y) => visited.add(`${x},${y}`));
        expect([...visited].sort()).toEqual(['-4,0', '0,0', '4,0']);
    });
});
