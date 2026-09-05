import { describe, it, expect } from 'vitest';
import { Image3 } from '../src/Image3.js';
import { ImageUtility3 } from '../src/ImageUtility3.js';
import { check, fc } from './helpers/arbitraries.js';

type Voxel = readonly [number, number, number];

// Build an Image3<number> of the given size whose voxels are 1 exactly at the
// listed (x,y,z) coordinates.
function makeImage(dim0: number, dim1: number, dim2: number,
    ones: readonly Voxel[]): Image3<number> {
    const image = new Image3<number>(dim0, dim1, dim2);
    for (const [x, y, z] of ones) {
        image.set(x, y, z, 1);
    }
    return image;
}

// The set of (x,y,z) with a nonzero voxel, as sorted "x,y,z" strings.
function nonzeroCoords(image: Image3<number>): string[] {
    const result: string[] = [];
    for (let z = 0; z < image.getDimension(2); ++z) {
        for (let y = 0; y < image.getDimension(1); ++y) {
            for (let x = 0; x < image.getDimension(0); ++x) {
                if (image.get(x, y, z) !== 0) {
                    result.push(`${x},${y},${z}`);
                }
            }
        }
    }
    return result.sort();
}

function coordKeys(voxels: readonly Voxel[]): string[] {
    return voxels.map(([x, y, z]) => `${x},${y},${z}`).sort();
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state ^= state << 13; state >>>= 0;
        state ^= state >>> 17;
        state ^= state << 5; state >>>= 0;
        return state / 4294967296;
    };
}

// Reference connected-component labeling over the 3D neighborhood.
function referenceComponents(voxels: readonly number[], dim0: number, dim1: number, dim2: number,
    offsets: readonly (readonly number[])[]): number[][] {
    const label = new Array<number>(dim0 * dim1 * dim2).fill(0);
    const components: number[][] = [[]];
    for (let i = 0; i < voxels.length; ++i) {
        if (voxels[i] !== 1 || label[i] !== 0) {
            continue;
        }
        const current = components.length;
        const member: number[] = [];
        const stack = [i];
        label[i] = current;
        while (stack.length > 0) {
            const v = stack.pop() as number;
            member.push(v);
            const x = v % dim0;
            const y = Math.floor(v / dim0) % dim1;
            const z = Math.floor(v / (dim0 * dim1));
            for (const offset of offsets) {
                const nx = x + offset[0];
                const ny = y + offset[1];
                const nz = z + offset[2];
                if (0 <= nx && nx < dim0 && 0 <= ny && ny < dim1 && 0 <= nz && nz < dim2) {
                    const w = nx + dim0 * (ny + dim1 * nz);
                    if (voxels[w] === 1 && label[w] === 0) {
                        label[w] = current;
                        stack.push(w);
                    }
                }
            }
        }
        member.sort((a, b) => a - b);
        components.push(member);
    }
    return components;
}

// The four probe voxels used by the connectivity tests, in a 6x6x6 image:
// A and B differ by a face diagonal (18-connected but not 6-connected),
// C and D differ by a body diagonal (26-connected only).
const A: Voxel = [1, 1, 1];
const B: Voxel = [2, 2, 1];
const C: Voxel = [1, 1, 3];
const D: Voxel = [2, 2, 4];
const indexA = 1 + 6 * (1 + 6 * 1);
const indexB = 2 + 6 * (2 + 6 * 1);
const indexC = 1 + 6 * (1 + 6 * 3);
const indexD = 2 + 6 * (2 + 6 * 4);

describe('ImageUtility3.getComponents', () => {
    it('separates face- and body-diagonal voxels with 6-connectivity', () => {
        const image = makeImage(6, 6, 6, [A, B, C, D]);
        const components = ImageUtility3.getComponents6(image);
        expect(components.length).toBe(5);
        expect(components[0]).toEqual([]);
        expect(components[1]).toEqual([indexA]);
        expect(components[2]).toEqual([indexB]);
        expect(components[3]).toEqual([indexC]);
        expect(components[4]).toEqual([indexD]);
        expect(image.get(1, 1, 1)).toBe(1);
        expect(image.get(2, 2, 4)).toBe(4);
    });

    it('merges face diagonals but not body diagonals with 18-connectivity', () => {
        const image = makeImage(6, 6, 6, [A, B, C, D]);
        const components = ImageUtility3.getComponents18(image);
        expect(components.length).toBe(4);
        expect(components[1]).toEqual([indexA, indexB]);
        expect(components[2]).toEqual([indexC]);
        expect(components[3]).toEqual([indexD]);
    });

    it('merges body diagonals with 26-connectivity', () => {
        const image = makeImage(6, 6, 6, [A, B, C, D]);
        const components = ImageUtility3.getComponents26(image);
        expect(components.length).toBe(3);
        expect(components[1]).toEqual([indexA, indexB]);
        expect(components[2]).toEqual([indexC, indexD]);
    });

    it('returns an empty array for an all-background image', () => {
        const image = new Image3<number>(5, 5, 5);
        expect(ImageUtility3.getComponents6(image)).toEqual([]);
    });

    it('rejects an empty neighborhood', () => {
        const image = makeImage(5, 5, 5, [[2, 2, 2]]);
        expect(() => ImageUtility3.getComponents(image, [])).toThrow('Invalid neighbors.');
    });

    it('labels a solid block as a single component', () => {
        const ones: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    ones.push([x, y, z]);
                }
            }
        }
        const image = makeImage(5, 5, 5, ones);
        const components = ImageUtility3.getComponents6(image);
        expect(components.length).toBe(2);
        expect(components[1].length).toBe(27);
    });

    it('agrees with a reference labeling on random images', () => {
        const random = makeRandom(4242);
        const dim = 7;
        for (let trial = 0; trial < 10; ++trial) {
            const original = new Array<number>(dim * dim * dim).fill(0);
            for (let z = 1; z < dim - 1; ++z) {
                for (let y = 1; y < dim - 1; ++y) {
                    for (let x = 1; x < dim - 1; ++x) {
                        original[x + dim * (y + dim * z)] = (random() < 0.3 ? 1 : 0);
                    }
                }
            }

            for (const connectivity of [6, 18, 26] as const) {
                const image = new Image3<number>(dim, dim, dim);
                for (let i = 0; i < original.length; ++i) {
                    image.set(i, original[i]);
                }
                const offsets = connectivity === 6 ? image.getNeighborhood6Coords()
                    : connectivity === 18 ? image.getNeighborhood18Coords()
                        : image.getNeighborhood26Coords();
                const expected = referenceComponents(original, dim, dim, dim, offsets);
                const actual = connectivity === 6 ? ImageUtility3.getComponents6(image)
                    : connectivity === 18 ? ImageUtility3.getComponents18(image)
                        : ImageUtility3.getComponents26(image);

                if (expected.length === 1) {
                    expect(actual).toEqual([]);
                } else {
                    expect(actual.length).toBe(expected.length);
                    for (let k = 1; k < expected.length; ++k) {
                        expect(actual[k]).toEqual(expected[k]);
                        for (const index of expected[k]) {
                            expect(image.get(index)).toBe(k);
                        }
                    }
                }
            }
        }
    });
});

describe('ImageUtility3.dilate and erode', () => {
    it('dilates a single voxel to the 6-, 18- and 26-neighborhoods', () => {
        const image = makeImage(6, 6, 6, [[2, 2, 2]]);
        expect(nonzeroCoords(ImageUtility3.dilate6(image)).length).toBe(7);
        expect(nonzeroCoords(ImageUtility3.dilate18(image)).length).toBe(19);
        expect(nonzeroCoords(ImageUtility3.dilate26(image)).length).toBe(27);

        expect(nonzeroCoords(ImageUtility3.dilate6(image))).toEqual(coordKeys([
            [2, 2, 2], [1, 2, 2], [3, 2, 2], [2, 1, 2], [2, 3, 2], [2, 2, 1], [2, 2, 3]
        ]));

        const block: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    block.push([x, y, z]);
                }
            }
        }
        expect(nonzeroCoords(ImageUtility3.dilate26(image))).toEqual(coordKeys(block));

        // The input image is not modified.
        expect(nonzeroCoords(image)).toEqual(['2,2,2']);
    });

    it('dilates a voxel in the x = 0 plane (upstream skips that column)', () => {
        // Upstream ImageUtility3::Dilate starts its innermost loop at i0 = 1,
        // so a foreground voxel at x = 0 would not dilate at all. The port
        // starts at i0 = 0 to match the 2D version and the documentation.
        const image = makeImage(4, 4, 4, [[0, 1, 1]]);
        expect(nonzeroCoords(ImageUtility3.dilate6(image))).toEqual(coordKeys([
            [0, 1, 1], [1, 1, 1], [0, 0, 1], [0, 2, 1], [0, 1, 0], [0, 1, 2]
        ]));
    });

    it('erodes a solid block to its center', () => {
        const ones: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    ones.push([x, y, z]);
                }
            }
        }
        const image = makeImage(5, 5, 5, ones);
        expect(nonzeroCoords(ImageUtility3.erode6(image, true))).toEqual(['2,2,2']);
        expect(nonzeroCoords(ImageUtility3.erode26(image, true))).toEqual(['2,2,2']);
    });

    it('honors zeroExterior when the foreground touches the border', () => {
        const ones: Voxel[] = [];
        for (let z = 0; z < 4; ++z) {
            for (let y = 0; y < 4; ++y) {
                for (let x = 0; x < 4; ++x) {
                    ones.push([x, y, z]);
                }
            }
        }
        const image = makeImage(4, 4, 4, ones);

        const interior: Voxel[] = [];
        for (let z = 1; z <= 2; ++z) {
            for (let y = 1; y <= 2; ++y) {
                for (let x = 1; x <= 2; ++x) {
                    interior.push([x, y, z]);
                }
            }
        }
        expect(nonzeroCoords(ImageUtility3.erode6(image, true))).toEqual(coordKeys(interior));
        expect(nonzeroCoords(ImageUtility3.erode6(image, false))).toEqual(coordKeys(ones));
    });

    it('rejects an empty structuring element', () => {
        const image = makeImage(4, 4, 4, [[1, 1, 1]]);
        expect(() => ImageUtility3.dilate(image, [])).toThrow('Invalid neighbors.');
        expect(() => ImageUtility3.erode(image, true, [])).toThrow('Invalid neighbors.');
    });

    it('dilation and erosion are dual under complementation', () => {
        const random = makeRandom(2468);
        const dim = 5;
        const image = new Image3<number>(dim, dim, dim);
        const complement = new Image3<number>(dim, dim, dim);
        for (let i = 0; i < dim * dim * dim; ++i) {
            const value = (random() < 0.5 ? 1 : 0);
            image.set(i, value);
            complement.set(i, 1 - value);
        }

        const eroded = ImageUtility3.erode26(image, false);
        const dilatedComplement = ImageUtility3.dilate26(complement);
        for (let i = 0; i < dim * dim * dim; ++i) {
            expect(eroded.get(i)).toBe(1 - dilatedComplement.get(i));
        }
    });
});

describe('ImageUtility3.open and close', () => {
    it('opening removes an isolated voxel and restores a solid block', () => {
        const ones: Voxel[] = [[5, 5, 5]];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    ones.push([x, y, z]);
                }
            }
        }
        const image = makeImage(8, 8, 8, ones);
        expect(nonzeroCoords(ImageUtility3.open26(image, true))).toEqual(
            coordKeys(ones.slice(1)));

        // With the 6-neighborhood the eroded center grows back into an
        // octahedron, not the full block.
        expect(nonzeroCoords(ImageUtility3.open6(image, true))).toEqual(coordKeys([
            [2, 2, 2], [1, 2, 2], [3, 2, 2], [2, 1, 2], [2, 3, 2], [2, 2, 1], [2, 2, 3]
        ]));
    });

    it('closing fills a one-voxel hole', () => {
        const ones: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    if (x !== 2 || y !== 2 || z !== 2) {
                        ones.push([x, y, z]);
                    }
                }
            }
        }
        const image = makeImage(5, 5, 5, ones);
        const closed = ImageUtility3.close6(image, true);
        expect(closed.get(2, 2, 2)).toBe(1);
        expect(nonzeroCoords(closed)).toEqual(coordKeys(ones.concat([[2, 2, 2]])));
    });

    it('caller-specified structuring elements match the 6-neighborhood forms', () => {
        const image = makeImage(7, 7, 7, [[2, 2, 2], [3, 2, 2], [2, 3, 2], [5, 5, 5]]);
        const nbrs = image.getNeighborhood6Coords();
        expect(nonzeroCoords(ImageUtility3.open(image, true, nbrs)))
            .toEqual(nonzeroCoords(ImageUtility3.open6(image, true)));
        expect(nonzeroCoords(ImageUtility3.close(image, true, nbrs)))
            .toEqual(nonzeroCoords(ImageUtility3.close6(image, true)));
    });
});

describe('ImageUtility3.computeCDConvex', () => {
    it('fills the interior hole of a solid block', () => {
        const ones: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    if (x !== 2 || y !== 2 || z !== 2) {
                        ones.push([x, y, z]);
                    }
                }
            }
        }
        const image = makeImage(5, 5, 5, ones);
        ImageUtility3.computeCDConvex(image);
        expect(nonzeroCoords(image)).toEqual(coordKeys(ones.concat([[2, 2, 2]])));
    });

    it('is idempotent on a solid block', () => {
        const ones: Voxel[] = [];
        for (let z = 1; z <= 3; ++z) {
            for (let y = 1; y <= 3; ++y) {
                for (let x = 1; x <= 3; ++x) {
                    ones.push([x, y, z]);
                }
            }
        }
        const image = makeImage(5, 5, 5, ones);
        ImageUtility3.computeCDConvex(image);
        expect(nonzeroCoords(image)).toEqual(coordKeys(ones));
        ImageUtility3.computeCDConvex(image);
        expect(nonzeroCoords(image)).toEqual(coordKeys(ones));
    });

    it('keeps a single voxel and clears everything else', () => {
        const image = makeImage(5, 5, 5, [[2, 2, 2]]);
        ImageUtility3.computeCDConvex(image);
        expect(nonzeroCoords(image)).toEqual(['2,2,2']);
    });

    it('clears an all-background image', () => {
        const image = new Image3<number>(4, 4, 4);
        ImageUtility3.computeCDConvex(image);
        expect(image.getPixels().every((v) => v === 0)).toBe(true);
    });
});

describe('ImageUtility3.floodFill6', () => {
    it('fills only the 6-connected region on one side of a wall', () => {
        const image = new Image3<number>(4, 4, 4);
        for (let z = 0; z < 4; ++z) {
            for (let y = 0; y < 4; ++y) {
                image.set(2, y, z, 5);
            }
        }
        ImageUtility3.floodFill6(image, 0, 0, 0, 1, 0);
        for (let z = 0; z < 4; ++z) {
            for (let y = 0; y < 4; ++y) {
                expect(image.get(0, y, z)).toBe(1);
                expect(image.get(1, y, z)).toBe(1);
                expect(image.get(2, y, z)).toBe(5);
                expect(image.get(3, y, z)).toBe(0);
            }
        }
    });

    it('does nothing for a seed outside the image', () => {
        const image = new Image3<number>(3, 3, 3);
        ImageUtility3.floodFill6(image, 0, 0, 3, 1, 0);
        ImageUtility3.floodFill6(image, -1, 0, 0, 1, 0);
        expect(image.getPixels().every((v) => v === 0)).toBe(true);
    });

    it('fills the whole image when there is no barrier', () => {
        const image = new Image3<number>(3, 4, 5);
        ImageUtility3.floodFill6(image, 1, 1, 1, 9, 0);
        expect(image.getPixels().every((v) => v === 9)).toBe(true);
    });
});

describe('ImageUtility3.drawLine', () => {
    it('reproduces the 3D Bresenham voxels', () => {
        const visited: Voxel[] = [];
        ImageUtility3.drawLine(0, 0, 0, 4, 2, 1, (x, y, z) => visited.push([x, y, z]));
        expect(visited).toEqual([[0, 0, 0], [1, 1, 0], [2, 1, 1], [3, 2, 1], [4, 2, 1]]);
    });

    it('handles the y- and z-dominant cases and a degenerate line', () => {
        const yDominant: Voxel[] = [];
        ImageUtility3.drawLine(0, 0, 0, 1, 3, 0, (x, y, z) => yDominant.push([x, y, z]));
        expect(yDominant.length).toBe(4);
        expect(yDominant[0]).toEqual([0, 0, 0]);
        expect(yDominant[3]).toEqual([1, 3, 0]);

        const zDominant: Voxel[] = [];
        ImageUtility3.drawLine(0, 0, 0, 1, 1, 5, (x, y, z) => zDominant.push([x, y, z]));
        expect(zDominant.length).toBe(6);
        expect(zDominant[0]).toEqual([0, 0, 0]);
        expect(zDominant[5]).toEqual([1, 1, 5]);

        const point: Voxel[] = [];
        ImageUtility3.drawLine(3, 4, 5, 3, 4, 5, (x, y, z) => point.push([x, y, z]));
        expect(point).toEqual([[3, 4, 5]]);
    });

    it('visits max(|dx|,|dy|,|dz|)+1 voxels and reaches the endpoint', () => {
        const random = makeRandom(13579);
        for (let trial = 0; trial < 40; ++trial) {
            const p: number[] = [];
            for (let k = 0; k < 6; ++k) {
                p.push(Math.floor(random() * 17) - 8);
            }
            const visited: Voxel[] = [];
            ImageUtility3.drawLine(p[0], p[1], p[2], p[3], p[4], p[5],
                (x, y, z) => visited.push([x, y, z]));
            const expectedCount = Math.max(Math.abs(p[3] - p[0]), Math.abs(p[4] - p[1]),
                Math.abs(p[5] - p[2])) + 1;
            expect(visited.length).toBe(expectedCount);
            expect(visited[0]).toEqual([p[0], p[1], p[2]]);
            expect(visited[visited.length - 1]).toEqual([p[3], p[4], p[5]]);
            // Consecutive voxels are 26-connected.
            for (let i = 1; i < visited.length; ++i) {
                for (let d = 0; d < 3; ++d) {
                    expect(Math.abs(visited[i][d] - visited[i - 1][d])).toBeLessThanOrEqual(1);
                }
            }
        }
    });

    it('is symmetric in the dominant direction for a diagonal line', () => {
        const forward: string[] = [];
        ImageUtility3.drawLine(0, 0, 0, 5, 5, 5, (x, y, z) => forward.push(`${x},${y},${z}`));
        expect(forward).toEqual(['0,0,0', '1,1,1', '2,2,2', '3,3,3', '4,4,4', '5,5,5']);
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V26): property-based cross-checks against the upstream
// header ImageUtility3.h.
// ---------------------------------------------------------------------------

describe('ImageUtility3 verification', () => {
    const D0 = 5;
    const D1 = 4;
    const D2 = 4;
    const N = D0 * D1 * D2;

    const bitsArb = fc.array(fc.integer({ min: 0, max: 1 }),
        { minLength: N, maxLength: N });

    // An image with no boundary constraint (dilate and erode accept these).
    const freeArb = bitsArb.map((bits) => {
        const image = new Image3<number>(D0, D1, D2);
        const voxels = image.getPixels();
        for (let i = 0; i < bits.length; ++i) {
            voxels[i] = bits[i];
        }
        return image;
    });

    // An image with zeros on the boundary, as the header requires for
    // getComponents, computeCDConvex and floodFill6.
    const borderedArb = bitsArb.map((bits) => {
        const image = new Image3<number>(D0, D1, D2);
        const voxels = image.getPixels();
        for (let z = 1; z < D2 - 1; ++z) {
            for (let y = 1; y < D1 - 1; ++y) {
                for (let x = 1; x < D0 - 1; ++x) {
                    const i = x + D0 * (y + D1 * z);
                    voxels[i] = bits[i];
                }
            }
        }
        return image;
    });

    const foreground = (image: Image3<number>): Set<number> => {
        const set = new Set<number>();
        const voxels = image.getPixels();
        for (let i = 0; i < voxels.length; ++i) {
            if (voxels[i] === 1) {
                set.add(i);
            }
        }
        return set;
    };

    const complement = (image: Image3<number>): Image3<number> => {
        const out = new Image3<number>(D0, D1, D2);
        const src = image.getPixels();
        const dst = out.getPixels();
        for (let i = 0; i < src.length; ++i) {
            dst[i] = 1 - src[i];
        }
        return out;
    };

    const clone = (image: Image3<number>): Image3<number> => {
        const out = new Image3<number>(D0, D1, D2);
        const src = image.getPixels();
        const dst = out.getPixels();
        for (let i = 0; i < src.length; ++i) {
            dst[i] = src[i];
        }
        return out;
    };

    it('dilation is extensive and erosion is anti-extensive', () => {
        check(freeArb, (image) => {
            const input = foreground(image);
            for (const dilated of [ImageUtility3.dilate6(image),
                ImageUtility3.dilate18(image), ImageUtility3.dilate26(image)]) {
                const out = foreground(dilated);
                for (const i of input) {
                    expect(out.has(i)).toBe(true);
                }
            }
            for (const zeroExterior of [false, true]) {
                for (const eroded of [
                    ImageUtility3.erode6(image, zeroExterior),
                    ImageUtility3.erode18(image, zeroExterior),
                    ImageUtility3.erode26(image, zeroExterior)]) {
                    for (const i of foreground(eroded)) {
                        expect(input.has(i)).toBe(true);
                    }
                }
            }
        }, 40);
    });

    it('the 6-, 18- and 26-neighborhoods are nested', () => {
        check(freeArb, (image) => {
            const d6 = foreground(ImageUtility3.dilate6(image));
            const d18 = foreground(ImageUtility3.dilate18(image));
            const d26 = foreground(ImageUtility3.dilate26(image));
            for (const i of d6) {
                expect(d18.has(i)).toBe(true);
            }
            for (const i of d18) {
                expect(d26.has(i)).toBe(true);
            }
            for (const zeroExterior of [false, true]) {
                const e6 = foreground(ImageUtility3.erode6(image,
                    zeroExterior));
                const e18 = foreground(ImageUtility3.erode18(image,
                    zeroExterior));
                const e26 = foreground(ImageUtility3.erode26(image,
                    zeroExterior));
                for (const i of e26) {
                    expect(e18.has(i)).toBe(true);
                }
                for (const i of e18) {
                    expect(e6.has(i)).toBe(true);
                }
            }
        }, 40);
    });

    it('dilation and erosion are dual under complementation', () => {
        check(freeArb, (image) => {
            for (const [dilate, erode] of [
                [ImageUtility3.dilate6, ImageUtility3.erode6],
                [ImageUtility3.dilate18, ImageUtility3.erode18],
                [ImageUtility3.dilate26, ImageUtility3.erode26]] as const) {
                const eroded = erode(image, false);
                const dual = complement(dilate(complement(image)));
                expect(nonzeroCoords(eroded)).toEqual(nonzeroCoords(dual));
            }
        }, 40);
    });

    it('dilation commutes with the x-reflection of the image', () => {
        // The 6-, 18- and 26-neighborhoods are symmetric under x -> -x, so
        // this must hold. Upstream starts the innermost loop at i0 = 1 and
        // never uses an x = 0 voxel as a dilation source, which breaks it;
        // the port starts at i0 = 0.
        check(freeArb, (image) => {
            const flip = (source: Image3<number>) => {
                const out = new Image3<number>(D0, D1, D2);
                for (let z = 0; z < D2; ++z) {
                    for (let y = 0; y < D1; ++y) {
                        for (let x = 0; x < D0; ++x) {
                            out.set(D0 - 1 - x, y, z, source.get(x, y, z));
                        }
                    }
                }
                return out;
            };
            for (const dilate of [ImageUtility3.dilate6,
                ImageUtility3.dilate18, ImageUtility3.dilate26]) {
                expect(nonzeroCoords(dilate(flip(image))))
                    .toEqual(nonzeroCoords(flip(dilate(image))));
            }
        }, 40);
    });

    it('getComponents partitions the foreground into maximal components',
        () => {
            const offsets6: number[][] = [];
            for (let d = 0; d < 3; ++d) {
                for (const s of [-1, 1]) {
                    const o = [0, 0, 0];
                    o[d] = s;
                    offsets6.push(o);
                }
            }
            const offsets18: number[][] = [];
            const offsets26: number[][] = [];
            for (let dz = -1; dz <= 1; ++dz) {
                for (let dy = -1; dy <= 1; ++dy) {
                    for (let dx = -1; dx <= 1; ++dx) {
                        const sum = Math.abs(dx) + Math.abs(dy) + Math.abs(dz);
                        if (sum === 0) {
                            continue;
                        }
                        offsets26.push([dx, dy, dz]);
                        if (sum <= 2) {
                            offsets18.push([dx, dy, dz]);
                        }
                    }
                }
            }

            check(borderedArb, (image) => {
                for (const [call, offsets] of [
                    [ImageUtility3.getComponents6, offsets6],
                    [ImageUtility3.getComponents18, offsets18],
                    [ImageUtility3.getComponents26, offsets26]] as const) {
                    const work = clone(image);
                    const source = [...image.getPixels()];
                    const expected = referenceComponents(source, D0, D1, D2,
                        offsets);
                    const actual = call(work);
                    if (expected.length === 1) {
                        expect(actual.length).toBe(0);
                        continue;
                    }
                    expect(actual.length).toBe(expected.length);
                    for (let k = 1; k < expected.length; ++k) {
                        expect([...actual[k]].sort((p, q) => p - q))
                            .toEqual(expected[k]);
                        for (const i of actual[k]) {
                            expect(work.get(i)).toBe(k);
                        }
                    }
                }
            }, 40);
        });

    it('computeCDConvex fills the coordinate-directional spans', () => {
        check(borderedArb, (image) => {
            const source = [...image.getPixels()];
            const work = clone(image);
            ImageUtility3.computeCDConvex(work);

            // A voxel is in the result exactly when, along each of the three
            // coordinate directions, its line carries a foreground voxel at
            // or below it and one at or above it.
            for (let z = 0; z < D2; ++z) {
                for (let y = 0; y < D1; ++y) {
                    for (let x = 0; x < D0; ++x) {
                        const p = [x, y, z];
                        let inside = true;
                        for (let d = 0; d < 3 && inside; ++d) {
                            const bound = [D0, D1, D2][d];
                            let below = false;
                            let above = false;
                            for (let t = 0; t < bound; ++t) {
                                const q = p.slice();
                                q[d] = t;
                                if (source[q[0] + D0 * (q[1] + D1 * q[2])]
                                    === 1) {
                                    if (t <= p[d]) {
                                        below = true;
                                    }
                                    if (t >= p[d]) {
                                        above = true;
                                    }
                                }
                            }
                            inside = below && above;
                        }
                        expect(work.get(x, y, z)).toBe(inside ? 1 : 0);
                    }
                }
            }

            // The result contains the input.
            for (let i = 0; i < source.length; ++i) {
                if (source[i] === 1) {
                    expect(work.getPixels()[i]).toBe(1);
                }
            }
        }, 40);
    });

    it('floodFill6 fills exactly the 6-connected region of the seed', () => {
        check(fc.tuple(borderedArb, fc.integer({ min: 1, max: D0 - 2 }),
            fc.integer({ min: 1, max: D1 - 2 }),
            fc.integer({ min: 1, max: D2 - 2 })), ([image, sx, sy, sz]) => {
            const source = [...image.getPixels()];
            const before = image.get(sx, sy, sz);
            const filled = clone(image);
            ImageUtility3.floodFill6(filled, sx, sy, sz, 7, before);

            const region = new Set<number>();
            const stack = [sx + D0 * (sy + D1 * sz)];
            region.add(stack[0]);
            while (stack.length > 0) {
                const v = stack.pop() as number;
                const x = v % D0;
                const y = Math.floor(v / D0) % D1;
                const z = Math.floor(v / (D0 * D1));
                for (let d = 0; d < 3; ++d) {
                    for (const s of [-1, 1]) {
                        const q = [x, y, z];
                        q[d] += s;
                        if (0 <= q[0] && q[0] < D0 && 0 <= q[1] && q[1] < D1
                            && 0 <= q[2] && q[2] < D2) {
                            const w = q[0] + D0 * (q[1] + D1 * q[2]);
                            if (source[w] === before && !region.has(w)) {
                                region.add(w);
                                stack.push(w);
                            }
                        }
                    }
                }
            }

            const out = filled.getPixels();
            for (let i = 0; i < out.length; ++i) {
                expect(out[i]).toBe(region.has(i) ? 7 : source[i]);
            }
        }, 40);
    });

    it('drawLine visits a connected chain from (x0,y0,z0) to (x1,y1,z1)',
        () => {
            const coord = fc.integer({ min: -4, max: 4 });
            check(fc.tuple(coord, coord, coord, coord, coord, coord),
                ([x0, y0, z0, x1, y1, z1]) => {
                    const visited: [number, number, number][] = [];
                    ImageUtility3.drawLine(x0, y0, z0, x1, y1, z1,
                        (x, y, z) => visited.push([x, y, z]));
                    expect(visited.length).toBe(1 + Math.max(
                        Math.abs(x1 - x0), Math.abs(y1 - y0),
                        Math.abs(z1 - z0)));
                    expect(visited[0]).toEqual([x0, y0, z0]);
                    expect(visited[visited.length - 1]).toEqual([x1, y1, z1]);
                    for (let i = 1; i < visited.length; ++i) {
                        const d = [0, 1, 2].map((k) =>
                            Math.abs(visited[i][k] - visited[i - 1][k]));
                        // A single voxel step in the dominant direction.
                        expect(Math.max(...d)).toBe(1);
                    }
                }, 100);
        });

    it('the morphology operations do not modify their input', () => {
        check(freeArb, (image) => {
            const source = [...image.getPixels()];
            ImageUtility3.dilate26(image);
            ImageUtility3.erode26(image, true);
            ImageUtility3.open18(image, false);
            ImageUtility3.close18(image, false);
            expect([...image.getPixels()]).toEqual(source);
        }, 40);
    });

    it('open and close are the documented compositions', () => {
        check(freeArb, (image) => {
            for (const zeroExterior of [false, true]) {
                expect(nonzeroCoords(ImageUtility3.open6(image, zeroExterior)))
                    .toEqual(nonzeroCoords(ImageUtility3.dilate6(
                        ImageUtility3.erode6(image, zeroExterior))));
                expect(nonzeroCoords(ImageUtility3.close26(image,
                    zeroExterior)))
                    .toEqual(nonzeroCoords(ImageUtility3.erode26(
                        ImageUtility3.dilate26(image), zeroExterior)));
                // The caller-specified overloads agree with the named ones.
                expect(nonzeroCoords(ImageUtility3.open(image, zeroExterior,
                    image.getNeighborhood18Coords())))
                    .toEqual(nonzeroCoords(ImageUtility3.open18(image,
                        zeroExterior)));
                expect(nonzeroCoords(ImageUtility3.close(image, zeroExterior,
                    image.getNeighborhood18Coords())))
                    .toEqual(nonzeroCoords(ImageUtility3.close18(image,
                        zeroExterior)));
            }
        }, 30);
    });
});
