import { describe, it, expect } from 'vitest';
import { Image } from '../src/Image.js';
import { Image3 } from '../src/Image3.js';
import { check, fc } from './helpers/arbitraries.js';

// Lexicographical index for an image with dimensions (d0,d1,d2).
function linear(d0: number, d1: number, x: number, y: number, z: number): number {
    return x + d0 * (y + d1 * z);
}

// Fill an image so that each voxel holds its own linear index times 10, which
// makes accessor mistakes obvious.
function makeImage(d0: number, d1: number, d2: number): Image3<number> {
    const image = new Image3<number>(d0, d1, d2);
    for (let i = 0; i < image.getNumPixels(); ++i) {
        image.set(i, 10 * i);
    }
    return image;
}

describe('Image3', () => {
    it('defaults to an empty image', () => {
        const image = new Image3<number>();
        expect(image.getNumDimensions()).toBe(0);
        expect(image.getNumPixels()).toBe(0);
        expect(image instanceof Image).toBe(true);
    });

    it('constructs with the dimensions, offsets and pixel count of the base', () => {
        const image = new Image3<number>(4, 3, 2);
        expect(image.getNumDimensions()).toBe(3);
        expect(image.getDimensions()).toEqual([4, 3, 2]);
        expect(image.getOffsets()).toEqual([1, 4, 12]);
        expect(image.getNumPixels()).toBe(24);
        expect(image.getPixels()).toEqual(new Array(24).fill(0));
    });

    it('leaves the image empty for a nonpositive dimension', () => {
        const image = new Image3<number>(4, 0, 2);
        expect(image.getNumDimensions()).toBe(0);
        expect(image.getNumPixels()).toBe(0);
    });

    it('uses the pixel factory so pixels do not alias', () => {
        const image = new Image3<number[]>(2, 2, 2, () => [0]);
        expect(image.getNumPixels()).toBe(8);
        image.get(0, 0, 0)[0] = 7;
        expect(image.get(1, 0, 0)[0]).toBe(0);
    });

    it('reconstructs from three dimensions or from a dimension array', () => {
        const image = new Image3<number>(2, 2, 2);
        image.set(0, 5);
        image.reconstruct(3, 4, 5);
        expect(image.getDimensions()).toEqual([3, 4, 5]);
        expect(image.getNumPixels()).toBe(60);
        expect(image.get(0)).toBe(0);

        image.reconstruct([2, 3, 4]);
        expect(image.getDimensions()).toEqual([2, 3, 4]);
        expect(image.getNumPixels()).toBe(24);
    });

    it('converts between 3-tuples and 1-dimensional indices', () => {
        const image = new Image3<number>(4, 3, 2);
        for (let z = 0; z < 2; ++z) {
            for (let y = 0; y < 3; ++y) {
                for (let x = 0; x < 4; ++x) {
                    const index = linear(4, 3, x, y, z);
                    expect(image.getIndex(x, y, z)).toBe(index);
                    expect(image.getIndex([x, y, z])).toBe(index);
                    expect(image.getCoordinates(index)).toEqual([x, y, z]);
                }
            }
        }
        // The 3D index agrees with the n-dimensional base formula.
        expect(image.getIndex(2, 1, 1)).toBe(18);
        expect(image.getCoordinates(18)).toEqual([2, 1, 1]);
    });

    it('reads and writes voxels through all accessor forms', () => {
        const image = makeImage(4, 3, 2);
        expect(image.get(2, 1, 1)).toBe(180);
        expect(image.get([2, 1, 1])).toBe(180);
        expect(image.get(18)).toBe(180);

        image.set(2, 1, 1, -1);
        expect(image.get(18)).toBe(-1);
        image.set([2, 1, 1], -2);
        expect(image.get(18)).toBe(-2);
        image.set(18, -3);
        expect(image.get(18)).toBe(-3);
    });

    it('clamps out-of-range coordinates in getClamped/setClamped', () => {
        const image = makeImage(3, 2, 2);
        // (-5,-5,-5) clamps to (0,0,0) and (100,100,100) clamps to (2,1,1).
        expect(image.getClamped(-5, -5, -5)).toBe(0);
        expect(image.getClamped(100, 100, 100)).toBe(10 * linear(3, 2, 2, 1, 1));
        expect(image.getClamped([-1, 3, 0])).toBe(10 * linear(3, 2, 0, 1, 0));
        // Clamping is per-component, so an in-range component is preserved.
        expect(image.getClamped(1, -1, 1)).toBe(10 * linear(3, 2, 1, 0, 1));

        image.setClamped(-5, 0, 0, 99);
        expect(image.get(0, 0, 0)).toBe(99);
        image.setClamped([100, 100, 100], 77);
        expect(image.get(2, 1, 1)).toBe(77);

        // The one-argument forms are the base-class linear-index versions,
        // which clamp an invalid index to element 0.
        expect(image.getClamped(-1)).toBe(99);
        image.setClamped(1000, 55);
        expect(image.get(0)).toBe(55);
    });

    it('does not modify the caller coordinate array when clamping', () => {
        const image = makeImage(3, 2, 2);
        const coord = [-4, 9, 0];
        image.getClamped(coord);
        expect(coord).toEqual([-4, 9, 0]);
    });

    it('builds the relative 1-dimensional neighborhoods from the dimensions', () => {
        const image = new Image3<number>(4, 3, 2);
        const dim0 = 4;
        const dim01 = 12;
        expect(image.getNeighborhood6()).toEqual([-1, 1, -dim0, dim0, -dim01, dim01]);
        expect(image.getNeighborhood18().length).toBe(18);
        expect(image.getNeighborhood26().length).toBe(26);
        expect(image.getCorners8()).toEqual([
            0, 1, dim0, dim0 + 1, dim01, dim01 + 1, dim01 + dim0, dim01 + dim0 + 1
        ]);
        expect(image.getFull27().length).toBe(27);
        // The full neighborhood is centered on the voxel itself.
        expect(image.getFull27()[13]).toBe(0);
        // The 18-connected list extends the 6-connected list, and the
        // 26-connected list extends the 18-connected list.
        expect(image.getNeighborhood18().slice(0, 6)).toEqual(image.getNeighborhood6());
        expect(image.getNeighborhood26().slice(0, 18)).toEqual(image.getNeighborhood18());
    });

    it('keeps the relative offsets and relative 3-tuples consistent', () => {
        const image = new Image3<number>(5, 4, 3);
        const dim0 = 5;
        const dim01 = 20;
        const toOffset = (c: number[]) => c[0] + dim0 * c[1] + dim01 * c[2];
        const pairs: Array<[number[], number[][]]> = [
            [image.getNeighborhood6(), image.getNeighborhood6Coords()],
            [image.getNeighborhood18(), image.getNeighborhood18Coords()],
            [image.getNeighborhood26(), image.getNeighborhood26Coords()],
            [image.getCorners8(), image.getCorners8Coords()],
            [image.getFull27(), image.getFull27Coords()]
        ];
        for (const [offsets, coords] of pairs) {
            expect(coords.length).toBe(offsets.length);
            expect(coords.map(toOffset)).toEqual(offsets);
        }
    });

    it('produces the expected relative 3-tuple sets', () => {
        const image = new Image3<number>(5, 4, 3);
        const key = (c: number[]) => c.join(',');

        // The 26-connected neighborhood is the 3x3x3 block minus the center.
        const n26 = new Set(image.getNeighborhood26Coords().map(key));
        expect(n26.size).toBe(26);
        expect(n26.has('0,0,0')).toBe(false);

        // The full neighborhood is the whole 3x3x3 block.
        const full = new Set(image.getFull27Coords().map(key));
        expect(full.size).toBe(27);
        for (let z = -1; z <= 1; ++z) {
            for (let y = -1; y <= 1; ++y) {
                for (let x = -1; x <= 1; ++x) {
                    expect(full.has(key([x, y, z]))).toBe(true);
                }
            }
        }

        // The 6-connected neighbors are the face neighbors, the 18-connected
        // neighbors add the edge neighbors, and the 26-connected neighbors add
        // the corner neighbors.
        const l1 = (c: number[]) => Math.abs(c[0]) + Math.abs(c[1]) + Math.abs(c[2]);
        expect(image.getNeighborhood6Coords().every(c => l1(c) === 1)).toBe(true);
        expect(image.getNeighborhood18Coords().every(c => l1(c) <= 2)).toBe(true);
        expect(image.getNeighborhood18Coords().filter(c => l1(c) === 2).length).toBe(12);
        expect(image.getNeighborhood26Coords().filter(c => l1(c) === 3).length).toBe(8);

        // The corners are the unit cube with the voxel as its minimum corner.
        expect(image.getCorners8Coords().every(c => c.every(v => v === 0 || v === 1)))
            .toBe(true);
    });

    it('converts the relative neighborhoods to absolute indices at a voxel', () => {
        const image = new Image3<number>(4, 3, 2);
        const index = image.getIndex(1, 1, 1);
        expect(index).toBe(17);
        expect(image.getNeighborhood6(1, 1, 1)).toEqual(
            image.getNeighborhood6().map(o => index + o));
        expect(image.getNeighborhood6(1, 1, 1)).toEqual([16, 18, 13, 21, 5, 29]);
        expect(image.getFull27(1, 1, 1)[13]).toBe(index);
        expect(image.getCorners8(0, 0, 0)).toEqual([0, 1, 4, 5, 12, 13, 16, 17]);

        // No clamping on the boundary: the 1-dimensional indices may leave the
        // image, and the 3-tuples may be negative.
        expect(image.getNeighborhood6(0, 0, 0)).toEqual([-1, 1, -4, 4, -12, 12]);
        expect(image.getNeighborhood6Coords(0, 0, 0)).toEqual([
            [-1, 0, 0], [1, 0, 0], [0, -1, 0], [0, 1, 0], [0, 0, -1], [0, 0, 1]
        ]);
    });

    it('keeps the absolute indices and absolute 3-tuples consistent', () => {
        const image = new Image3<number>(5, 4, 3);
        const [x, y, z] = [2, 2, 1];
        const toIndex = (c: number[]) => linear(5, 4, c[0], c[1], c[2]);
        const pairs: Array<[number[], number[][]]> = [
            [image.getNeighborhood6(x, y, z), image.getNeighborhood6Coords(x, y, z)],
            [image.getNeighborhood18(x, y, z), image.getNeighborhood18Coords(x, y, z)],
            [image.getNeighborhood26(x, y, z), image.getNeighborhood26Coords(x, y, z)],
            [image.getCorners8(x, y, z), image.getCorners8Coords(x, y, z)],
            [image.getFull27(x, y, z), image.getFull27Coords(x, y, z)]
        ];
        for (const [indices, coords] of pairs) {
            expect(coords.map(toIndex)).toEqual(indices);
        }
    });

    it('returns fresh arrays so callers cannot corrupt later results', () => {
        const image = new Image3<number>(4, 3, 2);
        const first = image.getNeighborhood6();
        first[0] = 12345;
        expect(image.getNeighborhood6()[0]).toBe(-1);

        const firstCoords = image.getFull27Coords();
        firstCoords[0][0] = 12345;
        expect(image.getFull27Coords()[0]).toEqual([-1, -1, -1]);
    });

    it('reads the 26 neighbors of an interior voxel correctly', () => {
        const d0 = 5, d1 = 4, d2 = 3;
        const image = makeImage(d0, d1, d2);
        const [x, y, z] = [2, 2, 1];
        const indices = image.getNeighborhood26(x, y, z);
        const coords = image.getNeighborhood26Coords(x, y, z);
        for (let i = 0; i < 26; ++i) {
            expect(image.get(indices[i])).toBe(10 * indices[i]);
            expect(image.get(coords[i][0], coords[i][1], coords[i][2]))
                .toBe(10 * linear(d0, d1, coords[i][0], coords[i][1], coords[i][2]));
        }
    });
});

describe('Image3 verification', () => {
    // Random dimensions plus an in-range voxel derived from them, so no
    // filtering is needed and shrinking stays deterministic.
    const dimsAndVoxel = fc.tuple(
        fc.integer({ min: 1, max: 9 }),
        fc.integer({ min: 1, max: 9 }),
        fc.integer({ min: 1, max: 9 }),
        fc.nat({ max: 1000 }), fc.nat({ max: 1000 }), fc.nat({ max: 1000 }))
        .map(([d0, d1, d2, rx, ry, rz]) => ({
            d0, d1, d2, x: rx % d0, y: ry % d1, z: rz % d2
        }));

    it('getIndex and getCoordinates are mutual inverses', () => {
        check(dimsAndVoxel, ({ d0, d1, d2, x, y, z }) => {
            const image = new Image3<number>(d0, d1, d2);

            // (x,y,z) -> index -> (x,y,z)
            const i = image.getIndex(x, y, z);
            expect(i).toBe(linear(d0, d1, x, y, z));
            expect(image.getCoordinates(i)).toEqual([x, y, z]);
            // The array overload agrees with the three-argument overload.
            expect(image.getIndex([x, y, z])).toBe(i);

            // index -> (x,y,z) -> index, over every voxel of the image.
            for (let k = 0; k < d0 * d1 * d2; ++k) {
                const c = image.getCoordinates(k);
                expect(c[0] >= 0 && c[0] < d0).toBe(true);
                expect(c[1] >= 0 && c[1] < d1).toBe(true);
                expect(c[2] >= 0 && c[2] < d2).toBe(true);
                expect(image.getIndex(c)).toBe(k);
            }
        });
    });

    it('the relative 1-dimensional offsets linearize the relative 3-tuples', () => {
        check(dimsAndVoxel, ({ d0, d1, d2 }) => {
            const image = new Image3<number>(d0, d1, d2);
            const cases: [number[], [number, number, number][]][] = [
                [image.getNeighborhood6(), image.getNeighborhood6Coords()],
                [image.getNeighborhood18(), image.getNeighborhood18Coords()],
                [image.getNeighborhood26(), image.getNeighborhood26Coords()],
                [image.getCorners8(), image.getCorners8Coords()],
                [image.getFull27(), image.getFull27Coords()]
            ];
            for (const [offsets, coords] of cases) {
                expect(offsets.length).toBe(coords.length);
                for (let k = 0; k < offsets.length; ++k) {
                    // The offset must be the linearization of the 3-tuple,
                    // which is what makes the two families interchangeable.
                    expect(offsets[k]).toBe(
                        coords[k][0] + d0 * (coords[k][1] + d1 * coords[k][2]));
                }
            }
        });
    });

    it('the relative 3-tuple sets are the documented neighborhoods', () => {
        check(dimsAndVoxel, ({ d0, d1, d2 }) => {
            const image = new Image3<number>(d0, d1, d2);
            const key = (c: readonly number[]) => c[0] + ',' + c[1] + ',' + c[2];
            const setOf = (cs: [number, number, number][]) => new Set(cs.map(key));

            // Brute-force the 3x3x3 stencil and classify by the number of
            // nonzero offsets, which is the upstream definition.
            const all: [number, number, number][] = [];
            for (let w = -1; w <= 1; ++w) {
                for (let v = -1; v <= 1; ++v) {
                    for (let u = -1; u <= 1; ++u) {
                        all.push([u, v, w]);
                    }
                }
            }
            const nnz = (c: readonly number[]) =>
                Math.abs(c[0]) + Math.abs(c[1]) + Math.abs(c[2]);

            const n6 = setOf(image.getNeighborhood6Coords());
            const n18 = setOf(image.getNeighborhood18Coords());
            const n26 = setOf(image.getNeighborhood26Coords());
            const full = setOf(image.getFull27Coords());
            expect(n6.size).toBe(6);
            expect(n18.size).toBe(18);
            expect(n26.size).toBe(26);
            expect(full.size).toBe(27);
            for (const c of all) {
                expect(n6.has(key(c))).toBe(nnz(c) === 1);
                expect(n18.has(key(c))).toBe(nnz(c) === 1 || nnz(c) === 2);
                expect(n26.has(key(c))).toBe(nnz(c) >= 1);
                expect(full.has(key(c))).toBe(true);
            }

            // The corners are the 2x2x2 cube of nonnegative offsets.
            const corners = setOf(image.getCorners8Coords());
            expect(corners.size).toBe(8);
            for (const c of all) {
                const isCorner = c[0] >= 0 && c[1] >= 0 && c[2] >= 0;
                expect(corners.has(key(c))).toBe(isCorner);
            }
        });
    });

    it('absolute neighborhoods translate the relative ones without clamping', () => {
        check(dimsAndVoxel, ({ d0, d1, d2, x, y, z }) => {
            const image = new Image3<number>(d0, d1, d2);
            const relative: [number[], [number, number, number][]][] = [
                [image.getNeighborhood6(x, y, z), image.getNeighborhood6Coords()],
                [image.getNeighborhood18(x, y, z), image.getNeighborhood18Coords()],
                [image.getNeighborhood26(x, y, z), image.getNeighborhood26Coords()],
                [image.getCorners8(x, y, z), image.getCorners8Coords()],
                [image.getFull27(x, y, z), image.getFull27Coords()]
            ];
            const absoluteCoords: [number, number, number][][] = [
                image.getNeighborhood6Coords(x, y, z),
                image.getNeighborhood18Coords(x, y, z),
                image.getNeighborhood26Coords(x, y, z),
                image.getCorners8Coords(x, y, z),
                image.getFull27Coords(x, y, z)
            ];
            const base = linear(d0, d1, x, y, z);
            for (let s = 0; s < relative.length; ++s) {
                const [absolute, relCoords] = relative[s];
                for (let k = 0; k < absolute.length; ++k) {
                    // 1-D indices are the base index plus the relative offset.
                    expect(absolute[k]).toBe(base
                        + relCoords[k][0] + d0 * (relCoords[k][1] + d1 * relCoords[k][2]));
                    // 3-tuples are plain translations: an out-of-range neighbor
                    // stays negative (or beyond the bound) instead of wrapping
                    // around as the upstream size_t arithmetic does.
                    expect(absoluteCoords[s][k]).toEqual([
                        x + relCoords[k][0], y + relCoords[k][1], z + relCoords[k][2]]);
                }
            }
        });
    });

    it('getClamped and setClamped agree with an explicit clamp', () => {
        const offset = fc.integer({ min: -2, max: 11 });
        check(fc.tuple(dimsAndVoxel, offset, offset, offset),
            ([{ d0, d1, d2 }, x, y, z]) => {
                const image = makeImage(d0, d1, d2);
                const cx = Math.min(Math.max(x, 0), d0 - 1);
                const cy = Math.min(Math.max(y, 0), d1 - 1);
                const cz = Math.min(Math.max(z, 0), d2 - 1);
                const expected = 10 * linear(d0, d1, cx, cy, cz);
                expect(image.getClamped(x, y, z)).toBe(expected);
                expect(image.getClamped([x, y, z])).toBe(expected);

                image.setClamped(x, y, z, -7);
                expect(image.get(cx, cy, cz)).toBe(-7);
                image.setClamped([x, y, z], -9);
                expect(image.get(linear(d0, d1, cx, cy, cz))).toBe(-9);
            });
    });

    it('each neighborhood call returns a fresh array', () => {
        check(dimsAndVoxel, ({ d0, d1, d2, x, y, z }) => {
            const image = new Image3<number>(d0, d1, d2);
            const a = image.getFull27(x, y, z);
            a[0] = 12345;
            expect(image.getFull27(x, y, z)[0]).not.toBe(12345);
            const c = image.getNeighborhood26Coords(x, y, z);
            c[0][0] = 12345;
            expect(image.getNeighborhood26Coords(x, y, z)[0][0]).not.toBe(12345);
        });
    });
});
