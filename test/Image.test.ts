import { describe, it, expect } from 'vitest';
import { Image } from '../src/Image.js';
import { check, fc } from './helpers/arbitraries.js';

describe('Image', () => {
    it('defaults to an empty image', () => {
        const image = new Image<number>();
        expect(image.getNumDimensions()).toBe(0);
        expect(image.getDimensions()).toEqual([]);
        expect(image.getOffsets()).toEqual([]);
        expect(image.getNumPixels()).toBe(0);
    });

    it('constructs from dimensions with correct offsets and pixel count', () => {
        const image = new Image<number>([4, 3, 2]);
        expect(image.getNumDimensions()).toBe(3);
        expect(image.getDimension(0)).toBe(4);
        expect(image.getDimension(1)).toBe(3);
        expect(image.getDimension(2)).toBe(2);
        // offsets[0] = 1, offsets[d] = dim[d-1] * offsets[d-1]
        expect(image.getOffsets()).toEqual([1, 4, 12]);
        expect(image.getOffset(2)).toBe(12);
        expect(image.getNumPixels()).toBe(24);
    });

    it('zero-fills pixels by default', () => {
        const image = new Image<number>([3, 3]);
        expect(image.getPixels()).toEqual(new Array(9).fill(0));
    });

    it('uses the pixel factory when provided', () => {
        const image = new Image<number[]>([2, 2], () => [0, 0]);
        expect(image.getNumPixels()).toBe(4);
        // Each pixel must be a distinct object (no aliasing).
        image.get(0)[0] = 5;
        expect(image.get(1)[0]).toBe(0);
    });

    it('reconstruct replaces dimensions and loses pixel data', () => {
        const image = new Image<number>([2, 2]);
        image.set(3, 7);
        image.reconstruct([5]);
        expect(image.getNumDimensions()).toBe(1);
        expect(image.getNumPixels()).toBe(5);
        expect(image.get(3)).toBe(0);
    });

    it('reconstruct with a nonpositive dimension leaves the image empty', () => {
        const image = new Image<number>([2, 2]);
        image.reconstruct([3, 0]);
        expect(image.getNumDimensions()).toBe(0);
        expect(image.getNumPixels()).toBe(0);
        const image2 = new Image<number>([3, -1]);
        expect(image2.getNumPixels()).toBe(0);
    });

    it('converts between coordinates and flat index (both directions)', () => {
        const image = new Image<number>([4, 3, 2]);
        // index = c0 + 4 * c1 + 12 * c2
        expect(image.getIndex([1, 2, 1])).toBe(1 + 8 + 12);
        expect(image.getCoordinates(21)).toEqual([1, 2, 1]);
        // Round trip over the whole image.
        for (let i = 0; i < image.getNumPixels(); ++i) {
            expect(image.getIndex(image.getCoordinates(i))).toBe(i);
        }
    });

    it('get/set access the flat pixel array', () => {
        const image = new Image<number>([2, 3]);
        image.set(5, 42);
        expect(image.get(5)).toBe(42);
        expect(image.getPixels()[5]).toBe(42);
    });

    it('clamped accessors clamp invalid indices to element 0', () => {
        const image = new Image<number>([2, 2]);
        image.set(0, 11);
        image.set(3, 33);
        expect(image.getClamped(3)).toBe(33);
        expect(image.getClamped(4)).toBe(11);
        expect(image.getClamped(-1)).toBe(11);
        image.setClamped(100, 77);
        expect(image.get(0)).toBe(77);
        image.setClamped(2, 55);
        expect(image.get(2)).toBe(55);
    });

    it('supports pixels of generic type', () => {
        const image = new Image<string>([2, 2], () => '');
        image.set(2, 'abc');
        expect(image.get(2)).toBe('abc');
        expect(image.get(1)).toBe('');
    });
});

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream Image.h semantics rather than restating the implementation.
// ---------------------------------------------------------------------------

describe('Image verification', () => {
    // Small positive dimension lists; the image is fully enumerated below, so
    // the total pixel count is kept modest.
    const dimensions = fc.array(fc.integer({ min: 1, max: 6 }),
        { minLength: 1, maxLength: 4 })
        .filter(d => d.reduce((a, b) => a * b, 1) <= 400);

    it('round trips every flat index through getCoordinates/getIndex', () => {
        check(dimensions, dims => {
            const image = new Image<number>(dims);
            const numPixels = image.getNumPixels();
            expect(numPixels).toBe(dims.reduce((a, b) => a * b, 1));
            for (let i = 0; i < numPixels; ++i) {
                const coord = image.getCoordinates(i);
                expect(coord.length).toBe(dims.length);
                for (let d = 0; d < dims.length; ++d) {
                    expect(coord[d]).toBeGreaterThanOrEqual(0);
                    expect(coord[d]).toBeLessThan(dims[d]);
                }
                expect(image.getIndex(coord)).toBe(i);
            }
        }, 60);
    });

    it('linearizes coordinates as sum_d coord[d] * prod_{e<d} dim[e]', () => {
        check(dimensions, dims => {
            const image = new Image<number>(dims);
            // Independent computation of the offsets from the dimensions.
            const offsets: number[] = [];
            let product = 1;
            for (let d = 0; d < dims.length; ++d) {
                offsets.push(product);
                product *= dims[d];
            }
            expect([...image.getOffsets()]).toEqual(offsets);

            const numPixels = image.getNumPixels();
            for (let i = 0; i < numPixels; ++i) {
                const coord = image.getCoordinates(i);
                let index = 0;
                for (let d = 0; d < dims.length; ++d) {
                    index += coord[d] * offsets[d];
                }
                expect(image.getIndex(coord)).toBe(index);
            }
        }, 60);
    });

    it('dimension 0 varies fastest (lexicographical storage order)', () => {
        check(dimensions.filter(d => d.length >= 2), dims => {
            const image = new Image<number>(dims);
            // Consecutive flat indices differ only in coordinate 0 unless
            // coordinate 0 wraps.
            for (let i = 0; i + 1 < image.getNumPixels(); ++i) {
                const c0 = image.getCoordinates(i);
                const c1 = image.getCoordinates(i + 1);
                if (c0[0] + 1 < dims[0]) {
                    expect(c1[0]).toBe(c0[0] + 1);
                    expect(c1.slice(1)).toEqual(c0.slice(1));
                } else {
                    expect(c1[0]).toBe(0);
                }
            }
        }, 60);
    });

    it('reconstruct leaves the image empty when any dimension is nonpositive', () => {
        check(fc.tuple(dimensions, fc.nat({ max: 3 }), fc.integer({ min: -4, max: 0 })),
            ([dims, position, bad]) => {
                const image = new Image<number>(dims);
                const spoiled = dims.slice();
                spoiled[position % spoiled.length] = bad;
                image.reconstruct(spoiled);
                expect(image.getNumDimensions()).toBe(0);
                expect(image.getNumPixels()).toBe(0);
                expect([...image.getOffsets()]).toEqual([]);
            });
    });

    it('clamped accessors fall back to element 0 exactly outside [0, numPixels)', () => {
        check(fc.tuple(dimensions, fc.integer({ min: -20, max: 420 })),
            ([dims, i]) => {
                const image = new Image<number>(dims);
                const numPixels = image.getNumPixels();
                for (let k = 0; k < numPixels; ++k) {
                    image.set(k, 100 + k);
                }
                const inRange = 0 <= i && i < numPixels;
                expect(image.getClamped(i)).toBe(inRange ? 100 + i : 100);

                image.setClamped(i, -7);
                expect(image.get(inRange ? i : 0)).toBe(-7);
            });
    });

    it('reconstruct discards all previous pixel data', () => {
        check(fc.tuple(dimensions, dimensions), ([dims0, dims1]) => {
            const image = new Image<number>(dims0);
            for (let k = 0; k < image.getNumPixels(); ++k) {
                image.set(k, 1 + k);
            }
            image.reconstruct(dims1);
            expect(image.getNumPixels()).toBe(dims1.reduce((a, b) => a * b, 1));
            expect(image.getPixels().every(p => p === 0)).toBe(true);
        }, 60);
    });

    it('does not alias the caller dimensions array', () => {
        check(dimensions, dims => {
            const input = dims.slice();
            const image = new Image<number>(input);
            input[0] = 999;
            expect(image.getDimension(0)).toBe(dims[0]);
        });
    });
});
