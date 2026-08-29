import { describe, it, expect } from 'vitest';
import { Image } from '../src/Image';

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
