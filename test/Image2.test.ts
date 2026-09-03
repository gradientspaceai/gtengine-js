import { describe, it, expect } from 'vitest';
import { Image2 } from '../src/Image2.js';
import { Image } from '../src/Image.js';

describe('Image2', () => {
    it('defaults to an empty image', () => {
        const image = new Image2<number>();
        expect(image.getNumDimensions()).toBe(0);
        expect(image.getNumPixels()).toBe(0);
    });

    it('constructs from two dimensions and extends Image', () => {
        const image = new Image2<number>(4, 3);
        expect(image).toBeInstanceOf(Image);
        expect(image.getNumDimensions()).toBe(2);
        expect(image.getDimension(0)).toBe(4);
        expect(image.getDimension(1)).toBe(3);
        expect(image.getNumPixels()).toBe(12);
        expect(image.getPixels()).toEqual(new Array(12).fill(0));
    });

    it('uses the pixel factory when provided', () => {
        const image = new Image2<number[]>(2, 2, () => [0]);
        image.get(0, 0)[0] = 5;
        expect(image.get(1, 0)[0]).toBe(0);
    });

    it('reconstructs from two dimensions or a dimensions array', () => {
        const image = new Image2<number>(2, 2);
        image.reconstruct(5, 3);
        expect(image.getDimension(0)).toBe(5);
        expect(image.getDimension(1)).toBe(3);
        expect(image.getNumPixels()).toBe(15);

        image.reconstruct([3, 2]);
        expect(image.getDimension(0)).toBe(3);
        expect(image.getDimension(1)).toBe(2);
        expect(image.getNumPixels()).toBe(6);
    });

    it('converts (x,y) coordinates to 1-dimensional indices', () => {
        const image = new Image2<number>(5, 4);
        expect(image.getIndex(0, 0)).toBe(0);
        expect(image.getIndex(3, 2)).toBe(3 + 5 * 2);
        expect(image.getIndex([3, 2])).toBe(13);
        // Base class array signature still works through the override.
        for (let y = 0; y < 4; ++y) {
            for (let x = 0; x < 5; ++x) {
                expect(image.getIndex(x, y)).toBe(x + 5 * y);
            }
        }
    });

    it('converts 1-dimensional indices to (x,y) coordinates', () => {
        const image = new Image2<number>(5, 4);
        expect(image.getCoordinates(0)).toEqual([0, 0]);
        expect(image.getCoordinates(13)).toEqual([3, 2]);
        for (let i = 0; i < image.getNumPixels(); ++i) {
            const [x, y] = image.getCoordinates(i);
            expect(image.getIndex(x, y)).toBe(i);
        }
    });

    it('accesses pixels through flat, (x,y) and coordinate overloads', () => {
        const image = new Image2<number>(4, 3);
        image.set(2, 1, 42);
        expect(image.get(2, 1)).toBe(42);
        expect(image.get(2 + 4 * 1)).toBe(42);
        expect(image.get([2, 1])).toBe(42);

        image.set([1, 2], 7);
        expect(image.get(1, 2)).toBe(7);

        image.set(5, 9);
        expect(image.get(1, 1)).toBe(9);
    });

    it('clamped accessors clamp coordinates to the image boundary', () => {
        const image = new Image2<number>(3, 3);
        for (let y = 0; y < 3; ++y) {
            for (let x = 0; x < 3; ++x) {
                image.set(x, y, 10 * y + x);
            }
        }

        expect(image.getClamped(-1, 1)).toBe(image.get(0, 1));
        expect(image.getClamped(5, 1)).toBe(image.get(2, 1));
        expect(image.getClamped(1, -2)).toBe(image.get(1, 0));
        expect(image.getClamped(1, 7)).toBe(image.get(1, 2));
        expect(image.getClamped(-1, -1)).toBe(image.get(0, 0));
        expect(image.getClamped([4, 4])).toBe(image.get(2, 2));
        expect(image.getClamped(1, 1)).toBe(11);

        image.setClamped(-5, 1, 77);
        expect(image.get(0, 1)).toBe(77);
        image.setClamped([3, -1], 88);
        expect(image.get(2, 0)).toBe(88);
        image.setClamped(1, 1, 99);
        expect(image.get(1, 1)).toBe(99);

        // Flat-index clamped access retains base class behavior (clamps to
        // element 0).
        expect(image.getClamped(100)).toBe(image.get(0));
    });

    it('computes relative 1-dimensional neighborhoods', () => {
        const image = new Image2<number>(7, 5);
        const dim0 = 7;
        expect(image.getNeighborhood4()).toEqual([-1, 1, -dim0, dim0]);
        expect(image.getNeighborhood8()).toEqual([
            -1, 1, -dim0, dim0, -1 - dim0, 1 - dim0, -1 + dim0, 1 + dim0]);
        expect(image.getCorners()).toEqual([0, 1, dim0, dim0 + 1]);
        expect(image.getFull()).toEqual([
            -1 - dim0, -dim0, 1 - dim0, -1, 0, 1, -1 + dim0, dim0, 1 + dim0]);
    });

    it('computes absolute 1-dimensional neighborhoods at (x,y)', () => {
        const image = new Image2<number>(7, 5);
        const index = image.getIndex(3, 2);
        expect(image.getNeighborhood4(3, 2)).toEqual(
            image.getNeighborhood4().map((offset) => index + offset));
        expect(image.getNeighborhood8(3, 2)).toEqual(
            image.getNeighborhood8().map((offset) => index + offset));
        expect(image.getCorners(3, 2)).toEqual(
            image.getCorners().map((offset) => index + offset));
        expect(image.getFull(3, 2)).toEqual(
            image.getFull().map((offset) => index + offset));
    });

    it('computes relative 2-tuple neighborhoods', () => {
        const image = new Image2<number>(7, 5);
        expect(image.getNeighborhood4Coords()).toEqual([
            [-1, 0], [1, 0], [0, -1], [0, 1]]);
        // The 2-tuple 8-neighborhood is in row-major order without the
        // center, unlike the 1-dimensional ordering (upstream behavior).
        expect(image.getNeighborhood8Coords()).toEqual([
            [-1, -1], [0, -1], [1, -1], [-1, 0], [1, 0], [-1, 1], [0, 1], [1, 1]]);
        expect(image.getCornersCoords()).toEqual([
            [0, 0], [1, 0], [0, 1], [1, 1]]);
        expect(image.getFullCoords()).toEqual([
            [-1, -1], [0, -1], [1, -1], [-1, 0], [0, 0], [1, 0], [-1, 1], [0, 1], [1, 1]]);
    });

    it('computes absolute 2-tuple neighborhoods at (x,y)', () => {
        const image = new Image2<number>(7, 5);
        expect(image.getNeighborhood4Coords(3, 2)).toEqual([
            [2, 2], [4, 2], [3, 1], [3, 3]]);
        expect(image.getNeighborhood8Coords(3, 2)).toEqual([
            [2, 1], [3, 1], [4, 1], [2, 2], [4, 2], [2, 3], [3, 3], [4, 3]]);
        expect(image.getCornersCoords(3, 2)).toEqual([
            [3, 2], [4, 2], [3, 3], [4, 3]]);
        expect(image.getFullCoords(3, 2)).toEqual([
            [2, 1], [3, 1], [4, 1], [2, 2], [3, 2], [4, 2], [2, 3], [3, 3], [4, 3]]);
    });

    it('has consistent 1-dimensional and 2-tuple neighborhoods', () => {
        const image = new Image2<number>(6, 4);
        const x = 2, y = 1;
        const nbr4 = image.getNeighborhood4(x, y);
        const nbr8 = image.getNeighborhood8(x, y);
        const corners = image.getCorners(x, y);
        const full = image.getFull(x, y);

        // Each 2-tuple neighborhood maps to the same set of 1-dimensional
        // indices (orderings differ for the 8-neighborhood).
        const toIndex = (coord: [number, number]): number => image.getIndex(coord);
        expect(image.getNeighborhood4Coords(x, y).map(toIndex)).toEqual(nbr4);
        expect(image.getNeighborhood8Coords(x, y).map(toIndex).sort((a, b) => a - b))
            .toEqual([...nbr8].sort((a, b) => a - b));
        expect(image.getCornersCoords(x, y).map(toIndex)).toEqual(corners);
        expect(image.getFullCoords(x, y).map(toIndex)).toEqual(full);
    });
});
