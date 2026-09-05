import { describe, it, expect } from 'vitest';
import { Image2 } from '../src/Image2.js';
import { check, fc } from './helpers/arbitraries.js';
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

// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream Image2.h semantics.
// ---------------------------------------------------------------------------

describe('Image2 verification', () => {
    const bounds = fc.tuple(fc.integer({ min: 1, max: 9 }), fc.integer({ min: 1, max: 9 }));

    function filled(dim0: number, dim1: number): Image2<number> {
        const image = new Image2<number>(dim0, dim1);
        for (let i = 0; i < image.getNumPixels(); ++i) {
            image.set(i, 1000 + i);
        }
        return image;
    }

    it('round trips every (x, y) through getIndex/getCoordinates', () => {
        check(bounds, ([dim0, dim1]) => {
            const image = new Image2<number>(dim0, dim1);
            for (let y = 0; y < dim1; ++y) {
                for (let x = 0; x < dim0; ++x) {
                    const index = image.getIndex(x, y);
                    expect(index).toBe(x + dim0 * y);
                    expect(image.getCoordinates(index)).toEqual([x, y]);
                    // The array overload and the inherited base-class
                    // linearization must agree with the (x, y) overload.
                    expect(image.getIndex([x, y])).toBe(index);
                }
            }
            for (let i = 0; i < image.getNumPixels(); ++i) {
                const [x, y] = image.getCoordinates(i);
                expect(image.getIndex(x, y)).toBe(i);
            }
        });
    });

    it('the three get/set overloads address the same pixel', () => {
        check(bounds, ([dim0, dim1]) => {
            const image = filled(dim0, dim1);
            for (let y = 0; y < dim1; ++y) {
                for (let x = 0; x < dim0; ++x) {
                    const i = x + dim0 * y;
                    expect(image.get(x, y)).toBe(image.get(i));
                    expect(image.get([x, y])).toBe(image.get(i));
                    image.set(x, y, -i);
                    expect(image.get(i)).toBe(-i);
                    image.set([x, y], i);
                    expect(image.get(i)).toBe(i);
                    image.set(i, 2 * i);
                    expect(image.get(x, y)).toBe(2 * i);
                }
            }
        });
    });

    it('clamped (x, y) accessors clamp each coordinate independently', () => {
        check(fc.tuple(bounds, fc.integer({ min: -4, max: 12 }), fc.integer({ min: -4, max: 12 })),
            ([[dim0, dim1], x, y]) => {
                const image = filled(dim0, dim1);
                const cx = Math.min(Math.max(x, 0), dim0 - 1);
                const cy = Math.min(Math.max(y, 0), dim1 - 1);
                const expected = 1000 + cx + dim0 * cy;
                expect(image.getClamped(x, y)).toBe(expected);
                expect(image.getClamped([x, y])).toBe(expected);

                image.setClamped(x, y, -1);
                expect(image.get(cx, cy)).toBe(-1);
                image.setClamped([x, y], -2);
                expect(image.get(cx, cy)).toBe(-2);
            });
    });

    it('relative neighborhood offsets equal the 2-tuple offsets linearized', () => {
        check(bounds, ([dim0]) => {
            const image = new Image2<number>(dim0, 5);
            const linearize = (c: [number, number]) => c[0] + dim0 * c[1];
            // The 4-neighborhood, corners and full neighborhood use the same
            // ordering in both forms; the 8-neighborhood does not (upstream
            // orders the 1-dimensional form with the 4-connected neighbors
            // first and the 2-tuple form in row-major order), so only the
            // sets are compared there.
            expect(image.getNeighborhood4Coords().map(linearize))
                .toEqual(image.getNeighborhood4());
            expect(image.getCornersCoords().map(linearize)).toEqual(image.getCorners());
            expect(image.getFullCoords().map(linearize)).toEqual(image.getFull());

            const sorted = (a: number[]) => [...a].sort((p, q) => p - q);
            expect(sorted(image.getNeighborhood8Coords().map(linearize)))
                .toEqual(sorted(image.getNeighborhood8()));
        });
    });

    it('absolute neighborhoods are the relative offsets added to getIndex(x, y)', () => {
        check(fc.tuple(bounds, fc.nat({ max: 8 }), fc.nat({ max: 8 })),
            ([[dim0, dim1], rx, ry]) => {
                const image = new Image2<number>(dim0, dim1);
                const x = rx % dim0, y = ry % dim1;
                const index = image.getIndex(x, y);
                const shift = (a: number[]) => a.map(v => v + index);
                expect(image.getNeighborhood4(x, y)).toEqual(shift(image.getNeighborhood4()));
                expect(image.getNeighborhood8(x, y)).toEqual(shift(image.getNeighborhood8()));
                expect(image.getCorners(x, y)).toEqual(shift(image.getCorners()));
                expect(image.getFull(x, y)).toEqual(shift(image.getFull()));

                const offset = (a: [number, number][]) =>
                    a.map(c => [c[0] + x, c[1] + y]);
                expect(image.getNeighborhood4Coords(x, y))
                    .toEqual(offset(image.getNeighborhood4Coords()));
                expect(image.getNeighborhood8Coords(x, y))
                    .toEqual(offset(image.getNeighborhood8Coords()));
                expect(image.getCornersCoords(x, y))
                    .toEqual(offset(image.getCornersCoords()));
                expect(image.getFullCoords(x, y))
                    .toEqual(offset(image.getFullCoords()));
            });
    });

    it('2-tuple neighborhoods stay signed on the boundary (upstream #64)', () => {
        // Upstream computes size_t(x) + inbr[i][0], which wraps to SIZE_MAX
        // for x == 0; the port keeps signed numbers so the documented
        // no-clamping behavior yields negative coordinates instead.
        check(bounds, ([dim0, dim1]) => {
            const image = new Image2<number>(dim0, dim1);
            expect(image.getFullCoords(0, 0)[0]).toEqual([-1, -1]);
            expect(image.getNeighborhood4Coords(0, 0)[0]).toEqual([-1, 0]);
            expect(image.getNeighborhood8Coords(0, 0)[0]).toEqual([-1, -1]);
            expect(image.getNeighborhood4(0, 0)[0]).toBe(-1);
            expect(image.getFull(0, 0)[0]).toBe(-1 - dim0);
        });
    });

    it('the full neighborhood is the 8-neighborhood plus the center', () => {
        check(bounds, ([dim0]) => {
            const image = new Image2<number>(dim0, 4);
            const sorted = (a: number[]) => [...a].sort((p, q) => p - q);
            expect(sorted(image.getFull()))
                .toEqual(sorted([...image.getNeighborhood8(), 0]));
            // Every offset in the full neighborhood is a distinct member of
            // the 3x3 block centered at the pixel, in row-major order.
            const brute: number[] = [];
            for (let dy = -1; dy <= 1; ++dy) {
                for (let dx = -1; dx <= 1; ++dx) {
                    brute.push(dx + dim0 * dy);
                }
            }
            expect(image.getFull()).toEqual(brute);
        });
    });

    it('reconstruct(dim0, dim1) and reconstruct([dim0, dim1]) agree', () => {
        check(fc.tuple(bounds, bounds), ([[a0, a1], [b0, b1]]) => {
            const first = new Image2<number>(a0, a1);
            const second = new Image2<number>(a0, a1);
            first.reconstruct(b0, b1);
            second.reconstruct([b0, b1]);
            expect([...first.getDimensions()]).toEqual([...second.getDimensions()]);
            expect(first.getNumPixels()).toBe(second.getNumPixels());
            expect([...first.getOffsets()]).toEqual([1, b0]);
        });
    });
});
