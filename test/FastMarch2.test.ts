import { describe, it, expect } from 'vitest';
import { FastMarch2 } from '../src/FastMarch2.js';

// Run the fast marching to completion (the heap empties, after which
// iterate() is a no-op).
function march(filter: FastMarch2): void {
    for (let i = 0; i < filter.getQuantity(); ++i) {
        filter.iterate();
    }
}

describe('FastMarch2', () => {
    it('exposes the bounds, spacings and lexicographic index', () => {
        const filter = new FastMarch2(5, 4, 0.5, 2, [12], 1);
        expect(filter.getXBound()).toBe(5);
        expect(filter.getYBound()).toBe(4);
        expect(filter.getXSpacing()).toBe(0.5);
        expect(filter.getYSpacing()).toBe(2);
        expect(filter.getQuantity()).toBe(20);
        expect(filter.index(0, 0)).toBe(0);
        expect(filter.index(2, 3)).toBe(2 + 5 * 3);
    });

    it('marks the whole image border as zero speed', () => {
        const bound = 5;
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(2, 2, bound)], 1);
        for (let y = 0; y < bound; ++y) {
            for (let x = 0; x < bound; ++x) {
                const i = filter2Index(x, y, bound);
                const onBorder = (x === 0 || y === 0 || x === bound - 1 || y === bound - 1);
                expect(filter.isZeroSpeed(i)).toBe(onBorder);
                if (onBorder) {
                    expect(filter.getTime(i)).toBe(-Number.MAX_VALUE);
                    expect(filter.isValid(i)).toBe(false);
                }
            }
        }
    });

    it('gives the seed time zero and its axis neighbors the reciprocal speed', () => {
        const bound = 7, c = 3;
        const seed = filter2Index(c, c, bound);
        const filter = new FastMarch2(bound, bound, 1, 1, [seed], 1);
        expect(filter.getTime(seed)).toBe(0);

        // Before any iteration the only trial pixels are the four 4-neighbors
        // of the seed, whose equation is linear: t = 1/speed.
        for (const [x, y] of [[c + 1, c], [c - 1, c], [c, c + 1], [c, c - 1]]) {
            const i = filter2Index(x, y, bound);
            expect(filter.isTrial(i)).toBe(true);
            expect(filter.getTime(i)).toBeCloseTo(1, 15);
        }
        // The diagonal neighbors have no *known* (valid, non-trial) neighbor
        // yet, so they are still far.
        for (const [x, y] of [[c + 1, c + 1], [c - 1, c + 1], [c + 1, c - 1],
            [c - 1, c - 1]]) {
            const i = filter2Index(x, y, bound);
            expect(filter.isTrial(i)).toBe(false);
            expect(filter.isFar(i)).toBe(true);
        }

        // With speed 2 the reciprocal speed is 1/2.
        const fast = new FastMarch2(bound, bound, 1, 1, [seed], 2);
        expect(fast.getTime(filter2Index(c + 1, c, bound))).toBeCloseTo(0.5, 15);
    });

    it('halves the crossing times when the speed doubles', () => {
        const bound = 9, c = 4;
        const seed = filter2Index(c, c, bound);
        const unit = new FastMarch2(bound, bound, 1, 1, [seed], 1);
        const fast = new FastMarch2(bound, bound, 1, 1, [seed], 2);
        march(unit);
        march(fast);
        for (let y = 1; y < bound - 1; ++y) {
            for (let x = 1; x < bound - 1; ++x) {
                const i = filter2Index(x, y, bound);
                expect(fast.getTime(i)).toBeCloseTo(0.5 * unit.getTime(i), 12);
            }
        }
    });

    it('propagates a point source to approximately the Euclidean distance', () => {
        const bound = 15, c = 7;
        const seed = filter2Index(c, c, bound);
        const filter = new FastMarch2(bound, bound, 1, 1, [seed], 1);
        march(filter);

        for (let y = 1; y < bound - 1; ++y) {
            for (let x = 1; x < bound - 1; ++x) {
                const i = filter2Index(x, y, bound);
                const distance = Math.hypot(x - c, y - c);
                expect(filter.isValid(i)).toBe(true);
                const time = filter.getTime(i);
                // The first-order upwind scheme never underestimates the
                // distance and overestimates it by at most ~21% (worst at the
                // diagonal immediately next to the source).
                expect(time).toBeGreaterThanOrEqual(distance - 1e-12);
                expect(time).toBeLessThanOrEqual(1.21 * distance + 1e-12);
                // Times along the axes are exact.
                if (y === c) {
                    expect(time).toBeCloseTo(distance, 10);
                }
            }
        }

        // The (1,1) diagonal neighbor has the closed-form value
        // 0.5*(sum + sqrt(2*invSpeed^2 - diff^2)) with sum = 2 and diff = 0.
        expect(filter.getTime(filter2Index(c + 1, c + 1, bound)))
            .toBeCloseTo(0.5 * (2 + Math.SQRT2), 12);
    });

    it('produces times that increase away from the seed', () => {
        const bound = 11, c = 5;
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(c, c, bound)], 1);
        march(filter);
        for (let x = c; x < bound - 2; ++x) {
            expect(filter.getTime(filter2Index(x + 1, c, bound)))
                .toBeGreaterThan(filter.getTime(filter2Index(x, c, bound)));
        }
        // The image is symmetric under the reflections that fix the seed.
        for (let y = 1; y < bound - 1; ++y) {
            for (let x = 1; x < bound - 1; ++x) {
                expect(filter.getTime(filter2Index(x, y, bound))).toBeCloseTo(
                    filter.getTime(filter2Index(bound - 1 - x, y, bound)), 12);
                expect(filter.getTime(filter2Index(x, y, bound))).toBeCloseTo(
                    filter.getTime(filter2Index(y, x, bound)), 12);
            }
        }
    });

    it('classifies interior, trial and boundary pixels during the march', () => {
        const bound = 9, c = 4;
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(c, c, bound)], 1);

        // The seed is the only known pixel initially, and it is a boundary
        // pixel because its neighbors are trial pixels.
        expect(filter.getInterior()).toEqual([filter2Index(c, c, bound)]);
        expect(filter.isBoundary(filter2Index(c, c, bound))).toBe(true);
        expect(filter.getBoundary()).toEqual([filter2Index(c, c, bound)]);

        for (let i = 0; i < 8; ++i) {
            filter.iterate();
        }
        const interior = filter.getInterior();
        expect(interior.length).toBe(9);
        for (const i of interior) {
            expect(filter.isInterior(i)).toBe(true);
            expect(filter.isTrial(i)).toBe(false);
        }
        // Every boundary pixel is interior and has a trial neighbor.
        const boundary = filter.getBoundary();
        expect(boundary.length).toBeGreaterThan(0);
        for (const i of boundary) {
            expect(interior).toContain(i);
        }

        march(filter);
        // Once the heap is empty nothing is a trial pixel any more.
        expect(filter.getBoundary()).toEqual([]);
    });

    it('never crosses a zero-speed barrier', () => {
        const bound = 11, c = 5;
        const speeds = new Array<number>(bound * bound).fill(1);
        // A wall at x = 5 with a gap at y = 1 (which lies on the image
        // border, so it is closed too): the right half is unreachable.
        for (let y = 0; y < bound; ++y) {
            speeds[filter2Index(5, y, bound)] = 0;
        }
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(2, c, bound)], speeds);
        march(filter);

        for (let y = 1; y < bound - 1; ++y) {
            expect(filter.isZeroSpeed(filter2Index(5, y, bound))).toBe(true);
            for (let x = 6; x < bound - 1; ++x) {
                const i = filter2Index(x, y, bound);
                expect(filter.isFar(i)).toBe(true);
                expect(filter.isValid(i)).toBe(false);
            }
            for (let x = 1; x < 5; ++x) {
                expect(filter.isValid(filter2Index(x, y, bound))).toBe(true);
            }
        }
    });

    it('reports the extremes of the valid crossing times', () => {
        const bound = 9, c = 4;
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(c, c, bound)], 1);
        march(filter);
        const { minValue, maxValue } = filter.getTimeExtremes();
        expect(minValue).toBe(0);
        // The far corner of the interior is the last pixel reached.
        expect(maxValue).toBeCloseTo(filter.getTime(filter2Index(1, 1, bound)), 12);
        expect(maxValue).toBeGreaterThan(Math.hypot(3, 3));
    });

    it('accepts multiple seeds and takes the minimum arrival time', () => {
        const bound = 11;
        const left = filter2Index(2, 5, bound);
        const right = filter2Index(8, 5, bound);
        const both = new FastMarch2(bound, bound, 1, 1, [left, right], 1);
        const single = new FastMarch2(bound, bound, 1, 1, [left], 1);
        march(both);
        march(single);
        expect(both.getTime(right)).toBe(0);
        for (let y = 1; y < bound - 1; ++y) {
            for (let x = 1; x < bound - 1; ++x) {
                const i = filter2Index(x, y, bound);
                expect(both.getTime(i)).toBeLessThanOrEqual(single.getTime(i) + 1e-12);
            }
        }
        // Times are symmetric about the midpoint between the two seeds.
        for (let y = 1; y < bound - 1; ++y) {
            for (let x = 1; x < bound - 1; ++x) {
                expect(both.getTime(filter2Index(x, y, bound))).toBeCloseTo(
                    both.getTime(filter2Index(bound - 1 - x, y, bound)), 12);
            }
        }
    });

    it('iterating an exhausted heap is a no-op', () => {
        const bound = 7, c = 3;
        const filter = new FastMarch2(bound, bound, 1, 1, [filter2Index(c, c, bound)], 1);
        march(filter);
        const before = filter.getInterior();
        filter.iterate();
        filter.iterate();
        expect(filter.getInterior()).toEqual(before);
    });
});

function filter2Index(x: number, y: number, xBound: number): number {
    return x + xBound * y;
}
