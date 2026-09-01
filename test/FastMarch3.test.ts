import { describe, it, expect } from 'vitest';
import { FastMarch3 } from '../src/FastMarch3';

function idx(x: number, y: number, z: number, xBound: number, yBound: number): number {
    return x + xBound * (y + yBound * z);
}

// Run the fast marching to completion (the heap empties, after which
// iterate() is a no-op).
function march(filter: FastMarch3): void {
    for (let i = 0; i < filter.getQuantity(); ++i) {
        filter.iterate();
    }
}

describe('FastMarch3', () => {
    it('exposes the bounds, spacings and lexicographic index', () => {
        const filter = new FastMarch3(5, 4, 3, 0.5, 2, 4, [30], 1);
        expect(filter.getXBound()).toBe(5);
        expect(filter.getYBound()).toBe(4);
        expect(filter.getZBound()).toBe(3);
        expect(filter.getXSpacing()).toBe(0.5);
        expect(filter.getYSpacing()).toBe(2);
        expect(filter.getZSpacing()).toBe(4);
        expect(filter.getQuantity()).toBe(60);
        expect(filter.index(0, 0, 0)).toBe(0);
        expect(filter.index(2, 3, 1)).toBe(2 + 5 * (3 + 4 * 1));
    });

    it('marks the six boundary faces as zero speed', () => {
        // Upstream marks only the box vertices and edges, which leaves the
        // face interiors reachable by the front and makes their off-grid
        // neighbor indices wrap. The port marks whole faces, as the 2D code
        // does for its whole border.
        const b = 5;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(2, 2, 2, b, b)], 1);
        for (let z = 0; z < b; ++z) {
            for (let y = 0; y < b; ++y) {
                for (let x = 0; x < b; ++x) {
                    const i = idx(x, y, z, b, b);
                    const onFace = (x === 0 || y === 0 || z === 0
                        || x === b - 1 || y === b - 1 || z === b - 1);
                    expect(filter.isZeroSpeed(i)).toBe(onFace);
                    if (onFace) {
                        expect(filter.getTime(i)).toBe(-Number.MAX_VALUE);
                        expect(filter.isValid(i)).toBe(false);
                        expect(filter.isTrial(i)).toBe(false);
                    }
                }
            }
        }
    });

    it('never lets the front reach a boundary face', () => {
        const b = 9, c = 4;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(c, c, c, b, b)], 1);
        march(filter);
        for (let y = 0; y < b; ++y) {
            for (let x = 0; x < b; ++x) {
                for (const face of [
                    idx(x, y, 0, b, b), idx(x, y, b - 1, b, b),
                    idx(x, 0, y, b, b), idx(x, b - 1, y, b, b),
                    idx(0, x, y, b, b), idx(b - 1, x, y, b, b)]) {
                    expect(filter.isValid(face)).toBe(false);
                    expect(filter.isTrial(face)).toBe(false);
                }
            }
        }
        // Every interior voxel is finalized and carries a finite time.
        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                for (let x = 1; x < b - 1; ++x) {
                    const i = idx(x, y, z, b, b);
                    expect(filter.isValid(i)).toBe(true);
                    expect(Number.isFinite(filter.getTime(i))).toBe(true);
                }
            }
        }
    });

    it('gives the seed time zero and its face neighbors the reciprocal speed', () => {
        const b = 7, c = 3;
        const seed = idx(c, c, c, b, b);
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [seed], 1);
        expect(filter.getTime(seed)).toBe(0);
        for (const [x, y, z] of [[c + 1, c, c], [c - 1, c, c], [c, c + 1, c],
            [c, c - 1, c], [c, c, c + 1], [c, c, c - 1]]) {
            const i = idx(x, y, z, b, b);
            expect(filter.isTrial(i)).toBe(true);
            expect(filter.getTime(i)).toBeCloseTo(1, 15);
        }
        // The edge and corner neighbors have no known neighbor yet.
        expect(filter.isFar(idx(c + 1, c + 1, c, b, b))).toBe(true);
        expect(filter.isFar(idx(c + 1, c + 1, c + 1, b, b))).toBe(true);

        const fast = new FastMarch3(b, b, b, 1, 1, 1, [seed], 4);
        expect(fast.getTime(idx(c + 1, c, c, b, b))).toBeCloseTo(0.25, 15);
    });

    it('propagates a point source to approximately the Euclidean distance', () => {
        const b = 13, c = 6;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(c, c, c, b, b)], 1);
        march(filter);

        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                for (let x = 1; x < b - 1; ++x) {
                    const i = idx(x, y, z, b, b);
                    const distance = Math.sqrt((x - c) ** 2 + (y - c) ** 2 + (z - c) ** 2);
                    const time = filter.getTime(i);
                    // The first-order upwind scheme never underestimates the
                    // distance and overestimates it by at most ~32% (worst at
                    // the body diagonal next to the source).
                    expect(time).toBeGreaterThanOrEqual(distance - 1e-12);
                    expect(time).toBeLessThanOrEqual(1.32 * distance + 1e-12);
                    if (y === c && z === c) {
                        // Times along the axes are exact.
                        expect(time).toBeCloseTo(distance, 10);
                    }
                }
            }
        }

        // Closed forms for the first two off-axis neighbors: a face diagonal
        // solves the 2-term quadratic with equal constants 1, and the body
        // diagonal solves the 3-term quadratic with equal constants
        // (2 + sqrt(2))/2.
        const face = 0.5 * (2 + Math.SQRT2);
        expect(filter.getTime(idx(c + 1, c + 1, c, b, b))).toBeCloseTo(face, 12);
        expect(filter.getTime(idx(c + 1, c + 1, c + 1, b, b)))
            .toBeCloseTo((3 * face + Math.sqrt(3)) / 3, 12);
    });

    it('halves the crossing times when the speed doubles', () => {
        const b = 9, c = 4;
        const seed = idx(c, c, c, b, b);
        const unit = new FastMarch3(b, b, b, 1, 1, 1, [seed], 1);
        const fast = new FastMarch3(b, b, b, 1, 1, 1, [seed], 2);
        march(unit);
        march(fast);
        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                for (let x = 1; x < b - 1; ++x) {
                    const i = idx(x, y, z, b, b);
                    expect(fast.getTime(i)).toBeCloseTo(0.5 * unit.getTime(i), 12);
                }
            }
        }
    });

    it('produces times that increase away from the seed and are symmetric', () => {
        const b = 9, c = 4;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(c, c, c, b, b)], 1);
        march(filter);
        for (let x = c; x < b - 2; ++x) {
            expect(filter.getTime(idx(x + 1, c, c, b, b)))
                .toBeGreaterThan(filter.getTime(idx(x, c, c, b, b)));
        }
        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                for (let x = 1; x < b - 1; ++x) {
                    const t = filter.getTime(idx(x, y, z, b, b));
                    expect(t).toBeCloseTo(filter.getTime(idx(b - 1 - x, y, z, b, b)), 12);
                    expect(t).toBeCloseTo(filter.getTime(idx(y, x, z, b, b)), 12);
                    expect(t).toBeCloseTo(filter.getTime(idx(x, z, y, b, b)), 12);
                }
            }
        }
    });

    it('classifies interior, trial and boundary voxels during the march', () => {
        const b = 9, c = 4;
        const seed = idx(c, c, c, b, b);
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [seed], 1);

        expect(filter.getInterior()).toEqual([seed]);
        expect(filter.isBoundary(seed)).toBe(true);
        expect(filter.getBoundary()).toEqual([seed]);

        for (let i = 0; i < 6; ++i) {
            filter.iterate();
        }
        const interior = filter.getInterior();
        expect(interior.length).toBe(7);
        for (const i of interior) {
            expect(filter.isInterior(i)).toBe(true);
        }
        for (const i of filter.getBoundary()) {
            expect(interior).toContain(i);
        }

        march(filter);
        expect(filter.getBoundary()).toEqual([]);
    });

    it('never crosses a zero-speed barrier', () => {
        const b = 9;
        const speeds = new Array<number>(b * b * b).fill(1);
        for (let z = 0; z < b; ++z) {
            for (let y = 0; y < b; ++y) {
                speeds[idx(4, y, z, b, b)] = 0;
            }
        }
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(2, 4, 4, b, b)], speeds);
        march(filter);
        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                expect(filter.isZeroSpeed(idx(4, y, z, b, b))).toBe(true);
                for (let x = 5; x < b - 1; ++x) {
                    expect(filter.isFar(idx(x, y, z, b, b))).toBe(true);
                }
                for (let x = 1; x < 4; ++x) {
                    expect(filter.isValid(idx(x, y, z, b, b))).toBe(true);
                }
            }
        }
    });

    it('accepts multiple seeds and takes the minimum arrival time', () => {
        const b = 9;
        const left = idx(2, 4, 4, b, b);
        const right = idx(6, 4, 4, b, b);
        const both = new FastMarch3(b, b, b, 1, 1, 1, [left, right], 1);
        const single = new FastMarch3(b, b, b, 1, 1, 1, [left], 1);
        march(both);
        march(single);
        expect(both.getTime(right)).toBe(0);
        for (let z = 1; z < b - 1; ++z) {
            for (let y = 1; y < b - 1; ++y) {
                for (let x = 1; x < b - 1; ++x) {
                    const i = idx(x, y, z, b, b);
                    expect(both.getTime(i)).toBeLessThanOrEqual(single.getTime(i) + 1e-12);
                    expect(both.getTime(i)).toBeCloseTo(
                        both.getTime(idx(b - 1 - x, y, z, b, b)), 12);
                }
            }
        }
    });

    it('reports the extremes of the valid crossing times', () => {
        const b = 9, c = 4;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(c, c, c, b, b)], 1);
        march(filter);
        const { minValue, maxValue } = filter.getTimeExtremes();
        expect(minValue).toBe(0);
        expect(maxValue).toBeCloseTo(filter.getTime(idx(1, 1, 1, b, b)), 12);
        expect(maxValue).toBeGreaterThan(Math.sqrt(27));
    });

    it('iterating an exhausted heap is a no-op', () => {
        const b = 7, c = 3;
        const filter = new FastMarch3(b, b, b, 1, 1, 1, [idx(c, c, c, b, b)], 1);
        march(filter);
        const before = filter.getInterior();
        filter.iterate();
        expect(filter.getInterior()).toEqual(before);
    });
});
