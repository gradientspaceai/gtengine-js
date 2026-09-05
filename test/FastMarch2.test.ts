import { describe, it, expect } from 'vitest';
import { FastMarch2 } from '../src/FastMarch2.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

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

describe('FastMarch2 verification', () => {
    const MAXREAL = Number.MAX_VALUE;

    // Bounds of at least 3 (so there is an interior) with an interior seed
    // derived from the bounds, which avoids a filtered generator.
    const grid = fc.tuple(
        fc.integer({ min: 3, max: 9 }),
        fc.integer({ min: 3, max: 9 }),
        fc.nat({ max: 1000 }), fc.nat({ max: 1000 }),
        fc.constantFrom(0.5, 1, 2, 4))
        .map(([xB, yB, rx, ry, speed]) => ({
            xB, yB, speed,
            sx: 1 + rx % (xB - 2),
            sy: 1 + ry % (yB - 2)
        }));

    it('marks the whole border zero speed and keeps the seed at time zero', () => {
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const seed = sx + xB * sy;
            const filter = new FastMarch2(xB, yB, 1, 1, [seed], speed);
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    const i = x + xB * y;
                    const onBorder = (x === 0 || x === xB - 1
                        || y === 0 || y === yB - 1);
                    expect(filter.isZeroSpeed(i)).toBe(onBorder);
                    if (onBorder) {
                        expect(filter.isValid(i)).toBe(false);
                        expect(filter.isFar(i)).toBe(false);
                    }
                }
            }
            expect(filter.getTime(seed)).toBe(0);
            expect(filter.isValid(seed)).toBe(true);
        });
    });

    it('accepts the trial pixels in nondecreasing time order', () => {
        // The defining property of the fast marching method: the front never
        // moves backwards. A heap-bookkeeping error shows up here.
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const seed = sx + xB * sy;
            const filter = new FastMarch2(xB, yB, 1, 1, [seed], speed);
            let previous = 0;
            for (let n = 0; n < filter.getQuantity(); ++n) {
                const trialsBefore: number[] = [];
                for (let i = 0; i < filter.getQuantity(); ++i) {
                    if (filter.isTrial(i)) {
                        trialsBefore.push(i);
                    }
                }
                if (trialsBefore.length === 0) {
                    break;
                }
                filter.iterate();
                // Exactly the pixels that stopped being trials this step.
                const accepted = trialsBefore.filter(i => !filter.isTrial(i));
                expect(accepted.length).toBe(1);
                const t = filter.getTime(accepted[0]);
                expect(t).toBeGreaterThanOrEqual(previous - 1e-12);
                previous = t;
            }
        });
    });

    it('leaves every interior pixel known and the border untouched', () => {
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const seed = sx + xB * sy;
            const filter = new FastMarch2(xB, yB, 1, 1, [seed], speed);
            march(filter);
            for (let y = 0; y < yB; ++y) {
                for (let x = 0; x < xB; ++x) {
                    const i = x + xB * y;
                    const interior = (0 < x && x < xB - 1 && 0 < y && y < yB - 1);
                    if (interior) {
                        expect(filter.isValid(i)).toBe(true);
                        expect(filter.isTrial(i)).toBe(false);
                        expect(filter.isInterior(i)).toBe(true);
                        expect(Number.isFinite(filter.getTime(i))).toBe(true);
                        expect(filter.getTime(i)).toBeGreaterThanOrEqual(0);
                    } else {
                        expect(filter.getTime(i)).toBe(-MAXREAL);
                    }
                }
            }
            // With no trials left nothing is a front boundary.
            expect(filter.getBoundary()).toEqual([]);
            expect(filter.getInterior().length)
                .toBe((xB - 2) * (yB - 2));
        });
    });

    it('scales the crossing times by the reciprocal of the speed', () => {
        // The update formula is homogeneous of degree one in 1/speed, so the
        // whole solution scales exactly.
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const seed = sx + xB * sy;
            const unit = new FastMarch2(xB, yB, 1, 1, [seed], 1);
            const scaled = new FastMarch2(xB, yB, 1, 1, [seed], speed);
            march(unit);
            march(scaled);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    const i = x + xB * y;
                    expectClose(scaled.getTime(i), unit.getTime(i) / speed,
                        1e-9, 1e-9);
                }
            }
        });
    });

    it('is equivariant under the grid symmetries', () => {
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const base = new FastMarch2(xB, yB, 1, 1, [sx + xB * sy], speed);
            march(base);

            // Reflect x.
            const mirrored = new FastMarch2(xB, yB, 1, 1,
                [(xB - 1 - sx) + xB * sy], speed);
            march(mirrored);
            // Transpose.
            const transposed = new FastMarch2(yB, xB, 1, 1, [sy + yB * sx], speed);
            march(transposed);

            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    expectClose(mirrored.getTime((xB - 1 - x) + xB * y),
                        base.getTime(x + xB * y), 1e-9, 1e-9);
                    expectClose(transposed.getTime(y + yB * x),
                        base.getTime(x + xB * y), 1e-9, 1e-9);
                }
            }
        });
    });

    it('brackets the arrival time between the Euclidean and L1 distances', () => {
        // The first-order upwind scheme on a unit grid overestimates the
        // Euclidean distance and never exceeds the grid (L1) distance: the
        // two-term update satisfies t <= min(a,b) + 1/speed, which is the L1
        // recursion.
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const filter = new FastMarch2(xB, yB, 1, 1, [sx + xB * sy], speed);
            march(filter);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    const dx = Math.abs(x - sx), dy = Math.abs(y - sy);
                    const t = filter.getTime(x + xB * y);
                    expect(t).toBeGreaterThanOrEqual(
                        Math.sqrt(dx * dx + dy * dy) / speed - 1e-9);
                    expect(t).toBeLessThanOrEqual((dx + dy) / speed + 1e-9);
                }
            }
        });
    });

    it('satisfies the local upwind bounds against its 4-neighbors', () => {
        // Every update is either invSpeed plus one neighbor time or a
        // two-term root that lies between the smallest neighbor time and that
        // time plus invSpeed.
        check(grid, ({ xB, yB, sx, sy, speed }) => {
            const seed = sx + xB * sy;
            const filter = new FastMarch2(xB, yB, 1, 1, [seed], speed);
            march(filter);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    const i = x + xB * y;
                    if (i === seed) {
                        continue;
                    }
                    const neighbors = [i - 1, i + 1, i - xB, i + xB]
                        .filter(j => filter.isValid(j))
                        .map(j => filter.getTime(j));
                    expect(neighbors.length).toBeGreaterThan(0);
                    const smallest = Math.min(...neighbors);
                    const t = filter.getTime(i);
                    expect(t).toBeGreaterThanOrEqual(smallest - 1e-9);
                    expect(t).toBeLessThanOrEqual(smallest + 1 / speed + 1e-9);
                }
            }
        });
    });
    it('reaches every interior pixel from any seed set', () => {
        // Adding seeds does NOT lower every arrival time: see the
        // upstream-quirk test below. What does hold is that every seed keeps
        // time 0 and the march still classifies the whole interior.
        const twoSeeds = fc.tuple(grid, fc.nat({ max: 1000 }), fc.nat({ max: 1000 }))
            .map(([g, rx, ry]) => ({
                ...g,
                tx: 1 + rx % (g.xB - 2),
                ty: 1 + ry % (g.yB - 2)
            }));
        check(twoSeeds, ({ xB, yB, sx, sy, tx, ty, speed }) => {
            const a = sx + xB * sy;
            const b = tx + xB * ty;
            const both = new FastMarch2(xB, yB, 1, 1, [a, b], speed);
            expect(both.getTime(a)).toBe(0);
            expect(both.getTime(b)).toBe(0);
            march(both);
            expect(both.getTime(a)).toBe(0);
            expect(both.getTime(b)).toBe(0);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    const i = x + xB * y;
                    expect(both.isInterior(i)).toBe(true);
                    expect(both.getTime(i)).toBeGreaterThanOrEqual(0);
                    expect(Number.isFinite(both.getTime(i))).toBe(true);
                }
            }
        });
    });

    it('upstream quirk: the negative-discriminant fallback raises times', () => {
        // ComputeTime falls back to "the maximum time of the neighbors" when
        // the two-term quadratic has no real root (upstream even asks "Is
        // there a better choice?"). The Godunov scheme would fall back to the
        // one-sided update min(xConst, yConst) + 1/speed instead, so the
        // upstream choice can RAISE a trial time above a one-term value that
        // was already computed. The quirk is preserved by the port; this test
        // pins it so a future change is deliberate.
        //
        // 6x4 grid, speed 1/2 (so 1/speed = 2), seeds at (4,2) and (1,2).
        const xB = 6, yB = 4, speed = 0.5;
        const a = 4 + xB * 2, b = 1 + xB * 2;
        const both = new FastMarch2(xB, yB, 1, 1, [a, b], speed);
        // The pixel (2,2) is a 4-neighbor of the seed (1,2), so its arrival
        // time starts (correctly) at 1/speed = 2.
        expect(both.getTime(2 + xB * 2)).toBeCloseTo(2, 12);
        march(both);
        // Accepting (3,2) recomputes (2,2) with xConst = 0 and yConst = 1 +
        // sqrt(2) times 1/speed; the discriminant 2/speed^2 - (xConst-yConst)^2
        // is negative, so the fallback takes the larger of the two and raises
        // (2,2) from 2 to 1 + sqrt(2) times 1/speed.
        expect(both.getTime(2 + xB * 2)).toBeCloseTo(2 * (1 + Math.SQRT1_2), 12);
        // A single seed at (1,2) gives the expected 1/speed there.
        const one = new FastMarch2(xB, yB, 1, 1, [b], speed);
        march(one);
        expect(one.getTime(2 + xB * 2)).toBeCloseTo(2, 12);
    });
});
