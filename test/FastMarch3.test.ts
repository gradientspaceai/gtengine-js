import { describe, it, expect } from 'vitest';
import { FastMarch3 } from '../src/FastMarch3.js';
import { FastMarch2 } from '../src/FastMarch2.js';
import { check, expectClose, fc } from './helpers/arbitraries.js';

// Run a 2D march to completion, for the slab cross-check below.
function marchFM2(filter: FastMarch2): void {
    for (let i = 0; i < filter.getQuantity(); ++i) {
        filter.iterate();
    }
}

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

describe('FastMarch3 verification', () => {
    const MAXREAL = Number.MAX_VALUE;

    // Bounds of at least 3 (so there is an interior) with an interior seed
    // derived from the bounds, which avoids a filtered generator.
    const grid = fc.tuple(
        fc.integer({ min: 3, max: 6 }),
        fc.integer({ min: 3, max: 6 }),
        fc.integer({ min: 3, max: 6 }),
        fc.nat({ max: 1000 }), fc.nat({ max: 1000 }), fc.nat({ max: 1000 }),
        fc.constantFrom(0.5, 1, 2, 4))
        .map(([xB, yB, zB, rx, ry, rz, speed]) => ({
            xB, yB, zB, speed,
            sx: 1 + rx % (xB - 2),
            sy: 1 + ry % (yB - 2),
            sz: 1 + rz % (zB - 2)
        }));

    it('marks all six boundary faces zero speed and keeps seeds at time zero', () => {
        // Upstream marks only the box vertices and edges, which lets a face
        // voxel join the front and then index off the grid; the port marks
        // the faces (upstream issue #121).
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const seed = idx(sx, sy, sz, xB, yB);
            const filter = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], speed);
            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        const i = idx(x, y, z, xB, yB);
                        const onFace = (x === 0 || x === xB - 1
                            || y === 0 || y === yB - 1
                            || z === 0 || z === zB - 1);
                        expect(filter.isZeroSpeed(i)).toBe(onFace);
                        if (onFace) {
                            expect(filter.isValid(i)).toBe(false);
                            expect(filter.isFar(i)).toBe(false);
                        }
                    }
                }
            }
            expect(filter.getTime(seed)).toBe(0);
        });
    });

    it('accepts the trial voxels in nondecreasing time order', () => {
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const seed = idx(sx, sy, sz, xB, yB);
            const filter = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], speed);
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
                const accepted = trialsBefore.filter(i => !filter.isTrial(i));
                expect(accepted.length).toBe(1);
                const t = filter.getTime(accepted[0]);
                expect(t).toBeGreaterThanOrEqual(previous - 1e-12);
                previous = t;
            }
        });
    });

    it('leaves every interior voxel known and the faces untouched', () => {
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const seed = idx(sx, sy, sz, xB, yB);
            const filter = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], speed);
            march(filter);
            for (let z = 0; z < zB; ++z) {
                for (let y = 0; y < yB; ++y) {
                    for (let x = 0; x < xB; ++x) {
                        const i = idx(x, y, z, xB, yB);
                        const interior = (0 < x && x < xB - 1 && 0 < y && y < yB - 1
                            && 0 < z && z < zB - 1);
                        if (interior) {
                            expect(filter.isInterior(i)).toBe(true);
                            expect(Number.isFinite(filter.getTime(i))).toBe(true);
                            expect(filter.getTime(i)).toBeGreaterThanOrEqual(0);
                        } else {
                            expect(filter.getTime(i)).toBe(-MAXREAL);
                        }
                    }
                }
            }
            expect(filter.getBoundary()).toEqual([]);
            expect(filter.getInterior().length)
                .toBe((xB - 2) * (yB - 2) * (zB - 2));
        });
    });

    it('scales the crossing times by the reciprocal of the speed', () => {
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const seed = idx(sx, sy, sz, xB, yB);
            const unit = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], 1);
            const scaled = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], speed);
            march(unit);
            march(scaled);
            for (let z = 1; z + 1 < zB; ++z) {
                for (let y = 1; y + 1 < yB; ++y) {
                    for (let x = 1; x + 1 < xB; ++x) {
                        const i = idx(x, y, z, xB, yB);
                        expectClose(scaled.getTime(i), unit.getTime(i) / speed,
                            1e-9, 1e-9);
                    }
                }
            }
        });
    });

    it('is equivariant under reflecting an axis and under permuting axes', () => {
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const base = new FastMarch3(xB, yB, zB, 1, 1, 1,
                [idx(sx, sy, sz, xB, yB)], speed);
            march(base);

            const mirrored = new FastMarch3(xB, yB, zB, 1, 1, 1,
                [idx(xB - 1 - sx, sy, sz, xB, yB)], speed);
            march(mirrored);

            // Relabel (x,y,z) as (y,z,x).
            const permuted = new FastMarch3(yB, zB, xB, 1, 1, 1,
                [idx(sy, sz, sx, yB, zB)], speed);
            march(permuted);

            for (let z = 1; z + 1 < zB; ++z) {
                for (let y = 1; y + 1 < yB; ++y) {
                    for (let x = 1; x + 1 < xB; ++x) {
                        const t = base.getTime(idx(x, y, z, xB, yB));
                        expectClose(mirrored.getTime(
                            idx(xB - 1 - x, y, z, xB, yB)), t, 1e-9, 1e-9);
                        expectClose(permuted.getTime(idx(y, z, x, yB, zB)), t,
                            1e-9, 1e-9);
                    }
                }
            }
        });
    });

    it('brackets the arrival time between the Euclidean and L1 distances', () => {
        // For a single interior seed on a unit grid the first-order upwind
        // scheme overestimates the Euclidean distance and never exceeds the
        // grid (L1) distance. (With several fronts meeting, the upstream
        // negative-discriminant fallback can break the upper bound; see the
        // FastMarch2 quirk test.)
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const filter = new FastMarch3(xB, yB, zB, 1, 1, 1,
                [idx(sx, sy, sz, xB, yB)], speed);
            march(filter);
            for (let z = 1; z + 1 < zB; ++z) {
                for (let y = 1; y + 1 < yB; ++y) {
                    for (let x = 1; x + 1 < xB; ++x) {
                        const dx = Math.abs(x - sx);
                        const dy = Math.abs(y - sy);
                        const dz = Math.abs(z - sz);
                        const t = filter.getTime(idx(x, y, z, xB, yB));
                        expect(t).toBeGreaterThanOrEqual(
                            Math.sqrt(dx * dx + dy * dy + dz * dz) / speed - 1e-9);
                        expect(t).toBeLessThanOrEqual((dx + dy + dz) / speed + 1e-9);
                    }
                }
            }
        });
    });

    it('satisfies the local upwind bounds against its 6-neighbors', () => {
        check(grid, ({ xB, yB, zB, sx, sy, sz, speed }) => {
            const seed = idx(sx, sy, sz, xB, yB);
            const filter = new FastMarch3(xB, yB, zB, 1, 1, 1, [seed], speed);
            march(filter);
            const xy = xB * yB;
            for (let z = 1; z + 1 < zB; ++z) {
                for (let y = 1; y + 1 < yB; ++y) {
                    for (let x = 1; x + 1 < xB; ++x) {
                        const i = idx(x, y, z, xB, yB);
                        if (i === seed) {
                            continue;
                        }
                        const neighbors = [i - 1, i + 1, i - xB, i + xB, i - xy, i + xy]
                            .filter(j => filter.isValid(j))
                            .map(j => filter.getTime(j));
                        expect(neighbors.length).toBeGreaterThan(0);
                        const smallest = Math.min(...neighbors);
                        const t = filter.getTime(i);
                        expect(t).toBeGreaterThanOrEqual(smallest - 1e-9);
                        expect(t).toBeLessThanOrEqual(smallest + 1 / speed + 1e-9);
                    }
                }
            }
        });
    });

    it('agrees with the 2D march on a slab that is one voxel deep inside', () => {
        // A 3-deep grid has exactly one interior z layer, so the z-neighbors
        // of every interior voxel are zero-speed faces: no z-term is ever
        // active and the update reduces to the 2D one.
        const slab = fc.tuple(
            fc.integer({ min: 3, max: 7 }),
            fc.integer({ min: 3, max: 7 }),
            fc.nat({ max: 1000 }), fc.nat({ max: 1000 }),
            fc.constantFrom(0.5, 1, 2))
            .map(([xB, yB, rx, ry, speed]) => ({
                xB, yB, speed, sx: 1 + rx % (xB - 2), sy: 1 + ry % (yB - 2)
            }));
        check(slab, ({ xB, yB, sx, sy, speed }) => {
            const three = new FastMarch3(xB, yB, 3, 1, 1, 1,
                [idx(sx, sy, 1, xB, yB)], speed);
            const two = new FastMarch2(xB, yB, 1, 1, [sx + xB * sy], speed);
            march(three);
            marchFM2(two);
            for (let y = 1; y + 1 < yB; ++y) {
                for (let x = 1; x + 1 < xB; ++x) {
                    expectClose(three.getTime(idx(x, y, 1, xB, yB)),
                        two.getTime(x + xB * y), 1e-9, 1e-9);
                }
            }
        });
    });
});
