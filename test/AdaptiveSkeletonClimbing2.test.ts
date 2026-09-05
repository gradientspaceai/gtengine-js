import { describe, it, expect } from 'vitest';
import { AdaptiveSkeletonClimbing2 } from '../src/AdaptiveSkeletonClimbing2.js';
import { check, fc } from './helpers/arbitraries.js';

// Build a (2^N+1)-by-(2^N+1) integer image from a function f(x, y).
function makeImage(N: number, f: (x: number, y: number) => number): number[] {
    const size = (1 << N) + 1;
    const pixels = new Array<number>(size * size);
    for (let y = 0; y < size; ++y) {
        for (let x = 0; x < size; ++x) {
            pixels[x + size * y] = f(x, y);
        }
    }
    return pixels;
}

// Every extracted vertex that lies on a grid line (one integer coordinate)
// must interpolate the image to the extraction level exactly. Branch-point
// vertices from the plus-sign configuration are skipped (neither coordinate
// is an integer).
function checkVerticesOnLevelSet(pixels: number[], size: number, level: number,
    vertices: [number, number][]): void {
    for (const [vx, vy] of vertices) {
        const xIsInt = Math.abs(vx - Math.round(vx)) < 1.0e-12;
        const yIsInt = Math.abs(vy - Math.round(vy)) < 1.0e-12;
        if (yIsInt && !xIsInt) {
            const y = Math.round(vy);
            const x0 = Math.floor(vx);
            const t = vx - x0;
            const f0 = pixels[x0 + size * y];
            const f1 = pixels[x0 + 1 + size * y];
            expect(f0 + t * (f1 - f0)).toBeCloseTo(level, 10);
        } else if (xIsInt && !yIsInt) {
            const x = Math.round(vx);
            const y0 = Math.floor(vy);
            const t = vy - y0;
            const f0 = pixels[x + size * y0];
            const f1 = pixels[x + size * (y0 + 1)];
            expect(f0 + t * (f1 - f0)).toBeCloseTo(level, 10);
        } else {
            // A vertex with both coordinates integer would mean the level
            // equals an image sample, which the non-integer level forbids.
            expect(xIsInt && yIsInt).toBe(false);
        }
    }
}

describe('AdaptiveSkeletonClimbing2', () => {
    it('throws on invalid construction input', () => {
        expect(() => new AdaptiveSkeletonClimbing2(0, [0, 0, 0, 0]))
            .toThrow('Invalid input.');
    });

    it('extracts a vertical line from a linear x-ramp', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const pixels = makeImage(N, (x) => x);
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = 3.5;

        // Full-resolution extraction (no merging): one unit segment per
        // row of cells.
        const fine = asc.extract(level, N);
        expect(fine.edges.length).toBe(size - 1);
        for (const [vx] of fine.vertices) {
            expect(vx).toBeCloseTo(3.5, 12);
        }
        checkVerticesOnLevelSet(pixels, size, level, fine.vertices);

        // Adaptive extraction merges the monotone rectangles into a single
        // segment spanning the whole image.
        const coarse = asc.extract(level, -1);
        expect(coarse.edges.length).toBe(1);
        expect(coarse.vertices.length).toBe(2);
        const ys = coarse.vertices.map((v) => v[1]).sort((a, b) => a - b);
        expect(coarse.vertices[0][0]).toBeCloseTo(3.5, 12);
        expect(coarse.vertices[1][0]).toBeCloseTo(3.5, 12);
        expect(ys[0]).toBe(0);
        expect(ys[1]).toBe(size - 1);
    });

    it('extracts a horizontal line from a linear y-ramp', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const pixels = makeImage(N, (_x, y) => y);
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = 5.5;

        const coarse = asc.extract(level, -1);
        expect(coarse.edges.length).toBe(1);
        const xs = coarse.vertices.map((v) => v[0]).sort((a, b) => a - b);
        expect(coarse.vertices[0][1]).toBeCloseTo(5.5, 12);
        expect(coarse.vertices[1][1]).toBeCloseTo(5.5, 12);
        expect(xs[0]).toBe(0);
        expect(xs[1]).toBe(size - 1);
    });

    it('extracts a diagonal line from a linear x+y ramp', () => {
        const N = 3;
        const size = (1 << N) + 1;
        const pixels = makeImage(N, (x, y) => x + y);
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = 3.5;

        // Full-resolution extraction exercises the corner cases (types 5,
        // 6, 9, 10): each vertex lies exactly on the line x + y = 3.5.
        const { vertices, edges } = asc.extract(level, N);
        expect(edges.length).toBeGreaterThan(0);
        for (const [vx, vy] of vertices) {
            expect(vx + vy).toBeCloseTo(3.5, 12);
        }
        checkVerticesOnLevelSet(pixels, size, level, vertices);

        // Every edge endpoint index is valid.
        for (const [v0, v1] of edges) {
            expect(v0).toBeGreaterThanOrEqual(0);
            expect(v1).toBeLessThan(vertices.length);
        }
    });

    it('extracts a closed curve near a circle', () => {
        const N = 4;
        const size = (1 << N) + 1;
        const cx = 8, cy = 8;
        const pixels = makeImage(N, (x, y) => (x - cx) * (x - cx) + (y - cy) * (y - cy));
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = 25.5;
        const radius = Math.sqrt(level);

        for (const depth of [N, -1]) {
            const { vertices, edges } = asc.extract(level, depth);
            expect(edges.length).toBeGreaterThan(0);
            checkVerticesOnLevelSet(pixels, size, level, vertices);

            // All vertices lie near the true iso-contour (the circle of
            // radius sqrt(level) centered at (8,8)).
            for (const [vx, vy] of vertices) {
                const distance = Math.hypot(vx - cx, vy - cy);
                expect(Math.abs(distance - radius)).toBeLessThan(0.1);
            }
        }

        // At full resolution the segments are contained in unit cells, so
        // segment midpoints are also near the contour.
        const fine = asc.extract(level, N);
        for (const [v0, v1] of fine.edges) {
            const mx = 0.5 * (fine.vertices[v0][0] + fine.vertices[v1][0]);
            const my = 0.5 * (fine.vertices[v0][1] + fine.vertices[v1][1]);
            const distance = Math.hypot(mx - cx, my - cy);
            expect(Math.abs(distance - radius)).toBeLessThan(0.35);
            // Unit-cell segments are short.
            const dx = fine.vertices[v0][0] - fine.vertices[v1][0];
            const dy = fine.vertices[v0][1] - fine.vertices[v1][1];
            expect(Math.hypot(dx, dy)).toBeLessThanOrEqual(Math.SQRT2);
        }
    });

    it('returns nothing when the level does not intersect the image range', () => {
        const N = 3;
        const pixels = makeImage(N, (x, y) => x + y);
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);

        const above = asc.extract(1000.5, -1);
        expect(above.vertices.length).toBe(0);
        expect(above.edges.length).toBe(0);

        const below = asc.extract(-0.5, -1);
        expect(below.vertices.length).toBe(0);
        expect(below.edges.length).toBe(0);
    });

    it('disambiguates the four-intersection saddle configuration', () => {
        // Checkerboard-like saddle: every unit cell has all four edges
        // crossed. det = i00 * i11 - i01 * i10 = -100 < 0 for cell (0,0),
        // giving disjoint hyperbolic segments.
        const N = 1;
        const size = (1 << N) + 1;
        const pixels = [
            0, 10, 0,
            10, 0, 10,
            0, 10, 0
        ];
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = 5.5;

        const { vertices, edges } = asc.extract(level, N);
        // Each of the 4 cells contributes 2 disjoint segments.
        expect(edges.length).toBe(8);
        expect(vertices.length).toBe(16);
        checkVerticesOnLevelSet(pixels, size, level, vertices);

        // Each segment stays within one unit cell.
        for (const [v0, v1] of edges) {
            const dx = vertices[v0][0] - vertices[v1][0];
            const dy = vertices[v0][1] - vertices[v1][1];
            expect(Math.hypot(dx, dy)).toBeLessThan(Math.SQRT2);
        }

        // The inverted checkerboard flips the determinant sign (det > 0)
        // and still yields 8 segments.
        const inverted = new AdaptiveSkeletonClimbing2(N, pixels.map((p) => 10 - p));
        const flipped = inverted.extract(level, N);
        expect(flipped.edges.length).toBe(8);
        checkVerticesOnLevelSet(pixels.map((p) => 10 - p), size, level, flipped.vertices);
    });

    it('adds a branch point for the plus-sign configuration (det = 0)', () => {
        // Cell (0,0) has corners i00 = -2, i10 = 4, i01 = 1, i11 = -2 with
        // level -0.5: all four edges are crossed and
        // det = (-2)(-2) - (1)(4) = 0, the plus-sign configuration.
        const N = 1;
        const size = (1 << N) + 1;
        const pixels = [
            -2, 4, 4,
            1, -2, -2,
            1, -2, -2
        ];
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const level = -0.5;

        const { vertices, edges } = asc.extract(level, N);
        checkVerticesOnLevelSet(pixels, size, level, vertices);

        // The branch point is (fx2, y0) = (0.25, 0.5): the interpolated
        // ymin crossing x and the interpolated xmin crossing y.
        const branchEdges = edges.filter(([v0, v1]) => {
            const isBranch = (v: [number, number]): boolean =>
                Math.abs(v[0] - 0.25) < 1.0e-12 && Math.abs(v[1] - 0.5) < 1.0e-12;
            return isBranch(vertices[v0]) || isBranch(vertices[v1]);
        });
        expect(branchEdges.length).toBe(4);
    });

    it('makeUnique removes duplicate vertices and edges', () => {
        const N = 1;
        const pixels = [
            0, 10, 0,
            10, 0, 10,
            0, 10, 0
        ];
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const { vertices, edges } = asc.extract(5.5, N);
        expect(vertices.length).toBe(16);
        expect(edges.length).toBe(8);

        const originalVertices = vertices.map((v) => [...v] as [number, number]);
        asc.makeUnique(vertices, edges);

        // The 3x3 grid has 12 unit grid edges, each crossed exactly once.
        expect(vertices.length).toBe(12);
        expect(edges.length).toBe(8);

        // No duplicate vertices remain.
        const keys = new Set(vertices.map((v) => v[0] + ',' + v[1]));
        expect(keys.size).toBe(vertices.length);

        // Every original vertex is still present.
        for (const v of originalVertices) {
            expect(keys.has(v[0] + ',' + v[1])).toBe(true);
        }

        // All edge indices are valid and every edge is unique.
        const edgeKeys = new Set<string>();
        for (const [v0, v1] of edges) {
            expect(v0).toBeGreaterThanOrEqual(0);
            expect(v0).toBeLessThan(vertices.length);
            expect(v1).toBeGreaterThanOrEqual(0);
            expect(v1).toBeLessThan(vertices.length);
            edgeKeys.add(v0 + ',' + v1);
        }
        expect(edgeKeys.size).toBe(edges.length);
    });

    it('makeUnique on the circle produces a closed curve', () => {
        const N = 4;
        const cx = 8, cy = 8;
        const pixels = makeImage(N, (x, y) => (x - cx) * (x - cx) + (y - cy) * (y - cy));
        const asc = new AdaptiveSkeletonClimbing2(N, pixels);
        const { vertices, edges } = asc.extract(25.5, N);

        const before = vertices.length;
        asc.makeUnique(vertices, edges);
        expect(vertices.length).toBeLessThan(before);

        // A closed curve at full resolution: every vertex has exactly two
        // incident edges.
        const valence = new Array<number>(vertices.length).fill(0);
        for (const [v0, v1] of edges) {
            ++valence[v0];
            ++valence[v1];
        }
        for (const v of valence) {
            expect(v).toBe(2);
        }
    });
});


// ---------------------------------------------------------------------------
// Verification wave (V24): properties cross-checking the port against the
// upstream AdaptiveSkeletonClimbing2.h algorithm.
// ---------------------------------------------------------------------------

// Independent marching squares on the same image. Each cell contributes one
// segment when exactly two of its four edges are crossed; cells with four
// crossings (the ambiguous saddle) are counted so the caller can skip them.
function marchingSquares(pixels: readonly number[], size: number, level: number) {
    const key = (x: number, y: number) => `${x.toFixed(9)},${y.toFixed(9)}`;
    const segments: string[] = [];
    const crossings = new Set<string>();
    let ambiguous = 0;
    for (let y = 0; y + 1 < size; ++y) {
        for (let x = 0; x + 1 < size; ++x) {
            const f00 = pixels[x + size * y];
            const f10 = pixels[x + 1 + size * y];
            const f01 = pixels[x + size * (y + 1)];
            const f11 = pixels[x + 1 + size * (y + 1)];
            const points: [number, number][] = [];
            if ((f00 < level) !== (f01 < level)) {
                points.push([x, y + (level - f00) / (f01 - f00)]);
            }
            if ((f10 < level) !== (f11 < level)) {
                points.push([x + 1, y + (level - f10) / (f11 - f10)]);
            }
            if ((f00 < level) !== (f10 < level)) {
                points.push([x + (level - f00) / (f10 - f00), y]);
            }
            if ((f01 < level) !== (f11 < level)) {
                points.push([x + (level - f01) / (f11 - f01), y + 1]);
            }
            for (const p of points) {
                crossings.add(key(p[0], p[1]));
            }
            if (points.length === 2) {
                const a = key(points[0][0], points[0][1]);
                const b = key(points[1][0], points[1][1]);
                segments.push(a < b ? `${a}|${b}` : `${b}|${a}`);
            } else if (points.length === 4) {
                ++ambiguous;
            }
        }
    }
    return { segments: segments.sort(), crossings, ambiguous };
}

function segmentKeys(vertices: readonly [number, number][],
    edges: readonly [number, number][]): string[] {
    const key = (v: readonly number[]) => `${v[0].toFixed(9)},${v[1].toFixed(9)}`;
    return edges.map(e => {
        const a = key(vertices[e[0]]), b = key(vertices[e[1]]);
        return a < b ? `${a}|${b}` : `${b}|${a}`;
    }).sort();
}

describe('AdaptiveSkeletonClimbing2 verification', () => {
    // Smooth integer images: a few low-frequency terms rounded to integers.
    // Saddle cells (the ambiguous four-crossing configuration) essentially
    // never occur for these, which lets the reference marching squares stay
    // an independent computation instead of replicating the upstream
    // determinant disambiguation.
    const coefficient = fc.integer({ min: -10, max: 10 });
    const smoothImage = (N: number) =>
        fc.tuple(coefficient, coefficient, coefficient, coefficient)
            .map(([a, b, c, d]) => {
                const size = (1 << N) + 1;
                const pixels = new Array<number>(size * size);
                for (let y = 0; y < size; ++y) {
                    for (let x = 0; x < size; ++x) {
                        const u = x / (size - 1), v = y / (size - 1);
                        pixels[x + size * y] = Math.round(
                            a * Math.sin(2 * Math.PI * u) + b * Math.cos(2 * Math.PI * v)
                            + c * u * v + d * Math.sin(Math.PI * (u + v)));
                    }
                }
                return pixels;
            });

    // An ellipse whose level curve is strictly inside the image, so the
    // extracted curve is closed.
    const ellipseImage = (N: number) =>
        fc.tuple(fc.integer({ min: 25, max: 40 }), fc.integer({ min: 25, max: 40 }))
            .map(([ra, rb]) => {
                const size = (1 << N) + 1;
                const center = (size - 1) / 2;
                const a = (ra / 100) * (size - 1), b = (rb / 100) * (size - 1);
                const pixels = new Array<number>(size * size);
                for (let y = 0; y < size; ++y) {
                    for (let x = 0; x < size; ++x) {
                        pixels[x + size * y] = Math.round(20 * (1
                            - ((x - center) / a) ** 2 - ((y - center) / b) ** 2));
                    }
                }
                return pixels;
            });

    it('full-resolution extraction reproduces marching squares', () => {
        // depth >= N prevents every merge, so each rectangle is a single
        // image cell and the result must be the marching-squares segments.
        const N = 3;
        const size = (1 << N) + 1;
        check(smoothImage(N), pixels => {
            const level = 0.5;
            const reference = marchingSquares(pixels, size, level);
            if (reference.ambiguous > 0) { return; }
            const asc = new AdaptiveSkeletonClimbing2(N, pixels);
            const { vertices, edges } = asc.extract(level, N);
            expect(segmentKeys(vertices, edges)).toEqual(reference.segments);
            // Every extracted edge has two distinct endpoints.
            for (const [v0, v1] of edges) {
                expect(v0).not.toBe(v1);
            }
        }, 80);
    });

    it('every vertex sits on a grid edge at the extraction level', () => {
        const N = 4;
        const size = (1 << N) + 1;
        check(fc.tuple(smoothImage(N), fc.constantFrom(-1, 0, 1, 2, N)),
            ([pixels, depth]) => {
                const asc = new AdaptiveSkeletonClimbing2(N, pixels);
                const { vertices } = asc.extract(0.5, depth);
                checkVerticesOnLevelSet(pixels, size, 0.5, vertices);
            }, 60);
    });

    it('adaptive vertices are a subset of the full-resolution crossings', () => {
        // Merging only removes crossings from the output; it never invents a
        // point that is not a zero crossing of a grid edge.
        const N = 4;
        const size = (1 << N) + 1;
        const key = (v: readonly number[]) => `${v[0].toFixed(9)},${v[1].toFixed(9)}`;
        check(smoothImage(N), pixels => {
            const reference = marchingSquares(pixels, size, 0.5);
            if (reference.ambiguous > 0) { return; }
            const asc = new AdaptiveSkeletonClimbing2(N, pixels);
            const adaptive = asc.extract(0.5, -1);
            for (const vertex of adaptive.vertices) {
                expect(reference.crossings.has(key(vertex)),
                    `vertex ${key(vertex)}`).toBe(true);
            }
            // The adaptive extraction never produces more segments than the
            // full-resolution one.
            const full = asc.extract(0.5, N);
            expect(adaptive.edges.length).toBeLessThanOrEqual(full.edges.length);
        }, 60);
    });

    it('a level outside the image range extracts nothing', () => {
        const N = 3;
        check(fc.tuple(smoothImage(N), fc.constantFrom(-1, 0, N)), ([pixels, depth]) => {
            const asc = new AdaptiveSkeletonClimbing2(N, pixels);
            const above = Math.max(...pixels) + 0.5;
            const below = Math.min(...pixels) - 0.5;
            for (const level of [above, below]) {
                const { vertices, edges } = asc.extract(level, depth);
                expect(vertices).toEqual([]);
                expect(edges).toEqual([]);
            }
        }, 60);
    });

    it('closed level curves have even vertex degree after makeUnique', () => {
        const N = 4;
        check(fc.tuple(ellipseImage(N), fc.constantFrom(-1, 0, 1, N)),
            ([pixels, depth]) => {
                const asc = new AdaptiveSkeletonClimbing2(N, pixels);
                const { vertices, edges } = asc.extract(0.5, depth);
                asc.makeUnique(vertices, edges);
                expect(vertices.length).toBeGreaterThan(0);
                const degree = new Array<number>(vertices.length).fill(0);
                for (const [v0, v1] of edges) {
                    ++degree[v0];
                    ++degree[v1];
                }
                for (let i = 0; i < vertices.length; ++i) {
                    expect(degree[i] % 2, `vertex ${i}`).toBe(0);
                    expect(degree[i]).toBeGreaterThan(0);
                }
            }, 60);
    });

    it('makeUnique preserves the segment geometry and packs first-encounter order', () => {
        const N = 3;
        check(fc.tuple(smoothImage(N), fc.constantFrom(-1, 0, N)), ([pixels, depth]) => {
            const asc = new AdaptiveSkeletonClimbing2(N, pixels);
            const { vertices, edges } = asc.extract(0.5, depth);
            if (edges.length === 0) { return; }
            const before = segmentKeys(vertices, edges);
            const key = (v: readonly number[]) => `${v[0].toFixed(9)},${v[1].toFixed(9)}`;

            // Independent expectation: distinct vertices in first-encounter
            // order (upstream packs by the std::map insertion index).
            const seen: string[] = [];
            const seenSet = new Set<string>();
            for (const vertex of vertices) {
                if (!seenSet.has(key(vertex))) {
                    seenSet.add(key(vertex));
                    seen.push(key(vertex));
                }
            }

            const original = vertices.map(v => [v[0], v[1]] as [number, number]);
            asc.makeUnique(vertices, edges);
            expect(vertices.map(key)).toEqual(seen);
            // Edges are unique and index the packed vertices.
            const asStrings = edges.map(e => `${e[0]}_${e[1]}`);
            expect(new Set(asStrings).size).toBe(asStrings.length);
            for (const [v0, v1] of edges) {
                expect(v0).toBeGreaterThanOrEqual(0);
                expect(v0).toBeLessThan(vertices.length);
                expect(v1).toBeGreaterThanOrEqual(0);
                expect(v1).toBeLessThan(vertices.length);
            }
            // The undirected segment set is unchanged apart from duplicates.
            expect(new Set(segmentKeys(vertices, edges)))
                .toEqual(new Set(before));
            // Every original vertex is still present.
            for (const vertex of original) {
                expect(seenSet.has(key(vertex))).toBe(true);
            }
        }, 60);
    });

    it('rejects a nonpositive N (upstream #52: the comment says N >= 0)', () => {
        check(fc.integer({ min: -4, max: 0 }), N => {
            const size = (1 << Math.max(N, 0)) + 1;
            const pixels = new Array<number>(size * size).fill(0);
            expect(() => new AdaptiveSkeletonClimbing2(N, pixels))
                .toThrow('Invalid input.');
        });
    });

    it('is idempotent across repeated extractions with the same arguments', () => {
        // The merge trees are rebuilt by every SetLevel call, so a second
        // extraction must not see stale state from the first.
        const N = 4;
        check(fc.tuple(smoothImage(N), fc.constantFrom(-1, 0, 1, N)),
            ([pixels, depth]) => {
                const asc = new AdaptiveSkeletonClimbing2(N, pixels);
                const first = asc.extract(0.5, depth);
                asc.extract(1.5, -1);
                asc.extract(-2.5, N);
                const second = asc.extract(0.5, depth);
                expect(segmentKeys(second.vertices, second.edges))
                    .toEqual(segmentKeys(first.vertices, first.edges));
            }, 60);
    });
});
