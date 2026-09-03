import { describe, it, expect } from 'vitest';
import { AdaptiveSkeletonClimbing2 } from '../src/AdaptiveSkeletonClimbing2.js';

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
