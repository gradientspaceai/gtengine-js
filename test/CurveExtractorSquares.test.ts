import { describe, it, expect } from 'vitest';
import { CurveExtractorSquares } from '../src/CurveExtractorSquares.js';
import { CurveExtractorVertex } from '../src/CurveExtractor.js';

// Build a pixel image from a function of the integer lattice.
function image(xBound: number, yBound: number,
    f: (x: number, y: number) => number): number[] {
    const pixels: number[] = [];
    for (let y = 0; y < yBound; ++y) {
        for (let x = 0; x < xBound; ++x) {
            pixels.push(f(x, y));
        }
    }
    return pixels;
}

function asPair(v: CurveExtractorVertex): [number, number] {
    return [v.xNumer / v.xDenom, v.yNumer / v.yDenom];
}

// Exact test that the rational vertex (xn/xd, yn/yd) satisfies
// a*x + b*y + c = 0 for integer a, b, c.
function onLine(v: CurveExtractorVertex, a: number, b: number, c: number): boolean {
    return a * v.xNumer * v.yDenom + b * v.yNumer * v.xDenom
        + c * v.xDenom * v.yDenom === 0;
}

describe('CurveExtractorSquares construction', () => {
    it('rejects images that are too small', () => {
        expect(() => new CurveExtractorSquares(1, 4, new Array(4).fill(0)))
            .toThrow('Invalid input.');
        expect(() => new CurveExtractorSquares(4, 1, new Array(4).fill(0)))
            .toThrow('Invalid input.');
        expect(() => new CurveExtractorSquares(4, 4, new Array(8).fill(0)))
            .toThrow('Invalid input.');
    });
});

describe('CurveExtractorSquares single-square sign patterns', () => {
    // The corners of the single square of a 2x2 image are
    //   f00 = pixels[0] at (0,0), f10 = pixels[1] at (1,0),
    //   f01 = pixels[2] at (0,1), f11 = pixels[3] at (1,1).
    function extract2x2(f00: number, f10: number, f01: number, f11: number) {
        const extractor = new CurveExtractorSquares(2, 2, [f00, f10, f01, f11]);
        return extractor.extract(0);
    }

    it('emits nothing when all corners are positive (++++)', () => {
        const r = extract2x2(1, 1, 1, 1);
        expect(r.vertices).toEqual([]);
        expect(r.edges).toEqual([]);
    });

    it('emits nothing when all corners are negative (----)', () => {
        const r = extract2x2(-1, -1, -1, -1);
        expect(r.vertices).toEqual([]);
        expect(r.edges).toEqual([]);
    });

    it('interpolates the two crossings for the case +++-', () => {
        // f01 = -1 at (0,1); the level curve crosses the top edge at x = 1/2
        // and the left edge at y = 1/2.
        const r = extract2x2(1, 1, -1, 1);
        expect(r.edges.length).toBe(1);
        expect(r.vertices.length).toBe(2);
        expect(asPair(r.vertices[0])).toEqual([0.5, 1]);
        expect(asPair(r.vertices[1])).toEqual([0, 0.5]);
        expect(r.edges[0].v).toEqual([0, 1]);
    });

    it('is sign-symmetric: negating the image gives the same curve', () => {
        const positive = extract2x2(1, 1, -1, 1);
        const negative = extract2x2(-1, -1, 1, -1);
        expect(negative.vertices.map(asPair)).toEqual(positive.vertices.map(asPair));
        expect(negative.edges.map(e => e.v)).toEqual(positive.edges.map(e => e.v));
    });

    it('produces a horizontal segment for the case ++-- (top row negative)', () => {
        // f00 = f10 = 1, f01 = f11 = -1: the curve is the horizontal line
        // y = 1/2.
        const r = extract2x2(1, 1, -1, -1);
        expect(r.edges.length).toBe(1);
        expect(r.vertices.map(asPair)).toEqual([[0, 0.5], [1, 0.5]]);
    });

    it('produces a vertical segment when the right column is negative', () => {
        const r = extract2x2(1, -1, 1, -1);
        expect(r.edges.length).toBe(1);
        const pairs = r.vertices.map(asPair).sort((p, q) => p[1] - q[1]);
        expect(pairs).toEqual([[0.5, 0], [0.5, 1]]);
    });

    it('emits a single vertex when exactly one corner is zero (+++0)', () => {
        const r = extract2x2(1, 1, 0, 1);
        expect(r.edges).toEqual([]);
        expect(r.vertices.length).toBe(1);
        expect(asPair(r.vertices[0])).toEqual([0, 1]);
    });

    it('emits two vertices for the diagonal-zero case +0+0', () => {
        const r = extract2x2(1, 0, 0, 1);
        expect(r.edges).toEqual([]);
        expect(r.vertices.map(asPair)).toEqual([[1, 0], [0, 1]]);
    });

    it('emits the four boundary edges for the all-zero square', () => {
        const r = extract2x2(0, 0, 0, 0);
        expect(r.edges.length).toBe(4);
        expect(r.vertices.length).toBe(8);
        expect(r.vertices.map(asPair)).toEqual([
            [0, 0], [1, 0],
            [1, 0], [1, 1],
            [1, 1], [0, 1],
            [0, 1], [0, 0]
        ]);
    });

    it('resolves the ambiguous saddle case +-+- by the determinant test', () => {
        // f00 = 4, f10 = -1, f11 = 4, f01 = -1 is a saddle. The upstream
        // determinant decides how to connect the four edge crossings, and
        // the result is two segments (four vertices, two edges).
        const r = extract2x2(4, -1, -1, 4);
        expect(r.edges.length).toBe(2);
        expect(r.vertices.length).toBe(4);
        for (const v of r.vertices) {
            const [px, py] = asPair(v);
            expect(px).toBeGreaterThanOrEqual(0);
            expect(px).toBeLessThanOrEqual(1);
            expect(py).toBeGreaterThanOrEqual(0);
            expect(py).toBeLessThanOrEqual(1);
        }
    });

    it('emits the degenerate fan when the saddle determinant is zero', () => {
        // f00 = 1, f10 = -1, f11 = 1, f01 = -1 is the symmetric saddle whose
        // crossings all meet at the center (1/2, 1/2); upstream emits four
        // edges radiating from that point.
        const r = extract2x2(1, -1, -1, 1);
        expect(r.edges.length).toBe(4);
        expect(r.vertices.length).toBe(8);
        // Every edge starts at the center.
        for (let e = 0; e < 4; ++e) {
            expect(asPair(r.vertices[2 * e])).toEqual([0.5, 0.5]);
        }
        expect(r.vertices.filter((_, i) => i % 2 === 1).map(asPair)).toEqual([
            [0.5, 0], [0.5, 1], [0, 0.5], [1, 0.5]
        ]);
    });
});

describe('CurveExtractorSquares on linear images', () => {
    // For f(x,y) = a*x + b*y + c the exact level curve is a line, and every
    // vertex that the extractor emits must lie on it exactly.
    const linearCases: [number, number, number][] = [
        [3, 5, -11],
        [1, 1, -4],
        [2, -3, 5],
        [-1, 4, -6],
        [1, 0, -3],
        [0, 1, -2]
    ];

    for (const [a, b, c] of linearCases) {
        it(`places every vertex on ${a}x + ${b}y + ${c} = 0`, () => {
            const xBound = 7;
            const yBound = 6;
            const pixels = image(xBound, yBound, (x, y) => a * x + b * y + c);
            const extractor = new CurveExtractorSquares(xBound, yBound, pixels);
            const r = extractor.extract(0);
            expect(r.vertices.length).toBeGreaterThan(0);
            for (const v of r.vertices) {
                expect(v.xDenom).toBeGreaterThan(0);
                expect(v.yDenom).toBeGreaterThan(0);
                expect(onLine(v, a, b, c)).toBe(true);
            }
            // Every edge index is in range.
            for (const e of r.edges) {
                expect(e.v[0]).toBeGreaterThanOrEqual(0);
                expect(e.v[1]).toBeLessThan(r.vertices.length);
                expect(e.v[0]).toBeLessThanOrEqual(e.v[1]);
            }
        });
    }

    it('keeps every vertex inside the image bounds', () => {
        const pixels = image(8, 8, (x, y) => 2 * x + 3 * y - 15);
        const extractor = new CurveExtractorSquares(8, 8, pixels);
        const r = extractor.extract(0);
        for (const v of r.vertices) {
            const [px, py] = asPair(v);
            expect(px).toBeGreaterThanOrEqual(0);
            expect(px).toBeLessThanOrEqual(7);
            expect(py).toBeGreaterThanOrEqual(0);
            expect(py).toBeLessThanOrEqual(7);
        }
    });
});

describe('CurveExtractorSquares level shifting', () => {
    it('extracting level L equals extracting level 0 of the shifted image', () => {
        const xBound = 6;
        const yBound = 5;
        const base = image(xBound, yBound, (x, y) => x * x + y * y);
        const level = 13;
        const shifted = base.map(p => p - level);

        const a = new CurveExtractorSquares(xBound, yBound, base).extract(level);
        const b = new CurveExtractorSquares(xBound, yBound, shifted).extract(0);
        expect(a.vertices.map(asPair)).toEqual(b.vertices.map(asPair));
        expect(a.edges.map(e => e.v)).toEqual(b.edges.map(e => e.v));
    });

    it('extracts a closed curve for a discretized circle', () => {
        const bound = 11;
        const center = 5;
        const pixels = image(bound, bound, (x, y) =>
            (x - center) ** 2 + (y - center) ** 2 - 16);
        const extractor = new CurveExtractorSquares(bound, bound, pixels);
        const r = extractor.extract(0);
        expect(r.edges.length).toBeGreaterThan(0);
        // Every vertex is close to the circle of radius 4.
        for (const v of r.vertices) {
            const [px, py] = asPair(v);
            const radius = Math.hypot(px - center, py - center);
            expect(radius).toBeGreaterThan(3.4);
            expect(radius).toBeLessThan(4.6);
        }
    });
});

describe('CurveExtractorSquares extractReal', () => {
    it('produces floating-point vertices and removes duplicates on request', () => {
        const xBound = 6;
        const yBound = 6;
        const pixels = image(xBound, yBound, (x, y) => 3 * x + 5 * y - 20);
        const extractor = new CurveExtractorSquares(xBound, yBound, pixels);

        const withDuplicates = extractor.extractReal(0, false);
        const withoutDuplicates = extractor.extractReal(0, true);
        expect(withoutDuplicates.vertices.length)
            .toBeLessThan(withDuplicates.vertices.length);

        // The unique vertices are exactly the distinct positions.
        const distinct = new Set(withDuplicates.vertices.map(p => `${p[0]},${p[1]}`));
        expect(withoutDuplicates.vertices.length).toBe(distinct.size);
        for (const p of withoutDuplicates.vertices) {
            expect(distinct.has(`${p[0]},${p[1]}`)).toBe(true);
        }
        // Every unique vertex still lies on the line.
        for (const [px, py] of withoutDuplicates.vertices) {
            expect(3 * px + 5 * py - 20).toBeCloseTo(0, 12);
        }
    });
});
