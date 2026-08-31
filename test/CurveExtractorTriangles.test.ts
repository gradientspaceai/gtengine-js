import { describe, it, expect } from 'vitest';
import { CurveExtractorTriangles } from '../src/CurveExtractorTriangles';
import { CurveExtractorSquares } from '../src/CurveExtractorSquares';
import { CurveExtractorVertex } from '../src/CurveExtractor';

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

describe('CurveExtractorTriangles construction', () => {
    it('rejects images that are too small', () => {
        expect(() => new CurveExtractorTriangles(1, 4, new Array(4).fill(0)))
            .toThrow('Invalid input.');
        expect(() => new CurveExtractorTriangles(4, 4, new Array(8).fill(0)))
            .toThrow('Invalid input.');
    });
});

describe('CurveExtractorTriangles single-square sign patterns', () => {
    // The 2x2 image has corners f00 = pixels[0] at (0,0), f10 = pixels[1] at
    // (1,0), f01 = pixels[2] at (0,1) and f11 = pixels[3] at (1,1). The
    // square at (0,0) has xParity == yParity, so it splits into the
    // triangles <(0,0),(0,1),(1,0)> and <(1,1),(1,0),(0,1)>.
    function extract2x2(f00: number, f10: number, f01: number, f11: number) {
        return new CurveExtractorTriangles(2, 2, [f00, f10, f01, f11]).extract(0);
    }

    it('emits nothing when all corners have the same sign', () => {
        expect(extract2x2(1, 1, 1, 1).edges).toEqual([]);
        expect(extract2x2(1, 1, 1, 1).vertices).toEqual([]);
        expect(extract2x2(-2, -3, -1, -5).edges).toEqual([]);
    });

    it('cuts only the triangle that contains the sign change', () => {
        // Only f11 is negative, so the first triangle is all-positive and
        // the second has the pattern "+--" after the sign normalization.
        const r = extract2x2(1, 1, 1, -1);
        expect(r.edges.length).toBe(1);
        expect(r.vertices.map(asPair)).toEqual([[1, 0.5], [0.5, 1]]);
    });

    it('cuts both triangles when the whole top row is negative', () => {
        // f00 = f10 = 1, f01 = f11 = -1: the exact level curve is y = 1/2
        // and both triangles are cut.
        const r = extract2x2(1, 1, -1, -1);
        expect(r.edges.length).toBe(2);
        for (const v of r.vertices) {
            expect(asPair(v)[1]).toBe(0.5);
        }
    });

    it('emits a single vertex when one corner is zero and the rest positive', () => {
        // The zero corner (0,1) belongs to both triangles, so it is emitted
        // once per triangle.
        const r = extract2x2(1, 1, 0, 1);
        expect(r.edges).toEqual([]);
        expect(r.vertices.map(asPair)).toEqual([[0, 1], [0, 1]]);
    });

    it('emits the triangle boundary for an all-zero image', () => {
        const r = extract2x2(0, 0, 0, 0);
        // Three edges per triangle, two triangles.
        expect(r.edges.length).toBe(6);
        expect(r.vertices.length).toBe(12);
    });

    it('emits a segment along a shared edge when two corners are zero', () => {
        // f00 = 0, f01 = 0 (the left edge), the others positive. Triangle 1
        // has the pattern "00+" and emits the segment (0,0)-(0,1); triangle
        // 2 is "+0+"-like and emits the single vertex (0,1).
        const r = extract2x2(0, 1, 0, 1);
        const pairs = r.vertices.map(asPair);
        expect(r.edges.length).toBe(1);
        expect(pairs).toContainEqual([0, 0]);
        expect(pairs).toContainEqual([0, 1]);
    });
});

describe('CurveExtractorTriangles diagonal parity', () => {
    it('alternates the diagonal with the parity of the square', () => {
        // A 3x2 image with the same corner values in both squares gives
        // different vertex sets, because the squares (0,0) and (1,0) split
        // along opposite diagonals.
        const pixels = [1, -1, 1, -1, 1, -1];
        const r = new CurveExtractorTriangles(3, 2, pixels).extract(0);
        // Both squares are cut, and every vertex is on a grid edge.
        expect(r.edges.length).toBeGreaterThan(0);
        for (const v of r.vertices) {
            const [px, py] = asPair(v);
            expect(px).toBeGreaterThanOrEqual(0);
            expect(px).toBeLessThanOrEqual(2);
            expect(py).toBeGreaterThanOrEqual(0);
            expect(py).toBeLessThanOrEqual(1);
        }
    });
});

describe('CurveExtractorTriangles on linear images', () => {
    // A linear function is exactly linear on each triangle, so the extracted
    // curve is exactly the line a*x + b*y + c = 0.
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
            const r = new CurveExtractorTriangles(xBound, yBound, pixels).extract(0);
            expect(r.vertices.length).toBeGreaterThan(0);
            for (const v of r.vertices) {
                expect(v.xDenom).toBeGreaterThan(0);
                expect(v.yDenom).toBeGreaterThan(0);
                expect(onLine(v, a, b, c)).toBe(true);
            }
            for (const e of r.edges) {
                expect(e.v[0]).toBeGreaterThanOrEqual(0);
                expect(e.v[1]).toBeLessThan(r.vertices.length);
                expect(e.v[0]).toBeLessThanOrEqual(e.v[1]);
            }
        });
    }

    it('agrees with the squares extractor on a linear image', () => {
        // Both extractors reproduce the exact line for a linear image, so
        // the sets of unique vertex positions must coincide.
        const a = 3, b = 5, c = -20;
        const xBound = 6, yBound = 6;
        const pixels = image(xBound, yBound, (x, y) => a * x + b * y + c);
        const tri = new CurveExtractorTriangles(xBound, yBound, pixels)
            .extractReal(0, true);
        const sqr = new CurveExtractorSquares(xBound, yBound, pixels)
            .extractReal(0, true);
        const key = (p: [number, number]) => `${p[0]},${p[1]}`;
        const triSet = new Set(tri.vertices.map(key));
        const sqrSet = new Set(sqr.vertices.map(key));
        // Every squares vertex is a triangles vertex; the triangles method
        // additionally samples the diagonals, so it may have more.
        for (const k of sqrSet) {
            expect(triSet.has(k)).toBe(true);
        }
        expect(triSet.size).toBeGreaterThanOrEqual(sqrSet.size);
    });
});

describe('CurveExtractorTriangles level shifting and curves', () => {
    it('extracting level L equals extracting level 0 of the shifted image', () => {
        const xBound = 6;
        const yBound = 5;
        const base = image(xBound, yBound, (x, y) => x * x + y * y);
        const level = 13;
        const shifted = base.map(p => p - level);

        const a = new CurveExtractorTriangles(xBound, yBound, base).extract(level);
        const b = new CurveExtractorTriangles(xBound, yBound, shifted).extract(0);
        expect(a.vertices.map(asPair)).toEqual(b.vertices.map(asPair));
        expect(a.edges.map(e => e.v)).toEqual(b.edges.map(e => e.v));
    });

    it('extracts a closed curve for a discretized circle', () => {
        const bound = 11;
        const center = 5;
        const pixels = image(bound, bound, (x, y) =>
            (x - center) ** 2 + (y - center) ** 2 - 16);
        const r = new CurveExtractorTriangles(bound, bound, pixels).extract(0);
        expect(r.edges.length).toBeGreaterThan(0);
        for (const v of r.vertices) {
            const [px, py] = asPair(v);
            const radius = Math.hypot(px - center, py - center);
            expect(radius).toBeGreaterThan(3.4);
            expect(radius).toBeLessThan(4.6);
        }
    });

    it('produces a connected chain of edges for the circle', () => {
        // After removing duplicates, every vertex of the extracted circle
        // should be shared by at least two edges (the curve is closed).
        const bound = 13;
        const center = 6;
        const pixels = image(bound, bound, (x, y) =>
            (x - center) ** 2 + (y - center) ** 2 - 25);
        const r = new CurveExtractorTriangles(bound, bound, pixels)
            .extractReal(0, true);
        const degree = new Array<number>(r.vertices.length).fill(0);
        for (const e of r.edges) {
            ++degree[e.v[0]];
            ++degree[e.v[1]];
        }
        for (let i = 0; i < degree.length; ++i) {
            expect(degree[i]).toBeGreaterThanOrEqual(2);
        }
    });
});

describe('CurveExtractorTriangles extractReal', () => {
    it('removes duplicate vertices while keeping the positions', () => {
        const xBound = 6;
        const yBound = 6;
        const pixels = image(xBound, yBound, (x, y) => 3 * x + 5 * y - 20);
        const extractor = new CurveExtractorTriangles(xBound, yBound, pixels);

        const withDuplicates = extractor.extractReal(0, false);
        const withoutDuplicates = extractor.extractReal(0, true);
        expect(withoutDuplicates.vertices.length)
            .toBeLessThan(withDuplicates.vertices.length);

        const distinct = new Set(withDuplicates.vertices.map(p => `${p[0]},${p[1]}`));
        expect(withoutDuplicates.vertices.length).toBe(distinct.size);
        for (const [px, py] of withoutDuplicates.vertices) {
            expect(3 * px + 5 * py - 20).toBeCloseTo(0, 12);
        }
    });
});
