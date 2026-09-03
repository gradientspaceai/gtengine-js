import { describe, it, expect } from 'vitest';
import {
    CurveExtractor, CurveExtractorEdge, CurveExtractorVertex
} from '../src/CurveExtractor.js';

// A concrete extractor used to exercise the protected helpers and the
// abstract extract() contract. It does not implement a real marching-squares
// algorithm; it emits a fixed set of edges so the base-class bookkeeping can
// be tested exactly.
class TestExtractor extends CurveExtractor {
    // Each entry is a pair of rational vertices (as the 8 numerator and
    // denominator values that addEdge takes).
    edgeSpecs: number[][] = [];

    constructor(xBound: number, yBound: number, pixels: number[]) {
        super(xBound, yBound, pixels);
    }

    override extract(level: number): {
        vertices: CurveExtractorVertex[];
        edges: CurveExtractorEdge[];
    } {
        void level;
        const vertices: CurveExtractorVertex[] = [];
        const edges: CurveExtractorEdge[] = [];
        for (const s of this.edgeSpecs) {
            this.addEdge(vertices, edges, s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
        }
        return { vertices, edges };
    }

    addVertexPublic(vertices: CurveExtractorVertex[],
        xNumer: number, xDenom: number, yNumer: number, yDenom: number): void {
        this.addVertex(vertices, xNumer, xDenom, yNumer, yDenom);
    }

    getPixelCount(): number {
        return this.mPixels.length;
    }
}

function makeExtractor(edgeSpecs: number[][]): TestExtractor {
    const extractor = new TestExtractor(3, 3, new Array<number>(9).fill(0));
    extractor.edgeSpecs = edgeSpecs;
    return extractor;
}

describe('CurveExtractorVertex', () => {
    it('normalizes the signs so both denominators are positive', () => {
        const v = new CurveExtractorVertex(-3, -4, 5, -6);
        expect(v.xNumer).toBe(3);
        expect(v.xDenom).toBe(4);
        expect(v.yNumer).toBe(-5);
        expect(v.yDenom).toBe(6);
    });

    it('leaves a positive denominator alone', () => {
        const v = new CurveExtractorVertex(1, 2, 3, 4);
        expect([v.xNumer, v.xDenom, v.yNumer, v.yDenom]).toEqual([1, 2, 3, 4]);
    });

    it('compares equivalent fractions as equal', () => {
        const a = new CurveExtractorVertex(1, 2, 3, 4);
        const b = new CurveExtractorVertex(2, 4, 6, 8);
        const c = new CurveExtractorVertex(-1, -2, -3, -4);
        expect(a.equals(b)).toBe(true);
        expect(b.equals(a)).toBe(true);
        expect(a.equals(c)).toBe(true);
        expect(a.lessThan(b)).toBe(false);
        expect(b.lessThan(a)).toBe(false);
    });

    it('orders lexicographically by x then y', () => {
        const a = new CurveExtractorVertex(1, 2, 1, 2);
        const b = new CurveExtractorVertex(2, 3, 0, 1);
        const c = new CurveExtractorVertex(1, 2, 3, 4);
        // 1/2 < 2/3, so a < b regardless of y.
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(a)).toBe(false);
        // Equal x, so the y components decide: 1/2 < 3/4.
        expect(a.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
    });

    it('is a strict weak ordering on a set of random rationals', () => {
        const vertices: CurveExtractorVertex[] = [];
        let state = 13579;
        for (let i = 0; i < 40; ++i) {
            const next = (m: number) => {
                state = (state * 1103515245 + 12345) & 0x7FFFFFFF;
                return 1 + (state % m);
            };
            vertices.push(new CurveExtractorVertex(next(7), next(5), next(7), next(5)));
        }
        for (const a of vertices) {
            expect(a.lessThan(a)).toBe(false);
            for (const b of vertices) {
                // Antisymmetry, and equality is exactly the incomparable case.
                expect(a.lessThan(b) && b.lessThan(a)).toBe(false);
                expect(a.equals(b)).toBe(!a.lessThan(b) && !b.lessThan(a));
                for (const c of vertices) {
                    if (a.lessThan(b) && b.lessThan(c)) {
                        expect(a.lessThan(c)).toBe(true);
                    }
                }
            }
        }
    });

    it('compares exactly for products beyond 2^53', () => {
        // 2^53 and 2^53 + 1 are distinct as int64 products, but the double
        // products are both 2^53. The bigint fallback keeps the comparison
        // exact.
        const a = new CurveExtractorVertex(2 ** 53, 1, 0, 1);
        const b = new CurveExtractorVertex(2 ** 53 + 2, 1, 0, 1);
        expect(a.lessThan(b)).toBe(true);
        expect(a.equals(b)).toBe(false);
    });
});

describe('CurveExtractorEdge', () => {
    it('sorts its two vertex indices', () => {
        expect(new CurveExtractorEdge(5, 2).v).toEqual([2, 5]);
        expect(new CurveExtractorEdge(2, 5).v).toEqual([2, 5]);
    });

    it('compares by index pair', () => {
        const a = new CurveExtractorEdge(0, 1);
        const b = new CurveExtractorEdge(0, 2);
        const c = new CurveExtractorEdge(1, 2);
        expect(a.equals(new CurveExtractorEdge(1, 0))).toBe(true);
        expect(a.lessThan(b)).toBe(true);
        expect(b.lessThan(c)).toBe(true);
        expect(c.lessThan(a)).toBe(false);
        expect(a.lessThan(a)).toBe(false);
    });
});

describe('CurveExtractor construction', () => {
    it('rejects bounds smaller than 2', () => {
        expect(() => new TestExtractor(1, 3, new Array<number>(3).fill(0)))
            .toThrow('Invalid input.');
        expect(() => new TestExtractor(3, 1, new Array<number>(3).fill(0)))
            .toThrow('Invalid input.');
    });

    it('rejects an undersized pixel array', () => {
        expect(() => new TestExtractor(3, 3, new Array<number>(8).fill(0)))
            .toThrow('Invalid input.');
    });

    it('allocates xBound * yBound working pixels', () => {
        const extractor = new TestExtractor(4, 5, new Array<number>(20).fill(0));
        expect(extractor.getPixelCount()).toBe(20);
    });
});

describe('CurveExtractor.addVertex / addEdge', () => {
    it('appends a vertex with normalized signs', () => {
        const extractor = makeExtractor([]);
        const vertices: CurveExtractorVertex[] = [];
        extractor.addVertexPublic(vertices, 1, -2, 3, 4);
        expect(vertices.length).toBe(1);
        expect(vertices[0].xNumer).toBe(-1);
        expect(vertices[0].xDenom).toBe(2);
    });

    it('appends an edge that references the two new vertices', () => {
        const extractor = makeExtractor([[0, 1, 1, 2, 1, 1, 1, 2]]);
        const result = extractor.extract(0);
        expect(result.vertices.length).toBe(2);
        expect(result.edges.length).toBe(1);
        expect(result.edges[0].v).toEqual([0, 1]);
        expect(result.vertices[0].xNumer / result.vertices[0].xDenom).toBe(0);
        expect(result.vertices[1].xNumer / result.vertices[1].xDenom).toBe(1);
    });
});

describe('CurveExtractor.convert', () => {
    it('divides the numerators by the denominators', () => {
        const extractor = makeExtractor([]);
        const converted = extractor.convert([
            new CurveExtractorVertex(1, 2, 3, 4),
            new CurveExtractorVertex(-5, 10, 0, 1)
        ]);
        expect(converted).toEqual([[0.5, 0.75], [-0.5, 0]]);
    });

    it('returns an empty array for empty input', () => {
        expect(makeExtractor([]).convert([])).toEqual([]);
    });
});

describe('CurveExtractor.makeUnique', () => {
    it('leaves empty inputs alone', () => {
        const extractor = makeExtractor([]);
        const vertices: CurveExtractorVertex[] = [];
        const edges: CurveExtractorEdge[] = [];
        extractor.makeUnique(vertices, edges);
        expect(vertices.length).toBe(0);
        expect(edges.length).toBe(0);
    });

    it('merges duplicated vertices of a shared endpoint', () => {
        // Two edges of a polyline that share the middle vertex (1/1, 0/1).
        // The extraction emits the shared vertex twice.
        const extractor = makeExtractor([
            [0, 1, 0, 1, 1, 1, 0, 1],
            [1, 1, 0, 1, 2, 1, 0, 1]
        ]);
        const result = extractor.extract(0);
        expect(result.vertices.length).toBe(4);
        expect(result.edges.length).toBe(2);

        extractor.makeUnique(result.vertices, result.edges);
        expect(result.vertices.length).toBe(3);
        expect(result.edges.length).toBe(2);

        const points = extractor.convert(result.vertices);
        // The unique vertices are numbered in first-occurrence order.
        expect(points).toEqual([[0, 0], [1, 0], [2, 0]]);
        expect(result.edges.map(e => e.v)).toEqual([[0, 1], [1, 2]]);
    });

    it('recognizes equivalent fractions as the same vertex', () => {
        // (2/2, 0/1) is the same point as (1/1, 0/1).
        const extractor = makeExtractor([
            [0, 1, 0, 1, 1, 1, 0, 1],
            [2, 2, 0, 3, 2, 1, 0, 1]
        ]);
        const result = extractor.extract(0);
        extractor.makeUnique(result.vertices, result.edges);
        expect(result.vertices.length).toBe(3);
        expect(extractor.convert(result.vertices)).toEqual([[0, 0], [1, 0], [2, 0]]);
    });

    it('removes duplicated edges', () => {
        // The same edge emitted twice with the same orientation collapses.
        const extractor = makeExtractor([
            [0, 1, 0, 1, 1, 1, 0, 1],
            [0, 1, 0, 1, 1, 1, 0, 1]
        ]);
        const result = extractor.extract(0);
        expect(result.edges.length).toBe(2);
        extractor.makeUnique(result.vertices, result.edges);
        expect(result.vertices.length).toBe(2);
        expect(result.edges.length).toBe(1);
        expect(result.edges[0].v).toEqual([0, 1]);
    });

    it('preserves the upstream quirk of not reordering remapped edges', () => {
        // Upstream replaces the edge endpoints in place and does not restore
        // the increasing order, so <1,0> and <0,1> are distinct map keys even
        // though they are the same undirected edge. Here the second edge is
        // emitted with its endpoints in the opposite order, which remaps to
        // (1,0) and therefore survives as a separate edge.
        const extractor = makeExtractor([
            [0, 1, 0, 1, 1, 1, 0, 1],
            [1, 1, 0, 1, 0, 1, 0, 1]
        ]);
        const result = extractor.extract(0);
        extractor.makeUnique(result.vertices, result.edges);
        expect(result.vertices.length).toBe(2);
        expect(result.edges.length).toBe(2);
        expect(result.edges.map(e => e.v)).toEqual([[0, 1], [1, 0]]);
    });
});

describe('CurveExtractor.extractReal', () => {
    it('optionally removes duplicate vertices', () => {
        const specs = [
            [0, 1, 0, 1, 1, 1, 0, 1],
            [1, 1, 0, 1, 2, 1, 0, 1]
        ];

        const withDuplicates = makeExtractor(specs).extractReal(0, false);
        expect(withDuplicates.vertices).toEqual([[0, 0], [1, 0], [1, 0], [2, 0]]);
        expect(withDuplicates.edges.map(e => e.v)).toEqual([[0, 1], [2, 3]]);

        const unique = makeExtractor(specs).extractReal(0, true);
        expect(unique.vertices).toEqual([[0, 0], [1, 0], [2, 0]]);
        expect(unique.edges.map(e => e.v)).toEqual([[0, 1], [1, 2]]);
    });

    it('keeps every edge endpoint a valid index into the vertices', () => {
        const extractor = makeExtractor([
            [0, 1, 0, 1, 1, 2, 1, 2],
            [1, 2, 1, 2, 1, 1, 1, 1],
            [1, 1, 1, 1, 0, 1, 0, 1]
        ]);
        const result = extractor.extractReal(0, true);
        expect(result.vertices.length).toBe(3);
        for (const edge of result.edges) {
            for (const v of edge.v) {
                expect(v).toBeGreaterThanOrEqual(0);
                expect(v).toBeLessThan(result.vertices.length);
            }
        }
    });
});
