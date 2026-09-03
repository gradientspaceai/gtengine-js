import { describe, it, expect } from 'vitest';
import { check, fc } from './helpers/arbitraries.js';
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

describe('CurveExtractor verification', () => {
    // ---- exact rational reference ----------------------------------------
    /** Compare n0/d0 with n1/d1 exactly (denominators strictly positive). */
    function cmpRational(n0: number, d0: number, n1: number, d1: number): number {
        const a = BigInt(n0) * BigInt(d1);
        const b = BigInt(n1) * BigInt(d0);
        return a < b ? -1 : (a > b ? 1 : 0);
    }
    function cmpVertex(a: CurveExtractorVertex, b: CurveExtractorVertex): number {
        const cx = cmpRational(a.xNumer, a.xDenom, b.xNumer, b.xDenom);
        if (cx !== 0) { return cx; }
        return cmpRational(a.yNumer, a.yDenom, b.yNumer, b.yDenom);
    }
    /** A canonical string for the rational pair, computed independently. */
    function canonical(v: CurveExtractorVertex): string {
        const red = (n: number, d: number) => {
            let x = BigInt(Math.abs(n)), y = BigInt(Math.abs(d));
            while (y !== 0n) { const t = x % y; x = y; y = t; }
            const g = x === 0n ? 1n : x;
            return `${BigInt(n) / g}/${BigInt(d) / g}`;
        };
        return `${red(v.xNumer, v.xDenom)},${red(v.yNumer, v.yDenom)}`;
    }

    const smallInt = fc.integer({ min: -8, max: 8 });
    const nonzeroInt = smallInt.filter(x => x !== 0);
    const vertexArb = fc.tuple(smallInt, nonzeroInt, smallInt, nonzeroInt)
        .map(([xn, xd, yn, yd]) => new CurveExtractorVertex(xn, xd, yn, yd));
    // Large components, to exercise the bigint fallback in compareProducts.
    const bigVertexArb = fc.tuple(
        fc.integer({ min: -(2 ** 40), max: 2 ** 40 }),
        fc.integer({ min: 1, max: 2 ** 40 }),
        fc.integer({ min: -(2 ** 40), max: 2 ** 40 }),
        fc.integer({ min: 1, max: 2 ** 40 }))
        .map(([xn, xd, yn, yd]) => new CurveExtractorVertex(xn, xd, yn, yd));

    // ---- CurveExtractorVertex --------------------------------------------
    it('the constructor makes both denominators positive without changing the value', () => {
        check(fc.tuple(smallInt, nonzeroInt, smallInt, nonzeroInt), ([xn, xd, yn, yd]) => {
            const v = new CurveExtractorVertex(xn, xd, yn, yd);
            return v.xDenom > 0 && v.yDenom > 0
                && cmpRational(v.xNumer, v.xDenom, xn * Math.sign(xd), Math.abs(xd)) === 0
                && cmpRational(v.yNumer, v.yDenom, yn * Math.sign(yd), Math.abs(yd)) === 0;
        });
    });

    it('the default constructor produces the upstream default (all components zero)', () => {
        const v = new CurveExtractorVertex();
        expect([v.xNumer, v.xDenom, v.yNumer, v.yDenom]).toEqual([0, 0, 0, 0]);
        expect(Object.is(v.xNumer, -0)).toBe(false);
        expect(Object.is(v.xDenom, -0)).toBe(false);
    });

    it('equals and lessThan match exact rational comparison', () => {
        check(fc.tuple(vertexArb, vertexArb), ([a, b]) => {
            const c = cmpVertex(a, b);
            return a.lessThan(b) === (c < 0)
                && b.lessThan(a) === (c > 0)
                && a.equals(b) === (c === 0);
        });
        // Same, with components that overflow the exact double product.
        check(fc.tuple(bigVertexArb, bigVertexArb), ([a, b]) => {
            const c = cmpVertex(a, b);
            return a.lessThan(b) === (c < 0) && a.equals(b) === (c === 0);
        });
    });

    it('lessThan is a strict weak ordering', () => {
        check(fc.array(vertexArb, { minLength: 3, maxLength: 3 }), ([a, b, c]) => {
            // Irreflexive.
            if (a!.lessThan(a!)) { return false; }
            // Asymmetric.
            if (a!.lessThan(b!) && b!.lessThan(a!)) { return false; }
            // Transitive.
            if (a!.lessThan(b!) && b!.lessThan(c!) && !a!.lessThan(c!)) { return false; }
            // Equivalence (neither less) is transitive too.
            const eq = (p: CurveExtractorVertex, q: CurveExtractorVertex) =>
                !p.lessThan(q) && !q.lessThan(p);
            if (eq(a!, b!) && eq(b!, c!) && !eq(a!, c!)) { return false; }
            // The equivalence induced by lessThan is exactly equals().
            return eq(a!, b!) === a!.equals(b!);
        });
    });

    it('sorting with lessThan orders the rational values', () => {
        check(fc.array(vertexArb, { minLength: 2, maxLength: 10 }), vs => {
            const sorted = [...vs].sort((p, q) => (p.lessThan(q) ? -1 : (q.lessThan(p) ? 1 : 0)));
            for (let i = 1; i < sorted.length; ++i) {
                if (cmpVertex(sorted[i - 1]!, sorted[i]!) > 0) { return false; }
            }
            return true;
        });
    });

    it('the map key used by makeUnique agrees with the upstream equivalence', () => {
        // Upstream keys a std::map by Vertex, so two vertices share a slot iff
        // neither is less than the other. The port keys by a reduced-fraction
        // string; the two must induce the same partition.
        check(fc.tuple(vertexArb, vertexArb), ([a, b]) =>
            (canonical(a) === canonical(b)) === a.equals(b));
    });

    // ---- CurveExtractorEdge ----------------------------------------------
    const edgeArb = fc.tuple(fc.integer({ min: 0, max: 6 }), fc.integer({ min: 0, max: 6 }))
        .map(([a, b]) => new CurveExtractorEdge(a, b));

    it('the edge constructor stores the indices in increasing order', () => {
        check(fc.tuple(fc.integer({ min: -5, max: 5 }), fc.integer({ min: -5, max: 5 })),
            ([a, b]) => {
                const e = new CurveExtractorEdge(a, b);
                return e.v[0] <= e.v[1] && e.v[0] === Math.min(a, b) && e.v[1] === Math.max(a, b);
            });
    });

    it('edge lessThan is the lexicographic order and equals is its equivalence', () => {
        check(fc.tuple(edgeArb, edgeArb), ([a, b]) => {
            const lex = a.v[0] !== b.v[0] ? (a.v[0] < b.v[0] ? -1 : 1)
                : (a.v[1] !== b.v[1] ? (a.v[1] < b.v[1] ? -1 : 1) : 0);
            return a.lessThan(b) === (lex < 0) && b.lessThan(a) === (lex > 0)
                && a.equals(b) === (lex === 0);
        });
    });

    // ---- addVertex / addEdge / convert ------------------------------------
    it('addEdge appends two vertices and one edge referencing them', () => {
        check(fc.array(fc.tuple(smallInt, nonzeroInt, smallInt, nonzeroInt,
            smallInt, nonzeroInt, smallInt, nonzeroInt), { minLength: 1, maxLength: 6 }),
            specs => {
                const extractor = makeExtractor(specs.map(s => [...s]));
                const { vertices, edges } = extractor.extract(0);
                if (vertices.length !== 2 * specs.length) { return false; }
                if (edges.length !== specs.length) { return false; }
                return edges.every((e, t) => e.v[0] === 2 * t && e.v[1] === 2 * t + 1);
            });
    });

    it('convert divides each numerator by its denominator', () => {
        check(fc.array(vertexArb, { minLength: 0, maxLength: 10 }), vs => {
            const extractor = makeExtractor([]);
            const out = extractor.convert(vs);
            return out.length === vs.length && out.every((p, i) =>
                Object.is(p[0], vs[i]!.xNumer / vs[i]!.xDenom)
                && Object.is(p[1], vs[i]!.yNumer / vs[i]!.yDenom));
        });
    });

    // ---- makeUnique -------------------------------------------------------
    /** A random extraction: a soup of edges over a small set of rationals. */
    const specsArb = fc.array(
        fc.tuple(fc.integer({ min: 0, max: 3 }), fc.integer({ min: 1, max: 2 }),
            fc.integer({ min: 0, max: 3 }), fc.integer({ min: 1, max: 2 }),
            fc.integer({ min: 0, max: 3 }), fc.integer({ min: 1, max: 2 }),
            fc.integer({ min: 0, max: 3 }), fc.integer({ min: 1, max: 2 })),
        { minLength: 1, maxLength: 12 }).map(a => a.map(s => [...s]));

    it('makeUnique produces distinct vertices and valid edge indices', () => {
        check(specsArb, specs => {
            const extractor = makeExtractor(specs);
            const { vertices, edges } = extractor.extract(0);
            const before = vertices.map(canonical);
            extractor.makeUnique(vertices, edges);
            // Postcondition: the vertices are pairwise distinct rationals and
            // are exactly the distinct rationals of the input.
            expect(new Set(vertices.map(canonical)).size).toBe(vertices.length);
            expect(new Set(vertices.map(canonical))).toEqual(new Set(before));
            // Every edge endpoint indexes the packed array.
            for (const e of edges) {
                for (const i of e.v) {
                    expect(i).toBeGreaterThanOrEqual(0);
                    expect(i).toBeLessThan(vertices.length);
                }
            }
            return true;
        });
    });

    it('makeUnique preserves the set of undirected geometric edges', () => {
        check(specsArb, specs => {
            const extractor = makeExtractor(specs);
            const { vertices, edges } = extractor.extract(0);
            const geoBefore = new Set(edges.map(e => {
                const a = canonical(vertices[e.v[0]]!), b = canonical(vertices[e.v[1]]!);
                return a <= b ? `${a}|${b}` : `${b}|${a}`;
            }));
            extractor.makeUnique(vertices, edges);
            const geoAfter = new Set(edges.map(e => {
                const a = canonical(vertices[e.v[0]]!), b = canonical(vertices[e.v[1]]!);
                return a <= b ? `${a}|${b}` : `${b}|${a}`;
            }));
            expect(geoAfter).toEqual(geoBefore);
            return true;
        });
    });

    it('makeUnique is idempotent', () => {
        check(specsArb, specs => {
            const extractor = makeExtractor(specs);
            const { vertices, edges } = extractor.extract(0);
            extractor.makeUnique(vertices, edges);
            const v1 = vertices.map(canonical);
            const e1 = edges.map(e => [...e.v]);
            extractor.makeUnique(vertices, edges);
            return JSON.stringify(vertices.map(canonical)) === JSON.stringify(v1)
                && JSON.stringify(edges.map(e => [...e.v])) === JSON.stringify(e1);
        });
    });

    it('makeUnique leaves the arrays untouched when either is empty', () => {
        check(fc.array(vertexArb, { minLength: 0, maxLength: 4 }), vs => {
            const extractor = makeExtractor([]);
            const vertices = [...vs];
            const edges: CurveExtractorEdge[] = [];
            extractor.makeUnique(vertices, edges);
            return vertices.length === vs.length && edges.length === 0;
        });
    });

    it('the surviving duplicate edges are exactly those hit by the upstream quirk', () => {
        // Upstream writes the new indices into the edge without restoring the
        // increasing order, so <a,b> and <b,a> are distinct map keys. The
        // number of output edges therefore equals the number of distinct
        // ORDERED remapped pairs, not the number of distinct undirected ones.
        check(specsArb, specs => {
            const extractor = makeExtractor(specs);
            const { vertices, edges } = extractor.extract(0);
            // Independent computation of the remapped ordered pairs.
            const index = new Map<string, number>();
            for (const v of vertices) {
                const k = canonical(v);
                if (!index.has(k)) { index.set(k, index.size); }
            }
            const ordered = new Set(edges.map(e =>
                `${index.get(canonical(vertices[e.v[0]]!))},`
                + `${index.get(canonical(vertices[e.v[1]]!))}`));
            extractor.makeUnique(vertices, edges);
            expect(edges.length).toBe(ordered.size);
            expect(new Set(edges.map(e => `${e.v[0]},${e.v[1]}`))).toEqual(ordered);
            return true;
        });
    });

    // ---- extractReal ------------------------------------------------------
    it('extractReal equals extract followed by makeUnique and convert', () => {
        check(fc.tuple(specsArb, fc.boolean()), ([specs, unique]) => {
            const a = makeExtractor(specs);
            const real = a.extractReal(0, unique);

            const b = makeExtractor(specs);
            const raw = b.extract(0);
            if (unique) { b.makeUnique(raw.vertices, raw.edges); }
            const expected = b.convert(raw.vertices);

            expect(real.vertices).toEqual(expected);
            expect(real.edges.map(e => [...e.v])).toEqual(raw.edges.map(e => [...e.v]));
            // Every endpoint indexes the converted array.
            for (const e of real.edges) {
                for (const i of e.v) {
                    expect(i).toBeGreaterThanOrEqual(0);
                    expect(i).toBeLessThan(real.vertices.length);
                }
            }
            return true;
        });
    });

    it('the constructor rejects degenerate bounds and short pixel arrays', () => {
        check(fc.tuple(fc.integer({ min: -2, max: 4 }), fc.integer({ min: -2, max: 4 })),
            ([x, y]) => {
                const pixels = new Array<number>(Math.max(0, x * y)).fill(0);
                const build = () => new TestExtractor(x, y, pixels);
                if (x > 1 && y > 1) {
                    expect(build).not.toThrow();
                    expect(build().getPixelCount()).toBe(x * y);
                } else {
                    expect(build).toThrow();
                }
                return true;
            });
    });
});
