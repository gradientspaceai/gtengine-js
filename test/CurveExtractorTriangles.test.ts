import { describe, it, expect } from 'vitest';
import { CurveExtractorTriangles } from '../src/CurveExtractorTriangles.js';
import { CurveExtractorSquares } from '../src/CurveExtractorSquares.js';
import { CurveExtractorVertex } from '../src/CurveExtractor.js';
import { check, fc } from './helpers/arbitraries.js';

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

// ---------------------------------------------------------------------------
// Verification block (V16).
// ---------------------------------------------------------------------------

// A canonical string for the reduced rational pair of a vertex, used to count
// distinct vertex positions.
function reducedKeyT(v: CurveExtractorVertex): string {
    const gcd = (a: number, b: number): number => {
        let x = Math.abs(a), y = Math.abs(b);
        while (y !== 0) { const t = x % y; x = y; y = t; }
        return x || 1;
    };
    const gx = gcd(v.xNumer, v.xDenom);
    const gy = gcd(v.yNumer, v.yDenom);
    return `${v.xNumer / gx}/${v.xDenom / gx},${v.yNumer / gy}/${v.yDenom / gy}`;
}

// Integer images of the form f(x,y) = a*x + b*y + c. All products used by the
// extractor stay far below 2^53, so the rational arithmetic is exact.
const linearImageT = fc.tuple(
    fc.integer({ min: -4, max: 4 }), fc.integer({ min: -4, max: 4 }),
    fc.integer({ min: -20, max: 20 }),
    fc.integer({ min: 2, max: 8 }), fc.integer({ min: 2, max: 8 }))
    .map(([a, b, c, xBound, yBound]) => ({
        a, b, c, xBound, yBound,
        pixels: image(xBound, yBound, (x, y) => a * x + b * y + c)
    }));

// An image of +-1 pixels with a +1 border, so the level set is a union of
// closed curves strictly inside the image.
const signImageT = fc.tuple(fc.integer({ min: 5, max: 8 }),
    fc.integer({ min: 5, max: 8 }),
    fc.array(fc.boolean(), { minLength: 64, maxLength: 64 }))
    .map(([xBound, yBound, flags]) => ({
        xBound, yBound,
        pixels: image(xBound, yBound, (x, y) => {
            if (x === 0 || y === 0 || x === xBound - 1 || y === yBound - 1) {
                return 1;
            }
            return flags[(x + xBound * y) % flags.length] ? 1 : -1;
        })
    }));

describe('CurveExtractorTriangles verification', () => {
    it('places every vertex of a linear image exactly on the level line', () => {
        check(linearImageT, ({ a, b, c, xBound, yBound, pixels }) => {
            const r = new CurveExtractorTriangles(xBound, yBound, pixels).extract(0);
            for (const v of r.vertices) {
                // The constructor normalizes the sign so both denominators
                // are positive.
                expect(v.xDenom).toBeGreaterThan(0);
                expect(v.yDenom).toBeGreaterThan(0);
                expect(onLine(v, a, b, c)).toBe(true);
                const [px, py] = asPair(v);
                expect(px).toBeGreaterThanOrEqual(0);
                expect(px).toBeLessThanOrEqual(xBound - 1);
                expect(py).toBeGreaterThanOrEqual(0);
                expect(py).toBeLessThanOrEqual(yBound - 1);
            }
            for (const e of r.edges) {
                expect(e.v[0]).toBeGreaterThanOrEqual(0);
                expect(e.v[1]).toBeLessThan(r.vertices.length);
                // The Edge constructor stores the pair in increasing order.
                expect(e.v[0]).toBeLessThanOrEqual(e.v[1]);
            }
        });
    });

    it('extracts level L from the image exactly as level 0 from the shifted image', () => {
        check(fc.tuple(linearImageT, fc.integer({ min: -20, max: 20 })),
            ([{ xBound, yBound, pixels }, level]) => {
                const a = new CurveExtractorTriangles(xBound, yBound, pixels)
                    .extract(level);
                const shifted = pixels.map(p => p - level);
                const b = new CurveExtractorTriangles(xBound, yBound, shifted)
                    .extract(0);
                expect(a.vertices.map(reducedKeyT))
                    .toEqual(b.vertices.map(reducedKeyT));
                expect(a.edges.map(e => e.v.slice()))
                    .toEqual(b.edges.map(e => e.v.slice()));
            });
    });

    it('is invariant when the image is negated', () => {
        // The case analysis converts every sign pattern to one starting with
        // a positive corner, so the extracted curve is the same set.
        check(signImageT, ({ xBound, yBound, pixels }) => {
            const a = new CurveExtractorTriangles(xBound, yBound, pixels).extract(0);
            const b = new CurveExtractorTriangles(xBound, yBound,
                pixels.map(p => -p)).extract(0);
            expect(new Set(b.vertices.map(reducedKeyT)))
                .toEqual(new Set(a.vertices.map(reducedKeyT)));
        });
    });

    it('produces closed level curves for an image with a positive border', () => {
        // Every pixel is nonzero and the border is positive, so the level set
        // is a union of closed curves in the interior of the image. In the
        // deduplicated graph every vertex must then have even degree.
        check(signImageT, ({ xBound, yBound, pixels }) => {
            const extractor = new CurveExtractorTriangles(xBound, yBound, pixels);
            const r = extractor.extract(0);
            extractor.makeUnique(r.vertices, r.edges);
            const degree = new Array<number>(r.vertices.length).fill(0);
            for (const e of r.edges) {
                ++degree[e.v[0]];
                ++degree[e.v[1]];
            }
            for (let i = 0; i < degree.length; ++i) {
                expect(degree[i]).toBeGreaterThan(0);
                expect(degree[i] % 2).toBe(0);
            }
        });
    });

    it('collapses exactly the coincident vertices in makeUnique', () => {
        check(signImageT, ({ xBound, yBound, pixels }) => {
            const extractor = new CurveExtractorTriangles(xBound, yBound, pixels);
            const r = extractor.extract(0);
            const before = r.vertices.map(reducedKeyT);
            const distinct = new Set(before);
            const oldEdges = r.edges.map(e => e.v.slice());
            extractor.makeUnique(r.vertices, r.edges);

            expect(r.vertices.length).toBe(distinct.size);
            // First-occurrence numbering: vertex i of the packed array is the
            // i-th distinct key in the original order.
            const order: string[] = [];
            const seen = new Set<string>();
            for (const key of before) {
                if (!seen.has(key)) { seen.add(key); order.push(key); }
            }
            expect(r.vertices.map(reducedKeyT)).toEqual(order);

            // Every remapped edge joins the images of its old endpoints, and
            // the packed edges are the distinct ordered remapped pairs (the
            // upstream MakeUnique does not reorder after remapping, which is
            // upstream issue "duplicate edges when remapped vertices reverse
            // order").
            const index = new Map<string, number>();
            order.forEach((key, i) => index.set(key, i));
            const remapped = oldEdges.map(e => [
                index.get(before[e[0]]) as number,
                index.get(before[e[1]]) as number
            ]);
            const uniqueOrdered: number[][] = [];
            const seenEdges = new Set<string>();
            for (const e of remapped) {
                const key = `${e[0]},${e[1]}`;
                if (!seenEdges.has(key)) { seenEdges.add(key); uniqueOrdered.push(e); }
            }
            expect(r.edges.map(e => e.v.slice())).toEqual(uniqueOrdered);
        });
    });

    it('converts the rational vertices to the quotient in extractReal', () => {
        check(linearImageT, ({ xBound, yBound, pixels }) => {
            const extractor = new CurveExtractorTriangles(xBound, yBound, pixels);
            const rational = extractor.extract(0);
            const real = new CurveExtractorTriangles(xBound, yBound, pixels)
                .extractReal(0, false);
            expect(real.vertices.length).toBe(rational.vertices.length);
            for (let i = 0; i < real.vertices.length; ++i) {
                expect(real.vertices[i][0])
                    .toBe(rational.vertices[i].xNumer / rational.vertices[i].xDenom);
                expect(real.vertices[i][1])
                    .toBe(rational.vertices[i].yNumer / rational.vertices[i].yDenom);
            }
        });
    });

    it('emits nothing when the level set misses the image', () => {
        check(fc.tuple(fc.integer({ min: 2, max: 8 }),
            fc.integer({ min: 2, max: 8 }), fc.integer({ min: 1, max: 9 })),
        ([xBound, yBound, value]) => {
            const pixels = image(xBound, yBound, () => value);
            const r = new CurveExtractorTriangles(xBound, yBound, pixels).extract(0);
            expect(r.vertices).toEqual([]);
            expect(r.edges).toEqual([]);
            // Extracting at that constant value gives the all-zero case: the
            // three boundary edges of each of the two triangles of a square.
            const all = new CurveExtractorTriangles(xBound, yBound, pixels)
                .extract(value);
            expect(all.edges.length)
                .toBe(6 * (xBound - 1) * (yBound - 1));
        });
    });

    it('rejects images smaller than a single square', () => {
        check(fc.tuple(fc.integer({ min: 0, max: 1 }),
            fc.integer({ min: 0, max: 8 })), ([small, other]) => {
            const size = Math.max(small * other, 1);
            expect(() => new CurveExtractorTriangles(small, other,
                new Array(size).fill(0))).toThrow('Invalid input.');
            expect(() => new CurveExtractorTriangles(other, small,
                new Array(size).fill(0))).toThrow('Invalid input.');
        });
    });
});

describe('CurveExtractorTriangles vs CurveExtractorSquares verification', () => {
    it('contains every vertex the squares extractor finds on a linear image', () => {
        // A linear image is linear on each of the two triangles of a square,
        // so both algorithms cross the same grid edges at the same points. The
        // triangles extractor additionally crosses the square diagonal, so its
        // vertex set is a superset of the squares extractor's.
        check(linearImageT, ({ a, b, c, xBound, yBound, pixels }) => {
            const t = new CurveExtractorTriangles(xBound, yBound, pixels);
            const rt = t.extract(0);
            t.makeUnique(rt.vertices, rt.edges);

            const s = new CurveExtractorSquares(xBound, yBound, pixels);
            const rs = s.extract(0);
            s.makeUnique(rs.vertices, rs.edges);

            const fromTriangles = new Set(rt.vertices.map(reducedKeyT));
            for (const v of rs.vertices) {
                expect(fromTriangles.has(reducedKeyT(v))).toBe(true);
            }
            // Both sets lie on the exact level line.
            for (const v of rt.vertices) { expect(onLine(v, a, b, c)).toBe(true); }
            for (const v of rs.vertices) { expect(onLine(v, a, b, c)).toBe(true); }
        });
    });
});
