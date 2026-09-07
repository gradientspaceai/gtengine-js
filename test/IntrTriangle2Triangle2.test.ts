import { describe, expect, it } from 'vitest';
import {
    IntrTriangle2Triangle2FI,
    IntrTriangle2Triangle2TI,
    defaultIntrTriangle2Triangle2FIResult,
    defaultIntrTriangle2Triangle2TIResult
} from '../src/IntrTriangle2Triangle2.js';
import { Triangle } from '../src/Triangle.js';
import { Vector, add, dot, mul, sub } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame,
    seededRandom, wellScaledVector
} from './helpers/arbitraries.js';
import { dotPerp } from '../src/Vector2.js';

const V2 = (x: number, y: number) => Vector.fromArray([x, y]);

function tri(a: readonly [number, number], b: readonly [number, number],
    c: readonly [number, number]): Triangle {
    return Triangle.fromVertices(V2(a[0], a[1]), V2(b[0], b[1]),
        V2(c[0], c[1]));
}

function translate(t: Triangle, d: Vector): Triangle {
    return Triangle.fromVertexArray(t.v.map(v => add(v, d)));
}

// Counterclockwise ordering is required by both queries.
function makeCcw(t: Triangle): Triangle {
    const area2 = dotPerp(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0]));
    return area2 >= 0 ? t
        : Triangle.fromVertices(t.v[0], t.v[2], t.v[1]);
}

// The signed area of a polygon, positive for counterclockwise ordering.
function polygonArea(poly: readonly Vector[]): number {
    let sum = 0;
    for (let i = poly.length - 1, j = 0; j < poly.length; i = j++) {
        sum += poly[i].values[0] * poly[j].values[1]
            - poly[j].values[0] * poly[i].values[1];
    }
    return sum / 2;
}

// The signed barycentric-style test: a point is inside (or on) a
// counterclockwise triangle when it is left of, or on, every edge.
function insideTriangle(p: Vector, t: Triangle, tol: number): boolean {
    for (let i = 0; i < 3; ++i) {
        const j = (i + 1) % 3;
        const e = sub(t.v[j], t.v[i]);
        if (dotPerp(e, sub(p, t.v[i])) < -tol) {
            return false;
        }
    }
    return true;
}

describe('IntrTriangle2Triangle2 default results', () => {
    it('match the upstream default constructors', () => {
        expect(defaultIntrTriangle2Triangle2TIResult())
            .toEqual({ intersect: false });
        expect(defaultIntrTriangle2Triangle2FIResult().intersection)
            .toEqual([]);
    });
});

describe('IntrTriangle2Triangle2TI known configurations', () => {
    const ti = new IntrTriangle2Triangle2TI();
    const unit = tri([0, 0], [1, 0], [0, 1]);

    it('reports a triangle intersecting itself', () => {
        expect(ti.test(unit, unit).intersect).toBe(true);
    });

    it('rejects separated triangles', () => {
        expect(ti.test(unit, translate(unit, V2(5, 0))).intersect).toBe(false);
        expect(ti.test(unit, translate(unit, V2(0, -5))).intersect).toBe(false);
        expect(ti.test(unit, translate(unit, V2(-3, 3))).intersect).toBe(false);
    });

    it('reports overlap for partially overlapping triangles', () => {
        expect(ti.test(unit, translate(unit, V2(0.5, 0.1))).intersect)
            .toBe(true);
    });

    it('reports containment', () => {
        const big = tri([-5, -5], [5, -5], [0, 5]);
        const small = tri([0, 0], [0.2, 0], [0, 0.2]);
        expect(ti.test(big, small).intersect).toBe(true);
        expect(ti.test(small, big).intersect).toBe(true);
    });

    it('treats edge-touching triangles as separated', () => {
        // The mirrored triangle shares the edge from (1,0) to (0,1). The
        // upstream WhichSide convention counts a projection interval [0,b]
        // with b > 0 as lying on the positive side, so a contact of measure
        // zero is reported as separation.
        const mirror = tri([1, 1], [0, 1], [1, 0]);
        expect(ti.test(unit, mirror).intersect).toBe(false);
        expect(ti.test(mirror, unit).intersect).toBe(false);
    });

    it('treats vertex-touching triangles as separated', () => {
        const other = tri([1, 0], [2, 0], [2, 1]);
        expect(ti.test(unit, other).intersect).toBe(false);
        expect(ti.test(other, unit).intersect).toBe(false);
    });

    it('reports intersection for an overlap of positive area', () => {
        // Push the mirrored triangle of the previous case slightly inward so
        // that the overlap has positive area.
        const mirror = tri([0.9, 0.9], [-0.1, 0.9], [0.9, -0.1]);
        expect(ti.test(unit, mirror).intersect).toBe(true);
        expect(ti.test(mirror, unit).intersect).toBe(true);
    });

    it('is symmetric in its arguments', () => {
        for (const d of [0.2, 0.9, 1.0, 1.1, 2.0]) {
            const t1 = translate(unit, V2(d, d * 0.3));
            expect(ti.test(unit, t1).intersect).toBe(ti.test(t1, unit)
                .intersect);
        }
    });
});

describe('IntrTriangle2Triangle2FI known configurations', () => {
    const fi = new IntrTriangle2Triangle2FI();
    const unit = tri([0, 0], [1, 0], [0, 1]);

    it('returns the triangle when intersected with itself', () => {
        const poly = fi.find(unit, unit).intersection;
        expect(poly.length).toBe(3);
        expect(polygonArea(poly)).toBeCloseTo(0.5, 12);
    });

    it('returns the contained triangle', () => {
        const big = tri([-5, -5], [5, -5], [0, 5]);
        const small = tri([0, 0], [0.2, 0], [0, 0.2]);
        const poly = fi.find(big, small).intersection;
        expect(poly.length).toBe(3);
        expect(polygonArea(poly)).toBeCloseTo(0.02, 12);
    });

    it('returns an empty set for separated triangles', () => {
        expect(fi.find(unit, translate(unit, V2(5, 0))).intersection)
            .toEqual([]);
    });

    it('computes a known quadrilateral overlap', () => {
        // Two axis-aligned right triangles whose overlap is the square
        // [0,1] x [0,1] minus nothing: the triangles are
        // {x>=0,y>=0,x+y<=3} and {x<=3,y<=3,x+y>=... }. Use a simpler pair:
        // the unit square split by two overlapping triangles.
        const t0 = tri([0, 0], [4, 0], [0, 4]);
        const t1 = tri([1, 1], [5, 1], [1, 5]);
        const poly = fi.find(t0, t1).intersection;
        // The overlap is the triangle with vertices (1,1), (3,1), (1,3).
        expect(poly.length).toBe(3);
        expect(polygonArea(poly)).toBeCloseTo(2, 12);
        const xs = poly.map(v => `${v.values[0]},${v.values[1]}`).sort();
        expect(xs).toEqual(['1,1', '1,3', '3,1']);
    });

    it('returns an empty set for edge-touching triangles', () => {
        // The clipping of the mirrored triangle against the line x + y = 1
        // leaves nothing on the positive side, so the query agrees with the
        // test-intersection query that a measure-zero contact is not an
        // intersection.
        const mirror = tri([1, 1], [0, 1], [1, 0]);
        expect(fi.find(unit, mirror).intersection).toEqual([]);
    });

    it('returns the sliver overlap of nearly edge-touching triangles', () => {
        const mirror = tri([0.99, 0.99], [-0.01, 0.99], [0.99, -0.01]);
        const poly = fi.find(unit, mirror).intersection;
        // The overlap is the hexagon with vertices (0.98,0), (0.99,0),
        // (0.99,0.01), (0.01,0.99), (0,0.99), (0,0.98): the sliver of the
        // unit triangle above x + y = 0.98, with the two corners beyond
        // x = 0.99 and y = 0.99 trimmed. Its area is
        // (0.5 - 0.98^2/2) - 2 * (0.01^2/2).
        expect(poly.length).toBe(6);
        expect(polygonArea(poly))
            .toBeCloseTo(0.5 - (0.98 * 0.98) / 2 - 0.01 * 0.01, 12);
        for (const v of poly) {
            expect(insideTriangle(v, unit, 1e-12)).toBe(true);
            expect(insideTriangle(v, mirror, 1e-12)).toBe(true);
        }
    });
});

describe('IntrTriangle2Triangle2 randomized cross-checks', () => {
    let seed = 13572468;
    const rand = () => {
        seed = (seed * 1103515245 + 12345) & 0x7fffffff;
        return seed / 0x7fffffff;
    };
    const rv = (s: number) => V2((rand() - 0.5) * s, (rand() - 0.5) * s);

    it('has TI and FI agreeing, with the clipped polygon inside both',
        () => {
            const ti = new IntrTriangle2Triangle2TI();
            const fi = new IntrTriangle2Triangle2FI();
            let hits = 0, misses = 0;
            for (let k = 0; k < 2000; ++k) {
                const t0 = makeCcw(Triangle.fromVertices(rv(4), rv(4), rv(4)));
                const t1 = makeCcw(Triangle.fromVertices(
                    add(rv(4), V2(1, 0)), add(rv(4), V2(1, 0)),
                    add(rv(4), V2(1, 0))));

                const tiHit = ti.test(t0, t1).intersect;
                const poly = fi.find(t0, t1).intersection;
                const fiHit = poly.length > 0;

                // The TI and FI queries agree except in configurations whose
                // intersection is a single point or a sliver, where the
                // clipping arithmetic can drop the degenerate result.
                if (fiHit) {
                    expect(tiHit).toBe(true);
                }

                if (fiHit) {
                    ++hits;
                    // Every clipped vertex must lie in both triangles.
                    for (const v of poly) {
                        expect(insideTriangle(v, t0, 1e-9)).toBe(true);
                        expect(insideTriangle(v, t1, 1e-9)).toBe(true);
                    }
                    // The clipped polygon is convex and counterclockwise.
                    expect(polygonArea(poly)).toBeGreaterThanOrEqual(-1e-9);
                } else {
                    ++misses;
                }
            }
            expect(hits).toBeGreaterThan(200);
            expect(misses).toBeGreaterThan(200);
        });

    it('reports intersection whenever the triangles share an interior point',
        () => {
            const ti = new IntrTriangle2Triangle2TI();
            const fi = new IntrTriangle2Triangle2FI();
            for (let k = 0; k < 800; ++k) {
                const q = rv(8);
                // Translate each triangle so that a strictly interior point
                // lands on the shared point q.
                const shift = (t: Triangle) => {
                    const w = [rand() + 0.3, rand() + 0.3, rand() + 0.3];
                    const sum = w[0] + w[1] + w[2];
                    let p = V2(0, 0);
                    for (let i = 0; i < 3; ++i) {
                        p = add(p, Vector.fromArray(
                            t.v[i].values.map(x => (x * w[i]) / sum)));
                    }
                    return Triangle.fromVertexArray(
                        t.v.map(v => add(v, sub(q, p))));
                };
                const t0 = shift(makeCcw(
                    Triangle.fromVertices(rv(3), rv(3), rv(3))));
                const t1 = shift(makeCcw(
                    Triangle.fromVertices(rv(3), rv(3), rv(3))));

                expect(ti.test(t0, t1).intersect).toBe(true);
                const poly = fi.find(t0, t1).intersection;
                expect(poly.length).toBeGreaterThanOrEqual(3);
                expect(insideTriangle(q, t0, 1e-9)).toBe(true);
                expect(insideTriangle(q, t1, 1e-9)).toBe(true);
            }
        });

    it('agrees with a Monte-Carlo area estimate of the intersection', () => {
        const fi = new IntrTriangle2Triangle2FI();
        for (let k = 0; k < 40; ++k) {
            const t0 = makeCcw(Triangle.fromVertices(rv(3), rv(3), rv(3)));
            const t1 = makeCcw(Triangle.fromVertices(rv(3), rv(3), rv(3)));
            const poly = fi.find(t0, t1).intersection;
            const exact = poly.length > 0 ? Math.abs(polygonArea(poly)) : 0;

            // Sample the bounding box of t0 and t1 and count the points that
            // are inside both triangles.
            let xmin = Number.MAX_VALUE, xmax = -Number.MAX_VALUE;
            let ymin = Number.MAX_VALUE, ymax = -Number.MAX_VALUE;
            for (const t of [t0, t1]) {
                for (const v of t.v) {
                    xmin = Math.min(xmin, v.values[0]);
                    xmax = Math.max(xmax, v.values[0]);
                    ymin = Math.min(ymin, v.values[1]);
                    ymax = Math.max(ymax, v.values[1]);
                }
            }
            const n = 4000;
            let inBoth = 0;
            for (let s = 0; s < n; ++s) {
                const p = V2(xmin + (xmax - xmin) * rand(),
                    ymin + (ymax - ymin) * rand());
                if (insideTriangle(p, t0, 0) && insideTriangle(p, t1, 0)) {
                    ++inBoth;
                }
            }
            const boxArea = (xmax - xmin) * (ymax - ymin);
            const estimate = (boxArea * inBoth) / n;
            // Monte-Carlo error scales like 1/sqrt(n); allow a generous
            // tolerance relative to the bounding-box area.
            expect(Math.abs(estimate - exact)).toBeLessThan(0.08 * boxArea);
        }
    });
});

// ---------------------------------------------------------------------------
// Verification (V35): the TI query is cross-checked against a brute-force
// method of separating axes over the six edge normals, and the FI polygon
// against an independent Sutherland-Hodgman clip of triangle1 by the three
// halfplanes of triangle0.
// ---------------------------------------------------------------------------

describe('IntrTriangle2Triangle2 verification', () => {
    const tiv = new IntrTriangle2Triangle2TI();
    const fiv = new IntrTriangle2Triangle2FI();

    function perp2(v: Vector): Vector {
        return Vector.fromArray([v.values[1], -v.values[0]]);
    }

    // The six candidate separating directions: the edge normals of both
    // triangles, unnormalized (so integer input stays exact).
    function candidateAxes(t0: Triangle, t1: Triangle): Vector[] {
        const axes: Vector[] = [];
        for (const t of [t0, t1]) {
            for (let i0 = 2, i1 = 0; i1 < 3; i0 = i1++) {
                axes.push(perp2(sub(t.v[i1], t.v[i0])));
            }
        }
        return axes;
    }

    function interval(t: Triangle, n: Vector): [number, number] {
        let lo = dot(n, t.v[0]);
        let hi = lo;
        for (let i = 1; i < 3; ++i) {
            const d = dot(n, t.v[i]);
            if (d < lo) { lo = d; } else if (d > hi) { hi = d; }
        }
        return [lo, hi];
    }

    // The two solid triangles have disjoint interiors exactly when some
    // candidate axis gives closed projection intervals that do not overlap;
    // touching intervals count as separated, matching the query, whose
    // WhichSide requires a strictly positive projection.
    function satSeparated(t0: Triangle, t1: Triangle): boolean {
        for (const n of candidateAxes(t0, t1)) {
            if (dot(n, n) === 0) {
                continue;
            }
            const [lo0, hi0] = interval(t0, n);
            const [lo1, hi1] = interval(t1, n);
            if (hi0 <= lo1 || hi1 <= lo0) {
                return true;
            }
        }
        return false;
    }

    function hasTie(t0: Triangle, t1: Triangle, tol: number): boolean {
        for (const n of candidateAxes(t0, t1)) {
            const nn = dot(n, n);
            if (nn === 0) {
                return true;
            }
            const scale = Math.sqrt(nn);
            const [lo0, hi0] = interval(t0, n);
            const [lo1, hi1] = interval(t1, n);
            if (Math.abs(lo1 - hi0) <= tol * scale
                || Math.abs(lo0 - hi1) <= tol * scale) {
                return true;
            }
        }
        return false;
    }

    function twiceArea(poly: readonly Vector[]): number {
        let a = 0;
        for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
            a += poly[j].values[0] * poly[i].values[1]
                - poly[i].values[0] * poly[j].values[1];
        }
        return a;
    }

    // Clip 'poly' by the halfplane { x : dotPerp(b - a, x - a) >= 0 }.
    function clipHalfplane(poly: readonly Vector[], a: Vector, b: Vector):
        Vector[] {
        const edge = sub(b, a);
        const side = (x: Vector): number => dotPerp(edge, sub(x, a));
        const out: Vector[] = [];
        for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
            const sj = side(poly[j]);
            const si = side(poly[i]);
            if (sj >= 0) {
                out.push(poly[j]);
            }
            if ((sj > 0 && si < 0) || (sj < 0 && si > 0)) {
                const t = sj / (sj - si);
                out.push(add(poly[j], mul(t, sub(poly[i], poly[j]))));
            }
        }
        return out;
    }

    // The convex intersection of two counterclockwise triangles.
    function clipTriangles(t0: Triangle, t1: Triangle): Vector[] {
        let poly: Vector[] = [t1.v[0], t1.v[1], t1.v[2]];
        for (let i0 = 0, i1 = 1; i0 < 3; i0 = i1++) {
            poly = clipHalfplane(poly, t0.v[i0], t0.v[(i0 + 1) % 3]);
            if (poly.length === 0) {
                return poly;
            }
        }
        return poly;
    }

    // Counterclockwise triangle with area bounded away from zero.
    const ccwTriangle = fc.array(wellScaledVector(2, -6, 6),
        { minLength: 3, maxLength: 3 })
        .map(vs => {
            const t = Triangle.fromVertices(vs[0], vs[1], vs[2]);
            if (dotPerp(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0])) < 0) {
                const swap = t.v[1];
                t.v[1] = t.v[2];
                t.v[2] = swap;
            }
            return t;
        })
        .filter(t => dotPerp(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0])) > 1);

    it('TI agrees with a brute-force separating-axis test on lattice input',
        () => {
            // Integer coordinates keep every projection onto an unnormalized
            // edge normal exact, so the reference decides separation exactly
            // and no tolerance is needed.
            const rng = seededRandom(0x7a11e5);
            const draw = (): Vector => Vector.fromArray([
                Math.round(rng() * 10) - 5, Math.round(rng() * 10) - 5]);
            let numIntersect = 0, numSeparate = 0;
            for (let iter = 0; iter < 4000; ++iter) {
                const mk = (): Triangle | null => {
                    const t = Triangle.fromVertices(draw(), draw(), draw());
                    const area = dotPerp(sub(t.v[1], t.v[0]),
                        sub(t.v[2], t.v[0]));
                    if (area === 0) {
                        return null;
                    }
                    if (area < 0) {
                        const swap = t.v[1];
                        t.v[1] = t.v[2];
                        t.v[2] = swap;
                    }
                    return t;
                };
                const t0 = mk();
                const t1 = mk();
                if (t0 === null || t1 === null) {
                    continue;
                }
                const separated = satSeparated(t0, t1);
                expect(tiv.test(t0, t1).intersect).toBe(!separated);
                if (separated) { ++numSeparate; } else { ++numIntersect; }
            }
            expect(numIntersect).toBeGreaterThan(100);
            expect(numSeparate).toBeGreaterThan(100);
        }, 30000);

    it('TI agrees with a brute-force separating-axis test on real input',
        () => {
            check(fc.tuple(ccwTriangle, ccwTriangle), ([t0, t1]) => {
                if (hasTie(t0, t1, 1e-9)) {
                    return;
                }
                expect(tiv.test(t0, t1).intersect).toBe(!satSeparated(t0, t1));
            });
        });

    it('TI is symmetric in its arguments', () => {
        check(fc.tuple(ccwTriangle, ccwTriangle), ([t0, t1]) => {
            if (hasTie(t0, t1, 1e-9)) {
                return;
            }
            expect(tiv.test(t1, t0).intersect).toBe(tiv.test(t0, t1).intersect);
        });
    });

    it('FI equals an independent convex clip', () => {
        check(fc.tuple(ccwTriangle, ccwTriangle), ([t0, t1]) => {
            const got = fiv.find(t0, t1).intersection;
            const want = clipTriangles(t0, t1);
            const areaGot = twiceArea(got);
            const areaWant = twiceArea(want);
            // Both polygons are formed from the same intersections of the
            // same lines, so the areas agree to rounding; the tolerance is
            // relative to the triangle sizes (coordinates below 6).
            expectClose(areaGot, areaWant, 1e-7, 1e-7);
            if (areaWant > 1e-6) {
                expect(got.length).toBeGreaterThanOrEqual(3);
            }
        });
    });

    it('FI returns a counterclockwise convex polygon inside both triangles',
        () => {
            check(fc.tuple(ccwTriangle, ccwTriangle), ([t0, t1]) => {
                const poly = fiv.find(t0, t1).intersection;
                if (poly.length < 3) {
                    return;
                }
                expect(twiceArea(poly)).toBeGreaterThanOrEqual(-1e-9);
                // Convexity: every turn is a left turn (or straight).
                for (let i = 0; i < poly.length; ++i) {
                    const a = poly[i];
                    const b = poly[(i + 1) % poly.length];
                    const c = poly[(i + 2) % poly.length];
                    expect(dotPerp(sub(b, a), sub(c, b)))
                        .toBeGreaterThanOrEqual(-1e-9);
                }
                // Containment in both triangles. The clipped vertices are
                // computed by linear interpolation along the clipping lines,
                // so the residual is proportional to the coordinate size.
                for (const x of poly) {
                    for (const t of [t0, t1]) {
                        for (let i0 = 0; i0 < 3; ++i0) {
                            const a = t.v[i0];
                            const b = t.v[(i0 + 1) % 3];
                            expect(dotPerp(sub(b, a), sub(x, a)))
                                .toBeGreaterThan(-1e-8);
                        }
                    }
                }
            });
        });

    it('a triangle intersected with itself returns its own area', () => {
        check(ccwTriangle, t => {
            const poly = fiv.find(t, t).intersection;
            expect(poly.length).toBeGreaterThanOrEqual(3);
            expectClose(twiceArea(poly),
                dotPerp(sub(t.v[1], t.v[0]), sub(t.v[2], t.v[0])),
                1e-9, 1e-9);
            expect(tiv.test(t, t).intersect).toBe(true);
        });
    });

    it('TI implies a nonempty FI polygon, and a positive-area FI polygon'
        + ' implies TI', () => {
            // TI reports separation for measure-zero contact (its WhichSide
            // requires a strictly positive projection), while FI still
            // returns the touching vertex or edge, so the two do not agree in
            // both directions. These implications do hold.
            check(fc.tuple(ccwTriangle, ccwTriangle), ([t0, t1]) => {
                const poly = fiv.find(t0, t1).intersection;
                if (tiv.test(t0, t1).intersect) {
                    expect(poly.length).toBeGreaterThan(0);
                }
                if (twiceArea(poly) > 1e-7) {
                    expect(tiv.test(t0, t1).intersect).toBe(true);
                }
            });
        });

    it('treats measure-zero contact as separation in TI', () => {
        // Two triangles meeting at the single point (1,0). The TI WhichSide
        // requires a strictly positive projection, so the shared-vertex axis
        // reports separation. The FI clip drops the contact as well here
        // (the third clipping line leaves the polygon on its negative side
        // with one zero-height vertex), but that is a knife edge: for other
        // touching configurations the clip returns a zero-area polygon
        // instead, which is why the TI/FI relation is only asserted in the
        // two safe directions above.
        const a = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([1, 0]), Vector.fromArray([0, 1]));
        const b = Triangle.fromVertices(Vector.fromArray([1, 0]),
            Vector.fromArray([2, 0]), Vector.fromArray([1, 1]));
        expect(tiv.test(a, b).intersect).toBe(false);
        expect(tiv.test(b, a).intersect).toBe(false);
        const poly = fiv.find(a, b).intersection;
        expect(Math.abs(twiceArea(poly))).toBeLessThan(1e-12);
        for (const x of poly) {
            expectVectorClose(x, Vector.fromArray([1, 0]), 1e-9, 1e-9);
        }
        // Two triangles sharing the edge from (0,0) to (0,1) likewise report
        // separation.
        const c = Triangle.fromVertices(Vector.fromArray([0, 0]),
            Vector.fromArray([-1, 0]), Vector.fromArray([0, 1]));
        expect(tiv.test(a, c).intersect).toBe(false);
        expect(Math.abs(twiceArea(fiv.find(a, c).intersection)))
            .toBeLessThan(1e-12);
    });

    it('is equivariant under rigid motions', () => {
        check(fc.tuple(ccwTriangle, ccwTriangle, rotationFrame(2),
            wellScaledVector(2, -5, 5)), ([t0, t1, frame, tr]) => {
                if (hasTie(t0, t1, 1e-7)) {
                    return;
                }
                const xf = (v: Vector): Vector => Vector.fromArray([
                    frame[0].values[0] * v.values[0]
                        + frame[1].values[0] * v.values[1] + tr.values[0],
                    frame[0].values[1] * v.values[0]
                        + frame[1].values[1] * v.values[1] + tr.values[1]]);
                const a = Triangle.fromVertices(xf(t0.v[0]), xf(t0.v[1]),
                    xf(t0.v[2]));
                const b = Triangle.fromVertices(xf(t1.v[0]), xf(t1.v[1]),
                    xf(t1.v[2]));
                expect(tiv.test(a, b).intersect)
                    .toBe(tiv.test(t0, t1).intersect);
                expectClose(twiceArea(fiv.find(a, b).intersection),
                    twiceArea(fiv.find(t0, t1).intersection), 1e-7, 1e-7);
            });
    });
});
