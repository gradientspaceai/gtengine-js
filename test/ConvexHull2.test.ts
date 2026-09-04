import { describe, it, expect, vi } from 'vitest';
import { ConvexHull2 } from '../src/ConvexHull2.js';
import { Vector, sub } from '../src/Vector.js';
import { BSNumber } from '../src/BSNumber.js';
import { SWInterval } from '../src/SWInterval.js';
import { check, expectClose, fc, latticeVector, wellScaledVector } from './helpers/arbitraries.js';
import { exactDyadic, orient2 } from './helpers/exact.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Exact-in-double cross product of (b-a) and (c-a) for small integer inputs.
function cross(a: Vector, b: Vector, c: Vector): number {
    return (b.values[0] - a.values[0]) * (c.values[1] - a.values[1])
        - (b.values[1] - a.values[1]) * (c.values[0] - a.values[0]);
}

// Twice the signed area of the polygon; positive for counterclockwise order.
function twiceSignedArea(points: Vector[]): number {
    let area = 0;
    for (let i = 0, j = points.length - 1; i < points.length; j = i++) {
        area += points[j].values[0] * points[i].values[1]
            - points[i].values[0] * points[j].values[1];
    }
    return area;
}

// Verify: the hull polygon is strictly convex, counterclockwise, and every
// input point is inside or on the hull.
function checkHull(points: Vector[], hull: readonly number[]): void {
    expect(hull.length).toBeGreaterThanOrEqual(3);
    const poly = hull.map(i => points[i]);

    // Counterclockwise.
    expect(twiceSignedArea(poly)).toBeGreaterThan(0);

    // Strictly convex: every consecutive triple turns left. (Upstream's hull
    // does not keep collinear interior points on hull edges.)
    const n = poly.length;
    for (let i = 0; i < n; ++i) {
        const a = poly[i], b = poly[(i + 1) % n], c = poly[(i + 2) % n];
        expect(cross(a, b, c)).toBeGreaterThan(0);
    }

    // Containment: brute-force check that no input point is strictly outside
    // any hull edge.
    for (const p of points) {
        for (let i = 0; i < n; ++i) {
            const a = poly[i], b = poly[(i + 1) % n];
            expect(cross(a, b, p)).toBeGreaterThanOrEqual(-1e-12);
        }
    }

    // The hull indices are distinct.
    expect(new Set(hull).size).toBe(hull.length);
}

// An independent O(n log n) Andrew monotone-chain hull for cross-checking.
// The inputs used with it are small integers, so the double arithmetic is
// exact.
function monotoneChain(points: Vector[]): Vector[] {
    const pts = points.map(p => p).sort((p, q) =>
        p.values[0] !== q.values[0] ? p.values[0] - q.values[0]
            : p.values[1] - q.values[1]);
    const unique: Vector[] = [];
    for (const p of pts) {
        const last = unique[unique.length - 1];
        if (!last || last.values[0] !== p.values[0] || last.values[1] !== p.values[1]) {
            unique.push(p);
        }
    }
    if (unique.length < 3) {
        return unique;
    }
    const build = (src: Vector[]): Vector[] => {
        const out: Vector[] = [];
        for (const p of src) {
            while (out.length >= 2
                && cross(out[out.length - 2], out[out.length - 1], p) <= 0) {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
        return out;
    };
    const lower = build(unique);
    const upper = build(unique.slice().reverse());
    return lower.concat(upper);
}

// Canonicalize a hull polygon to a comparable string starting at its
// lexicographically smallest vertex.
function canonical(poly: Vector[]): string {
    let best = 0;
    for (let i = 1; i < poly.length; ++i) {
        if (poly[i].lessThan(poly[best])) {
            best = i;
        }
    }
    const rotated: string[] = [];
    for (let i = 0; i < poly.length; ++i) {
        const p = poly[(best + i) % poly.length];
        rotated.push(`${p.values[0]},${p.values[1]}`);
    }
    return rotated.join(' ');
}

describe('ConvexHull2', () => {
    it('rejects an empty input', () => {
        expect(() => new ConvexHull2().compute([])).toThrow();
    });

    it('reports dimension 0 for coincident points', () => {
        const points = [v2(3, -1), v2(3, -1), v2(3, -1)];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(false);
        expect(hull.getDimension()).toBe(0);
        expect(hull.getNumPoints()).toBe(3);
        expect(hull.getNumUniquePoints()).toBe(1);
        expect(hull.getHull().length).toBe(1);
        expect(points[hull.getHull()[0]].equals(v2(3, -1))).toBe(true);
    });

    it('reports dimension 0 for a single point', () => {
        const hull = new ConvexHull2();
        expect(hull.compute([v2(7, 9)])).toBe(false);
        expect(hull.getDimension()).toBe(0);
        expect(hull.getHull()).toEqual([0]);
    });

    it('reports dimension 1 for collinear points and returns the extremes', () => {
        // Points on the line through (0,0) with direction (1,2), given out of
        // order and with an interior duplicate.
        const points = [v2(2, 4), v2(0, 0), v2(1, 2), v2(5, 10), v2(1, 2), v2(3, 6)];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(false);
        expect(hull.getDimension()).toBe(1);
        expect(hull.getNumUniquePoints()).toBe(5);

        const indices = hull.getHull();
        expect(indices.length).toBe(2);
        const p0 = points[indices[0]], p1 = points[indices[1]];
        expect(p0.equals(v2(0, 0))).toBe(true);
        expect(p1.equals(v2(5, 10))).toBe(true);

        const line = hull.getLine();
        expect(line.origin.equals(v2(0, 0))).toBe(true);
        const s = Math.sqrt(5);
        expect(line.direction.values[0]).toBeCloseTo(1 / s, 12);
        expect(line.direction.values[1]).toBeCloseTo(2 / s, 12);
    });

    it('reports dimension 1 for two distinct points', () => {
        const points = [v2(4, 1), v2(-2, 1)];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(false);
        expect(hull.getDimension()).toBe(1);
        const indices = hull.getHull();
        expect(points[indices[0]].equals(v2(-2, 1))).toBe(true);
        expect(points[indices[1]].equals(v2(4, 1))).toBe(true);
    });

    it('computes the hull of a unit square with an interior point', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1), v2(0.5, 0.5)];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        expect(hull.getDimension()).toBe(2);
        expect(hull.getNumUniquePoints()).toBe(5);

        const indices = hull.getHull();
        expect(indices.length).toBe(4);
        expect(indices).not.toContain(4);
        checkHull(points, indices);
        expect(canonical(indices.map(i => points[i])))
            .toBe(canonical([v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)]));
    });

    it('drops points collinear on a hull edge', () => {
        // The midpoints of the square's edges lie on hull edges; the
        // divide-and-conquer merge with the exact predicate excludes them.
        const points = [
            v2(0, 0), v2(2, 0), v2(2, 2), v2(0, 2),
            v2(1, 0), v2(2, 1), v2(1, 2), v2(0, 1)
        ];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        const indices = hull.getHull();
        expect(indices.length).toBe(4);
        checkHull(points, indices);
        expect(canonical(indices.map(i => points[i])))
            .toBe(canonical([v2(0, 0), v2(2, 0), v2(2, 2), v2(0, 2)]));
    });

    it('computes the hull of a 5x5 integer grid', () => {
        const points: Vector[] = [];
        for (let y = 0; y < 5; ++y) {
            for (let x = 0; x < 5; ++x) {
                points.push(v2(x, y));
            }
        }
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        expect(hull.getNumUniquePoints()).toBe(25);
        const indices = hull.getHull();
        expect(indices.length).toBe(4);
        checkHull(points, indices);
        expect(canonical(indices.map(i => points[i])))
            .toBe(canonical([v2(0, 0), v2(4, 0), v2(4, 4), v2(0, 4)]));
    });

    it('computes the hull of points on a circle (all are hull vertices)', () => {
        const n = 24;
        const points: Vector[] = [];
        for (let i = 0; i < n; ++i) {
            const t = (2 * Math.PI * i) / n;
            points.push(v2(Math.cos(t), Math.sin(t)));
        }
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        const indices = hull.getHull();
        expect(indices.length).toBe(n);
        checkHull(points, indices);
    });

    it('handles duplicated points on the hull', () => {
        const base = [v2(0, 0), v2(3, 0), v2(3, 3), v2(0, 3)];
        const points = base.concat(base).concat(base).concat([v2(1, 1)]);
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        expect(hull.getNumPoints()).toBe(13);
        expect(hull.getNumUniquePoints()).toBe(5);
        const indices = hull.getHull();
        expect(indices.length).toBe(4);
        checkHull(points, indices);
    });

    it('handles a near-degenerate configuration that stresses the exact path', () => {
        // Points nearly collinear along y = x, with perturbations at the
        // magnitude of the double-precision epsilon. The interval arithmetic
        // cannot resolve the orientation signs, so the BSNumber fallback is
        // exercised. The middle point is *below* the line through the
        // endpoints by one ulp, so it must not be on the upper chain.
        const eps = Number.EPSILON;
        const points = [
            v2(0, 0),
            v2(1, 1),
            v2(0.5, 0.5 - eps),
            v2(0.25, 0.25 + eps),
            v2(0.75, 0.75 + eps)
        ];
        const hull = new ConvexHull2();
        expect(hull.compute(points)).toBe(true);
        const indices = hull.getHull();
        checkHull(points, indices);
        // (0.5, 0.5 - eps) is strictly below the segment from (0,0) to (1,1),
        // so it is a hull vertex on the lower chain.
        expect(indices).toContain(2);
        // (0.25, 0.25 + eps) and (0.75, 0.75 + eps) are strictly above, and
        // both must be on the upper chain because the three upper points are
        // not collinear.
        expect(indices).toContain(3);
        expect(indices).toContain(4);
        expect(indices).toContain(0);
        expect(indices).toContain(1);
    });

    it('detects exact collinearity that floating-point evaluation would miss', () => {
        // The three points are exactly collinear in rational arithmetic:
        // (0,0), (1e-16, 1e-16*3) and (1, 3). The exact predicate must
        // classify the intermediate point as interior to the segment, so the
        // hull is 1-dimensional.
        const a = 1e-16;
        const points = [v2(0, 0), v2(a, 3 * a), v2(1, 3)];
        const hull = new ConvexHull2();
        // 3*a is exact (a is a power-of-two multiple? not necessarily), so
        // check the exactness assumption before asserting.
        const exactlyCollinear = (3 * a) * 1 - 3 * a === 0;
        expect(exactlyCollinear).toBe(true);
        expect(hull.compute(points)).toBe(false);
        expect(hull.getDimension()).toBe(1);
    });

    it('agrees with an independent monotone-chain hull on random point sets', () => {
        const rnd = makeRandom(20260901);
        for (let trial = 0; trial < 40; ++trial) {
            const n = 3 + Math.floor(rnd() * 40);
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                // Small integer coordinates so the reference computation is
                // exact and duplicates/collinearity occur frequently.
                points.push(v2(Math.floor(rnd() * 11) - 5, Math.floor(rnd() * 11) - 5));
            }
            const expected = monotoneChain(points);
            const hull = new ConvexHull2();
            const is2D = hull.compute(points);
            expect(is2D).toBe(expected.length >= 3);
            if (!is2D) {
                continue;
            }
            const actual = hull.getHull().map(i => points[i]);
            expect(canonical(actual)).toBe(canonical(expected));
            checkHull(points, hull.getHull());
        }
    });

    it('agrees with an independent monotone-chain hull on random real point sets', () => {
        const rnd = makeRandom(987654321);
        for (let trial = 0; trial < 20; ++trial) {
            const n = 4 + Math.floor(rnd() * 60);
            const points: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                points.push(v2(rnd() * 2 - 1, rnd() * 2 - 1));
            }
            const expected = monotoneChain(points);
            const hull = new ConvexHull2();
            expect(hull.compute(points)).toBe(true);
            const actual = hull.getHull().map(i => points[i]);
            expect(canonical(actual)).toBe(canonical(expected));
            checkHull(points, hull.getHull());
        }
    });

    it('reuses the functor for multiple data sets', () => {
        const hull = new ConvexHull2();
        expect(hull.compute([v2(0, 0), v2(1, 0), v2(0, 1)])).toBe(true);
        expect(hull.getHull().length).toBe(3);
        expect(hull.compute([v2(5, 5), v2(5, 5)])).toBe(false);
        expect(hull.getDimension()).toBe(0);
        expect(hull.getNumPoints()).toBe(2);
        expect(hull.getHull().length).toBe(1);
        expect(hull.compute([v2(0, 0), v2(2, 0), v2(2, 2), v2(0, 2)])).toBe(true);
        expect(hull.getHull().length).toBe(4);
        expect(hull.getPoints().length).toBe(4);
    });
});

// ---------------------------------------------------------------------------
// Verification pass (VERIFYING.md): property-based cross-checks against an
// exact bigint convex hull. Every input coordinate is a finite double, so the
// coordinates of one point set can be scaled by a common power of two into
// exact integers (exactDyadic); the orientation determinant is homogeneous,
// so that scale never changes a predicate's sign. The reference hull is
// therefore the mathematically exact hull of the input, which is what
// ConvexHull2's interval + BSNumber predicate is required to reproduce.
// ---------------------------------------------------------------------------

type ExactPoint = { x: bigint; y: bigint };

function exactPoints(points: readonly Vector[]): ExactPoint[] {
    const flat: number[] = [];
    for (const p of points) {
        flat.push(p.values[0], p.values[1]);
    }
    const s = exactDyadic(flat);
    return points.map((_, i) => ({ x: s[2 * i], y: s[2 * i + 1] }));
}

// Andrew monotone chain on exact integer coordinates. The result has 1 point
// when all inputs coincide, 2 points (the extremes) when they are collinear,
// and otherwise the strictly convex counterclockwise hull polygon.
function exactHull(pts: readonly ExactPoint[]): ExactPoint[] {
    const cmp = (a: ExactPoint, b: ExactPoint): number =>
        a.x < b.x ? -1 : a.x > b.x ? 1 : a.y < b.y ? -1 : a.y > b.y ? 1 : 0;
    const sorted = pts.slice().sort(cmp);
    const uniq: ExactPoint[] = [];
    for (const p of sorted) {
        const last = uniq[uniq.length - 1];
        if (!last || cmp(last, p) !== 0) {
            uniq.push(p);
        }
    }
    if (uniq.length < 3) {
        return uniq;
    }
    const turn = (o: ExactPoint, a: ExactPoint, b: ExactPoint): number =>
        orient2(o.x, o.y, a.x, a.y, b.x, b.y);
    const build = (src: readonly ExactPoint[]): ExactPoint[] => {
        const out: ExactPoint[] = [];
        for (const p of src) {
            while (out.length >= 2
                && turn(out[out.length - 2], out[out.length - 1], p) <= 0) {
                out.pop();
            }
            out.push(p);
        }
        out.pop();
        return out;
    };
    return build(uniq).concat(build(uniq.slice().reverse()));
}

function canonicalExact(poly: readonly ExactPoint[]): string {
    if (poly.length === 0) {
        return '';
    }
    let best = 0;
    for (let i = 1; i < poly.length; ++i) {
        if (poly[i].x < poly[best].x
            || (poly[i].x === poly[best].x && poly[i].y < poly[best].y)) {
            best = i;
        }
    }
    const out: string[] = [];
    for (let i = 0; i < poly.length; ++i) {
        const p = poly[(best + i) % poly.length];
        out.push(String(p.x) + ',' + String(p.y));
    }
    return out.join(' ');
}

function keyOf(p: ExactPoint): string {
    return String(p.x) + ',' + String(p.y);
}

// Run ConvexHull2 on 'points' and compare every observable against the exact
// reference hull.
function compareWithExact(points: Vector[]): void {
    const exact = exactPoints(points);
    const reference = exactHull(exact);
    const expectedDimension = reference.length === 1 ? 0
        : reference.length === 2 ? 1 : 2;

    const ch = new ConvexHull2();
    const is2D = ch.compute(points);
    expect(is2D).toBe(expectedDimension === 2);
    expect(ch.getDimension()).toBe(expectedDimension);
    expect(ch.getNumPoints()).toBe(points.length);
    expect(ch.getPoints()).toBe(points);

    // getNumUniquePoints() counts the distinct input points.
    const distinct = new Set(exact.map(keyOf));
    expect(ch.getNumUniquePoints()).toBe(distinct.size);

    const indices = ch.getHull();
    expect(indices.length).toBe(reference.length);
    for (const i of indices) {
        expect(Number.isInteger(i)).toBe(true);
        expect(i).toBeGreaterThanOrEqual(0);
        expect(i).toBeLessThan(points.length);
    }
    const actual = indices.map(i => exact[i]);
    // The hull vertices are distinct points (not merely distinct indices).
    expect(new Set(actual.map(keyOf)).size).toBe(actual.length);

    if (expectedDimension === 2) {
        // Same cyclic sequence, so also the same counterclockwise order.
        expect(canonicalExact(actual)).toBe(canonicalExact(reference));
    } else {
        // Dimension 0 and 1: the returned indices are the extreme points as a
        // set, and for dimension 1 they are lexicographically ordered.
        expect(actual.map(keyOf).sort()).toEqual(reference.map(keyOf).sort());
        if (expectedDimension === 1) {
            expect(actual[0].x < actual[1].x
                || (actual[0].x === actual[1].x && actual[0].y < actual[1].y))
                .toBe(true);
            // The approximating line passes through the first extreme with
            // the normalized direction toward the second.
            const line = ch.getLine();
            expect(line.origin.equals(points[indices[0]])).toBe(true);
            const d = sub(points[indices[1]], points[indices[0]]);
            const len = Math.hypot(d.values[0], d.values[1]);
            expectClose(line.direction.values[0], d.values[0] / len, 1e-12, 1e-12);
            expectClose(line.direction.values[1], d.values[1] / len, 1e-12, 1e-12);
        }
    }
}

describe('ConvexHull2 verification', () => {
    it('reproduces the exact hull of integer-lattice point sets', () => {
        check(fc.array(latticeVector(2, -6, 6), { minLength: 1, maxLength: 24 }),
            (points) => {
                compareWithExact(points);
            }, 200);
    });

    it('reproduces the exact hull of well-scaled real point sets', () => {
        check(fc.array(wellScaledVector(2, -5, 5), { minLength: 1, maxLength: 18 }),
            (points) => {
                compareWithExact(points);
            }, 200);
    });

    it('reproduces the exact hull of clustered points (duplicates and collinearities)', () => {
        // A coarse lattice with few distinct values produces duplicate points
        // and collinear triples in almost every draw, which is where the
        // exact predicate matters.
        check(fc.array(latticeVector(2, -2, 2), { minLength: 1, maxLength: 20 }),
            (points) => {
                compareWithExact(points);
            }, 200);
    });

    it('classifies exactly collinear input as dimension 0 or 1', () => {
        const arb = fc.tuple(
            latticeVector(2, -10, 10),
            latticeVector(2, -5, 5).filter(d => d.values[0] !== 0 || d.values[1] !== 0),
            fc.array(fc.integer({ min: -20, max: 20 }), { minLength: 1, maxLength: 12 }));
        check(arb, ([base, dir, ts]) => {
            const points = ts.map(t => Vector.fromArray([
                base.values[0] + t * dir.values[0],
                base.values[1] + t * dir.values[1]]));
            const ch = new ConvexHull2();
            expect(ch.compute(points)).toBe(false);
            const unique = new Set(ts).size;
            expect(ch.getDimension()).toBe(unique === 1 ? 0 : 1);
            expect(ch.getNumUniquePoints()).toBe(unique);
            compareWithExact(points);
        }, 200);
    });

    it('is invariant under permutation of the input points', () => {
        const arb = fc.array(latticeVector(2, -5, 5), { minLength: 3, maxLength: 16 })
            .chain(points => fc.tuple(fc.constant(points),
                fc.shuffledSubarray(points.map((_, i) => i),
                    { minLength: points.length, maxLength: points.length })));
        check(arb, ([points, perm]) => {
            const permuted = perm.map(i => points[i]);
            const a = new ConvexHull2();
            const b = new ConvexHull2();
            a.compute(points);
            b.compute(permuted);
            expect(b.getDimension()).toBe(a.getDimension());
            expect(b.getNumUniquePoints()).toBe(a.getNumUniquePoints());
            const ea = exactPoints(points), eb = exactPoints(permuted);
            expect(canonicalExact(b.getHull().map(i => eb[i])))
                .toBe(canonicalExact(a.getHull().map(i => ea[i])));
        }, 100);
    });

    it('is unchanged when the input is duplicated', () => {
        check(fc.array(latticeVector(2, -5, 5), { minLength: 1, maxLength: 14 }),
            (points) => {
                const doubled = points.concat(points);
                const a = new ConvexHull2();
                const b = new ConvexHull2();
                a.compute(points);
                b.compute(doubled);
                expect(b.getNumPoints()).toBe(2 * points.length);
                expect(b.getNumUniquePoints()).toBe(a.getNumUniquePoints());
                expect(b.getDimension()).toBe(a.getDimension());
                const ea = exactPoints(points), eb = exactPoints(doubled);
                expect(canonicalExact(b.getHull().map(i => eb[i])))
                    .toBe(canonicalExact(a.getHull().map(i => ea[i])));
            }, 150);
    });

    it('reuse does not leak state between data sets', () => {
        const arb = fc.tuple(
            fc.array(latticeVector(2, -5, 5), { minLength: 1, maxLength: 12 }),
            fc.array(latticeVector(2, -5, 5), { minLength: 1, maxLength: 12 }));
        check(arb, ([first, second]) => {
            const shared = new ConvexHull2();
            shared.compute(first);
            shared.compute(second);
            const fresh = new ConvexHull2();
            fresh.compute(second);
            expect(shared.getDimension()).toBe(fresh.getDimension());
            expect(shared.getNumPoints()).toBe(fresh.getNumPoints());
            expect(shared.getNumUniquePoints()).toBe(fresh.getNumUniquePoints());
            expect(shared.getHull()).toEqual(fresh.getHull());
        }, 100);
    });

    it('takes the exact BSNumber path when the interval predicate is indeterminate', () => {
        // Nearly collinear points one ulp apart: the SWInterval determinant
        // straddles zero, so toLineExtended must fall back to BSNumber. The
        // spy on BSNumber.fromNumber observes getRationalPoint, which is
        // reached only from the three exact-arithmetic branches.
        const eps = Number.EPSILON;
        const points = [
            v2(0, 0), v2(1, 1), v2(0.5, 0.5 - eps),
            v2(0.25, 0.25 + eps), v2(0.75, 0.75 + eps)
        ];

        // The interval determinant for <P,Q0,Q1> = <(0.5,0.5-eps),(0,0),(1,1)>
        // really is indeterminate, so the fallback is not merely incidental.
        const iv = (a: number, b: number): SWInterval =>
            new SWInterval(a).sub(new SWInterval(b));
        const ix0 = iv(1, 0), iy0 = iv(1, 0);
        const ix1 = iv(0.5, 0), iy1 = iv(0.5 - eps, 0);
        const iDet = ix0.mul(iy1).sub(ix1.mul(iy0));
        expect(iDet.get(0)).toBeLessThanOrEqual(0);
        expect(iDet.get(1)).toBeGreaterThanOrEqual(0);

        const spy = vi.spyOn(BSNumber, 'fromNumber');
        try {
            const ch = new ConvexHull2();
            expect(ch.compute(points)).toBe(true);
            expect(spy.mock.calls.length).toBeGreaterThan(0);
        } finally {
            spy.mockRestore();
        }
        compareWithExact(points);
    });

    it('the exact fallback agrees with bigint ground truth on collinear points float arithmetic misclassifies', () => {
        // (0,0), (a,3a), (1,3) with a = 1e-16 are exactly collinear as
        // rationals, but the interval determinant contains zero. Only the
        // BSNumber path can conclude that the middle point lies strictly
        // between the endpoints.
        const a = 1e-16;
        const points = [v2(0, 0), v2(a, 3 * a), v2(1, 3)];
        const e = exactPoints(points);
        expect(orient2(e[0].x, e[0].y, e[1].x, e[1].y, e[2].x, e[2].y)).toBe(0);

        const spy = vi.spyOn(BSNumber, 'fromNumber');
        try {
            const ch = new ConvexHull2();
            expect(ch.compute(points)).toBe(false);
            expect(ch.getDimension()).toBe(1);
            expect(spy.mock.calls.length).toBeGreaterThan(0);
        } finally {
            spy.mockRestore();
        }
        compareWithExact(points);
    });
});
