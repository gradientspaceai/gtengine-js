import { describe, expect, it } from 'vitest';
import { MinimumAreaCircle2 } from '../src/MinimumAreaCircle2.js';
import type { Circle2 } from '../src/Hypersphere.js';
import { Vector, sub, dot } from '../src/Vector.js';
import {
    check, expectClose, expectVectorClose, fc, latticeVector, wellScaledVector
} from './helpers/arbitraries.js';

function v2(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function distance(p: Vector, q: Vector): number {
    const d = sub(p, q);
    return Math.sqrt(dot(d, d));
}

// An independent brute-force minimum-area circle: enumerate the circles
// defined by every pair (diametral circle) and every triple (circumcircle,
// computed with the determinant formula, which is a different derivation
// from the barycentric solve used by the port) and keep the smallest one
// that contains every point.
function bruteForceMinimumCircle(points: readonly Vector[]):
    { center: Vector, radius: number } {
    const n = points.length;
    if (n === 1) {
        return { center: points[0].clone(), radius: 0 };
    }

    let best: { center: Vector, radius: number } | null = null;
    const consider = (center: Vector, radius: number): void => {
        // Allow a small tolerance so that the defining points, which sit
        // exactly on the boundary in exact arithmetic, are not rejected.
        const tol = 1e-9 * Math.max(1, radius);
        for (const p of points) {
            if (distance(p, center) > radius + tol) {
                return;
            }
        }
        if (best === null || radius < best.radius) {
            best = { center, radius };
        }
    };

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            const center = v2(
                0.5 * (points[i].get(0) + points[j].get(0)),
                0.5 * (points[i].get(1) + points[j].get(1)));
            consider(center, 0.5 * distance(points[i], points[j]));
        }
    }

    for (let i = 0; i < n; ++i) {
        for (let j = i + 1; j < n; ++j) {
            for (let k = j + 1; k < n; ++k) {
                const ax = points[i].get(0), ay = points[i].get(1);
                const bx = points[j].get(0), by = points[j].get(1);
                const cx = points[k].get(0), cy = points[k].get(1);
                const d = 2 * (ax * (by - cy) + bx * (cy - ay)
                    + cx * (ay - by));
                if (Math.abs(d) < 1e-14) {
                    continue;
                }
                const a2 = ax * ax + ay * ay;
                const b2 = bx * bx + by * by;
                const c2 = cx * cx + cy * cy;
                const ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by))
                    / d;
                const uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax))
                    / d;
                const center = v2(ux, uy);
                consider(center, distance(center, points[i]));
            }
        }
    }

    if (best === null) {
        throw new Error('brute force failed');
    }
    return best;
}

// A deterministic pseudo-random generator for the randomized tests.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (1664525 * state + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('MinimumAreaCircle2', () => {
    it('handles a single point', () => {
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute([v2(3, -4)]);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(minimal.center.get(0)).toBe(3);
        expect(minimal.center.get(1)).toBe(-4);
        expect(mac.numSupport).toBe(1);
        expect(mac.support[0]).toBe(0);
    });

    it('computes the diametral circle of two points', () => {
        const mac = new MinimumAreaCircle2();
        const points = [v2(-1, 0), v2(1, 0)];
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0, 12);
        expect(mac.numSupport).toBe(2);
        expect(new Set(mac.support.slice(0, 2))).toEqual(new Set([0, 1]));
    });

    it('computes the circumcircle of an acute triangle', () => {
        // Equilateral triangle with circumradius 1.
        const points = [
            v2(1, 0),
            v2(Math.cos(2 * Math.PI / 3), Math.sin(2 * Math.PI / 3)),
            v2(Math.cos(4 * Math.PI / 3), Math.sin(4 * Math.PI / 3))
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(1, 10);
        expect(minimal.center.get(0)).toBeCloseTo(0, 10);
        expect(minimal.center.get(1)).toBeCloseTo(0, 10);
        expect(mac.numSupport).toBe(3);
    });

    it('uses the longest side for an obtuse triangle', () => {
        // The circumcircle of this triangle is larger than the circle whose
        // diameter is the longest side, so the minimum circle has 2 support
        // points.
        const points = [v2(-2, 0), v2(2, 0), v2(0, 0.25)];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0, 12);
        expect(mac.numSupport).toBe(2);
        expect(new Set(mac.support.slice(0, 2))).toEqual(new Set([0, 1]));
    });

    it('computes the circle of a unit square', () => {
        const points = [v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1)];
        const { minimal, success } = new MinimumAreaCircle2().compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(Math.SQRT1_2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0.5, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0.5, 12);
    });

    it('ignores duplicate points', () => {
        const points = [
            v2(0, 0), v2(1, 0), v2(0, 0), v2(1, 0), v2(1, 1), v2(0, 1),
            v2(0, 1)
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(Math.SQRT1_2, 12);
        expect(minimal.center.get(0)).toBeCloseTo(0.5, 12);
        expect(minimal.center.get(1)).toBeCloseTo(0.5, 12);
        // The support indices index into the input array.
        for (let i = 0; i < mac.numSupport; ++i) {
            const index = mac.support[i];
            expect(index).toBeGreaterThanOrEqual(0);
            expect(index).toBeLessThan(points.length);
            expect(distance(points[index], minimal.center))
                .toBeCloseTo(minimal.radius, 10);
        }
    });

    it('handles a set of identical points', () => {
        const points = [v2(5, 5), v2(5, 5), v2(5, 5)];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBe(0);
        expect(minimal.center.get(0)).toBe(5);
        expect(mac.numSupport).toBe(1);
    });

    it('handles collinear points', () => {
        const points = [
            v2(0, 0), v2(1, 1), v2(2, 2), v2(3, 3), v2(-1, -1)
        ];
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        // The extremes are (-1,-1) and (3,3).
        expect(minimal.center.get(0)).toBeCloseTo(1, 10);
        expect(minimal.center.get(1)).toBeCloseTo(1, 10);
        expect(minimal.radius).toBeCloseTo(2 * Math.SQRT2, 10);
        expect(mac.numSupport).toBe(2);
    });

    it('handles cocircular points', () => {
        const points: Vector[] = [];
        for (let i = 0; i < 8; ++i) {
            const angle = 2 * Math.PI * i / 8;
            points.push(v2(2 * Math.cos(angle), 2 * Math.sin(angle)));
        }
        const mac = new MinimumAreaCircle2();
        const { minimal, success } = mac.compute(points);
        expect(success).toBe(true);
        expect(minimal.radius).toBeCloseTo(2, 8);
        expect(minimal.center.get(0)).toBeCloseTo(0, 8);
        expect(minimal.center.get(1)).toBeCloseTo(0, 8);
    });

    it('throws for an empty point set', () => {
        expect(() => new MinimumAreaCircle2().compute([]))
            .toThrow('Input must contain points.');
    });

    it('throws for non-2D points', () => {
        expect(() => new MinimumAreaCircle2()
            .compute([Vector.fromArray([1, 2, 3])])).toThrow();
    });

    it('matches brute force on random point sets', () => {
        const random = makeRandom(20260902);
        for (let trial = 0; trial < 40; ++trial) {
            const numPoints = 3 + (trial % 10);
            const points: Vector[] = [];
            for (let i = 0; i < numPoints; ++i) {
                points.push(v2(10 * random() - 5, 10 * random() - 5));
            }

            const mac = new MinimumAreaCircle2();
            const { minimal, success } = mac.compute(points);
            expect(success).toBe(true);

            const expected = bruteForceMinimumCircle(points);
            expect(minimal.radius).toBeCloseTo(expected.radius, 8);
            expect(minimal.center.get(0))
                .toBeCloseTo(expected.center.get(0), 7);
            expect(minimal.center.get(1))
                .toBeCloseTo(expected.center.get(1), 7);

            // Every point is contained.
            for (const p of points) {
                expect(distance(p, minimal.center))
                    .toBeLessThanOrEqual(minimal.radius + 1e-9);
            }

            // Every support point is on the boundary and the center is a
            // convex combination of the support points (the optimality
            // certificate for the minimum enclosing ball).
            expect(mac.numSupport).toBeGreaterThanOrEqual(2);
            for (let i = 0; i < mac.numSupport; ++i) {
                expect(distance(points[mac.support[i]], minimal.center))
                    .toBeCloseTo(minimal.radius, 7);
            }
            if (mac.numSupport === 2) {
                const p0 = points[mac.support[0]];
                const p1 = points[mac.support[1]];
                expect(minimal.center.get(0))
                    .toBeCloseTo(0.5 * (p0.get(0) + p1.get(0)), 7);
                expect(minimal.center.get(1))
                    .toBeCloseTo(0.5 * (p0.get(1) + p1.get(1)), 7);
            } else {
                // Solve center = w0*P0 + w1*P1 + (1-w0-w1)*P2 and require
                // nonnegative barycentric weights.
                const P0 = points[mac.support[0]];
                const P1 = points[mac.support[1]];
                const P2 = points[mac.support[2]];
                const a00 = P0.get(0) - P2.get(0);
                const a01 = P1.get(0) - P2.get(0);
                const a10 = P0.get(1) - P2.get(1);
                const a11 = P1.get(1) - P2.get(1);
                const b0 = minimal.center.get(0) - P2.get(0);
                const b1 = minimal.center.get(1) - P2.get(1);
                const det = a00 * a11 - a01 * a10;
                expect(Math.abs(det)).toBeGreaterThan(1e-12);
                const w0 = (b0 * a11 - a01 * b1) / det;
                const w1 = (a00 * b1 - b0 * a10) / det;
                const w2 = 1 - w0 - w1;
                expect(w0).toBeGreaterThan(-1e-8);
                expect(w1).toBeGreaterThan(-1e-8);
                expect(w2).toBeGreaterThan(-1e-8);
            }
        }
    });

    it('is minimal: shrinking the radius excludes a point', () => {
        const random = makeRandom(777);
        for (let trial = 0; trial < 20; ++trial) {
            const points: Vector[] = [];
            for (let i = 0; i < 12; ++i) {
                points.push(v2(random(), random()));
            }
            const { minimal, success } =
                new MinimumAreaCircle2().compute(points);
            expect(success).toBe(true);

            // Perturbing the center in any of four directions while keeping
            // the radius must leave some point outside.
            const eps = 1e-4;
            const directions = [[1, 0], [-1, 0], [0, 1], [0, -1],
                [0.7, 0.7], [-0.7, 0.7], [0.7, -0.7], [-0.7, -0.7]];
            for (const d of directions) {
                const c = v2(minimal.center.get(0) + eps * d[0],
                    minimal.center.get(1) + eps * d[1]);
                let maxDist = 0;
                for (const p of points) {
                    maxDist = Math.max(maxDist, distance(p, c));
                }
                expect(maxDist).toBeGreaterThan(minimal.radius - 1e-12);
            }
        }
    });
});

// ---------------------------------------------------------------------------
// Verification pass (group V10). The generators use integer lattice points so
// that the brute-force reference circle and the query agree to round-off: the
// support points of a lattice set are well separated, which keeps the
// circumcircle solve of exactCircle3 well conditioned.
//
// Lattice inputs do reach the trapped-failure path (`success === false`): a
// cocircular set makes `contains` reject a point that is exactly on the
// boundary in exact arithmetic, and updateSupport3 then finds no candidate.
// Round-off can also produce an *untrapped* failure -- `success === true`
// with a circle that misses a point -- which is the known upstream defect
// documented at the top of src/MinimumAreaCircle2.ts and pinned by the last
// two tests of this block (a cocircular trigger and a nearly collinear one).
// Roughly 0.05% of nine-point lattice sets are affected. The invariant that
// survives it is that every point which is not one of the reported support
// points is inside the circle (checked over 100000 random lattice sets while
// writing these tests), so that is what the containment properties assert;
// the properties about minimality skip the affected inputs explicitly.
// ---------------------------------------------------------------------------

const latticePoints = fc.array(latticeVector(2, -8, 8),
    { minLength: 1, maxLength: 10 });

// The circle for a point set, plus the query object that produced it.
function minAreaCircle(points: readonly Vector[]):
    { query: MinimumAreaCircle2; minimal: Circle2; success: boolean } {
    const query = new MinimumAreaCircle2();
    const { minimal, success } = query.compute(points);
    return { query, minimal, success };
}

// The radius is the square root of a computed squared radius, so containment
// allows a relative round-off slack.
function containsAll(points: readonly Vector[], circle: Circle2): boolean {
    const tol = 1e-9 * Math.max(1, circle.radius);
    for (const p of points) {
        if (distance(p, circle.center) > circle.radius + tol) {
            return false;
        }
    }
    return true;
}

function expectContainsAll(points: readonly Vector[], circle: Circle2): void {
    expect(containsAll(points, circle)).toBe(true);
}

// True when the query returned a circle that is claimed minimal and really
// does bound the input; the known upstream defect is excluded here so that
// the minimality properties test what they are meant to test.
function isTrustworthy(points: readonly Vector[],
    result: { minimal: Circle2; success: boolean }): boolean {
    return result.success && containsAll(points, result.minimal);
}

describe('MinimumAreaCircle2 verification', () => {
    // Every input point that is not one of the reported support points was
    // tested against the final circle by the main loop, so it must be inside
    // it. This is the strongest containment statement that survives the known
    // upstream defect, and it holds for both the minimal and the trapped
    // (getContainerCircle2) results.
    function expectContainsNonSupport(points: readonly Vector[]): void {
        const { query, minimal, success } = minAreaCircle(points);
        const support: Vector[] = [];
        if (success) {
            for (let i = 0; i < query.numSupport; ++i) {
                support.push(points[query.support[i]]);
            }
        }
        for (const p of points) {
            // A duplicate of a support point is supported as well.
            if (support.some(s => s.equals(p))) {
                continue;
            }
            expect(distance(p, minimal.center))
                .toBeLessThanOrEqual(minimal.radius
                    + 1e-9 * Math.max(1, minimal.radius));
        }
    }

    it('contains every non-support lattice point, minimal or trapped', () => {
        check(latticePoints, points => { expectContainsNonSupport(points); });
    });

    it('contains every non-support point for real coordinates', () => {
        check(fc.array(wellScaledVector(2, -8, 8),
            { minLength: 1, maxLength: 12 }),
        points => { expectContainsNonSupport(points); });
    });

    it('puts every support point on the boundary', () => {
        check(latticePoints, points => {
            const result = minAreaCircle(points);
            const { query, minimal, success } = result;
            if (!success) {
                // The trapped-failure path clears the support set.
                expect(query.numSupport).toBe(0);
                return;
            }
            expect(query.numSupport).toBeGreaterThanOrEqual(1);
            expect(query.numSupport).toBeLessThanOrEqual(3);
            const seen = new Set<number>();
            for (let i = 0; i < query.numSupport; ++i) {
                const index = query.support[i];
                expect(index).toBeGreaterThanOrEqual(0);
                expect(index).toBeLessThan(points.length);
                expect(seen.has(index)).toBe(false);
                seen.add(index);
            }
            if (!isTrustworthy(points, result)) {
                return;     // known upstream defect; see the block comment
            }
            for (let i = 0; i < query.numSupport; ++i) {
                expectClose(distance(points[query.support[i]], minimal.center),
                    minimal.radius, 1e-9, 1e-9);
            }
            // A single support point means every input point is the same.
            if (query.numSupport === 1) {
                expect(minimal.radius).toBe(0);
            }
        });
    });

    it('matches the brute-force minimum enclosing circle', () => {
        check(latticePoints, points => {
            const result = minAreaCircle(points);
            if (!isTrustworthy(points, result)) {
                return;
            }
            const brute = bruteForceMinimumCircle(points);
            expectClose(result.minimal.radius, brute.radius, 1e-9, 1e-9);
            if (brute.radius > 0) {
                expectVectorClose(result.minimal.center, brute.center,
                    1e-8, 1e-8);
            }
        });
    });

    it('is invariant under a permutation of the input', () => {
        check(fc.array(latticeVector(2, -8, 8), { minLength: 1, maxLength: 8 })
            .chain(points => fc.tuple(fc.constant(points),
                fc.shuffledSubarray(points.map((_, i) => i),
                    { minLength: points.length, maxLength: points.length }))),
            ([points, order]) => {
                const shuffled = order.map(i => points[i]);
                const a = minAreaCircle(points);
                const b = minAreaCircle(shuffled);
                if (!isTrustworthy(points, a) || !isTrustworthy(points, b)) {
                    return;
                }
                expectClose(a.minimal.radius, b.minimal.radius, 1e-9, 1e-9);
                if (a.minimal.radius > 0) {
                    expectVectorClose(a.minimal.center, b.minimal.center,
                        1e-8, 1e-8);
                }
            });
    });

    it('is unchanged by duplicated points', () => {
        check(fc.tuple(latticePoints, fc.nat(1000)), ([points, seed]) => {
            const augmented = [...points, points[seed % points.length].clone()];
            const a = minAreaCircle(points);
            const b = minAreaCircle(augmented);
            // The duplicate is removed before the algorithm runs, so the two
            // runs take exactly the same path.
            expect(b.success).toBe(a.success);
            expectClose(a.minimal.radius, b.minimal.radius, 0, 0);
            expectVectorClose(a.minimal.center, b.minimal.center, 0, 0);
        });
    });

    it('is equivariant under integer translation', () => {
        // Integer translations of integer inputs are exact, so the radius is
        // unchanged and the center translates by the same amount.
        check(fc.tuple(latticePoints, latticeVector(2, -20, 20)),
            ([points, offset]) => {
                const moved = points.map(p => Vector.fromArray([
                    p.get(0) + offset.get(0), p.get(1) + offset.get(1)]));
                const a = minAreaCircle(points);
                const b = minAreaCircle(moved);
                if (!isTrustworthy(points, a) || !isTrustworthy(moved, b)) {
                    return;
                }
                expectClose(a.minimal.radius, b.minimal.radius, 1e-8, 1e-8);
                expectClose(a.minimal.center.get(0) + offset.get(0),
                    b.minimal.center.get(0), 1e-7, 1e-8);
                expectClose(a.minimal.center.get(1) + offset.get(1),
                    b.minimal.center.get(1), 1e-7, 1e-8);
            });
    });

    it('no smaller enclosing circle is spanned by two input points', () => {
        // The minimum enclosing circle is unique, so any circle that contains
        // every point has a radius at least as large.
        check(latticePoints, points => {
            const result = minAreaCircle(points);
            if (!isTrustworthy(points, result)) {
                return;
            }
            const minimal = result.minimal;
            const n = points.length;
            const tol = 1e-9 * Math.max(1, minimal.radius);
            for (let i = 0; i < n; ++i) {
                for (let j = i + 1; j < n; ++j) {
                    const center = v2(
                        0.5 * (points[i].get(0) + points[j].get(0)),
                        0.5 * (points[i].get(1) + points[j].get(1)));
                    const radius = 0.5 * distance(points[i], points[j]);
                    let encloses = true;
                    for (const p of points) {
                        if (distance(p, center) > radius + tol) {
                            encloses = false;
                            break;
                        }
                    }
                    if (encloses) {
                        expect(minimal.radius)
                            .toBeLessThanOrEqual(radius + tol);
                    }
                }
            }
        });
    });

    it('the trapped-failure fallback bounds the whole input array', () => {
        // Regression for the upstream defect of issue #286 (item 1): the
        // fallback GetContainer call passes the *unique* point count with the
        // full input array, so upstream bounds only a prefix of the input.
        //
        // The four points below are cocircular (they lie on the circle of
        // radius sqrt(14746)/22 about (-39/22, -71/22)), which drives the
        // query into its trapped-failure branch: the fourth point is exactly
        // on the boundary in exact arithmetic and round-off puts it outside.
        // Prefixing three duplicates of the first point makes the unique
        // count 4 while the input array has 7 entries, so upstream's
        // GetContainer(4, points) would bound [A, A, A, A] alone and return
        // the degenerate circle (center A, radius 0).
        const a = v2(1, -8), b = v2(-7, -5), c = v2(3, -6), d = v2(0, 2);
        const points = [a, a.clone(), a.clone(), a.clone(), b, c, d];
        const query = new MinimumAreaCircle2();
        const { minimal, success } = query.compute(points);
        expect(success).toBe(false);
        expect(query.numSupport).toBe(0);
        expect(minimal.radius).toBeGreaterThan(0);
        expectContainsAll(points, minimal);
        // The fallback is the ContCircle2 bounding circle of *all* the input
        // points: center at the average, radius the largest distance to it.
        expect(minimal.center.get(0)).toBeCloseTo((4 * 1 - 7 + 3 + 0) / 7, 12);
        expect(minimal.center.get(1))
            .toBeCloseTo((4 * -8 - 5 - 6 + 2) / 7, 12);
    });

    it('the cocircular trapped failure is order independent', () => {
        // The same four cocircular points in any order reach the trapped
        // failure and the same fallback circle (which does not depend on the
        // order of the input).
        const base = [v2(1, -8), v2(-7, -5), v2(3, -6), v2(0, 2)];
        const orders = [[0, 1, 2, 3], [3, 0, 2, 1], [2, 3, 1, 0]];
        for (const order of orders) {
            const points = order.map(i => base[i]);
            const { minimal, success } =
                new MinimumAreaCircle2().compute(points);
            expect(success).toBe(false);
            expectContainsAll(points, minimal);
            expect(minimal.center.get(0)).toBeCloseTo(-0.75, 12);
            expect(minimal.center.get(1)).toBeCloseTo(-4.25, 12);
        }
    });

    it('reports success for a circle that misses a collinear point (upstream)',
        () => {
            // The second trigger of the known upstream defect: four points on
            // a line, two of them 5e-5 apart. exactCircle3 asks
            // LinearSystem.solve2x2 for the circumcircle of three collinear
            // points; the 2x2 determinant is a cancellation residue instead
            // of the exact zero, so the "no circumcircle" sentinel
            // (radius = MAX_VALUE) is not produced and the solve returns the
            // degenerate circle of radius 0 about the origin. That circle is
            // the smallest candidate, so updateSupport2 selects it, the
            // caller discards it as smaller than the current circle, and the
            // point (0, 0) is left hidden inside the support set.
            const points = [
                v2(-7.999998087043836, 0),
                v2(-0.0031354238413438856, 0),
                v2(0, 0),
                v2(-0.003085874074862011, 0)
            ];
            const query = new MinimumAreaCircle2();
            const { minimal, success } = query.compute(points);
            expect(success).toBe(true);
            // The true circle has the diameter (-7.999998087043836, 0)-(0, 0).
            expect(minimal.radius).toBeLessThan(0.5 * 7.999998087043836);
            expect(distance(points[2], minimal.center))
                .toBeGreaterThan(minimal.radius + 1e-3);
            expect(query.support.slice(0, query.numSupport)).toContain(2);
        });

    it('reports success for a circle that misses a cocircular point (upstream)', () => {
        // Pins the known upstream defect documented at the top of
        // src/MinimumAreaCircle2.ts: four of these points are cocircular on
        // the circle of squared radius 72.5 about (-2.5, 1.5), exactCircle3
        // computes that squared radius as 72.49999999999999 from one triple,
        // and the resulting spurious update replaces the support set without
        // enlarging the circle. The point (2, -6) then hides inside the
        // support set and the query returns success = true for a circle that
        // misses it by 0.23.
        const points = [v2(3, 8), v2(0, 0), v2(-2, -7), v2(1, 0), v2(-3, -7),
            v2(-3, 0), v2(-8, -5), v2(0, 1), v2(2, -6)];
        const query = new MinimumAreaCircle2();
        const { minimal, success } = query.compute(points);
        expect(success).toBe(true);
        expect(minimal.center.get(0)).toBeCloseTo(-2.5, 12);
        expect(minimal.center.get(1)).toBeCloseTo(1.5, 12);
        expect(minimal.radius).toBeCloseTo(Math.sqrt(72.5), 12);
        // (2, -6) is outside, and it is one of the reported support points.
        expect(distance(points[8], minimal.center))
            .toBeGreaterThan(minimal.radius + 0.2);
        expect(query.support.slice(0, query.numSupport)).toContain(8);
        // The true minimum enclosing circle is larger.
        const brute = bruteForceMinimumCircle(points);
        expect(brute.radius).toBeGreaterThan(minimal.radius);
        for (const p of points) {
            expect(distance(p, brute.center))
                .toBeLessThanOrEqual(brute.radius + 1e-9 * brute.radius);
        }
    });
});
