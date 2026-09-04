import { describe, it, expect } from 'vitest';
import { InscribedFixedAspectRectInQuad } from '../src/InscribedFixedAspectRectInQuad.js';
import type { InscribedFixedAspectRectInQuadResult } from '../src/InscribedFixedAspectRectInQuad.js';
import { Vector } from '../src/Vector.js';
import { check, expectClose, fc, scaled } from './helpers/arbitraries.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

// Test whether the point is inside or on the counterclockwise convex quad,
// allowing a small tolerance for floating-point roundoff.
function insideQuad(quad: Vector[], x: number, y: number, epsilon = 1e-9): boolean {
    for (let i = 0; i < 4; ++i) {
        const p0 = quad[i];
        const p1 = quad[(i + 1) % 4];
        const ex = p1.values[0] - p0.values[0];
        const ey = p1.values[1] - p0.values[1];
        const px = x - p0.values[0];
        const py = y - p0.values[1];
        if (ex * py - ey * px < -epsilon) {
            return false;
        }
    }
    return true;
}

function rectCorners(origin: Vector, width: number, height: number):
    Array<[number, number]> {
    const u = origin.values[0];
    const v = origin.values[1];
    return [[u, v], [u + width, v], [u + width, v + height], [u, v + height]];
}

describe('InscribedFixedAspectRectInQuad', () => {
    it('fills an axis-aligned rectangle when the aspect ratios match', () => {
        // The quad is the rectangle [0,4]x[0,2] whose aspect ratio is 2, so
        // the inscribed rectangle is the quad itself.
        const quad = [v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 2);
        expect(result.isUnique).toBe(true);
        expect(result.rectWidth).toBeCloseTo(4, 12);
        expect(result.rectHeight).toBeCloseTo(2, 12);
        expect(result.rectOrigin.values[0]).toBeCloseTo(0, 12);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 12);
    });

    it('reports infinitely many solutions when the rectangle can slide', () => {
        // The largest square inside [0,4]x[0,2] is 2x2 and can be translated
        // horizontally, so the solution is not unique.
        const quad = [v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 1);
        expect(result.isUnique).toBe(false);
        expect(result.rectWidth).toBeCloseTo(2, 12);
        expect(result.rectHeight).toBeCloseTo(2, 12);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 12);
    });

    it('inscribes the maximum square in a symmetric trapezoid', () => {
        // The trapezoid has the bottom edge from (0,0) to (4,0) and the top
        // edge from (3,2) to (1,2). The 2x2 square with origin (1,0) has its
        // top corners exactly at the top vertices, so it is the maximum.
        const quad = [v2(0, 0), v2(4, 0), v2(3, 2), v2(1, 2)];
        const result = InscribedFixedAspectRectInQuad.execute(quad, 1);
        expect(result.isUnique).toBe(true);
        expect(result.rectWidth).toBeCloseTo(2, 10);
        expect(result.rectHeight).toBeCloseTo(2, 10);
        expect(result.rectOrigin.values[0]).toBeCloseTo(1, 10);
        expect(result.rectOrigin.values[1]).toBeCloseTo(0, 10);
    });

    it('honors the aspect ratio w/h', () => {
        const quad = [v2(0, 0), v2(6, 0), v2(6, 6), v2(0, 6)];
        for (const aspectRatio of [0.5, 1, 2, 3]) {
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);
            expect(result.rectWidth / result.rectHeight).toBeCloseTo(aspectRatio, 10);
            // The limiting dimension is 6.
            expect(Math.max(result.rectWidth, result.rectHeight)).toBeCloseTo(6, 10);
        }
    });

    it('is translation equivariant', () => {
        const quad = [v2(0, 0), v2(5, 1), v2(4, 4), v2(-1, 3)];
        const shift = v2(10, -7);
        const shifted = quad.map(p =>
            v2(p.values[0] + shift.values[0], p.values[1] + shift.values[1]));
        const base = InscribedFixedAspectRectInQuad.execute(quad, 1.5);
        const moved = InscribedFixedAspectRectInQuad.execute(shifted, 1.5);
        expect(moved.rectWidth).toBeCloseTo(base.rectWidth, 9);
        expect(moved.rectHeight).toBeCloseTo(base.rectHeight, 9);
        expect(moved.rectOrigin.values[0])
            .toBeCloseTo(base.rectOrigin.values[0] + shift.values[0], 9);
        expect(moved.rectOrigin.values[1])
            .toBeCloseTo(base.rectOrigin.values[1] + shift.values[1], 9);
    });

    it('scales the rectangle with a uniform scaling of the quad', () => {
        const quad = [v2(0, 0), v2(5, 1), v2(4, 4), v2(-1, 3)];
        const scale = 3;
        const scaled = quad.map(p => v2(scale * p.values[0], scale * p.values[1]));
        const base = InscribedFixedAspectRectInQuad.execute(quad, 2);
        const bigger = InscribedFixedAspectRectInQuad.execute(scaled, 2);
        expect(bigger.rectWidth).toBeCloseTo(scale * base.rectWidth, 9);
        expect(bigger.rectHeight).toBeCloseTo(scale * base.rectHeight, 9);
        expect(bigger.rectOrigin.values[0])
            .toBeCloseTo(scale * base.rectOrigin.values[0], 9);
    });

    it('throws when the quad does not have 4 vertices', () => {
        expect(() => InscribedFixedAspectRectInQuad.execute(
            [v2(0, 0), v2(1, 0), v2(1, 1)], 1)).toThrow();
    });

    it('produces an inscribed rectangle that no larger one can beat', () => {
        const random = makeRandom(160934);
        for (let trial = 0; trial < 100; ++trial) {
            // Build a convex counterclockwise quadrilateral by perturbing the
            // corners of the unit square outward.
            const quad = [
                v2(-1 - random(), -1 - random()),
                v2(1 + random(), -1 - random()),
                v2(1 + random(), 1 + random()),
                v2(-1 - random(), 1 + random())
            ];
            const aspectRatio = 0.4 + 2 * random();
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);

            expect(result.rectWidth).toBeGreaterThan(0);
            expect(result.rectHeight).toBeCloseTo(result.rectWidth / aspectRatio, 9);

            // The rectangle is inscribed: all four corners are in the quad.
            for (const [x, y] of rectCorners(result.rectOrigin, result.rectWidth,
                result.rectHeight)) {
                expect(insideQuad(quad, x, y, 1e-8)).toBe(true);
            }

            // Growing the rectangle by 1% about its origin pushes at least one
            // corner outside, so the returned width is maximal.
            const grownWidth = 1.01 * result.rectWidth;
            const grownHeight = grownWidth / aspectRatio;
            let allInside = true;
            for (const [x, y] of rectCorners(result.rectOrigin, grownWidth, grownHeight)) {
                if (!insideQuad(quad, x, y, 1e-8)) {
                    allInside = false;
                    break;
                }
            }
            expect(allInside).toBe(false);
        }
    });

    it('never beats the maximum width found by a brute-force search', () => {
        const random = makeRandom(6022140);
        for (let trial = 0; trial < 20; ++trial) {
            const quad = [
                v2(-1 - random(), -1 - random()),
                v2(1 + random(), -1 - random()),
                v2(1 + random(), 1 + random()),
                v2(-1 - random(), 1 + random())
            ];
            const aspectRatio = 0.5 + random();
            const result = InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);

            // Sample origins on a grid and, for each, bisect on the width to
            // find the largest inscribed rectangle with that origin.
            let best = 0;
            const lo = -2.5;
            const hi = 2.5;
            const steps = 40;
            for (let i = 0; i <= steps; ++i) {
                const u = lo + (hi - lo) * i / steps;
                for (let j = 0; j <= steps; ++j) {
                    const v = lo + (hi - lo) * j / steps;
                    if (!insideQuad(quad, u, v)) {
                        continue;
                    }
                    let low = 0;
                    let high = 6;
                    for (let k = 0; k < 40; ++k) {
                        const mid = 0.5 * (low + high);
                        const corners = rectCorners(v2(u, v), mid, mid / aspectRatio);
                        const fits = corners.every(([x, y]) => insideQuad(quad, x, y));
                        if (fits) {
                            low = mid;
                        }
                        else {
                            high = mid;
                        }
                    }
                    best = Math.max(best, low);
                }
            }
            // The grid search is coarse, so it can only give a lower bound on
            // the optimum. The query result must be at least as large.
            expect(result.rectWidth).toBeGreaterThan(best - 1e-6);
        }
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// The upstream algorithm reduces the containment of the rectangle to four
// linear constraints by choosing, for each edge normal, the rectangle corner
// that is extreme in the -N direction (the floor(2*angle/pi) case analysis).
// The reference below does not make that reduction: it writes all sixteen
// constraints "corner k is on the inner side of edge i" as linear inequalities
// in (u, v, w) and solves the linear program "maximize w" by enumerating the
// vertices of the feasible polyhedron (every triple of constraint planes).
// That is an independent derivation of both the constraint set and the
// optimum, so a mistake in the case analysis or in the interval intersection
// would show up as a different maximum width.
// ---------------------------------------------------------------------------

interface LinearConstraint { a: [number, number, number]; c: number }

// The sixteen constraints Dot(N_i, R_k - V_i) >= 0 written as a.x >= c with
// x = (u, v, w) and the corners
//   R0 = (u, v), R1 = (u+w, v), R2 = (u+w, v+w/r), R3 = (u, v+w/r).
function allConstraints(quad: Vector[], r: number): LinearConstraint[] {
    const constraints: LinearConstraint[] = [];
    for (let i = 0; i < 4; ++i) {
        const V = quad[i];
        const W = quad[(i + 1) % 4];
        // Inner-pointing normal of the edge <V, W> of a counterclockwise quad.
        const n0 = -(W.values[1] - V.values[1]);
        const n1 = W.values[0] - V.values[0];
        const c = n0 * V.values[0] + n1 * V.values[1];
        const wCoefficients = [0, n0, n0 + n1 / r, n1 / r];
        for (const wc of wCoefficients) {
            constraints.push({ a: [n0, n1, wc], c });
        }
    }
    return constraints;
}

// Maximize w over the feasible polyhedron by enumerating its vertices. The
// returned 'origins' are the (u,v) of every feasible vertex whose w is within
// a whisker of the optimum, so a caller can tell whether the optimal origin is
// determined or whether the rectangle can slide.
function referenceOptimum(quad: Vector[], r: number):
    { width: number, origins: Array<[number, number]> } {
    const vertices = feasibleVertices(quad, r);
    let width = Number.NEGATIVE_INFINITY;
    for (const vertex of vertices) { width = Math.max(width, vertex[2]); }
    const origins: Array<[number, number]> = [];
    for (const vertex of vertices) {
        if (vertex[2] >= width - 1e-7 * Math.max(1, Math.abs(width))) {
            origins.push([vertex[0], vertex[1]]);
        }
    }
    return { width, origins };
}

// True when every near-optimal vertex has the same (u,v), i.e. the maximizing
// rectangle origin is determined rather than sliding.
function optimalOriginIsDetermined(quad: Vector[], r: number): boolean {
    const { origins } = referenceOptimum(quad, r);
    for (const origin of origins) {
        if (Math.abs(origin[0] - origins[0][0]) > 1e-7
            || Math.abs(origin[1] - origins[0][1]) > 1e-7) {
            return false;
        }
    }
    return origins.length > 0;
}

function referenceMaxWidth(quad: Vector[], r: number): number {
    return referenceOptimum(quad, r).width;
}

// The vertices of the feasible polyhedron: every triple of constraint planes
// whose intersection satisfies all sixteen constraints.
function feasibleVertices(quad: Vector[], r: number): Array<[number, number, number]> {
    const constraints = allConstraints(quad, r);
    const scale = Math.max(1, ...quad.map(p =>
        Math.max(Math.abs(p.values[0]), Math.abs(p.values[1]))));
    const vertices: Array<[number, number, number]> = [];
    for (let i = 0; i < constraints.length; ++i) {
        for (let j = i + 1; j < constraints.length; ++j) {
            for (let k = j + 1; k < constraints.length; ++k) {
                const m = [constraints[i], constraints[j], constraints[k]];
                const A = m.map(t => t.a);
                const det =
                    A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
                    - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
                    + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
                if (Math.abs(det) < 1e-9 * scale * scale * scale) { continue; }
                const b = m.map(t => t.c);
                const solve = (column: number): number => {
                    const M = A.map((row, rowIndex) =>
                        row.map((value, col) => (col === column ? b[rowIndex] : value)));
                    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
                        - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                        + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
                };
                const x = [solve(0), solve(1), solve(2)];
                let feasible = true;
                for (const constraint of constraints) {
                    const value = constraint.a[0] * x[0] + constraint.a[1] * x[1]
                        + constraint.a[2] * x[2] - constraint.c;
                    const magnitude = Math.max(1, Math.abs(constraint.c),
                        Math.abs(constraint.a[0] * x[0]) + Math.abs(constraint.a[1] * x[1])
                        + Math.abs(constraint.a[2] * x[2]));
                    if (value < -1e-8 * magnitude) { feasible = false; break; }
                }
                if (feasible) { vertices.push([x[0], x[1], x[2]]); }
            }
        }
    }
    return vertices;
}

// A convex counterclockwise quadrilateral: four points at increasing angles on
// an ellipse, mapped by a rotation and a translation.
// The coordinates are drawn from a fine grid rather than a coarse one so that
// the quads are in general position: symmetric or axis-parallel quads make the
// optimum touch all four edges at once, which is exactly the degenerate case
// the upstream solver mishandles (pinned separately below).
const convexQuad = fc.record({
    angles: fc.uniqueArray(fc.integer({ min: 0, max: 359 }),
        { minLength: 4, maxLength: 4 }),
    radii: fc.array(scaled(2, 6, 4096), { minLength: 4, maxLength: 4 }),
    tx: scaled(-3, 3, 512),
    ty: scaled(-3, 3, 512),
    aspectRatio: scaled(0.3, 3, 512)
}).map(({ angles, radii, tx, ty, aspectRatio }) => {
    const sorted = [...angles].sort((a, b) => a - b);
    const quad = sorted.map((degrees, k) => {
        const t = (Math.PI * degrees) / 180;
        return v2(radii[k] * Math.cos(t) + tx, radii[k] * Math.sin(t) + ty);
    });
    return { quad, aspectRatio };
}).filter(({ quad }) => {
    // Strictly convex, counterclockwise, and not too thin (the reference LP
    // and the algorithm both lose accuracy on slivers).
    for (let i = 0; i < 4; ++i) {
        const p0 = quad[i], p1 = quad[(i + 1) % 4], p2 = quad[(i + 2) % 4];
        const cross = (p1.values[0] - p0.values[0]) * (p2.values[1] - p0.values[1])
            - (p1.values[1] - p0.values[1]) * (p2.values[0] - p0.values[0]);
        if (cross < 1) { return false; }
    }
    // Skip the configurations whose constraint planes are coplanar; that
    // limitation is pinned by its own test below.
    return isSolvableByUpstream(quad);
});

// Run the query, tolerating the two configurations the upstream solver cannot
// handle (both pinned by their own tests below). The helper asserts that a
// thrown error really is one of them, so a port defect that threw on ordinary
// input would still fail the property.
function executeOrDegenerate(quad: Vector[], aspectRatio: number):
    InscribedFixedAspectRectInQuadResult | null {
    try {
        return InscribedFixedAspectRectInQuad.execute(quad, aspectRatio);
    } catch (error) {
        expect(degenerateForUpstream(quad, aspectRatio),
            `unexpected throw: ${(error as Error).message}`).toBe(true);
        return null;
    }
}

// The two known upstream failure modes:
//  (a) three of the four constraint normals share a case index j, so they are
//      coplanar and the alpha assert fires ("Unexpected condition");
//  (b) the optimum makes all four upstream constraints tight, so the feasible
//      interval on the line of two of the planes degenerates to a point and
//      round-off can report it as empty for both pairings ("Unexpected
//      interval intersection type").
function degenerateForUpstream(quad: Vector[], r: number): boolean {
    if (!isSolvableByUpstream(quad)) { return true; }
    const { width, origins } = referenceOptimum(quad, r);
    if (!Number.isFinite(width) || origins.length === 0) { return true; }
    const x = [origins[0][0], origins[0][1], width];
    let tight = 0;
    for (const constraint of upstreamConstraints(quad, r)) {
        const value = constraint.a[0] * x[0] + constraint.a[1] * x[1]
            + constraint.a[2] * x[2] - constraint.c;
        if (Math.abs(value) <= 1e-9 * Math.max(1, Math.abs(constraint.c))) { ++tight; }
    }
    return tight >= 4;
}

// The four constraints upstream builds: one per edge, using the rectangle
// corner selected by the quadrant of the edge normal.
function upstreamConstraints(quad: Vector[], r: number): LinearConstraint[] {
    const quadrants = normalQuadrants(quad);
    const constraints: LinearConstraint[] = [];
    for (let i = 0; i < 4; ++i) {
        const V = quad[i], W = quad[(i + 1) % 4];
        const n0 = V.values[1] - W.values[1];
        const n1 = W.values[0] - V.values[0];
        const c = n0 * V.values[0] + n1 * V.values[1];
        const wCoefficients = [0, n0, n0 + n1 / r, n1 / r];
        constraints.push({ a: [n0, n1, wCoefficients[quadrants[i]]], c });
    }
    return constraints;
}

// Whether the upstream solver can handle the quad: constraint normals sharing
// a case index j lie in a common plane of (u,v,w)-space, so if normals 0, 2
// and one of 1, 3 share a quadrant the corresponding alpha is zero and Execute
// asserts "Unexpected condition". That limitation is pinned by its own test
// below. Note that a quad with an axis-parallel edge sits exactly on a
// quadrant boundary, so a translation or a scaling can round the normal to the
// other side and change the answer of this predicate; callers re-test the
// transformed quad.
function isSolvableByUpstream(quad: Vector[]): boolean {
    const q = normalQuadrants(quad);
    return !((q[0] === q[1] && q[1] === q[2]) || (q[0] === q[2] && q[2] === q[3]));
}

// The quadrant index floor(2*angle/pi) of each inner edge normal, the value
// upstream uses to pick the binding rectangle corner for that edge.
function normalQuadrants(quad: Vector[]): number[] {
    const quadrants: number[] = [];
    for (let i = 0; i < 4; ++i) {
        const V = quad[i], W = quad[(i + 1) % 4];
        // Perp(V - W) is the inner-pointing normal of the edge <V, W>.
        const n0 = V.values[1] - W.values[1];
        const n1 = W.values[0] - V.values[0];
        let angle = Math.atan2(n1, n0);
        if (angle < 0) { angle += 2 * Math.PI; }
        quadrants.push(Math.floor((2 * angle) / Math.PI) + 0);   // avoid -0
    }
    return quadrants;
}

describe('InscribedFixedAspectRectInQuad verification', () => {
    it('the reported rectangle is inscribed in the quad', () => {
        check(convexQuad, ({ quad, aspectRatio }) => {
            const result = executeOrDegenerate(quad, aspectRatio);
            if (result === null) { return; }
            expect(result.rectWidth).toBeGreaterThan(0);
            expect(result.rectHeight).toBe(result.rectWidth / aspectRatio);
            const scale = Math.max(1, ...quad.map(p =>
                Math.max(Math.abs(p.values[0]), Math.abs(p.values[1]))));
            for (const [x, y] of rectCorners(result.rectOrigin, result.rectWidth,
                result.rectHeight)) {
                expect(insideQuad(quad, x, y, 1e-9 * scale * scale)).toBe(true);
            }
        });
    });

    it('the width is the linear-programming optimum', () => {
        check(convexQuad, ({ quad, aspectRatio }) => {
            const result = executeOrDegenerate(quad, aspectRatio);
            if (result === null) { return; }
            const reference = referenceMaxWidth(quad, aspectRatio);
            expect(Number.isFinite(reference)).toBe(true);
            // Both computations solve the same LP by different routes; the
            // tolerance covers the conditioning of the 3x3 solves.
            expectClose(result.rectWidth, reference, 1e-7, 1e-7);
        });
    });

    it('no wider rectangle of the same aspect ratio fits anywhere', () => {
        check(convexQuad, ({ quad, aspectRatio }) => {
            const result = executeOrDegenerate(quad, aspectRatio);
            if (result === null) { return; }
            const width = result.rectWidth * 1.02;
            const height = width / aspectRatio;
            // Scan candidate origins over the bounding box of the quad.
            const xs = quad.map(p => p.values[0]);
            const ys = quad.map(p => p.values[1]);
            const x0 = Math.min(...xs), x1 = Math.max(...xs);
            const y0 = Math.min(...ys), y1 = Math.max(...ys);
            const steps = 24;
            for (let i = 0; i <= steps; ++i) {
                const u = x0 + ((x1 - x0) * i) / steps;
                for (let j = 0; j <= steps; ++j) {
                    const v = y0 + ((y1 - y0) * j) / steps;
                    const fits = rectCorners(v2(u, v), width, height)
                        .every(([x, y]) => insideQuad(quad, x, y));
                    expect(fits).toBe(false);
                }
            }
        }, 40);
    });

    it('is equivariant under translation and uniform scaling', () => {
        check(fc.tuple(convexQuad, scaled(-6, 6, 24), scaled(-6, 6, 24),
            scaled(0.25, 4, 16)), ([{ quad, aspectRatio }, tx, ty, s]) => {
            const base = executeOrDegenerate(quad, aspectRatio);
            if (base === null) { return; }
            const translated = quad.map(p => v2(p.values[0] + tx, p.values[1] + ty));
            const scaledQuad = quad.map(p => v2(s * p.values[0], s * p.values[1]));
            // An axis-parallel edge can round to the other side of a quadrant
            // boundary under the transform and land in a configuration the
            // upstream solver cannot handle.
            const moved = executeOrDegenerate(translated, aspectRatio);
            if (moved === null) { return; }
            const scaleTol = 1e-9 * Math.max(1, Math.abs(tx), Math.abs(ty));
            expectClose(moved.rectWidth, base.rectWidth, 1e-9, 1e-9);
            if (base.isUnique && optimalOriginIsDetermined(quad, aspectRatio)) {
                // The origin is comparable only when the optimum is attained
                // at a single (u,v): a sliding rectangle may be reported at
                // either end, and a near-tie between two distinct optimal
                // vertices can switch under a translation.
                expectClose(moved.rectOrigin.values[0], base.rectOrigin.values[0] + tx,
                    1e-8 + scaleTol, 1e-9);
                expectClose(moved.rectOrigin.values[1], base.rectOrigin.values[1] + ty,
                    1e-8 + scaleTol, 1e-9);
            }
            const bigger = executeOrDegenerate(scaledQuad, aspectRatio);
            if (bigger === null) { return; }
            expectClose(bigger.rectWidth, s * base.rectWidth, 1e-8, 1e-9);
            expectClose(bigger.rectHeight, s * base.rectHeight, 1e-8, 1e-9);
        });
    });

    it('sliding solutions are reported as non-unique', () => {
        // A 4-by-2 axis-aligned quad and aspect ratio 1: the maximal square is
        // 2-by-2 and can slide horizontally, so the solution is not unique.
        const quad = [v2(0, 0), v2(4, 0), v2(4, 2), v2(0, 2)];
        const sliding = InscribedFixedAspectRectInQuad.execute(quad, 1);
        expect(sliding.isUnique).toBe(false);
        expect(sliding.rectWidth).toBeCloseTo(2, 12);
        expect(sliding.rectHeight).toBeCloseTo(2, 12);
        // Matching aspect ratio: the rectangle is the quad itself and unique.
        const exact = InscribedFixedAspectRectInQuad.execute(quad, 2);
        expect(exact.rectWidth).toBeCloseTo(4, 12);
        expect(exact.rectHeight).toBeCloseTo(2, 12);
        // A trapezoid pins the square: unique solution.
        const trapezoid = [v2(0, 0), v2(4, 0), v2(3, 2), v2(1, 2)];
        const pinned = InscribedFixedAspectRectInQuad.execute(trapezoid, 1);
        expect(pinned.isUnique).toBe(true);
        expect(pinned.rectWidth)
            .toBeCloseTo(referenceMaxWidth(trapezoid, 1), 9);
    });

    it('documented upstream limitation: three normals in one quadrant throw', () => {
        // <(0,3),(0,2),(1,1),(3,0)> is a strictly convex counterclockwise quad
        // (area 2), so a maximum inscribed rectangle exists for every aspect
        // ratio. Its normals 0, 1 and 2 share the quadrant index 0, so the
        // three constraint planes have a zero third component, the line of
        // planes 0 and 2 is parallel to plane 1 (alpha1 == 0), and upstream's
        // LogAssert fires. The port preserves that behavior.
        const quad = [v2(0, 3), v2(0, 2), v2(1, 1), v2(3, 0)];
        for (let i = 0; i < 4; ++i) {
            const p0 = quad[i], p1 = quad[(i + 1) % 4], p2 = quad[(i + 2) % 4];
            const cross = (p1.values[0] - p0.values[0]) * (p2.values[1] - p0.values[1])
                - (p1.values[1] - p0.values[1]) * (p2.values[0] - p0.values[0]);
            expect(cross).toBeGreaterThan(0);
        }
        expect(normalQuadrants(quad).slice(0, 3)).toEqual([0, 0, 0]);
        for (const aspectRatio of [0.5, 1, 2]) {
            expect(() => InscribedFixedAspectRectInQuad.execute(quad, aspectRatio))
                .toThrow('Unexpected condition.');
        }
        // A rectangle of positive width does fit, so the failure is a solver
        // limitation and not an empty feasible set.
        expect(referenceMaxWidth(quad, 1)).toBeGreaterThan(0);
    });

    it('documented upstream limitation: a four-tight optimum can be lost', () => {
        // For this convex counterclockwise quad the maximum rectangle touches
        // all four edges, so the feasible interval along the line of planes 0
        // and 2 degenerates to a single point t*. Upstream compares the two
        // half-line endpoints without any tolerance, and here they straddle t*
        // by one ulp in the wrong order (end1 - end3 = 4e-16 > 0), so the
        // intersection is reported as empty; the complementary pairing is
        // empty for the same reason and Execute reports "Unexpected interval
        // intersection type". The port preserves that behavior.
        const quad = [
            v2(-2.2928932188134525, 0.4393398282201786),
            v2(-2.034074173710932, 1.888228567653781),
            v2(-2.2928932188134525, 2.560660171779821),
            v2(-3.965925826289068, 1.8882285676537809)
        ];
        for (let i = 0; i < 4; ++i) {
            const p0 = quad[i], p1 = quad[(i + 1) % 4], p2 = quad[(i + 2) % 4];
            const cross = (p1.values[0] - p0.values[0]) * (p2.values[1] - p0.values[1])
                - (p1.values[1] - p0.values[1]) * (p2.values[0] - p0.values[0]);
            expect(cross).toBeGreaterThan(0);
        }
        // The constraint planes are independent, so this is not the coplanar
        // failure mode pinned above.
        expect(isSolvableByUpstream(quad)).toBe(true);
        expect(() => InscribedFixedAspectRectInQuad.execute(quad, 0.46875))
            .toThrow('Unexpected interval intersection type.');
        // The linear program is feasible and bounded: the optimum exists and
        // makes all four upstream constraints tight.
        const optimum = referenceOptimum(quad, 0.46875);
        expect(optimum.width).toBeGreaterThan(0.65);
        expect(degenerateForUpstream(quad, 0.46875)).toBe(true);
    });

    it('rejects a quad that does not have four vertices', () => {
        check(fc.integer({ min: 0, max: 6 }).filter(n => n !== 4), (n) => {
            const quad: Vector[] = [];
            for (let i = 0; i < n; ++i) { quad.push(v2(i, i * i)); }
            expect(() => InscribedFixedAspectRectInQuad.execute(quad, 1))
                .toThrow('The quadrilateral must have 4 vertices.');
        });
    });
});
