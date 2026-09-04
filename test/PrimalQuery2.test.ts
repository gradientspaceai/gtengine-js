import { describe, it, expect } from 'vitest';
import { PrimalQuery2, PrimalQuery2OrderType } from '../src/PrimalQuery2.js';
import { Vector } from '../src/Vector.js';
import { check, fc } from './helpers/arbitraries.js';
import { inCircle2, orient2 } from './helpers/exact.js';

const v2 = (x: number, y: number): Vector => Vector.fromArray([x, y]);

// Deterministic LCG so the randomized cross-checks are reproducible.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

function sign(value: number): number {
    return (value > 0 ? +1 : (value < 0 ? -1 : 0));
}

describe('PrimalQuery2', () => {
    it('has a default constructor with no vertices', () => {
        const query = new PrimalQuery2();
        expect(query.getNumVertices()).toBe(0);
        expect(query.getVertices()).toEqual([]);
    });

    it('supports set() member access', () => {
        const vertices = [v2(0, 0), v2(1, 0)];
        const query = new PrimalQuery2();
        query.set(vertices.length, vertices);
        expect(query.getNumVertices()).toBe(2);
        expect(query.getVertices()).toBe(vertices);
    });

    describe('toLine', () => {
        // The line has origin V0 = (0,0) and direction <V0,V1> = (1,0), the
        // positive x-axis. Points with y < 0 are on the right of the line.
        const vertices = [
            v2(0, 0),   // 0: V0
            v2(1, 0),   // 1: V1
            v2(2, -3),  // 2: right of the line
            v2(2, 3),   // 3: left of the line
            v2(5, 0),   // 4: on the line, beyond V1
            v2(-5, 0)   // 5: on the line, before V0
        ];
        const query = new PrimalQuery2(vertices.length, vertices);

        it('classifies points relative to the directed line', () => {
            expect(query.toLine(2, 0, 1)).toBe(+1);
            expect(query.toLine(3, 0, 1)).toBe(-1);
            expect(query.toLine(4, 0, 1)).toBe(0);
            expect(query.toLine(5, 0, 1)).toBe(0);
            expect(query.toLine(0, 0, 1)).toBe(0);
            expect(query.toLine(1, 0, 1)).toBe(0);
        });

        it('accepts a test point as well as a vertex index', () => {
            expect(query.toLine(v2(2, -3), 0, 1)).toBe(+1);
            expect(query.toLine(v2(2, 3), 0, 1)).toBe(-1);
            expect(query.toLine(v2(0.5, 0), 0, 1)).toBe(0);
        });

        it('flips sign when the line direction is reversed', () => {
            expect(query.toLine(2, 1, 0)).toBe(-1);
            expect(query.toLine(3, 1, 0)).toBe(+1);
        });

        it('agrees with the sign of the cross product (randomized)', () => {
            const random = makeRandom(112233);
            for (let trial = 0; trial < 500; ++trial) {
                // Small integer coordinates keep the determinant exact in
                // double precision, matching the upstream exact-arithmetic
                // instantiation.
                const points = [
                    v2(Math.floor(21 * random()) - 10, Math.floor(21 * random()) - 10),
                    v2(Math.floor(21 * random()) - 10, Math.floor(21 * random()) - 10),
                    v2(Math.floor(21 * random()) - 10, Math.floor(21 * random()) - 10)
                ];
                const q = new PrimalQuery2(3, points);
                const x0 = points[2].values[0] - points[0].values[0];
                const y0 = points[2].values[1] - points[0].values[1];
                const x1 = points[1].values[0] - points[0].values[0];
                const y1 = points[1].values[1] - points[0].values[1];
                expect(q.toLine(2, 0, 1)).toBe(sign(x0 * y1 - x1 * y0));
            }
        });
    });

    describe('toLineWithOrder', () => {
        // V0 = (0,0), V1 = (4,0).
        const vertices = [v2(0, 0), v2(4, 0)];
        const query = new PrimalQuery2(vertices.length, vertices);

        it('returns +1/+3 for a point on the right of the line', () => {
            const result = query.toLineWithOrder(v2(2, -1), 0, 1);
            expect(result.sign).toBe(+1);
            expect(result.order).toBe(+3);
        });

        it('returns -1/-3 for a point on the left of the line', () => {
            const result = query.toLineWithOrder(v2(2, 1), 0, 1);
            expect(result.sign).toBe(-1);
            expect(result.order).toBe(-3);
        });

        it('returns order -1 for P = V0 and +1 for P = V1', () => {
            expect(query.toLineWithOrder(v2(0, 0), 0, 1)).toEqual({ sign: 0, order: -1 });
            expect(query.toLineWithOrder(v2(4, 0), 0, 1)).toEqual({ sign: 0, order: +1 });
        });

        it('returns order -2 for a collinear point before V0', () => {
            expect(query.toLineWithOrder(v2(-3, 0), 0, 1)).toEqual({ sign: 0, order: -2 });
        });

        it('preserves the upstream swap of the collinear orders 0 and +2', () => {
            // Upstream compares Dot(P-V0,V1-V0) against |P-V0|^2 rather than
            // |V1-V0|^2, which swaps the interior (0) and beyond-V1 (+2)
            // classifications. The port keeps that behavior; see the port
            // note on PrimalQuery2ToLineOrderResult.
            //
            // P = (2,0) is interior to [V0,V1] but yields order +2.
            expect(query.toLineWithOrder(v2(2, 0), 0, 1)).toEqual({ sign: 0, order: +2 });
            // P = (9,0) is beyond V1 but yields order 0.
            expect(query.toLineWithOrder(v2(9, 0), 0, 1)).toEqual({ sign: 0, order: 0 });
        });

        it('agrees with toLine on the sign', () => {
            const random = makeRandom(445566);
            for (let trial = 0; trial < 300; ++trial) {
                const test = v2(Math.floor(11 * random()) - 5, Math.floor(11 * random()) - 5);
                expect(query.toLineWithOrder(test, 0, 1).sign).toBe(query.toLine(test, 0, 1));
            }
        });
    });

    describe('toTriangle', () => {
        // Counterclockwise triangle V0 = (0,0), V1 = (4,0), V2 = (0,4).
        const vertices = [
            v2(0, 0), v2(4, 0), v2(0, 4),
            v2(1, 1),    // 3: strictly inside
            v2(5, 5),    // 4: strictly outside
            v2(2, 0),    // 5: on the edge <V0,V1>
            v2(2, 2),    // 6: on the edge <V1,V2>
            v2(-1, 2)    // 7: outside, left of the edge <V2,V0>
        ];
        const query = new PrimalQuery2(vertices.length, vertices);

        it('returns -1 for interior points', () => {
            expect(query.toTriangle(3, 0, 1, 2)).toBe(-1);
            expect(query.toTriangle(v2(0.5, 0.5), 0, 1, 2)).toBe(-1);
        });

        it('returns +1 for exterior points', () => {
            expect(query.toTriangle(4, 0, 1, 2)).toBe(+1);
            expect(query.toTriangle(7, 0, 1, 2)).toBe(+1);
            expect(query.toTriangle(v2(-1, -1), 0, 1, 2)).toBe(+1);
        });

        it('returns 0 for points on the triangle boundary', () => {
            expect(query.toTriangle(5, 0, 1, 2)).toBe(0);
            expect(query.toTriangle(6, 0, 1, 2)).toBe(0);
            for (const i of [0, 1, 2]) {
                expect(query.toTriangle(i, 0, 1, 2)).toBe(0);
            }
        });

        it('agrees with a barycentric containment test (randomized)', () => {
            const random = makeRandom(987654);
            const triangle = [v2(0, 0), v2(8, 0), v2(3, 7)];
            const query2 = new PrimalQuery2(3, triangle);
            const area = (triangle[1].values[0] - triangle[0].values[0])
                * (triangle[2].values[1] - triangle[0].values[1])
                - (triangle[1].values[1] - triangle[0].values[1])
                * (triangle[2].values[0] - triangle[0].values[0]);
            expect(area).toBeGreaterThan(0);

            for (let trial = 0; trial < 500; ++trial) {
                const p = v2(Math.floor(13 * random()) - 2, Math.floor(11 * random()) - 2);
                const sub = (a: Vector, b: Vector): number[] =>
                    [a.values[0] - b.values[0], a.values[1] - b.values[1]];
                const cross = (a: number[], b: number[]): number => a[0] * b[1] - a[1] * b[0];
                const s0 = cross(sub(triangle[1], triangle[0]), sub(p, triangle[0]));
                const s1 = cross(sub(triangle[2], triangle[1]), sub(p, triangle[1]));
                const s2 = cross(sub(triangle[0], triangle[2]), sub(p, triangle[2]));
                let expected: number;
                if (s0 > 0 && s1 > 0 && s2 > 0) {
                    expected = -1;
                }
                else if (s0 < 0 || s1 < 0 || s2 < 0) {
                    expected = +1;
                }
                else {
                    expected = 0;
                }
                expect(query2.toTriangle(p, 0, 1, 2)).toBe(expected);
            }
        });
    });

    describe('toCircumcircle', () => {
        // Counterclockwise triangle inscribed in the unit circle.
        const vertices = [
            v2(1, 0), v2(0, 1), v2(-1, 0),
            v2(0, 0),     // 3: the circumcenter, strictly inside
            v2(5, 5),     // 4: far outside
            v2(0, -1),    // 5: exactly on the circumcircle
            v2(0, 0.5)    // 6: inside
        ];
        const query = new PrimalQuery2(vertices.length, vertices);

        it('returns -1 inside, +1 outside and 0 on the circumcircle', () => {
            expect(query.toCircumcircle(3, 0, 1, 2)).toBe(-1);
            expect(query.toCircumcircle(6, 0, 1, 2)).toBe(-1);
            expect(query.toCircumcircle(4, 0, 1, 2)).toBe(+1);
            expect(query.toCircumcircle(5, 0, 1, 2)).toBe(0);
        });

        it('returns 0 for the triangle vertices themselves', () => {
            for (const i of [0, 1, 2]) {
                expect(query.toCircumcircle(i, 0, 1, 2)).toBe(0);
            }
        });

        it('agrees with an explicit circumcenter/radius test (randomized)', () => {
            const random = makeRandom(24681012);
            let numInside = 0;
            let numOutside = 0;
            for (let trial = 0; trial < 400; ++trial) {
                // A counterclockwise triangle with integer coordinates.
                const t0 = v2(Math.floor(9 * random()) - 4, Math.floor(9 * random()) - 4);
                const t1 = v2(Math.floor(9 * random()) - 4, Math.floor(9 * random()) - 4);
                const t2 = v2(Math.floor(9 * random()) - 4, Math.floor(9 * random()) - 4);
                const ax = t1.values[0] - t0.values[0];
                const ay = t1.values[1] - t0.values[1];
                const bx = t2.values[0] - t0.values[0];
                const by = t2.values[1] - t0.values[1];
                const twiceArea = ax * by - ay * bx;
                if (twiceArea <= 0) {
                    // Skip degenerate and clockwise triangles; the query
                    // documents counterclockwise ordering.
                    continue;
                }

                // Circumcenter of the triangle, computed independently.
                const d = 2 * twiceArea;
                const sa = ax * ax + ay * ay;
                const sb = bx * bx + by * by;
                const cx = t0.values[0] + (by * sa - ay * sb) / d;
                const cy = t0.values[1] + (ax * sb - bx * sa) / d;
                const radius = Math.hypot(t0.values[0] - cx, t0.values[1] - cy);

                const p = v2(Math.floor(13 * random()) - 6, Math.floor(13 * random()) - 6);
                const distance = Math.hypot(p.values[0] - cx, p.values[1] - cy);
                if (Math.abs(distance - radius) < 1e-9) {
                    // Skip the on-circle cases where the independent test is
                    // itself ambiguous.
                    continue;
                }

                const points = [t0, t1, t2];
                const q = new PrimalQuery2(3, points);
                const expected = (distance < radius ? -1 : +1);
                if (expected < 0) {
                    ++numInside;
                }
                else {
                    ++numOutside;
                }
                expect(q.toCircumcircle(p, 0, 1, 2)).toBe(expected);
            }
            expect(numInside).toBeGreaterThan(0);
            expect(numOutside).toBeGreaterThan(0);
        });
    });

    describe('toLineExtended', () => {
        const query = new PrimalQuery2();
        const Q0 = v2(0, 0);
        const Q1 = v2(4, 0);

        it('detects the degenerate and coincident cases', () => {
            expect(query.toLineExtended(v2(1, 1), Q0, Q0))
                .toBe(PrimalQuery2OrderType.Q0_EQUALS_Q1);
            expect(query.toLineExtended(Q0, Q0, Q1))
                .toBe(PrimalQuery2OrderType.P_EQUALS_Q0);
            expect(query.toLineExtended(Q1, Q0, Q1))
                .toBe(PrimalQuery2OrderType.P_EQUALS_Q1);
        });

        it('classifies noncollinear points by triangle orientation', () => {
            // <P,Q0,Q1> is counterclockwise when P is above the x-axis: the
            // vertices (2,1) -> (0,0) -> (4,0) wind counterclockwise, while
            // (2,-1) -> (0,0) -> (4,0) wind clockwise.
            expect(query.toLineExtended(v2(2, 1), Q0, Q1))
                .toBe(PrimalQuery2OrderType.POSITIVE);
            expect(query.toLineExtended(v2(2, -1), Q0, Q1))
                .toBe(PrimalQuery2OrderType.NEGATIVE);
        });

        it('classifies collinear points by their ordering on the line', () => {
            expect(query.toLineExtended(v2(-2, 0), Q0, Q1))
                .toBe(PrimalQuery2OrderType.COLLINEAR_LEFT);
            expect(query.toLineExtended(v2(9, 0), Q0, Q1))
                .toBe(PrimalQuery2OrderType.COLLINEAR_RIGHT);
            expect(query.toLineExtended(v2(2, 0), Q0, Q1))
                .toBe(PrimalQuery2OrderType.COLLINEAR_CONTAIN);
        });

        it('agrees with toLine on the noncollinear sign (randomized)', () => {
            const random = makeRandom(13579);
            for (let trial = 0; trial < 300; ++trial) {
                const P = v2(Math.floor(11 * random()) - 5, Math.floor(11 * random()) - 5);
                const q0 = v2(Math.floor(11 * random()) - 5, Math.floor(11 * random()) - 5);
                const q1 = v2(Math.floor(11 * random()) - 5, Math.floor(11 * random()) - 5);
                const points = [P, q0, q1];
                const q = new PrimalQuery2(3, points);
                const order = q.toLineExtended(P, q0, q1);
                if (order === PrimalQuery2OrderType.POSITIVE) {
                    // toLineExtended forms the determinant with the opposite
                    // sign convention of toLine: its det is
                    // Cross(Q1-Q0, P-Q0) while toLine computes
                    // Cross(P-V0, V1-V0). So POSITIVE (P left of the line)
                    // corresponds to toLine returning -1.
                    expect(q.toLine(0, 1, 2)).toBe(-1);
                }
                else if (order === PrimalQuery2OrderType.NEGATIVE) {
                    expect(q.toLine(0, 1, 2)).toBe(+1);
                }
                else if (order === PrimalQuery2OrderType.COLLINEAR_LEFT
                    || order === PrimalQuery2OrderType.COLLINEAR_RIGHT
                    || order === PrimalQuery2OrderType.COLLINEAR_CONTAIN) {
                    expect(q.toLine(0, 1, 2)).toBe(0);
                }
            }
        });
    });
});

// ---------------------------------------------------------------------------
// Independent verification pass (VERIFYING.md).
//
// PrimalQuery2 was ported for T = number only (no BSNumber/BSRational
// instantiation), so the port is exact exactly when every intermediate value
// is representable in binary64. The properties below therefore use integer
// coordinates in [-60, 60]: the largest intermediate in any of the queries is
// the circumcircle determinant, bounded by 3 * 120 * (2 * 120 * 2 * 120^2)
// < 2^40, so the double evaluation is exact and every sign the queries report
// must match an exact bigint evaluation of the same predicate. The references
// are the exact predicates of test/helpers/exact.ts (orient2, inCircle2) plus
// barycentric signs built from them, not transcriptions of the port, so a
// mis-transcribed term would show up.
// ---------------------------------------------------------------------------

// Negation that keeps 0 as +0, so toBe() (Object.is) does not see -0.
const negateSign = (s: number): number => (s === 0 ? 0 : -s);

// The exact predicates of test/helpers/exact.ts applied to integer-coordinate
// points: orient(a, b, c) is the sign of cross(b - a, c - a) and
// inCircle(a, b, c, d) is +1 when d is inside the circumcircle of the
// counterclockwise triangle <a,b,c>.
const orient = (a: Vector, b: Vector, c: Vector): number =>
    orient2(BigInt(a.values[0]), BigInt(a.values[1]),
        BigInt(b.values[0]), BigInt(b.values[1]),
        BigInt(c.values[0]), BigInt(c.values[1]));

const inCircle = (a: Vector, b: Vector, c: Vector, d: Vector): number =>
    inCircle2(BigInt(a.values[0]), BigInt(a.values[1]),
        BigInt(b.values[0]), BigInt(b.values[1]),
        BigInt(c.values[0]), BigInt(c.values[1]),
        BigInt(d.values[0]), BigInt(d.values[1]));

// Integer point coordinates small enough that every product below is exact.
const ipoint = fc.tuple(fc.integer({ min: -60, max: 60 }),
    fc.integer({ min: -60, max: 60 })).map(([x, y]) => v2(x, y));

// A generator that produces collinear and coincident configurations often:
// each point is either free or built on the lattice line through two anchors.
const degenerateTriple = fc.tuple(ipoint, ipoint,
    fc.integer({ min: -3, max: 4 }), fc.integer({ min: -2, max: 2 }),
    fc.boolean())
    .map(([A, B, num, perp, useLattice]) => {
        if (!useLattice) { return [A, B, v2(num * 7, perp * 5)] as Vector[]; }
        // P = A + num*(B - A) + perp*Perp(B - A)/1, all integer.
        const dx = B.values[0] - A.values[0];
        const dy = B.values[1] - A.values[1];
        const P = v2(A.values[0] + num * dx - perp * dy,
            A.values[1] + num * dy + perp * dx);
        return [A, B, P] as Vector[];
    });

describe('PrimalQuery2 verification', () => {
    it('toLine returns the exact sign of the cross-product determinant', () => {
        check(fc.tuple(ipoint, ipoint, ipoint), ([P, V0, V1]) => {
            const query = new PrimalQuery2(3, [P, V0, V1]);
            // det = cross(P - V0, V1 - V0) = -cross(V1 - V0, P - V0).
            const expected = negateSign(orient(V0, V1, P));
            expect(query.toLine(P, 1, 2)).toBe(expected);
            // The index overload uses mVertices[i] as the test point.
            expect(query.toLine(0, 1, 2)).toBe(expected);
            // Reversing the line direction negates the sign.
            expect(query.toLine(0, 2, 1)).toBe(negateSign(expected));
        });
    });

    it('toLine is exact on collinear and coincident configurations', () => {
        check(degenerateTriple, ([A, B, P]) => {
            const query = new PrimalQuery2(3, [P, A, B]);
            expect(query.toLine(0, 1, 2)).toBe(negateSign(orient(A, B, P)));
        });
    });

    it('toLineWithOrder agrees with toLine and with the upstream order rule', () => {
        check(degenerateTriple, ([A, B, P]) => {
            const query = new PrimalQuery2(3, [P, A, B]);
            const { sign: s, order } = query.toLineWithOrder(0, 1, 2);
            expect(s).toBe(query.toLine(0, 1, 2));
            const det = negateSign(orient(A, B, P));
            if (det !== 0) {
                expect(order).toBe(det > 0 ? +3 : -3);
                return;
            }
            // Collinear. Upstream compares dot = Dot(P-V0, V1-V0) against
            // |P-V0|^2 instead of |V1-V0|^2 (issue #100); reproduce that rule
            // exactly, including the swap it induces.
            const x0 = BigInt(P.values[0] - A.values[0]);
            const y0 = BigInt(P.values[1] - A.values[1]);
            const x1 = BigInt(B.values[0] - A.values[0]);
            const y1 = BigInt(B.values[1] - A.values[1]);
            const d = x0 * x1 + y0 * y1;
            const sqrLength = x0 * x0 + y0 * y0;
            const expected = d === 0n ? -1
                : d < 0n ? -2
                    : d === sqrLength ? +1
                        : d > sqrLength ? +2 : 0;
            expect(order).toBe(expected);
        });
    });

    it('the upstream order swap: interior points report +2, beyond-V1 report 0', () => {
        // Regression pin for issue #100 on exact input: P = V0 + c*(V1-V0).
        const V0 = v2(0, 0), V1 = v2(4, 2);
        const query = new PrimalQuery2(2, [V0, V1]);
        // c = 1/2 (interior): documented order 0, upstream reports +2.
        expect(query.toLineWithOrder(v2(2, 1), 0, 1).order).toBe(+2);
        // c = 2 (beyond V1): documented order +2, upstream reports 0.
        expect(query.toLineWithOrder(v2(8, 4), 0, 1).order).toBe(0);
        // The unaffected cases.
        expect(query.toLineWithOrder(v2(0, 0), 0, 1).order).toBe(-1);
        expect(query.toLineWithOrder(v2(4, 2), 0, 1).order).toBe(+1);
        expect(query.toLineWithOrder(v2(-4, -2), 0, 1).order).toBe(-2);
        // toLineExtended performs the same test correctly.
        expect(query.toLineExtended(v2(2, 1), V0, V1))
            .toBe(PrimalQuery2OrderType.COLLINEAR_CONTAIN);
        expect(query.toLineExtended(v2(8, 4), V0, V1))
            .toBe(PrimalQuery2OrderType.COLLINEAR_RIGHT);
    });

    it('toTriangle matches the exact barycentric sign classification', () => {
        check(fc.tuple(ipoint, ipoint, ipoint, ipoint), ([P, A, B, C]) => {
            if (orient(A, B, C) <= 0) { return; }   // the query needs a CCW triangle
            const query = new PrimalQuery2(4, [P, A, B, C]);
            // Exact barycentric signs of P with respect to <A,B,C>.
            const b0 = orient(P, B, C);
            const b1 = orient(P, C, A);
            const b2 = orient(P, A, B);
            const expected = (b0 < 0 || b1 < 0 || b2 < 0) ? +1
                : (b0 > 0 && b1 > 0 && b2 > 0) ? -1 : 0;
            expect(query.toTriangle(0, 1, 2, 3)).toBe(expected);
            expect(query.toTriangle(P, 1, 2, 3)).toBe(expected);
            // The classification does not depend on which vertex starts the
            // counterclockwise cycle.
            expect(query.toTriangle(0, 2, 3, 1)).toBe(expected);
            expect(query.toTriangle(0, 3, 1, 2)).toBe(expected);
        });
    });

    it('toCircumcircle matches the exact in-circle determinant', () => {
        check(fc.tuple(ipoint, ipoint, ipoint, ipoint), ([P, A, B, C]) => {
            if (orient(A, B, C) <= 0) { return; }
            const query = new PrimalQuery2(4, [P, A, B, C]);
            // inCircle2 is +1 when P is strictly inside the circumcircle of
            // the counterclockwise triangle; the query returns -1 for inside.
            const expected = negateSign(inCircle(A, B, C, P));
            expect(query.toCircumcircle(0, 1, 2, 3)).toBe(expected);
            expect(query.toCircumcircle(P, 1, 2, 3)).toBe(expected);
            // Cyclic permutations keep the orientation and the answer.
            expect(query.toCircumcircle(0, 2, 3, 1)).toBe(expected);
            // A transposition reverses the orientation and negates the answer.
            expect(query.toCircumcircle(0, 2, 1, 3)).toBe(negateSign(expected));
        });
    });

    it('points on a lattice circumcircle are reported as on the circle', () => {
        // (0,0), (10,0), (0,10) has circumcenter (5,5) and radius^2 = 50; the
        // lattice points at distance^2 = 50 from (5,5) all lie on the circle.
        const A = v2(0, 0), B = v2(10, 0), C = v2(0, 10);
        const query = new PrimalQuery2(3, [A, B, C]);
        for (const [dx, dy] of [[5, 5], [-5, 5], [5, -5], [-5, -5],
            [1, 7], [7, 1], [-1, 7], [7, -1], [1, -7], [-7, 1], [-1, -7], [-7, -1]]) {
            expect(query.toCircumcircle(v2(5 + dx, 5 + dy), 0, 1, 2)).toBe(0);
        }
        expect(query.toCircumcircle(v2(5, 5), 0, 1, 2)).toBe(-1);
        expect(query.toCircumcircle(v2(20, 20), 0, 1, 2)).toBe(+1);
    });

    it('toLineExtended matches its documented classification exactly', () => {
        check(degenerateTriple, ([Q0, Q1, P]) => {
            const query = new PrimalQuery2();
            const equal = (a: Vector, b: Vector): boolean =>
                a.values[0] === b.values[0] && a.values[1] === b.values[1];
            let expected: PrimalQuery2OrderType;
            if (equal(Q0, Q1)) {
                expected = PrimalQuery2OrderType.Q0_EQUALS_Q1;
            } else if (equal(P, Q0)) {
                expected = PrimalQuery2OrderType.P_EQUALS_Q0;
            } else if (equal(P, Q1)) {
                expected = PrimalQuery2OrderType.P_EQUALS_Q1;
            } else {
                const det = orient(Q0, Q1, P);
                if (det > 0) {
                    expected = PrimalQuery2OrderType.POSITIVE;
                } else if (det < 0) {
                    expected = PrimalQuery2OrderType.NEGATIVE;
                } else {
                    const dx = BigInt(Q1.values[0] - Q0.values[0]);
                    const dy = BigInt(Q1.values[1] - Q0.values[1]);
                    const d = dx * BigInt(P.values[0] - Q0.values[0])
                        + dy * BigInt(P.values[1] - Q0.values[1]);
                    const sqrLength = dx * dx + dy * dy;
                    expected = d < 0n ? PrimalQuery2OrderType.COLLINEAR_LEFT
                        : d > sqrLength ? PrimalQuery2OrderType.COLLINEAR_RIGHT
                            : PrimalQuery2OrderType.COLLINEAR_CONTAIN;
                }
            }
            expect(query.toLineExtended(P, Q0, Q1)).toBe(expected);
        });
    });

    it('toLineExtended checks the equality cases in upstream order', () => {
        const query = new PrimalQuery2();
        const A = v2(3, 4);
        // Q0 == Q1 wins over P == Q0.
        expect(query.toLineExtended(A, A, A)).toBe(PrimalQuery2OrderType.Q0_EQUALS_Q1);
        // P == Q0 is tested before P == Q1.
        expect(query.toLineExtended(A, A, v2(0, 0)))
            .toBe(PrimalQuery2OrderType.P_EQUALS_Q0);
        expect(query.toLineExtended(A, v2(0, 0), A))
            .toBe(PrimalQuery2OrderType.P_EQUALS_Q1);
    });

    it('set() and the vertex array are held by reference, not copied', () => {
        check(fc.tuple(ipoint, ipoint, ipoint), ([P, V0, V1]) => {
            const vertices = [P, V0, V1];
            const query = new PrimalQuery2(3, vertices);
            expect(query.getVertices()).toBe(vertices);
            expect(query.getNumVertices()).toBe(3);
            const other = [V1, P, V0];
            query.set(3, other);
            expect(query.getVertices()).toBe(other);
        });
    });
});
