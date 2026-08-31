import { describe, it, expect } from 'vitest';
import { PrimalQuery2, PrimalQuery2OrderType } from '../src/PrimalQuery2';
import { Vector } from '../src/Vector';

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
