import { describe, expect, it } from 'vitest';
import { DistPointTriangle } from '../src/DistPointTriangle';
import { DistSegment2Triangle2 } from '../src/DistSegment2Triangle2';
import { Segment } from '../src/Segment';
import { Triangle } from '../src/Triangle';
import { Vector, add, mul, sub } from '../src/Vector';

function v(...values: number[]): Vector {
    return Vector.fromArray(values);
}

function segment(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(v(...p0), v(...p1));
}

function triangle(v0: number[], v1: number[], v2: number[]): Triangle {
    return Triangle.fromVertices(v(...v0), v(...v1), v(...v2));
}

const ptQuery = new DistPointTriangle();

function pointTriangleSqrDistance(p: Vector, t: Triangle): number {
    return ptQuery.compute(p, t).sqrDistance;
}

// The squared distance from segment(t) to the solid triangle is a convex
// function of t on [0,1], so a ternary search finds its minimum.
function bruteForceSqrDistance(s: Segment, t: Triangle): number {
    const f = (u: number) => pointTriangleSqrDistance(
        add(s.p[0], mul(u, sub(s.p[1], s.p[0]))), t);

    let lo = 0;
    let hi = 1;
    for (let i = 0; i < 200; ++i) {
        const m0 = lo + (hi - lo) / 3;
        const m1 = hi - (hi - lo) / 3;
        if (f(m0) < f(m1)) {
            hi = m1;
        }
        else {
            lo = m0;
        }
    }

    // Also sample densely, guarding against a nonconvexity in the search.
    let best = f(0.5 * (lo + hi));
    for (let i = 0; i <= 2000; ++i) {
        const value = f(i / 2000);
        if (value < best) {
            best = value;
        }
    }
    return best;
}

// Verify the internal consistency of a result: the closest points lie on
// their primitives and realize the reported distance.
function expectConsistent(result: ReturnType<DistSegment2Triangle2['compute']>,
    s: Segment, t: Triangle): void {
    expect(result.parameter).toBeGreaterThanOrEqual(-1e-12);
    expect(result.parameter).toBeLessThanOrEqual(1 + 1e-12);

    const onSegment = add(s.p[0],
        mul(result.parameter, sub(s.p[1], s.p[0])));
    for (let i = 0; i < 2; ++i) {
        expect(result.closest[0].values[i]).toBeCloseTo(onSegment.values[i], 6);
    }

    const b = result.barycentric;
    expect(b[0] + b[1] + b[2]).toBeCloseTo(1, 8);
    for (let i = 0; i < 3; ++i) {
        expect(b[i]).toBeGreaterThanOrEqual(-1e-12);
    }
    const onTriangle = add(add(mul(b[0], t.v[0]), mul(b[1], t.v[1])),
        mul(b[2], t.v[2]));
    for (let i = 0; i < 2; ++i) {
        expect(result.closest[1].values[i]).toBeCloseTo(onTriangle.values[i], 6);
    }

    const diff = sub(result.closest[0], result.closest[1]);
    const length = Math.sqrt(diff.values[0] * diff.values[0]
        + diff.values[1] * diff.values[1]);
    expect(length).toBeCloseTo(result.distance, 8);
    expect(result.sqrDistance).toBeCloseTo(result.distance * result.distance, 8);
}

// A small deterministic linear congruential generator.
function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (state * 1664525 + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

describe('DistSegment2Triangle2', () => {
    const query = new DistSegment2Triangle2();
    const unitTriangle = triangle([0, 0], [1, 0], [0, 1]);

    it('reports zero distance for a segment crossing the triangle', () => {
        const result = query.compute(segment([-1, 0.25], [2, 0.25]),
            unitTriangle);
        expect(result.distance).toBe(0);
        expect(result.sqrDistance).toBe(0);
        expectConsistent(result, segment([-1, 0.25], [2, 0.25]),
            unitTriangle);
    });

    it('reports zero distance for a segment strictly inside the triangle',
        () => {
            const s = segment([0.2, 0.2], [0.3, 0.3]);
            const result = query.compute(s, unitTriangle);
            expect(result.distance).toBeCloseTo(0, 12);
            expectConsistent(result, s, unitTriangle);
        });

    it('measures a segment whose interior point is closest', () => {
        // The segment y = 2 runs parallel to the edge <V[0],V[1]> (y = 0).
        const s = segment([-1, 2], [1, 2]);
        const result = query.compute(s, unitTriangle);
        expect(result.distance).toBeCloseTo(1, 12);
        expectConsistent(result, s, unitTriangle);
    });

    it('clamps to the first endpoint when the line minimum is behind it',
        () => {
            // The closest line point is at t < 0, so the closest segment
            // point is p[0] = (3,0).
            const s = segment([3, 0], [5, 0]);
            const result = query.compute(s, unitTriangle);
            expect(result.parameter).toBe(0);
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
            expect(result.closest[1].values[0]).toBeCloseTo(1, 12);
            expectConsistent(result, s, unitTriangle);
        });

    it('clamps to the second endpoint when the line minimum is beyond it',
        () => {
            const s = segment([5, 0], [3, 0]);
            const result = query.compute(s, unitTriangle);
            expect(result.parameter).toBe(1);
            expect(result.distance).toBeCloseTo(2, 12);
            expect(result.closest[0].values[0]).toBeCloseTo(3, 12);
            expectConsistent(result, s, unitTriangle);
        });

    it('gives the same distance for both segment orientations', () => {
        const cases: Array<[number[], number[]]> = [
            [[3, 0], [5, 0]],
            [[-2, -2], [2, 3]],
            [[0.5, 0.5], [4, 4]],
            [[-1, 0.25], [2, 0.25]],
            [[-3, 5], [7, 5]]
        ];
        for (const [p0, p1] of cases) {
            const forward = query.compute(segment(p0, p1), unitTriangle);
            const backward = query.compute(segment(p1, p0), unitTriangle);
            expect(forward.distance).toBeCloseTo(backward.distance, 10);
        }
    });

    it('reports zero distance when an endpoint touches a vertex', () => {
        const s = segment([1, 0], [3, 2]);
        const result = query.compute(s, unitTriangle);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitTriangle);
    });

    it('handles a degenerate zero-length segment outside the triangle', () => {
        const s = segment([3, 4], [3, 4]);
        const result = query.compute(s, unitTriangle);
        const expected = ptQuery.compute(v(3, 4), unitTriangle);
        expect(result.distance).toBeCloseTo(expected.distance, 12);
        expect(result.parameter).toBe(0);
        expect(result.closest[0].values[0]).toBe(3);
        expect(result.closest[0].values[1]).toBe(4);
        expectConsistent(result, s, unitTriangle);
    });

    it('handles a degenerate zero-length segment inside the triangle', () => {
        const s = segment([0.25, 0.25], [0.25, 0.25]);
        const result = query.compute(s, unitTriangle);
        expect(result.distance).toBeCloseTo(0, 12);
        expectConsistent(result, s, unitTriangle);
    });

    it('handles a degenerate (single point) triangle', () => {
        const t = triangle([2, 2], [2, 2], [2, 2]);
        const s = segment([0, 0], [0, 4]);
        const result = query.compute(s, t);
        expect(result.distance).toBeCloseTo(2, 10);
    });

    it('matches a brute-force minimization on analytic configurations', () => {
        const cases: Array<[number[], number[]]> = [
            [[-4, 3], [4, 3]],
            [[2, 2], [3, 5]],
            [[-1, -1], [-1, 4]],
            [[0.1, 0.1], [0.2, 0.6]],
            [[-5, 0.5], [-0.5, 0.5]],
            [[1, 1], [2, 2]],
            [[0, -3], [0, -1]]
        ];
        for (const [p0, p1] of cases) {
            const s = segment(p0, p1);
            const result = query.compute(s, unitTriangle);
            const expected =
                Math.sqrt(bruteForceSqrDistance(s, unitTriangle));
            expect(result.distance).toBeCloseTo(expected, 6);
            expectConsistent(result, s, unitTriangle);
        }
    });

    it('matches a brute-force minimization on random configurations', () => {
        const random = makeRandom(20250510);
        const coord = () => 6 * random() - 3;
        for (let trial = 0; trial < 400; ++trial) {
            const t = triangle([coord(), coord()], [coord(), coord()],
                [coord(), coord()]);
            const s = segment([coord(), coord()], [coord(), coord()]);
            const result = query.compute(s, t);
            const expected = Math.sqrt(bruteForceSqrDistance(s, t));
            expect(result.distance).toBeCloseTo(expected, 6);
            expectConsistent(result, s, t);

            // Reversing the segment must not change the distance.
            const reversed = query.compute(
                Segment.fromEndpoints(s.p[1], s.p[0]), t);
            expect(reversed.distance).toBeCloseTo(result.distance, 6);
        }
    });
});
