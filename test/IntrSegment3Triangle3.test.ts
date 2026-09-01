import { describe, it, expect } from 'vitest';
import {
    IntrSegment3Triangle3TI,
    IntrSegment3Triangle3FI
} from '../src/IntrSegment3Triangle3';
import {
    IntrLine3Triangle3FI
} from '../src/IntrLine3Triangle3';
import { Line } from '../src/Line';
import { Segment } from '../src/Segment';
import { Triangle } from '../src/Triangle';
import { Vector, add, dot, mul, sub } from '../src/Vector';
import { cross } from '../src/Vector3';

const ti = new IntrSegment3Triangle3TI();
const fi = new IntrSegment3Triangle3FI();
const lineFI = new IntrLine3Triangle3FI();

function v3(x: number, y: number, z: number): Vector {
    return Vector.fromArray([x, y, z]);
}

function makeRandom(seed: number): () => number {
    let state = seed >>> 0;
    return () => {
        state = (Math.imul(state, 1664525) + 1013904223) >>> 0;
        return state / 4294967296;
    };
}

const unitTriangle = Triangle.fromVertices(
    v3(0, 0, 0), v3(1, 0, 0), v3(0, 1, 0));

// An independent oracle: intersect the segment with the triangle plane and
// compute the barycentric coordinates from signed sub-triangle areas. The
// parameter u is in [0,1] along p0 -> p1.
function oracle(segment: Segment, triangle: Triangle):
    { hit: boolean, u: number, bary: number[] } | null {
    const e1 = sub(triangle.v[1], triangle.v[0]);
    const e2 = sub(triangle.v[2], triangle.v[0]);
    const n = cross(e1, e2);
    const dir = sub(segment.p[1], segment.p[0]);
    const den = dot(dir, n);
    if (den === 0) {
        return null;
    }
    const u = dot(sub(triangle.v[0], segment.p[0]), n) / den;
    const p = add(segment.p[0], mul(u, dir));
    const nn = dot(n, n);
    const b1 = dot(cross(sub(p, triangle.v[0]), e2), n) / nn;
    const b2 = dot(cross(e1, sub(p, triangle.v[0])), n) / nn;
    const b0 = 1 - b1 - b2;
    return {
        hit: u >= 0 && u <= 1 && b0 >= 0 && b1 >= 0 && b2 >= 0,
        u,
        bary: [b0, b1, b2]
    };
}

describe('IntrSegment3Triangle3TI', () => {
    it('reports a segment crossing the triangle interior', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -1), v3(0.25, 0.25, 1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(true);
    });

    it('rejects a segment that stops short of the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -3), v3(0.25, 0.25, -1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('accepts a segment whose endpoint touches the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -2), v3(0.25, 0.25, 0));
        expect(ti.test(segment, unitTriangle).intersect).toBe(true);
    });

    it('rejects a segment that starts past the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, 1), v3(0.25, 0.25, 3));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('rejects a segment crossing the plane outside the triangle', () => {
        const segment = Segment.fromEndpoints(
            v3(0.75, 0.75, -1), v3(0.75, 0.75, 1));
        expect(ti.test(segment, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection when the segment is parallel to the plane', () => {
        const coplanar = Segment.fromEndpoints(v3(-1, 0.25, 0), v3(1, 0.25, 0));
        expect(ti.test(coplanar, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate (zero-length) segment', () => {
        // The centered direction is the zero vector, so Dot(D,N) is zero and
        // the query takes the parallel branch.
        const degenerate = Segment.fromEndpoints(
            v3(0.25, 0.25, 0), v3(0.25, 0.25, 0));
        expect(ti.test(degenerate, unitTriangle).intersect).toBe(false);
    });

    it('reports no intersection for a degenerate triangle', () => {
        const degenerate = Triangle.fromVertices(
            v3(0, 0, 0), v3(1, 0, 0), v3(2, 0, 0));
        const segment = Segment.fromEndpoints(
            v3(0.5, 0, -1), v3(0.5, 0, 1));
        expect(ti.test(segment, degenerate).intersect).toBe(false);
    });
});

describe('IntrSegment3Triangle3FI', () => {
    it('computes the centered parameter, barycentrics and point', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -1), v3(0.25, 0.25, 3));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(true);
        // The centered form has center (0.25,0.25,1), extent 2 and direction
        // (0,0,1), so the intersection is at s = -1.
        expect(result.parameter).toBeCloseTo(-1, 12);
        expect(result.triangleBary[0]).toBeCloseTo(0.5, 12);
        expect(result.triangleBary[1]).toBeCloseTo(0.25, 12);
        expect(result.triangleBary[2]).toBeCloseTo(0.25, 12);
        expect(result.point.values[0]).toBeCloseTo(0.25, 12);
        expect(result.point.values[1]).toBeCloseTo(0.25, 12);
        expect(result.point.values[2]).toBeCloseTo(0, 12);
    });

    it('reports the endpoint parameter for a touching endpoint', () => {
        const segment = Segment.fromEndpoints(
            v3(0.25, 0.25, -2), v3(0.25, 0.25, 0));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(true);
        const { extent } = segment.getCenteredForm();
        expect(result.parameter).toBeCloseTo(extent, 12);
    });

    it('leaves the default result when there is no intersection', () => {
        const segment = Segment.fromEndpoints(
            v3(5, 5, -1), v3(5, 5, 1));
        const result = fi.find(segment, unitTriangle);
        expect(result.intersect).toBe(false);
        expect(result.parameter).toBe(0);
        expect(result.triangleBary).toEqual([0, 0, 0]);
        expect(result.point.values).toEqual([0, 0, 0]);
    });

    it('is watertight along an edge shared by two triangles', () => {
        const t0 = unitTriangle;
        const t1 = Triangle.fromVertices(
            v3(1, 0, 0), v3(1, 1, 0), v3(0, 1, 0));
        for (let i = 0; i <= 10; ++i) {
            const s = i / 10;
            const p = add(mul(1 - s, v3(1, 0, 0)), mul(s, v3(0, 1, 0)));
            const segment = Segment.fromEndpoints(
                add(p, v3(0, 0, -1)), add(p, v3(0, 0, 1)));
            expect(fi.find(segment, t0).intersect
                || fi.find(segment, t1).intersect).toBe(true);
        }
    });
});

describe('IntrSegment3Triangle3 consistency', () => {
    it('agrees between TI, FI, the line query and the oracle', () => {
        const rand = makeRandom(97531);
        const triangle = Triangle.fromVertices(
            v3(0.5, -1, 0.25), v3(2, 0.5, -0.5), v3(-0.75, 1.5, 1));
        let hits = 0;
        for (let trial = 0; trial < 500; ++trial) {
            const p0 = v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5);
            const p1 = v3(3 * rand() - 1.5, 3 * rand() - 1.5, 3 * rand() - 1.5);
            const segment = Segment.fromEndpoints(p0, p1);
            const { center, direction, extent } = segment.getCenteredForm();
            if (extent < 1e-6) {
                continue;
            }

            const tiResult = ti.test(segment, triangle);
            const fiResult = fi.find(segment, triangle);
            expect(tiResult.intersect).toBe(fiResult.intersect);

            const expected = oracle(segment, triangle);
            expect(expected).not.toBeNull();
            expect(tiResult.intersect).toBe(expected!.hit);

            // The segment hits exactly when the line hits with |s| <= extent.
            const line = Line.fromOriginDirection(center, direction);
            const lineResult = lineFI.find(line, triangle);
            expect(tiResult.intersect).toBe(lineResult.intersect
                && Math.abs(lineResult.parameter) <= extent);

            if (fiResult.intersect) {
                ++hits;
                expect(Math.abs(fiResult.parameter))
                    .toBeLessThanOrEqual(extent);
                for (let i = 0; i < 3; ++i) {
                    expect(fiResult.triangleBary[i])
                        .toBeCloseTo(expected!.bary[i], 8);
                }
                const combo = add(add(
                    mul(fiResult.triangleBary[0], triangle.v[0]),
                    mul(fiResult.triangleBary[1], triangle.v[1])),
                    mul(fiResult.triangleBary[2], triangle.v[2]));
                for (let i = 0; i < 3; ++i) {
                    expect(combo.values[i])
                        .toBeCloseTo(fiResult.point.values[i], 8);
                }
                // The centered parameter maps to the [0,1] parameter of the
                // endpoint form.
                const u = (fiResult.parameter + extent) / (2 * extent);
                expect(u).toBeCloseTo(expected!.u, 8);
            }
        }
        expect(hits).toBeGreaterThan(10);
    });
});
