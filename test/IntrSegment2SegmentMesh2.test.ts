import { describe, it, expect } from 'vitest';
import { IntrLine2SegmentMesh2FI } from '../src/IntrLine2SegmentMesh2.js';
import {
    IntrSegment2SegmentMesh2FI,
    defaultIntrSegment2SegmentMesh2FIResult
} from '../src/IntrSegment2SegmentMesh2.js';
import { Line } from '../src/Line.js';
import { Segment } from '../src/Segment.js';
import { SegmentMesh } from '../src/SegmentMesh.js';
import { Vector, add, dot, length, mul, sub } from '../src/Vector.js';
import { IntrRay2SegmentMesh2FI } from '../src/IntrRay2SegmentMesh2.js';
import { Ray } from '../src/Ray.js';
import {
    check, expectClose, expectVectorClose, fc, rotationFrame,
    segment as arbSegment, unitVector, wellScaledVector
} from './helpers/arbitraries.js';

function vec(x: number, y: number): Vector {
    return Vector.fromArray([x, y]);
}

function seg(p0: number[], p1: number[]): Segment {
    return Segment.fromEndpoints(Vector.fromArray(p0), Vector.fromArray(p1));
}

// A unit square as a closed contiguous mesh.
function squareMesh(): SegmentMesh {
    return SegmentMesh.fromContiguous(
        [vec(0, 0), vec(1, 0), vec(1, 1), vec(0, 1)], false);
}

const fi = new IntrSegment2SegmentMesh2FI();

describe('IntrSegment2SegmentMesh2', () => {
    it('has an empty default result', () => {
        expect(defaultIntrSegment2SegmentMesh2FIResult().intersections)
            .toHaveLength(0);
    });

    it('finds both crossings of a square and reports [0,1] parameters', () => {
        // The segment spans x in [-1, 2] at y = 0.5, so the crossings at
        // x = 0 and x = 1 are at parameters 1/3 and 2/3.
        const result = fi.find(seg([-1, 0.5], [2, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(1 / 3, 12);
        expect(result.intersections[1].segmentParameter).toBeCloseTo(2 / 3, 12);
        expect(result.intersections[0].point.values[0]).toBeCloseTo(0, 12);
        expect(result.intersections[1].point.values[0]).toBeCloseTo(1, 12);
    });

    it('discards crossings beyond the segment endpoints', () => {
        // The segment stops inside the square, so only the x = 0 edge is hit.
        const result = fi.find(seg([-1, 0.5], [0.5, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(1);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(2 / 3, 12);
        // A segment fully inside the square hits nothing.
        expect(fi.find(seg([0.2, 0.5], [0.8, 0.5]), squareMesh()).intersections)
            .toHaveLength(0);
        // A segment fully outside on the line of the square hits nothing.
        expect(fi.find(seg([-3, 0.5], [-2, 0.5]), squareMesh()).intersections)
            .toHaveLength(0);
    });

    it('keeps crossings exactly at the segment endpoints', () => {
        const result = fi.find(seg([0, 0.5], [1, 0.5]), squareMesh());
        expect(result.intersections).toHaveLength(2);
        expect(result.intersections[0].segmentParameter).toBeCloseTo(0, 12);
        expect(result.intersections[1].segmentParameter).toBeCloseTo(1, 12);
    });

    it('reports nothing for a segment that misses the mesh', () => {
        expect(fi.find(seg([-1, 5], [2, 5]), squareMesh()).intersections)
            .toHaveLength(0);
    });

    it('matches the line query filtered by parameter on random inputs', () => {
        let state = 909090;
        const rand = () => {
            state = (1103515245 * state + 12345) % 2147483648;
            return state / 2147483648 * 2 - 1;
        };

        const lsQuery = new IntrLine2SegmentMesh2FI();
        const vertices: Vector[] = [];
        for (let i = 0; i < 10; ++i) {
            vertices.push(vec(i - 5, (i % 2 === 0 ? -1 : 1) * 1.5));
        }
        const mesh = SegmentMesh.fromContiguous(vertices, true);

        let numHits = 0;
        for (let trial = 0; trial < 300; ++trial) {
            const p0 = [rand() * 6, rand() * 3];
            const p1 = [rand() * 6, rand() * 3];
            if (p0[0] === p1[0] && p0[1] === p1[1]) {
                continue;
            }
            const S = seg(p0, p1);
            const result = fi.find(S, mesh);

            const direction = sub(S.p[1], S.p[0]);
            const lineResult = lsQuery.find(
                Line.fromOriginDirection(S.p[0], direction), mesh);
            const expected = lineResult.intersections.filter(
                (object) => object.lineParameter >= 0
                    && object.lineParameter <= 1);
            expect(result.intersections).toHaveLength(expected.length);
            numHits += result.intersections.length;

            for (let i = 0; i < expected.length; ++i) {
                const actual = result.intersections[i];
                expect(actual.indexPair).toEqual(expected[i].indexPair);
                expect(actual.segmentParameter)
                    .toBeCloseTo(expected[i].lineParameter, 12);
                expect(actual.segmentParameter).toBeGreaterThanOrEqual(0);
                expect(actual.segmentParameter).toBeLessThanOrEqual(1);
                // The point is on the segment and on the named mesh segment.
                const onSegment = add(S.p[0],
                    mul(actual.segmentParameter, direction));
                expect(length(sub(actual.point, onSegment)))
                    .toBeCloseTo(0, 10);
                const V0 = mesh.getVertices()[actual.indexPair[0]];
                const V1 = mesh.getVertices()[actual.indexPair[1]];
                const onMeshSegment = add(
                    mul(1 - actual.meshSegmentParameter, V0),
                    mul(actual.meshSegmentParameter, V1));
                expect(length(sub(actual.point, onMeshSegment)))
                    .toBeCloseTo(0, 10);
            }
        }
        expect(numHits).toBeGreaterThan(50);
    });
});

// ---------------------------------------------------------------------------
// Verification group V34: property-based re-verification against
// GTE/Mathematics/IntrSegment2SegmentMesh2.h at commit d29e7758ae26.
// ---------------------------------------------------------------------------

describe('IntrSegment2SegmentMesh2 verification', () => {
    const fiQ = new IntrSegment2SegmentMesh2FI();
    const lsQ = new IntrLine2SegmentMesh2FI();

    // The line-segment tests are knife-edge when the query line passes
    // through a mesh vertex: which of the two incident mesh segments
    // reports the crossing then depends on the rounding of the line
    // origin and direction, so a reversed or moved copy of the same line
    // can report a different number of hits. Properties that compare two
    // such copies skip those configurations.
    function clearsMeshVertices(origin: Vector, direction: Vector,
        mesh: SegmentMesh, tolerance = 1e-6): boolean {
        const dd = dot(direction, direction);
        if (dd === 0) { return false; }
        for (const v of mesh.getVertices()) {
            const w = sub(v, origin);
            const perp = sub(w, mul(dot(w, direction) / dd, direction));
            if (Math.sqrt(dot(perp, perp)) < tolerance) { return false; }
        }
        return true;
    }

    const arbMesh = fc.tuple(fc.integer({ min: 3, max: 9 }),
        fc.array(fc.double({ min: 0.5, max: 2, noNaN: true }),
            { minLength: 9, maxLength: 9 }))
        .map(([n, radii]) => {
            const vs: Vector[] = [];
            for (let i = 0; i < n; ++i) {
                const a = (2 * Math.PI * i) / n;
                const r = radii[i];
                vs.push(Vector.fromArray([r * Math.cos(a), r * Math.sin(a)]));
            }
            return SegmentMesh.fromContiguous(vs, false);
        });

    it('is the line-mesh query filtered to parameters in [0,1]', () => {
        check(fc.tuple(arbSegment(2), arbMesh), ([S, M]) => {
            // The query builds the line from p[0] with the NON-unit direction
            // p[1] - p[0], so the reported parameter is in [0,1].
            const line = Line.fromOriginDirection(S.p[0],
                sub(S.p[1], S.p[0]));
            const ls = lsQ.find(line, M);
            const kept = ls.intersections.filter(
                o => o.lineParameter >= 0 && o.lineParameter <= 1);
            const r = fiQ.find(S, M);
            expect(r.intersections.length).toBe(kept.length);
            for (let i = 0; i < kept.length; ++i) {
                expect(r.intersections[i].segmentParameter)
                    .toBe(kept[i].lineParameter);
                expect(r.intersections[i].meshSegmentParameter)
                    .toBe(kept[i].meshSegmentParameter);
                expect(r.intersections[i].indexPair)
                    .toEqual(kept[i].indexPair);
            }
        });
    });

    it('the reported points are on the segment and on the mesh segment', () => {
        check(fc.tuple(arbSegment(2), arbMesh), ([S, M]) => {
            const vertices = M.getVertices();
            const r = fiQ.find(S, M);
            for (const o of r.intersections) {
                const s = o.segmentParameter;
                expect(s).toBeGreaterThanOrEqual(0);
                expect(s).toBeLessThanOrEqual(1);
                expect(o.meshSegmentParameter).toBeGreaterThanOrEqual(0);
                expect(o.meshSegmentParameter).toBeLessThanOrEqual(1);
                const scale = 1 + length(S.p[0]) + length(S.p[1]);
                // The documented convention: (1-s)*p[0] + s*p[1].
                expectVectorClose(o.point,
                    add(mul(1 - s, S.p[0]), mul(s, S.p[1])),
                    1e-9 * scale, 1e-9);
                const V0 = vertices[o.indexPair[0]];
                const V1 = vertices[o.indexPair[1]];
                const t = o.meshSegmentParameter;
                expectVectorClose(o.point,
                    add(mul(1 - t, V0), mul(t, V1)), 1e-9 * scale, 1e-9);
            }
        });
    });

    it('reports the intersections in increasing segment parameter', () => {
        check(fc.tuple(arbSegment(2), arbMesh), ([S, M]) => {
            const r = fiQ.find(S, M);
            for (let i = 1; i < r.intersections.length; ++i) {
                expect(r.intersections[i].segmentParameter)
                    .toBeGreaterThanOrEqual(
                        r.intersections[i - 1].segmentParameter);
            }
        });
    });

    it('reverses consistently: the parameters mirror about 1/2', () => {
        // The reported parameter is for (1-s)*p[0] + s*p[1], so reversing the
        // segment must give the same points with s -> 1-s and the order
        // reversed. Hits within rounding of an endpoint or of a mesh vertex
        // are skipped, where the inclusive filters can flip.
        check(fc.tuple(arbSegment(2), arbMesh), ([S, M]) => {
            if (!clearsMeshVertices(S.p[0], sub(S.p[1], S.p[0]), M)) {
                return;
            }
            const r0 = fiQ.find(S, M);
            for (const o of r0.intersections) {
                if (o.segmentParameter < 1e-6 || o.segmentParameter > 1 - 1e-6) {
                    return;
                }
            }
            const S2 = Segment.fromEndpoints(S.p[1], S.p[0]);
            const r1 = fiQ.find(S2, M);
            expect(r1.intersections.length).toBe(r0.intersections.length);
            const n = r0.intersections.length;
            const scale = 1 + length(S.p[0]) + length(S.p[1]);
            for (let i = 0; i < n; ++i) {
                const a = r0.intersections[i];
                const b = r1.intersections[n - 1 - i];
                expectClose(b.segmentParameter, 1 - a.segmentParameter,
                    1e-9, 1e-9);
                expectVectorClose(b.point, a.point, 1e-9 * scale, 1e-9);
                expect(b.indexPair).toEqual(a.indexPair);
            }
        });
    });

    it('a segment with both endpoints inside a closed mesh crosses it twice', () => {
        check(fc.tuple(arbMesh, unitVector(2)), ([M, d]) => {
            // A long segment through the (star-shaped) polygon centre.
            const S = Segment.fromEndpoints(mul(-20, d), mul(20, d));
            const r = fiQ.find(S, M);
            for (const o of r.intersections) {
                if (o.meshSegmentParameter < 1e-9
                    || o.meshSegmentParameter > 1 - 1e-9) {
                    return;
                }
            }
            expect(r.intersections.length).toBe(2);
        });
    });

    it('does not alias the input mesh or the line-query result', () => {
        const M = SegmentMesh.fromContiguous([
            Vector.fromArray([-1, -1]), Vector.fromArray([1, -1]),
            Vector.fromArray([1, 1]), Vector.fromArray([-1, 1])], false);
        const S = Segment.fromEndpoints(Vector.fromArray([-5, 0]),
            Vector.fromArray([5, 0]));
        const r = fiQ.find(S, M);
        expect(r.intersections.length).toBe(2);
        const before = M.getVertices().map(v => v.values.slice());
        r.intersections[0].point.values[0] = 999;
        r.intersections[0].indexPair[1] = 7;
        expect(M.getVertices().map(v => v.values.slice())).toEqual(before);
        const r2 = fiQ.find(S, M);
        expect(r2.intersections[0].point.values[0]).toBeCloseTo(-1, 12);
    });

    it('is equivariant under a rigid motion', () => {
        check(fc.tuple(arbSegment(2), arbMesh, rotationFrame(2),
            wellScaledVector(2)),
            ([S, M, Rot, t]) => {
                const rot = (v: Vector) => Vector.fromArray([
                    dot(Rot[0], v), dot(Rot[1], v)]);
                const map = (v: Vector) => add(rot(v), t);
                if (!clearsMeshVertices(S.p[0], sub(S.p[1], S.p[0]),
                    M)) {
                    return;
                }
                const r0 = fiQ.find(S, M);
                for (const o of r0.intersections) {
                    if (o.segmentParameter < 1e-6
                        || o.segmentParameter > 1 - 1e-6) {
                        return;
                    }
                }
                const M2 = SegmentMesh.fromContiguous(
                    M.getVertices().map(map), false);
                const S2 = Segment.fromEndpoints(map(S.p[0]), map(S.p[1]));
                const r1 = fiQ.find(S2, M2);
                expect(r1.intersections.length).toBe(r0.intersections.length);
                for (let i = 0; i < r0.intersections.length; ++i) {
                    expectClose(r1.intersections[i].segmentParameter,
                        r0.intersections[i].segmentParameter, 1e-7, 1e-7);
                }
            });
    });

    it('reports nothing for a segment far from the mesh', () => {
        check(fc.tuple(arbMesh, unitVector(2)), ([M, d]) => {
            const S = Segment.fromEndpoints(mul(20, d), mul(30, d));
            expect(fiQ.find(S, M).intersections.length).toBe(0);
        });
    });

    it('has an empty default result that is not shared', () => {
        const a = defaultIntrSegment2SegmentMesh2FIResult();
        const b = defaultIntrSegment2SegmentMesh2FIResult();
        expect(a.intersections).toEqual([]);
        expect(a.intersections).not.toBe(b.intersections);
    });

    it('documents the degenerate (point) segment behaviour', () => {
        // The query builds the line direction p[1] - p[0], which is the zero
        // vector for a degenerate segment. Upstream keeps whatever the
        // line-mesh query returns (its remove_if predicate is false for NaN);
        // the port's inclusive filter '0 <= t <= 1' drops NaN parameters. The
        // two agree here because the line-segment query reports no
        // intersection for a zero-length line direction.
        const M = SegmentMesh.fromContiguous([
            Vector.fromArray([-1, -1]), Vector.fromArray([1, -1]),
            Vector.fromArray([1, 1]), Vector.fromArray([-1, 1])], false);
        const onEdge = Vector.fromArray([0, -1]);
        const S = Segment.fromEndpoints(onEdge, onEdge.clone());
        const r = fiQ.find(S, M);
        for (const o of r.intersections) {
            expect(Number.isFinite(o.segmentParameter)).toBe(true);
            expect(Number.isFinite(o.meshSegmentParameter)).toBe(true);
        }
    });
});
